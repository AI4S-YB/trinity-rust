//! 模型引导组装（M5）: reads 招募到模型转录本的基因分组 → 每组定向
//! 迷你重组装（复用 [`crate::recursive::run_mini`]）→ 合并去冗余
//! （[`crate::recursive::harvest_mini`]）。
//!
//! 三个阶段:
//! 1. [`cluster_transcripts_to_genes`]: minimap2 self-alignment（asm5）→
//!    identity ≥80% 且覆盖 ≥50% 的转录本对建边 → 并查集弱连通分量 = "基因";
//! 2. [`recruit_reads`]: minimap2 `-x sr` reads vs 模型 → 每 read best hit
//!    （identity × 覆盖 最高）→ 转录本 → 基因; PE 对（read 名去 /1 /2 后缀）
//!    取并集——任一端命中即整对入组, 两端不同组时取 R1 的命中;
//! 3. [`assemble_groups`]: 每组（+ orphan 兜底组）`run_mini` → 拼合 →
//!    harvest_mini（min_len 过滤 + 全序列去重）。
//!
//! minimap2 以外部进程调用（写临时 fasta 到 `std::env::temp_dir()` →
//! Command → 解析 PAF stdout）。

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use trinity_common::error::CommonError;

use crate::recursive::{harvest_mini, run_mini, MiniParams};

/// minimap2 可执行文件（环境变量 `TRINITY_MINIMAP2` 可覆盖）。
pub(crate) fn minimap2_bin() -> &'static str {
    static DEFAULT: &str = "/public/home/senior007/miniconda3/envs/cactus/bin/minimap2";
    match std::env::var("TRINITY_MINIMAP2") {
        Ok(v) if !v.is_empty() => Box::leak(v.into_boxed_str()),
        _ => DEFAULT,
    }
}

/// 逐调用唯一临时目录。
fn temp_dir(tag: &str) -> PathBuf {
    static SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let n = SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let d = std::env::temp_dir().join(format!("model_guided_{tag}_{}_{n}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    d
}

/// 运行 minimap2（query vs target）→ PAF 全文。`secondary_ok` = 保留次级
/// 比对（self-alignment 聚类时 T0↔T1 交叉命中是自身命中的 secondary,
/// 必须保留; reads 招募时只要 best hit, 关掉省输出）。
fn run_minimap2(
    query: &Path,
    target: &Path,
    preset: &str,
    secondary_ok: bool,
) -> Result<String, CommonError> {
    let mut cmd = std::process::Command::new(minimap2_bin());
    cmd.args(["-x", preset, "-t", "8"]);
    if !secondary_ok {
        cmd.arg("--secondary=no");
    }
    let out = cmd
        .args([&target.display().to_string(), &query.display().to_string()])
        .stderr(std::process::Stdio::null())
        .output()
        .map_err(|e| {
            CommonError::Parse(format!("cannot spawn minimap2 ({}): {e}", minimap2_bin()))
        })?;
    if !out.status.success() {
        return Err(CommonError::Parse(format!(
            "minimap2 failed with status {}",
            out.status
        )));
    }
    Ok(String::from_utf8_lossy(&out.stdout).into_owned())
}

/// PAF 记录（只取本模块需要的列）。
struct PafHit {
    qname: String,
    /// query 覆盖率（0-1）。
    qcov: f64,
    tname: String,
    /// target 覆盖率（0-1）。
    tcov: f64,
    /// identity（nmatch/alnlen, 0-1）。
    identity: f64,
}

/// PAF: qname qlen qstart qend strand tname tlen tstart tend nmatch alnlen mapq
fn parse_paf(text: &str) -> Vec<PafHit> {
    let mut out = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        let (Ok(qlen), Ok(qs), Ok(qe)) = (
            f[1].parse::<f64>(),
            f[2].parse::<f64>(),
            f[3].parse::<f64>(),
        ) else {
            continue;
        };
        let (Ok(tlen), Ok(ts), Ok(te)) = (
            f[6].parse::<f64>(),
            f[7].parse::<f64>(),
            f[8].parse::<f64>(),
        ) else {
            continue;
        };
        let (Ok(nmatch), Ok(alen)) = (f[9].parse::<f64>(), f[10].parse::<f64>()) else {
            continue;
        };
        if qlen <= 0.0 || tlen <= 0.0 || alen <= 0.0 {
            continue;
        }
        out.push(PafHit {
            qname: f[0].to_string(),
            qcov: (qe - qs) / qlen,
            tname: f[5].to_string(),
            tcov: (te - ts) / tlen,
            identity: nmatch / alen,
        });
    }
    out
}

/// 简单 FASTA 记录（accession 首空白前 token, 序列已折行合并）。
fn parse_fasta(data: &[u8]) -> Result<Vec<(String, String)>, CommonError> {
    let text = std::str::from_utf8(data)
        .map_err(|e| CommonError::FastaFormat(format!("non-UTF8 fasta: {e}")))?;
    let mut out = Vec::new();
    let mut header: Option<String> = None;
    let mut seq = String::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if let Some(h0) = header.take() {
                out.push((h0, std::mem::take(&mut seq)));
            }
            header = Some(h.split_whitespace().next().unwrap_or_default().to_string());
        } else if header.is_some() {
            seq.push_str(line.trim());
        }
    }
    if let Some(h0) = header {
        out.push((h0, seq));
    }
    Ok(out)
}

/// 原始 fasta 记录字节（header 行 + 序列行照抄）, 供组内重组装复用。
struct RawRecord {
    /// 原始 accession（含 /1 /2 后缀, PAF qname 匹配用）。
    acc: String,
    /// 配对 key（accession 去 /1 /2 后缀）。
    key: String,
    bytes: Vec<u8>,
}

fn parse_raw_records(data: &[u8]) -> Result<Vec<RawRecord>, CommonError> {
    let text = std::str::from_utf8(data)
        .map_err(|e| CommonError::FastaFormat(format!("non-UTF8 fasta: {e}")))?;
    let mut out = Vec::new();
    let mut cur: Vec<&str> = Vec::new();
    for line in text.lines() {
        if line.starts_with('>') && !cur.is_empty() {
            out.push(take_record(&cur)?);
            cur.clear();
        }
        cur.push(line);
    }
    if !cur.is_empty() {
        out.push(take_record(&cur)?);
    }
    Ok(out)
}

fn take_record(lines: &[&str]) -> Result<RawRecord, CommonError> {
    let acc = lines[0]
        .strip_prefix('>')
        .unwrap_or(lines[0])
        .split_whitespace()
        .next()
        .unwrap_or_default();
    let key = acc
        .strip_suffix("/1")
        .or_else(|| acc.strip_suffix("/2"))
        .unwrap_or(acc)
        .to_string();
    let mut bytes = String::new();
    for l in lines {
        bytes.push_str(l);
        bytes.push('\n');
    }
    Ok(RawRecord {
        acc: acc.to_string(),
        key,
        bytes: bytes.into_bytes(),
    })
}

// ---------------------------------------------------------------------------
// 并查集
// ---------------------------------------------------------------------------

struct Dsu {
    parent: Vec<usize>,
}

impl Dsu {
    fn new(n: usize) -> Self {
        Dsu {
            parent: (0..n).collect(),
        }
    }
    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            let r = self.find(self.parent[x]);
            self.parent[x] = r;
        }
        self.parent[x]
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[rb] = ra;
        }
    }
    /// 弱连通分量（root 升序, 组内索引升序）。
    fn groups(&mut self) -> Vec<Vec<usize>> {
        let mut map: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..self.parent.len() {
            map.entry(self.find(i)).or_default().push(i);
        }
        let mut out: Vec<Vec<usize>> = map.into_values().collect();
        for g in &mut out {
            g.sort_unstable();
        }
        out.sort_by_key(|g| g[0]);
        out
    }
}

// ---------------------------------------------------------------------------
// 阶段 1: 基因聚类
// ---------------------------------------------------------------------------

/// 模型转录本按 ≥50% 覆盖（identity ≥80%）聚类成"基因"（弱连通分量）。
/// 输入: 模型 fasta 字节; 输出: `Vec<Vec<usize>>`（每组是 transcript 索引列表）。
pub fn cluster_transcripts_to_genes(model_fa: &[u8]) -> Result<Vec<Vec<usize>>, CommonError> {
    let records = parse_fasta(model_fa)?;
    let n = records.len();
    if n == 0 {
        return Err(CommonError::FastaFormat(
            "model fasta has no records".to_string(),
        ));
    }
    let mut dsu = Dsu::new(n);
    if n == 1 {
        return Ok(dsu.groups());
    }
    let tmp = temp_dir("cluster");
    fs::create_dir_all(&tmp)?;
    let model_path = tmp.join("model.fa");
    fs::write(&model_path, model_fa)?;
    let paf = run_minimap2(&model_path, &model_path, "asm5", true)?;
    let idx: HashMap<&str, usize> = records
        .iter()
        .map(|(h, _)| h.as_str())
        .enumerate()
        .map(|(i, h)| (h, i))
        .collect();
    for h in parse_paf(&paf) {
        if h.qname == h.tname {
            continue;
        }
        let (Some(&q), Some(&t)) = (idx.get(h.qname.as_str()), idx.get(h.tname.as_str())) else {
            continue;
        };
        // 覆盖取 query/target 较大者（self-alignment 双向对称, 保守取并）。
        if h.identity >= 0.8 && h.qcov.max(h.tcov) >= 0.5 {
            dsu.union(q, t);
        }
    }
    let _ = fs::remove_dir_all(&tmp);
    Ok(dsu.groups())
}

// ---------------------------------------------------------------------------
// 阶段 2: reads 招募
// ---------------------------------------------------------------------------

/// reads（">name/1" PE 或单端格式 fasta 字节）招募到基因组。
/// 返回 `(组 reads fasta 字节, orphan reads fasta 字节)`; 组数与
/// `gene_groups` 一一对应（无 reads 的组为空字节）。
pub fn recruit_reads(
    reads_fa: &[u8],
    model_fa: &[u8],
    gene_groups: &[Vec<usize>],
) -> Result<(Vec<Vec<u8>>, Vec<u8>), CommonError> {
    let records = parse_raw_records(reads_fa)?;
    let transcripts = parse_fasta(model_fa)?;
    let tmp = temp_dir("recruit");
    fs::create_dir_all(&tmp)?;
    let reads_path = tmp.join("reads.fa");
    let model_path = tmp.join("model.fa");
    fs::write(&reads_path, reads_fa)?;
    fs::write(&model_path, model_fa)?;
    let paf = run_minimap2(&reads_path, &model_path, "sr", false)?;
    let _ = fs::remove_dir_all(&tmp);

    // transcript index → gene group index。
    let mut trans_group: HashMap<String, usize> = HashMap::new();
    for (gi, group) in gene_groups.iter().enumerate() {
        for &ti in group {
            if let Some((h, _)) = transcripts.get(ti) {
                trans_group.insert(h.clone(), gi);
            }
        }
    }

    // best hit per 原始 read 名（score = identity × qcov）。
    // minimap2 `-x sr` 会把 /1 /2 配对后缀剥掉再回填（qname 可能变成
    // "name/1/1"）——查询时剥一次尾缀还原原始 accession。
    let mut best: HashMap<String, (f64, usize)> = HashMap::new();
    for h in parse_paf(&paf) {
        let Some(&gi) = trans_group.get(&h.tname) else {
            continue;
        };
        let q = h
            .qname
            .strip_suffix("/1")
            .or_else(|| h.qname.strip_suffix("/2"))
            .unwrap_or(&h.qname)
            .to_string();
        let score = h.identity * h.qcov;
        // 同时以剥后缀名与原始 qname 为键（minimap2 对配对后缀的处理
        // 因输入排布而异: 有的场景回填成 "x/2/2", 有的保持 "x/1" 原样）。
        for key in [q, h.qname] {
            let e = best.entry(key).or_insert((f64::MIN, gi));
            if score > e.0 {
                *e = (score, gi);
            }
        }
    }

    // PE 并集: key（去 /1 /2）的组 = R1 命中优先, 否则 R2 命中。
    let mut key_group: HashMap<&str, usize> = HashMap::new();
    let mut key_r1: HashMap<&str, usize> = HashMap::new();
    for r in &records {
        if let Some(&(_, gi)) = best.get(&r.acc) {
            if r.acc.ends_with("/1") {
                key_r1.insert(r.key.as_str(), gi);
            }
            key_group.entry(r.key.as_str()).or_insert(gi);
        }
    }
    // R1 命中覆盖（两端不同组时取 R1 的命中）。
    for (k, gi) in key_r1 {
        key_group.insert(k, gi);
    }

    let mut groups: Vec<Vec<u8>> = vec![Vec::new(); gene_groups.len()];
    let mut orphan: Vec<u8> = Vec::new();
    let mut assigned = 0usize;
    for r in &records {
        match key_group.get(r.key.as_str()) {
            Some(&gi) => {
                groups[gi].extend_from_slice(&r.bytes);
                assigned += 1;
            }
            None => orphan.extend_from_slice(&r.bytes),
        }
    }
    eprintln!(
        "-- model-guided recruitment: {} reads records assigned, {} orphan",
        assigned,
        records.len() - assigned
    );
    Ok((groups, orphan))
}

// ---------------------------------------------------------------------------
// 阶段 3: 分组组装
// ---------------------------------------------------------------------------

/// 每组（+ orphan 兜底组）跑 mini-assembly → 合并 + 去冗余。
/// 返回合并后的 Trinity fasta 字节。组产物落 `outdir/g{i}.mini.fa`
/// （存在即跳过 = 断点）; `outdir` 为 None 时用临时目录（无断点）。
pub fn assemble_groups(
    groups: &[Vec<u8>],
    orphan: &[u8],
    params: &MiniParams,
    outdir: Option<&Path>,
    threads: usize,
) -> Result<Vec<u8>, CommonError> {
    let tmp_owned;
    let dir = match outdir {
        Some(d) => d.to_path_buf(),
        None => {
            tmp_owned = temp_dir("assemble");
            fs::create_dir_all(&tmp_owned)?;
            tmp_owned
        }
    };
    fs::create_dir_all(&dir)?;

    // jobs: (组名, reads 字节); 空组跳过。
    let mut jobs: Vec<(String, Vec<u8>)> = groups
        .iter()
        .enumerate()
        .filter(|(_, g)| !g.is_empty())
        .map(|(i, g)| (format!("g{i}"), g.clone()))
        .collect();
    if !orphan.is_empty() {
        jobs.push(("orphan".to_string(), orphan.to_vec()));
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .map_err(|e| CommonError::Parse(format!("cannot build model-guided pool: {e}")))?;
    let failed: std::sync::Mutex<Vec<String>> = std::sync::Mutex::new(Vec::new());
    let mini_files: Vec<PathBuf> = pool.install(|| {
        jobs.par_iter()
            .map(|(name, reads)| {
                let mini_fa = dir.join(format!("{name}.mini.fa"));
                if !mini_fa.exists() {
                    // 大栈线程（对齐 run_partition_pool: 递归会溢出默认 2MB 栈）。
                    let mut p = params.clone();
                    p.name = name.clone();
                    let reads = reads.clone();
                    let r = std::thread::Builder::new()
                        .stack_size(256 * 1024 * 1024)
                        .spawn(move || run_mini(&reads, &p))
                        .and_then(|h| h.join().map_err(|_| std::io::Error::other("join")))
                        .and_then(|r| {
                            r.map(|o| o.fasta)
                                .map_err(|e| std::io::Error::other(e.to_string()))
                        });
                    match r {
                        Ok(fasta) => fs::write(&mini_fa, &fasta)?,
                        Err(e) => {
                            eprintln!("-model-guided group {name} failed: {e}");
                            failed.lock().unwrap().push(format!("{name}: {e}"));
                        }
                    }
                }
                Ok(mini_fa)
            })
            .collect::<Result<Vec<_>, CommonError>>()
    })?;

    let failed = failed.into_inner().unwrap();
    if !failed.is_empty() {
        return Err(CommonError::Parse(format!(
            "Error, model-guided mini pipeline failed for {} group(s) (first: {})",
            failed.len(),
            failed[0]
        )));
    }

    let final_path = dir.join("Trinity.tmp.fasta");
    let n = harvest_mini(&mini_files, params.min_contig_length, &final_path)?;
    if n == 0 {
        return Err(CommonError::Parse(
            "WARNING: no transcripts harvested from model-guided group output".to_string(),
        ));
    }
    let bytes = fs::read(&final_path)?;
    if outdir.is_none() {
        let _ = fs::remove_dir_all(&dir);
    }
    Ok(bytes)
}

// ---------------------------------------------------------------------------
// orchestrate 接线（M5 主线）
// ---------------------------------------------------------------------------

/// 归一化 + prep 之后的 model-guided 主线。
/// `target_fa` = both.fa/single.fa（reads 来源）; 返回最终 fasta 路径。
pub fn run_model_guided(
    args: &crate::args::TrinityArgs,
    outdir: &Path,
    target_fa: &Path,
) -> Result<PathBuf, CommonError> {
    let mg_ok = outdir.join(".model_guided.ok");
    let final_fa = PathBuf::from(format!("{}.Trinity.fasta", outdir.display()));
    if crate::checkpoint::checkpoint_exists(&mg_ok) && final_fa.exists() {
        eprintln!(
            "---- Trinity (rust) checkpoint found, skipping: {}",
            mg_ok.display()
        );
        return Ok(final_fa);
    }

    // ---- 1. 模型获取: --model 外部文件; 否则自引导（first_pass/ 全流程）----
    let model_bytes: Vec<u8> = match &args.model {
        Some(p) => fs::read(p)?,
        None => {
            let fp_dir = outdir.join("first_pass");
            let fp_fa = PathBuf::from(format!("{}.Trinity.fasta", fp_dir.display()));
            let fp_ok = outdir.join(".model_first_pass.ok");
            crate::checkpoint::run_with_checkpoint(
                &fp_ok,
                "Model first pass (self-guided)",
                || {
                    let mut a = args.clone();
                    a.model_guided = false;
                    a.recursive = false;
                    a.super_reads = false;
                    a.path_engine = crate::args::PathEngine::PogFlow;
                    a.output = fp_dir.clone();
                    if fp_fa.exists() {
                        return Ok(());
                    }
                    let p = crate::orchestrate::run_trinity(&a)?;
                    if p != fp_fa {
                        return Err(CommonError::Parse(format!(
                            "first pass output path mismatch: {} != {}",
                            p.display(),
                            fp_fa.display()
                        )));
                    }
                    Ok(())
                },
            )?;
            fs::read(&fp_fa)?
        }
    };

    // ---- 2. 基因聚类 + reads 招募 ----
    let gene_groups = cluster_transcripts_to_genes(&model_bytes)?;
    eprintln!(
        "-- model-guided: {} model transcripts → {} gene groups",
        parse_fasta(&model_bytes)?.len(),
        gene_groups.len()
    );
    let reads_bytes = fs::read(target_fa)?;
    let (groups, orphan) = recruit_reads(&reads_bytes, &model_bytes, &gene_groups)?;

    // ---- 3. 分组组装（断点: model_guided_groups/g{i}.mini.fa 存在即跳过）----
    let params = MiniParams {
        k: args.kmer_size,
        strand: args.ss_lib_type.is_some(),
        min_contig_length: args.min_contig_length,
        path_reinforcement_distance: 25,
        max_pair_distance: args.group_pairs_distance,
        min_kmer_count: args.min_kmer_count,
        stack_size_mb: args.bfly_stack_mb,
        name: String::new(),
    };
    let mg_dir = outdir.join("model_guided_groups");
    let merged = assemble_groups(&groups, &orphan, &params, Some(&mg_dir), args.cpu)?;

    // ---- 4. 输出（命名/gene_trans_map 与 classic 尾段一致）----
    let trinity_tmp = outdir.join("Trinity.tmp.fasta");
    fs::write(&trinity_tmp, &merged)?;
    fs::rename(&trinity_tmp, &final_fa)?;
    crate::orchestrate::write_gene_trans_map(
        &final_fa,
        &final_fa.with_file_name(format!(
            "{}.gene_trans_map",
            final_fa.file_name().unwrap_or_default().to_string_lossy()
        )),
    )?;

    // ---- 5. 收尾清理（与 recursive 相同集合）----
    if !args.no_cleanup {
        for p in [
            target_fa.to_path_buf(),
            outdir.join("left.fa"),
            outdir.join("right.fa"),
        ] {
            let _ = fs::remove_file(p);
        }
    }
    fs::write(&mg_ok, b"")?;
    Ok(final_fa)
}

// ---------------------------------------------------------------------------
// 测试（合成集手推, 计划 Task 1 Step 4）
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn xorshift_nt(x: &mut u64) -> char {
        *x ^= *x << 13;
        *x ^= *x >> 7;
        *x ^= *x << 17;
        match *x % 4 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        }
    }

    fn randseq(seed: u64, n: usize) -> String {
        let mut x = seed | 1;
        (0..n).map(|_| xorshift_nt(&mut x)).collect()
    }

    fn revcomp(s: &str) -> String {
        s.bytes()
            .rev()
            .map(|b| match b {
                b'A' => 'T',
                b'C' => 'G',
                b'G' => 'C',
                _ => 'A',
            })
            .collect()
    }

    /// 3 条转录本: T0/T1 共享 500nt 前缀（83% 覆盖, 100% identity）→ 同基因;
    /// T2 独立。
    fn synth_model() -> (String, String, String) {
        let common = randseq(0x243F6A8885A308D3, 500);
        let t0 = format!("{}{}", common, randseq(11, 100));
        let t1 = format!("{}{}", common, randseq(22, 100));
        let t2 = randseq(0x9E3779B97F4A7C15, 600);
        (t0, t1, t2)
    }

    /// 每条转录本平铺 n_pairs 对 PE reads（100bp, 步进 = 全长/n/2）。
    fn tiled_pe_reads(name_prefix: &str, seq: &str, n_pairs: usize) -> String {
        let b = seq.as_bytes();
        let step = (b.len() - 250) / n_pairs.max(1);
        let mut text = String::new();
        for i in 0..n_pairs {
            let s = i * step;
            text.push_str(&format!(">{}{}_1\n{}\n", name_prefix, i, &seq[s..s + 100]));
            text.push_str(&format!(
                ">{}{}_2\n{}\n",
                name_prefix,
                i,
                revcomp(std::str::from_utf8(&b[s + 150..s + 250]).unwrap())
            ));
        }
        text
    }

    fn model_fa(t0: &str, t1: &str, t2: &str) -> Vec<u8> {
        format!(">T0\n{t0}\n>T1\n{t1}\n>T2\n{t2}\n").into_bytes()
    }

    #[test]
    fn cluster_and_recruit_synthetic() {
        let (t0, t1, t2) = synth_model();
        let model = model_fa(&t0, &t1, &t2);
        let groups = cluster_transcripts_to_genes(&model).unwrap();
        assert_eq!(groups.len(), 2, "T0/T1 同基因, T2 独立 → 2 组: {groups:?}");
        assert!(groups.iter().any(|g| g == &vec![0, 1]));
        assert!(groups.iter().any(|g| g == &vec![2]));

        let reads = format!(
            "{}{}{}",
            tiled_pe_reads("a", &t0, 10),
            tiled_pe_reads("b", &t1, 10),
            tiled_pe_reads("c", &t2, 10)
        );
        let (recruited, orphan) = recruit_reads(reads.as_bytes(), &model, &groups).unwrap();
        assert_eq!(recruited.len(), 2);
        assert!(orphan.is_empty(), "全部 reads 都应命中模型");
        // a/b reads 同组（T0/T1 同基因）; c reads 另一组。
        let n_records = |b: &[u8]| b.iter().filter(|&&c| c == b'>').count();
        let (g01, g2) = if groups[0] == vec![0, 1] {
            (0, 1)
        } else {
            (1, 0)
        };
        assert_eq!(n_records(&recruited[g01]), 40, "a+b 各 10 对（每对 2 条）");
        assert_eq!(n_records(&recruited[g2]), 20, "c 10 对");
    }

    #[test]
    fn pe_union_uses_either_mate() {
        // R1 不命中（来自新序列）, R2 命中 T2 → 整对应入 T2 的组。
        let (t0, _t1, t2) = synth_model();
        // 模型不能有空记录——重建只含 T0/T2。
        let model = format!(">T0\n{t0}\n>T2\n{t2}\n").into_bytes();
        let groups = cluster_transcripts_to_genes(&model).unwrap();
        let novel = randseq(0xDEADBEEFCAFEBABE, 300);
        let reads = format!(">x/1\n{novel}\n>x/2\n{}\n", revcomp(&t2[..100]));
        let (recruited, orphan) = recruit_reads(reads.as_bytes(), &model, &groups).unwrap();
        let n_records = |b: &[u8]| b.iter().filter(|&&c| c == b'>').count();
        let total: usize = recruited.iter().map(|g| n_records(g)).sum();
        assert_eq!(total, 2, "R2 命中 → 整对入组（并集）");
        assert!(orphan.is_empty());
    }

    #[test]
    fn orphan_reads_assembled() {
        // 2 条独立模型转录本 + 一组 reads 来自第三条"新"转录本 →
        // orphan 组产出该 reads 的组装。
        let (t0, _t1, t2) = synth_model();
        let model = format!(">T0\n{t0}\n>T2\n{t2}\n").into_bytes();
        let groups_idx = cluster_transcripts_to_genes(&model).unwrap();
        let novel_ref = randseq(0x1234ABCD5678EF01, 600);
        let novel_reads = tiled_pe_reads("n", &novel_ref, 25);
        // 模型命中的 reads（给两组各一点, 保证非空组路径也被走到）。
        let m0 = tiled_pe_reads("m0", &t0, 5);
        let m2 = tiled_pe_reads("m2", &t2, 5);
        let reads = format!("{m0}{m2}{novel_reads}");
        let (recruited, orphan) = recruit_reads(reads.as_bytes(), &model, &groups_idx).unwrap();
        assert!(
            orphan.windows(2).any(|w| w == b">n"),
            "新转录本 reads 应落入 orphan: {}",
            String::from_utf8_lossy(&orphan[..orphan.len().min(200)])
        );
        assert!(recruited.iter().any(|g| g.windows(3).any(|w| w == b">m0")));
        assert!(recruited.iter().any(|g| g.windows(3).any(|w| w == b">m2")));

        let params = MiniParams {
            k: 25,
            strand: false,
            min_contig_length: 100,
            path_reinforcement_distance: 25,
            max_pair_distance: 500,
            min_kmer_count: 1,
            stack_size_mb: 64,
            name: String::new(),
        };
        let merged = assemble_groups(&recruited, &orphan, &params, None, 2).unwrap();
        let s = String::from_utf8(merged.clone()).unwrap();
        assert!(s.contains('>'), "合并产物非空: {} bytes", merged.len());
        // 最长记录 ≥400bp（新转录本 600nt 被 25 对 reads 覆盖, 参考 recursive 测试）。
        let longest = parse_fasta(&merged)
            .unwrap()
            .iter()
            .map(|(_, seq)| seq.len())
            .max()
            .unwrap_or(0);
        assert!(longest >= 400, "最长转录本 {longest}bp < 400");
    }
}
