//! M6a 蛋白引导组装（探路版）: DIAMOND blastx 模板选择 + 缺口轮廓 +
//! 引导延长 stage（插在 butterfly 池与 harvest 之间）。

use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

use trinity_common::error::CommonError;

/// DIAMOND 可执行文件（env `TRINITY_DIAMOND` 可覆盖）。
fn diamond_bin() -> &'static str {
    static DEFAULT: &str = "/public/home/senior007/miniconda3/envs/diamond2go/bin/diamond";
    match std::env::var("TRINITY_DIAMOND") {
        Ok(v) if !v.is_empty() => Box::leak(v.into_boxed_str()),
        _ => DEFAULT,
    }
}

pub struct DiamondHit {
    pub qid: String,
    pub sid: String,
    pub pident: f64,
    pub qlen: i64,
    pub slen: i64,
    pub qstart: i64,
    pub qend: i64,
    pub sstart: i64,
    pub send: i64,
    pub evalue: f64,
    pub bitscore: f64,
    pub qframe: i32,
}

/// outfmt 6: qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qframe
pub fn parse_diamond_tab(text: &str) -> Vec<DiamondHit> {
    let mut out = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 13 {
            continue;
        }
        let (Ok(qlen), Ok(slen)) = (f[4].parse::<i64>(), f[5].parse::<i64>()) else {
            continue;
        };
        let (Ok(qs), Ok(qe), Ok(ss), Ok(se)) = (
            f[6].parse::<i64>(),
            f[7].parse::<i64>(),
            f[8].parse::<i64>(),
            f[9].parse::<i64>(),
        ) else {
            continue;
        };
        out.push(DiamondHit {
            qid: f[0].to_string(),
            sid: f[1].to_string(),
            pident: f[2].parse().unwrap_or(0.0),
            qlen,
            slen,
            qstart: qs,
            qend: qe,
            sstart: ss,
            send: se,
            evalue: f[10].parse().unwrap_or(f64::MAX),
            bitscore: f[11].parse().unwrap_or(0.0),
            qframe: f[12].parse().unwrap_or(0),
        });
    }
    out
}

/// 每查询最佳命中（bitscore 最高; 平局取先出现者）。
pub fn best_hit_per_query(hits: &[DiamondHit]) -> HashMap<String, &DiamondHit> {
    let mut best: HashMap<String, &DiamondHit> = HashMap::new();
    for h in hits {
        match best.get(&h.qid) {
            Some(b) if b.bitscore >= h.bitscore => {}
            _ => {
                best.insert(h.qid.clone(), h);
            }
        }
    }
    best
}

/// makedb（存在即跳过）→ 返回 dmnd 路径。
pub fn diamond_makedb_if_needed(
    db_fa: &Path,
    workdir: &Path,
    threads: usize,
) -> Result<PathBuf, CommonError> {
    fs::create_dir_all(workdir)?;
    let dmnd = workdir.join("protein.dmnd");
    if dmnd.exists() {
        return Ok(dmnd);
    }
    let out = std::process::Command::new(diamond_bin())
        .arg("makedb")
        .arg("--in")
        .arg(db_fa)
        .arg("-d")
        .arg(&dmnd)
        .arg("--threads")
        .arg(threads.to_string())
        .output()
        .map_err(|e| {
            CommonError::Parse(format!("cannot spawn diamond ({}): {e}", diamond_bin()))
        })?;
    if !out.status.success() {
        return Err(CommonError::Parse(format!(
            "diamond makedb failed: {}",
            String::from_utf8_lossy(&out.stderr)
        )));
    }
    Ok(dmnd)
}

/// blastx（DNA 查询 vs 蛋白库）→ tabular 全文（临时输出文件读回）。
pub fn diamond_blastx(
    query_fa: &Path,
    dmnd: &Path,
    threads: usize,
    tag: &str,
) -> Result<String, CommonError> {
    let out_path = std::env::temp_dir().join(format!("diamond_{tag}_{}", std::process::id()));
    let out = std::process::Command::new(diamond_bin())
        .arg("blastx")
        .arg("-q")
        .arg(query_fa)
        .arg("-d")
        .arg(dmnd)
        .arg("--threads")
        .arg(threads.to_string())
        .arg("-e")
        .arg("1e-5")
        .arg("--max-target-seqs")
        .arg("5")
        .arg("--outfmt")
        .arg("6")
        .arg("qseqid")
        .arg("sseqid")
        .arg("pident")
        .arg("length")
        .arg("qlen")
        .arg("slen")
        .arg("qstart")
        .arg("qend")
        .arg("sstart")
        .arg("send")
        .arg("evalue")
        .arg("bitscore")
        .arg("qframe")
        .arg("-o")
        .arg(&out_path)
        .output()
        .map_err(|e| {
            CommonError::Parse(format!("cannot spawn diamond ({}): {e}", diamond_bin()))
        })?;
    if !out.status.success() {
        return Err(CommonError::Parse(format!(
            "diamond blastx failed: {}",
            String::from_utf8_lossy(&out.stderr)
        )));
    }
    let text = fs::read_to_string(&out_path).map_err(|e| {
        CommonError::Parse(format!(
            "cannot read diamond output {}: {e}",
            out_path.display()
        ))
    })?;
    let _ = fs::remove_file(&out_path);
    Ok(text)
}

/// junction 连接的最小 PE 支持对数（对齐管线 MIN_SEED_COVERAGE=2 惯例）。
const MIN_JUNCTION_PAIRS: usize = 2;

/// junction PE 支持门: 对每个候选延长, 取原转录本端侧 100nt (a) 与延长区
/// 端侧 100nt (b), minimap2 -x sr 比对 both.fa, 统计一端命中 a、另一端命中
/// b 的 read 对数（mate 名去 /1 /2 后配对）。≥ min_pairs 视为分子层面
/// 支持该连接（旁系游荡不会有跨越对）。
/// 返回每个 j{i} 的支持对数。
pub fn junction_pair_support(
    junction_fa: &Path, // 记录: j{i}_a / j{i}_b
    reads_fa: &Path,    // both.fa
    min_len: usize,     // 忽略过短命中 (<30)
    threads: usize,
) -> Result<HashMap<String, usize>, CommonError> {
    // minimap2 -x sr <junction_fa=target> <reads_fa=query> → PAF stdout;
    // PAF 列: qname=read 名, tname=junction 记录名。
    let out = std::process::Command::new(crate::model_guided::minimap2_bin())
        .args(["-x", "sr", "-t", &threads.to_string()])
        .arg(junction_fa)
        .arg(reads_fa)
        .stderr(std::process::Stdio::null())
        .output()
        .map_err(|e| {
            CommonError::Parse(format!(
                "cannot spawn minimap2 ({}): {e}",
                crate::model_guided::minimap2_bin()
            ))
        })?;
    if !out.status.success() {
        return Err(CommonError::Parse(format!(
            "minimap2 failed with status {}",
            out.status
        )));
    }
    Ok(junction_support_from_paf(
        &String::from_utf8_lossy(&out.stdout),
        min_len,
    ))
}

/// PAF 解析核心（纯函数, 便于单测）: qname = read（剥 /1 /2 配对）,
/// tname = `j{i}_a` / `j{i}_b`。同一 read 名的两条命中覆盖同一组的
/// a+b 则该组 +1（每对只计一次）。alnlen < min_len 的行忽略。
pub fn junction_support_from_paf(text: &str, min_len: usize) -> HashMap<String, usize> {
    // read 名(剥后) → 命中的 tname 集合
    let mut hits: HashMap<String, HashSet<String>> = HashMap::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 11 {
            continue;
        }
        let Ok(alnlen) = f[10].parse::<usize>() else {
            continue;
        };
        if alnlen < min_len {
            continue;
        }
        let qname = f[0]
            .strip_suffix("/1")
            .or_else(|| f[0].strip_suffix("/2"))
            .unwrap_or(f[0]);
        hits.entry(qname.to_string())
            .or_default()
            .insert(f[5].to_string());
    }
    // 组 i 的 a/b 名 → 组名
    let mut groups: HashMap<String, String> = HashMap::new();
    for names in hits.values() {
        for t in names {
            // t 形如 j{i}_a / j{i}_b; 逆序剥 _a/_b + `_`
            let Some((group, side)) = t.rsplit_once('_') else {
                continue;
            };
            if side != "a" && side != "b" {
                continue;
            }
            groups.insert(t.clone(), group.to_string());
        }
    }
    let mut out: HashMap<String, usize> = HashMap::new();
    for names in hits.values() {
        // 本对覆盖的 (组, 侧) 集合
        let mut sides: HashMap<&str, HashSet<&str>> = HashMap::new();
        for t in names {
            if let Some(g) = groups.get(t) {
                sides.entry(g).or_default().insert(&t[t.len() - 1..]);
            }
        }
        for (g, s) in sides {
            if s.contains("a") && s.contains("b") {
                *out.entry(g.to_string()).or_insert(0) += 1;
            }
        }
    }
    out
}

/// 蛋白库 fasta → sseqid → 氨基酸序列（header 取首空白 token）。
pub fn index_protein_db(path: &Path) -> Result<HashMap<String, Vec<u8>>, CommonError> {
    let text = fs::read_to_string(path).map_err(|e| {
        CommonError::Parse(format!("cannot read protein db {}: {e}", path.display()))
    })?;
    let mut out = HashMap::new();
    let mut cur: Option<(String, Vec<u8>)> = None;
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if let Some((k, v)) = cur.take() {
                out.insert(k, v);
            }
            cur = Some((
                h.split_whitespace().next().unwrap_or_default().to_string(),
                Vec::new(),
            ));
        } else if let Some((_, seq)) = cur.as_mut() {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if let Some((k, v)) = cur.take() {
        out.insert(k, v);
    }
    Ok(out)
}

/// 蛋白库净化: 丢弃零长度序列记录（DIAMOND makedb 拒绝空序列）。
/// 返回 (保留数, 丢弃数)。dst 已存在则直接读取计数返回（断点复用）。
pub fn sanitize_protein_db(src: &Path, dst: &Path) -> Result<(usize, usize), CommonError> {
    if let Some(parent) = dst.parent() {
        fs::create_dir_all(parent)?;
    }
    let text = fs::read_to_string(src).map_err(|e| {
        CommonError::Parse(format!("cannot read protein db {}: {e}", src.display()))
    })?;
    if dst.exists() {
        let (_, kept, _) = sanitize_protein_db_text(
            &fs::read_to_string(dst)
                .map_err(|e| CommonError::Parse(format!("cannot read {}: {e}", dst.display())))?,
        );
        return Ok((kept, 0));
    }
    let (out, kept, dropped) = sanitize_protein_db_text(&text);
    let tmp = dst.with_extension("tmp");
    fs::write(&tmp, &out)?;
    fs::rename(&tmp, dst)?;
    Ok((kept, dropped))
}

/// 内部: 过滤文本, 返回 (新文本, 保留数, 丢弃数)。
fn sanitize_protein_db_text(text: &str) -> (String, usize, usize) {
    let mut out = String::new();
    let mut kept = 0usize;
    let mut dropped = 0usize;
    let mut header: Option<String> = None;
    let mut seq = String::new();
    let mut flush = |header: &Option<String>, seq: &mut String| {
        if let Some(h) = header {
            if seq.is_empty() {
                dropped += 1;
            } else {
                kept += 1;
                out.push_str(&format!(">{h}\n"));
                for chunk in seq.as_bytes().chunks(60) {
                    out.push_str(std::str::from_utf8(chunk).unwrap_or_default());
                    out.push('\n');
                }
            }
        }
        seq.clear();
    };
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            flush(&header, &mut seq);
            header = Some(h.to_string());
        } else if header.is_some() {
            seq.push_str(line.trim());
        }
    }
    flush(&header, &mut seq);
    (out, kept, dropped)
}

/// 缺口轮廓: 延长方向 + 预算（单 HSP 版）。
#[derive(Debug, Clone, Copy)]
pub struct GapProfile {
    pub plus: bool,
    /// true = 向转录本 3' 延长。
    pub dir_right: bool,
    n_missing: usize,
    c_missing: usize,
    budget_nt: usize,
}

impl GapProfile {
    pub fn missing_nterm(&self) -> bool {
        self.n_missing > 0
    }
    pub fn missing_cterm(&self) -> bool {
        self.c_missing > 0
    }
    pub fn budget_nt(&self) -> usize {
        self.budget_nt
    }
}

/// 单 HSP: 覆盖区间 [min(sstart,send), max(sstart,send)]。
/// 补缺失跨度大的一侧（平局取 C 侧）; plus: C→Right / N→Left; minus 反转。
/// budget = 该侧缺失 aa × 3 nt × 3/2（容错余量）, 上限 5000。
pub fn gap_profile(h: &DiamondHit) -> GapProfile {
    let (lo, hi) = (h.sstart.min(h.send), h.sstart.max(h.send));
    let n_missing = (lo - 1).max(0) as usize;
    let c_missing = (h.slen - hi).max(0) as usize;
    let plus = h.qframe > 0;
    let fill_c = c_missing >= n_missing;
    let (dir_right, missing) = if fill_c {
        (plus, c_missing)
    } else {
        (!plus, n_missing)
    };
    GapProfile {
        plus,
        dir_right,
        n_missing,
        c_missing,
        budget_nt: (missing * 3 * 3 / 2).min(5000),
    }
}

use trinity_butterfly::graph_io::build_new_graph_use_kmers;
use trinity_butterfly::protein_extend::{extend_transcript, AnchorHsp, Direction};

/// 延长后验收门槛: 新 bitscore 须超旧值 × (1 + 该系数)。
const ACCEPT_MARGIN: f64 = 0.05;

/// stage 入口: 插在 butterfly 池之后、harvest 之前。
/// 逐组件改写 `<base>.graph.allProbPaths.fasta`（原件备份 `.orig`, 一次性）。
pub fn run_protein_guided(
    listing: &[(u64, PathBuf)],
    args: &crate::args::TrinityArgs,
    outdir: &Path,
    reads_fa: &Path,
) -> Result<(), CommonError> {
    use crate::checkpoint::checkpoint_exists;
    let ok = outdir.join(".protein_guided.ok");
    if checkpoint_exists(&ok) {
        eprintln!(
            "---- Trinity (rust) checkpoint found, skipping: {}",
            ok.display()
        );
        return Ok(());
    }
    let Some(db_fa) = args.protein_db.clone() else {
        return Err(CommonError::Parse(
            "--protein-guided requires --protein-db".to_string(),
        ));
    };
    let work = outdir.join("protein_guided");
    let threads = args.cpu.max(1);

    // ---- 1. 收集全部初步转录本 → queries.fa ----
    let mut queries: Vec<(String, PathBuf, String, String)> = Vec::new(); // (qid, base, header, seq)
    let mut qfa = String::new();
    for (i, (_id, base)) in listing.iter().enumerate() {
        let fa = PathBuf::from(format!("{}.graph.allProbPaths.fasta", base.display()));
        let Ok(text) = fs::read_to_string(&fa) else {
            continue;
        };
        for (header, seq) in parse_fasta_pairs(&text) {
            let acc = header.split_whitespace().next().unwrap_or("t");
            let qid = format!("c{i}_{acc}");
            qfa.push_str(&format!(">{qid}\n{seq}\n"));
            queries.push((qid, base.clone(), header, seq));
        }
    }
    if queries.is_empty() {
        eprintln!("-- protein-guided: no transcripts found, skipping");
        fs::write(&ok, b"")?;
        return Ok(());
    }
    let clean_fa = work.join("protein.clean.fa");
    let (kept, dropped) = sanitize_protein_db(&db_fa, &clean_fa)?;
    if dropped > 0 {
        eprintln!(
            "-- protein-guided: dropped {dropped} empty record(s) from protein db ({kept} kept)"
        );
    }
    let dmnd = diamond_makedb_if_needed(&clean_fa, &work, threads)?;
    let qfa_path = work.join("queries.fa");
    fs::write(&qfa_path, &qfa)?;

    // ---- 2. blastx + 模板 ----
    let hits = parse_diamond_tab(&diamond_blastx(&qfa_path, &dmnd, threads, "q")?);
    let best = best_hit_per_query(&hits);
    let proteins = index_protein_db(&clean_fa)?;

    // ---- 3. 逐转录本延长（有缺口才动）----
    struct ExtOut {
        base: PathBuf,
        header: String,
        new_seq: String,
        /// 原始序列（junction 端侧取样用）。
        orig_seq: String,
        /// true = Right（3' 延长）; false = Left。
        dir_right: bool,
        old_sid: String,
        old_bitscore: f64,
    }
    let mut extended: Vec<ExtOut> = Vec::new();
    let mut with_t = 0usize;
    // 组件分组: base → 有最佳命中且预算非零的查询下标（图按组件只读/解析一次）。
    let mut by_base: HashMap<PathBuf, Vec<usize>> = HashMap::new();
    for (i, (qid, base, _header, _seq)) in queries.iter().enumerate() {
        let Some(h) = best.get(qid) else { continue };
        with_t += 1;
        if gap_profile(h).budget_nt() == 0 {
            continue;
        }
        by_base.entry(base.clone()).or_default().push(i);
    }
    for (base, idxs) in &by_base {
        let gout = PathBuf::from(format!("{}.graph.out", base.display()));
        let Ok(gtext) = fs::read_to_string(&gout) else {
            eprintln!("-protein-guided: cannot reload {}", gout.display());
            continue;
        };
        let Ok(br) = build_new_graph_use_kmers(&gtext) else {
            eprintln!("-protein-guided: cannot parse {}", gout.display());
            continue;
        };
        for &i in idxs {
            let (qid, _base, header, seq) = &queries[i];
            let Some(h) = best.get(qid) else { continue };
            let gp = gap_profile(h);
            let anchor = AnchorHsp {
                qstart: if h.qframe > 0 {
                    h.qstart.min(h.qend)
                } else {
                    h.qstart.max(h.qend)
                },
                sstart: h.sstart,
                plus: h.qframe > 0,
            };
            let dir = if gp.dir_right {
                Direction::Right
            } else {
                Direction::Left
            };
            let template = proteins.get(&h.sid).cloned().unwrap_or_default();
            // 顶点名长度取图解析器自带的 ctx.kmer_size（chrysalis graph.out 的
            // 顶点是 K−1-mer, 与 --KMER_SIZE 不同——用 args 值会锚定/校验双失败）。
            let vlen = if br.ctx.kmer_size > 0 {
                br.ctx.kmer_size
            } else {
                args.kmer_size
            };
            let ext = extend_transcript(
                &br.graph,
                vlen,
                seq.as_bytes(),
                dir,
                &anchor,
                &template,
                gp.budget_nt(),
            );
            if ext.is_empty() {
                continue;
            }
            let Ok(ext_s) = std::str::from_utf8(&ext) else {
                eprintln!("-protein-guided: non-UTF8 extension bytes for {qid}, skipping");
                continue;
            };
            let new_seq = match dir {
                Direction::Right => format!("{seq}{ext_s}"),
                Direction::Left => format!("{ext_s}{seq}"),
            };
            extended.push(ExtOut {
                base: base.clone(),
                header: header.clone(),
                new_seq,
                orig_seq: seq.clone(),
                dir_right: gp.dir_right,
                old_sid: h.sid.clone(),
                old_bitscore: h.bitscore,
            });
        }
    }

    // ---- 4. 验收: junction PE 支持门 + 延长产物重 blastx ----
    let mut jpass = 0usize;
    let accepted: Vec<&ExtOut> = if extended.is_empty() {
        Vec::new()
    } else {
        // junction.fa: j{i}_a = 原转录本端侧 100nt, j{i}_b = 延长区连接侧
        // 100nt（Right: 延长区头部; Left: 延长区尾部——前插序的尾最贴近连接）。
        // 任一侧 <30nt 无法配对 → 该条目不支持（计数 0, 门会拒绝）。
        let mut jfa = String::new();
        for (i, e) in extended.iter().enumerate() {
            let ext_len = e.new_seq.len() - e.orig_seq.len();
            let (a, b) = if e.dir_right {
                (
                    &e.orig_seq[e.orig_seq.len().saturating_sub(100)..],
                    &e.new_seq[e.orig_seq.len()..e.orig_seq.len() + ext_len.min(100)],
                )
            } else {
                (
                    &e.orig_seq[..e.orig_seq.len().min(100)],
                    &e.new_seq[ext_len.saturating_sub(100)..ext_len],
                )
            };
            if a.len() >= 30 && b.len() >= 30 {
                jfa.push_str(&format!(">j{i}_a\n{a}\n>j{i}_b\n{b}\n"));
            }
        }
        let jfa_path = work.join("junction.fa");
        fs::write(&jfa_path, &jfa)?;
        let support = junction_pair_support(&jfa_path, reads_fa, 30, threads)?;
        jpass = extended
            .iter()
            .enumerate()
            .filter(|(i, _e)| {
                support.get(&format!("j{i}")).copied().unwrap_or(0) >= MIN_JUNCTION_PAIRS
            })
            .count();
        let mut xfa = String::new();
        for (i, e) in extended.iter().enumerate() {
            xfa.push_str(&format!(">x{i}\n{}\n", e.new_seq));
        }
        let xfa_path = work.join("extended.fa");
        fs::write(&xfa_path, &xfa)?;
        let xhits = parse_diamond_tab(&diamond_blastx(&xfa_path, &dmnd, threads, "x")?);
        let xb = best_hit_per_query(&xhits);
        let accepted: Vec<&ExtOut> = extended
            .iter()
            .enumerate()
            .filter(|(i, e)| {
                support.get(&format!("j{i}")).copied().unwrap_or(0) >= MIN_JUNCTION_PAIRS
                    && accepted(
                        &e.old_sid,
                        e.old_bitscore,
                        xb.get(&format!("x{i}")).map(|v| &**v),
                    )
            })
            .map(|(_, e)| e)
            .collect();
        eprintln!(
            "-- protein-guided: junction gate: {jpass}/{} with >={MIN_JUNCTION_PAIRS} junction pairs",
            extended.len()
        );
        accepted
    };

    // ---- 5. 改写组件 fasta（原件备份一次, 临时文件 + 原子改名）----
    let mut by_comp: HashMap<PathBuf, Vec<&ExtOut>> = HashMap::new();
    for e in &accepted {
        by_comp.entry(e.base.clone()).or_default().push(e);
    }
    let mut total_nt = 0usize;
    for (base, exts) in &by_comp {
        let fa = PathBuf::from(format!("{}.graph.allProbPaths.fasta", base.display()));
        let backup = PathBuf::from(format!("{}.graph.allProbPaths.orig", base.display()));
        let text = fs::read_to_string(&fa)?;
        if !backup.exists() {
            fs::write(&backup, &text)?;
        }
        let pairs: Vec<(String, String)> = exts
            .iter()
            .map(|e| (e.header.clone(), e.new_seq.clone()))
            .collect();
        let (out, _n, delta) = rewrite_component(&text, &pairs);
        total_nt += delta;
        let tmp = PathBuf::from(format!("{}.pg_tmp", fa.display()));
        fs::write(&tmp, &out)?;
        fs::rename(&tmp, &fa)?;
    }
    eprintln!(
        "-- protein-guided: {with_t} transcripts with template, {} extended, {jpass} with >={MIN_JUNCTION_PAIRS} junction pairs, {} accepted, {} reverted, +{total_nt} nt",
        extended.len(),
        accepted.len(),
        extended.len() - accepted.len()
    );
    fs::write(&ok, b"")?;
    Ok(())
}

/// 用通过的延长替换匹配 header 的记录, 其余原样返回; 返回 (新文本, 替换数, 净增碱基)。
fn rewrite_component(text: &str, exts: &[(String, String)]) -> (String, usize, usize) {
    let mut out = String::new();
    let mut n = 0usize;
    let mut delta = 0usize;
    for (header, seq) in parse_fasta_pairs(text) {
        match exts.iter().find(|(h, _)| *h == header) {
            Some((_, new_seq)) => {
                n += 1;
                delta += new_seq.len().saturating_sub(seq.len());
                out.push_str(&format!(">{header}\n{new_seq}\n"));
            }
            None => out.push_str(&format!(">{header}\n{seq}\n")),
        }
    }
    (out, n, delta)
}

/// 延长验收: 新最佳命中存在、sid 与原模板一致, 且 bitscore 须超旧值 × (1 + ACCEPT_MARGIN)。
/// 换目标蛋白（旁系同源）视为未通过验收——延长必须忠于原模板。
fn accepted(old_sid: &str, old_bitscore: f64, new_hit: Option<&DiamondHit>) -> bool {
    match new_hit {
        Some(h) => h.sid == old_sid && h.bitscore > old_bitscore * (1.0 + ACCEPT_MARGIN),
        None => false,
    }
}

fn parse_fasta_pairs(text: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut header: Option<String> = None;
    let mut seq = String::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if let Some(h0) = header.take() {
                out.push((h0, std::mem::take(&mut seq)));
            }
            header = Some(h.to_string());
        } else if header.is_some() {
            seq.push_str(line.trim());
        }
    }
    if let Some(h0) = header {
        out.push((h0, seq));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn rewrite_replaces_matching_header_and_keeps_others() {
        let text = ">a x y\nAAAA\n>b z\nCCCCCC\n";
        let exts = vec![("a x y".to_string(), "AAAAAA".to_string())];
        let (out, n, delta) = rewrite_component(text, &exts);
        assert_eq!(n, 1);
        assert_eq!(delta, 2);
        assert_eq!(out, ">a x y\nAAAAAA\n>b z\nCCCCCC\n");
    }

    #[test]
    fn rewrite_no_match_returns_verbatim() {
        let text = ">a x y\nAAAA\n>b z\nCCCCCC\n";
        let exts = vec![("absent".to_string(), "TTTT".to_string())];
        let (out, n, delta) = rewrite_component(text, &exts);
        assert_eq!(n, 0);
        assert_eq!(delta, 0);
        assert_eq!(out, ">a x y\nAAAA\n>b z\nCCCCCC\n");
    }

    #[test]
    fn rewrite_duplicate_header_replaces_all_copies() {
        let text = ">a\nAAAA\n>b\nCCCCCC\n>a\nAAAA\n";
        let exts = vec![("a".to_string(), "TT".to_string())];
        let (out, n, delta) = rewrite_component(text, &exts);
        assert_eq!(n, 2, "重复 header 全部被替换（当前语义）");
        assert_eq!(delta, 0); // 4→2 两次, 净 0
        assert_eq!(out, ">a\nTT\n>b\nCCCCCC\n>a\nTT\n");
    }

    #[test]
    fn acceptance_margin() {
        let hit = |sid: &str, bitscore: f64| DiamondHit {
            qid: "x0".into(),
            sid: sid.into(),
            pident: 90.0,
            qlen: 300,
            slen: 200,
            qstart: 1,
            qend: 300,
            sstart: 1,
            send: 100,
            evalue: 1e-40,
            bitscore,
            qframe: 1,
        };
        assert!(accepted("P1", 100.0, Some(&hit("P1", 106.0))));
        assert!(!accepted("P1", 100.0, Some(&hit("P1", 104.0))));
        assert!(
            !accepted("P1", 100.0, Some(&hit("P1", 105.0))),
            "105 不严格大于 100×1.05"
        );
        assert!(
            !accepted("P1", 100.0, Some(&hit("P2", 9999.0))),
            "换目标蛋白（旁系同源）拒绝, 无论 bitscore 多高"
        );
        assert!(!accepted("P1", 100.0, None), "无新命中拒绝");
    }

    #[test]
    fn parse_tab_roundtrip() {
        let text = "t1\tP1\t90.1\t100\t300\t200\t1\t300\t1\t100\t1e-40\t185.2\t1\n\
                    t1\tP2\t50.0\t50\t300\t80\t10\t160\t5\t55\t1e-9\t45.0\t-1\n";
        let hits = parse_diamond_tab(text);
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].qid, "t1");
        assert_eq!(hits[0].sid, "P1");
        assert_eq!(hits[0].qlen, 300);
        assert_eq!(hits[0].slen, 200);
        assert_eq!(hits[0].qframe, 1);
        assert_eq!(hits[1].qframe, -1);
        assert!((hits[0].bitscore - 185.2).abs() < 1e-9);
    }

    /// DIAMOND `--outfmt 6` 的默认 12 列（无 qlen/slen/qframe）必须被拒绝，
    /// 否则字段错位产生垃圾比对坐标。diamond_blastx 因此显式指定 13 列。
    #[test]
    fn parse_tab_rejects_default_12_columns() {
        // 默认 12 列: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        let default12 = "t1\tP1\t90.1\t100\t3\t0\t1\t300\t1\t100\t1e-40\t185.2\n";
        assert!(parse_diamond_tab(default12).is_empty());
        // 13 列（diamond_blastx 的实际输出布局）正常解析。
        let full13 = "t1\tP1\t90.1\t100\t300\t200\t1\t300\t1\t100\t1e-40\t185.2\t1\n";
        assert_eq!(parse_diamond_tab(full13).len(), 1);
    }

    #[test]
    fn best_hit_per_query_max_bitscore() {
        let text = "a\tP1\t90\t100\t300\t200\t1\t300\t1\t100\t1e-40\t185.2\t1\n\
                    a\tP2\t99\t300\t300\t310\t1\t300\t1\t300\t1e-99\t520.0\t1\n\
                    b\tP3\t88\t90\t100\t120\t5\t95\t8\t98\t1e-20\t95.0\t1\n";
        let hits = parse_diamond_tab(text);
        let best = best_hit_per_query(&hits);
        assert_eq!(best["a"].sid, "P2");
        assert_eq!(best["b"].sid, "P3");
    }

    #[test]
    fn sanitize_drops_empty_records() {
        let dir = std::env::temp_dir().join(format!("pg_san_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let src = dir.join("src.fa");
        fs::write(&src, ">P1 desc\nMKVL\nMK\n>P2\n\n>P3 x\nAAA\n").unwrap();
        let dst = dir.join("clean.fa");
        let (kept, dropped) = sanitize_protein_db(&src, &dst).unwrap();
        assert_eq!((kept, dropped), (2, 1));
        let pairs = parse_fasta_pairs(&fs::read_to_string(&dst).unwrap());
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].0, "P1 desc");
        assert_eq!(pairs[0].1, "MKVLMK");
        assert_eq!(pairs[1].0, "P3 x");
        // dst 已存在 → 断点复用, 计数从 dst 读回。
        let (kept2, dropped2) = sanitize_protein_db(&src, &dst).unwrap();
        assert_eq!((kept2, dropped2), (2, 0));
    }

    #[test]
    fn protein_db_index() {
        let dir = std::env::temp_dir().join(format!("pg_idx_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("p.fa");
        fs::write(&fa, ">P1\nMKVL\n>P2\nAAA\n").unwrap();
        let idx = index_protein_db(&fa).unwrap();
        assert_eq!(idx["P1"], b"MKVL");
        assert_eq!(idx.len(), 2);
    }

    /// 转录本 300nt 命中 P1(200aa) 的 aa 60..160（plus 框）:
    /// N 缺 59 > C 缺 39 → 补 N → plus 框向左; budget = 59 × 3 × 3/2。
    #[test]
    fn gap_profile_missing_nterm() {
        let h = DiamondHit {
            qid: "t".into(),
            sid: "P1".into(),
            pident: 90.0,
            qlen: 300,
            slen: 200,
            qstart: 180,
            qend: 300,
            sstart: 60,
            send: 160,
            evalue: 1e-40,
            bitscore: 200.0,
            qframe: 1,
        };
        let gp = gap_profile(&h);
        assert!(gp.plus);
        assert!(!gp.dir_right, "缺 N 端更大 → plus 框向左");
        assert!(gp.missing_nterm());
        assert!(gp.missing_cterm());
        assert_eq!(gp.budget_nt(), 59 * 3 * 3 / 2);
    }

    /// 只缺 C 端（覆盖 aa 4..483 ↔ 1..160）→ 向右。
    #[test]
    fn gap_profile_missing_cterm_only() {
        let h = DiamondHit {
            qid: "t".into(),
            sid: "P1".into(),
            pident: 90.0,
            qlen: 483,
            slen: 200,
            qstart: 4,
            qend: 483,
            sstart: 1,
            send: 160,
            evalue: 1e-40,
            bitscore: 200.0,
            qframe: 1,
        };
        let gp = gap_profile(&h);
        assert!(gp.dir_right);
        assert_eq!(gp.budget_nt(), 40 * 3 * 3 / 2);
    }

    /// PAF 行（12 标准列够用; 解析只用 qname/tname/alnlen=f[10]）。
    fn paf(qname: &str, tname: &str, alnlen: usize) -> String {
        format!("{qname}\t100\t0\t50\t+\t{tname}\t100\t0\t50\t50\t{alnlen}\t60\n")
    }

    #[test]
    fn junction_support_counts_pairs() {
        // r1 一对: /1 命中 j0_a, /2 命中 j0_b → j0 +1
        let mut text = String::new();
        text.push_str(&paf("r1/1", "j0_a", 50));
        text.push_str(&paf("r1/2", "j0_b", 50));
        // r2 一对: 两端都只命中 j0_a → 不计
        text.push_str(&paf("r2/1", "j0_a", 50));
        text.push_str(&paf("r2/2", "j0_a", 50));
        // r3 一对命中 j1 的 a+b → j1 +1
        text.push_str(&paf("r3/1", "j1_b", 40));
        text.push_str(&paf("r3/2", "j1_a", 40));
        // 无关组名（非 j{i}_a/_b 形态）忽略
        text.push_str(&paf("r4/1", "chr1", 50));
        text.push_str(&paf("r4/2", "chr1", 50));
        let m = junction_support_from_paf(&text, 30);
        assert_eq!(m.get("j0"), Some(&1), "r1 跨越对计 1, r2 不计: {m:?}");
        assert_eq!(m.get("j1"), Some(&1));
        assert!(!m.contains_key("chr1"), "无 a/b 侧的 tname 不成组");
    }

    #[test]
    fn junction_support_strips_suffix_and_filters_short() {
        let mut text = String::new();
        // /1 /2 剥离后配对成功
        text.push_str(&paf("readA/1", "j5_a", 35));
        text.push_str(&paf("readA/2", "j5_b", 35));
        // alnlen < 30 的行忽略 → 未配对
        text.push_str(&paf("readB/1", "j5_a", 29));
        text.push_str(&paf("readB/2", "j5_b", 50));
        let m = junction_support_from_paf(&text, 30);
        assert_eq!(m.get("j5"), Some(&1), "readA 计 1, readB 短命中不计: {m:?}");
    }

    /// minus 框方向翻转: 只缺 C 端 → 向左。
    #[test]
    fn gap_profile_minus_directions() {
        let h = DiamondHit {
            qid: "t".into(),
            sid: "P1".into(),
            pident: 90.0,
            qlen: 300,
            slen: 200,
            qstart: 480,
            qend: 1,
            sstart: 160,
            send: 1,
            evalue: 1e-40,
            bitscore: 200.0,
            qframe: -1,
        };
        let gp = gap_profile(&h);
        assert!(!gp.plus);
        assert!(!gp.dir_right, "minus + 只缺 C → 向左");
    }
}
