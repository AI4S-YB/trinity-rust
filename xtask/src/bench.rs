//! xtask bench —— F1 基准: 组装(缓存) → minimap2 → R/P/F1。
//! 口径锁定见 M0 计划; 纯函数 compute_metrics 便于单测。

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use sha2::{Digest, Sha256};

pub const MINIMAP2: &str = "/public/home/senior007/miniconda3/envs/cactus/bin/minimap2";

#[derive(Debug, Default, PartialEq)]
pub struct Metrics {
    pub refs_total: usize,
    pub r_full: f64,    // ≥95%cov & ≥95%id
    pub r70: f64,       // ≥70%cov & ≥90%id
    pub precision: f64, // 碱基级
    pub f1: f64,
    pub n50: usize,
    pub denovo_total: usize,
}

/// PAF 行: q qlen qs qe strand t tlen ts te mch aln ...
/// ref_lens: 参照名→长度; denovo_lens: de novo 名→长度
pub fn compute_metrics(
    paf: &str,
    ref_lens: &HashMap<String, usize>,
    denovo_lens: &HashMap<String, usize>,
) -> Metrics {
    let mut best: HashMap<&str, Option<(f64, f64)>> = HashMap::new(); // ref -> best (cov, ident)（按 cov*ident 取）; None=未见命中
    let mut qhit: HashMap<&str, usize> = HashMap::new();
    for line in paf.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 11 {
            continue;
        }
        let (q, t) = (f[0], f[5]);
        let qs: usize = f[2].parse().unwrap_or(0);
        let qe: usize = f[3].parse().unwrap_or(0);
        let ts: usize = f[7].parse().unwrap_or(0);
        let te: usize = f[8].parse().unwrap_or(0);
        let mch: usize = f[9].parse().unwrap_or(0);
        let aln: usize = f[10].parse().unwrap_or(0);
        let ident = if aln > 0 {
            mch as f64 / aln as f64
        } else {
            0.0
        };
        if let Some(&tl) = ref_lens.get(t) {
            let cov = te.saturating_sub(ts) as f64 / tl as f64;
            let score = cov * ident;
            let e = best.entry(t).or_insert(None);
            #[allow(clippy::unnecessary_map_or)]
            if e.map_or(true, |(c, i)| c * i < score) {
                *e = Some((cov, ident));
            }
        }
        if ident >= 0.95 {
            let span = qe.saturating_sub(qs);
            let e = qhit.entry(q).or_insert(0);
            if span > *e {
                *e = span;
            }
        }
    }
    let refs_total = ref_lens.len();
    let mut r_full = 0usize;
    let mut r70 = 0usize;
    for &(cov, ident) in best.values().flatten() {
        if cov >= 0.95 && ident >= 0.95 {
            r_full += 1;
        }
        if cov >= 0.70 && ident >= 0.90 {
            r70 += 1;
        }
    }
    let total_bases: usize = denovo_lens.values().sum();
    let hit_bases: usize = qhit.values().sum();
    let precision = if total_bases > 0 {
        hit_bases as f64 / total_bases as f64
    } else {
        0.0
    };
    let r = if refs_total > 0 {
        r_full as f64 / refs_total as f64
    } else {
        0.0
    };
    let f1 = if r + precision > 0.0 {
        2.0 * r * precision / (r + precision)
    } else {
        0.0
    };
    // N50
    let mut lens: Vec<usize> = denovo_lens.values().copied().collect();
    lens.sort_unstable_by(|a, b| b.cmp(a));
    let half = total_bases as f64 / 2.0;
    let (mut acc, mut n50) = (0usize, 0usize);
    for &l in &lens {
        acc += l;
        if acc as f64 >= half {
            n50 = l;
            break;
        }
    }
    Metrics {
        refs_total,
        r_full: r,
        r70: if refs_total > 0 {
            r70 as f64 / refs_total as f64
        } else {
            0.0
        },
        precision,
        f1,
        n50,
        denovo_total: denovo_lens.len(),
    }
}

fn read_fasta_lens(path: &std::path::Path) -> std::io::Result<HashMap<String, usize>> {
    let mut m = HashMap::new();
    let mut name: Option<String> = None;
    let mut len = 0usize;
    for line in std::io::BufRead::lines(&mut std::io::BufReader::new(std::fs::File::open(path)?)) {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(n) = name {
                m.insert(n, len);
            }
            name = Some(rest.split_whitespace().next().unwrap_or("").to_string());
            len = 0;
        } else {
            len += line.trim().len();
        }
    }
    if let Some(n) = name {
        m.insert(n, len);
    }
    Ok(m)
}

/// 缓存键 = sha256(reads 文件内容 ‖ 0x00 ‖ 二进制内容 ‖ 0x00 ‖ 参数) 前 16 hex
fn cache_key(args: &[&str], reads: &[Vec<u8>], bin: &[u8]) -> String {
    let mut h = Sha256::new();
    for r in reads {
        h.update(r);
        h.update([0u8]);
    }
    h.update(bin);
    h.update([0u8]);
    for a in args {
        h.update(a.as_bytes());
        h.update([0u8]);
    }
    let out = h.finalize();
    out.iter().take(8).map(|b| format!("{b:02x}")).collect()
}

fn bench_root() -> PathBuf {
    std::env::var_os("ATHAL")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("/storage/home/senior007/test/athal"))
        .join("bench")
}

pub fn run(argv: &[String]) -> i32 {
    // 用法: bench --tag dev1 --reads L.fq.gz,R.fq.gz --ref refs/dev1.fa [-- passthrough...]
    let mut tag = String::new();
    let mut reads = String::new();
    let mut reff = String::new();
    let mut prot_ref = String::new();
    let mut i = 0;
    let mut passthrough: Vec<String> = Vec::new();
    while i < argv.len() {
        match argv[i].as_str() {
            "--tag" => {
                tag = argv.get(i + 1).cloned().unwrap_or_default();
                i += 2;
            }
            "--reads" => {
                reads = argv.get(i + 1).cloned().unwrap_or_default();
                i += 2;
            }
            "--ref" => {
                reff = argv.get(i + 1).cloned().unwrap_or_default();
                i += 2;
            }
            "--prot-ref" => {
                prot_ref = argv.get(i + 1).cloned().unwrap_or_default();
                i += 2;
            }
            "--" => {
                passthrough = argv[i + 1..].to_vec();
                break;
            }
            _ => {
                eprintln!("未知 bench 参数: {}", argv[i]);
                return 2;
            }
        }
    }
    if tag.is_empty() || reads.is_empty() || reff.is_empty() {
        eprintln!("用法: cargo xtask bench --tag <名> --reads <L.fq.gz,R.fq.gz> --ref <ref.fa> [-- <trinity-cli 参数...>]");
        return 2;
    }
    let read_files: Vec<&str> = reads.split(',').collect();
    // 默认参数注入（未显式给定时）
    fn has(pt: &[String], name: &str) -> bool {
        pt.iter().any(|a| {
            a == name || a.starts_with(&format!("{name}=")) || a.starts_with(&format!("{name} "))
        })
    }
    if !has(&passthrough, "--CPU") {
        passthrough.push("--CPU".into());
        passthrough.push("8".into());
    }
    if !has(&passthrough, "--max_memory") {
        passthrough.push("--max_memory".into());
        passthrough.push("8G".into());
    }
    // 键与缓存
    let bin_path = concat!(env!("CARGO_MANIFEST_DIR"), "/../target/release/trinity-cli");
    if !std::path::Path::new(bin_path).exists() {
        eprintln!("Error: 未找到 {bin_path} —— 先 cargo build --release -p trinity-cli");
        return 1;
    }
    let bin = std::fs::read(bin_path).unwrap_or_default();
    let rcontents: Vec<Vec<u8>> = read_files
        .iter()
        .map(|p| std::fs::read(p).unwrap_or_default())
        .collect();
    let pt: Vec<&str> = passthrough.iter().map(|s| s.as_str()).collect();
    let key = cache_key(&pt, &rcontents, &bin);
    let dir = bench_root().join("cache").join(&key);
    let final_fa = dir.join("out.Trinity.fasta");
    let bin_path_s = std::path::PathBuf::from(bin_path);
    if !final_fa.exists() {
        let _ = std::fs::remove_dir_all(dir.join("out"));
        std::fs::create_dir_all(&dir).unwrap();
        let mut cmd = std::process::Command::new(&bin_path_s);
        cmd.args([
            "--seqType",
            "fq",
            "--left",
            read_files[0],
            "--right",
            read_files.get(1).unwrap_or(&read_files[0]),
        ]);
        cmd.args(["--output"]).arg(dir.join("out"));
        cmd.args(&passthrough);
        let st = cmd
            .status()
            .unwrap_or_else(|e| panic!("无法启动 trinity-cli: {e}"));
        if !st.success() {
            eprintln!("组装失败 ({st})");
            return 1;
        }
    } else {
        println!("[bench] 组装缓存命中: {}", dir.display());
    }
    // minimap2
    let paf = std::process::Command::new(MINIMAP2)
        .args(["-x", "asm5", "-t", "8", "--secondary=no"])
        .arg(&reff)
        .arg(&final_fa)
        .output()
        .expect("无法启动 minimap2");
    if !paf.status.success() {
        eprintln!(
            "minimap2 失败 ({})\n{}",
            paf.status,
            String::from_utf8_lossy(&paf.stderr)
        );
        return 1;
    }
    let paf = String::from_utf8_lossy(&paf.stdout).to_string();
    let ref_lens = read_fasta_lens(std::path::Path::new(&reff)).unwrap();
    let dn = read_fasta_lens(&final_fa).unwrap();
    let m = compute_metrics(&paf, &ref_lens, &dn);
    // 蛋白级回收: DIAMOND blastx vs --prot-ref（缺省: --ref 最长 ORF 翻译, 缓存于 cache 目录）
    let (p_full, p70, prot_total) = {
        let pm_dir = dir.join("protein_metric");
        let _ = std::fs::create_dir_all(&pm_dir);
        let prot_fa = if prot_ref.is_empty() {
            let p = pm_dir.join("ref_orf.fa");
            if !p.exists() {
                translate_longest_orf(std::path::Path::new(&reff), &p)
                    .unwrap_or_else(|e| panic!("translate_longest_orf: {e}"));
            }
            p
        } else {
            std::path::PathBuf::from(&prot_ref)
        };
        let tab = diamond_blastx_cached(&final_fa, &prot_fa, 8).unwrap_or_default();
        // 蛋白参考长度表
        let mut lens: HashMap<String, usize> = HashMap::new();
        let mut cur: Option<String> = None;
        for line in std::fs::read_to_string(&prot_fa)
            .unwrap_or_default()
            .lines()
        {
            if let Some(h) = line.strip_prefix('>') {
                cur = Some(h.split_whitespace().next().unwrap_or("").to_string());
                lens.entry(cur.clone().unwrap()).or_insert(0);
            } else if let Some(nm) = &cur {
                *lens.entry(nm.clone()).or_insert(0) += line.trim().len();
            }
        }
        protein_recall(&tab, &lens)
    };
    // 指标落盘 cache/<key>/metrics.json（手写 JSON, 与 stdout 同口径; sweep 汇总用）
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(
        dir.join("metrics.json"),
        format!(
            "{{\"tag\":\"{tag}\",\"date\":\"{}\",\"r_full\":{:.4},\"r70\":{:.4},\"precision\":{:.4},\"f1\":{:.4},\"n50\":{},\"n\":{},\"p_full\":{:.4},\"p70\":{:.4},\"prot_total\":{},\"args\":\"{}\"}}\n",
            chrono_free_date(),
            m.r_full,
            m.r70,
            m.precision,
            m.f1,
            m.n50,
            m.denovo_total,
            p_full,
            p70,
            prot_total,
            passthrough.join(" ")
        ),
    )
    .unwrap();
    // 报告 + 日志
    println!("== bench[{tag}] {key} ==\n  转录本 {}/{} | R={:.1}% R70={:.1}% P={:.1}% F1={:.1}% | N50={} | P_full={:.1}% P70={:.1}%",
        m.denovo_total, m.refs_total, 100.0*m.r_full, 100.0*m.r70, 100.0*m.precision, 100.0*m.f1, m.n50, 100.0 * p_full, 100.0 * p70);
    let log = concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/opt-log.md");
    // 文件头并发注释（幂等: 已存在则不重复）
    {
        let marker = "<!-- 并发注意: 禁止同 key 并发 bench（缓存目录互踩） -->\n";
        let cur = std::fs::read_to_string(log).unwrap_or_default();
        if !cur.contains("禁止同 key 并发 bench") {
            std::fs::write(log, format!("{marker}{cur}")).unwrap();
        }
    }
    let row = format!(
        "| {} | {} | R={:.1}% P={:.1}% **F1={:.1}%** | N50={} n={} | {} |\n",
        chrono_free_date(),
        tag,
        100.0 * m.r_full,
        100.0 * m.precision,
        100.0 * m.f1,
        m.n50,
        m.denovo_total,
        passthrough.join(" ")
    );
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(log)
        .unwrap();
    use std::io::Write;
    f.write_all(row.as_bytes()).unwrap();

    // M7 探针: top-5 耗时段（两份 .stage_timings 合并）
    {
        let mut stages: Vec<(String, f64)> = Vec::new();
        for p in [
            dir.join("out/.stage_timings"),
            dir.join("out/chrysalis/.stage_timings"),
        ] {
            if let Ok(text) = std::fs::read_to_string(&p) {
                for line in text.lines() {
                    let f: Vec<&str> = line.split('\t').collect();
                    if f.len() == 2 {
                        if let Ok(s) = f[1].parse::<f64>() {
                            stages.push((f[0].to_string(), s));
                        }
                    }
                }
            }
        }
        stages.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        if !stages.is_empty() {
            let top: Vec<String> = stages
                .iter()
                .take(5)
                .map(|(n, s)| format!("{n} {s:.0}s"))
                .collect();
            println!("  Top stages: {}", top.join(" | "));
        }
    }
    0
}

fn chrono_free_date() -> String {
    std::process::Command::new("date")
        .arg("+%Y-%m-%d")
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_else(|_| "unknown".into())
}

// ---------------------------------------------------------------------------
// 蛋白级回收指标: DIAMOND blastx 组装产物 vs 蛋白参考
// ---------------------------------------------------------------------------

fn diamond_bin() -> &'static str {
    static DEFAULT: &str = "/public/home/senior007/miniconda3/envs/diamond2go/bin/diamond";
    match std::env::var("TRINITY_DIAMOND") {
        Ok(v) if !v.is_empty() => Box::leak(v.into_boxed_str()),
        _ => DEFAULT,
    }
}

/// DIAMOND blastx → tabular（12 列默认）。dmnd 按参考内容 sha256 缓存于 bench cache。
fn diamond_blastx_cached(
    query_fa: &Path,
    prot_ref: &Path,
    threads: usize,
) -> std::io::Result<String> {
    let bytes = std::fs::read(prot_ref)?;
    let mut h = Sha256::new();
    h.update(&bytes);
    let kh: String = h
        .finalize()
        .iter()
        .take(6)
        .map(|b| format!("{b:02x}"))
        .collect();
    std::fs::create_dir_all(bench_root().join("cache"))?;
    let dmnd = bench_root()
        .join("cache")
        .join(format!("protref_{kh}.dmnd"));
    if !dmnd.exists() {
        let st = std::process::Command::new(diamond_bin())
            .arg("makedb")
            .arg("--in")
            .arg(prot_ref)
            .arg("-d")
            .arg(&dmnd)
            .arg("--threads")
            .arg(threads.to_string())
            .stdout(std::process::Stdio::null())
            .status()
            .map_err(|e| std::io::Error::other(format!("spawn diamond: {e}")))?;
        if !st.success() {
            return Err(std::io::Error::other("diamond makedb failed"));
        }
    }
    let out_path = std::env::temp_dir().join(format!("xt_blastx_{}_{}", kh, std::process::id()));
    let st = std::process::Command::new(diamond_bin())
        .arg("blastx")
        .arg("-q")
        .arg(query_fa)
        .arg("-d")
        .arg(&dmnd)
        .arg("--threads")
        .arg(threads.to_string())
        .arg("-e")
        .arg("1e-5")
        .arg("--max-target-seqs")
        .arg("20")
        .arg("--outfmt")
        .arg("6")
        .arg("-o")
        .arg(&out_path)
        .stdout(std::process::Stdio::null())
        .status()
        .map_err(|e| std::io::Error::other(format!("spawn diamond: {e}")))?;
    if !st.success() {
        return Err(std::io::Error::other("diamond blastx failed"));
    }
    let text = std::fs::read_to_string(&out_path)?;
    let _ = std::fs::remove_file(&out_path);
    Ok(text)
}

/// 标准密码子表翻译（translate_longest_orf 用; 3 字母以外/含简并碱基 → X）。
fn translate_codon(codon: &[u8]) -> u8 {
    let n = |b: u8| match b {
        b'T' => 0,
        b'C' => 1,
        b'A' => 2,
        b'G' => 3,
        _ => 4,
    };
    match codon.len() {
        3 => match (n(codon[0]), n(codon[1]), n(codon[2])) {
            (0, 0, 0) => b'F',
            (0, 0, 1) => b'F',
            (0, 0, 2) => b'L',
            (0, 0, 3) => b'L',
            (0, 1, 0) => b'S',
            (0, 1, 1) => b'S',
            (0, 1, 2) => b'S',
            (0, 1, 3) => b'S',
            (0, 2, 0) => b'Y',
            (0, 2, 1) => b'Y',
            (0, 2, 2) => b'*',
            (0, 2, 3) => b'*',
            (0, 3, 0) => b'C',
            (0, 3, 1) => b'C',
            (0, 3, 2) => b'*',
            (0, 3, 3) => b'W',
            (1, 0, 0) => b'L',
            (1, 0, 1) => b'L',
            (1, 0, 2) => b'L',
            (1, 0, 3) => b'L',
            (1, 1, 0) => b'P',
            (1, 1, 1) => b'P',
            (1, 1, 2) => b'P',
            (1, 1, 3) => b'P',
            (1, 2, 0) => b'H',
            (1, 2, 1) => b'H',
            (1, 2, 2) => b'Q',
            (1, 2, 3) => b'Q',
            (1, 3, 0) => b'R',
            (1, 3, 1) => b'R',
            (1, 3, 2) => b'R',
            (1, 3, 3) => b'R',
            (2, 0, 0) => b'I',
            (2, 0, 1) => b'I',
            (2, 0, 2) => b'I',
            (2, 0, 3) => b'M',
            (2, 1, 0) => b'T',
            (2, 1, 1) => b'T',
            (2, 1, 2) => b'T',
            (2, 1, 3) => b'T',
            (2, 2, 0) => b'N',
            (2, 2, 1) => b'N',
            (2, 2, 2) => b'K',
            (2, 2, 3) => b'K',
            (2, 3, 0) => b'S',
            (2, 3, 1) => b'S',
            (2, 3, 2) => b'R',
            (2, 3, 3) => b'R',
            (3, 0, 0) => b'V',
            (3, 0, 1) => b'V',
            (3, 0, 2) => b'V',
            (3, 0, 3) => b'V',
            (3, 1, 0) => b'A',
            (3, 1, 1) => b'A',
            (3, 1, 2) => b'A',
            (3, 1, 3) => b'A',
            (3, 2, 0) => b'D',
            (3, 2, 1) => b'D',
            (3, 2, 2) => b'E',
            (3, 2, 3) => b'E',
            (3, 3, 0) => b'G',
            (3, 3, 1) => b'G',
            (3, 3, 2) => b'G',
            (3, 3, 3) => b'G',
            _ => b'X',
        },
        _ => b'X',
    }
}

/// 参考转录本 → 最长 ORF 蛋白集（≥50aa; 3 框按 stop 分段, M 锚定优先, 否则整段）。
/// 返回写出的蛋白记录数。写 60 列折行。
fn translate_longest_orf(cdna_fa: &Path, out_fa: &Path) -> std::io::Result<usize> {
    let text = std::fs::read_to_string(cdna_fa)?;
    let mut out = String::new();
    let mut n = 0usize;
    let mut name: Option<String> = None;
    let mut seq = String::new();
    let flush = |name: &Option<String>, seq: &mut String, out: &mut String, n: &mut usize| {
        if let Some(nm) = name {
            let mut best: Vec<u8> = Vec::new();
            // 输入序列先大写化（参考 fasta 常为 soft-masked）。
            let su = seq.to_ascii_uppercase();
            for frame in 0..3 {
                let b = su.as_bytes();
                let mut prot: Vec<u8> = Vec::new();
                let mut i = frame;
                while i + 3 <= b.len() {
                    prot.push(translate_codon(&b[i..i + 3]));
                    i += 3;
                }
                for seg in prot.split(|&c| c == b'*') {
                    let cand: &[u8] = match seg.iter().position(|&c| c == b'M') {
                        Some(p) => &seg[p..],
                        None => seg,
                    };
                    if cand.len() > best.len() {
                        best = cand.to_vec();
                    }
                }
            }
            if best.len() >= 50 {
                *n += 1;
                out.push('>');
                out.push_str(nm);
                out.push('\n');
                for ch in best.chunks(60) {
                    out.push_str(&String::from_utf8_lossy(ch));
                    out.push('\n');
                }
            }
        }
        seq.clear();
    };
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            flush(&name.take(), &mut seq, &mut out, &mut n);
            name = Some(h.split_whitespace().next().unwrap_or("").to_string());
        } else {
            seq.push_str(line.trim());
        }
    }
    flush(&name, &mut seq, &mut out, &mut n);
    std::fs::write(out_fa, out)?;
    Ok(n)
}

/// blastx tabular → (P_full, P70, 蛋白总数)。
/// 同一 (q,s) 的 HSP 求 subject 并集跨度; 每蛋白取 c×pid 最高的 query。
fn protein_recall(tab: &str, prot_lens: &HashMap<String, usize>) -> (f64, f64, usize) {
    use std::collections::HashMap as Map;
    /// (q,s) → HSP subject 区间集 + 最高 pident
    type HspAcc = (Vec<(i64, i64)>, f64);
    let mut pair: Map<(String, String), HspAcc> = Map::new();
    for line in tab.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        let (Ok(pid), Ok(ss), Ok(se)) = (
            f[2].parse::<f64>(),
            f[8].parse::<i64>(),
            f[9].parse::<i64>(),
        ) else {
            continue;
        };
        let pid = pid / 100.0; // tabular pident 为百分数 → 归一化
        let e = pair
            .entry((f[0].to_string(), f[1].to_string()))
            .or_insert((Vec::new(), 0.0));
        e.0.push((ss.min(se), ss.max(se)));
        e.1 = e.1.max(pid);
    }
    let mut best: Map<String, (f64, f64)> = Map::new();
    for ((_q, s), (ivs, pid)) in &pair {
        let Some(&len) = prot_lens.get(s) else {
            continue;
        };
        let mut ivs = ivs.clone();
        ivs.sort_unstable();
        let mut cov = 0i64;
        let mut cur: Option<(i64, i64)> = None;
        for &(a, b) in &ivs {
            match cur {
                None => cur = Some((a, b)),
                Some((cs, ce)) if a <= ce + 1 => cur = Some((cs, ce.max(b))),
                Some((cs, ce)) => {
                    cov += ce - cs + 1;
                    cur = Some((a, b));
                }
            }
        }
        if let Some((cs, ce)) = cur {
            cov += ce - cs + 1;
        }
        let c = cov as f64 / len.max(1) as f64;
        let score = c * pid;
        match best.get(s) {
            Some((bc, bp)) if *bc * *bp >= score => {}
            _ => {
                best.insert(s.clone(), (c, *pid));
            }
        }
    }
    let total = prot_lens.len();
    if total == 0 {
        return (0.0, 0.0, 0);
    }
    let full = best
        .values()
        .filter(|&&(c, p)| c >= 0.95 && p >= 0.95)
        .count();
    let p70 = best
        .values()
        .filter(|&&(c, p)| c >= 0.70 && p >= 0.90)
        .count();
    (full as f64 / total as f64, p70 as f64 / total as f64, total)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hm(pairs: &[(&str, usize)]) -> HashMap<String, usize> {
        pairs.iter().map(|(a, b)| (a.to_string(), *b)).collect()
    }

    #[test]
    fn metrics_hand_case() {
        // 参照 r1=1000bp, r2=500bp; de novo d1=900bp, d2=400bp, d3=100bp(垃圾)
        let ref_l = hm(&[("r1", 1000), ("r2", 500)]);
        let dn_l = hm(&[("d1", 900), ("d2", 400), ("d3", 100)]);
        // d1 覆盖 r1 960/1000=0.96, id=0.97 → r1 全长
        // d2 覆盖 r2 400/500=0.80, id=0.95 → R70 计入
        let paf = "d1\t900\t0\t960\t+\tr1\t1000\t40\t1000\t941\t970\n\
                   d2\t400\t0\t400\t+\tr2\t500\t100\t500\t390\t400\n\
                   d3\t100\t0\t50\t+\tr1\t1000\t0\t50\t30\t50\n";
        // d3 id=30/50=0.6 <0.95 不计精确度; d1/d2 全量计入
        let m = compute_metrics(paf, &ref_l, &dn_l);
        assert_eq!(m.refs_total, 2);
        assert!((m.r_full - 0.5).abs() < 1e-9); // r1 全长
        assert!((m.r70 - 1.0).abs() < 1e-9); // r1+r2
        assert!((m.precision - (960.0 + 400.0) / 1400.0).abs() < 1e-9);
        assert!((m.f1 - 2.0 * 0.5 * (1360.0 / 1400.0) / (0.5 + 1360.0 / 1400.0)).abs() < 1e-9);
    }

    #[test]
    fn n50_case() {
        let dn_l = hm(&[("a", 100), ("b", 900), ("c", 400)]);
        let m = compute_metrics("", &HashMap::new(), &dn_l);
        assert_eq!(m.n50, 900); // total 1400, half 700, 900 先跨过
        assert!((m.precision - 0.0).abs() < 1e-9);
        assert_eq!(m.denovo_total, 3);
    }

    #[test]
    fn cache_key_stable_and_sensitive() {
        let a = cache_key(
            &["--CPU", "8"],
            &[b"reads1".to_vec(), b"reads2".to_vec()],
            b"BIN",
        );
        let b = cache_key(
            &["--CPU", "8"],
            &[b"reads1".to_vec(), b"reads2".to_vec()],
            b"BIN",
        );
        let c = cache_key(
            &["--CPU", "4"],
            &[b"reads1".to_vec(), b"reads2".to_vec()],
            b"BIN",
        );
        let d = cache_key(
            &["--CPU", "8"],
            &[b"readsX".to_vec(), b"reads2".to_vec()],
            b"BIN",
        );
        assert_eq!(a, b);
        assert_ne!(a, c); // 参数变 → 键变
        assert_ne!(a, d); // reads 变 → 键变
    }

    /// 手算用例: P1=100aa, P2=200aa。
    /// q1↔P1 pident 97, HSP [1,60]+[61,100] → 并集 cov=1.0
    /// q2↔P1 pident 99, HSP [1,70] → cov=0.70
    /// P1 最佳 = q1 (1.0×0.97 > 0.70×0.99) → full ✓ p70 ✓
    /// q3↔P2 pident 80, HSP [1,200] → cov=1.0 但 pid<90 → 都不满足
    #[test]
    fn protein_recall_hand_case() {
        let tab = "q1\tP1\t97.0\t60\t200\t100\t1\t180\t1\t60\t0.0\t100\n\
                   q1\tP1\t97.0\t40\t200\t100\t181\t300\t61\t100\t0.0\t70\n\
                   q2\tP1\t99.0\t70\t200\t100\t1\t210\t1\t70\t0.0\t110\n\
                   q3\tP2\t80.0\t200\t600\t200\t1\t600\t1\t200\t0.0\t300\n";
        let lens: HashMap<String, usize> =
            [("P1".to_string(), 100usize), ("P2".to_string(), 200usize)]
                .into_iter()
                .collect();
        let (p_full, p70, total) = protein_recall(tab, &lens);
        assert_eq!(total, 2);
        assert!((p_full - 0.5).abs() < 1e-9);
        assert!((p70 - 0.5).abs() < 1e-9);
    }

    /// 最长 ORF: 记录1 = ATG + 50×GCT + TAA → M+50A（51aa ≥50 保留）;
    /// 记录2 = 短 ORF → 过滤。
    #[test]
    fn longest_orf_translation() {
        let dir = std::env::temp_dir().join(format!("xt_orf_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        // 第 2 个密码子小写: 验证 soft-masked 输入也被大写化翻译（不得出 X）。
        let cds = format!(">T1\nATGgct{}TAA\n>T2\nATGAAATAA\n", "GCT".repeat(49));
        let in_fa = dir.join("in.fa");
        std::fs::write(&in_fa, cds).unwrap();
        let out_fa = dir.join("out.fa");
        let n = translate_longest_orf(&in_fa, &out_fa).unwrap();
        assert_eq!(n, 1);
        let text = std::fs::read_to_string(&out_fa).unwrap();
        assert!(text.contains(">T1"));
        let seq: String = text.lines().skip(1).collect();
        assert_eq!(seq.len(), 51);
        assert!(seq.starts_with('M') && seq.ends_with('A'));
        assert_eq!(seq.matches('A').count(), 50);
        assert!(!seq.contains('X'));
    }
}
