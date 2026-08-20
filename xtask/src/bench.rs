//! xtask bench —— F1 基准: 组装(缓存) → minimap2 → R/P/F1。
//! 口径锁定见 M0 计划; 纯函数 compute_metrics 便于单测。

use std::collections::HashMap;
use std::path::PathBuf;

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
            if (*e).map_or(true, |(c, i)| score > c * i) {
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
    // 报告 + 日志
    println!("== bench[{tag}] {key} ==\n  转录本 {}/{} | R={:.1}% R70={:.1}% P={:.1}% F1={:.1}% | N50={}",
        m.denovo_total, m.refs_total, 100.0*m.r_full, 100.0*m.r70, 100.0*m.precision, 100.0*m.f1, m.n50);
    let log = concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/opt-log.md");
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
    0
}

fn chrono_free_date() -> String {
    std::process::Command::new("date")
        .arg("+%Y-%m-%d")
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_else(|_| "unknown".into())
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
}
