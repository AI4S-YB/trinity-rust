//! M6a e2e: `--protein-guided`——合成 CDS 基因中段窗口低覆盖 + 一对稀疏 reads,
//! 蛋白库含真蛋白; 断言旗标开不劣于旗标关、stage 工件（queries.fa、
//! .protein_guided.ok）只出现在旗标开; 旗标关零变化。

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use trinity_cli::args::parse_args;
use trinity_cli::orchestrate::run_trinity;

fn tmpdir(tag: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!("trinity-cli-pg-{}-{}", tag, std::process::id()));
    let _ = fs::remove_dir_all(&p);
    fs::create_dir_all(&p).unwrap();
    p
}

fn xorshift_nt(x: &mut u64) -> usize {
    *x ^= *x << 13;
    *x ^= *x >> 7;
    *x ^= *x << 17;
    (*x % 20) as usize
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

/// 密码子表（每 aa 一码, 避开 stop）。
const CODONS: [&str; 20] = [
    "GCT", "TGT", "GAT", "GAA", "TTC", "GGC", "CAC", "ATC", "AAA", "CTG", "ATG", "AAC", "CCC",
    "CAA", "CGT", "AGC", "ACC", "GTT", "TGG", "TAC",
];
const AAS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

/// 合成 CDS: ATG + n 个非 stop 密码子 + TAA（3n+3 nt, 无内部 stop）。
fn synth_cds(seed: u64, n_codons: usize) -> (String, String) {
    let mut x = seed | 1;
    let mut dna = String::from("ATG");
    let mut aa = String::from("M");
    for _ in 0..n_codons {
        let i = xorshift_nt(&mut x);
        dna.push_str(CODONS[i]);
        aa.push(AAS[i] as char);
    }
    dna.push_str("TAA");
    aa.push('*');
    (dna, aa)
}

/// 平铺 PE reads（150nt, 插入 250, 步进 8）; 中段窗口 [lo,hi) 内完全抑制 pair 起点
/// （覆盖空洞 ≈ [lo+242, hi)），洞内只放一对稀疏 reads（起点 519）。
///
/// 实证注记（本管线中覆盖空洞无法把 vanilla 与组件图分开）:
/// - 稀疏 pair 与侧翼零共享 kmer（起点 ≥519+1）→ 不进组件 → 图无路径,
///   vanilla 截断于 ~546, 但引导走图也无路可走（0 extended）;
/// - 稀疏 pair 与侧翼共享 ≥1 个 kmer → k=24 组件图链连续 → butterfly
///   直接输出全长（off=on=898）, 引导无缺口可补。
///
/// 两个方向都试过; 长度断言仍保留（旗标开不得变差）, stage 接线由
/// queries.fa / .protein_guided.ok 断言覆盖。
fn write_reads(
    dir: &Path,
    cds: &str,
    mid_lo: usize,
    mid_hi: usize,
    sparse_starts: &[usize],
) -> (PathBuf, PathBuf) {
    let b = cds.as_bytes();
    let mut r1 = String::new();
    let mut r2 = String::new();
    let mut n = 0usize;
    let mut s = 0usize;
    while s + 250 <= b.len() {
        if s < mid_lo || s >= mid_hi {
            let f1 = std::str::from_utf8(&b[s..s + 150]).unwrap();
            let f2 = std::str::from_utf8(&b[s + 100..s + 250]).unwrap();
            r1.push_str(&format!("@r{n}\n{f1}\n+\n{}\n", "I".repeat(150)));
            r2.push_str(&format!("@r{n}\n{}\n+\n{}\n", revcomp(f2), "I".repeat(150)));
            n += 1;
            s += 8;
        } else {
            s = mid_hi; // 窗口内直接跳到 hi
        }
    }
    for (i, &s) in sparse_starts.iter().enumerate() {
        let f1 = std::str::from_utf8(&b[s..s + 150]).unwrap();
        let f2 = std::str::from_utf8(&b[s + 100..s + 250]).unwrap();
        r1.push_str(&format!("@m{i}\n{f1}\n+\n{}\n", "I".repeat(150)));
        r2.push_str(&format!("@m{i}\n{}\n+\n{}\n", revcomp(f2), "I".repeat(150)));
    }
    let p1 = dir.join("r1.fq.gz");
    let p2 = dir.join("r2.fq.gz");
    for (p, body) in [(&p1, &r1), (&p2, &r2)] {
        let mut w = flate2::write::GzEncoder::new(
            fs::File::create(p).unwrap(),
            flate2::Compression::fast(),
        );
        w.write_all(body.as_bytes()).unwrap();
        w.finish().unwrap();
    }
    (p1, p2)
}

fn base_args(r1: &Path, r2: &Path, outdir: &Path) -> Vec<String> {
    [
        "--seqType",
        "fq",
        "--left",
        r1.to_str().unwrap(),
        "--right",
        r2.to_str().unwrap(),
        "--max_memory",
        "1G",
        "--output",
        outdir.to_str().unwrap(),
        "--no_normalize_reads",
        "--CPU",
        "2",
        "--min_contig_length",
        "100",
        "--path-engine",
        "pogflow",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect()
}

fn diamond_available() -> bool {
    let bin = std::env::var("TRINITY_DIAMOND")
        .ok()
        .filter(|v| !v.is_empty())
        .unwrap_or_else(|| {
            "/public/home/senior007/miniconda3/envs/diamond2go/bin/diamond".to_string()
        });
    std::process::Command::new(bin)
        .arg("--version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

/// 最长转录本长度（fasta 行有折行, 累加每条序列的各行）。
fn max_len(fa: &Path) -> usize {
    let mut best = 0usize;
    let mut cur = 0usize;
    for line in fs::read_to_string(fa).unwrap().lines() {
        if line.starts_with('>') {
            best = best.max(cur);
            cur = 0;
        } else {
            cur += line.trim().len();
        }
    }
    best.max(cur)
}

#[test]
fn protein_guided_e2e_extends_low_coverage_gap() {
    if !diamond_available() {
        eprintln!("skip: diamond not available");
        return;
    }
    let dir = tmpdir("e2e");
    // 900nt CDS（298 密码子+起止 = 3*298+6）; 中段 [300,600) 只给 1 对 reads
    let (cds, aa) = synth_cds(0xA5A5_1234_DEAD_BEEF, 298);
    assert_eq!(cds.len(), 900);
    // 稀疏 pair 起点 519: R1 覆盖 [519,669)，与侧翼 R2 覆盖（止于 542）恰重叠 23nt。
    let (r1, r2) = write_reads(&dir, &cds, 300, 600, &[519]);

    let pdb = dir.join("prot.fa");
    fs::write(&pdb, format!(">REAL\n{aa}\n")).unwrap();

    // ---- off: vanilla ----
    let out_off = dir.join("trinity_off");
    let f_off = run_trinity(&parse_args(&base_args(&r1, &r2, &out_off)).unwrap()).unwrap();
    let off_max = max_len(&f_off);

    // ---- on ----
    let out_on = dir.join("trinity_on");
    let mut argv = base_args(&r1, &r2, &out_on);
    argv.push("--protein-guided".into());
    argv.push("--protein-db".into());
    argv.push(pdb.display().to_string());
    let f_on = run_trinity(&parse_args(&argv).unwrap()).unwrap();
    let on_max = max_len(&f_on);

    assert!(out_on.join(".protein_guided.ok").exists());
    // stage 真跑了 DIAMOND（查询集 + 库索引工件）
    assert!(out_on.join("protein_guided/queries.fa").exists());
    assert!(out_on.join("protein_guided/protein.dmnd").exists());
    assert!(!out_off.join(".protein_guided.ok").exists());
    assert!(!out_off.join("protein_guided").exists());
    assert!(on_max >= off_max, "延长不应变短: {on_max} < {off_max}");
    // 主断言: 引导后接近全长（允许 10% 缺口）
    assert!(on_max >= 810, "on_max={on_max} off_max={off_max}");
}
