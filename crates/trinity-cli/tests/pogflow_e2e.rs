//! M2.1 e2e: `--path-engine pogflow` 全流程（合成 PE reads → butterfly POG 图
//! → 流分解 → Trinity.fasta）。复用 flow_e2e 的合成集手法。

use std::fs;
use std::path::{Path, PathBuf};

use trinity_cli::args::parse_args;
use trinity_cli::orchestrate::run_trinity;

fn tmpdir(tag: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!(
        "trinity-cli-pogflow-{}-{}",
        tag,
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&p);
    fs::create_dir_all(&p).unwrap();
    p
}

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

/// 合成 PE reads: 600nt 参考（xorshift）× 25 对, FR 插入 250。
fn write_synthetic_pe(dir: &Path) -> (PathBuf, PathBuf) {
    let mut x: u64 = 0x243F6A8885A308D3;
    let reference: String = (0..600).map(|_| xorshift_nt(&mut x)).collect();
    let rb = reference.as_bytes();
    let mut r1 = String::new();
    let mut r2 = String::new();
    for i in 0..25usize {
        let s = i * 14;
        let f1 = std::str::from_utf8(&rb[s..s + 100]).unwrap();
        let f2 = std::str::from_utf8(&rb[s + 150..s + 250]).unwrap();
        r1.push_str(&format!("@r{i}\n{f1}\n+\n{}\n", "I".repeat(100)));
        r2.push_str(&format!("@r{i}\n{}\n+\n{}\n", revcomp(f2), "I".repeat(100)));
    }
    let p1 = dir.join("r1.fq");
    let p2 = dir.join("r2.fq");
    fs::write(&p1, &r1).unwrap();
    fs::write(&p2, &r2).unwrap();
    (p1, p2)
}

fn pogflow_args(r1: &Path, r2: &Path, outdir: &Path) -> trinity_cli::TrinityArgs {
    let argv: Vec<String> = [
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
    .collect();
    parse_args(&argv).unwrap()
}

#[test]
fn pogflow_e2e_produces_trinity_fasta() {
    let dir = tmpdir("e2e");
    let (r1, r2) = write_synthetic_pe(&dir);
    let outdir = dir.join("trinity_out");

    let final_fa = run_trinity(&pogflow_args(&r1, &r2, &outdir)).unwrap();
    let text = fs::read_to_string(&final_fa).unwrap();
    assert!(!text.is_empty() && text.contains('>'), "Trinity.fasta 空");
    assert!(
        text.lines()
            .any(|l| l.starts_with('>') && l.contains("_g1_i")),
        "缺 pogflow 命名 header: {text}"
    );
    assert!(outdir.join(".butterfly.ok").exists());
    let map = final_fa.with_file_name(format!(
        "{}.gene_trans_map",
        final_fa.file_name().unwrap().to_string_lossy()
    ));
    let map_text = fs::read_to_string(&map).unwrap();
    let n: usize = text.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(map_text.lines().filter(|l| !l.is_empty()).count(), n);
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn path_engine_pogflow_parse() {
    let argv: Vec<String> = [
        "--seqType",
        "fa",
        "--single",
        "s.fa",
        "--max_memory",
        "1G",
        "--output",
        "out",
        "--path-engine",
        "pogflow",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    let a = parse_args(&argv).unwrap();
    assert_eq!(a.path_engine, trinity_cli::args::PathEngine::PogFlow);
}
