//! M1-T3 e2e: `--recursive` 递归分区模式全流程（合成 PE reads → 分区池 →
//! Trinity.fasta）+ 断点续跑（同 outdir 重跑秒回且输出不变）。

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use trinity_cli::args::parse_args;
use trinity_cli::orchestrate::run_trinity;

fn tmpdir(tag: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!(
        "trinity-cli-recursive-{}-{}",
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

/// 合成 PE reads（T2 recursive.rs 同手法）: 600nt 参考（xorshift）× 25 对,
/// R1 = ref[s..s+100], R2 = revcomp(ref[s+150..s+250])（FR 文库插入 250）。
/// 写临时 fq.gz（gzip level 1）。
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
    let p1 = dir.join("r1.fq.gz");
    let p2 = dir.join("r2.fq.gz");
    let mut w1 =
        flate2::write::GzEncoder::new(fs::File::create(&p1).unwrap(), flate2::Compression::fast());
    w1.write_all(r1.as_bytes()).unwrap();
    w1.finish().unwrap();
    let mut w2 =
        flate2::write::GzEncoder::new(fs::File::create(&p2).unwrap(), flate2::Compression::fast());
    w2.write_all(r2.as_bytes()).unwrap();
    w2.finish().unwrap();
    (p1, p2)
}

fn recursive_args(r1: &Path, r2: &Path, outdir: &Path) -> trinity_cli::TrinityArgs {
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
        "--recursive",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    parse_args(&argv).unwrap()
}

#[test]
fn recursive_e2e_with_checkpoint_resume() {
    let dir = tmpdir("e2e");
    let (r1, r2) = write_synthetic_pe(&dir);
    let outdir = dir.join("trinity_out");

    let final_fa = run_trinity(&recursive_args(&r1, &r2, &outdir)).unwrap();
    // 断言: Trinity.fasta 存在非空
    let text = fs::read_to_string(&final_fa).unwrap();
    assert!(!text.is_empty() && text.contains('>'), "Trinity.fasta 空");
    // 默认 scratch: 中间产物与 checkpoint 在 default_scratch_dir 内
    let work = trinity_cli::orchestrate::default_scratch_dir(&outdir);
    // 断言: 分区产物 ≥1
    let minis: Vec<_> = fs::read_dir(work.join("read_partitions"))
        .unwrap()
        .filter_map(|e| {
            let p = e.unwrap().path();
            (p.extension().map(|x| x == "fa").unwrap_or(false)
                && p.file_name()?.to_string_lossy().ends_with(".mini.fa"))
            .then_some(p)
        })
        .collect();
    assert!(!minis.is_empty(), "read_partitions/*.mini.fa 应 ≥1");
    // 断言: 断点标记存在（scratch 内）
    assert!(work.join(".recursive.ok").exists(), ".recursive.ok 缺失");
    // classic 的 butterfly 路径不应触发; output 目录无中间目录
    assert!(!work.join(".butterfly.ok").exists());
    assert!(!outdir.join("read_partitions").exists());
    // gene_trans_map 与 classic 路径一致
    let map = final_fa.with_file_name(format!(
        "{}.gene_trans_map",
        final_fa.file_name().unwrap().to_string_lossy()
    ));
    let map_text = fs::read_to_string(&map).unwrap();
    let n: usize = text.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(map_text.lines().filter(|l| !l.is_empty()).count(), n);

    // ---- 断点续跑: 同 outdir 重跑 → 秒回（<5s）且 Trinity.fasta 内容不变 ----
    let before = fs::read_to_string(&final_fa).unwrap();
    let mini_mtimes: Vec<_> = minis
        .iter()
        .map(|p| fs::metadata(p).unwrap().modified().unwrap())
        .collect();
    let t0 = Instant::now();
    let final_fa2 = run_trinity(&recursive_args(&r1, &r2, &outdir)).unwrap();
    let elapsed = t0.elapsed();
    assert!(
        elapsed.as_secs_f64() < 5.0,
        "断点续跑应 <5s, 实际 {elapsed:?}"
    );
    assert_eq!(final_fa2, final_fa);
    assert_eq!(fs::read_to_string(&final_fa).unwrap(), before, "输出应不变");
    // mini.fa 未被重算（mtime 不变）
    for (p, m) in minis.iter().zip(&mini_mtimes) {
        assert_eq!(
            fs::metadata(p).unwrap().modified().unwrap(),
            *m,
            "mini.fa {p:?} 不应重算"
        );
    }
}
