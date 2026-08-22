//! M4 e2e: `--super_reads` K19 预组装段（合成 PE reads → super_reads.fa
//! header `>sr_N`）+ 关闭旗标无产物 + classic 不受影响。

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use trinity_cli::args::parse_args;
use trinity_cli::orchestrate::run_trinity;

fn tmpdir(tag: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!(
        "trinity-cli-superreads-{}-{}",
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

/// 800nt 参考 × 40 对, 150nt reads（FR 插入 250）——足够深的 K19 覆盖
/// 让 inchworm 预组装产出 ≥200bp contig。
fn write_synthetic_pe(dir: &Path) -> (PathBuf, PathBuf) {
    let mut x: u64 = 0x9E3779B97F4A7C15;
    let reference: String = (0..800).map(|_| xorshift_nt(&mut x)).collect();
    let rb = reference.as_bytes();
    let mut r1 = String::new();
    let mut r2 = String::new();
    for i in 0..40usize {
        let s = i * 8;
        let f1 = std::str::from_utf8(&rb[s..s + 150]).unwrap();
        let f2 = std::str::from_utf8(&rb[s + 100..s + 250]).unwrap();
        r1.push_str(&format!("@r{i}\n{f1}\n+\n{}\n", "I".repeat(150)));
        r2.push_str(&format!("@r{i}\n{}\n+\n{}\n", revcomp(f2), "I".repeat(150)));
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
    ]
    .iter()
    .map(|s| s.to_string())
    .collect()
}

fn with_flag(mut argv: Vec<String>, flag: &str) -> trinity_cli::TrinityArgs {
    argv.push(flag.to_string());
    parse_args(&argv).unwrap()
}

#[test]
fn super_reads_e2e_on_off_and_classic_unaffected() {
    let dir = tmpdir("e2e");
    let (r1, r2) = write_synthetic_pe(&dir);

    // ---- on: super_reads.fa 存在, header 含 >sr_ ----
    let out_on = dir.join("trinity_on");
    let final_on = run_trinity(&with_flag(base_args(&r1, &r2, &out_on), "--super_reads")).unwrap();
    let sr = out_on.join("super_reads.fa");
    assert!(sr.exists(), "super_reads.fa 应存在");
    assert!(
        out_on.join(".k19_preasm.ok").exists(),
        ".k19_preasm.ok 缺失"
    );
    let sr_text = fs::read_to_string(&sr).unwrap();
    assert!(sr_text.contains(">sr_"), "header 应含 >sr_: {sr_text}");
    let text_on = fs::read_to_string(&final_on).unwrap();
    assert!(!text_on.is_empty() && text_on.contains('>'));

    // ---- off: 无 super_reads.fa / .k19_preasm.ok ----
    let out_off = dir.join("trinity_off");
    run_trinity(&parse_args(&base_args(&r1, &r2, &out_off)).unwrap()).unwrap();
    assert!(!out_off.join("super_reads.fa").exists());
    assert!(!out_off.join(".k19_preasm.ok").exists());
}
