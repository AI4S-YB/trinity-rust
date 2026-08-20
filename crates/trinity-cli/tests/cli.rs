//! P5-T2 集成测试: args 解析/校验、checkpoint 跳过语义、prep（fq PE 三形态
//! revcomp / fa SE / 多文件拼接 / both.fa 顺序与字节 / 断点跳过）。

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use trinity_cli::args::{parse_args, SeqType};
use trinity_cli::checkpoint::run_with_checkpoint;
use trinity_cli::prep::prep_reads;

fn argv(items: &[&str]) -> Vec<String> {
    items.iter().map(|s| s.to_string()).collect()
}

fn tmpdir(tag: &str) -> PathBuf {
    static N: AtomicU64 = AtomicU64::new(0);
    let p = std::env::temp_dir().join(format!(
        "trinity-cli-t2-{}-{}-{}",
        tag,
        std::process::id(),
        N.fetch_add(1, Ordering::SeqCst)
    ));
    fs::create_dir_all(&p).unwrap();
    p
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/p1/diginorm")
        .join(name)
}

fn base_args() -> Vec<String> {
    argv(&[
        "--seqType",
        "fq",
        "--left",
        "a.fq",
        "--right",
        "b.fq",
        "--max_memory",
        "1G",
        "--output",
        "out",
    ])
}

// ---------------------------------------------------------------- args

#[test]
fn args_defaults_and_names() {
    let a = parse_args(&base_args()).unwrap();
    assert_eq!(a.seq_type, SeqType::Fq);
    assert_eq!(a.left, vec![PathBuf::from("a.fq")]);
    assert_eq!(a.cpu, 2);
    assert_eq!(a.max_memory, 1024u64 * 1024 * 1024);
    assert_eq!(a.output, PathBuf::from("out"));
    assert_eq!(a.kmer_size, 25);
    assert_eq!(a.min_contig_length, 200);
    assert_eq!(a.min_kmer_count, 1);
    assert_eq!(a.group_pairs_distance, 500);
    assert_eq!(a.normalize_max_read_cov, 200);
    assert_eq!(a.inchworm_cpu, 2); // min(6, cpu)
    assert_eq!(a.max_reads_per_graph, 200000);
    assert_eq!(a.bfly_stack_mb, 256);
    assert!(!a.no_normalize_reads);
    assert!(!a.no_cleanup);
    assert!(a.ss_lib_type.is_none());
}

#[test]
fn args_explicit_values_and_comma_lists() {
    let a = parse_args(&argv(&[
        "--seqType",
        "fa",
        "--single",
        "s1.fa,s2.fa,s3.fa",
        "--SS_lib_type",
        "F",
        "--CPU",
        "8",
        "--max_memory",
        "10G",
        "--output",
        "o",
        "--KMER_SIZE",
        "31",
        "--min_contig_length",
        "300",
        "--min_kmer_count",
        "3",
        "--group_pairs_distance",
        "900",
        "--no_normalize_reads",
        "--normalize_max_read_cov",
        "50",
        "--inchworm_cpu",
        "4",
        "--bfly_stack_mb",
        "512",
        "--max_reads_per_graph",
        "1000",
        "--no_cleanup",
    ]))
    .unwrap();
    assert_eq!(a.seq_type, SeqType::Fa);
    assert_eq!(a.single.len(), 3);
    assert_eq!(a.single[2], PathBuf::from("s3.fa"));
    assert_eq!(a.cpu, 8);
    assert_eq!(a.max_memory, 10 * 1024 * 1024 * 1024);
    assert_eq!(a.kmer_size, 31);
    assert_eq!(a.min_contig_length, 300);
    assert_eq!(a.min_kmer_count, 3);
    assert_eq!(a.group_pairs_distance, 900);
    assert!(a.no_normalize_reads);
    assert_eq!(a.normalize_max_read_cov, 50);
    assert_eq!(a.inchworm_cpu, 4);
    assert_eq!(a.bfly_stack_mb, 512);
    assert_eq!(a.max_reads_per_graph, 1000);
    assert!(a.no_cleanup);
}

#[test]
fn args_inchworm_cpu_capped_by_cpu() {
    let a = parse_args(&argv(&[
        "--seqType",
        "fq",
        "--single",
        "s.fq",
        "--CPU",
        "3",
        "--max_memory",
        "2G",
        "--output",
        "o",
    ]))
    .unwrap();
    assert_eq!(a.inchworm_cpu, 3);
}

#[test]
fn args_max_memory_required_and_format() {
    let e = parse_args(&argv(&[
        "--seqType",
        "fq",
        "--single",
        "s.fq",
        "--output",
        "o",
    ]))
    .unwrap_err();
    assert!(e.usage);
    for bad in ["10", "10M", "G", "1xG"] {
        let e = parse_args(&argv(&[
            "--seqType",
            "fq",
            "--single",
            "s.fq",
            "--max_memory",
            bad,
            "--output",
            "o",
        ]))
        .unwrap_err();
        assert!(e.usage, "max_memory {bad} should be rejected");
    }
}

#[test]
fn args_ss_lib_type_validation() {
    // PE: 只 FR|RF
    for good in ["FR", "RF"] {
        let mut v = base_args();
        v.extend(["--SS_lib_type".into(), good.into()]);
        assert!(parse_args(&v).is_ok(), "PE {good} should pass");
    }
    for bad in ["F", "R", "XF"] {
        let mut v = base_args();
        v.extend(["--SS_lib_type".into(), bad.into()]);
        assert!(parse_args(&v).is_err(), "PE {bad} should fail");
    }
    // SE: 只 F|R
    for good in ["F", "R"] {
        assert!(parse_args(&argv(&[
            "--seqType",
            "fq",
            "--single",
            "s.fq",
            "--SS_lib_type",
            good,
            "--max_memory",
            "1G",
            "--output",
            "o",
        ]))
        .is_ok());
    }
    for bad in ["FR", "RF"] {
        assert!(parse_args(&argv(&[
            "--seqType",
            "fq",
            "--single",
            "s.fq",
            "--SS_lib_type",
            bad,
            "--max_memory",
            "1G",
            "--output",
            "o",
        ]))
        .is_err());
    }
}

#[test]
fn args_mixing_and_missing_reads_rejected() {
    // PE 与 SE 混用（原版 L1093-1095）
    let mut v = base_args();
    v.extend(["--single".into(), "s.fq".into()]);
    assert!(parse_args(&v).is_err());
    // 只有 left 无 right
    assert!(parse_args(&argv(&[
        "--seqType",
        "fq",
        "--left",
        "a.fq",
        "--max_memory",
        "1G",
        "--output",
        "o",
    ]))
    .is_err());
    // 无 reads
    assert!(parse_args(&argv(&[
        "--seqType",
        "fq",
        "--max_memory",
        "1G",
        "--output",
        "o",
    ]))
    .is_err());
    // seqType 必填/取值
    assert!(parse_args(&argv(&[
        "--single",
        "s.fq",
        "--max_memory",
        "1G",
        "--output",
        "o",
    ]))
    .is_err());
    assert!(parse_args(&argv(&[
        "--seqType",
        "cfa",
        "--single",
        "s.fa",
        "--max_memory",
        "1G",
        "--output",
        "o",
    ]))
    .is_err());
}

#[test]
fn args_ranges_and_unknown_option() {
    let mut v = base_args();
    v.extend(["--KMER_SIZE".into(), "33".into()]);
    assert!(parse_args(&v).is_err());
    let mut v = base_args();
    v.extend(["--min_contig_length".into(), "99".into()]);
    let e = parse_args(&v).unwrap_err();
    assert!(e.msg.contains("imposed threshold of 100"));
    let mut v = base_args();
    v.extend(["--trimmomatic".into()]);
    let e = parse_args(&v).unwrap_err();
    assert!(e.usage && e.msg.contains("do not understand option"));
}

// ---------------------------------------------------------------- checkpoint

#[test]
fn checkpoint_runs_and_skips() {
    let dir = tmpdir("ckpt");
    let ckpt = dir.join("step.ok");
    let mut count = 0usize;
    let did = run_with_checkpoint(&ckpt, "step one", || {
        count += 1;
        Ok(())
    })
    .unwrap();
    assert!(did && ckpt.exists() && count == 1);

    let did = run_with_checkpoint(&ckpt, "step one again", || {
        count += 1;
        Ok(())
    })
    .unwrap();
    assert!(!did && count == 1, "existing checkpoint must skip");
}

#[test]
fn checkpoint_failure_leaves_no_marker() {
    let dir = tmpdir("ckptfail");
    let ckpt = dir.join("bad.ok");
    let r = run_with_checkpoint(&ckpt, "bad step", || {
        Err(trinity_common::error::CommonError::Parse("boom".into()))
    });
    assert!(r.is_err() && !ckpt.exists());
}

// ---------------------------------------------------------------- prep

fn pe_args(ss: Option<&str>) -> trinity_cli::TrinityArgs {
    let mut v = argv(&["--seqType", "fq", "--max_memory", "1G", "--output", "out"]);
    let l = fixture("pe.l.fq");
    let r = fixture("pe.r.fq");
    v.extend(["--left".into(), l.to_str().unwrap().into()]);
    v.extend(["--right".into(), r.to_str().unwrap().into()]);
    if let Some(ss) = ss {
        v.extend(["--SS_lib_type".into(), ss.into()]);
    }
    parse_args(&v).unwrap()
}

fn read(p: &Path) -> String {
    fs::read_to_string(p).unwrap()
}

#[test]
fn prep_fq_pe_no_ss() {
    let out = tmpdir("prep-pe");
    let args = pe_args(None);
    let (target, n) = prep_reads(&args, &out).unwrap();
    assert_eq!(target, out.join("both.fa"));
    let l = read(&out.join("left.fa"));
    let r = read(&out.join("right.fa"));
    // seqtk -R 1/-R 2: 名字追加 /1、/2（fixtures 名已带 /1、/2 → 原样）
    assert!(l.starts_with(">pair1/1\n"));
    assert!(r.starts_with(">pair1/2\n"));
    // both.fa = left ++ right（先左后右、逐字节）
    let both = read(&target);
    assert_eq!(both, format!("{l}{r}"));
    assert!(out.join("left.fa.ok").exists());
    assert!(out.join("right.fa.ok").exists());
    assert!(out.join("both.fa.ok").exists());
    // read count = 两侧记录数之和
    let nl = l.lines().filter(|x| x.starts_with('>')).count() as u64;
    let nr = r.lines().filter(|x| x.starts_with('>')).count() as u64;
    assert_eq!(n, nl + nr);
}

#[test]
fn prep_fq_pe_ss_rf_revcomp_left() {
    let out = tmpdir("prep-rf");
    let args = pe_args(Some("RF"));
    let (_, _) = prep_reads(&args, &out).unwrap();
    // RF → left 端 'R' → revcomp; right 端 'F' → 原样
    let l = read(&out.join("left.fa"));
    assert!(
        l.contains("GTGAAATTAATACAATCGTCCCTTAGATGTATTCATCATT"),
        "RF left must be reverse-complemented:\n{l}"
    );
    let r = read(&out.join("right.fa"));
    assert!(
        r.contains("TTCTTCTCTGCGTAGAGGGTATGTTGACCTAAGCGAGGGCGAGTCCAACGTAGGAGGATG"),
        "RF right (F) must stay forward:\n{r}"
    );
    assert!(!l.contains("AATGATGAATACATCTAAGGGACGATTGTATTAATTTCAC"));
}

#[test]
fn prep_fq_pe_ss_fr_no_revcomp() {
    let out = tmpdir("prep-fr");
    let args = pe_args(Some("FR"));
    prep_reads(&args, &out).unwrap();
    let l = read(&out.join("left.fa"));
    assert!(l.contains("AATGATGAATACATCTAAGGGACGATTGTATTAATTTCAC"));
}

#[test]
fn prep_fq_multi_file_concat_order() {
    let out = tmpdir("prep-multi");
    let f1 = fixture("pe.l.fq");
    let f2 = fixture("ss.pe.l.fq");
    let v = argv(&[
        "--seqType",
        "fq",
        "--single",
        &format!("{},{}", f1.display(), f2.display()),
        "--max_memory",
        "1G",
        "--output",
        "out",
    ]);
    let args = parse_args(&v).unwrap();
    let (target, n) = prep_reads(&args, &out).unwrap();
    assert_eq!(target, out.join("single.fa"));
    let s = read(&target);
    assert!(s.starts_with(">pair1/1\n"));
    assert!(s.contains(">ssA1/1\n"), "second file must follow in order");
    let n_expected = s.lines().filter(|x| x.starts_with('>')).count() as u64;
    assert_eq!(n, n_expected);
}

#[test]
fn prep_fa_se_ss_r_revcomp() {
    let dir = tmpdir("prep-fa");
    let inp = dir.join("s.fa");
    fs::write(&inp, ">a\nACGTacgtNRY\n>b\nTTTT\n").unwrap();
    let v = argv(&[
        "--seqType",
        "fa",
        "--single",
        inp.to_str().unwrap(),
        "--SS_lib_type",
        "R",
        "--max_memory",
        "1G",
        "--output",
        "out",
    ]);
    let args = parse_args(&v).unwrap();
    let (target, n) = prep_reads(&args, &dir.join("out")).unwrap();
    let s = read(&target);
    // revcomp_fasta.pl 镜像（tr 表）+ 60 列折行
    assert!(s.contains(">a\n"));
    assert!(s.contains("RYNacgtACGT"), "revcomp via tr table:\n{s}");
    assert!(s.contains(">b\nAAAA"));
    assert_eq!(n, 2);
}

#[test]
fn prep_checkpoint_skip() {
    let out = tmpdir("prep-skip");
    // 预置 left.fa（哨兵内容）+ left.fa.ok → 该侧跳过
    fs::write(out.join("left.fa"), ">sentinel\nACGT\n").unwrap();
    fs::write(out.join("left.fa.ok"), b"").unwrap();
    fs::write(out.join("right.fa"), ">r\nTTTT\n").unwrap();
    fs::write(out.join("right.fa.ok"), b"").unwrap();
    fs::write(out.join("both.fa"), ">sentinel\nACGT\n>r\nTTTT\n").unwrap();
    fs::write(out.join("both.fa.ok"), b"").unwrap();
    let args = pe_args(None);
    let (target, n) = prep_reads(&args, &out).unwrap();
    assert_eq!(read(&target), ">sentinel\nACGT\n>r\nTTTT\n");
    assert_eq!(n, 2);
    // left fq 未被重新生成（哨兵保留）
    assert_eq!(read(&out.join("left.fa")), ">sentinel\nACGT\n");
}

#[test]
fn prep_both_fa_size_mismatch_rebuilds() {
    let out = tmpdir("prep-rebuild");
    // left/right 无 .ok → 重新生成; both.fa 预置为空且无 .ok → 重建为字节精确拼接
    fs::write(out.join("both.fa"), b"").unwrap();
    let args = pe_args(None);
    let (target, _) = prep_reads(&args, &out).unwrap();
    let l = read(&out.join("left.fa"));
    let r = read(&out.join("right.fa"));
    assert_eq!(read(&target), format!("{l}{r}"));
}
