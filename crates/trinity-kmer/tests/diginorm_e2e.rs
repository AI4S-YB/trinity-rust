//! DigiNorm together 模式端到端（fixtures/p1/diginorm，手推名单锁定）。
//!
//! fixture 手推（K=25、DS、dump -L 2、maxC200/minC1/maxCV10000）:
//! - left.fa: pair1/1=S_A、pair2/1=S_A（完全同序列）、pair3/1=S_N（29 位 N）、
//!   pair4/1=20bp、pair5/1、pair6/1 唯一序列;
//! - right.fa: pair1/2..pair3/2=S_B（3 份）、pair4/2=20bp、pair5/2、pair6/2 唯一。
//! - 计数表（≥2）: S_A 的 25-mer = 2（N 不含窗）或 3（S_N 贡献 +1）; S_B = 3。
//! - left median: pair1/2/1→2、pair3/1→1、pair4/1→0（<25bp 空向量）、pair5/6→1;
//!   right median: pair1-3/2→3、pair4/2→0、pair5/6→1。
//! - PE merge median: pair1/2→2.5、pair3→2.0、pair4→0.0、pair5/6→1.0。
//! - select: ratio=200/median ≥ 80 ≥ 1 恒保留（不依赖 drand48）; pair4 median 0 < min_cov 1 → 丢。
//!   ⇒ 名单 = [pair1, pair2, pair3, pair5, pair6]（字节序），输出 fq 每侧 5×4=20 行。

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::process::Command;

use trinity_kmer::counter::{CountMap, KmerCountTable};
use trinity_kmer::coverage_stats::{coverage_stats_rows, write_stats_tsv};
use trinity_kmer::diginorm::{run, DigiNormParams, ReadsInput};
use trinity_kmer::nbkc::{merge_pairs, StatsRow};
use trinity_kmer::read_names::{fq_records_to_fa, ReadType};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/p1/diginorm")
        .join(name)
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("trinity_kmer_diginorm_{}", tag));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// fq/fa 文件的 (首 token, 整条记录字节) 列表（fq 4 行/条、fa 2 行/条——fixture 均单行序列）。
fn fq_records(data: &[u8]) -> Vec<(String, Vec<u8>)> {
    let mut out = Vec::new();
    let lines: Vec<&[u8]> = split_lines(data);
    let per = if lines.first().is_some_and(|l| l.starts_with(b"@")) {
        4
    } else {
        2
    };
    for chunk in lines.chunks(per) {
        let header = chunk[0].strip_suffix(b"\n").unwrap_or(chunk[0]);
        let body = header
            .strip_prefix(b"@")
            .or_else(|| header.strip_prefix(b">"))
            .unwrap_or(header);
        let token = body
            .split(|&b| b == b' ' || b == b'\t')
            .next()
            .unwrap()
            .to_vec();
        out.push((String::from_utf8(token).unwrap(), chunk.concat()));
    }
    out
}

fn split_lines(data: &[u8]) -> Vec<&[u8]> {
    let mut out = Vec::new();
    let mut start = 0;
    while start < data.len() {
        let end = data[start..]
            .iter()
            .position(|&b| b == b'\n')
            .map_or(data.len(), |i| start + i + 1);
        out.push(&data[start..end]);
        start = end;
    }
    // 末行无换行时也计入; 尾部空行剔除
    while out.last().is_some_and(|l| l == b"\n" || l.is_empty()) {
        out.pop();
    }
    out
}

const KEPT: [&str; 5] = ["pair1", "pair2", "pair3", "pair5", "pair6"];

fn kept_fq_bytes(fixture_path: &Path) -> Vec<u8> {
    let data = std::fs::read(fixture_path).unwrap();
    let kept: HashSet<&str> = KEPT.iter().copied().collect();
    let mut out = Vec::new();
    for (token, rec) in fq_records(&data) {
        let core = token
            .strip_suffix("/1")
            .or_else(|| token.strip_suffix("/2"))
            .unwrap_or(&token);
        if kept.contains(core) {
            out.extend_from_slice(&rec);
        }
    }
    out
}

#[test]
fn paired_fq_end_to_end_hand_derived() {
    let out_dir = tmpdir("pe_fq");
    let params = DigiNormParams::default();
    let reads = ReadsInput::Paired(vec![fixture("pe.l.fq")], vec![fixture("pe.r.fq")]);
    let outs = run(&params, &reads, &out_dir).unwrap();

    // 输出命名（镜像 L398-410）
    let left = outs.left.unwrap();
    assert_eq!(
        left.file_name().unwrap().to_str().unwrap(),
        "pe.l.fq.normalized_K25_maxC200_minC1_maxCV10000.fq"
    );
    let right = outs.right.unwrap();
    assert!(right
        .file_name()
        .unwrap()
        .to_str()
        .unwrap()
        .ends_with("pe.r.fq.normalized_K25_maxC200_minC1_maxCV10000.fq"));
    assert!(outs.single.is_none());

    // 抽取内容 = 原始记录逐字节（fq 4 行透传）; 短 read pair4 不保留
    let left_out = std::fs::read(&left).unwrap();
    assert_eq!(left_out, kept_fq_bytes(&fixture("pe.l.fq")));
    let right_out = std::fs::read(&right).unwrap();
    assert_eq!(right_out, kept_fq_bytes(&fixture("pe.r.fq")));
    // 行数 = 名单数 × 4
    assert_eq!(split_lines(&left_out).len(), 5 * 4);
    assert_eq!(split_lines(&right_out).len(), 5 * 4);
    assert!(!String::from_utf8_lossy(&left_out).contains("pair4"));
    assert!(!String::from_utf8_lossy(&right_out).contains("pair4"));

    // 便捷名双写（原版 symlink 的位置与内容）
    let conv_left = std::fs::read(out_dir.join("left.norm.fq")).unwrap();
    assert_eq!(conv_left, left_out);
    let conv_right = std::fs::read(out_dir.join("right.norm.fq")).unwrap();
    assert_eq!(conv_right, right_out);
}

#[test]
fn single_fq_end_to_end() {
    let out_dir = tmpdir("se_fq");
    let params = DigiNormParams::default();
    let outs = run(
        &params,
        &ReadsInput::Single(vec![fixture("pe.l.fq")]),
        &out_dir,
    )
    .unwrap();
    let single = outs.single.unwrap();
    assert!(single
        .file_name()
        .unwrap()
        .to_str()
        .unwrap()
        .starts_with("pe.l.fq.normalized_K25_"));
    // 单端不经 merge: pair3/1 的 median 是整数 1（不是 merge 后的 2.0），仍 ≥ min_cov → 保留;
    // 名单与 PE 相同（本 fixture 两侧过滤行为一致）
    let out = std::fs::read(&single).unwrap();
    assert_eq!(out, kept_fq_bytes(&fixture("pe.l.fq")));
    assert_eq!(std::fs::read(out_dir.join("single.norm.fq")).unwrap(), out);
    assert!(outs.left.is_none() && outs.right.is_none());
}

/// 多文件（逗号列表展开）：按列表序拼接 = 单文件合并流的同一归一化结果
/// （原版 prep_list_of_seqs 把每侧文件列表串成一个 single.fa，全部一起
/// 归一化）；输出 basename = 第一个文件名 + `_ext_all_reads`（L399）。
#[test]
fn single_multi_file_matches_concat() {
    let src = std::fs::read(fixture("pe.l.fq")).unwrap();
    let lines = split_lines(&src);
    // 6 条记录 → 前 3 条 + 后 3 条两个文件
    assert_eq!(lines.len(), 24);
    let dir = tmpdir("se_multi_in");
    let a = dir.join("readsA.fq");
    let b = dir.join("readsB.fq");
    std::fs::write(&a, lines[..12].concat()).unwrap();
    std::fs::write(&b, lines[12..].concat()).unwrap();

    // 单文件（两半拼回）基线
    let base_dir = tmpdir("se_multi_base");
    let whole = dir.join("whole.fq");
    std::fs::write(&whole, src).unwrap();
    let base = run(
        &DigiNormParams::default(),
        &ReadsInput::Single(vec![whole]),
        &base_dir,
    )
    .unwrap()
    .single
    .unwrap();

    let out_dir = tmpdir("se_multi");
    let outs = run(
        &DigiNormParams::default(),
        &ReadsInput::Single(vec![a, b]),
        &out_dir,
    )
    .unwrap()
    .single
    .unwrap();
    let name = outs.file_name().unwrap().to_str().unwrap();
    assert!(
        name.starts_with("readsA.fq_ext_all_reads.normalized_K25_"),
        "多文件 basename 应为 first + _ext_all_reads: {name}"
    );
    assert_eq!(std::fs::read(&outs).unwrap(), std::fs::read(&base).unwrap());
    assert_eq!(
        std::fs::read(out_dir.join("single.norm.fq")).unwrap(),
        std::fs::read(&base).unwrap()
    );
}

#[test]
fn paired_fa_end_to_end() {
    let out_dir = tmpdir("pe_fa");
    let params = DigiNormParams::default();
    let outs = run(
        &params,
        &ReadsInput::Paired(vec![fixture("pe.l.fa")], vec![fixture("pe.r.fa")]),
        &out_dir,
    )
    .unwrap();
    let left = std::fs::read(outs.left.unwrap()).unwrap();
    // fa 抽取: ">header\nseq\n" 重排（单行序列原样）; pair4 短 read 不保留
    let expect = kept_fq_bytes(&fixture("pe.l.fa"));
    assert_eq!(left, expect);
    let text = String::from_utf8(left.clone()).unwrap();
    assert!(text.contains(">pair1/1\n"));
    assert!(!text.contains("pair4"));
    assert!(std::fs::read(out_dir.join("left.norm.fa")).unwrap() == left);
}

/// SS_lib_type RF: left 侧 fq→fa 走 -r（revcomp），名单不变（本 fixture 的计数
/// 在 SS 计数下无跨链碰撞，left.fa revcomp 后计数表同步）——锁定 SS 路径可运行且
/// 选择一致。
#[test]
fn paired_fq_ss_rf() {
    let out_dir = tmpdir("pe_ss");
    let params = DigiNormParams {
        ss_revcomp_left: true,
        ds: false,
        ..Default::default()
    };
    let outs = run(
        &params,
        &ReadsInput::Paired(vec![fixture("pe.l.fq")], vec![fixture("pe.r.fq")]),
        &out_dir,
    )
    .unwrap();
    let left_out = std::fs::read(outs.left.unwrap()).unwrap();
    // 抽取回原始 fq（revcomp 只影响计数用的 left.fa，不影响输出格式）
    assert_eq!(left_out, kept_fq_bytes(&fixture("pe.l.fq")));
}

// ---------------------------------------------------------------- SS F oracle（原版管线 golden）

/// SS F 互补链 fixture（左 X / 右 revcomp(X)，典型 dUTP 文库）。K=25、maxC200/minC3，
/// 计数(literal)→-L2→**canonical 合并**后:
/// - ssA(×2/侧): canonical 4 → merged median 4.0 ≥ 3 → **保留**（buggy 无合并: 2.0 < 3 → 丢）;
/// - ssB(左1/右2): X_B literal 1 被 dump -L 丢弃 → canonical 2（仅右链贡献）→ merged
///   2.0 < 3 → 丢（**锁定先过滤后合并**: 错误的先合并得 3.0 → 保留）;
/// - ssC(×3/侧): canonical 6 → merged 6.0 → 保留（两版一致，对照）。
///
/// golden 由原版管线生成（见 [`ss_pe_f_complementary_oracle`]）。
fn ss_fixture_params() -> DigiNormParams {
    DigiNormParams {
        ds: false, // SS_lib_type F: 无 revcomp、jellyfish 无 --canonical
        min_cov: 3.0,
        ..Default::default()
    }
}

/// 原版管线对拍（trinityrnaseq-v2.15.2 + /usr/bin/perl，--SS_lib_type F --pairs_together
/// --max_cov 200 --min_cov 3）: ss.pe.{l,r}.norm.golden.fq 即其 left/right.norm.fq 产物。
/// 本测试断言两侧输出与 golden **逐字节一致**。
#[test]
fn ss_pe_f_complementary_oracle() {
    let out_dir = tmpdir("ss_pe");
    let params = ss_fixture_params();
    let outs = run(
        &params,
        &ReadsInput::Paired(vec![fixture("ss.pe.l.fq")], vec![fixture("ss.pe.r.fq")]),
        &out_dir,
    )
    .unwrap();
    let left = std::fs::read(outs.left.unwrap()).unwrap();
    let right = std::fs::read(outs.right.unwrap()).unwrap();
    assert_eq!(
        left,
        std::fs::read(fixture("ss.pe.left.norm.golden.fq")).unwrap()
    );
    assert_eq!(
        right,
        std::fs::read(fixture("ss.pe.right.norm.golden.fq")).unwrap()
    );
    // 名单语义（golden 的内容）: ssA×2 + ssC×3 保留、ssB 丢弃——buggy 版丢 ssA（锁定修复）
    let text = String::from_utf8_lossy(&left).into_owned();
    let names: Vec<&str> = text
        .lines()
        .step_by(4)
        .map(|l| l.trim_start_matches('@'))
        .collect();
    assert_eq!(names, ["ssA1/1", "ssA2/1", "ssC1/1", "ssC2/1", "ssC3/1"]);
}

/// 单端混链（同一文件含 X 与 revcomp(X)）原版对拍（同上参数，--single）:
/// ssA(X×2+rc×2) median 4 → 保留（buggy: 2 → 丢）、ssB(X×1+rc×2) median 2 → 丢、
/// ssC(×3+×3) median 6 → 保留。输出与 golden 逐字节一致。
#[test]
fn ss_single_mixed_strands_oracle() {
    let out_dir = tmpdir("ss_se");
    let params = ss_fixture_params();
    let outs = run(
        &params,
        &ReadsInput::Single(vec![fixture("ss.single.fq")]),
        &out_dir,
    )
    .unwrap();
    let out = std::fs::read(outs.single.unwrap()).unwrap();
    assert_eq!(
        out,
        std::fs::read(fixture("ss.single.norm.golden.fq")).unwrap()
    );
    let text = String::from_utf8_lossy(&out).into_owned();
    let names: Vec<&str> = text
        .lines()
        .step_by(4)
        .map(|l| l.trim_start_matches('@'))
        .collect();
    assert_eq!(names.len(), 10); // ssA×4 + ssC×6
    assert!(names
        .iter()
        .all(|n| n.starts_with("ssA") || n.starts_with("ssC")));
}

// ---------------------------------------------------------------- 中间产物 oracle（原版管线产物 golden）

/// 中间产物 golden（由原版 insilico_read_normalization.pl 在本 fixture 上的运行产物切列生成，
/// tid 列因原版多线程为 thread:0/1 而本移植恒 thread:0，golden 切去该列）:
/// - pe.{l,r}.stats.sort.golden.tsv: left/right.fa.K25.stats.sort 的前 4 列
/// - pairs.stats.golden.tsv: pairs.K25.stats 的 acc + 合成 3 列（cut -f1,10,11,12）
///
/// 本测试用公开 API 重走 diginorm 内部路径（含双次文本往返），逐字节对照 golden。
#[test]
fn intermediates_match_original_pipeline_golden() {
    let k = 25;
    let l_fq = std::fs::read(fixture("pe.l.fq")).unwrap();
    let r_fq = std::fs::read(fixture("pe.r.fq")).unwrap();
    let left_fa = fq_records_to_fa(&l_fq, ReadType::R1, false).unwrap();
    let right_fa = fq_records_to_fa(&r_fq, ReadType::R2, false).unwrap();
    let mut both = left_fa.clone();
    both.extend_from_slice(&right_fa);
    // dump -L 2 语义
    let table: CountMap = KmerCountTable::count_fasta_data(&both, k, true)
        .into_iter()
        .filter(|&(_, c)| c >= 2)
        .collect();

    let render_sorted = |fa: &[u8]| -> Vec<String> {
        let rows = coverage_stats_rows(fa, &table, k, true).unwrap();
        let mut buf = Vec::new();
        write_stats_tsv(&mut buf, &rows).unwrap();
        let mut full: Vec<String> = String::from_utf8(buf)
            .unwrap()
            .lines()
            .skip(1) // 表头
            .map(str::to_string)
            .collect();
        // 按 acc 字节序（镜像原版 sort -k1,1）
        full.sort_by(|a, b| {
            a.split('\t')
                .next()
                .unwrap()
                .cmp(b.split('\t').next().unwrap())
        });
        full
    };

    let parse_row = |line: &str| -> StatsRow {
        // 与 diginorm::parse_num 同一 perl 数值化语义（-0 归一 +0）
        let num = |s: &str| -> f64 {
            let v: f64 = s.parse().unwrap();
            if v == 0.0 {
                0.0
            } else {
                v
            }
        };
        let f: Vec<&str> = line.split('\t').collect();
        StatsRow {
            acc: f[0].to_string(),
            median: num(f[1]),
            mean: num(f[2]),
            stdev: num(f[3]),
        }
    };

    let left_lines = render_sorted(&left_fa);
    let right_lines = render_sorted(&right_fa);
    // golden 对照（切掉 tid 列 = 每行 split('\t') 前 4 列）
    let cut4 = |l: &str| -> String { l.split('\t').take(4).collect::<Vec<_>>().join("\t") };
    for (name, lines) in [
        ("pe.l.stats.sort.golden.tsv", &left_lines),
        ("pe.r.stats.sort.golden.tsv", &right_lines),
    ] {
        let golden = std::fs::read_to_string(fixture(name)).unwrap();
        let golden_rows: Vec<String> = golden.lines().skip(1).map(cut4).collect();
        let ours: Vec<String> = lines.iter().map(|l| cut4(l)).collect();
        assert_eq!(ours, golden_rows, "{name}");
    }

    // PE merge（第一次文本往返: 行文本 parse → merge → %.1f 列再 parse，均由库完成）
    let lrows: Vec<StatsRow> = left_lines.iter().map(|l| parse_row(l)).collect();
    let rrows: Vec<StatsRow> = right_lines.iter().map(|l| parse_row(l)).collect();
    let merged = merge_pairs(&lrows, &rrows);
    let rendered: Vec<String> = merged
        .iter()
        .map(|m| format!("{}\t{}\t{}\t{}", m.core, m.median, m.mean, m.stdev))
        .collect();
    let golden = std::fs::read_to_string(fixture("pairs.stats.golden.tsv")).unwrap();
    let golden_rows: Vec<String> = golden.lines().skip(1).map(str::to_string).collect();
    assert_eq!(rendered, golden_rows);
    // 第二次文本往返: %.1f 列 parse 回 f64（nbkc 输入侧），median 值域 sanity
    let medians: Vec<f64> = merged
        .iter()
        .map(|m| m.median.parse::<f64>().unwrap())
        .collect();
    assert_eq!(medians, vec![2.5, 2.5, 2.0, 0.0, 1.0, 1.0]);
}

// ---------------------------------------------------------------- CLI

fn bin_path() -> &'static str {
    env!("CARGO_BIN_EXE_trinity-kmer")
}

fn run_cli(args: &[&str]) -> (i32, String, String) {
    let out = Command::new(bin_path())
        .args(args)
        .output()
        .expect("spawn trinity-kmer");
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

/// T2 审查 Important（CLI 校验，count/diginorm 一并）: 未知 flag / 缺值 / -K 越界 → exit 2。
#[test]
fn cli_validation_exit_2() {
    let reads = fixture("pe.l.fa");
    let reads_str = reads.to_str().unwrap();
    // 未知 flag
    let (code, _, err) = run_cli(&["count", "--reads", reads_str, "--bogus"]);
    assert_eq!(code, 2, "unknown flag: {err}");
    assert!(err.contains("bogus"), "{err}");
    // 值 flag 缺值
    let (code, _, err) = run_cli(&["count", "--reads"]);
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("expects a value"), "{err}");
    // -K 非法
    let (code, _, err) = run_cli(&["count", "--reads", reads_str, "-K", "abc"]);
    assert_eq!(code, 2, "{err}");
    // -K 越界（>32）
    let (code, _, err) = run_cli(&["count", "--reads", reads_str, "-K", "33"]);
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("32"), "{err}");
    let (code, _, _) = run_cli(&["count", "--reads", reads_str, "-K", "0"]);
    assert_eq!(code, 2);
    // coverage-stats 缺 --kmers
    let (code, _, err) = run_cli(&["coverage-stats", "--reads", reads_str]);
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("--kmers required"), "{err}");
    // diginorm: left 无 right
    let (code, _, err) = run_cli(&["diginorm", "--left", reads_str, "-o", "/tmp/x_dn"]);
    assert_eq!(code, 2, "{err}");
    // diginorm: SS_lib_type 非法（原版消息逐字风格）
    let (code, _, err) = run_cli(&[
        "diginorm",
        "--single",
        reads_str,
        "-o",
        "/tmp/x_dn",
        "--SS_lib_type",
        "RR",
    ]);
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("unrecognized SS_lib_type"), "{err}");
    // diginorm: max_cov < 2（原版 L231）
    let (code, _, err) = run_cli(&[
        "diginorm",
        "--single",
        reads_str,
        "-o",
        "/tmp/x_dn",
        "--max_cov",
        "1",
    ]);
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("at least 2"), "{err}");
}

/// coverage-stats CLI: count 产 dump → stats median 列与手推一致（2,2,1,0,1,1）。
#[test]
fn cli_coverage_stats_medians_hand_derived() {
    let out_dir = tmpdir("cli_stats");
    let reads = fixture("pe.l.fa");
    let dump = out_dir.join("kmers.fa");
    let (code, _, err) = run_cli(&[
        "count",
        "--reads",
        reads.to_str().unwrap(),
        "--canonical",
        "--min-count",
        "2",
        "-K",
        "25",
        "-o",
        dump.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "{err}");
    let (code, tsv, err) = run_cli(&[
        "coverage-stats",
        "--reads",
        reads.to_str().unwrap(),
        "--kmers",
        dump.to_str().unwrap(),
        "-K",
        "25",
    ]);
    assert_eq!(code, 0, "{err}");
    let lines: Vec<&str> = tsv.lines().collect();
    assert_eq!(lines[0], "acc\tmedian_cov\tmean_cov\tstdev\ttid");
    let expect = [
        ("pair1/1", "2"),
        ("pair2/1", "2"),
        ("pair3/1", "1"),
        ("pair4/1", "0"),
        ("pair5/1", "1"),
        ("pair6/1", "1"),
    ];
    for (line, (acc, med)) in lines.iter().skip(1).zip(expect) {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f[0], acc, "{line}");
        assert_eq!(f[1], med, "{line}");
    }
    assert_eq!(lines.len(), 7);
}

/// diginorm CLI 快乐路径: exit 0 + 输出文件与库层 run() 一致。
#[test]
fn cli_diginorm_happy_path() {
    let out_dir = tmpdir("cli_dn");
    let (code, stdout, err) = run_cli(&[
        "diginorm",
        "--left",
        fixture("pe.l.fq").to_str().unwrap(),
        "--right",
        fixture("pe.r.fq").to_str().unwrap(),
        "-o",
        out_dir.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "stdout={stdout} err={err}");
    let left = out_dir.join("pe.l.fq.normalized_K25_maxC200_minC1_maxCV10000.fq");
    assert!(left.exists());
    assert_eq!(
        std::fs::read(&left).unwrap(),
        kept_fq_bytes(&fixture("pe.l.fq"))
    );
    assert!(out_dir.join("left.norm.fq").exists());
    assert!(out_dir.join("right.norm.fq").exists());
}
