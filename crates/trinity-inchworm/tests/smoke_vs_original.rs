//! P2-T5 冒烟对拍——原版 inchworm 二进制（trinityrnaseq-v2.15.2）锁定黄金。
//!
//! fixture 生成（原版二进制一次性运行，产物入库）:
//! ```sh
//! TRINITY_SRC=.../trinityrnaseq-v2.15.2
//! $TRINITY_SRC/Inchworm/bin/inchworm --kmers fixtures/p1/smoke.kmers.fa \
//!     --run_inchworm -K 25 --monitor 2 --DS --num_threads 1 \
//!     > fixtures/p2/smoke.orig.fa 2> smoke.orig.log
//! grep 'SEED kmer' smoke.orig.log | sed 's/SEED kmer: //; s/, count: /\t/' \
//!     > fixtures/p2/smoke.seed_order.orig.tsv
//! ```
//!
//! 已证结论（2026-08-16，对拍三组）:
//! 1. **核心逐字节一致**: 以原版 --monitor 2 抓到的种子序重放本库
//!    （populate → 默认剪枝 → compute_sequence_assemblies_from_seeds），三组输入
//!    （本 smoke fixture、sample_data 全量 jellyfish dump 288,470 kmers、
//!    sample_data --reads 模式）全部 **BYTE-MATCH**——贪心延伸/tie 打破/glibc
//!    rand/记录清零/去重/输出格式全链路位级复刻。
//! 2. 默认种子序（本移植: count 降序 + 平局 kmer 值降序的全序）与原版
//!    （__gnu_cxx::hash_map 迭代序 + libstdc++ std::sort 平局不稳定序）不同:
//!    - smoke fixture: 6 vs 6 contig，header 多重集（去 aN）相等，**rc 不变多重集
//!      相等**（DS 语义下 contig 与其 revcomp 同义——3 条仅链方向不同）
//!    - sample_data: 822 vs 822，rc 不变多重集差 15/822（kmers 模式）与
//!      18/822（reads 模式）——分支点平局的划分差异，已由 (1) 归因于种子序
//!      （同种子序即逐字节一致）。

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;

use trinity_common::kmer::kmer_to_intval;
use trinity_inchworm::irke::{
    compute_sequence_assemblies_from_seeds, populate_from_kmers, prune_some_kmers, AssemblyParams,
    IrkeParams, Monitor,
};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/p2")
        .join(name)
}

fn p1_fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/p1")
        .join(name)
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("trinity_iworm_smoke_{tag}"));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// FASTA → Vec<(header, sequence)>（单行/折行序列均可）。
fn read_fasta(path: &Path) -> Vec<(String, String)> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut recs = Vec::new();
    let mut header: Option<String> = None;
    let mut seq = String::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if let Some(h) = header.take() {
                recs.push((h, std::mem::take(&mut seq)));
            }
            header = Some(h.to_string());
        } else if !line.is_empty() {
            seq.push_str(line);
        }
    }
    if let Some(h) = header {
        recs.push((h, seq));
    }
    recs
}

fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'G' => b'C',
            b'A' => b'T',
            b'T' => b'A',
            _ => b'G',
        })
        .map(|b| b as char)
        .collect()
}

/// 多重集（std Counter 未稳定——BTreeMap 计数）。
fn multiset<T: Ord + Clone>(items: impl Iterator<Item = T>) -> BTreeMap<T, usize> {
    let mut m = BTreeMap::new();
    for i in items {
        *m.entry(i).or_insert(0) += 1;
    }
    m
}

/// rc 不变 key: (s, revcomp(s)) 的字典序小者（DS 语义下二者同义分子）。
fn rc_key(s: &str) -> String {
    std::cmp::min(s.to_string(), revcomp(s))
}

// ---------------------------------------------------------------- 核心逐字节重放

/// **贪心核心位级锁定**: 原版种子序（--monitor 2 抓取）+ 同输入 → 逐字节相等。
/// 该测试通过即证明: 与原版的全部残余分歧只能来自种子平局序（原版哈希迭代序
/// + 不稳定 sort），而非贪心延伸/tie 打破/rand/清零/去重/格式化逻辑。
#[test]
fn replay_with_original_seed_order_is_byte_identical() {
    let data = std::fs::read(p1_fixture("smoke.kmers.fa")).unwrap();
    let (mut counter, _parsed) = populate_from_kmers(&data, 25, true).unwrap();
    // CLI 默认剪枝路径（--kmers 无 --no_prune_error_kmers 时 prune_error_kmers=true）
    prune_some_kmers(&mut counter, 1, 0.0, true, 0.005);

    let seeds_text = std::fs::read_to_string(fixture("smoke.seed_order.orig.tsv")).unwrap();
    let seeds: Vec<(u64, u32)> = seeds_text
        .lines()
        .map(|l| {
            let (kmer, count) = l.split_once('\t').unwrap();
            (
                kmer_to_intval(kmer.as_bytes()).unwrap(),
                count.parse().unwrap(),
            )
        })
        .collect();
    assert!(!seeds.is_empty());

    let mut out = Vec::new();
    let n = compute_sequence_assemblies_from_seeds(
        &mut counter,
        &IrkeParams::default(),
        &AssemblyParams::default(),
        &Monitor::default(),
        &seeds,
        &mut out,
    )
    .unwrap();
    assert_eq!(n, 6);
    // 与原版二进制产物逐字节一致（含 60 列折行与 header 全字段）
    assert_eq!(
        out,
        std::fs::read(fixture("smoke.orig.fa")).unwrap(),
        "同种子序重放必须逐字节复刻原版——分歧即核心移植缺陷"
    );
}

// ---------------------------------------------------------------- CLI 对拍

fn bin_path() -> &'static str {
    env!("CARGO_BIN_EXE_inchworm")
}

fn run_cli_in(cwd: &Path, args: &[&str]) -> (i32, String, String) {
    let out = Command::new(bin_path())
        .args(args)
        .current_dir(cwd)
        .output()
        .expect("spawn inchworm");
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

/// CLI 冒烟对拍（默认种子序）: contig 数 / header 多重集（去 aN）/ **rc 不变序列
/// 多重集**（DS 语义）与原版黄金一致;残余 3 条仅链方向（原版哈希序任意选链，
/// 本移植确定性取 canonical 大值链）——由上方重放测试归因。
#[test]
fn cli_smoke_matches_original_rc_invariant_multiset() {
    let cwd = tmpdir("cli");
    let kmers = p1_fixture("smoke.kmers.fa").to_str().unwrap().to_string();
    let (code, stdout, stderr) = run_cli_in(
        &cwd,
        &[
            "--kmers",
            &kmers,
            "--run_inchworm",
            "-K",
            "25",
            "--monitor",
            "1",
            "--DS",
            "--num_threads",
            "1",
        ],
    );
    assert_eq!(code, 0, "stderr:\n{stderr}");

    let golden = read_fasta(&fixture("smoke.orig.fa"));
    let ours = read_fasta_txt(&stdout, "cli_single");

    // contig 数
    assert_eq!(ours.len(), golden.len(), "contig 数与原版一致");
    // header 多重集（去 aN 序号——aN 依产出序，平局序不同即不同，不断言）
    let strip = |h: &String| h.split_once(';').map(|(_, r)| r.to_string()).unwrap();
    assert_eq!(
        multiset(ours.iter().map(|(h, _)| strip(h))),
        multiset(golden.iter().map(|(h, _)| strip(h))),
        "header（avg_cov/total_counts/Seed/K/length）多重集须一致"
    );
    // rc 不变序列多重集（DS: contig 与 revcomp 同义分子）
    assert_eq!(
        multiset(ours.iter().map(|(_, s)| rc_key(s))),
        multiset(golden.iter().map(|(_, s)| rc_key(s))),
        "rc 不变 contig 多重集须一致（残余差异应仅为链方向，由重放测试归因于种子序）"
    );

    // stderr 关键 monitor 行（供 xcheck 抓取的文案在此锁定）
    for line in [
        "Kmer length set to: 25",
        "Monitor turned on, set to: 1",
        "double stranded mode set",
        "setting number of threads to: 1",
        "-reading Kmer occurrences...",
        " done parsing 1887 Kmers, 1887 added, taking 0 seconds.",
        "TIMING KMER_DB_BUILDING 0 s.",
        "Pruning kmers (min_kmer_count=1 min_any_entropy=0 min_ratio_non_error=0.005)",
        "Pruned 0 kmers from catalog.",
        "TIMING PRUNING 0 s.",
        "-populating the kmer seed candidate list.",
        "Kcounter hash size: 1887",
        "Processed 1887 non-zero abundance kmers in kcounter.",
        "Total kcounter hash size: 1887 vs. sorted list size: 1887",
        "num threads set to: 1",
        "TIMING CONTIG_BUILDING 0 s.",
        "TIMING PROG_RUNTIME ",
    ] {
        assert!(
            stderr.contains(line),
            "stderr 缺关键 monitor 行: {line:?}\n--- stderr ---\n{stderr}"
        );
    }

    // inchworm.kmer_count 写在 CWD（原版行为，IRKE.cpp:144-147）
    assert_eq!(
        std::fs::read_to_string(cwd.join("inchworm.kmer_count")).unwrap(),
        "1887\n"
    );
}

fn read_fasta_txt(text: &str, tag: &str) -> Vec<(String, String)> {
    // tag 须调用方唯一——本测试二进制内各 CLI 测试并发运行，共享目录会互删
    let p = tmpdir(&format!("txt_{tag}"));
    let f = p.join("out.fa");
    std::fs::write(&f, text).unwrap();
    read_fasta(&f)
}

/// 参数校验: 未知参数 / 缺 --run_inchworm / 缺输入 / -K 越界一律 exit 2;
/// --PARALLEL_IWORM 自 T6 起启用（exit 0）。
#[test]
fn cli_validation_exit_2() {
    let cwd = tmpdir("validation");
    let kmers = p1_fixture("smoke.kmers.fa").to_str().unwrap().to_string();
    let base = |extra: &[&str]| -> Vec<String> {
        ["--kmers", &kmers, "--run_inchworm"]
            .iter()
            .map(|s| s.to_string())
            .chain(extra.iter().map(|s| s.to_string()))
            .collect()
    };
    let run = |args: &[String]| {
        let refs: Vec<&str> = args.iter().map(|s| s.as_str()).collect();
        run_cli_in(&cwd, &refs)
    };

    // 快乐路径 exit 0
    assert_eq!(run(&base(&[])).0, 0);

    let (code, _, err) = run(&base(&["--bogus"]));
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("bogus"), "{err}");

    let (code, _, err) = run(&["--kmers".to_string(), kmers.clone()]);
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("--run_inchworm"), "{err}");

    let (code, _, _) = run(&["--run_inchworm".to_string()]);
    assert_eq!(code, 2, "缺 --reads/--kmers");

    for bad in ["0", "33", "abc"] {
        let (code, _, err) = run(&base(&["-K", bad]));
        assert_eq!(code, 2, "-K {bad}: {err}");
    }

    let (code, _, err) = run(&base(&["--PARALLEL_IWORM"]));
    assert_eq!(code, 0, "--PARALLEL_IWORM 自 T6 起启用: {err}");

    // --SINGLE_PHASE 单独给出可收下（原版同样仅在 PARALLEL 下读取）
    assert_eq!(run(&base(&["--SINGLE_PHASE"])).0, 0);
    // --keep_tmp_files / --no_prune_error_kmers 收下
    assert_eq!(
        run(&base(&["--keep_tmp_files", "--no_prune_error_kmers"])).0,
        0
    );
}

// ---------------------------------------------------------------- PARALLEL 对拍

/// PARALLEL_IWORM CLI 冒烟: rayon chunk 并行 + dashmap 弱一致清零 + TWO_PHASE。
/// 断言（PARALLEL 语义本就 nondeterministic——与原版同为多线程竞态，故锁弱不变量）:
/// - exit 0、产出非空 FASTA
/// - 两次运行（6 线程 vs 2 线程）的 **rc 不变序列多重集**一致——chunk 划分由
///   全局种子序确定 + 每 chunk 独立 srand(1)，竞态窗口外的行为完全确定;
///   窗口内差异属原版同款 nondeterminism（本 fixture 实测稳定，见 T6 对拍报告）
/// - stderr 含 PARALLEL 特征行（不排序 + 线程数）
/// - --SINGLE_PHASE（关闭 TWO_PHASE）正常运行
#[test]
fn cli_parallel_iworm_smoke() {
    let cwd = tmpdir("parallel");
    let kmers = p1_fixture("smoke.kmers.fa").to_str().unwrap().to_string();

    let run_parallel = |nthreads: &str| {
        run_cli_in(
            &cwd,
            &[
                "--kmers",
                &kmers,
                "--run_inchworm",
                "-K",
                "25",
                "--monitor",
                "1",
                "--DS",
                "--num_threads",
                nthreads,
                "--PARALLEL_IWORM",
                "-L",
                "25",
            ],
        )
    };

    let (code, stdout, stderr) = run_parallel("6");
    assert_eq!(code, 0, "stderr:\n{stderr}");

    let ours = read_fasta_txt(&stdout, "par6");
    assert!(!ours.is_empty(), "PARALLEL 须产出 contig");

    // PARALLEL 特征 stderr 行（供 xcheck 抓取）
    for line in [
        "-setting parallel iworm mode.",
        "-Not sorting list of kmers, given parallel mode in effect.",
        "num threads set to: 6",
    ] {
        assert!(
            stderr.contains(line),
            "stderr 缺 PARALLEL 特征行: {line:?}\n--- stderr ---\n{stderr}"
        );
    }

    // 二次运行（不同线程数）: rc 不变多重集**高度重合**——锁"竞态窗口外确定"的
    // 弱不变量。PARALLEL 语义本就 nondeterministic（原版同为多线程竞态）:
    // chunk 划分由全局种子序确定 + 每 chunk 独立 srand(1) → 窗口外完全确定;
    // 窗口内（两 chunk 并发走/清同一图区域）偶发增补短 contig——本 fixture
    // 实测 100 次运行 98 次逐字节相同、2 次各多 1 条竞态 contig（6 条稳定
    // contig 始终齐备）。断言交集 ≥ 6-1（允许一条被竞态分裂）。
    let (code2, stdout2, _) = run_parallel("2");
    assert_eq!(code2, 0);
    let ours2 = read_fasta_txt(&stdout2, "par2");
    let m1 = multiset(ours.iter().map(|(_, s)| rc_key(s)));
    let m2 = multiset(ours2.iter().map(|(_, s)| rc_key(s)));
    let intersection: usize = m1
        .keys()
        .map(|k| (*m1.get(k).unwrap_or(&0)).min(*m2.get(k).unwrap_or(&0)))
        .sum();
    assert!(
        intersection + 1 >= m1.len().max(m2.len()),
        "PARALLEL 两次运行的 rc 不变 contig 多重集须高度重合（竞态窗口外确定）:\
         \n  6 线程: {m1:?}\n  2 线程: {m2:?}"
    );

    // --SINGLE_PHASE（关闭 TWO_PHASE）: 正常运行，contig 仍产出
    let (code_sp, stdout_sp, _) = run_cli_in(
        &cwd,
        &[
            "--kmers",
            &kmers,
            "--run_inchworm",
            "-K",
            "25",
            "--DS",
            "--num_threads",
            "4",
            "--PARALLEL_IWORM",
            "--SINGLE_PHASE",
        ],
    );
    assert_eq!(code_sp, 0);
    assert!(!read_fasta_txt(&stdout_sp, "par_sp").is_empty());
}
