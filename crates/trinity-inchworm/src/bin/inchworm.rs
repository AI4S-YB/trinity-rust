//! inchworm CLI — IRKE_run.cpp run_IRKE 的单线程移植（P2-T5）。
//!
//! ```text
//! inchworm --kmers <kmer_dump.fa> --run_inchworm [-K 25] [--DS] [--monitor 1] ...
//! inchworm --reads <reads.fa>    --run_inchworm [...]
//! ```
//!
//! 与原版的行为对应（IRKE_run.cpp 行号）:
//! - 必填: (--reads | --kmers) 之一 + --run_inchworm，否则打印 usage 退出
//!   （原版 return 1，本 CLI 按本仓库约定 exit 2——已知偏差，见下方 usage 注）
//! - 未知参数报错（原版 ArgProcessor 也报 "do not understand option"）
//! - --PARALLEL_IWORM: P2-T6 起启用 rayon chunk 并行组装（T5 之前 exit 2 拒绝）
//! - --SINGLE_PHASE: 仅在 --PARALLEL_IWORM 下生效（关闭 TWO_PHASE 二段重建，
//!   IRKE_run.cpp:436-438;单独给出本就是 no-op——原版同样）
//! - --keep_tmp_files: 收下但 no-op——本实现不落 tmp.iworm.fa（原版写后即删）
//! - --num_threads: PARALLEL 下建对应大小的 rayon 池;未给出时用 RAYON_NUM_THREADS
//!   或核数（对齐原版 omp_get_max_threads() 默认）。非 PARALLEL 下仅读解析级
//!   参数，组装恒单线程（原版 omp_set_num_threads(1)，IRKE.cpp:467）
//! - gzip 魔数嗅探读入（io_util）——原版 Fasta_reader 只读裸文件，此处是超集
//!
//! stderr 文案镜像原版（关键行供 xcheck 抓取）: -reading Kmer occurrences... /
//! done parsing N Kmers, M added... / TIMING KMER_DB_BUILDING / Pruning kmers... /
//! Pruned N kmers from catalog. / TIMING PRUNING / -populating the kmer seed
//! candidate list. / Kcounter hash size / TIMING CONTIG_BUILDING / TIMING
//! PROG_RUNTIME。

use std::io::{BufWriter, Read, Write};
use std::path::Path;
use std::time::Instant;

use trinity_common::cli::{Cli, CliError};
use trinity_common::io_util::open_maybe_gz;
use trinity_inchworm::irke::{
    compute_sequence_assemblies, compute_sequence_assemblies_parallel, populate_from_kmers,
    populate_from_reads, prune_some_kmers, write_kmer_count_report, AssemblyParams, IrkeParams,
    Monitor,
};
use trinity_inchworm::kmer_counter::KmerCounter;

const USAGE: &str = "\
Usage:
  inchworm --reads <filename.fasta> --run_inchworm [opts]
  inchworm --kmers <filename.fasta> --run_inchworm [opts]

Required
  --reads  <str>   :fasta file containing reads
  --kmers  <str>   :fasta file containing kmers (jellyfish dump)
  --run_inchworm   :run inchworm assembly

Common options
  -K <int>          :kmer length (default 25, max 32)
  --DS              :double stranded mode (default: strand-specific)
  --monitor <int>   :monitoring level (default 0)
  --num_threads <int>
  --minKmerCount <int>
  -L <int>          :minimum length of an inchworm assembly (default 25)
  --min_assembly_coverage <int>  (default 2)
  --min_seed_entropy <float>     (default 1.5)
  --min_seed_coverage <int>      (default 2)
  --min_any_entropy <float>      (default 0)
  --min_con_ratio <float>        (default 0)
  --min_ratio_non_error <float>  (default 0.005)
  --no_prune_error_kmers
  --keep_tmp_files  (accepted; implementation writes no tmp files)
  --PARALLEL_IWORM  (parallel contig building via rayon chunked seeds)
  --SINGLE_PHASE    (with --PARALLEL_IWORM: skip the two-phase seed re-selection)";

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    std::process::exit(run(&args));
}

// 解析器: trinity-common::cli（P3-T0 起两 bin 共用;`CliError` 亦来自该模块——
// usage=true → exit 2（参数问题，本仓库 CLI 约定;原版 usage 路径 return 1）;
// false → exit 1（运行失败））。
fn read_file(path: &str) -> Result<Vec<u8>, CliError> {
    let mut data = Vec::new();
    match open_maybe_gz(Path::new(path)) {
        Ok(mut r) => {
            r.read_to_end(&mut data)
                .map_err(|e| CliError::run(format!("cannot read {path}: {e}")))?;
            Ok(data)
        }
        Err(e) => Err(CliError::run(format!("cannot read {path}: {e}"))),
    }
}

/// time(NULL) 级（秒）计时的移植载体——原版全部 TIMING 行都是整秒。
fn secs(elapsed: std::time::Duration) -> u64 {
    elapsed.as_secs()
}

fn run(args: &[String]) -> i32 {
    match run_inner(args) {
        Ok(code) => code,
        Err(e) => {
            if e.usage {
                eprintln!("{USAGE}");
                eprintln!("\nError, {}", e.msg);
                2
            } else {
                eprintln!("Error, {}", e.msg);
                1
            }
        }
    }
}

fn run_inner(args: &[String]) -> Result<i32, CliError> {
    let prog_start = Instant::now();

    let mut cli = Cli::new(args);

    // ---- 必填校验（IRKE_run.cpp:172-189） ----
    if cli.bool_flag("--help") || cli.bool_flag("-h") {
        eprintln!("{USAGE}");
        return Ok(2);
    }
    let reads_fasta = cli.value_flag("--reads")?;
    let kmers_fasta = cli.value_flag("--kmers")?;
    if reads_fasta.is_none() && kmers_fasta.is_none() {
        return Err(CliError::usage(
            "must set --reads or --kmers (and --run_inchworm)",
        ));
    }
    if !cli.bool_flag("--run_inchworm") {
        return Err(CliError::usage("--run_inchworm required"));
    }
    // IRKE_run.cpp:371-378: PARALLEL_IWORM / SINGLE_PHASE（后者仅在前者下读取）
    let parallel_iworm = cli.bool_flag("--PARALLEL_IWORM");
    let single_phase = cli.bool_flag("--SINGLE_PHASE");
    // --keep_tmp_files: 收下但 no-op——本实现不落 tmp.iworm.fa（下方回显处消费）

    // ---- 参数提取（输出顺序镜像 IRKE_run.cpp:191-406 的处理序） ----
    let kmer_length: usize = match cli.value_flag("-K")? {
        None => 25, // IRKE_run.cpp:89
        Some(v) => {
            let k = v
                .parse::<usize>()
                .map_err(|_| CliError::usage(format!("invalid -K value: {v}")))?;
            if !(1..=32).contains(&k) {
                return Err(CliError::usage(format!(
                    "kmer length exceeds 32 or is zero: {k}"
                )));
            }
            eprintln!("Kmer length set to: {k}"); // IRKE_run.cpp:193
            k
        }
    };

    let min_kmer_count = cli.uint_flag("--minKmerCount", 1u32)?; // IRKE_run.cpp:90
    if cli.was_given("--minKmerCount") {
        eprintln!("Min coverage set to: {min_kmer_count}"); // IRKE_run.cpp:199
    }

    // IRKE_run.cpp:90: `MIN_ASSEMBLY_LENGTH = kmer_length`（=25）在 -K 解析前
    // 定值——默认恒 25，不随 -K 变
    let min_assembly_length = cli.uint_flag::<u32>("-L", 25)? as usize;
    if cli.was_given("-L") {
        eprintln!("Min assembly length set to: {min_assembly_length}"); // IRKE_run.cpp:204
    }

    let min_assembly_coverage = cli.uint_flag("--min_assembly_coverage", 2u32)?; // IRKE_run.cpp:91
    if cli.was_given("--min_assembly_coverage") {
        eprintln!("Min assembly coverage set to: {min_assembly_coverage}"); // IRKE_run.cpp:208
    }

    let monitor = Monitor::new(cli.uint_flag("--monitor", 0u32)?);
    if cli.was_given("--monitor") {
        eprintln!("Monitor turned on, set to: {}", monitor.level); // IRKE_run.cpp:231
    }

    if cli.bool_flag("--keep_tmp_files") {
        eprintln!("-retaining tmp files"); // IRKE_run.cpp:237（本实现无 tmp 文件，仅文案）
    }

    let min_con_ratio = cli.float_flag("--min_con_ratio", 0.0f32)?; // IRKE_run.cpp:92

    let double_stranded = cli.bool_flag("--DS");
    if double_stranded {
        eprintln!("double stranded mode set"); // IRKE_run.cpp:244
    }

    let min_seed_entropy = cli.float_flag("--min_seed_entropy", 1.5f32)?; // IRKE_run.cpp:101
    if cli.was_given("--min_seed_entropy") {
        eprintln!("Min seed entropy set to: {min_seed_entropy}"); // IRKE_run.cpp:248
    }

    let min_seed_coverage = cli.uint_flag("--min_seed_coverage", 2u32)?; // IRKE_run.cpp:102
    if cli.was_given("--min_seed_coverage") {
        eprintln!("min seed coverage set to: {min_seed_coverage}"); // IRKE_run.cpp:253
    }

    let min_any_entropy = cli.float_flag("--min_any_entropy", 0.0f32)?; // IRKE_run.cpp:93
    if cli.was_given("--min_any_entropy") {
        eprintln!("min entropy set to: {min_any_entropy}"); // IRKE_run.cpp:258
    }

    let num_threads = cli.uint_flag("--num_threads", 1u32)?;
    let num_threads_given = cli.was_given("--num_threads");
    if num_threads_given {
        eprintln!("setting number of threads to: {num_threads}"); // IRKE_run.cpp:283
        if num_threads > 1 && !parallel_iworm {
            // 原版此参数只影响读解析与 PARALLEL 组装;非 PARALLEL 的组装恒单线程
            // （IRKE.cpp:467 omp_set_num_threads(1)），本实现读解析也是单线程
            eprintln!(
                "Error, --num_threads > 1 requested but contig building remains single-threaded without --PARALLEL_IWORM; running with 1 thread"
            );
        }
    }
    // PARALLEL 组装的线程数: 显式给出用之;否则 None → rayon 默认
    // （RAYON_NUM_THREADS 或核数，对齐原版未给 --num_threads 时取
    //  omp_get_max_threads() 的行为）
    let parallel_num_threads = if num_threads_given {
        Some(num_threads as usize)
    } else {
        None
    };

    // IRKE_run.cpp:354-359: 默认 prune_error_kmers=true
    let prune_error_kmers = !cli.bool_flag("--no_prune_error_kmers");
    let min_ratio_non_error = cli.float_flag("--min_ratio_non_error", 0.005f32)?; // IRKE_run.cpp:107
    if prune_error_kmers && cli.was_given("--min_ratio_non_error") {
        // IRKE_run.cpp:362
        eprintln!("Set to prune kmers below min ratio non-erro: {min_ratio_non_error}");
    }

    // IRKE_run.cpp:371-380: PARALLEL_IWORM/SINGLE_PHASE 的模式行（种子不排序
    // 标志 __DEVEL_no_kmer_sort 的生效点在 compute_sequence_assemblies_parallel）
    if parallel_iworm {
        eprintln!("-setting parallel iworm mode."); // IRKE_run.cpp:374
        if single_phase {
            eprintln!("-setting single phase parallel iworm build."); // IRKE_run.cpp:379
        }
    }

    // 未知参数校验（原版 ArgProcessor "do not understand option"）——
    // 须在任何副作用（读文件/写 inchworm.kmer_count）之前
    cli.finish()?;

    // ---- IRKE 对象等价参数（IRKE_run.cpp:448） ----
    let params = IrkeParams {
        min_connectivity: min_con_ratio,
        min_seed_entropy,
        min_seed_coverage,
    };
    let aparams = AssemblyParams {
        min_connectivity: min_con_ratio,
        min_assembly_length,
        min_assembly_coverage,
    };

    // ---- Kmer catalog construction（IRKE_run.cpp:454-472） ----
    let (mut counter, _parsed): (KmerCounter, usize) = if let Some(reads) = &reads_fasta {
        // reads 优先（原版 if/else if 序）——原版 IRKE.cpp:157-280 "-storing Kmers..."
        eprintln!("-storing Kmers...");
        let t0 = Instant::now();
        let data = read_file(reads)?;
        let r = populate_from_reads(&data, kmer_length, double_stranded)?;
        let t = secs(t0.elapsed());
        // IRKE.cpp:276: done parsing N **sequences**, extracted M kmers（reads 模式
        // 无 TIMING 行、不写 inchworm.kmer_count——均为原版行为）
        eprintln!(
            "\n done parsing {} sequences, extracted {} kmers, taking {t} seconds.",
            r.1,
            r.0.size()
        );
        r
    } else {
        let kmers = kmers_fasta.as_deref().unwrap();
        eprintln!("-reading Kmer occurrences..."); // IRKE.cpp:96
        let t0 = Instant::now();
        let data = read_file(kmers)?;
        let r = populate_from_kmers(&data, kmer_length, double_stranded)?;
        let t = secs(t0.elapsed());
        eprintln!(
            "\n done parsing {} Kmers, {} added, taking {t} seconds.",
            r.1,
            r.0.size()
        ); // IRKE.cpp:140
        eprintln!("\nTIMING KMER_DB_BUILDING {t} s."); // IRKE.cpp:142
                                                       // IRKE.cpp:144-147: 写 CWD 的 inchworm.kmer_count（原版固定相对路径）
        write_kmer_count_report(Path::new("inchworm.kmer_count"), r.0.size())?;
        r
    };

    // ---- 剪枝（IRKE_run.cpp:493-506） ----
    if min_kmer_count > 1 || min_any_entropy > 0.0 || prune_error_kmers {
        eprintln!(
            "Pruning kmers (min_kmer_count={min_kmer_count} min_any_entropy={min_any_entropy} min_ratio_non_error={min_ratio_non_error})"
        );
        let t0 = Instant::now();
        let pruned = prune_some_kmers(
            &mut counter,
            min_kmer_count,
            min_any_entropy,
            prune_error_kmers,
            min_ratio_non_error,
        );
        let t = secs(t0.elapsed());
        eprintln!("Pruned {pruned} kmers from catalog."); // KC:266（原版无条件打印）
        eprintln!("\tPruning time: {t} seconds = {} minutes.", t as f32 / 60.0);
        eprintln!("\nTIMING PRUNING {t} s.");
    }

    // ---- 组装（IRKE_run.cpp:520-536） ----
    // 注: 原版 "-populating the kmer seed candidate list."（IRKE_run.cpp:521）先于
    // 本行，本移植中种子列表由 compute_sequence_assemblies* 内建、其文案随函数
    // 输出——与 "-beginning" 行交换次序（已证文案逐条一致）
    eprintln!("-beginning inchworm contig assembly.");
    let t0 = Instant::now();
    let stdout = std::io::stdout();
    let mut writer = BufWriter::new(stdout.lock());
    let n_contigs = if parallel_iworm {
        // IRKE.cpp:504-620 并行分支: rayon chunk（dynamic,1000 镜像）+ dashmap
        // 弱一致目录 + TWO_PHASE（--SINGLE_PHASE 关闭，IRKE_run.cpp:436-438）
        compute_sequence_assemblies_parallel(
            counter,
            &params,
            &aparams,
            &monitor,
            !single_phase, // TWO_PHASE 默认 true（IRKE.cpp:59）
            parallel_num_threads,
            &mut writer,
        )?
    } else {
        compute_sequence_assemblies(
            &mut counter,
            &params,
            &aparams,
            &monitor,
            true, // 单线程: 种子按 count 降序（PARALLEL 走 --PARALLEL_IWORM）
            &mut writer,
        )?
    };
    writer
        .flush()
        .map_err(|e| CliError::run(format!("writing stdout: {e}")))?;
    let t = secs(t0.elapsed());
    eprintln!(
        "\tIworm contig assembly time: {t} seconds = {} minutes.",
        t as f32 / 60.0
    );
    eprintln!("\nTIMING CONTIG_BUILDING {t} s.");

    eprintln!("\nTIMING PROG_RUNTIME {} s.", secs(prog_start.elapsed())); // IRKE_run.cpp:541-545

    let _ = n_contigs; // （产物数由调用方从 stdout 统计;保留返回值便于 T6 记日志）
    Ok(0)
}
