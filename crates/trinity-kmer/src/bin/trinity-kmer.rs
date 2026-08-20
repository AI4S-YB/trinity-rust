//! trinity-kmer CLI — jellyfish count+dump 的一体化替代 + DigiNorm 编排。
//!
//! count:          trinity-kmer count --reads in.fa [-K 25] [--canonical] [--min-count 1] [-o out.fa]
//! coverage-stats: trinity-kmer coverage-stats --reads reads.fa --kmers dump.fa [-K 25] [--SS] [-o out.tsv]
//! diginorm:       trinity-kmer diginorm (--left l.fq --right r.fq | --single s.fq) -o outdir
//!                 [--SS_lib_type F|R|RF|FR] [--max_cov 200] [--min_cov 1] [--max_CV 10000] [-K 25]
//!
//! 参数校验（T2 审查 Important 修复）: 未知 flag 报错、值 flag 缺值报错、
//! -K 非 1..=32 报错——一律 exit 2（不再静默忽略/采用默认值）。
//! 解析器: trinity-common::cli（P3-T0 起两 bin 共用）。

use std::io::{BufWriter, Read};
use std::path::PathBuf;

use trinity_common::cli::{Cli, CliError};
use trinity_common::io_util::open_maybe_gz;
use trinity_kmer::diginorm::{run, DigiNormParams, ReadsInput};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let cmd = args.first().cloned().unwrap_or_default();
    // 空参数时 get(1..) 返回 None，不能直接切片 [1..]（会 panic 而非打印用法）
    let rest = args.get(1..).unwrap_or(&[]);
    let code = match cmd.as_str() {
        "count" => run_cli(cmd_count(rest)),
        "coverage-stats" => run_cli(cmd_coverage_stats(rest)),
        "diginorm" => run_cli(cmd_diginorm(rest)),
        _ => {
            eprintln!("用法: trinity-kmer <count|coverage-stats|diginorm> [参数]");
            eprintln!("  count --reads in.fa [-K 25] [--canonical] [--min-count 1] [-o out.fa]");
            eprintln!(
                "  coverage-stats --reads reads.fa --kmers dump.fa [-K 25] [--SS] [-o out.tsv]"
            );
            eprintln!("  diginorm (--left l.fq --right r.fq | --single s.fq) -o outdir \\");
            eprintln!("          [--SS_lib_type F|R|RF|FR] [--max_cov 200] [--min_cov 1] [--max_CV 10000] [-K 25]");
            2
        }
    };
    std::process::exit(code);
}

/// 统一出口: 参数类错误（unknown/缺值/-K 越界/缺必填）→ exit 2; 运行类错误 → 1。
fn run_cli(r: Result<(), CliError>) -> i32 {
    match r {
        Ok(()) => 0,
        Err(e) => {
            eprintln!("Error, {}", e.msg);
            if e.usage {
                2
            } else {
                1
            }
        }
    }
}

/// -K 解析与 1..=32 校验（sequenceUtil "kmer length exceeds 32"）。
fn parse_kmer(cli: &mut Cli) -> Result<usize, CliError> {
    match cli.value_flag("-K")? {
        None => Ok(25),
        Some(v) => {
            let k: usize = v.parse().map_err(|_| {
                CliError::usage(format!("invalid -K value: {v} (expected integer 1..=32)"))
            })?;
            if !(1..=32).contains(&k) {
                return Err(CliError::usage(format!(
                    "kmer length exceeds 32 or is zero: {k}"
                )));
            }
            Ok(k)
        }
    }
}

/// f64 flag + 有限性检查（"inf"/"nan" 拒绝——原 parse_f64 语义）。
fn parse_finite(cli: &mut Cli, name: &str, default: f64) -> Result<f64, CliError> {
    match cli.value_flag(name)? {
        None => Ok(default),
        Some(v) => {
            let x = v
                .parse::<f64>()
                .map_err(|_| CliError::usage(format!("invalid {name} value: {v}")))?;
            if x.is_finite() {
                Ok(x)
            } else {
                Err(CliError::usage(format!("invalid {name} value: {v}")))
            }
        }
    }
}

fn read_file(path: &str) -> Result<Vec<u8>, CliError> {
    let mut data = Vec::new();
    match open_maybe_gz(std::path::Path::new(path)) {
        Ok(mut r) => {
            r.read_to_end(&mut data)
                .map_err(|e| CliError::run(format!("cannot read {path}: {e}")))?;
            Ok(data)
        }
        Err(e) => Err(CliError::run(format!("cannot read {path}: {e}"))),
    }
}

fn writer_for(out: Option<String>) -> Result<Box<dyn std::io::Write>, CliError> {
    match out {
        Some(path) => match std::fs::File::create(&path) {
            Ok(f) => Ok(Box::new(BufWriter::new(f))),
            Err(e) => Err(CliError::run(format!("cannot write {path}: {e}"))),
        },
        None => Ok(Box::new(BufWriter::new(std::io::stdout().lock()))),
    }
}

// ---------------------------------------------------------------- count

fn cmd_count(rest: &[String]) -> Result<(), CliError> {
    let mut cli = Cli::new(rest);
    let reads = cli.req_flag("--reads")?;
    let k = parse_kmer(&mut cli)?;
    let canonical = cli.bool_flag("--canonical");
    let min_count = cli.uint_flag("--min-count", 1u32)?;
    let out = cli.value_flag("-o")?;
    cli.finish()?;
    let data = read_file(&reads)?;
    let counts = trinity_kmer::counter::KmerCountTable::count_fasta_data(&data, k, canonical);
    // -o 指定输出文件，缺省 stdout（jellyfish count -o 写 .jf，我们直接写 dump FASTA）
    let mut w = writer_for(out)?;
    trinity_kmer::dump::write_dump(&mut w, &counts, k, min_count, canonical)
        .map_err(|e| CliError::run(format!("writing dump: {e}")))?;
    Ok(())
}

// ---------------------------------------------------------------- coverage-stats

/// fastaToKmerCoverageStats 镜像: 默认 DS，--SS 关闭 canonical。
fn cmd_coverage_stats(rest: &[String]) -> Result<(), CliError> {
    let mut cli = Cli::new(rest);
    let reads = cli.req_flag("--reads")?;
    let kmers = cli.req_flag("--kmers")?;
    let k = parse_kmer(&mut cli)?;
    let ds = !cli.bool_flag("--SS");
    let out = cli.value_flag("-o")?;
    cli.finish()?;
    let reads_data = read_file(&reads)?;
    let kmer_data = read_file(&kmers)?;
    let counts = trinity_kmer::coverage_stats::load_kmer_dump(&kmer_data, k, ds)?;
    let rows = trinity_kmer::coverage_stats::coverage_stats_rows(&reads_data, &counts, k, ds)?;
    let mut w = writer_for(out)?;
    trinity_kmer::coverage_stats::write_stats_tsv(&mut w, &rows)
        .map_err(|e| CliError::run(format!("writing stats: {e}")))?;
    Ok(())
}

// ---------------------------------------------------------------- diginorm

/// SS_lib_type → (revcomp_left, revcomp_right, ds)，镜像原版:
/// - 空 = DS（jellyfish --canonical + stats --DS）;
/// - 单端 'R' 先按 L200-202 改写为 'F' → **单端恒不 revcomp**;
/// - PE `split(//, $ss)` 逐字符: 'RF' → left revcomp; 'FR' → **right** revcomp。
fn ss_flags(ss: &str, single: bool) -> Result<(bool, bool, bool), CliError> {
    if !["", "F", "R", "RF", "FR"].contains(&ss) {
        return Err(CliError::usage(format!(
            "unrecognized SS_lib_type value of {ss}. Should be: F, R, RF, or FR"
        )));
    }
    let ss = if single && ss == "R" { "F" } else { ss };
    let bytes = ss.as_bytes();
    let left = bytes.first() == Some(&b'R');
    let right = bytes.get(1) == Some(&b'R');
    Ok((left, right, ss.is_empty()))
}

/// 逗号列表展开（原版 L219-226 `split(",", join(",", @files))`——`--left a.fq,b.fq`）。
fn split_list(v: &str) -> Vec<PathBuf> {
    v.split(',').map(PathBuf::from).collect()
}

fn cmd_diginorm(rest: &[String]) -> Result<(), CliError> {
    let mut cli = Cli::new(rest);
    let reads = match (
        cli.value_flag("--left")?,
        cli.value_flag("--right")?,
        cli.value_flag("--single")?,
    ) {
        (Some(l), Some(r), None) => ReadsInput::Paired(split_list(&l), split_list(&r)),
        (None, None, Some(s)) => ReadsInput::Single(split_list(&s)),
        _ => {
            return Err(CliError::usage(
                "need either options 'left' and 'right' or option 'single'",
            ))
        }
    };
    let out_dir = cli.req_flag("-o")?;
    let ss_lib = cli.str_flag("--SS_lib_type", "")?;
    let single = matches!(reads, ReadsInput::Single(_));
    let (rev_l, rev_r, ds) = ss_flags(&ss_lib, single)?;
    let max_cov = parse_finite(&mut cli, "--max_cov", 200.0)?;
    let min_cov = parse_finite(&mut cli, "--min_cov", 1.0)?;
    let max_cv = parse_finite(&mut cli, "--max_CV", 10000.0)?;
    // 原版 L231: max_cov >= 2
    if max_cov < 2.0 {
        return Err(CliError::usage("need to set --max_cov at least 2"));
    }
    let k = parse_kmer(&mut cli)?;
    cli.finish()?;
    let params = DigiNormParams {
        k,
        max_cov,
        min_cov,
        max_cv,
        ss_revcomp_left: rev_l,
        ss_revcomp_right: rev_r,
        ds,
        ..Default::default()
    };
    let out = run(&params, &reads, std::path::Path::new(&out_dir))?;
    let files: Vec<String> = [out.left, out.right, out.single]
        .iter()
        .flatten()
        .map(|p| p.display().to_string())
        .collect();
    eprintln!(
        "\nNormalization complete. See outputs: \n\t{}",
        files.join("\n\t")
    );
    Ok(())
}
