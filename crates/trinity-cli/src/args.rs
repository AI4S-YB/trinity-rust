//! 原版 `Trinity` 主程序参数面移植（GetOptions L660-760 主线子集 +
//! 校验 L1088-1130 / L969-985 镜像）。参数名与原版逐一同名;
//! 未知参数经 [`Cli::finish`] 报 "do not understand option"。

use std::path::PathBuf;

use trinity_common::cli::{Cli, CliError};

/// `--seqType fq|fa`（原版还有 cfa/cfq——主线不支持，prep_seqs 直接 confess）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqType {
    Fq,
    Fa,
}

/// 主线参数集（原版 17 个同名 flag + Rust 版扩展 `--bfly_stack_mb`）。
///
/// 继续性: workdir 既有 outdir 恢复（`.ok` 断点）天然支持，无 `--resume`
/// 专门旗标（原版即如此）。
#[derive(Debug, Clone)]
pub struct TrinityArgs {
    pub seq_type: SeqType,
    /// 逗号列表（原版 `left=s{,}` 逗号展开，L1127）。
    pub left: Vec<PathBuf>,
    pub right: Vec<PathBuf>,
    pub single: Vec<PathBuf>,
    /// F/R/FR/RF（PE 只 FR|RF; SE 只 F|R; 与 PE/SE 混用报错——原版 L1098-1117）。
    pub ss_lib_type: Option<String>,
    /// 原版默认 2。
    pub cpu: usize,
    /// `--max_memory xG` 字节数（必填，原版 L975-984 die）。
    pub max_memory: u64,
    pub output: PathBuf,
    /// 原版 `__KMER_SIZE` 默认 25（本实现暴露为 `--KMER_SIZE`）; 1..=32。
    pub kmer_size: usize,
    /// 默认 200; >= ABSOLUTE_MIN_CONTIG_LENGTH（=100，原版 L29/L969-970）。
    pub min_contig_length: usize,
    /// 原版 `min_kmer_cov` 默认 1。
    pub min_kmer_count: u32,
    /// 默认 500（原版 L96）。
    pub group_pairs_distance: usize,
    pub no_normalize_reads: bool,
    /// 默认 200（原版 L214）。
    pub normalize_max_read_cov: u32,
    /// 原版默认 6（L119; 后续 min(CPU,6) 收敛——原版 L1036/L1233 语义）。
    pub inchworm_cpu: usize,
    /// **Rust 版扩展**（原版 java -Xmx 侧; 这里是 butterfly 线程栈 MB）: 默认 256。
    pub bfly_stack_mb: usize,
    /// 默认 200000（原版 L144）。
    pub max_reads_per_graph: u64,
    /// 收尾不删中间文件。
    pub no_cleanup: bool,
}

/// `--max_memory 10G` → 字节。镜像原版 L1156 `^(\d+)G`（大小写不敏感）;
/// 其余格式报错（L1174 "must specify max memory ... eg. --max_memory 10G"）。
fn parse_max_memory(v: &str) -> Result<u64, CliError> {
    let body = v.strip_suffix(['G', 'g']).filter(|s| !s.is_empty());
    match body.and_then(|s| s.parse::<u64>().ok()) {
        Some(g) => Ok(g * 1024 * 1024 * 1024),
        _ => Err(CliError::usage(
            "Error, must specify max memory for jellyfish to use, eg.  --max_memory 10G",
        )),
    }
}

/// 逗号列表 → PathBuf（原版 split(",", ...)）。
fn path_list(v: &str) -> Vec<PathBuf> {
    v.split(',').map(PathBuf::from).collect()
}

/// 双名旗标（主名/别名）取整数值; 均未给出 → default。重复名后取先生效者。
fn aliased_uint<T: std::str::FromStr>(
    cli: &mut Cli,
    primary: &str,
    alias: &str,
    default: T,
) -> Result<T, CliError> {
    match cli.value_flag(primary)?.or(cli.value_flag(alias)?) {
        None => Ok(default),
        Some(v) => v
            .parse::<T>()
            .map_err(|_| CliError::usage(format!("invalid {primary} value: {v}"))),
    }
}

pub fn parse_args(argv: &[String]) -> Result<TrinityArgs, CliError> {
    let mut cli = Cli::new(argv);

    let seq_type = match cli.req_flag("--seqType")?.as_str() {
        "fq" => SeqType::Fq,
        "fa" => SeqType::Fa,
        other => {
            return Err(CliError::usage(format!(
                "Error, unrecognized seqType value of {other}. Should be: fq or fa"
            )))
        }
    };
    let left = cli
        .value_flag("--left")?
        .map(|v| path_list(&v))
        .unwrap_or_default();
    let right = cli
        .value_flag("--right")?
        .map(|v| path_list(&v))
        .unwrap_or_default();
    let single = cli
        .value_flag("--single")?
        .map(|v| path_list(&v))
        .unwrap_or_default();
    let ss_lib_type = cli.value_flag("--SS_lib_type")?;
    let cpu: usize = cli.uint_flag("--CPU", 2)?;
    let max_memory = parse_max_memory(&cli.req_flag("--max_memory")?)?;
    let output = PathBuf::from(cli.req_flag("--output")?);
    // T3 别名: `__KMER_SIZE`（原版主程序名）与 `--min_kmer_cov`（原版名）。
    let kmer_size: usize = aliased_uint(&mut cli, "--KMER_SIZE", "__KMER_SIZE", 25)?;
    let min_contig_length: usize = cli.uint_flag("--min_contig_length", 200)?;
    let min_kmer_count: u32 = aliased_uint(&mut cli, "--min_kmer_count", "--min_kmer_cov", 1)?;
    let group_pairs_distance: usize = cli.uint_flag("--group_pairs_distance", 500)?;
    let no_normalize_reads = cli.bool_flag("--no_normalize_reads");
    let normalize_max_read_cov: u32 = cli.uint_flag("--normalize_max_read_cov", 200)?;
    let inchworm_cpu_in: Option<usize> = match cli.value_flag("--inchworm_cpu")? {
        None => None,
        Some(v) => Some(
            v.parse::<usize>()
                .map_err(|_| CliError::usage(format!("invalid --inchworm_cpu value: {v}")))?,
        ),
    };
    let bfly_stack_mb: usize = cli.uint_flag("--bfly_stack_mb", 256)?;
    let max_reads_per_graph: u64 = cli.uint_flag("--max_reads_per_graph", 200000)?;
    let no_cleanup = cli.bool_flag("--no_cleanup");

    cli.finish()?;

    // ---- 校验（镜像原版 L1088-1130 / L969-985）----
    if !left.is_empty() && !single.is_empty() {
        return Err(CliError::usage(
            "Error, cannot mix PE and SE reads by using --left, --right, and --single.  See Trinity FAQ for how to combine SE and PE data",
        ));
    }
    if let Some(ss) = &ss_lib_type {
        if !matches!(ss.as_str(), "R" | "F" | "RF" | "FR") {
            return Err(CliError::usage(format!(
                "Error, unrecognized SS_lib_type value of {ss}. Should be: F, R, RF, or FR"
            )));
        }
        if !single.is_empty() && ss.len() != 1 {
            return Err(CliError::usage(
                "Error, with --single reads, the --SS_lib_type can be 'F' or 'R' only.",
            ));
        }
        if !left.is_empty() && ss.len() != 2 {
            return Err(CliError::usage(
                "Error, with paired end reads, the --SS_lib_type can be 'RF' or 'FR' only.",
            ));
        }
    }
    if left.is_empty() != right.is_empty() {
        return Err(CliError::usage(
            "Error, need either options 'left' and 'right' or option 'single'",
        ));
    }
    if left.is_empty() && single.is_empty() {
        return Err(CliError::usage(
            "Error, need either options 'left' and 'right' or option 'single'",
        ));
    }
    if !(1..=32).contains(&kmer_size) {
        return Err(CliError::usage(format!(
            "invalid --KMER_SIZE value: {kmer_size} (must be 1..=32)"
        )));
    }
    // ABSOLUTE_MIN_CONTIG_LENGTH = 100（原版 L29）。
    if min_contig_length < 100 {
        return Err(CliError::usage(format!(
            "sorry, min contig length set at {min_contig_length} is below our imposed threshold of 100 and might lead to undesirably long runtimes and numbers of transcript clusters to pursue (and number of intermediate files generated)."
        )));
    }
    // 原版: inchworm_cpu 默认 6，但 CPU 更小时收敛到 CPU（L1233 语义）。
    let inchworm_cpu = inchworm_cpu_in.unwrap_or(6).min(cpu.max(1));

    Ok(TrinityArgs {
        seq_type,
        left,
        right,
        single,
        ss_lib_type,
        cpu,
        max_memory,
        output,
        kmer_size,
        min_contig_length,
        min_kmer_count,
        group_pairs_distance,
        no_normalize_reads,
        normalize_max_read_cov,
        inchworm_cpu,
        bfly_stack_mb,
        max_reads_per_graph,
        no_cleanup,
    })
}
