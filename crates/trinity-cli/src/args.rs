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

/// `--path-engine classic|gtasm|rnaspades`（**Rust 版扩展**, 原版无）。
///
/// - `classic`（默认）: butterfly allProbPaths 全链, 行为与现状逐字节一致;
/// - `gtasm`: GT 边打分（GTW1 权重, `--gt-weights`）× GTasm 式解码 × 原始
///   reads 边支持（python 原型 gt_decode2.py C 变体的产品化, 详见
///   trinity-butterfly `gtmodel`/`gtdecode` 模块）;
/// - `rnaspades`: rnaSPAdes PathExtend MultiExtender 核心规则复刻（覆盖
///   一致性贪心延伸 + 简单 bulge 双 isoform 分叉 + 其余歧义停止 +
///   path_cleaning 110bp/cov≥2, 详见 trinity-butterfly `rnaspades` 模块）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PathEngine {
    Classic,
    Gtasm,
    RnaSpades,
}

/// `--multi-k`（**Rust 版扩展**, rnaSPAdes 式多 k 迭代 + trusted contig 注入）。
/// 缺省 = None（单 k classic, 行为逐字节不变）。
///
/// - `Auto`: 按 reads 读长（前 1000 条抽样）走 rnaSPAdes `rna_k_values` 公式
///   （multik::rna_k_values）, 再钳入支持域 [15, 32];
/// - `List`: 显式递升 k 列表（逗号分隔, 每个 k ∈ [15, 32]）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MultiK {
    Auto,
    List(Vec<usize>),
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
    /// **Rust 版扩展**: 中间产物（归一化目录、both.fa、mer_counts、inchworm、
    /// chrysalis 整个目录、butterfly 中间物）全部写到 `<scratch>/<output 目录名>/`
    /// 而非 --output 下; 最终 `<output>.Trinity.fasta` 与 gene_trans_map 仍落在
    /// --output 旁。checkpoint（.ok）随产物放在 scratch 内——scratch 是易失的,
    /// 换 scratch 目录 = 全部重算。默认 None = 原行为（全落 --output）。
    pub scratch: Option<PathBuf>,
    /// `--recursive`（M1）: 递归分区模式——phase-1 chrysalis 止于分区
    /// （skip_quantify）, 逐分区进程内迷你 Trinity, 产物落
    /// `read_partitions/c<id>.mini.fa`（存在即跳过 = 断点）。
    pub recursive: bool,
    /// **Rust 版扩展**: butterfly 组件的路径引擎（默认 classic）。
    pub path_engine: PathEngine,
    /// **Rust 版扩展**: GTasm 引擎的 GTW1 权重路径（gtasm 必填）。
    pub gt_weights: Option<PathBuf>,
    /// **Rust 版扩展**: rnaSPAdes 式多 k 迭代（`--multi-k auto|k1,k2,...`）。
    /// None = 关闭（单 k classic, 现有行为逐字节不变）。与 --KMER_SIZE /
    /// --recursive 互斥（防歧义, 沿 gtasm×recursive 先例）。
    pub multi_k: Option<MultiK>,
    /// **Rust 版扩展**: rnaSPAdes 式相对覆盖图清理（chrysalis QuantifyGraph
    /// 后对 c*.graph.out 做 tip/EC/rcc 删除, 见 trinity-chrysalis::graph_clean）。
    /// 默认关 = 现有行为逐字节不变。--recursive 下暂不生效（迷你管线未接线）。
    pub graph_clean_rnaspades: bool,
    /// **Rust 版扩展**: SPAdes 式 gap closer（chrysalis QuantifyGraph/清理后
    /// 对 c*.graph.out 做 tip 对 tip overlap 缝合, PE 证据排序候选, 见
    /// trinity-chrysalis::gap_close）。默认关 = 现有行为逐字节不变。
    /// 仅 PE 数据生效（单端 = SPAdes "no paired-end libraries" 跳过语义）;
    /// --recursive 下暂不生效（迷你管线未接线）。
    pub gap_close_spades: bool,
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

/// `--multi-k k1,k2,...` 显式列表解析: 逗号分隔整数, 每个 k ∈
/// [multik::MULTI_K_MIN, multik::MULTI_K_MAX]（15..=32, KmerId 位宽上限）,
/// 严格递升。空段/非数字/越界/乱序/重复一律 usage 报错。
fn parse_multi_k_list(v: &str) -> Result<Vec<usize>, CliError> {
    use crate::multik::{MULTI_K_MAX, MULTI_K_MIN};
    let mut ks = Vec::new();
    for tok in v.split(',') {
        let k: usize = tok
            .parse()
            .map_err(|_| CliError::usage(format!("invalid --multi-k value: {v}")))?;
        if !(MULTI_K_MIN..=MULTI_K_MAX).contains(&k) {
            return Err(CliError::usage(format!(
                "invalid --multi-k value: {k} (must be {MULTI_K_MIN}..={MULTI_K_MAX})"
            )));
        }
        ks.push(k);
    }
    if ks.is_empty() {
        return Err(CliError::usage("invalid --multi-k value: empty list"));
    }
    if !ks.windows(2).all(|w| w[0] < w[1]) {
        return Err(CliError::usage(format!(
            "invalid --multi-k value: {v} (k 列表须严格递升)"
        )));
    }
    Ok(ks)
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
    let recursive = cli.bool_flag("--recursive");
    let scratch = cli.value_flag("--scratch")?.map(PathBuf::from);
    let path_engine = match cli.value_flag("--path-engine")?.as_deref() {
        None | Some("classic") => PathEngine::Classic,
        Some("gtasm") => PathEngine::Gtasm,
        Some("rnaspades") => PathEngine::RnaSpades,
        Some(other) => {
            return Err(CliError::usage(format!(
                "invalid --path-engine value: {other} (must be classic|gtasm|rnaspades)"
            )))
        }
    };
    let gt_weights = cli.value_flag("--gt-weights")?.map(PathBuf::from);
    // --multi-k auto | k1,k2,...（Rust 版扩展, rnaSPAdes 式多 k 迭代）
    let multi_k = match cli.value_flag("--multi-k")?.as_deref() {
        None => None,
        Some("auto") => Some(MultiK::Auto),
        Some(v) => Some(MultiK::List(parse_multi_k_list(v)?)),
    };
    // --graph-clean-rnaspades（Rust 版扩展, 默认关）
    let graph_clean_rnaspades = cli.bool_flag("--graph-clean-rnaspades");
    // --gap-close-spades（Rust 版扩展, 默认关）
    let gap_close_spades = cli.bool_flag("--gap-close-spades");

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

    // gtasm 必须给权重; 非 gtasm 给了权重也只警告（不误判为 gtasm）。
    if path_engine == PathEngine::Gtasm && gt_weights.is_none() {
        return Err(CliError::usage(
            "Error, --path-engine gtasm requires --gt-weights <GTW1 权重文件路径>",
        ));
    }
    if path_engine != PathEngine::Gtasm && gt_weights.is_some() {
        eprintln!("WARNING: --gt-weights 仅在 --path-engine gtasm 下生效, 忽略");
    }
    // gtasm/rnaspades × recursive 未实现（递归模式的分区池不走 butterfly 组件池）
    if path_engine != PathEngine::Classic && recursive {
        return Err(CliError::usage(format!(
            "Error, --path-engine {} 与 --recursive 暂不支持同时使用",
            match path_engine {
                PathEngine::Gtasm => "gtasm",
                PathEngine::RnaSpades => "rnaspades",
                PathEngine::Classic => unreachable!(),
            }
        )));
    }
    // multi-k 互斥校验（防歧义）: 多 k 轮次自带 k 序列, 与 --KMER_SIZE 同给
    // 语义冲突（哪个生效?）; 与 --recursive 的组合未实现（沿 gtasm 先例）。
    if multi_k.is_some() {
        if cli.was_given("--KMER_SIZE") || cli.was_given("__KMER_SIZE") {
            return Err(CliError::usage(
                "Error, --multi-k 与 --KMER_SIZE 互斥: 多 k 迭代的 k 序列由 --multi-k 给定, 请勿再指定 --KMER_SIZE",
            ));
        }
        if recursive {
            return Err(CliError::usage(
                "Error, --multi-k 与 --recursive 暂不支持同时使用",
            ));
        }
    }
    // --graph-clean-rnaspades 在 --recursive 下未接线（迷你管线走默认参数）→ 明告
    if graph_clean_rnaspades && recursive {
        eprintln!("WARNING: --graph-clean-rnaspades 在 --recursive 下暂不生效, 忽略");
    }
    // --gap-close-spades: 同清理的 recursive 处置; 单端数据无 PE 证据,
    // 对应 SPAdes "No paired-end libraries exist, skipping gap closer"
    // （gap_closer.cpp execute 开头）→ 明告并置关。
    if gap_close_spades && recursive {
        eprintln!("WARNING: --gap-close-spades 在 --recursive 下暂不生效, 忽略");
    }
    let gap_close_spades = if gap_close_spades && left.is_empty() {
        eprintln!("WARNING: --gap-close-spades 需要 paired-end reads（--left/--right）, 单端数据下忽略");
        false
    } else {
        gap_close_spades
    };

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
        scratch,
        recursive,
        path_engine,
        gt_weights,
        multi_k,
        graph_clean_rnaspades,
        gap_close_spades,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(extra: &[&str]) -> Vec<String> {
        let mut v: Vec<String> = [
            "--seqType",
            "fq",
            "--left",
            "l.fq",
            "--right",
            "r.fq",
            "--max_memory",
            "2G",
            "--output",
            "out",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        v.extend(extra.iter().map(|s| s.to_string()));
        v
    }

    #[test]
    fn scratch_defaults_to_none() {
        let a = parse_args(&argv(&[])).unwrap();
        assert_eq!(a.scratch, None);
    }

    #[test]
    fn gap_close_spades_flag_and_se_downgrade() {
        // 默认关
        let a = parse_args(&argv(&[])).unwrap();
        assert!(!a.gap_close_spades);
        // PE + flag → 开
        let a = parse_args(&argv(&["--gap-close-spades"])).unwrap();
        assert!(a.gap_close_spades);
        // 单端 + flag → 警告并置关（SPAdes "no paired-end libraries" 语义）
        let se: Vec<String> = [
            "--seqType", "fq", "--single", "s.fq", "--max_memory", "2G", "--output", "out",
            "--gap-close-spades",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        let a = parse_args(&se).unwrap();
        assert!(!a.gap_close_spades);
    }

    #[test]
    fn scratch_parses_dir() {
        let a = parse_args(&argv(&["--scratch", "/tmp/scratch"])).unwrap();
        assert_eq!(a.scratch, Some(PathBuf::from("/tmp/scratch")));
    }

    #[test]
    fn path_engine_defaults_classic() {
        let a = parse_args(&argv(&[])).unwrap();
        assert_eq!(a.path_engine, PathEngine::Classic);
        assert_eq!(a.gt_weights, None);
    }

    #[test]
    fn path_engine_gtasm_requires_weights() {
        // gtasm 无 --gt-weights → 报错
        assert!(parse_args(&argv(&["--path-engine", "gtasm"])).is_err());
        // 带上权重 → 通过
        let a = parse_args(&argv(&[
            "--path-engine",
            "gtasm",
            "--gt-weights",
            "/tmp/m.gtw",
        ]))
        .unwrap();
        assert_eq!(a.path_engine, PathEngine::Gtasm);
        assert_eq!(a.gt_weights, Some(PathBuf::from("/tmp/m.gtw")));
        // 非法值 → 报错
        assert!(parse_args(&argv(&["--path-engine", "xgb"])).is_err());
        // classic + 权重 → 通过（忽略并警告）
        let a = parse_args(&argv(&["--gt-weights", "/tmp/m.gtw"])).unwrap();
        assert_eq!(a.path_engine, PathEngine::Classic);
    }

    #[test]
    fn path_engine_rnaspades_parse() {
        // rnaspades 无需额外 flag
        let a = parse_args(&argv(&["--path-engine", "rnaspades"])).unwrap();
        assert_eq!(a.path_engine, PathEngine::RnaSpades);
        // 与 --recursive 互斥（沿 gtasm 先例）
        assert!(parse_args(&argv(&["--path-engine", "rnaspades", "--recursive"])).is_err());
        // rnaspades + 权重 → 通过（忽略并警告）
        let a = parse_args(&argv(&[
            "--path-engine",
            "rnaspades",
            "--gt-weights",
            "/tmp/m.gtw",
        ]))
        .unwrap();
        assert_eq!(a.path_engine, PathEngine::RnaSpades);
    }

    #[test]
    fn multi_k_defaults_off() {
        // 缺省 = None（单 k classic 路径, 行为不变）
        let a = parse_args(&argv(&[])).unwrap();
        assert_eq!(a.multi_k, None);
    }

    #[test]
    fn multi_k_auto_and_list_parse() {
        let a = parse_args(&argv(&["--multi-k", "auto"])).unwrap();
        assert_eq!(a.multi_k, Some(MultiK::Auto));
        let a = parse_args(&argv(&["--multi-k", "21,25,29"])).unwrap();
        assert_eq!(a.multi_k, Some(MultiK::List(vec![21, 25, 29])));
        // 单元素列表合法（退化为该 k 的单轮）
        let a = parse_args(&argv(&["--multi-k", "29"])).unwrap();
        assert_eq!(a.multi_k, Some(MultiK::List(vec![29])));
    }

    #[test]
    fn multi_k_list_validation() {
        // 越界: <15 或 >32 报错
        assert!(parse_args(&argv(&["--multi-k", "14,25"])).is_err());
        assert!(parse_args(&argv(&["--multi-k", "25,33"])).is_err());
        assert!(parse_args(&argv(&["--multi-k", "15,32"])).is_ok());
        // 非数字 / 空段
        assert!(parse_args(&argv(&["--multi-k", "x"])).is_err());
        assert!(parse_args(&argv(&["--multi-k", "21,,25"])).is_err());
        // 非递升 / 重复
        assert!(parse_args(&argv(&["--multi-k", "25,21"])).is_err());
        assert!(parse_args(&argv(&["--multi-k", "25,25"])).is_err());
    }

    #[test]
    fn multi_k_mutex_with_kmer_size_and_recursive() {
        // 与 --KMER_SIZE（含原版别名 __KMER_SIZE）同给 → 报错
        assert!(parse_args(&argv(&["--multi-k", "auto", "--KMER_SIZE", "25"])).is_err());
        assert!(parse_args(&argv(&["--multi-k", "21,25", "__KMER_SIZE", "25"])).is_err());
        // 与 --recursive 同给 → 报错（auto / list 同规则）
        assert!(parse_args(&argv(&["--multi-k", "auto", "--recursive"])).is_err());
        assert!(parse_args(&argv(&["--multi-k", "21,25", "--recursive"])).is_err());
        // 单独 --recursive / 单独 --multi-k 均不受影响
        assert!(parse_args(&argv(&["--recursive"])).is_ok());
        assert!(parse_args(&argv(&["--multi-k", "21,25"])).is_ok());
    }
}
