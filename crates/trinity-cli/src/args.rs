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

/// `--multi-k`（**Rust 版扩展**, rnaSPAdes 式多 k 迭代 + trusted contig 注入）。
/// **缺省（无任何相关 flag）= Auto**（2026-08-30 转正为默认: 三数据集交叉验证
/// 为净增益, 见 docs/opt-log.md「multi-k 转正为默认」）; `--no-multi-k` /
/// 单独 `--KMER_SIZE` / recursive 特例回退为 None（单 k, 见
/// [`resolve_multi_k`]）。
///
/// - `Auto`: 按 reads 读长（前 1000 条抽样）走 rnaSPAdes `rna_k_values` 公式
///   （multik::rna_k_values）, 再钳入支持域 [15, 63]（KmerId=u128 扩容后,
///   RL=76 的公式原值 [29,37] 直接落地）;
/// - `List`: 显式递升 k 列表（逗号分隔, 每个 k ∈ [15, 63]）。
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
    /// 原版 `__KMER_SIZE` 默认 25（本实现暴露为 `--KMER_SIZE`）; 1..=63
    /// （KmerId=u128 扩容, 上限 63 = MAX_KMER_LENGTH-1: `(1 << 2k)` 族表达式
    /// 恒落在位宽内, 不依赖满宽边界的自限回绕; 原版 u64 时代为 1..=32）。
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
    /// **Rust 版扩展**: rnaSPAdes 式多 k 迭代的生效值（`--multi-k auto|k1,...` /
    /// `--no-multi-k` 经 [`resolve_multi_k`] 解析后的结果）。None = 单 k classic
    /// （k = --KMER_SIZE）; 缺省（无任何相关 flag）= Some(Auto)（2026-08-30 转正）。
    /// 显式 --multi-k 与 --KMER_SIZE / --no-multi-k / --recursive 互斥（防歧义）。
    pub multi_k: Option<MultiK>,
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
/// [multik::MULTI_K_MIN, multik::MULTI_K_MAX]（15..=63, KmerId=u128 支持域）,
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

/// multi-k 生效值解析（2026-08-30 转正为默认, 三数据集净增益见
/// docs/opt-log.md「多 k 增益三数据集交叉验证」）。返回（生效值, 警告列表）——
/// 纯函数不打印, 警告由 parse_args 统一 eprintln（便于测试断言文案）。
///
/// 优先级（explicit = 显式 --multi-k 值, 已解析）:
/// 1. 显式 `--multi-k auto|list` → 照用。
/// 2. `--no-multi-k` → 单 k（与显式 --multi-k 同给 = 矛盾, usage 报错）。
/// 3. 单独 `--KMER_SIZE`（无 --multi-k/--no-multi-k）→ 隐式单 k + 提示警告。
/// 4. 缺省: --recursive → 回退单 k + 警告（迷你管线未接线; 显式组合仍是报错）;
///    其余 → Some(Auto)。
fn resolve_multi_k(
    explicit: Option<MultiK>,
    no_multi_k: bool,
    kmer_given: bool,
    kmer_size: usize,
    recursive: bool,
) -> Result<(Option<MultiK>, Vec<String>), CliError> {
    if no_multi_k && explicit.is_some() {
        return Err(CliError::usage(
            "Error, --no-multi-k 与 --multi-k 互斥: 前者关闭多 k 迭代, 后者显式开启, 语义矛盾",
        ));
    }
    match explicit {
        Some(mk) => Ok((Some(mk), Vec::new())),
        None => {
            if no_multi_k {
                return Ok((None, Vec::new()));
            }
            if kmer_given {
                return Ok((
                    None,
                    vec![format!(
                        "WARNING: 已给出 --KMER_SIZE {kmer_size} 但未给 --multi-k——按隐式单 k 模式运行（多 k 已是默认行为）; 如需多 k 请显式 --multi-k auto|k1,k2,..."
                    )],
                ));
            }
            if recursive {
                return Ok((
                    None,
                    vec![
                        "WARNING: --multi-k 与 --recursive 暂不支持组合（迷你管线未接线）, 已回退单 k（--KMER_SIZE 可选定 k 值）".to_string(),
                    ],
                ));
            }
            Ok((Some(MultiK::Auto), Vec::new()))
        }
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
    // --multi-k auto | k1,k2,...（Rust 版扩展, rnaSPAdes 式多 k 迭代; 显式值）
    let multi_k_explicit = match cli.value_flag("--multi-k")?.as_deref() {
        None => None,
        Some("auto") => Some(MultiK::Auto),
        Some(v) => Some(MultiK::List(parse_multi_k_list(v)?)),
    };
    // --no-multi-k（沿 --no_normalize_reads 命名先例）: 回退 classic 单 k
    let no_multi_k = cli.bool_flag("--no-multi-k");
    let kmer_given = cli.was_given("--KMER_SIZE") || cli.was_given("__KMER_SIZE");

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
    // u128 扩容: 支持域上限与 multik::MULTI_K_MAX 一致（63 = MAX_KMER_LENGTH-1,
    // 见 multik.rs 注记; 原版 u64 时代为 1..=32）
    if !(1..=crate::multik::MULTI_K_MAX).contains(&kmer_size) {
        return Err(CliError::usage(format!(
            "invalid --KMER_SIZE value: {kmer_size} (must be 1..={})",
            crate::multik::MULTI_K_MAX
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

    // multi-k 生效值解析（默认化, 2026-08-30 转正）: 缺省 = auto; --no-multi-k /
    // 单独 --KMER_SIZE / recursive 特例 → 单 k（回退均带 WARNING, 让用户对生效
    // 模式无歧义）。警告文案由 resolve_multi_k 以返回值给出, 此处统一打印
    // （纯函数便于测试断言文案）。
    let (multi_k, mk_warns) = resolve_multi_k(
        multi_k_explicit,
        no_multi_k,
        kmer_given,
        kmer_size,
        recursive,
    )?;
    for w in &mk_warns {
        eprintln!("{w}");
    }
    // multi-k 互斥校验（防歧义）: 多 k 轮次自带 k 序列, 与 --KMER_SIZE 同给
    // 语义冲突（哪个生效?）; 与 --recursive 的组合未实现。
    // 解析后 multi_k 仍 Some ⇒ 用户显式给了 --multi-k（缺省 auto 在 recursive
    // 下已被 resolve 回退为单 k, 不误伤纯 --recursive 用法）。
    if multi_k.is_some() {
        if kmer_given {
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
        multi_k,
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
    fn scratch_parses_dir() {
        let a = parse_args(&argv(&["--scratch", "/tmp/scratch"])).unwrap();
        assert_eq!(a.scratch, Some(PathBuf::from("/tmp/scratch")));
    }

    #[test]
    fn multi_k_defaults_to_auto() {
        // 缺省（无任何相关 flag）= --multi-k auto（2026-08-30 转正）
        let a = parse_args(&argv(&[])).unwrap();
        assert_eq!(a.multi_k, Some(MultiK::Auto));
    }

    #[test]
    fn no_multi_k_forces_single_k() {
        // --no-multi-k → 单 k classic（k = --KMER_SIZE 值, 默认 25）
        let a = parse_args(&argv(&["--no-multi-k"])).unwrap();
        assert_eq!(a.multi_k, None);
        assert_eq!(a.kmer_size, 25);
        // --no-multi-k + --KMER_SIZE → 单 k 取给定值, 无警告
        let a = parse_args(&argv(&["--no-multi-k", "--KMER_SIZE", "31"])).unwrap();
        assert_eq!(a.multi_k, None);
        assert_eq!(a.kmer_size, 31);
        // 与显式 --multi-k 同给 = 矛盾 → usage 报错
        assert!(parse_args(&argv(&["--no-multi-k", "--multi-k", "auto"])).is_err());
        assert!(parse_args(&argv(&["--no-multi-k", "--multi-k", "21,25"])).is_err());
    }

    #[test]
    fn kmer_size_alone_implicit_single_k_with_warning() {
        // 单独 --KMER_SIZE → 隐式单 k（既有互斥只挡显式 --multi-k 组合）
        let a = parse_args(&argv(&["--KMER_SIZE", "31"])).unwrap();
        assert_eq!(a.multi_k, None);
        assert_eq!(a.kmer_size, 31);
        assert!(parse_args(&argv(&["__KMER_SIZE", "31"])).unwrap().multi_k.is_none());
        // 警告文案（resolve_multi_k 纯函数断言）: 提示可用 --multi-k
        let (mk, warns) = resolve_multi_k(None, false, true, 31, false).unwrap();
        assert_eq!(mk, None);
        assert_eq!(warns.len(), 1);
        assert!(
            warns[0].starts_with("WARNING") && warns[0].contains("--multi-k"),
            "{}",
            warns[0]
        );
    }

    #[test]
    fn recursive_default_downgrades_to_single_k() {
        // --recursive 未显式 --multi-k → 回退单 k（不报错, 保既有用法）;
        // 显式组合仍是报错（见 multi_k_mutex_with_kmer_size_and_recursive）
        let a = parse_args(&argv(&["--recursive"])).unwrap();
        assert_eq!(a.multi_k, None);
        let (_, warns) = resolve_multi_k(None, false, false, 25, true).unwrap();
        assert_eq!(warns.len(), 1);
        assert!(warns[0].contains("--recursive"), "{}", warns[0]);
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
        // 越界: <15 或 >63 报错（u128 支持域 [15,63]; 33 扩容后合法）
        assert!(parse_args(&argv(&["--multi-k", "14,25"])).is_err());
        assert!(parse_args(&argv(&["--multi-k", "25,64"])).is_err());
        assert!(parse_args(&argv(&["--multi-k", "15,32"])).is_ok());
        assert!(parse_args(&argv(&["--multi-k", "15,63"])).is_ok());
        // RL=76 公式原值 [29,37] 作为显式列表合法（扩容动机）
        assert!(parse_args(&argv(&["--multi-k", "29,37"])).is_ok());
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
