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
///   trinity-butterfly `gtmodel`/`gtdecode` 模块）; 只在 k=25 图上验证过——
///   未显式 `--multi-k` 时强制单 k（[`resolve_multi_k`] 特例）;
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
/// **缺省（无任何相关 flag）= Auto**（2026-08-30 转正为默认: 三数据集交叉验证
/// 为净增益, 见 docs/opt-log.md「multi-k 转正为默认」）; `--no-multi-k` /
/// 单独 `--KMER_SIZE` / gtasm·recursive 特例回退为 None（单 k, 见
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
    /// （KmerId=u128 扩容, 上限与 multik::MULTI_K_MAX 同一支持域）。
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
    /// **Rust 版扩展**: gtasm 解码的 visited 穿越（walk 可穿越已 visited 节点,
    /// 两级候选过滤 + 完全重复产出丢弃, 机制见 trinity-butterfly `gtdecode`
    /// 模块文档）。**默认开**（实验默认, 修复共享外显子被强势 isoform 先占
    /// 导致弱势 isoform 截断的缺口）; `--no-gtasm-visited-walk` 回退原型语义
    /// （沿 --no-multi-k 命名先例）。仅 gtasm 引擎生效。
    pub gtasm_visited_walk: bool,
    /// **Rust 版扩展**: gtasm 解码的简单 bulge 双臂产出（前向延伸遇简单 bulge
    /// 分叉点双臂各走一条、胜出起点的全部 fork 参与产出, 机制见 gtdecode
    /// 模块文档）。**默认关**; `--gtasm-bulge-dual` 开启。仅 gtasm 引擎生效。
    pub gtasm_bulge_dual: bool,
    /// **Rust 版扩展**: gtasm 分支决策的 read pair 证据（乘性: gsc ×
    /// (1+log1p(pair)); 支持 pass 产物 gtasm_pair_support.bin, 机制见
    /// gt_support/gtdecode 模块文档）。**默认关**; `--gtasm-pair-support`
    /// 开启。仅 gtasm 引擎、双端 reads（--left/--right）生效。
    pub gtasm_pair_support: bool,
    /// **Rust 版扩展**: butterfly classic 引擎的 cd-hit 式去冗余开关
    /// （**默认关** = 原版语义; `--bfly-no-path-merging` 开启则跳过两处
    /// twoPathsAreTooSimilar 去冗余——同基因高相似 isoform 不再被删）。
    pub bfly_no_path_merging: bool,
    /// **Rust 版扩展**: butterfly 搜索 0 产出分量的 R=1 定向豁免
    /// （**默认关**; `--bfly-empty-retry` 开启——381 档定向修复）。
    pub bfly_empty_retry: bool,
    /// **Rust 版扩展**: gtasm 组件产出的碎片 pair-guided 合并（组件产出
    /// contigs 后、入 harvest 前, 同组件内 contig 末端 × 末端的 mate 链接
    /// ≥2 对 → OLC/N 保守拼接, 机制见 frag_merge 模块文档）。**默认关**;
    /// `--gtasm-frag-merge` 开启。仅 gtasm 引擎、双端 reads（--left/--right）
    /// 生效。
    pub gtasm_frag_merge: bool,
    /// **Rust 版扩展**: 输出层 split-read 缝合（harvest 之后、Trinity.fasta
    /// 落盘前: contig 末端锚 × 首端锚在同一条 read 内共现 ≥2 → 链式缝合,
    /// 锚唯一性 + 端唯一化 + 并查集禁环防误连, 机制见 stitch 模块文档）。
    /// **默认关**; `--stitch` 开启。
    pub stitch: bool,
    /// **Rust 版扩展**: rnaSPAdes 式多 k 迭代的生效值（`--multi-k auto|k1,...` /
    /// `--no-multi-k` 经 [`resolve_multi_k`] 解析后的结果）。None = 单 k classic
    /// （k = --KMER_SIZE）; 缺省（无任何相关 flag）= Some(Auto)（2026-08-30 转正）。
    /// 显式 --multi-k 与 --KMER_SIZE / --no-multi-k / --recursive 互斥（防歧义,
    /// 沿 gtasm×recursive 先例）。
    pub multi_k: Option<MultiK>,
    /// **Rust 版扩展**: rnaSPAdes 式相对覆盖图清理（chrysalis QuantifyGraph
    /// 后对 c*.graph.out 做 tip/EC/rcc 删除, 见 trinity-chrysalis::graph_clean）。
    /// 默认关 = 现有行为逐字节不变。--recursive 下暂不生效（迷你管线未接线）。
    pub graph_clean_rnaspades: bool,
    /// **Rust 版扩展**: ReadsToTranscripts 的 25-mer 索引多值化（共享 kmer
    /// 记全部组件, 原版为后写者赢单值语义）。救跨 weld 口 read 被其他组件
    /// 抢占 → 组件图断口无桥（QM4 缺失解剖 2026-09-03）。**2026-09-04 转正
    /// 默认开**（dev1/dev3 零回退 + QM4 R+0.3/召回+99/230 组件进图）;
    /// `--no-rtt-multi-kmer` 回退原版单值语义。
    pub rtt_multi_kmer: bool,
    /// **Rust 版扩展**: GraphFromFasta weld 证据因子（--weld-glue-factor,
    /// 默认 0.05 = 原版）。聚类割裂档放宽方向（min_glue 随 cov 比例下调）。
    pub weld_glue_factor: f64,
    /// **Rust 版扩展**: weld 最小 read 支持（--weld-min-glue）。**2026-09-04
    /// 默认 auto（哨兵 0, orchestrate 处 resolve: 读长 ≥150 → 1/mg1, 短读 →
    /// 2/原版——dev1 76bp 上单 read 证据 P−1.7 系噪声）**; 显式值覆盖。
    pub weld_min_glue: u32,
    /// **Rust 版扩展**: BubbleUpClustering 最大簇大小（--cluster-max-size,
    /// 默认 25 = 原版 -max_cluster_size）。
    pub cluster_max_size: usize,
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
/// 1. 显式 `--multi-k auto|list` → 照用; gtasm 引擎下附「组合未验证」警告
///    （显式指定 = 用户强制, 放行）。
/// 2. `--no-multi-k` → 单 k（与显式 --multi-k 同给 = 矛盾, usage 报错）。
/// 3. 单独 `--KMER_SIZE`（无 --multi-k/--no-multi-k）→ 隐式单 k + 提示警告。
/// 4. 缺省: gtasm 引擎 → 强制单 k + 回退警告（多 k 组合仅在 k=25 图上验证过）;
///    --recursive → 回退单 k + 警告（迷你管线未接线; 显式组合仍是报错）;
///    其余 → Some(Auto)。
fn resolve_multi_k(
    explicit: Option<MultiK>,
    no_multi_k: bool,
    kmer_given: bool,
    kmer_size: usize,
    engine: PathEngine,
    recursive: bool,
) -> Result<(Option<MultiK>, Vec<String>), CliError> {
    if no_multi_k && explicit.is_some() {
        return Err(CliError::usage(
            "Error, --no-multi-k 与 --multi-k 互斥: 前者关闭多 k 迭代, 后者显式开启, 语义矛盾",
        ));
    }
    match explicit {
        Some(mk) => {
            let warns = if engine == PathEngine::Gtasm {
                vec![
                    "WARNING: gtasm 引擎的多 k 组合未验证（gtasm 只在 k=25 图上验证过）——按显式 --multi-k 强制执行".to_string(),
                ]
            } else {
                Vec::new()
            };
            Ok((Some(mk), warns))
        }
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
            if engine == PathEngine::Gtasm {
                return Ok((
                    None,
                    vec![
                        "WARNING: gtasm 引擎的多 k 组合未验证, 已回退单 k; 显式 --multi-k 可强制组合".to_string(),
                    ],
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
    // --no-gtasm-visited-walk（沿 --no-multi-k 命名先例）: 回退原型 visited 硬排除
    let no_gtasm_visited_walk = cli.bool_flag("--no-gtasm-visited-walk");
    // --gtasm-bulge-dual（Rust 版扩展, 默认关）: 简单 bulge 分叉点双臂各走一条
    let gtasm_bulge_dual = cli.bool_flag("--gtasm-bulge-dual");
    // --gtasm-pair-support（Rust 版扩展, 默认关）: read pair 分支证据乘性消费
    let gtasm_pair_support = cli.bool_flag("--gtasm-pair-support");
    // --bfly-no-path-merging（Rust 版扩展, 默认关）: 跳过 butterfly 两处
    // cd-hit 式去冗余（twoPathsAreTooSimilar）
    let bfly_no_path_merging = cli.bool_flag("--bfly-no-path-merging");
    // --bfly-empty-retry（Rust 版扩展, 默认关）: 0 产出分量 R=1 定向豁免
    let bfly_empty_retry = cli.bool_flag("--bfly-empty-retry");
    // --gtasm-frag-merge（Rust 版扩展, 默认关）: 组件产出碎片 pair-guided 合并
    let gtasm_frag_merge = cli.bool_flag("--gtasm-frag-merge");
    // --stitch（Rust 版扩展, 默认关）: 输出层 split-read 缝合
    let stitch = cli.bool_flag("--stitch");
    // --multi-k auto | k1,k2,...（Rust 版扩展, rnaSPAdes 式多 k 迭代; 显式值）
    let multi_k_explicit = match cli.value_flag("--multi-k")?.as_deref() {
        None => None,
        Some("auto") => Some(MultiK::Auto),
        Some(v) => Some(MultiK::List(parse_multi_k_list(v)?)),
    };
    // --no-multi-k（沿 --no_normalize_reads 命名先例）: 回退 classic 单 k
    let no_multi_k = cli.bool_flag("--no-multi-k");
    let kmer_given = cli.was_given("--KMER_SIZE") || cli.was_given("__KMER_SIZE");
    // --graph-clean-rnaspades（Rust 版扩展, 默认关）
    let graph_clean_rnaspades = cli.bool_flag("--graph-clean-rnaspades");
    // --rtt-multi-kmer（Rust 版扩展, 2026-09-04 转正默认开）: RTT 25-mer 索引
    // 多值化; 旧旗标保留为幂等兼容, --no-rtt-multi-kmer 回退原版单值语义
    let _rtt_multi_kmer_on = cli.bool_flag("--rtt-multi-kmer");
    let rtt_multi_kmer = !cli.bool_flag("--no-rtt-multi-kmer");
    // 聚类参数（Rust 版扩展; weld-min-glue 默认 auto 哨兵 0, orchestrate 按
    // 读长 resolve: RL≥150 → 1/mg1, 短读 → 2/原版; 显式值覆盖）
    let weld_glue_factor = cli.float_flag("--weld-glue-factor", 0.05f64)?;
    let weld_min_glue = cli.uint_flag("--weld-min-glue", 0u32)?;
    let cluster_max_size = cli.uint_flag("--cluster-max-size", 25usize)?;
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
    // U256 扩容: 支持域上限与 multik::MULTI_K_MAX 一致（127 = MAX_KMER_LENGTH-1,
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

    // gtasm 必须给权重; 非 gtasm 给了权重也只警告（不误判为 gtasm）。
    if path_engine == PathEngine::Gtasm && gt_weights.is_none() {
        return Err(CliError::usage(
            "Error, --path-engine gtasm requires --gt-weights <GTW1 权重文件路径>",
        ));
    }
    if path_engine != PathEngine::Gtasm && gt_weights.is_some() {
        eprintln!("WARNING: --gt-weights 仅在 --path-engine gtasm 下生效, 忽略");
    }
    // visited 穿越默认开（实验默认）; 显式 --no-gtasm-visited-walk 回退原型语义。
    // 非 gtasm 引擎给了关闭旗标也只警告（沿 --gt-weights 先例）。
    if path_engine != PathEngine::Gtasm && no_gtasm_visited_walk {
        eprintln!("WARNING: --no-gtasm-visited-walk 仅在 --path-engine gtasm 下生效, 忽略");
    }
    let gtasm_visited_walk = !no_gtasm_visited_walk;
    // bulge 双臂默认关; 非 gtasm 引擎给了旗标也只警告（沿 --gt-weights 先例）。
    if path_engine != PathEngine::Gtasm && gtasm_bulge_dual {
        eprintln!("WARNING: --gtasm-bulge-dual 仅在 --path-engine gtasm 下生效, 忽略");
    }
    // pair 支持同理（单端 reads 的运行期跳过警告在 orchestrate 发出）。
    if path_engine != PathEngine::Gtasm && gtasm_pair_support {
        eprintln!("WARNING: --gtasm-pair-support 仅在 --path-engine gtasm 下生效, 忽略");
    }
    // 碎片合并同理（单端 reads 的运行期跳过警告在 orchestrate 发出）。
    if path_engine != PathEngine::Gtasm && gtasm_frag_merge {
        eprintln!("WARNING: --gtasm-frag-merge 仅在 --path-engine gtasm 下生效, 忽略");
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
    // multi-k 生效值解析（默认化, 2026-08-30 转正）: 缺省 = auto; --no-multi-k /
    // 单独 --KMER_SIZE / gtasm·recursive 特例 → 单 k（回退均带 WARNING,
    // 让用户对生效模式无歧义）。警告文案由 resolve_multi_k 以返回值给出,
    // 此处统一打印（纯函数便于测试断言文案）。
    let (multi_k, mk_warns) = resolve_multi_k(
        multi_k_explicit,
        no_multi_k,
        kmer_given,
        kmer_size,
        path_engine,
        recursive,
    )?;
    for w in &mk_warns {
        eprintln!("{w}");
    }
    // multi-k 互斥校验（防歧义）: 多 k 轮次自带 k 序列, 与 --KMER_SIZE 同给
    // 语义冲突（哪个生效?）; 与 --recursive 的组合未实现（沿 gtasm 先例）。
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
        gtasm_visited_walk,
        gtasm_bulge_dual,
        gtasm_pair_support,
        bfly_no_path_merging,
        bfly_empty_retry,
        gtasm_frag_merge,
        stitch,
        multi_k,
        graph_clean_rnaspades,
        rtt_multi_kmer,
        weld_glue_factor,
        weld_min_glue,
        cluster_max_size,
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

    /// 2026-09-04 转正: rtt-multi-kmer 默认开 + weld-min-glue 默认 auto（哨兵 0）。
    #[test]
    fn rtt_multi_kmer_and_weld_min_glue_promoted_defaults() {
        // 默认: rtt-multi-kmer 开, weld-min-glue = 0（auto, orchestrate 按读长 resolve）
        let a = parse_args(&argv(&[])).unwrap();
        assert!(a.rtt_multi_kmer);
        assert_eq!(a.weld_min_glue, 0);
        // 旧旗标保留为幂等兼容（不报 "do not understand option"）
        let a = parse_args(&argv(&["--rtt-multi-kmer"])).unwrap();
        assert!(a.rtt_multi_kmer);
        // --no-rtt-multi-kmer 回退原版单值语义
        let a = parse_args(&argv(&["--no-rtt-multi-kmer"])).unwrap();
        assert!(!a.rtt_multi_kmer);
        // 显式 weld-min-glue 覆盖 auto（含显式 2 = 原版）
        let a = parse_args(&argv(&["--weld-min-glue", "1"])).unwrap();
        assert_eq!(a.weld_min_glue, 1);
        let a = parse_args(&argv(&["--weld-min-glue", "2"])).unwrap();
        assert_eq!(a.weld_min_glue, 2);
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
    fn gtasm_visited_walk_defaults_on() {
        // 默认开（gtasm 与 classic 皆然——字段恒有值, 仅 gtasm 引擎消费）
        let a = parse_args(&argv(&[
            "--path-engine",
            "gtasm",
            "--gt-weights",
            "/tmp/m.gtw",
        ]))
        .unwrap();
        assert!(a.gtasm_visited_walk);
        // --no-gtasm-visited-walk → 回退原型语义
        let a = parse_args(&argv(&[
            "--path-engine",
            "gtasm",
            "--gt-weights",
            "/tmp/m.gtw",
            "--no-gtasm-visited-walk",
        ]))
        .unwrap();
        assert!(!a.gtasm_visited_walk);
        // classic + 关闭旗标 → 通过（忽略并警告, 沿 --gt-weights 先例）
        let a = parse_args(&argv(&["--no-gtasm-visited-walk"])).unwrap();
        assert!(!a.gtasm_visited_walk);
    }

    #[test]
    fn gtasm_bulge_dual_flag() {
        // 默认关; --gtasm-bulge-dual 开; 非 gtasm 引擎给旗标 → 通过（忽略并警告）
        let a = parse_args(&argv(&[
            "--path-engine",
            "gtasm",
            "--gt-weights",
            "/tmp/m.gtw",
        ]))
        .unwrap();
        assert!(!a.gtasm_bulge_dual);
        let a = parse_args(&argv(&[
            "--path-engine",
            "gtasm",
            "--gt-weights",
            "/tmp/m.gtw",
            "--gtasm-bulge-dual",
        ]))
        .unwrap();
        assert!(a.gtasm_bulge_dual);
        let a = parse_args(&argv(&["--gtasm-bulge-dual"])).unwrap();
        assert!(a.gtasm_bulge_dual); // 字段恒有值, 仅 gtasm 引擎消费
    }

    #[test]
    fn gtasm_pair_support_flag() {
        // 默认关; --gtasm-pair-support 开; 非 gtasm 引擎给旗标 → 通过（忽略并警告）
        let a = parse_args(&argv(&[
            "--path-engine",
            "gtasm",
            "--gt-weights",
            "/tmp/m.gtw",
        ]))
        .unwrap();
        assert!(!a.gtasm_pair_support);
        let a = parse_args(&argv(&[
            "--path-engine",
            "gtasm",
            "--gt-weights",
            "/tmp/m.gtw",
            "--gtasm-pair-support",
        ]))
        .unwrap();
        assert!(a.gtasm_pair_support);
        let a = parse_args(&argv(&["--gtasm-pair-support"])).unwrap();
        assert!(a.gtasm_pair_support);
    }

    #[test]
    fn gtasm_frag_merge_flag() {
        // 默认关; --gtasm-frag-merge 开; 非 gtasm 引擎给旗标 → 通过（忽略并警告）
        let a = parse_args(&argv(&[
            "--path-engine",
            "gtasm",
            "--gt-weights",
            "/tmp/m.gtw",
        ]))
        .unwrap();
        assert!(!a.gtasm_frag_merge);
        let a = parse_args(&argv(&[
            "--path-engine",
            "gtasm",
            "--gt-weights",
            "/tmp/m.gtw",
            "--gtasm-frag-merge",
        ]))
        .unwrap();
        assert!(a.gtasm_frag_merge);
        let a = parse_args(&argv(&["--gtasm-frag-merge"])).unwrap();
        assert!(a.gtasm_frag_merge); // 字段恒有值, 仅 gtasm 引擎消费
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
        let (mk, warns) =
            resolve_multi_k(None, false, true, 31, PathEngine::Classic, false).unwrap();
        assert_eq!(mk, None);
        assert_eq!(warns.len(), 1);
        assert!(
            warns[0].starts_with("WARNING") && warns[0].contains("--multi-k"),
            "{}",
            warns[0]
        );
    }

    #[test]
    fn gtasm_engine_special_case() {
        // gtasm + 未显式 --multi-k → 强制单 k + 回退警告
        let a = parse_args(&argv(&["--path-engine", "gtasm", "--gt-weights", "/tmp/m.gtw"]))
            .unwrap();
        assert_eq!(a.multi_k, None);
        let (mk, warns) =
            resolve_multi_k(None, false, false, 25, PathEngine::Gtasm, false).unwrap();
        assert_eq!(mk, None);
        assert!(
            warns.iter().any(|w| w.contains("gtasm") && w.contains("回退")),
            "{warns:?}"
        );
        // 显式 --multi-k + gtasm → 允许 + 未验证警告
        let a = parse_args(&argv(&[
            "--path-engine",
            "gtasm",
            "--gt-weights",
            "/tmp/m.gtw",
            "--multi-k",
            "29,37",
        ]))
        .unwrap();
        assert_eq!(a.multi_k, Some(MultiK::List(vec![29, 37])));
        let (mk, warns) =
            resolve_multi_k(Some(MultiK::Auto), false, false, 25, PathEngine::Gtasm, false)
                .unwrap();
        assert_eq!(mk, Some(MultiK::Auto));
        assert!(
            warns.iter().any(|w| w.contains("gtasm") && w.contains("未验证")),
            "{warns:?}"
        );
        // classic / rnaspades 引擎不受此限: 缺省仍 auto
        let a = parse_args(&argv(&["--path-engine", "rnaspades"])).unwrap();
        assert_eq!(a.multi_k, Some(MultiK::Auto));
        // gtasm + --no-multi-k → 单 k 无警告（用户显式选择, 非回退）
        let (mk, warns) =
            resolve_multi_k(None, true, false, 25, PathEngine::Gtasm, false).unwrap();
        assert_eq!(mk, None);
        assert!(warns.is_empty());
    }

    #[test]
    fn recursive_default_downgrades_to_single_k() {
        // --recursive 未显式 --multi-k → 回退单 k（不报错, 保既有用法）;
        // 显式组合仍是报错（见 multi_k_mutex_with_kmer_size_and_recursive）
        let a = parse_args(&argv(&["--recursive"])).unwrap();
        assert_eq!(a.multi_k, None);
        let (_, warns) =
            resolve_multi_k(None, false, false, 25, PathEngine::Classic, true).unwrap();
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
        assert!(parse_args(&argv(&["--multi-k", "25,128"])).is_err());
        assert!(parse_args(&argv(&["--multi-k", "15,32"])).is_ok());
        assert!(parse_args(&argv(&["--multi-k", "15,127"])).is_ok());
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
