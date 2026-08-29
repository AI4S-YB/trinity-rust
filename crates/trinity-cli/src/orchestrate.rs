//! 编排主线（P5-T3）——原版 `Trinity` 主流程（L1376-1990）上移:
//! 归一化 → prep → k-mer 计数 → inchworm → chrysalis → butterfly 组件池
//! → 汇总输出。每阶段 checkpoint 沿用原版同名文件（chdir outdir 语义,
//! 这里统一以绝对路径拼出、落在 outdir 下）。
//!
//! 与原版的粒度差异（详见 P5 计划「移植契约」）:
//! - chrysalis 整段用**一个** checkpoint `outdir/.quantify_graph.ok` 包裹
//!   （T1 pipeline 内部六子阶段未落原版各自的 `.ok`——welds/sorted/
//!   components/bundle/rtc/debruijn/partition 粒度在此不可分）;
//! - jellyfish histo 与其 `.ok` 未复刻（无下游消费, 记录于 concerns）;
//! - jellyfish 计数落盘为自定义二进制 `mer_counts.<k>.asm.rs.bin`
//!   （count 与 dump 两个 checkpoint 之间传递, 原版 .jf 在 dump 后即删）。

use std::fs;
use std::path::{Path, PathBuf};

use trinity_chrysalis::pipeline::{run_chrysalis_pipeline, ChrysalisParams};
use trinity_common::error::CommonError;
use trinity_inchworm::irke::{
    compute_sequence_assemblies_parallel, populate_from_table, prune_some_kmers,
    write_kmer_count_report, AssemblyParams, IrkeParams, Monitor,
};
use trinity_kmer::counter::{estimate_hash_size, KmerCountTable};
use trinity_kmer::diginorm::{self, DigiNormParams, ReadsInput};

use crate::args::{MultiK, PathEngine, SeqType, TrinityArgs};
use crate::butterfly_pool::{run_butterfly_pool, ButterflyPoolParams, GtasmShared};
use crate::checkpoint::{checkpoint_exists, run_with_checkpoint};
use crate::harvest::harvest;
use crate::prep::prep_reads;
use crate::recursive::{harvest_mini, run_partition_pool};

/// 原版 L2556 `SR_flag`（super_read_mode 关闭时为 "asm"）。
const SR_FLAG: &str = "asm";

// ---------------------------------------------------------------------------
// --max_memory 护栏（P1 审查 I2）
// ---------------------------------------------------------------------------

/// k-mer 哈希单条目的保守字节下限（u64 key + u32 count + 槽位开销）。
pub const KMER_ENTRY_BYTES: u64 = 18;

/// 归一化 / k-mer 计数阶段的**硬护栏**估算（保守代理公式）:
/// 输入文件需整读进内存（`input_bytes`），计数哈希按"文件量级 distinct
/// k-mer × 18B 下限"估计——去重后的表通常远小于文件，故取
/// `input_bytes + input_bytes/5`（即 1.2× 输入）作为最低需求代理。
/// 返回 `Some(错误信息)` 表示 `max_memory` 过小（应报错而非 OOM）。
pub fn memory_guard_error(max_memory: u64, input_bytes: u64, stage: &str) -> Option<String> {
    let need = input_bytes + input_bytes / 5;
    if max_memory >= need {
        return None;
    }
    let gib = |b: u64| format!("{:.2}G", b as f64 / (1024.0 * 1024.0 * 1024.0));
    let suggest = (need.next_power_of_two() / (1024 * 1024 * 1024)).max(1);
    Some(format!(
        "Error, --max_memory ({}) is too small for the {stage} stage: input is {} bytes ({}) \
         and the k-mer hash must hold on the order of {} distinct k-mers x {KMER_ENTRY_BYTES} B \
         (conservative floor; after dedup the guard requires 1.2x input = {} = {}). \
         Please increase --max_memory (e.g. --max_memory {suggest}G).",
        gib(max_memory),
        input_bytes,
        gib(input_bytes),
        input_bytes / 2,
        need,
        gib(need),
    ))
}

/// chrysalis / butterfly 阶段的**软护栏**（粗检）: 文件量级 × 1.2 超过
/// max_memory 时向 stderr 警告但继续（决策: diginorm/计数硬、后段软——
/// 后段逐组件处理, 峰值由最大组件决定, 文件量级只是上界代理）。
pub fn warn_memory(max_memory: u64, input_bytes: u64, stage: &str) {
    let need = input_bytes + input_bytes / 5;
    if max_memory < need {
        eprintln!(
            "WARNING: {stage} stage input (~{} bytes) may exceed --max_memory ({max_memory} \
             bytes) at current settings; continuing (per-component processing usually fits). \
             If the run is OOM-killed, increase --max_memory.",
            need
        );
    }
}

/// 输入 reads 文件总字节数（归一化前的原始输入; 元数据缺失按 0 计）。
fn input_files_bytes(args: &TrinityArgs) -> u64 {
    args.left
        .iter()
        .chain(&args.right)
        .chain(&args.single)
        .filter_map(|p| fs::metadata(p).ok())
        .map(|m| m.len())
        .sum()
}

fn file_len(p: &Path) -> u64 {
    fs::metadata(p).map(|m| m.len()).unwrap_or(0)
}

/// 相对 outdir → 绝对路径（原版 create_full_path + chdir 语义）。
fn absolutize(p: &Path) -> PathBuf {
    if p.is_absolute() {
        // 原版 `s/\/+$//`：剥尾部斜杠（`<outdir>.Trinity.fasta` 拼接用）。
        let s = p.to_string_lossy();
        return PathBuf::from(s.trim_end_matches('/'));
    }
    let cwd = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
    cwd.join(p)
}

/// outdir 绝对路径的 FNV-1a 64bit 哈希（默认 scratch 目录名用; 确定性,
/// 同一 --output 重跑命中同一目录 → checkpoint 续跑语义不变）。
fn fnv1a64(s: &str) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for b in s.as_bytes() {
        h = (h ^ *b as u64).wrapping_mul(0x100000001b3);
    }
    h
}

/// 默认 scratch 目录: `${TMPDIR:-/tmp}/trinity-rust-<outdir 路径哈希>`。
pub fn default_scratch_dir(outdir: &Path) -> PathBuf {
    let tmp = std::env::var_os("TMPDIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("/tmp"));
    tmp.join(format!(
        "trinity-rust-{:016x}",
        fnv1a64(&outdir.to_string_lossy())
    ))
}

/// statvfs 可用字节数（bavail × frsize）; 调用失败 → None。
fn avail_bytes(dir: &Path) -> Option<u64> {
    use std::os::unix::ffi::OsStrExt;
    let c = std::ffi::CString::new(dir.as_os_str().as_bytes()).ok()?;
    let mut st: libc::statvfs = unsafe { std::mem::zeroed() };
    // SAFETY: c 为合法 C 字符串, st 指向已分配的 statvfs 结构。
    if unsafe { libc::statvfs(c.as_ptr(), &mut st) } == 0 {
        Some(st.f_bavail as u64 * st.f_frsize.max(1) as u64)
    } else {
        None
    }
}

/// 默认 scratch 的空间护栏: 需求 = 压缩输入总大小 × 6, 下限 4GB。
/// TMPDIR 不可创建/不可写或空间不足 → Err（调用方打印警告并回落 outdir）。
fn check_default_scratch(cand: &Path, need: u64) -> Result<(), String> {
    if let Err(e) = fs::create_dir_all(cand) {
        return Err(format!(
            "默认 scratch 目录 {} 不可创建/不可写（{e}）",
            cand.display()
        ));
    }
    match avail_bytes(cand) {
        None => Err(format!("无法 statvfs 默认 scratch {}", cand.display())),
        Some(avail) if avail < need => {
            let _ = fs::remove_dir(cand); // 刚建的空目录顺手清掉（非空则保留）
            Err(format!(
                "默认 scratch {} 可用空间 {:.1}GB < 需求 {:.1}GB",
                cand.display(),
                avail as f64 / (1 << 30) as f64,
                need as f64 / (1 << 30) as f64
            ))
        }
        Some(_) => Ok(()),
    }
}

fn ext_of(t: SeqType) -> &'static str {
    match t {
        SeqType::Fq => "fq",
        SeqType::Fa => "fa",
    }
}

/// 归一化的 SS 方向规则（insilico_read_normalization.pl L195-204 镜像）:
/// PE RF → left revcomp、FR → right revcomp; 单端 'R' 原版改写为 'F'
/// （不 revcomp）。
fn diginorm_params(args: &TrinityArgs, paired: bool) -> DigiNormParams {
    let ss = args.ss_lib_type.as_deref();
    let (rl, rr) = match ss {
        Some("RF") if paired => (true, false),
        Some("FR") if paired => (false, true),
        // 单端（或无 SS）恒不 revcomp（'R'→'F' 改写）。
        _ => (false, false),
    };
    DigiNormParams {
        ss_revcomp_left: rl,
        ss_revcomp_right: rr,
        ds: ss.is_none(),
        max_cov: args.normalize_max_read_cov as f64,
        min_cov: args.min_kmer_count as f64,
        // IO 优化: fq 输入直出归一化 fa（跳过 norm.fq 物化）; --no_cleanup
        // 时额外保留 norm.fq（原版把归一化 fq 留给用户——既定差异见 README）。
        fa_direct: args.seq_type == SeqType::Fq,
        keep_fq: args.no_cleanup,
        ..Default::default()
    }
}

/// mer_counts 二进制存取已收编至 `trinity_kmer::counter::{save_counts_bin,
/// load_counts_bin}`; 保留本地薄封装以稳调用点。
fn save_counts(
    path: &Path,
    counts: &trinity_kmer::counter::CountMap,
    k: usize,
    ds: bool,
) -> Result<(), CommonError> {
    trinity_kmer::counter::save_counts_bin(path, counts, k, ds)
}

fn load_counts(
    path: &Path,
    k: usize,
    ds: bool,
) -> Result<trinity_kmer::counter::CountMap, CommonError> {
    trinity_kmer::counter::load_counts_bin(path, k, ds)
}

/// component_base_listing.txt（T1 partition 落盘）→ (comp id, base 前缀) 升序。
pub fn read_component_base_listing(p: &Path) -> Result<Vec<(u64, PathBuf)>, CommonError> {
    let text = fs::read_to_string(p)?;
    let mut out = Vec::new();
    for line in text.lines() {
        let line = line.trim_end_matches('\n');
        if line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let id: u64 = f.next().unwrap_or("").parse().map_err(|_| {
            CommonError::Parse(format!("bad component id in {}: {line}", p.display()))
        })?;
        let base = PathBuf::from(f.next().unwrap_or(""));
        out.push((id, base));
    }
    out.sort_by_key(|(id, _)| *id);
    Ok(out)
}

/// Trinity 主线全流程。返回最终 fasta 路径（`<outdir>.Trinity.fasta`）。
pub fn run_trinity(args: &TrinityArgs) -> Result<PathBuf, CommonError> {
    let outdir = absolutize(&args.output);
    // --scratch: 中间产物（归一化目录/both.fa/mer_counts/inchworm/chrysalis/
    // butterfly）与 checkpoint（.ok）全部落 `<scratch>/<output 目录名>/`,
    // 最终 `<output>.Trinity.fasta` 与 gene_trans_map 拷回 --output 旁。
    // checkpoint 随产物走: scratch 易失, 换 scratch 目录 = 全部重算。
    // **默认化**: 不传 --scratch 时中间产物落 `${TMPDIR:-/tmp}/trinity-rust-
    // <outdir 路径哈希>/`（哈希确定 → 同 --output 重跑续跑语义不变）;
    // TMPDIR 不可写或空间不足（粗估: 压缩输入 ×6, 下限 4GB）时警告并回落
    // 旧行为（中间产物直接放 --output）。显式 --scratch 不检查不回落。
    let workdir = match &args.scratch {
        Some(s) => absolutize(s).join(outdir.file_name().unwrap_or_default()),
        None => {
            let cand = default_scratch_dir(&outdir);
            let need = input_files_bytes(args).saturating_mul(6).max(4u64 << 30);
            match check_default_scratch(&cand, need) {
                Ok(()) => cand,
                Err(msg) => {
                    eprintln!(
                        "WARNING: {msg}——中间产物回落到 --output 目录（可用 --scratch 显式指定跳过检查）"
                    );
                    outdir.clone()
                }
            }
        }
    };
    // scratch 模式下也建 outdir: 保证最终产物（outdir 的同级文件）父链存在,
    // 与默认模式"outdir 必被创建"的行为对齐。
    fs::create_dir_all(&outdir)?;
    fs::create_dir_all(&workdir)?;

    let paired = !args.left.is_empty();
    let ss = args.ss_lib_type.as_deref();
    let ds = ss.is_none();
    let mut a = args.clone();

    // rayon 全局池按 --CPU 设一次（局部需不同线程数的入口——inchworm/
    // chrysalis GFF/RTT/quantify/butterfly 池——各自显式 scoped pool 覆盖;
    // `-t 1`/`--CPU 1` 时全局默认即单线程）。
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(args.cpu.max(1))
        .build_global();

    // ---- multi-k 轮次解析（rnaSPAdes 式多 k 迭代 + trusted contig 注入,
    // 机制调研见 docs/opt-log.md「rnaSPAdes 构图机制调研」）:
    // 缺省 None → 单 k classic（k = --KMER_SIZE, 行为逐字节不变）;
    // auto → 读长抽样 + rna_k_values 公式 + 支持域钳制（multik.rs）。
    // k = 最后一轮的 k（单 k 时即唯一轮）。----
    let k_list: Vec<usize> = match &args.multi_k {
        None => vec![args.kmer_size],
        Some(MultiK::List(ks)) => ks.clone(),
        Some(MultiK::Auto) => {
            let rl = crate::multik::detect_read_length(args)?;
            let raw = crate::multik::rna_k_values(rl);
            let ks = crate::multik::clamp_k_list(raw.clone());
            eprintln!(
                "--multi-k auto: 读长 RL={rl}（各文件前 1000 条抽样取 max, 跨文件取 min）→ rnaSPAdes 公式 {raw:?} → 钳制后 {ks:?}"
            );
            ks
        }
    };
    // 恒非空: None→[kmer_size]; List 解析校验非空; auto 公式恒非空。
    let k = *k_list.last().unwrap();

    // ---- 1. 归一化（L1441-1453; checkpoint insilico_read_normalization/normalization.ok）----
    // fa_direct（fq 输入直出归一化 fa）时: 下游按 fa 消费, 且 SS revcomp 已在
    // diginorm 内部 prep 完成 → prep 阶段不得再 revcomp（防双重反转）。
    let norm_fa_direct = !args.no_normalize_reads && args.seq_type == SeqType::Fq;
    if !args.no_normalize_reads {
        let norm: Result<(), CommonError> =
            trinity_common::timer::timed(&workdir, "normalize", || {
                // 硬护栏: 归一化需把两侧读入整读进内存。
                if let Some(msg) = memory_guard_error(
                    args.max_memory,
                    input_files_bytes(args),
                    "read normalization",
                ) {
                    return Err(CommonError::Parse(msg));
                }
                let norm_dir = workdir.join("insilico_read_normalization");
                let norm_ok = norm_dir.join("normalization.ok");
                if !run_with_checkpoint(&norm_ok, "In silico Read Normalization", || {
                    let params = diginorm_params(args, paired);
                    let reads = if paired {
                        ReadsInput::Paired(a.left.clone(), a.right.clone())
                    } else {
                        ReadsInput::Single(a.single.clone())
                    };
                    diginorm::run(&params, &reads, &norm_dir)?;
                    Ok(())
                })? {
                    // skipped: 输出文件应已在（resume 语义）。
                }
                let ext = if norm_fa_direct {
                    "fa"
                } else {
                    ext_of(args.seq_type)
                };
                if paired {
                    a.left = vec![norm_dir.join(format!("left.norm.{ext}"))];
                    a.right = vec![norm_dir.join(format!("right.norm.{ext}"))];
                } else {
                    a.single = vec![norm_dir.join(format!("single.norm.{ext}"))];
                }
                if norm_fa_direct {
                    a.seq_type = SeqType::Fa;
                }
                Ok(())
            });
        norm?;
    }

    // ---- 2. prep → both.fa/single.fa + .read_count（L1604-1705）----
    let (target_fa, read_count) = trinity_common::timer::timed(&workdir, "prep", || {
        prep_reads(&a, &workdir, norm_fa_direct)
    })?;
    let read_count_file = workdir.join(format!(
        "{}.read_count",
        target_fa.file_name().unwrap_or_default().to_string_lossy()
    ));
    if fs::metadata(&read_count_file).map(|m| m.len()).unwrap_or(0) == 0 {
        fs::write(&read_count_file, format!("{read_count}\n"))?;
    }
    let _ = read_count;

    // ---- 多 k 中间构图轮（i < last）: 只跑 kmer 计数 + inchworm, 跳过
    // chrysalis/butterfly（SPAdes 中间轮不做路径解析）; 产物 inchworm.*.fa
    // 为 trusted contig 源。各轮独立子目录 iter_k<k>/（断点/计时与主线隔离）。
    // 链式注入: 第 i 轮（i>0）计数装载时注入第 i-1 轮产物（SPAdes
    // construction.cpp:233 语义——trusted contigs 贡献序列、不贡献覆盖）。----
    let mut trusted: Option<PathBuf> = None;
    let n_rounds = k_list.len();
    for (i, &ki) in k_list.iter().enumerate().take(n_rounds.saturating_sub(1)) {
        let iter_dir = workdir.join(format!("iter_k{ki}"));
        fs::create_dir_all(&iter_dir)?;
        eprintln!(
            "==== multi-k 构图轮 {}/{n_rounds}: K={ki}（跳过 chrysalis/butterfly）====",
            i + 1
        );
        let iworm = run_count_and_inchworm(&iter_dir, &target_fa, ki, ds, args, trusted.as_deref())?;
        // 中间轮产物为空 = 数据对该 k 太稀疏: 警告并降级（本轮无 trusted 可注）,
        // 不走主线的报错中止。
        trusted = if file_len(&iworm) == 0 {
            eprintln!(
                "WARNING: multi-k 轮 K={ki} 的 inchworm 产物为空——本轮无 trusted contigs, 链式注入降级"
            );
            None
        } else {
            Some(iworm)
        };
        // 中间轮的 counts_bin 使命已毕（inchworm 断点已立）, 清掉省盘;
        // resume 时本轮各 .ok 在 → 整轮跳过, 无需 counts_bin。
        if !args.no_cleanup {
            let _ = fs::remove_file(iter_dir.join(format!("mer_counts.{ki}.{SR_FLAG}.rs.bin")));
        }
    }

    // 最后一轮（单 k classic = 唯一轮）: 多 k 时落 iter_k<k>/ 子目录——
    // .quantify_graph.ok / .butterfly.ok 等断点名不含 k, 必须与 classic 主
    // workdir 的同名断点隔离; 单 k 时 round_dir == workdir, 行为逐字节不变。
    let multi = n_rounds > 1;
    let round_dir = if multi {
        let d = workdir.join(format!("iter_k{k}"));
        fs::create_dir_all(&d)?;
        eprintln!("==== multi-k 最终轮 {n_rounds}/{n_rounds}: K={k}（完整管线 + trusted 注入）====");
        d
    } else {
        workdir.clone()
    };
    let chrysalis_dir = round_dir.join("chrysalis");
    fs::create_dir_all(&chrysalis_dir)?;

    // ---- 3/4. jellyfish 计数 → inchworm（L2584-2727）----
    // trusted 注入在 kmer_load 的最终计数表上、prune_some_kmers 之前统一进行
    // （详见 run_count_and_inchworm 文档）。
    let iworm_fa = run_count_and_inchworm(&round_dir, &target_fa, k, ds, args, trusted.as_deref())?;
    // 命名镜像 L1576-1580: SS → inchworm.fa, 无 SS → inchworm.DS.fa。
    let iworm_name = if ss.is_some() {
        "inchworm.fa"
    } else {
        "inchworm.DS.fa"
    };
    if fs::metadata(&iworm_fa).map(|m| m.len()).unwrap_or(0) == 0 {
        // 原版 L1841-1843 NON_FATAL_EXCEPTION（稀疏数据）——这里直接报错,
        // 由调用方决定是否容忍。
        return Err(CommonError::Parse(format!(
            "WARNING, no Inchworm output is detected at: {}",
            iworm_fa.display()
        )));
    }

    // ---- 4.5 min100 过滤（Trinity:2068-2086: filter_iworm_by_min_length_or_cov.pl 100 10）----
    // 主线（phase-1）在 chrysalis 聚类前只保留 len>=100 || cov>=10 的 contig
    // （cov 取自 accession `aN;cov` 的分号第二段）。产物落 chrysalis/ 下,
    // checkpoint 名 = `<basename>.min100.ok`（镜像原版）。
    let iworm_min100 = chrysalis_dir.join(format!("{iworm_name}.min100"));
    let iworm_min100_ok = chrysalis_dir.join(format!("{iworm_name}.min100.ok"));
    run_with_checkpoint(&iworm_min100_ok, "filter iworm min100", || {
        if !checkpoint_exists(&iworm_min100) {
            let data = fs::read(&iworm_fa)?;
            let filtered = filter_iworm_min_len_or_cov(&data, 100, 10)?;
            fs::write(&iworm_min100, filtered)?;
        }
        Ok(())
    })?;

    // ---- 5. chrysalis（整段 checkpoint .quantify_graph.ok; 粒度差异见模块文档）----
    // 软护栏（粗检, 继续执行）。
    warn_memory(
        args.max_memory,
        file_len(&iworm_min100) + file_len(&target_fa),
        "Chrysalis",
    );
    let iworm_data = fs::read(&iworm_min100)?;
    let reads_data = fs::read(&target_fa)?;

    // ---- 5/6/7. 递归模式（--recursive, M1）: phase-1 止于分区 → 迷你分区池
    // ---- → 汇总（与 classic 的 chrysalis+butterfly+harvest 段并列）。 ----
    // 注: multi-k 与 --recursive 在参数层互斥, 此处 round_dir 恒 == workdir。
    if args.recursive {
        return run_recursive(
            a,
            args,
            &outdir,
            &round_dir,
            &chrysalis_dir,
            &iworm_data,
            &reads_data,
        );
    }

    let mut listing_opt: Option<Vec<(u64, PathBuf)>> = None;
    let quantify_ok = round_dir.join(".quantify_graph.ok");
    trinity_common::timer::timed(&round_dir, "chrysalis_total", || {
        run_with_checkpoint(
            &quantify_ok,
            "Chrysalis (clustering & de Bruijn graph)",
            || {
                let listing = run_chrysalis_pipeline(
                    &iworm_data,
                    &reads_data,
                    &chrysalis_dir,
                    &ChrysalisParams {
                        strand: ss.is_some(),
                        min_contig_length: args.min_contig_length,
                        max_reads: args.max_reads_per_graph,
                        threads: args.cpu.max(1),
                        graph_clean_rnaspades: args.graph_clean_rnaspades,
                        ..Default::default()
                    },
                )?;
                listing_opt = Some(listing);
                Ok(())
            },
        )
    })?;
    let listing = match listing_opt {
        Some(l) => l,
        None => read_component_base_listing(
            &chrysalis_dir
                .join("Component_bins")
                .join("component_base_listing.txt"),
        )?,
    };
    if listing.is_empty() {
        return Err(CommonError::Parse(
            "WARNING, component base listing file is empty - likely sparse data".to_string(),
        ));
    }

    // ---- 5.5 gtasm 边支持 pass（仅 --path-engine gtasm; 读 Component_bins
    // 建 k2v, 扫**原始** reads（归一化前）, 产物 gtasm_edge_support.bin;
    // checkpoint .gtasm_support.ok 保护, 重跑直接加载不重扫）----
    let gtasm_shared: Option<std::sync::Arc<GtasmShared>> = if args.path_engine
        == PathEngine::Gtasm
    {
        let weights = args.gt_weights.as_ref().ok_or_else(|| {
            CommonError::Parse("--path-engine gtasm 需要 --gt-weights".to_string())
        })?;
        // 模型先加载（廉价）——坏权重在昂贵的 reads 扫描之前报错
        let model = trinity_butterfly::GtModel::load(weights)
            .map_err(|e| CommonError::Parse(format!("GT 权重加载失败: {e}")))?;
        eprintln!(
            "-gtasm: 模型 {} 加载完成（hidden={}, {} 层）",
            weights.display(),
            model.hidden,
            model.n_layers()
        );
        let sup_bin = round_dir.join("gtasm_edge_support.bin");
        let sup_ok = round_dir.join(".gtasm_support.ok");
        trinity_common::timer::timed(&round_dir, "gtasm_support", || {
            run_with_checkpoint(&sup_ok, "GTasm edge support (raw reads)", || {
                let raw_reads: Vec<PathBuf> = if paired {
                    args.left.iter().chain(args.right.iter()).cloned().collect()
                } else {
                    args.single.clone()
                };
                let st = crate::gt_support::run_support_pass(
                    &listing,
                    &raw_reads,
                    args.seq_type,
                    &sup_bin,
                )?;
                eprintln!(
                    "-gtasm support: reads {}, 有支持转移 {}, 有支持边 {}",
                    st.reads, st.pairs, st.edges_with_support
                );
                Ok(())
            })
        })?;
        let support = trinity_butterfly::EdgeSupport::load(&sup_bin)?;
        Some(std::sync::Arc::new(GtasmShared { model, support }))
    } else {
        None
    };

    // ---- 6. butterfly 组件池（checkpoint .butterfly.ok; L1934-1936）----
    warn_memory(args.max_memory, file_len(&iworm_fa), "Butterfly");
    trinity_common::timer::timed(&round_dir, "butterfly", || {
        run_with_checkpoint(
            &round_dir.join(".butterfly.ok"),
            "Butterfly (component pool)",
            || {
                run_butterfly_pool(
                    &listing,
                    &round_dir,
                    &ButterflyPoolParams {
                        cpu: args.cpu.max(1),
                        min_contig_length: args.min_contig_length,
                        group_pairs_distance: args.group_pairs_distance,
                        stack_size_mb: args.bfly_stack_mb,
                        path_engine: args.path_engine,
                        gtasm: gtasm_shared.clone(),
                    },
                )
            },
        )
    })?;

    // ---- 7. 汇总（无 checkpoint; 每次重写 Trinity.tmp.fasta, L1916/L1943）----
    let trinity_tmp = round_dir.join("Trinity.tmp.fasta");
    let n_transcripts = trinity_common::timer::timed(&round_dir, "harvest", || {
        harvest(&listing, args.min_contig_length, &trinity_tmp)
    })?;
    if n_transcripts == 0 {
        return Err(CommonError::Parse(
            "WARNING: no transcripts harvested from butterfly output".to_string(),
        ));
    }

    // 最终命名: `<outdir 绝对路径>.Trinity.fasta`（L1512）+ gene_trans_map。
    // scratch 模式下 tmp 在 workdir（可能跨设备）→ rename 失败回落 copy+remove。
    let final_fa = PathBuf::from(format!("{}.Trinity.fasta", outdir.display()));
    move_file(&trinity_tmp, &final_fa)?;
    write_gene_trans_map(
        &final_fa,
        &final_fa.with_file_name(format!(
            "{}.gene_trans_map",
            final_fa.file_name().unwrap_or_default().to_string_lossy()
        )),
    )?;

    // ---- 8. 收尾清理（原版 unlink both.fa/left.fa/right.fa 与 jellyfish
    // 中间物; `--no_cleanup` 保留全部）。保守集合: 只删明确知道的中间物,
    // chrysalis/butterfly 产物保留（原版也保留大部分）。.ok 断点保留
    // （prep 侧对"产物缺失但 .ok 在"会重建, resume 语义不受影响）。
    // both.fa/left.fa/right.fa 是各 k 轮共享的 prep 产物（主 workdir）;
    // mer_counts 是本轮的（round_dir; classic 时与 workdir 相同）。
    if !args.no_cleanup {
        let counts_bin = round_dir.join(format!("mer_counts.{k}.{SR_FLAG}.rs.bin"));
        for p in [
            &target_fa,
            &counts_bin,
            &workdir.join("left.fa"),
            &workdir.join("right.fa"),
        ] {
            let _ = fs::remove_file(p);
        }
    }
    Ok(final_fa)
}

/// k-mer 计数 + inchworm 组装（主线阶段 3-4 的抽出; classic 唯一轮与 multi-k
/// 各轮共用同一实现）。断点/产物全部落 `round_dir`（classic: 主 workdir;
/// multi-k: `iter_k<k>/` 子目录）。
///
/// `trusted` = 上一轮 inchworm 产物（multi-k 链式注入; classic 恒 None）:
/// 在 kmer_load 产出的**最终计数表**上注入其全部 k-mer——表中已存在的键不动
/// 计数, 缺席的以 count=1 插入（presence-only, 见 counter::inject_trusted_kmers）。
/// 注入点统一在 `populate_from_table` 之前: 计数表的两条构建路径（counts_bin
/// 断点装载 / target_fa 重算）在此汇合, 注入因而先于 prune_some_kmers 生效,
/// 且对两条路径全覆盖（counts_bin 断点内容保持 reads-only; 注入幂等, resume
/// 走装载路径再注一次结果不变）。
///
/// 返回 inchworm 产物 fasta 路径（可能指向空文件——判定归调用方: 主线报错,
/// multi-k 中间轮降级）。
fn run_count_and_inchworm(
    round_dir: &Path,
    target_fa: &Path,
    k: usize,
    ds: bool,
    args: &TrinityArgs,
    trusted: Option<&Path>,
) -> Result<PathBuf, CommonError> {
    // ---- jellyfish 计数（L2584-2632）----
    let count_ok = round_dir.join(format!(
        ".jellyfish_count.{k}.mincov{}.{SR_FLAG}.ok",
        args.min_kmer_count
    ));
    let counts_bin = round_dir.join(format!("mer_counts.{k}.{SR_FLAG}.rs.bin"));
    trinity_common::timer::timed(round_dir, "kmer_count", || -> Result<(), CommonError> {
        if !checkpoint_exists(&count_ok) {
            eprintln!("-- Trinity (rust) | Jellyfish count (K={k})");
            // 硬护栏: target_fa 整读进内存 + 计数哈希。
            if let Some(msg) =
                memory_guard_error(args.max_memory, file_len(target_fa), "k-mer counting")
            {
                return Err(CommonError::Parse(msg));
            }
            let data = fs::read(target_fa)?;
            // estimate_hash_size（Trinity:2598-2604 容量提示语义）→ 计数
            // HashMap 预 reserve（每线程 estimate/18/CPU，函数内部再截上限）,
            // 避免增长期 rehash 内存峰值。
            let estimate = estimate_hash_size(args.max_memory, data.len() as u64);
            let per_thread = (estimate / KMER_ENTRY_BYTES / args.cpu.max(1) as u64) as usize;
            let counts = KmerCountTable::count_fasta_data_with_capacity(&data, k, ds, per_thread);
            save_counts(&counts_bin, &counts, k, ds)?;
            fs::write(&count_ok, b"")?;
        }
        Ok(())
    })?;

    // ---- inchworm（L2654-2727）----
    // 命名镜像 L1576-1580: SS → inchworm.fa, 无 SS → inchworm.DS.fa。
    let iworm_name = if args.ss_lib_type.is_some() {
        "inchworm.fa"
    } else {
        "inchworm.DS.fa"
    };
    let iworm_fa = round_dir.join(iworm_name);
    let iworm_tmp = round_dir.join(format!("{iworm_name}.tmp"));
    let iworm_ok = round_dir.join(format!(".iworm.{k}.{SR_FLAG}.ok"));
    trinity_common::timer::timed(round_dir, "inchworm", || {
        run_with_checkpoint(&iworm_ok, "Inchworm (linear contig construction)", || {
            // M7 直桥: 计数表（checkpoint 二进制或重算）直接构建 counter,
            // 跳过 jellyfish dump 文本往返; 原 dump 的 -L min_cov 过滤由
            // 下方 prune_some_kmers 等价承担。
            let (mut counter, _parsed) = trinity_common::timer::timed(
                round_dir,
                "kmer_load",
                || -> Result<(trinity_inchworm::kmer_counter::KmerCounter, usize), CommonError> {
                    let mut counts = if checkpoint_exists(&counts_bin) {
                        load_counts(&counts_bin, k, ds)?
                    } else {
                        let data = fs::read(target_fa)?;
                        let estimate = estimate_hash_size(args.max_memory, data.len() as u64);
                        let per_thread =
                            (estimate / KMER_ENTRY_BYTES / args.cpu.max(1) as u64) as usize;
                        KmerCountTable::count_fasta_data_with_capacity(&data, k, ds, per_thread)
                    };
                    // multi-k trusted contig 注入（presence-only）: 计数表两条
                    // 构建路径汇合后的统一注入点, 先于 populate/剪枝。
                    if let Some(t) = trusted {
                        let tdata = fs::read(t)?;
                        let inserted =
                            trinity_kmer::counter::inject_trusted_kmers(&mut counts, &tdata, k, ds);
                        eprintln!(
                            "-multi-k: trusted contigs {} → 新注入 {inserted} 个 k-mer（count=1, presence-only）",
                            t.display()
                        );
                    }
                    populate_from_table(&counts, k, ds)
                },
            )?;
            write_kmer_count_report(&round_dir.join("inchworm.kmer_count"), counter.size())?;
            // 主线实参: --min_any_entropy 1.0 --no_prune_error_kmers（L1050/L2707）
            // + dump -L 已按 min_kmer_cov 过滤（minKmerCount 此处为等价二次过滤）。
            eprintln!(
                "Pruning kmers (min_kmer_count={} min_any_entropy=1 min_ratio_non_error=0.005)",
                args.min_kmer_count
            );
            prune_some_kmers(&mut counter, args.min_kmer_count, 1.0, false, 0.005);
            let params = IrkeParams::default();
            // -L = MIN_IWORM_LEN = KMER_SIZE（L1064/L2698）
            let aparams = AssemblyParams {
                min_assembly_length: k,
                ..Default::default()
            };
            let monitor = Monitor::new(1);
            let mut buf = Vec::new();
            compute_sequence_assemblies_parallel(
                counter,
                &params,
                &aparams,
                &monitor,
                true, // TWO_PHASE 默认开（IRKE.cpp:59）
                Some(args.inchworm_cpu),
                &mut buf,
            )?;
            fs::write(&iworm_tmp, &buf)?;
            Ok(())
        })
    })?;
    let iworm_renamed_ok = round_dir.join(format!("iworm_renamed.{k}.{SR_FLAG}.ok"));
    run_with_checkpoint(&iworm_renamed_ok, "mv inchworm.fa.tmp", || {
        if checkpoint_exists(&iworm_tmp) {
            fs::rename(&iworm_tmp, &iworm_fa)?;
        } else if !checkpoint_exists(&iworm_fa) {
            return Err(CommonError::Parse(format!(
                "inchworm output missing: {}",
                iworm_fa.display()
            )));
        }
        Ok(())
    })?;
    // L1835-1837: touch <inchworm>.finished
    let iworm_finished = round_dir.join(format!("{iworm_name}.finished"));
    if !checkpoint_exists(&iworm_finished) {
        fs::write(&iworm_finished, b"")?;
    }
    Ok(iworm_fa)
}

/// 跨设备安全移动（--scratch 与 --output 常在不同文件系统, rename 直接失败）。
pub(crate) fn move_file(src: &Path, dst: &Path) -> Result<(), CommonError> {
    match fs::rename(src, dst) {
        Ok(()) => Ok(()),
        Err(_) => {
            fs::copy(src, dst)?;
            fs::remove_file(src)?;
            Ok(())
        }
    }
}

/// 递归模式主线（`--recursive`, M1）: 取代 classic 的 chrysalis-quantify +
/// butterfly 池 + harvest 三段。
///
/// - phase-1: chrysalis 止于分区（`skip_quantify: true`, pct 50 钳制,
///   与原版 TRINITY_COMPLETE 语义一致）, checkpoint 用 `.recursive_phase1.ok`（防跨模式污染, M1 终审 Issue 2）;
/// - phase-2: rayon 分区池逐分区 `run_mini` → `read_partitions/c<id>.mini.fa`
///   （存在即跳过 = 断点）, 整段 checkpoint `.recursive.ok`;
/// - 汇总: mini.fa 拼接 → min_len 过滤 + 全序列去重 → `Trinity.tmp.fasta`
///   → 最终命名/gene_trans_map 与 classic 尾段一致（无 checkpoint, 每次重汇总）。
fn run_recursive(
    a: TrinityArgs,
    args: &TrinityArgs,
    outdir: &Path,
    workdir: &Path,
    chrysalis_dir: &Path,
    iworm_data: &[u8],
    reads_data: &[u8],
) -> Result<PathBuf, CommonError> {
    let ss = a.ss_lib_type.as_deref();
    warn_memory(
        args.max_memory,
        iworm_data.len() as u64 + reads_data.len() as u64,
        "Chrysalis (recursive phase-1)",
    );
    // ---- phase-1: chrysalis 止于分区（.tmp 保留, 供分区池消费）----
    let mut listing_opt: Option<Vec<(u64, PathBuf)>> = None;
    let quantify_ok = workdir.join(".recursive_phase1.ok");
    run_with_checkpoint(
        &quantify_ok,
        "Chrysalis (recursive phase-1: clustering & partitioning)",
        || {
            let listing = run_chrysalis_pipeline(
                iworm_data,
                reads_data,
                chrysalis_dir,
                &ChrysalisParams {
                    skip_quantify: true,
                    strand: ss.is_some(),
                    min_contig_length: args.min_contig_length,
                    pct_required: 50,
                    max_reads: args.max_reads_per_graph,
                    threads: args.cpu.max(1),
                    graphs_per_partition: 1000,
                    graph_clean_rnaspades: false, // 递归 phase-1 跳过 quantify, 清理不生效
                },
            )?;
            listing_opt = Some(listing);
            Ok(())
        },
    )?;
    let listing = match listing_opt {
        Some(l) => l,
        None => read_component_base_listing(
            &chrysalis_dir
                .join("Component_bins")
                .join("component_base_listing.txt"),
        )?,
    };
    if listing.is_empty() {
        return Err(CommonError::Parse(
            "WARNING, component base listing file is empty - likely sparse data".to_string(),
        ));
    }

    // ---- phase-2: 分区池（checkpoint .recursive.ok; 命中时直接从 mini.fa 重汇总）----
    warn_memory(
        args.max_memory,
        reads_data.len() as u64,
        "Recursive mini pool",
    );
    let mini_dir = workdir.join("read_partitions");
    run_with_checkpoint(
        &workdir.join(".recursive.ok"),
        "Recursive mini partition pool",
        || run_partition_pool(&listing, chrysalis_dir, args, &mini_dir).map(|_| ()),
    )?;
    let mini_files: Vec<PathBuf> = listing
        .iter()
        .map(|(id, _)| mini_dir.join(format!("c{id}.mini.fa")))
        .collect();

    // ---- 汇总: min_len 过滤 + 全局去重（print_butterfly_assemblies 语义）----
    let trinity_tmp = workdir.join("Trinity.tmp.fasta");
    let n_transcripts = harvest_mini(&mini_files, args.min_contig_length, &trinity_tmp)?;
    if n_transcripts == 0 {
        return Err(CommonError::Parse(
            "WARNING: no transcripts harvested from recursive partition output".to_string(),
        ));
    }

    // 最终命名与 gene_trans_map: 与 classic 尾段一致（scratch 跨设备走 move_file）。
    let final_fa = PathBuf::from(format!("{}.Trinity.fasta", outdir.display()));
    move_file(&trinity_tmp, &final_fa)?;
    write_gene_trans_map(
        &final_fa,
        &final_fa.with_file_name(format!(
            "{}.gene_trans_map",
            final_fa.file_name().unwrap_or_default().to_string_lossy()
        )),
    )?;

    // 收尾清理与 classic 相同集合（mini.fa 与 chrysalis 产物保留——断点依赖）。
    if !args.no_cleanup {
        let target_fa = workdir.join("both.fa");
        for p in [
            &target_fa,
            &workdir.join(format!("mer_counts.{}.{}.rs.bin", args.kmer_size, SR_FLAG)),
            &workdir.join(format!("jellyfish.kmers.{}.{}.fa", args.kmer_size, SR_FLAG)),
            &workdir.join("left.fa"),
            &workdir.join("right.fa"),
        ] {
            let _ = fs::remove_file(p);
        }
    }
    Ok(final_fa)
}

/// get_Trinity_gene_to_trans_map.pl: header 正则 `^(.*c\d+_g\d+)(_i\d+)` →
/// `gene\ttrans`; 不匹配的 header 警告跳过。
pub(crate) fn write_gene_trans_map(fasta: &Path, map_path: &Path) -> Result<(), CommonError> {
    let text = fs::read_to_string(fasta)?;
    let mut out = String::new();
    for line in text.lines() {
        if !line.starts_with('>') {
            continue;
        }
        let acc = line[1..].split_whitespace().next().unwrap_or("");
        match split_gene_trans(acc) {
            Some((gene, trans)) => {
                out.push_str(&gene);
                out.push('\t');
                out.push_str(&trans);
                out.push('\n');
            }
            None => eprintln!("Error, could not parse transcript name: {acc}"),
        }
    }
    fs::write(map_path, out)?;
    Ok(())
}

/// 正则 `^(.*c\d+_g\d+)(_i\d+)`（贪婪 gene 前缀; trans = 整体匹配）的手工展开。
fn split_gene_trans(acc: &str) -> Option<(String, String)> {
    let b = acc.as_bytes();
    // 边界 i 从大到小（`.*` 贪婪）: acc[i..] = "_i<digits>"（≥1 位, 尽量长）,
    // acc[..i] 以 `c<digits>_g<digits>` 结尾。
    for i in (1..b.len()).rev() {
        if b[i] != b'_' || i + 2 >= b.len() || b[i + 1] != b'i' || !b[i + 2].is_ascii_digit() {
            continue;
        }
        let head = &acc[..i];
        let hb = head.as_bytes();
        let mut j = hb.len();
        while j > 0 && hb[j - 1].is_ascii_digit() {
            j -= 1;
        }
        if j == hb.len() || j == 0 || hb[j - 1] != b'g' || hb[j - 2] != b'_' {
            continue;
        }
        let m = j - 2; // `_g<digits>` 的 `_`
        if m == 0 {
            continue;
        }
        let mut p = m;
        while p > 0 && hb[p - 1].is_ascii_digit() {
            p -= 1;
        }
        if p == m + 1 || p == 0 || hb[p - 1] != b'c' {
            continue;
        }
        let mut e = i + 2;
        while e < b.len() && b[e].is_ascii_digit() {
            e += 1;
        }
        return Some((head.to_string(), acc[..e].to_string()));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn guard_too_small_reports_estimate_basis() {
        // 1G 内存 vs 1G 输入 → 1.2G 需求 → 报错且含估算依据与建议
        let g = 1024 * 1024 * 1024u64;
        let msg = memory_guard_error(g, g, "k-mer counting").expect("should error");
        assert!(msg.contains("1.2x input"), "{msg}");
        assert!(msg.contains("distinct k-mers x 18 B"), "{msg}");
        assert!(msg.contains("increase --max_memory"), "{msg}");
        assert!(msg.contains("k-mer counting"), "{msg}");
    }

    #[test]
    fn guard_exactly_enough_passes() {
        let g = 6_000_000_000u64;
        // 输入 = 5G → need = 5G + 1G = 6G == max_memory → 恰好通过（>= 判定）。
        let input = 5_000_000_000u64;
        assert_eq!(memory_guard_error(g, input, "x"), None);
        assert!(memory_guard_error(g, input + 1, "x").is_some());
    }

    #[test]
    fn warn_memory_only_prints() {
        // 软护栏无返回值; 极端入参不 panic 即可。
        warn_memory(1, u64::MAX / 2, "Chrysalis");
        warn_memory(u64::MAX, 0, "Butterfly");
    }

    #[test]
    fn default_scratch_dir_deterministic_per_outdir() {
        let a = default_scratch_dir(Path::new("/x/out"));
        assert_eq!(a, default_scratch_dir(Path::new("/x/out")));
        assert_ne!(a, default_scratch_dir(Path::new("/x/out2")));
        assert!(a
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .starts_with("trinity-rust-"));
    }

    #[test]
    fn default_scratch_check_falls_back_when_unwritable() {
        // 不可创建/不可写 → Err（编排层回落 --output）
        assert!(check_default_scratch(Path::new("/proc/nonexistent/tr-rs-x"), 1).is_err());
        // 需求超过可用空间 → Err（刚建的空目录被清掉）
        assert!(check_default_scratch(
            &std::env::temp_dir().join(format!("tr_rs_huge_{}", std::process::id())),
            u64::MAX / 2
        )
        .is_err());
        // 可写且空间充足 → Ok
        let d = std::env::temp_dir().join(format!("tr_rs_ok_{}", std::process::id()));
        assert!(check_default_scratch(&d, 1).is_ok());
        let _ = fs::remove_dir(&d);
    }
}

/// filter_iworm_by_min_length_or_cov.pl 镜像（Trinity:2073 调用, 参数 100 10）:
/// 保留 `len(seq) >= min_len || cov >= min_cov` 的记录原样字节（header 行 + 序列行
/// 照抄——下游 GraphFromFasta 重读时格式细节无关语义）。cov 解析失败按 0（Perl 数值上下文）。
pub fn filter_iworm_min_len_or_cov(
    data: &[u8],
    min_len: usize,
    min_cov: u32,
) -> Result<Vec<u8>, CommonError> {
    let mut out = Vec::with_capacity(data.len());
    let mut i = 0usize;
    while i < data.len() {
        // header 行
        let hdr_end = data[i..]
            .iter()
            .position(|&b| b == b'\n')
            .map(|p| i + p)
            .unwrap_or(data.len());
        let header = &data[i..hdr_end];
        // 记录体到下一个 '>' 行首或文件尾
        let mut j = hdr_end + 1;
        let body_end = loop {
            if j >= data.len() {
                break data.len();
            }
            if data[j] == b'>' {
                break j;
            }
            let nl = data[j..]
                .iter()
                .position(|&b| b == b'\n')
                .map(|p| j + p + 1)
                .unwrap_or(data.len());
            j = nl;
        };
        let body = &data[hdr_end + 1..body_end];
        // accession = header 首个空白分隔 token（去 '>'）; cov = split(';') 第二段的数字前缀
        let acc = header
            .strip_prefix(b">")
            .unwrap_or(header)
            .split(|&b| b == b' ' || b == b'\t')
            .find(|t| !t.is_empty())
            .unwrap_or(&[]);
        let cov: u32 = acc
            .split(|&b| b == b';')
            .nth(1)
            .map(|seg| {
                let digits: Vec<u8> = seg
                    .iter()
                    .copied()
                    .skip_while(|&b| b == b' ' || b == b'\t')
                    .take_while(|b| b.is_ascii_digit())
                    .collect();
                String::from_utf8_lossy(&digits).parse().unwrap_or(0)
            })
            .unwrap_or(0);
        let seq_len = body.iter().filter(|&&b| b != b'\n' && b != b'\r').count();
        if seq_len >= min_len || cov >= min_cov {
            out.extend_from_slice(header);
            out.push(b'\n');
            out.extend_from_slice(body);
        }
        i = body_end;
    }
    Ok(out)
}

#[cfg(test)]
mod min100_tests {
    use super::filter_iworm_min_len_or_cov;

    #[test]
    fn filter_keeps_len_or_cov() {
        // 注意: 过滤按实际序列长度（length($sequence)）, 非 header 的 length 字段
        let seq300 = "ACGT".repeat(75);
        let fa = format!(">a1;5 total_counts: 10 K: 25 length: 300\n{seq300}\n>a2;50 total_counts: 99 K: 25 length: 50\nAAAA\n>a3;1 total_counts: 2 K: 25 length: 40\nTTTT\n").into_bytes();
        let out = filter_iworm_min_len_or_cov(&fa, 100, 10).unwrap();
        let s = String::from_utf8(out).unwrap();
        // a1: len 300>=100 保留; a2: cov 50>=10 保留; a3: 都不满足 剔除
        assert!(s.contains(">a1;5"));
        assert!(s.contains(">a2;50"));
        assert!(!s.contains(">a3"));
    }

    #[test]
    fn filter_boundary_exact() {
        // 恰好 100bp 保留; 恰好 cov=10 保留
        let seq100 = "A".repeat(100);
        let fa = format!(">x;0 length: 100\n{seq100}\n>y;10 length: 5\nACGTA\n");
        let out = filter_iworm_min_len_or_cov(fa.as_bytes(), 100, 10).unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(s.contains(">x;") && s.contains(">y;"));
    }
}
