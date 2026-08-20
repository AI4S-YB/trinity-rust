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
use std::io::Write;
use std::path::{Path, PathBuf};

use trinity_chrysalis::pipeline::{run_chrysalis_pipeline, ChrysalisParams};
use trinity_common::error::CommonError;
use trinity_inchworm::irke::{
    compute_sequence_assemblies_parallel, populate_from_kmers, prune_some_kmers,
    write_kmer_count_report, AssemblyParams, IrkeParams, Monitor,
};
use trinity_kmer::counter::{estimate_hash_size, KmerCountTable};
use trinity_kmer::diginorm::{self, DigiNormParams, ReadsInput};
use trinity_kmer::dump::write_dump;

use crate::args::{SeqType, TrinityArgs};
use crate::butterfly_pool::{run_butterfly_pool, ButterflyPoolParams};
use crate::checkpoint::{checkpoint_exists, run_with_checkpoint};
use crate::harvest::harvest;
use crate::prep::prep_reads;

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
        ..Default::default()
    }
}

/// mer_counts 二进制记录布局: `u64 key, u32 count`（小端）, 头部魔数
/// `TRMC` + u32 k + u8 ds（读回时校验）。
fn save_counts(
    path: &Path,
    counts: &trinity_kmer::counter::CountMap,
    k: usize,
    ds: bool,
) -> Result<(), CommonError> {
    let mut buf: Vec<u8> = Vec::new();
    buf.extend_from_slice(b"TRMC");
    buf.extend_from_slice(&(k as u32).to_le_bytes());
    buf.push(ds as u8);
    for (&key, &c) in counts.iter() {
        buf.extend_from_slice(&key.to_le_bytes());
        buf.extend_from_slice(&c.to_le_bytes());
    }
    fs::write(path, buf)?;
    Ok(())
}

fn load_counts(
    path: &Path,
    k: usize,
    ds: bool,
) -> Result<trinity_kmer::counter::CountMap, CommonError> {
    let data = fs::read(path)?;
    if data.len() < 9 || &data[..4] != b"TRMC" {
        return Err(CommonError::Parse(format!(
            "corrupt k-mer table {}: bad header",
            path.display()
        )));
    }
    let k_in = u32::from_le_bytes([data[4], data[5], data[6], data[7]]) as usize;
    let ds_in = data[8] != 0;
    if k_in != k || ds_in != ds {
        return Err(CommonError::Parse(format!(
            "k-mer table {} was built with different k/DS settings",
            path.display()
        )));
    }
    let mut counts = trinity_kmer::counter::CountMap::default();
    let mut i = 9;
    while i + 12 <= data.len() {
        let key = u64::from_le_bytes(data[i..i + 8].try_into().unwrap());
        let c = u32::from_le_bytes(data[i + 8..i + 12].try_into().unwrap());
        counts.insert(key, c);
        i += 12;
    }
    Ok(counts)
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
    let chrysalis_dir = outdir.join("chrysalis");
    fs::create_dir_all(&chrysalis_dir)?;

    let paired = !args.left.is_empty();
    let ss = args.ss_lib_type.as_deref();
    let k = args.kmer_size;
    let ds = ss.is_none();
    let mut a = args.clone();

    // rayon 全局池按 --CPU 设一次（局部需不同线程数的入口——inchworm/
    // chrysalis GFF/RTT/quantify/butterfly 池——各自显式 scoped pool 覆盖;
    // `-t 1`/`--CPU 1` 时全局默认即单线程）。
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(args.cpu.max(1))
        .build_global();

    // ---- 1. 归一化（L1441-1453; checkpoint insilico_read_normalization/normalization.ok）----
    if !args.no_normalize_reads {
        // 硬护栏: 归一化需把两侧读入整读进内存。
        if let Some(msg) = memory_guard_error(
            args.max_memory,
            input_files_bytes(args),
            "read normalization",
        ) {
            return Err(CommonError::Parse(msg));
        }
        let norm_dir = outdir.join("insilico_read_normalization");
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
        let ext = ext_of(args.seq_type);
        if paired {
            a.left = vec![norm_dir.join(format!("left.norm.{ext}"))];
            a.right = vec![norm_dir.join(format!("right.norm.{ext}"))];
        } else {
            a.single = vec![norm_dir.join(format!("single.norm.{ext}"))];
        }
    }

    // ---- 2. prep → both.fa/single.fa + .read_count（L1604-1705）----
    let (target_fa, read_count) = prep_reads(&a, &outdir)?;
    let read_count_file = outdir.join(format!(
        "{}.read_count",
        target_fa.file_name().unwrap_or_default().to_string_lossy()
    ));
    if fs::metadata(&read_count_file).map(|m| m.len()).unwrap_or(0) == 0 {
        fs::write(&read_count_file, format!("{read_count}\n"))?;
    }
    let _ = read_count;

    // ---- 3. jellyfish 计数 + dump（L2584-2632）----
    let count_ok = outdir.join(format!(
        ".jellyfish_count.{k}.mincov{}.{SR_FLAG}.ok",
        args.min_kmer_count
    ));
    let dump_ok = outdir.join(format!(
        ".jellyfish_dump.{k}.mincov{}.{SR_FLAG}.ok",
        args.min_kmer_count
    ));
    let kmer_fa = outdir.join(format!("jellyfish.kmers.{k}.{SR_FLAG}.fa"));
    let counts_bin = outdir.join(format!("mer_counts.{k}.{SR_FLAG}.rs.bin"));
    if !checkpoint_exists(&count_ok) {
        eprintln!("-- Trinity (rust) | Jellyfish count (K={k})");
        // 硬护栏: target_fa 整读进内存 + 计数哈希。
        if let Some(msg) =
            memory_guard_error(args.max_memory, file_len(&target_fa), "k-mer counting")
        {
            return Err(CommonError::Parse(msg));
        }
        let data = fs::read(&target_fa)?;
        // estimate_hash_size（Trinity:2598-2604 容量提示语义）→ 计数
        // HashMap 预 reserve（每线程 estimate/18/CPU，函数内部再截上限）,
        // 避免增长期 rehash 内存峰值。
        let estimate = estimate_hash_size(args.max_memory, data.len() as u64);
        let per_thread = (estimate / KMER_ENTRY_BYTES / args.cpu.max(1) as u64) as usize;
        let counts = KmerCountTable::count_fasta_data_with_capacity(&data, k, ds, per_thread);
        save_counts(&counts_bin, &counts, k, ds)?;
        fs::write(&count_ok, b"")?;
    }
    if !checkpoint_exists(&dump_ok) {
        eprintln!(
            "-- Trinity (rust) | Jellyfish dump -L {}",
            args.min_kmer_count
        );
        let counts = if checkpoint_exists(&counts_bin) {
            load_counts(&counts_bin, k, ds)?
        } else {
            let data = fs::read(&target_fa)?;
            let estimate = estimate_hash_size(args.max_memory, data.len() as u64);
            let per_thread = (estimate / KMER_ENTRY_BYTES / args.cpu.max(1) as u64) as usize;
            KmerCountTable::count_fasta_data_with_capacity(&data, k, ds, per_thread)
        };
        let mut out = fs::File::create(&kmer_fa)?;
        write_dump(&mut out, &counts, k, args.min_kmer_count.max(1), ds)?;
        out.flush()?;
        fs::write(&dump_ok, b"")?;
    }

    // ---- 4. inchworm（L2654-2727）----
    // 命名镜像 L1576-1580: SS → inchworm.fa, 无 SS → inchworm.DS.fa。
    let iworm_name = if ss.is_some() {
        "inchworm.fa"
    } else {
        "inchworm.DS.fa"
    };
    let iworm_fa = outdir.join(iworm_name);
    let iworm_tmp = outdir.join(format!("{iworm_name}.tmp"));
    let iworm_ok = outdir.join(format!(".iworm.{k}.{SR_FLAG}.ok"));
    run_with_checkpoint(&iworm_ok, "Inchworm (linear contig construction)", || {
        let data = fs::read(&kmer_fa)?;
        let (mut counter, _parsed) = populate_from_kmers(&data, k, ds)?;
        write_kmer_count_report(&outdir.join("inchworm.kmer_count"), counter.size())?;
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
    })?;
    let iworm_renamed_ok = outdir.join(format!("iworm_renamed.{k}.{SR_FLAG}.ok"));
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
    let _ = iworm_fa;
    // L1835-1837: touch <inchworm>.finished
    let iworm_finished = outdir.join(format!("{iworm_name}.finished"));
    if !checkpoint_exists(&iworm_finished) {
        fs::write(&iworm_finished, b"")?;
    }
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
    let mut listing_opt: Option<Vec<(u64, PathBuf)>> = None;
    let quantify_ok = outdir.join(".quantify_graph.ok");
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
                    ..Default::default()
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

    // ---- 6. butterfly 组件池（checkpoint .butterfly.ok; L1934-1936）----
    warn_memory(args.max_memory, file_len(&iworm_fa), "Butterfly");
    run_with_checkpoint(
        &outdir.join(".butterfly.ok"),
        "Butterfly (component pool)",
        || {
            run_butterfly_pool(
                &listing,
                &outdir,
                &ButterflyPoolParams {
                    cpu: args.cpu.max(1),
                    min_contig_length: args.min_contig_length,
                    group_pairs_distance: args.group_pairs_distance,
                    stack_size_mb: args.bfly_stack_mb,
                },
            )
        },
    )?;

    // ---- 7. 汇总（无 checkpoint; 每次重写 Trinity.tmp.fasta, L1916/L1943）----
    let trinity_tmp = outdir.join("Trinity.tmp.fasta");
    let n_transcripts = harvest(&listing, args.min_contig_length, &trinity_tmp)?;
    if n_transcripts == 0 {
        return Err(CommonError::Parse(
            "WARNING: no transcripts harvested from butterfly output".to_string(),
        ));
    }

    // 最终命名: `<outdir 绝对路径>.Trinity.fasta`（L1512）+ gene_trans_map。
    let final_fa = PathBuf::from(format!("{}.Trinity.fasta", outdir.display()));
    fs::rename(&trinity_tmp, &final_fa)?;
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
    if !args.no_cleanup {
        for p in [
            &target_fa,
            &counts_bin,
            &kmer_fa,
            &outdir.join("left.fa"),
            &outdir.join("right.fa"),
        ] {
            let _ = fs::remove_file(p);
        }
    }
    Ok(final_fa)
}

/// get_Trinity_gene_to_trans_map.pl: header 正则 `^(.*c\d+_g\d+)(_i\d+)` →
/// `gene\ttrans`; 不匹配的 header 警告跳过。
fn write_gene_trans_map(fasta: &Path, map_path: &Path) -> Result<(), CommonError> {
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
