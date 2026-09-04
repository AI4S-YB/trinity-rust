//! chrysalis-all 全链管线（库层）——原版 Trinity 管线里 shell/perl 串联
//! Chrysalis 六子阶段的编排上移（`chrysalis-all` 子命令与库层同入口）。
//!
//! 输入是**内存字节**（iworm fasta 全文 + reads fasta 全文），输出布局镜像
//! 原版 Chrysalis 目录（bundled_iworm_contigs.fasta / Component_bins/Cbin_*/）。
//!
//! 链路：GraphFromFasta → 焊接图 sort → BubbleUpClustering →
//! CreateIwormFastaBundle → ReadsToTranscripts → sort → FastaToDeBruijn
//! (-K 24) → partition → QuantifyGraph（rayon 线程池，`-t` 语义保持）。

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use rayon::prelude::*;
use trinity_common::error::CommonError;
use trinity_common::fasta::FastaReader;
use trinity_common::timer::timed;
use trinity_inchworm::debruijn::graph_per_record;

use crate::bubble_up::{bubble_up_clustering, BubbleParams};
use crate::bundle::create_iworm_fasta_bundle;
use crate::dna_vector::{read_fasta_bytes, stream_fasta_records};
use crate::graph_from_fasta::{graph_from_fasta, sort_weld_graph, GffParams};
use crate::pack::{
    PackIndex, PackWriter, GRAPHS_PACK, KIND_GRAPH_OUT, KIND_GRAPH_READS, KIND_GRAPH_TMP,
    KIND_READS_TMP, TMP_PACK,
};
use crate::partition::{partition, PartParams};
use crate::quantify::{quantify_graph_mem, QgParams};
use crate::reads_to_transcripts::{reads_to_transcripts, sort_reads_to_components, RttParams};

/// chrysalis-all 管线参数（CLI `-o` 之外的全部旗标面）。
#[derive(Debug, Clone)]
pub struct ChrysalisParams {
    /// `--SS`（strand-specific）。
    pub strand: bool,
    /// `-L`：进入组件的最短 contig 长度（bubble-up 与 partition 共用）。
    pub min_contig_length: usize,
    /// `-p`：ReadsToTranscripts 的 pct_required。
    pub pct_required: u32,
    /// `-N`：每个分区（Cbin_*）的图数上限。
    pub graphs_per_partition: usize,
    /// `--max_reads`：QuantifyGraph 每组件读入上限。
    pub max_reads: u64,
    /// `-t`：QuantifyGraph 阶段 rayon 线程池大小。
    pub threads: usize,
    /// 跳过 QuantifyGraph（phase-2 分割点）：止于分区，保留 c*.graph.tmp /
    /// c*.reads.tmp 供递归模式逐分区消费（默认 false，行为不变）。
    pub skip_quantify: bool,
    /// `--graph-clean-rnaspades`：QuantifyGraph 之后对 graph.out 做 rnaSPAdes
    /// 式相对覆盖清理（tip/EC/rcc, 见 graph_clean 模块）。默认 false = 逐字节不变。
    pub graph_clean_rnaspades: bool,
    /// `--rtt-multi-kmer`（Rust 版扩展, 默认 false = 原版语义）：ReadsToTranscripts
    /// 的 25-mer 索引多值化——共享 kmer 记全部组件。救跨 weld 口 read 被
    /// 其他组件抢占导致的断口无桥（见 RttParams::multi_kmer 文档）。
    pub rtt_multi_kmer: bool,
    /// GraphFromFasta 的 weld 证据门槛 `glue_factor`（默认 0.05 = 原版）:
    /// min_cov = clamp(max(cov_i,cov_c) * glue_factor, min_glue, max_glue)。
    /// 实验: 聚类割裂档（missing tx 物理相邻组件被拆）的放宽方向。
    pub weld_glue_factor: f64,
    /// weld 最小 read 支持数 `min_glue_required`（默认 2 = 原版）。
    pub weld_min_glue: u32,
    /// BubbleUpClustering 的 `-max_cluster_size`（默认 25 = 原版）。
    pub cluster_max_size: usize,
    /// **Rust 版扩展**: weld（GraphFromFasta）阶段专用 reads 文件路径——
    /// None = 用管线 reads_fa（原版语义）。分级归一化实验: weld 用
    /// max_cov 更高的归一化版本（低/中表达区保留更多桥 read 证据）,
    /// RTT/quantify 仍用默认归一化（控图规模）。
    pub weld_reads_path: Option<PathBuf>,
    /// `--gap-close-spades`：清理之后对 graph.out 做 SPAdes 式 tip 对 tip
    /// overlap 缝合（PE 证据排序候选, 见 gap_close 模块）。默认 false = 逐字节不变。
    pub gap_close_spades: bool,
}

impl Default for ChrysalisParams {
    fn default() -> Self {
        ChrysalisParams {
            strand: false,
            min_contig_length: 200,
            pct_required: 50,
            graphs_per_partition: 1000,
            max_reads: 200_000,
            threads: 1,
            skip_quantify: false,
            graph_clean_rnaspades: false,
            rtt_multi_kmer: false,
            weld_glue_factor: 0.05,
            weld_min_glue: 2,
            cluster_max_size: 25,
            weld_reads_path: None,
            gap_close_spades: false,
        }
    }
}

/// chrysalis-all 全链：iworm fasta 字节 + reads fasta 字节 → out_dir 下的
/// Chrysalis 目录布局。返回 partition listing（(comp id, 路径前缀) 升序）。
///
/// QuantifyGraph 逐组件在 `threads` 大小的 rayon 池上并行（原版
/// `perl -P`/xargs 语义）；成功后按 QuantifyGraph.cc:489-493 删 .tmp 输入。
pub fn run_chrysalis_pipeline(
    iworm_fa: &[u8],
    reads_fa: &[u8],
    out_dir: &Path,
    p: &ChrysalisParams,
) -> Result<Vec<(u64, PathBuf)>, CommonError> {
    std::fs::create_dir_all(out_dir)?;

    // 1. GraphFromFasta → weld 图 → sort
    println!("-running GraphFromFasta...");
    let (iworm_seqs, welds_sorted) = timed(out_dir, "gff", || {
        let reads_text: String = match &p.weld_reads_path {
            Some(path) => {
                let bytes = std::fs::read(path).map_err(|e| {
                    CommonError::Parse(format!("weld reads {}: {e}", path.display()))
                })?;
                String::from_utf8_lossy(&bytes).into_owned()
            }
            None => String::from_utf8_lossy(reads_fa).into_owned(),
        };
        let reads_streamed = stream_fasta_records(&reads_text);
        let iworm_seqs = read_fasta_bytes(iworm_fa, false);
        let welds = graph_from_fasta(
            &iworm_seqs,
            &reads_streamed,
            &GffParams {
                strand: p.strand,
                threads: p.threads.max(1),
                glue_factor: p.weld_glue_factor,
                min_glue_required: p.weld_min_glue,
                ..Default::default()
            },
        )?;
        let welds_sorted = sort_weld_graph(&welds);
        Ok::<_, CommonError>((iworm_seqs, welds_sorted))
    })?;

    // 2. BubbleUpClustering → component.out → bundle
    let comps = timed(out_dir, "bubble", || {
        println!("-running BubbleUpClustering...");
        bubble_up_clustering(
            &iworm_seqs,
            &welds_sorted,
            &BubbleParams {
                min_contig_length: p.min_contig_length,
                max_cluster_size: p.cluster_max_size,
                ..Default::default()
            },
        )
    })?;
    let (bundle, bundle_path) = timed(out_dir, "bundle", || {
        println!("-running CreateIwormFastaBundle...");
        // 2026-09-03: 原版此处硬编码 200; 改为随 --min_contig_length 联动——
        // QM4 端部截断解剖显示孤立 100-200bp 组件在此被整体丢弃（TPM≥20 截断
        // 缺口的 42% 归因于该档, 见 opt-log「端部截断全景」）。
        let bundle = create_iworm_fasta_bundle(&comps, p.min_contig_length)?;
        let bundle_path = out_dir.join("bundled_iworm_contigs.fasta");
        std::fs::write(&bundle_path, &bundle)?;
        Ok::<_, CommonError>((bundle, bundle_path))
    })?;

    // 3. ReadsToTranscripts → sort
    let rtt_sorted = timed(out_dir, "rtt", || {
        println!("-running ReadsToTranscripts...");
        let bundles = read_fasta_bytes(bundle.as_bytes(), false);
        let reads_seqs = read_fasta_bytes(reads_fa, true);
        let rtt = reads_to_transcripts(
            &reads_seqs,
            &bundles,
            &RttParams {
                strand: p.strand,
                pct_required: p.pct_required,
                threads: p.threads.max(1),
                multi_kmer: p.rtt_multi_kmer,
                ..Default::default()
            },
        )?;
        Ok::<_, CommonError>(sort_reads_to_components(&rtt.text))
    })?;

    // 4. FastaToDeBruijn（Trinity 实参 -K 24）
    let debruijn = timed(out_dir, "f2db", || {
        println!("-running FastaToDeBruijn...");
        let f = std::fs::File::open(&bundle_path)?;
        let mut reader = FastaReader::new(std::io::BufReader::new(f));
        let mut recs = Vec::new();
        while let Some(r) = reader.next_record()? {
            recs.push(r);
        }
        graph_per_record(&recs, 24, p.strand)
    })?;

    // 5. partition → Component_bins/Cbin_*/c*.{graph.tmp,reads.tmp}
    let (bins, listing) = timed(out_dir, "partition", || {
        println!("-running partition...");
        let bins = out_dir.join("Component_bins");
        let listing = partition(
            &debruijn,
            &rtt_sorted,
            &bins,
            &PartParams {
                graphs_per_partition: p.graphs_per_partition,
                min_contig_length: p.min_contig_length,
            },
        )?;
        Ok::<_, CommonError>((bins, listing))
    })?;

    // 6. QuantifyGraph 逐组件（rayon 并行，组件独立）
    //    skip_quantify: 递归模式分割点——止于分区, 保留 tmp 包供递归分区池消费。
    //    输入从每 Cbin 的 `tmp.pack` 按偏移读; 产物写每 Cbin 一个
    //    `graphs.pack`+`graphs.idx`（包在阶段开始 truncate 重建 = 阶段级
    //    checkpoint 重跑的幂等语义）。整阶段成功后删 tmp 包（镜像原版
    //    逐组件删 .tmp 输入的空间回收语义, 粗一档: 整包删, 失败保留供排查）。
    if !p.skip_quantify {
        println!(
            "-running QuantifyGraph on {} components ({} threads)...",
            listing.len(),
            p.threads
        );
        let qg = QgParams {
            strand: p.strand,
            max_reads: p.max_reads as i64,
            ..Default::default()
        };
        // rnaSPAdes 图清理的跨组件统计与耗时累计（各组件并行, 锁合并）
        let clean_stats: Mutex<crate::graph_clean::CleanStats> = Mutex::new(Default::default());
        let clean_nanos = std::sync::atomic::AtomicU64::new(0);
        // SPAdes 式 gap 缝合的跨组件统计与耗时累计（同清理的合并模式）
        let gap_stats: Mutex<crate::gap_close::GapCloseStats> = Mutex::new(Default::default());
        let gap_nanos = std::sync::atomic::AtomicU64::new(0);
        let mut writers: HashMap<PathBuf, Mutex<PackWriter>> = HashMap::new();
        let mut tmp_indexes: HashMap<PathBuf, PackIndex> = HashMap::new();
        for (_, base) in &listing {
            let dir = base.parent().map(|d| d.to_path_buf()).unwrap_or_default();
            if !writers.contains_key(&dir) {
                writers.insert(
                    dir.clone(),
                    Mutex::new(PackWriter::create(&dir, GRAPHS_PACK)?),
                );
                tmp_indexes.insert(
                    dir.clone(),
                    PackIndex::load(&dir.join(format!("{TMP_PACK}.idx")))?,
                );
            }
        }
        let failures: Vec<String> = timed(
            out_dir,
            "quantify",
            || -> Result<Vec<String>, CommonError> {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(p.threads.max(1))
                    .build()
                    .map_err(|e| CommonError::Parse(format!("cannot build thread pool: {e}")))?;
                Ok(pool.install(|| {
                    listing
                        .par_iter()
                        .filter_map(|(id, base)| {
                            let dir = base.parent().map(|d| d.to_path_buf()).unwrap_or_default();
                            let w = writers.get(&dir).expect("writer per Cbin pre-created");
                            let idx = tmp_indexes.get(&dir).expect("tmp index per Cbin loaded");
                            match run_quantify_packed(
                                *id,
                                &dir,
                                idx,
                                w,
                                &qg,
                                p.graph_clean_rnaspades,
                                &clean_stats,
                                &clean_nanos,
                                p.gap_close_spades,
                                &gap_stats,
                                &gap_nanos,
                            ) {
                                Ok(()) => None,
                                Err(e) => Some(format!("{}: {e}", base.display())),
                            }
                        })
                        .collect()
                }))
            },
        )?;
        for w in writers.values() {
            w.lock().unwrap().finish()?;
        }
        // 清理统计行 + 落 .stage_timings（CPU 时间累计: 各组件并行清理的耗时和）
        if p.graph_clean_rnaspades {
            let st = clean_stats.lock().unwrap();
            println!(
                "-graph-clean-rnaspades: {} components, tip chains {} (g1 {} + g2 {}, {} nodes), ec chains {} ({} nodes), rcc components {} ({} nodes), lines removed {}, dangling repaired {}, reverted {}",
                st.components,
                st.tip_g1_chains + st.tip_g2_chains,
                st.tip_g1_chains,
                st.tip_g2_chains,
                st.tip_nodes,
                st.ec_chains,
                st.ec_nodes,
                st.rcc_components,
                st.rcc_nodes,
                st.lines_removed,
                st.repaired_dangling,
                st.reverted_components,
            );
            drop(st);
            let secs = clean_nanos.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9;
            if let Ok(mut fh) = std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(out_dir.join(".stage_timings"))
            {
                use std::io::Write;
                let _ = writeln!(fh, "graph_clean\t{secs:.3}");
            }
        }
        // 缝合统计行 + 落 .stage_timings（CPU 时间累计: 各组件并行缝合的耗时和）
        if p.gap_close_spades {
            let st = gap_stats.lock().unwrap();
            println!(
                "-gap-close-spades: {} components, tips (out {} + in {}), tip reads {}, pair votes {}, candidates {}, checked {}, closed {} (perfect {} + corrected {}), rejected (low-complexity {} + hamming {} + kmer-clash {}), bridge nodes {}",
                st.components,
                st.out_tips,
                st.in_tips,
                st.tip_reads,
                st.pair_votes,
                st.candidate_pairs,
                st.pairs_checked,
                st.gaps_closed,
                st.closed_perfect,
                st.closed_corrected,
                st.rejected_low_complexity,
                st.rejected_hamming,
                st.rejected_kmer_clash,
                st.bridge_nodes,
            );
            drop(st);
            let secs = gap_nanos.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9;
            if let Ok(mut fh) = std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(out_dir.join(".stage_timings"))
            {
                use std::io::Write;
                let _ = writeln!(fh, "gap_close\t{secs:.3}");
            }
        }
        if !failures.is_empty() {
            return Err(CommonError::Parse(format!(
                "QuantifyGraph failed for {} components: {}",
                failures.len(),
                failures.join("; ")
            )));
        }
        // 整阶段成功 → 删 tmp 包（重跑由 partition 重建）
        for dir in tmp_indexes.keys() {
            let _ = std::fs::remove_file(dir.join(format!("{TMP_PACK}.pack")));
            let _ = std::fs::remove_file(dir.join(format!("{TMP_PACK}.idx")));
        }
    }
    println!(
        "-chrysalis-all complete: {} components in {}",
        listing.len(),
        bins.display()
    );
    Ok(listing)
}

// tmp.pack 的 per-thread 打开句柄缓存（同一线程连续命中同 Cbin 时省
// open/close; NFS 上每次 open 都是一次元数据往返）。
thread_local! {
    static TMP_FILES: std::cell::RefCell<HashMap<PathBuf, std::fs::File>> =
        std::cell::RefCell::new(HashMap::new());
}

fn read_tmp_record(dir: &Path, idx: &PackIndex, id: u64, kind: u64) -> Result<Vec<u8>, String> {
    let (offset, len) = idx.find(id, kind).ok_or_else(|| {
        format!(
            "tmp pack {} has no record (id={id}, kind={kind})",
            dir.display()
        )
    })?;
    let path = dir.join(format!("{TMP_PACK}.pack"));
    TMP_FILES.with(|m| {
        let mut m = m.borrow_mut();
        if !m.contains_key(&path) {
            let f = std::fs::File::open(&path).map_err(|e| e.to_string())?;
            m.insert(path.clone(), f);
        }
        crate::pack::read_at(&m[&path], offset, len).map_err(|e| e.to_string())
    })
}

/// 逐组件 QuantifyGraph: 从 Cbin 的 tmp 包读 graph.tmp/reads.tmp 记录,
/// 产物追加进 graphs 包。输入的删除由调用方在整阶段成功后整包完成。
/// `graph_clean`：开启时对 graph.out 文本做 rnaSPAdes 式清理（tip→EC→rcc,
/// 见 graph_clean 模块）; `gap_close`：清理之后做 SPAdes 式 tip 对 tip
/// overlap 缝合（PE 证据来自 reads.tmp, 见 gap_close 模块）。统计/耗时
/// 分别并入对应 Mutex/原子量。
#[allow(clippy::too_many_arguments)]
fn run_quantify_packed(
    id: u64,
    dir: &Path,
    tmp_idx: &PackIndex,
    writer: &Mutex<PackWriter>,
    p: &QgParams,
    graph_clean: bool,
    clean_stats: &Mutex<crate::graph_clean::CleanStats>,
    clean_nanos: &std::sync::atomic::AtomicU64,
    gap_close: bool,
    gap_stats: &Mutex<crate::gap_close::GapCloseStats>,
    gap_nanos: &std::sync::atomic::AtomicU64,
) -> Result<(), String> {
    let graph = read_tmp_record(dir, tmp_idx, id, KIND_GRAPH_TMP)?;
    let reads = read_tmp_record(dir, tmp_idx, id, KIND_READS_TMP)?;
    let graph = String::from_utf8(graph).map_err(|e| e.to_string())?;
    let (gout, rout) = quantify_graph_mem(&graph, &reads, p).map_err(|e| e.to_string())?;
    let gout = if graph_clean {
        let t0 = std::time::Instant::now();
        let (cleaned, st) = crate::graph_clean::clean_graph_rnaspades(&gout);
        clean_nanos
            .fetch_add(t0.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
        clean_stats.lock().unwrap().add_assign(&st);
        cleaned
    } else {
        gout
    };
    // SPAdes 流程里 gap closer 在 simplification 之后跑（config.info
    // after_simplify=true）——此处清理之后缝合, 顺序一致。
    let gout = if gap_close {
        let t0 = std::time::Instant::now();
        let (closed, st) = crate::gap_close::close_gaps_spades(&gout, &reads, p.strand);
        gap_nanos
            .fetch_add(t0.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
        gap_stats.lock().unwrap().add_assign(&st);
        closed
    } else {
        gout
    };
    {
        let mut w = writer.lock().unwrap();
        w.append(id, KIND_GRAPH_OUT, gout.as_bytes())
            .map_err(|e| e.to_string())?;
        w.append(id, KIND_GRAPH_READS, rout.as_bytes())
            .map_err(|e| e.to_string())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn skip_quantify_stops_after_partition() {
        // 用 fixtures/p3 的小输入（gff.iworm.fa + gff.reads.fa）
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../fixtures/p3/");
        let iworm = std::fs::read(format!("{root}gff.iworm.fa")).unwrap();
        let reads = std::fs::read(format!("{root}gff.reads.fa")).unwrap();
        let dir = std::env::temp_dir().join("chry_skip_quantify_test");
        let _ = std::fs::remove_dir_all(&dir);
        let listing = run_chrysalis_pipeline(
            &iworm,
            &reads,
            &dir,
            &ChrysalisParams {
                skip_quantify: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(!listing.is_empty());
        // tmp 包在且含 c0 的 graph.tmp/reads.tmp 记录（递归模式要读）;
        // graphs 包不在（quantify 未跑）; tmp 包未被清理
        let cbin0 = dir.join("Component_bins/Cbin0");
        let g = crate::pack::read_entry(
            &cbin0,
            crate::pack::TMP_PACK,
            0,
            crate::pack::KIND_GRAPH_TMP,
        )
        .expect("c0 graph.tmp 记录应在 tmp 包内");
        assert!(!g.is_empty());
        let r = crate::pack::read_entry(
            &cbin0,
            crate::pack::TMP_PACK,
            0,
            crate::pack::KIND_READS_TMP,
        )
        .expect("c0 reads.tmp 记录应在 tmp 包内");
        assert!(!r.is_empty());
        assert!(!cbin0
            .join(format!("{}.pack", crate::pack::GRAPHS_PACK))
            .exists());
        // per-component 小文件不再产生
        let c0 = dir.join("Component_bins/Cbin0/c0");
        assert!(!c0.with_extension("graph.tmp").exists());
        assert!(!c0.with_extension("graph.out").exists());
    }
}
