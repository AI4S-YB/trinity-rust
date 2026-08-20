//! chrysalis-all 全链管线（库层）——原版 Trinity 管线里 shell/perl 串联
//! Chrysalis 六子阶段的编排上移（`chrysalis-all` 子命令与库层同入口）。
//!
//! 输入是**内存字节**（iworm fasta 全文 + reads fasta 全文），输出布局镜像
//! 原版 Chrysalis 目录（bundled_iworm_contigs.fasta / Component_bins/Cbin_*/）。
//!
//! 链路：GraphFromFasta → 焊接图 sort → BubbleUpClustering →
//! CreateIwormFastaBundle → ReadsToTranscripts → sort → FastaToDeBruijn
//! (-K 24) → partition → QuantifyGraph（rayon 线程池，`-t` 语义保持）。

use std::path::{Path, PathBuf};

use rayon::prelude::*;
use trinity_common::error::CommonError;
use trinity_common::fasta::FastaReader;
use trinity_inchworm::debruijn::graph_per_record;

use crate::bubble_up::{bubble_up_clustering, BubbleParams};
use crate::bundle::create_iworm_fasta_bundle;
use crate::dna_vector::{read_fasta_bytes, stream_fasta_records};
use crate::graph_from_fasta::{graph_from_fasta, sort_weld_graph, GffParams};
use crate::partition::{partition, PartParams};
use crate::quantify::{quantify_graph, reads_ext_filename, QgParams};
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
    let iworm_seqs = read_fasta_bytes(iworm_fa, false);
    let reads_text = String::from_utf8_lossy(reads_fa).into_owned();
    let reads_streamed = stream_fasta_records(&reads_text);
    let welds = graph_from_fasta(
        &iworm_seqs,
        &reads_streamed,
        &GffParams {
            strand: p.strand,
            threads: p.threads.max(1),
            ..Default::default()
        },
    )?;
    let welds_sorted = sort_weld_graph(&welds);

    // 2. BubbleUpClustering → component.out → bundle
    println!("-running BubbleUpClustering...");
    let comps = bubble_up_clustering(
        &iworm_seqs,
        &welds_sorted,
        &BubbleParams {
            min_contig_length: p.min_contig_length,
            ..Default::default()
        },
    )?;
    println!("-running CreateIwormFastaBundle...");
    let bundle = create_iworm_fasta_bundle(&comps, 200)?;
    let bundle_path = out_dir.join("bundled_iworm_contigs.fasta");
    std::fs::write(&bundle_path, &bundle)?;

    // 3. ReadsToTranscripts → sort
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
            ..Default::default()
        },
    )?;
    let rtt_sorted = sort_reads_to_components(&rtt.text);

    // 4. FastaToDeBruijn（Trinity 实参 -K 24）
    println!("-running FastaToDeBruijn...");
    let f = std::fs::File::open(&bundle_path)?;
    let mut reader = FastaReader::new(std::io::BufReader::new(f));
    let mut recs = Vec::new();
    while let Some(r) = reader.next_record()? {
        recs.push(r);
    }
    let debruijn = graph_per_record(&recs, 24, p.strand)?;

    // 5. partition → Component_bins/Cbin_*/c*.{graph.tmp,reads.tmp}
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

    // 6. QuantifyGraph 逐组件（rayon 并行，组件独立）
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
    let failures: Vec<String> = {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(p.threads.max(1))
            .build()
            .map_err(|e| CommonError::Parse(format!("cannot build thread pool: {e}")))?;
        pool.install(|| {
            listing
                .par_iter()
                .filter_map(|(_, base)| {
                    let g = base.with_extension("graph.tmp");
                    let r = base.with_extension("reads.tmp");
                    let out = base.with_extension("graph.out");
                    match run_quantify_for_all(&g, &r, &out, &qg) {
                        Ok(()) => None,
                        Err(e) => Some(format!("{}: {e}", base.display())),
                    }
                })
                .collect()
        })
    };
    if !failures.is_empty() {
        return Err(CommonError::Parse(format!(
            "QuantifyGraph failed for {} components: {}",
            failures.len(),
            failures.join("; ")
        )));
    }
    println!(
        "-chrysalis-all complete: {} components in {}",
        listing.len(),
        bins.display()
    );
    Ok(listing)
}

fn run_quantify_for_all(g: &Path, r: &Path, out: &Path, p: &QgParams) -> Result<(), String> {
    let o =
        quantify_graph(&g.to_string_lossy(), &r.to_string_lossy(), p).map_err(|e| e.to_string())?;
    let want = out.to_string_lossy().to_string();
    let derived = PathBuf::from(&o.graph_out);
    if derived != out {
        std::fs::rename(&derived, out).map_err(|e| e.to_string())?;
        std::fs::rename(o.reads_out, reads_ext_filename(&want)).map_err(|e| e.to_string())?;
    }
    // 原版语义：成功后删输入
    let _ = std::fs::remove_file(g);
    let _ = std::fs::remove_file(r);
    Ok(())
}
