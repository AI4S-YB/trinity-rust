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
                ..Default::default()
            },
        )
    })?;
    let (bundle, bundle_path) = timed(out_dir, "bundle", || {
        println!("-running CreateIwormFastaBundle...");
        let bundle = create_iworm_fasta_bundle(&comps, 200)?;
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
                            match run_quantify_packed(*id, &dir, idx, w, &qg) {
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
fn run_quantify_packed(
    id: u64,
    dir: &Path,
    tmp_idx: &PackIndex,
    writer: &Mutex<PackWriter>,
    p: &QgParams,
) -> Result<(), String> {
    let graph = read_tmp_record(dir, tmp_idx, id, KIND_GRAPH_TMP)?;
    let reads = read_tmp_record(dir, tmp_idx, id, KIND_READS_TMP)?;
    let graph = String::from_utf8(graph).map_err(|e| e.to_string())?;
    let (gout, rout) = quantify_graph_mem(&graph, &reads, p).map_err(|e| e.to_string())?;
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
