//! 顶层流分解引擎: graph.out + graph.reads → FASTA。

use trinity_common::error::CommonError;

use crate::graph::FlowGraph;
use crate::paths::{dedup_similar, path_to_seq, to_fasta, PathOutput};
use crate::peel::{peel_component, PeelParams};

pub struct FlowParams {
    /// 最短输出（bp）。
    pub min_len: usize,
    /// 路径瓶颈下限（边剩余容量 ≥ 此值才产出）。
    pub min_support: u64,
    /// 去冗余同一率。
    pub dedup_threshold: f64,
    /// kmer 长（原始 24-mer 图 = 24; 相邻节点重叠 K-1, 每步追加 1 碱基）。
    pub k: usize,
}

impl Default for FlowParams {
    fn default() -> Self {
        FlowParams {
            min_len: 200,
            min_support: 2,
            dedup_threshold: 0.95,
            k: 24,
        }
    }
}

/// 单组件流分解。`graph_out`/`graph_reads` 为 Chrysalis 量化产物内容;
/// 返回 allProbPaths 等价格式 FASTA（header `{comp}_g1_i{n}`）。
///
/// M2 首版不消费 reads（纯流; PE 约束为 M2.1 迭代项）, 参数保留以稳定接口。
pub fn run_flow_component(
    graph_out: &str,
    _graph_reads: &str,
    p: &FlowParams,
) -> Result<String, CommonError> {
    let g = FlowGraph::parse(graph_out.as_bytes())?;
    let kmers: Vec<(i64, Vec<u8>)> = g
        .nodes
        .iter()
        .filter(|(_, n)| !n.kmer.is_empty())
        .map(|(&id, n)| (id, n.kmer.clone()))
        .collect();

    let mut all_outputs: Vec<PathOutput> = Vec::new();
    for comp in g.weak_components() {
        // 跳过太小的分量（无法构成路径）。
        if comp.len() < 2 {
            continue;
        }
        let paths = peel_component(
            &g,
            &comp,
            &PeelParams {
                min_support: p.min_support,
                min_path_nodes: 2,
            },
        );
        for path in paths {
            let seq = path_to_seq(&path, &kmers, p.k);
            if seq.len() < p.min_len {
                continue;
            }
            let support = path
                .windows(2)
                .map(|w| g.edge_weight(w[0], w[1]).unwrap_or(0))
                .min()
                .unwrap_or(0);
            all_outputs.push(PathOutput { seq, support, path });
        }
    }
    let deduped = dedup_similar(all_outputs, p.dedup_threshold);
    Ok(to_fasta(&deduped, &format!("c{}", g.component_id), 0, 60))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flow_on_c0_fixture() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../fixtures/p3/quantify/c0/"
        );
        let g = std::fs::read_to_string(format!("{root}orig.graph.out")).unwrap();
        let r = std::fs::read_to_string(format!("{root}orig.reads.out")).unwrap();
        let out = run_flow_component(&g, &r, &FlowParams::default()).unwrap();
        // 产出非空 FASTA
        assert!(out.contains('>') && out.len() > 200, "c0 输出过短");
        // 转录本数 1-20（c0 的 classic butterfly 产出 5 条, 流分解应同量级）
        let n = out.matches('>').count();
        assert!((1..=20).contains(&n), "c0 流分解产出 {n} 条");
        // 全部 ≥ min_len
        for line in out.lines() {
            if line.starts_with('>') {
                let len: usize = line
                    .split("len=")
                    .nth(1)
                    .unwrap()
                    .split(' ')
                    .next()
                    .unwrap()
                    .parse()
                    .unwrap();
                assert!(len >= FlowParams::default().min_len, "len={len}");
            }
        }
    }

    #[test]
    fn empty_input_yields_empty_output() {
        let out = run_flow_component("Component 7\n", "", &FlowParams::default()).unwrap();
        assert!(out.is_empty());
    }
}
