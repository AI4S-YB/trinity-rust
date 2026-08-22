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
            // support 上报口径: 路径上**原始**边权的最小值（非剥离时的剩余
            // 容量/瓶颈）——同一物理转录本的多条残余路径各自报满 support,
            // 由 dedup 取高者; 保守高估而非低估。
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

    /// dev1 Cbin0/c11 实数据截取的 ~100 节点有环子图
    /// （fixtures/p3/quantify/c11cyc/, 源分量 6155 节点含环——曾触发
    /// 递归剪环的栈溢出/O(E²) 内存）。回归: SCC 缩点后不死锁、
    /// 有环分量被识别且剥离终止。
    #[test]
    fn flow_on_c11_cyclic_fixture() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../fixtures/p3/quantify/c11cyc/"
        );
        let g_text = std::fs::read_to_string(format!("{root}orig.graph.out")).unwrap();
        let g = crate::graph::FlowGraph::parse(g_text.as_bytes()).unwrap();
        let mut n_cyclic = 0;
        for comp in g.weak_components() {
            if g.toposort_component(&comp).is_none() {
                n_cyclic += 1;
            }
        }
        assert!(n_cyclic > 0, "fixture 应含环");
        // min_len=1: 该子图的 DAG 尾巴路径短（24bp/节点）, 默认 200bp 会滤空;
        // 回归点是「不死锁不崩溃 + DAG 部分产出」, 非全长过滤。
        let out = run_flow_component(
            &g_text,
            "",
            &FlowParams {
                min_len: 1,
                ..FlowParams::default()
            },
        )
        .unwrap();
        assert!(!out.is_empty(), "c11 有环子图应产出非空 FASTA");
    }

    /// dev1 全量 Cbin0/c11（6155 节点, 含 6155 节点环分量——曾触发栈溢出）。
    /// 需 chrysalis 缓存: `TRINITY_C11_FULL=<path>/Cbin0/c11.graph.out cargo test -- --ignored`
    #[test]
    #[ignore]
    fn flow_on_c11_full_component() {
        let path = std::env::var("TRINITY_C11_FULL").expect("set TRINITY_C11_FULL");
        let g_text = std::fs::read_to_string(&path).unwrap();
        let out = run_flow_component(
            &g_text,
            "",
            &FlowParams {
                min_len: 1,
                ..FlowParams::default()
            },
        )
        .unwrap();
        assert!(!out.is_empty());
    }
}
