//! 最重路径剥离（DAG 动态规划 + 瓶颈减流）。

use crate::graph::FlowGraph;
use std::collections::HashMap;

pub struct PeelParams {
    pub min_support: u64,
    pub min_path_nodes: usize,
}

/// 分量内剥离。返回路径列表（node id 序列, 源→汇）。
/// 注：graph 自身不被修改（在副本的 capacity 上操作）。
pub fn peel_component(g: &FlowGraph, comp: &[i64], p: &PeelParams) -> Vec<Vec<i64>> {
    let order = match g.toposort_component(comp) {
        Some(o) => o,
        None => return handle_cycle(g, comp, p), // 环 → 剪最小边重试
    };
    let set: std::collections::HashSet<i64> = comp.iter().copied().collect();
    // 剩余容量
    let mut cap: HashMap<(i64, i64), u64> = g
        .edges
        .iter()
        .filter(|((u, v), _)| set.contains(u) && set.contains(v))
        .map(|(&k, &w)| (k, w))
        .collect();

    let mut paths = Vec::new();
    loop {
        // DAG DP: best[v] = (最大总权, 前驱)
        let mut best: HashMap<i64, (u64, Option<i64>)> = HashMap::new();
        for &v in &order {
            let mut bv = (0u64, None);
            for pred in g.predecessors(v) {
                if let Some(&c) = cap.get(&(pred, v)) {
                    if c == 0 {
                        continue;
                    }
                    if let Some(&(bp, _)) = best.get(&pred) {
                        let tot = bp.saturating_add(c);
                        if tot > bv.0 {
                            bv = (tot, Some(pred));
                        }
                    }
                }
            }
            best.insert(v, bv);
        }
        // 汇 = 分量中无后继（或后继容量全 0）且 best > 0 的节点中 best 最大者
        let mut sink = None;
        let mut sink_w = 0u64;
        for &v in &order {
            let succs: Vec<i64> = g
                .successors(v)
                .into_iter()
                .filter(|s| cap.get(&(v, *s)).copied().unwrap_or(0) > 0)
                .collect();
            if succs.is_empty() {
                if let Some(&(bv, _)) = best.get(&v) {
                    if bv > sink_w {
                        sink_w = bv;
                        sink = Some(v);
                    }
                }
            }
        }
        let sink = match sink {
            Some(s) if sink_w > 0 => s,
            _ => break,
        };
        // 回溯路径
        let mut path = vec![sink];
        let mut cur = sink;
        while let Some(Some(pred)) = best.get(&cur).map(|(_, p)| *p) {
            path.push(pred);
            cur = pred;
        }
        path.reverse();
        // 瓶颈 = 路径相邻对的最小剩余容量
        let mut bottleneck = u64::MAX;
        for w in path.windows(2) {
            let c = cap.get(&(w[0], w[1])).copied().unwrap_or(0);
            if c < bottleneck {
                bottleneck = c;
            }
        }
        if bottleneck == u64::MAX || bottleneck < p.min_support || path.len() < p.min_path_nodes {
            break;
        }
        // 减流
        for w in path.windows(2) {
            if let Some(c) = cap.get_mut(&(w[0], w[1])) {
                *c = c.saturating_sub(bottleneck);
            }
        }
        paths.push(path);
    }
    paths
}

/// 有环分量的保守处理: 全分量最小权边剪断 → 递归 peel
fn handle_cycle(g: &FlowGraph, comp: &[i64], p: &PeelParams) -> Vec<Vec<i64>> {
    let set: std::collections::HashSet<i64> = comp.iter().copied().collect();
    let min_edge = g
        .edges
        .iter()
        .filter(|((u, v), _)| set.contains(u) && set.contains(v))
        .min_by_key(|(_, &w)| w);
    match min_edge {
        Some((&(u, v), _)) => {
            let mut g2 = FlowGraph {
                nodes: g.nodes.clone(),
                edges: g.edges.clone(),
                component_id: g.component_id,
            };
            g2.edges.remove(&(u, v));
            peel_component(&g2, comp, p)
        }
        None => Vec::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::FlowGraph;

    fn mini() -> FlowGraph {
        // 菱形: 1→{2,3}→4, 权重 1→2=8, 1→3=3, 2→4=8, 3→4=3
        // 剥离应先走 1→2→4 (瓶颈 8), 再走 1→3→4 (瓶颈 3)
        let mut g = FlowGraph::default();
        for (id, kmer) in [
            (1i64, b"AAAA".to_vec()),
            (2, b"CCCC".to_vec()),
            (3, b"GGGG".to_vec()),
            (4, b"TTTT".to_vec()),
        ] {
            g.nodes.insert(id, crate::graph::FlowNode { kmer });
        }
        g.edges.insert((1, 2), 8);
        g.edges.insert((1, 3), 3);
        g.edges.insert((2, 4), 8);
        g.edges.insert((3, 4), 3);
        g
    }

    #[test]
    fn peel_diamond() {
        let g = mini();
        let params = PeelParams {
            min_support: 1,
            min_path_nodes: 2,
        };
        let paths = peel_component(&g, &[1, 2, 3, 4], &params);
        // 第一条走 1→2→4（总权 16 > 6）; 第二条走 1→3→4（剩余 3）
        assert_eq!(paths.len(), 2);
        assert_eq!(paths[0], vec![1, 2, 4]);
        assert_eq!(paths[1], vec![1, 3, 4]);
    }

    #[test]
    fn peel_respects_min_support() {
        let g = mini();
        let params = PeelParams {
            min_support: 5,
            min_path_nodes: 2,
        };
        let paths = peel_component(&g, &[1, 2, 3, 4], &params);
        // 1→3→4 瓶颈 3 < 5 → 不产出
        assert_eq!(paths.len(), 1);
        assert_eq!(paths[0], vec![1, 2, 4]);
    }

    #[test]
    fn peel_stops_when_exhausted() {
        let g = mini();
        let params = PeelParams {
            min_support: 1,
            min_path_nodes: 2,
        };
        // 同参数再剥——全消费完不无限循环
        let p1 = peel_component(&g, &[1, 2, 3, 4], &params);
        let p2 = peel_component(&g, &[1, 2, 3, 4], &params);
        assert_eq!(p1.len(), 2);
        assert_eq!(p2.len(), 2); // 确定性
    }
}
