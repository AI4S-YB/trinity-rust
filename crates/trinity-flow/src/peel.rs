//! 最重路径剥离（DAG 动态规划 + 瓶颈减流）。

use crate::graph::FlowGraph;
use std::collections::HashMap;

pub struct PeelParams {
    pub min_support: u64,
    pub min_path_nodes: usize,
}

/// 分量内剥离。返回路径列表（node id 序列, 源→汇）。
/// 注：graph 自身不被修改（在副本的 capacity 上操作）。
///
/// 性能: 邻接表 + 边容量用紧凑索引数组（一次构建, 剥离循环内零哈希查找——
/// 逐迭代 O(V+E) 而非原 HashMap 边表全扫描的 O(V*E)）。
/// 语义: 瓶颈 < min_support 的路径仍消耗其流量（流分解完整收敛）但不产出
/// （否则单条低权边会阻断整条分量, 见 c0 fixture: 起点首边权重 1）。
pub fn peel_component(g: &FlowGraph, comp: &[i64], p: &PeelParams) -> Vec<Vec<i64>> {
    let order = match g.toposort_component(comp) {
        Some(o) => o,
        None => return handle_cycle(g, comp, p), // 环 → 剪最小边重试
    };
    let set: std::collections::HashSet<i64> = comp.iter().copied().collect();

    // 索引空间: 节点 → usize, 边 → usize。边按 (u, v) 排序保证确定性。
    let idx: HashMap<i64, usize> = comp.iter().enumerate().map(|(i, &id)| (id, i)).collect();
    let n = comp.len();
    let mut edge_list: Vec<(i64, i64)> = g
        .edges
        .keys()
        .filter(|&&(u, v)| set.contains(&u) && set.contains(&v))
        .copied()
        .collect();
    edge_list.sort_unstable();
    let mut cap: Vec<u64> = edge_list
        .iter()
        .map(|&(u, v)| g.edges.get(&(u, v)).copied().unwrap_or(0))
        .collect();
    let mut edge_pos: HashMap<(usize, usize), usize> = HashMap::with_capacity(edge_list.len());
    let mut out_edges: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut in_edges: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (i, &(u, v)) in edge_list.iter().enumerate() {
        let (ui, vi) = (idx[&u], idx[&v]);
        edge_pos.insert((ui, vi), i);
        out_edges[ui].push(i);
        in_edges[vi].push(i);
    }

    let mut paths: Vec<Vec<i64>> = Vec::new();
    loop {
        // DAG DP: best[v] = (最大总权, 前驱边)
        let mut best: Vec<(u64, Option<usize>)> = vec![(0, None); n];
        for &v in &order {
            let vi = idx[&v];
            let mut bv = (0u64, None);
            for &e in &in_edges[vi] {
                let c = cap[e];
                if c == 0 {
                    continue;
                }
                let pred = edge_list[e].0;
                let tot = best[idx[&pred]].0.saturating_add(c);
                if tot > bv.0 {
                    bv = (tot, Some(e));
                }
            }
            best[vi] = bv;
        }
        // 汇 = 分量中无正容量后继且 best > 0 的节点中 best 最大者
        let mut sink = None;
        let mut sink_w = 0u64;
        for &v in &order {
            let vi = idx[&v];
            let has_succ = out_edges[vi].iter().any(|&e| cap[e] > 0);
            if !has_succ && best[vi].0 > sink_w {
                sink_w = best[vi].0;
                sink = Some(vi);
            }
        }
        let sink = match sink {
            Some(s) if sink_w > 0 => s,
            _ => break,
        };
        // 回溯路径（边 → 节点, 索引空间）
        let mut path = vec![sink];
        let mut cur = sink;
        while let Some(e) = best[cur].1 {
            cur = idx[&edge_list[e].0];
            path.push(cur);
        }
        path.reverse();
        // 瓶颈 = 路径相邻对的最小剩余容量
        let mut bottleneck = u64::MAX;
        for w in path.windows(2) {
            let c = cap[edge_pos[&(w[0], w[1])]];
            if c < bottleneck {
                bottleneck = c;
            }
        }
        if bottleneck == u64::MAX || path.len() < p.min_path_nodes {
            break;
        }
        // 减流（< min_support 时仍消耗但不产出——低权边不再阻断分量）
        for w in path.windows(2) {
            let e = edge_pos[&(w[0], w[1])];
            cap[e] = cap[e].saturating_sub(bottleneck);
        }
        if bottleneck >= p.min_support {
            paths.push(path.iter().map(|&i| comp[i]).collect());
        }
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
