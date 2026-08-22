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
        None => return handle_cycle(g, comp, p), // 环 → SCC 缩点成 DAG 再剥
    };
    peel_dag(g, comp, &order, p, &std::collections::HashSet::new())
}

/// DAG 分量剥离主体。`silent` 为缩点代表节点（多节点 SCC）——路径穿过时
/// 仍消耗流量但不产出（首版保守: 不产出 SCC 内部路径, 见 handle_cycle）。
fn peel_dag(
    g: &FlowGraph,
    comp: &[i64],
    order: &[i64],
    p: &PeelParams,
    silent: &std::collections::HashSet<i64>,
) -> Vec<Vec<i64>> {
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
        for &v in order {
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
        for &v in order {
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
        if bottleneck == u64::MAX {
            break;
        }
        // 减流（短于 min_path_nodes / < min_support 时仍消耗但不产出——
        // 否则短路径耗尽后剩余流被整体丢弃; bottleneck > 0 保证循环收敛）
        for w in path.windows(2) {
            let e = edge_pos[&(w[0], w[1])];
            cap[e] = cap[e].saturating_sub(bottleneck);
        }
        if bottleneck >= p.min_support
            && path.len() >= p.min_path_nodes
            && !path.iter().any(|&i| silent.contains(&comp[i]))
        {
            paths.push(path.iter().map(|&i| comp[i]).collect());
        }
    }
    paths
}

/// 有环分量的处理: Tarjan SCC 缩点 → DAG → 常规 peel。
///
/// - 每个 SCC 缩为单节点, 代表取 SCC 内最小节点 id（确定性）;
/// - 桥接权重: 缩点后 A→B 的边权 = 原 A→B 所有边权**之和**（总流量守恒,
///   多条平行桥接边的流不会在缩点时丢失/重复）;
/// - SCC 内部边（含自环）丢弃——缩点图必为 DAG;
/// - 首版保守: **不产出 SCC 内部路径**——穿过 >1 节点 SCC 的路径只消耗
///   流量不产出（SCC 内部回溯路径的序列续接无意义, 留待后续版本展开）,
///   因此只产 DAG 部分的路径。
///
/// 取代原「剪最小权边 → 递归重剥」: dev1 实数据 129/6387 分量含环
/// （最大 5398 环节点）, 逐边剪断的递归克隆为 O(E²) 内存且栈溢出;
/// Tarjan 为显式栈实现, 缩点后一次 DAG peel, 全程无深递归。
fn handle_cycle(g: &FlowGraph, comp: &[i64], p: &PeelParams) -> Vec<Vec<i64>> {
    let sccs = g.scc_component(comp);
    let mut rep_of: HashMap<i64, i64> = HashMap::with_capacity(comp.len());
    let mut condensed = FlowGraph {
        component_id: g.component_id,
        ..Default::default()
    };
    // silent = 多节点 SCC 的代表节点: 穿过它的路径不产出
    let mut silent = std::collections::HashSet::new();
    for scc in &sccs {
        let rep = scc[0]; // scc 内已排序 → 最小 id
        if scc.len() > 1 {
            silent.insert(rep);
        }
        rep_of.extend(scc.iter().map(|&id| (id, rep)));
        let kmer = g
            .nodes
            .get(&rep)
            .map(|n| n.kmer.clone())
            .unwrap_or_default();
        condensed.nodes.insert(rep, crate::graph::FlowNode { kmer });
    }
    let set: std::collections::HashSet<i64> = comp.iter().copied().collect();
    for (&(u, v), &w) in &g.edges {
        if !set.contains(&u) || !set.contains(&v) {
            continue;
        }
        let (ru, rv) = (rep_of[&u], rep_of[&v]);
        if ru != rv {
            // 桥接权重求和: 流量守恒
            let e = condensed.edges.entry((ru, rv)).or_insert(0);
            *e = e.saturating_add(w);
        }
    }
    let mut comp2: Vec<i64> = condensed.nodes.keys().copied().collect();
    comp2.sort_unstable();
    let order = match condensed.toposort_component(comp2.as_slice()) {
        Some(o) => o,
        // 缩点图必为 DAG; 到达此处说明逻辑错误, 保守返回空而非死循环
        None => return Vec::new(),
    };
    peel_dag(&condensed, &comp2, &order, p, &silent)
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

    /// 3 节点环 {1,2,3} + DAG 尾巴 3→4→5 与独立入度 6→4。
    /// SCC 缩点后: 只产出 DAG 尾巴的路径（6→4→5, 6→4）,
    /// 不产出穿过环的路径, 不死循环。
    #[test]
    fn peel_cycle_component_via_scc() {
        let mut g = FlowGraph::default();
        for id in 1..=6i64 {
            g.nodes
                .insert(id, crate::graph::FlowNode { kmer: Vec::new() });
        }
        // 环: 1→2→3→1（权重 5）
        g.edges.insert((1, 2), 5);
        g.edges.insert((2, 3), 5);
        g.edges.insert((3, 1), 5);
        // 尾巴: 3→4→5（权重 4）, 独立入口 6→4（权重 9, 保证 DP 先走 DAG 路）
        g.edges.insert((3, 4), 4);
        g.edges.insert((4, 5), 4);
        g.edges.insert((6, 4), 9);
        let params = PeelParams {
            min_support: 1,
            min_path_nodes: 2,
        };
        let paths = peel_component(&g, &[1, 2, 3, 4, 5, 6], &params);
        // 不死循环（测试本身会超时）+ 只含 DAG 节点
        assert!(paths
            .iter()
            .all(|p| !p.contains(&1) && !p.contains(&2) && !p.contains(&3)));
        // DAG 尾巴的两条路径都产出
        assert!(paths.contains(&vec![6, 4, 5]));
        assert!(paths.contains(&vec![6, 4]));
    }
}
