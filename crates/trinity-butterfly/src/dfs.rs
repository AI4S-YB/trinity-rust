//! 镜像 `My_DFS.java`：`runDFS2()`（双向 visitVertex2 + down-only 两阶段）与
//! finish time 降序拓扑序。
//!
//! 与 Java 的差异点：
//! - 顶点迭代序：JUNG `getVertices()` 是 HashSet 任意序，我们用插入序（确定性）。
//!   discovery/finish time 因此可能与 Java 不同（同一棵 DFS 森林的形状一致时才相同）；
//!   depth 结果对迭代序不敏感（down-only 阶段的平局用 (depth, 插入序) 稳定排序）。
//! - 递归实现（Java 同为递归）。c0 类组件小；大组件（数万节点链）存在栈深风险，
//!   T10 端到端若溢出再转显式栈（转显式栈时须保持 finish time 的"全部邻居
//!   递归完成后"语义 = 出栈序）。

use rustc_hash::FxHashMap;

use crate::graph::DiGraph;

const WHITE: u8 = 0;
const BLACK: u8 = 1;
const GRAY: u8 = 2;

/// `My_DFS.runDFS2()`：initDFS + 源点双向 DFS + min-depth 归一化 +
/// employ_regular_down_only_DFS_search。
///
/// 填充每个顶点的 `node_depth` / `depth` / `dfs_discovery_time` / `dfs_finish_time`。
pub fn run_dfs2(graph: &mut DiGraph) {
    // ---- initDFS() ----
    let mut colors = FxHashMap::default();
    let ids: Vec<i32> = graph.vertex_ids().to_vec();
    for &id in &ids {
        colors.insert(id, WHITE);
        let v = graph.get_vertex_mut(id).unwrap();
        v.node_depth = -1;
        v.depth = 0;
        v.dfs_discovery_time = -1;
        v.dfs_finish_time = -1;
    }
    let mut time: i32 = 0;

    // ---- 源点（in_degree==0）双向 DFS ----
    for &id in &ids {
        if graph.in_degree(id) == 0 && colors[&id] == WHITE {
            visit_vertex2(graph, &mut colors, &mut time, id, 0);
        }
    }

    // ---- min-depth 归一化（平移到 0）----
    let min_depth = ids
        .iter()
        .map(|&id| graph.get_vertex(id).unwrap().node_depth)
        .min()
        .unwrap_or(0);
    for &id in &ids {
        let v = graph.get_vertex_mut(id).unwrap();
        let adj = v.node_depth - min_depth;
        v.depth = adj;
        v.node_depth = adj;
    }

    // ---- employ_regular_down_only_DFS_search() ----
    employ_regular_down_only_dfs_search(graph);
}

/// `employ_regular_down_only_DFS_search()`：
/// 按 node_depth 升序重算 depth（get_depth_based_on_parents + down_only_visit_vertex）。
///
/// Java 用 `Collections.sort`（稳定排序，平局序 = HashSet 任意序）；这里用
/// (node_depth, 插入序) 的确定性序——depth 结果对平局序不敏感（见测试）。
fn employ_regular_down_only_dfs_search(graph: &mut DiGraph) {
    let mut order: Vec<i32> = graph.vertex_ids().to_vec();
    order.sort_by_key(|&id| graph.get_vertex(id).unwrap().node_depth);

    let mut colors = FxHashMap::default();
    for &id in &order {
        colors.insert(id, WHITE);
        let v = graph.get_vertex_mut(id).unwrap();
        v.node_depth = -1;
        v.depth = -1;
    }

    for &id in &order {
        if colors[&id] != WHITE {
            continue;
        }
        let depth = if graph.get_predecessors(id).is_empty() {
            0
        } else {
            get_depth_based_on_parents(graph, id)
        };
        {
            let v = graph.get_vertex_mut(id).unwrap();
            v.depth = depth;
            v.node_depth = depth;
        }
        colors.insert(id, GRAY);
        for s in graph.get_successors(id).to_vec() {
            down_only_visit_vertex(graph, &mut colors, s, 1);
        }
        colors.insert(id, BLACK);
    }
}

/// `get_depth_based_on_parents(v)`：max(前驱 nodeDepth) + 1。
fn get_depth_based_on_parents(graph: &DiGraph, id: i32) -> i32 {
    let mut depth = -1;
    for &p in graph.get_predecessors(id) {
        depth = depth.max(graph.get_vertex(p).unwrap().node_depth);
    }
    depth + 1
}

/// `down_only_visit_vertex(v, d)`：非 WHITE 或有 WHITE 前驱则延后；否则沿后继下行。
fn down_only_visit_vertex(graph: &mut DiGraph, colors: &mut FxHashMap<i32, u8>, id: i32, d: i32) {
    if colors[&id] != WHITE {
        return;
    }
    if has_white_parent(graph, colors, id) {
        // delay
        return;
    }
    colors.insert(id, GRAY);
    {
        let v = graph.get_vertex_mut(id).unwrap();
        v.depth = d;
        v.node_depth = d;
    }
    let succs = graph.get_successors(id).to_vec();
    for s in succs {
        down_only_visit_vertex(graph, colors, s, d + 1);
    }
    colors.insert(id, BLACK);
}

/// `has_white_parent(v)`：任一前驱为 WHITE。
fn has_white_parent(graph: &DiGraph, colors: &FxHashMap<i32, u8>, id: i32) -> bool {
    graph
        .get_predecessors(id)
        .iter()
        .any(|&p| colors[&p] == WHITE)
}

/// `visitVertex2(v, rec)`：双向递归（先后继 rec+1，再前驱 rec−1）。
///
/// GRAY 时记 discovery，BLACK 时记 finish；前后邻接递归全部完成后才出递归。
fn visit_vertex2(
    graph: &mut DiGraph,
    colors: &mut FxHashMap<i32, u8>,
    time: &mut i32,
    id: i32,
    rec: i32,
) {
    if colors[&id] != WHITE {
        return; // already visited
    }

    {
        let v = graph.get_vertex_mut(id).unwrap();
        v.depth = rec;
        v.node_depth = rec;
    }
    colors.insert(id, GRAY);
    *time += 1;
    graph.get_vertex_mut(id).unwrap().dfs_discovery_time = *time;

    // 先访问全部后继（rec+1）
    let succs = graph.get_successors(id).to_vec();
    for s in succs {
        visit_vertex2(graph, colors, time, s, rec + 1);
    }
    // 再访问全部前驱（rec−1）——双向
    let preds = graph.get_predecessors(id).to_vec();
    for p in preds {
        visit_vertex2(graph, colors, time, p, rec - 1);
    }

    colors.insert(id, BLACK);
    *time += 1;
    graph.get_vertex_mut(id).unwrap().dfs_finish_time = *time;
}

/// finish time 降序的顶点 id 序（DAG 拓扑序；平局按 id 升序保证确定性）。
///
/// 注意：visitVertex2 是双向 DFS，经前驱递归访问到的分支 finish 嵌套在内层，
/// finish 降序在含这类分支的图上不是严格拓扑序——与 Java 语义一致，调用方
/// （路径搜索）自行依赖 finish 相对序。
pub fn get_topological_order(graph: &DiGraph) -> Vec<i32> {
    let mut ids: Vec<i32> = graph.vertex_ids().to_vec();
    ids.sort_by(|&a, &b| {
        let fa = graph.get_vertex(a).unwrap().dfs_finish_time;
        let fb = graph.get_vertex(b).unwrap().dfs_finish_time;
        fb.cmp(&fa).then(a.cmp(&b))
    });
    ids
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::SeqVertex;

    fn build(edges: &[(i32, i32)]) -> DiGraph {
        let mut g = DiGraph::new();
        let mut ids: Vec<i32> = edges.iter().flat_map(|&(a, b)| [a, b]).collect();
        ids.sort();
        ids.dedup();
        for id in ids {
            g.add_vertex(SeqVertex::new(id, "ACGT"));
        }
        for &(a, b) in edges {
            g.add_edge(a, b, crate::graph::SimpleEdge::new(1.0, a, b));
        }
        g
    }

    fn d(g: &DiGraph, id: i32) -> i32 {
        g.get_vertex(id).unwrap().node_depth
    }
    fn disc(g: &DiGraph, id: i32) -> i32 {
        g.get_vertex(id).unwrap().dfs_discovery_time
    }
    fn fin(g: &DiGraph, id: i32) -> i32 {
        g.get_vertex(id).unwrap().dfs_finish_time
    }

    #[test]
    fn diamond_graph_times_and_depths() {
        // A(1)->B(2), A(1)->C(3), B->D(4), C->D
        let mut g = build(&[(1, 2), (1, 3), (2, 4), (3, 4)]);
        run_dfs2(&mut g);

        // 手推：visit(1,0) disc1 → succ 2 disc2 → succ 4 disc3
        //   → preds [2(BLACK 跳过? 尚 GRAY), 3]：2 GRAY 跳过；3 WHITE disc4
        //     → succ 4 非白跳过; preds [1 GRAY] 跳过; fin3=5
        //   → 4 fin=6; 2 preds [1] 跳过 fin=7; 1 preds 无 fin=8
        assert_eq!((disc(&g, 1), fin(&g, 1)), (1, 8));
        assert_eq!((disc(&g, 2), fin(&g, 2)), (2, 7));
        assert_eq!((disc(&g, 3), fin(&g, 3)), (4, 5));
        assert_eq!((disc(&g, 4), fin(&g, 4)), (3, 6));

        // depth 归一化（min=0 无平移）+ down-only 重算
        assert_eq!(
            (
                d(&g, 1),
                d(&g, 2),
                d(&g, 3),
                d(&g, 4),
                g.get_vertex(4).unwrap().depth
            ),
            (0, 1, 1, 2, 2)
        );

        // finish 降序拓扑序
        assert_eq!(get_topological_order(&g), vec![1, 2, 4, 3]);
    }

    #[test]
    fn bidirectional_recursion_reaches_upstream_of_visited() {
        // 1(S)->4(M)；P(2)→M；Q(3)→P、P→Q（P,Q 成环）。唯一 in_degree==0 顶点是 1。
        // P、Q 只能通过 M 的前驱递归被编号（双向递归的作用）。
        let mut g = build(&[(1, 4), (2, 4), (3, 2), (2, 3)]);
        run_dfs2(&mut g);

        // phase1: 1 disc1 → 4 disc2 → preds[1 GRAY,2 WHITE]: 2 disc3
        //   → succs[4 GRAY,3]: 3 disc4 → succs[2 GRAY] pred[2] fin3=5;
        //   2 fin=6; 4 fin=7; 1 fin=8
        assert_eq!((disc(&g, 2), disc(&g, 3)), (3, 4));
        assert_eq!(fin(&g, 1), 8);

        // down-only：归一化深度 1:0 2:0 3:1 4:1；序 (depth,id) = 1,2,3,4
        // v=1: depth0（无前驱）；down(4,1) 被前驱 2 的 WHITE 阻断
        // v=2: preds=[3] 但 3 已重置 nodeDepth=-1 → depth = -1+1 = 0；
        //      down(4,1)（前驱 1 BLACK、2 GRAY → 通过）、down(3,1)
        assert_eq!((d(&g, 1), d(&g, 2), d(&g, 3), d(&g, 4)), (0, 0, 1, 1));
    }

    #[test]
    fn has_white_parent_blocks_and_delays() {
        // S(1)->T(3)，S->B(2)，B->T：T 有两个前驱；S 处理时 B 仍 WHITE → T 阻断，
        // 之后按序 B 先于 T 得到深度，T 由 B 的 down_only 访问补上。
        let mut g = build(&[(1, 3), (1, 2), (2, 3)]);
        run_dfs2(&mut g);
        // phase1: succ 序 [3,2] → 3 先访问（rec1），3 的前驱递归发现 2（rec0）
        // down-only 序 (depth,id) = 1,2,3：v=1 down(3,1) 被 2 的 WHITE 阻断;
        // down(2,1) 通过 → 2 depth1 → down(3,2)（前驱全非白）
        assert_eq!((d(&g, 1), d(&g, 2), d(&g, 3)), (0, 1, 2));
        assert_eq!(g.get_vertex(3).unwrap().depth, 2);
    }

    #[test]
    fn cycle_without_source_still_numbered_via_down_only() {
        // 无 in_degree==0 顶点：A(1)→B(2)→C(3)→A，C→D(4)。
        // phase1 不访问任何点（node_depth 全 -1 → 归一化全 0），编号完全来自
        // down-only 阶段。
        let mut g = build(&[(1, 2), (2, 3), (3, 1), (3, 4)]);
        run_dfs2(&mut g);
        // 序（全 depth0, id 升序）: 1,2,3,4
        // v=1: preds=[3] 但 3 已重置 nodeDepth=-1 → depth 0；down(2,1)→down(3,2)
        //   → down(1,3) 跳过 / down(4,3)；其余 BLACK 跳过
        assert_eq!((d(&g, 1), d(&g, 2), d(&g, 3), d(&g, 4)), (0, 1, 2, 3));
        // discovery/finish 未被记录（phase1 无源点）
        assert_eq!(disc(&g, 1), -1);
        assert_eq!(fin(&g, 1), -1);
    }

    #[test]
    fn min_depth_normalization_shifts_to_zero() {
        // 前驱递归链产生负深度：1→2、3→2、4→3（4 是源但插入序靠后，
        // 1 先被访问；2(rec1) → 前驱 3(rec0) → 前驱 4(rec-1)）。
        let mut g = build(&[(1, 2), (4, 3), (3, 2)]);
        run_dfs2(&mut g);
        // phase1 深度 1:0 2:1 3:0 4:-1，min=-1 → 平移后 1:1 2:2 3:1 4:0
        assert_eq!((disc(&g, 1), disc(&g, 3), disc(&g, 4)), (1, 3, 4));
        assert_eq!(
            (fin(&g, 1), fin(&g, 2), fin(&g, 3), fin(&g, 4)),
            (8, 7, 6, 5)
        );
        // down-only 重算（重置后序 4,1,3,2）：
        //   v=4: 无前驱 depth0 → down(3,1) → down(2,2) 被前驱 1 的 WHITE 阻断
        //   v=1: 无前驱 depth0 → down(2,1)（前驱 1 GRAY、3 BLACK）→ 通过
        assert_eq!((d(&g, 1), d(&g, 2), d(&g, 3), d(&g, 4)), (0, 1, 1, 0));
        assert_eq!(g.get_vertex(4).unwrap().depth, 0);
    }

    #[test]
    fn topological_order_finish_descending() {
        let mut g = build(&[(1, 2), (2, 3), (1, 3)]);
        run_dfs2(&mut g);
        // 1 disc1 fin6; 2 disc2 fin3; 3 disc3 fin... visit(1,0)→2(disc2)→3(disc3,fin4)→2 fin5→1 fin6
        assert_eq!(
            get_topological_order(&g),
            vec![1, 2, 3] // fin: 1=6, 2=5, 3=4
        );
    }
}
