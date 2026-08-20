//! P4 Butterfly：read-path 重叠图（POG）——构建 / PE links / 破环 / 转 SeqVertex DAG / zipping。
//!
//! 镜像 trinityrnaseq v2.15.2：
//! - `TransAssembly_allProbPaths.java`：`create_DAG_from_OverlapLayout`（L1617）、
//!   `populate_pairpaths_and_readsupport`（L8501）、`remove_containments`（L3460）、
//!   `find_dispersed_repeat_nodes`（L1865）、`construct_path_overlap_graph`（L3321）、
//!   `addPairPathsToOverlapGraph`（L1771）、`break_cycles_in_path_overlap_graph`（L3147）、
//!   `convert_path_DAG_to_SeqVertex_DAG`（L2146）、`DFS_add_path_to_graph`（L2842）、
//!   zipping 调度（L2220-2330）、`zip_up`/`zip_down`（L2573/2621）、
//!   `attempt_zip_merge_SeqVertices`（L2671）、`destroy_unzipped_duplicates_above`（L2393）、
//!   `update_PairPaths_using_overlapDAG_refined_paths`（L1935）、
//!   `get_all_possible_updated_path_mappings`（L2116）、POG `writeDotFile`（L3280）。
//! - `Path.java`：`pathB_extends_pathA_allowRepeats`（L553）、
//!   `pathA_contains_pathB_allowRepeats`（L498）、`getRepeatNodesAndCounts`（L158）。
//! - `PathWithOrig.java`：`align_path_by_orig_id`。
//! - `TopologicalSort.java`：Kahn 拓扑排序 + 深度赋值。
//!
//! **迭代序说明**：Java 的 POG 节点编号（PN#）、PairPath 遍历序与 SeqVertex 新 id
//! 都取决于 JUNG/HashMap 迭代序，跨实现不可复现。本模块用确定的插入序；对拍以
//! **结构比较**为准（POG：路径内容多重集 + 内容对边集；SeqVertex 图：orig-id
//! 节点多重集 + orig-id 边多重集），见 `tests/pog_c0.rs`。

use std::collections::VecDeque;

use rustc_hash::{FxHashMap, FxHashSet};

use crate::context::BflyContext;
use crate::graph::{DiGraph, SeqVertex, SimpleEdge};
use crate::pair_paths::{PairPath, SuffStats};

// ---------------------------------------------------------------------------
// 基础结构（PathOverlap / Path / PathWithOrig / SimplePathNodeEdge）
// ---------------------------------------------------------------------------

/// `PathOverlap.java`：match_score = 非重复节点匹配数，match_length = 比对步数。
/// 发布版 jar 额外携带 `path_A/path_B/idx_start_A/idx_start_B/A_contains_B`
/// （javap 反编译证实），其中 `idx_start_A` = 匹配在 pathA 的起点、
/// `idx_start_B` 恒 0（B 从自身开头延伸 A）；DFS 跨路径边依赖这两个下标。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PathOverlap {
    pub match_score: i32,
    pub match_length: usize,
    pub idx_start_a: usize,
    pub idx_start_b: usize,
}

/// `Path.java`：POG 节点 = 一条 read 路径。`pn_id` 从 1 起（`pathNodeCounter`）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PathNode {
    pub pn_id: i32,
    pub vertices: Vec<i32>,
}

/// `PathWithOrig.java`：新图 id 链 + 原图 id 链（`align_path_by_orig_id` 的模板）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PathWithOrig {
    pub vertex_id_list: Vec<i32>,
    pub orig_vertex_id_list: Vec<i32>,
}

impl PathWithOrig {
    /// `PathWithOrig(path)`：按 orig 查询函数把 id 链翻译成 orig id 链。
    pub fn from_path(path: &[i32], orig_of: &dyn Fn(i32) -> i32) -> Self {
        Self {
            vertex_id_list: path.to_vec(),
            orig_vertex_id_list: path.iter().map(|&id| orig_of(id)).collect(),
        }
    }

    /// `align_path_by_orig_id(template)`：本路径能否完整贴到模板的某个起点上；
    /// 能则返回按模板 id 重写的路径（Java 在 mismatch 处也 append 模板 id —— 原样）。
    pub fn align_path_by_orig_id(&self, template: &PathWithOrig) -> Option<PathWithOrig> {
        for i in 0..template.orig_vertex_id_list.len() {
            if self.orig_vertex_id_list.first() != template.orig_vertex_id_list.get(i) {
                continue;
            }
            let mut restructured: Vec<i32> = Vec::new();
            let mut matches = 0usize;
            let (mut tp, mut sp) = (i, 0usize);
            while tp < template.orig_vertex_id_list.len() && sp < self.orig_vertex_id_list.len() {
                restructured.push(template.vertex_id_list[tp]);
                if self.orig_vertex_id_list[sp] == template.orig_vertex_id_list[tp] {
                    matches += 1;
                } else {
                    break;
                }
                tp += 1;
                sp += 1;
            }
            if matches == self.orig_vertex_id_list.len() {
                return Some(PathWithOrig {
                    vertex_id_list: restructured,
                    orig_vertex_id_list: self.orig_vertex_id_list.clone(),
                });
            }
        }
        None
    }
}

/// `SimplePathNodeEdge`：POG 边（weight = match_score；num_loops 破环用）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SimplePathNodeEdge {
    pub weight: i32,
    pub num_loops: i32,
}

// ---------------------------------------------------------------------------
// POG 图
// ---------------------------------------------------------------------------

/// 路径重叠图：节点索引（= `path_list` 插入序）+ 边表（保留插入序）。
#[derive(Debug, Clone, Default)]
pub struct PathOverlapGraph {
    pub nodes: Vec<PathNode>,
    /// (from_idx, to_idx) → 边；from = 被延伸（precursor），to = 延伸者。
    pub edges: FxHashMap<(usize, usize), SimplePathNodeEdge>,
    /// 边插入序（迭代确定性）。
    pub edge_order: Vec<(usize, usize)>,
}

impl PathOverlapGraph {
    fn add_edge(&mut self, from: usize, to: usize, weight: i32) {
        if !self.edges.contains_key(&(from, to)) {
            self.edge_order.push((from, to));
        }
        self.edges.insert(
            (from, to),
            SimplePathNodeEdge {
                weight,
                num_loops: 0,
            },
        );
    }

    /// 后继（边插入序）。
    pub fn successors(&self, from: usize) -> Vec<usize> {
        self.edge_order
            .iter()
            .filter(|(f, _)| *f == from)
            .map(|(_, t)| *t)
            .collect()
    }

    /// 入度。
    pub fn in_degree(&self, idx: usize) -> usize {
        self.edge_order.iter().filter(|(_, t)| *t == idx).count()
    }

    /// `_POG.dot`（`writeDotFile` 的 POG 变体，L3280）：`PN#` label + 裸边。
    pub fn write_dot(&self) -> String {
        let mut s = String::from("digraph G {\n");
        for n in &self.nodes {
            s.push_str(&format!("\tPN{} [label=\"PN{}\"]\n", n.pn_id, n.pn_id));
            for t in self.successors(n.pn_id as usize - 1) {
                s.push_str(&format!("\tPN{}->PN{}\n", n.pn_id, self.nodes[t].pn_id));
            }
        }
        s.push_str("}\n");
        s
    }

    /// BFS 可达性：返回一条最短路径的边序列（Java DijkstraShortestPath 仅用于
    /// null 判定 + 任一最短路；环边集判定等价）。
    fn shortest_path(&self, from: usize, to: usize) -> Option<Vec<(usize, usize)>> {
        let n = self.nodes.len();
        let mut prev: Vec<Option<usize>> = vec![None; n];
        let mut visited = vec![false; n];
        let mut queue = VecDeque::new();
        visited[from] = true;
        queue.push_back(from);
        while let Some(u) = queue.pop_front() {
            if u == to {
                break;
            }
            for v in self.successors(u) {
                if !visited[v] {
                    visited[v] = true;
                    prev[v] = Some(u);
                    queue.push_back(v);
                }
            }
        }
        if !visited[to] {
            return None;
        }
        let mut path = Vec::new();
        let mut cur = to;
        while let Some(p) = prev[cur] {
            path.push((p, cur));
            cur = p;
        }
        path.reverse();
        Some(path)
    }
}

// ---------------------------------------------------------------------------
// Path.java 移植
// ---------------------------------------------------------------------------

/// `Path.pathB_extends_pathA_allowRepeats`（L553）：把 pathB 的前缀贴到 pathA 的
/// 某个后缀起点；repeat 节点不计分但仍要求相等。取 match_score 最高的起点
/// （严格 `>`，平局保留先遇到者——i 从 pathA 末尾向前扫）。
pub fn path_b_extends_path_a_allow_repeats(
    path_b: &[i32],
    path_a: &[i32],
    repeat_node_ids: &FxHashSet<i32>,
) -> PathOverlap {
    let mut best_match_count = -1;
    let mut path_overlap = PathOverlap {
        match_score: 0,
        match_length: 0,
        idx_start_a: 0,
        idx_start_b: 0,
    };

    if path_b.is_empty() || path_a.is_empty() {
        return path_overlap;
    }

    // Java: for (int i = pathA.size()-1; i >= pathA.size()-pathB.size() && i >= 0; i--)
    let start = path_a.len() as isize - path_b.len() as isize;
    let mut i = path_a.len() as isize - 1;
    while i >= start && i >= 0 {
        let iu = i as usize;
        if path_b[0] == path_a[iu] {
            let mut matches = 0i32;
            let mut overlap_len = 0usize;
            let (mut bp, mut ap) = (0usize, iu);
            while ap < path_a.len() && bp < path_b.len() {
                overlap_len += 1;
                if path_b[bp] == path_a[ap] {
                    if !repeat_node_ids.contains(&path_a[ap]) {
                        matches += 1;
                    }
                } else {
                    matches = -1;
                    break;
                }
                ap += 1;
                bp += 1;
            }
            if matches > best_match_count {
                best_match_count = matches;
                path_overlap.match_length = overlap_len;
                path_overlap.match_score = matches;
                // jar：new PathOverlap(pathA, pathB, i, 0, matches, overlap_len)
                path_overlap.idx_start_a = iu;
                path_overlap.idx_start_b = 0;
            }
        }
        i -= 1;
    }
    path_overlap
}

/// `Path.pathA_contains_pathB_allowRepeats`（L498）：pathB 是 pathA 的连续子串。
pub fn path_a_contains_path_b_allow_repeats(path_a: &[i32], path_b: &[i32]) -> bool {
    if path_b.is_empty() {
        return false; // Java 会对空 pathB 的 get(0) 抛错；调用方保证非空
    }
    for i in 0..path_a.len() {
        if path_b[0] != path_a[i] {
            continue;
        }
        let mut matches = 0usize;
        let (mut bp, mut ap) = (0usize, i);
        while ap < path_a.len() && bp < path_b.len() {
            if path_b[bp] == path_a[ap] {
                matches += 1;
            } else {
                break;
            }
            ap += 1;
            bp += 1;
        }
        if matches == path_b.len() {
            return true;
        }
    }
    false
}

/// `Path.getRepeatNodesAndCounts`（L158）：路径内出现 >1 次的节点 → 次数。
pub fn get_repeat_nodes_and_counts(path: &[i32]) -> FxHashMap<i32, i32> {
    let mut counts: FxHashMap<i32, i32> = FxHashMap::default();
    for &id in path {
        *counts.entry(id).or_insert(0) += 1;
    }
    counts.retain(|_, c| *c > 1);
    counts
}

// ---------------------------------------------------------------------------
// create_DAG_from_OverlapLayout 子步骤
// ---------------------------------------------------------------------------

/// `populate_pairpaths_and_readsupport`（L8501，非 MAKE_PE_SE 分支）：
/// trimSinkNodes 后收进 set；**重复 pp 的 support 覆盖而非累加**（Java quirk）。
/// 返回确定序（start_id 升序，再按 path 内容序）的 (pp, support) 列表。
pub fn populate_pairpaths_and_readsupport(
    combined_read_hash: &FxHashMap<i32, FxHashMap<PairPath, i64>>,
) -> Vec<(PairPath, i64)> {
    let mut keys: Vec<i32> = combined_read_hash.keys().copied().collect();
    keys.sort_unstable();
    let mut support: FxHashMap<PairPath, i64> = FxHashMap::default();
    let mut order: Vec<PairPath> = Vec::new();
    for k in keys {
        let inner = &combined_read_hash[&k];
        let mut pps: Vec<&PairPath> = inner.keys().collect();
        pps.sort_by(|a, b| (&a.path1, &a.path2).cmp(&(&b.path1, &b.path2)));
        for pp in pps {
            let s = inner[pp];
            let trimmed = pp.trim_sink_nodes();
            if !support.contains_key(&trimmed) {
                order.push(trimmed.clone());
            }
            support.insert(trimmed, s); // Java put() —— 覆盖
        }
    }
    order
        .into_iter()
        .map(|pp| (pp.clone(), support[&pp]))
        .collect()
}

/// `remove_containments`（L3460）：降序输入；先选入的路径若完全包含当前路径则
/// 剔除当前路径；**所有包含者逐个记录**（一个路径可有多个 container）。
pub fn remove_containments(
    paths: &[Vec<i32>],
    contained_path_to_containers: &mut FxHashMap<Vec<i32>, Vec<Vec<i32>>>,
) -> Vec<Vec<i32>> {
    let mut noncontained: Vec<Vec<i32>> = Vec::new();
    for path in paths {
        let mut contained = false;
        for chosen in &noncontained {
            if path_a_contains_path_b_allow_repeats(chosen, path) {
                contained = true;
                contained_path_to_containers
                    .entry(path.clone())
                    .or_default()
                    .push(chosen.clone());
            }
        }
        if !contained {
            noncontained.push(path.clone());
        }
    }
    noncontained
}

/// `find_dispersed_repeat_nodes`（L1865）：出现在 >= `MIN_OCCURRENCE_REPEAT_NODE`(10)
/// 条（非包含）路径中的节点（每条路径至多计一次）。
pub fn find_dispersed_repeat_nodes(paths: &[Vec<i32>]) -> FxHashSet<i32> {
    const MIN_OCCURRENCE_REPEAT_NODE: i32 = 10;
    let mut counter: FxHashMap<i32, i32> = FxHashMap::default();
    for path in paths {
        for id in path.iter().copied().collect::<FxHashSet<i32>>() {
            *counter.entry(id).or_insert(0) += 1;
        }
    }
    counter
        .into_iter()
        .filter(|(_, c)| *c >= MIN_OCCURRENCE_REPEAT_NODE)
        .map(|(id, _)| id)
        .collect()
}

/// `construct_path_overlap_graph`（L3321）：全对 (i, j≠i) 判定
/// pathB=paths[i] extends pathA=paths[j]；score>0 → 边 j→i（weight=score）。
/// `store_best_extension_match_only = false`：**所有**延伸关系都建边。
/// repeat 集 = dispersed ∪ 局部（路径内出现 >2 次的节点）。
pub fn construct_path_overlap_graph(
    path_list: &[Vec<i32>],
    path_matches: &mut FxHashMap<(usize, usize), PathOverlap>,
    dispersed_repeat_nodes: &FxHashSet<i32>,
) -> PathOverlapGraph {
    let mut pog = PathOverlapGraph {
        nodes: path_list
            .iter()
            .enumerate()
            .map(|(idx, p)| PathNode {
                pn_id: idx as i32 + 1,
                vertices: p.clone(),
            })
            .collect(),
        ..Default::default()
    };

    let mut repeat_node_ids: FxHashSet<i32> = dispersed_repeat_nodes.clone();
    let max_internal_repeat_count = 2;
    for path in path_list {
        for (id, count) in get_repeat_nodes_and_counts(path) {
            if count > max_internal_repeat_count {
                repeat_node_ids.insert(id);
            }
        }
    }

    for i in 0..path_list.len() {
        let mut best_precursor_j_indices: Vec<usize> = Vec::new();
        for j in 0..path_list.len() {
            if i == j {
                continue;
            }
            let po =
                path_b_extends_path_a_allow_repeats(&path_list[i], &path_list[j], &repeat_node_ids);
            if po.match_score <= 0 {
                continue;
            }
            path_matches.insert((j, i), po); // i extends j
            best_precursor_j_indices.push(j);
        }
        for j in best_precursor_j_indices {
            let w = path_matches[&(j, i)].match_score;
            pog.add_edge(j, i, w);
        }
    }
    pog
}

/// `addPairPathsToOverlapGraph`（L1771）：**Java v2.15.2 中被 `if (true) return`
/// 整体关闭**——恒返回空集（`_POG.dot` 与 `_POG.PE_links_added.dot` 恒相同）。
pub fn add_pair_paths_to_overlap_graph(
    _pog: &mut PathOverlapGraph,
    _pair_path_to_read_support: &FxHashMap<PairPath, i64>,
    _contained_path_to_containers: &FxHashMap<Vec<i32>, Vec<Vec<i32>>>,
) -> FxHashSet<(usize, usize)> {
    FxHashSet::default()
}

/// `break_cycles_in_path_overlap_graph`（L3147）：逐边找回到起点的路径构成环
/// 边集；统计每条边涉及的环数；按 **环数降序、weight 升序** 优先删环数最多的边，
/// 每删一条同步消解包含它的环。返回是否删过边（true → 下一轮）。
pub fn break_cycles_in_path_overlap_graph(pog: &mut PathOverlapGraph) -> bool {
    // 收集环（边集合去重）
    let mut cur_loops: FxHashSet<Vec<(usize, usize)>> = FxHashSet::default();
    for p in 0..pog.nodes.len() {
        for s in pog.successors(p) {
            if let Some(mut loop_path) = pog.shortest_path(s, p) {
                loop_path.insert(0, (p, s)); // 补 p->s 边构成完整环
                let mut edge_set: Vec<(usize, usize)> = loop_path;
                edge_set.sort_unstable();
                edge_set.dedup();
                cur_loops.insert(edge_set);
            }
        }
    }
    if cur_loops.is_empty() {
        return false;
    }

    let mut loop_sets: Vec<Vec<(usize, usize)>> = cur_loops.into_iter().collect();
    loop_sets.sort();
    for lp in &loop_sets {
        for e in lp {
            pog.edges.get_mut(e).unwrap().num_loops += 1;
        }
    }

    let mut res = false;
    loop {
        if loop_sets.is_empty() || pog.edge_order.is_empty() {
            break;
        }
        // 优先队列弹出：num_loops 降序，weight 升序，平局按边插入序
        let mut best: Option<(usize, usize)> = None;
        let mut best_key: Option<(i32, i32, usize)> = None;
        for (pos, &e) in pog.edge_order.iter().enumerate() {
            let edge = &pog.edges[&e];
            if edge.num_loops <= 0 {
                continue;
            }
            let key = (edge.num_loops, edge.weight, pos);
            if best_key.is_none_or(|b| {
                key.0 > b.0 || (key.0 == b.0 && (key.1 < b.1 || (key.1 == b.1 && key.2 < b.2)))
            }) {
                best_key = Some(key);
                best = Some(e);
            }
        }
        let Some(next_e) = best else { break };

        // 消解包含 next_e 的环（递减其中每条边的环计数）
        let mut remaining: Vec<Vec<(usize, usize)>> = Vec::new();
        for lp in loop_sets.iter() {
            if lp.contains(&next_e) {
                for e2 in lp {
                    pog.edges.get_mut(e2).unwrap().num_loops -= 1;
                }
            } else {
                remaining.push(lp.clone());
            }
        }
        loop_sets = remaining;

        pog.edges.remove(&next_e);
        pog.edge_order.retain(|e| *e != next_e);
        res = true;
    }
    res
}

// ---------------------------------------------------------------------------
// 拓扑排序（TopologicalSort.java）+ DOT（SeqVertex 变体）
// ---------------------------------------------------------------------------

/// Kahn 拓扑排序（FIFO，起始集为插入序）。
/// 返回拓扑序 id 列表；有环返回 None（Java 抛异常）。
pub fn topo_sort_seq_vertices_dag(graph: &DiGraph) -> Option<Vec<i32>> {
    let mut in_deg: FxHashMap<i32, usize> = FxHashMap::default();
    for &id in graph.vertex_ids() {
        in_deg.insert(id, graph.in_degree(id));
    }
    let mut queue: VecDeque<i32> = graph
        .vertex_ids()
        .iter()
        .copied()
        .filter(|&id| in_deg[&id] == 0)
        .collect();
    let mut order: Vec<i32> = Vec::new();
    while let Some(n) = queue.pop_front() {
        order.push(n);
        for &m in graph.get_successors(n) {
            let d = in_deg.get_mut(&m).unwrap();
            *d -= 1;
            if *d == 0 {
                queue.push_back(m);
            }
        }
    }
    if order.len() != graph.vertex_count() {
        return None;
    }
    Some(order)
}

/// 排序后为节点赋 `node_depth = 序号`（Java depth 从 -1 自增 → 首节点 0）。
pub fn assign_depths(graph: &mut DiGraph, order: &[i32]) {
    // Java TopologicalSort L119-120：`v.setDepth(depth); v.setNodeDepth(depth);`
    // —— 两个深度字段同步为拓扑序号（BflyQueue 的 SeqVertexNodeDepthComparator
    // 读的是 getDepth()/_depth）。
    for id in graph.vertex_ids().to_vec() {
        let v = graph.get_vertex_mut(id).unwrap();
        v.node_depth = -1;
        v.depth = -1;
    }
    for (d, &id) in order.iter().enumerate() {
        let v = graph.get_vertex_mut(id).unwrap();
        v.node_depth = d as i32;
        v.depth = d as i32;
    }
}

fn short_seq(v: &SeqVertex) -> String {
    let name = v.name.as_str();
    let res = if name.len() > 30 {
        format!("{}...{}", &name[..10], &name[name.len() - 10..])
    } else {
        name.to_string()
    };
    format!("{}:W{}", res, v.get_weight_avg())
}

/// `writeDotFile(graph, file, name, /*printFullSeq=*/false)`（L12868）：
/// 短序列 label（getShortSeqWID + [L:len][T:discoveryTime]）+ 权重边。
pub fn write_seqvertex_dot(graph: &DiGraph) -> String {
    let mut s = String::from("digraph G {\n");
    for &id in graph.vertex_ids() {
        let v = graph.get_vertex(id).unwrap();
        let mut label = format!("{}(V{}", short_seq(v), v.id);
        if v.orig_butterfly_id != v.id {
            label.push_str(&format!("_{}", v.orig_butterfly_id));
        }
        label.push_str(&format!(
            "_D{})[L:{}][T:{}]",
            v.node_depth,
            v.name.len(),
            v.dfs_discovery_time
        ));
        if v.get_weight_avg() > 25 {
            label.push_str(" ,style=bold,color=\"#AF0000\"");
        }
        s.push_str(&format!("\t{} [label=\"{}\"]\n", v.id, label));
        for &to in graph.get_successors(id) {
            let weight = crate::graph_io::java_math_round(graph.find_edge(id, to).unwrap().weight);
            let edge_style = if weight > 20 {
                format!("[style=bold,label={},color=\"#AF0000\"]", weight)
            } else {
                format!("[label={}]", weight)
            };
            s.push_str(&format!("\t{}->{}{}\n", id, to, edge_style));
        }
    }
    s.push_str("}\n");
    s
}

// ---------------------------------------------------------------------------
// zip 状态（__tmp_compressed_vertices + is_replacement_vertex）
// ---------------------------------------------------------------------------

/// SeqVertex 的 POG 附着状态（Java `__tmp_compressed_vertices` /
/// `is_replacement_vertex`），独立于 `DiGraph` 存放。
#[derive(Debug, Clone, Default)]
pub struct ZipState {
    pub tmp_compressed_vertices: FxHashMap<i32, Vec<i32>>,
    pub is_replacement: FxHashSet<i32>,
}

impl ZipState {
    fn compressed(&self, id: i32) -> &[i32] {
        self.tmp_compressed_vertices
            .get(&id)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// 读取某节点的 compressed 列表（测试 / 下游用）。
    pub fn compressed_public(&self, id: i32) -> &[i32] {
        self.compressed(id)
    }

    fn add_compressed(&mut self, id: i32, vals: Vec<i32>) {
        self.tmp_compressed_vertices
            .entry(id)
            .or_default()
            .extend(vals);
    }
}

/// `init_replacement_vertices`：每轮 zip 前清 replacement 标记。
pub fn init_replacement_vertices(zip: &mut ZipState) {
    zip.is_replacement.clear();
}

/// `attempt_zip_merge_SeqVertices`（L2671）：同 orig id 的一组节点合并成新节点。
/// 约束：parent/child 均非空时要求 `max(parent_depths) < min(child_depths)`；
/// 任一 parent/child 已是 replacement → 推迟（返回 0）。新节点深度 = min/max。
pub fn attempt_zip_merge_seq_vertices(
    targets: &[i32],
    graph: &mut DiGraph,
    zip: &mut ZipState,
    ctx: &mut BflyContext,
    dir_min: bool,
) -> usize {
    let replacement_vertex_id = ctx.get_next_id();

    let mut parent_vertices: Vec<i32> = Vec::new(); // 去重（保插入序；Java HashSet）
    let mut child_vertices: Vec<i32> = Vec::new();
    let mut edges_to_delete: Vec<(i32, i32)> = Vec::new();
    let mut parent_depths: Vec<i32> = Vec::new();
    let mut child_depths: Vec<i32> = Vec::new();
    let mut target_depths: Vec<i32> = Vec::new();

    for &v in targets {
        target_depths.push(graph.get_vertex(v).unwrap().node_depth);
        for &p in graph.get_predecessors(v) {
            if zip.is_replacement.contains(&p) {
                return 0; // delay till next round
            }
            if !parent_vertices.contains(&p) {
                parent_vertices.push(p);
            }
            parent_depths.push(graph.get_vertex(p).unwrap().node_depth);
            edges_to_delete.push((p, v));
        }
        for &c in graph.get_successors(v) {
            if zip.is_replacement.contains(&c) {
                return 0;
            }
            if !child_vertices.contains(&c) {
                child_vertices.push(c);
            }
            child_depths.push(graph.get_vertex(c).unwrap().node_depth);
            edges_to_delete.push((v, c));
        }
    }

    // 深度顺序约束（parent/child 都非空时才检查 —— Java 同）
    if !parent_depths.is_empty() && !child_depths.is_empty() {
        let max_p = *parent_depths.iter().max().unwrap();
        let min_c = *child_depths.iter().min().unwrap();
        if max_p >= min_c {
            return 0; // 会破坏相对顺序
        }
    }

    for e in &edges_to_delete {
        graph.remove_edge(e.0, e.1);
    }

    // 模板字段先拷出（同 orig → 同序列），再删目标节点
    let t0 = graph.get_vertex(targets[0]).unwrap().clone();
    let mut merged_vertex_ids: Vec<i32> = Vec::new();
    for &v in targets {
        merged_vertex_ids.push(v);
        merged_vertex_ids.extend_from_slice(zip.compressed(v));
        graph.remove_vertex(v);
        zip.tmp_compressed_vertices.remove(&v);
    }

    // new SeqVertex(id, v)：拷贝 name/weights/orig；随后 L2786-2787
    // `setDepth(replacement_vertex_depth); setNodeDepth(...)` 把两个深度字段
    // 都覆盖为 min/max 目标深度（构造器的 depth+len 被丢弃）。
    let mut new_v = SeqVertex::new(replacement_vertex_id, t0.name.clone());
    new_v.weights = t0.weights.clone();
    new_v.orig_butterfly_id = t0.orig_butterfly_id;
    let depth = if dir_min {
        *target_depths.iter().min().unwrap()
    } else {
        *target_depths.iter().max().unwrap()
    };
    new_v.node_depth = depth;
    new_v.depth = depth;
    graph.add_vertex(new_v);

    for &p in &parent_vertices {
        graph.add_edge(
            p,
            replacement_vertex_id,
            SimpleEdge::new(1.0, p, replacement_vertex_id),
        );
    }
    for &c in &child_vertices {
        graph.add_edge(
            replacement_vertex_id,
            c,
            SimpleEdge::new(1.0, replacement_vertex_id, c),
        );
    }

    // 本轮局部环境冻结：邻接 replacement 节点不再参与合并
    for &p in &parent_vertices {
        zip.is_replacement.insert(p);
    }
    for &c in &child_vertices {
        zip.is_replacement.insert(c);
    }
    zip.is_replacement.insert(replacement_vertex_id);
    zip.tmp_compressed_vertices
        .insert(replacement_vertex_id, merged_vertex_ids);

    targets.len()
}

/// `zip_up`（L2573）：v 的多个**前驱**若同 orig id → 合并（depth 取 min）。
/// 任一前驱是 replacement → 整个 v 推迟到下一轮（Java return 0）。
pub fn zip_up(graph: &mut DiGraph, zip: &mut ZipState, ctx: &mut BflyContext, v: i32) -> usize {
    let pred_list: Vec<i32> = graph.get_predecessors(v).to_vec();
    if pred_list.len() <= 1 {
        return 0;
    }
    // orig id → 首现序分组
    let mut groups: Vec<(i32, Vec<i32>)> = Vec::new();
    for &p in &pred_list {
        if zip.is_replacement.contains(&p) {
            return 0; // Java：整个 v 推迟到下一轮
        }
        let orig = graph.get_vertex(p).unwrap().orig_butterfly_id;
        match groups.iter_mut().find(|(o, _)| *o == orig) {
            Some((_, g)) => g.push(p),
            None => groups.push((orig, vec![p])),
        }
    }
    let mut count = 0;
    for (_, mut g) in groups {
        if g.len() == 1 {
            continue;
        }
        g.sort_unstable(); // Java HashSet 任意序；取确定序
        count += attempt_zip_merge_seq_vertices(&g, graph, zip, ctx, true);
    }
    count
}

/// `zip_down`（L2621）：v 的多个**后继**若同 orig id → 合并（depth 取 max）。
pub fn zip_down(graph: &mut DiGraph, zip: &mut ZipState, ctx: &mut BflyContext, v: i32) -> usize {
    let child_list: Vec<i32> = graph.get_successors(v).to_vec();
    if child_list.len() <= 1 {
        return 0;
    }
    let mut groups: Vec<(i32, Vec<i32>)> = Vec::new();
    for &c in &child_list {
        if zip.is_replacement.contains(&c) {
            return 0;
        }
        let orig = graph.get_vertex(c).unwrap().orig_butterfly_id;
        match groups.iter_mut().find(|(o, _)| *o == orig) {
            Some((_, g)) => g.push(c),
            None => groups.push((orig, vec![c])),
        }
    }
    let mut count = 0;
    for (_, mut g) in groups {
        if g.len() == 1 {
            continue;
        }
        g.sort_unstable();
        count += attempt_zip_merge_seq_vertices(&g, graph, zip, ctx, false);
    }
    count
}

/// `zipper_collapse_DAG_zip_up`（L2490）：拓扑序**逆序**（bottom-up）扫。
pub fn zipper_collapse_dag_zip_up(
    graph: &mut DiGraph,
    zip: &mut ZipState,
    ctx: &mut BflyContext,
) -> usize {
    let order = topo_sort_seq_vertices_dag(graph).expect("seqvertex graph 非 DAG（zip_up 前）");
    assign_depths(graph, &order);
    let mut total = 0;
    for &v in order.iter().rev() {
        if zip.is_replacement.contains(&v) || !graph.contains_vertex(v) {
            continue;
        }
        total += zip_up(graph, zip, ctx, v);
    }
    total
}

/// `zipper_collapse_DAG_zip_down`（L2529）：拓扑序正序（top-down）扫。
pub fn zipper_collapse_dag_zip_down(
    graph: &mut DiGraph,
    zip: &mut ZipState,
    ctx: &mut BflyContext,
) -> usize {
    let order = topo_sort_seq_vertices_dag(graph).expect("seqvertex graph 非 DAG（zip_down 前）");
    assign_depths(graph, &order);
    let mut total = 0;
    for &v in order.iter() {
        if zip.is_replacement.contains(&v) || !graph.contains_vertex(v) {
            continue;
        }
        total += zip_down(graph, zip, ctx, v);
    }
    total
}

/// `destroy_unzipped_duplicates_above`（L2393）：无前驱的 v 若与某个**有前驱**的
/// O 同 orig id 且共享孩子 c → 删 v（其 id 并入 O 的 compressed 列表）。
pub fn destroy_unzipped_duplicates_above(graph: &mut DiGraph, zip: &mut ZipState) {
    let start_vertices: Vec<i32> = graph
        .vertex_ids()
        .iter()
        .copied()
        .filter(|&v| graph.in_degree(v) == 0)
        .collect();
    if start_vertices.is_empty() {
        return;
    }
    let mut edges_to_delete: Vec<(i32, i32)> = Vec::new();
    let mut vertices_to_delete: Vec<i32> = Vec::new();
    for &v in &start_vertices {
        let v_orig = graph.get_vertex(v).unwrap().orig_butterfly_id;
        let mut target_merge_vertex: Option<i32> = None;
        'outer: for &c in graph.get_successors(v) {
            for &o in graph.get_predecessors(c) {
                if o != v
                    && graph.get_vertex(o).unwrap().orig_butterfly_id == v_orig
                    && graph.in_degree(o) > 0
                {
                    target_merge_vertex = Some(o);
                    edges_to_delete.push((v, c));
                    break 'outer;
                }
            }
        }
        if let Some(o) = target_merge_vertex {
            vertices_to_delete.push(v);
            let mut merged = vec![v, o]; // Java 先 add v 再 add O 自身
            merged.extend_from_slice(zip.compressed(v));
            zip.add_compressed(o, merged);
        }
    }
    for e in edges_to_delete {
        graph.remove_edge(e.0, e.1);
    }
    for v in vertices_to_delete {
        graph.remove_vertex(v);
    }
}

// ---------------------------------------------------------------------------
// update_PairPaths_using_overlapDAG_refined_paths
// ---------------------------------------------------------------------------

/// `get_all_possible_updated_path_mappings`（L2116）：把 read 路径按 orig id 对齐到
/// 每条 revised path 上，收集所有可行映射（去重，保序）。
pub fn get_all_possible_updated_path_mappings(
    p1: &[i32],
    orig_of: &dyn Fn(i32) -> i32,
    revised_paths: &[PathWithOrig],
) -> Vec<Vec<i32>> {
    let needs_updating = PathWithOrig::from_path(p1, orig_of);
    let mut all: Vec<Vec<i32>> = Vec::new();
    for pwo in revised_paths {
        if let Some(updated) = needs_updating.align_path_by_orig_id(pwo) {
            if !all.contains(&updated.vertex_id_list) {
                all.push(updated.vertex_id_list);
            }
        }
    }
    assert!(
        !all.is_empty(),
        "Unable to remap read: {p1:?} given revised paths"
    );
    all
}

/// `construct_combinedReadhHash_from_PairPath_list`（L3786）。
fn construct_combined_read_hash_from_pair_path_list(
    pairpath_hmap: &FxHashMap<PairPath, i64>,
) -> FxHashMap<i32, FxHashMap<PairPath, i64>> {
    let mut combined: FxHashMap<i32, FxHashMap<PairPath, i64>> = FxHashMap::default();
    for (pp, &support) in pairpath_hmap {
        combined
            .entry(pp.get_first_id())
            .or_default()
            .insert(pp.clone(), support);
    }
    combined
}

/// `update_PairPaths_using_overlapDAG_refined_paths`（L1935）：按 POG 修订路径
/// 重映射每条 PairPath；双端都唯一映射 → 保留配对（原 support），否则拆成单端
/// support=1（Java `containsKey(List)` 恒 false 的 quirk 一并镜像：重复覆盖）。
pub fn update_pair_paths_using_overlap_dag_refined_paths(
    revised_paths: &[PathWithOrig],
    orig_paths: &[Vec<i32>],
    pair_path_to_read_support: &[(PairPath, i64)],
    orig_of: &dyn Fn(i32) -> i32,
) -> FxHashMap<i32, FxHashMap<PairPath, i64>> {
    let mut old_to_new_path: FxHashMap<Vec<i32>, Vec<i32>> = FxHashMap::default();
    for (orig, pwo) in orig_paths.iter().zip(revised_paths.iter()) {
        old_to_new_path.insert(orig.clone(), pwo.vertex_id_list.clone());
    }

    let mut updated_pair_paths: FxHashMap<PairPath, i64> = FxHashMap::default();
    for (pp, support) in pair_path_to_read_support {
        let p1 = &pp.path1;
        let p1_list: Vec<Vec<i32>> = match old_to_new_path.get(p1) {
            Some(m) => vec![m.clone()],
            None => get_all_possible_updated_path_mappings(p1, orig_of, revised_paths),
        };

        if pp.has_second_path() {
            let p2 = &pp.path2;
            let p2_list: Vec<Vec<i32>> = match old_to_new_path.get(p2) {
                Some(m) => vec![m.clone()],
                None => get_all_possible_updated_path_mappings(p2, orig_of, revised_paths),
            };
            if p1_list.len() == 1 && p2_list.len() == 1 {
                let new_pp = PairPath::with_pair(p1_list[0].clone(), p2_list[0].clone());
                updated_pair_paths.insert(new_pp, *support);
            } else {
                // Java quirk：containsKey(List) 恒 false → 每条都插入（覆盖）
                for p in p1_list {
                    updated_pair_paths.insert(PairPath::new(p), 1);
                }
                for p in p2_list {
                    updated_pair_paths.insert(PairPath::new(p), 1);
                }
            }
        } else {
            for p in p1_list {
                updated_pair_paths.insert(PairPath::new(p), *support);
            }
        }
    }
    construct_combined_read_hash_from_pair_path_list(&updated_pair_paths)
}

// ---------------------------------------------------------------------------
// 编排：create_DAG_from_OverlapLayout（L1617）
// ---------------------------------------------------------------------------

/// `create_DAG_from_OverlapLayout` 的产物（图 + 检查点 + 重映射后的 read hash）。
#[derive(Debug)]
pub struct OverlapLayoutResult {
    pub pog: PathOverlapGraph,
    pub path_matches: FxHashMap<(usize, usize), PathOverlap>,
    pub contained_path_to_containers: FxHashMap<Vec<i32>, Vec<Vec<i32>>>,
    /// 非 contained 路径（降序），与 pog.nodes 一致。
    pub noncontained_paths: Vec<Vec<i32>>,
    pub seqvertex_graph: DiGraph,
    pub zip_state: ZipState,
    /// POG 节点 idx → zipping 后的修订路径。
    pub revised_paths: Vec<PathWithOrig>,
    pub combined_read_hash: FxHashMap<i32, FxHashMap<PairPath, i64>>,
    /// `_POG.dot`。
    pub pog_dot: String,
    /// `_POG.PE_links_added.dot`（Java 已关闭 PE → 与 pog_dot 恒相同）。
    pub pe_links_dot: String,
    /// `_POG.cyclesRemoved.r{N}.dot`（含最后一轮无环轮）。
    pub cycle_round_dots: Vec<String>,
    /// `_before_zippingUpSeqVertexGraph.dot`。
    pub before_zipping_dot: String,
    /// `_before_zippingUpSeqVertexGraph.TopoSort.dot`。
    pub before_zipping_toposort_dot: String,
    /// `_zip_round_{N}_{zip_up|zip_down}.dot`（按轮次序）。
    pub zip_round_dots: Vec<String>,
    /// (round, "zip_up"/"zip_down", 本轮合并数)。
    pub zip_round_merges: Vec<(usize, &'static str, usize)>,
}

/// `create_DAG_from_OverlapLayout`（L1617）全链：
/// populate → 降序 → remove_containments → dispersed repeats → POG 构建 →
/// PE links（关闭）→ 破环（rN）→ convert（展开 + DFS 跨路径边）→ topo →
/// zipping 各轮 → destroy_unzipped → PairPath 重映射。
///
/// 输入 `orig_graph` = 剪枝/压缩后的 de Bruijn 图；`suff` = `getSuffStats_wPairs`
/// 的结果；`ctx` 提供 `getNextID`（新 SeqVertex id 延续全局计数器）。
pub fn create_dag_from_overlap_layout(
    orig_graph: &DiGraph,
    ctx: &mut BflyContext,
    suff: &SuffStats,
) -> OverlapLayoutResult {
    // populate_pairpaths_and_readsupport
    let pair_path_to_read_support = populate_pairpaths_and_readsupport(&suff.combined_read_hash);
    let support_map: FxHashMap<PairPath, i64> = pair_path_to_read_support.iter().cloned().collect();

    // 收集路径并按长度降序（Java 稳定排序 + reverse → 长度相同保持原序）
    let mut paths: Vec<Vec<i32>> = Vec::new();
    for (pp, _) in &pair_path_to_read_support {
        paths.push(pp.path1.clone());
        if pp.has_second_path() {
            paths.push(pp.path2.clone());
        }
    }
    paths.sort_by_key(|p| p.len());
    paths.reverse();

    // remove_containments
    let mut contained_path_to_containers: FxHashMap<Vec<i32>, Vec<Vec<i32>>> = FxHashMap::default();
    let noncontained_paths = remove_containments(&paths, &mut contained_path_to_containers);

    // find_dispersed_repeat_nodes
    let dispersed_repeat_nodes = find_dispersed_repeat_nodes(&noncontained_paths);

    // construct_path_overlap_graph
    let mut path_matches: FxHashMap<(usize, usize), PathOverlap> = FxHashMap::default();
    let mut pog = construct_path_overlap_graph(
        &noncontained_paths,
        &mut path_matches,
        &dispersed_repeat_nodes,
    );
    let pog_dot = pog.write_dot();

    // PE links（Java 已关闭 → 恒空；转换前移除 pair-link 边也是 no-op）
    let _pair_links =
        add_pair_paths_to_overlap_graph(&mut pog, &support_map, &contained_path_to_containers);
    let pe_links_dot = pog.write_dot();

    // 破环（迭代，直到一轮无环可破；最后一轮仍写 rN 检查点）
    let mut cycle_round_dots: Vec<String> = Vec::new();
    loop {
        let breaking = break_cycles_in_path_overlap_graph(&mut pog);
        cycle_round_dots.push(pog.write_dot());
        if !breaking {
            break;
        }
    }

    // ---- convert_path_DAG_to_SeqVertex_DAG（L2146，内联以持有图本体）----

    // 每条 POG 节点展开成新 SeqVertex 链 + 链内边（weight=1）
    let mut seqgraph = DiGraph::new();
    let mut zip = ZipState::default();
    let mut node_to_vertex_ids: Vec<Vec<i32>> = Vec::new();
    let mut revised_paths: Vec<PathWithOrig> = Vec::new();

    for node in &pog.nodes {
        let mut vertex_listing: Vec<i32> = Vec::new();
        for &node_id in &node.vertices {
            let orig_v = orig_graph.get_vertex(node_id).unwrap();
            let new_id = ctx.get_next_id();
            let mut new_v = SeqVertex::new(new_id, orig_v.name.clone());
            new_v.weights = orig_v.weights.clone();
            new_v.orig_butterfly_id = orig_v.orig_butterfly_id;
            new_v.depth = orig_v.depth + orig_v.name.len() as i32;
            vertex_listing.push(new_id);
            seqgraph.add_vertex(new_v);
        }
        for w in vertex_listing.windows(2) {
            seqgraph.add_edge(w[0], w[1], SimpleEdge::new(1.0, w[0], w[1]));
        }
        node_to_vertex_ids.push(vertex_listing.clone());
        revised_paths.push(PathWithOrig {
            vertex_id_list: vertex_listing,
            orig_vertex_id_list: node.vertices.clone(),
        });
    }

    // DFS_add_path_to_graph（发布版 jar 的实际实现——由 Butterfly.jar 字节码
    // 反编译还原；源码树此函数是旧版，行为不同）：
    // * 驱动（convert_path_DAG_to_SeqVertex_DAG 尾部）：按 POG 拓扑序遍历，
    //   只对**入度 0** 的节点发起 DFS。
    // * 每帧：后继按 PN id **字符串序**（`"PN10" < "PN9"`，Comparator 为
    //   String.compareTo）；
    //   `if succs_seen.contains(succ) && parents_seen.contains(p) { return }`
    //   （整个帧返回，丢弃剩余后继）；
    //   跨边（idxA = idx_start_A+match_length-1，idxB = idx_start_B+同式）：
    //     idxA>=1 → A[idxA-1] → B[idxB]
    //     idxB>=1 → B[idxB-1] → A[idxA]
    //   然后 parents_seen += p、succs_seen += succ，递归 succ。
    //（典型 ml=1 场景退化为单边 A[len-2]→B[0]。）
    #[allow(clippy::too_many_arguments)]
    fn dfs_add_path_to_graph(
        p: usize,
        seqgraph: &mut DiGraph,
        pog: &PathOverlapGraph,
        path_matches: &FxHashMap<(usize, usize), PathOverlap>,
        node_to_vertex_ids: &[Vec<i32>],
        visited: &mut FxHashSet<usize>,
        succs_seen: &mut FxHashSet<usize>,
        parents_seen: &mut FxHashSet<usize>,
    ) {
        if visited.contains(&p) {
            return;
        }
        visited.insert(p);
        let mut succs = pog.successors(p);
        succs.sort_by_key(|&s| format!("PN{}", pog.nodes[s].pn_id));
        for s in succs {
            if succs_seen.contains(&s) && parents_seen.contains(&p) {
                return;
            }
            let po = path_matches[&(p, s)];
            let (a, b) = (&node_to_vertex_ids[p], &node_to_vertex_ids[s]);
            let idx_a = po.idx_start_a + po.match_length.saturating_sub(1);
            let idx_b = po.idx_start_b + po.match_length.saturating_sub(1);
            if idx_a >= 1 {
                let (u, v) = (a[idx_a - 1], b[idx_b]);
                seqgraph.add_edge(u, v, SimpleEdge::new(1.0, u, v));
            }
            if idx_b >= 1 {
                let (u, v) = (b[idx_b - 1], a[idx_a]);
                seqgraph.add_edge(u, v, SimpleEdge::new(1.0, u, v));
            }
            parents_seen.insert(p);
            succs_seen.insert(s);
            dfs_add_path_to_graph(
                s,
                seqgraph,
                pog,
                path_matches,
                node_to_vertex_ids,
                visited,
                succs_seen,
                parents_seen,
            );
        }
    }

    // 驱动拓扑序（Kahn，节点下标序——对拍以结构比较，序仅影响平局）
    let mut in_deg: Vec<usize> = (0..pog.nodes.len()).map(|i| pog.in_degree(i)).collect();
    let mut topo_order: Vec<usize> = (0..pog.nodes.len()).filter(|&i| in_deg[i] == 0).collect();
    let mut qi = 0;
    while qi < topo_order.len() {
        let u = topo_order[qi];
        qi += 1;
        for v in pog.successors(u) {
            in_deg[v] -= 1;
            if in_deg[v] == 0 {
                topo_order.push(v);
            }
        }
    }

    let mut visited: FxHashSet<usize> = FxHashSet::default();
    let mut succs_seen: FxHashSet<usize> = FxHashSet::default();
    let mut parents_seen: FxHashSet<usize> = FxHashSet::default();
    for &p in &topo_order {
        if pog.in_degree(p) != 0 {
            continue;
        }
        dfs_add_path_to_graph(
            p,
            &mut seqgraph,
            &pog,
            &path_matches,
            &node_to_vertex_ids,
            &mut visited,
            &mut succs_seen,
            &mut parents_seen,
        );
    }

    let before_zipping_dot = write_seqvertex_dot(&seqgraph);
    let topo = topo_sort_seq_vertices_dag(&seqgraph).expect("seqvertex_graph 含环");
    assign_depths(&mut seqgraph, &topo);
    let before_zipping_toposort_dot = write_seqvertex_dot(&seqgraph);

    // ---- zipping 轮次（L2220-2330：外层 sum>0；内层先 zip_up 不动点再 zip_down）----
    let mut zip_round_dots: Vec<String> = Vec::new();
    let mut zip_round_merges: Vec<(usize, &'static str, usize)> = Vec::new();
    let mut zip_round = 0usize;
    let mut sum_merged = 1;
    while sum_merged > 0 {
        sum_merged = 0;
        let mut count = 1;
        while count > 0 {
            zip_round += 1;
            init_replacement_vertices(&mut zip);
            count = zipper_collapse_dag_zip_up(&mut seqgraph, &mut zip, ctx);
            sum_merged += count;
            zip_round_merges.push((zip_round, "zip_up", count));
            zip_round_dots.push(write_seqvertex_dot(&seqgraph));
        }
        let mut count = 1;
        while count > 0 {
            zip_round += 1;
            init_replacement_vertices(&mut zip);
            count = zipper_collapse_dag_zip_down(&mut seqgraph, &mut zip, ctx);
            sum_merged += count;
            zip_round_merges.push((zip_round, "zip_down", count));
            zip_round_dots.push(write_seqvertex_dot(&seqgraph));
        }
    }

    destroy_unzipped_duplicates_above(&mut seqgraph, &mut zip);
    let topo = topo_sort_seq_vertices_dag(&seqgraph).expect("seqvertex_graph 含环（zipping 后）");
    assign_depths(&mut seqgraph, &topo);

    // 旧→新 id 映射（compressed 列表 → 保留节点），更新修订路径
    let mut old_to_new: FxHashMap<i32, i32> = FxHashMap::default();
    for &vid in &topo {
        let compressed = zip.compressed(vid).to_vec();
        if !compressed.is_empty() {
            for old in compressed {
                old_to_new.insert(old, vid);
            }
        } else {
            old_to_new.insert(vid, vid);
        }
    }
    for pwo in revised_paths.iter_mut() {
        pwo.vertex_id_list = pwo
            .vertex_id_list
            .iter()
            .map(|&id| *old_to_new.get(&id).expect("旧→新 id 映射缺失"))
            .collect();
    }

    // update_PairPaths_using_overlapDAG_refined_paths
    let orig_of = |id: i32| -> i32 {
        orig_graph
            .get_vertex(id)
            .map(|v| v.orig_butterfly_id)
            .unwrap_or(id)
    };
    let combined_read_hash = update_pair_paths_using_overlap_dag_refined_paths(
        &revised_paths,
        &noncontained_paths,
        &pair_path_to_read_support,
        &orig_of,
    );

    OverlapLayoutResult {
        pog,
        path_matches,
        contained_path_to_containers,
        noncontained_paths,
        seqvertex_graph: seqgraph,
        zip_state: zip,
        revised_paths,
        combined_read_hash,
        pog_dot,
        pe_links_dot,
        cycle_round_dots,
        before_zipping_dot,
        before_zipping_toposort_dot,
        zip_round_dots,
        zip_round_merges,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hs(ids: &[i32]) -> FxHashSet<i32> {
        ids.iter().copied().collect()
    }

    // ---------------- Path.java 移植 ----------------

    #[test]
    fn extends_basic_suffix_overlap() {
        let po = path_b_extends_path_a_allow_repeats(&[3, 4, 5], &[1, 2, 3, 4, 5], &hs(&[]));
        assert_eq!(po.match_score, 3);
        assert_eq!(po.match_length, 3);
    }

    #[test]
    fn extends_repeat_nodes_match_but_not_scored() {
        // pathA=[7,1], pathB=[1,7]：1 是 repeat 节点 → 匹配但不计分（score 0）
        let po = path_b_extends_path_a_allow_repeats(&[1, 7], &[7, 1], &hs(&[1]));
        assert_eq!(po.match_score, 0);
        assert_eq!(po.match_length, 1);
    }

    #[test]
    fn extends_path_b_longer_than_path_a_scans_all_positions() {
        let po = path_b_extends_path_a_allow_repeats(&[2, 3, 9], &[2, 3], &hs(&[]));
        // start = pathA.len - pathB.len 为负 → pathA 从尾向前全扫（i=1 不匹配，i=0 匹配）
        assert_eq!(po.match_score, 2);
        assert_eq!(po.match_length, 2);
    }

    #[test]
    fn extends_best_start_wins() {
        let po = path_b_extends_path_a_allow_repeats(&[1, 2], &[1, 2, 1, 2], &hs(&[]));
        assert_eq!(po.match_score, 2);
        assert_eq!(po.match_length, 2);
        // 无匹配 → (0,0)
        let po = path_b_extends_path_a_allow_repeats(&[5], &[1, 2], &hs(&[]));
        assert_eq!(po.match_score, 0);
        assert_eq!(po.match_length, 0);
    }

    #[test]
    fn contains_contiguous_substring() {
        assert!(path_a_contains_path_b_allow_repeats(&[1, 2, 3, 4], &[2, 3]));
        assert!(path_a_contains_path_b_allow_repeats(&[1, 2, 3], &[1, 2, 3]));
        assert!(!path_a_contains_path_b_allow_repeats(&[1, 3, 2], &[1, 2]));
        assert!(path_a_contains_path_b_allow_repeats(&[5, 5, 5], &[5, 5]));
    }

    #[test]
    fn repeat_nodes_and_counts() {
        let m = get_repeat_nodes_and_counts(&[1, 2, 2, 3, 2]);
        assert_eq!(m.len(), 1);
        assert_eq!(m[&2], 3);
    }

    // ---------------- remove_containments / dispersed ----------------

    #[test]
    fn remove_containments_keeps_first_and_records_all_containers() {
        let paths = vec![vec![1, 2, 3, 4], vec![2, 3], vec![2, 3, 4], vec![9]];
        let mut containers = FxHashMap::default();
        let kept = remove_containments(&paths, &mut containers);
        assert_eq!(kept, vec![vec![1, 2, 3, 4], vec![9]]);
        assert_eq!(containers[&vec![2, 3]], vec![vec![1, 2, 3, 4]]);
        assert_eq!(containers[&vec![2, 3, 4]], vec![vec![1, 2, 3, 4]]);
    }

    #[test]
    fn dispersed_repeats_need_ten_paths() {
        let mut paths: Vec<Vec<i32>> = (0..9).map(|i| vec![i, 100]).collect();
        assert!(find_dispersed_repeat_nodes(&paths).is_empty());
        paths.push(vec![50, 100]); // 第 10 条含 100
        assert_eq!(find_dispersed_repeat_nodes(&paths), hs(&[100]));
        // 每条路径只计一次
        let paths2: Vec<Vec<i32>> = (0..3).map(|_| vec![100, 100, 100, 7]).collect();
        assert!(find_dispersed_repeat_nodes(&paths2).is_empty());
    }

    // ---------------- POG 构建 / 破环 ----------------

    fn make_pog(paths: &[Vec<i32>]) -> (PathOverlapGraph, FxHashMap<(usize, usize), PathOverlap>) {
        let mut pm = FxHashMap::default();
        let dispersed = find_dispersed_repeat_nodes(paths);
        let pog = construct_path_overlap_graph(paths, &mut pm, &dispersed);
        (pog, pm)
    }

    #[test]
    fn pog_edges_for_suffix_extensions() {
        let (pog, pm) = make_pog(&[vec![1, 2, 3], vec![2, 3, 4], vec![3, 4, 5]]);
        // [2,3,4] extends [1,2,3] → 边 0→1；[3,4,5] extends [2,3,4] → 1→2；
        // [3,4,5] 也 extends [1,2,3]（后缀 [3] 重叠）→ 0→2（store_best=false 全建边）
        assert_eq!(pog.edge_order, vec![(0, 1), (0, 2), (1, 2)]);
        assert_eq!(pm[&(0, 1)].match_score, 2);
        assert_eq!(pm[&(0, 1)].match_length, 2);
        assert_eq!(pm[&(0, 2)].match_score, 1);
    }

    #[test]
    fn break_cycles_removes_loop_edges() {
        // 三条路径两两成环：A=[1,2], B=[2,3], C=[3,1]
        let (mut pog, _) = make_pog(&[vec![1, 2], vec![2, 3], vec![3, 1]]);
        assert_eq!(pog.edge_order.len(), 3, "应构成 3 边环");
        assert!(break_cycles_in_path_overlap_graph(&mut pog));
        assert_eq!(pog.edge_order.len(), 2, "删 1 条边破环");
        // 再跑一轮：无环 → false
        assert!(!break_cycles_in_path_overlap_graph(&mut pog));
    }

    #[test]
    fn pog_dot_format_matches_java() {
        let (pog, _) = make_pog(&[vec![1, 2, 3], vec![2, 3, 4]]);
        let dot = pog.write_dot();
        assert_eq!(
            dot,
            "digraph G {\n\tPN1 [label=\"PN1\"]\n\tPN1->PN2\n\tPN2 [label=\"PN2\"]\n}\n"
        );
    }

    // ---------------- PathWithOrig / 重映射 ----------------

    #[test]
    fn align_path_by_orig_id_full_match() {
        let template = PathWithOrig {
            vertex_id_list: vec![10, 11, 12, 13],
            orig_vertex_id_list: vec![1, 2, 3, 4],
        };
        let read = PathWithOrig {
            vertex_id_list: vec![0, 0],
            orig_vertex_id_list: vec![2, 3],
        };
        let aligned = read.align_path_by_orig_id(&template).unwrap();
        assert_eq!(aligned.vertex_id_list, vec![11, 12]);
        // 起点不匹配 → None
        let bad = PathWithOrig {
            vertex_id_list: vec![0, 0],
            orig_vertex_id_list: vec![9, 3],
        };
        assert!(bad.align_path_by_orig_id(&template).is_none());
    }

    #[test]
    fn all_possible_mappings_dedupe() {
        let orig_of = |id: i32| id;
        let t1 = PathWithOrig {
            vertex_id_list: vec![10, 11],
            orig_vertex_id_list: vec![1, 2],
        };
        let t2 = PathWithOrig {
            vertex_id_list: vec![20, 21],
            orig_vertex_id_list: vec![1, 2],
        };
        let mappings = get_all_possible_updated_path_mappings(&[1, 2], &orig_of, &[t1, t2]);
        assert_eq!(mappings, vec![vec![10, 11], vec![20, 21]]);
    }

    // ---------------- 拓扑排序 / zipping ----------------

    fn toy_graph() -> DiGraph {
        let mut g = DiGraph::new();
        // 1 → {2,3}(同 orig 100) → 4
        for (id, orig, name) in [(1, 1, "A"), (2, 100, "B"), (3, 100, "B"), (4, 4, "D")] {
            let mut v = SeqVertex::new(id, name);
            v.orig_butterfly_id = orig;
            g.add_vertex(v);
        }
        for (a, b) in [(1, 2), (1, 3), (2, 4), (3, 4)] {
            g.add_edge(a, b, SimpleEdge::new(1.0, a, b));
        }
        g
    }

    #[test]
    fn topo_sort_assigns_depths() {
        let mut g = toy_graph();
        let order = topo_sort_seq_vertices_dag(&g).unwrap();
        assert_eq!(order.first(), Some(&1));
        assign_depths(&mut g, &order);
        assert_eq!(g.get_vertex(1).unwrap().node_depth, 0);
        assert_eq!(g.get_vertex(4).unwrap().node_depth, 3);
    }

    #[test]
    fn topo_sort_detects_cycle() {
        let mut g = DiGraph::new();
        g.add_vertex(SeqVertex::new(1, "A"));
        g.add_vertex(SeqVertex::new(2, "B"));
        g.add_edge(1, 2, SimpleEdge::new(1.0, 1, 2));
        g.add_edge(2, 1, SimpleEdge::new(1.0, 2, 1));
        assert!(topo_sort_seq_vertices_dag(&g).is_none());
    }

    #[test]
    fn zip_up_merges_same_orig_preds() {
        let (mut g, mut ctx, mut zip) = (toy_graph(), BflyContext::new(), ZipState::default());
        ctx.last_id = 100; // 新 id 不与既有节点冲突
        let order = topo_sort_seq_vertices_dag(&g).unwrap();
        assign_depths(&mut g, &order);
        let merged = zip_up(&mut g, &mut zip, &mut ctx, 4);
        assert_eq!(merged, 2);
        // 图：1 → R(101) → 4
        assert_eq!(g.vertex_count(), 3);
        assert_eq!(g.edge_count(), 2);
        assert_eq!(zip.compressed_public(101), &[2, 3]);
        assert_eq!(g.get_vertex(101).unwrap().orig_butterfly_id, 100);
    }

    #[test]
    fn zip_up_blocked_by_depth_ordering() {
        let mut g = toy_graph();
        // 手工设深度：parent(1)=0、targets 2/3 = 1、child(4)=0 → max(parent)=0 不小于
        // min(child)=0 → 拒绝合并
        for (id, d) in [(1, 0), (2, 1), (3, 1), (4, 0)] {
            g.get_vertex_mut(id).unwrap().node_depth = d;
        }
        let mut ctx = BflyContext::new();
        let mut zip = ZipState::default();
        assert_eq!(zip_up(&mut g, &mut zip, &mut ctx, 4), 0);
        assert_eq!(g.vertex_count(), 4, "约束不满足 → 不合并");
    }

    #[test]
    fn destroy_unzipped_duplicates_removes_parentless_dup() {
        let mut g = DiGraph::new();
        // v(1, orig 7, 无前驱) → c(4)；O(3, orig 7, 有前驱 2) → 4
        for (id, orig) in [(1, 7), (2, 8), (3, 7), (4, 9)] {
            let mut v = SeqVertex::new(id, "X");
            v.orig_butterfly_id = orig;
            g.add_vertex(v);
        }
        for (a, b) in [(1, 4), (3, 4), (2, 3)] {
            g.add_edge(a, b, SimpleEdge::new(1.0, a, b));
        }
        let mut zip = ZipState::default();
        destroy_unzipped_duplicates_above(&mut g, &mut zip);
        assert!(!g.contains_vertex(1), "无前驱重复节点 v 被删");
        assert_eq!(g.in_degree(4), 1);
        let mut merged = zip.compressed_public(3).to_vec();
        merged.sort_unstable();
        assert_eq!(merged, vec![1, 3]);
    }

    // ---------------- populate 覆盖 quirk ----------------

    #[test]
    fn populate_support_overwrite_quirk() {
        let mut crh: FxHashMap<i32, FxHashMap<PairPath, i64>> = FxHashMap::default();
        // 两条不同 pp trimSinkNodes 后相同 → support 覆盖（Java put 语义）
        crh.insert(
            1,
            [
                (PairPath::with_pair(vec![-1, 1, 2], vec![3]), 5),
                (PairPath::with_pair(vec![1, 2], vec![-1, 3]), 9),
            ]
            .into_iter()
            .collect(),
        );
        let out = populate_pairpaths_and_readsupport(&crh);
        assert_eq!(out.len(), 1, "trimSinkNodes 后两条 pp 相同 → 合一");
        assert_eq!(out[0].0, PairPath::with_pair(vec![1, 2], vec![3]));
        assert_eq!(out[0].1, 9, "后者覆盖前者（Java HashMap.put quirk）");
    }
}
