//! `getAllProbablePaths` 主算法（TransAssembly_allProbPaths.java L9601-10084）与
//! 其前后的编排（reorganizeReadPairings / 分量 / triplet 提取 / addSandT /
//! remove_identical_subseqs / remove_short_seqs）及 allProbPaths.fasta 雏形输出。
//!
//! 默认参数对应 `Butterfly.jar -L 200 -F 10000 -R 2`：
//! MIN_OUTPUT_SEQ=200、MAX_PAIR_DISTANCE=10000 →
//! PATH_REINFORCEMENT_DISTANCE = 25% * 10000 = 2500、MIN_READ_SUPPORT_THR=2。
//!
//! 默认开关（Java 静态默认，c0 主线）：ALL_POSSIBLE_PATHS=false、
//! LENIENT_PATH_CHECKING=false、USE_TRIPLETS=false、ORIGINAL_PATH_EXTENSIONS=false、
//! FRACTURE_UNRESOLVED_XSTRUCTURE=false、MISO_OUTPUT=true（T9：path=[id:j-k] 输出格式）。

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use rustc_hash::{FxHashMap, FxHashSet};

use crate::dfs::run_dfs2;
use crate::graph::{DiGraph, SeqVertex, SimpleEdge, T_VERTEX_ID, VERTEX_ROOT_ID};
use crate::pair_paths::PairPath;

/// 路径搜索参数（命令行可调项；Default = c0 黄金运行配置）。
#[derive(Debug, Clone)]
pub struct PathSearchParams {
    /// `-L`：输出序列最小长度（remove_short_seqs 用 `>=`）。
    pub min_output_seq: usize,
    /// `-R`：MIN_READ_SUPPORT_THR。
    pub min_read_support_thr: i64,
    /// PATH_REINFORCEMENT_DISTANCE = 25% * `-F`。
    pub path_reinforcement_distance: usize,
    /// `--max_number_of_paths_per_node_extend`（默认 25）。
    pub max_paths_per_node_extend: usize,
    /// `--max_number_of_paths_per_node_init`（默认 100）。
    pub max_paths_per_node_init: usize,
    /// KMER_SIZE（S/T PairPath、getPathSeq 的 k-1 截断）。
    pub kmer_size: usize,
}

impl Default for PathSearchParams {
    fn default() -> Self {
        Self {
            min_output_seq: 200,
            min_read_support_thr: 2,
            path_reinforcement_distance: 2500,
            max_paths_per_node_extend: 25,
            max_paths_per_node_init: 100,
            kmer_size: 25,
        }
    }
}

/// `getAllProbablePaths` 的产物：最终路径列表（按捕获序；Java 为 HashMap 无序）。
#[derive(Debug, Clone, Default)]
pub struct PathSearchResult {
    pub final_paths: Vec<Vec<i32>>,
}

pub type ReadHash = FxHashMap<i32, FxHashMap<PairPath, i64>>;

// ---------------------------------------------------------------------------
// Dijkstra / 可达性（edu.uci.ics.jung 补丁 + SeqVertex.isAncestral）
// ---------------------------------------------------------------------------

/// `SeqVertex.isAncestral(v1, v2)`（L618）：a→b 可达返回 1，b→a 返回 -1，都不可达 0。
/// Dijkstra 距离只判 null/非 null → 纯可达性（BFS）。
pub fn is_ancestral(graph: &DiGraph, a: i32, b: i32) -> i32 {
    if reachable(graph, a, b) {
        return 1;
    }
    if reachable(graph, b, a) {
        return -1;
    }
    0
}

/// 正向可达性（BFS）。
pub fn reachable(graph: &DiGraph, src: i32, dst: i32) -> bool {
    if !graph.contains_vertex(src) || !graph.contains_vertex(dst) {
        return false;
    }
    if src == dst {
        return true;
    }
    let mut seen: FxHashSet<i32> = FxHashSet::default();
    let mut queue = std::collections::VecDeque::from([src]);
    seen.insert(src);
    while let Some(v) = queue.pop_front() {
        for &u in graph.get_successors(v) {
            if u == dst {
                return true;
            }
            if seen.insert(u) {
                queue.push_back(u);
            }
        }
    }
    false
}

/// `DijkstraDistanceWoVer.getDistanceWoVer(source, target, verToExclude)`
/// （edu 补丁 125 行）：最短路 Dijkstra，扩展时跳过 `ver_to_exclude` 顶点
///（含出边端点检查 `!w.equals(verToExclude)`）。不可达返回 None。
/// 边权用 `SimpleEdge.weight`（非负；f64 用 to_bits 进堆——非负浮点下单调）。
pub fn get_distance_wo_ver(
    graph: &DiGraph,
    source: i32,
    target: i32,
    ver_to_exclude: i32,
) -> Option<f64> {
    if !graph.contains_vertex(source) || !graph.contains_vertex(target) {
        return None;
    }
    if source == target {
        return Some(0.0);
    }
    let mut dist: FxHashMap<i32, f64> = FxHashMap::default();
    let mut heap: BinaryHeap<Reverse<(u64, i32)>> = BinaryHeap::new();
    dist.insert(source, 0.0);
    heap.push(Reverse((0.0f64.to_bits(), source)));
    while let Some(Reverse((dbits, v))) = heap.pop() {
        let d = f64::from_bits(dbits);
        if v == target {
            return Some(d);
        }
        if d > *dist.get(&v).unwrap_or(&f64::INFINITY) {
            continue; // 过期记录
        }
        if v == ver_to_exclude {
            continue;
        }
        for &w in graph.get_successors(v) {
            if w == ver_to_exclude {
                continue;
            }
            let weight = graph.find_edge(v, w).map(|e| e.weight).unwrap_or(1.0);
            let nd = d + weight;
            if nd < *dist.get(&w).unwrap_or(&f64::INFINITY) {
                dist.insert(w, nd);
                heap.push(Reverse((nd.to_bits(), w)));
            }
        }
    }
    None
}

// ---------------------------------------------------------------------------
// 序列 / 长度（getPathSeq / getSeqPathLength）
// ---------------------------------------------------------------------------

/// `getPathSeq`（L9258）：首节点取全名，其余去前 K-1 碱基；id<0（S/T）跳过。
pub fn get_path_seq(graph: &DiGraph, path: &[i32], kmer_size: usize) -> String {
    let mut seq = String::new();
    let mut first_node = true;
    for &node_id in path {
        if node_id >= 0 {
            let name = &graph.get_vertex(node_id).unwrap().name;
            if !first_node && name.len() >= kmer_size.saturating_sub(1) {
                seq.push_str(&name[kmer_size - 1..]);
            } else {
                seq.push_str(name);
            }
            first_node = false;
        }
    }
    seq
}

fn get_seq_path_length(graph: &DiGraph, path: &[i32], kmer_size: usize) -> usize {
    get_path_seq(graph, path, kmer_size).len()
}

/// `getNameKmerAdj()` 长度（pathHasEnoughReadSupport 的 look-back 计量）。
fn name_kmer_adj_len(graph: &DiGraph, id: i32, kmer_size: usize) -> usize {
    if id < 0 {
        return 0;
    }
    let v = graph.get_vertex(id).unwrap();
    if graph.in_degree(id) > 0 {
        v.name.len().saturating_sub(kmer_size - 1)
    } else {
        v.name.len()
    }
}

// ---------------------------------------------------------------------------
// 前置编排：reorganizeReadPairings（combinePaths）/ 权重重标定 / 分量
// ---------------------------------------------------------------------------

/// 驱动 L1000-1010：把 SeqVertex DAG 的边权按原始 de Bruijn 图的 kmer 对重标定
///（zip 产生的边权为 1；键 = `源节点全名_目标节点全名`）。
pub fn relabel_edge_weights_using_orig_kmers(graph: &mut DiGraph, orig_graph: &DiGraph) {
    let k = orig_graph
        .vertex_ids()
        .first()
        .map(|&id| orig_graph.get_vertex(id).unwrap().name.len())
        .unwrap_or(25);
    let mut orig_weights: FxHashMap<(String, String), f64> = FxHashMap::default();
    for &id in orig_graph.vertex_ids() {
        let v = orig_graph.get_vertex(id).unwrap();
        for &succ in orig_graph.get_successors(id) {
            let w = orig_graph.get_vertex(succ).unwrap();
            if let Some(e) = orig_graph.find_edge(id, succ) {
                orig_weights.insert((v.name.clone(), w.name.clone()), e.weight);
            }
        }
    }
    for u in graph.vertex_ids().to_vec() {
        for v in graph.get_successors(u).to_vec() {
            let (fk, tk) = {
                let a = graph.get_vertex(u).unwrap();
                let b = graph.get_vertex(v).unwrap();
                (
                    a.get_last_kmer(k).to_string(),
                    b.get_first_kmer(k).to_string(),
                )
            };
            if let Some(&w) = orig_weights.get(&(fk, tk)) {
                if let Some(e) = graph.find_edge_mut(u, v) {
                    e.weight = w;
                }
            }
        }
    }
}

/// `reorganizeReadPairings`（L4206）：配对 PairPath 在（zip 后的）SeqVertex 图上
/// 重新组合成单条路径；不可组合则拆成两条单端。
pub fn reorganize_read_pairings(graph: &DiGraph, combined_read_hash: &ReadHash) -> ReadHash {
    let mut new_hash: ReadHash = FxHashMap::default();
    for pp_map in combined_read_hash.values() {
        for (pp, &read_support) in pp_map {
            if pp.has_second_path() {
                if let Some(combined) = combine_paths(graph, &pp.path1, &pp.path2) {
                    store_pair_path_by_first_vertex(combined, &mut new_hash, read_support);
                } else {
                    store_pair_path_by_first_vertex(
                        PairPath::new(pp.path1.clone()),
                        &mut new_hash,
                        read_support,
                    );
                    store_pair_path_by_first_vertex(
                        PairPath::new(pp.path2.clone()),
                        &mut new_hash,
                        read_support,
                    );
                }
            } else {
                store_pair_path_by_first_vertex(pp.clone(), &mut new_hash, read_support);
            }
        }
    }
    new_hash
}

fn store_pair_path_by_first_vertex(pp: PairPath, hash: &mut ReadHash, read_support: i64) {
    let first = pp.get_first_id();
    hash.entry(first)
        .or_default()
        .entry(pp)
        .and_modify(|c| *c += read_support)
        .or_insert(read_support);
}

/// `combinePaths`（L9417）：两半 read 在 DAG 上的方向判定 + 前缀重叠拼接 +
/// 中间路径 imputation（唯一可连通后继逐跳外推）。空路径返回 None。
pub fn combine_paths(graph: &DiGraph, path1: &[i32], path2: &[i32]) -> Option<PairPath> {
    let first_v1 = path1[0];
    let last_v1 = path1[path1.len() - 1];
    let first_v2 = path2[0];
    let last_v2 = path2[path2.len() - 1];

    let mut p1: Vec<i32> = Vec::new();
    let mut p2: Vec<i32> = Vec::new();

    let contains_all = |a: &[i32], b: &[i32]| b.iter().all(|x| a.contains(x));

    if contains_all(path1, path2) {
        p1 = path1.to_vec();
    } else if contains_all(path2, path1) {
        p2 = path2.to_vec();
    }
    // path1 --> path2
    else if is_ancestral(graph, last_v1, first_v2) > 0 && last_v1 != first_v2 {
        p1 = path1.to_vec();
        p2 = path2.to_vec();
    }
    // path2 --> path1
    else if is_ancestral(graph, last_v2, first_v1) > 0 && last_v2 != first_v1 {
        p1 = path2.to_vec();
        p2 = path1.to_vec();
    } else if is_ancestral(graph, first_v2, first_v1) == 0
        && is_ancestral(graph, last_v2, last_v1) == 0
    {
        // 无一致方向：留空
    }
    // path1(部分) -> path2
    else if is_ancestral(graph, first_v1, first_v2) > 0 && path1.contains(&first_v2) {
        let i = path1.iter().position(|&x| x == first_v2).unwrap();
        p1 = path1[..i].to_vec();
        p1.extend_from_slice(path2);
    }
    // path2(部分) -> path1
    else if is_ancestral(graph, first_v2, first_v1) > 0 && path2.contains(&first_v1) {
        let i = path2.iter().position(|&x| x == first_v1).unwrap();
        p1 = path2[..i].to_vec();
        p1.extend_from_slice(path1);
    }

    if p1.is_empty() && !p2.is_empty() {
        std::mem::swap(&mut p1, &mut p2);
    }

    // imputation：从 p1 末尾向 p2 开头外推（唯一 isAncestral>0 的后继）
    if !p1.is_empty() && !p2.is_empty() {
        let l1 = p1[p1.len() - 1];
        let f2 = p2[0];
        if is_ancestral(graph, l1, f2) > 0 {
            let mut impute = true;
            let mut v = l1;
            let mut intervening: Vec<i32> = Vec::new();
            loop {
                let mut next: Option<i32> = None;
                let mut count_connectable = 0;
                for &succ in graph.get_successors(v) {
                    if is_ancestral(graph, succ, f2) > 0 {
                        count_connectable += 1;
                        next = Some(succ);
                    }
                }
                match next {
                    Some(n) if count_connectable == 1 => {
                        if n == f2 {
                            break;
                        }
                        intervening.push(n);
                        v = n;
                    }
                    _ => {
                        impute = false;
                        break;
                    }
                }
            }
            if impute {
                p1.extend(intervening);
                p1.append(&mut p2);
            }
        }
    }

    if p1.is_empty() {
        None
    } else {
        Some(PairPath::with_pair(p1, p2))
    }
}

/// `divideIntoComponents`（L13377，WeakComponentClusterer）：弱连通分量。
pub fn divide_into_components(graph: &DiGraph) -> Vec<FxHashSet<i32>> {
    let mut comps: Vec<FxHashSet<i32>> = Vec::new();
    let mut assigned: FxHashSet<i32> = FxHashSet::default();
    for start in graph.vertex_ids().to_vec() {
        if assigned.contains(&start) {
            continue;
        }
        let mut comp: FxHashSet<i32> = FxHashSet::default();
        let mut queue = std::collections::VecDeque::from([start]);
        comp.insert(start);
        while let Some(v) = queue.pop_front() {
            for &u in graph.get_successors(v) {
                if comp.insert(u) {
                    queue.push_back(u);
                }
            }
            for &u in graph.get_predecessors(v) {
                if comp.insert(u) {
                    queue.push_back(u);
                }
            }
        }
        assigned.extend(comp.iter().copied());
        comps.push(comp);
    }
    comps
}

/// `getComponentReads`（L15759）：与分量共享任一节点的 PairPath（保留起始节点键）。
pub fn get_component_reads(comp: &FxHashSet<i32>, combined_read_hash: &ReadHash) -> ReadHash {
    let mut res: ReadHash = FxHashMap::default();
    for (node_id, pp_map) in combined_read_hash {
        for (pp, &count) in pp_map {
            let any = pp
                .path1
                .iter()
                .chain(pp.path2.iter())
                .any(|id| comp.contains(id));
            if any {
                res.entry(*node_id).or_default().insert(pp.clone(), count);
            }
        }
    }
    res
}

/// `reduce_to_max_paths_per_node`（L5806）：起始节点 PairPath 数超限时按
/// read 支持降序保留前 max 条（Java 稳定排序 + subList 截断）。
pub fn reduce_to_max_paths_per_node(read_hash: &mut ReadHash, max: usize) {
    for pp_map in read_hash.values_mut() {
        if pp_map.len() > max {
            let mut list: Vec<(PairPath, i64)> =
                pp_map.iter().map(|(k, v)| (k.clone(), *v)).collect();
            // 稳定排序：等支持保持插入序（Java 保持 HashMap keySet 序）
            list.sort_by_key(|(_, supp)| std::cmp::Reverse(*supp));
            for (pp, _) in list.into_iter().skip(max) {
                pp_map.remove(&pp);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// triplet 提取族
// ---------------------------------------------------------------------------

pub type TripletMapper = FxHashMap<i32, Vec<Vec<i32>>>;

/// `extractTripletsFromReads`（L15206）。
pub fn extract_triplets_from_reads(read_hash: &ReadHash) -> TripletMapper {
    let mut mapper: TripletMapper = FxHashMap::default();
    for pp_map in read_hash.values() {
        for pp in pp_map.keys() {
            for read_path in [&pp.path1, &pp.path2] {
                if read_path.len() < 3 {
                    continue;
                }
                for i in 1..read_path.len() - 1 {
                    let triplet = vec![read_path[i - 1], read_path[i], read_path[i + 1]];
                    let list = mapper.entry(read_path[i]).or_default();
                    if !list.contains(&triplet) {
                        list.push(triplet);
                    }
                }
            }
        }
    }
    mapper
}

/// `getXstructuresResolvedByTriplets`（L8874）。
pub fn get_xstructures_resolved_by_triplets(
    graph: &DiGraph,
    comp: &FxHashSet<i32>,
    triplet_mapper: &TripletMapper,
) -> FxHashMap<i32, bool> {
    let mut res: FxHashMap<i32, bool> = FxHashMap::default();
    for &vid in comp {
        if graph.in_degree(vid) > 1 && graph.out_degree(vid) > 1 {
            res.insert(vid, triplet_mapper.contains_key(&vid));
        }
    }
    res
}

/// `Path.share_suffix_fully_contained`（L617）—— **镜像 Java Integer 引用比较
/// quirk**：`path_i.get(i) != path_j.get(j)` 是对象同一性比较，只有 [-128,127]
/// 缓存区内的数值才可能跨列表相等。节点 id 大于 127 时首对比较即失败 →
/// 该函数实际恒 false（复杂前缀 purge 从不发生）。
fn share_suffix_fully_contained_java(a: &[i32], b: &[i32]) -> bool {
    let (mut i, mut j) = (a.len() as isize - 1, b.len() as isize - 1);
    while i >= 0 && j >= 0 {
        let (x, y) = (a[i as usize], b[j as usize]);
        // 数值相等且落在 Integer 缓存区内 → 同一对象 → "=="
        if !(x == y && (-128..=127).contains(&x)) {
            return false;
        }
        i -= 1;
        j -= 1;
    }
    true
}

/// `extractComplexPathPrefixesFromReads`（L15270）：每条 read 路径的每个
/// 长度 >= 3 的前缀（尾节点为键），再做（quirk 化的）子前缀 purge。
pub fn extract_complex_path_prefixes_from_reads(read_hash: &ReadHash) -> TripletMapper {
    let mut mapper: TripletMapper = FxHashMap::default();
    for pp_map in read_hash.values() {
        for pp in pp_map.keys() {
            for read_path in [&pp.path1, &pp.path2] {
                if read_path.len() < 3 {
                    continue;
                }
                for i in (2..read_path.len()).rev() {
                    let node_id = read_path[i];
                    let prefix = read_path[..i + 1].to_vec();
                    let list = mapper.entry(node_id).or_default();
                    if !list.contains(&prefix) {
                        list.push(prefix);
                    }
                }
            }
        }
    }
    // purge：更长前缀包含的短前缀删除（Java 引用比较 quirk → 见上）
    for prefixes in mapper.values_mut() {
        let mut to_purge: Vec<Vec<i32>> = Vec::new();
        for (ia, prefix) in prefixes.iter().enumerate() {
            for (ib, prefix2) in prefixes.iter().enumerate() {
                if ia != ib
                    && prefix2.len() > prefix.len()
                    && share_suffix_fully_contained_java(prefix, prefix2)
                {
                    to_purge.push(prefix.clone());
                    break;
                }
            }
        }
        prefixes.retain(|p| !to_purge.contains(p));
    }
    mapper
}

/// `tripletSupported`（L15355）。
pub fn triplet_supported(triplet_list: &[Vec<i32>], triplet: &[i32]) -> bool {
    triplet_list.iter().any(|t| t == triplet)
}

// ---------------------------------------------------------------------------
// addSandT / removeAllEdgesOfSandT
// ---------------------------------------------------------------------------

/// `addSandT`（L13405）：S(-1) → 每个 0 入度节点、每个 0 出度节点 → T(-2)，
/// 并向 read hash 注入 [S,v] / [v,T] PairPath（support = MIN_READ_SUPPORT_THR）。
/// S/T 的 _depth 设为 -1 / MAX（BflyQueue 序）。
pub fn add_s_and_t(
    graph: &mut DiGraph,
    comp: &FxHashSet<i32>,
    read_hash: &mut ReadHash,
    min_read_support_thr: i64,
) {
    let mut root = SeqVertex::new(VERTEX_ROOT_ID, "S");
    root.depth = -1;
    let mut t = SeqVertex::new(T_VERTEX_ID, "E");
    t.depth = i32::MAX;
    graph.add_vertex(root);
    graph.add_vertex(t);

    // 先收集（增边会改变度数）
    let mut s_targets: Vec<i32> = Vec::new();
    let mut t_targets: Vec<i32> = Vec::new();
    for &vid in comp {
        if graph.in_degree(vid) == 0 && vid != VERTEX_ROOT_ID && vid != T_VERTEX_ID {
            s_targets.push(vid);
        }
        if graph.out_degree(vid) == 0 && vid != T_VERTEX_ID && vid != VERTEX_ROOT_ID {
            t_targets.push(vid);
        }
    }
    for vid in s_targets {
        let w = graph
            .get_vertex(vid)
            .map(|v| v.weights.first().copied().unwrap_or(-1.0))
            .unwrap_or(-1.0);
        let w = if w == -1.0 { 1.0 } else { w };
        graph.add_edge(VERTEX_ROOT_ID, vid, SimpleEdge::new(w, VERTEX_ROOT_ID, vid));
        let pp = PairPath::new(vec![VERTEX_ROOT_ID, vid]);
        read_hash
            .entry(VERTEX_ROOT_ID)
            .or_default()
            .insert(pp, min_read_support_thr);
    }
    for vid in t_targets {
        let w = graph
            .get_vertex(vid)
            .map(|v| v.weights.last().copied().unwrap_or(-1.0))
            .unwrap_or(-1.0);
        let w = if w == -1.0 { 1.0 } else { w };
        graph.add_edge(vid, T_VERTEX_ID, SimpleEdge::new(w, vid, T_VERTEX_ID));
        let pp = PairPath::new(vec![vid, T_VERTEX_ID]);
        read_hash
            .entry(vid)
            .or_default()
            .insert(pp, min_read_support_thr);
    }
}

/// `removeAllEdgesOfSandT`（L13492）。
pub fn remove_all_edges_of_s_and_t(graph: &mut DiGraph) {
    if graph.contains_vertex(VERTEX_ROOT_ID) {
        for v in graph.get_successors(VERTEX_ROOT_ID).to_vec() {
            graph.remove_edge(VERTEX_ROOT_ID, v);
        }
    }
    if graph.contains_vertex(T_VERTEX_ID) {
        for u in graph.get_predecessors(T_VERTEX_ID).to_vec() {
            graph.remove_edge(u, T_VERTEX_ID);
        }
    }
}

// ---------------------------------------------------------------------------
// getAllProbablePaths 主体
// ---------------------------------------------------------------------------

/// 优先队列元素：node depth（升序 poll），同深度按入队序 FIFO
///（Java PriorityQueue 对相等元素不保证次序；此处取确定性的 FIFO）。
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct QueueItem {
    depth: i32,
    seq: u64,
    id: i32,
}

/// `parents_all_visited`（L10088）：所有 id>0 的前驱都已访问。
fn parents_all_visited(graph: &DiGraph, v: i32, visited: &FxHashSet<i32>) -> bool {
    graph
        .get_predecessors(v)
        .iter()
        .all(|&p| p <= 0 || visited.contains(&p))
}

/// `updateReadsOfPath`（L11675）：把 readsOfPathUntilV 中与 pathWu 兼容的
/// reads 传递给 pathWu（pathMinusU 已"完整包含"的直接传播）。
fn update_reads_of_path(
    path_reads: &mut FxHashMap<Vec<i32>, FxHashMap<PairPath, i64>>,
    path_contained_reads: &mut FxHashMap<Vec<i32>, FxHashSet<PairPath>>,
    path_wu: &[i32],
    reads_until_v: &FxHashMap<PairPath, i64>,
) {
    let path_minus_u = &path_wu[..path_wu.len() - 1];
    path_reads.entry(path_wu.to_vec()).or_default();
    path_contained_reads.entry(path_wu.to_vec()).or_default();

    for (pp, &count) in reads_until_v {
        let wu_reads = path_reads.get_mut(path_wu).unwrap();
        if wu_reads.contains_key(pp) {
            continue;
        }
        // Java：PathContainedReads.get(pathMinusU).contains(pPath)
        //（种子路径场景下 [v] 键可能缺失 → Java NPE；实际只有
        // readsStartingAtV 为空时才播种，故此处空集语义等价）
        let contained_in_minus = path_contained_reads
            .get(path_minus_u)
            .map(|s| s.contains(pp))
            .unwrap_or(false);
        if contained_in_minus {
            path_contained_reads
                .get_mut(path_wu)
                .unwrap()
                .insert(pp.clone());
            wu_reads.insert(pp.clone(), count);
        } else if pp.is_compatible_path(path_wu) {
            if pp.is_compatible_and_contained_by_single_path(path_wu) {
                path_contained_reads
                    .get_mut(path_wu)
                    .unwrap()
                    .insert(pp.clone());
            }
            wu_reads.insert(pp.clone(), count);
        }
    }
}

/// `subPathHasEnoughReadSupport`（L11271，COMPATIBLE_PATH_EXTENSIONS 默认分支）。
fn sub_path_has_enough_read_support(
    graph: &DiGraph,
    full_path_wu: &[i32],
    reads_until_v: &FxHashMap<PairPath, i64>,
    sub_path: &[i32],
    min_read_support_thr: i64,
) -> bool {
    let last_subpath_id = sub_path[sub_path.len() - 1];
    let mut first_subpath_id = sub_path[0];
    if first_subpath_id < 0 {
        first_subpath_id = sub_path[1]; // 不含汇点
    }
    let mut number_reads_supporting: i64 = 0;
    for (pp, &count) in reads_until_v {
        let ok = pp.contains_id(last_subpath_id)
            && pp.is_compatible_path(full_path_wu)
            && pp.node_is_contained_or_possibly_in_gap(first_subpath_id, &|a, b| {
                is_ancestral(graph, a, b)
            });
        if ok {
            number_reads_supporting += count;
            if number_reads_supporting >= min_read_support_thr {
                break;
            }
        }
    }
    number_reads_supporting >= min_read_support_thr
}

/// `pathHasEnoughReadSupport`（L11208，默认分支：look-back PATH_REINFORCEMENT_DISTANCE）。
fn path_has_enough_read_support(
    graph: &DiGraph,
    reads_until_v: &FxHashMap<PairPath, i64>,
    path: &[i32],
    u: i32,
    params: &PathSearchParams,
) -> bool {
    let mut path_wu: Vec<i32> = path.to_vec();
    path_wu.push(u);

    let mut sub_path: Vec<i32> = vec![u];
    let look_back = params.path_reinforcement_distance;
    let mut len_so_far = name_kmer_adj_len(graph, u, params.kmer_size);
    let mut j = path.len() as isize - 1;
    while j >= 0 && len_so_far < look_back {
        let v_last = path[j as usize];
        sub_path.insert(0, v_last);
        len_so_far += name_kmer_adj_len(graph, v_last, params.kmer_size);
        j -= 1;
    }
    sub_path_has_enough_read_support(
        graph,
        &path_wu,
        reads_until_v,
        &sub_path,
        params.min_read_support_thr,
    )
}

/// `getSuppCalculation`（L11655）。
pub fn get_supp_calculation(reads: &FxHashMap<PairPath, i64>) -> i64 {
    reads.values().sum()
}

/// `getAllProbablePaths`（L9601）。graph 需已含 S/T（addSandT）；
/// `component_read_hash` 为该分量的 read hash（reduce 后）。
#[allow(clippy::too_many_arguments)]
pub fn get_all_probable_paths(
    graph: &DiGraph,
    comp: &FxHashSet<i32>,
    component_read_hash: &ReadHash,
    triplet_mapper: &TripletMapper,
    extended_triplet_mapper: &TripletMapper,
    _x_structures: &FxHashMap<i32, bool>,
    params: &PathSearchParams,
) -> PathSearchResult {
    let kmer = params.kmer_size;

    // 顶点 → 已构造路径；路径 → PairPath 计数；路径 → 完整包含的 reads
    let mut paths: FxHashMap<i32, Vec<Vec<i32>>> = FxHashMap::default();
    let mut path_reads: FxHashMap<Vec<i32>, FxHashMap<PairPath, i64>> = FxHashMap::default();
    let mut path_contained_reads: FxHashMap<Vec<i32>, FxHashSet<PairPath>> = FxHashMap::default();
    let mut extensions: FxHashMap<Vec<i32>, bool> = FxHashMap::default();

    // FinalPaths_all：Java HashMap（按 path 去重）
    let mut final_all: Vec<Vec<i32>> = Vec::new();
    let mut final_all_set: FxHashSet<Vec<i32>> = FxHashSet::default();

    let mut queue: BinaryHeap<Reverse<QueueItem>> = BinaryHeap::new();
    let mut in_queue: FxHashSet<i32> = FxHashSet::default();
    let mut seq_counter: u64 = 0;
    let mut node_visited: FxHashSet<i32> = FxHashSet::default();

    macro_rules! enqueue {
        ($id:expr) => {{
            queue.push(Reverse(QueueItem {
                depth: graph.get_vertex($id).map(|v| v.depth).unwrap_or(-1),
                seq: seq_counter,
                id: $id,
            }));
            seq_counter += 1;
            in_queue.insert($id);
        }};
    }

    enqueue!(VERTEX_ROOT_ID);
    paths.insert(VERTEX_ROOT_ID, vec![vec![VERTEX_ROOT_ID]]);

    while !queue.is_empty() {
        let mut item = queue.pop().unwrap().0;
        in_queue.remove(&item.id);

        // 延迟处理：id>0 且父未全访问 → 轮询其他节点，全部受阻则报错
        //（受阻项最后原样放回队列，Java addAll）
        let mut delayed: Vec<QueueItem> = Vec::new();
        let mut blocked = item.id > 0 && !parents_all_visited(graph, item.id, &node_visited);
        while blocked {
            delayed.push(item);
            let Some(Reverse(next)) = queue.pop() else {
                break;
            };
            in_queue.remove(&next.id);
            item = next;
            blocked = item.id > 0 && !parents_all_visited(graph, item.id, &node_visited);
        }
        assert!(
            !blocked,
            "queue ran out of nodes and current node has unvisited parents"
        );
        for d in delayed {
            queue.push(Reverse(QueueItem {
                depth: d.depth,
                seq: seq_counter,
                id: d.id,
            }));
            seq_counter += 1;
            in_queue.insert(d.id);
        }

        let v = item.id;
        if node_visited.contains(&v) {
            continue; // 已访问（loop 防护）
        }
        node_visited.insert(v);

        let reads_starting_at_v = component_read_hash.get(&v);

        // 初始化 v 处各路径的 reads / contained / extensions
        for path in paths.get(&v).cloned().unwrap_or_default() {
            path_reads.entry(path.clone()).or_default();
            path_contained_reads.entry(path.clone()).or_default();
            if let Some(reads) = reads_starting_at_v {
                for (pp, &c) in reads {
                    path_reads.get_mut(&path).unwrap().insert(pp.clone(), c);
                }
            }
            extensions.insert(path, false);
        }

        // 逐后继扩展
        for u in graph.get_successors(v).to_vec() {
            if !(comp.contains(&u) || u == T_VERTEX_ID) {
                continue; // 只在本分量内扩展（汇点除外）
            }

            let mut path_counter = 0usize;
            let mut v_extended_to_u = false;

            let mut paths_ending_at_v = paths.get(&v).cloned().unwrap_or_default();
            // Java：PathReadSupportComparator 升序稳定排序后整体 reverse
            // → 支持数降序、等支持路径相对次序反转
            paths_ending_at_v
                .sort_by_key(|p| path_reads.get(p).map(get_supp_calculation).unwrap_or(0));
            paths_ending_at_v.reverse();

            for path in paths_ending_at_v {
                // ---- triplet 锁定 ----
                let mut path_wvu_acceptable = true;
                let mut extended_triplet_path_compatible = false;

                if path.len() >= 3 {
                    let w = path[path.len() - 2]; // w-v-u 三元组
                    if triplet_mapper.get(&v).map(|l| l.len() > 1).unwrap_or(false) {
                        let triplet = [w, v, u];
                        let triplet_list = &triplet_mapper[&v];
                        if triplet_supported(triplet_list, &triplet) {
                            // 扩展三元组：pathWu 与某复杂前缀兼容
                            let mut path_wu: Vec<i32> = path.clone();
                            path_wu.push(u);
                            for prefix_path in extended_triplet_mapper.get(&u).into_iter().flatten()
                            {
                                let ppath = PairPath::new(prefix_path.clone());
                                if ppath.is_compatible_and_contained_by_single_path(&path_wu) {
                                    extended_triplet_path_compatible = true;
                                    break;
                                }
                            }
                        } else {
                            path_wvu_acceptable = false; // 锁定
                        }
                    }
                    // FRACTURE_UNRESOLVED_XSTRUCTURE = false → 无动作
                }

                let reads_until_v = path_reads.get(&path).cloned().unwrap_or_default();

                if path_wvu_acceptable
                    && (extended_triplet_path_compatible
                        || path_counter <= params.max_paths_per_node_extend)
                    && (path_has_enough_read_support(graph, &reads_until_v, &path, u, params)
                        || u < 0)
                {
                    path_counter += 1;

                    let mut path_wu = path.clone();
                    path_wu.push(u);
                    let list = paths.entry(u).or_default();
                    if !list.contains(&path_wu) {
                        list.push(path_wu.clone());
                    }
                    update_reads_of_path(
                        &mut path_reads,
                        &mut path_contained_reads,
                        &path_wu,
                        &reads_until_v,
                    );
                    extensions.insert(path.clone(), true);
                    v_extended_to_u = true;
                }
            }

            if !in_queue.contains(&u) {
                enqueue!(u);
            }

            // 裸边保底：v 的任何路径都没用到边 (v,u) → 播种 [v,u]
            if !v_extended_to_u {
                let vu_path = vec![v, u];
                paths.entry(u).or_default().push(vu_path.clone());
                // Java：PathReads[vuPath].putAll(readsStartingAtV) 后调
                // updateReadsOfPath —— 所有 reads 已在 map 中 → 实际只做
                // putAll（contained 集保持空）
                let vu_reads = path_reads.entry(vu_path).or_default();
                if let Some(reads) = reads_starting_at_v {
                    for (pp, &c) in reads {
                        vu_reads.insert(pp.clone(), c);
                    }
                }
            }
        }

        // v 处理完：未扩展且未到 T 的路径收进 FinalPaths（len > MIN_OUTPUT_SEQ）
        let mut remove_paths: Vec<Vec<i32>> = Vec::new();
        for path in paths.get(&v).cloned().unwrap_or_default() {
            let last = path[path.len() - 1];
            if last != T_VERTEX_ID && extensions.get(&path).map(|e| !e).unwrap_or(false) {
                if get_seq_path_length(graph, &path, kmer) > params.min_output_seq
                    && final_all_set.insert(path.clone())
                {
                    final_all.push(path.clone());
                }
                remove_paths.push(path);
            }
        }
        if let Some(list) = paths.get_mut(&v) {
            for path in &remove_paths {
                list.retain(|p| p != path);
                extensions.remove(path);
            }
        }
    }

    // 到 T 的路径
    for path in paths.get(&T_VERTEX_ID).cloned().unwrap_or_default() {
        if get_seq_path_length(graph, &path, kmer) > params.min_output_seq
            && final_all_set.insert(path.clone())
        {
            final_all.push(path);
        }
    }

    if final_all.len() > 1 {
        final_all = remove_identical_subseqs(graph, final_all, kmer);
    }

    PathSearchResult {
        final_paths: final_all,
    }
}

/// `remove_identical_subseqs`（L15506）：序列长度降序（稳定），互为包含
///（indexOf）的路径删除被包含者；被过滤的路径不能再作为过滤证据。
pub fn remove_identical_subseqs(
    graph: &DiGraph,
    final_paths: Vec<Vec<i32>>,
    kmer_size: usize,
) -> Vec<Vec<i32>> {
    let mut path_vec: Vec<(Vec<i32>, String)> = final_paths
        .into_iter()
        .map(|p| {
            let s = get_path_seq(graph, &p, kmer_size);
            (p, s)
        })
        .collect();
    // FinalPaths.compareTo：长度降序（稳定排序保持捕获序）
    path_vec.sort_by_key(|(_, seq)| std::cmp::Reverse(seq.len()));

    let mut filtered: FxHashSet<usize> = FxHashSet::default();
    for i in 0..path_vec.len().saturating_sub(1) {
        if filtered.contains(&i) {
            continue;
        }
        for j in (i + 1)..path_vec.len() {
            if filtered.contains(&j) {
                continue;
            }
            let seq_i = &path_vec[i].1;
            let seq_j = &path_vec[j].1;
            if seq_i.contains(seq_j.as_str()) {
                filtered.insert(j);
            } else if seq_j.contains(seq_i.as_str()) {
                filtered.insert(i);
            }
        }
    }
    path_vec
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !filtered.contains(i))
        .map(|(_, p)| p.0)
        .collect()
}

/// `remove_short_seqs`（L1529）：保留 seq.len() >= MIN_OUTPUT_SEQ。
pub fn remove_short_seqs(
    graph: &DiGraph,
    final_paths: Vec<Vec<i32>>,
    min_output_seq: usize,
    kmer_size: usize,
) -> Vec<Vec<i32>> {
    final_paths
        .into_iter()
        .filter(|p| get_path_seq(graph, p, kmer_size).len() >= min_output_seq)
        .collect()
}

// ---------------------------------------------------------------------------
// 组件级编排（驱动 L960-1209）+ fasta 输出
// ---------------------------------------------------------------------------

/// `run_butterfly`（驱动 L940-1420）的完整结果。
#[derive(Debug, Default)]
pub struct ButterflyResult {
    /// 处理后的图（S/T 边已移除）
    pub graph: DiGraph,
    /// 最终路径（orig id，Java HashMap 迭代序——决定 isoform 编号）
    pub final_paths: Vec<Vec<i32>>,
    /// 与 final_paths 对齐的基因编号（1 起）
    pub gene_ids: Vec<usize>,
    /// final path → 兼容包含的 reads（orig id）
    pub contained_reads: FxHashMap<Vec<i32>, FxHashMap<PairPath, i64>>,
}

/// 全链（T9 完成）：My_DFS(runDFS2) → reorganizeReadPairings → 权重重标定 →
/// 分量循环（getComponentReads → reduce_to_max(100) → triplet/extendedTriplet
/// → addSandT → getAllProbablePaths → remove_short_seqs →
/// assignCompatibleReadsToPaths（删未分配 read 的路径）→ convert_to_orig_ids
/// → reduce_cdhit_like）→ 全集 reduce_cdhit_like → group_paths_into_genes →
/// EM stub（NO_EM_REDUCE）→ removeAllEdgesOfSandT。
pub fn run_butterfly_all_prob_paths(
    seqvertex_graph: &DiGraph,
    orig_kmer_graph: &DiGraph,
    orig_bfly_graph: &DiGraph,
    combined_read_hash: &ReadHash,
    params: &PathSearchParams,
) -> ButterflyResult {
    run_butterfly_all_prob_paths_with(
        seqvertex_graph,
        orig_kmer_graph,
        orig_bfly_graph,
        combined_read_hash,
        params,
        &crate::postprocess::PostProcessParams::default(),
    )
}

/// 同上，带后处理参数（cd-hit 阈值等）。
///
/// * `seqvertex_graph`：POG/zip 产物（路径搜索用）
/// * `orig_kmer_graph`：原始 kmer 图（边权重重标定用）
/// * `orig_bfly_graph`：剪枝链后的压缩图（Java 驱动里的 `graph`——
///   T9 后处理与 printFinalPaths 都以 orig id 在它上面取序列/finish time）
pub fn run_butterfly_all_prob_paths_with(
    seqvertex_graph: &DiGraph,
    orig_kmer_graph: &DiGraph,
    orig_bfly_graph: &DiGraph,
    combined_read_hash: &ReadHash,
    params: &PathSearchParams,
    post: &crate::postprocess::PostProcessParams,
) -> ButterflyResult {
    use crate::postprocess::{
        assign_compatible_reads_to_paths, convert_to_orig_ids, group_paths_into_genes,
        java_hashmap_order, java_hashmap_order_cap, java_putall_cap, reduce_cdhit_like,
        run_em_reduce,
    };

    // Java L946：POG 之后、分量处理之前，对 SeqVertex DAG 跑 My_DFS(runDFS2)
    //（findLastSharedNode 依赖 finish time；addSandT 的 S/T 保持默认 -1）
    let mut graph = seqvertex_graph.clone();
    run_dfs2(&mut graph);

    // Java 顺序：reorganize（L961，全 1 边权）→ relabel（L1000）→ 分量（L973）
    let read_hash = reorganize_read_pairings(&graph, combined_read_hash);
    relabel_edge_weights_using_orig_kmers(&mut graph, orig_kmer_graph);

    // NUM_MISMATCHES_HASH：全程共享的比对缓存
    let mut mismatch_cache = FxHashMap::default();

    let mut collection: Vec<Vec<i32>> = Vec::new();
    let mut collection_reads: FxHashMap<Vec<i32>, FxHashMap<PairPath, i64>> = FxHashMap::default();

    for comp in divide_into_components(&graph) {
        let mut comp_reads = get_component_reads(&comp, &read_hash);
        if comp_reads.is_empty() {
            continue;
        }
        reduce_to_max_paths_per_node(&mut comp_reads, params.max_paths_per_node_init);
        let triplets = extract_triplets_from_reads(&comp_reads);
        let x_structures = get_xstructures_resolved_by_triplets(&graph, &comp, &triplets);
        let extended = extract_complex_path_prefixes_from_reads(&comp_reads);

        add_s_and_t(
            &mut graph,
            &comp,
            &mut comp_reads,
            params.min_read_support_thr,
        );

        let result = get_all_probable_paths(
            &graph,
            &comp,
            &comp_reads,
            &triplets,
            &extended,
            &x_structures,
            params,
        );
        let final_paths = remove_short_seqs(
            &graph,
            result.final_paths,
            params.min_output_seq,
            params.kmer_size,
        );
        if final_paths.is_empty() {
            continue;
        }

        // ---- T9 后处理（驱动 L1224-1300）----
        // assignCompatibleReadsToPaths（含 addSandT 注入的 sink 读）
        let contained = assign_compatible_reads_to_paths(&final_paths, &comp_reads);
        // 删除未被分配任何 read 的路径
        let final_paths: Vec<Vec<i32>> = final_paths
            .into_iter()
            .filter(|p| contained.contains_key(p))
            .collect();

        // orig id 转换（路径 + contained reads 键）
        let (orig_paths, orig_contained) = convert_to_orig_ids(&graph, &final_paths, &contained);

        // 分量内 cd-hit 式去冗余
        let mut orig_contained = orig_contained;
        let orig_paths = if !post.no_path_merging && orig_paths.len() > 1 {
            reduce_cdhit_like(
                orig_bfly_graph,
                orig_paths,
                &mut orig_contained,
                params.kmer_size,
                post,
                &mut mismatch_cache,
            )
        } else {
            orig_paths
        };

        collection.extend(orig_paths);
        for (k, v) in orig_contained {
            collection_reads.insert(k, v);
        }
    }

    // 全集 cd-hit（L1328）
    let collection = if !post.no_path_merging && collection.len() > 1 {
        reduce_cdhit_like(
            orig_bfly_graph,
            collection,
            &mut collection_reads,
            params.kmer_size,
            post,
            &mut mismatch_cache,
        )
    } else {
        collection
    };

    // 基因分组（L1337）——先按 Java HashMap 序重排（printFinalPaths 的
    // isoform 编号与 golden 对齐依赖此序）
    let order = java_hashmap_order(&collection);
    let collection: Vec<Vec<i32>> = order.into_iter().map(|i| collection[i].clone()).collect();
    let gene_ids = group_paths_into_genes(orig_bfly_graph, &collection, post);

    // EM 削减（L1365；Java 默认运行，Trinity 主脚本传 --NO_EM_REDUCE 跳过）
    let final_paths = run_em_reduce(
        orig_bfly_graph,
        collection.clone(),
        &collection_reads,
        &gene_ids,
        params.kmer_size,
        post,
    );

    // printFinalPaths 的迭代序：EM 路径下最终 map 是
    // `filtered_paths_to_keep = new HashMap(); putAll(EM 结果)` 构建的——
    // 空表 putAll 的容量 = tableSizeFor(n/0.75+1)（非默认 16），c0 的
    // 黄金 i1/i2 正是同桶（cap 8 下 1&7 == 9&7 == 1）相邻序。
    let gene_by_path: FxHashMap<Vec<i32>, usize> = collection
        .iter()
        .zip(gene_ids.iter())
        .map(|(p, &g)| (p.clone(), g))
        .collect();
    let order = if post.no_em_reduce {
        java_hashmap_order(&final_paths)
    } else {
        java_hashmap_order_cap(&final_paths, java_putall_cap(final_paths.len()))
    };
    let final_paths: Vec<Vec<i32>> = order.iter().map(|&i| final_paths[i].clone()).collect();
    let gene_ids: Vec<usize> = final_paths.iter().map(|p| gene_by_path[p]).collect();

    remove_all_edges_of_s_and_t(&mut graph);
    ButterflyResult {
        graph,
        final_paths,
        gene_ids,
        contained_reads: collection_reads,
    }
}

/// allProbPaths.fasta 雏形（printFinalPaths L8960 的简化版）：
/// `>{comp}_g{i}_i{j} len=... path=[<id>'<seqlen>'...]`，60 列折行。
/// T9 完成基因分组前先用序号（每条路径独立 g{i}_i1）。
pub fn paths_to_fasta(
    graph: &DiGraph,
    paths: &[Vec<i32>],
    comp_name: &str,
    kmer_size: usize,
) -> String {
    let mut out = String::new();
    for (i, path) in paths.iter().enumerate() {
        let seq = get_path_seq(graph, path, kmer_size);
        let path_desc: Vec<String> = path
            .iter()
            .filter(|&&id| id >= 0)
            .map(|&id| {
                let len = name_kmer_adj_len(graph, id, kmer_size);
                format!("{id}'{len}'")
            })
            .collect();
        out.push_str(&format!(
            ">{}_g{}_i1 len={} path=[{}]\n",
            comp_name,
            i + 1,
            seq.len(),
            path_desc.join("")
        ));
        for chunk in seq.as_bytes().chunks(60) {
            out.push_str(std::str::from_utf8(chunk).unwrap());
            out.push('\n');
        }
    }
    out
}

// ---------------------------------------------------------------------------
// tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn v(id: i32, name: &str) -> SeqVertex {
        SeqVertex::new(id, name.to_string())
    }

    /// 链图 1→2→3→4 + 旁路 2→5→4。
    fn chain_graph() -> DiGraph {
        let mut g = DiGraph::new();
        for (id, name) in [
            (1, "AAAAACCCCCGGGGGTTTTT"),
            (2, "CCCCCGGGGGTTTTTAAAAA"),
            (3, "GGGGGTTTTTAAAAACCCCC"),
            (4, "TTTTTAAAAACCCCCGGGGG"),
            (5, "GGGGGTTTTTCCCCCAAAAA"),
        ] {
            g.add_vertex(v(id, name));
        }
        for (a, b, w) in [
            (1, 2, 1.0),
            (2, 3, 2.0),
            (3, 4, 3.0),
            (2, 5, 1.0),
            (5, 4, 1.0),
        ] {
            g.add_edge(a, b, SimpleEdge::new(w, a, b));
        }
        g
    }

    #[test]
    fn wo_ver_dijkstra_skips_excluded_vertex() {
        let g = chain_graph();
        // 2→4 最短经 5：1+1=2
        assert_eq!(get_distance_wo_ver(&g, 2, 4, 99), Some(2.0));
        // 删 5 → 只能经 3：2+3=5
        assert_eq!(get_distance_wo_ver(&g, 2, 4, 5), Some(5.0));
        // 删 3 不影响 2→4
        assert_eq!(get_distance_wo_ver(&g, 2, 4, 3), Some(2.0));
        // 删 4 → 不可达
        assert_eq!(get_distance_wo_ver(&g, 2, 4, 4), None);
        // 逆向不可达（DAG）
        assert_eq!(get_distance_wo_ver(&g, 4, 1, 99), None);
        assert_eq!(get_distance_wo_ver(&g, 3, 3, 3), Some(0.0));
    }

    #[test]
    fn is_ancestral_directionality() {
        let g = chain_graph();
        assert_eq!(is_ancestral(&g, 1, 4), 1);
        assert_eq!(is_ancestral(&g, 4, 1), -1);
        assert_eq!(is_ancestral(&g, 3, 5), 0); // 兄弟分支
    }

    #[test]
    fn share_suffix_java_integer_cache_quirk() {
        // 小 id（缓存区）：数值比较生效——重叠后缀逐位相等（一侧先走完）即 true
        assert!(share_suffix_fully_contained_java(&[2, 3], &[7, 9, 2, 3]));
        // 走到 (1,9)：数值不等 → false
        assert!(!share_suffix_fully_contained_java(&[1, 2, 3], &[9, 2, 3]));
        // 大 id（>127）：跨列表 Integer 引用比较 → 首对即失败（Java quirk）
        assert!(!share_suffix_fully_contained_java(&[200, 201], &[200, 201]));
    }

    #[test]
    fn triplet_extraction_and_support() {
        let mut hash: ReadHash = FxHashMap::default();
        let mut m = FxHashMap::default();
        m.insert(PairPath::new(vec![1, 2, 3, 4]), 3);
        m.insert(PairPath::new(vec![5, 2, 7]), 1);
        hash.insert(1, m);
        let t = extract_triplets_from_reads(&hash);
        assert!(triplet_supported(&t[&2], &[1, 2, 3]));
        assert!(triplet_supported(&t[&3], &[2, 3, 4]));
        assert!(triplet_supported(&t[&2], &[5, 2, 7]));
        assert!(!triplet_supported(&t[&2], &[9, 2, 3]));
        // 复杂前缀：节点 4 有 [1,2,3,4]，节点 3 有 [1,2,3]（purge 因 quirk 不触发）
        let ext = extract_complex_path_prefixes_from_reads(&hash);
        assert!(ext[&3].contains(&vec![1, 2, 3]));
        assert!(ext[&4].contains(&vec![1, 2, 3, 4]));
    }

    #[test]
    fn add_s_and_t_connects_terminals_and_reads() {
        let mut g = chain_graph();
        for id in 1..=5 {
            let vert = g.get_vertex_mut(id).unwrap();
            vert.weights = vec![2.0; vert.name.len()];
        }
        let comp: FxHashSet<i32> = [1, 2, 3, 4, 5].into_iter().collect();
        let mut reads: ReadHash = FxHashMap::default();
        add_s_and_t(&mut g, &comp, &mut reads, 2);

        assert!(g.contains_vertex(VERTEX_ROOT_ID) && g.contains_vertex(T_VERTEX_ID));
        assert_eq!(g.get_successors(VERTEX_ROOT_ID), &[1]); // 只有 1 入度 0
        assert_eq!(g.get_predecessors(T_VERTEX_ID), &[4]); // 只有 4 出度 0
                                                           // ROOT read hash: [-1,1]→2
        assert_eq!(reads[&VERTEX_ROOT_ID][&PairPath::new(vec![-1, 1])], 2);
        assert_eq!(reads[&4][&PairPath::new(vec![4, -2])], 2);

        remove_all_edges_of_s_and_t(&mut g);
        assert_eq!(g.out_degree(VERTEX_ROOT_ID), 0);
        assert_eq!(g.in_degree(T_VERTEX_ID), 0);
    }

    #[test]
    fn combine_paths_orders_and_imputes() {
        let g = chain_graph();
        // p1 在 p2 上游（1→3 存在路径经 2）
        let combined = combine_paths(&g, &[1], &[3, 4]).unwrap();
        // imputation：1 的唯一可连通后继是 2，2 的是 3 → [1,2,3,4]
        assert_eq!(combined.path1, vec![1, 2, 3, 4]);
        assert!(!combined.has_second_path());
        // 逆向也可组合（1→4 可达 → path2 在前）；imputation 因 2 有两个
        // 可达 4 的后继（3/5）而放弃 → 保持配对
        let c = combine_paths(&g, &[4], &[1]).unwrap();
        assert_eq!(c.path1, vec![1]);
        assert_eq!(c.path2, vec![4]);
        // 兄弟分支（3 与 5 互不可达）不可组合
        assert!(combine_paths(&g, &[3], &[5]).is_none());
        // 包含关系
        let c = combine_paths(&g, &[1, 2, 3], &[2, 3]).unwrap();
        assert_eq!(c.path1, vec![1, 2, 3]);
    }

    #[test]
    fn reorganize_splits_unconnectable_pairs() {
        let g = chain_graph();
        let mut hash: ReadHash = FxHashMap::default();
        let mut m = FxHashMap::default();
        m.insert(PairPath::with_pair(vec![3], vec![5]), 5); // 兄弟分支，不可组合
        hash.insert(3, m);
        let out = reorganize_read_pairings(&g, &hash);
        // 拆成两条单端，各 support 5，按首节点存
        assert_eq!(out[&3][&PairPath::new(vec![3])], 5);
        assert_eq!(out[&5][&PairPath::new(vec![5])], 5);
    }

    #[test]
    fn reduce_to_max_keeps_top_supported() {
        let mut m: FxHashMap<PairPath, i64> = FxHashMap::default();
        m.insert(PairPath::new(vec![1, 2]), 5);
        m.insert(PairPath::new(vec![1, 3]), 9);
        m.insert(PairPath::new(vec![1, 4]), 7);
        let mut hash: ReadHash = FxHashMap::default();
        hash.insert(1, m);
        reduce_to_max_paths_per_node(&mut hash, 2);
        let kept = &hash[&1];
        assert_eq!(kept.len(), 2);
        assert!(kept.contains_key(&PairPath::new(vec![1, 3])));
        assert!(kept.contains_key(&PairPath::new(vec![1, 4])));
    }

    #[test]
    fn divide_into_components_weak() {
        let mut g = chain_graph();
        g.add_vertex(v(6, "TTTTTGGGGGCCCCCAAAAA"));
        g.add_edge(6, 1, SimpleEdge::new(1.0, 6, 1));
        let comps = divide_into_components(&g);
        assert_eq!(comps.len(), 1); // 6→1 弱连通
        g.remove_edge(6, 1);
        assert_eq!(divide_into_components(&g).len(), 2);
    }

    #[test]
    fn remove_identical_subseqs_drops_contained() {
        let g = chain_graph();
        // [1,2] 序列是 [1,2,3] 序列的前缀（kmer 邻接拼接）→ 被包含删除
        let full = vec![1, 2, 3];
        let sub = vec![1, 2];
        let res = remove_identical_subseqs(&g, vec![sub, full.clone()], 5);
        assert_eq!(res, vec![full.clone()]);
        // 非包含的两条都保留（2-5 与 1-2-3 起点不同）
        let res = remove_identical_subseqs(&g, vec![vec![2, 5], full], 5);
        assert_eq!(res.len(), 2);
    }

    #[test]
    fn get_path_seq_kmer_adj_trimming() {
        let g = chain_graph(); // kmer len 20 → 截 19
        let seq = get_path_seq(&g, &[VERTEX_ROOT_ID, 1, 2, T_VERTEX_ID], 20);
        let n1 = &g.get_vertex(1).unwrap().name;
        let n2 = &g.get_vertex(2).unwrap().name;
        assert_eq!(seq, format!("{}{}", n1, &n2[19..]));
    }

    #[test]
    fn simple_bubble_search() {
        // S→1, 1→2, 1→3, 2→T, 3→T；read [1,2]、[1,3] 各 2 条
        let mut g = DiGraph::new();
        g.add_vertex(v(1, &"A".repeat(30)));
        g.add_vertex(v(2, &"C".repeat(30)));
        g.add_vertex(v(3, &"G".repeat(30)));
        g.add_vertex(v(VERTEX_ROOT_ID, "S"));
        g.add_vertex(v(T_VERTEX_ID, "E"));
        for (a, b) in [
            (VERTEX_ROOT_ID, 1),
            (1, 2),
            (1, 3),
            (2, T_VERTEX_ID),
            (3, T_VERTEX_ID),
        ] {
            g.add_edge(a, b, SimpleEdge::new(1.0, a, b));
        }
        for id in [1, 2, 3] {
            g.get_vertex_mut(id).unwrap().depth = id;
        }
        g.get_vertex_mut(VERTEX_ROOT_ID).unwrap().depth = -1;
        g.get_vertex_mut(T_VERTEX_ID).unwrap().depth = 99;

        let mut reads: ReadHash = FxHashMap::default();
        let mut m = FxHashMap::default();
        m.insert(PairPath::new(vec![1, 2]), 2);
        m.insert(PairPath::new(vec![1, 3]), 2);
        reads.insert(1, m);
        let mut root_reads = FxHashMap::default();
        root_reads.insert(PairPath::new(vec![VERTEX_ROOT_ID, 1]), 2);
        reads.insert(VERTEX_ROOT_ID, root_reads);
        let mut t_reads = FxHashMap::default();
        t_reads.insert(PairPath::new(vec![2, T_VERTEX_ID]), 2);
        t_reads.insert(PairPath::new(vec![3, T_VERTEX_ID]), 2);
        reads.insert(2, t_reads);

        let comp: FxHashSet<i32> = [1, 2, 3].into_iter().collect();
        let triplets = extract_triplets_from_reads(&reads);
        let ext = extract_complex_path_prefixes_from_reads(&reads);
        let params = PathSearchParams {
            kmer_size: 25,
            min_output_seq: 0, // 无长度门槛
            min_read_support_thr: 2,
            path_reinforcement_distance: 30,
            ..PathSearchParams::default()
        };
        let res = get_all_probable_paths(
            &g,
            &comp,
            &reads,
            &triplets,
            &ext,
            &FxHashMap::default(),
            &params,
        );
        // 两条路径 [-1,1,2,-2] / [-1,1,3,-2]
        assert_eq!(res.final_paths.len(), 2, "{:?}", res.final_paths);
        assert!(res
            .final_paths
            .contains(&vec![VERTEX_ROOT_ID, 1, 2, T_VERTEX_ID]));
        assert!(res
            .final_paths
            .contains(&vec![VERTEX_ROOT_ID, 1, 3, T_VERTEX_ID]));
    }
}
