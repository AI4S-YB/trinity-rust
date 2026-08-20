//! T9 后处理：assignCompatibleReadsToPaths / convert_to_orig_ids /
//! reduce_cdhit_like（CD-HIT 式去冗余：twoPathsAreTooSimilar +
//! getPrevCalcNumMismatches + findLastSharedNode）/ group_paths_into_genes /
//! printFinalPaths（`_g{i}_i{j}` 命名 + `path=[...]` MISO 格式）。
//!
//! 镜像 TransAssembly_allProbPaths.java L15392（reduce_cdhit_like）、
//! L10608（twoPathsAreTooSimilar）、L10748（getPrevCalcNumMismatches）、
//! L10683（findLastSharedNode）、L9129（group_paths_into_genes）、
//! L8954（printFinalPaths）、L1438（convert_to_orig_ids）、L15663
//!（get_pathName_string）、L7000（assignCompatibleReadsToPaths）、
//! L11130（removeTheLesserSupportedPath）。
//!
//! 保留的 Java quirk（逐条对应源码）：
//!  1. `DIFFS_WINDOW_SIZE`/`MAX_FINAL_DIFFS_IN_WINDOW` 在 v2.15.2 的
//!     isThisTooSimilar 路径中**未被使用**（参数存在但不参与判定）——
//!     这里同样只作为字段保留。
//!  2. getPrevCalcNumMismatches 的"单侧空路径"分支**不写缓存**（Java 漏了
//!     NUM_MISMATCHES_HASH.put）。
//!  3. reduce_cdhit_like 中若 path_i 被删除，内层 j 循环**继续**用已删除的
//!     path_i 作为过滤证据（Java 不 break）。
//!  4. 去冗余期间的 NW 比对长度上限被临时替换为 ALL_VS_ALL_MAX_DP_LEN(1000)。
//!  5. isThisTooSimilar 只用 `mismatches`（不含 gaps）算 numMM 与 per-id。
//!  6. 输出/排序的迭代顺序按 Java HashMap 桶序仿真（List<Integer> 的
//!     31 进制 hashCode + spread；见 `java_hashmap_order`）。

use std::collections::BTreeMap;

use rustc_hash::{FxHashMap, FxHashSet};

use crate::align::{increment_alignment_stats, run_nw_alignment, zipper_align, AlignmentStats};
use crate::graph::{DiGraph, T_VERTEX_ID, VERTEX_ROOT_ID};
use crate::pair_paths::{have_any_node_in_common_paths, trim_sink_nodes, PairPath};
use crate::paths::{get_path_seq, ReadHash};

/// 后处理参数（CLI 默认值，TransAssembly_allProbPaths.java L65-69 / L133 / L155）。
#[derive(Debug, Clone)]
pub struct PostProcessParams {
    /// `MAX_DIFFS_SAME_PATH`（默认 2）
    pub max_diffs_same_path: i32,
    /// `MIN_PERCENT_IDENTITY_SAME_PATH`（默认 98.0）
    pub min_per_id_same_path: f32,
    /// `MAX_INTERNAL_GAP_SAME_PATH`（默认 10）
    pub max_internal_gap_same_path: usize,
    /// `DIFFS_WINDOW_SIZE`（默认 100；v2.15.2 中本路径未使用，QUIRK 1）
    pub diffs_window_size: usize,
    /// `MAX_FINAL_DIFFS_IN_WINDOW`（默认 5；同上未使用）
    pub max_final_diffs_in_window: usize,
    /// `ALL_VS_ALL_MAX_DP_LEN`（默认 1000）：去冗余期间替换 MAX_SEQ_LEN_DP_ALIGN
    pub all_vs_all_max_dp_len: usize,
    /// `MIN_ISOFORM_PCT_LEN_OVERLAP`（默认 30）
    pub min_isoform_pct_len_overlap: f32,
    /// `MAX_SEQ_LEN_DP_ALIGN`（默认 10000）：非去冗余场景的比对长度上限
    pub max_seq_len_dp_align: usize,
    /// Java 默认 false（跑 EM）；Trinity 主脚本传 `--NO_EM_REDUCE`（L2387）
    pub no_em_reduce: bool,
    /// `--no_path_merging`（L389）：跳过两处 cd-hit 式去冗余（L1295/L1322 的
    /// `!NO_PATH_MERGING && size > 1` 门）
    pub no_path_merging: bool,
    /// `MIN_TOTAL_ISOFORM_EXPRESSION`（默认 0 = 关闭）
    pub min_total_isoform_expression: f32,
    /// `MIN_RELATIVE_ISOFORM_EXPRESSION`（默认 5%）
    pub min_relative_isoform_expression: f32,
}

impl Default for PostProcessParams {
    fn default() -> Self {
        Self {
            max_diffs_same_path: 2,
            min_per_id_same_path: 98.0,
            max_internal_gap_same_path: 10,
            diffs_window_size: 100,
            max_final_diffs_in_window: 5,
            all_vs_all_max_dp_len: 1000,
            min_isoform_pct_len_overlap: 30.0,
            max_seq_len_dp_align: 10000,
            no_em_reduce: false,
            no_path_merging: false,
            min_total_isoform_expression: 0.0,
            min_relative_isoform_expression: 5.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Java HashMap 迭代序仿真（QUIRK 6）
// ---------------------------------------------------------------------------

/// Java `List<Integer>.hashCode()`：`h = 31*h + Integer.hashCode(x)`（i32 回绕）。
pub fn java_list_hash(path: &[i32]) -> i32 {
    path.iter()
        .fold(1i32, |h, &x| h.wrapping_mul(31).wrapping_add(x))
}

/// Java `HashMap` 的 spread：`(h ^ (h >>> 16)) & (cap-1)`。
fn java_spread(h: i32) -> u32 {
    (h as u32) ^ ((h as u32) >> 16)
}

/// 仿真 Java HashMap（默认初始容量 16、负载因子 0.75、尾插链表）在
/// `keySet()` 迭代时给出的条目顺序：按桶号升序、同桶按插入序。
/// 返回原顺序的下标排列。路径数 >8 落同一桶的 treeify 场景未仿真
///（转录本路径长度互异，实际不会触发）。
pub fn java_hashmap_order(paths: &[Vec<i32>]) -> Vec<usize> {
    // 模拟逐条插入时的扩容（size > 0.75*cap 即翻倍）
    let mut cap = 16usize;
    while paths.len() > cap * 3 / 4 {
        cap *= 2;
    }
    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); cap];
    for (i, p) in paths.iter().enumerate() {
        let b = (java_spread(java_list_hash(p)) & (cap as u32 - 1)) as usize;
        buckets[b].push(i);
    }
    let mut order = Vec::with_capacity(paths.len());
    for b in buckets {
        order.extend(b);
    }
    order
}

/// `HashMap.putMapEntries` 的容量推导（putAll 到 table 尚未分配的空 map）：
/// `t = (int)(n/0.75 + 1)`，`cap = tableSizeFor(t)`（≥ t 的最小 2 幂）。
pub fn java_putall_cap(n: usize) -> usize {
    let t = (n as f32 / 0.75 + 1.0) as usize;
    t.max(1).next_power_of_two()
}

/// 同上，但指定表容量（putAll 构建的 map，容量由 `java_putall_cap` 决定
/// 而非默认 16——见 c0：5 条路径 → cap 8，黄金 i1/i2 同桶相邻）。
pub fn java_hashmap_order_cap(paths: &[Vec<i32>], cap: usize) -> Vec<usize> {
    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); cap];
    for (i, p) in paths.iter().enumerate() {
        let b = (java_spread(java_list_hash(p)) & (cap as u32 - 1)) as usize;
        buckets[b].push(i);
    }
    let mut order = Vec::with_capacity(paths.len());
    for b in buckets {
        order.extend(b);
    }
    order
}

/// Java `List.toString()`：`"[1, -1, 2]"`。
pub fn java_list_string(path: &[i32]) -> String {
    let inner: Vec<String> = path.iter().map(|x| x.to_string()).collect();
    format!("[{}]", inner.join(", "))
}

// ---------------------------------------------------------------------------
// assignCompatibleReadsToPaths（L7000）
// ---------------------------------------------------------------------------

/// 每条 final path → 兼容且被包含的 reads（PairPath → 计数）。
pub type ContainedReads = FxHashMap<Vec<i32>, FxHashMap<PairPath, i64>>;

/// Java `assignCompatibleReadsToPaths`：对每条 path，扫全部 read，
/// `p.isCompatibleAndContainedBySinglePath(path)` 者计入。
pub fn assign_compatible_reads_to_paths(
    final_paths: &[Vec<i32>],
    combined_read_hash: &ReadHash,
) -> ContainedReads {
    let mut path_to_contained: ContainedReads = FxHashMap::default();
    for path in final_paths {
        for read_map in combined_read_hash.values() {
            for (p, &count) in read_map {
                if p.is_compatible_and_contained_by_single_path(path) {
                    path_to_contained
                        .entry(path.clone())
                        .or_default()
                        .insert(p.clone(), count);
                }
            }
        }
    }
    path_to_contained
}

// ---------------------------------------------------------------------------
// convert_to_orig_ids（L1438）
// ---------------------------------------------------------------------------

/// Java `PairPath.setOrigIds()`：两条路径逐 id 映射到 origButterflyID。
fn pair_path_set_orig_ids(pp: &PairPath, graph: &DiGraph) -> PairPath {
    let map = |p: &[i32]| -> Vec<i32> {
        p.iter()
            .map(|&id| {
                graph
                    .get_vertex(id)
                    .map(|v| v.orig_butterfly_id)
                    .unwrap_or(id)
            })
            .collect()
    };
    let mut out = PairPath::new(map(&pp.path1));
    if !pp.path2.is_empty() {
        out = PairPath::with_pair(map(&pp.path1), map(&pp.path2));
    }
    out
}

/// Java `convert_to_orig_ids`：path 与 contained-reads 键都换成 orig id
///（计数合并——不同 zipped id 可能映射到同一 orig id 路径）。
pub fn convert_to_orig_ids(
    graph: &DiGraph,
    final_paths: &[Vec<i32>],
    contained: &ContainedReads,
) -> (Vec<Vec<i32>>, ContainedReads) {
    let mut out_paths: Vec<Vec<i32>> = Vec::new();
    let mut seen: FxHashSet<Vec<i32>> = FxHashSet::default();
    let mut out_contained: ContainedReads = FxHashMap::default();

    for path in final_paths {
        let revised: Vec<i32> = path
            .iter()
            .map(|&id| {
                graph
                    .get_vertex(id)
                    .map(|v| v.orig_butterfly_id)
                    .unwrap_or(id)
            })
            .collect();
        if seen.insert(revised.clone()) {
            out_paths.push(revised.clone());
        }
        if let Some(reads) = contained.get(path) {
            let entry = out_contained.entry(revised).or_default();
            for (pp, &count) in reads {
                *entry.entry(pair_path_set_orig_ids(pp, graph)).or_default() += count;
            }
        }
    }
    (out_paths, out_contained)
}

// ---------------------------------------------------------------------------
// getPrevCalcNumMismatches（L10748）+ findLastSharedNode（L10683）
// ---------------------------------------------------------------------------

/// `SeqVertexFinishTimeComparator`：finish time 降序为"升序"。
fn finish_time_cmp(graph: &DiGraph, a: i32, b: i32) -> std::cmp::Ordering {
    let fa = graph.get_vertex(a).map(|v| v.dfs_finish_time).unwrap_or(-1);
    let fb = graph.get_vertex(b).map(|v| v.dfs_finish_time).unwrap_or(-1);
    // Java: f1<f2 → 1, f1>f2 → -1, == → 0
    fb.cmp(&fa)
}

/// Java `findLastSharedNode`：sink 修剪后，从两条路径尾部各持一指针，
/// 按 finish time 双指针后退，相遇即最后共享节点。
fn find_last_shared_node(graph: &DiGraph, path1: &[i32], path2: &[i32]) -> i32 {
    let p1 = trim_sink_nodes(path1);
    let p2 = trim_sink_nodes(path2);
    if p1.is_empty() || p2.is_empty() {
        return -1;
    }
    // 反向遍历（getReverseSeqVertexPath）
    let mut i1 = 0usize; // p1 反向游标
    let mut i2 = 0usize;
    let mut v1 = p1[p1.len() - 1 - i1];
    let mut v2 = p2[p2.len() - 1 - i2];
    while v1 != v2 {
        if finish_time_cmp(graph, v1, v2) != std::cmp::Ordering::Less {
            // compare >= 0 → 前进 p1
            if i1 + 1 < p1.len() {
                i1 += 1;
                v1 = p1[p1.len() - 1 - i1];
            } else {
                break;
            }
        } else if i2 + 1 < p2.len() {
            i2 += 1;
            v2 = p2[p2.len() - 1 - i2];
        } else {
            break;
        }
    }
    if v1 == v2 {
        v1
    } else {
        -1
    }
}

/// Java `String.CASE_INSENSITIVE_ORDER`（路径串只含数字/`[`/`]`/`,`/` `/`-`，
/// 等价于逐字符比较，再比长度）。
fn case_insensitive_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    let ab = a.as_bytes();
    let bb = b.as_bytes();
    for i in 0..ab.len().min(bb.len()) {
        let (ca, cb) = (ab[i].to_ascii_lowercase(), bb[i].to_ascii_lowercase());
        match ca.cmp(&cb) {
            std::cmp::Ordering::Equal => continue,
            o => return o,
        }
    }
    ab.len().cmp(&bb.len())
}

/// getPrevCalcNumMismatches 的比对上下文：当前生效的 DP 长度上限
///（cd-hit 期间被临时替换为 ALL_VS_ALL_MAX_DP_LEN，QUIRK 4）。
#[derive(Clone, Copy)]
struct AlignCtx {
    kmer_size: usize,
    max_dp_len: usize,
}

fn path_is_effectively_empty(path: &[i32]) -> bool {
    path.is_empty() || (path.len() == 1 && path[0] < 0)
}

/// Java `getPrevCalcNumMismatches`：缓存 + 三段分解（prefix / shared /
/// suffix）+ 无共享节点时的序列比对（含 gap 归属三规则）。
fn get_prev_calc_num_mismatches(
    graph: &DiGraph,
    path1: &[i32],
    path2: &[i32],
    ctx: AlignCtx,
    cache: &mut FxHashMap<String, AlignmentStats>,
) -> AlignmentStats {
    let is_at_start =
        (path1.first() == Some(&VERTEX_ROOT_ID)) || (path2.first() == Some(&VERTEX_ROOT_ID));
    let is_at_end = (path1.last() == Some(&T_VERTEX_ID)) || (path2.last() == Some(&T_VERTEX_ID));

    let s1 = java_list_string(path1);
    let s2 = java_list_string(path2);
    let key = if case_insensitive_cmp(&s1, &s2) != std::cmp::Ordering::Less {
        format!("{s1};{s2}")
    } else {
        format!("{s2};{s1}")
    };
    if let Some(a) = cache.get(&key) {
        return a.clone();
    }

    // 两侧都（有效地）为空 → 全零
    if path_is_effectively_empty(path1) && path_is_effectively_empty(path2) {
        let a = AlignmentStats::default();
        cache.insert(key, a.clone());
        return a;
    }

    // 路径相同 → 完美匹配
    if path1 == path2 {
        let mut a = AlignmentStats::default();
        let len = get_path_seq(graph, path1, ctx.kmer_size).len();
        a.matches = len;
        a.alignment_length = len;
        cache.insert(key, a.clone());
        return a;
    }

    // 单侧为空：另一侧全长（减 k-1）计为 gap（QUIRK 2：不写缓存）
    if path_is_effectively_empty(path1) {
        let mut a = AlignmentStats::default();
        let l = get_path_seq(graph, path2, ctx.kmer_size).len();
        let l = l.saturating_sub(ctx.kmer_size - 1);
        a.max_internal_gap_length = l;
        a.gaps = l;
        return a;
    }
    if path_is_effectively_empty(path2) {
        let mut a = AlignmentStats::default();
        let l = get_path_seq(graph, path1, ctx.kmer_size).len();
        let l = l.saturating_sub(ctx.kmer_size - 1);
        a.max_internal_gap_length = l;
        a.gaps = l;
        return a;
    }

    let last_shared = find_last_shared_node(graph, path1, path2);
    if last_shared != -1 {
        // 三段分解：prefix（不含共享节点，递归）+ shared（单节点，恒完美）
        // + suffix（递归），stats 相加（increment_alignment_stats）。
        let p1_idx = path1.iter().position(|&x| x == last_shared).unwrap();
        let p2_idx = path2.iter().position(|&x| x == last_shared).unwrap();

        let prefix_stats =
            get_prev_calc_num_mismatches(graph, &path1[..p1_idx], &path2[..p2_idx], ctx, cache);
        let shared_stats = get_prev_calc_num_mismatches(
            graph,
            &path1[p1_idx..p1_idx + 1],
            &path2[p2_idx..p2_idx + 1],
            ctx,
            cache,
        );
        let suffix_1: &[i32] = if p1_idx < path1.len() - 1 {
            &path1[p1_idx + 1..]
        } else {
            &[]
        };
        let suffix_2: &[i32] = if p2_idx < path2.len() - 1 {
            &path2[p2_idx + 1..]
        } else {
            &[]
        };
        let suffix_stats = get_prev_calc_num_mismatches(graph, suffix_1, suffix_2, ctx, cache);

        let combined = increment_alignment_stats(
            &increment_alignment_stats(&suffix_stats, &shared_stats),
            &prefix_stats,
        );
        cache.insert(key, combined.clone());
        return combined;
    }

    // 无共享节点：真正比对
    let seq1 = get_path_seq(graph, path1, ctx.kmer_size);
    let seq2 = get_path_seq(graph, path2, ctx.kmer_size);
    let (l1, l2) = (seq1.len(), seq2.len());

    let mut stats;
    if (l1 > ctx.max_dp_len && l2 > ctx.max_dp_len) || l1 > 100_000 || l2 > 100_000 {
        // 超长：拉链启发式
        stats = zipper_align(seq1.as_bytes(), seq2.as_bytes());
    } else if (l1 < 10 && l2 > 20) || (l1 > 20 && l2 < 10) {
        // 一侧极短：图首锚右 / 图尾锚左 / 否则拉链
        stats = if is_at_start {
            crate::align::zipper_anchor_right(seq1.as_bytes(), seq2.as_bytes())
        } else if is_at_end {
            crate::align::zipper_anchor_left(seq1.as_bytes(), seq2.as_bytes())
        } else {
            zipper_align(seq1.as_bytes(), seq2.as_bytes())
        };
    } else {
        // Needleman-Wunsch（Gotoh）全局比对，4/-5/10/1
        let (_aln, st) = run_nw_alignment(seq1.as_bytes(), seq2.as_bytes(), 4.0, -5.0, 10.0, 1.0);
        stats = st;
    }

    // ---- gap 归属三规则 ----
    let mut total_significant_diffs;
    if is_at_start || is_at_end {
        total_significant_diffs = stats.mismatches + stats.gaps;
        if is_at_start {
            stats.left_gap_length = 0;
            if !is_at_end {
                stats.max_internal_gap_length =
                    stats.max_internal_gap_length.max(stats.right_gap_length);
                total_significant_diffs += stats.right_gap_length;
                stats.gaps += stats.right_gap_length;
            }
        }
        if is_at_end {
            stats.right_gap_length = 0;
            if !is_at_start {
                stats.max_internal_gap_length =
                    stats.max_internal_gap_length.max(stats.left_gap_length);
                total_significant_diffs += stats.left_gap_length;
                stats.gaps += stats.left_gap_length;
            }
        }
    } else {
        // 图内部：所有 gap 都算
        total_significant_diffs =
            stats.mismatches + stats.gaps + stats.left_gap_length + stats.right_gap_length;
        stats.max_internal_gap_length = stats
            .max_internal_gap_length
            .max(stats.left_gap_length)
            .max(stats.right_gap_length);
    }
    stats.total_not_matched = total_significant_diffs;

    cache.insert(key, stats.clone());
    stats
}

/// Java `isThisTooSimilar`（L11098）：`gap<=10 && (numMM<=2 || per_id>=98)`。
/// 注意 numMM 只取 mismatches（不含 gaps，QUIRK 5）；窗口参数未参与。
fn is_this_too_similar(
    num_mm: usize,
    max_internal_gap_length: usize,
    percent_identity: f32,
    params: &PostProcessParams,
) -> bool {
    max_internal_gap_length <= params.max_internal_gap_same_path
        && ((num_mm as i32) <= params.max_diffs_same_path
            || percent_identity >= params.min_per_id_same_path)
}

/// Java `twoPathsAreTooSimilar`（L10608）。
pub fn two_paths_are_too_similar(
    graph: &DiGraph,
    path1: &[i32],
    path2: &[i32],
    kmer_size: usize,
    params: &PostProcessParams,
    cache: &mut FxHashMap<String, AlignmentStats>,
) -> bool {
    if !have_any_node_in_common_paths(path1, path2) {
        return false; // 无公共节点不可能相似
    }
    let ctx = AlignCtx {
        kmer_size,
        max_dp_len: params.all_vs_all_max_dp_len, // cd-hit 期间的临时上限
    };
    let stats = get_prev_calc_num_mismatches(graph, path1, path2, ctx, cache);
    let len1 = get_path_seq(graph, path1, kmer_size).len();
    let len2 = get_path_seq(graph, path2, kmer_size).len();
    let shorter_len = len1.min(len2);
    let path_per_id = 100.0 - stats.mismatches as f32 / shorter_len as f32 * 100.0;
    is_this_too_similar(
        stats.mismatches,
        stats.max_internal_gap_length,
        path_per_id,
        params,
    )
}

// ---------------------------------------------------------------------------
// reduce_cdhit_like（L15392）+ removeTheLesserSupportedPath（L11130）
// ---------------------------------------------------------------------------

/// Java `ensure_path_has_sinks`（L15359）：补 -1 头 / -2 尾。
fn ensure_path_has_sinks(path: &[i32]) -> Vec<i32> {
    let mut p: Vec<i32> = path.to_vec();
    if p.first() != Some(&VERTEX_ROOT_ID) {
        p.insert(0, VERTEX_ROOT_ID);
    }
    if p.last() != Some(&T_VERTEX_ID) {
        p.push(T_VERTEX_ID);
    }
    p
}

/// Java `removeTheLesserSupportedPath`：低 read 支持者删；平手保长者
///（仍平保前者）。返回 1 = 删 path1、2 = 删 path2；并从 PathReads 移除被删者。
fn remove_the_lesser_supported_path(
    seq1_len: usize,
    seq2_len: usize,
    path1: &[i32],
    path2: &[i32],
    path_reads: &mut ContainedReads,
) -> u8 {
    let sum = |p: &[i32]| -> i64 { path_reads.get(p).map(|m| m.values().sum()).unwrap_or(0) };
    let (sum1, sum2) = (sum(path1), sum(path2));
    let remove_first = if sum1 < sum2 {
        true
    } else if sum1 > sum2 {
        false
    } else {
        // 平手：保长者（>= 保 path1）
        seq1_len < seq2_len
    };
    let path2remove = if remove_first { path1 } else { path2 };
    path_reads.remove(path2remove);
    if remove_first {
        1
    } else {
        2
    }
}

/// Java `reduce_cdhit_like`：按序列长度降序（稳定，平序按 Java HashMap 桶序），
/// 两两 twoPathsAreTooSimilar 者删低支持。返回保留路径（输入序）。
pub fn reduce_cdhit_like(
    graph: &DiGraph,
    final_paths: Vec<Vec<i32>>,
    path_reads: &mut ContainedReads,
    kmer_size: usize,
    params: &PostProcessParams,
    cache: &mut FxHashMap<String, AlignmentStats>,
) -> Vec<Vec<i32>> {
    // HashMap 迭代序 → FinalPaths 序（稳定排序的平序基础）
    let hm_order = java_hashmap_order(&final_paths);
    let mut path_vec: Vec<(Vec<i32>, String)> = hm_order
        .into_iter()
        .map(|i| {
            let p = final_paths[i].clone();
            let s = get_path_seq(graph, &p, kmer_size);
            (p, s)
        })
        .collect();
    // FinalPaths.compareTo：长度降序（稳定）
    path_vec.sort_by_key(|(_, seq)| std::cmp::Reverse(seq.len()));

    let mut filtered: FxHashSet<usize> = FxHashSet::default();
    for i in 0..path_vec.len().saturating_sub(1) {
        if filtered.contains(&i) {
            continue; // 被删者不能作为过滤证据（但见 QUIRK 3：仍继续用）
        }
        let path_i_sinks = ensure_path_has_sinks(&path_vec[i].0);
        for j in (i + 1)..path_vec.len() {
            if filtered.contains(&j) {
                continue;
            }
            let path_j_sinks = ensure_path_has_sinks(&path_vec[j].0);
            if two_paths_are_too_similar(
                graph,
                &path_i_sinks,
                &path_j_sinks,
                kmer_size,
                params,
                cache,
            ) {
                let r = remove_the_lesser_supported_path(
                    path_vec[i].1.len(),
                    path_vec[j].1.len(),
                    &path_vec[i].0,
                    &path_vec[j].0,
                    path_reads,
                );
                if r == 1 {
                    // QUIRK 3：i 被删后 j 循环继续，path_i 仍作证据
                    filtered.insert(i);
                } else {
                    filtered.insert(j);
                }
            }
        }
    }

    (0..path_vec.len())
        .filter(|&i| !filtered.contains(&i))
        .map(|i| path_vec[i].0.clone())
        .collect()
}

// ---------------------------------------------------------------------------
// group_paths_into_genes（L9129）
// ---------------------------------------------------------------------------

/// Java `group_paths_into_genes`：节点长度重叠率（双向取 max）≥30% 建边，
/// 弱连通分量 = 基因。返回与输入对齐的基因编号（1 起）。
///
/// 分量编号次序的偏差说明：Java 用 jung WeakComponentClusterer 对 HashSet
/// 聚类迭代（HashSet 次序不可复现），这里按分量成员在 Java HashMap 序中的
/// 最小出现次序编号——单分量（如 c0）与黄金一致；多分量时只影响 g 编号
/// 分配次序，不影响分组本身。
pub fn group_paths_into_genes(
    graph: &DiGraph,
    final_paths: &[Vec<i32>],
    params: &PostProcessParams,
) -> Vec<usize> {
    // 节点长度表（nodeID>0 取 name 长，sink 记 0）
    let mut node_length: FxHashMap<i32, usize> = FxHashMap::default();
    for path in final_paths {
        for &id in path {
            let l = if id > 0 {
                graph.get_vertex(id).map(|v| v.name.len()).unwrap_or(0)
            } else {
                0
            };
            node_length.insert(id, l);
        }
    }

    let hm_order = java_hashmap_order(final_paths);
    let paths: Vec<&Vec<i32>> = hm_order.iter().map(|&i| &final_paths[i]).collect();

    // 并查集（无向边 = 弱连通）
    let mut parent: Vec<usize> = (0..paths.len()).collect();
    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        if parent[x] != x {
            let p = parent[x];
            let root = find(parent, p);
            parent[x] = root;
            root
        } else {
            x
        }
    }

    for i in 0..paths.len() {
        let path_i_len: usize = paths[i].iter().map(|&n| node_length[&n]).sum();
        for j in (i + 1)..paths.len() {
            let mut path_j_len = 0usize;
            let mut nodes_same_length = 0usize;
            for &node in paths[j].iter() {
                path_j_len += node_length[&node];
                if paths[i].contains(&node) {
                    nodes_same_length += node_length[&node];
                }
            }
            let iso_pct_overlap = if path_i_len == 0 || path_j_len == 0 {
                0.0
            } else {
                ((nodes_same_length as f32 / path_i_len as f32)
                    .max(nodes_same_length as f32 / path_j_len as f32))
                    * 100.0
            };
            if iso_pct_overlap >= params.min_isoform_pct_len_overlap {
                let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                if ri != rj {
                    parent[ri] = rj;
                }
            }
        }
    }

    // 分量编号：按成员最小下标出现序
    let mut root_to_cluster: BTreeMap<usize, usize> = BTreeMap::new();
    let mut gene_of_slot: Vec<usize> = vec![0; final_paths.len()];
    for i in 0..paths.len() {
        let r = find(&mut parent, i);
        let next = root_to_cluster.len() + 1;
        let cluster = *root_to_cluster.entry(r).or_insert(next);
        gene_of_slot[hm_order[i]] = cluster;
    }
    gene_of_slot
}

// ---------------------------------------------------------------------------
// printFinalPaths（L8954）+ get_pathName_string（L15663）
// ---------------------------------------------------------------------------

/// Java `get_pathName_string`（MISO 格式）：
/// `path=[id:j-k id:j-k ...]`（首节点全长、后续减 k-1），随后 `" " + path`
/// 追加 List.toString 形式的完整路径（含 -1/-2 sink）。
pub fn get_path_name_string(graph: &DiGraph, path: &[i32], kmer_size: usize) -> String {
    let mut start = 0usize;
    let mut end = path.len();
    if path.first() == Some(&VERTEX_ROOT_ID) {
        start = 1;
    }
    if path.last() == Some(&T_VERTEX_ID) {
        end -= 1;
    }
    let mut s = String::from("[");
    let mut j = 0usize;
    let n = end.saturating_sub(start);
    for (idx, &id) in path[start..end].iter().enumerate() {
        let name_len = graph.get_vertex(id).map(|v| v.name.len()).unwrap_or(0);
        let seg_len = if idx == 0 {
            name_len
        } else {
            name_len.saturating_sub(kmer_size - 1)
        };
        s.push_str(&format!("{id}:{j}-{}", j + seg_len.saturating_sub(1)));
        if idx + 1 < n {
            s.push(' ');
        }
        j += seg_len;
    }
    s.push(']');
    s
}

/// Java `printFinalPaths`：`>{name}_g{g}_i{i} len=L path=[...] [完整路径]`
/// + 60 列折行 FASTA。`paths` 必须按 Java HashMap 迭代序给出（isoform
///   编号取决于此序）；`gene_ids` 与 paths 对齐。
pub fn print_final_paths(
    graph: &DiGraph,
    paths: &[Vec<i32>],
    gene_ids: &[usize],
    comp_name: &str,
    kmer_size: usize,
) -> String {
    let name = comp_name.replace(".graph", "");
    let mut local_gene_id_mapping: FxHashMap<usize, usize> = FxHashMap::default();
    let mut local_seq_counter: FxHashMap<usize, usize> = FxHashMap::default();
    let mut gene_counter = 0usize;

    let mut out = String::new();
    for (path, &local_gene) in paths.iter().zip(gene_ids) {
        let seq = get_path_seq(graph, path, kmer_size);

        let gene_id;
        let iso_id;
        if let Some(&g) = local_gene_id_mapping.get(&local_gene) {
            gene_id = g;
            let next = local_seq_counter.get(&gene_id).copied().unwrap_or(0) + 1;
            local_seq_counter.insert(gene_id, next);
            iso_id = next;
        } else {
            gene_counter += 1;
            gene_id = gene_counter;
            local_gene_id_mapping.insert(local_gene, gene_id);
            local_seq_counter.insert(gene_id, 1);
            iso_id = 1;
        }

        let path_name = get_path_name_string(graph, path, kmer_size);
        out.push_str(&format!(
            ">{}_g{}_i{} len={} path={} {}\n",
            name,
            gene_id,
            iso_id,
            seq.len(),
            path_name,
            java_list_string(path)
        ));
        for chunk in seq.as_bytes().chunks(60) {
            out.push_str(std::str::from_utf8(chunk).unwrap());
            out.push('\n');
        }
    }
    out
}

// ---------------------------------------------------------------------------
// run_EM_REDUCE（L1556）+ PathExpressionComparator（EM 表达估计）
// ---------------------------------------------------------------------------

/// PathExpressionComparator 的核心：多轮 EM 估计各 transcript 的相对表达。
///
/// * `path_reads`：final path →（PairPath → read 计数）（orig id）
/// * `seq_lengths`：各 path 的序列长度
///
/// 全程 f32（Java float）。返回 path → 相对表达值（fractional expr）；
/// 未出现在 path_reads 中的 path 不在 EM 中（get_expr = 0，Java 同）。
pub fn run_path_expression_em(
    all_paths: &[Vec<i32>],
    path_reads: &ContainedReads,
    seq_lengths: &FxHashMap<Vec<i32>, usize>,
) -> FxHashMap<Vec<i32>, f32> {
    // pp → 支持它的 transcript 集；pp → read 支持数（putAll：后写覆盖）
    let mut pp_to_transcripts: FxHashMap<PairPath, Vec<Vec<i32>>> = FxHashMap::default();
    let mut pp_support: FxHashMap<PairPath, i64> = FxHashMap::default();
    for path in all_paths {
        let Some(pp_map) = path_reads.get(path) else {
            continue;
        };
        for (pp, &count) in pp_map {
            pp_support.insert(pp.clone(), count);
            pp_to_transcripts
                .entry(pp.clone())
                .or_default()
                .push(path.clone());
        }
    }

    let em_paths: Vec<Vec<i32>> = all_paths
        .iter()
        .filter(|p| path_reads.contains_key(*p))
        .cloned()
        .collect();
    let len_of = |p: &Vec<i32>| -> f32 { seq_lengths.get(p).copied().unwrap_or(1).max(1) as f32 };

    // 确定性迭代序（Java HashMap 序不可复现；求和顺序只影响 f32 末位）
    let mut pps: Vec<&PairPath> = pp_to_transcripts.keys().collect();
    pps.sort_by_key(|a| a.to_key_string());

    let mut sum_frags: FxHashMap<Vec<i32>, f32> = FxHashMap::default();
    let mut expr: FxHashMap<Vec<i32>, f32> = FxHashMap::default();
    let mut mass: FxHashMap<Vec<i32>, f32> = FxHashMap::default();
    for p in &em_paths {
        sum_frags.insert(p.clone(), 0.0);
    }

    // ---- init_transcript_expr ----
    for pp in &pps {
        let transcripts = &pp_to_transcripts[*pp];
        let n = transcripts.len() as f32;
        let support = *pp_support.get(*pp).unwrap_or(&0) as f32;
        let frac = support / n;
        for t in transcripts {
            *sum_frags.entry(t.clone()).or_default() += frac;
        }
    }
    let mut sum_expr: f32 = em_paths.iter().map(|p| sum_frags[p] / len_of(p)).sum();
    for p in &em_paths {
        let e = sum_frags[p] / len_of(p) / sum_expr;
        expr.insert(p.clone(), e);
    }
    let mut sum_mass: f32 = 0.0;
    for p in &em_paths {
        let m = expr[p] * len_of(p);
        mass.insert(p.clone(), m);
        sum_mass += m;
    }
    for p in &em_paths {
        *mass.get_mut(p).unwrap() /= sum_mass;
    }

    let likelihood = |_expr: &FxHashMap<Vec<i32>, f32>, mass: &FxHashMap<Vec<i32>, f32>| -> f32 {
        let mut ll = 0f32;
        for pp in &pps {
            let count = *pp_support.get(*pp).unwrap_or(&0) as f32;
            let mut sum_pp = 0f32;
            for t in &pp_to_transcripts[*pp] {
                sum_pp += mass.get(t).copied().unwrap_or(0.0) * (1.0 / len_of(t));
            }
            if sum_pp > 0.0 {
                ll += count * sum_pp.ln();
            }
        }
        ll
    };

    // ---- EM 主循环（|Δ|<0.01 或 100 轮）----
    let mut prev = likelihood(&expr, &mass);
    for _round in 0..100 {
        // E-step
        for p in &em_paths {
            sum_frags.insert(p.clone(), 0.0);
        }
        for pp in &pps {
            let count = *pp_support.get(*pp).unwrap_or(&0) as f32;
            let transcripts = &pp_to_transcripts[*pp];
            let mut sum_probs = 0f32;
            for t in transcripts {
                sum_probs += mass.get(t).copied().unwrap_or(0.0) * (1.0 / len_of(t));
            }
            if sum_probs == 0.0 {
                continue;
            }
            for t in transcripts {
                let prob = mass.get(t).copied().unwrap_or(0.0) * (1.0 / len_of(t));
                *sum_frags.entry(t.clone()).or_default() += count * (prob / sum_probs);
            }
        }
        // M-step
        sum_expr = em_paths.iter().map(|p| sum_frags[p] / len_of(p)).sum();
        for p in &em_paths {
            expr.insert(p.clone(), sum_frags[p] / len_of(p) / sum_expr);
        }
        sum_mass = 0.0;
        for p in &em_paths {
            let m = expr[p] * len_of(p);
            *mass.get_mut(p).unwrap() = m;
            sum_mass += m;
        }
        for p in &em_paths {
            *mass.get_mut(p).unwrap() /= sum_mass;
        }

        let curr = likelihood(&expr, &mass);
        let delta = curr - prev;
        prev = curr;
        if delta.abs() < 0.01 {
            break;
        }
    }

    expr
}

/// `remove_lesser_supported_paths_EM`（L10104）：按基因保留表达 ≥ 主导
/// isoform 5%（`MIN_RELATIVE_ISOFORM_EXPRESSION`；`MIN_TOTAL_ISOFORM_EXPRESSION`
/// 默认 0 = 关闭）的 isoform；每基因最高表达者恒保留。
fn remove_lesser_supported_paths_em(
    all_paths: &[Vec<i32>],
    expr: &FxHashMap<Vec<i32>, f32>,
    gene_ids: &[usize],
    params: &PostProcessParams,
) -> Vec<Vec<i32>> {
    // 每基因最大表达（严格 >，先到先得——Java 语义）
    let mut gene_max: FxHashMap<usize, f32> = FxHashMap::default();
    let mut gene_top: FxHashMap<usize, Vec<i32>> = FxHashMap::default();
    for (i, path) in all_paths.iter().enumerate() {
        let gene = gene_ids[i];
        let e = expr.get(path).copied().unwrap_or(0.0);
        match gene_max.get(&gene) {
            Some(&m) if m >= e => {}
            _ => {
                gene_max.insert(gene, e);
                gene_top.insert(gene, path.clone());
            }
        }
    }
    let mut keep = Vec::new();
    for (i, path) in all_paths.iter().enumerate() {
        let gene = gene_ids[i];
        let e = expr.get(path).copied().unwrap_or(0.0);
        let max = gene_max[&gene];
        let pct = e / max * 100.0;
        let is_top = gene_top.get(&gene) == Some(path);
        if is_top
            || (e * 100.0 >= params.min_total_isoform_expression
                && pct >= params.min_relative_isoform_expression)
        {
            keep.push(path.clone());
        }
    }
    keep
}

/// `run_EM_REDUCE`（L1556）：EM 表达估计 + 低表达 isoform 削减。
/// `no_em_reduce == true`（Trinity 主脚本 `--NO_EM_REDUCE`）时恒等返回。
pub fn run_em_reduce(
    graph: &DiGraph,
    final_paths: Vec<Vec<i32>>,
    path_reads: &ContainedReads,
    gene_ids: &[usize],
    kmer_size: usize,
    params: &PostProcessParams,
) -> Vec<Vec<i32>> {
    if params.no_em_reduce || final_paths.len() <= 1 {
        return final_paths;
    }
    let mut seq_lengths: FxHashMap<Vec<i32>, usize> = FxHashMap::default();
    for p in &final_paths {
        seq_lengths.insert(p.clone(), get_path_seq(graph, p, kmer_size).len());
    }
    let expr = run_path_expression_em(&final_paths, path_reads, &seq_lengths);
    remove_lesser_supported_paths_em(&final_paths, &expr, gene_ids, params)
}

// ---------------------------------------------------------------------------
// tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{DiGraph, SeqVertex, SimpleEdge};

    fn chain_graph() -> DiGraph {
        // 1(A*50) → 2(C*50) → 3(G*50)，kmer=25 时边拼接各 -24
        let mut g = DiGraph::new();
        g.add_vertex(SeqVertex::new(1, "A".repeat(50)));
        g.add_vertex(SeqVertex::new(2, "C".repeat(50)));
        g.add_vertex(SeqVertex::new(3, "G".repeat(50)));
        g.add_edge(1, 2, SimpleEdge::new(1.0, 1, 2));
        g.add_edge(2, 3, SimpleEdge::new(1.0, 2, 3));
        for v in [1, 2, 3] {
            g.get_vertex_mut(v).unwrap().dfs_finish_time = v;
        }
        g
    }

    #[test]
    fn putall_cap_order_matches_c0_golden_bucket_layout() {
        // c0 黄金 5 条（EM 后）：putAll 构建 → cap = tableSizeFor(5/0.75+1) = 8。
        // cap16 下 1103(bucket 4) 会排在 1287(bucket 9) 前；cap8 下 1287 与
        // 1034 同桶（1），随后 1078(3)/1103(4)/914(6)——恰为黄金 i1..i5 序。
        let paths: Vec<Vec<i32>> = vec![
            vec![1599, 3956, 3934],
            vec![1599, 5646, 3934],
            vec![-1, 2341, 4651, 1599],
            vec![3934, 3467, -2],
            vec![3934, 2342, -2],
        ];
        assert_eq!(java_putall_cap(5), 8);
        let order = java_hashmap_order_cap(&paths, java_putall_cap(5));
        assert_eq!(order, vec![0, 2, 3, 1, 4]); // 1034,1287,1078,1103,914
    }

    #[test]
    fn java_hash_order_matches_computed_buckets() {
        // 黄金 c0 的 5 条路径（orig id，含 sink）——手算 Java HashMap 桶序
        // 应恰为 allprobpaths.orig.fasta 的 i1..i5 顺序。
        let paths: Vec<Vec<i32>> = vec![
            vec![1599, 3956, 3934],
            vec![-1, 2341, 4651, 1599],
            vec![3934, 3467, -2],
            vec![1599, 5646, 3934],
            vec![3934, 2342, -2],
        ];
        let order = java_hashmap_order(&paths);
        // cap16 默认桶序：1034(1) 1103(4) 1287(9) 1078(11) 914(14)——与黄金
        // 打印序（putAll cap8，见下一测试）不同，属预期。
        assert_eq!(order, vec![0, 3, 1, 2, 4]);
    }

    #[test]
    fn identical_paths_are_too_similar() {
        let g = chain_graph();
        let params = PostProcessParams::default();
        let mut cache = FxHashMap::default();
        let p = vec![1, 2, 3];
        assert!(two_paths_are_too_similar(
            &g, &p, &p, 25, &params, &mut cache
        ));
    }

    #[test]
    fn disjoint_paths_are_not_too_similar() {
        let g = chain_graph();
        let mut g = g;
        g.add_vertex(SeqVertex::new(4, "T".repeat(50)));
        g.add_vertex(SeqVertex::new(5, "A".repeat(50)));
        g.add_edge(4, 5, SimpleEdge::new(1.0, 4, 5));
        let params = PostProcessParams::default();
        let mut cache = FxHashMap::default();
        assert!(!two_paths_are_too_similar(
            &g,
            &[1, 2, 3],
            &[4, 5],
            25,
            &params,
            &mut cache
        ));
    }

    #[test]
    fn find_last_shared_node_basic() {
        let g = chain_graph();
        // QUIRK（Java 同）：[1,2] vs [1,2,3] 双指针从尾部后退——p1 的 2 先于
        // p2 的 3（finish time 序），p1 率先耗尽 → -1（"最后共享节点"并不可靠）
        assert_eq!(find_last_shared_node(&g, &[1, 2], &[1, 2, 3]), -1);
        // 尾节点相同 → 直接命中
        assert_eq!(find_last_shared_node(&g, &[1, 2, 3], &[1, 4, 3]), 3);
        // 无共享
        assert_eq!(find_last_shared_node(&g, &[1], &[3]), -1);
    }

    #[test]
    fn empty_side_counts_full_gap() {
        // 单侧空：max_internal_gap = seqLen-(k-1)
        let g = chain_graph();
        let mut cache = FxHashMap::default();
        let ctx = AlignCtx {
            kmer_size: 25,
            max_dp_len: 10000,
        };
        let st = get_prev_calc_num_mismatches(&g, &[], &[1, 2, 3], ctx, &mut cache);
        // seq([1,2,3]) = 50 + 26 + 26 = 102 → gap = 102 - (25-1) = 78
        assert_eq!(st.max_internal_gap_length, 78);
        assert_eq!(st.gaps, 78);
    }

    #[test]
    fn three_segment_decomposition_is_exact_for_shared_tail() {
        // [1,2,3] vs [1,2,3]（共享尾）→ 分解为完美匹配
        let g = chain_graph();
        let mut cache = FxHashMap::default();
        let ctx = AlignCtx {
            kmer_size: 25,
            max_dp_len: 10000,
        };
        let st = get_prev_calc_num_mismatches(&g, &[1, 2, 3], &[1, 2, 3], ctx, &mut cache);
        assert_eq!(st.matches, 102); // 50 + 26 + 26
        assert_eq!(st.mismatches, 0);
        // [-1,1,2,3] vs [1,2,3]：prefix（[-1] vs []）触发单侧空分支——
        // Java 的 gap = 0-(k-1) = -24（负数！≤10 成立）；我们 saturate 到 0，
        // 判定等价（同样"gap 小"）
        let st2 = get_prev_calc_num_mismatches(&g, &[-1, 1, 2, 3], &[1, 2, 3], ctx, &mut cache);
        assert_eq!(st2.max_internal_gap_length, 0);
    }

    #[test]
    fn reduce_removes_lesser_supported_near_duplicate() {
        // 两条仅 1bp 差异的路径：高支持者保留
        let mut g = DiGraph::new();
        let mut name2 = "C".repeat(50);
        name2.replace_range(10..11, "G"); // 与 path A 差 1bp
        g.add_vertex(SeqVertex::new(1, "A".repeat(50)));
        g.add_vertex(SeqVertex::new(2, name2));
        g.add_vertex(SeqVertex::new(3, "G".repeat(50)));
        g.add_vertex(SeqVertex::new(4, "C".repeat(50)));
        g.add_edge(1, 2, SimpleEdge::new(1.0, 1, 2));
        g.add_edge(2, 3, SimpleEdge::new(1.0, 2, 3));
        g.add_edge(1, 4, SimpleEdge::new(1.0, 1, 4));
        g.add_edge(4, 3, SimpleEdge::new(1.0, 4, 3));
        for v in [1, 2, 3, 4] {
            g.get_vertex_mut(v).unwrap().dfs_finish_time = v;
        }
        let params = PostProcessParams::default();
        let mut cache = FxHashMap::default();
        let mut path_reads: ContainedReads = FxHashMap::default();
        path_reads.insert(vec![1, 2, 3], FxHashMap::default());
        path_reads
            .get_mut(&vec![1, 2, 3])
            .unwrap()
            .insert(PairPath::new(vec![1, 2, 3]), 10);
        let kept = reduce_cdhit_like(
            &g,
            vec![vec![1, 2, 3], vec![1, 4, 3]],
            &mut path_reads,
            25,
            &params,
            &mut cache,
        );
        assert_eq!(kept, vec![vec![1, 2, 3]]);
        assert!(!path_reads.contains_key(&vec![1, 4, 3]));
    }

    #[test]
    fn group_genes_by_overlap() {
        let mut g = DiGraph::new();
        g.add_vertex(SeqVertex::new(1, "A".repeat(100)));
        g.add_vertex(SeqVertex::new(2, "C".repeat(100)));
        g.add_vertex(SeqVertex::new(3, "G".repeat(100)));
        g.add_vertex(SeqVertex::new(9, "T".repeat(100)));
        let paths = vec![
            vec![1, 2, 3], // 与下条共享 1+2（200/300 ≥30%）
            vec![1, 2],
            vec![9], // 完全独立
        ];
        let genes = group_paths_into_genes(&g, &paths, &PostProcessParams::default());
        assert_eq!(genes[0], genes[1]);
        assert_ne!(genes[0], genes[2]);
    }

    #[test]
    fn miso_path_string_format() {
        let g = chain_graph();
        // c0 黄金 i1 形态：首节点全长、后续 -24
        let s = get_path_name_string(&g, &[1, 2, 3], 25);
        assert_eq!(s, "[1:0-49 2:50-75 3:76-101]");
        // sink 修剪
        let s2 = get_path_name_string(&g, &[-1, 1, 2, 3, -2], 25);
        assert_eq!(s2, "[1:0-49 2:50-75 3:76-101]");
    }

    #[test]
    fn em_expression_and_removal() {
        // 两条共享读路径的 transcript：t1 专有读多 → 高表达保留；
        // t2 共享读为主且短 → 低表达（<5%）被删。
        let t1 = vec![1, 2, 3];
        let t2 = vec![1, 4, 3];
        let mut reads: ContainedReads = FxHashMap::default();
        let mut m1 = FxHashMap::default();
        m1.insert(PairPath::new(vec![1, 2]), 100);
        m1.insert(PairPath::new(vec![2, 3]), 100);
        reads.insert(t1.clone(), m1);
        let mut m2 = FxHashMap::default();
        m2.insert(PairPath::new(vec![1, 4]), 1);
        m2.insert(PairPath::new(vec![4, 3]), 1);
        reads.insert(t2.clone(), m2);
        let mut lens = FxHashMap::default();
        lens.insert(t1.clone(), 100);
        lens.insert(t2.clone(), 100);
        let expr = run_path_expression_em(&[t1.clone(), t2.clone()], &reads, &lens);
        assert!(expr[&t1] > 0.9);
        assert!(expr[&t2] < 0.1);
        let params = PostProcessParams::default();
        let kept =
            remove_lesser_supported_paths_em(&[t1.clone(), t2.clone()], &expr, &[1, 1], &params);
        assert_eq!(kept, vec![t1.clone()]);
        // no_em_reduce：恒等
        let g = DiGraph::new();
        let kept_all = run_em_reduce(
            &g,
            vec![t1.clone(), t2.clone()],
            &reads,
            &[1, 1],
            25,
            &PostProcessParams {
                no_em_reduce: true,
                ..Default::default()
            },
        );
        assert_eq!(kept_all.len(), 2);
    }

    #[test]
    fn fasta_line_wrap_60() {
        let mut g = DiGraph::new();
        g.add_vertex(SeqVertex::new(1, "ACGT".repeat(20))); // 80bp
        let paths = vec![vec![1]];
        let genes = vec![1usize];
        let out = print_final_paths(&g, &paths, &genes, "c0", 25);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines[0], ">c0_g1_i1 len=80 path=[1:0-79] [1]");
        assert_eq!(lines[1].len(), 60);
        assert_eq!(lines[2].len(), 20);
    }
}
