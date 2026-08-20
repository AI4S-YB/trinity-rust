//! read 穿线全链（TransAssembly_allProbPaths.java L11766-12687）：
//! `getReadStarts` → `readAndMapSingleRead` → `findPathInGraph` →
//! `updatePathRecursively`（递归 + memoization + zipper/短NW预检/带状NW 三级比对切换）。
//!
//! `loc_in_node` 精读结论（`getOriginalVerIDsMappingHash` + 递归调用点）：
//! 它是 **verSeq（节点全长 name，非 kmerAdj）的 0-based 碱基下标**，zipper 的 i 直接
//! 从它开始。取值来源有三种：
//! 1. 原始节点 id 自身（vid ≤ LAST_REAL_ID）→ `LocInGraph(vid, 0)`：read 的锚 kmer
//!    与 name 的前 K-1 重叠区对齐，即从 name 下标 0 起比对；
//! 2. prevID 段 → `LocInGraph(合并节点id, loc)`，loc 是 prevID 段计数（首个段 loc=1，
//!    若合并节点 vid > LAST_REAL_ID 则自身不注册、首个段 loc=0）——**quirk：段计数被
//!    直接当作碱基下标用**（Java 原文如此，段序号≠真实碱基偏移，但保守地落在 name
//!    内部即被当作 zipper 起点）；
//! 3. 后继递归 → `KMER_SIZE-1`：read 消费完前驱末碱基后，与前驱共享 K-1 重叠区的
//!    后继 name 的下标 K-1 处继续。

use rustc_hash::FxHashMap;

use crate::align::{
    get_num_left_end_gaps, get_num_right_end_gaps, run_nw_alignment, run_nw_banded_alignment,
};
use crate::context::BflyContext;
use crate::graph::DiGraph;
use crate::graph_io::RawRead;

/// Java `LocInGraph`：图中位置 = 节点 id + 节点内下标（见模块文档的 loc_in_node 语义）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocInGraph {
    pub node_id: i32,
    pub index_in_node: i32,
}

impl LocInGraph {
    pub fn new(node_id: i32, index_in_node: i32) -> Self {
        Self {
            node_id,
            index_in_node,
        }
    }
}

/// Java `Path_n_MM_count`。
///
/// `path`/`positions` 由 `add` **头插**构建（从右向左），最终为左→右序。
/// `positions` 元素 = (vertex_id, vertex_base_end i, read_base_end j)。
#[derive(Debug, Clone, PartialEq)]
pub struct ReadPath {
    pub mismatch_count: i32,
    pub path: Vec<i32>,
    pub positions: Vec<(i32, i32, i32)>,
}

impl ReadPath {
    /// 基例：单节点路径（Java `Path_n_MM_count(node, mm, i, j)`）。
    fn from_node_end(node: i32, mm: i32, end_vertex_base: i32, end_read_base: i32) -> Self {
        Self {
            mismatch_count: mm,
            path: vec![node],
            positions: vec![(node, end_vertex_base, end_read_base)],
        }
    }

    /// Java `add_path_n_mm`：头插当前节点并累加 mm。
    fn add_node_front(&mut self, node: i32, mm: i32, end_vertex_base: i32, end_read_base: i32) {
        self.path.insert(0, node);
        self.mismatch_count += mm;
        self.positions
            .insert(0, (node, end_vertex_base, end_read_base));
    }

    /// Java `get_trimmed_path(min_end_length)`：
    /// 前端 `read_base_end < min` 逐段删除；末端 `vertex_base_end < min` 逐段删除。
    pub fn get_trimmed_path(&self, min_end_length: i32) -> Vec<i32> {
        let mut vpps: Vec<(i32, i32, i32)> = self.positions.clone();
        while let Some(&(_, _, read_base_end)) = vpps.first() {
            if read_base_end < min_end_length {
                vpps.remove(0);
            } else {
                break;
            }
        }
        while let Some(&(_, vertex_base_end, _)) = vpps.last() {
            if vertex_base_end < min_end_length {
                vpps.pop();
            } else {
                break;
            }
        }
        vpps.into_iter().map(|(vid, _, _)| vid).collect()
    }
}

/// 穿线上下文（graph 只读；`max_mm_allowed` 与 `memo` 均**每 read 重置**——
/// Java 中它们是静态全局，靠 findPathInGraph 每次重赋值）。
pub struct ThreadingCtx<'a> {
    pub graph: &'a DiGraph,
    pub kmer_size: usize,
    pub max_read_seq_divergence: f64,
    pub max_read_local_seq_divergence: f64,
    pub use_dp: bool,
    /// MAX_MM_ALLOWED = ceil(seq.len * MAX_READ_SEQ_DIVERGENCE)，每 read 重置。
    max_mm_allowed: i32,
    /// best_path_memoization（value None = 死端缓存，Java 存 null）。
    memo: FxHashMap<String, Option<ReadPath>>,
}

impl<'a> ThreadingCtx<'a> {
    pub fn new(graph: &'a DiGraph, ctx: &BflyContext) -> Self {
        Self {
            graph,
            kmer_size: ctx.kmer_size,
            max_read_seq_divergence: ctx.max_read_seq_divergence,
            max_read_local_seq_divergence: ctx.max_read_local_seq_divergence,
            use_dp: ctx.use_dp_read_to_vertex_align,
            max_mm_allowed: 0,
            memo: FxHashMap::default(),
        }
    }
}

/// `areTwoNucleotidesEqual`（USE_DEGENERATE_CODE 默认 false → 纯等值比较）。
/// 注意与比对器打分不同：zipper **不**把 N 当匹配（quirk 保留）。
fn are_two_nucleotides_equal(a: u8, b: u8) -> bool {
    a == b
}

/// `getOriginalVerIDsMappingHash`（L12647-12687）。
///
/// 先对每个节点 `clearDoubleEntriesToPrevIDs`（相邻完全相等 prevID 段去重）；
/// `vid > LAST_REAL_ID` 不注册自身（loc 从 -1 起，首个 prev 段 loc=0）；
/// 否则自身注册为 (vid, 0)，随后每个 prevID 段 loc++ 逐 id 注册
/// (合并节点id, loc)。后写的条目覆盖先写的。
pub fn get_original_ver_ids_mapping_hash(
    graph: &mut DiGraph,
    last_real_id: i32,
) -> FxHashMap<i32, LocInGraph> {
    let ids: Vec<i32> = graph.vertex_ids().to_vec();
    for vid in ids {
        graph
            .get_vertex_mut(vid)
            .expect("in vertex_ids")
            .clear_double_entries_to_prev_ids();
    }
    let mut hash: FxHashMap<i32, LocInGraph> = FxHashMap::default();
    let ids: Vec<i32> = graph.vertex_ids().to_vec();
    for vid in ids {
        let v = graph.get_vertex(vid).expect("in vertex_ids");
        let mut loc: i32 = 0;
        if vid > last_real_id {
            loc = -1;
        } else {
            hash.insert(vid, LocInGraph::new(vid, 0));
        }
        for vec in &v.prev_vertices_id {
            loc += 1;
            for &id in vec {
                hash.insert(id, LocInGraph::new(vid, loc));
            }
        }
    }
    hash
}

/// `getReadStarts`（L11766）。
///
/// 输入为已解析的 raw reads（组件 reads 文件去头）。返回 read 名 → 该名下所有
/// 成功穿线的路径（顺序 = 文件出现序；Read 的 path 存**未 trim** 的 best_path.path，
/// 镜像 Java：Read 对象在 trim 前入 hash）。
pub fn get_read_starts(
    graph: &DiGraph,
    ctx: &BflyContext,
    orig_id_map: &FxHashMap<i32, LocInGraph>,
    orig_kmer_to_node: &FxHashMap<String, i32>,
    raw_reads: &[RawRead],
) -> FxHashMap<String, Vec<ReadPath>> {
    let mut read_name_hash: FxHashMap<String, Vec<ReadPath>> = FxHashMap::default();
    for (name, p) in get_read_starts_ordered(graph, ctx, orig_id_map, orig_kmer_to_node, raw_reads)
    {
        read_name_hash.entry(name).or_default().push(p);
    }
    read_name_hash
}

/// 同 `get_read_starts`，但保留**文件出现序**的 `(read 名, 路径)` 列表。
/// Java readNameHash 是名 → List<Read>，物理 occurrence（/1 vs /2 行）谁成功
/// 穿线只有顺序信息能区分——debug dump / 后续层自回归用。
pub fn get_read_starts_ordered(
    graph: &DiGraph,
    ctx: &BflyContext,
    orig_id_map: &FxHashMap<i32, LocInGraph>,
    orig_kmer_to_node: &FxHashMap<String, i32>,
    raw_reads: &[RawRead],
) -> Vec<(String, ReadPath)> {
    let mut result: Vec<(String, ReadPath)> = Vec::new();
    let k = ctx.kmer_size;

    for raw in raw_reads {
        // readAndMapSingleRead 字段解析
        let mut seq: Vec<u8> = raw.seq.clone().into_bytes();
        let end_in_read = raw.end_in_read;

        if end_in_read >= seq.len() as i64 {
            continue; // Java: 长度短于 endInRead 标记 → 返回空 path
        }
        let start = raw.start_in_read as usize;
        let end = end_in_read as usize;
        // Java substring(start, end+1)（endInRead = f3 + K 的 off-by-one 保留）
        seq = seq[start..=end].to_vec();

        // 锚点失效再锚定：滑窗 i 从 2 起，找第一个仍存活的 kmer 并截前缀
        let mut from_v: Option<&LocInGraph> = orig_id_map.get(&raw.from_orig_v);
        if from_v.is_none() && seq.len() >= k {
            for i in 2..=(seq.len() - k) {
                let kmer = String::from_utf8_lossy(&seq[i..i + k]).into_owned();
                if let Some(&id) = orig_kmer_to_node.get(&kmer) {
                    if let Some(loc) = orig_id_map.get(&id) {
                        from_v = Some(loc);
                        seq = seq[i..].to_vec();
                        break;
                    }
                }
            }
        }

        let Some(from_v) = from_v else {
            continue;
        };

        // findPathInGraph
        let mut tctx = ThreadingCtx::new(graph, ctx);
        if let Some(best_path) = find_path_in_graph(&mut tctx, &seq, from_v) {
            result.push((raw.name.clone(), best_path));
        }
    }
    result
}

/// `findPathInGraph`（L11992）：每 read 重置 MAX_MM_ALLOWED 与 memo。
pub fn find_path_in_graph(
    ctx: &mut ThreadingCtx,
    seq: &[u8],
    from_v: &LocInGraph,
) -> Option<ReadPath> {
    // MAX_MM_ALLOWED_CAP = (int) Math.ceil(seq.length() * MAX_READ_SEQ_DIVERGENCE)
    ctx.max_mm_allowed = (seq.len() as f64 * ctx.max_read_seq_divergence).ceil() as i32;
    ctx.memo = FxHashMap::default();
    let total_num_mm = 0;
    update_path_recursively(
        ctx,
        from_v.node_id,
        seq,
        0,
        from_v.index_in_node,
        total_num_mm,
    )
}

/// `updatePathRecursively`（L12053-12466）。返回深拷贝（Java 逐字镜像）。
fn update_path_recursively(
    ctx: &mut ThreadingCtx,
    from_v_id: i32,
    seq: &[u8],
    loc_in_seq: i32,
    loc_in_node: i32,
    total_num_mm: i32,
) -> Option<ReadPath> {
    const MIN_SEQ_LENGTH_TEST_DIVERGENCE: i32 = 20;
    const MAX_LEFT_END_GAPS: i32 = 5;
    const MIN_LENGTH_TEST_DP: usize = 100;

    let ver_seq: &[u8] = ctx.graph.get_vertex(from_v_id)?.name.as_bytes();
    let graph = ctx.graph;
    let kmer_size = ctx.kmer_size;

    let mut num_mm = total_num_mm;

    let start_i = loc_in_node;
    let mut j = loc_in_seq;
    let mut i = start_i;

    let token = format!("{from_v_id}_{loc_in_node}_{loc_in_seq}");
    if let Some(cached) = ctx.memo.get(&token) {
        // 命中：None = 死端；Some = 返回深拷贝
        return cached.clone();
    }

    let length_to_align = ((ver_seq.len() as i32 - i).min(seq.len() as i32 - j)).max(0) as usize;

    // ---- zipper align ----
    let mut failed_alignment = false;
    let mut mm_encountered_here: i32 = 0;
    while i >= 0 && (i as usize) < ver_seq.len() && (j as usize) < seq.len() {
        let read_letter = seq[j as usize];
        let ver_letter = ver_seq[i as usize];
        if !are_two_nucleotides_equal(read_letter, ver_letter) {
            num_mm += 1;
            mm_encountered_here += 1;
            if num_mm > ctx.max_mm_allowed
                || (i >= MIN_SEQ_LENGTH_TEST_DIVERGENCE
                    && (mm_encountered_here as f32 / i as f32)
                        > ctx.max_read_local_seq_divergence as f32)
            {
                failed_alignment = true;
                break;
            }
        }
        i += 1;
        j += 1;
    }

    let zipper_i = i;
    let zipper_j = j;
    let zipper_mm = mm_encountered_here;

    // ---- 100bp 短 NW 预检 ----
    let mut short_dp_test_passes = true;
    if ctx.use_dp && length_to_align > MIN_LENGTH_TEST_DP && mm_encountered_here > 1 {
        // 重置 i/j 做短 NW；无论结果如何随后恢复 zipper 状态
        let (_aln, stats) = run_nw_alignment(
            &ver_seq[start_i as usize..start_i as usize + MIN_LENGTH_TEST_DP],
            &seq[loc_in_seq as usize..loc_in_seq as usize + MIN_LENGTH_TEST_DP],
            4.0,
            -5.0,
            10.0,
            1.0,
        );
        mm_encountered_here = (stats.mismatches + stats.gaps + stats.left_gap_length) as i32;
        let pct_divergence = mm_encountered_here as f32 / MIN_LENGTH_TEST_DP as f32;
        if pct_divergence as f64 > ctx.max_read_local_seq_divergence {
            short_dp_test_passes = false;
            // failed_alignment 状态保持
        }
        // 恢复 zipper 统计
        i = zipper_i;
        j = zipper_j;
        mm_encountered_here = zipper_mm;
    }

    // ---- banded NW ----
    let mut max_left_gaps: i32 = 0;

    if ctx.use_dp && ver_seq.len() > 2 && mm_encountered_here > 1 && short_dp_test_passes {
        j = loc_in_seq;
        i = start_i;

        // read_length_to_align = (int)(verSeq.length() * 1.05f)
        let mut read_length_to_align = (ver_seq.len() as f32 * 1.05f32) as i32; // Java int 截断
        if read_length_to_align + j > seq.len() as i32 {
            read_length_to_align = seq.len() as i32 - j;
        }
        let bandwidth =
            (ctx.max_read_local_seq_divergence as f32 * read_length_to_align as f32) as i32;

        let ver_part = &ver_seq[i as usize..];
        let read_part = &seq[j as usize..(j + read_length_to_align) as usize];
        let (aln, stats) = run_nw_banded_alignment(
            ver_part,
            read_part,
            4.0,
            -5.0,
            10.0,
            1.0,
            bandwidth.max(0) as usize,
        );

        mm_encountered_here =
            (stats.mismatches + stats.gaps + stats.left_gap_length + stats.right_gap_length) as i32;

        // Java：name1=="Read" 时交换，使 aligned1 = vertex 行。
        // 我们的端口约定 aligned1 对应较长的输入（等长时对 s1=ver_part）。
        let (vertex_align, read_align): (&[u8], &[u8]) = if ver_part.len() >= read_part.len() {
            (&aln.aligned1, &aln.aligned2)
        } else {
            (&aln.aligned2, &aln.aligned1)
        };

        let vertex_num_right_end_gaps = get_num_right_end_gaps(vertex_align) as i32;
        let read_num_right_end_gaps = get_num_right_end_gaps(read_align) as i32;
        max_left_gaps = (get_num_left_end_gaps(vertex_align) as i32)
            .max(get_num_left_end_gaps(read_align) as i32);

        i = ver_seq.len() as i32;
        j += read_length_to_align;

        if vertex_num_right_end_gaps > 0 {
            // read 超出 vertex 末端：读回退
            j -= vertex_num_right_end_gaps;
            mm_encountered_here -= vertex_num_right_end_gaps;
        } else if read_num_right_end_gaps > 0 {
            // vertex 超出 read 末端
            mm_encountered_here -= read_num_right_end_gaps;
        }

        if mm_encountered_here >= zipper_mm && zipper_i == ver_seq.len() as i32 {
            // 回退 zipper
            i = zipper_i;
            j = zipper_j;
            mm_encountered_here = zipper_mm;
            max_left_gaps = 0;
            // failed_alignment 状态保持
        } else {
            failed_alignment = false;
        }

        num_mm = mm_encountered_here + total_num_mm;
    }

    // ---- 失败判定 ----
    let current_alignment_divergence = num_mm as f32 / j as f32;
    let local_vertex_alignment_divergence = mm_encountered_here as f32 / i as f32;

    if i >= MIN_SEQ_LENGTH_TEST_DIVERGENCE
        && local_vertex_alignment_divergence as f64 >= ctx.max_read_local_seq_divergence
    {
        failed_alignment = true;
    }
    if max_left_gaps > MAX_LEFT_END_GAPS {
        failed_alignment = true;
    }

    let successors_empty = graph.get_successors(from_v_id).is_empty();

    if (current_alignment_divergence as f64 > ctx.max_read_seq_divergence
        || num_mm > ctx.max_mm_allowed)
        || failed_alignment
    {
        if failed_alignment {
            ctx.memo.insert(token, None);
        }
        None // go back and try alternative vertex if available
    } else if j == seq.len() as i32 || successors_empty {
        // 基例：read 结束（或无后继——未对齐末端计为 mismatch）
        if successors_empty {
            mm_encountered_here += seq.len() as i32 - j;
        }
        let best_path = ReadPath::from_node_end(from_v_id, mm_encountered_here, i, j);
        ctx.memo.insert(token, Some(best_path.clone()));
        Some(best_path)
    } else if i == ver_seq.len() as i32 {
        // 到达节点末端：对每个后继递归（loc_in_node = KMER_SIZE-1）
        let continue_vers_ids: Vec<i32> = graph.get_successors(from_v_id).to_vec();

        let mut best_path: Option<ReadPath> = None;

        for &successor_vertex_id in &continue_vers_ids {
            let best_extension = update_path_recursively(
                ctx,
                successor_vertex_id,
                seq,
                j,
                kmer_size as i32 - 1,
                num_mm,
            );

            if let Some(best_extension) = best_extension {
                // best==null || ext.mm <= best.mm（<= 平局替换——取迭代序最后一条）
                if best_path.is_none()
                    || best_extension.mismatch_count
                        <= best_path.as_ref().expect("checked").mismatch_count
                {
                    // Java 此处是引用赋值（平局后 add 直接改子对象；非平局先拷贝再
                    // add）。我们独占所有权，两者等价——统一原地头插。tie 的区分在
                    // 结果上无差异（TRUNCATE_TIED_PATH 恒为 false）。
                    best_path = Some(best_extension);
                }
            }
        }

        if let Some(mut best) = best_path {
            best.add_node_front(from_v_id, mm_encountered_here, i, j);
            ctx.memo.insert(token, Some(best.clone()));
            Some(best)
        } else {
            ctx.memo.insert(token, None);
            None
        }
    } else {
        // Java: throw RuntimeException("should never end up here ...")
        panic!(
            "should never end up here, supposedly.  i={i}, j={j} ver length = {}  and readSeq length = {} ",
            ver_seq.len(),
            seq.len()
        );
    }
}

/// debug dump：read 名 → path 节点序列（T8 自回归用）。
/// 输入为 `get_read_starts_ordered` 的结果；每行 `name\tv1,v2,...`（文件出现序）。
pub fn format_read_paths_dump(ordered: &[(String, ReadPath)]) -> String {
    let mut out = String::new();
    for (name, p) in ordered {
        let ids: Vec<String> = p.path.iter().map(|v| v.to_string()).collect();
        out.push_str(name);
        out.push('\t');
        out.push_str(&ids.join(","));
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::BflyContext;
    use crate::graph::{DiGraph, SeqVertex, SimpleEdge};

    fn chain_graph() -> DiGraph {
        // K=4 三节点链：v2/v3 的 name 前 K-1 与前驱末 K-1 重叠
        let mut g = DiGraph::new();
        g.add_vertex(SeqVertex::new(1, "ACGTAAAA"));
        g.add_vertex(SeqVertex::new(2, "AAACCCC"));
        g.add_vertex(SeqVertex::new(3, "CCCTTTT"));
        g.add_edge(1, 2, SimpleEdge::new(1.0, 1, 2));
        g.add_edge(2, 3, SimpleEdge::new(1.0, 2, 3));
        g
    }

    fn ctx_k4() -> BflyContext {
        let mut c = BflyContext::new();
        c.kmer_size = 4;
        c
    }

    fn map1() -> FxHashMap<i32, LocInGraph> {
        let mut m = FxHashMap::default();
        m.insert(1, LocInGraph::new(1, 0));
        m
    }

    #[test]
    fn perfect_match_three_node_chain() {
        let g = chain_graph();
        let ctx = ctx_k4();
        let seq = b"ACGTAAAACCCCTTTT";
        let mut t = ThreadingCtx::new(&g, &ctx);
        let p = find_path_in_graph(&mut t, seq, &LocInGraph::new(1, 0)).unwrap();
        assert_eq!(p.path, vec![1, 2, 3]);
        assert_eq!(p.mismatch_count, 0);
        assert_eq!(p.positions, vec![(1, 8, 8), (2, 7, 12), (3, 7, 16)]);
        // max_mm = ceil(16*0.05) = ceil(0.8) = 1
        assert_eq!(t.max_mm_allowed, 1);
    }

    #[test]
    fn single_mismatch_within_thresholds() {
        // 24bp read（ceil(24*0.05)=2 允许 2 个 mm；1/24=0.042 < 0.05）
        let mut g = DiGraph::new();
        g.add_vertex(SeqVertex::new(1, "ACGTAAAACCCC"));
        g.add_vertex(SeqVertex::new(2, "CCCTTTTGGGGAAAA"));
        g.add_edge(1, 2, SimpleEdge::new(1.0, 1, 2));
        let ctx = ctx_k4();
        // 读序列与路径一致但末位 1 个 mismatch
        let seq = b"ACGTAAAACCCCTTTTGGGGAAAX";
        let mut t = ThreadingCtx::new(&g, &ctx);
        assert_eq!(t.max_mm_allowed, 0); // 未跑前为 0
        let p = find_path_in_graph(&mut t, seq, &LocInGraph::new(1, 0)).unwrap();
        assert_eq!(p.path, vec![1, 2]);
        assert_eq!(p.mismatch_count, 1);
        assert_eq!(t.max_mm_allowed, 2);
    }

    #[test]
    fn f32_quirk_one_mismatch_of_twenty_fails() {
        // Java quirk 保留：numMM/(float)j 在 j=20、mm=1 时 = 0.05f32，
        // 而 MAX_READ_SEQ_DIVERGENCE 是 double 0.05——0.05f32 > 0.05 成立 → 失败。
        let mut g = DiGraph::new();
        g.add_vertex(SeqVertex::new(1, "ACGTAAAACCCC"));
        g.add_vertex(SeqVertex::new(2, "CCCTTTTGGGG"));
        g.add_edge(1, 2, SimpleEdge::new(1.0, 1, 2));
        let ctx = ctx_k4();
        let seq = b"ACGTAAAACCCCTTTTGGGX";
        let mut t = ThreadingCtx::new(&g, &ctx);
        assert!(find_path_in_graph(&mut t, seq, &LocInGraph::new(1, 0)).is_none());
    }

    #[test]
    fn divergence_failure_memoizes_dead_end() {
        // 短 read 1 mismatch：numMM/j = 1/16 > 0.05 → 失败且 memo 死端
        let g = chain_graph();
        let ctx = ctx_k4();
        let mut seq = b"ACGTAAAACCCCTTTT".to_vec();
        seq[15] = b'G'; // 1 mismatch，j=16 时 div=0.0625>0.05
        let mut t = ThreadingCtx::new(&g, &ctx);
        assert!(find_path_in_graph(&mut t, &seq, &LocInGraph::new(1, 0)).is_none());
        // 死端已缓存（vertex1@0_0）
        assert!(t.memo.contains_key("1_0_0"));
        assert!(t.memo.get("1_0_0").unwrap().is_none());
    }

    #[test]
    fn tie_breaks_to_last_successor() {
        // 平局：两个后继 mm 相等 → 取迭代序（插入序）最后一条
        let mut g = DiGraph::new();
        g.add_vertex(SeqVertex::new(1, "ACGTAAAA"));
        // 两个分支 name 完全一致（read 后缀走 4bp 后无后继，各计 8 个尾部 mm → 平局）
        g.add_vertex(SeqVertex::new(2, "AAACCCC"));
        g.add_vertex(SeqVertex::new(3, "AAACCCC"));
        g.add_edge(1, 2, SimpleEdge::new(1.0, 1, 2));
        g.add_edge(1, 3, SimpleEdge::new(1.0, 1, 3));
        let ctx = ctx_k4();
        let seq = b"ACGTAAAACCCCTTTTTTTT";
        let mut t = ThreadingCtx::new(&g, &ctx);
        let p = find_path_in_graph(&mut t, seq, &LocInGraph::new(1, 0)).unwrap();
        assert_eq!(p.path, vec![1, 3]);
        // 两个分支尾部未对齐 8bp 均在基例中累加（失败判定发生在累加前，Java 语义）
        assert_eq!(p.mismatch_count, 8);
    }

    #[test]
    fn zipper_local_divergence_shortcircuit() {
        // i>=20 后 mm/i > 0.1 短路 + 超 max_mm：i=21..23 密集 mismatch
        let mut g = DiGraph::new();
        let name = "ACGTAAAACCCCTTTTGGGGAAAAAAACCC".to_string(); // 30bp
        g.add_vertex(SeqVertex::new(1, name.clone()));
        let ctx = ctx_k4();
        let mut seq = name.clone().into_bytes();
        seq[21] = b'X';
        seq[22] = b'X';
        seq[23] = b'X';
        // len 30，max_mm=ceil(1.5)=2；i=23 时 mm=3：3 > 2 → zipper 短路 failed
        let mut t = ThreadingCtx::new(&g, &ctx);
        assert!(find_path_in_graph(&mut t, &seq, &LocInGraph::new(1, 0)).is_none());
        assert!(t.memo.get("1_0_0").unwrap().is_none()); // 死端缓存
    }

    #[test]
    fn banded_nw_triggered_by_two_mismatches() {
        // 150bp 节点，3 个稀疏 mismatch（60/80/100）：zipper mm=3 > 1 →
        // 触发 100bp 短 NW 预检（divergence 小，通过）→ banded NW 全长比对
        let mut g = DiGraph::new();
        let mut name: Vec<u8> = Vec::new();
        for i in 0..150 {
            name.push(match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            });
        }
        g.add_vertex(SeqVertex::new(1, String::from_utf8(name.clone()).unwrap()));
        let ctx = ctx_k4();
        let mut seq = name.clone();
        seq[60] = if seq[60] == b'A' { b'T' } else { b'A' };
        seq[80] = if seq[80] == b'A' { b'T' } else { b'A' };
        seq[100] = if seq[100] == b'A' { b'T' } else { b'A' };
        let mut t = ThreadingCtx::new(&g, &ctx);
        let p = find_path_in_graph(&mut t, &seq, &LocInGraph::new(1, 0));
        // DP 对齐应成功（3 mismatch 分散，div 3/150=0.02）
        let p = p.unwrap();
        assert_eq!(p.path, vec![1]);
        assert_eq!(p.mismatch_count, 3);
        assert_eq!(p.positions[0].1, 150); // i = verSeq.len
        assert_eq!(p.positions[0].2, 150); // j = seq.len
    }

    #[test]
    fn read_and_map_full_flow_and_reanchor() {
        // get_read_starts：锚点失效 → 再锚定滑窗 i 从 2 起
        let mut g = DiGraph::new();
        g.add_vertex(SeqVertex::new(5, "GGGGCCCC"));
        let ctx = ctx_k4();
        let mut orig_map = FxHashMap::default();
        orig_map.insert(5, LocInGraph::new(5, 0));
        let mut kmer_map = FxHashMap::default();
        kmer_map.insert("GGGG".to_string(), 5);
        // from_orig_v=999 不存在；i=2..4 的 kmer 不在图；i=5 命中 "GGGG"
        let raw = RawRead {
            name: "r".to_string(),
            seq: "TTTTTGGGGCCCC".to_string(),
            start_in_read: 0,
            end_in_read: 12, // f3=8 + K=4；end<len(13)；substring(0,13) 全长
            from_orig_v: 999,
        };
        let reads = get_read_starts(&g, &ctx, &orig_map, &kmer_map, &[raw]);
        let paths = reads.get("r").unwrap();
        assert_eq!(paths.len(), 1);
        assert_eq!(paths[0].path, vec![5]);
        assert_eq!(paths[0].mismatch_count, 0);
    }

    #[test]
    fn end_in_read_past_seq_returns_unmapped() {
        let g = chain_graph();
        let ctx = ctx_k4();
        let raw = RawRead {
            name: "r".to_string(),
            seq: "ACGT".to_string(),
            start_in_read: 0,
            end_in_read: 4, // >= len(4)
            from_orig_v: 1,
        };
        let reads = get_read_starts(&g, &ctx, &map1(), &FxHashMap::default(), &[raw]);
        assert!(reads.is_empty());
    }

    #[test]
    fn substring_off_by_one_is_preserved() {
        // endInRead = f3+K，substring(start, end+1)：f3=7,K=4 → end=11，取 [0,11]
        let g = chain_graph();
        let ctx = ctx_k4();
        let raw = RawRead {
            name: "r".to_string(),
            seq: "ACGTAAAACCCCTTTT".to_string(),
            start_in_read: 0,
            end_in_read: 11,
            from_orig_v: 1,
        };
        let reads = get_read_starts(&g, &ctx, &map1(), &FxHashMap::default(), &[raw]);
        let p = &reads["r"][0];
        // 12bp：v1 完整 + v2 从 K-1 起的 4 个 → path [1,2]
        assert_eq!(p.path, vec![1, 2]);
        assert_eq!(p.mismatch_count, 0);
    }

    #[test]
    fn orig_ver_ids_mapping_hash_segments() {
        // 合并节点 10（>last_real_id? 设 last_real_id=9 → 10 为新节点）prev 段
        let mut g = DiGraph::new();
        let mut v = SeqVertex::new(10, "ACGTAAAACCCC");
        v.prev_vertices_id = vec![vec![1, 2], vec![1, 2], vec![3]];
        g.add_vertex(v);
        let h = get_original_ver_ids_mapping_hash(&mut g, 9);
        // 相邻相等段去重 → [[1,2],[3]]；vid=10>last_real_id → loc 从 -1 起：
        // 1,2 段 → loc=0；3 段 → loc=1
        assert_eq!(h.get(&1).unwrap(), &LocInGraph::new(10, 0));
        assert_eq!(h.get(&2).unwrap(), &LocInGraph::new(10, 0));
        assert_eq!(h.get(&3).unwrap(), &LocInGraph::new(10, 1));
        assert!(!h.contains_key(&10)); // vid > last_real_id 不注册自身
                                       // 自身注册路径
        let mut g2 = DiGraph::new();
        let mut v2 = SeqVertex::new(4, "ACGT");
        v2.prev_vertices_id = vec![vec![7]];
        g2.add_vertex(v2);
        let h2 = get_original_ver_ids_mapping_hash(&mut g2, 9);
        assert_eq!(h2.get(&4).unwrap(), &LocInGraph::new(4, 0));
        assert_eq!(h2.get(&7).unwrap(), &LocInGraph::new(4, 1));
    }

    #[test]
    fn trimmed_path_behavior() {
        let mut p = ReadPath::from_node_end(3, 1, 50, 60);
        p.add_node_front(2, 0, 30, 40);
        p.add_node_front(1, 0, 10, 20);
        assert_eq!(p.path, vec![1, 2, 3]);
        assert_eq!(p.mismatch_count, 1);
        // min=25：前端 read_base_end(20)<25 删 → [2,3]（末端 vertex_base_end 50≥25 保留）
        assert_eq!(p.get_trimmed_path(25), vec![2, 3]);
        // min=55：前端 20、40 均 <55 连删，末端 50<55 再删 → []
        assert_eq!(p.get_trimmed_path(55), Vec::<i32>::new());
        assert_eq!(p.get_trimmed_path(0), vec![1, 2, 3]);
    }

    #[test]
    fn memo_hit_returns_deep_copy() {
        // 同一 read 二次穿线：命中 memo 返回拷贝（改一份不影响缓存）
        let g = chain_graph();
        let ctx = ctx_k4();
        let seq = b"ACGTAAAACCCCTTTT";
        let mut t = ThreadingCtx::new(&g, &ctx);
        let p1 = find_path_in_graph(&mut t, seq, &LocInGraph::new(1, 0)).unwrap();
        let key = "3_3_12".to_string(); // v3 起点（memo 存从该节点起的子路径）
        let cached = t.memo.get(&key).unwrap().as_ref().unwrap().clone();
        assert_eq!(cached, ReadPath::from_node_end(3, 0, 7, 16));
        // 深拷贝断言：修改返回值不影响缓存
        let mut p1m = p1.clone();
        p1m.path.push(99);
        assert_ne!(t.memo.get(&key).unwrap().as_ref().unwrap(), &p1m);
    }
}
