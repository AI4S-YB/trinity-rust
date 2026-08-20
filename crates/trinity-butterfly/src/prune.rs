//! 剪枝链：镜像 `TransAssembly_allProbPaths.java` 的
//! fixExtremelyHighSingleEdges (L8787) / removeLightEdges 系 (L13010-13222) /
//! compactLinearPaths (L12936) / removeSingleNtBubbles (L8603) /
//! calcSubComponentsStats (L13320)。
//!
//! main 调用序（L780-890）：
//! `fix → removeLight → [检查点 B] → compact → (My_DFS, T3) → [检查点 C]
//!  → (COLLAPSE_SNPs 时 bubbles → [检查点 H]) → calcSubComponents → [检查点 D]`。
//!
//! 比较符 quirk（逐字保留）：
//! - light in/out：`weight <= total*EDGE_THR` 删（**含等号**）
//! - light flow：`weight < weightAvg*FLOW_THR` 删（**严格小于**）
//! - fix：`supp > in*200 && supp > out*200`（严格大于两侧）
//! - SNP bubbles：`w(v,v1) > w(v,v2)` 保 v1 否则保 v2（平局保 v2）

use rustc_hash::FxHashMap;

use crate::context::BflyContext;
use crate::graph::{DiGraph, SeqVertex, SimpleEdge, T_VERTEX_ID, VERTEX_ROOT_ID};

/// EXTREME_EDGE_FLOW_FACTOR（L102）。
pub const EXTREME_EDGE_FLOW_FACTOR: f64 = 200.0;
/// COMP_AVG_COV_THR（L57）= 1；移除条件 `avgCov < COMP_AVG_COV_THR - 0.5`。
pub const COMP_AVG_COV_THR: f32 = 1.0;

// ---------------------------------------------------------------------------
// fixExtremelyHighSingleEdges (L8787)
// ---------------------------------------------------------------------------

/// 把支撑度远超两端流量的边压回 max(in, out)。
///
/// 条件：`inFlow(source)` 与 `outFlow(target)` 都存在，且
/// `supp > in*200 && supp > out*200`（严格 > 两侧）→ `weight = max(in, out)`。
pub fn fix_extremely_high_single_edges(
    graph: &mut DiGraph,
    out_flow: &FxHashMap<i32, f64>,
    in_flow: &FxHashMap<i32, f64>,
) {
    let edge_keys: Vec<(i32, i32)> = graph.edges_sorted().into_iter().map(|(k, _)| k).collect();
    for (u, v) in edge_keys {
        let (Some(&in_f), Some(&out_f)) = (in_flow.get(&u), out_flow.get(&v)) else {
            continue;
        };
        let supp = graph.find_edge(u, v).expect("edge exists").weight;
        if supp > in_f * EXTREME_EDGE_FLOW_FACTOR && supp > out_f * EXTREME_EDGE_FLOW_FACTOR {
            let new_supp = in_f.max(out_f);
            graph.find_edge_mut(u, v).expect("edge exists").weight = new_supp;
        }
    }
}

// ---------------------------------------------------------------------------
// removeLightEdges 系 (L13010-13222)
// ---------------------------------------------------------------------------

/// `removeLightEdges` = in + out + flow 三连（comp 部分已注释掉）。
pub fn remove_light_edges(graph: &mut DiGraph, ctx: &BflyContext) -> bool {
    let in_c = remove_light_in_edges(graph, ctx.edge_thr);
    let out_c = remove_light_out_edges(graph, ctx.edge_thr);
    let flow_c = remove_light_flow_edges(graph, ctx.flow_thr);
    in_c || out_c || flow_c
}

/// `atSimpleCycle`：任一入边源有 v→源的回边，或任一出边目标有 目标→v 的回边。
fn at_simple_cycle(graph: &DiGraph, v: i32) -> bool {
    graph
        .get_predecessors(v)
        .iter()
        .any(|&p| graph.find_edge(v, p).is_some())
        || graph
            .get_successors(v)
            .iter()
            .any(|&s| graph.find_edge(s, v).is_some())
}

/// `removeLightInEdges`：inDegree<=1 跳过；atSimpleCycle 跳过；
/// `totalIn` 按**删除前**的入边一次算清；`weight <= totalIn*EDGE_THR` 删（含等号）。
pub fn remove_light_in_edges(graph: &mut DiGraph, edge_thr: f64) -> bool {
    let queue: Vec<i32> = graph.vertex_ids().to_vec();
    let mut changed = false;
    for v in queue {
        if graph.in_degree(v) <= 1 {
            continue;
        }
        if at_simple_cycle(graph, v) {
            continue;
        }
        let preds = graph.get_predecessors(v).to_vec();
        let total_in: f64 = preds
            .iter()
            .map(|&p| graph.find_edge(p, v).expect("edge exists").weight)
            .sum();
        let thr = total_in * edge_thr;
        for p in preds {
            if graph.find_edge(p, v).expect("edge exists").weight <= thr {
                graph.remove_edge(p, v);
                changed = true;
            }
        }
    }
    changed
}

/// `removeLightOutEdges`：outDegree<=1 跳过；atSimpleCycle 跳过；`<=` 含等号。
pub fn remove_light_out_edges(graph: &mut DiGraph, edge_thr: f64) -> bool {
    let queue: Vec<i32> = graph.vertex_ids().to_vec();
    let mut changed = false;
    for v in queue {
        if graph.out_degree(v) <= 1 {
            continue;
        }
        if at_simple_cycle(graph, v) {
            continue;
        }
        let succs = graph.get_successors(v).to_vec();
        let total_out: f64 = succs
            .iter()
            .map(|&s| graph.find_edge(v, s).expect("edge exists").weight)
            .sum();
        let thr = total_out * edge_thr;
        for s in succs {
            if graph.find_edge(v, s).expect("edge exists").weight <= thr {
                graph.remove_edge(v, s);
                changed = true;
            }
        }
    }
    changed
}

/// `removeLightFlowEdges`：in==0 && out==0 跳过；
/// `weight < weightAvg*FLOW_THR` 删（**严格 <**，weightAvg 为 Java round 后的 int）。
pub fn remove_light_flow_edges(graph: &mut DiGraph, flow_thr: f64) -> bool {
    let all: Vec<i32> = graph.vertex_ids().to_vec();
    let mut changed = false;
    for v in all {
        if graph.in_degree(v) == 0 && graph.out_degree(v) == 0 {
            continue;
        }
        let thr = graph.get_vertex(v).expect("vertex exists").get_weight_avg() as f64 * flow_thr;
        // 快照收集后统一删（Java 用 HashSet removeEdges）
        let mut to_remove: Vec<(i32, i32)> = Vec::new();
        for &s in graph.get_successors(v) {
            if graph.find_edge(v, s).expect("edge exists").weight < thr {
                to_remove.push((v, s));
            }
        }
        for &p in graph.get_predecessors(v) {
            if graph.find_edge(p, v).expect("edge exists").weight < thr {
                to_remove.push((p, v));
            }
        }
        for (u, w) in to_remove {
            graph.remove_edge(u, w);
            changed = true;
        }
    }
    changed
}

// ---------------------------------------------------------------------------
// compactLinearPaths (L12936)
// ---------------------------------------------------------------------------

/// 线性链压缩。快照遍历顶点；`while (v1 != ROOT && outDegree(v1)==1)` 语义。
///
/// break 条件：`inDegree(v2)!=1 || v2.toBeDeleted || v2==T || v1==v2`。
/// 吸收时 `v1.concat_vertex(v2, e.w, lastRealID)`，v2 出边整体搬到 v1
/// （SimpleEdge 拷贝构造：weight/isInCircle/numLoops 保留），被吸收节点延迟到末尾统一删除。
pub fn compact_linear_paths(graph: &mut DiGraph, ctx: &BflyContext) -> bool {
    let snapshot: Vec<i32> = graph.vertex_ids().to_vec();
    let mut remove_vertices: Vec<i32> = Vec::new();
    let mut changed = false;

    for &v1 in &snapshot {
        while v1 != VERTEX_ROOT_ID && graph.out_degree(v1) == 1 {
            let v2 = graph.get_successors(v1)[0];
            {
                let v2v = graph.get_vertex(v2).expect("vertex exists");
                if graph.in_degree(v2) != 1 || v2v.to_be_deleted || v2 == T_VERTEX_ID || v1 == v2 {
                    break;
                }
            }
            let e_w = graph.find_edge(v1, v2).expect("edge exists").weight;
            let v2_clone = graph.get_vertex(v2).expect("vertex exists").clone();
            graph
                .get_vertex_mut(v1)
                .expect("vertex exists")
                .concat_vertex(v2_clone, e_w, ctx.last_real_id, ctx.kmer_size);
            graph
                .get_vertex_mut(v2)
                .expect("vertex exists")
                .to_be_deleted = true;
            remove_vertices.push(v2);
            changed = true;

            // v2 出边搬到 v1（先加新边再删旧边，Java 同序）
            let out_snapshot: Vec<i32> = graph.get_successors(v2).to_vec();
            for &v3 in &out_snapshot {
                let e2 = graph.find_edge(v2, v3).expect("edge exists").clone();
                let new_edge = SimpleEdge {
                    from_vertex_id: v1,
                    to_vertex_id: v3,
                    ..e2
                };
                graph.add_edge(v1, v3, new_edge);
                graph.remove_edge(v2, v3);
            }
            graph.remove_edge(v1, v2);
        }
    }
    for v in remove_vertices {
        graph.remove_vertex(v);
    }
    changed
}

// ---------------------------------------------------------------------------
// removeSingleNtBubbles (L8603)
// ---------------------------------------------------------------------------

/// 单核苷酸 bubble 坍缩（USE_DEGENERATE_CODE=false 的默认变体）。
///
/// 条件：succ(v)==2，两分支 `getNameKmerAdj().len()==K`，分支各自 succ==1 且汇合同一 vend。
/// `w(v,v1) > w(v,v2)` 保 v1 否则保 v2（**平局保 v2**）。新节点 `getNextID`，
/// `copyTheRest(vToKeep)` + `addToPrevIDs(keep, remove, lastRealID)`（quirk 见 SeqVertex），
/// 新边权 = 两平行边之和。removeV 延迟统一删除。
pub fn remove_single_nt_bubbles(graph: &mut DiGraph, ctx: &mut BflyContext) {
    let all_v: Vec<i32> = graph.vertex_ids().to_vec();
    let k = ctx.kmer_size;
    let mut remove_v: Vec<i32> = Vec::new();

    for &v in &all_v {
        if remove_v.contains(&v) {
            continue;
        }
        if graph.get_successors(v).len() != 2 {
            continue;
        }
        let (v1, v2) = (graph.get_successors(v)[0], graph.get_successors(v)[1]);
        let len1 = graph
            .get_vertex(v1)
            .expect("vertex exists")
            .get_name_kmer_adj(k, graph.in_degree(v1) > 0)
            .len();
        let len2 = graph
            .get_vertex(v2)
            .expect("vertex exists")
            .get_name_kmer_adj(k, graph.in_degree(v2) > 0)
            .len();
        if len1 != k || len2 != k {
            continue;
        }
        if graph.get_successors(v1).len() != 1 || graph.get_successors(v2).len() != 1 {
            continue;
        }
        let (s1, s2) = (graph.get_successors(v1)[0], graph.get_successors(v2)[0]);
        if s1 != s2 {
            continue;
        }
        let vend = s1;
        let w_v1 = graph.find_edge(v, v1).expect("edge exists").weight;
        let w_v2 = graph.find_edge(v, v2).expect("edge exists").weight;
        let (v_keep, v_remove) = if w_v1 > w_v2 { (v1, v2) } else { (v2, v1) };
        let e1_to_keep = graph.find_edge(v, v_keep).expect("edge exists").weight;
        let e2_to_keep = graph.find_edge(v_keep, vend).expect("edge exists").weight;
        let e1_to_remove = graph.find_edge(v, v_remove).expect("edge exists").weight;
        let e2_to_remove = graph.find_edge(v_remove, vend).expect("edge exists").weight;

        let new_id = ctx.get_next_id();
        let mut new_v = SeqVertex::new(
            new_id,
            graph
                .get_vertex(v_keep)
                .expect("vertex exists")
                .name
                .clone(),
        );
        new_v.copy_the_rest(graph.get_vertex(v_keep).expect("vertex exists"));
        new_v.add_to_prev_ids(
            graph.get_vertex(v_keep).expect("vertex exists"),
            graph.get_vertex(v_remove).expect("vertex exists"),
            ctx.last_real_id,
        );
        graph.add_vertex(new_v);
        graph.add_edge(
            v,
            new_id,
            SimpleEdge::new(e1_to_keep + e1_to_remove, v, new_id),
        );
        graph.add_edge(
            new_id,
            vend,
            SimpleEdge::new(e2_to_keep + e2_to_remove, new_id, vend),
        );

        remove_v.push(v_remove);
        remove_v.push(v_keep);
    }
    for rv in remove_v {
        graph.remove_vertex(rv);
    }
}

// ---------------------------------------------------------------------------
// calcSubComponentsStats (L13320)
// ---------------------------------------------------------------------------

/// 无向弱连通分量（JUNG WeakComponentClusterer；含 ROOT/T，按 vertex_ids 插入序输出）。
fn divide_into_components(graph: &DiGraph) -> Vec<Vec<i32>> {
    fn find(mut x: i32, parent: &mut FxHashMap<i32, i32>) -> i32 {
        loop {
            let p = *parent.entry(x).or_insert(x);
            if p == x {
                break x;
            }
            x = p;
        }
    }
    let mut parent: FxHashMap<i32, i32> = FxHashMap::default();
    for &v in graph.vertex_ids() {
        parent.insert(v, v);
    }
    for &u in graph.vertex_ids() {
        for &w in graph.get_successors(u) {
            let (ru, rw) = (find(u, &mut parent), find(w, &mut parent));
            if ru != rw {
                parent.insert(ru, rw);
            }
        }
    }
    let mut comps: FxHashMap<i32, Vec<i32>> = FxHashMap::default();
    for &v in graph.vertex_ids() {
        comps.entry(find(v, &mut parent)).or_default().push(v);
    }
    let mut order: Vec<i32> = graph.vertex_ids().to_vec();
    // 保持与插入序对应的稳定输出（Java HashSet 顺序任意，但各分量处理互不影响）
    order.sort_by_key(|&v| find(v, &mut parent));
    let mut out: Vec<Vec<i32>> = Vec::new();
    let mut last_root = i32::MIN;
    for v in order {
        let r = find(v, &mut parent);
        if r != last_root {
            out.push(comps.remove(&r).expect("root present"));
            last_root = r;
        }
    }
    out
}

/// `calcSubComponentsStats`：小组件移除。
///
/// 每分量 allW = 所有节点 weights + 所有出边权重；`allW 空` 或
/// `单节点且 name.len < min_output_seq`（-L，Java MIN_OUTPUT_SEQ）→ 删该节点；
/// 否则 `avgCov = (int 累加 t)/size`（Java `int t += Double` 每步截断），
/// `avgCov < COMP_AVG_COV_THR - 0.5`（= 0.5）→ 删整分量。
pub fn calc_sub_components_stats(graph: &mut DiGraph, min_output_seq: usize) {
    let comps = divide_into_components(graph);
    for comp in comps {
        let mut all_w: Vec<f64> = Vec::new();
        for &v in &comp {
            all_w.extend_from_slice(&graph.get_vertex(v).expect("vertex exists").weights);
            for &s in graph.get_successors(v) {
                all_w.push(graph.find_edge(v, s).expect("edge exists").weight);
            }
        }
        let comp_id = comp[0];
        if all_w.is_empty()
            || (comp.len() == 1
                && graph.get_vertex(comp_id).expect("vertex exists").name.len() < min_output_seq)
        {
            graph.remove_vertex(comp_id);
            continue;
        }
        // Java: int t; for (Double w : allW) t += w;（复合赋值每步向零截断）
        let mut t: i32 = 0;
        for &w in &all_w {
            t = (t as f64 + w) as i32;
        }
        let avg_cov = t as f32 / all_w.len() as f32;
        if avg_cov < COMP_AVG_COV_THR - 0.5 {
            for &v in &comp {
                graph.remove_vertex(v);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// main 剪枝链（L780-890 的本任务覆盖段；My_DFS 属 T3）
// ---------------------------------------------------------------------------

/// 依次执行 fix → removeLight → compact → bubbles → calcSubComponents。
/// 检查点 B/C/H/D 由调用方在各步之间自行 `write_dot_string`。
pub fn run_pruning_chain(
    graph: &mut DiGraph,
    ctx: &mut BflyContext,
    in_flow: &FxHashMap<i32, f64>,
    out_flow: &FxHashMap<i32, f64>,
    min_output_seq: usize,
) {
    fix_extremely_high_single_edges(graph, out_flow, in_flow);
    remove_light_edges(graph, ctx);
    compact_linear_paths(graph, ctx);
    remove_single_nt_bubbles(graph, ctx);
    calc_sub_components_stats(graph, min_output_seq);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn edge(graph: &DiGraph, u: i32, v: i32) -> f64 {
        graph.find_edge(u, v).expect("edge").weight
    }

    // ------------------------------------------------------------------
    // fixExtremelyHighSingleEdges
    // ------------------------------------------------------------------

    #[test]
    fn fix_extreme_requires_both_sides_strictly() {
        let mut g = DiGraph::new();
        g.add_vertex(SeqVertex::new(1, "AAAA"));
        g.add_vertex(SeqVertex::new(2, "CCCC"));
        g.add_vertex(SeqVertex::new(3, "GGGG"));
        // 200*5=1000：supp=1000 不严格大于任一侧 → 不动
        g.add_edge(1, 2, SimpleEdge::new(1000.0, 1, 2));
        // 两侧都 1 → 200 阈：supp=201 严格大于两侧 → 压回 max(1,1)=1
        g.add_edge(2, 3, SimpleEdge::new(201.0, 2, 3));
        // 只一侧超：in(3)=1001? out(3) 缺失 → 跳过
        let mut in_flow = FxHashMap::default();
        let mut out_flow = FxHashMap::default();
        in_flow.insert(1, 5.0);
        out_flow.insert(2, 5.0);
        in_flow.insert(2, 1.0);
        out_flow.insert(3, 1.0);
        fix_extremely_high_single_edges(&mut g, &out_flow, &in_flow);
        assert_eq!(edge(&g, 1, 2), 1000.0); // 等号不触发
        assert_eq!(edge(&g, 2, 3), 1.0); // max(in,out)
    }

    // ------------------------------------------------------------------
    // removeLight in/out：<= 含等号；flow：< 严格
    // ------------------------------------------------------------------

    #[test]
    fn light_in_edges_inclusive_threshold() {
        let mut g = DiGraph::new();
        for id in [1, 2, 3, 4] {
            g.add_vertex(SeqVertex::new(id, "ACGT"));
        }
        // totalIn = 21；thr = 21*0.05 = 1.05 → 1 <= 1.05 删，10 保留
        g.add_edge(1, 4, SimpleEdge::new(10.0, 1, 4));
        g.add_edge(2, 4, SimpleEdge::new(10.0, 2, 4));
        g.add_edge(3, 4, SimpleEdge::new(1.0, 3, 4)); // total=21, thr=1.05
        assert!(remove_light_in_edges(&mut g, 0.05));
        assert!(g.find_edge(3, 4).is_none());
        assert!(g.find_edge(1, 4).is_some());
        assert!(g.find_edge(2, 4).is_some());
        // 恰等：total = 19+1 = 20, thr(0.05) = 1，边权恰 1 → 删（<= 含等号）
        let mut g2 = DiGraph::new();
        for id in [1, 2, 3] {
            g2.add_vertex(SeqVertex::new(id, "ACGT"));
        }
        g2.add_edge(1, 3, SimpleEdge::new(19.0, 1, 3));
        g2.add_edge(2, 3, SimpleEdge::new(1.0, 2, 3));
        assert!(remove_light_in_edges(&mut g2, 0.05));
        assert!(g2.find_edge(2, 3).is_none());
    }

    #[test]
    fn light_flow_edges_strict_threshold() {
        let mut g = DiGraph::new();
        for id in [1, 2, 3] {
            g.add_vertex(SeqVertex::new(id, "ACGT"));
        }
        // v=2：weightAvg=round(50)=50；thr = 50*0.02 = 1；边权恰 1 → 严格 < 不删
        let v2 = g.get_vertex_mut(2).unwrap();
        v2.weights = vec![50.0];
        g.add_edge(1, 2, SimpleEdge::new(1.0, 1, 2));
        g.add_edge(2, 3, SimpleEdge::new(0.99, 2, 3));
        assert!(remove_light_flow_edges(&mut g, 0.02));
        assert!(g.find_edge(1, 2).is_some(), "恰等不删（严格 <）");
        assert!(g.find_edge(2, 3).is_none(), "0.99 < 1 删");
    }

    #[test]
    fn at_simple_cycle_skips_pruning() {
        // 1→2, 2→1（简单环）；2 的出度 2：2→1, 2→3
        let mut g = DiGraph::new();
        for id in [1, 2, 3] {
            g.add_vertex(SeqVertex::new(id, "ACGT"));
        }
        g.add_edge(1, 2, SimpleEdge::new(10.0, 1, 2));
        g.add_edge(2, 1, SimpleEdge::new(1.0, 2, 1)); // 回边
        g.add_edge(2, 3, SimpleEdge::new(1.0, 2, 3));
        assert!(!remove_light_out_edges(&mut g, 0.02), "atSimpleCycle 跳过");
        assert_eq!(g.edge_count(), 3);
        // 去掉回边后：total=2, thr=0.04 → 1.0 > 0.04 都不删
        let mut g2 = g.clone();
        g2.remove_edge(2, 1);
        assert!(!remove_light_out_edges(&mut g2, 0.02));
        // 出度 2 且无回边：total=2, thr(1.0)=2 → 两条 1.0 边都 <= 2 删
        let mut g3 = DiGraph::new();
        for id in [1, 2, 3, 4] {
            g3.add_vertex(SeqVertex::new(id, "ACGT"));
        }
        g3.add_edge(2, 3, SimpleEdge::new(1.0, 2, 3));
        g3.add_edge(2, 4, SimpleEdge::new(1.0, 2, 4));
        assert!(remove_light_out_edges(&mut g3, 1.0));
        assert_eq!(g3.edge_count(), 0);
    }

    // ------------------------------------------------------------------
    // compactLinearPaths
    // ------------------------------------------------------------------

    #[test]
    fn compact_merges_linear_chain() {
        // ROOT→1→2→3→T，K=4；名字首尾相接
        let mut g = DiGraph::new();
        let mut ctx = BflyContext::new();
        ctx.kmer_size = 4;
        ctx.last_real_id = 3;
        g.add_vertex(SeqVertex::new(VERTEX_ROOT_ID, "S"));
        g.add_vertex(SeqVertex::with_per_base_weight(1, "ACGT", 2.0));
        g.add_vertex(SeqVertex::new(2, "CGTA"));
        g.add_vertex(SeqVertex::new(3, "GTAC"));
        g.add_vertex(SeqVertex::new(T_VERTEX_ID, "E"));
        g.add_edge(VERTEX_ROOT_ID, 1, SimpleEdge::new(2.0, VERTEX_ROOT_ID, 1));
        g.add_edge(1, 2, SimpleEdge::new(5.0, 1, 2));
        g.add_edge(2, 3, SimpleEdge::new(7.0, 2, 3));
        g.add_edge(3, T_VERTEX_ID, SimpleEdge::new(9.0, 3, T_VERTEX_ID));
        assert!(compact_linear_paths(&mut g, &ctx));
        // 1 吸收 2、3（3 的后继是 T，但吸收发生在检查 v2==T 之前——检查的是 v2 的
        // 后继无要求；T 不可作为被吸收者）→ 1→T
        assert_eq!(g.vertex_count(), 3); // ROOT, 1, T
        let v1 = g.get_vertex(1).unwrap();
        assert_eq!(v1.name, "ACGTAC"); // +A +C
        assert_eq!(v1.weights, vec![2.0, 2.0, 2.0, 2.0, 5.0, 7.0]);
        assert_eq!(v1.prev_vertices_id, vec![vec![2], vec![3]]);
        assert_eq!(edge(&g, 1, T_VERTEX_ID), 9.0);
    }

    #[test]
    fn compact_breaks_on_branch_and_t() {
        // 1 出度 2 → 不压缩；2→T：v2==T break
        let mut g = DiGraph::new();
        let mut ctx = BflyContext::new();
        ctx.kmer_size = 4;
        g.add_vertex(SeqVertex::new(1, "ACGT"));
        g.add_vertex(SeqVertex::new(2, "CGTA"));
        g.add_vertex(SeqVertex::new(3, "GTAC"));
        g.add_vertex(SeqVertex::new(T_VERTEX_ID, "E"));
        g.add_edge(1, 2, SimpleEdge::new(1.0, 1, 2));
        g.add_edge(1, 3, SimpleEdge::new(1.0, 1, 3));
        g.add_edge(2, T_VERTEX_ID, SimpleEdge::new(1.0, 2, T_VERTEX_ID));
        assert!(!compact_linear_paths(&mut g, &ctx));
        assert_eq!(g.vertex_count(), 4);
        assert_eq!(g.edge_count(), 3);
    }

    // ------------------------------------------------------------------
    // removeSingleNtBubbles
    // ------------------------------------------------------------------

    fn bubble_graph(k: usize) -> (DiGraph, BflyContext) {
        // v=1 → (2,3) → vend=4；K=4 时 kmerAdj len==4 → 分支名长 2K-1=7
        let mut g = DiGraph::new();
        let mut ctx = BflyContext::new();
        ctx.kmer_size = k;
        ctx.last_id = 4;
        ctx.last_real_id = 4;
        g.add_vertex(SeqVertex::new(1, "ACGTA"));
        g.add_vertex(SeqVertex::new(2, "CGTGAGC")); // kmerAdj = GTGA? len(name)-(k-1)=3 ≠ 4
        g.add_vertex(SeqVertex::new(3, "CGTGACT"));
        g.add_vertex(SeqVertex::new(4, "CTACGT"));
        g.add_edge(1, 2, SimpleEdge::new(10.0, 1, 2));
        g.add_edge(1, 3, SimpleEdge::new(10.0, 1, 3));
        g.add_edge(2, 4, SimpleEdge::new(10.0, 2, 4));
        g.add_edge(3, 4, SimpleEdge::new(10.0, 3, 4));
        (g, ctx)
    }

    #[test]
    fn bubbles_tie_keeps_v2_and_sums_weights() {
        let (mut g, mut ctx) = bubble_graph(4);
        // 平局（10 == 10）→ 保 v2（第二个后继）；新 id = 5
        remove_single_nt_bubbles(&mut g, &mut ctx);
        assert_eq!(g.vertex_count(), 3); // 1, 4, new(5)
        assert!(!g.contains_vertex(2));
        assert!(!g.contains_vertex(3));
        let nv = g.get_vertex(5).unwrap();
        // copyTheRest(v2)：name 用 vToKeep=v2 的名
        assert_eq!(nv.name, "CGTGACT");
        // 新边权求和
        assert_eq!(edge(&g, 1, 5), 20.0);
        assert_eq!(edge(&g, 5, 4), 20.0);
        // addToPrevIDs quirk：copyTheRest 后 prev 为空（v2 无历史）→ 执行；
        // id 2/3 >= lastRealID(4)? 2<4、3<4 → 都不加；v2/v3 prev 空 → 空 Vec
        assert_eq!(nv.prev_vertices_id, vec![Vec::<i32>::new()]);
    }

    #[test]
    fn bubbles_strict_greater_keeps_v1_and_prev_quirk() {
        let (mut g, mut ctx) = bubble_graph(4);
        g.find_edge_mut(1, 2).unwrap().weight = 11.0; // 严格 > → 保 v1
                                                      // 给 v1(prev keep) 一段 prev 历史 → copyTheRest 拷过来 → addToPrevIDs no-op
        g.get_vertex_mut(2).unwrap().prev_vertices_id = vec![vec![9]];
        remove_single_nt_bubbles(&mut g, &mut ctx);
        let nv = g.get_vertex(5).unwrap();
        assert_eq!(nv.name, "CGTGAGC");
        assert_eq!(
            nv.prev_vertices_id,
            vec![vec![9]],
            "quirk：prev 非空 → no-op"
        );
    }

    #[test]
    fn bubbles_wrong_kmer_length_ignored() {
        let (mut g, mut ctx) = bubble_graph(7); // kmerAdj len 3 != 7
        let before = g.edge_count();
        remove_single_nt_bubbles(&mut g, &mut ctx);
        assert_eq!(g.edge_count(), before);
        assert_eq!(g.vertex_count(), 4);
        assert_eq!(ctx.last_id, 4, "无新 id 分配");
    }

    // ------------------------------------------------------------------
    // calcSubComponentsStats
    // ------------------------------------------------------------------

    #[test]
    fn sub_components_avg_cov_removal() {
        let mut g = DiGraph::new();
        // 分量 A：1→2，权重均值 10 → 保留
        g.add_vertex(SeqVertex::new(1, "ACGT"));
        g.add_vertex(SeqVertex::new(2, "CGTA"));
        g.get_vertex_mut(1).unwrap().weights = vec![10.0; 4];
        g.add_edge(1, 2, SimpleEdge::new(10.0, 1, 2));
        // 分量 B：3→4，权重全 0.2 → avg=0.2 < 0.5 → 删
        g.add_vertex(SeqVertex::new(3, "GGGG"));
        g.add_vertex(SeqVertex::new(4, "TTTT"));
        g.get_vertex_mut(3).unwrap().weights = vec![0.2; 4];
        g.add_edge(3, 4, SimpleEdge::new(0.2, 3, 4));
        calc_sub_components_stats(&mut g, 200);
        assert!(g.contains_vertex(1) && g.contains_vertex(2));
        assert!(!g.contains_vertex(3) && !g.contains_vertex(4));
    }

    #[test]
    fn sub_components_single_short_node_removed() {
        let mut g = DiGraph::new();
        // 单节点分量：无权重 → allW 只有出边权重 5；min_output_seq=200 → 单节点 name < 200 → 删
        g.add_vertex(SeqVertex::new(1, "ACGT"));
        g.add_vertex(SeqVertex::new(2, "TTTT"));
        g.add_edge(1, 2, SimpleEdge::new(5.0, 1, 2));
        // 让 2 成为独立单节点分量：删 1→2 边
        g.remove_edge(1, 2);
        calc_sub_components_stats(&mut g, 200);
        assert!(!g.contains_vertex(1), "单节点短序列被删");
        assert!(
            !g.contains_vertex(2),
            "无权重（allW 只有出边？2 无出边→allW 空）删"
        );
    }

    #[test]
    fn weak_components_via_undirected_connectivity() {
        // 1→2, 3→2：弱连通为一个分量（3 与 1 无向相通）
        let mut g = DiGraph::new();
        for id in [1, 2, 3] {
            g.add_vertex(SeqVertex::new(id, "ACGT"));
        }
        g.add_edge(1, 2, SimpleEdge::new(1.0, 1, 2));
        g.add_edge(3, 2, SimpleEdge::new(1.0, 3, 2));
        let comps = divide_into_components(&g);
        assert_eq!(comps.len(), 1);
        assert_eq!(comps[0], vec![1, 2, 3]);
        // 独立分量
        g.add_vertex(SeqVertex::new(9, "GGGG"));
        let comps = divide_into_components(&g);
        assert_eq!(comps.len(), 2);
    }

    #[test]
    fn java_int_truncating_sum() {
        // Java: int t; for (Double w : allW) t += w;
        // t 是整数，trunc(t + w) == t + trunc(w)，逐步截断等价于逐元素截断求和；
        // 负数向零截断（Java (int) cast 语义，Rust as i32 同）
        let all_w = vec![0.6_f64, 0.6, -0.5];
        let mut t: i32 = 0;
        for &w in &all_w {
            t = (t as f64 + w) as i32;
        }
        assert_eq!(t, 0);
        let all_w = vec![1.9_f64, 1.9];
        let mut t: i32 = 0;
        for &w in &all_w {
            t = (t as f64 + w) as i32;
        }
        assert_eq!(t, 2);
    }
}
