//! 单组件全链驱动（库层）——原版 `java -jar Butterfly.jar -C <prefix>` 的
//! main 编排（TransAssembly_allProbPaths.java L726-912）上移。
//!
//! 输入是 **文本**（graph.out / graph.reads 内容），输出是 allProbPaths FASTA
//! 文本——CLI、集成测试（c0/c1/c2）与 xtask 对拍共用同一入口，防漂移。
//!
//! 递归保护：POG 的 DFS_add_path 与 threading 递归在深组件上可能爆默认
//! 栈（Java 默认虚拟机栈可增长到数百 MB），因此在 `stack_size` 配置的
//! 专用线程上运行（默认 256MB，Trinity 主线 JVM 经验值）。

use std::io::Write;
use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::postprocess::PostProcessParams;
use crate::{DiGraph, PathSearchParams, SuffStats};
use trinity_common::error::CommonError;
use trinity_common::kmer::{get_ds_kmer_val, kmer_to_intval};

use crate::{
    build_new_graph_use_kmers, calc_sub_components_stats, compact_linear_paths,
    create_dag_from_overlap_layout, edge_read_coverage, fix_extremely_high_single_edges,
    get_original_ver_ids_mapping_hash, get_read_starts_ordered, get_suff_stats_w_pairs,
    parse_graph_reads, path_read_support, pre_process_graph_file, print_final_paths,
    remove_light_edges, remove_single_nt_bubbles, run_butterfly_all_prob_paths_with, run_dfs2,
};

/// 组件驱动参数（镜像 Butterfly.jar getopt 面的库层子集）。
#[derive(Debug, Clone)]
pub struct ComponentParams {
    /// `-N`：total reads / fragment pairs（Java printUsage 只要求非零）。
    pub n: u64,
    /// `-L`：最短报告序列长度。
    pub min_len: usize,
    /// `-F`：配对端最大跨度（max_pair_distance）。
    pub max_pair_distance: usize,
    /// `--NO_EM_REDUCE`（Trinity 主线传 true；getopt 默认 false）。
    pub no_em_reduce: bool,
    /// `--no_path_merging`。
    pub no_path_merging: bool,
    /// `-R`：最小 read 支持（默认 2）。
    pub min_read_support: i64,
    /// 显式 reinforcement 距离（覆盖 `-O` 百分比推算；None → 25% * max_pair_distance）。
    pub path_reinforcement_distance: Option<usize>,
    /// FASTA header 里的组件名（原版 = -C 值最后一段；printFinalPaths 去后缀）。
    pub name: String,
    /// 专用线程栈大小（递归逃生口；默认 256MB）。
    pub stack_size: usize,
    /// 输出层路径引擎（M2.1；默认 Classic 不变）。
    pub path_engine: PathEngine,
    /// PogFlow 流分解的去冗余阈值（M2.1 终审微调: 0.95→0.98 找回 ~120 条近缘 isoform）。
    pub flow_dedup_threshold: f64,
    /// `--flow-cov-fusion`（M8）: 无 read 证据边的 25-mer 覆盖度兜底权重
    /// （计数表, CAP）；None = 现状行为（floor 1）。仅 PogFlow 引擎消费。
    pub flow_cov: Option<(Arc<trinity_kmer::counter::CountMap>, u64)>,
}

impl Default for ComponentParams {
    fn default() -> Self {
        ComponentParams {
            n: 100_000,
            min_len: 200,
            max_pair_distance: 10_000,
            no_em_reduce: false,
            no_path_merging: false,
            min_read_support: 2,
            path_reinforcement_distance: None,
            name: "component".to_string(),
            stack_size: 256 * 1024 * 1024,
            flow_dedup_threshold: 0.98,
            path_engine: PathEngine::Classic,
            flow_cov: None,
        }
    }
}

/// 输出层路径引擎（M2.1）。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PathEngine {
    /// 原 `run_butterfly_all_prob_paths` 路径（不变）。
    #[default]
    Classic,
    /// POG 图（zipping 后 SeqVertex DAG）+ 最重路径剥离流分解（M2.1）。
    PogFlow,
}

/// POG 边 (u→v) 的 junction 25-mer（M8）: u 名末 ≤24 碱基 + v 经最长
/// 后缀-前缀重叠后的延伸碱基补足 25（重叠语义同 `path_to_seq`, 逐对节点
/// 计算——ov > len(u) 的角落情形与逐路径拼接的接缝可能差 1 碱基, 仅影响
/// 覆盖率查表键的一个窗口, 有界; 见 M8 spec）。k = --KMER_SIZE（计数表键长,
/// 默认 25; graph.out 顶点名长 k−1）。
/// v 无延伸（包含关系）或拼接不足 k → None。
pub fn junction_kmer(u_name: &[u8], v_name: &[u8], k: usize) -> Option<Vec<u8>> {
    let max_ov = u_name.len().min(v_name.len());
    let mut ov = max_ov;
    while ov > 0 && u_name[u_name.len() - ov..] != v_name[..ov] {
        ov -= 1;
    }
    let ext = &v_name[ov..];
    if ext.is_empty() {
        return None; // v 被 u 后缀包含, 无 junction
    }
    let take_u = u_name.len().min(k - 1);
    let need_v = k - take_u;
    if ext.len() < need_v {
        return None;
    }
    let mut kj = Vec::with_capacity(k);
    kj.extend_from_slice(&u_name[u_name.len() - take_u..]);
    kj.extend_from_slice(&ext[..need_v]);
    Some(kj)
}

/// DS canonical k-mer 计数查询（表缺/编码失败 → 0）。
fn lookup_cov(table: &trinity_kmer::counter::CountMap, kj: &[u8], k: usize) -> u64 {
    let Ok(val) = kmer_to_intval(kj) else {
        return 0;
    };
    let key = get_ds_kmer_val(val, k);
    table.get(&key).map(|&c| c as u64).unwrap_or(0)
}

/// 边权融合（M8 spec）: read 证据（POG 边权 pog_w / PairPath 覆盖 read_w）
/// 优先原样; 无证据边（两者皆 0）用 min(cov25, cap) 兜底, 下限 1。
/// flow_cov=None 时与 M3 行为逐边一致（floor 1）。
fn fuse_edge_weight(pog_w: u64, read_w: u64, cov25: u64, cap: u64) -> u64 {
    debug_assert!(cap > 0);
    let w = pog_w.max(read_w);
    if w > 0 {
        w
    } else {
        cov25.min(cap).max(1)
    }
}

/// butterfly DiGraph → flow FlowGraph。
/// 节点 kmer = SeqVertex.name（zipping 后全长序列，不等长）；边权 = SimpleEdge.weight
/// （f64 read 支持，max(0) round 防 NaN/负）与 PairPath read 边覆盖取 max 融合
/// （M3：read 覆盖高的边不被 POG 低权截断），下限 1。
/// M8: flow_cov 给出时, 无 read 证据的边改以 junction 25-mer 覆盖度（封顶）兜底。
fn digraph_to_flow(
    g: &DiGraph,
    read_cov: &FxHashMap<(i32, i32), i64>,
    flow_cov: Option<(&trinity_kmer::counter::CountMap, u64)>,
    k: usize,
) -> trinity_flow::graph::FlowGraph {
    let mut fg = trinity_flow::graph::FlowGraph::default();
    for &vid in g.vertex_ids() {
        if let Some(v) = g.get_vertex(vid) {
            fg.nodes.insert(
                vid as i64,
                trinity_flow::graph::FlowNode {
                    kmer: v.name.as_bytes().to_vec(),
                },
            );
        }
    }
    for ((u, v), e) in g.edges_sorted() {
        let pog_w = e.weight.max(0.0).round() as u64;
        let read_w = read_cov.get(&(u, v)).copied().unwrap_or(0) as u64;
        // 仅无证据边查 k-mer 表（省编码开销）; None 时 cov=0 → fuse 退化为 floor 1
        let cov25 = if pog_w.max(read_w) == 0 {
            match flow_cov {
                Some((table, _)) => junction_kmer(
                    g.get_vertex(u).unwrap().name.as_bytes(),
                    g.get_vertex(v).unwrap().name.as_bytes(),
                    k,
                )
                .map(|kj| lookup_cov(table, &kj, k))
                .unwrap_or(0),
                None => 0,
            }
        } else {
            0
        };
        let cap = flow_cov.map(|(_, c)| c).unwrap_or(1);
        let w = fuse_edge_weight(pog_w, read_w, cov25, cap);
        fg.edges.insert((u as i64, v as i64), w);
    }
    fg
}

/// 单组件运行结果。
pub struct ComponentResult {
    /// `<prefix>.allProbPaths.fasta` 全文（调用方决定写哪/打 stderr）。
    pub all_prob_paths_fasta: String,
    /// 最终路径（转录本）数。
    pub num_paths: usize,
}

/// 单组件全链：graph.out 文本 + graph.reads 文本 → allProbPaths FASTA 文本。
///
/// 链路 = graph 解析 → 剪枝链（B/C/H/D + dfs2 + bubble + stats）→ 穿线
/// → SuffStats → POG → 路径搜索 → T9 后处理 → printFinalPaths。
/// 在 `stack_size` 配置的专用线程上运行（递归逃生口）；失败返回 Err
/// （组件级隔离由调用方线程池做）。
pub fn run_component(
    graph_out: &str,
    graph_reads: &str,
    params: &ComponentParams,
) -> Result<ComponentResult, CommonError> {
    if params.n == 0 || params.max_pair_distance == 0 || params.min_read_support < 1 {
        return Err(CommonError::Parse(
            "缺少必选参数（-N/-F/-R 无效值）".to_string(),
        ));
    }
    // Java L700+：显式 --path_reinforcement_distance 优先，否则 O% * F（O 默认 25）
    let prd = params
        .path_reinforcement_distance
        .unwrap_or(params.max_pair_distance / 4);

    let graph_out = graph_out.to_string();
    let graph_reads = graph_reads.to_string();
    let params = params.clone();

    // 递归逃生口：POG DFS_add_path / threading 递归在大组件上会爆默认 8MB 栈。
    let handle = std::thread::Builder::new()
        .stack_size(params.stack_size)
        .spawn(move || run_component_inner(&graph_out, &graph_reads, &params, prd))
        .map_err(|e| CommonError::Parse(format!("cannot spawn component thread: {e}")))?;
    match handle.join() {
        Ok(r) => r,
        Err(_) => Err(CommonError::Parse(
            "butterfly 组件线程 panic（深递归/数据异常）".to_string(),
        )),
    }
}

fn run_component_inner(
    graph_text: &str,
    reads_text: &str,
    params: &ComponentParams,
    path_reinforcement_distance: usize,
) -> Result<ComponentResult, CommonError> {
    // 原始 kmer 图（边权重重标定基准）
    let br0 = build_new_graph_use_kmers(graph_text)
        .map_err(|e| CommonError::Parse(format!("graph 解析失败: {e}")))?;
    let mut orig_kmer_to_node: FxHashMap<String, i32> = FxHashMap::default();
    for &vid in br0.graph.vertex_ids() {
        orig_kmer_to_node.insert(br0.graph.get_vertex(vid).unwrap().name.clone(), vid);
    }

    // 主链 B/C/H/D + 穿线 + SuffStats + POG
    let br = build_new_graph_use_kmers(graph_text)
        .map_err(|e| CommonError::Parse(format!("graph 解析失败: {e}")))?;
    let (mut graph, mut ctx) = (br.graph, br.ctx);
    let (in_flow, out_flow, _) = pre_process_graph_file(graph_text);
    fix_extremely_high_single_edges(&mut graph, &out_flow, &in_flow);
    remove_light_edges(&mut graph, &ctx);
    compact_linear_paths(&mut graph, &ctx);
    run_dfs2(&mut graph);
    remove_single_nt_bubbles(&mut graph, &mut ctx);
    calc_sub_components_stats(&mut graph, params.min_len);

    let orig_id_map = get_original_ver_ids_mapping_hash(&mut graph, ctx.last_real_id);
    let raw_reads = parse_graph_reads(reads_text, ctx.kmer_size, false)
        .map_err(|e| CommonError::Parse(format!("reads 解析失败: {e}")))?;
    let ordered =
        get_read_starts_ordered(&graph, &ctx, &orig_id_map, &orig_kmer_to_node, &raw_reads);
    let suff = get_suff_stats_w_pairs(&graph, &ordered);
    let mut result = create_dag_from_overlap_layout(&graph, &mut ctx, &suff);

    // M2.1 PogFlow：POG 图上流分解，不进 run_butterfly_all_prob_paths。
    // M3：PairPath read 边覆盖注入流权重 + 产出路径 read 支持过滤。
    if params.path_engine == PathEngine::PogFlow {
        // PairPath 已重映射到 seqvertex 新 id（result.combined_read_hash），
        // 与流图节点 id 同域；suff 本身的路径是旧图 id，不能直接用。
        let suff_pog = SuffStats {
            combined_read_hash: std::mem::take(&mut result.combined_read_hash),
            ..Default::default()
        };
        let read_cov = edge_read_coverage(&suff_pog);
        let flow_g = digraph_to_flow(
            &result.seqvertex_graph,
            &read_cov,
            params.flow_cov.as_ref().map(|(t, c)| (t.as_ref(), *c)),
            ctx.kmer_size,
        );
        let kmers: Vec<(i64, Vec<u8>)> = flow_g
            .nodes
            .iter()
            .map(|(&id, n)| (id, n.kmer.clone()))
            .collect();
        let mut all: Vec<trinity_flow::paths::PathOutput> = Vec::new();
        for comp in flow_g.weak_components() {
            if comp.is_empty() {
                continue;
            }
            // 单节点分量（zipping 完全线性化）：节点自身即转录本候选。
            if comp.len() == 1 {
                let seq = trinity_flow::paths::path_to_seq(&comp, &kmers, ctx.kmer_size);
                if seq.len() >= params.min_len {
                    // M3：单节点分量同样做 read 支持过滤（替代 u64::MAX 哨兵）
                    let comp_i32: Vec<i32> = comp.iter().map(|&x| x as i32).collect();
                    let support = path_read_support(&suff_pog, &comp_i32);
                    if support >= params.min_read_support {
                        all.push(trinity_flow::paths::PathOutput {
                            seq,
                            support: support as u64,
                            path: comp,
                        });
                    }
                }
                continue;
            }
            let paths = trinity_flow::peel::peel_component(
                &flow_g,
                &comp,
                &trinity_flow::peel::PeelParams {
                    min_support: 1,
                    min_path_nodes: 2,
                },
            );
            for path in paths {
                let seq = trinity_flow::paths::path_to_seq(&path, &kmers, ctx.kmer_size);
                if seq.len() < params.min_len {
                    continue;
                }
                // M3 核心变更：PairPath 兼容计数替代边权最小值
                let path_i32: Vec<i32> = path.iter().map(|&x| x as i32).collect();
                let support = path_read_support(&suff_pog, &path_i32);
                if support < params.min_read_support {
                    continue;
                }
                all.push(trinity_flow::paths::PathOutput {
                    seq,
                    support: support as u64,
                    path,
                });
            }
        }
        let deduped = trinity_flow::paths::dedup_similar(all, params.flow_dedup_threshold);
        let fasta = trinity_flow::paths::to_fasta(&deduped, &params.name, 0, 60);
        return Ok(ComponentResult {
            all_prob_paths_fasta: fasta,
            num_paths: deduped.len(),
        });
    }

    let sp = PathSearchParams {
        kmer_size: ctx.kmer_size,
        min_output_seq: params.min_len,
        min_read_support_thr: params.min_read_support,
        path_reinforcement_distance,
        ..PathSearchParams::default()
    };
    let post = PostProcessParams {
        no_em_reduce: params.no_em_reduce,
        no_path_merging: params.no_path_merging,
        ..PostProcessParams::default()
    };
    let bfly = run_butterfly_all_prob_paths_with(
        &result.seqvertex_graph,
        &br0.graph,
        &graph,
        &result.combined_read_hash,
        &sp,
        &post,
    );

    let fasta = print_final_paths(
        &graph,
        &bfly.final_paths,
        &bfly.gene_ids,
        &params.name,
        ctx.kmer_size,
    );
    let _ = std::io::stderr().flush();
    Ok(ComponentResult {
        all_prob_paths_fasta: fasta,
        num_paths: bfly.final_paths.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{DiGraph, SeqVertex, SimpleEdge};

    #[test]
    fn junction_k25_normal_24mer_shift() {
        // 等长 24-mer 图, 重叠 23, 延伸 1 碱基
        let u = b"TGCAGTCCGATGACAATTCCGGTA";
        let v = b"GCAGTCCGATGACAATTCCGGTAA";
        assert_eq!(
            junction_kmer(u, v, 25).unwrap(),
            b"TGCAGTCCGATGACAATTCCGGTAA".to_vec()
        );
    }

    #[test]
    fn junction_k25_contained_v_is_none() {
        let u = b"AAAACCCCGGGGTTTT";
        let v = b"GGGGTTTT"; // 被 u 后缀包含
        assert!(junction_kmer(u, v, 25).is_none());
    }

    #[test]
    fn junction_k25_short_u_borrows_v() {
        // u 仅 8 碱基 → 借 v 延伸补足 25
        let u = b"ACGTACGT";
        let v = b"CGTAAAAAAAAAAAAAAAAAAAAA"; // ov=3, ext=21 个 A
        let k25 = junction_kmer(u, v, 25).unwrap();
        assert_eq!(k25.len(), 25);
        assert_eq!(&k25[..8], b"ACGTACGT");
        assert_eq!(&k25[8..], b"AAAAAAAAAAAAAAAAA"); // 17 个 A
    }

    #[test]
    fn junction_k25_insufficient_ext_is_none() {
        let u = b"ACGT";
        let v = b"CGTAA"; // ext 仅 2, 补不满 25
        assert!(junction_kmer(u, v, 25).is_none());
    }

    #[test]
    fn junction_k25_long_zipped_nodes() {
        // zipping 后变长节点: u 30nt, 重叠 25, 延伸 "TT"
        let u = b"GGGGAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let v = b"AAAAAAAAAAAAAAAAAAAAAAAAATT";
        let k25 = junction_kmer(u, v, 25).unwrap();
        assert_eq!(k25.len(), 25);
        assert_eq!(&k25[..24], b"AAAAAAAAAAAAAAAAAAAAAAAA");
        assert_eq!(&k25[24..], b"T");
    }

    #[test]
    fn fuse_read_evidence_untouched() {
        assert_eq!(fuse_edge_weight(5, 0, 100, 10), 5); // pog 证据
        assert_eq!(fuse_edge_weight(0, 3, 100, 10), 3); // PairPath 证据
        assert_eq!(fuse_edge_weight(1, 0, 100, 10), 1); // 证据=1 不被覆盖度顶掉
    }

    #[test]
    fn fuse_floor_by_cov_capped() {
        assert_eq!(fuse_edge_weight(0, 0, 100, 10), 10); // cap 生效
        assert_eq!(fuse_edge_weight(0, 0, 7, 10), 7); // 未达 cap
        assert_eq!(fuse_edge_weight(0, 0, 0, 10), 1); // 表缺 → floor 1（现状）
    }

    #[test]
    fn lookup_cov25_canonical() {
        let k25 = b"TGCAGTCCGATGACAATTCCGGTAA";
        let key = get_ds_kmer_val(kmer_to_intval(k25).unwrap(), 25);
        let mut table = trinity_kmer::counter::CountMap::default();
        table.insert(key, 9);
        assert_eq!(lookup_cov(&table, k25, 25), 9);
        // 反义查同一 canonical 键
        let rc: Vec<u8> = k25
            .iter()
            .rev()
            .map(|b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                x => *x,
            })
            .collect();
        assert_eq!(lookup_cov(&table, &rc, 25), 9);
        // 缺失 → 0; 非 gatc → 0
        assert_eq!(lookup_cov(&table, b"C", 25), 0);
        assert_eq!(lookup_cov(&table, b"N", 25), 0);
    }

    fn two_branch_graph() -> DiGraph {
        // 1→2 / 1→3 分叉, 三条 24-mer 名（ov=23 各延伸 1 碱基）
        let mut g = DiGraph::new();
        g.add_vertex(SeqVertex::new(1, "TGCAGTCCGATGACAATTCCGGTA"));
        g.add_vertex(SeqVertex::new(2, "GCAGTCCGATGACAATTCCGGTAA"));
        g.add_vertex(SeqVertex::new(3, "GCAGTCCGATGACAATTCCGGTAC"));
        g.add_edge(1, 2, SimpleEdge::new(0.0, 1, 2));
        g.add_edge(1, 3, SimpleEdge::new(0.0, 1, 3));
        g
    }

    fn cov_table() -> trinity_kmer::counter::CountMap {
        // 仅 (1→2) junction 25-mer 计数 100
        let k25 = b"TGCAGTCCGATGACAATTCCGGTAA";
        let key = get_ds_kmer_val(kmer_to_intval(k25).unwrap(), 25);
        let mut t = trinity_kmer::counter::CountMap::default();
        t.insert(key, 100);
        t
    }

    #[test]
    fn digraph_to_flow_none_matches_legacy() {
        let g = two_branch_graph();
        let fg = digraph_to_flow(&g, &FxHashMap::default(), None, 25);
        assert_eq!(fg.edges[&(1, 2)], 1);
        assert_eq!(fg.edges[&(1, 3)], 1);
    }

    #[test]
    fn digraph_to_flow_cov_fusion_fills_floor() {
        let g = two_branch_graph();
        let t = cov_table();
        let fg = digraph_to_flow(&g, &FxHashMap::default(), Some((&t, 10)), 25);
        assert_eq!(fg.edges[&(1, 2)], 10); // 覆盖度 100 → cap 10
        assert_eq!(fg.edges[&(1, 3)], 1); // 表缺 → floor 1
    }

    #[test]
    fn digraph_to_flow_read_evidence_wins_over_cov() {
        let mut g = DiGraph::new();
        g.add_vertex(SeqVertex::new(1, "TGCAGTCCGATGACAATTCCGGTA"));
        g.add_vertex(SeqVertex::new(2, "GCAGTCCGATGACAATTCCGGTAA"));
        g.add_edge(1, 2, SimpleEdge::new(5.0, 1, 2));
        let t = cov_table();
        let fg = digraph_to_flow(&g, &FxHashMap::default(), Some((&t, 10)), 25);
        assert_eq!(fg.edges[&(1, 2)], 5); // read 证据不被覆盖度压过
    }
}
