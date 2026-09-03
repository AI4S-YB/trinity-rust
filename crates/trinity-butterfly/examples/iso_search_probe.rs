//! 单案例深挖：孤立链分量的完整搜索实验。
//!
//! 用法: iso_search_probe <graph.out> <reads>
//! 全链跑到分量, 对每分量（含孤立链的）跑 add_s_and_t + 搜索,
//! dump 搜索返回的 final_paths 与被过滤原因的粗分类。

use std::collections::HashSet;

use rustc_hash::FxHashMap;

use trinity_butterfly::graph_io::{
    build_new_graph_use_kmers, parse_graph_reads, pre_process_graph_file,
};
use trinity_butterfly::pair_paths::get_suff_stats_w_pairs;
use trinity_butterfly::paths::{
    add_s_and_t, divide_into_components, extract_complex_path_prefixes_from_reads,
    extract_triplets_from_reads, get_all_probable_paths,
    get_component_reads, get_xstructures_resolved_by_triplets,
    reduce_to_max_paths_per_node, reorganize_read_pairings,
    relabel_edge_weights_using_orig_kmers, PathSearchParams,
};
use trinity_butterfly::pog::create_dag_from_overlap_layout;
use trinity_butterfly::prune::{
    calc_sub_components_stats, compact_linear_paths, fix_extremely_high_single_edges,
    remove_light_edges, remove_single_nt_bubbles,
};
use trinity_butterfly::threading::{
    get_original_ver_ids_mapping_hash, get_read_starts_ordered,
};
use trinity_butterfly::{run_dfs2, BflyContext, DiGraph};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let graph_text = std::fs::read_to_string(&args[1]).unwrap();
    let reads_raw = std::fs::read(&args[2]).unwrap();

    let br0 = build_new_graph_use_kmers(&graph_text).unwrap();
    let mut orig_kmer_to_node: FxHashMap<String, i32> = FxHashMap::default();
    for &vid in br0.graph.vertex_ids() {
        orig_kmer_to_node.insert(br0.graph.get_vertex(vid).unwrap().name.clone(), vid);
    }
    let br = build_new_graph_use_kmers(&graph_text).unwrap();
    let (mut graph, mut ctx): (DiGraph, BflyContext) = (br.graph, br.ctx);
    let (in_flow, out_flow, _) = pre_process_graph_file(&graph_text);
    fix_extremely_high_single_edges(&mut graph, &out_flow, &in_flow);
    remove_light_edges(&mut graph, &ctx);
    compact_linear_paths(&mut graph, &ctx);
    run_dfs2(&mut graph);
    remove_single_nt_bubbles(&mut graph, &mut ctx);
    calc_sub_components_stats(&mut graph, 200);

    let orig_id_map = get_original_ver_ids_mapping_hash(&mut graph, ctx.last_real_id);
    let raw_reads =
        parse_graph_reads(String::from_utf8_lossy(&reads_raw).as_ref(), ctx.kmer_size, false)
            .unwrap();
    let ordered =
        get_read_starts_ordered(&graph, &ctx, &orig_id_map, &orig_kmer_to_node, &raw_reads);
    let suff = get_suff_stats_w_pairs(&graph, &ordered);
    let mut ctx2 = ctx.clone();
    let result = create_dag_from_overlap_layout(&graph, &mut ctx2, &suff);

    let mut sg = result.seqvertex_graph.clone();
    run_dfs2(&mut sg);
    let read_hash = reorganize_read_pairings(&sg, &result.combined_read_hash);
    relabel_edge_weights_using_orig_kmers(&mut sg, &br0.graph);

    let sp = PathSearchParams {
        kmer_size: ctx.kmer_size,
        min_output_seq: std::env::var("PROBE_L").ok().and_then(|v| v.parse().ok()).unwrap_or(200),
        min_read_support_thr: std::env::var("PROBE_R").ok().and_then(|v| v.parse().ok()).unwrap_or(2),
        path_reinforcement_distance: std::env::var("PROBE_PRD").ok().and_then(|v| v.parse().ok()).unwrap_or(25),
        max_paths_per_node_init: 100_000,
        ..PathSearchParams::default()
    };

    for comp in divide_into_components(&sg) {
        let mut comp_reads = get_component_reads(&comp, &read_hash);
        let bases: usize = comp
            .iter()
            .filter_map(|v| sg.get_vertex(*v).map(|x| x.name.len()))
            .sum();
        eprintln!(
            "comp: {} vertices (bases ~{}), {} reads-paths",
            comp.len(),
            bases,
            comp_reads.values().map(|m| m.len()).sum::<usize>()
        );
        if comp_reads.is_empty() {
            eprintln!("  -> skipped (empty reads)");
            continue;
        }
        reduce_to_max_paths_per_node(&mut comp_reads, sp.max_paths_per_node_init);
        let triplets = extract_triplets_from_reads(&comp_reads);
        let x_structures = get_xstructures_resolved_by_triplets(&sg, &comp, &triplets);
        let extended = extract_complex_path_prefixes_from_reads(&comp_reads);
        let mut graph2 = sg.clone();
        add_s_and_t(&mut graph2, &comp, &mut comp_reads, sp.min_read_support_thr);
        let r = get_all_probable_paths(
            &graph2,
            &comp,
            &comp_reads,
            &triplets,
            &extended,
            &x_structures,
            &sp,
        );
        eprintln!(
            "  search -> {} final paths (lens {:?})",
            r.final_paths.len(),
            r.final_paths.iter().map(|p| p.len()).take(8).collect::<Vec<_>>()
        );
        // 序列长度（粗算: 节点数 × ~1 + k）
        for p in r.final_paths.iter().take(5) {
            let total: usize = p
                .iter()
                .filter_map(|v| graph2.get_vertex(*v).map(|x| x.name.len()))
                .sum();
            eprintln!("    path len nodes={} bases~{}", p.len(), total);
        }
        // comp 顶点集合是否全被某 final path 覆盖
        let covered: HashSet<i32> = r.final_paths.iter().flatten().copied().collect();
        let inter = comp.iter().filter(|v| covered.contains(v)).count();
        eprintln!(
            "  coverage: {}/{} comp vertices in final paths",
            inter,
            comp.len()
        );
    }
}
