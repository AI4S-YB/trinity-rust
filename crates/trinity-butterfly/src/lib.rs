//! P4 Butterfly: de Bruijn 图基础（DiGraph / SeqVertex / SimpleEdge / BflyContext / graph IO / DOT 写出）。
//!
//! 镜像 trinityrnaseq v2.15.2 `Butterfly/Butterfly/src/src/{SeqVertex,SimpleEdge}.java`
//! 与 `TransAssembly_allProbPaths.java` 的 preProcessGraphFile / buildNewGraphUseKmers /
//! writeDotFile / getReadStarts(readAndMapSingleRead 的解析部分)。

pub mod align;
pub mod component;
pub mod context;
pub mod dfs;
pub mod graph;
pub mod graph_io;
pub mod pair_paths;
pub mod paths;
pub mod pog;
pub mod postprocess;
pub mod prune;
pub mod threading;

pub use component::{run_component, ComponentParams, ComponentResult};
pub use context::BflyContext;
pub use dfs::{get_topological_order, run_dfs2};
pub use graph::{DiGraph, SeqVertex, SimpleEdge, T_VERTEX_ID, VERTEX_ROOT_ID};
pub use graph_io::{
    build_new_graph_use_kmers, parse_graph_reads, pre_process_graph_file, write_dot_file,
    write_dot_string, BuildResult, RawRead, ReadParseError,
};
pub use pair_paths::{
    get_suff_stats_w_pairs, individual_paths_are_compatible, trim_sink_nodes, PairPath, SuffStats,
};
pub use paths::{
    add_s_and_t, divide_into_components, extract_complex_path_prefixes_from_reads,
    extract_triplets_from_reads, get_all_probable_paths, get_component_reads, get_distance_wo_ver,
    get_xstructures_resolved_by_triplets, is_ancestral, paths_to_fasta,
    reduce_to_max_paths_per_node, relabel_edge_weights_using_orig_kmers,
    remove_all_edges_of_s_and_t, reorganize_read_pairings, run_butterfly_all_prob_paths,
    run_butterfly_all_prob_paths_with, ButterflyResult, PathSearchParams, PathSearchResult,
};
pub use pog::{
    break_cycles_in_path_overlap_graph, create_dag_from_overlap_layout,
    find_dispersed_repeat_nodes, get_all_possible_updated_path_mappings,
    path_a_contains_path_b_allow_repeats, path_b_extends_path_a_allow_repeats,
    populate_pairpaths_and_readsupport, remove_containments, topo_sort_seq_vertices_dag,
    OverlapLayoutResult, PathNode, PathOverlap, PathOverlapGraph, PathWithOrig, ZipState,
};
pub use postprocess::{
    assign_compatible_reads_to_paths, convert_to_orig_ids, get_path_name_string,
    group_paths_into_genes, java_hashmap_order, java_list_string, print_final_paths,
    reduce_cdhit_like, run_em_reduce, run_path_expression_em, two_paths_are_too_similar,
    ContainedReads, PostProcessParams,
};
pub use prune::{
    calc_sub_components_stats, compact_linear_paths, fix_extremely_high_single_edges,
    remove_light_edges, remove_light_flow_edges, remove_light_in_edges, remove_light_out_edges,
    remove_single_nt_bubbles, run_pruning_chain,
};
pub use threading::{
    find_path_in_graph, format_read_paths_dump, get_original_ver_ids_mapping_hash, get_read_starts,
    get_read_starts_ordered, LocInGraph, ReadPath, ThreadingCtx,
};

/// Butterfly 级别的错误（graph IO 等）。
#[derive(Debug, thiserror::Error)]
pub enum BflyError {
    #[error("graph 文件格式错误: {0}")]
    GraphFileFormat(String),
}

pub type Result<T> = std::result::Result<T, BflyError>;
