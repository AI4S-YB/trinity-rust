//! 探针：单组件穿线 dump——验证断口横跨 read 的路径形态。
//!
//! 用法: split_bridge_probe <graph.out> <graph.reads>
//! 输出: 每条成功穿线 read 的 path + mismatch_count（stdout）。

use rustc_hash::FxHashMap;

use trinity_butterfly::graph_io::{build_new_graph_use_kmers, parse_graph_reads, pre_process_graph_file};
use trinity_butterfly::prune::{
    calc_sub_components_stats, compact_linear_paths, fix_extremely_high_single_edges,
    remove_light_edges, remove_single_nt_bubbles,
};
use trinity_butterfly::threading::{
    format_read_paths_dump, get_original_ver_ids_mapping_hash, get_read_starts_ordered,
};
use trinity_butterfly::{run_dfs2, BflyContext, DiGraph};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let graph_text = std::fs::read_to_string(&args[1]).unwrap();
    let reads_text = std::fs::read(&args[2]).unwrap();

    // 镜像 run_component_inner 的剪枝链（至穿线前）
    let br0 = build_new_graph_use_kmers(&graph_text).unwrap();
    let mut orig_kmer_to_node: FxHashMap<String, i32> = FxHashMap::default();
    for &vid in br0.graph.vertex_ids() {
        orig_kmer_to_node.insert(
            br0.graph.get_vertex(vid).unwrap().name.clone(),
            vid,
        );
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
    let raw_reads = parse_graph_reads(
        String::from_utf8_lossy(&reads_text).as_ref(),
        ctx.kmer_size,
        false,
    )
    .unwrap();
    eprintln!(
        "graph: {} vertices (last_real_id {}), raw reads: {}",
        graph.vertex_ids().len(),
        ctx.last_real_id,
        raw_reads.len()
    );
    let ordered =
        get_read_starts_ordered(&graph, &ctx, &orig_id_map, &orig_kmer_to_node, &raw_reads);
    eprintln!("threaded read paths: {}", ordered.len());
    print!("{}", format_read_paths_dump(&ordered));

    // 剪枝后图节点组成: 合并节点 (id>last_real_id 的 orig 映射来源) + name 长度
    let mut merged_of: FxHashMap<i32, Vec<i32>> = FxHashMap::default();
    for (orig, loc) in &orig_id_map {
        merged_of.entry(loc.node_id).or_default().push(*orig);
    }
    let mut ids = graph.vertex_ids().to_vec();
    ids.sort_unstable();
    eprintln!("--- pruned graph vertices ---");
    for vid in ids {
        let v = graph.get_vertex(vid).unwrap();
        let origs = merged_of.get(&vid).cloned().unwrap_or_default();
        let mut origs = origs;
        origs.sort_unstable();
        eprintln!(
            "v{} name_len={} origs={:?}",
            vid,
            v.name.len(),
            if vid <= ctx.last_real_id { vec![vid] } else { origs }
        );
    }

    // 横跨 read 检查: 单节点链内路径的 mm 分布（尾部 mm 巨大 = 横跨特征）
    let big_mm = ordered
        .iter()
        .filter(|(_, p)| p.mismatch_count > 20)
        .count();
    eprintln!("paths with mismatch_count > 20: {}", big_mm);
}
