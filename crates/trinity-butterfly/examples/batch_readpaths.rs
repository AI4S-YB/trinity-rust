//! 批量探针：B4 组件逐个跑剪枝链+穿线，dump read paths + orig id 映射。
//!
//! 用法: batch_readpaths <probe_dir>
//!   <probe_dir>/comps.txt = 组件号清单; <c>.graph/<c>.reads = 输入
//!   输出 <c>.paths（name \t v1,v2,...）与 <c>.map（orig_id \t merged_id）。
//! rayon 并行（组件独立）。

use std::path::Path;

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use trinity_butterfly::graph_io::{
    build_new_graph_use_kmers, parse_graph_reads, pre_process_graph_file,
};
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
    let dir = Path::new(&args[1]);
    let comps: Vec<u64> = std::fs::read_to_string(dir.join("comps.txt"))
        .unwrap()
        .lines()
        .filter_map(|l| l.trim().parse().ok())
        .collect();
    eprintln!("{} comps", comps.len());

    let fail: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    comps.par_iter().for_each(|c| {
        let graph_text = match std::fs::read_to_string(dir.join(format!("{c}.graph"))) {
            Ok(t) => t,
            Err(_) => return,
        };
        let reads_raw = match std::fs::read(dir.join(format!("{c}.reads"))) {
            Ok(t) => t,
            Err(_) => return,
        };
        let res = run_one(&graph_text, &reads_raw);
        let (paths, map_text) = match res {
            Some(v) => v,
            None => {
                fail.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                return;
            }
        };
        let _ = std::fs::write(dir.join(format!("{c}.paths")), paths);
        let _ = std::fs::write(dir.join(format!("{c}.map")), map_text);
    });
    eprintln!(
        "done, failures: {}",
        fail.load(std::sync::atomic::Ordering::Relaxed)
    );
}

fn run_one(graph_text: &str, reads_raw: &[u8]) -> Option<(String, String)> {
    let br0 = build_new_graph_use_kmers(graph_text).ok()?;
    let mut orig_kmer_to_node: FxHashMap<String, i32> = FxHashMap::default();
    for &vid in br0.graph.vertex_ids() {
        orig_kmer_to_node.insert(br0.graph.get_vertex(vid).unwrap().name.clone(), vid);
    }
    let br = build_new_graph_use_kmers(graph_text).ok()?;
    let (mut graph, mut ctx): (DiGraph, BflyContext) = (br.graph, br.ctx);
    let (in_flow, out_flow, _) = pre_process_graph_file(graph_text);
    fix_extremely_high_single_edges(&mut graph, &out_flow, &in_flow);
    remove_light_edges(&mut graph, &ctx);
    compact_linear_paths(&mut graph, &ctx);
    run_dfs2(&mut graph);
    remove_single_nt_bubbles(&mut graph, &mut ctx);
    calc_sub_components_stats(&mut graph, 200);

    let orig_id_map = get_original_ver_ids_mapping_hash(&mut graph, ctx.last_real_id);
    let raw_reads =
        parse_graph_reads(String::from_utf8_lossy(reads_raw).as_ref(), ctx.kmer_size, false)
            .ok()?;
    let ordered =
        get_read_starts_ordered(&graph, &ctx, &orig_id_map, &orig_kmer_to_node, &raw_reads);
    let paths = format_read_paths_dump(&ordered);
    let mut map_text = String::new();
    for (orig, loc) in &orig_id_map {
        map_text.push_str(&format!("{}\t{}\n", orig, loc.node_id));
    }
    Some((paths, map_text))
}
