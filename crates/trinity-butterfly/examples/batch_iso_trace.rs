//! 终版探针：完整跑到分量判定, 定位孤立链消失环节。
//!
//! 用法: batch_iso_trace <probe_dir>
//! 输出:
//!   <c>.comp  每行: 分量节点 id（逗号） \t reads 数
//!   <c>.iso   每行: POG 孤立节点的 revised new ids（逗号） \t orig ids（逗号）
//!   <c>.readhash 每行: read 路径（重映射后 new ids 逗号） \t support

use std::path::Path;

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use trinity_butterfly::graph_io::{
    build_new_graph_use_kmers, parse_graph_reads, pre_process_graph_file,
};
use trinity_butterfly::pair_paths::get_suff_stats_w_pairs;
use trinity_butterfly::paths::{
    divide_into_components, get_component_reads, reorganize_read_pairings,
    relabel_edge_weights_using_orig_kmers,
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
    let dir = Path::new(&args[1]);
    let comps: Vec<u64> = std::fs::read_to_string(dir.join("comps.txt"))
        .unwrap()
        .lines()
        .filter_map(|l| l.trim().parse().ok())
        .collect();
    eprintln!("{} comps", comps.len());
    let fail = std::sync::atomic::AtomicUsize::new(0);
    comps.par_iter().for_each(|c| {
        let graph_text = match std::fs::read_to_string(dir.join(format!("{c}.graph"))) {
            Ok(t) => t,
            Err(_) => return,
        };
        let reads_raw = match std::fs::read(dir.join(format!("{c}.reads"))) {
            Ok(t) => t,
            Err(_) => return,
        };
        match run_one(&graph_text, &reads_raw) {
            Some((comp_t, iso_t, rh_t)) => {
                let _ = std::fs::write(dir.join(format!("{c}.comp")), comp_t);
                let _ = std::fs::write(dir.join(format!("{c}.iso")), iso_t);
                let _ = std::fs::write(dir.join(format!("{c}.readhash")), rh_t);
            }
            None => {
                fail.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
        }
    });
    eprintln!("done, failures: {}", fail.load(std::sync::atomic::Ordering::Relaxed));
}

fn run_one(graph_text: &str, reads_raw: &[u8]) -> Option<(String, String, String)> {
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
    let suff = get_suff_stats_w_pairs(&graph, &ordered);

    let mut ctx2 = ctx.clone();
    let result = create_dag_from_overlap_layout(&graph, &mut ctx2, &suff);

    // paths.rs 前半（分量判定）
    let mut sg = result.seqvertex_graph.clone();
    run_dfs2(&mut sg);
    let read_hash = reorganize_read_pairings(&sg, &result.combined_read_hash);
    relabel_edge_weights_using_orig_kmers(&mut sg, &br0.graph);

    let mut comp_text = String::new();
    for comp in divide_into_components(&sg) {
        let cr = get_component_reads(&comp, &read_hash);
        let ids: Vec<String> = comp.iter().map(|x| x.to_string()).collect();
        let nreads: usize = cr.values().map(|m| m.len()).sum();
        comp_text.push_str(&ids.join(","));
        comp_text.push('\t');
        comp_text.push_str(&nreads.to_string());
        comp_text.push('\n');
    }

    // 孤立 POG 节点（.pogedges 已知; 这里直接按 POG 度重新算更稳）
    let mut iso_text = String::new();
    for pwo in &result.revised_paths {
        iso_text.push_str(
            &pwo.vertex_id_list
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(","),
        );
        iso_text.push('\t');
        iso_text.push_str(
            &pwo.orig_vertex_id_list
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(","),
        );
        iso_text.push('\n');
    }

    // 重映射后 read 路径
    let mut rh_text = String::new();
    for (_k, inner) in read_hash.iter() {
        for (pp, s) in inner {
            let mut ids: Vec<i32> = pp.path1.clone();
            if pp.has_second_path() {
                ids.extend(pp.path2.iter().copied());
            }
            rh_text.push_str(
                &ids.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","),
            );
            rh_text.push('\t');
            rh_text.push_str(&s.to_string());
            rh_text.push('\n');
        }
    }
    Some((comp_text, iso_text, rh_text))
}
