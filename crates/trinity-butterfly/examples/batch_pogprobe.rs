//! 批量 POG 探针：B4 组件跑至 POG 破环层，dump 拼接证据。
//!
//! 用法: batch_pogprobe <probe_dir>
//! 输出（每组件）:
//!   <c>.pog       非包含 read path 集（POG 节点）: path \t support
//!   <c>.cont      containment: 被含 path \t container;...;
//!   <c>.pogedges  overlap 边 j->i（j,i 为 .pog 行号 0-based）+ 破环删除标记
//!                 行格式: j \t i \t weight \t [del] （del=破环删除）

use std::path::Path;

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use trinity_butterfly::graph_io::{
    build_new_graph_use_kmers, parse_graph_reads, pre_process_graph_file,
};
use trinity_butterfly::pair_paths::get_suff_stats_w_pairs;
use trinity_butterfly::pog::{
    break_cycles_in_path_overlap_graph, construct_path_overlap_graph,
    find_dispersed_repeat_nodes, populate_pairpaths_and_readsupport, remove_containments,
};
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
            Some((pog, cont, edges)) => {
                let _ = std::fs::write(dir.join(format!("{c}.pog")), pog);
                let _ = std::fs::write(dir.join(format!("{c}.cont")), cont);
                let _ = std::fs::write(dir.join(format!("{c}.pogedges")), edges);
            }
            None => {
                fail.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
        }
    });
    eprintln!("done, failures: {}", fail.load(std::sync::atomic::Ordering::Relaxed));
}

fn run_one(graph_text: &str, reads_raw: &[u8]) -> Option<(String, String, String)> {
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
    // 锚点 kmer 表（穿线再锚定用）
    let br0 = build_new_graph_use_kmers(graph_text).ok()?;
    let mut orig_kmer_to_node: FxHashMap<String, i32> = FxHashMap::default();
    for &vid in br0.graph.vertex_ids() {
        orig_kmer_to_node.insert(br0.graph.get_vertex(vid).unwrap().name.clone(), vid);
    }
    let raw_reads =
        parse_graph_reads(String::from_utf8_lossy(reads_raw).as_ref(), ctx.kmer_size, false)
            .ok()?;
    let ordered =
        get_read_starts_ordered(&graph, &ctx, &orig_id_map, &orig_kmer_to_node, &raw_reads);
    let suff = get_suff_stats_w_pairs(&graph, &ordered);

    // POG 前半（镜像 create_dag_from_overlap_layout）
    let pair_path_to_read_support = populate_pairpaths_and_readsupport(&suff.combined_read_hash);
    let support_map: FxHashMap<_, _> = pair_path_to_read_support.iter().cloned().collect();
    let mut paths: Vec<Vec<i32>> = Vec::new();
    for (pp, _) in &pair_path_to_read_support {
        paths.push(pp.path1.clone());
        if pp.has_second_path() {
            paths.push(pp.path2.clone());
        }
    }
    paths.sort_by_key(|p| p.len());
    paths.reverse();

    let mut contained_map: FxHashMap<Vec<i32>, Vec<Vec<i32>>> = FxHashMap::default();
    let noncontained = remove_containments(&paths, &mut contained_map);

    let mut pog_text = String::new();
    let mut path_support = FxHashMap::default();
    for (pp, s) in &pair_path_to_read_support {
        path_support.insert(pp.path1.clone(), *s);
        if pp.has_second_path() {
            path_support.insert(pp.path2.clone(), *s);
        }
    }
    for p in &noncontained {
        let sup = path_support.get(p).copied().unwrap_or(0);
        let ids: Vec<String> = p.iter().map(|x| x.to_string()).collect();
        pog_text.push_str(&ids.join(","));
        pog_text.push('\t');
        pog_text.push_str(&sup.to_string());
        pog_text.push('\n');
    }
    let mut cont_text = String::new();
    let mut cont_keys: Vec<&Vec<i32>> = contained_map.keys().collect();
    cont_keys.sort_by_key(|p| p.len());
    for k in cont_keys {
        let ids: Vec<String> = k.iter().map(|x| x.to_string()).collect();
        cont_text.push_str(&ids.join(","));
        cont_text.push('\t');
        let cs: Vec<String> = contained_map[k]
            .iter()
            .map(|c| c.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","))
            .collect();
        cont_text.push_str(&cs.join(";"));
        cont_text.push('\n');
    }

    let dispersed = find_dispersed_repeat_nodes(&noncontained);
    let mut path_matches: FxHashMap<(usize, usize), _> = FxHashMap::default();
    let mut pog = construct_path_overlap_graph(&noncontained, &mut path_matches, &dispersed);
    let edges_before: Vec<(usize, usize)> = {
        let mut v: Vec<_> = pog.edges.keys().copied().collect();
        v.sort_unstable();
        v
    };
    loop {
        if !break_cycles_in_path_overlap_graph(&mut pog) {
            break;
        }
    }
    let mut edge_text = String::new();
    for (j, i) in edges_before {
        let w = path_matches
            .get(&(j, i))
            .map(|po| po.match_score)
            .unwrap_or(0);
        let del = if pog.edges.contains_key(&(j, i)) { "" } else { "del" };
        edge_text.push_str(&format!("{j}\t{i}\t{w}\t{del}\n"));
    }
    let _ = support_map;
    Some((pog_text, cont_text, edge_text))
}
