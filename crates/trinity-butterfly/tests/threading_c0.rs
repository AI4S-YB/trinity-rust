//! c0 组件全链穿线对拍：build → 剪枝链(B/C/H/D) → getOriginalVerIDsMappingHash →
//! getReadStarts，与原版 Butterfly.jar 的 `-V 17 --stderr` 逐 read 路径黄金比对。
//!
//! 黄金生成见 fixtures/p4/c0/threading_README.md。

use std::fs;
use std::path::Path;

use rustc_hash::FxHashMap;
use trinity_butterfly::{
    build_new_graph_use_kmers, calc_sub_components_stats, compact_linear_paths,
    fix_extremely_high_single_edges, format_read_paths_dump, get_original_ver_ids_mapping_hash,
    get_read_starts_ordered, parse_graph_reads, pre_process_graph_file, remove_light_edges,
    remove_single_nt_bubbles, run_dfs2,
};

fn repo_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
}

#[test]
fn c0_threading_matches_java_per_read() {
    let graph_text =
        fs::read_to_string(repo_root().join("fixtures/p3/quantify/c0/orig.graph.out")).unwrap();
    let reads_text =
        fs::read_to_string(repo_root().join("fixtures/p3/quantify/c0/orig.reads.out")).unwrap();

    // 初始图（剪枝前）建 kmer → node id 映射（originalGraphKmerToNodeID）
    let br0 = build_new_graph_use_kmers(&graph_text).unwrap();
    let mut orig_kmer_to_node: FxHashMap<String, i32> = FxHashMap::default();
    for &vid in br0.graph.vertex_ids() {
        orig_kmer_to_node.insert(br0.graph.get_vertex(vid).unwrap().name.clone(), vid);
    }

    // 主链：B/C/H/D（与 checkpoints_bchd 相同的管道）
    let br = build_new_graph_use_kmers(&graph_text).unwrap();
    let (mut graph, mut ctx) = (br.graph, br.ctx);
    let (in_flow, out_flow, _) = pre_process_graph_file(&graph_text);
    fix_extremely_high_single_edges(&mut graph, &out_flow, &in_flow);
    remove_light_edges(&mut graph, &ctx);
    compact_linear_paths(&mut graph, &ctx);
    run_dfs2(&mut graph);
    remove_single_nt_bubbles(&mut graph, &mut ctx);
    calc_sub_components_stats(&mut graph, 200);

    let orig_id_map = get_original_ver_ids_mapping_hash(&mut graph, ctx.last_real_id);

    let raw_reads = parse_graph_reads(&reads_text, ctx.kmer_size, false).unwrap();
    let ordered =
        get_read_starts_ordered(&graph, &ctx, &orig_id_map, &orig_kmer_to_node, &raw_reads);

    // 统计对拍：4336 mapped / 4342 total（threading_stats.txt）
    assert_eq!(ordered.len(), 4336, "穿线成功 read 数应与 jar 一致");
    let names: std::collections::HashSet<&str> = ordered.iter().map(|(n, _)| n.as_str()).collect();
    assert_eq!(names.len(), 2423, "不同 read 名数应与 jar pairs 统计一致");

    // 逐 read 路径对拍（文件出现序）
    let dump = format_read_paths_dump(&ordered);
    let golden = fs::read_to_string(repo_root().join("fixtures/p4/c0/threading_golden.tsv"))
        .unwrap()
        .replace("\r", "");
    let dump_lines: Vec<&str> = dump.lines().collect();
    let golden_lines: Vec<&str> = golden.lines().collect();
    assert_eq!(dump_lines.len(), golden_lines.len());

    let mut diff = 0;
    let mut first_diffs: Vec<String> = Vec::new();
    for (i, (ours, theirs)) in dump_lines.iter().zip(golden_lines.iter()).enumerate() {
        if ours != theirs {
            diff += 1;
            if first_diffs.len() < 10 {
                first_diffs.push(format!("line {i}: rust={ours:?} java={theirs:?}"));
            }
        }
    }
    // JUNG 后继迭代序为 HashSet 任意序：平局(<=)取末的分支选择可能不同。
    // 允许极小比例差异，但要求 >99.5% 逐字节一致。
    let pct = diff as f64 / golden_lines.len() as f64;
    assert!(
        pct <= 0.005,
        "路径差异 {diff}/{} ({pct:.4}) 超过 0.5%：\n{}",
        golden_lines.len(),
        first_diffs.join("\n")
    );
    eprintln!("threading diff: {diff}/{} ({pct:.4})", golden_lines.len());
}
