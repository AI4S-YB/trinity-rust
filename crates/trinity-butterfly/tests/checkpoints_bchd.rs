//! 检查点 B/C/H/D 对拍：剪枝链各阶段 DOT 与 Java Butterfly 参照比较。
//!
//! 参照由内层 Butterfly.jar 生成（--generate_intermediate_dot_files）：
//! - B = removeLightEdges 后（DFS 未跑）
//! - C = compactLinearPaths + My_DFS2 后
//! - H = removeSingleNtBubbles 后（COLLAPSE_SNPs 默认 true）
//! - D = calcSubComponentsStats 后
//!
//! My_DFS 已移植（T3）：node depth（_D）全字段对拍一致；只有 discovery/finish
//! time（[T:..]）与 Java 不同——JUNG getVertices() 是 HashSet 任意序，我们的
//! 迭代序是插入序，DFS 森林访问顺序不同导致 time 编号不同（depth 不受影响）。
//! 因此跨实现比较仍剥 [T:]；另存我们自己的 C 快照（含 [T:]）做全字段自回归。

use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

use trinity_butterfly::{
    build_new_graph_use_kmers, calc_sub_components_stats, compact_linear_paths,
    fix_extremely_high_single_edges, pre_process_graph_file, remove_light_edges,
    remove_single_nt_bubbles, run_dfs2, write_dot_string,
};

fn repo_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
}

/// 剥离 [T:...]（discovery time，迭代序相关，与 Java 不可对拍；_D depth 不剥）。
fn strip_dfs_fields(line: &str) -> String {
    let mut s = line.to_string();
    if let Some(pos) = s.find("[T:") {
        if let Some(end) = s[pos..].find(']') {
            s = format!("{}{}", &s[..pos], &s[pos + end + 1..]);
        }
    }
    s
}

fn normalize_dot(dot: &str) -> (BTreeSet<String>, BTreeSet<String>) {
    let mut nodes = BTreeSet::new();
    let mut edges = BTreeSet::new();
    for line in dot.lines() {
        let line = line.trim();
        if line.is_empty() || line == "digraph G {" || line == "}" {
            continue;
        }
        if line.contains("->") {
            edges.insert(line.to_string());
        } else {
            nodes.insert(strip_dfs_fields(line));
        }
    }
    (nodes, edges)
}

struct Stage {
    name: &'static str,
    fixture: &'static str,
}

fn compare_checkpoint(stage: Stage, our_dot: &str) {
    let orig = fs::read_to_string(repo_root().join(format!("fixtures/p4/c0/{}", stage.fixture)))
        .unwrap_or_else(|e| panic!("读取 {} 参照失败: {}", stage.fixture, e));
    let (our_n, our_e) = normalize_dot(our_dot);
    let (orig_n, orig_e) = normalize_dot(&orig);
    let nd: Vec<_> = our_n.symmetric_difference(&orig_n).collect();
    let ed: Vec<_> = our_e.symmetric_difference(&orig_e).collect();
    assert!(
        nd.is_empty() && ed.is_empty(),
        "[{}] 节点差异({} vs {}): {:#?}\n边差异({} vs {}): {:#?}",
        stage.name,
        our_n.len(),
        orig_n.len(),
        nd,
        our_e.len(),
        orig_e.len(),
        ed,
    );
}

#[test]
fn checkpoints_b_c_h_d_match_java() {
    let graph_text =
        fs::read_to_string(repo_root().join("fixtures/p3/quantify/c0/orig.graph.out")).unwrap();
    let br = build_new_graph_use_kmers(&graph_text).unwrap();
    let (mut graph, mut ctx) = (br.graph, br.ctx);

    let (in_flow, out_flow, _) = pre_process_graph_file(&graph_text);

    // B：fix → removeLight
    fix_extremely_high_single_edges(&mut graph, &out_flow, &in_flow);
    remove_light_edges(&mut graph, &ctx);
    compare_checkpoint(
        Stage {
            name: "B",
            fixture: "c0.graph_removeLightEdges_init.B.dot",
        },
        &write_dot_string(&graph, ctx.kmer_size, false),
    );

    // C：compact + My_DFS2（runDFS2：双向 DFS + down-only 重算 depth）
    compact_linear_paths(&mut graph, &ctx);
    run_dfs2(&mut graph);
    let c_dot = write_dot_string(&graph, ctx.kmer_size, false);
    compare_checkpoint(
        Stage {
            name: "C",
            fixture: "c0.graph_compactLinearPaths_init.C.dot",
        },
        &c_dot,
    );

    // H：SNP bubbles（COLLAPSE_SNPs 默认 true）
    remove_single_nt_bubbles(&mut graph, &mut ctx);
    compare_checkpoint(
        Stage {
            name: "H",
            fixture: "c0.graph_SNPs_removed.H.dot",
        },
        &write_dot_string(&graph, ctx.kmer_size, false),
    );

    // D：calcSubComponents（-L 200 → min_output_seq=200）
    calc_sub_components_stats(&mut graph, 200);
    compare_checkpoint(
        Stage {
            name: "D",
            fixture: "c0.graph_compactLinearPaths_removeSmallComp.D.dot",
        },
        &write_dot_string(&graph, ctx.kmer_size, false),
    );
}

#[test]
fn strip_dfs_fields_examples() {
    assert_eq!(
        strip_dfs_fields("4651 [label=\"X:W110(V4651_D1)[L:173][T:4]\" ,style=bold]"),
        "4651 [label=\"X:W110(V4651_D1)[L:173]\" ,style=bold]"
    );
    assert_eq!(
        strip_dfs_fields("5645 [label=\"ABC:W-1(V5645_D0)[L:24][T:1]\"]"),
        "5645 [label=\"ABC:W-1(V5645_D0)[L:24]\"]"
    );
    // B 检查点行（DFS 未跑）只剥 [T:-1]，_D-1 保留
    let stripped = strip_dfs_fields("1 [label=\"A:W5(V1_D-1)[L:1][T:-1]\"]");
    assert_eq!(stripped, "1 [label=\"A:W5(V1_D-1)[L:1]\"]");
}

/// T3 复验（全字段自回归）：我们自己的 post-My_DFS2 C 快照（含 [T:]）。
/// 防 DFS/depth 逻辑未来回归；快照本身已与 Java 的 _D 值一致（见上）。
#[test]
fn checkpoint_c_full_field_self_snapshot() {
    let graph_text =
        fs::read_to_string(repo_root().join("fixtures/p3/quantify/c0/orig.graph.out")).unwrap();
    let br = build_new_graph_use_kmers(&graph_text).unwrap();
    let (mut graph, ctx) = (br.graph, br.ctx);
    let (in_flow, out_flow, _) = pre_process_graph_file(&graph_text);
    fix_extremely_high_single_edges(&mut graph, &out_flow, &in_flow);
    remove_light_edges(&mut graph, &ctx);
    compact_linear_paths(&mut graph, &ctx);
    run_dfs2(&mut graph);
    let our_dot = write_dot_string(&graph, ctx.kmer_size, false);
    let snapshot =
        fs::read_to_string(repo_root().join("fixtures/p4/c0/c0.rust_post_dfs.C.dot")).unwrap();
    assert_eq!(our_dot, snapshot);
}
