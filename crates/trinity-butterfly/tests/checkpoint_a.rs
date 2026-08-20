//! 检查点 A 对拍：buildNewGraphUseKmers 后的 DOT 与 Java Butterfly 的
//! `c0.graph_deBruijn.A.dot` 比较（节点/边集合；剥离 [T:...] discovery time 字段，
//! 因为此阶段 DFS 尚未运行，Java 侧是构造器默认值 -1）。

use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

use trinity_butterfly::{build_new_graph_use_kmers, parse_graph_reads, write_dot_string};

fn repo_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
}

/// 把 DOT 行分成 (节点行, 边行)，并对节点行剥离 [T:..] 字段。
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
            // 剥离 "[T:-1]" / "[T:12]"（DFS 未跑时无意义）
            let stripped = strip_discovery_time(line);
            nodes.insert(stripped);
        }
    }
    (nodes, edges)
}

fn strip_discovery_time(line: &str) -> String {
    if let Some(pos) = line.find("[T:") {
        if let Some(end) = line[pos..].find(']') {
            return format!("{}{}", &line[..pos], &line[pos + end + 1..]);
        }
    }
    line.to_string()
}

#[test]
fn checkpoint_a_matches_java_debruijn_dot() {
    let graph_text = fs::read_to_string(repo_root().join("fixtures/p3/quantify/c0/orig.graph.out"))
        .expect("需要 P3 的 c0 graph.out fixture");
    let orig_dot = fs::read_to_string(repo_root().join("fixtures/p4/c0/checkpoint_A.orig.dot"))
        .expect("需要 Java 生成的检查点 A 参照 DOT");

    let br = build_new_graph_use_kmers(&graph_text).unwrap();
    let our_dot = write_dot_string(&br.graph, br.ctx.kmer_size, false);

    let (our_nodes, our_edges) = normalize_dot(&our_dot);
    let (orig_nodes, orig_edges) = normalize_dot(&orig_dot);

    assert_eq!(our_nodes.len(), orig_nodes.len(), "节点数不一致");
    assert_eq!(our_edges.len(), orig_edges.len(), "边数不一致");
    let node_diff: Vec<_> = our_nodes.symmetric_difference(&orig_nodes).collect();
    let edge_diff: Vec<_> = our_edges.symmetric_difference(&orig_edges).collect();
    assert!(
        node_diff.is_empty() && edge_diff.is_empty(),
        "节点差异: {node_diff:#?}\n边差异: {edge_diff:#?}"
    );
}

#[test]
fn checkpoint_a_graph_invariants() {
    let graph_text = fs::read_to_string(repo_root().join("fixtures/p3/quantify/c0/orig.graph.out"))
        .expect("需要 P3 的 c0 graph.out fixture");
    let br = build_new_graph_use_kmers(&graph_text).unwrap();

    // KMER_SIZE 推断 = 24（err.log 中 Java 输出 "KMER_SIZE=24"）
    assert_eq!(br.ctx.kmer_size, 24);
    // 每个 root 节点：无前驱、每碱基权重 = supp
    assert!(!br.root_ids.is_empty());
    for &rid in &br.root_ids {
        assert_eq!(br.graph.in_degree(rid), 0, "root {rid} 应无入边");
        let v = br.graph.get_vertex(rid).unwrap();
        assert_eq!(v.weights.len(), v.name.len());
        assert!(v.weights.iter().all(|&w| w == v.weights[0]));
    }
    // 非 root 节点：无权重（weightAvg = -1）
    for &id in br.graph.vertex_ids() {
        if !br.root_ids.contains(&id) && br.graph.in_degree(id) > 0 {
            let v = br.graph.get_vertex(id).unwrap();
            if !br.root_ids.contains(&id) {
                assert!(v.weights.is_empty(), "非 root 节点 {id} 不应有权重");
            }
        }
    }
    // 边数 + root 数 <= 数据行数（root 行不产生边；toV 已存在的 root 行什么都不做）
    let n_data_lines = graph_text.lines().count() - 1;
    assert!(br.graph.edge_count() + br.root_ids.len() <= n_data_lines);
    // c0 有 4 个 from=-1 的 root 行，但 5645 更早作为普通节点入图（L1267 作为 from
    // 被创建），root 行到来时 toV 已存在 → 不重复建、无每碱基权重（Java 同语义）
    assert_eq!(br.root_ids.len(), 3);
    assert!(!br.root_ids.contains(&5645));
    // Java 检查点 A：2881 节点行 / 2880 边行
    assert_eq!(br.graph.vertex_count(), 2881);
    assert_eq!(br.graph.edge_count(), 2880);
}

#[test]
fn checkpoint_a_reads_parse() {
    let reads_text = fs::read_to_string(repo_root().join("fixtures/p3/quantify/c0/orig.reads.out"))
        .expect("需要 P3 的 c0 reads.out fixture");
    let reads = parse_graph_reads(&reads_text, 24, false).unwrap();
    assert!(!reads.is_empty());
    // 名字后缀已剥离
    assert!(reads
        .iter()
        .all(|r| !r.name.ends_with("/1") && !r.name.ends_with("/2")));
    // off-by-one：end = f3 + 24 > f3
    assert!(reads.iter().all(|r| r.end_in_read > r.start_in_read));
}
