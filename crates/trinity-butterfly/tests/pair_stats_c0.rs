//! c0 组件配对组装对拍：复用 T5 的 threading_golden.tsv（jar `-V 17 --stderr`
//! "Threaded Read as" 黄金，顺序 = readNameHash occurrence 序），跑
//! getSuffStats_wPairs，与已知统计（4336 reads / 2423 pairs）自洽对拍。
//!
//! Java 侧 getSuffStats_wPairs 的 "## Read PathPair results: ..." 统计行需 -V 10，
//! 当时 err.txt 未留存；按任务规格退化为总数自洽断言（计数与分桶
//! start-vertex 语义由单元测试锁定）。

use std::collections::HashSet;
use std::fs;
use std::path::Path;

use trinity_butterfly::{
    build_new_graph_use_kmers, calc_sub_components_stats, compact_linear_paths,
    fix_extremely_high_single_edges, get_suff_stats_w_pairs, pre_process_graph_file,
    remove_light_edges, remove_single_nt_bubbles, run_dfs2, ReadPath,
};

/// c0 剪枝链产物（combinePaths 的可达性判定需要）。
fn c0_pruned_graph() -> trinity_butterfly::DiGraph {
    let graph_text =
        fs::read_to_string(repo_root().join("fixtures/p3/quantify/c0/orig.graph.out")).unwrap();
    let br = build_new_graph_use_kmers(&graph_text).unwrap();
    let (mut graph, mut ctx) = (br.graph, br.ctx);
    let (in_flow, out_flow, _) = pre_process_graph_file(&graph_text);
    fix_extremely_high_single_edges(&mut graph, &out_flow, &in_flow);
    remove_light_edges(&mut graph, &ctx);
    compact_linear_paths(&mut graph, &ctx);
    run_dfs2(&mut graph);
    remove_single_nt_bubbles(&mut graph, &mut ctx);
    calc_sub_components_stats(&mut graph, 200);
    graph
}

fn repo_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
}

/// 解析 threading_golden.tsv（name \t 逗号分隔节点 id）为 (name, ReadPath) 序。
fn load_golden() -> Vec<(String, ReadPath)> {
    let text = fs::read_to_string(repo_root().join("fixtures/p4/c0/threading_golden.tsv")).unwrap();
    text.lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let (name, ids) = l.split_once('\t').expect("name\\tpath");
            let path: Vec<i32> = if ids.is_empty() {
                Vec::new()
            } else {
                ids.split(',').map(|x| x.trim().parse().unwrap()).collect()
            };
            assert!(!path.is_empty(), "成功穿线的 read 路径非空");
            (
                name.to_string(),
                ReadPath {
                    mismatch_count: 0,
                    path,
                    positions: Vec::new(),
                },
            )
        })
        .collect()
}

#[test]
fn c0_pair_stats_self_consistent() {
    let ordered = load_golden();
    assert_eq!(
        ordered.len(),
        4336,
        "穿线成功 read 数（threading_stats.txt）"
    );

    let names: HashSet<&str> = ordered.iter().map(|(n, _)| n.as_str()).collect();
    assert_eq!(names.len(), 2423, "不同 read 名数 = jar pairs 统计");

    // 从黄金直接推导 singleton/pair 划分（Java 只取同名前两条）
    let mut occ: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
    for (n, _) in &ordered {
        *occ.entry(n.as_str()).or_default() += 1;
    }
    let exp_singletons = occ.values().filter(|&&c| c == 1).count() as i64;
    let exp_pairs = occ.values().filter(|&&c| c > 1).count() as i64;

    let s = get_suff_stats_w_pairs(&c0_pruned_graph(), &ordered);
    assert_eq!(s.num_singletons, exp_singletons);
    assert_eq!(s.num_pairs, exp_pairs);
    assert_eq!(s.num_reads_used, 2423, "每名恰好一条 PairPath");
    assert_eq!(
        s.num_pairs_discarded, 0,
        "c0 配对全部可合并（jar 实际执行 combinePaths，见 pair_paths.rs 文档）"
    );
    assert_eq!(s.total_count(), 2423, "combined_read_hash 总计数自洽");
    // 每个桶的键 firstV 与桶 id 一致
    for (first_v, m) in &s.combined_read_hash {
        for pp in m.keys() {
            assert_eq!(pp.path1.first().copied(), Some(*first_v));
        }
    }
    // c0 数据无 LR$| 长读
    assert!(s.long_read_name_to_ppath.is_empty());
    assert!(s.long_read_path_map.is_empty());
    println!(
        "c0 pair stats: {} start vertices, {} singletons, {} pairs",
        s.num_start_vertices(),
        s.num_singletons,
        s.num_pairs
    );
}
