//! gtasm 边概率导出工具——数值对拍用（验收 #2）。
//!
//! 用法: `gtasm_score_dump <graph.out> <weights.gtw>`
//! 输出: 每边一行 `eid<TAB>src_vid<TAB>dst_vid<TAB>prob`（prob 为 %.8e f32 文本）。
//! python 侧 `scripts/gt_score_compare.py` 以 torch 重载同一 gtw 在同一图上
//! 重算, 逐边对比（容差 1e-4）。

use trinity_butterfly::{build_gt_features, parse_component_graph, GtModel};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() != 2 {
        eprintln!("用法: gtasm_score_dump <graph.out> <weights.gtw>");
        std::process::exit(2);
    }
    let graph_text = std::fs::read_to_string(&args[0])
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", args[0]));
    let model = GtModel::load(std::path::Path::new(&args[1]))
        .unwrap_or_else(|e| panic!("cannot load {}: {e}", args[1]));
    let g = parse_component_graph(&graph_text).unwrap_or_else(|e| panic!("parse: {e}"));
    let (x, e) = build_gt_features(&g);
    let prob = model.forward(&x, &e, &g.src, &g.dst);
    for (ei, p) in prob.iter().enumerate() {
        println!("{}\t{}\t{}\t{:.8e}", ei, g.vids[g.src[ei]], g.vids[g.dst[ei]], p);
    }
}
