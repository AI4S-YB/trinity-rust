//! M2.1 PogFlow 引擎 c0 fixture：POG 图 + 流分解产出非空 FASTA。

use std::fs;
use std::path::Path;

use trinity_butterfly::{run_component, ComponentParams, PathEngine};

fn repo_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
}

#[test]
fn pog_flow_produces_fasta() {
    let graph_text =
        fs::read_to_string(repo_root().join("fixtures/p3/quantify/c0/orig.graph.out")).unwrap();
    let reads_text =
        fs::read_to_string(repo_root().join("fixtures/p3/quantify/c0/orig.reads.out")).unwrap();
    let out = run_component(
        &graph_text,
        &reads_text,
        &ComponentParams {
            path_engine: PathEngine::PogFlow,
            min_len: 200,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(
        out.all_prob_paths_fasta.contains('>'),
        "PogFlow 输出无 FASTA header"
    );
    assert!(
        out.all_prob_paths_fasta.len() > 200,
        "PogFlow 输出过短: {}",
        out.all_prob_paths_fasta.len()
    );
    assert!(out.num_paths >= 1, "PogFlow 路径数为 0");
}
