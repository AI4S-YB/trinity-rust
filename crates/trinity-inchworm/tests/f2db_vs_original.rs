//! P3-T5 对拍——`fixtures/p3/bundle.orig.fa`（真实 Chrysalis 捆绑输出）经
//! 原版 `Inchworm/bin/FastaToDeBruijn --fasta bundle.orig.fa -K 24
//! --graph_per_record`（默认 DS）产物 `f2db.orig.txt` 比对。
//!
//! 比较契约：**每 Component 块的行多重集**（块序不管——原版 omp 并行 +
//! critical 输出块序不定；行内 5 列全比，包括节点 id 与 '-' 定向 id+N 偏移
//! ——id 由首次插入序决定，两侧 add 顺序相同则应一致）。

use std::collections::BTreeMap;
use std::io::BufReader;
use std::path::Path;

use trinity_common::fasta::FastaReader;
use trinity_inchworm::debruijn::graph_per_record;

fn fixture(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/p3")
        .join(name)
}

fn read_bundle(path: &std::path::PathBuf) -> Vec<trinity_common::fasta::FastaRecord> {
    let f = std::fs::File::open(path).unwrap();
    let mut reader = FastaReader::new(BufReader::new(f));
    let mut recs = Vec::new();
    while let Some(r) = reader.next_record().unwrap() {
        recs.push(r);
    }
    recs
}

/// "Component <id>" 起始的块 → (component_id, 块内行多重集)
fn blocks(text: &str) -> BTreeMap<i64, Vec<String>> {
    let mut out: BTreeMap<i64, Vec<String>> = BTreeMap::new();
    let mut cur: Option<i64> = None;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("Component ") {
            cur = Some(rest.parse().unwrap());
        } else if let Some(id) = cur {
            out.entry(id).or_default().push(line.to_string());
        }
    }
    out
}

#[test]
fn graph_per_record_matches_original_on_real_bundle() {
    let bundles = read_bundle(&fixture("bundle.orig.fa"));
    assert_eq!(bundles.len(), 55);

    let ours = graph_per_record(&bundles, 24, false).unwrap();
    let orig = std::fs::read_to_string(fixture("f2db.orig.txt")).unwrap();

    let ours_blocks = blocks(&ours);
    let orig_blocks = blocks(&orig);

    assert_eq!(
        ours_blocks.len(),
        orig_blocks.len(),
        "component 块数不一致: ours={} orig={}",
        ours_blocks.len(),
        orig_blocks.len()
    );

    let mut block_diffs = 0;
    let mut total_lines = 0;
    for (id, orig_lines) in &orig_blocks {
        total_lines += orig_lines.len();
        let our_lines = ours_blocks.get(id).expect("缺 component 块");
        let mut a = our_lines.clone();
        let mut b = orig_lines.clone();
        a.sort();
        b.sort();
        if a.len() != b.len() {
            block_diffs += 1;
            continue;
        }
        if a != b {
            // id 应一致（首次插入序相同）；若只有行序差异已在 sort 消除
            let only_ours: Vec<_> = a.iter().filter(|l| !b.contains(l)).take(3).collect();
            let only_orig: Vec<_> = b.iter().filter(|l| !a.contains(l)).take(3).collect();
            panic!("component {id} 行多重集不一致\n仅我方: {only_ours:?}\n仅原版: {only_orig:?}");
        }
    }
    assert_eq!(block_diffs, 0, "存在行数不一致的块");
    assert!(total_lines > 1000, "总行数异常: {total_lines}");
}
