//! P3-T2 对拍——原版 GraphFromFasta 二进制（trinityrnaseq-v2.15.2）锁定黄金。
//!
//! fixture 生成命令与结论见 `fixtures/p3/README.md`。比较契约：**边多重集**
//! （有向 (A,B) 对 + weldmers/total/min_len 全字段）相等——行序非契约
//! （原版 report 非稳定 sort + OMP 并行插入序 vs 本库确定性序，
//! 见 `graph_from_fasta` 模块文档）。

use std::collections::BTreeMap;
use std::path::Path;

use trinity_chrysalis::dna_vector::{read_fasta, stream_fasta_records};
use trinity_chrysalis::graph_from_fasta::{graph_from_fasta, sort_weld_graph, GffParams};

fn fixture(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/p3")
        .join(name)
}

/// 行 → (A, B, weldmers, total, min_len) 全字段元组（缺字段行跳过）。
fn parse_edges(text: &str) -> BTreeMap<(usize, usize, u32, u32, i64), usize> {
    let mut m = BTreeMap::new();
    for l in text.lines() {
        let f: Vec<&str> = l.split_whitespace().collect();
        if f.len() < 11 {
            continue;
        }
        let a: usize = f[0].parse().unwrap();
        let b: usize = f[2].parse().unwrap();
        let w: u32 = f[4].parse().unwrap();
        let t: u32 = f[8].parse().unwrap();
        let mlen: i64 = f[10].parse().unwrap();
        *m.entry((a, b, w, t, mlen)).or_insert(0) += 1;
    }
    m
}

#[test]
fn edge_multiset_matches_original() {
    let iworm = read_fasta(fixture("gff.iworm.fa")).unwrap();
    assert_eq!(iworm.len(), 219);
    let reads = stream_fasta_records(&std::fs::read_to_string(fixture("gff.reads.fa")).unwrap());
    assert_eq!(reads.len(), 61150);

    let ours = graph_from_fasta(&iworm, &reads, &GffParams::default()).unwrap();
    let orig = std::fs::read_to_string(fixture("gff.welds.orig.txt")).unwrap();

    let (e_ours, e_orig) = (parse_edges(&ours), parse_edges(&orig));
    assert!(
        !e_orig.is_empty() && e_orig.len() >= 20,
        "黄金边数异常: {}",
        e_orig.len()
    );
    assert_eq!(
        e_ours, e_orig,
        "边多重集（含 weldmers/total/min_len）应相等"
    );
}

/// sort_weld_graph 与 GNU `sort -k9,9gr` 在黄金 weld 图上的逐行一致
/// （total 值两两不同时无 tie，稳定/非稳定无差别）。
#[test]
fn sort_weld_graph_matches_gnu_sort_on_fixture() {
    let orig = std::fs::read_to_string(fixture("gff.welds.orig.txt")).unwrap();
    let sorted = sort_weld_graph(&orig);
    let totals: Vec<f64> = sorted
        .lines()
        .map(|l| l.split_whitespace().nth(8).unwrap().parse().unwrap())
        .collect();
    let mut expect = totals.clone();
    expect.sort_by(|a, b| b.partial_cmp(a).unwrap());
    assert_eq!(totals, expect);
    assert_eq!(sorted.lines().count(), orig.lines().count());
}
