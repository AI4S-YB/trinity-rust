//! 对拍 fixture: 原版 Inchworm/bin/fastaToKmerCoverageStats 输出（--num_threads 1，
//! tid 恒为 thread:0）vs 本移植。5 列全等（acc/median/mean/stdev/tid）。
//!
//! fixture 生成（k=25，DS）:
//!   jellyfish count -m 25 --canonical <reads> && jellyfish dump -L 1 → {name}.kmers.fa
//!   fastaToKmerCoverageStats --reads <reads> --kmers <kmers> --kmer_size 25 --num_threads 1
//!     → {name}.stats.orig.tsv

use std::fs;

use trinity_kmer::coverage_stats::{coverage_stats_rows, load_kmer_dump, write_stats_tsv};

fn fixture(name: &str) -> String {
    format!("{}/../../fixtures/p1/{name}", env!("CARGO_MANIFEST_DIR"))
}

fn assert_matches_original(name: &str) {
    let reads = fs::read(fixture(&format!("{name}.fa"))).unwrap();
    let kmers = fs::read(fixture(&format!("{name}.kmers.fa"))).unwrap();
    let orig = fs::read_to_string(fixture(&format!("{name}.stats.orig.tsv"))).unwrap();

    let counts = load_kmer_dump(&kmers, 25, true).unwrap();
    let rows = coverage_stats_rows(&reads, &counts, 25, true).unwrap();
    let mut ours = Vec::new();
    write_stats_tsv(&mut ours, &rows).unwrap();
    let ours = String::from_utf8(ours).unwrap();

    let orig_lines: Vec<&str> = orig.lines().collect();
    let our_lines: Vec<&str> = ours.lines().collect();
    assert_eq!(orig_lines.len(), our_lines.len(), "{name}: 行数不一致");
    // 表头逐字（cpp L116-120）
    assert_eq!(orig_lines[0], "acc\tmedian_cov\tmean_cov\tstdev\ttid");
    assert_eq!(our_lines[0], orig_lines[0], "{name}: 表头不一致");

    for (i, (o, u)) in orig_lines.iter().zip(&our_lines).enumerate().skip(1) {
        // 5 列全等; 失败信息带双方完整行便于定位
        assert_eq!(
            o, u,
            "{name}: 第 {i} 行不一致\n  orig: {o:?}\n  ours: {u:?}"
        );
    }
}

#[test]
fn smoke_stats_match_original() {
    assert_matches_original("smoke");
}

#[test]
fn edge_stats_match_original() {
    // 覆盖: <25bp 短 read（0/0/-0）、含 N read、全同 read、小写、折行、恰 25bp（-nan）
    assert_matches_original("edge");
}

/// 短 read / 单窗口 read 的原版特殊形态逐字锁定（防止 fixture 被重生成后悄悄漂移）。
#[test]
fn edge_special_rows_locked() {
    let orig = fs::read_to_string(fixture("edge.stats.orig.tsv")).unwrap();
    let lines: Vec<&str> = orig.lines().collect();
    let by_acc = |acc: &str| *lines.iter().find(|l| l.starts_with(acc)).unwrap();
    assert_eq!(
        by_acc("edge_short_20bp"),
        "edge_short_20bp\t0\t0\t-0\tthread:0"
    );
    assert_eq!(by_acc("edge_exact25"), "edge_exact25\t1\t1\t-nan\tthread:0");
    assert_eq!(by_acc("edge_dup_a"), "edge_dup_a\t3\t3.25\t0.5\tthread:0");
    // 含 N read: 表头带描述但 acc 只取首 token
    assert_eq!(
        by_acc("edge_with_N"),
        "edge_with_N\t7\t5\t3.31662\tthread:0"
    );
}
