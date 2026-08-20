//! Golden comparison against Butterfly's jaligner jar (NeedlemanWunschGotoh and
//! NeedlemanWunschGotohBanded, f32): aligned strings, score (bit-level f32) and
//! traceback start positions must match Java exactly.

use std::fs;
use std::path::PathBuf;

use trinity_butterfly::align::{nw_gotoh, nw_gotoh_banded};

fn repo_root() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.pop();
    p.pop();
    p
}

#[test]
fn golden_nw_and_banded() {
    let root = repo_root();
    let input = fs::read_to_string(root.join("fixtures/p4/align/align_golden_input.tsv")).unwrap();
    let golden = fs::read_to_string(root.join("fixtures/p4/align/align_golden.tsv")).unwrap();

    let inputs: Vec<Vec<&str>> = input
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').collect())
        .collect();
    let golds: Vec<Vec<&str>> = golden
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(inputs.len(), golds.len(), "input/golden line count");

    let mut checked = 0;
    for (inp, gold) in inputs.iter().zip(golds.iter()) {
        let mode = inp[0];
        let s1 = inp[1].as_bytes();
        let s2 = inp[2].as_bytes();
        let aln = if mode == "B" {
            let bw: usize = inp[3].parse().unwrap();
            nw_gotoh_banded(s1, s2, 4.0, -5.0, 10.0, 1.0, bw)
        } else {
            nw_gotoh(s1, s2, 4.0, -5.0, 10.0, 1.0)
        };

        let ctx = format!(
            "mode={} s1={} s2={} bw={}",
            mode,
            inp[1],
            inp[2],
            inp.get(3).copied().unwrap_or("-")
        );
        assert_eq!(
            String::from_utf8_lossy(&aln.aligned1),
            gold[0],
            "aligned1: {ctx}"
        );
        assert_eq!(
            String::from_utf8_lossy(&aln.aligned2),
            gold[1],
            "aligned2: {ctx}"
        );
        let gscore: f32 = gold[2].parse().unwrap();
        assert_eq!(
            aln.score.to_bits(),
            gscore.to_bits(),
            "score bits: {ctx} (rust {} vs java {})",
            aln.score,
            gscore
        );
        let g_start1: usize = gold[3].parse().unwrap();
        let g_start2: usize = gold[4].parse().unwrap();
        assert_eq!(aln.start1, g_start1, "start1: {ctx}");
        assert_eq!(aln.start2, g_start2, "start2: {ctx}");
        checked += 1;
    }
    assert!(checked >= 30, "expected a decent case count, got {checked}");
}
