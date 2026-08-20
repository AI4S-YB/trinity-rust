//! P3-T4 对拍——`fixtures/p3` 真实链（bundle.orig.fa + gff.reads.fa →
//! ReadsToTranscripts）与原版二进制产物（`rtt.orig.out` /
//! `rtt.orig.out.rcts.out`，命令 `-i reads -f bundle -o out -t 1 -p 50`）
//! 比对。
//!
//! 比较契约：两侧各自 `sort_reads_to_components` / `sort -k1,1n -k3,3nr
//! -k2,2` 后**逐行相等**（原版行序 = multimap 序，排序后消除序差）；
//! readCount 逐值相等。

use std::path::Path;

use trinity_chrysalis::dna_vector::{read_fasta, read_fasta_short_names};
use trinity_chrysalis::reads_to_transcripts::{reads_to_transcripts, RttParams};

fn fixture(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/p3")
        .join(name)
}

#[test]
fn rtt_matches_original_on_real_bundle() {
    let bundles = read_fasta(fixture("bundle.orig.fa")).unwrap();
    let reads = read_fasta_short_names(fixture("gff.reads.fa")).unwrap();
    assert!(!bundles.is_empty() && !reads.is_empty());

    let out = reads_to_transcripts(
        &reads,
        &bundles,
        &RttParams {
            strand: false,
            pct_required: 50,
            min_kmer_entropy: 1.5,
            max_mem_reads: usize::MAX,
            threads: 4,
        },
    )
    .unwrap();

    // readCount 契约
    let orig_count: u64 = std::fs::read_to_string(fixture("rtt.orig.out.rcts.out"))
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    assert_eq!(out.mapped_count, orig_count, "成功映射 read 总数");

    let ours = trinity_chrysalis::reads_to_transcripts::sort_reads_to_components(&out.text);
    let orig_raw = std::fs::read_to_string(fixture("rtt.orig.out")).unwrap();
    let orig = trinity_chrysalis::reads_to_transcripts::sort_reads_to_components(&orig_raw);

    let a: Vec<&str> = ours.lines().collect();
    let b: Vec<&str> = orig.lines().collect();
    assert_eq!(a.len(), b.len(), "行数");
    let mut diffs = 0;
    for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
        if x != y {
            if diffs < 5 {
                eprintln!("line {i}\n  ours: {x}\n  orig: {y}");
            }
            diffs += 1;
        }
    }
    assert_eq!(diffs, 0, "排序后逐行差异数");
}
