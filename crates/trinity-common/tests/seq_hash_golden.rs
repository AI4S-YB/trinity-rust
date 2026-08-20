//! generateHash 黄金对拍: fixtures/p2/hash_golden.tsv 由 xtask gen-fixtures
//! 经 C harness（链原版 sequenceUtil.cpp）生成，逐行断言 u64 相等。

use trinity_common::seq_hash::generate_hash;

#[test]
fn hash_golden_tsv() {
    let path = format!(
        "{}/../../fixtures/p2/hash_golden.tsv",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&path).unwrap();
    let mut n = 0;
    for line in text.lines() {
        let (seq, val) = line
            .split_once('\t')
            .unwrap_or_else(|| panic!("行无 TAB: {line:?}"));
        let expected: u64 = val.parse().unwrap();
        assert_eq!(
            generate_hash(seq.as_bytes()),
            expected,
            "seq={seq:?} 不匹配"
        );
        n += 1;
    }
    assert_eq!(n, 20, "黄金行数应为 20（含 1 行空串）");
}
