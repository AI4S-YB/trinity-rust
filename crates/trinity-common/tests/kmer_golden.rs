//! 黄金向量: 与原版 C++（sequenceUtil.cpp 直链 harness）位级一致。

use trinity_common::kmer::{
    compute_entropy, decode_kmer_from_intval, get_ds_kmer_val, kmer_to_intval, revcomp_val,
};

#[test]
fn kmer_ops_match_original_cpp_bit_for_bit() {
    let tsv = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../fixtures/kmer_golden.tsv"
    ))
    .expect("缺少 fixtures/kmer_golden.tsv — 先跑 cargo xtask gen-fixtures");
    let mut n = 0;
    for line in tsv.lines() {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 5, "行格式异常: {line}");
        let kmer = cols[0].as_bytes();
        let intval: u64 = cols[1].parse().unwrap();
        let revcomp: u64 = cols[2].parse().unwrap();
        let dsval: u64 = cols[3].parse().unwrap();
        // %.9g 十进制位足以唯一确定 f32 → 解析回的 f32 与原版逐位一致
        let entropy: f32 = cols[4].parse().unwrap();
        let k = kmer.len();

        assert_eq!(
            kmer_to_intval(kmer).unwrap(),
            intval,
            "intval 不一致: {line}"
        );
        assert_eq!(revcomp_val(intval, k), revcomp, "revcomp 不一致: {line}");
        assert_eq!(get_ds_kmer_val(intval, k), dsval, "dsval 不一致: {line}");
        assert_eq!(
            compute_entropy(intval, k),
            entropy,
            "entropy 不一致: {line}"
        );
        // decode 与输入大写形式互逆（小写输入编码相同、解码为大写）
        assert_eq!(
            decode_kmer_from_intval(intval, k),
            kmer.to_ascii_uppercase()
        );

        n += 1;
    }
    assert!(n >= 26, "黄金向量行数异常: {n}");
}
