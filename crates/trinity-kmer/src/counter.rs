//! 并行 k-mer 计数 — jellyfish count 的 Rust 替代（多重集等价，见 xcheck-kmer）
//! 语义: 小写 gatc 视同大写; 窗口遇非 gatc 跳过该 kmer 并继续滑 1 碱基;
//! DS 模式表键 = max(kmer, revcomp(kmer))（sequenceUtil.cpp:376）。
//! 记录切分对齐 jellyfish: 序列段内滤除 `\n`/`\r` 跨行拼接（T2 实测确认）;
//! 行内空白等其他非 gatc 字符仍按非法字符断窗（Trinity 管道产物无行内空白）。

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use trinity_common::kmer::{get_ds_kmer_val, KmerId};

pub type CountMap = FxHashMap<KmerId, u32>;

/// 对一条序列（原始字节）枚举有效 kmer，回调规范化后的键。
/// 镜像 jellyfish: 非 gatc 字符使包含它的窗口失效，窗口继续前滑（不是重启）。
pub fn for_each_kmer(seq: &[u8], k: usize, ds: bool, mut f: impl FnMut(KmerId)) {
    let mask: KmerId = if k >= 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k as u32)) - 1
    };
    let mut kmer: KmerId = 0;
    let mut valid: usize = 0;
    for &c in seq {
        let v = match c {
            b'G' | b'g' => 0u64,
            b'A' | b'a' => 1,
            b'T' | b't' => 2,
            b'C' | b'c' => 3,
            _ => {
                valid = 0;
                kmer = 0;
                continue;
            }
        };
        kmer = ((kmer << 2) | v) & mask;
        if valid < k {
            valid += 1;
        }
        if valid == k {
            let key = if ds { get_ds_kmer_val(kmer, k) } else { kmer };
            f(key);
        }
    }
}

/// 按记录切分 FASTA 字节流: 每条记录 = `>` header 行之后的序列字节段
/// （到下一个行首 `>` 或文件尾），段内滤除 `\n`/`\r`（jellyfish 跨行拼接，
/// 实测对 ">r1\nACGT\nACGTACGT\n" 按连续 12 字符滑窗）。
/// 其余非 gatc 字符不过滤，由 for_each_kmer 按非法字符断窗。
/// 返回 owned 记录（滤行需脱离原缓冲区）。
pub fn split_fasta_sequences(data: &[u8]) -> Vec<Vec<u8>> {
    let mut records: Vec<Vec<u8>> = Vec::new();
    let mut pos = 0;
    while pos < data.len() && data[pos] == b'>' {
        // header 行结束于 \n; 未闭合（无 \n）则该记录无序列
        let Some(hdr_nl) = data[pos..].iter().position(|&b| b == b'\n') else {
            break;
        };
        let seq_start = pos + hdr_nl + 1;
        // 序列段终点 = 下一个行首 '>'（即 (\n, '>') 字节对）或文件尾
        let end = data
            .windows(2)
            .skip(seq_start - 1)
            .position(|w| w[0] == b'\n' && w[1] == b'>')
            .map_or(data.len(), |p| seq_start + p);
        // 段内滤除 \n 与 \r（对齐 jellyfish 跨行拼接），其余字节保留
        let rec: Vec<u8> = data[seq_start..end]
            .iter()
            .copied()
            .filter(|&b| b != b'\n' && b != b'\r')
            .collect();
        records.push(rec);
        pos = end;
    }
    records
}

pub struct KmerCountTable;

impl KmerCountTable {
    /// 并行统计 FASTA 字节流的 k-mer 计数。
    /// rayon 按记录分片 → 每线程局部 FxHashMap → reduce 合并（jellyfish 同型策略）。
    pub fn count_fasta_data(data: &[u8], k: usize, ds: bool) -> CountMap {
        Self::count_fasta_data_with_capacity(data, k, ds, 0)
    }

    /// 带每线程预 reserve 的计数（`estimate_hash_size` 接线）:
    /// `per_thread_capacity` 为每线程局部 HashMap 的初始容量（0 = 不预分配,
    /// 即旧行为）。调用方需保证 capacity 总量（× 线程数 × 每项 ~20B）在
    /// max_memory 预算内——orchestrate 侧按 estimate/18/threads 截断。
    pub fn count_fasta_data_with_capacity(
        data: &[u8],
        k: usize,
        ds: bool,
        per_thread_capacity: usize,
    ) -> CountMap {
        let cap = per_thread_capacity.min(1 << 22); // 单线程上限 ~4M 槽，防 rehash 峰值即可
        split_fasta_sequences(data)
            .par_iter()
            .fold(
                || CountMap::with_capacity_and_hasher(cap, Default::default()),
                |mut local, rec| {
                    for_each_kmer(rec, k, ds, |key| {
                        *local.entry(key).or_insert(0) += 1;
                    });
                    local
                },
            )
            .reduce(CountMap::default, merge_counts)
    }
}

/// reduce 合并: b 并入 a（计数累加）。
fn merge_counts(mut a: CountMap, b: CountMap) -> CountMap {
    for (key, cnt) in b {
        *a.entry(key).or_insert(0) += cnt;
    }
    a
}

/// Trinity:2598-2604 hash 大小护栏（信息性容量提示用，语义与原版一致）
pub fn estimate_hash_size(max_memory: u64, read_file_size: u64) -> u64 {
    let mut h = if max_memory > read_file_size {
        (max_memory - read_file_size) / 7
    } else {
        0
    };
    if h < 100_000_000 || read_file_size < 5_000_000_000 {
        h = 100_000_000;
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;
    use trinity_common::kmer::{decode_kmer_from_intval, kmer_to_intval};

    fn counts_of(seq: &str, k: usize, ds: bool) -> Vec<(String, u32)> {
        let mut m: CountMap = FxHashMap::default();
        for_each_kmer(seq.as_bytes(), k, ds, |key| *m.entry(key).or_insert(0) += 1);
        let mut v: Vec<_> = m
            .into_iter()
            .map(|(key, c)| {
                (
                    String::from_utf8(decode_kmer_from_intval(key, k)).unwrap(),
                    c,
                )
            })
            .collect();
        v.sort();
        v
    }

    #[test]
    fn basic_ss_counting() {
        assert_eq!(
            counts_of("ACGT", 2, false),
            vec![("AC".into(), 1), ("CG".into(), 1), ("GT".into(), 1)]
        );
    }

    #[test]
    fn ds_collapses_revcomp() {
        assert_eq!(counts_of("AA TT", 2, true), vec![("TT".into(), 2)]);
        assert_eq!(counts_of("GT", 2, true), vec![("AC".into(), 1)]);
    }

    #[test]
    fn lowercase_and_n_window_slide() {
        assert_eq!(counts_of("acgt", 2, false)[0].1, 1);
        assert_eq!(counts_of("ACNT", 2, false), vec![("AC".into(), 1)]);
        // 手推 ACGNNT: 无 N 窗口 = AC(0,1), CG(1,2) → 计两个（计划稿漏了 CG，已修正;
        // 只有包含非法字符的窗口失效，G 在 N 前的 CG 窗口有效）
        assert_eq!(
            counts_of("ACGNNT", 2, false),
            vec![("AC".into(), 1), ("CG".into(), 1)]
        );
        // 滑动 vs 重启判别: N 后窗口恢复即计 CA（重启实现会漏）
        // 手推 AANCAA: 无 N 的 2-mer 窗口 = AA(0,1), CA(3,4), AA(4,5) → AA×2, CA×1
        assert_eq!(
            counts_of("AANCAA", 2, false),
            vec![("AA".into(), 2), ("CA".into(), 1)]
        );
    }

    #[test]
    fn short_seq_no_kmer() {
        assert_eq!(counts_of("AC", 2, false).len(), 1);
        assert!(counts_of("A", 2, false).is_empty());
    }

    #[test]
    fn estimate_hash_size_mirrors_trinity() {
        // Trinity:2601 护栏是 `||`: read_file_size < 5e9 即触发 100e6 下限
        // 手推: 1e9 < 5e9 → 护栏触发（计划稿的 9e9/7 期望未过护栏，已修正）
        assert_eq!(
            estimate_hash_size(10_000_000_000, 1_000_000_000),
            100_000_000
        );
        assert_eq!(
            estimate_hash_size(10_000_000_000, 4_000_000_000),
            100_000_000
        );
        // 大输入文件（≥5e9）且估计值 ≥100e6 时 /7 估计保留
        assert_eq!(
            estimate_hash_size(10_000_000_000, 5_000_000_000),
            5_000_000_000 / 7
        );
        // 下溢防护分支: max_memory ≤ read_file_size → 0 → 护栏
        assert_eq!(estimate_hash_size(1_000, 2_000), 100_000_000);
    }

    #[test]
    fn count_fasta_data_end_to_end() {
        let fa = b">r1\nACGTACGT\n>r2\nACGTACGT\n>r3\nTTTT\n";
        let m = KmerCountTable::count_fasta_data(fa, 4, false);
        let key = |s: &str| kmer_to_intval(s.as_bytes()).unwrap();
        assert_eq!(m[&key("ACGT")], 4);
        assert_eq!(m[&key("CGTA")], 2);
        assert_eq!(m[&key("TTTT")], 1);
        let mds = KmerCountTable::count_fasta_data(fa, 4, true);
        assert_eq!(mds[&key("TTTT")], 1);
    }

    #[test]
    fn count_fasta_multiline_and_crlf() {
        // jellyfish 跨行拼接: 序列段内 \n/\r 滤除后滑窗（实测 jellyfish 2.3.x 对
        // ">r1\nACGT\nACGTACGT\n" 输出 ACGT=3 CGTA=2 GTAC=2 TACG=2，即 12 字符连续滑窗）。
        // 其余非 gatc 字符（N、行内空格等）仍断窗（T1 滑动语义测试锁定）。
        let fa = b">r1 desc\nACGT\nACGT\r\nTT\n";
        let m = KmerCountTable::count_fasta_data(fa, 4, false);
        let key = |s: &str| kmer_to_intval(s.as_bytes()).unwrap();
        // 滤 \n\r 后序列 = ACGTACGTTT（10 字符）→ 7 个窗口:
        // ACGT CGTA GTAC TACG ACGT CGTT GTTT
        assert_eq!(m[&key("ACGT")], 2);
        assert_eq!(m[&key("CGTA")], 1);
        assert_eq!(m[&key("GTAC")], 1);
        assert_eq!(m[&key("TACG")], 1);
        assert_eq!(m[&key("CGTT")], 1);
        assert_eq!(m[&key("GTTT")], 1);
        assert_eq!(m.len(), 6);
    }

    #[test]
    fn count_fasta_joins_lines_like_jellyfish() {
        // 与 jellyfish 对拍原例（T2 冒烟固化）: 跨行窗口 CGTA/GTAC/TACG 只在拼接语义下出现
        let fa = b">r1\nACGT\nACGTACGT\n";
        let m = KmerCountTable::count_fasta_data(fa, 4, false);
        let key = |s: &str| kmer_to_intval(s.as_bytes()).unwrap();
        assert_eq!(m[&key("ACGT")], 3);
        assert_eq!(m[&key("CGTA")], 2);
        assert_eq!(m[&key("GTAC")], 2);
        assert_eq!(m[&key("TACG")], 2);
        assert_eq!(m.len(), 4);
    }
}

#[cfg(test)]
mod capacity_tests {
    use super::*;

    /// 预 reserve 版本与默认版本多重集等价。
    #[test]
    fn with_capacity_matches_default() {
        let fa = b">r1\nACGTACGTACGTACGTACGTACGTACGTACGT\n>r2\nTTTTACGTACGTACGTACGTACGTACGTACG\n";
        let a = KmerCountTable::count_fasta_data(fa, 25, false);
        let b = KmerCountTable::count_fasta_data_with_capacity(fa, 25, false, 1 << 20);
        assert_eq!(a.len(), b.len());
        assert!(a.iter().all(|(k, v)| b.get(k) == Some(v)));
    }
}
