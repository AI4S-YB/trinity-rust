//! 并行 k-mer 计数 — jellyfish count 的 Rust 替代（多重集等价，见 xcheck-kmer）
//! 语义: 小写 gatc 视同大写; 窗口遇非 gatc 跳过该 kmer 并继续滑 1 碱基;
//! DS 模式表键 = max(kmer, revcomp(kmer))（sequenceUtil.cpp:376）。
//! 记录切分对齐 jellyfish: 序列段内滤除 `\n`/`\r` 跨行拼接（T2 实测确认）;
//! 行内空白等其他非 gatc 字符仍按非法字符断窗（Trinity 管道产物无行内空白）。

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use trinity_common::kmer::{KmerId, MAX_KMER_LENGTH, RollingKey};

pub type CountMap = FxHashMap<KmerId, u32>;

/// 对一条序列（原始字节）枚举有效 kmer，回调规范化后的键。
/// 镜像 jellyfish: 非 gatc 字符使包含它的窗口失效，窗口继续前滑（不是重启）。
///
/// 性能分层（2026-09-02）: 泛型核心 [`for_each_kmer_of`] 按键类型实例化——
/// u128 紧凑键（k≤64, diginorm 热路径, 24B/槽 vs U256 40B/槽）与 KmerId 全宽。
/// 本函数 = KmerId 实例化（管线对外签名不变）。
pub fn for_each_kmer(seq: &[u8], k: usize, ds: bool, f: impl FnMut(KmerId)) {
    for_each_kmer_of::<KmerId>(seq, k, ds, f)
}

/// 滚动枚举泛型核心; 语义与 jellyfish/原 for_each_kmer 全域一致
/// （DS canonical = max(kmer, revcomp), 非滚动兜底走 K::revcomp 等价式）。
pub fn for_each_kmer_of<K: RollingKey>(seq: &[u8], k: usize, ds: bool, mut f: impl FnMut(K)) {
    // 满宽兜底自限正确（超窗高位随左移自然掉出）; u128 实例的 low_mask 自带 k≤64 断言。
    let mask: K = K::low_mask(k);
    let mut kmer: K = K::zero();
    let mut valid: usize = 0;
    // M7-norm: 滚动 revcomp（revcomp(km<<2|v) == (revcomp(km)>>2) | (v^3)<<(2k-2)），
    // 免去每 kmer 的 k 迭代 revcomp_val 循环; 与 canonical max 语义数值恒等。
    // rev_shift = 2k-2 ≤ 254（k ≤ MAX_KMER_LENGTH=128）恒安全。
    if ds && k <= MAX_KMER_LENGTH {
        let rev_shift: u32 = (2 * k - 2) as u32;
        let mut rev: K = K::zero();
        for &c in seq {
            let v = match c {
                b'G' | b'g' => 0u8,
                b'A' | b'a' => 1,
                b'T' | b't' => 2,
                b'C' | b'c' => 3,
                _ => {
                    valid = 0;
                    kmer = K::zero();
                    rev = K::zero();
                    continue;
                }
            };
            kmer = kmer.roll_forward(v, mask);
            rev = rev.roll_rev(v, rev_shift, mask);
            if valid < k {
                valid += 1;
            }
            if valid == k {
                f(if rev > kmer { rev } else { kmer });
            }
        }
        return;
    }
    for &c in seq {
        let v = match c {
            b'G' | b'g' => 0u8,
            b'A' | b'a' => 1,
            b'T' | b't' => 2,
            b'C' | b'c' => 3,
            _ => {
                valid = 0;
                kmer = K::zero();
                continue;
            }
        };
        kmer = kmer.roll_forward(v, mask);
        if valid < k {
            valid += 1;
        }
        if valid == k {
            let key = if ds {
                // get_ds_kmer_val 的 K 泛型等价式（非滚动兜底路径）
                let rc = kmer.revcomp(k);
                if rc > kmer { rc } else { kmer }
            } else {
                kmer
            };
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

/// 借用视图版切分（性能）: 返回各记录的 (start, end) 区间 + "是否存在多行
/// 序列段"标记。单行序列段（段内 \n 仅在末尾）下, 对区间切片直接喂
/// for_each_kmer 与滤除 \n 后滑窗完全等价（段内无 \n 可断窗）;
/// 多行段存在时调用方须回退 [`split_fasta_sequences`]（跨行拼接语义）。
pub fn split_fasta_ranges(data: &[u8]) -> (Vec<(usize, usize)>, bool) {
    let mut ranges = Vec::new();
    let mut multiline = false;
    let mut pos = 0;
    while pos < data.len() && data[pos] == b'>' {
        let Some(hdr_nl) = data[pos..].iter().position(|&b| b == b'\n') else {
            break;
        };
        let seq_start = pos + hdr_nl + 1;
        let end = data
            .windows(2)
            .skip(seq_start - 1)
            .position(|w| w[0] == b'\n' && w[1] == b'>')
            .map_or(data.len(), |p| seq_start + p);
        // 段内换行仅允许出现在末尾（单行序列）; 段内 \r 一并滤除语义
        // 由调用方判断——此处只报多行
        let body = &data[seq_start..end];
        let trimmed = body.strip_suffix(b"\n").unwrap_or(body);
        if trimmed.contains(&b'\n') || trimmed.contains(&b'\r') {
            multiline = true;
        }
        ranges.push((seq_start, end));
        pos = end;
    }
    (ranges, multiline)
}

pub struct KmerCountTable;

impl KmerCountTable {
    /// 并行统计 FASTA 字节流的 k-mer 计数（KmerId 口径; 泛型核心见
    /// [`count_fasta_data_of`]——diginorm 热路径按 k 分层用 u64/u128
    /// 紧凑键实例化）。
    pub fn count_fasta_data(data: &[u8], k: usize, ds: bool) -> CountMap {
        Self::count_fasta_data_with_capacity(data, k, ds, 0)
    }

    /// 带每线程预 reserve 的计数（`estimate_hash_size` 接线）:
    /// `per_thread_capacity` 为每线程局部 HashMap 的初始容量（0 = 不预分配,
    /// 即旧行为）。调用方需保证 capacity 总量（× 线程数 × 每项开销）在
    /// max_memory 预算内——orchestrate 侧按 estimate/threads 截断。
    ///
    /// 性能（2026-09-01）: 单行序列走借用 range 视图（免逐记录拷贝,
    /// 切片无内部 \n 时 for_each_kmer 断窗行为与滤除版逐字节等价）;
    /// 检测到多行序列段自动回退 owned 拷贝路径（跨行拼接语义）。
    pub fn count_fasta_data_with_capacity(
        data: &[u8],
        k: usize,
        ds: bool,
        per_thread_capacity: usize,
    ) -> CountMap {
        count_fasta_data_of::<KmerId>(data, k, ds, per_thread_capacity)
    }
}

/// 并行计数泛型核心: rayon 按记录分片 → 每线程局部 FxHashMap<K,u32> →
/// reduce 合并（jellyfish 同型策略）。K = u64（k≤32）/u128（k≤64 紧凑键）
/// 或 KmerId。
pub fn count_fasta_data_of<K>(
    data: &[u8],
    k: usize,
    ds: bool,
    per_thread_capacity: usize,
) -> FxHashMap<K, u32>
where
    K: RollingKey + Send,
{
    let cap = per_thread_capacity.min(1 << 22); // 单线程上限 ~4M 槽，防 rehash 峰值即可
    let (ranges, multiline) = split_fasta_ranges(data);
    if multiline {
        return split_fasta_sequences(data)
            .par_iter()
            .fold(
                || FxHashMap::<K, u32>::with_capacity_and_hasher(cap, Default::default()),
                |mut local, rec| {
                    for_each_kmer_of::<K>(rec, k, ds, |key| {
                        *local.entry(key).or_insert(0) += 1;
                    });
                    local
                },
            )
            .reduce(FxHashMap::default, merge_counts);
    }
    ranges
        .par_iter()
        .fold(
            || FxHashMap::<K, u32>::with_capacity_and_hasher(cap, Default::default()),
            |mut local, &(s, e)| {
                for_each_kmer_of::<K>(&data[s..e], k, ds, |key| {
                    *local.entry(key).or_insert(0) += 1;
                });
                local
            },
        )
        .reduce(FxHashMap::default, merge_counts)
}

/// 窄键表 → CountMap（KmerId）一次性 widen（压时间战役 2026-09-04 刀1）:
/// 精确预留容量（无生长 rehash）, 键值多重集原样平移——窄键域内
/// `K::to_wide()` 与 KmerId 低 limb 数值恒等（kmer.rs 数值一致性单测）。
fn widen_to_count_map<K: RollingKey>(m: FxHashMap<K, u32>) -> CountMap {
    let mut out = FxHashMap::with_capacity_and_hasher(m.len(), Default::default());
    for (key, cnt) in m {
        out.insert(key.to_wide(), cnt);
    }
    out
}

/// inchworm 轮计数的分层入口: k≤32 用 u64 键、k≤64 用 u128 键在窄键空间
/// 完成 per-thread 计数与 reduce 合并（键带宽与槽位减半~四分之一, 大表
/// cache miss 与 FxHash 输入同步缩）, 终点一次 widen 回 CountMap——管线
/// 既有消费方（trusted 注入 / counts_bin 落盘 / KmerCounter::from_count_map）
/// 接口与语义不变。等价性: 各窄键层与 KmerId 层走同一泛型代码, 键数值
/// 域内恒等 → 计数多重集相等（单测 count_layered_matches_wide 全 k 对拍）。
/// 迭代序差异不构成回归: 并行 reduce 的合并序本就非确定（inchworm 并行
/// 平局序非确定, backlog 已立案）。
pub fn count_fasta_data_layered(data: &[u8], k: usize, ds: bool, per_thread_capacity: usize) -> CountMap {
    if k <= 32 {
        widen_to_count_map(count_fasta_data_of::<u64>(
            data,
            k,
            ds,
            per_thread_capacity,
        ))
    } else if k <= 64 {
        widen_to_count_map(count_fasta_data_of::<u128>(
            data,
            k,
            ds,
            per_thread_capacity,
        ))
    } else {
        count_fasta_data_of::<KmerId>(data, k, ds, per_thread_capacity)
    }
}

/// reduce 合并: b 并入 a（计数累加）。小表并入大表（rehash 次数最少化）。
pub fn merge_counts<K: RollingKey>(
    mut a: FxHashMap<K, u32>,
    mut b: FxHashMap<K, u32>,
) -> FxHashMap<K, u32> {
    if a.len() < b.len() {
        std::mem::swap(&mut a, &mut b);
    }
    for (key, cnt) in b {
        *a.entry(key).or_insert(0) += cnt;
    }
    a
}

/// trusted contig 注入（rnaSPAdes construction.cpp:233 的 Trinity 移植——
/// 上一轮组装的简化 contig 并入本轮 k-mer 计数流: **贡献序列、不贡献覆盖**）。
/// 表中已存在的键**不动其计数**; 缺席的键以 count=1 插入（presence-only:
/// 能被 min_count=1 的剪枝保留, 且在 inchworm 贪心延伸中优先级最低——
/// 恰是"只补序列不扰覆盖"）。注入 k-mer 不获任何豁免: 后续熵过滤等剪枝
/// 对它们与 reads k-mer 同等处理（忠实接受低复杂度注入 kmer 被剪掉）。
/// 调用方须保证 `k`/`ds` 与计数表的构建参数一致（注入发生在
/// prune_some_kmers 之前的最终表上）。返回新插入的键数（日志用）。
pub fn inject_trusted_kmers(counts: &mut CountMap, trusted_fasta: &[u8], k: usize, ds: bool) -> usize {
    let mut inserted = 0usize;
    for rec in split_fasta_sequences(trusted_fasta) {
        for_each_kmer(&rec, k, ds, |key| {
            // entry().or_insert(1): 已存在（reads 已覆盖或上轮已注入）→ 不动;
            // 缺席 → count=1 占位（presence-only）。
            if let std::collections::hash_map::Entry::Vacant(e) = counts.entry(key) {
                e.insert(1);
                inserted += 1;
            }
        });
    }
    inserted
}

/// mer_counts 二进制记录布局: `U256 key, u32 count`（小端, 36B/条）, 头部魔数
/// `TRM3` + u32 k + u8 ds（读回时校验）。M7 消 dump 直桥的落盘格式,
/// 管线（orchestrate）与离线工具（xtask export-edges）共用。
///
/// 格式沿革: u64 时代 `TRMC` + 8B key（12B/条）; u128 扩容（2026-08-30）升
/// `TRM2` + 16B key（20B/条）; U256 扩容（2026-09-02）升 `TRM3` + 32B key。
/// 键宽变化后旧文件静默错位解析的风险不可接受 → 每次魔数升版硬拒绝
/// （"bad header" 报错）, 不做双格式读（counts_bin 是易失缓存, 重算即可）。
pub fn save_counts_bin(
    path: &std::path::Path,
    counts: &CountMap,
    k: usize,
    ds: bool,
) -> Result<(), trinity_common::error::CommonError> {
    let mut buf: Vec<u8> = Vec::new();
    buf.extend_from_slice(b"TRM3");
    buf.extend_from_slice(&(k as u32).to_le_bytes());
    buf.push(ds as u8);
    for (&key, &c) in counts.iter() {
        buf.extend_from_slice(&key.to_le_bytes());
        buf.extend_from_slice(&c.to_le_bytes());
    }
    std::fs::write(path, buf)?;
    Ok(())
}

/// save_counts_bin 的读回侧; k/ds 不匹配或头损坏报错（`TRMC`/`TRM2` 等
/// 旧格式同走 bad header 硬拒绝——见 save_counts_bin 的格式沿革注记）。
pub fn load_counts_bin(
    path: &std::path::Path,
    k: usize,
    ds: bool,
) -> Result<CountMap, trinity_common::error::CommonError> {
    use trinity_common::error::CommonError;
    let data = std::fs::read(path)?;
    if data.len() < 9 || &data[..4] != b"TRM3" {
        return Err(CommonError::Parse(format!(
            "corrupt k-mer table {}: bad header",
            path.display()
        )));
    }
    let k_in = u32::from_le_bytes([data[4], data[5], data[6], data[7]]) as usize;
    let ds_in = data[8] != 0;
    if k_in != k || ds_in != ds {
        return Err(CommonError::Parse(format!(
            "k-mer table {} was built with different k/DS settings",
            path.display()
        )));
    }
    let mut counts = CountMap::default();
    let mut i = 9;
    while i + 36 <= data.len() {
        let key = KmerId::from_le_bytes(data[i..i + 32].try_into().unwrap());
        let c = u32::from_le_bytes(data[i + 32..i + 36].try_into().unwrap());
        counts.insert(key, c);
        i += 36;
    }
    Ok(counts)
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
    use trinity_common::kmer::{decode_kmer_from_intval, get_ds_kmer_val, kmer_to_intval};

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

    /// u128 紧凑键计数与 KmerId 全宽计数多重集等价（k≤64 域, hi 恒 0;
    /// 性能分层 2026-09-02 的数值正确性锁）。
    #[test]
    fn u128_compact_counting_matches_kmerid() {
        let seq = b">x\nACGTACGTNNAGGCACGTACGTACGTagct\n>y\nAAAAACCCCCGGGGGTTTTTACGTACGTA";
        for k in [15usize, 25, 31, 63, 64] {
            for ds in [true, false] {
                let wide = KmerCountTable::count_fasta_data(seq, k, ds);
                let narrow = count_fasta_data_of::<u128>(seq, k, ds, 0);
                assert_eq!(wide.len(), narrow.len(), "k={k} ds={ds}");
                for (key, c) in &narrow {
                    assert_eq!(wide.get(&key.to_wide()), Some(c), "k={k} ds={ds}");
                }
            }
        }
    }

    /// 分层计数入口（count_fasta_data_layered, 压时间战役 2026-09-04 刀1）与
    /// KmerId 直计在全部键层（u64 k≤32 / u128 k≤64 / KmerId 直通）多重集等价。
    /// 多行记录 fixture 额外覆盖 split_fasta_ranges 的 multiline 回退分支。
    #[test]
    fn count_layered_matches_wide() {
        let single = b">x\nACGTACGTNNAGGCACGTACGTACGTagct\n>y\nAAAAACCCCCGGGGGTTTTTACGTACGTA";
        let multiline: Vec<u8> = {
            let mut v = b">m\n".to_vec();
            for _ in 0..10 {
                v.extend_from_slice(b"ACGTACGT"); // 80 base, 分两行
                v.push(b'\n');
            }
            v.extend_from_slice(b"ACGTACGTNNAGG");
            v
        };
        for (tag, seq) in [("single", &single[..]), ("multiline", &multiline[..])] {
            for k in [7usize, 25, 31, 32, 33, 49, 63, 64, 65] {
                for ds in [true, false] {
                    let wide = KmerCountTable::count_fasta_data(seq, k, ds);
                    let layered = count_fasta_data_layered(seq, k, ds, 0);
                    assert_eq!(layered.len(), wide.len(), "{tag} k={k} ds={ds}");
                    for (key, c) in &wide {
                        assert_eq!(layered.get(key), Some(c), "{tag} k={k} ds={ds}");
                    }
                }
            }
        }
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

    /// M7-norm oracle: 滚动 revcomp 路径产出的 DS 键与逐 kmer get_ds_kmer_val 恒等
    /// （含 N 断窗、小写、滑窗恢复、k=1/k=2 边界; u128 扩容后加 k=33/37——
    /// 滚动路径在 u64 时代对 k>32 回落, 现对 k≤64 全程生效）。
    #[test]
    fn rolling_revcomp_matches_get_ds_kmer_val() {
        let seqs = [
            "ACGTACGTACGTACGTACGTACGTACG",
            "acgtNNNNacgtacgtacgtacgtacg",
            "GGGGG",
            "TTTTGGGGCCCCAAAATTTTGGGGCACGTN",
            // 70bp 长序列覆盖 k=33/37 的多窗滑窗（含首尾窗口的高位碱基）
            "ACGTACGTTTGCAGTCGATCGTACGATCGATGCATCGATCGTACGATGCTAGCTAGCATCGA",
        ];
        for &seq in &seqs {
            for k in [1usize, 2, 4, 25, 33, 37] {
                if k > seq.len() {
                    continue;
                }
                let mut rolling: Vec<KmerId> = Vec::new();
                for_each_kmer(seq.as_bytes(), k, true, |key| rolling.push(key));
                // naive: 逐窗口 kmer_to_intval + get_ds_kmer_val
                let b = seq.as_bytes();
                let mut naive: Vec<KmerId> = Vec::new();
                let mut i = 0;
                while i + k <= b.len() {
                    if let Ok(v) = trinity_common::kmer::kmer_to_intval(&b[i..i + k]) {
                        naive.push(get_ds_kmer_val(v, k));
                    }
                    i += 1;
                }
                assert_eq!(rolling, naive, "seq={seq} k={k}");
            }
        }
    }

    /// counts_bin TRM2（u128 键）落盘/读回往返, 含 k=37 的高位键;
    /// u64 时代 TRMC 旧格式被硬拒绝（防静默错位解析）。
    #[test]
    fn counts_bin_trm2_roundtrip_and_legacy_rejected() {
        let dir = std::env::temp_dir().join(format!("trmc2_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bin = dir.join("m.bin");
        let mut m: CountMap = FxHashMap::default();
        let hi = kmer_to_intval(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTA").unwrap(); // 37-mer
        assert!(hi > KmerId::from(u64::MAX as u128)); // 确为 u64 放不下的键
        m.insert(hi, 7);
        m.insert(kmer_to_intval(b"AAAA").unwrap(), 3);
        save_counts_bin(&bin, &m, 37, true).unwrap();
        let back = load_counts_bin(&bin, 37, true).unwrap();
        assert_eq!(back.len(), 2);
        assert_eq!(back[&hi], 7);
        // k/ds 不匹配 → 报错
        assert!(load_counts_bin(&bin, 25, true).is_err());
        assert!(load_counts_bin(&bin, 37, false).is_err());
        // 旧 TRMC 魔数（u64 版）→ bad header 硬拒绝
        let mut legacy = b"TRMC".to_vec();
        legacy.extend_from_slice(&37u32.to_le_bytes());
        legacy.push(1);
        std::fs::write(&bin, &legacy).unwrap();
        assert!(load_counts_bin(&bin, 37, true).is_err());
        let _ = std::fs::remove_dir_all(&dir);
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

#[cfg(test)]
mod trusted_inject_tests {
    use super::*;
    use trinity_common::kmer::kmer_to_intval;

    fn key(s: &str) -> KmerId {
        kmer_to_intval(s.as_bytes()).unwrap()
    }

    /// 核心语义: 已存在的键不动计数, 缺席的键以 count=1 插入。
    #[test]
    fn existing_untouched_absent_inserted_one() {
        // reads 表: ACG=7; trusted 携带 ACG/CGT/GTA/TAC（k=3, SS）
        let mut m: CountMap = FxHashMap::default();
        m.insert(key("ACG"), 7);
        let n = inject_trusted_kmers(&mut m, b">c1\nACGTAC\n", 3, false);
        assert_eq!(n, 3); // CGT/GTA/TAC 新插入
        assert_eq!(m[&key("ACG")], 7); // 已存在 → 不动（不累加、不覆盖）
        assert_eq!(m[&key("CGT")], 1);
        assert_eq!(m[&key("GTA")], 1);
        assert_eq!(m[&key("TAC")], 1);
        assert_eq!(m.len(), 4);
    }

    /// 幂等: 重复注入（resume 两条构建路径汇合后）不改变表。
    #[test]
    fn injection_is_idempotent() {
        let mut m: CountMap = FxHashMap::default();
        let fa = b">c1\nACGTAC\n";
        assert_eq!(inject_trusted_kmers(&mut m, fa, 3, false), 4);
        // 第二次注入: 全部已存在（含上轮注入的 count=1）→ 一个都不动
        assert_eq!(inject_trusted_kmers(&mut m, fa, 3, false), 0);
        assert_eq!(m.len(), 4);
        assert!(m.values().all(|&c| c == 1));
    }

    /// DS canonical 折叠与计数路径同键: reads 里的 CTT 与 trusted 里的 AAG
    /// 折叠为同一 DS 键 → 已存在不动。
    #[test]
    fn ds_canonical_folding_matches_reads_table() {
        let mut m: CountMap = FxHashMap::default();
        for_each_kmer(b"CTT", 3, true, |kk| *m.entry(kk).or_insert(0) += 5);
        assert_eq!(m.len(), 1);
        let n = inject_trusted_kmers(&mut m, b">c1\nAAG\n", 3, true);
        assert_eq!(n, 0); // AAG 的 DS canonical == CTT 的 canonical → 已存在
        assert_eq!(m.values().copied().collect::<Vec<_>>(), vec![5]);
    }

    /// 注入按**本轮** k 枚举（trusted 由上一轮较小 k 产出）;
    /// 记录短于 k → 无 k-mer 可注入。
    #[test]
    fn inject_uses_current_k_not_contig_build_k() {
        let mut m: CountMap = FxHashMap::default();
        assert_eq!(inject_trusted_kmers(&mut m, b">c1\nACGT\n", 5, false), 0);
        assert!(m.is_empty());
        assert_eq!(inject_trusted_kmers(&mut m, b">c1\nACGT\n", 3, false), 2);
    }

    /// 多记录 + N 断窗语义与计数路径一致（同一 for_each_kmer）。
    #[test]
    fn multi_record_and_n_breaks() {
        let mut m: CountMap = FxHashMap::default();
        let fa = b">c1\nACG\n>c2 desc\nTTNAC\n";
        let n = inject_trusted_kmers(&mut m, fa, 2, false);
        // c1: AC CG; c2 "TTNAC": TT 有效, TN/NA 含 N 失效, N 后窗口恢复 → AC 有效
        assert_eq!(n, 3);
        assert_eq!(m[&key("AC")], 1);
        assert_eq!(m[&key("CG")], 1);
        assert_eq!(m[&key("TT")], 1);
    }
}
