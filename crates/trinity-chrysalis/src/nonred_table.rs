//! NonRedKmerTable 移植（P3-T1 基础层）。
//! 镜像 `Chrysalis/analysis/{NonRedKmerTable.h,NonRedKmerTable.cc}`。
//!
//! 结构 = 排序去重的 k-mer 字符串数组（`UniqueSort` 语义）+ 平行 counts；
//! 查询走二分。双语义复用（原版两处调用点）：
//! - **weldmer 计数**（GraphFromFasta）：templates = 候选 48-mer 的两半 24-mer，
//!   `add_counts_from_reads` 累加 read 支持；
//! - **k-mer → 组件索引**（CreateIwormFastaBundle / ReadsToTranscripts）：
//!   `set_count` 直接写入组件号，`get_count` 读回。
//!
//! 已证差异（有意为之）：模板/读长度 < k 的窗口，原版因 `size_t` 下溢
//! `d.size()-m_k` 越界读（UB），本版安全跳过。

/// 排序字符串数组 + 二分的非冗余 k-mer 表。
use rustc_hash::FxHashMap;
use trinity_common::kmer::{base_to_int, kmer_to_intval, KmerId};

/// NonRedKmerTable（NonRedKmerTable.h/.cc 直译）。
/// 表 = 排序去重的 k-mer 集（字节序）+ 同长计数; 查询 = 二分（`index`）。
/// 性能扩展（2026-09-01）: 纯 GATC k-mer 懒建 2bit 哈希索引（`hash_index`），
/// `add_counts_from_reads_par` 的逐窗查询从 O(log n)+memcmp 降为 O(1)；
/// 含非 GATC 的 kmer 不进索引（编码无定义）, 查询侧含非 GATC 的窗走
/// 字节二分回退——与原版逐字节语义完全一致（原版对 N 也做字节相等匹配）。
pub struct NonRedKmerTable {
    k: usize,
    data: Vec<Vec<u8>>,
    counts: Vec<i32>,
    /// 懒建哈希索引: 纯 GATC kmer 的 2bit 值 → data 下标（OnceLock 一次建成）。
    hash_idx: std::sync::OnceLock<FxHashMap<KmerId, u32>>,
}

impl NonRedKmerTable {
    /// 镜像 `SetUp(templ, noNs)`（NonRedKmerTable.cc:12-96）：
    /// 收集全部 k-mer 窗口（noNs=true 时跳过含非**大写** ACGT 的窗口——
    /// `Regular()` 只认大写），字典序排序 + 原地去重，counts 清零。
    ///
    /// 性能（2026-09-01）: 窗口平铺收集（大 Vec + 偏移表, 免逐窗堆分配）,
    /// 对偏移表按窗口字节序排序（移动 u32 而非 Vec）, dedup 后重建——
    /// 与逐 Vec 排序的产出逐字节一致（同排序键, 同去重）。
    pub fn set_up_templates(templ: &[Vec<u8>], k: usize, no_ns: bool) -> Self {
        let mut flat: Vec<u8> = Vec::new();
        let mut offs: Vec<u32> = Vec::new();
        for d in templ {
            if d.len() < k {
                continue; // 原版 size_t 下溢越界读（见模块文档）
            }
            for j in 0..=d.len() - k {
                let w = &d[j..j + k];
                if no_ns && w.iter().any(|&c| !matches!(c, b'A' | b'C' | b'G' | b'T')) {
                    continue;
                }
                offs.push(flat.len() as u32);
                flat.extend_from_slice(w);
            }
        }
        offs.sort_unstable_by(|&a, &b| {
            flat[a as usize..a as usize + k].cmp(&flat[b as usize..b as usize + k])
        });
        let mut data: Vec<Vec<u8>> = Vec::with_capacity(offs.len());
        let mut prev: Option<&[u8]> = None;
        for &o in &offs {
            let w = &flat[o as usize..o as usize + k];
            if Some(w) == prev {
                continue;
            }
            data.push(w.to_vec());
            prev = Some(w);
        }
        let counts = vec![0i32; data.len()];
        NonRedKmerTable {
            k,
            data,
            counts,
            hash_idx: std::sync::OnceLock::new(),
        }
    }

    /// 懒建 2bit 哈希索引（纯 GATC kmer; 含非 GATC 者不进索引, 由 `index`
    /// 字节二分兜底）。
    fn hash_index(&self) -> &FxHashMap<KmerId, u32> {
        self.hash_idx.get_or_init(|| {
            self.data
                .iter()
                .enumerate()
                .filter_map(|(i, w)| kmer_to_intval(w).ok().map(|v| (v, i as u32)))
                .collect()
        })
    }

    pub fn k(&self) -> usize {
        self.k
    }

    /// 表内 k-mer 数。
    pub fn len(&self) -> usize {
        self.data.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// 二分定位（`Index`/`BinSearch` 语义）；窗口越界（pos+k > len）视为 miss
    /// （原版 substr 截短后必然搜不到）。
    fn index(&self, seq: &[u8], pos: usize) -> Option<usize> {
        if pos + self.k > seq.len() {
            return None;
        }
        let w = &seq[pos..pos + self.k];
        let i = self.data.partition_point(|x| x.as_slice() < w);
        if i < self.data.len() && self.data[i].as_slice() == w {
            Some(i)
        } else {
            None
        }
    }

    /// 镜像 `AddData(DNAStringStreamFast &)`（cc:162-200，omp 并行版本——
    /// 本版先单线程逐窗口 +1，rayon 并行 + 原子累加在 T2 接线）：
    /// 每 read **toupper 后**逐窗口二分，命中 +1，miss 跳过。
    pub fn add_counts_from_reads(&mut self, reads: &[Vec<u8>]) {
        for r in reads {
            if r.len() < self.k {
                continue;
            }
            let upper: Vec<u8> = r.iter().map(|&c| c.to_ascii_uppercase()).collect();
            for j in 0..=upper.len() - self.k {
                if let Some(i) = self.index(&upper, j) {
                    self.counts[i] += 1;
                }
            }
        }
    }

    /// **并行版** `add_counts_from_reads`（镜像 NonRedKmerTable.cc:161-200
    /// `AddData(DNAStringStreamFast&)` 的 omp atomic 语义）：分块局部计数 +
    /// 按块序合并——整数累加，与原版原子累加/与单线程版结果完全一致。
    /// 消费方（`graph_from_fasta`）等价于原版 `#pragma omp parallel` 抢读。
    ///
    /// 性能（2026-09-01）: 滚动 2-bit 编码 + 哈希查询（`hash_index`）替换
    /// 逐窗二分+整段大写拷贝; 含非 GATC 的窗走字节二分回退（原版 toupper
    /// 后字节匹配语义）。read 长度 < k 的跳过语义不变。
    pub fn add_counts_from_reads_par(&mut self, reads: &[Vec<u8>]) {
        use rayon::prelude::*;
        let n = self.data.len();
        if n == 0 || reads.is_empty() {
            return;
        }
        let hidx = self.hash_index(); // 先建索引（OnceLock, 一次）
        let k = self.k;
        let mask: KmerId = KmerId::low_bits(2 * k as u32);
        let data = &self.data;
        let chunk = (reads.len() / (rayon::current_num_threads() * 4)).max(1);
        let partials: Vec<Vec<i32>> = reads
            .par_chunks(chunk)
            .map(|ch| {
                let mut local = vec![0i32; n];
                for r in ch {
                    if r.len() < k {
                        continue;
                    }
                    let mut val: KmerId = KmerId::ZERO;
                    let mut run = 0usize;
                    for j in 0..r.len() {
                        match base_to_int(r[j]) {
                            Some(v) => {
                                val = ((val << 2) | KmerId::from(v)) & mask;
                                run += 1;
                            }
                            None => {
                                val = KmerId::ZERO;
                                run = 0;
                            }
                        }
                        if j + 1 < k {
                            continue;
                        }
                        if run >= k {
                            if let Some(&i) = hidx.get(&val) {
                                local[i as usize] += 1;
                            }
                        } else {
                            // 含非 GATC 窗: 大写后字节二分回退（原版语义）
                            let ws = j + 1 - k;
                            let w: Vec<u8> =
                                r[ws..ws + k].iter().map(|c| c.to_ascii_uppercase()).collect();
                            let i = data.partition_point(|x| x.as_slice() < w.as_slice());
                            if i < data.len() && data[i].as_slice() == w.as_slice() {
                                local[i] += 1;
                            }
                        }
                    }
                }
                local
            })
            .collect();
        for p in &partials {
            for (i, &c) in p.iter().enumerate() {
                self.counts[i] += c;
            }
        }
    }

    /// builder 形态：`set_up_templates(...).with_read_counts(reads)`
    /// （GraphFromFasta 的 `SetUp` + `AddData` 序列）。
    pub fn with_read_counts(mut self, reads: &[Vec<u8>]) -> Self {
        self.add_counts_from_reads_par(reads);
        self
    }

    /// `GetCount`（NonRedKmerTable.h:40-45）：miss 返回 **0**。
    pub fn get_count(&self, seq: &[u8], pos: usize) -> i32 {
        match self.index(seq, pos) {
            Some(i) => self.counts[i],
            None => 0,
        }
    }

    /// `GetCountReal`（h:47-52）：miss 返回 **-1**（可区分"真零计数"）。
    pub fn get_count_real(&self, seq: &[u8], pos: usize) -> i32 {
        match self.index(seq, pos) {
            Some(i) => self.counts[i],
            None => -1,
        }
    }

    /// `SetAllCounts`（h:68-71）。
    pub fn set_all_counts(&mut self, v: i32) {
        for c in self.counts.iter_mut() {
            *c = v;
        }
    }

    /// `SetCount`（h:33-38）：miss 是 **no-op**。
    pub fn set_count(&mut self, seq: &[u8], pos: usize, v: i32) {
        if let Some(i) = self.index(seq, pos) {
            self.counts[i] = v;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// setUp：两模板窗口并集 → 字典序排序去重；counts 清零。
    /// ACGTTGCA → {ACGT,CGTT,GTTG,TTGC,TGCA}；ACGTTGC → {ACGT,CGTT,GTTG,TTGC}
    /// → 去重后 5 个，字典序 ACGT < CGTT < GTTG < TGCA < TTGC。
    #[test]
    fn set_up_sorts_and_dedups() {
        let t = NonRedKmerTable::set_up_templates(
            &[b"ACGTTGCA".to_vec(), b"ACGTTGC".to_vec()],
            4,
            false,
        );
        assert_eq!(t.len(), 5);
        assert_eq!(t.k(), 4);
        // 全部计数为 0（get_count miss 语义 = 0，用 real 区分）
        for w in [b"ACGT", b"CGTT", b"GTTG", b"TGCA", b"TTGC"] {
            assert_eq!(t.get_count_real(w, 0), 0, "{:?}", w as &[u8]);
        }
        // 不在表内
        assert_eq!(t.get_count_real(b"AAAA", 0), -1);
    }

    /// 同模板内部重复窗口也去重（AAAAAA k=3 → 只 1 个 AAA）。
    #[test]
    fn set_up_dedup_within_template() {
        let t = NonRedKmerTable::set_up_templates(&[b"AAAAAA".to_vec()], 3, false);
        assert_eq!(t.len(), 1);
    }

    /// add_counts：toupper 后滑窗命中 +1；miss 窗口不计数。
    #[test]
    fn add_counts_uppercases_and_counts() {
        let mut t = NonRedKmerTable::set_up_templates(&[b"ACGTTGCA".to_vec()], 4, false);
        t.add_counts_from_reads(&[b"acgttgca".to_vec()]); // 小写 read → toupper
        for (w, expect) in [
            (&b"ACGT"[..], 1),
            (&b"CGTT"[..], 1),
            (&b"GTTG"[..], 1),
            (&b"TTGC"[..], 1),
            (&b"TGCA"[..], 1),
        ] {
            assert_eq!(t.get_count(w, 0), expect, "{:?}", w);
        }
        t.add_counts_from_reads(&[b"ACGTACGT".to_vec()]);
        // "ACGTACGT" 的窗口（len-k+1=5 个）里 ACGT 出现在 pos 0 与 pos 4 → 1+2=3
        assert_eq!(t.get_count(b"ACGT", 0), 3);
        assert_eq!(t.get_count(b"CGTT", 0), 1);
        assert_eq!(t.get_count(b"CGTA", 0), 0); // 读内窗口但不在模板 → miss → 0
    }

    /// pos 偏移查询（Index 用 seq[pos..pos+k]）与 pos 越界 miss。
    #[test]
    fn get_count_with_pos_offset() {
        let mut t = NonRedKmerTable::set_up_templates(&[b"ACGTTGCA".to_vec()], 4, false);
        t.add_counts_from_reads(&[b"ACGTTGCA".to_vec(), b"TTGCAAAA".to_vec()]);
        // "TTGC" 出现在 read1[3..7) 与 read2[0..4)；查询串 pos=3 起的窗口恰为 TTGC
        assert_eq!(t.get_count(b"AAATTGCAA", 3), 2); // pos 偏移命中
        assert_eq!(t.get_count(b"AAATTGCAA", 3 + 4), 0); // pos=7 越界 → 0
        assert_eq!(t.get_count_real(b"AAATTGCAA", 7), -1); // 越界 → -1
    }

    /// get_count 与 get_count_real 的 miss 语义差（0 vs -1）。
    #[test]
    fn get_count_vs_real_miss_semantics() {
        let t = NonRedKmerTable::set_up_templates(&[b"ACGT".to_vec()], 4, false);
        assert_eq!(t.get_count(b"TTTT", 0), 0);
        assert_eq!(t.get_count_real(b"TTTT", 0), -1);
        // 真零计数：real 也是 0
        assert_eq!(t.get_count_real(b"ACGT", 0), 0);
    }

    /// set_count：命中覆写；miss no-op；set_all_counts 全量覆写。
    #[test]
    fn set_count_and_set_all() {
        let mut t = NonRedKmerTable::set_up_templates(&[b"ACGTTG".to_vec()], 3, false);
        // ACGTTG k=3 → {ACG, CGT, GTT, TTG}
        t.set_count(b"ACGT", 0, 42); // 命中 ACG
        assert_eq!(t.get_count(b"ACGT", 0), 42);
        t.set_count(b"AAAA", 0, 99); // miss → no-op
        assert_eq!(t.get_count(b"AAAA", 0), 0);
        t.set_all_counts(7);
        for w in [&b"ACG"[..], &b"CGT"[..], &b"GTT"[..], &b"TTG"[..]] {
            assert_eq!(t.get_count(w, 0), 7);
        }
    }

    /// no_ns=true：跳过含非大写 ACGT 的窗口（N 与小写都算"不规则"）。
    #[test]
    fn no_ns_skips_irregular_windows() {
        let t = NonRedKmerTable::set_up_templates(&[b"ACGTNAAA".to_vec()], 4, true);
        // 窗口 ACGT(留) CGTN GTNA TNAA NAAA(全跳)
        assert_eq!(t.len(), 1);
        assert_eq!(t.get_count_real(b"ACGT", 0), 0);
        assert_eq!(t.get_count_real(b"CGTN", 0), -1);

        // 小写在模板中也算不规则（Regular() 只认大写）
        let t2 = NonRedKmerTable::set_up_templates(&[b"acgtACGT".to_vec()], 4, true);
        // acgt cgtA gtAC tACG(跳) ACGT(留)
        assert_eq!(t2.len(), 1);
        assert_eq!(t2.get_count_real(b"ACGT", 0), 0);
    }

    /// no_ns=false：含 N 窗口照常入表（计数语义用）。
    #[test]
    fn no_ns_false_keeps_n_windows() {
        let t = NonRedKmerTable::set_up_templates(&[b"ACGTN".to_vec()], 2, false);
        // len-k+1 = 4 个窗口：AC CG GT TN
        assert_eq!(t.len(), 4);
        assert_eq!(t.get_count_real(b"TN", 0), 0);
        assert_eq!(t.get_count_real(b"NT", 0), -1);
    }

    /// 短于 k 的模板/read 安全跳过（原版 size_t 下溢 UB）。
    #[test]
    fn shorter_than_k_is_skipped() {
        let t = NonRedKmerTable::set_up_templates(&[b"ACG".to_vec()], 4, false);
        assert!(t.is_empty());
        let mut t2 = NonRedKmerTable::set_up_templates(&[b"ACGT".to_vec()], 4, false);
        t2.add_counts_from_reads(&[b"AC".to_vec(), b"".to_vec()]); // 不 panic
        assert_eq!(t2.get_count_real(b"ACGT", 0), 0);
    }

    /// 并行版与单线程版计数逐项一致（整数累加无序无关）。
    #[test]
    fn par_counts_equal_serial_counts() {
        let templ: Vec<Vec<u8>> = (0..20)
            .map(|i| {
                let mut v = Vec::new();
                let mut s = i as u8;
                for _ in 0..30 {
                    s = s.wrapping_mul(31).wrapping_add(7);
                    v.push(b"ACGT"[s as usize % 4]);
                }
                v
            })
            .collect();
        let reads: Vec<Vec<u8>> = (0..200)
            .map(|i| {
                // 模板拼接成 read——窗口必与模板 6-mer 重合，保证有非零计数
                [
                    templ[i % templ.len()].clone(),
                    templ[(i * 7) % templ.len()].clone(),
                ]
                .concat()
            })
            .collect();
        let mut a = NonRedKmerTable::set_up_templates(&templ, 6, false);
        let mut b = NonRedKmerTable::set_up_templates(&templ, 6, false);
        a.add_counts_from_reads(&reads);
        b.add_counts_from_reads_par(&reads);
        for w in &a.data {
            assert_eq!(a.get_count(w, 0), b.get_count(w, 0));
        }
        assert!(a.data.iter().any(|w| a.get_count(w, 0) > 0));
    }
}
