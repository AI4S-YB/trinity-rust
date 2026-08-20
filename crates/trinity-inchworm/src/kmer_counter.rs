//! KmerCounter — 直译 Inchworm/src/KmerCounter.cpp 的 kmer 计数目录（单线程版）。
//!
//! 语义要点（对照原版行号）:
//! - 键 = kmer_int_type_t(u64)，值 = unsigned int(u32) 计数（KC:13）
//! - DS 模式下 add/find/clear 全部先 canonical: max(k, revcomp(k))
//!   （KC:367-375 find_kmer / KC:476-480 add_kmer / KC:420-436 clear_kmer）
//! - 删除是惰性的: 置 0 不 erase（KC:430 `it->second = 0`），size() 含 0 值键（KC:51-53）
//! - clear_kmer 对不存在的键不做任何事（KC:428 `if (it != end)` 才置 0）——不插入
//! - 候选生成是纯位运算（KC:518-575）:
//!   forward 前缀 = (kmer << (33-K)*2) >> (32-K)*2，候选 = prefix | i；
//!   reverse 后缀 = kmer >> 2，候选 = (i << (2K-2)) | suffix。
//!   i = 0..3 恰为 G,A,T,C 序（编码 G=0,A=1,T=2,C=3）；过滤 count == 0（KC:528）
//! - 候选排序（KC:823-831 发布版比较器只有 count 降序；std::sort 不稳定）:
//!   此处显式化为稳定排序，平局保持 G,A,T,C 收集序（计划 Task 1 规格选择）

use rustc_hash::FxHashMap;
use trinity_common::error::CommonError;
use trinity_common::kmer::{
    decode_kmer_from_intval, get_ds_kmer_val, kmer_to_intval, KmerId, MAX_KMER_LENGTH,
};

/// kmer 计数目录的**只读**视图 —— KmerCounter（单线程）与 SyncKmerCounter
/// （dashmap，PARALLEL_IWORM 组装期，counter_sync.rs）共同实现，贪心延伸核心
/// （is_good_seed_kmer / build_inchworm_contig_from_seed / extract_best_seed /
/// reconstruct_path_sequence）以 `&impl KmerCatalog` 泛型在两实现上复用同一
/// 直译逻辑（原版 C++ 只有一份 KmerCounter& 代码，多线程直接共享同一对象）。
///
/// 只含读侧方法（组装期贪心核心从不写目录）:
/// - clear_kmer 不入 trait——单线程版 `&mut self`（FxHashMap get_mut）、
///   PARALLEL 版 `&self`（dashmap 原子置 0），接收者不同型，各自保留固有方法，
///   清零发生在主循环而非贪心核心；
/// - add_kmer / get_kmer_intval / iter_* 只在单线程装载/剪枝阶段使用，留在
///   KmerCounter 固有方法；
/// - 候选生成是纯位运算 + 4 次 get_kmer_count + 稳定排序（与存储结构无关），
///   作为默认方法两实现共用一份逻辑。
pub trait KmerCatalog {
    /// 键总数（含 count=0 的键——惰性删除不缩表，KC:51-53）
    fn size(&self) -> usize;

    fn get_kmer_length(&self) -> usize;

    /// DOUBLE_STRANDED_MODE 标志——IRKE 构造 Kmer_visitor 时需要（IRKE.cpp:775）
    fn is_double_stranded(&self) -> bool;

    /// KC:448-457 get_kmer_count: DS canonical 后查，缺省 0
    fn get_kmer_count(&self, kmer: KmerId) -> u32;

    /// KC:460-465 get_kmer_string: 按存储值解码（不做 canonical）
    fn get_kmer_string(&self, kmer: KmerId) -> Vec<u8>;

    /// KC:507-533 forward 候选（count>0 过滤 + count 降序稳定排序——默认实现）
    fn get_forward_kmer_candidates(&self, kmer: KmerId) -> Vec<(KmerId, u32)> {
        let kmer_length = self.get_kmer_length();
        // K=32 时 (33-K)*2=2、(32-K)*2=0；K=1 时左移 64 在原版是 UB（x86 硬件按
        // 移位量 mod 64 回绕），wrapping_shl 的回绕语义与之相同。
        let prefix = kmer.wrapping_shl(((33 - kmer_length) * 2) as u32) >> ((32 - kmer_length) * 2);
        let mut candidates = Vec::with_capacity(4);
        for i in 0..4 {
            let candidate = prefix | i;
            let count = self.get_kmer_count(candidate);
            if count != 0 {
                candidates.push((candidate, count));
            }
        }
        candidates.sort_by_key(|&(_, count)| std::cmp::Reverse(count)); // 稳定: 平局保持 G,A,T,C 收集序
        candidates
    }

    /// KC:549-575 reverse 候选（同上过滤与排序——默认实现）
    fn get_reverse_kmer_candidates(&self, kmer: KmerId) -> Vec<(KmerId, u32)> {
        let kmer_length = self.get_kmer_length();
        let suffix = kmer >> 2;
        let top_shift = kmer_length * 2 - 2;
        let mut candidates = Vec::with_capacity(4);
        for i in 0..4 {
            let candidate = (i << top_shift) | suffix;
            let count = self.get_kmer_count(candidate);
            if count != 0 {
                candidates.push((candidate, count));
            }
        }
        candidates.sort_by_key(|&(_, count)| std::cmp::Reverse(count)); // 稳定: 平局保持 G,A,T,C 收集序
        candidates
    }
}

/// 单线程 kmer 目录。PARALLEL 版见 counter_sync.rs（Task 6）。
#[derive(Debug)]
pub struct KmerCounter {
    kmer_length: usize,
    ds_mode: bool,
    counter: FxHashMap<KmerId, u32>,
}

impl KmerCounter {
    /// KC:13-17: kmer_length > 32 抛错（64bit / 2bit 编码上限），消息原样保留。
    pub fn new(kmer_length: usize, ds_mode: bool) -> Self {
        assert!(
            kmer_length <= MAX_KMER_LENGTH,
            "Kmer length exceeds max of 32"
        );
        KmerCounter {
            kmer_length,
            ds_mode,
            counter: FxHashMap::default(),
        }
    }

    /// DS canonical（KC:369-372 的 find_kmer 前置折叠，add/clear 同）
    fn canonical(&self, kmer: KmerId) -> KmerId {
        if self.ds_mode {
            get_ds_kmer_val(kmer, self.kmer_length)
        } else {
            kmer
        }
    }

    /// 键总数（含 count=0 的键——惰性删除不缩表，KC:51-53）
    pub fn size(&self) -> usize {
        self.counter.len()
    }

    pub fn get_kmer_length(&self) -> usize {
        self.kmer_length
    }

    /// DOUBLE_STRANDED_MODE 标志——IRKE 构造 Kmer_visitor 时需要（IRKE.cpp:775）
    pub fn is_double_stranded(&self) -> bool {
        self.ds_mode
    }

    /// KC:476-488 add_kmer(intval, count): DS canonical 后 `counter[key] += count`
    /// （unsigned int 溢出回绕 → wrapping_add）
    pub fn add_kmer(&mut self, kmer: KmerId, count: u32) {
        let entry = self.counter.entry(self.canonical(kmer)).or_insert(0);
        *entry = entry.wrapping_add(count);
    }

    /// KC:448-457 get_kmer_count: DS canonical 后查，缺省 0
    pub fn get_kmer_count(&self, kmer: KmerId) -> u32 {
        self.counter
            .get(&self.canonical(kmer))
            .copied()
            .unwrap_or(0)
    }

    /// KC:420-436 clear_kmer: DS canonical; 键存在才置 0（不存在不插入，size 不变）
    pub fn clear_kmer(&mut self, kmer: KmerId) {
        if let Some(count) = self.counter.get_mut(&self.canonical(kmer)) {
            *count = 0;
        }
    }

    /// KC:460-465 get_kmer_string: 按存储值解码（不做 canonical）
    pub fn get_kmer_string(&self, kmer: KmerId) -> Vec<u8> {
        decode_kmer_from_intval(kmer, self.kmer_length)
    }

    /// KC:467-473 get_kmer_intval(string): **仅编码**（kmer_to_intval），不做 DS
    /// 折叠——canonical 发生在 add_kmer（KC:479）。非 gatc → Err（原版 throw）。
    pub fn get_kmer_intval(&self, seq: &[u8]) -> Result<KmerId, CommonError> {
        kmer_to_intval(seq)
    }

    /// count > 0 的 (kmer, count)（get_kmers_sort_descending_counts 收集阶段，KC:719-728）
    pub fn iter_nonzero(&self) -> impl Iterator<Item = (KmerId, u32)> + '_ {
        self.counter
            .iter()
            .filter(|&(_, &count)| count > 0)
            .map(|(&kmer, &count)| (kmer, count))
    }

    /// 全键迭代（**含 count=0 的键**）——C++ hash_map 迭代域，map 迭代序。
    /// prune_some_kmers（irke.rs）用它镜像 KC:143 的全表遍历。
    pub fn iter_all(&self) -> impl Iterator<Item = (KmerId, u32)> + '_ {
        self.counter.iter().map(|(&kmer, &count)| (kmer, count))
    }

    /// 拆出内部存储（kmer_length, ds_mode, 整表）——SyncKmerCounter::from_counter
    /// （counter_sync.rs）把装载+剪枝完成的单线程目录整表转入并发结构时使用。
    pub(crate) fn into_parts(self) -> (usize, bool, FxHashMap<KmerId, u32>) {
        (self.kmer_length, self.ds_mode, self.counter)
    }
}

/// 只读视图委托（候选两方法用 trait 默认实现——位运算 + get_kmer_count 共用一份）
impl KmerCatalog for KmerCounter {
    fn size(&self) -> usize {
        self.counter.len()
    }

    fn get_kmer_length(&self) -> usize {
        self.kmer_length
    }

    fn is_double_stranded(&self) -> bool {
        self.ds_mode
    }

    fn get_kmer_count(&self, kmer: KmerId) -> u32 {
        self.counter
            .get(&self.canonical(kmer))
            .copied()
            .unwrap_or(0)
    }

    fn get_kmer_string(&self, kmer: KmerId) -> Vec<u8> {
        decode_kmer_from_intval(kmer, self.kmer_length)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 测试辅助: 替换向量末位碱基（K=32 断言用，避免到处 concat）。
    trait WithLast {
        fn with_last(self, b: u8) -> Vec<u8>;
    }
    impl WithLast for Vec<u8> {
        fn with_last(mut self, b: u8) -> Vec<u8> {
            let n = self.len();
            self[n - 1] = b;
            self
        }
    }

    /// 便捷编码（测试内手算黄金值都用裸数字并注释推导）
    fn enc(s: &[u8]) -> KmerId {
        kmer_to_intval(s).unwrap()
    }

    // ---- 候选位运算手算（K=2） -------------------------------------------

    #[test]
    fn forward_candidates_k2_hand_derived() {
        // seed "AC" = (1<<2)|3 = 7。
        // prefix = (7 << (33-2)*2) >> (32-2)*2 = (7 << 62) >> 60
        //        = 0xC000_0000_0000_0000 >> 60 = 12 = "CG"
        // 候选 = 12|i（i=0..3）→ CG(12), CA(13), CT(14), CC(15)，恰为 G,A,T,C 序
        // 语义: "AC" 分别向 3' 延伸 G/A/T/C 后取最后 2 碱基
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(12, 5); // CG
        kc.add_kmer(13, 5); // CA
        kc.add_kmer(14, 1); // CT
                            // CC 不存在 → count 0 → 被过滤
        let cands = kc.get_forward_kmer_candidates(7);
        assert_eq!(cands, vec![(12, 5), (13, 5), (14, 1)]);
        // 平局 5:5 保持 GATC 收集序（CG 的 i=0 在 CA 的 i=1 前）
        assert_eq!(kc.get_kmer_string(cands[0].0), b"CG".to_vec());
        assert_eq!(kc.get_kmer_string(cands[1].0), b"CA".to_vec());
        assert_eq!(kc.get_kmer_string(cands[2].0), b"CT".to_vec());
    }

    #[test]
    fn reverse_candidates_k2_hand_derived() {
        // seed "AC" = 7。suffix = 7 >> 2 = 1。
        // 候选 = (i << (2*2-2)) | 1 = (i<<2)|1 → GA(1), AA(5), TA(9), CA(13)
        // 语义: "AC" 分别向 5' 前置 G/A/T/C 后取最前 2 碱基
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(1, 2); // GA
        kc.add_kmer(5, 7); // AA
        kc.add_kmer(9, 7); // TA
                           // CA 不存在 → 过滤
        let cands = kc.get_reverse_kmer_candidates(7);
        // 收集序 [GA:2, AA:7, TA:7] → count 降序稳定排序: AA, TA（平局 AA 的 i=1 在 TA 的 i=2 前）, GA
        assert_eq!(cands, vec![(5, 7), (9, 7), (1, 2)]);
        assert_eq!(kc.get_kmer_string(cands[0].0), b"AA".to_vec());
        assert_eq!(kc.get_kmer_string(cands[1].0), b"TA".to_vec());
    }

    #[test]
    fn ds_mode_candidates_lookup_canonicalizes() {
        // DS 模式: 候选 id 是原始位运算值，但计数查询经 canonical（KC:527 → get_kmer_count）
        // seed "GT" = 2（canonical 为 "AC"=7，但前缀按传入的原始 2 计算）
        // prefix = (2 << 62) >> 60 = 8 → 候选 TG(8), TA(9), TT(10), TC(11)
        // canonical: TG(8)→CA(13)（revcomp("TG")="CA"）; TA(9) 回文自映; TT(10)→TT(10)… max(10,5)=10
        let mut kc = KmerCounter::new(2, true);
        kc.add_kmer(enc(b"CA"), 7); // canonical 13
        kc.add_kmer(enc(b"TA"), 2); // canonical 9（回文）
        let cands = kc.get_forward_kmer_candidates(2);
        // TG(8) 查到 canonical CA 的 7；TA(9) 查到 2；TT/TC 为 0 被过滤
        assert_eq!(cands, vec![(8, 7), (9, 2)]);
    }

    // ---- K=32 边界 -------------------------------------------------------

    #[test]
    fn k32_forward_boundary_no_panic() {
        // K=32: (33-32)*2 = 2，(32-32)*2 = 0 —— 退化为 (kmer << 2) | i
        let seed = enc(&[b'A'; 32]); // 0x5555_5555_5555_5555
        assert_eq!(seed, 0x5555_5555_5555_5555);
        let prefix = seed << 2; // 顶 2 位（首个 A）被移出
        assert_eq!(prefix, 0x5555_5555_5555_5554);
        let mut kc = KmerCounter::new(32, false);
        kc.add_kmer(prefix, 4); // A*31 + G
        kc.add_kmer(prefix | 1, 9); // A*31 + A
        kc.add_kmer(prefix | 2, 4); // A*31 + T
                                    // prefix|3（A*31+C）不存在 → 过滤
        let fwd = kc.get_forward_kmer_candidates(seed);
        assert_eq!(fwd, vec![(prefix | 1, 9), (prefix, 4), (prefix | 2, 4)]);
        // 平局 4:4 → GATC 序（|0 在 |2 前）
        assert_eq!(kc.get_kmer_string(fwd[1].0), vec![b'A'; 32].with_last(b'G'));
    }

    #[test]
    fn k32_reverse_boundary_no_panic() {
        // K=32: 候选 = (i << 62) | (seed >> 2)。前置 A 时 = seed 本身（A*32 自指）
        let seed = enc(&[b'A'; 32]);
        let suffix = seed >> 2;
        let mut kc = KmerCounter::new(32, false);
        kc.add_kmer(seed, 9); // A*32（前置 A）
        kc.add_kmer(suffix, 2); // G + A*31（前置 G，i=0）
        kc.add_kmer((2u64 << 62) | suffix, 2); // T + A*31（前置 T，i=2）
        let rev = kc.get_reverse_kmer_candidates(seed);
        assert_eq!(
            rev,
            vec![(seed, 9), (suffix, 2), ((2u64 << 62) | suffix, 2)]
        );
        // 平局 2:2 → GATC 序（前置 G 在前置 T 前）
        let mut expected = vec![b'G'];
        expected.extend(vec![b'A'; 31]);
        assert_eq!(kc.get_kmer_string(rev[1].0), expected);
    }

    // ---- DS canonical / 累加 / 回绕 --------------------------------------

    #[test]
    fn ds_mode_canonicalizes_add_and_lookup() {
        // "GT"=2 → canonical = max(2, revcomp("GT")=7) = 7（"AC"）
        let mut kc = KmerCounter::new(2, true);
        kc.add_kmer(enc(b"GT"), 3);
        assert_eq!(kc.size(), 1);
        assert_eq!(kc.get_kmer_count(enc(b"GT")), 3);
        assert_eq!(kc.get_kmer_count(enc(b"AC")), 3);
        // 互补链累加到同一键
        kc.add_kmer(enc(b"AC"), 2);
        assert_eq!(kc.get_kmer_count(enc(b"GT")), 5);
        assert_eq!(kc.size(), 1);
        // clear 经 canonical 生效（从 GT 清，AC 键被置 0 但仍在表中）
        kc.clear_kmer(enc(b"GT"));
        assert_eq!(kc.get_kmer_count(enc(b"AC")), 0);
        assert_eq!(kc.size(), 1);
    }

    #[test]
    fn add_accumulates_and_wraps() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(7, u32::MAX);
        kc.add_kmer(7, 1); // 原版 unsigned int += 溢出回绕 → 0
        assert_eq!(kc.get_kmer_count(7), 0);
        assert_eq!(kc.size(), 1); // 回绕后键仍占位
    }

    // ---- 惰性删除 / size 语义 --------------------------------------------

    #[test]
    fn clear_kmer_is_lazy_and_size_stable() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(7, 3);
        // 不存在的键: KC:428 find==end → 什么都不做（不插入！）
        kc.clear_kmer(9);
        assert_eq!(kc.size(), 1);
        assert_eq!(kc.get_kmer_count(9), 0);
        // 存在的键: 置 0，size 不变（不 erase）
        kc.clear_kmer(7);
        assert_eq!(kc.get_kmer_count(7), 0);
        assert_eq!(kc.size(), 1);
        assert_eq!(kc.iter_nonzero().count(), 0);
    }

    #[test]
    fn iter_nonzero_skips_zero_counts() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(7, 3);
        kc.add_kmer(9, 0); // add 0 → 键存在值为 0
        kc.add_kmer(12, 2);
        kc.clear_kmer(12); // 置 0 后不再出现
        kc.add_kmer(13, 1);
        let mut v: Vec<_> = kc.iter_nonzero().collect();
        v.sort_unstable();
        assert_eq!(v, vec![(7, 3), (13, 1)]);
        assert_eq!(kc.size(), 4); // 7, 9, 12, 13 全部在表（含 0 值）
    }

    // ---- 编码 / 解码 / 构造护栏 ------------------------------------------

    #[test]
    fn get_kmer_intval_encodes_without_ds_fold() {
        // KC:467-473: get_kmer_intval 只编码，不 canonical——DS 折叠在 add_kmer。
        // （IRKE.cpp:128 的调用序正是 get_kmer_intval → add_kmer，折叠由后者完成）
        let ss = KmerCounter::new(2, false);
        assert_eq!(ss.get_kmer_intval(b"AC").unwrap(), 7);
        let ds = KmerCounter::new(2, true);
        assert_eq!(ds.get_kmer_intval(b"GT").unwrap(), 2); // 原值，不折叠
        assert_eq!(ds.get_kmer_intval(b"AC").unwrap(), 7);
        assert!(ds.get_kmer_intval(b"GN").is_err());
        assert_eq!(ds.get_kmer_length(), 2);
    }

    #[test]
    fn iter_all_includes_zero_counts() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(7, 3);
        kc.add_kmer(9, 0);
        kc.clear_kmer(7);
        assert_eq!(kc.size(), 2);
        let mut v: Vec<_> = kc.iter_all().collect();
        v.sort_unstable();
        assert_eq!(v, vec![(7, 0), (9, 0)]); // 含 0 值键
        assert_eq!(kc.iter_nonzero().count(), 0);
    }

    #[test]
    #[should_panic(expected = "Kmer length exceeds max of 32")]
    fn kmer_length_over_32_rejected() {
        let _ = KmerCounter::new(33, false);
    }
}
