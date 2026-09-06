//! SyncKmerCounter — PARALLEL_IWORM 组装期的并发 kmer 计数目录（Task 6）。
//!
//! 原版 IRKE.cpp:504-620 的 `#pragma omp parallel for schedule(dynamic, 1000)`
//! 循环里，多个 OpenMP 线程对**同一个** KmerCounter（__gnu_cxx::hash_map）
//! 并发 get_kmer_count / clear_kmer——C++ 侧没有任何同步（读-写竞态按标准是
//! UB，实践上单键 u32 置 0/读"碰巧"可用，Trinity 就这么跑了十几年）。
//!
//! 此移植用 dashmap（分片 RwLock 哈希表）承载组装期目录：
//! - **单键操作原子**（get / 置 0），**跨键无任何一致性**——另一线程可能同时
//!   读到已被置 0 前的非零值（或反之）。这正是原版无锁竞态的语义，只是消除了
//!   UB 与段错误风险（"弱一致"，PARALLEL 模式本就 nondeterministic，by-design）；
//! - 键/值/DS canonical/惰性删除语义与 KmerCounter 完全一致——组装期贪心核心
//!   只经 [`KmerCatalog`] 只读视图访问，两实现共用 irke.rs 的同一份直译逻辑；
//! - 装载与剪枝仍在单线程 KmerCounter 上完成，之后**整表转入**本结构
//!   （组装期不增不删键 → size 恒为转入时键数，与惰性删除的 size 语义一致）。

use dashmap::DashMap;
use rustc_hash::{FxBuildHasher, FxHashMap};

use trinity_common::kmer::{decode_kmer_from_intval, get_ds_kmer_val, KmerId};

use crate::kmer_counter::{
    ConcurrentClearCatalog, ExclusiveClearCatalog, KmerCatalog, KmerCounter,
};

/// 并发 kmer 目录。构造点唯一：[`SyncKmerCounter::from_counter`]
/// （单线程装载+剪枝完成的 KmerCounter 整表转入）。
pub struct SyncKmerCounter {
    kmer_length: usize,
    ds_mode: bool,
    counter: DashMap<KmerId, u32, FxBuildHasher>,
}

impl SyncKmerCounter {
    /// 单线程 KmerCounter（装载+剪枝完成态）整表转入——转移后原目录不可再写，
    /// 此结构即 PARALLEL 组装期目录的唯一定义点。
    ///
    /// 迭代序：种子快照走 [`SyncKmerCounter::nonzero_keys`]（DashMap 分片序，
    /// 2026-09-06 刀M6 起转入先于种子收集执行——from_counter 整表重建的
    /// 双表并存瞬时需避开亿级种子列表在场）, 分片/桶布局只影响性能与
    /// 种子收集序（非确定源之一, 平局序 backlog 已立案）不影响键值语义。
    pub fn from_counter(counter: KmerCounter) -> Self {
        let (kmer_length, ds_mode, map) = counter.into_parts();
        let counter = DashMap::with_capacity_and_hasher(map.len(), FxBuildHasher);
        for (kmer, count) in map {
            counter.insert(kmer, count);
        }
        SyncKmerCounter {
            kmer_length,
            ds_mode,
            counter,
        }
    }

    /// count > 0 的键快照（PARALLEL 种子列表, 刀M6）。只取键——count 主循环
    /// 不参与（实时重查, IRKE.cpp:546）; `with_capacity` 精确预留免倍增
    /// （末次翻倍 +11.7GB 瞬时, 见 collect_nonzero_seeds 注记）。
    pub fn nonzero_keys(&self) -> Vec<KmerId> {
        let mut out = Vec::with_capacity(self.counter.len());
        self.counter.iter().for_each(|item| {
            if *item.value() > 0 {
                out.push(*item.key());
            }
        });
        out
    }

    /// DS canonical（KC:369-372 find_kmer 前置折叠的镜像）
    fn canonical(&self, kmer: KmerId) -> KmerId {
        if self.ds_mode {
            get_ds_kmer_val(kmer, self.kmer_length)
        } else {
            kmer
        }
    }

    /// KC:420-436 clear_kmer 的并发版: DS canonical 后，键存在才**原子**置 0
    /// （不插入、不删键——get_mut 缺键即无操作，镜像 KC:428 的 find==end 分支）。
    ///
    /// 与并发的 [`KmerCatalog::get_kmer_count`] 之间只保证单键线性一致，
    /// 跨键不保证——镜像原版无锁竞态（另一线程可能同时读到置 0 前的旧值）。
    pub fn clear_kmer(&self, kmer: KmerId) {
        let key = self.canonical(kmer);
        if let Some(mut count) = self.counter.get_mut(&key) {
            *count = 0;
        }
    }
}

/// 只读视图：与 KmerCounter 同一查询语义（候选两方法用 trait 默认实现）。
/// 组装期不增不删键 → size() 恒等于转入时键数（含 count=0 键，惰性删除语义）。
impl KmerCatalog for SyncKmerCounter {
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
            .map(|count| *count)
            .unwrap_or(0)
    }

    fn get_kmer_string(&self, kmer: KmerId) -> Vec<u8> {
        decode_kmer_from_intval(kmer, self.kmer_length)
    }
}

impl ConcurrentClearCatalog for SyncKmerCounter {
    fn clear_kmer_concurrent(&self, kmer: KmerId) {
        self.clear_kmer(kmer)
    }
}

/// u128 并发 kmer 目录（k≤64 域; 刀M8 窄链 PARALLEL 组装期目录）。
/// 键 = KmerId 低 limb（k≤64 域内数值恒等, hi 恒 0）——DashMap 槽位
/// 40B→20B（QM4 ~2 亿键 ≈ 10.7→5.4GB）。查询侧 KmerId 参数经 canonical
/// 折叠后取 lo（与 KmerCounter.canonical 同式）, debug_assert 守域。
pub struct SyncKmerCounter128 {
    kmer_length: usize,
    ds_mode: bool,
    counter: DashMap<u128, u32, FxBuildHasher>,
}

impl SyncKmerCounter128 {
    /// u128 计数表（trusted 注入 + 剪枝完成态）整表转入。
    /// 转入瞬态 = 窄表 + DashMap 两份（宽链为宽 Fx 表 + 宽 DashMap 两份
    /// 大表再叠 KmerId 宽表物化——刀M8 前的全链峰值构成）。
    pub fn from_narrow(map: FxHashMap<u128, u32>, kmer_length: usize, ds_mode: bool) -> Self {
        let counter = DashMap::with_capacity_and_hasher(map.len(), FxBuildHasher);
        for (kmer, count) in map {
            counter.insert(kmer, count);
        }
        SyncKmerCounter128 {
            kmer_length,
            ds_mode,
            counter,
        }
    }

    /// canonical 折叠 + 收窄（k≤64 域内与 KmerCounter.canonical 数值恒等）。
    fn narrow_key(&self, kmer: KmerId) -> u128 {
        let c = if self.ds_mode {
            get_ds_kmer_val(kmer, self.kmer_length)
        } else {
            kmer
        };
        debug_assert_eq!(c.hi, 0, "k≤64 域外键进入窄目录");
        c.lo
    }

    /// count > 0 的键快照（PARALLEL 种子列表, 刀M6/8: 只取键 + 精确容量）。
    pub fn nonzero_keys(&self) -> Vec<u128> {
        let mut out = Vec::with_capacity(self.counter.len());
        self.counter.iter().for_each(|item| {
            if *item.value() > 0 {
                out.push(*item.key());
            }
        });
        out
    }
}

impl KmerCatalog for SyncKmerCounter128 {
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
            .get(&self.narrow_key(kmer))
            .map(|count| *count)
            .unwrap_or(0)
    }

    fn get_kmer_string(&self, kmer: KmerId) -> Vec<u8> {
        decode_kmer_from_intval(kmer, self.kmer_length)
    }
}

impl ConcurrentClearCatalog for SyncKmerCounter128 {
    fn clear_kmer_concurrent(&self, kmer: KmerId) {
        let key = self.narrow_key(kmer);
        if let Some(mut count) = self.counter.get_mut(&key) {
            *count = 0;
        }
    }
}

/// u128 计数表的独占借用适配器（刀M8）: 让泛型 [`crate::irke::prune_some_kmers`]
/// 直接作用于窄表, 免物化 KmerId 宽表。canonical/收窄语义与
/// SyncKmerCounter128 同式。
pub struct NarrowCounter128<'a> {
    map: &'a mut FxHashMap<u128, u32>,
    kmer_length: usize,
    ds_mode: bool,
}

impl<'a> NarrowCounter128<'a> {
    pub fn new(map: &'a mut FxHashMap<u128, u32>, kmer_length: usize, ds_mode: bool) -> Self {
        NarrowCounter128 {
            map,
            kmer_length,
            ds_mode,
        }
    }

    fn narrow_key(&self, kmer: KmerId) -> u128 {
        let c = if self.ds_mode {
            get_ds_kmer_val(kmer, self.kmer_length)
        } else {
            kmer
        };
        debug_assert_eq!(c.hi, 0, "k≤64 域外键进入窄适配器");
        c.lo
    }
}

impl KmerCatalog for NarrowCounter128<'_> {
    fn size(&self) -> usize {
        self.map.len()
    }

    fn get_kmer_length(&self) -> usize {
        self.kmer_length
    }

    fn is_double_stranded(&self) -> bool {
        self.ds_mode
    }

    fn get_kmer_count(&self, kmer: KmerId) -> u32 {
        self.map.get(&self.narrow_key(kmer)).copied().unwrap_or(0)
    }

    fn get_kmer_string(&self, kmer: KmerId) -> Vec<u8> {
        decode_kmer_from_intval(kmer, self.kmer_length)
    }
}

impl ExclusiveClearCatalog for NarrowCounter128<'_> {
    fn clear_kmer_exclusive(&mut self, kmer: KmerId) {
        let key = self.narrow_key(kmer);
        if let Some(count) = self.map.get_mut(&key) {
            *count = 0;
        }
    }

    fn snapshot_all_keys(&self) -> Vec<KmerId> {
        self.map
            .keys()
            .map(|&k| KmerId { hi: 0, lo: k })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer_counter::KmerCounter;
    use trinity_common::kmer::kmer_to_intval;

    fn enc(s: &[u8]) -> KmerId {
        kmer_to_intval(s).unwrap()
    }

    /// 与 KmerCounter 的查询逐项一致（SS 与 DS 两模式，含 DS canonical 互补链查询）。
    #[test]
    fn query_parity_with_kmer_counter() {
        for ds in [false, true] {
            let mut kc = KmerCounter::new(4, ds);
            kc.add_kmer(enc(b"ACGT"), 9); // DS: canonical 折叠
            kc.add_kmer(enc(b"ACGA"), 3);
            kc.add_kmer(enc(b"ACGG"), 5);
            kc.clear_kmer(enc(b"ACGG")); // 惰性删除: 键在、count=0

            let sync = SyncKmerCounter::from_counter(kc);
            assert_eq!(sync.size(), 3); // 含 0 值键
            assert_eq!(sync.get_kmer_length(), 4);
            assert_eq!(sync.is_double_stranded(), ds);
            assert_eq!(sync.get_kmer_count(enc(b"ACGT")), 9);
            assert_eq!(sync.get_kmer_count(enc(b"ACGG")), 0); // 已清
            assert_eq!(sync.get_kmer_count(enc(b"TTTT")), 0); // 缺键 → 0
            if ds {
                // 互补链同键: revcomp("ACGT") = "ACGT" 回文自映;
                // revcomp("ACGA") = "TCGT" 折到同键 → 查到 3
                assert_eq!(sync.get_kmer_count(enc(b"TCGT")), 3);
            }
            // 候选（trait 默认实现）: "ACGT" 前向候选 TGT?/TGA?/... 逐项对照
            // KmerCounter 的同键查询
            let mut kc2 = KmerCounter::new(4, ds);
            kc2.add_kmer(enc(b"ACGT"), 9);
            kc2.add_kmer(enc(b"ACGA"), 3);
            kc2.add_kmer(enc(b"ACGG"), 5);
            kc2.clear_kmer(enc(b"ACGG"));
            for probe in [enc(b"ACGT"), enc(b"ACGA"), enc(b"TTTT")] {
                assert_eq!(
                    sync.get_forward_kmer_candidates(probe),
                    kc2.get_forward_kmer_candidates(probe)
                );
                assert_eq!(
                    sync.get_reverse_kmer_candidates(probe),
                    kc2.get_reverse_kmer_candidates(probe)
                );
            }
            assert_eq!(
                sync.get_kmer_string(enc(b"ACGT")),
                b"ACGT".to_vec() // 解码不做 canonical
            );
        }
    }

    /// clear_kmer: 存在键原子置 0 且不删键；不存在键不插入（size 不变）。
    #[test]
    fn clear_kmer_zeroes_existing_and_never_inserts() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(KmerId::from(7), 5);
        let sync = SyncKmerCounter::from_counter(kc);
        sync.clear_kmer(KmerId::from(9)); // 缺键 → 无操作（不插入！）
        assert_eq!(sync.size(), 1);
        assert_eq!(sync.get_kmer_count(KmerId::from(9)), 0);
        sync.clear_kmer(KmerId::from(7));
        assert_eq!(sync.get_kmer_count(KmerId::from(7)), 0);
        assert_eq!(sync.size(), 1); // 惰性删除: 键仍在
    }

    /// DS 下从互补链方向 clear 经 canonical 折到同一键。
    #[test]
    fn clear_kmer_ds_canonicalizes() {
        let mut kc = KmerCounter::new(2, true);
        kc.add_kmer(enc(b"AC"), 4); // canonical 7
        let sync = SyncKmerCounter::from_counter(kc);
        sync.clear_kmer(enc(b"GT")); // revcomp("GT") = "AC"
        assert_eq!(sync.get_kmer_count(enc(b"AC")), 0);
        assert_eq!(sync.get_kmer_count(enc(b"GT")), 0);
    }

    /// 并发弱一致语义: 多线程同时 clear + get 无 panic、无撕裂；汇合后被清键
    /// 恒为 0。读到 0 或旧值皆合法（单键原子），断言只锁最终态。
    #[test]
    fn concurrent_clear_and_get_no_panic_final_state_zero() {
        let mut kc = KmerCounter::new(4, false);
        let keys: Vec<KmerId> = (0..200).map(|i| KmerId::from(1000 + i)).collect();
        for &k in &keys {
            kc.add_kmer(k, 7);
        }
        let sync = std::sync::Arc::new(SyncKmerCounter::from_counter(kc));

        std::thread::scope(|s| {
            for t in 0..4 {
                let sync = &sync;
                let keys = &keys;
                s.spawn(move || {
                    for (i, &k) in keys.iter().enumerate() {
                        if (i + t) % 2 == 0 {
                            sync.clear_kmer(k);
                        } else {
                            let _ = sync.get_kmer_count(k); // 0 或 7 皆合法
                        }
                    }
                });
            }
        });
        // 多轮把所有键清干净（各轮部分线程 clear 部分线程读）
        for _ in 0..3 {
            std::thread::scope(|s| {
                let sync = &sync;
                let keys = &keys;
                s.spawn(move || {
                    for &k in keys {
                        sync.clear_kmer(k);
                    }
                });
            });
        }
        for &k in &keys {
            assert_eq!(sync.get_kmer_count(k), 0);
        }
        assert_eq!(sync.size(), keys.len()); // 从未删键
    }

    /// 刀M8: 窄目录 SyncKmerCounter128 与宽目录 SyncKmerCounter 的
    /// 查询/候选/并发置零逐项一致（SS 与 DS 两模式, 含互补链探测）。
    #[test]
    fn sync128_query_parity_with_wide() {
        use crate::counter_sync::SyncKmerCounter128;
        for ds in [false, true] {
            let pairs = [
                (enc(b"ACGT"), 9u32),
                (enc(b"ACGA"), 3),
                (enc(b"ACGG"), 5),
                (enc(b"TTTT"), 7),
                (enc(b"ACGC"), 2),
            ];
            let mut kc = KmerCounter::new(4, ds);
            for &(k, c) in &pairs {
                kc.add_kmer(k, c); // DS: canonical 折叠 +=（同窄表构造口径）
            }
            kc.clear_kmer(enc(b"ACGG")); // 惰性置零
            let wide = SyncKmerCounter::from_counter(kc);

            let canon = |k: KmerId| if ds { get_ds_kmer_val(k, 4) } else { k };
            let mut narrow_map: FxHashMap<u128, u32> = FxHashMap::default();
            for &(k, c) in &pairs {
                *narrow_map.entry(canon(k).lo).or_insert(0) += c;
            }
            narrow_map.insert(canon(enc(b"ACGG")).lo, 0);
            let narrow = SyncKmerCounter128::from_narrow(narrow_map, 4, ds);

            assert_eq!(wide.size(), narrow.size());
            // 探测集: 表键 + 互补链 + 缺席键（含 1/2/3 碱基前缀邻域）
            let probes = [
                b"ACGT".as_ref(),
                b"ACGA".as_ref(),
                b"ACGG".as_ref(),
                b"TTTT".as_ref(),
                b"ACGC".as_ref(),
                b"TCGT".as_ref(), // ACGA 互补链（DS 折叠命中）
                b"GCGT".as_ref(), // ACGC 互补链
                b"AAAA".as_ref(),
                b"CCGT".as_ref(),
            ];
            for p in probes {
                let k = enc(p);
                assert_eq!(
                    wide.get_kmer_count(k),
                    narrow.get_kmer_count(k),
                    "ds={ds} probe={p:?}"
                );
                assert_eq!(
                    wide.get_forward_kmer_candidates(k),
                    narrow.get_forward_kmer_candidates(k),
                    "ds={ds} fwd {p:?}"
                );
                assert_eq!(
                    wide.get_reverse_kmer_candidates(k),
                    narrow.get_reverse_kmer_candidates(k),
                    "ds={ds} rev {p:?}"
                );
            }
            // 并发置零一致
            wide.clear_kmer_concurrent(enc(b"ACGA"));
            narrow.clear_kmer_concurrent(enc(b"ACGA"));
            for p in [b"ACGA".as_ref(), b"TCGT".as_ref()] {
                assert_eq!(wide.get_kmer_count(enc(p)), narrow.get_kmer_count(enc(p)), "ds={ds}");
                assert_eq!(narrow.get_kmer_count(enc(p)), 0, "置零后");
            }
            assert_eq!(wide.size(), narrow.size()); // 未删键
        }
    }

    /// 刀M8: 窄链 prune（NarrowCounter128 适配泛型 prune）与宽链
    /// prune_some_kmers(KmerCounter) 的非零计数多重集一致。
    #[test]
    fn prune_narrow_parity_with_wide() {
        use crate::counter_sync::NarrowCounter128;
        use crate::irke::{prune_some_kmers, prune_some_kmers_narrow};
        for ds in [false, true] {
            let seqs: &[&[u8]] = &[
                b"ACGTACGTAC", b"ACGTACGTAG", b"ACGTACGTAT", b"TTTTTTTTTT",
                b"ACGTACGTCG", b"GGGGGGGGGG", b"ACGACGACGT",
            ];
            let mut kc = KmerCounter::new(5, ds);
            let mut narrow_map: FxHashMap<u128, u32> = FxHashMap::default();
            let canon = |k: KmerId| if ds { get_ds_kmer_val(k, 5) } else { k };
            for (i, s) in seqs.iter().enumerate() {
                for j in 0..=s.len() - 5 {
                    let k = kmer_to_intval(&s[j..j + 5]).unwrap();
                    kc.add_kmer(k, 1);
                    let c = canon(k);
                    *narrow_map.entry(c.lo).or_insert(0) += 1;
                }
                let _ = i;
            }
            let mut narrow = narrow_map;
            let pruned_wide = prune_some_kmers(&mut kc, 2, 1.0, true, 0.005);
            let pruned_narrow = prune_some_kmers_narrow(&mut narrow, 5, ds, 2, 1.0, true, 0.005);
            assert_eq!(pruned_wide, pruned_narrow, "ds={ds}");
            // 非零多重集一致（宽侧收集 / 窄侧收窄比较）
            let mut w: Vec<(KmerId, u32)> = kc
                .iter_nonzero()
                .map(|(k, c)| (k, c))
                .collect();
            w.sort();
            let mut n: Vec<(KmerId, u32)> = narrow
                .iter()
                .filter(|&(_, &c)| c > 0)
                .map(|(&k, &c)| (KmerId { hi: 0, lo: k }, c))
                .collect();
            n.sort();
            assert_eq!(w, n, "ds={ds}");
        }
    }
}
