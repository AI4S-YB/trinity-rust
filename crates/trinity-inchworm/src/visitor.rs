//! Kmer_visitor — 直译 KmerCounter.cpp:756-799。
//!
//! 原版 std::set<kmer_int_type_t> → std::HashSet<u64>：集合语义一致；
//! 4^32 键空间下两容器都不发生行为级碰撞差异（计划 Task 1 规格选择）。
//! DS 模式下 add/exists/erase 全部先 canonical（KC:766-789）。

use std::collections::HashSet;

use trinity_common::kmer::{get_ds_kmer_val, KmerId};

pub struct KmerVisitor {
    kmer_length: usize,
    ds_mode: bool,
    visited: HashSet<KmerId>,
}

impl KmerVisitor {
    /// KC:760-764
    pub fn new(kmer_length: usize, ds_mode: bool) -> Self {
        KmerVisitor {
            kmer_length,
            ds_mode,
            visited: HashSet::new(),
        }
    }

    /// DS canonical 后的键（KC:768-770）
    fn canonical(&self, kmer: KmerId) -> KmerId {
        if self.ds_mode {
            get_ds_kmer_val(kmer, self.kmer_length)
        } else {
            kmer
        }
    }

    /// KC:766-772: DS canonical 后插入
    pub fn add(&mut self, kmer: KmerId) {
        self.visited.insert(self.canonical(kmer));
    }

    /// KC:774-780: DS canonical 后查询
    pub fn exists(&self, kmer: KmerId) -> bool {
        self.visited.contains(&self.canonical(kmer))
    }

    /// KC:782-789: DS canonical 后移除（不存在则无操作）
    pub fn erase(&mut self, kmer: KmerId) {
        self.visited.remove(&self.canonical(kmer));
    }

    /// KC:791-794
    pub fn clear(&mut self) {
        self.visited.clear();
    }

    /// KC:796-799
    pub fn size(&self) -> usize {
        self.visited.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use trinity_common::kmer::kmer_to_intval;

    fn enc(s: &[u8]) -> KmerId {
        kmer_to_intval(s).unwrap()
    }

    #[test]
    fn ds_visitor_canonicalizes() {
        // "GT"=2 → canonical 7（"AC"）: GT 与 AC 互斥共存于同一键
        let mut v = KmerVisitor::new(2, true);
        v.add(enc(b"GT"));
        assert!(v.exists(enc(b"GT")));
        assert!(v.exists(enc(b"AC")));
        assert_eq!(v.size(), 1);
        assert!(!v.exists(enc(b"TT"))); // canonical("TT") = "TT"=10，未加入
                                        // 经任一互补链 erase，两条链同时失效
        v.erase(enc(b"AC"));
        assert!(!v.exists(enc(b"GT")));
        assert!(!v.exists(enc(b"AC")));
        assert_eq!(v.size(), 0);
    }

    #[test]
    fn ss_visitor_raw_keys() {
        // 非 DS: 原始键，互补链互不干扰
        let mut v = KmerVisitor::new(2, false);
        v.add(enc(b"GT"));
        assert!(v.exists(enc(b"GT")));
        assert!(!v.exists(enc(b"AC")));
        v.clear();
        assert_eq!(v.size(), 0);
        assert!(!v.exists(enc(b"GT")));
    }

    #[test]
    fn visitor_backtrack_semantics() {
        // inchworm_step 递归回溯依赖: erase 后同键可再次 add（KC:782-789 erase + 重入）
        let mut v = KmerVisitor::new(2, false);
        v.add(7);
        assert!(v.exists(7));
        v.erase(7);
        assert!(!v.exists(7));
        v.add(7);
        assert_eq!(v.size(), 1);
        // erase 不存在的键无副作用
        v.erase(9);
        assert_eq!(v.size(), 1);
    }
}
