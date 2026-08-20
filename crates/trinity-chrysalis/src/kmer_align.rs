//! KmerAlignCore 移植 —— 2×12-mer 倒排索引求交得到 24-mer 精确匹配
//! （P3-T1 基础层）。镜像 `Chrysalis/aligns/{KmerAlignCore.h,KmerAlignCore.cc}`。
//!
//! 原版 `m_table` 是 4^12 ≈ 1.68e7 个 `KmerAlignCoreRecordStore`（svec 对象）
//! 的数组；本版用 **CSR**（`offsets` 桶前缀和 + `records` 紧凑数组）等价替换：
//! - `AddData(vecDNAVector)` 的"计数遍 → 定容 → 填充遍"两遍式布局完全保留；
//! - 桶内记录序 = 填充序 = (contig 升序, pos 升序)，另做显式桶内排序兜底；
//! - `GetMatches`（m_numTables=2、m_lookAhead=0、m_max12=∞ 的默认配置）：
//!   对查询序列前/后 12-mer 各查一次表，后者位置归一 `pos - 12`，
//!   再以两指针线性交（MergeSortFilter）取公共 (contig, pos)。
//!
//! 查询 12-mer 含 N 或非法字符（小写等）→ `bases_to_number` = -1 → 该表 hits
//! 为空 → 交集自然为空（原版同型）。

use crate::dna_vector::{bases_to_number, DnaSeq};

/// TranslateBasesToNumber 默认 m_range = 12（KmerAlignCore.h:16）
pub const KMER_SIZE: usize = 12;

/// GetWordSize = m_numTables(2) × 12 = 24（KmerAlignCore.h:293）
pub const WORD_SIZE: usize = 2 * KMER_SIZE;

/// GetBoundValue = 4^12 桶数（KmerAlignCore.h:26-34）
const NUM_BUCKETS: usize = 1 << (2 * KMER_SIZE);

/// KmerAlignCoreRecord（KmerAlignCore.h:156-204）：比较只看 (contig, pos)。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MatchRecord {
    pub contig: i32,
    pub pos: i32,
}

impl MatchRecord {
    #[inline]
    fn key(&self) -> (i32, i32) {
        (self.contig, self.pos)
    }
}

/// CSR 形态的 KmerAlignCore。
pub struct KmerAlignCore {
    /// 桶前缀和，长度 NUM_BUCKETS+1；桶 b = records[offsets[b]..offsets[b+1]]
    offsets: Vec<i32>,
    /// 紧凑记录（桶内 (contig,pos) 升序）
    records: Vec<MatchRecord>,
}

impl KmerAlignCore {
    /// 镜像 `AddData(const vecDNAVector &)`（KmerAlignCore.cc:18-24 → 三参重载
    /// min=1 且 tags 为空 → IsRepeat 恒 false → 所有窗口入表）。
    pub fn build(seqs: &[DnaSeq]) -> Self {
        // ---- 第一遍：计数（cc:54-86）----
        let mut counts = vec![0i32; NUM_BUCKETS];
        for d in seqs {
            let len = d.seq.len();
            if len < KMER_SIZE {
                continue;
            }
            for k in 0..=len - KMER_SIZE {
                let n = bases_to_number(&d.seq, k, KMER_SIZE);
                if n >= 0 {
                    counts[n as usize] += 1;
                }
            }
        }

        // ---- 前缀和 → CSR offsets（cc:91-97 的 Resize 定容）----
        let mut offsets = Vec::with_capacity(NUM_BUCKETS + 1);
        let mut running: i32 = 0;
        offsets.push(0);
        for &c in &counts {
            running += c;
            offsets.push(running);
        }

        // ---- 第二遍：填充（cc:100-120）；counts 复用为桶内游标 ----
        let mut records = vec![MatchRecord { contig: 0, pos: 0 }; running as usize];
        for c in counts.iter_mut() {
            *c = 0;
        }
        for (j, d) in seqs.iter().enumerate() {
            let len = d.seq.len();
            if len < KMER_SIZE {
                continue;
            }
            for k in 0..=len - KMER_SIZE {
                let n = bases_to_number(&d.seq, k, KMER_SIZE);
                if n >= 0 {
                    let b = n as usize;
                    let idx = (offsets[b] + counts[b]) as usize;
                    records[idx] = MatchRecord {
                        contig: j as i32,
                        pos: k as i32,
                    };
                    counts[b] += 1;
                }
            }
        }
        drop(counts);

        // ---- 桶内 (contig,pos) 升序（填充序已天然升序，显式排序兜底）----
        for b in 0..NUM_BUCKETS {
            let s = offsets[b] as usize;
            let e = offsets[b + 1] as usize;
            if e - s > 1 {
                records[s..e].sort_unstable_by_key(|r| r.key());
            }
        }

        KmerAlignCore { offsets, records }
    }

    #[inline]
    fn bucket(&self, n: i64) -> &[MatchRecord] {
        let b = n as usize;
        let s = self.offsets[b] as usize;
        let e = self.offsets[b + 1] as usize;
        &self.records[s..e]
    }

    /// 镜像 `GetMatches(matches, b, start=0)` 默认配置（numTables=2、
    /// lookAhead=0、max12=∞）。查询取 seq 的**前 24bp**；len < 24 → 空
    /// （原版 cerr "Error: sequence length=..." 后返回 false）。
    ///
    /// 返回值 = 与查询 24-mer 精确相等的位置集合（(contig, pos) 升序），
    /// 含查询来源 contig 的自匹配（原版无自排除）。
    pub fn get_matches(&self, seq: &[u8]) -> Vec<MatchRecord> {
        if seq.len() < WORD_SIZE {
            return Vec::new();
        }
        let n0 = bases_to_number(seq, 0, KMER_SIZE);
        let n1 = bases_to_number(seq, KMER_SIZE, KMER_SIZE);

        let one: Vec<MatchRecord> = if n0 >= 0 {
            self.bucket(n0).to_vec()
        } else {
            Vec::new()
        };
        let two: Vec<MatchRecord> = if n1 >= 0 {
            // i=1 表的位置归一：pos - i*size（cc:221）
            self.bucket(n1)
                .iter()
                .map(|&r| MatchRecord {
                    contig: r.contig,
                    pos: r.pos - KMER_SIZE as i32,
                })
                .collect()
        } else {
            Vec::new()
        };
        merge_sort_filter(&one, &two)
    }
}

/// MergeSortFilter（KmerAlignCore.cc:287-340）：两指针线性交，保 one 的顺序。
/// 任一侧为空 → 空结果。
fn merge_sort_filter(one: &[MatchRecord], two: &[MatchRecord]) -> Vec<MatchRecord> {
    let mut result = Vec::new();
    if one.is_empty() || two.is_empty() {
        return result;
    }
    let mut y = 0usize;
    for x in one {
        while y < two.len() && two[y].key() < x.key() {
            y += 1;
        }
        if y >= two.len() {
            break;
        }
        if two[y].key() == x.key() {
            result.push(*x);
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dna(name: &str, seq: &[u8]) -> DnaSeq {
        DnaSeq {
            name: name.to_string(),
            seq: seq.to_vec(),
        }
    }

    /// 周期 8 序列 X=ACGTTGCA，S = X×3（24bp，首/尾 12-mer 不同）。
    /// contig0 = S+ACG（S 出现于 0）；contig1 = TTTTT+S+G（S 出现于 5）；
    /// contig2 = S[..12]+N+S[12..]+AAAA（S 被 N 阻断——q0 命中 (2,0)、
    /// q1 命中 (2,13)→归一 (2,1)，两者不等 → 不产匹配）；
    /// contig3 = GG+q1+q1（q1 = S[12..24)，构造 q1q1 24-mer 的匹配）。
    fn fixture() -> (Vec<u8>, Vec<u8>, KmerAlignCore) {
        let x = b"ACGTTGCA";
        let mut s = x.to_vec();
        s.extend_from_slice(x);
        s.extend_from_slice(x);
        let q1 = s[12..].to_vec();

        let contigs = vec![
            dna("c0", &[s.clone(), b"ACG".to_vec()].concat()),
            dna(
                "c1",
                &[b"TTTTT".to_vec(), s.clone(), b"G".to_vec()].concat(),
            ),
            dna(
                "c2",
                &[
                    s[..12].to_vec(),
                    b"N".to_vec(),
                    s[12..].to_vec(),
                    b"AAAA".to_vec(),
                ]
                .concat(),
            ),
            dna("c3", &[b"GG".to_vec(), q1.clone(), q1.clone()].concat()),
        ];
        let q1q1 = [q1.clone(), q1.clone()].concat();
        (s, q1q1, KmerAlignCore::build(&contigs))
    }

    /// 手推桶内容：S 的首 12-mer "ACGTTGCAACGT" 命中 (0,0),(0,8),(1,5),(1,13),(2,0)；
    /// 交集后只剩真正的 24-mer 匹配 (0,0)（自匹配包含）与 (1,5)（跨 contig 共享）。
    #[test]
    fn get_matches_cross_contig_and_self_match() {
        let (s, _, core) = fixture();
        let m = core.get_matches(&s);
        assert_eq!(m.len(), 2);
        assert_eq!((m[0].contig, m[0].pos), (0, 0)); // 自匹配无排除
        assert_eq!((m[1].contig, m[1].pos), (1, 5)); // contig1 中偏移 5
    }

    /// N 阻断的 contig2：q0 与 q1 各自命中但归一位置差 1（(2,0) vs (2,1)）
    /// → 交集空（上文 fixture() 推导），验证位置归一参与判等。
    #[test]
    fn get_matches_n_broken_contig_excluded() {
        let (s, _, core) = fixture();
        assert!(
            !core.get_matches(&s).iter().any(|r| r.contig == 2),
            "被 N 阻断的 contig2 不应产生 24-mer 匹配"
        );
    }

    /// 位置归一直接验证：查询 q1q1（首/尾 12-mer 相同），q1 在 contig3 出现于
    /// 2 与 14；仅 (3,2) 处两半同时成立 → 唯一匹配 (3,2)，(3,14) 被归一过滤。
    #[test]
    fn get_matches_position_normalization() {
        let (_, q1q1, core) = fixture();
        let m = core.get_matches(&q1q1);
        assert_eq!(m.len(), 1);
        assert_eq!((m[0].contig, m[0].pos), (3, 2));
    }

    /// 查询含 N → 首 12-mer 编码 -1 → 该表 hits 空 → 交集空。
    #[test]
    fn get_matches_query_with_n_is_empty() {
        let (s, _, core) = fixture();
        let mut q = vec![b'N'];
        q.extend_from_slice(&s[..23]);
        assert!(core.get_matches(&q).is_empty());
        // 后半含 N 同理
        let mut q2 = s[..12].to_vec();
        q2.push(b'N');
        q2.extend_from_slice(&s[12..23]);
        assert!(core.get_matches(&q2).is_empty());
    }

    /// len < 24 → 空（原版 Error 分支）。
    #[test]
    fn get_matches_too_short_is_empty() {
        let (s, _, core) = fixture();
        assert!(core.get_matches(&s[..20]).is_empty());
    }

    /// 无共享 24-mer 的查询 → 空；交集保序（按 (contig,pos) 升序）。
    #[test]
    fn get_matches_no_hit_and_ordering() {
        let (_, _, core) = fixture();
        let q = b"GGGGGGGGGGGGCCCCCCCCCCCC";
        assert!(core.get_matches(q).is_empty());
        // 保序：S 的匹配本就以 (0,0) < (1,5) 输出（见第一测试）；
        // 再验多个命中时输出与 one 侧顺序一致——构造不了多命中场景（表已定），
        // 由 merge_sort_filter 的两指针保序实现保证。
    }

    #[test]
    fn merge_filter_two_pointer_semantics() {
        use MatchRecord as R;
        let one = [
            R { contig: 0, pos: 0 },
            R { contig: 0, pos: 8 },
            R { contig: 1, pos: 5 },
        ];
        let two = [
            R {
                contig: 0,
                pos: -12,
            },
            R { contig: 0, pos: 0 },
            R { contig: 1, pos: 5 },
        ];
        let m = merge_sort_filter(&one, &two);
        assert_eq!(m.len(), 2);
        assert_eq!(m[0].key(), (0, 0));
        assert_eq!(m[1].key(), (1, 5));
        // 空 → 空
        assert!(merge_sort_filter(&one, &[]).is_empty());
        assert!(merge_sort_filter(&[], &two).is_empty());
    }

    /// 空 build（无 contig）可用；任意查询 → 空。
    #[test]
    fn build_empty() {
        let core = KmerAlignCore::build(&[]);
        assert!(core.get_matches(b"ACGTTGCAACGTTGCAACGTTGCA").is_empty());
    }

    /// 恰 24bp 的 contig 只产 1 个窗口；查询其自身 → 自匹配 (0,0)。
    #[test]
    fn build_exactly_24bp_contig() {
        let core = KmerAlignCore::build(&[dna("c", b"ACGTTGCAACGTTGCAACGTTGCA")]);
        let m = core.get_matches(b"ACGTTGCAACGTTGCAACGTTGCA");
        assert_eq!(m.len(), 1);
        assert_eq!((m[0].contig, m[0].pos), (0, 0));
    }
}
