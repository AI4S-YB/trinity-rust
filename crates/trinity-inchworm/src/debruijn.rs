//! `Inchworm/src/FastaToDeBruijn.cpp` + `Inchworm/src/DeBruijnGraph.cpp` 的移植。
//!
//! 物理位置放 trinity-inchworm（原版在 Inchworm 二进制里），供 Chrysalis
//! GraphFromFasta 流程使用：
//! - `DeBruijnGraph`：kmer 节点表 + 4-bit prev/next 掩码邻接（DS 模式按
//!   canonical 键存储，方向翻转时前驱/后继角色互换）。
//! - `to_chrysalis_format`：Chrysalis `Component` 图格式（优先队列按
//!   kmer_count 降序遍历，'-' 定向节点 id 加 N 偏移）。
//! - `graph_per_record`：捆绑 FASTA（`>s_<n> <cov>...` + `X` 连接序列）逐
//!   记录建图输出。

use std::collections::BinaryHeap;

use trinity_common::error::CommonError;
use trinity_common::fasta::FastaRecord;
use trinity_common::kmer::{
    decode_kmer_from_intval, get_ds_kmer_val, kmer_to_intval, revcomp_val, KmerId,
};

/// DeBruijnGraph.cpp:11-14 — 4-bit 邻接掩码。
/// 编码序（base_to_int：G=0,A=1,T=2,C=3）对应位 8/4/2/1。
const G_MASK: u8 = 8;
const A_MASK: u8 = 4;
const T_MASK: u8 = 2;
const C_MASK: u8 = 1;

fn base_mask(code: u64) -> u8 {
    match code {
        0 => G_MASK, // G
        1 => A_MASK, // A
        2 => T_MASK, // T
        3 => C_MASK, // C
        _ => 0,
    }
}

/// DeBruijnGraph.hpp DeBruijnKmer — 节点 = kmer 值 + 1-based id + 覆盖计数
/// + prev/next 掩码。annotations 字段（仅调试用）不移植。
#[derive(Debug, Clone)]
pub struct DeBruijnKmer {
    kmer: KmerId,
    /// 1-based，由首次插入序决定（++_kmer_id_counter）
    id: u64,
    /// increment_kmer_count：每次出现 += 整条 contig 覆盖度（cov_val，非 +1）
    kmer_count: i64,
    prev: u8,
    next: u8,
}

impl DeBruijnKmer {
    fn new(kmer: KmerId, id: u64) -> Self {
        DeBruijnKmer {
            kmer,
            id,
            kmer_count: 0,
            prev: 0,
            next: 0,
        }
    }

    pub fn id(&self) -> u64 {
        self.id
    }

    pub fn kmer_count(&self) -> i64 {
        self.kmer_count
    }

    pub fn prev_mask(&self) -> u8 {
        self.prev
    }

    pub fn next_mask(&self) -> u8 {
        self.next
    }

    /// DeBruijnGraph.cpp:176 add_prev_kmer — 取 k 的**首碱基**置位。
    fn add_prev_kmer(&mut self, k: KmerId, kmer_length: usize) {
        let first = k >> (kmer_length * 2 - 2);
        self.prev |= base_mask(first);
    }

    /// DeBruijnGraph.cpp:214 add_next_kmer — 取 k 的**末碱基**置位。
    fn add_next_kmer(&mut self, k: KmerId) {
        let last = k & 3;
        self.next |= base_mask(last);
    }

    /// DeBruijnGraph.cpp:105 get_prev_kmers — 候选 = 首碱基拼上 k 的后缀
    /// （k >> 2）；序 G,A,T,C（只产掩码位非零者）。
    pub fn get_prev_kmers(&self, kmer_length: usize) -> Vec<KmerId> {
        let reverse_suffix = self.kmer >> 2;
        let mut prev_kmers = Vec::new();
        for (mask, code) in [(G_MASK, 0u64), (A_MASK, 1), (T_MASK, 2), (C_MASK, 3)] {
            if self.prev & mask != 0 {
                prev_kmers.push((code << (kmer_length * 2 - 2)) | reverse_suffix);
            }
        }
        prev_kmers
    }

    /// DeBruijnGraph.cpp:145 get_next_kmers — 候选 = k 左移一位碱基（截去
    /// 首碱基）拼上末碱基；序 G,A,T,C。
    pub fn get_next_kmers(&self, kmer_length: usize) -> Vec<KmerId> {
        // 原版 (k << (33-K)*2) >> (32-K)*2 == (k & 低 2K-2 位) << 2
        // （截掉首碱基的高 2 位后整体左移一个碱基）
        let forward_prefix = (self.kmer & ((1u64 << (kmer_length * 2 - 2)) - 1)) << 2;
        let mut next_kmers = Vec::new();
        for (mask, code) in [(G_MASK, 0u64), (A_MASK, 1), (T_MASK, 2), (C_MASK, 3)] {
            if self.next & mask != 0 {
                next_kmers.push(forward_prefix | code);
            }
        }
        next_kmers
    }
}

/// 优先队列元素（DeBruijnGraph.cpp:545 kmer_count_comparer 的 max-heap）。
/// 平局弹出序在原版由 std::priority_queue 内部堆序决定、无语义保证——
/// 这里用（count 降序，kmer 值升序）做确定性 tie-break；由此产生的**输出
/// 行序非契约**（对拍按行多重集比较）。
#[derive(Debug, Clone, PartialEq, Eq)]
struct QItem {
    kmer: KmerId,
    count: i64,
}

impl Ord for QItem {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.count
            .cmp(&other.count)
            .then_with(|| other.kmer.cmp(&self.kmer)) // kmer 升序 tie-break
    }
}

impl PartialOrd for QItem {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// DeBruijnGraph.hpp DeBruijnGraph — std::map → BTreeMap（kmer 升序遍历，
/// 与原版一致，影响 root 拾取与"全空取首键"的确定性）。
#[derive(Debug)]
pub struct DeBruijnGraph {
    kmer_length: usize,
    kmer_map: std::collections::BTreeMap<KmerId, DeBruijnKmer>,
    kmer_id_counter: u64,
}

impl DeBruijnGraph {
    pub fn new(kmer_length: usize) -> Self {
        DeBruijnGraph {
            kmer_length,
            kmer_map: std::collections::BTreeMap::new(),
            kmer_id_counter: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.kmer_map.len()
    }

    pub fn is_empty(&self) -> bool {
        self.kmer_map.is_empty()
    }

    pub fn get_node(&self, kmer: KmerId) -> Option<&DeBruijnKmer> {
        self.kmer_map.get(&kmer)
    }

    /// DeBruijnGraph.cpp:271 add_sequence — 逐 kmer 入图。
    ///
    /// - DS：键 = canonical(k) = max(k, rc(k))；k 非 canonical 时该 kmer 视为
    ///   '-' 定向，其前驱/后继经 revcomp 后**角色互换**（前驱变后继边）。
    /// - kmer_count += cov_val（整条 contig 的覆盖度，非 +1）。
    pub fn add_sequence(
        &mut self,
        _accession: &str,
        sequence: &[u8],
        s_strand: bool,
        cov_val: u32,
    ) -> Result<(), CommonError> {
        // sequence_string_to_kmer_int_type_vector（sequenceUtil.cpp:387）。
        // 原版 `i <= seq.len() - k` 用无符号下溢——调用方（graph_per_record）
        // 已跳过 len < k 的区域；这里显式守卫为空向量。
        if sequence.len() < self.kmer_length {
            return Ok(());
        }
        let kit_vec: Vec<KmerId> = sequence
            .windows(self.kmer_length)
            .map(kmer_to_intval)
            .collect::<Result<_, _>>()?;

        let k_len = self.kmer_length;
        for i in 0..kit_vec.len() {
            let orig_sequence_k = kit_vec[i];
            let mut k = orig_sequence_k;

            if !s_strand {
                k = get_ds_kmer_val(k, k_len);
            }

            // get_kmer_node（DeBruijnGraph.cpp:341）：miss 时新建，id = ++counter
            let next_id = self.kmer_id_counter + 1;
            let dk = self
                .kmer_map
                .entry(k)
                .or_insert_with(|| DeBruijnKmer::new(k, next_id));
            if dk.id == next_id {
                self.kmer_id_counter = next_id;
            }
            dk.kmer_count += cov_val as i64;

            let flipped = k != orig_sequence_k; // kmer_orient == "-"

            if i != 0 {
                // SS 原样使用；DS 未翻转同向；DS 翻转 → 前驱 revcomp 后变后继边
                let mut pk = kit_vec[i - 1];
                if !s_strand && flipped {
                    pk = revcomp_val(pk, k_len);
                    self.kmer_map.get_mut(&k).unwrap().add_next_kmer(pk);
                } else {
                    self.kmer_map.get_mut(&k).unwrap().add_prev_kmer(pk, k_len);
                }
            }
            if i != kit_vec.len() - 1 {
                // 镜像：DS 翻转 → 后继 revcomp 后变前驱边
                let mut nk = kit_vec[i + 1];
                if !s_strand && flipped {
                    nk = revcomp_val(nk, k_len);
                    self.kmer_map.get_mut(&k).unwrap().add_prev_kmer(nk, k_len);
                } else {
                    self.kmer_map.get_mut(&k).unwrap().add_next_kmer(nk);
                }
            }
        }
        Ok(())
    }

    /// DeBruijnGraph.cpp:639 get_root_kmers — 起点 = 无 prev 者；DS 下无
    /// next 者也算（可反向走）。全空（环）时由调用方取 map 首键。
    fn get_root_kmers(&self, s_strand: bool) -> Vec<DeBruijnKmer> {
        self.kmer_map
            .values()
            .filter(|dk| {
                dk.get_prev_kmers(self.kmer_length).is_empty()
                    || (!s_strand && dk.get_next_kmers(self.kmer_length).is_empty())
            })
            .cloned()
            .collect()
    }

    /// DeBruijnGraph.cpp:556 toChrysalisFormat — Chrysalis 图格式输出。
    ///
    /// 输出行序依赖优先队列平局弹出序（非契约）；id 语义：存储节点 id 为
    /// 首次插入序，'-' 定向输出 id+N（N = 图节点总数）。
    pub fn to_chrysalis_format(&self, component_id: i64, s_strand: bool) -> String {
        let k_len = self.kmer_length;
        let mut s = format!("Component {component_id}\n");

        let mut seen: std::collections::HashMap<KmerId, bool> = std::collections::HashMap::new();
        let mut reported: std::collections::HashMap<KmerId, bool> =
            std::collections::HashMap::new();

        let mut root_kmers = self.get_root_kmers(s_strand);
        if root_kmers.is_empty() {
            // circular, initiate from the first kmer
            if let Some(first) = self.kmer_map.values().next() {
                root_kmers.push(first.clone());
            }
        }

        let total_kmer_count = self.kmer_map.len() as u64;

        // kmer_sorter_by_count_desc（std::sort 不稳定，平局序无保证）——
        // 用 count 降序 + kmer 升序 tie-break 保证确定性。
        root_kmers.sort_by(|a, b| {
            b.kmer_count
                .cmp(&a.kmer_count)
                .then_with(|| a.kmer.cmp(&b.kmer))
        });

        let mut collected_kmers: Vec<KmerId> = Vec::new();

        let mut kmer_queue: BinaryHeap<QItem> = BinaryHeap::new();

        for root in &root_kmers {
            // root 存储键恒视为 '+' 起向
            kmer_queue.push(QItem {
                kmer: root.kmer,
                count: root.kmer_count,
            });

            while let Some(pk_item) = kmer_queue.pop() {
                let k = pk_item.kmer; // 已按当前定向

                let kmer_seq = decode_kmer_from_intval(k, k_len);

                // 存储键：SS 用 k / DS 用 canonical(k)
                let dk_stored_kmer_val = if s_strand {
                    k
                } else {
                    get_ds_kmer_val(k, k_len)
                };
                let dk_orient = if dk_stored_kmer_val == k { '+' } else { '-' };

                let dk = match self.kmer_map.get(&dk_stored_kmer_val) {
                    Some(dk) => dk,
                    None => continue, // 原版 find()->second 为 UB；不发生
                };
                let mut dk_id = dk.id;
                if dk_orient == '-' {
                    dk_id += total_kmer_count;
                }

                if seen.get(&k).copied().unwrap_or(false) {
                    // already seen it
                    if !s_strand && !reported.get(&k).copied().unwrap_or(false) {
                        s += &format!(
                            "{dk_id}\t-1\t1\t{}\t1\n",
                            String::from_utf8_lossy(&kmer_seq)
                        );
                        reported.insert(k, true);
                    }
                    continue;
                }
                seen.insert(k, true);

                if !s_strand {
                    collected_kmers.push(k);
                }
                // palindromic kmer 警告（原版仅 cerr）不移植

                // 左行 prev：'+' 用 prev 掩码；'-' 用 next 候选的 revcomp
                let prev_kmers: Vec<KmerId> = if dk_orient == '+' {
                    dk.get_prev_kmers(k_len)
                } else {
                    dk.get_next_kmers(k_len)
                        .into_iter()
                        .map(|nk| revcomp_val(nk, k_len))
                        .collect()
                };

                let mut reached_terminal_extension = false;
                let kmer_seq_str = String::from_utf8_lossy(&kmer_seq).into_owned();

                if !prev_kmers.is_empty() {
                    for pk_kmer in prev_kmers {
                        let pk_stored_kmer = if s_strand {
                            pk_kmer
                        } else {
                            get_ds_kmer_val(pk_kmer, k_len)
                        };
                        let pk_orient = if pk_kmer == pk_stored_kmer { '+' } else { '-' };
                        let pkd = match self.kmer_map.get(&pk_stored_kmer) {
                            Some(pkd) => pkd,
                            None => continue,
                        };
                        let mut pkd_id = pkd.id;
                        if pk_orient == '-' {
                            pkd_id += total_kmer_count;
                        }
                        let pk_kmer_count = pkd.kmer_count;

                        s += &format!("{dk_id}\t{pkd_id}\t1\t{kmer_seq_str}\t1\n");

                        kmer_queue.push(QItem {
                            kmer: pk_kmer,
                            count: pk_kmer_count,
                        });
                    }
                } else {
                    // hit left end, no prev extension
                    s += &format!("{dk_id}\t-1\t1\t{kmer_seq_str}\t1\n");
                    reached_terminal_extension = true;
                }
                reported.insert(k, true);

                // 右行 next：只入队不打印；'-' 用 prev 候选的 revcomp
                let next_kmers: Vec<KmerId> = if dk_orient == '+' {
                    dk.get_next_kmers(k_len)
                } else {
                    dk.get_prev_kmers(k_len)
                        .into_iter()
                        .map(|pk| revcomp_val(pk, k_len))
                        .collect()
                };

                if next_kmers.is_empty() {
                    reached_terminal_extension = true;
                } else {
                    for nk_kmer in next_kmers {
                        let nk_stored_kmer = if s_strand {
                            nk_kmer
                        } else {
                            get_ds_kmer_val(nk_kmer, k_len)
                        };
                        let nk_count = match self.kmer_map.get(&nk_stored_kmer) {
                            Some(nkd) => nkd.kmer_count,
                            None => continue,
                        };
                        kmer_queue.push(QItem {
                            kmer: nk_kmer,
                            count: nk_count,
                        });
                    }
                }

                if !s_strand && reached_terminal_extension {
                    // 终端清理：collected 的 revcomp 标 seen
                    for kmer in collected_kmers.drain(..) {
                        seen.insert(revcomp_val(kmer, k_len), true);
                    }
                }
            }
        }

        s
    }
}

/// string_util.cpp tokenize — 跳过前导分隔符；连续分隔符不产空 token。
fn tokenize(s: &str, delim: char) -> Vec<&str> {
    s.split(delim).filter(|t| !t.is_empty()).collect()
}

/// atoi 等价：解析失败按 0。
fn atoi(s: &str) -> i64 {
    // atoi 停在首个非数字处（"12x" → 12）；取前导数字（含符号）部分
    let t = s.trim_start();
    let bytes = t.as_bytes();
    let mut end = 0;
    if end < bytes.len() && (bytes[end] == b'-' || bytes[end] == b'+') {
        end += 1;
    }
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        end += 1;
    }
    if end == 0 {
        return 0;
    }
    t[..end].parse().unwrap_or(0)
}

/// sequenceUtil.cpp contains_non_gatc / replace_nonGATC_chars_with_A —
/// 非 gatc（含小写 gatc 之外的）字符替换为 'A'。
fn replace_non_gatc_with_a(seq: &[u8]) -> Vec<u8> {
    use trinity_common::kmer::base_to_int;
    seq.iter()
        .map(|&c| if base_to_int(c).is_some() { c } else { b'A' })
        .collect()
}

/// FastaToDeBruijn.cpp:152 createGraphPerRecord — 逐捆绑记录建图并输出
/// Chrysalis 图格式（`--graph_per_record` 默认路径）。
///
/// 原版 omp 并行 + critical 输出（块序不定）；这里顺序处理、按输入序输出
/// ——块序非契约，对拍按每 Component 块的行多重集比较。
pub fn graph_per_record(
    bundles: &[FastaRecord],
    kmer_length: usize,
    s_strand: bool,
) -> Result<String, CommonError> {
    let mut out = String::new();

    for fe in bundles {
        let accession = fe.accession.as_str();

        // component value: tokenize(accession, "_")[1] → atoi
        let acc_pts = tokenize(accession, '_');
        let component_id = atoi(acc_pts.get(1).copied().unwrap_or(""));

        // iworm coverage values: tokenize(header, " ") 去首 token（accession）
        let mut iworm_cov_vals_str_vec = tokenize(&fe.header, ' ');
        if !iworm_cov_vals_str_vec.is_empty() {
            iworm_cov_vals_str_vec.remove(0);
        }

        // inchworm bundles concatenated with 'X' delimiters by Chrysalis
        let seq_regions = tokenize(&fe.sequence, 'X');

        if seq_regions.len() != iworm_cov_vals_str_vec.len() {
            return Err(CommonError::Parse(
                "Error, number of seqs and number of cov vals don't match".to_string(),
            ));
        }

        let mut g = DeBruijnGraph::new(kmer_length);

        for (s_idx, seq_region) in seq_regions.iter().enumerate() {
            let cov_val = atoi(iworm_cov_vals_str_vec[s_idx]);

            // < kmer_length 跳过（jaccard-clip 模式偶发）；cov 对应关系
            // 保持原位（错位行为与原版一致）
            if seq_region.len() < kmer_length {
                continue;
            }

            let seq_region = replace_non_gatc_with_a(seq_region.as_bytes());
            // accession^s 作为注解键（annotations 不移植，仅保留调用形参）
            let s_acc_val = format!("{accession}^{s_idx}");
            g.add_sequence(&s_acc_val, &seq_region, s_strand, cov_val as u32)?;
        }

        out += &g.to_chrysalis_format(component_id, s_strand);
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use trinity_common::kmer::{get_ds_kmer_val, kmer_to_intval};

    const A24: &[u8] = b"AAAAAAAAAAAAAAAAAAAAAAAA";
    const T24: &[u8] = b"TTTTTTTTTTTTTTTTTTTTTTTT";
    // A + 23T（rc(23A+T)）
    const AT23: &[u8] = b"ATTTTTTTTTTTTTTTTTTTTTTT";
    // 23A + T
    const A23T: &[u8] = b"AAAAAAAAAAAAAAAAAAAAAAAT";

    fn key(seq: &[u8]) -> KmerId {
        kmer_to_intval(seq).unwrap()
    }

    #[test]
    fn ss_single_chain_nodes_edges_counts() {
        // 26-mer → 3 个 24-mer 链
        let seq = b"AAAAAAAAAAAAAAAAAAAAAAAATT"; // A^24 + T + T（26-mer → 3 个 24-mer）
        assert_eq!(seq.len(), 26);
        let mut g = DeBruijnGraph::new(24);
        g.add_sequence("acc", seq, true, 5).unwrap();
        g.add_sequence("acc", seq, true, 7).unwrap();

        assert_eq!(g.len(), 3);
        // cov 累加（整条 contig 覆盖度，非 +1）
        let ids: Vec<(KmerId, i64)> = g
            .kmer_map
            .values()
            .map(|d| (d.kmer, d.kmer_count))
            .collect();
        for &(_, c) in &ids {
            assert_eq!(c, 12);
        }
        // 中间节点：prev = 首碱基 A 位、next = 末碱基 T 位
        let mid = key(A23T);
        let dk = g.get_node(mid).unwrap();
        assert_eq!(dk.prev_mask(), A_MASK);
        assert_eq!(dk.next_mask(), T_MASK);
        // 首/尾节点 id 为首次插入序
        let first = g.get_node(key(b"AAAAAAAAAAAAAAAAAAAAAAAA")).unwrap();
        let last = g.get_node(key(b"AAAAAAAAAAAAAAAAAAAAAATT")).unwrap();
        assert_eq!((first.id, mid_key_id(&g, mid), last.id), (1, 2, 3));
    }

    fn mid_key_id(g: &DeBruijnGraph, k: KmerId) -> u64 {
        g.get_node(k).unwrap().id()
    }

    #[test]
    fn ds_flipped_chain_direction_swap() {
        // S = A^24 + T（25-mer，2 个 kmer）：k1 = A^24（非 canonical，翻转），
        // k2 = A^23T（rc = A T^23 更大，也翻转）。
        let seq = b"AAAAAAAAAAAAAAAAAAAAAAAAT";
        assert_eq!(seq.len(), 25);
        let mut g = DeBruijnGraph::new(24);
        g.add_sequence("acc", seq, false, 1).unwrap();

        assert_eq!(g.len(), 2);
        // 手推（见任务报告）：节点 rc(k1)=T^24 → prev=A 位、next=0；
        // 节点 rc(k2)=A T^23 → prev=0、next=T 位（前驱 revcomp 后变后继边）。
        let n1 = g.get_node(get_ds_kmer_val(key(A24), 24)).unwrap();
        assert_eq!((n1.prev_mask(), n1.next_mask()), (A_MASK, 0));
        let n2 = g.get_node(get_ds_kmer_val(key(A23T), 24)).unwrap();
        assert_eq!((n2.prev_mask(), n2.next_mask()), (0, T_MASK));
        assert_eq!(get_ds_kmer_val(key(A23T), 24), key(AT23));
    }

    #[test]
    fn revcomp_of_a24_is_t24() {
        assert_eq!(trinity_common::kmer::revcomp_val(key(A24), 24), key(T24));
    }

    #[test]
    fn masks_and_candidates() {
        let mut dk = DeBruijnKmer::new(key(b"AC"), 1); // K=2, intval = 1*4+3
        dk.prev = G_MASK;
        dk.next = T_MASK;
        // prev 候选：G 拼上 k>>2 → "GA"
        assert_eq!(dk.get_prev_kmers(2), vec![key(b"GA")]);
        // next 候选：(k<<2)&mask | T → "CT"
        assert_eq!(dk.get_next_kmers(2), vec![key(b"CT")]);
        // add_prev 取首碱基；add_next 取末碱基
        let mut dk2 = DeBruijnKmer::new(key(b"GA"), 2);
        dk2.add_prev_kmer(key(b"TC"), 2); // 首碱基 T → T_MASK
        assert_eq!(dk2.prev_mask(), T_MASK);
        dk2.add_next_kmer(key(b"AC")); // 末碱基 C → C_MASK
        assert_eq!(dk2.next_mask(), C_MASK);
        // 候选序 G,A,T,C（多掩码位）
        let mut dk3 = DeBruijnKmer::new(key(b"AC"), 3);
        dk3.prev = G_MASK | C_MASK;
        assert_eq!(dk3.get_prev_kmers(2), vec![key(b"GA"), key(b"CA")]);
    }

    #[test]
    fn root_detection_ss_and_ds() {
        let mut g = DeBruijnGraph::new(24);
        g.add_sequence("acc", b"AAAAAAAAAAAAAAAAAAAAAAAATT", true, 1)
            .unwrap();
        let roots = g.get_root_kmers(true);
        assert_eq!(roots.len(), 1); // 仅链首（无 prev）
        assert_eq!(roots[0].kmer, key(b"AAAAAAAAAAAAAAAAAAAAAAAA"));

        // DS：无 next 者也算 root
        let mut g2 = DeBruijnGraph::new(24);
        g2.add_sequence("acc", b"AAAAAAAAAAAAAAAAAAAAAAAAT", false, 1)
            .unwrap();
        let roots2 = g2.get_root_kmers(false);
        assert_eq!(roots2.len(), 2); // prev 空者 + next 空者
    }

    #[test]
    fn to_chrysalis_ss_single_chain_hand_trace() {
        let mut g = DeBruijnGraph::new(24);
        g.add_sequence("acc", b"AAAAAAAAAAAAAAAAAAAAAAAATT", true, 1)
            .unwrap();
        let out = g.to_chrysalis_format(42, true);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 4); // Component 头 + 3 节点行
        assert_eq!(lines[0], "Component 42");
        // 根节点无 prev → -1 终端行；其后每行指向左邻 id
        assert_eq!(lines[1], "1\t-1\t1\tAAAAAAAAAAAAAAAAAAAAAAAA\t1");
        assert_eq!(lines[2], "2\t1\t1\tAAAAAAAAAAAAAAAAAAAAAAAT\t1");
        assert_eq!(lines[3], "3\t2\t1\tAAAAAAAAAAAAAAAAAAAAAATT\t1");
    }

    #[test]
    fn to_chrysalis_ds_negative_orientation_id_offset_and_terminal_line() {
        // 沿用 ds_flipped 链：root1 = A T^23（无 prev）→ 行 "2 -1 …"；
        // T^24（'+'，id 1）左行指向 A T^23（'+' id 2）；
        // 终端清理把 rc 标 seen → A T^23 再次弹出时输出 -1 行。
        let mut g = DeBruijnGraph::new(24);
        g.add_sequence("acc", b"AAAAAAAAAAAAAAAAAAAAAAAAT", false, 1)
            .unwrap();
        let out = g.to_chrysalis_format(0, false);
        let lines: Vec<&str> = out.lines().collect();
        // 手推（3 行）：root1 = A T^23（无 prev，id2）→ "-1" 终端行；终端清理
        // 把 rc(A T^23) = A^23 T 标 seen。弹出 T^24（'+'，id1）左行指向
        // A T^23（'+'，id2）。A T^23 再弹出时已 seen 且 reported（首访已输出），
        // 故无第二条 -1 行——"seen 未 reported" 的 -1 行只在 rc 定向被再弹出
        // 时出现（本图不触发）。
        let expect = vec![
            "Component 0",
            "2\t-1\t1\tATTTTTTTTTTTTTTTTTTTTTTT\t1",
            "1\t2\t1\tTTTTTTTTTTTTTTTTTTTTTTTT\t1",
        ];
        assert_eq!(lines, expect);

        // 负定向 id 偏移：K=2，"GTA" → k1=GT（翻转，存 AC id1）、k2=TA
        //（自 rc，canonical，id2）。手推边：k1 翻转 → next 边经 revcomp(TA)=TA
        // 落到 **prev** 位（首碱基 T）；k2 未翻转 → prev 边 add_prev(GT)（首碱基 G）。
        let mut g2 = DeBruijnGraph::new(2);
        g2.add_sequence("a", b"GTA", false, 3).unwrap();
        assert_eq!(g2.len(), 2);
        let ac = g2.get_node(key(b"AC")).unwrap();
        let ta = g2.get_node(key(b"TA")).unwrap();
        assert_eq!((ac.id(), ac.prev_mask(), ac.next_mask()), (1, T_MASK, 0));
        assert_eq!((ta.id(), ta.prev_mask(), ta.next_mask()), (2, G_MASK, 0));
        // 遍历（count 全等 → root 按 kmer 升序：AC 先）：
        // AC('+',id1) 左行 → TA('+',id2)；TA('+',id2) 左行候选 GT **非 canonical**
        // → '-' 定向 pkd_id = 1 + N(=2) = 3；GT 再弹出时 seen 未 reported（rc 清理）
        // → '-' 定向 dk_id = 3 的 -1 行。
        let out2 = g2.to_chrysalis_format(1, false);
        let lines2: Vec<&str> = out2.lines().collect();
        assert_eq!(
            lines2,
            vec![
                "Component 1",
                "1\t2\t1\tAC\t1",
                "2\t3\t1\tTA\t1",
                "3\t-1\t1\tGT\t1",
            ]
        );
    }

    #[test]
    fn graph_per_record_regions_covs_and_filters() {
        use trinity_common::fasta::FastaRecord;
        // X 切分 + cov 对应
        let rec = FastaRecord::new(
            ">s_7 5 9",
            format!(
                "{}X{}",
                String::from_utf8_lossy(b"AAAAAAAAAAAAAAAAAAAAAAAA"),
                String::from_utf8_lossy(b"CCCCCCCCCCCCCCCCCCCCCCCC")
            ),
        );
        let out = graph_per_record(&[rec], 24, false).unwrap();
        assert!(out.starts_with("Component 7\n"));
        let mut g = DeBruijnGraph::new(24);
        g.add_sequence("x", b"AAAAAAAAAAAAAAAAAAAAAAAA", false, 5)
            .unwrap();
        g.add_sequence("x", b"CCCCCCCCCCCCCCCCCCCCCCCC", false, 9)
            .unwrap();
        assert_eq!(out, g.to_chrysalis_format(7, false));

        // <24 区域跳过；cov 对应错位保留（第二个 cov 落到跳过位之后的区域不变）
        let rec2 = FastaRecord::new(
            ">s_1 5 3",
            "ACXAAAAAAAAAAAAAAAAAAAAAAAA".to_string(), // 区域1 "AC" 跳过
        );
        let out2 = graph_per_record(&[rec2], 24, true).unwrap();
        let mut g2 = DeBruijnGraph::new(24);
        g2.add_sequence("x", b"AAAAAAAAAAAAAAAAAAAAAAAA", true, 3)
            .unwrap();
        assert_eq!(out2, g2.to_chrysalis_format(1, true)); // cov 5 随 "AC" 丢弃

        // 非 gatc 替换为 A
        let rec3 = FastaRecord::new(">s_2 4", "NNNNNNNNNNNNNNNNNNNNNNNN".to_string());
        let out3 = graph_per_record(&[rec3], 24, true).unwrap();
        let mut g3 = DeBruijnGraph::new(24);
        g3.add_sequence("x", b"AAAAAAAAAAAAAAAAAAAAAAAA", true, 4)
            .unwrap();
        assert_eq!(out3, g3.to_chrysalis_format(2, true));

        // 数量不一致 → Err
        let rec4 = FastaRecord::new(">s_3 1 2", "AAAAAAAAAAAAAAAAAAAAAAAAX".to_string());
        assert!(graph_per_record(&[rec4], 24, true).is_err());
    }

    #[test]
    fn atoi_behaviour() {
        assert_eq!(atoi("s_12"), 0);
        assert_eq!(atoi("12"), 12);
        assert_eq!(atoi("12x9"), 12); // atoi 停在首个非数字
        assert_eq!(atoi(""), 0);
        assert_eq!(atoi("-5"), -5);
    }
}
