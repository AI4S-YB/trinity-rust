//! `Chrysalis/analysis/BubbleUpClustering.cc` 的移植——单链接 weld 聚类
//! （grow_prioritized_clusters）+ 未聚类补齐 + 长度和过滤 + COMPONENT 块输出。
//!
//! 复刻的三个 quirk（均有测试锁定）：
//!
//! 1. **+2 合并判据**：两簇合并条件是 `sizeA + sizeB + 2 <= MAX_CLUSTER_SIZE`
//!    （不是 `<= MAX`）——簇内 iworm 数可达 MAX-2+... 组合，与单端加入的
//!    `size < MAX` 判据不一致，按原样保留。
//! 2. **EOF 伪边**：原版 `while (!in.eof()) { getline(...) }` 在文件以换行结尾时
//!    会多跑一轮空行解析，得到 `0 -> 0` 伪边。若 iworm 0 在此前未被聚类，
//!    伪边会建出 `Pool[0,0]`（Pool::add 不去重！）——iworm 0 最终输出**两次**
//!    （成簇、长度和、#POOL_INFO 均重复计）。本移植按 `text.split('\n')`
//!    天然包含尾部空串（等价于那轮空 getline），原样复刻。
//! 3. **PrintSeq 空行**：80 整倍长序列在 80 列折行后再来一个 `"\n"`，
//!    产生一条空行。

use std::collections::HashMap;

use trinity_common::error::CommonError;

use crate::dna_vector::DnaSeq;

/// BubbleUpClustering.cc:26-27 静态默认；Trinity 管线实际传
/// `-min_contig_length 200 -max_cluster_size 25`。
#[derive(Debug, Clone)]
pub struct BubbleParams {
    pub min_contig_length: usize,
    pub max_cluster_size: usize,
    #[allow(dead_code)]
    pub debug_weld_all: bool,
}

impl Default for BubbleParams {
    fn default() -> Self {
        Self {
            min_contig_length: 24,
            max_cluster_size: 25,
            debug_weld_all: false,
        }
    }
}

/// `atoi` 语义：跳过前导空白、可选正负号，取连续数字，其余忽略，失败得 0。
fn atoi(s: &str) -> usize {
    let t = s.trim_start();
    let bytes = t.as_bytes();
    let mut i = 0;
    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        i += 1;
    }
    let mut val: i64 = 0;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        val = val
            .saturating_mul(10)
            .saturating_add((bytes[i] - b'0') as i64);
        i += 1;
    }
    if t.starts_with('-') {
        0 // node id 为 unsigned/usize，负数绕回语义无意义，取 0
    } else {
        val as usize
    }
}

/// grow_prioritized_clusters（BubbleUpClustering.cc:117-214）。
///
/// 输入 `weld_graph_sorted` **假定已按 total 降序排序**（管线责任，
/// 见 [`crate::graph_from_fasta::sort_weld_graph`]）。每行只读前 3 个 token
/// （tok0 = 左端点、tok1 = "->" 丢弃、tok2 = 右端点）；空行/缺 token 行
/// 按 atoi("") = 0 解析为 `0 -> 0`（含 EOF 伪边，见模块文档）。
fn grow_prioritized_clusters(weld_graph_sorted: &str, max_cluster_size: usize) -> Vec<Vec<usize>> {
    let mut clustered: Vec<Vec<usize>> = Vec::new();
    let mut id_to_cluster: HashMap<usize, usize> = HashMap::new();

    // split('\n') 的尾部空串 == 原版 eof 误判多出的那轮空 getline。
    for line in weld_graph_sorted.split('\n') {
        let mut tok = line.split_whitespace();
        let node_id = tok.next().map_or(0, atoi);
        let _arrow = tok.next();
        let linked = tok.next().map_or(0, atoi);

        let a = id_to_cluster.get(&node_id).copied();
        let b = id_to_cluster.get(&linked).copied();
        match (a, b) {
            (Some(ca), Some(cb)) => {
                if ca == cb {
                    // 已同簇（传递闭合）→ 跳过
                    continue;
                }
                // quirk：+2 判据
                if clustered[ca].len() + clustered[cb].len() + 2 <= max_cluster_size {
                    let migrate: Vec<usize> = clustered[cb].clone();
                    for n in migrate {
                        clustered[ca].push(n); // Pool::add 不去重
                        id_to_cluster.insert(n, ca);
                    }
                    clustered[cb].clear(); // 槽清空但保留索引
                }
                // 否则忽略该连接（过度聚合）
            }
            (Some(ca), None) => {
                if clustered[ca].len() < max_cluster_size {
                    clustered[ca].push(linked);
                    id_to_cluster.insert(linked, ca);
                }
            }
            (None, Some(cb)) => {
                if clustered[cb].len() < max_cluster_size {
                    clustered[cb].push(node_id);
                    id_to_cluster.insert(node_id, cb);
                }
            }
            (None, None) => {
                let idx = clustered.len();
                clustered.push(vec![node_id, linked]);
                id_to_cluster.insert(node_id, idx);
                id_to_cluster.insert(linked, idx);
            }
        }
    }

    // 剔除被合并清空的槽（保序）
    clustered.into_iter().filter(|p| !p.is_empty()).collect()
}

/// add_unclustered_iworm_contigs（:500-535）：未出现在任何簇中的 iworm
/// 按下标升序各成单元素簇，追加尾部。
fn add_unclustered(clustered: &mut Vec<Vec<usize>>, n_seqs: usize) {
    let mut found = std::collections::HashSet::new();
    for p in clustered.iter() {
        found.extend(p.iter().copied());
    }
    for i in 0..n_seqs {
        if !found.contains(&i) {
            clustered.push(vec![i]);
        }
    }
}

/// PrintSeq（:38-46）：80 列折行 + 结尾 "\n"——80 整倍长度产生额外空行（quirk）。
fn print_seq(seq: &[u8], out: &mut String) {
    for (i, &c) in seq.iter().enumerate() {
        out.push(c as char);
        if (i + 1) % 80 == 0 {
            out.push('\n');
        }
    }
    out.push('\n');
}

/// BubbleUpClustering 主流程：返回 COMPONENT 块全文（stdout 等价）。
pub fn bubble_up_clustering(
    iworm_seqs: &[DnaSeq],
    weld_graph_sorted: &str,
    p: &BubbleParams,
) -> Result<String, CommonError> {
    let mut clustered = if p.debug_weld_all {
        vec![(0..iworm_seqs.len()).collect::<Vec<usize>>()]
    } else {
        grow_prioritized_clusters(weld_graph_sorted, p.max_cluster_size)
    };
    if !p.debug_weld_all {
        add_unclustered(&mut clustered, iworm_seqs.len());
    }

    let mut out = String::new();
    let mut component_count = 0usize;
    for members in &mut clustered {
        if members.is_empty() {
            continue;
        }
        members.sort_unstable(); // sortvec

        // 长度和过滤（重复成员重复计数，复刻原版 quirk）
        let sum: usize = members.iter().map(|&z| iworm_seqs[z].seq.len()).sum();
        if sum < p.min_contig_length {
            continue; // 整簇丢弃，不占 component 号
        }

        let mut pool_info = format!("#POOL_INFO\t{component_count}:\t");
        out.push_str(&format!("COMPONENT {component_count}\t{}\n", members.len()));
        for &z in members.iter() {
            pool_info.push_str(&format!("{z} "));
            out.push_str(&format!(
                ">Component_{component_count} {} {z} [iworm{}]\n",
                members.len(),
                iworm_seqs[z].name
            ));
            print_seq(&iworm_seqs[z].seq, &mut out);
        }
        pool_info.push('\n');
        out.push_str(&pool_info);
        out.push_str("END\n");
        component_count += 1;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dna(name: &str, seq: &str) -> DnaSeq {
        DnaSeq {
            name: name.into(),
            seq: seq.as_bytes().to_vec(),
        }
    }

    /// 不带尾换行（getline 读完末行即置 eof，无伪边）；
    /// EOF 伪边场景由测试各自显式补 '\n'。
    fn weld(lines: &[&str]) -> String {
        lines.join("\n")
    }

    fn clusters(weld_text: &str, max: usize) -> Vec<Vec<usize>> {
        grow_prioritized_clusters(weld_text, max)
    }

    #[test]
    fn atoi_semantics() {
        assert_eq!(atoi("42"), 42);
        assert_eq!(atoi("7x9"), 7);
        assert_eq!(atoi(""), 0);
        assert_eq!(atoi("->"), 0);
    }

    #[test]
    fn both_unclustered_new_pool() {
        let c = clusters(&weld(&["5 -> 6 total: 9"]), 25);
        assert_eq!(c, vec![vec![5, 6]]);
    }

    #[test]
    fn same_cluster_continue() {
        let c = clusters(&weld(&["5 -> 6 total: 9", "6 -> 5 total: 8"]), 25);
        assert_eq!(c, vec![vec![5, 6]]);
    }

    #[test]
    fn single_end_add_bounded_by_max() {
        // 簇 [0,1]，MAX=3：size 2 < 3 → 加入 2 → [0,1,2]
        let c = clusters(&weld(&["0 -> 1 total: 9", "0 -> 2 total: 5"]), 3);
        assert_eq!(c, vec![vec![0, 1, 2]]);
        // 再加 3：size 3 == MAX → 拒绝
        let c = clusters(
            &weld(&["0 -> 1 total: 9", "0 -> 2 total: 5", "1 -> 3 total: 1"]),
            3,
        );
        assert_eq!(c, vec![vec![0, 1, 2]]);
    }

    #[test]
    fn merge_plus_two_quirk() {
        // MAX=4：A={0,1} B={2,3} → 2+2+2=6 > 4 → 不合并（若无 +2 quirk，4<=4 会合并）
        let c = clusters(
            &weld(&["0 -> 1 total: 9", "2 -> 3 total: 8", "0 -> 2 total: 5"]),
            4,
        );
        assert_eq!(c, vec![vec![0, 1], vec![2, 3]]);
        // MAX=5：2+2+2=6 > 5 → 仍不合并（quirk-free 4<=5 会合并）
        let c = clusters(
            &weld(&["0 -> 1 total: 9", "2 -> 3 total: 8", "0 -> 2 total: 5"]),
            5,
        );
        assert_eq!(c, vec![vec![0, 1], vec![2, 3]]);
        // MAX=6：2+2+2=6 <= 6 → 合并（保持 A 的槽位，B 槽清空后剔除）
        let c = clusters(
            &weld(&["0 -> 1 total: 9", "2 -> 3 total: 8", "0 -> 2 total: 5"]),
            6,
        );
        assert_eq!(c, vec![vec![0, 1, 2, 3]]);
    }

    #[test]
    fn eof_pseudo_edge_duplicates_iworm_zero() {
        // iworm 0 全图无边 → 尾换行触发 EOF 伪边 0->0 建 Pool[0,0]（重复成员）
        let c = clusters(&(weld(&["1 -> 2 total: 9"]) + "\n"), 25);
        assert_eq!(c, vec![vec![1, 2], vec![0, 0]]);
    }

    #[test]
    fn eof_pseudo_edge_no_trailing_newline_absent() {
        // 无尾换行（getline 后 eof 置位）→ 不产生伪边
        let c = clusters("1 -> 2 total: 9", 25);
        assert_eq!(c, vec![vec![1, 2]]);
    }

    #[test]
    fn add_unclustered_appends_singletons() {
        let mut c = vec![vec![3, 1]];
        add_unclustered(&mut c, 5);
        assert_eq!(c, vec![vec![3, 1], vec![0], vec![2], vec![4]]);
    }

    #[test]
    fn print_seq_eighty_multiple_extra_blank_line() {
        let mut s = String::new();
        print_seq(&[b'A'; 160], &mut s);
        // 两行 80 + 折行后的结尾 '\n' = 一条空行
        assert_eq!(s, format!("{}\n{}\n\n", "A".repeat(80), "A".repeat(80)));
        let mut s = String::new();
        print_seq(&[b'A'; 81], &mut s);
        assert_eq!(s, format!("{}\nA\n", "A".repeat(80)));
    }

    #[test]
    fn length_sum_filter_and_component_numbering() {
        let seqs = vec![dna(">a1;1_x", "AAAA"), dna(">a2;1_x", "CCCCCCCC")];
        let p = BubbleParams {
            min_contig_length: 10,
            ..Default::default()
        };
        // 簇 {0,1} 和 = 12 >= 10 → 保留；簇 {1} 单元素长 8 < 10 → 滤
        let seqs2 = vec![
            dna(">a1;1_x", "AAAA"),
            dna(">a2;1_x", "CCCCCCCC"),
            dna(">a3;1_x", "GGGG"),
        ];
        let out = bubble_up_clustering(&seqs2, &weld(&["0 -> 1 total: 9"]), &p).unwrap();
        assert!(out.starts_with("COMPONENT 0\t2\n"));
        assert!(!out.contains("COMPONENT 1"));
        // 全滤 → 空输出，component 号不占
        let p2 = BubbleParams {
            min_contig_length: 100,
            ..Default::default()
        };
        assert_eq!(
            bubble_up_clustering(&seqs, &weld(&["0 -> 1 total: 9"]), &p2).unwrap(),
            ""
        );
    }

    #[test]
    fn component_block_format() {
        let seqs = vec![
            dna(">a1;43_total_counts:_9", "ACGT"),
            dna(">a2;7_total_counts:_9", "TTTT"),
        ];
        let p = BubbleParams {
            min_contig_length: 1,
            ..Default::default()
        };
        let out = bubble_up_clustering(&seqs, &weld(&["1 -> 0 total: 9"]), &p).unwrap();
        let expected = "\
COMPONENT 0\t2
>Component_0 2 0 [iworm>a1;43_total_counts:_9]
ACGT
>Component_0 2 1 [iworm>a2;7_total_counts:_9]
TTTT
#POOL_INFO\t0:\t0 1 \n\
END
";
        assert_eq!(out, expected);
    }

    #[test]
    fn debug_weld_all_single_pool() {
        let seqs = vec![dna(">a", "AAAA"), dna(">b", "CCCC")];
        let p = BubbleParams {
            debug_weld_all: true,
            min_contig_length: 1,
            ..Default::default()
        };
        let out = bubble_up_clustering(&seqs, "", &p).unwrap();
        assert!(out.contains("COMPONENT 0\t2\n"));
        assert!(out.contains("#POOL_INFO\t0:\t0 1 \n"));
    }
}
