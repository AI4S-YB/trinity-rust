//! rnaSPAdes 式路径引擎（`--path-engine rnaspades`）——PathExtend MultiExtender
//! 核心规则的组件级复刻（B4, 收窄 scope 版）。
//!
//! 源码依据（spades-main）:
//! - `path_extenders.cpp:303-347` `MultiExtender::AddCandidates`: 0 候选 → 停;
//!   1 候选 → 直接延伸; **恰好 2 候选且同起点同终点（简单 bulge）→ 复制路径
//!   各走一支**（产 2 个 isoform）; 其余（2 候选非 bulge / ≥3 候选）→ 停止
//!   延伸（宁可断不可错连）。
//! - `rna_mode.info:155-197` path_cleaning default 档: `min_length 110` +
//!   `min_coverage 2`（isolated_/rel_ 各档未复刻）。
//!
//! 图语义映射: SPAdes 图经压缩, "边"是 unitig, bulge = 同起点同终点的两条
//! 平行边; 本组件图（graph.out, [`parse_component_graph`] 解析）**节点 =
//! kmer(24/25)、边 = (k-1) 重叠**, 未压缩——"简单 bulge"在此对应两条内部
//! 简单（出度 1）的分叉链在下游重汇合到同一节点 C（SNP bubble, 支长 ≈ k）。
//! 见 [`bulge_merge`]。
//!
//! **相对 rnaSPAdes 的简化点（明示）**:
//! - 不做 PE pair 权重（RNAExtensionChooser/PathCoverWeightCounter,
//!   extension_chooser.hpp:427-473/573-591 未复刻）; 候选排序用**覆盖一致
//!   性**: 候选得分 = 1/(1+|cnt−路径均覆盖|), 保留 score ≥ max/1.5
//!   （priority_coeff, pe_params.info 默认）者, 并列全部保留以触发 bulge
//!   判定。组件 reads 因此完全不消费。
//! - 无 LoopDetector/used_unique_storage; 以"路径内节点不复走"防环。
//! - path_cleaning 只实现 default 档两条硬阈值; isolated/relative 过滤与
//!   包含路径删除未复刻（harvest 层另有全序列去重兜底）。

use rustc_hash::{FxHashMap, FxHashSet};

use crate::gtdecode::{parse_component_graph, GtGraph};
use crate::Result;

/// 引擎参数（默认 = rna_mode.info path_cleaning default 档 + pe_params.info
/// extension_options 的 priority_coeff）。
#[derive(Debug, Clone)]
pub struct RnaSpadesParams {
    /// path_cleaning: 最短输出长度 bp（rna_mode.info default `min_length 110`）。
    pub min_len: usize,
    /// path_cleaning: 路径平均 kmer 覆盖下限（default `min_coverage 2`）。
    pub min_cov: f64,
    /// 候选并列保留系数（pe_params.info `priority_coeff 1.5`）: 保留
    /// score ≥ max_score / priority_coeff 的全部候选。
    pub priority_coeff: f64,
    /// 每组件路径总数安全阀（SPAdes 无此帽; 防串联 bulge 2^n 爆炸）。
    pub max_paths: usize,
    /// bulge 支链长度上限（节点数; SNP bubble 支长 ≈ k, 50 富余）。
    pub bulge_max_chain: usize,
}

impl Default for RnaSpadesParams {
    fn default() -> Self {
        Self {
            min_len: 110,
            min_cov: 2.0,
            priority_coeff: 1.5,
            max_paths: 4096,
            bulge_max_chain: 50,
        }
    }
}

/// 组件图邻接（节点级, 平行边去重, 哨兵 -1 剔除）。
struct Adj {
    /// u → 后继节点（首现序, 去重, 不含 -1）。
    succ: FxHashMap<i32, Vec<i32>>,
    /// 真实前驱数（不含来自 <0 哨兵的 root 边）。
    indeg: FxHashMap<i32, u32>,
}

fn build_adj(g: &GtGraph) -> Adj {
    let mut succ: FxHashMap<i32, Vec<i32>> = FxHashMap::default();
    let mut indeg: FxHashMap<i32, u32> = FxHashMap::default();
    let mut pred_seen: FxHashSet<(i32, i32)> = FxHashSet::default();
    for ei in 0..g.n_edges() {
        let (s, d) = (g.vids[g.src[ei]], g.vids[g.dst[ei]]);
        if d < 0 {
            continue; // 防御: 哨兵不作路径节点
        }
        let e = succ.entry(s).or_default();
        if !e.contains(&d) {
            e.push(d);
        }
        if s >= 0 && pred_seen.insert((s, d)) {
            *indeg.entry(d).or_insert(0) += 1;
        }
    }
    Adj { succ, indeg }
}

/// 生长中的路径（节点列 + 成员集 + 覆盖累计）。
#[derive(Clone)]
struct Growing {
    nodes: Vec<i32>,
    set: FxHashSet<i32>,
    cov_sum: f64,
}

impl Growing {
    fn new(seed: i32, cov: f64) -> Self {
        let mut set = FxHashSet::default();
        set.insert(seed);
        Self {
            nodes: vec![seed],
            set,
            cov_sum: cov,
        }
    }

    fn push(&mut self, v: i32, cov: f64) {
        self.nodes.push(v);
        self.set.insert(v);
        self.cov_sum += cov;
    }

    /// 当前路径平均 kmer 覆盖。
    fn avg_cov(&self) -> f64 {
        self.cov_sum / self.nodes.len() as f64
    }
}

/// 覆盖一致性得分: 候选覆盖越贴近路径均覆盖得分越高（1/(1+|差|), ∈(0,1]）。
fn cov_score(cand_cov: f64, path_avg: f64) -> f64 {
    1.0 / (1.0 + (cand_cov - path_avg).abs())
}

/// 简单 bulge 判定（`MultiExtender::AddCandidates` 2 候选情形在未压缩 kmer
/// 图上的翻译）: 从两支各自沿**唯一后继**前进（支内节点出度恒 1, 否则非
/// 简单 bulge）, 两支在 `max_chain` 节点内首次到达**同一**节点 C 且
/// C 的入度恰为 2（仅这两支汇合）→ `Some(C)`; 死端/支内再分叉/超上限/
/// C 已在当前路径上 → `None`。
fn bulge_merge(
    adj: &Adj,
    b1: i32,
    b2: i32,
    path: &FxHashSet<i32>,
    max_chain: usize,
) -> Option<i32> {
    // 返回支链终点下游首个"汇合候选"节点（入度 != 1 者）; 支内死端/分叉/超限 → None
    let follow = |mut x: i32| -> Option<i32> {
        for _ in 0..max_chain {
            let nexts = adj.succ.get(&x)?;
            if nexts.len() != 1 {
                return None; // 死端或支内再分叉 → 非简单 bulge
            }
            let n = nexts[0];
            if adj.indeg.get(&n).copied().unwrap_or(0) != 1 {
                return Some(n); // n 为汇合候选
            }
            x = n;
        }
        None // 超支长上限
    };
    let (c1, c2) = (follow(b1)?, follow(b2)?);
    if c1 == c2 && adj.indeg.get(&c1).copied().unwrap_or(0) == 2 && !path.contains(&c1) {
        Some(c1)
    } else {
        None
    }
}

/// 路径生长主循环（`MultiExtender` 语义）: 从种子出发, 逐步延伸;
/// 每步候选 = 末节点的未走过后继, 按覆盖一致性过滤（保留 score ≥
/// max/priority_coeff 者）, 再按 0/1/2(bulge)/其他 分派。分叉路径入
/// worklist 独立生长。返回全部完成路径（节点列）。
pub fn grow_paths(g: &GtGraph, params: &RnaSpadesParams) -> Vec<Vec<i32>> {
    let adj = build_adj(g);
    let cnt = |v: i32| g.cnt.get(&v).copied().unwrap_or(0) as f64;

    // 种子: 无真实前驱的源点, 按 (-cnt, vid) 排序（确定性, 高覆盖先行）;
    // 纯环组件无源点 → 退化为最高覆盖单种子（取舍, SPAdes 同以源边起步）。
    let mut seeds: Vec<i32> = g
        .vids
        .iter()
        .copied()
        .filter(|&v| v >= 0 && adj.indeg.get(&v).copied().unwrap_or(0) == 0)
        .collect();
    seeds.sort_by_key(|&v| (std::cmp::Reverse(g.cnt.get(&v).copied().unwrap_or(0)), v));
    if seeds.is_empty() {
        if let Some(&top) = g
            .vids
            .iter()
            .filter(|&&v| v >= 0)
            .max_by_key(|&&v| g.cnt.get(&v).copied().unwrap_or(0))
        {
            seeds.push(top);
        }
    }

    let mut finished: Vec<Vec<i32>> = Vec::new();
    for seed in seeds {
        let mut worklist: Vec<Growing> = vec![Growing::new(seed, cnt(seed))];
        while let Some(mut p) = worklist.pop() {
            loop {
                let last = *p.nodes.last().expect("non-empty path");
                // 候选 = 未走过（防环）的后继
                let cands: Vec<i32> = adj
                    .succ
                    .get(&last)
                    .map(|v| {
                        v.iter()
                            .copied()
                            .filter(|c| !p.set.contains(c))
                            .collect()
                    })
                    .unwrap_or_default();
                if cands.is_empty() {
                    break;
                }
                // 覆盖一致性过滤（简化点: 无 PE 权重）
                let avg = p.avg_cov();
                let max_score = cands
                    .iter()
                    .map(|&c| cov_score(cnt(c), avg))
                    .fold(0.0f64, f64::max);
                let cutoff = max_score / params.priority_coeff;
                let top: Vec<i32> = cands
                    .iter()
                    .copied()
                    .filter(|&c| cov_score(cnt(c), avg) >= cutoff)
                    .collect();
                match top.len() {
                    0 => break, // max ≥ cutoff 恒成立, 防御
                    1 => {
                        let c = top[0];
                        p.push(c, cnt(c));
                    }
                    2 => {
                        if bulge_merge(&adj, top[0], top[1], &p.set, params.bulge_max_chain)
                            .is_some()
                        {
                            // 简单 bulge: 复制路径走 top[1], 当前路径走 top[0]
                            if finished.len() + worklist.len() < params.max_paths {
                                let mut q = p.clone();
                                let c1 = top[1];
                                q.push(c1, cnt(c1));
                                worklist.push(q);
                            } // 超安全阀: 放弃分叉只走主支（SPAdes 无此帽）
                            let c0 = top[0];
                            p.push(c0, cnt(c0));
                        } else {
                            break; // 非 bulge 双歧义 → 停（宁可断）
                        }
                    }
                    _ => break, // ≥3 歧义 → 停
                }
            }
            finished.push(p.nodes);
        }
    }
    finished
}

/// 拼序列（gtasm 同规）: 首个有名节点给全 kmer, 后续节点各补末字符;
/// 无名节点跳过（哨兵不入路径, 防御保留）。
fn path_sequence(g: &GtGraph, path: &[i32]) -> String {
    let mut seq = String::new();
    for &v in path {
        if let Some(km) = g.names.get(&v) {
            if seq.is_empty() {
                seq.push_str(km);
            } else {
                seq.push(km.as_bytes()[km.len() - 1] as char);
            }
        }
    }
    seq
}

/// 路径平均 kmer 覆盖（path_cleaning 的 min_coverage 判定基准）。
fn path_cov(g: &GtGraph, path: &[i32]) -> f64 {
    let named: Vec<i64> = path
        .iter()
        .filter(|v| g.names.contains_key(v))
        .map(|v| g.cnt.get(v).copied().unwrap_or(0))
        .collect();
    if named.is_empty() {
        return 0.0;
    }
    named.iter().sum::<i64>() as f64 / named.len() as f64
}

/// 单组件 rnaspades 全链（池调用入口）: graph.out 文本 → allProbPaths
/// FASTA 文本。FASTA 头与 classic/gtasm 同族（`>{name}_g{j}_i1 len=L
/// path=[rnaspades]`）——harvest 长度过滤/去重与 gene_trans_map 兼容。
///
/// path_cleaning（rna_mode.info default 档）: 序列长 ≥ 110 且平均 kmer
/// 覆盖 ≥ 2; 另做全序列去重（harvest 之前先收敛, 免重复警告刷屏）。
pub fn run_component_rnaspades(graph_text: &str, name: &str) -> Result<String> {
    run_component_rnaspades_with(graph_text, name, &RnaSpadesParams::default())
}

/// 参数化变体（单测用）。
pub fn run_component_rnaspades_with(
    graph_text: &str,
    name: &str,
    params: &RnaSpadesParams,
) -> Result<String> {
    let g = parse_component_graph(graph_text)?;
    if g.n_edges() == 0 {
        return Ok(String::new());
    }
    let paths = grow_paths(&g, params);
    let mut seen: FxHashSet<String> = FxHashSet::default();
    let mut out = String::new();
    let mut j = 0usize;
    for path in &paths {
        let seq = path_sequence(&g, path);
        // path_cleaning default 档硬阈值
        if seq.len() < params.min_len || path_cov(&g, path) < params.min_cov {
            continue;
        }
        if !seen.insert(seq.clone()) {
            continue;
        }
        j += 1;
        out.push_str(&format!(">{name}_g{j}_i1 len={} path=[rnaspades]\n", seq.len()));
        for chunk in seq.as_bytes().chunks(60) {
            // 序列由图 kmer 拼成, 恒为 ASCII
            out.push_str(std::str::from_utf8(chunk).unwrap_or(""));
            out.push('\n');
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const K: usize = 5;

    /// 由若干全长序列建 graph.out 文本: 每条序列贡献连续 kmer 行
    /// （节点 id 按 kmer 首现分配; cnt 由 cov_of 给出）。
    fn text_from_seqs(seqs: &[&str], cov_of: &dyn Fn(&str) -> i64) -> String {
        let mut id_of: FxHashMap<String, i32> = FxHashMap::default();
        let mut next = 1i32;
        let mut text = String::from("Component 0\n");
        for s in seqs {
            assert!(s.len() >= K);
            let mut prev: Option<i32> = None;
            for i in 0..=s.len() - K {
                let km = &s[i..i + K];
                let id = *id_of.entry(km.to_string()).or_insert_with(|| {
                    let id = next;
                    next += 1;
                    id
                });
                let from = prev.unwrap_or(-1);
                text.push_str(&format!("{id}\t{from}\t{}\t{km}\t1\n", cov_of(km)));
                prev = Some(id);
            }
        }
        text
    }

    fn const_cov(c: i64) -> impl Fn(&str) -> i64 {
        move |_| c
    }

    fn test_params() -> RnaSpadesParams {
        RnaSpadesParams {
            min_len: 5,
            min_cov: 0.0,
            ..Default::default()
        }
    }

    /// SNP bubble 测试图（k=5）: 共享前缀 "ATGCCGTA" + alt A/C + 共享后缀
    /// "TGGAACCTT"。支长 = k = 5 节点, 汇合点 = "TGGAA"（入度 2）。
    const SX: &str = "ATGCCGTAATGGAACCTT";
    const SY: &str = "ATGCCGTACTGGAACCTT";

    #[test]
    fn bulge_fork_yields_two_isoforms() {
        let text = text_from_seqs(&[SX, SY], &const_cov(10));
        let g = parse_component_graph(&text).unwrap();
        let paths = grow_paths(&g, &test_params());
        assert_eq!(paths.len(), 2, "简单 bulge 应分叉产出 2 条路径: {paths:?}");
        let mut seqs: Vec<String> = paths.iter().map(|p| path_sequence(&g, p)).collect();
        seqs.sort();
        assert_eq!(seqs, vec![SX.to_string(), SY.to_string()]);
    }

    #[test]
    fn non_bulge_ambiguity_stops() {
        // 两支分叉后不再汇合 → 在分歧点停止, 只报共享前缀
        let text = text_from_seqs(&["ATGCCGTAATTGGGCCC", "ATGCCGTACTTTAACG"], &const_cov(10));
        let g = parse_component_graph(&text).unwrap();
        let paths = grow_paths(&g, &test_params());
        assert_eq!(paths.len(), 1);
        assert_eq!(path_sequence(&g, &paths[0]), "ATGCCGTA");
    }

    #[test]
    fn three_way_ambiguity_stops() {
        // 三岔口（≥3 候选）→ 停
        let text = text_from_seqs(
            &["ATGCCGTAATTGG", "ATGCCGTACTTGG", "ATGCCGTAGTTGG"],
            &const_cov(10),
        );
        let g = parse_component_graph(&text).unwrap();
        let paths = grow_paths(&g, &test_params());
        assert_eq!(paths.len(), 1);
        assert_eq!(path_sequence(&g, &paths[0]), "ATGCCGTA");
    }

    #[test]
    fn coverage_consistency_picks_closest() {
        // 无 PE（简化点）: 覆盖一致性在双歧义中选出贴近路径覆盖的一支。
        // 主干 cnt=100, 干扰支 cnt=10 → 干扰支 score≈0.11 < 1.0/1.5 → 被滤,
        // 剩 1 候选继续延伸（若两支都保留则非 bulge 会停在分歧点）。
        let main = "ATGCCGTAATTGGGCCC";
        let decoy = "ATGCCGTACTTTAACG";
        let decoy_kmers: FxHashSet<String> = (0..=decoy.len() - K)
            .map(|i| decoy[i..i + K].to_string())
            .collect();
        let main_kmers: FxHashSet<String> =
            (0..=main.len() - K).map(|i| main[i..i + K].to_string()).collect();
        let cov = move |km: &str| {
            if decoy_kmers.contains(km) && !main_kmers.contains(km) {
                10
            } else {
                100
            }
        };
        let text = text_from_seqs(&[main, decoy], &cov);
        let g = parse_component_graph(&text).unwrap();
        let paths = grow_paths(&g, &test_params());
        assert_eq!(paths.len(), 1);
        assert_eq!(path_sequence(&g, &paths[0]), main);
    }

    #[test]
    fn cycle_component_uses_high_cov_seed_and_stops_at_loop() {
        // 纯环（8 节点 k=5 环）: 无源点 → 高覆盖种子兜底; 防环止于一周
        let ring = "ACGTTGCAACGTTGC"; // 线性化: 第 8 个 kmer 之后回到节点 0
        let text = text_from_seqs(&[ring], &const_cov(5));
        // 去掉 root 行（from=-1）: 重建为纯环——首行 from 改为末节点
        let mut lines: Vec<&str> = text.lines().collect();
        // 行 1 = 头; 行 2 = 首 kmer 的 root 行, 改成 from = 第 8 个节点 id (=8)
        let first_data = lines[1].to_string();
        let mut f = first_data.split('\t');
        let to = f.next().unwrap();
        assert_eq!(f.next(), Some("-1"));
        lines[1] = Box::leak(format!("{to}\t8\t5\t{}\t1", &ring[0..K]).into_boxed_str());
        let text = lines.join("\n") + "\n";
        let g = parse_component_graph(&text).unwrap();
        assert_eq!(g.vids.len(), 8);
        let paths = grow_paths(&g, &test_params());
        assert_eq!(paths.len(), 1);
        assert_eq!(paths[0].len(), 8, "环上每节点恰好走一次: {:?}", paths[0]);
    }

    /// 生成 n 节点无重复 kmer 的链序列（贪心选碱基, 保证图是单链）。
    fn uniq_chain(n_nodes: usize, mut x: u64) -> String {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut s: Vec<u8> = Vec::new();
        let mut seen: FxHashSet<String> = FxHashSet::default();
        for i in 0..(n_nodes + K - 1) {
            for attempt in 0..4u8 {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                let b = bases[((x + attempt as u64) % 4) as usize];
                s.push(b);
                if i + 1 >= K {
                    let km = std::str::from_utf8(&s[i + 1 - K..=i]).unwrap().to_string();
                    if seen.contains(&km) {
                        s.pop();
                        continue; // 该碱基造成重复 kmer, 换下一个
                    }
                    seen.insert(km);
                }
                break;
            }
        }
        String::from_utf8(s).unwrap()
    }

    #[test]
    fn cleaning_length_boundary() {
        // len 110 保留 / 109 丢弃（default min_length=110, 边界含）
        let s110 = uniq_chain(106, 0x9E3779B97F4A7C15);
        assert_eq!(s110.len(), 110);
        let text = text_from_seqs(&[&s110], &const_cov(5));
        let out = run_component_rnaspades(&text, "c0").unwrap();
        assert!(out.contains("len=110"), "{out}");
        assert!(out.contains("path=[rnaspades]"));
        // 109bp 丢弃
        let s109 = uniq_chain(105, 0x9E3779B97F4A7C15);
        assert_eq!(s109.len(), 109);
        let text = text_from_seqs(&[&s109], &const_cov(5));
        let out = run_component_rnaspades(&text, "c0").unwrap();
        assert!(out.is_empty(), "109bp 应被 min_length=110 滤掉: {out}");
    }

    #[test]
    fn cleaning_coverage_boundary() {
        let s = uniq_chain(116, 0xBF58476D1CE4E5B9);
        assert_eq!(s.len(), 120);
        // 平均覆盖恰 2.0 → 保留（>= 边界含）
        let text = text_from_seqs(&[&s], &const_cov(2));
        let out = run_component_rnaspades(&text, "c0").unwrap();
        assert!(out.contains("len=120"), "avg cov=2.0 边界应保留: {out}");
        // 平均覆盖 1 < 2 → 丢弃
        let text = text_from_seqs(&[&s], &const_cov(1));
        let out = run_component_rnaspades(&text, "c0").unwrap();
        assert!(out.is_empty(), "avg cov=1 应被 min_coverage=2 滤掉");
        // 混合 3/1 交替 → 均值在 [1,3] 区间, 边界路径保留判定由 cov 决定
        let cov = |km: &str| {
            if km.as_bytes()[0].is_multiple_of(2) {
                3
            } else {
                1
            }
        };
        let text = text_from_seqs(&[&s], &cov);
        let g = parse_component_graph(&text).unwrap();
        let paths = grow_paths(&g, &test_params());
        assert_eq!(paths.len(), 1);
        let c = path_cov(&g, &paths[0]);
        assert!((1.0..=3.0).contains(&c));
    }

    #[test]
    fn empty_graph_yields_empty() {
        assert_eq!(run_component_rnaspades("Component 0\n", "c0").unwrap(), "");
    }

    #[test]
    fn fasta_header_and_wrapping() {
        let s = uniq_chain(126, 0x2545F4914F6CDD1D);
        assert_eq!(s.len(), 130);
        let text = text_from_seqs(&[&s], &const_cov(7));
        let out = run_component_rnaspades(&text, "c3").unwrap();
        assert!(out.starts_with(">c3_g1_i1 len=130 path=[rnaspades]\n"), "{out}");
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 4); // 头 + 60 + 60 + 10
        assert_eq!(lines[1].len(), 60);
        assert_eq!(lines[3].len(), 10);
    }
}
