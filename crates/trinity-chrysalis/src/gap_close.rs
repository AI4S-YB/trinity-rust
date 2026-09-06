//! SPAdes 式 gap closer（`--gap-close-spades`，默认关 = 逐字节不变）。
//!
//! 镜像 `spades-main/src/projects/spades/gap_closer.cpp`（GapCloserPairedIndexFiller
//! L28-178 / GapCloser L180-469）与 `configs/config.info` 的 `gap_closer` 段
//! （`minimal_intersection 10`、`weight_threshold 2.0`、`max_dist_to_tip 5000`），
//! 在 chrysalis QuantifyGraph（与可选的 graph_clean）之后逐组件缝合
//! `c*.graph.out`：覆盖山谷把一条 contig 断成的两段，其 tip 序列往往（近）
//! 重叠——用 PE pair 证据给候选排序，再做 overlap 缝合。**不做局部重组装**
//! （python 原型已证伪多路径歧义），只加桥接边，假阳性低。
//!
//! **映射取舍**（我们的图是未压缩 de Bruijn 图，节点=24mer，边=1bp 漂移）：
//! - SPAdes 的 out-tip（`IsDeadEnd(EdgeEnd)`）≈ 本图的汇向死端链（末端节点
//!   outdeg==0 且有真实入边，沿唯一前驱回溯到分支点(不含)/源头）；in-tip
//!   ≈ 源向死端链（首节点 indeg==0，沿唯一后继走到汇合点(不含)）。SPAdes 的
//!   tip 邻域含汇合点处的全部入边，我们只取唯一主链（保守，少收票源）。
//! - `max_dist_to_tip 5000` 为压缩图长边设计；本图组件小，取链端 200bp
//!   （177 个 24mer）作 PE 命中窗口（任务书指定取舍，注明）。
//! - PE 证据：reads.tmp（该组件 reads，short name 含 `/1`、`/2` 后缀）中同一
//!   pair（name 去尾 `/1`/`/2`）的两个 mate 分落 out-tip A 末端区与 in-tip B
//!   首端区（任一 mate 在哪侧均可；非链特异模式正/反向都试，镜像
//!   QuantifyGraph 的双向计数）→ (A,B) 权重+1（SPAdes RawPoint weight 1.0）。
//!   权重 < 2.0 丢弃（`weight_threshold 2.0`）。
//! - 缝合（SPAdes `ProcessPair`）：gap 1..=k-min_intersection（overlap 自 k-1
//!   递减到 10），重叠区 = out-tip 末端节点 kmer 后缀 vs in-tip 首节点 kmer
//!   前缀。汉明 > bound 试下一 gap；低复杂度（单碱基计数 > ratio×overlap,
//!   ratio = 0.8+0.2×(gap-1)/(k-min-1) 线性, 严格 >）→ 整对放弃（SPAdes
//!   return false 语义）; 否则成败都不再试更小 overlap。
//! - `hamming_dist_bound`：config.info 未显式给，GapCloser 构造默认 0，构造
//!   注释建议 `min_intersection/5`；取 2（启用错配修正路径的小常数）。
//! - 完美匹配 → 机械桥接：gap=24-overlap，从 a 漂移 gap 步正好到 b，插入
//!   gap-1 个新节点 + gap 条新边（行 `to frm count kmer\t1`），桥接覆盖
//!   = max(1, min(两链覆盖))（保守下限, 表达 PE 证据支持）。
//! - 错配修正（`HandlePositiveHammingDistanceCase` 语义，任务书许可的简化：
//!   以覆盖高的一侧为权威；tie 右侧优先，贴近 SPAdes 先 CanCorrectLeft 的
//!   顺序）。未压缩图节点不可分割（SPAdes 的 SplitEdge 不可移植），可修正
//!   窗口只剩贴端单位置：右权威=错配唯一且在 overlap 末位（改写 a 的 kmer
//!   末位，前驱兼容不受影响）; 左权威=错配唯一且在 overlap 首位（改写 b 的
//!   kmer 首位，后继兼容不受影响）。窗口外 → 放弃。
//! - 防假阳性：同链/共节点对不缝（近似 SPAdes avoid_rc_connections 的
//!   first==second/conjugate 防护）；改写 kmer 与新节点 kmer 均查重（撞组件
//!   既有 kmer 或互撞 → 放弃，防 butterfly 按 kmer 并节点造成意外短路）；
//!   每 out-tip 只接受第一个可行候选（SPAdes break 语义），缝合前复查两端
//!   仍是死端/源头（SPAdes L434-437）。
//! - 缝合表示 = **新增桥接边行**（不改既有行的 to/frm；修正路径才改写末端
//!   节点定义行的 kmer 列），对下游 butterfly 就是普通图边——决策依据：
//!   butterfly 重建图后能看到 a→…→b 的完整路径，两段 contig 可连成一条
//!   转录本; 合并链序列则会丢掉图内覆盖结构，不自然。

use std::collections::{HashMap, HashSet};

use crate::dna_vector::{read_fasta_bytes, revcomp};
use crate::graph_clean::{build_adj, split_blocks, Adj, NodeLine};

/// 图节点 kmer 长（FastaToDeBruijn 固定 -K 24）。
pub const K: usize = 24;
/// config.info `gap_closer.minimal_intersection`（最小重叠 bp）。
const MIN_INTERSECTION: usize = 10;
/// config.info `gap_closer.weight_threshold`（PE 票数下限, 权重为整数对数）。
const WEIGHT_THRESHOLD: u32 = 2;
/// 重叠区汉明差上限（构造注释 `min_intersection/5`；config.info 未给, 默认 0）。
const HAMMING_DIST_BOUND: usize = 2;
/// tip 末端区长度（bp）——SPAdes max_dist_to_tip 5000 的小图取舍（任务书指定）。
const TIP_REGION_BP: usize = 200;
/// 末端区覆盖的节点数上限（200bp 窗口含 177 个重叠 24mer）。
const TIP_REGION_NODES: usize = TIP_REGION_BP - K + 1;

/// 缝合统计（跨组件累加，管线层汇总打印）。
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct GapCloseStats {
    /// 处理的组件数。
    pub components: usize,
    /// 汇向死端链（out-tip）数。
    pub out_tips: usize,
    /// 源向死端链（in-tip）数。
    pub in_tips: usize,
    /// 命中任一 tip 末端区的 read 数。
    pub tip_reads: usize,
    /// (pair, tip 对) 投票总数（权重和）。
    pub pair_votes: usize,
    /// 权重 ≥ 阈值的 (out,in) 候选对数。
    pub candidate_pairs: usize,
    /// 进入 overlap 检验的对数（复查死端/源头通过后）。
    pub pairs_checked: usize,
    /// 缝合成功总数。
    pub gaps_closed: usize,
    /// 其中完美 overlap。
    pub closed_perfect: usize,
    /// 其中错配修正后缝合。
    pub closed_corrected: usize,
    /// 低复杂度拒绝对数（SPAdes return false 语义, 整对放弃）。
    pub rejected_low_complexity: usize,
    /// 错配不可修正（窗口外）拒绝对数。
    pub rejected_hamming: usize,
    /// kmer 查重冲突拒绝对数。
    pub rejected_kmer_clash: usize,
    /// 新增桥接中间节点总数。
    pub bridge_nodes: usize,
}

impl GapCloseStats {
    /// 逐字段累加（rayon 并行下由管线层合并）。
    pub fn add_assign(&mut self, o: &GapCloseStats) {
        self.components += o.components;
        self.out_tips += o.out_tips;
        self.in_tips += o.in_tips;
        self.tip_reads += o.tip_reads;
        self.pair_votes += o.pair_votes;
        self.candidate_pairs += o.candidate_pairs;
        self.pairs_checked += o.pairs_checked;
        self.gaps_closed += o.gaps_closed;
        self.closed_perfect += o.closed_perfect;
        self.closed_corrected += o.closed_corrected;
        self.rejected_low_complexity += o.rejected_low_complexity;
        self.rejected_hamming += o.rejected_hamming;
        self.rejected_kmer_clash += o.rejected_kmer_clash;
        self.bridge_nodes += o.bridge_nodes;
    }
}

/// 一条死端链（out-tip: 汇向; in-tip: 源向）。
struct Tip {
    /// 缝合锚点：out-tip 的末端节点 / in-tip 的首节点。
    anchor: i64,
    /// 链覆盖 = 链内真实入边 count max（graph_clean 惯例; 无真实边为 0）。
    cov: i64,
    /// 链全部节点（防同链自连）。
    nodes: HashSet<i64>,
}

/// 汇向死端链回溯：从死端 v 沿唯一前驱走，含源头、停在分支点（不含）。
/// 返回图序节点列 [v1..vn=v]。环防护：节点不重复。
fn collect_out_chain(adj: &Adj, lines: &[(NodeLine, &str)], v: i64) -> Vec<i64> {
    let mut rev = vec![v];
    let mut cur = v;
    loop {
        if adj.indeg(cur) != 1 {
            break; // 源头（indeg==0）或汇合点（≥2）→ 链止于此（含 cur）
        }
        let p = lines[adj.in_idx[&cur][0]].0.frm;
        if p < 0 || rev.contains(&p) {
            break;
        }
        if adj.outdeg(p) == 1 {
            rev.push(p);
            cur = p;
        } else {
            break; // p 出度≥2 = 分支点，不含
        }
    }
    rev.reverse();
    rev
}

/// 源向死端链前走：从源头 s 沿唯一后继走，停在汇合点（不含）。
/// 返回图序节点列 [s=w1..wn]。
fn collect_in_chain(adj: &Adj, lines: &[(NodeLine, &str)], s: i64) -> Vec<i64> {
    let mut nodes = vec![s];
    let mut cur = s;
    loop {
        if adj.outdeg(cur) != 1 {
            break;
        }
        let x = lines[adj.out_idx[&cur][0]].0.to;
        if nodes.contains(&x) {
            break;
        }
        if adj.indeg(x) == 1 {
            nodes.push(x);
            cur = x;
        } else {
            break; // x 入度≥2 = 汇合点，不含
        }
    }
    nodes
}

/// 链覆盖：链内真实入边行 count 的 max（根行 frm<0 不在 in_idx, 自动跳过）。
fn chain_cov(adj: &Adj, lines: &[(NodeLine, &str)], nodes: &[i64]) -> i64 {
    nodes
        .iter()
        .filter_map(|v| adj.in_idx.get(v))
        .flatten()
        .map(|&i| lines[i].0.count)
        .max()
        .unwrap_or(0)
}

/// 节点 kmer 的定长拷贝（非 24 字节的畸形行跳过, 返回 None）。
fn kmer24(adj: &Adj, lines: &[(NodeLine, &str)], v: i64) -> Option<[u8; K]> {
    let k = adj.kmer(lines, v)?;
    if k.len() != K {
        return None;
    }
    let mut out = [0u8; K];
    out.copy_from_slice(k);
    Some(out)
}

/// read 的 tip 命中集（(out-tip 序号集, in-tip 序号集)）。
type TipHits = (HashSet<u32>, HashSet<u32>);

/// read name → pair key：去 `>` 前缀与尾部 `/1`/`/2`（Trinity both.fa 惯例;
/// 新 Illumina 格式 pair 标记在第二 token, short name 天然同 key）。
fn pair_key(name: &str) -> &str {
    let n = name.strip_prefix('>').unwrap_or(name);
    n.strip_suffix("/1").or(n.strip_suffix("/2")).unwrap_or(n)
}

/// 单条 read 命中的 tip 集合（(out 集, in 集)）。非链特异模式正/反向都试
/// （镜像 QuantifyGraph 对 25-mer 的双向计数）。
fn read_tip_hits(
    seq: &[u8],
    index: &HashMap<[u8; K], Vec<(bool, u32)>>,
    strand: bool,
) -> TipHits {
    let mut out_hits: HashSet<u32> = HashSet::new();
    let mut in_hits: HashSet<u32> = HashSet::new();
    let mut scan = |s: &[u8]| {
        if s.len() < K {
            return;
        }
        for w in s.windows(K) {
            let mut kmer = [0u8; K];
            kmer.copy_from_slice(w);
            if let Some(marks) = index.get(&kmer) {
                for &(is_out, ti) in marks {
                    if is_out {
                        out_hits.insert(ti);
                    } else {
                        in_hits.insert(ti);
                    }
                }
            }
        }
    };
    scan(seq);
    if !strand {
        scan(&revcomp(seq));
    }
    (out_hits, in_hits)
}

/// 缝合计划（gap 确定后即定案; 修正路径至多改一侧锚点 kmer）。
struct BridgePlan {
    /// 缺口碱基数（1..=K-MIN_INTERSECTION）。
    gap: usize,
    /// 重叠长度 = K - gap。
    overlap: usize,
    /// 右权威修正：改写 a 节点 kmer（错配唯一且贴 a 末位）。
    fix_a: Option<[u8; K]>,
    /// 左权威修正：改写 b 节点 kmer（错配唯一且贴 b 首位）。
    fix_b: Option<[u8; K]>,
    /// 是否错配修正缝合（统计用）。
    corrected: bool,
}

/// SPAdes `ProcessPair`：overlap 自 K-1 递减到 MIN_INTERSECTION 逐一试。
/// 汉明 > bound → 下一 gap; 低复杂度 → 整对放弃; 否则定案（成败不继续）。
fn plan_pair(
    ka: &[u8; K],
    kb: &[u8; K],
    cov_a: i64,
    cov_b: i64,
    stats: &mut GapCloseStats,
) -> Option<BridgePlan> {
    for gap in 1..=(K - MIN_INTERSECTION) {
        let overlap = K - gap;
        let s1 = &ka[gap..]; // out-tip 末端 kmer 后缀
        let s2 = &kb[..overlap]; // in-tip 首端 kmer 前缀
        // 限定汉明（超过 bound 提前退出, SPAdes LimitedHammingDistance）
        let mut diff: Vec<usize> = Vec::new();
        for (i, (x, y)) in s1.iter().zip(s2.iter()).enumerate() {
            if x != y {
                diff.push(i);
                if diff.len() > HAMMING_DIST_BOUND {
                    break;
                }
            }
        }
        if diff.len() > HAMMING_DIST_BOUND {
            continue;
        }
        // 低复杂度检查（SPAdes L341-352: ratio 随 gap 自 0.8 线性升到 1.0,
        // 单碱基计数严格 > ratio×overlap → return false 整对放弃）
        let mut counts = [0usize; 256];
        for &c in s1 {
            counts[c as usize] += 1;
        }
        let curm = *counts.iter().max().unwrap_or(&0);
        let ratio = 0.8 + 0.2 * (gap - 1) as f64 / (K - MIN_INTERSECTION - 1) as f64;
        if (curm as f64) > ratio * overlap as f64 {
            stats.rejected_low_complexity += 1;
            return None;
        }
        if diff.is_empty() {
            return Some(BridgePlan {
                gap,
                overlap,
                fix_a: None,
                fix_b: None,
                corrected: false,
            });
        }
        // 错配修正（HandlePositiveHammingDistanceCase 简化语义: 以覆盖高的
        // 一侧为权威; tie 右侧优先 ≈ SPAdes 先试 CanCorrectLeft）。未压缩图
        // 节点不可分割, 只有贴端单位置可修: 右权威=错配在 overlap 末位
        // （改 a 的 kmer 末位, 不动与前驱兼容的前 23 位）; 左权威=错配在
        // overlap 首位（改 b 的 kmer 首位, 不动与后继兼容的后 23 位）。
        let right_ok = diff.len() == 1 && diff[0] == overlap - 1;
        let left_ok = diff.len() == 1 && diff[0] == 0;
        let right_first = cov_b >= cov_a;
        let order = if right_first {
            [true, false]
        } else {
            [false, true]
        };
        for right in order {
            if right && right_ok {
                let mut ka2 = *ka;
                ka2[K - 1] = kb[overlap - 1];
                return Some(BridgePlan {
                    gap,
                    overlap,
                    fix_a: Some(ka2),
                    fix_b: None,
                    corrected: true,
                });
            }
            if !right && left_ok {
                let mut kb2 = *kb;
                kb2[0] = ka[gap];
                return Some(BridgePlan {
                    gap,
                    overlap,
                    fix_a: None,
                    fix_b: Some(kb2),
                    corrected: true,
                });
            }
        }
        stats.rejected_hamming += 1;
        return None;
    }
    None
}

/// 改写节点 v 的全部定义行（to==v）的 kmer 列（第 4 列, 尾部列原样保留）。
fn rewrite_node_kmer(lines: &mut [(NodeLine, &str)], v: i64, new_kmer: &[u8; K]) {
    let kmer_str = std::str::from_utf8(new_kmer).unwrap_or("");
    for (l, raw) in lines.iter_mut() {
        if l.to != v {
            continue;
        }
        let tail: Vec<&str> = raw.split('\t').skip(4).collect();
        let mut text = format!("{}\t{}\t{}\t{kmer_str}", l.to, l.frm, l.count);
        for t in tail {
            text.push('\t');
            text.push_str(t);
        }
        l.rewrite = Some(text);
    }
}

/// 追加一条桥接行（格式镜像 QuantifyGraph 数据行: `to frm count kmer\t1`）。
fn push_bridge_line(lines: &mut Vec<(NodeLine, &str)>, to: i64, frm: i64, cov: i64, kmer: &[u8; K]) {
    let kmer_str = std::str::from_utf8(kmer).unwrap_or("");
    let text = format!("{to}\t{frm}\t{cov}\t{kmer_str}\t1");
    lines.push((
        NodeLine {
            to,
            frm,
            count: cov,
            alive: true,
            rewrite: Some(text),
        },
        "",
    ));
}

/// 应用缝合：改写锚点 kmer（如有）+ 追加桥接行 + 更新度数/kmer 占用。
/// 任一 kmer 撞组件既有集合或互撞 → 不应用, 返回 false（保守拒缝）。
#[allow(clippy::too_many_arguments)] // 缝合所需上下文一次带齐
fn apply_bridge(
    lines: &mut Vec<(NodeLine, &str)>,
    kmer_used: &mut HashSet<[u8; K]>,
    outdeg: &mut HashMap<i64, usize>,
    indeg: &mut HashMap<i64, usize>,
    next_id: &mut i64,
    a: i64,
    b: i64,
    ka: [u8; K],
    kb: [u8; K],
    plan: &BridgePlan,
    cov: i64,
    stats: &mut GapCloseStats,
) -> bool {
    let ka_eff = plan.fix_a.unwrap_or(ka);
    let kb_eff = plan.fix_b.unwrap_or(kb);
    // 桥接中间节点 kmer（机械漂移: 第 i 节点 = ka_eff[i..24] ++
    // kb_eff[overlap..overlap+i]; 完美时 i=gap 恰为 kb_eff, 每步都是合法
    // 1bp 漂移——修正路径改写锚点后同样成立）
    let mut mids: Vec<[u8; K]> = Vec::new();
    for i in 1..plan.gap {
        let mut kmer = [0u8; K];
        kmer[..K - i].copy_from_slice(&ka_eff[i..]);
        kmer[K - i..].copy_from_slice(&kb_eff[plan.overlap..plan.overlap + i]);
        mids.push(kmer);
    }
    // kmer 查重（防撞并既有节点; 周期序列互撞也拒）
    fn kmer_fresh(kmer: [u8; K], used: &HashSet<[u8; K]>, fresh: &mut HashSet<[u8; K]>) -> bool {
        !used.contains(&kmer) && fresh.insert(kmer)
    }
    let mut fresh: HashSet<[u8; K]> = HashSet::new();
    if let Some(k2) = plan.fix_a {
        if k2 != ka && !kmer_fresh(k2, kmer_used, &mut fresh) {
            stats.rejected_kmer_clash += 1;
            return false;
        }
    }
    if let Some(k2) = plan.fix_b {
        if k2 != kb && !kmer_fresh(k2, kmer_used, &mut fresh) {
            stats.rejected_kmer_clash += 1;
            return false;
        }
    }
    for m in &mids {
        if !kmer_fresh(*m, kmer_used, &mut fresh) {
            stats.rejected_kmer_clash += 1;
            return false;
        }
    }
    // 应用改写与桥接行
    if let Some(k2) = plan.fix_a {
        rewrite_node_kmer(lines, a, &k2);
    }
    if let Some(k2) = plan.fix_b {
        rewrite_node_kmer(lines, b, &k2);
    }
    let mut prev = a;
    for m in &mids {
        let id = *next_id;
        *next_id += 1;
        push_bridge_line(lines, id, prev, cov, m);
        *outdeg.entry(prev).or_insert(0) += 1;
        *indeg.entry(id).or_insert(0) += 1;
        kmer_used.insert(*m);
        prev = id;
    }
    push_bridge_line(lines, b, prev, cov, &kb_eff);
    *outdeg.entry(prev).or_insert(0) += 1;
    *indeg.entry(b).or_insert(0) += 1;
    kmer_used.insert(ka_eff);
    kmer_used.insert(kb_eff);
    stats.bridge_nodes += mids.len();
    stats.gaps_closed += 1;
    if plan.corrected {
        stats.closed_corrected += 1;
    } else {
        stats.closed_perfect += 1;
    }
    true
}

/// 单组件处理：tip 收集 → PE 配对索引 → 按权重降序缝合。
fn process_block(
    block: &mut crate::graph_clean::Block<'_>,
    reads: &[crate::dna_vector::DnaSeq],
    strand: bool,
    stats: &mut GapCloseStats,
) {
    let lines = &block.lines;
    let adj = build_adj(lines);

    // ---- tip 收集（SPAdes PrepareTipMap 的链式简化）----
    let mut out_anchors: Vec<i64> = adj
        .in_idx
        .keys()
        .copied()
        .filter(|&v| adj.outdeg(v) == 0)
        .collect();
    out_anchors.sort_unstable();
    let mut in_anchors: Vec<i64> = adj
        .out_idx
        .keys()
        .copied()
        .filter(|&v| adj.indeg(v) == 0)
        .collect();
    in_anchors.sort_unstable();

    let mut out_tips: Vec<Tip> = Vec::new();
    let mut in_tips: Vec<Tip> = Vec::new();
    let mut index: HashMap<[u8; K], Vec<(bool, u32)>> = HashMap::new();
    for &v in &out_anchors {
        let chain = collect_out_chain(&adj, lines, v);
        let cov = chain_cov(&adj, lines, &chain);
        // 末端区 = 链尾 ≤177 节点的 24mer
        let region = &chain[chain.len().saturating_sub(TIP_REGION_NODES)..];
        let ti = out_tips.len() as u32;
        for &u in region {
            if let Some(kmer) = kmer24(&adj, lines, u) {
                index.entry(kmer).or_default().push((true, ti));
            }
        }
        out_tips.push(Tip {
            anchor: v,
            cov,
            nodes: chain.into_iter().collect(),
        });
    }
    for &s in &in_anchors {
        let chain = collect_in_chain(&adj, lines, s);
        let cov = chain_cov(&adj, lines, &chain);
        // 首端区 = 链首 ≤177 节点的 24mer
        let region = &chain[..chain.len().min(TIP_REGION_NODES)];
        let ti = in_tips.len() as u32;
        for &u in region {
            if let Some(kmer) = kmer24(&adj, lines, u) {
                index.entry(kmer).or_default().push((false, ti));
            }
        }
        in_tips.push(Tip {
            anchor: s,
            cov,
            nodes: chain.into_iter().collect(),
        });
    }
    stats.out_tips += out_tips.len();
    stats.in_tips += in_tips.len();
    if out_tips.is_empty() || in_tips.is_empty() {
        return;
    }

    // ---- PE 配对索引（SPAdes ProcessPairedRead 的逐组件简化）----
    let mut pairs: HashMap<&str, Vec<TipHits>> = HashMap::new();
    for r in reads {
        let (oh, ih) = read_tip_hits(&r.seq, &index, strand);
        if oh.is_empty() && ih.is_empty() {
            continue;
        }
        stats.tip_reads += 1;
        pairs.entry(pair_key(&r.name)).or_default().push((oh, ih));
    }
    // 同一 pair 的两个 mate 分落 out-tip A / in-tip B（任一 mate 在哪侧均可;
    // 同一条 read 不自配, SPAdes path1×path2 语义）
    let mut votes: HashMap<(u32, u32), u32> = HashMap::new();
    for mates in pairs.values() {
        for (i, (oh_i, _)) in mates.iter().enumerate() {
            for (j, (_, ih_j)) in mates.iter().enumerate() {
                if i == j {
                    continue;
                }
                for &a in oh_i {
                    for &b in ih_j {
                        *votes.entry((a, b)).or_insert(0) += 1;
                    }
                }
            }
        }
    }
    stats.pair_votes += votes.values().map(|&w| w as usize).sum::<usize>();

    // ---- 候选收集（SPAdes CollectCandidates: 权重阈值 + 降序, tie 锚点 id 升序）----
    let mut by_a: HashMap<u32, Vec<(u32, u32)>> = HashMap::new(); // out → [(in, weight)]
    for ((a, b), w) in votes {
        if w < WEIGHT_THRESHOLD {
            continue;
        }
        // 防自连: 两链共节点（近似 avoid_rc_connections 的 first==second/
        // conjugate 防护——无链特异图没有显式共轭表示）
        if out_tips[a as usize]
            .nodes
            .iter()
            .any(|u| in_tips[b as usize].nodes.contains(u))
        {
            continue;
        }
        by_a.entry(a).or_default().push((b, w));
        stats.candidate_pairs += 1;
    }
    if by_a.is_empty() {
        return;
    }
    for cands in by_a.values_mut() {
        cands.sort_by(|&(b1, w1), &(b2, w2)| {
            w2.cmp(&w1)
                .then(in_tips[b1 as usize].anchor.cmp(&in_tips[b2 as usize].anchor))
        });
    }
    let mut a_order: Vec<u32> = by_a.keys().copied().collect();
    a_order.sort_by_key(|&a| out_tips[a as usize].anchor);

    // ---- 缝合（SPAdes CloseShortGaps: 复查死端/源头, 每 tip 首个可行候选即 break）----
    let lines = &mut block.lines;
    let mut outdeg: HashMap<i64, usize> = HashMap::new();
    let mut indeg: HashMap<i64, usize> = HashMap::new();
    let mut kmer_used: HashSet<[u8; K]> = HashSet::new();
    let mut max_id: i64 = 0;
    for (l, raw) in lines.iter() {
        outdeg.entry(l.frm).or_insert(0);
        if l.frm >= 0 {
            *outdeg.entry(l.frm).or_insert(0) += 1;
            *indeg.entry(l.to).or_insert(0) += 1;
            max_id = max_id.max(l.frm);
        }
        max_id = max_id.max(l.to);
        if let Some(k) = raw.split('\t').nth(3) {
            if k.len() == K {
                let mut kmer = [0u8; K];
                kmer.copy_from_slice(k.as_bytes());
                kmer_used.insert(kmer);
            }
        }
    }
    let mut next_id = max_id + 1;
    let deg = |m: &HashMap<i64, usize>, v: i64| m.get(&v).copied().unwrap_or(0);

    for a in a_order {
        let a_tip = &out_tips[a as usize];
        for &(b, _w) in &by_a[&a] {
            let b_tip = &in_tips[b as usize];
            // SPAdes L434-437 复查: 先前的缝合可能已把端点接走
            if deg(&outdeg, a_tip.anchor) != 0 || deg(&indeg, b_tip.anchor) != 0 {
                continue;
            }
            let (Some(ka), Some(kb)) = (
                kmer24(&adj, lines.as_slice(), a_tip.anchor),
                kmer24(&adj, lines.as_slice(), b_tip.anchor),
            ) else {
                continue;
            };
            stats.pairs_checked += 1;
            let Some(plan) = plan_pair(&ka, &kb, a_tip.cov, b_tip.cov, stats) else {
                continue;
            };
            let cov = a_tip.cov.min(b_tip.cov).max(1);
            if apply_bridge(
                lines,
                &mut kmer_used,
                &mut outdeg,
                &mut indeg,
                &mut next_id,
                a_tip.anchor,
                b_tip.anchor,
                ka,
                kb,
                &plan,
                cov,
                stats,
            ) {
                break; // 每 out-tip 只接受第一个可行候选
            }
        }
    }
}

/// 对 graph.out 文本做 SPAdes 式 tip 对 tip 缝合, 返回 (新文本, 统计)。
/// `reads_fasta` = 该组件 reads.tmp 字节（`>name pct\nseq\n` FASTA）。
/// 不做任何修改的组件（无候选）逐字节透传。
pub fn close_gaps_spades(text: &str, reads_fasta: &[u8], strand: bool) -> (String, GapCloseStats) {
    let mut stats = GapCloseStats::default();
    let reads = read_fasta_bytes(reads_fasta, true);
    let mut blocks = split_blocks(text);
    for b in &mut blocks {
        if b.lines.is_empty() {
            continue;
        }
        stats.components += 1;
        process_block(b, &reads, strand, &mut stats);
    }
    let mut out = String::with_capacity(text.len());
    for b in &blocks {
        for h in &b.head {
            out.push_str(h);
            out.push('\n');
        }
        for (l, raw) in &b.lines {
            // gap closer 只加不改删: alive 恒 true, rewrite 承载新行/改写行
            out.push_str(l.rewrite.as_deref().unwrap_or(raw));
            out.push('\n');
        }
    }
    (out, stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 后继 kmer = 前 23 碱基右移 + 新碱基（de Bruijn 漂移 1bp）。
    fn sh(k: &str, c: char) -> String {
        format!("{}{}", &k[1..], c)
    }

    fn mkline(to: i64, frm: i64, count: i64, kmer: &str) -> String {
        format!("{to}\t{frm}\t{count}\t{kmer}\t1")
    }

    fn graph(lines: &[String]) -> String {
        format!("Component 0\n{}\n", lines.join("\n"))
    }

    /// 一对 PE reads（reads.tmp 形态: `>name/1 pct\nseq\n`）。
    fn pe_pair(name: &str, s1: &str, s2: &str) -> String {
        format!(">{name}/1 96%\n{s1}\n>{name}/2 96%\n{s2}\n")
    }

    /// 链 X 的三节点 kmer（漂移合法）: k1 → k2 → k3。
    fn seed_kmers() -> (String, String, String) {
        let k1 = "TTGCACGATCGTACGGCATGCTAG".to_string();
        let k2 = sh(&k1, 'A');
        let k3 = sh(&k2, 'C');
        (k1, k2, k3)
    }

    /// 基础断裂图: 链 X 根 1→2→3（3 死端, out-tip）, 链 Y 根 50→51
    /// （50 源头, in-tip）。`kb` = 节点 50 的 kmer（按期望 overlap 构造）。
    fn broken_graph(kb: &str, cov_x: i64, cov_y: i64) -> String {
        let (k1, k2, k3) = seed_kmers();
        let k51 = sh(kb, 'G');
        graph(&[
            mkline(1, -1, 0, &k1),
            mkline(2, 1, cov_x, &k2),
            mkline(3, 2, cov_x, &k3),
            mkline(50, -1, 0, kb),
            mkline(51, 50, cov_y, &k51),
        ])
    }

    /// n 对 PE: mate1 命中链 X 末端区（含 k3 的 read）, mate2 命中链 Y 首端区
    /// （含 k50=kb 的 read）。
    fn pe_reads(k3: &str, kb: &str, n_pairs: usize) -> String {
        let mut s = String::new();
        for i in 0..n_pairs {
            s.push_str(&pe_pair(&format!("r{i}"), k3, kb));
        }
        s
    }

    /// 数据行的 (to, frm, count, kmer) 四元组集合（顺序无关断言用）。
    fn rows(text: &str) -> Vec<(i64, i64, i64, String)> {
        text.lines()
            .filter_map(|l| {
                let t: Vec<&str> = l.split('\t').collect();
                if t.len() >= 4 {
                    Some((
                        t[0].parse().ok()?,
                        t[1].parse().ok()?,
                        t[2].parse().ok()?,
                        t[3].to_string(),
                    ))
                } else {
                    None
                }
            })
            .collect()
    }

    #[test]
    fn pair_key_strips_suffix() {
        assert_eq!(pair_key(">61DF:1:100:1:1/1"), "61DF:1:100:1:1");
        assert_eq!(pair_key(">61DF:1:100:1:1/2"), "61DF:1:100:1:1");
        assert_eq!(pair_key(">M001:2:3:4"), "M001:2:3:4"); // 新格式天然同 key
    }

    #[test]
    fn perfect_overlap_bridged_with_two_pairs() {
        // overlap=20（gap=4）: kb 前 20 = k3 后 20。
        // 注意小测试图里两链的另一端也是死端/源头（51 是 out-tip, 1 是
        // in-tip）, 反向对 (51,1) 也会得票进候选, 但 k51 后缀与 k1 前缀无
        // 重叠 → 只缝 (3,50)。
        let (_, _, k3) = seed_kmers();
        let kb = format!("{}GATT", &k3[4..]);
        let g = broken_graph(&kb, 10, 8);
        let reads = pe_reads(&k3, &kb, 2);
        let (out, st) = close_gaps_spades(&g, reads.as_bytes(), false);
        assert_eq!(st.out_tips, 2, "3 与 51 都是死端");
        assert_eq!(st.in_tips, 2, "1 与 50 都是源头");
        assert_eq!(st.tip_reads, 4, "两条 mate 各命中一次 × 2 对");
        assert_eq!(st.pair_votes, 4, "正向 (3,50) 与反向 (51,1) 各 2 票");
        assert_eq!(st.candidate_pairs, 2);
        assert_eq!(st.pairs_checked, 2);
        assert_eq!(st.gaps_closed, 1);
        assert_eq!(st.closed_perfect, 1);
        assert_eq!(st.bridge_nodes, 3, "gap=4 → 3 个中间节点");
        // 桥接行: 52←3, 53←52, 54←53, 50←54; 覆盖 = min(10,8)=8
        let r = rows(&out);
        for (to, frm) in [(52, 3), (53, 52), (54, 53), (50, 54)] {
            let row = r
                .iter()
                .find(|(t, f, _, _)| *t == to && *f == frm)
                .unwrap_or_else(|| panic!("缺桥接行 {to}←{frm}"));
            assert_eq!(row.2, 8, "桥接覆盖 = min(两链覆盖)");
        }
        // 漂移合法: 52 的 kmer = k3[1..24] + kb[20]
        let k52 = &r.iter().find(|(t, f, _, _)| *t == 52 && *f == 3).unwrap().3;
        assert_eq!(*k52, format!("{}{}", &k3[1..], &kb[20..21]));
        // b=50 既有根行保留, 新增真实入边
        assert!(r.iter().any(|(t, f, _, _)| *t == 50 && *f == -1));
        // 原行不丢: 总行数 = 原 5 + 新 4
        assert_eq!(r.len(), 9);
    }

    #[test]
    fn single_pair_below_threshold_not_bridged() {
        let (_, _, k3) = seed_kmers();
        let kb = format!("{}GATT", &k3[4..]);
        let g = broken_graph(&kb, 10, 8);
        let reads = pe_reads(&k3, &kb, 1); // 权重 1 < 2.0
        let (out, st) = close_gaps_spades(&g, reads.as_bytes(), false);
        assert_eq!(st.gaps_closed, 0);
        assert_eq!(st.candidate_pairs, 0);
        assert_eq!(st.pair_votes, 2, "正反两方向各 1 票, 均 < 2.0");
        assert_eq!(out, g, "未达阈值 → 逐字节不变");
    }

    #[test]
    fn overlap_below_min_intersection_not_bridged() {
        // 最大共享 9bp（< min_intersection 10）: kb 前 9 = k3 后 9, 其余独特
        let (_, _, k3) = seed_kmers();
        let kb = format!("{}GATTCCGGAATTCCA", &k3[15..]);
        assert_eq!(kb.len(), 24);
        let g = broken_graph(&kb, 10, 8);
        let reads = pe_reads(&k3, &kb, 3);
        let (out, st) = close_gaps_spades(&g, reads.as_bytes(), false);
        assert_eq!(st.gaps_closed, 0);
        assert_eq!(st.pairs_checked, 2, "正反两对都进了 ProcessPair 但 overlap 全不匹配");
        assert_eq!(out, g);
    }

    #[test]
    fn min_intersection_boundary_bridged() {
        // 恰好 overlap=10（gap=14）: kb 前 10 = k3 后 10
        let (_, _, k3) = seed_kmers();
        let kb = format!("{}GATTCCGGAATTCC", &k3[14..]);
        assert_eq!(kb.len(), 24);
        let g = broken_graph(&kb, 10, 8);
        let reads = pe_reads(&k3, &kb, 2);
        let (_out, st) = close_gaps_spades(&g, reads.as_bytes(), false);
        assert_eq!(st.gaps_closed, 1);
        assert_eq!(st.bridge_nodes, 13, "gap=14 → 13 个中间节点");
    }

    #[test]
    fn low_complexity_overlap_rejected() {
        // overlap 区 poly-A（k3 末 23 位全 A）: gap=1 时 hamming=0 但
        // curm=23 > 0.8×23 → 整对放弃（SPAdes return false 语义）。
        // 链节点漂移合法: TCG+A*21 → CG+A*22 → G+A*23。
        let k1 = format!("TCG{}", "A".repeat(21));
        let k2 = sh(&k1, 'A');
        let k3 = sh(&k2, 'A');
        assert_eq!(&k3[1..], "A".repeat(23));
        let kb = format!("{}G", "A".repeat(23)); // 前 23 = k3 后 23（全 A）
        let k51 = sh(&kb, 'G');
        let g = graph(&[
            mkline(1, -1, 0, &k1),
            mkline(2, 1, 10, &k2),
            mkline(3, 2, 10, &k3),
            mkline(50, -1, 0, &kb),
            mkline(51, 50, 8, &k51),
        ]);
        let reads = pe_reads(&k3, &kb, 2);
        let (out, st) = close_gaps_spades(&g, reads.as_bytes(), false);
        assert_eq!(st.gaps_closed, 0);
        assert_eq!(st.rejected_low_complexity, 1);
        assert_eq!(out, g);
    }

    #[test]
    fn mismatch_at_a_end_corrected_right_authority() {
        // overlap=20（gap=4）, 唯一错配贴 a 末位（diff=[19]）: kb[19] = flip(k3[23])。
        // 覆盖 10 vs 8（左权威优先）但 left_ok 不成立（错配不在首位）→
        // 右权威: 改 a(3) kmer 末位。
        let (_, _, k3) = seed_kmers();
        let mut kb: String = k3[4..].to_string();
        let flip = if kb.as_bytes()[19] == b'A' { 'C' } else { 'A' };
        kb.replace_range(19..20, &flip.to_string());
        kb.push_str("GATT");
        assert_eq!(kb.len(), 24);
        let g = broken_graph(&kb, 10, 8);
        let reads = pe_reads(&k3, &kb, 2);
        let (out, st) = close_gaps_spades(&g, reads.as_bytes(), false);
        assert_eq!(st.gaps_closed, 1);
        assert_eq!(st.closed_corrected, 1);
        // a(3) 的定义行 kmer 末位被改写为 kb[19], 前 23 位不动（前驱兼容）
        let r = rows(&out);
        let kmer3 = &r.iter().find(|(t, _, _, _)| *t == 3).unwrap().3;
        assert_eq!(&kmer3[..23], &k3[..23]);
        assert_eq!(kmer3.as_bytes()[23], kb.as_bytes()[19]);
    }

    #[test]
    fn mismatch_at_b_start_corrected_left_authority() {
        // overlap=20（gap=4）, 唯一错配贴 b 首位（diff=[0]）且 a 覆盖更高
        // （左权威优先）: kb[0] = flip(k3[4]), kb[1..20] = k3[5..24]。
        let (_, _, k3) = seed_kmers();
        let mut kb: String = k3[4..].to_string();
        let flip = if kb.as_bytes()[0] == b'A' { 'C' } else { 'A' };
        kb.replace_range(0..1, &flip.to_string());
        kb.push_str("GATT");
        assert_eq!(kb.len(), 24);
        let g = broken_graph(&kb, 10, 8); // cov_x=10 > cov_y=8 → 左权威
        let reads = pe_reads(&k3, &kb, 2);
        let (out, st) = close_gaps_spades(&g, reads.as_bytes(), false);
        assert_eq!(st.gaps_closed, 1);
        assert_eq!(st.closed_corrected, 1);
        // b(50) 的根行 kmer 首位被改写为 k3[4], 后 23 位不动（后继兼容）
        let r = rows(&out);
        let kmer50 = &r
            .iter()
            .find(|(t, f, _, _)| *t == 50 && *f == -1)
            .unwrap()
            .3;
        assert_eq!(kmer50.as_bytes()[0], k3.as_bytes()[4]);
        assert_eq!(&kmer50[1..], &kb[1..]);
    }

    #[test]
    fn mismatch_in_middle_rejected() {
        // overlap=20, 唯一错配在中间（diff=[10]）→ 两侧窗口都不满足 → 拒
        let (_, _, k3) = seed_kmers();
        let mut kb: String = k3[4..].to_string();
        let flip = if kb.as_bytes()[10] == b'A' { 'C' } else { 'A' };
        kb.replace_range(10..11, &flip.to_string());
        kb.push_str("GATT");
        assert_eq!(kb.len(), 24);
        let g = broken_graph(&kb, 10, 8);
        let reads = pe_reads(&k3, &kb, 2);
        let (out, st) = close_gaps_spades(&g, reads.as_bytes(), false);
        assert_eq!(st.gaps_closed, 0);
        assert_eq!(st.rejected_hamming, 1);
        assert_eq!(out, g);
    }

    #[test]
    fn self_chain_pair_not_bridged() {
        // 单条孤立链同时是 out-tip 和 in-tip（两端）: 共节点 → 防自连
        let k1 = "TTGCACGATCGTACGGCATGCTAG".to_string();
        let k2 = sh(&k1, 'A');
        let g = graph(&[mkline(1, -1, 0, &k1), mkline(2, 1, 10, &k2)]);
        // mate1 命中 k2（out 末端区）, mate2 命中 k1（in 首端区）
        let reads = pe_reads(&k2, &k1, 3);
        let (out, st) = close_gaps_spades(&g, reads.as_bytes(), false);
        assert_eq!(st.out_tips, 1);
        assert_eq!(st.in_tips, 1);
        assert_eq!(st.candidate_pairs, 0, "同链共节点 → 候选被防自连过滤");
        assert_eq!(st.gaps_closed, 0);
        assert_eq!(out, g);
    }

    #[test]
    fn no_pe_evidence_byte_identical() {
        let kb = "AAAAAAAAAAAAAAAAAAAAAAAG".to_string();
        let g = broken_graph(&kb, 10, 8);
        let (out, st) = close_gaps_spades(&g, b"", false);
        assert_eq!(st.gaps_closed, 0);
        assert_eq!(out, g, "无 reads → 逐字节不变");
    }
}
