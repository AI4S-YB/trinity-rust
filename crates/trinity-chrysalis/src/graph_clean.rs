//! rnaSPAdes 式相对覆盖图清理（`--graph-clean-rnaspades`，默认关 = 逐字节不变）。
//!
//! 在 chrysalis QuantifyGraph 之后、butterfly 之前对单组件 `c*.graph.out`
//! （行: `to\tfrm\tcount\tkmer24\tflag`）做三类清理，镜像
//! `spades-main/src/projects/spades/configs/rna_mode.info` 的 simp 段条件
//! （判定结构见 `common/modules/simplification/{tip_clipper,erroneous_connection_remover,
//! relative_coverage_remover}.hpp`）：
//!
//! - **tip（死端链剪除）**，两组满足其一（rna_mode.info `simp.tc`）：
//!   1. `{mmm 3, tc_lb 4, rctc 0.5}`：链长 ≤4 节点，链覆盖 ≤ 0.5×(竞争分支最大
//!      覆盖+1)，且与竞争分支序列汉明差 ≤3（SPAdes `RelativeCoverageTipCondition`
//!      的 `le(cov, rctc*(max+1))` 与 `MismatchTipCondition` 的近似）；
//!   2. `{tc_lb 2, cb 1}`：链长 ≤2 节点且链覆盖 ≤1（`rctc 10000` 恒真，略）。
//! - **EC（错误连接删除）** `{ec_lb 30, icb 200, rcec_cb 1.0}`：两端都是
//!   junction（起点出度≥2 且终点入度≥2）的短链（内部节点 ≤2，ec_lb 30bp 的
//!   节点数近似），链覆盖 ≤200 且 ≤ 两端邻接分支平均覆盖+1
//!   （`RelativeCoverageECCondition`：`le(cov, ratio*(avg+1))`）。
//! - **rcc（相对低覆盖组件）** `coverage_gap 20`：从"覆盖 < 邻域最大覆盖/20"
//!   的节点 flood-fill 低覆盖子图；内部顶点 ≤10、边界顶点个个都有 ≥20 倍于
//!   内部基线的高覆盖外边 → 整片删（任务书简化版，省略 SPAdes 的长度/最长
//!   连接路径上界）。
//!
//! **映射取舍**（我们的图是未压缩 de Bruijn 图，每边=1 节点=24-mer 漂移 1bp）：
//! - SPAdes 的"边" ≈ 我们的死端/junction 间**链**（condensed edge 的对应物）；
//!   `tc_lb`/`ec_lb` 的 bp 阈值按节点数近似（任务书指定：4→4 节点、2→2 节点、
//!   30bp→内部 ≤2 节点）。
//! - `count` 列是 (frm 首碱基 + to 的 24-mer) 25-mer 的 read 计数 = **边覆盖**
//!   （SPAdes `g.coverage(e)` 的直接对应）。链覆盖取链内各边 max（保守）。
//! - 根行（frm<0）count 恒为 0（QuantifyGraph 不计数），只作顶点定义，不作
//!   覆盖依据；rcc 邻域/判定均只看真实边（frm≥0）。
//!
//! **安全性**：
//! - 每次删除以"行"为单位；同一轮内所有判定基于轮初图（同步语义，近似
//!   SPAdes 的 chunk 并行删除），链间天然不相交（内部节点度数 (1,1) 唯一）。
//! - 删除后修复悬空顶点：被存活边引用（作 frm）但自身定义行（to==v）全灭
//!   的节点，把其最高 count 的死行改写为根行（frm=-1, count=0，与
//!   QuantifyGraph 根行惯例一致）复活——保证 butterfly 建图不会遇到
//!   无 kmer 的空顶点。
//! - 兜底阀：组件数据行若被清空（蝴蝶组件池不容忍空图），整组件回滚。

use std::collections::{HashMap, HashSet};

/// rna_mode.info `simp.tc` 第一组：tc_lb 4 → 链长上限（节点数近似）。
const TIP_G1_MAX_CHAIN: usize = 4;
/// 第一组：mmm 3（与竞争分支的汉明差上限）。
const TIP_G1_MAX_MISMATCH: usize = 3;
/// 第一组：rctc 0.5（tip 覆盖 ≤ 0.5×(竞争最大覆盖+1)）。
const TIP_G1_RCTC: f64 = 0.5;
/// 第二组：tc_lb 2 → 链长上限 2 节点。
const TIP_G2_MAX_CHAIN: usize = 2;
/// 第二组：cb 1 → 链覆盖上限 1。
const TIP_G2_MAX_COV: i64 = 1;
/// rna_mode.info `simp.ec`：ec_lb 30 ≈ junction 间内部节点 ≤2（节点数近似）。
const EC_MAX_INTERNAL: usize = 2;
/// icb 200（迭代覆盖上界, 单轮 = 200）。
const EC_MAX_COV: i64 = 200;
/// rcec_cb 1.0（链覆盖 ≤ 1.0×(邻接平均覆盖+1)）。
const EC_RCEC_CB: f64 = 1.0;
/// rna_mode.info `simp.rcc`：coverage_gap 20。
const RCC_GAP: f64 = 20.0;
/// 低覆盖子图内部顶点数上限（任务书版; SPAdes 默认 max_vertex_cnt 30）。
const RCC_MAX_INNER: usize = 10;

/// 清理统计（跨组件累加，管线层汇总打印）。
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct CleanStats {
    /// 处理的组件数。
    pub components: usize,
    /// tip 第一组（mmm/rctc）命中链数。
    pub tip_g1_chains: usize,
    /// tip 第二组（短链低覆盖）命中链数。
    pub tip_g2_chains: usize,
    /// tip 删除的节点总数。
    pub tip_nodes: usize,
    /// EC 删除链数。
    pub ec_chains: usize,
    /// EC 删除的内部节点数。
    pub ec_nodes: usize,
    /// rcc 删除的低覆盖子图数。
    pub rcc_components: usize,
    /// rcc 删除的内部顶点总数。
    pub rcc_nodes: usize,
    /// 悬空修复（改根行复活）的节点数。
    pub repaired_dangling: usize,
    /// 兜底回滚（会被清空的组件）数。
    pub reverted_components: usize,
    /// 删除的数据行总数（修复前）。
    pub lines_removed: usize,
}

impl CleanStats {
    /// 逐字段累加（rayon 并行下由管线层合并）。
    pub fn add_assign(&mut self, o: &CleanStats) {
        self.components += o.components;
        self.tip_g1_chains += o.tip_g1_chains;
        self.tip_g2_chains += o.tip_g2_chains;
        self.tip_nodes += o.tip_nodes;
        self.ec_chains += o.ec_chains;
        self.ec_nodes += o.ec_nodes;
        self.rcc_components += o.rcc_components;
        self.rcc_nodes += o.rcc_nodes;
        self.repaired_dangling += o.repaired_dangling;
        self.reverted_components += o.reverted_components;
        self.lines_removed += o.lines_removed;
    }
}

/// 单行数据（`to frm count kmer flag...`）。`alive=false` 即被某轮删除。
struct NodeLine {
    to: i64,
    frm: i64, // <0 = 根行（源点挂接, 无真实前驱）
    count: i64,
    alive: bool,
    /// 悬空修复时改写为根行的替代文本；None = 输出原文。
    rewrite: Option<String>,
}

/// 一个组件块：头部行（`Component N` 等 token<4 的行，原样透传）+ 数据行。
struct Block<'a> {
    head: Vec<&'a str>,
    lines: Vec<(NodeLine, &'a str)>,
}

/// 把 graph.out 文本切成组件块（数据行 = ≥4 列且前三列可解析）。
fn split_blocks(text: &str) -> Vec<Block<'_>> {
    let mut blocks: Vec<Block> = Vec::new();
    let mut cur = Block {
        head: Vec::new(),
        lines: Vec::new(),
    };
    for line in text.lines() {
        let toks: Vec<&str> = line.split('\t').collect();
        let parsed = if toks.len() >= 4 {
            match (
                toks[0].parse::<i64>(),
                toks[1].parse::<i64>(),
                toks[2].parse::<i64>(),
            ) {
                (Ok(to), Ok(frm), Ok(count)) => Some(NodeLine {
                    to,
                    frm,
                    count,
                    alive: true,
                    rewrite: None,
                }),
                _ => None,
            }
        } else {
            None
        };
        match parsed {
            Some(nl) => cur.lines.push((nl, line)),
            None => {
                // 头部行：若当前块已有数据行则先封存, 新起一块
                if !cur.lines.is_empty() {
                    blocks.push(std::mem::replace(
                        &mut cur,
                        Block {
                            head: Vec::new(),
                            lines: Vec::new(),
                        },
                    ));
                }
                cur.head.push(line);
            }
        }
    }
    if !cur.head.is_empty() || !cur.lines.is_empty() {
        blocks.push(cur);
    }
    blocks
}

/// 单组件邻接索引（只看 alive 行）。
struct Adj {
    /// frm(≥0) → 行号列表（出边）。
    out_idx: HashMap<i64, Vec<usize>>,
    /// to → 行号列表（真实入边, 仅 frm≥0）。
    in_idx: HashMap<i64, Vec<usize>>,
    /// to → 该节点首个 alive 行号（取 kmer/自有行用; 同节点各行 kmer 相同）。
    first_line: HashMap<i64, usize>,
}

fn build_adj(lines: &[(NodeLine, &str)]) -> Adj {
    let mut adj = Adj {
        out_idx: HashMap::new(),
        in_idx: HashMap::new(),
        first_line: HashMap::new(),
    };
    for (i, (nl, _raw)) in lines.iter().enumerate() {
        if !nl.alive {
            continue;
        }
        if nl.frm >= 0 {
            adj.out_idx.entry(nl.frm).or_default().push(i);
            adj.in_idx.entry(nl.to).or_default().push(i);
        }
        adj.first_line.entry(nl.to).or_insert(i);
    }
    adj
}

impl Adj {
    fn outdeg(&self, v: i64) -> usize {
        self.out_idx.get(&v).map_or(0, |v| v.len())
    }
    fn indeg(&self, v: i64) -> usize {
        self.in_idx.get(&v).map_or(0, |v| v.len())
    }
    /// 节点 kmer 字节串（取自其首个 alive 行的第 4 列）。
    fn kmer<'a>(&self, lines: &'a [(NodeLine, &'a str)], v: i64) -> Option<&'a [u8]> {
        self.first_line.get(&v).map(|&i| {
            lines[i]
                .1
                .split('\t')
                .nth(3)
                .unwrap_or("")
                .as_bytes()
        })
    }
}

/// 链序列 = 首节点 kmer + 后续每节点 kmer 的末字符（漂移 1bp 的 de Bruijn 链）。
fn chain_seq(adj: &Adj, lines: &[(NodeLine, &str)], nodes: &[i64]) -> Vec<u8> {
    let mut seq = Vec::new();
    for (i, &v) in nodes.iter().enumerate() {
        let k = adj.kmer(lines, v).unwrap_or(b"");
        if i == 0 {
            seq.extend_from_slice(k);
        } else if let Some(&last) = k.last() {
            seq.push(last);
        }
    }
    seq
}

/// 汉明差（等长前缀逐位比较）。
fn hamming(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
}

// ---------------------------------------------------------------------------
// pass 1: tip（死端链剪除）
// ---------------------------------------------------------------------------

/// 死端链判定结果（kill 集 = 链节点自有行 + 源向链的连接行）。
struct TipHit {
    kill: Vec<usize>,
    group2: bool,
    nodes: usize,
}

/// 在单个组件上做 tip 判定；返回是否删了任何行。
fn tip_pass(lines: &mut [(NodeLine, &str)], stats: &mut CleanStats) -> bool {
    let adj = build_adj(lines);
    let mut hits: Vec<TipHit> = Vec::new();

    // ---- 汇向死端链（终点 outdeg==0）----
    let mut sinks: Vec<i64> = adj
        .in_idx
        .keys()
        .copied()
        .filter(|&v| adj.outdeg(v) == 0 && adj.indeg(v) == 1)
        .collect();
    sinks.sort_unstable();
    for t in sinks {
        // 从死端向回走：穿过 indeg==1&&outdeg==1 的内部节点, 停在出度≥2 的分支点 B
        let mut nodes_rev: Vec<i64> = vec![t]; // [vn, vn-1, ..., v1]
        let mut branch: Option<i64> = None;
        let mut cur = t;
        let mut ok = true;
        loop {
            let l = &lines[adj.in_idx[&cur][0]].0;
            let p = l.frm; // 唯一真实前驱（indeg==1 已保证）
            if p < 0 || nodes_rev.contains(&p) {
                ok = false; // 链上游无分支点（孤立路径/环）→ 非 tip
                break;
            }
            if adj.outdeg(p) >= 2 {
                branch = Some(p);
                break;
            }
            if adj.outdeg(p) == 1 && adj.indeg(p) == 1 {
                nodes_rev.push(p);
                if nodes_rev.len() > TIP_G1_MAX_CHAIN {
                    ok = false; // 超过两组条件的链长上限
                    break;
                }
                cur = p;
                continue;
            }
            ok = false; // p 是汇合点（outdeg==1 && indeg≥2）→ SPAdes TipCondition 不成立
            break;
        }
        if !ok {
            continue;
        }
        let b = branch.unwrap();
        // 链节点按图序 = nodes_rev 反转 [v1..vn]
        let nodes: Vec<i64> = nodes_rev.iter().rev().copied().collect();
        let v1 = nodes[0];
        let own: Vec<usize> = nodes.iter().map(|&v| adj.in_idx[&v][0]).collect();
        let chain_cov = own.iter().map(|&i| lines[i].0.count).max().unwrap_or(0);
        // 竞争分支 = B 的其它出边（SPAdes MaxCompetitorCoverage）
        let comps: Vec<usize> = adj.out_idx[&b]
            .iter()
            .copied()
            .filter(|&i| lines[i].0.to != v1)
            .collect();
        if comps.is_empty() {
            continue;
        }
        let max_comp = comps.iter().map(|&i| lines[i].0.count).max().unwrap_or(0);
        if let Some(group2) = tip_check(
            lines, &adj, &nodes, chain_cov, max_comp, &comps, true, b,
        ) {
            hits.push(TipHit {
                kill: own,
                group2,
                nodes: nodes.len(),
            });
        }
    }

    // ---- 源向死端链（起点无真实入边: indeg==0 且 outdeg==1）----
    let mut sources: Vec<i64> = adj
        .out_idx
        .keys()
        .copied()
        .filter(|&v| adj.indeg(v) == 0 && adj.outdeg(v) == 1)
        .collect();
    sources.sort_unstable();
    for s in sources {
        // 向前走：穿过 (1,1) 内部节点, 停在入度≥2 的汇合点 E
        let mut nodes: Vec<i64> = vec![s]; // [w1..wn]
        let mut end: Option<(i64, usize)> = None; // (E, 连接行行号)
        let mut cur = s;
        let mut ok = true;
        loop {
            let li = adj.out_idx[&cur][0];
            let x = lines[li].0.to;
            if nodes.contains(&x) {
                ok = false;
                break;
            }
            if adj.indeg(x) >= 2 {
                end = Some((x, li));
                break;
            }
            if adj.indeg(x) == 1 && adj.outdeg(x) == 1 {
                nodes.push(x);
                if nodes.len() > TIP_G1_MAX_CHAIN {
                    ok = false;
                    break;
                }
                cur = x;
                continue;
            }
            ok = false; // x 是另一死端或分支点 → 非 tip
            break;
        }
        if !ok {
            continue;
        }
        let (e, connector) = end.unwrap();
        let wn = *nodes.last().unwrap();
        // 链节点自有行 = 各自 to==v 的行（w1 是根行, 其余为真实入边行）
        let own: Vec<usize> = nodes
            .iter()
            .map(|&v| adj.first_line[&v]) // first_line 存的是该节点首个 alive 行号
            .collect();
        let chain_cov = own.iter().map(|&i| lines[i].0.count).max().unwrap_or(0);
        // 竞争分支 = E 的其它入边
        let comps: Vec<usize> = adj.in_idx[&e]
            .iter()
            .copied()
            .filter(|&i| lines[i].0.frm != wn)
            .collect();
        if comps.is_empty() {
            continue;
        }
        let max_comp = comps.iter().map(|&i| lines[i].0.count).max().unwrap_or(0);
        if let Some(group2) = tip_check(
            lines, &adj, &nodes, chain_cov, max_comp, &comps, false, e,
        ) {
            let mut kill = own;
            kill.push(connector); // 连接行 (E, wn) 一并删
            hits.push(TipHit {
                kill,
                group2,
                nodes: nodes.len(),
            });
        }
    }

    // ---- 应用（同步语义: 全部判定基于轮初图）----
    let mut changed = false;
    for h in hits {
        for i in h.kill {
            if lines[i].0.alive {
                lines[i].0.alive = false;
                stats.lines_removed += 1;
                changed = true;
            }
        }
        if h.group2 {
            stats.tip_g2_chains += 1;
        } else {
            stats.tip_g1_chains += 1;
        }
        stats.tip_nodes += h.nodes;
    }
    changed
}

/// tip 两组条件判定。`sink`：汇向链（竞争分支在分支点 B 的出侧）；
/// 否则源向链（竞争分支在汇合点 E 的入侧, 序列按 3' 端对齐比后缀）。
/// 返回 Some(是否第二组) / None（不删）。
#[allow(clippy::too_many_arguments)] // 两组条件所需的判定上下文一次带齐
fn tip_check(
    lines: &[(NodeLine, &str)],
    adj: &Adj,
    nodes: &[i64],
    chain_cov: i64,
    max_comp: i64,
    comps: &[usize],
    sink: bool,
    junction: i64,
) -> Option<bool> {
    let n = nodes.len();
    // 组 2 {tc_lb 2, cb 1}: 链 ≤2 节点且链覆盖 ≤1（rctc 10000 恒真, 略）
    if n <= TIP_G2_MAX_CHAIN && chain_cov <= TIP_G2_MAX_COV {
        return Some(true);
    }
    // 组 1 {mmm 3, tc_lb 4, rctc 0.5}
    if n > TIP_G1_MAX_CHAIN {
        return None;
    }
    if (chain_cov as f64) > TIP_G1_RCTC * (max_comp as f64 + 1.0) {
        return None;
    }
    // 竞争链: 从最高覆盖竞争者向外延伸（确定性 tie-break: count 降序, id 升序）
    let mut comp_nodes: Vec<i64> = {
        let mut cs: Vec<(i64, i64)> = comps
            .iter()
            .map(|&i| {
                let l = &lines[i].0;
                let id = if sink { l.to } else { l.frm };
                (-l.count, id)
            })
            .collect();
        cs.sort_unstable();
        vec![cs[0].1]
    };
    while comp_nodes.len() <= n {
        let cur = *comp_nodes.last().unwrap();
        let next_lines = if sink {
            adj.out_idx.get(&cur)
        } else {
            adj.in_idx.get(&cur)
        };
        let Some(cands) = next_lines else { break };
        // 选 (count 降序, 邻居 id 升序) 的最优后继/前驱
        let mut best: Option<(i64, i64)> = None;
        for &i in cands {
            let l = &lines[i].0;
            let id = if sink { l.to } else { l.frm };
            if id < 0 {
                continue; // 根行不作延伸
            }
            let key = (-l.count, id);
            if best.is_none_or(|b| key < b) {
                best = Some(key);
            }
        }
        let Some((_, nxt)) = best else { break };
        if nxt == junction || comp_nodes.contains(&nxt) || nodes.contains(&nxt) {
            break; // 回环/绕回链内 → 止步
        }
        comp_nodes.push(nxt);
    }
    // SPAdes MismatchTipCondition 要求竞争边严格更长（len < length(alt)）
    if comp_nodes.len() <= n {
        return None;
    }
    let chain = chain_seq(adj, lines, nodes);
    let comp = chain_seq(adj, lines, &comp_nodes);
    let diff = if sink {
        // 两侧均从分支点出发（5' 对齐）
        hamming(&chain, &comp[..chain.len()])
    } else {
        // 源向链: 两侧都止于汇合点 E（3' 对齐, 比后缀）
        hamming(&chain, &comp[comp.len() - chain.len()..])
    };
    if diff <= TIP_G1_MAX_MISMATCH {
        Some(false)
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// pass 2: EC（junction 间短链错误连接删除）
// ---------------------------------------------------------------------------

/// EC 候选: junction J1(出度≥2) → [≤2 个 (1,1) 内部节点] → junction J2(入度≥2)。
/// 命中则删链内全部行（内部节点随之消失; J1/J2 由其它边行保持定义）。
fn ec_pass(lines: &mut [(NodeLine, &str)], stats: &mut CleanStats) -> bool {
    let adj = build_adj(lines);
    // 候选首边: frm 为 junction 的 alive 真实边, 按 (frm, to) 排序保证确定性
    let mut starts: Vec<usize> = (0..lines.len())
        .filter(|&i| {
            let l = &lines[i].0;
            l.alive && l.frm >= 0 && adj.outdeg(l.frm) >= 2
        })
        .collect();
    starts.sort_by_key(|&i| (lines[i].0.frm, lines[i].0.to));

    let mut kills: Vec<(Vec<usize>, usize)> = Vec::new(); // (链行集, 内部节点数)
    for l0 in starts {
        let f0 = lines[l0].0.frm;
        if lines[l0].0.to == f0 {
            continue; // 自环不处理
        }
        // 沿唯一路径前行, 收集内部 (1,1) 节点, 止于入度≥2 的 J2
        let mut chain_edges: Vec<usize> = vec![l0];
        let mut internal: Vec<i64> = Vec::new();
        let mut cur = lines[l0].0.to;
        let j2 = loop {
            if cur == f0 || internal.contains(&cur) {
                break None; // 回环
            }
            let (din, dout) = (adj.indeg(cur), adj.outdeg(cur));
            if din >= 2 {
                break Some(cur);
            }
            if din == 1 && dout == 1 {
                internal.push(cur);
                if internal.len() > EC_MAX_INTERNAL {
                    break None; // 超长
                }
                let li = adj.out_idx[&cur][0];
                chain_edges.push(li);
                cur = lines[li].0.to;
                continue;
            }
            break None; // 死端(dout==0)或纯分支(din==1&&dout≥2) → 非 junction-junction
        };
        let Some(j2) = j2 else { continue };
        // 链覆盖 = 链内各边 max（保守）
        let cov = chain_edges
            .iter()
            .map(|&i| lines[i].0.count)
            .max()
            .unwrap_or(0);
        if cov > EC_MAX_COV {
            continue; // icb 200
        }
        // 邻接平均覆盖（RelativeCoverageECCondition：排除链自身两端边）
        let (first, last) = (chain_edges[0], *chain_edges.last().unwrap());
        let sum_out: i64 = adj.out_idx[&f0]
            .iter()
            .filter(|&&i| i != first)
            .map(|&i| lines[i].0.count)
            .sum();
        let sum_in: i64 = adj.in_idx[&j2]
            .iter()
            .filter(|&&i| i != last)
            .map(|&i| lines[i].0.count)
            .sum();
        let denom = adj.outdeg(f0) + adj.indeg(j2) - 2; // ≥2（两侧 junction）
        let avg = (sum_out + sum_in) as f64 / denom as f64;
        if (cov as f64) <= EC_RCEC_CB * (avg + 1.0) {
            kills.push((chain_edges, internal.len()));
        }
    }

    let mut changed = false;
    for (edges, n_internal) in kills {
        for i in edges {
            if lines[i].0.alive {
                lines[i].0.alive = false;
                stats.lines_removed += 1;
                changed = true;
            }
        }
        stats.ec_chains += 1;
        stats.ec_nodes += n_internal;
    }
    changed
}

// ---------------------------------------------------------------------------
// pass 3: rcc（相对低覆盖组件 flood-fill 删除）
// ---------------------------------------------------------------------------

/// 从低覆盖种子节点 flood-fill：邻居若不终止（无 ≥20× 基线的高覆盖外边）则
/// 并入组件; 终止者留作边界。内部顶点 ≤10、边界非空且各边界顶点都有
/// ≥20× 内部基线的高覆盖外边 → 删除组件内部与内部边（边界顶点及其外边保留）。
fn rcc_pass(lines: &mut [(NodeLine, &str)], stats: &mut CleanStats) -> bool {
    let adj = build_adj(lines);
    // 节点覆盖 = 关联真实边 count 的 max
    let mut ncov: HashMap<i64, i64> = HashMap::new();
    for (l, _) in lines.iter() {
        if !l.alive || l.frm < 0 {
            continue;
        }
        for v in [l.to, l.frm] {
            let e = ncov.entry(v).or_insert(l.count);
            if l.count > *e {
                *e = l.count;
            }
        }
    }
    // 节点级邻接（仅真实边）
    let neighbors = |v: i64| -> Vec<i64> {
        let mut ns = Vec::new();
        if let Some(os) = adj.out_idx.get(&v) {
            ns.extend(os.iter().map(|&i| lines[i].0.to));
        }
        if let Some(is) = adj.in_idx.get(&v) {
            ns.extend(is.iter().map(|&i| lines[i].0.frm));
        }
        ns
    };
    // 种子: 覆盖 < 邻域最大覆盖/20（严格 <）, 按 (覆盖, id) 升序处理保证确定性
    let mut order: Vec<(i64, i64)> = ncov.iter().map(|(&v, &c)| (c, v)).collect();
    order.sort_unstable();
    let mut tried: HashSet<i64> = HashSet::new();
    let mut kill: HashSet<usize> = HashSet::new();

    for (_, s) in order {
        if tried.contains(&s) {
            continue;
        }
        let c_s = ncov[&s] as f64;
        let max_nbr = neighbors(s)
            .iter()
            .map(|&u| ncov.get(&u).copied().unwrap_or(0))
            .max()
            .unwrap_or(0) as f64;
        if c_s * RCC_GAP >= max_nbr {
            continue; // 不满足 "覆盖 < 邻域/20"（严格 <）→ 非种子
        }
        // flood-fill
        let mut in_comp: HashSet<i64> = HashSet::from([s]);
        let mut inner: HashSet<i64> = HashSet::new();
        let mut border: HashSet<i64> = HashSet::new();
        let mut queue: std::collections::VecDeque<i64> = std::collections::VecDeque::from([s]);
        while let Some(u) = queue.pop_front() {
            // u 的基线 = 与**内部**顶点（不含边界）间真实边的 count max。
            // （边界顶点不算内部：边界边是"低覆盖侧挂接", 不应抬高基线——
            // 否则主路高覆盖边会把汇点拉进内部, 泡就永远删不掉。）
            let mut base: i64 = 0;
            let mut has_internal = false;
            for &i in adj
                .out_idx
                .get(&u)
                .into_iter()
                .flatten()
                .chain(adj.in_idx.get(&u).into_iter().flatten())
            {
                let l = &lines[i].0;
                let other = if l.frm == u { l.to } else { l.frm };
                if inner.contains(&other) && other != u {
                    has_internal = true;
                    base = base.max(l.count);
                }
            }
            // 终止判定（SPAdes IsTerminateVertex）: 有组件外边覆盖 > 20× 内边基线
            // （math::gr 严格 >）。无内边（种子首访）必不终止。
            let has_high_out = has_internal
                && adj
                    .out_idx
                    .get(&u)
                    .into_iter()
                    .flatten()
                    .chain(adj.in_idx.get(&u).into_iter().flatten())
                    .any(|&i| {
                        let l = &lines[i].0;
                        let other = if l.frm == u { l.to } else { l.frm };
                        !in_comp.contains(&other) && (l.count as f64) > RCC_GAP * base as f64
                    });
            if !has_high_out {
                // 并入内部, 吸收全部邻居
                inner.insert(u);
                for w in neighbors(u) {
                    if in_comp.insert(w) {
                        queue.push_back(w);
                    }
                }
            } else {
                border.insert(u);
            }
        }
        tried.extend(in_comp.iter().copied());
        // 检查: 内部顶点数 / 边界非空 / 边界四面高覆盖
        if inner.is_empty() || inner.len() > RCC_MAX_INNER || border.is_empty() {
            continue;
        }
        // 内部基线 = 两端都是内部顶点的真实边 count max
        let mut base: i64 = 0;
        for (l, _) in lines.iter() {
            if !l.alive || l.frm < 0 {
                continue;
            }
            if inner.contains(&l.to) && inner.contains(&l.frm) {
                base = base.max(l.count);
            }
        }
        let thr = RCC_GAP * base as f64;
        let all_high = border.iter().all(|&t| {
            adj.out_idx
                .get(&t)
                .into_iter()
                .flatten()
                .chain(adj.in_idx.get(&t).into_iter().flatten())
                .any(|&i| {
                    let l = &lines[i].0;
                    let other = if l.frm == t { l.to } else { l.frm };
                    !in_comp.contains(&other) && (l.count as f64) >= thr
                })
        });
        if !all_high {
            continue;
        }
        // 删除: 触及内部顶点的所有行（含其根行）。边界顶点及其外边保留;
        // 边界–边界边（如主路捷径）不是低覆盖泡的一部分, 不动。
        for (i, (l, _)) in lines.iter().enumerate() {
            if !l.alive {
                continue;
            }
            if inner.contains(&l.to) || inner.contains(&l.frm) {
                kill.insert(i);
            }
        }
        stats.rcc_components += 1;
        stats.rcc_nodes += inner.len();
    }

    let mut changed = false;
    for i in kill {
        if lines[i].0.alive {
            lines[i].0.alive = false;
            stats.lines_removed += 1;
            changed = true;
        }
    }
    changed
}

// ---------------------------------------------------------------------------
// 悬空修复
// ---------------------------------------------------------------------------

/// 删除后修复：被存活边引用（作 frm）但自身定义行（to==v）全灭的节点,
/// 把其最高 count 的死行改写为根行（`v\t-1\t0\tkmer\tflag...`）复活——
/// 避免 butterfly 建图时遇到无 kmer 的空顶点。改根行不引入新引用, 无需迭代。
fn repair_dangling(lines: &mut [(NodeLine, &str)], stats: &mut CleanStats) {
    let mut own_alive: HashMap<i64, usize> = HashMap::new();
    let mut ref_alive: HashSet<i64> = HashSet::new();
    for (l, _) in lines.iter() {
        if !l.alive {
            continue;
        }
        *own_alive.entry(l.to).or_insert(0) += 1;
        if l.frm >= 0 {
            ref_alive.insert(l.frm);
        }
    }
    let mut dangling: Vec<i64> = ref_alive
        .iter()
        .copied()
        .filter(|v| own_alive.get(v).copied().unwrap_or(0) == 0)
        .collect();
    dangling.sort_unstable();
    for v in dangling {
        // 最高 count（tie: 行号最小）的死行 → 改根行
        let mut best: Option<usize> = None;
        for (i, (l, _)) in lines.iter().enumerate() {
            if l.to == v && !l.alive && best.is_none_or(|b| lines[b].0.count < l.count) {
                best = Some(i);
            }
        }
        if let Some(i) = best {
            let kmer = lines[i]
                .1
                .split('\t')
                .nth(3)
                .unwrap_or("")
                .to_string();
            let tail: Vec<&str> = lines[i].1.split('\t').skip(4).collect();
            let mut text = format!("{v}\t-1\t0\t{kmer}");
            for t in tail {
                text.push('\t');
                text.push_str(t);
            }
            lines[i].0.alive = true;
            lines[i].0.frm = -1;
            lines[i].0.count = 0;
            lines[i].0.rewrite = Some(text);
            stats.repaired_dangling += 1;
        }
    }
}

// ---------------------------------------------------------------------------
// 入口
// ---------------------------------------------------------------------------

/// 对 graph.out 文本做 rnaSPAdes 式清理, 返回 (新文本, 统计)。
/// 不做任何文件的组件（无数据行）原样透传。
pub fn clean_graph_rnaspades(text: &str) -> (String, CleanStats) {
    let mut stats = CleanStats::default();
    let mut blocks = split_blocks(text);
    for b in &mut blocks {
        if b.lines.is_empty() {
            continue;
        }
        stats.components += 1;
        tip_pass(&mut b.lines, &mut stats);
        ec_pass(&mut b.lines, &mut stats);
        rcc_pass(&mut b.lines, &mut stats);
        repair_dangling(&mut b.lines, &mut stats);
        // 兜底阀: 组件被清空 → 整组件回滚（蝴蝶组件池不容忍空图）
        if b.lines.iter().all(|(l, _)| !l.alive) {
            for (l, _) in &mut b.lines {
                l.alive = true;
                l.rewrite = None;
            }
            stats.reverted_components += 1;
        }
    }
    let mut out = String::with_capacity(text.len());
    for b in &blocks {
        for h in &b.head {
            out.push_str(h);
            out.push('\n');
        }
        for (l, raw) in &b.lines {
            if l.alive {
                out.push_str(l.rewrite.as_deref().unwrap_or(raw));
                out.push('\n');
            }
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

    fn alive_nodes(out: &str) -> Vec<i64> {
        out.lines()
            .skip(1)
            .map(|l| l.split('\t').next().unwrap().parse::<i64>().unwrap())
            .collect()
    }

    const K0: &str = "ACGTGCACTGACCTGGCCTAAGTC"; // 分支点 B 的 kmer

    /// 基础图：根 1 → B(2) → 竞争支 3 → 4（高覆盖 cov_hi），
    /// B(2) → tip 10 → 11（cov_tip）。
    fn tip_fixture(cov_hi: i64, cov_tip: i64, tip_tail: (char, char)) -> String {
        let k1 = format!("T{}", &K0[..23]); // 根 kmer
        let k3 = sh(K0, 'A'); // 竞争者首节点
        let k4 = sh(&k3, 'G');
        let k10 = sh(K0, tip_tail.0);
        let k11 = sh(&k10, tip_tail.1);
        graph(&[
            mkline(1, -1, 0, &k1),
            mkline(2, 1, cov_hi, K0),
            mkline(3, 2, cov_hi, &k3),
            mkline(4, 3, cov_hi, &k4),
            mkline(10, 2, cov_tip, &k10),
            mkline(11, 10, cov_tip, &k11),
        ])
    }

    #[test]
    fn tip_group2_short_low_cov_removed() {
        // 2 节点链, 覆盖 1 ≤ 1 → 组 2 命中（注意组 1 不命中: 竞争链 2 节点不严格更长）
        let g = tip_fixture(50, 1, ('T', 'T'));
        let (out, st) = clean_graph_rnaspades(&g);
        assert_eq!(st.tip_g2_chains, 1);
        assert_eq!(st.tip_g1_chains, 0);
        assert_eq!(st.tip_nodes, 2);
        assert_eq!(alive_nodes(&out), vec![1, 2, 3, 4]);
    }

    #[test]
    fn tip_group2_cov2_kept() {
        // 链覆盖 3 > 1 → 组 2 失败; 组 1: 竞争链不长 → 失败; 覆盖 3 不满足
        // rcc 种子（20*3=60 > 邻域 50）→ 整体保留
        let g = tip_fixture(50, 3, ('T', 'T'));
        let (out, st) = clean_graph_rnaspades(&g);
        assert_eq!(st.tip_g2_chains, 0);
        assert_eq!(st.tip_g1_chains, 0);
        assert_eq!(st.rcc_components, 0);
        assert_eq!(out, g, "不满足任一条件 → 逐字节不变");
    }

    #[test]
    fn tip_group1_mmm_rctc_removed() {
        // 竞争支延长到 4 节点（严格更长）; tip 覆盖 5 ≤ 0.5*101; 汉明差 = 1 ≤ 3
        let k1 = format!("T{}", &K0[..23]);
        let k3 = sh(K0, 'A');
        let k4 = sh(&k3, 'G');
        let k5 = sh(&k4, 'C');
        let k6 = sh(&k5, 'A');
        let k10 = sh(K0, 'T'); // 与 k3 仅末位不同（1 diff）
        let k11 = sh(&k10, 'G'); // 尾字符与竞争支同位相同（不再新增 diff）
        let g = graph(&[
            mkline(1, -1, 0, &k1),
            mkline(2, 1, 100, K0),
            mkline(3, 2, 100, &k3),
            mkline(4, 3, 100, &k4),
            mkline(5, 4, 100, &k5),
            mkline(6, 5, 100, &k6),
            mkline(10, 2, 5, &k10),
            mkline(11, 10, 5, &k11),
        ]);
        let (out, st) = clean_graph_rnaspades(&g);
        assert_eq!(st.tip_g1_chains, 1);
        assert_eq!(st.tip_nodes, 2);
        assert_eq!(alive_nodes(&out), vec![1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn tip_group1_rctc_fail_kept() {
        // tip 覆盖 40 > 0.5*(50+1)=25.5 → rctc 失败（即便汉明差 0）
        let k1 = format!("T{}", &K0[..23]);
        let k3 = sh(K0, 'A');
        let k4 = sh(&k3, 'G');
        let k5 = sh(&k4, 'C');
        let k10 = sh(K0, 'A'); // 与 k3 相同
        let k11 = sh(&k10, 'G'); // 与 k4 相同 → 汉明差 0
        let g = graph(&[
            mkline(1, -1, 0, &k1),
            mkline(2, 1, 50, K0),
            mkline(3, 2, 50, &k3),
            mkline(4, 3, 50, &k4),
            mkline(5, 4, 50, &k5),
            mkline(10, 2, 40, &k10),
            mkline(11, 10, 40, &k11),
        ]);
        let (out, st) = clean_graph_rnaspades(&g);
        assert_eq!(st.tip_g1_chains, 0);
        assert_eq!(out, g);
    }

    #[test]
    fn tip_group1_mmm_fail_kept() {
        // 4 节点链, 每节点与竞争支同位皆不同 → 汉明差 ≥4 > 3 → 保留
        // （tip 覆盖 6: 20*6=120 !< 100, 保证不被后续 rcc 轮误删——rcc 未实现时亦无碍）
        let k1 = format!("T{}", &K0[..23]);
        let mut comp = vec![sh(K0, 'A')];
        for c in ['G', 'C', 'A', 'T'] {
            comp.push(sh(comp.last().unwrap(), c));
        }
        let mut tip = vec![sh(K0, 'T')];
        for c in ['C', 'T', 'G', 'C'] {
            tip.push(sh(tip.last().unwrap(), c));
        }
        let mut lines = vec![mkline(1, -1, 0, &k1), mkline(2, 1, 100, K0)];
        let mut prev = 2;
        for (i, k) in comp.iter().enumerate() {
            lines.push(mkline(100 + i as i64, prev, 100, k));
            prev = 100 + i as i64;
        }
        let mut prev = 2;
        for (i, k) in tip.iter().enumerate() {
            lines.push(mkline(200 + i as i64, prev, 6, k));
            prev = 200 + i as i64;
        }
        let g = graph(&lines);
        let (out, st) = clean_graph_rnaspades(&g);
        assert_eq!(st.tip_g1_chains, 0, "汉明差 4 > 3 不删");
        assert_eq!(out, g);
    }

    #[test]
    fn tip_chain_too_long_kept() {
        // 5 节点链: 组 2 (n>2) 与组 1 (n>4) 都不及; 覆盖 3 不满足 rcc 种子
        // （20*3=60 > 邻域 50）→ 整体保留
        let k1 = format!("T{}", &K0[..23]);
        let k3 = sh(K0, 'A');
        let k4 = sh(&k3, 'G');
        let mut tip = vec![sh(K0, 'T')];
        for c in ['C', 'T', 'G', 'C'] {
            tip.push(sh(tip.last().unwrap(), c));
        }
        let mut lines = vec![
            mkline(1, -1, 0, &k1),
            mkline(2, 1, 50, K0),
            mkline(3, 2, 50, &k3),
            mkline(4, 3, 50, &k4),
        ];
        let mut prev = 2;
        for (i, k) in tip.iter().enumerate() {
            lines.push(mkline(200 + i as i64, prev, 3, k));
            prev = 200 + i as i64;
        }
        let g = graph(&lines);
        let (out, st) = clean_graph_rnaspades(&g);
        assert_eq!(st.tip_g1_chains + st.tip_g2_chains, 0);
        assert_eq!(st.rcc_components, 0);
        assert_eq!(out, g);
    }

    #[test]
    fn tip_source_chain_removed_with_connector() {
        // 源向链: 30(根, cov 1) → 31(cov 1) → E(2); E 另有入边 41←40←…（cov 50）
        let k_e = K0;
        let k31 = format!("A{}", &k_e[..23]); // 31 的出边 → E: kmer(E)[..23] = kmer(31)[1..]+...
        // 即 kmer(31) = X + kmer(E)[..23]
        let k30 = format!("C{}", &k31[..23]);
        let k41 = format!("G{}", &k_e[..23]);
        let k40 = format!("T{}", &k41[..23]);
        let g = graph(&[
            mkline(30, -1, 0, &k30),
            mkline(31, 30, 1, &k31),
            mkline(40, -1, 0, &k40),
            mkline(41, 40, 50, &k41),
            mkline(2, 31, 1, k_e),
            mkline(2, 41, 50, k_e),
            mkline(3, 2, 50, &sh(k_e, 'A')),
        ]);
        let (out, st) = clean_graph_rnaspades(&g);
        assert_eq!(st.tip_g2_chains, 1);
        assert_eq!(st.tip_nodes, 2);
        // 30/31 消失; 连接行 (2,31) 删除, (2,41) 保留
        assert!(!out.contains("31\t"));
        assert!(out.contains("2\t41\t50"));
        assert_eq!(alive_nodes(&out), vec![40, 41, 2, 3]);
    }

    // ---------------- EC ----------------

    /// EC 基础夹具: 1(根)→2(J1) → 10(cov_low) → 20(J2); 竞争: 2→5→6,
    /// 40(根)→30→20（均 cov_hi）; 20→7。
    fn ec_fixture(cov_low: i64, cov_hi: i64) -> String {
        let k = |i: i64| format!("{i:0>24}");
        graph(&[
            mkline(1, -1, 0, &k(1)),
            mkline(2, 1, cov_hi, &k(2)),
            mkline(5, 2, cov_hi, &k(5)),
            mkline(6, 5, cov_hi, &k(6)),
            mkline(10, 2, cov_low, &k(10)),
            mkline(20, 10, cov_low, &k(20)),
            mkline(40, -1, 0, &k(40)),
            mkline(30, 40, cov_hi, &k(30)),
            mkline(20, 30, cov_hi, &k(20)),
            mkline(7, 20, cov_hi, &k(7)),
        ])
    }

    #[test]
    fn ec_single_internal_node_removed() {
        // 链 2→10→20: 内部 1 节点, cov 3 ≤ (100+100)/2 + 1 = 101 且 ≤200 → 删
        let (out, st) = clean_graph_rnaspades(&ec_fixture(3, 100));
        assert_eq!(st.ec_chains, 1);
        assert_eq!(st.ec_nodes, 1);
        assert_eq!(alive_nodes(&out), vec![1, 2, 5, 6, 40, 30, 20, 7]);
        assert!(out.contains("20\t30\t100"), "J2 由竞争入边保持定义");
    }

    #[test]
    fn ec_relative_coverage_fail_kept() {
        // 链覆盖 150 > 平均 100+1 → 不删
        let g = ec_fixture(150, 100);
        let (out, st) = clean_graph_rnaspades(&g);
        assert_eq!(st.ec_chains, 0);
        assert_eq!(out, g);
    }

    #[test]
    fn ec_icb_bound_kept() {
        // 链覆盖 250 ≤ 平均 300+1 但 > icb 200 → 不删
        let g = ec_fixture(250, 300);
        let (out, st) = clean_graph_rnaspades(&g);
        assert_eq!(st.ec_chains, 0);
        assert_eq!(out, g);
    }

    #[test]
    fn ec_direct_junction_edge_removed() {
        // 无内部节点的直连 junction 边 (20,2): cov 5 ≤ 101 → 删该行, J2 存活
        let k = |i: i64| format!("{i:0>24}");
        let g = graph(&[
            mkline(1, -1, 0, &k(1)),
            mkline(2, 1, 100, &k(2)),
            mkline(5, 2, 100, &k(5)),
            mkline(6, 5, 100, &k(6)),
            mkline(20, 2, 5, &k(20)),
            mkline(40, -1, 0, &k(40)),
            mkline(30, 40, 100, &k(30)),
            mkline(20, 30, 100, &k(20)),
            mkline(7, 20, 100, &k(7)),
        ]);
        let (out, st) = clean_graph_rnaspades(&g);
        assert_eq!(st.ec_chains, 1);
        assert_eq!(st.ec_nodes, 0);
        assert!(!out.contains("20\t2\t5"));
        assert!(out.contains("20\t30\t100"));
    }

    #[test]
    fn ec_two_internal_removed_three_kept() {
        // 内部 2 节点 → 删; 内部 3 节点（>2）→ 不删（且覆盖 5 不触发 rcc:
        // 20*5=100 不严格小于邻域 100）
        let k = |i: i64| format!("{i:0>24}");
        let base: Vec<String> = vec![
            mkline(1, -1, 0, &k(1)),
            mkline(2, 1, 100, &k(2)),
            mkline(5, 2, 100, &k(5)),
            mkline(6, 5, 100, &k(6)),
            mkline(40, -1, 0, &k(40)),
            mkline(30, 40, 100, &k(30)),
            mkline(20, 30, 100, &k(20)),
            mkline(7, 20, 100, &k(7)),
        ];
        // 2 内部节点: 2→10→11→20
        let mut lines2 = base.clone();
        lines2.push(mkline(10, 2, 5, &k(10)));
        lines2.push(mkline(11, 10, 5, &k(11)));
        lines2.push(mkline(20, 11, 5, &k(20)));
        let (out, st) = clean_graph_rnaspades(&graph(&lines2));
        assert_eq!(st.ec_chains, 1);
        assert_eq!(st.ec_nodes, 2);
        assert!(!out.contains("10\t"));
        // 3 内部节点: 2→10→11→12→20 → 超 ec_lb → 保留
        let mut lines3 = base.clone();
        lines3.push(mkline(10, 2, 5, &k(10)));
        lines3.push(mkline(11, 10, 5, &k(11)));
        lines3.push(mkline(12, 11, 5, &k(12)));
        lines3.push(mkline(20, 12, 5, &k(20)));
        let g3 = graph(&lines3);
        let (out3, st3) = clean_graph_rnaspades(&g3);
        assert_eq!(st3.ec_chains, 0);
        assert_eq!(st3.rcc_components, 0);
        assert_eq!(out3, g3);
    }

    // ---------------- rcc ----------------

    /// rcc 气泡夹具: 主路 1→2→3→4；低覆盖泡 2→10→11→…→3 (cov_low)。
    /// (2,1)=cov_hi; (3,2) 与 (4,3)=cov_rejoin_main。
    /// （泡 ≥3 内部节点使 EC 不及; 首尾衔接处无死端故 tip 不及。）
    fn rcc_fixture(n_blob: usize, cov_low: i64, cov_hi: i64, cov_rejoin_main: i64) -> String {
        let k = |i: i64| format!("{i:0>24}");
        let mut lines = vec![
            mkline(1, -1, 0, &k(1)),
            mkline(2, 1, cov_hi, &k(2)),
            mkline(3, 2, cov_rejoin_main, &k(3)),
            mkline(4, 3, cov_rejoin_main, &k(4)),
        ];
        let mut prev = 2;
        for i in 0..n_blob {
            let id = 10 + i as i64;
            lines.push(mkline(id, prev, cov_low, &k(id)));
            prev = id;
        }
        lines.push(mkline(3, prev, cov_low, &k(3))); // 泡重新汇入 3
        graph(&lines)
    }

    #[test]
    fn rcc_low_cov_bubble_removed() {
        // 3 节点泡 cov 2, 主路 100: 种子 10 (40<100), flood 得 inner={10,11,12},
        // 边界 {2,3} 各有 ≥40 的高覆盖外边 → 整片删
        let (out, st) = clean_graph_rnaspades(&rcc_fixture(3, 2, 100, 100));
        assert_eq!(st.rcc_components, 1);
        assert_eq!(st.rcc_nodes, 3);
        assert_eq!(alive_nodes(&out), vec![1, 2, 3, 4]);
        assert!(out.contains("3\t2\t100"));
        assert!(!out.contains("\t10\t") && !out.contains("10\t2"));
    }

    #[test]
    fn rcc_oversize_kept() {
        // 11 节点泡 > 10 上限 → 保留（逐字节不变）
        let g = rcc_fixture(11, 2, 100, 100);
        let (out, st) = clean_graph_rnaspades(&g);
        assert_eq!(st.rcc_components, 0);
        assert_eq!(out, g);
    }

    #[test]
    fn rcc_boundary_not_high_kept() {
        // 汇点下游边 (3,2)/(4,3) 只有 30（< 20×2=40 不足以让 3 终止）: flood 把
        // 3、4 一并吞入 → 内部基线抬到 30, 阈值 600, 边界 2 的外边 100 < 600
        // → 拒绝, 整体保留（非均匀低覆盖, 主路并入区不删）
        let g = rcc_fixture(3, 2, 100, 30);
        let (out, st) = clean_graph_rnaspades(&g);
        assert_eq!(st.rcc_components, 0);
        assert_eq!(out, g);
    }

    #[test]
    fn rcc_uniform_low_component_kept() {
        // 整组件均匀低覆盖: 无 20× 邻域 → 无种子 → 不删（真实低表达基因保护）
        let k = |i: i64| format!("{i:0>24}");
        let g = graph(&[
            mkline(1, -1, 0, &k(1)),
            mkline(2, 1, 2, &k(2)),
            mkline(3, 2, 2, &k(3)),
        ]);
        let (out, st) = clean_graph_rnaspades(&g);
        assert_eq!(st.rcc_components, 0);
        assert_eq!(out, g);
    }

    // ---------------- 悬空修复 ----------------

    #[test]
    fn dangling_vertex_rerooted() {
        // 两条源向 tip 链汇入同一合并点 E(5)（入度恰 2）: 两条连接行都删后
        // E 失去全部定义行但仍被出边 (6,5) 引用 → 修复为根行
        let k = |i: i64| format!("{i:0>24}");
        let g = graph(&[
            mkline(30, -1, 0, &k(30)),
            mkline(31, 30, 1, &k(31)),
            mkline(40, -1, 0, &k(40)),
            mkline(41, 40, 1, &k(41)),
            mkline(5, 31, 1, &k(5)),
            mkline(5, 41, 1, &k(5)),
            mkline(6, 5, 50, &k(6)),
            mkline(7, 6, 50, &k(7)),
        ]);
        let (out, st) = clean_graph_rnaspades(&g);
        assert_eq!(st.tip_g2_chains, 2);
        assert_eq!(st.repaired_dangling, 1);
        // E 以根行形式复活（count 归零, 与 QuantifyGraph 根行惯例一致）
        assert!(out.contains(&format!("5\t-1\t0\t{}\t1", k(5))), "{out}");
        assert!(out.contains("6\t5\t50"));
        assert_eq!(alive_nodes(&out), vec![5, 6, 7]);
    }
}
