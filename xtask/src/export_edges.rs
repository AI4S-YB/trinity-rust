//! 边表导出器（解图模型的特征/标签数据集构建, M8 探针）。
//!
//! 输入: 一个跑完的 trinity out 目录（或其上级 cache 目录 / chrysalis 目录）
//! 以及参照转录本 FASTA。遍历 chrysalis/Component_bins/Cbin*/c*.graph.out,
//! 每条边吐一行特征 + 标签到 TSV。
//!
//! 标签语义: 参照转录本滑 K−1=24-mer 窗（graph.out 顶点名是 24-mer, 见
//! opt-log M6a 实现坑 3）, 相邻两窗 (k_i, k_{i+1}) 命中同一组件同方向的两
//! 顶点且边存在 → 该边标 1。反义命中按 revcomp 解析, 图游走方向反转。
//! 转换级（v2）: 三连 (k_i, k_{i+1}, k_{i+2}) 标注 (入边,出边) 转换——
//! 解图决策单元的真数据集（边级 95% 正例无训练价值, 探针已证）。
//!
//! 已知取舍（v1, 如实记录）:
//! - 无名源点（只以 from 身份出现的顶点）的 kmer 首位不可知, 记 'N'+出邻居
//!   前缀 23 碱基; 其出边经 **tip 规则**推断标注（label_edges 函数注释）,
//!   启发式噪声量级由 pairs_tip 统计单列监控;
//! - 回文 24-mer / 正反义同名冲突 → 归入 ambiguous, 解析跳过;
//! - 边权 supp 恒 1（opt-log M6a 已证）, 仍输出占位。

use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use trinity_common::kmer::{get_ds_kmer_val, kmer_to_intval};
use trinity_kmer::counter::CountMap;

const K: usize = 24; // graph.out 顶点名长 = KMER_SIZE(25) - 1

// ---------------------------------------------------------------------------
// k-mer 覆盖度（选择切片的灵魂特征）
// ---------------------------------------------------------------------------
//
// 边 (u→v) 恰对应一个 25-mer = u 的 24-mer + v 的末碱基（de Bruijn 邻接）。
// 计数表是管线 jellyfish 阶段的 DS canonical 25-mer 表（mer_counts.*.rs.bin,
// 收尾清理会删——缓存里没有时回退从 insilico_read_normalization/{left,right}
// .norm.fq 重数, 与 both.fa 序列多重集相同, 计数表恒等）。
// 无名源点（kmer 首碱基未知, 'N' 占位）无法编码 → -1（缺失标记, 与 0 区分:
// 0 = 25-mer 在表中不存在）。

/// 边的 25-mer 覆盖度: -1 = 不可知（无名源点）, ≥0 = 表中计数。
pub fn edge_cov(counts: &CountMap, kmer_from: &str, kmer_to: &str) -> i64 {
    if kmer_from.len() != K || kmer_to.len() != K || !clean(kmer_from) {
        return -1;
    }
    let mut k25 = String::with_capacity(K + 1);
    k25.push_str(kmer_from);
    k25.push(kmer_to.as_bytes()[K - 1] as char);
    let Ok(val) = kmer_to_intval(k25.as_bytes()) else {
        return -1;
    };
    let key = get_ds_kmer_val(val, K + 1);
    counts.get(&key).map(|&c| c as i64).unwrap_or(0)
}

/// 定位管线 asm out 目录（含 chrysalis/ 的那层）: 接受 cache 根或 out 本身。
pub fn locate_asm_out(p: &Path) -> Option<PathBuf> {
    [p.to_path_buf(), p.join("out")]
        .into_iter()
        .find(|cand| cand.join("chrysalis").is_dir())
}

/// 取 25-mer DS 计数表: override 路径 > mer_counts bin > norm.fq 重数。
/// 返回 (表, 来源说明)。v1 只支持 K=25 DS（dev/real 全部数据集的形态）。
pub fn counts_for_out_dir(asm_out: &Path, override_path: &str) -> (CountMap, String) {
    let bin = if !override_path.is_empty() {
        PathBuf::from(override_path)
    } else {
        asm_out.join("mer_counts.25.asm.rs.bin")
    };
    if bin.is_file() {
        let m = trinity_kmer::counter::load_counts_bin(&bin, K + 1, true)
            .unwrap_or_else(|e| panic!("load {}: {e}", bin.display()));
        return (m, format!("bin {}", bin.display()));
    }
    let norm = asm_out.join("insilico_read_normalization");
    let mut merged = CountMap::default();
    let mut n_files = 0;
    for name in ["left.norm.fq", "right.norm.fq"] {
        let p = norm.join(name);
        if !p.is_file() {
            continue;
        }
        let data = std::fs::read(&p).unwrap_or_else(|e| panic!("read {}: {e}", p.display()));
        // fq → fa（4 行记录取序列行; both.fa 经 seqtk -A 转换, 序列多重集相同）
        let mut fa = Vec::with_capacity(data.len() / 2 + 64);
        for (i, line) in data.split(|&b| b == b'\n').enumerate() {
            if i % 4 == 1 && !line.is_empty() {
                fa.extend_from_slice(b">r\n");
                fa.extend_from_slice(line);
                fa.push(b'\n');
            }
        }
        let m = trinity_kmer::counter::KmerCountTable::count_fasta_data_with_capacity(
            &fa,
            K + 1,
            true,
            1 << 20,
        );
        for (k, v) in m {
            *merged.entry(k).or_insert(0) += v;
        }
        n_files += 1;
    }
    assert!(
        n_files > 0,
        "{} 下既无 mer_counts.25.asm.rs.bin 也无 left/right.norm.fq",
        asm_out.display()
    );
    (merged, format!("norm.fq 重数 ×{n_files}"))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ori {
    Sense,
    Anti,
}

/// 单组件图: 边按文件行序, names 只含以 to 身份出现过的顶点。
#[derive(Debug, Default)]
pub struct Component {
    pub id: i32,
    /// (from, to, supp)
    pub edges: Vec<(i32, i32, f64)>,
    /// to 顶点 id → 24-mer 名
    pub names: HashMap<i32, String>,
}

impl Component {
    pub fn vertices(&self) -> HashSet<i32> {
        let mut s = HashSet::new();
        for &(f, t, _) in &self.edges {
            s.insert(f);
            s.insert(t);
        }
        s
    }

    /// 顶点 kmer: 命名点直接取; 无名源点取 'N' + 任一出邻居名前缀 23。
    pub fn vertex_kmer(&self, v: i32) -> Option<String> {
        if let Some(k) = self.names.get(&v) {
            return Some(k.clone());
        }
        self.edges
            .iter()
            .find(|&&(f, _, _)| f == v)
            .and_then(|&(_, t, _)| self.names.get(&t))
            .map(|k| format!("N{}", &k[..K - 1]))
    }

    pub fn in_deg(&self, v: i32) -> usize {
        self.edges.iter().filter(|&&(_, t, _)| t == v).count()
    }
    pub fn out_deg(&self, v: i32) -> usize {
        self.edges.iter().filter(|&&(f, _, _)| f == v).count()
    }
}

/// 解析单个 graph.out 文件文本。组件 id 由调用方（文件名）赋予。
pub fn parse_graph_out(id: i32, text: &str) -> Component {
    let mut c = Component {
        id,
        ..Default::default()
    };
    for (i, line) in text.lines().enumerate() {
        if i == 0 || line.is_empty() {
            continue; // "Component N" 头
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 4 {
            continue;
        }
        let to: i32 = f[0].parse().unwrap_or_else(|_| panic!("bad to: {line}"));
        let from: i32 = f[1].parse().unwrap_or_else(|_| panic!("bad from: {line}"));
        let supp: f64 = f[2].parse().unwrap_or_else(|_| panic!("bad supp: {line}"));
        c.edges.push((from, to, supp));
        c.names.insert(to, f[3].to_string());
    }
    c
}

pub fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// 全局 24-mer 索引: kmer → (组件下标, 顶点 id, 方向)。
/// 每个顶点名以正/反义双写; 同名冲突（含跨组件）归入 ambiguous。
pub struct KmerIndex {
    map: HashMap<String, (usize, i32, Ori)>,
    ambiguous: HashSet<String>,
}

impl KmerIndex {
    pub fn build(comps: &[Component]) -> Self {
        let mut idx = KmerIndex {
            map: HashMap::new(),
            ambiguous: HashSet::new(),
        };
        for (ci, c) in comps.iter().enumerate() {
            for (&v, kmer) in &c.names {
                for (key, ori) in [(kmer.clone(), Ori::Sense), (revcomp(kmer), Ori::Anti)] {
                    if idx.ambiguous.contains(&key) {
                        continue;
                    }
                    match idx.map.get(&key) {
                        Some(&(pci, pv, _)) if pci != ci || pv != v => {
                            idx.map.remove(&key);
                            idx.ambiguous.insert(key);
                        }
                        _ => {
                            idx.map.insert(key, (ci, v, ori));
                        }
                    }
                }
            }
        }
        idx
    }

    /// 解析参照 24-mer: 正向命中 → Sense; revcomp 命中 → Anti（顶点为反义走过）。
    pub fn resolve(&self, k: &str) -> Option<(usize, i32, Ori)> {
        if self.ambiguous.contains(k) {
            return None;
        }
        if let Some(&hit) = self.map.get(k) {
            return Some(hit);
        }
        let rc = revcomp(k);
        if self.ambiguous.contains(&rc) {
            return None;
        }
        self.map.get(&rc).copied()
    }

    pub fn n_ambiguous(&self) -> usize {
        self.ambiguous.len()
    }
}

/// 干净 ACGT 窗口才参与（含 N 的窗跳过, 与 kmer 计数语义一致）。
fn clean(k: &str) -> bool {
    k.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

#[derive(Debug, Default)]
pub struct LabelStats {
    pub ref_txs: usize,
    pub ref_windows: usize,
    /// 单窗命中顶点（方向不限）
    pub window_hits: usize,
    /// 成对解析成功且边存在 → 标 1（含 tip 推断）
    pub pairs_labeled: usize,
    /// 其中经 tip 规则标上的（源点出边, 见下）
    pub pairs_tip: usize,
    /// tip 候选但源点入度分叉 ≥2 无名前驱, 无法断定, 跳过
    pub pairs_tip_ambiguous: usize,
    /// 成对解析成功但边不存在（参照声称的邻接图里没有）
    pub pairs_broken: usize,
    /// 成对解析成功但跨组件/跨方向（理论上不该有, 单列观察）
    pub pairs_mismatch: usize,
}

/// 参照走图标边: 返回 (被标 1 的边集合（comp 下标, 边行号）, 统计)。
///
/// tip 规则: 对 (k_i, k_{i+1}), 若 k_{i+1} 正向命中顶点 v 而 k_i 解析失败,
/// 且 v 恰好有 1 个无名前驱 u（源点 kmer 不在文件中）, 则推断 u 的 kmer 即
/// k_i（de Bruijn 邻接强制后缀 23 相等, 仅首碱基未知）, 标边 u→v。启发式:
/// 若参照在 contig 起点之外继续延伸且首碱基不同会误标——转录本边界边恰是
/// 模型最需要学的样本, 收益大于噪声, 量级由 pairs_tip 单列监控。反义方向
/// 不适用（反义对的 from 端必为 to 顶点, 恒有名; 参照与 contig 异向时源点
/// tip 不可见, 如实放弃）。
pub fn label_edges(
    comps: &[Component],
    idx: &KmerIndex,
    refs: &[(String, String)],
) -> (HashSet<(usize, usize)>, LabelStats) {
    let mut labeled = HashSet::new();
    let mut st = LabelStats::default();
    for (_, seq) in refs {
        st.ref_txs += 1;
        let seq = seq.to_uppercase();
        if seq.len() < K + 1 {
            continue;
        }
        for i in 0..=(seq.len() - K - 1) {
            let k1 = &seq[i..i + K];
            let k2 = &seq[i + 1..i + 1 + K];
            if !clean(k1) || !clean(k2) {
                continue;
            }
            st.ref_windows += 1;
            let h1 = idx.resolve(k1);
            let h2 = idx.resolve(k2);
            if h1.is_some() {
                st.window_hits += 1;
            }
            let Some((c2, v2, o2)) = h2 else {
                continue;
            };
            let Some((c1, v1, o1)) = h1 else {
                // tip 规则（仅正向, 见函数注释）
                if o2 == Ori::Sense {
                    let c = &comps[c2];
                    let mut unnamed = c
                        .edges
                        .iter()
                        .enumerate()
                        .filter(|&(_, (f, t, _))| *t == v2 && !c.names.contains_key(f));
                    match (unnamed.next(), unnamed.next()) {
                        (Some((ei, _)), None) => {
                            labeled.insert((c2, ei));
                            st.pairs_labeled += 1;
                            st.pairs_tip += 1;
                        }
                        (Some(_), Some(_)) => st.pairs_tip_ambiguous += 1,
                        _ => {}
                    }
                }
                continue;
            };
            if c1 != c2 || o1 != o2 {
                st.pairs_mismatch += 1;
                continue;
            }
            let (from, to) = match o1 {
                Ori::Sense => (v1, v2),
                Ori::Anti => (v2, v1), // 反义走过: 物理边为两端顶点的正向边
            };
            let c = &comps[c1];
            match c.edges.iter().position(|&(f, t, _)| f == from && t == to) {
                Some(ei) => {
                    labeled.insert((c1, ei));
                    st.pairs_labeled += 1;
                }
                None => st.pairs_broken += 1,
            }
        }
    }
    (labeled, st)
}

// ---------------------------------------------------------------------------
// 序列特征（纯函数, 分母只计 A/C/G/T——无名源点的 'N' 位自然剔除）
// ---------------------------------------------------------------------------

pub fn gc_frac(s: &str) -> f64 {
    let mut gc = 0usize;
    let mut n = 0usize;
    for b in s.bytes() {
        match b {
            b'G' | b'C' => {
                gc += 1;
                n += 1;
            }
            b'A' | b'T' => n += 1,
            _ => {}
        }
    }
    if n == 0 {
        0.0
    } else {
        gc as f64 / n as f64
    }
}

pub fn entropy(s: &str) -> f64 {
    let mut cnt = [0usize; 4];
    let mut n = 0usize;
    for b in s.bytes() {
        let i = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => continue,
        };
        cnt[i] += 1;
        n += 1;
    }
    if n == 0 {
        return 0.0;
    }
    let mut h = 0.0f64;
    for c in cnt {
        if c > 0 {
            let p = c as f64 / n as f64;
            h -= p * p.log2();
        }
    }
    h
}

pub fn max_homopolymer(s: &str) -> usize {
    let mut best = 0;
    let mut run = 0;
    let mut prev = 0u8;
    for b in s.bytes() {
        if b == prev && matches!(b, b'A' | b'C' | b'G' | b'T') {
            run += 1;
        } else {
            run = if matches!(b, b'A' | b'C' | b'G' | b'T') {
                1
            } else {
                0
            };
        }
        best = best.max(run);
        prev = b;
    }
    best
}

// ---------------------------------------------------------------------------
// TSV 输出
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// v2: 转换级（transition）标注——解图决策单元的真数据集
// ---------------------------------------------------------------------------
//
// 探针结论（dev1, 2026-08-26）: 全边级标签 95.2% 正例无训练价值; 决策点出边
// 正例率 38.1%; 全正分叉的配对（phasing）问题边级标签表达不了。转换级:
// 参照三连 (k_i, k_{i+1}, k_{i+2}) 全部解析成功且两边存在 → 标 (入边, 出边)。
// 反义三连: 中间顶点 = v(k_{i+1}), 入边 = (v(k_{i+2})→v(k_{i+1})),
// 出边 = (v(k_{i+1})→v(k_i))——与正向同一条物理转换。
// tip 扩展: 首窗解析失败但中/末窗正向成功且中间顶点恰有 1 个无名前驱 →
// 推断源点入边（与边级 tip 规则同启发式, 量由 triples_tip 单列）。

#[derive(Debug, Default)]
pub struct TransStats {
    pub triples: usize,
    pub triples_labeled: usize,
    pub triples_tip: usize,
    pub triples_tip_ambiguous: usize,
    /// 三窗全解析但至少一条边不存在
    pub triples_broken_edge: usize,
    pub triples_mismatch: usize,
}

/// 返回 (被标 1 的转换集合（comp 下标, 入边行号, 出边行号）, 统计)。
pub fn label_transitions(
    comps: &[Component],
    idx: &KmerIndex,
    refs: &[(String, String)],
) -> (HashSet<(usize, usize, usize)>, TransStats) {
    let mut labeled = HashSet::new();
    let mut st = TransStats::default();
    for (_, seq) in refs {
        let seq = seq.to_uppercase();
        if seq.len() < K + 2 {
            continue;
        }
        for i in 0..=(seq.len() - K - 2) {
            let k1 = &seq[i..i + K];
            let k2 = &seq[i + 1..i + 1 + K];
            let k3 = &seq[i + 2..i + 2 + K];
            if !clean(k1) || !clean(k2) || !clean(k3) {
                continue;
            }
            st.triples += 1;
            let (h1, h2, h3) = (idx.resolve(k1), idx.resolve(k2), idx.resolve(k3));
            let (Some((c2, v2, o2)), Some((c3, v3, o3))) = (h2, h3) else {
                continue;
            };
            if c2 != c3 || o2 != o3 {
                st.triples_mismatch += 1;
                continue;
            }
            let c = &comps[c2];
            let Some((c1, v1, o1)) = h1 else {
                // tip 扩展（仅正向）: v2 恰 1 个无名前驱 → (源点入边, v2→v3)
                if o2 == Ori::Sense {
                    let Some(ei_out) = c.edges.iter().position(|&(f, t, _)| f == v2 && t == v3)
                    else {
                        st.triples_broken_edge += 1;
                        continue;
                    };
                    let mut unnamed = c
                        .edges
                        .iter()
                        .enumerate()
                        .filter(|&(_, (f, t, _))| *t == v2 && !c.names.contains_key(f));
                    match (unnamed.next(), unnamed.next()) {
                        (Some((ei_in, _)), None) => {
                            labeled.insert((c2, ei_in, ei_out));
                            st.triples_labeled += 1;
                            st.triples_tip += 1;
                        }
                        (Some(_), Some(_)) => st.triples_tip_ambiguous += 1,
                        _ => {}
                    }
                }
                continue;
            };
            if c1 != c2 || o1 != o2 {
                st.triples_mismatch += 1;
                continue;
            }
            // 反义: 中间顶点 v(k2) 不变, 入边来自 v(k3), 出边去向 v(k1)
            let (in_from, mid, out_to) = match o1 {
                Ori::Sense => (v1, v2, v3),
                Ori::Anti => (v3, v2, v1),
            };
            let ei_in = c
                .edges
                .iter()
                .position(|&(f, t, _)| f == in_from && t == mid);
            let ei_out = c
                .edges
                .iter()
                .position(|&(f, t, _)| f == mid && t == out_to);
            match (ei_in, ei_out) {
                (Some(a), Some(b)) => {
                    labeled.insert((c1, a, b));
                    st.triples_labeled += 1;
                }
                _ => st.triples_broken_edge += 1,
            }
        }
    }
    (labeled, st)
}

pub const TRANS_TSV_HEADER: &str = "comp\tv\tfrom_u\tto_w\tin_deg_v\tout_deg_v\tout_deg_from_u\tin_deg_to_w\tu_named\tgc_u\tgc_v\tgc_w\tent_u\tent_v\tent_w\thom_u\thom_v\thom_w\te1_label\te2_label\te1_cov\te2_cov\tcomp_vertices\tcomp_edges\tlabel";

/// 候选转换行: 内部决策点（in≥1, out≥2）的 in×out 全组合。
/// e1_label/e2_label 来自边级标注（v1, oracle 列——仅切片分析用, 训练禁用）;
/// e1_cov/e2_cov 为边 25-mer 覆盖度（推理时可用）。
pub fn transition_rows(
    comps: &[Component],
    edge_labeled: &HashSet<(usize, usize)>,
    trans_labeled: &HashSet<(usize, usize, usize)>,
    counts: Option<&CountMap>,
) -> String {
    let mut out = String::new();
    for (ci, c) in comps.iter().enumerate() {
        let nv = c.vertices().len();
        let ne = c.edges.len();
        // 顶点 → 入/出边行号（行号升序 = 文件序, 确定）
        let mut ins: HashMap<i32, Vec<usize>> = HashMap::new();
        let mut outs: HashMap<i32, Vec<usize>> = HashMap::new();
        for (ei, &(f, t, _)) in c.edges.iter().enumerate() {
            outs.entry(f).or_default().push(ei);
            ins.entry(t).or_default().push(ei);
        }
        let mut vs: Vec<i32> = outs
            .iter()
            .filter(|(_, l)| l.len() >= 2)
            .map(|(&v, _)| v)
            .filter(|v| ins.contains_key(v))
            .collect();
        vs.sort_unstable();
        for v in vs {
            let kv = c.vertex_kmer(v).unwrap_or_default();
            for &ei1 in &ins[&v] {
                let (fu, _, _) = c.edges[ei1];
                let ku = c.vertex_kmer(fu).unwrap_or_default();
                for &ei2 in &outs[&v] {
                    let (_, tw, _) = c.edges[ei2];
                    let kw = c.vertex_kmer(tw).unwrap_or_default();
                    let lab = if trans_labeled.contains(&(ci, ei1, ei2)) {
                        1
                    } else {
                        0
                    };
                    let e1l = if edge_labeled.contains(&(ci, ei1)) {
                        1
                    } else {
                        0
                    };
                    let e2l = if edge_labeled.contains(&(ci, ei2)) {
                        1
                    } else {
                        0
                    };
                    let c1v = counts.map_or(-1, |m| edge_cov(m, &ku, &kv));
                    let c2v = counts.map_or(-1, |m| edge_cov(m, &kv, &kw));
                    writeln!(
                        out,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        c.id,
                        v,
                        fu,
                        tw,
                        ins[&v].len(),
                        outs[&v].len(),
                        c.out_deg(fu),
                        c.in_deg(tw),
                        if c.names.contains_key(&fu) { 1 } else { 0 },
                        gc_frac(&ku),
                        gc_frac(&kv),
                        gc_frac(&kw),
                        entropy(&ku),
                        entropy(&kv),
                        entropy(&kw),
                        max_homopolymer(&ku),
                        max_homopolymer(&kv),
                        max_homopolymer(&kw),
                        e1l,
                        e2l,
                        c1v,
                        c2v,
                        nv,
                        ne,
                        lab,
                    )
                    .unwrap();
                }
            }
        }
    }
    out
}

pub const TSV_HEADER: &str = "comp\tfrom\tto\tsupp\tin_deg_from\tout_deg_from\tin_deg_to\tout_deg_to\tfrom_named\tgc_from\tgc_to\tent_from\tent_to\thom_from\thom_to\tcomp_vertices\tcomp_edges\tcov\tlabel";

pub fn edge_rows(
    comps: &[Component],
    labeled: &HashSet<(usize, usize)>,
    counts: Option<&CountMap>,
) -> String {
    let mut out = String::new();
    for (ci, c) in comps.iter().enumerate() {
        let nv = c.vertices().len();
        let ne = c.edges.len();
        for (ei, &(from, to, supp)) in c.edges.iter().enumerate() {
            let kf = c.vertex_kmer(from).unwrap_or_default();
            let kt = c.vertex_kmer(to).unwrap_or_default();
            let lab = if labeled.contains(&(ci, ei)) { 1 } else { 0 };
            let cov = counts.map_or(-1, |m| edge_cov(m, &kf, &kt));
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}",
                c.id,
                from,
                to,
                supp,
                c.in_deg(from),
                c.out_deg(from),
                c.in_deg(to),
                c.out_deg(to),
                if c.names.contains_key(&from) { 1 } else { 0 },
                gc_frac(&kf),
                gc_frac(&kt),
                entropy(&kf),
                entropy(&kt),
                max_homopolymer(&kf),
                max_homopolymer(&kt),
                nv,
                ne,
                cov,
                lab,
            )
            .unwrap();
        }
    }
    out
}

// ---------------------------------------------------------------------------
// 文件定位与 FASTA 读取
// ---------------------------------------------------------------------------

/// 接受 cache 根（内含 out/）、out 目录、或 chrysalis 目录, 返回 Component_bins。
pub fn locate_component_bins(p: &Path) -> Option<PathBuf> {
    [
        p.join("chrysalis/Component_bins"),
        p.join("out/chrysalis/Component_bins"),
        p.join("Component_bins"),
    ]
    .into_iter()
    .find(|cand| cand.is_dir())
}

pub fn read_fasta(path: &Path) -> Vec<(String, String)> {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let text: String = if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut s = String::new();
        use std::io::Read;
        flate2::read::GzDecoder::new(&bytes[..])
            .read_to_string(&mut s)
            .unwrap_or_else(|e| panic!("gunzip {}: {e}", path.display()));
        s
    } else {
        String::from_utf8(bytes).unwrap_or_else(|e| panic!("utf8 {}: {e}", path.display()))
    };
    let mut out = Vec::new();
    let mut name = String::new();
    let mut seq = String::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if !name.is_empty() {
                out.push((std::mem::take(&mut name), std::mem::take(&mut seq)));
            }
            name = h.split_whitespace().next().unwrap_or("").to_string();
        } else {
            seq.push_str(line.trim());
        }
    }
    if !name.is_empty() {
        out.push((name, seq));
    }
    out
}

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

/// 用法: cargo xtask export-edges --out-dir <cache根|out|chrysalis> --ref <ref.fa[.gz]> --tsv <out.tsv>
pub fn run(argv: &[String]) -> i32 {
    let mut out_dir = String::new();
    let mut reff = String::new();
    let mut tsv = String::new();
    let mut trans_tsv = String::new();
    let mut counts_override = String::new();
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "--out-dir" => {
                out_dir = argv.get(i + 1).cloned().unwrap_or_default();
                i += 1;
            }
            "--ref" => {
                reff = argv.get(i + 1).cloned().unwrap_or_default();
                i += 1;
            }
            "--tsv" => {
                tsv = argv.get(i + 1).cloned().unwrap_or_default();
                i += 1;
            }
            "--trans-tsv" => {
                trans_tsv = argv.get(i + 1).cloned().unwrap_or_default();
                i += 1;
            }
            "--counts" => {
                counts_override = argv.get(i + 1).cloned().unwrap_or_default();
                i += 1;
            }
            other => {
                eprintln!("未知参数: {other}");
                return 2;
            }
        }
        i += 1;
    }
    if out_dir.is_empty() || reff.is_empty() || (tsv.is_empty() && trans_tsv.is_empty()) {
        eprintln!("用法: cargo xtask export-edges --out-dir <cache根|out|chrysalis> --ref <ref.fa[.gz]> [--tsv <边表.tsv>] [--trans-tsv <转换表.tsv>] [--counts <mer_counts.rs.bin>]");
        return 2;
    }
    let bins = locate_component_bins(Path::new(&out_dir))
        .unwrap_or_else(|| panic!("在 {out_dir} 下找不到 Component_bins"));
    eprintln!("[export-edges] Component_bins: {}", bins.display());

    // 收集并按组件号排序（文件名 c<N>.graph.out）
    let mut files: Vec<(i32, PathBuf)> = Vec::new();
    for cb in std::fs::read_dir(&bins).unwrap() {
        let cb = cb.unwrap().path();
        if !cb.is_dir() {
            continue;
        }
        for f in std::fs::read_dir(&cb).unwrap() {
            let f = f.unwrap().path();
            let name = f.file_name().unwrap().to_string_lossy().to_string();
            if let Some(rest) = name.strip_suffix(".graph.out") {
                if let Some(idstr) = rest.strip_prefix('c') {
                    if let Ok(id) = idstr.parse::<i32>() {
                        files.push((id, f));
                    }
                }
            }
        }
    }
    files.sort_by_key(|(id, _)| *id);

    let mut comps = Vec::with_capacity(files.len());
    for (id, f) in &files {
        let text =
            std::fs::read_to_string(f).unwrap_or_else(|e| panic!("read {}: {e}", f.display()));
        comps.push(parse_graph_out(*id, &text));
    }
    let n_edges: usize = comps.iter().map(|c| c.edges.len()).sum();
    eprintln!("[export-edges] 组件 {} 个, 边 {} 条", comps.len(), n_edges);

    let idx = KmerIndex::build(&comps);
    eprintln!(
        "[export-edges] 顶点名索引 {} 键（含反义）, ambiguous {} 个",
        idx.map.len(),
        idx.n_ambiguous()
    );

    let refs = read_fasta(Path::new(&reff));
    eprintln!("[export-edges] 参照转录本 {} 条", refs.len());

    let asm_out = locate_asm_out(Path::new(&out_dir))
        .unwrap_or_else(|| panic!("在 {out_dir} 下找不到 asm out（含 chrysalis/）"));
    eprintln!("[export-edges] 构建/装载 25-mer 计数表…");
    let t0 = std::time::Instant::now();
    let (counts, provenance) = counts_for_out_dir(&asm_out, &counts_override);
    eprintln!(
        "[export-edges] 计数表 {} 键（来源: {}, {:.1?}）",
        counts.len(),
        provenance,
        t0.elapsed()
    );
    let counts = Some(&counts);

    let (labeled, st) = label_edges(&comps, &idx, &refs);
    let hit_rate = if st.ref_windows > 0 {
        100.0 * st.window_hits as f64 / st.ref_windows as f64
    } else {
        0.0
    };
    eprintln!(
        "[export-edges] 标签: 窗口命中率 {:.1}% ({}/{}), 成对 labeled {}（其中 tip {} / tip 歧义跳过 {}）/ broken {} / mismatch {}",
        hit_rate, st.window_hits, st.ref_windows, st.pairs_labeled, st.pairs_tip, st.pairs_tip_ambiguous, st.pairs_broken, st.pairs_mismatch
    );
    eprintln!(
        "[export-edges] 标 1 边 {} 条 / 总 {} 条（{:.2}%）",
        labeled.len(),
        n_edges,
        100.0 * labeled.len() as f64 / n_edges.max(1) as f64
    );

    if !tsv.is_empty() {
        let mut out = String::with_capacity(1 << 20);
        out.push_str(TSV_HEADER);
        out.push('\n');
        out.push_str(&edge_rows(&comps, &labeled, counts));
        std::fs::write(&tsv, &out).unwrap_or_else(|e| panic!("write {tsv}: {e}"));
        eprintln!("[export-edges] 写出 {tsv}");
    }

    if !trans_tsv.is_empty() {
        let (tlab, tst) = label_transitions(&comps, &idx, &refs);
        let rows = transition_rows(&comps, &labeled, &tlab, counts);
        let n_rows = rows.lines().count();
        let n_pos = rows.lines().filter(|l| l.ends_with("\t1")).count();
        eprintln!(
            "[export-edges] 转换: 三连 labeled {}（其中 tip {} / tip 歧义 {}）/ broken_edge {} / mismatch {}",
            tst.triples_labeled, tst.triples_tip, tst.triples_tip_ambiguous, tst.triples_broken_edge, tst.triples_mismatch
        );
        eprintln!(
            "[export-edges] 候选转换 {} 行, 正例 {}（{:.1}%）",
            n_rows,
            n_pos,
            100.0 * n_pos as f64 / n_rows.max(1) as f64
        );
        let mut t = String::with_capacity(1 << 20);
        t.push_str(TRANS_TSV_HEADER);
        t.push('\n');
        t.push_str(&rows);
        std::fs::write(&trans_tsv, t).unwrap_or_else(|e| panic!("write {trans_tsv}: {e}"));
        eprintln!("[export-edges] 写出 {trans_tsv}");
    }
    0
}

// ---------------------------------------------------------------------------
// 测试（手算用例锁定）
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// 合成组件: 4 顶点直线 1->2->3->4（v1 无名源点）。
    /// 底层路径 S 长 27nt → 4 个 24-mer 滑窗构成真实 de Bruijn 邻接;
    /// S 特意选为非回文（ACGT 周期序列 24-mer 是自身 revcomp, 会干扰方向测试）。
    fn synth() -> (Component, String, String, String) {
        let s = "TGCAGTCCGATGACAATTCCGGTAAGC";
        let k1 = &s[0..24];
        let k2 = &s[1..25];
        let k3 = &s[2..26];
        let k4 = &s[3..27];
        let text = format!("Component 7\n2\t1\t1\t{k2}\t1\n3\t2\t1\t{k3}\t1\n4\t3\t1\t{k4}\t1\n");
        (
            parse_graph_out(7, &text),
            k1.to_string(),
            k2.to_string(),
            k3.to_string(),
        )
    }

    #[test]
    fn parse_edges_and_names() {
        let (c, _, k2, k3) = synth();
        assert_eq!(c.edges.len(), 3);
        assert_eq!(c.edges[0], (1, 2, 1.0));
        assert_eq!(c.names.get(&2).unwrap(), &k2);
        assert_eq!(c.names.get(&3).unwrap(), &k3);
        assert!(!c.names.contains_key(&1)); // 源点无名
        assert_eq!(c.vertices().len(), 4);
    }

    #[test]
    fn source_kmer_derived_with_n() {
        let (c, k1, _, _) = synth();
        let got = c.vertex_kmer(1).unwrap();
        // 'N' + 出邻居 k2 前缀23 = 'N' + k1[1..24]
        assert_eq!(got.len(), 24);
        assert!(got.starts_with('N'));
        assert_eq!(got[1..], k1[1..]);
    }

    #[test]
    fn label_sense_walk() {
        let (c, _, _, _) = synth();
        let comps = vec![c];
        let idx = KmerIndex::build(&comps);
        assert_eq!(idx.n_ambiguous(), 0);
        // 参照 = 完整底层路径 → 3 条边全标 1
        let refs = vec![("r1".into(), "TGCAGTCCGATGACAATTCCGGTAAGC".to_string())];
        let (labeled, st) = label_edges(&comps, &idx, &refs);
        assert_eq!(labeled.len(), 3);
        assert_eq!(st.pairs_labeled, 3);
        assert_eq!(st.pairs_broken, 0);
        assert_eq!(st.ref_windows, 3);
    }

    #[test]
    fn label_anti_walk() {
        let (c, _, _, _) = synth();
        let comps = vec![c];
        let idx = KmerIndex::build(&comps);
        // 参照 = 底层路径的 revcomp → 反义解析。
        // 边 (2,3),(3,4) 标 1; 源点 tip 边 (1,2) 反义不可见（见 label_edges 注释）。
        let refs = vec![("r1".into(), revcomp("TGCAGTCCGATGACAATTCCGGTAAGC"))];
        let (labeled, st) = label_edges(&comps, &idx, &refs);
        assert_eq!(labeled.len(), 2);
        assert_eq!(st.pairs_labeled, 2);
        assert_eq!(st.pairs_tip, 0);
    }

    #[test]
    fn tip_rule_labels_source_edge() {
        let (c, _, _, _) = synth();
        let comps = vec![c];
        let idx = KmerIndex::build(&comps);
        let refs = vec![("r1".into(), "TGCAGTCCGATGACAATTCCGGTAAGC".to_string())];
        let (labeled, st) = label_edges(&comps, &idx, &refs);
        // 3 条边全标 1, 其中源点边 (1,2) 来自 tip 规则
        assert_eq!(labeled.len(), 3);
        assert_eq!(st.pairs_labeled, 3);
        assert_eq!(st.pairs_tip, 1);
    }

    #[test]
    fn broken_pair_counted() {
        // 分叉组件: 0->1(k_a) ->2(k_b) ->8(k_y); 另有孤立边 6->7(k_x)。
        // k_a="G"+A*23, k_b="A"*23+"C", k_y="A"*22+"CA", k_x="A"*22+"CG"。
        // 参照窗对 (k_b, k_x) 是合法 de Bruijn 邻接（后缀/前缀 = "A"*22+"C"）,
        // 两顶点同组件同方向但边 2->7 不存在 → broken。反义同构。
        let k_a = "GAAAAAAAAAAAAAAAAAAAAAAA";
        let k_b = "AAAAAAAAAAAAAAAAAAAAAAAC";
        let k_y = "AAAAAAAAAAAAAAAAAAAAAACA";
        let k_x = "AAAAAAAAAAAAAAAAAAAAAACG";
        let text = format!(
            "Component 1\n1\t0\t1\t{k_a}\t1\n2\t1\t1\t{k_b}\t1\n8\t2\t1\t{k_y}\t1\n7\t6\t1\t{k_x}\t1\n"
        );
        let c = parse_graph_out(1, &text);
        let comps = vec![c];
        let idx = KmerIndex::build(&comps);
        let refs = vec![("r".into(), format!("{k_b}G"))]; // 窗对 (k_b, k_x)
        let (labeled, st) = label_edges(&comps, &idx, &refs);
        assert_eq!(st.pairs_broken, 1);
        assert_eq!(st.pairs_labeled, 0);
        assert!(labeled.is_empty());
        // 反义: from=v(第二窗)=v2, to=v(第一窗)=v7, 边 2->7 不存在 → 同样 broken
        let refs2 = vec![("r".into(), revcomp(&format!("{k_b}G")))];
        let (_, st2) = label_edges(&comps, &idx, &refs2);
        assert_eq!(st2.pairs_broken, 1);
    }

    #[test]
    fn features_hand_computed() {
        assert_eq!(gc_frac("ACGT"), 0.5);
        assert_eq!(gc_frac("NNNN"), 0.0);
        assert_eq!(entropy("AAAA"), 0.0);
        assert!((entropy("ACGT") - 2.0).abs() < 1e-9);
        assert!((entropy("AACC") - 1.0).abs() < 1e-9);
        assert_eq!(max_homopolymer("AACCCGT"), 3);
        assert_eq!(max_homopolymer("ACGT"), 1);
    }

    #[test]
    fn rows_schema_and_label() {
        let (c, _, _, _) = synth();
        let comps = vec![c];
        let idx = KmerIndex::build(&comps);
        let refs = vec![("r1".into(), "TGCAGTCCGATGACAATTCCGGTAAGC".to_string())];
        let (labeled, _) = label_edges(&comps, &idx, &refs);
        let rows = edge_rows(&comps, &labeled, None);
        let lines: Vec<&str> = rows.lines().collect();
        assert_eq!(lines.len(), 3);
        let cols: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(cols.len(), TSV_HEADER.split('\t').count());
        assert_eq!(cols[0], "7"); // comp id
        assert!(lines.iter().all(|l| l.ends_with("\t1"))); // 全标 1
        assert_eq!(TSV_HEADER.split('\t').count(), 19);
    }

    #[test]
    fn ambiguous_revcomp_collision_skipped() {
        // 两个顶点互为 revcomp（非回文）→ 该 24-mer 及其 rc 均 ambiguous, 解析跳过
        let k = "TGCAGTCCGATGACAATTCCGGTA";
        let kr = revcomp(k);
        assert_ne!(k, kr);
        let text = format!("Component 1\n2\t1\t1\t{k}\t1\n3\t2\t1\t{kr}\t1\n");
        let c = parse_graph_out(1, &text);
        let comps = vec![c];
        let idx = KmerIndex::build(&comps);
        assert!(idx.resolve(k).is_none());
        assert!(idx.resolve(&kr).is_none());
        assert_eq!(idx.n_ambiguous(), 2);
    }

    // -----------------------------------------------------------------------
    // v2 转换级测试
    // -----------------------------------------------------------------------

    /// 分叉组件（comp 7）: v2 有两个入边（1→2 无名源点, 10→2 命名 m1）
    /// 与两个出边（2→3 命名 k3, 2→5 命名 m3）。
    /// k1..k3 同 synth; m1="A"+k2[0..23], m3=k2[1..24]+"T"（de Bruijn 一致）。
    /// 边行号: 0=(1→2) 1=(2→3) 2=(9→10) 3=(10→2) 4=(2→5)。
    fn fork() -> (Component, String, String, String, String) {
        let k2 = "GCAGTCCGATGACAATTCCGGTAA";
        let k3 = "CAGTCCGATGACAATTCCGGTAAG";
        let m1 = "AGCAGTCCGATGACAATTCCGGTA";
        let m3 = "CAGTCCGATGACAATTCCGGTAAT";
        let text = format!(
            "Component 7\n2\t1\t1\t{k2}\t1\n3\t2\t1\t{k3}\t1\n10\t9\t1\t{m1}\t1\n2\t10\t1\t{k2}\t1\n5\t2\t1\t{m3}\t1\n"
        );
        (
            parse_graph_out(7, &text),
            k2.to_string(),
            k3.to_string(),
            m1.to_string(),
            m3.to_string(),
        )
    }

    #[test]
    fn trans_sense_fork() {
        let (c, _k2, k3, m1, m3) = fork();
        let comps = vec![c];
        let idx = KmerIndex::build(&comps);
        assert_eq!(idx.n_ambiguous(), 0);
        // seqA: 窗 m1,k2,k3 → 转换 (10→2, 2→3)=(edge3,edge1)
        // seqB: 窗 m1,k2,m3 → 转换 (10→2, 2→5)=(edge3,edge4)
        let seq_a = format!("{m1}AG");
        let seq_b = format!("{m1}AT");
        assert_eq!(seq_a.len(), 26);
        assert_eq!(&seq_a[2..26], &k3[..]);
        assert_eq!(&seq_b[2..26], &m3[..]);
        let refs = vec![("a".into(), seq_a), ("b".into(), seq_b)];
        let (tlab, st) = label_transitions(&comps, &idx, &refs);
        assert_eq!(st.triples, 2);
        assert_eq!(st.triples_labeled, 2);
        assert!(tlab.contains(&(0, 3, 1)));
        assert!(tlab.contains(&(0, 3, 4)));
        assert_eq!(tlab.len(), 2);
    }

    #[test]
    fn trans_anti_labels_same_physical() {
        let (c, _k2, k3, m1, _m3) = fork();
        let comps = vec![c];
        let idx = KmerIndex::build(&comps);
        let seq_a = format!("{m1}AG");
        assert_eq!(&seq_a[2..26], &k3[..]);
        let refs = vec![("a".into(), revcomp(&seq_a))];
        let (tlab, st) = label_transitions(&comps, &idx, &refs);
        assert_eq!(st.triples_labeled, 1);
        assert!(tlab.contains(&(0, 3, 1))); // 与正向同一条物理转换
    }

    #[test]
    fn trans_tip_extension() {
        let (c, _, _, _) = synth();
        let comps = vec![c];
        let idx = KmerIndex::build(&comps);
        // 完整路径 27nt → 2 个三连: (k1,k2,k3) 走 tip（k1 不可解析,
        // v2 恰 1 个无名前驱 v1）, (k2,k3,k4) 正常。
        let refs = vec![("r".into(), "TGCAGTCCGATGACAATTCCGGTAAGC".to_string())];
        let (tlab, st) = label_transitions(&comps, &idx, &refs);
        assert_eq!(st.triples, 2);
        assert_eq!(st.triples_labeled, 2);
        assert_eq!(st.triples_tip, 1);
        assert!(tlab.contains(&(0, 0, 1))); // (1→2, 2→3) tip 推断
        assert!(tlab.contains(&(0, 1, 2))); // (2→3, 3→4) 正常
    }

    #[test]
    fn trans_rows_candidates_and_labels() {
        let (c, _k2, _k3, m1, _m3) = fork();
        let comps = vec![c];
        let idx = KmerIndex::build(&comps);
        let refs = vec![
            ("a".into(), format!("{m1}AG")),
            ("b".into(), format!("{m1}AT")),
        ];
        let (elab, _) = label_edges(&comps, &idx, &refs);
        let (tlab, _) = label_transitions(&comps, &idx, &refs);
        let rows = transition_rows(&comps, &elab, &tlab, None);
        let lines: Vec<&str> = rows.lines().collect();
        // 仅 v2 是决策点: 入边 {edge0,edge3} × 出边 {edge1,edge4} = 4 行
        assert_eq!(lines.len(), 4);
        let ncol = TRANS_TSV_HEADER.split('\t').count();
        assert_eq!(ncol, 25);
        let mut pos = 0;
        for l in &lines {
            let cols: Vec<&str> = l.split('\t').collect();
            assert_eq!(cols.len(), ncol);
            let from_u: i32 = cols[2].parse().unwrap();
            let lab = cols[ncol - 1];
            let e1l = cols[ncol - 7];
            if from_u == 10 {
                // (10→2) 的两条转换都被参照走过
                assert_eq!(lab, "1");
                assert_eq!(e1l, "1"); // 边 (10→2) 边级标签也是 1
                pos += 1;
            } else {
                assert_eq!(from_u, 1);
                assert_eq!(lab, "0"); // 无名源点入边: 转换不可标注（tip 不含）
                assert_eq!(e1l, "0");
            }
        }
        assert_eq!(pos, 2);
    }

    #[test]
    fn edge_cov_canonical_and_missing() {
        // 25-mer "TGCAGTCCGATGACAATTCCGGTAA"（= k1+"A"）计数 7, DS canonical 存表。
        let k25 = "TGCAGTCCGATGACAATTCCGGTAA";
        let key = get_ds_kmer_val(kmer_to_intval(k25.as_bytes()).unwrap(), 25);
        let mut counts = CountMap::default();
        counts.insert(key, 7u32);
        // from=k1(24) to=k2(24, 末碱基 'A') → 25-mer 恰好是 k25
        let k1 = "TGCAGTCCGATGACAATTCCGGTA";
        let k2 = "GCAGTCCGATGACAATTCCGGTAA";
        assert_eq!(edge_cov(&counts, k1, k2), 7);
        // 反义 25-mer 查同一 canonical 键: from=rc(k25)[0..24], to 末碱基=rc(k25)[24]
        let rc = revcomp(k25);
        let (rf, rt) = (&rc[..24], &rc[1..25]);
        assert_eq!(edge_cov(&counts, rf, rt), 7);
        // 表中不存在的 25-mer → 0（与 -1 区分）
        let absent_from = "CCCCCCCCCCCCCCCCCCCCCCCC";
        let absent_to = "CCCCCCCCCCCCCCCCCCCCCCCA";
        assert_eq!(edge_cov(&counts, absent_from, absent_to), 0);
        // 无名源点（'N' 占位）→ -1
        let n_from = format!("N{}", &k2[..23]);
        assert_eq!(edge_cov(&counts, &n_from, k2), -1);
    }
}
