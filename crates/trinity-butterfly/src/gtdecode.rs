//! GTasm 式解码头（`--path-engine gtasm`）——python 原型 `gt_decode2.py` 的
//! C 变体（greedy_mode="gt_sup" × best_mode="sup_sum"）逐语义复刻。
//!
//! 链路: graph.out 文本 → [`parse_component_graph`]（`gt_train3.py load_comps`
//! 语义）→ [`build_gt_features`]（`probe_e2e.py build_all` 语义）→
//! [`GtModel::forward`] 得边概率 → [`decode_component`]。
//!
//! 解码（逐条对应 python `decode()`）:
//! - gsc_e = prob_e × (1 + ln(1 + sup_e))（sup = 原始 reads 跨边计数, 无记录 = 0）;
//! - 循环: 从 alive 边按 ∝(gsc+1e-9) 无放回采 [`NB_STARTS`] 个起点（固定种子
//!   xorshift64*, 采逐 comp 确定性; 不复刻 numpy RNG——采样噪声 ±0.2pp 已验收）;
//! - 每起点边 (s0,d0): 从 d0 前向贪心延伸（候选 = 出边中未 visited/未走过者,
//!   argmax (gsc, vid, eid)——平局 vid 大者优先, 再平 eid 大者）; 从 s0 反向同理;
//! - walk = reverse(back) + fwd（≥2 节点才有效）;
//! - 竞争: 各 walk 按延伸边 sup 总和取最优（首个达到最大者胜, 与 python
//!   `key > best_key` 严格大于一致）;
//! - 无合法 walk → 终止;
//! - 拼序列: 首个有 kmer 名的节点给全 24-mer, 后续节点各补末字符;
//!   **无名节点跳过**（原型 quirk 兼容, 仅 -1 哨兵/纯源节点触发）;
//!   长度 ≥ min_len 才输出;
//! - visited |= walk 节点; alive = 至少一端未 visited 的边。
//!
//! 与 python 的已知差异（验收允许）: RNG 不复刻 numpy; sup 累计用 u64 精确整数
//! （python 是 float32 值的 float64 累加, 计数 << 2^24 时序等价）。

use std::path::Path;

use rustc_hash::{FxHashMap, FxHashSet};

use trinity_common::error::CommonError;

use crate::gtmodel::GtModel;
use crate::{BflyError, Result};

/// 每轮起点采样数（python gt_decode2.NB_STARTS）。
pub const NB_STARTS: usize = 50;

/// 默认解码随机种子（逐 comp 再以 cid 混合——与线程调度无关, 跨运行确定）。
pub const GTASM_DEFAULT_SEED: u64 = 0x67DA_9A11_3C5E_7B2D;

// ---------------------------------------------------------------------------
// 图解析（gt_train3.py load_comps 语义）
// ---------------------------------------------------------------------------

/// load_comps 语义的组件图（注意: **不是** butterfly DiGraph——
/// from=-1 的 root 行也成边, -1 作为普通节点出现且无 kmer 名）。
#[derive(Debug, Clone)]
pub struct GtGraph {
    /// 全部节点 id（升序）= sorted(cnt 键 ∪ out 键)。
    pub vids: Vec<i32>,
    /// vid → 24-mer（仅在 "to" 列出现过的节点有名）。
    pub names: FxHashMap<i32, String>,
    /// vid → 该节点作为 "to" 的 max(count)。
    pub cnt: FxHashMap<i32, i64>,
    /// 边 src/dst（vids 下标）。顺序 = python `out.items()` 的插入序
    /// （frm 首现序）× 组内行序——eid 平局规则依赖此序, 不可改动。
    pub src: Vec<usize>,
    pub dst: Vec<usize>,
}

impl GtGraph {
    pub fn n_edges(&self) -> usize {
        self.src.len()
    }
}

/// 逐行解析 graph.out 文本: 首行 "Component N" 跳过; 数据行 <5 列跳过;
/// `p[0]=to, p[1]=from, p[2]=count, p[3]=kmer`; `names[to]=kmer`,
/// `cnt[to]=max(累计, count)`, `out[from].append(to)`。
pub fn parse_component_graph(text: &str) -> Result<GtGraph> {
    let mut names: FxHashMap<i32, String> = FxHashMap::default();
    let mut cnt: FxHashMap<i32, i64> = FxHashMap::default();
    // python dict 插入序仿真: out_order 记录 frm 键首现序
    let mut out: FxHashMap<i32, Vec<i32>> = FxHashMap::default();
    let mut out_order: Vec<i32> = Vec::new();

    for (i, line) in text.lines().enumerate() {
        if i == 0 {
            continue; // "Component N" 头
        }
        let mut fields = line.split('\t');
        let (Some(f_to), Some(f_from), Some(f_cnt), Some(f_kmer)) = (
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
        ) else {
            continue; // python: len(p) < 5 → skip（含空行）
        };
        if fields.next().is_none() {
            continue; // 恰 4 列 → python len(p) < 5 → skip
        }
        let to: i32 = f_to
            .parse()
            .map_err(|_| BflyError::GraphFileFormat(format!("bad to at line {}: {line}", i + 1)))?;
        let from: i32 = f_from.parse().map_err(|_| {
            BflyError::GraphFileFormat(format!("bad from at line {}: {line}", i + 1))
        })?;
        let c: i64 = f_cnt.parse().map_err(|_| {
            BflyError::GraphFileFormat(format!("bad count at line {}: {line}", i + 1))
        })?;
        names.insert(to, f_kmer.to_string());
        let e = cnt.entry(to).or_insert(0);
        if c > *e {
            *e = c;
        }
        out.entry(from).or_insert_with(|| {
            out_order.push(from);
            Vec::new()
        });
        out.get_mut(&from).unwrap().push(to);
    }

    // vids = sorted(cnt 键 ∪ out 键)
    let mut vids: Vec<i32> = cnt
        .keys()
        .chain(out.keys())
        .copied()
        .collect::<FxHashSet<i32>>()
        .into_iter()
        .collect();
    vids.sort_unstable();
    let remap: FxHashMap<i32, usize> = vids
        .iter()
        .enumerate()
        .map(|(i, &v)| (v, i))
        .collect();

    let mut src = Vec::new();
    let mut dst = Vec::new();
    for f in &out_order {
        for t in &out[f] {
            // 两端必在 vids（构造保证）, 与 python `if f in remap and t in remap` 等价
            src.push(remap[f]);
            dst.push(remap[t]);
        }
    }
    Ok(GtGraph {
        vids,
        names,
        cnt,
        src,
        dst,
    })
}

/// 解码期图特征（probe_e2e.py build_all 语义, 注意与训练期 build_graph 不同——
/// 边特征第 4 维是 `log1p(cnt[dst] 缺省 1)`, 不是训练期的 e2cov）。
///
/// - x[i] = [log1p(cnt), log1p(in_deg), log1p(out_deg)]（log1p(cnt) 以 f64 算后
///   转 f32, 与 `math.log1p → np.float32` 一致; 度数 log1p 直接 f32, 与
///   `np.log1p(float32)` 一致）;
/// - e[ei] = [logc_src, logc_dst, logc_src-logc_dst, log1p(cnt_dst 缺省 1)]。
pub fn build_gt_features(g: &GtGraph) -> (Vec<f32>, Vec<f32>) {
    let n = g.vids.len();
    let logc: Vec<f32> = g
        .vids
        .iter()
        .map(|v| (g.cnt.get(v).copied().unwrap_or(0) as f64).ln_1p() as f32)
        .collect();
    let mut in_deg = vec![0u32; n];
    let mut out_deg = vec![0u32; n];
    for ei in 0..g.n_edges() {
        out_deg[g.src[ei]] += 1;
        in_deg[g.dst[ei]] += 1;
    }
    let mut x = vec![0f32; n * 3];
    for i in 0..n {
        x[i * 3] = logc[i];
        x[i * 3 + 1] = (in_deg[i] as f32).ln_1p();
        x[i * 3 + 2] = (out_deg[i] as f32).ln_1p();
    }
    let ne = g.n_edges();
    let mut e = vec![0f32; ne * 4];
    for ei in 0..ne {
        let (s, d) = (g.src[ei], g.dst[ei]);
        let cnt_dst = g.cnt.get(&g.vids[d]).copied().unwrap_or(1);
        e[ei * 4] = logc[s];
        e[ei * 4 + 1] = logc[d];
        e[ei * 4 + 2] = logc[s] - logc[d];
        e[ei * 4 + 3] = (cnt_dst as f64).ln_1p() as f32;
    }
    (x, e)
}

// ---------------------------------------------------------------------------
// 原始 reads 边支持表（calc_edge_support2.py 的产物语义）
// ---------------------------------------------------------------------------

/// (u, v) 两个 i32 节点 id 打包成 u64 键（-1 哨兵亦可打包——按位拼接）。
#[inline]
pub fn pack_uv(u: i32, v: i32) -> u64 {
    ((u as u32 as u64) << 32) | (v as u32 as u64)
}

/// 边支持表: (cid, pack(u,v)) → 原始 reads 跨边计数。
/// 只存图中真实存在的边（非边转移解码期永不查询, 计数时即丢弃——
/// 与 python Counter 全量记录对解码结果等价, 内存小一个量级）。
#[derive(Debug, Default)]
pub struct EdgeSupport {
    map: FxHashMap<(u64, u64), u32>,
}

impl EdgeSupport {
    pub fn new() -> Self {
        Self::default()
    }

    /// 计数 +1（支持 pass 用）。
    pub fn add(&mut self, cid: u64, u: i32, v: i32) {
        *self.map.entry((cid, pack_uv(u, v))).or_insert(0) += 1;
    }

    /// 查询（无记录 = 0, 与 python `support.get(key, 0)` 一致）。
    pub fn get(&self, cid: u64, u: i32, v: i32) -> u32 {
        self.map
            .get(&(cid, pack_uv(u, v)))
            .copied()
            .unwrap_or(0)
    }

    /// 并入另一张表（多线程分片合并用）。
    pub fn merge(&mut self, other: &EdgeSupport) {
        for (&k, &c) in &other.map {
            *self.map.entry(k).or_insert(0) += c;
        }
    }

    pub fn len(&self) -> usize {
        self.map.len()
    }

    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// 落盘（checkpoint 语义: 池阶段重跑时直接加载, 不重扫 reads）。
    /// 格式（小端）: magic "GTS1" | u64 n | n × (u64 cid, i32 u, i32 v, u32 cnt),
    /// 记录按 (cid, u, v) 升序（确定性）。
    pub fn save(&self, path: &Path) -> std::result::Result<(), CommonError> {
        let mut recs: Vec<(u64, i32, i32, u32)> = self
            .map
            .iter()
            .map(|(&(cid, uv), &c)| {
                (
                    cid,
                    (uv >> 32) as u32 as i32,
                    (uv & 0xFFFF_FFFF) as u32 as i32,
                    c,
                )
            })
            .collect();
        recs.sort_unstable();
        let mut buf = Vec::with_capacity(8 + recs.len() * 20);
        buf.extend_from_slice(b"GTS1");
        buf.extend_from_slice(&(recs.len() as u64).to_le_bytes());
        for (cid, u, v, c) in recs {
            buf.extend_from_slice(&cid.to_le_bytes());
            buf.extend_from_slice(&u.to_le_bytes());
            buf.extend_from_slice(&v.to_le_bytes());
            buf.extend_from_slice(&c.to_le_bytes());
        }
        std::fs::write(path, buf)?;
        Ok(())
    }

    /// 从落盘文件加载（[`EdgeSupport::save`] 的逆）。
    pub fn load(path: &Path) -> std::result::Result<EdgeSupport, CommonError> {
        let data = std::fs::read(path)?;
        if data.len() < 8 || &data[..4] != b"GTS1" {
            return Err(CommonError::Parse(format!(
                "{}: 不是 GTS1 边支持文件",
                path.display()
            )));
        }
        let n = u64::from_le_bytes(data[4..12].try_into().unwrap()) as usize;
        let rec_len = 20usize;
        if data.len() != 12 + n * rec_len {
            return Err(CommonError::Parse(format!(
                "{}: GTS1 长度不符（头 {} 字节 + {n}×{rec_len}, 实际 {}）",
                path.display(),
                12,
                data.len()
            )));
        }
        let mut map = FxHashMap::with_capacity_and_hasher(n, Default::default());
        for i in 0..n {
            let b = &data[12 + i * rec_len..12 + (i + 1) * rec_len];
            let cid = u64::from_le_bytes(b[0..8].try_into().unwrap());
            let u = i32::from_le_bytes(b[8..12].try_into().unwrap());
            let v = i32::from_le_bytes(b[12..16].try_into().unwrap());
            let c = u32::from_le_bytes(b[16..20].try_into().unwrap());
            map.insert((cid, pack_uv(u, v)), c);
        }
        Ok(EdgeSupport { map })
    }
}

// ---------------------------------------------------------------------------
// 解码器（C 变体）
// ---------------------------------------------------------------------------

/// 解码参数。
#[derive(Debug, Clone)]
pub struct GtasmParams {
    /// 最短输出长度（原型恒 200; 产品化接到 --min_contig_length）。
    pub min_len: usize,
    /// 采样种子（逐 comp 再混 cid）。
    pub seed: u64,
}

/// xorshift64* PRNG（采样用; 确定性, 与线程调度无关）。
pub struct XorShift64(u64);

impl XorShift64 {
    pub fn new(seed: u64) -> Self {
        XorShift64(seed.max(1))
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// [0,1) 均匀（取高 53 位）。
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// python `max((gsc, vid, eid))` 元组语义: 先比分, 平则 vid 大者, 再平 eid 大者。
#[inline]
fn better_cand(cand: (f32, i32, usize), best: Option<(f32, i32, usize)>) -> bool {
    match best {
        None => true,
        Some(b) => {
            cand.0 > b.0 || (cand.0 == b.0 && (cand.1 > b.1 || (cand.1 == b.1 && cand.2 > b.2)))
        }
    }
}

/// 单向贪心延伸（前向传 succ / 反向传 pred）。
///
/// 逐 python while 循环: 当前节点不在 visited / frozen / 自身 seen 才入列;
/// 候选 = 邻接中对端不在任何排除集者, 取 argmax（gsc, vid, eid 元组序）;
/// 每迈出一步累计该边 sup。返回 (walk 节点列, 延伸边 sup 和)。
///
/// `frozen` = 另一向已走过的节点（反向延伸时排除前向 walk, python 的
/// `cur not in seen` 交叉检查）; 前向调用传空集。
fn greedy_extend(
    adj: &FxHashMap<i32, Vec<(f32, i32, usize)>>,
    start: i32,
    visited: &FxHashSet<i32>,
    frozen: &FxHashSet<i32>,
    sup: &[u32],
) -> (Vec<i32>, u64) {
    let mut walk = Vec::new();
    let mut seen: FxHashSet<i32> = FxHashSet::default();
    let mut sup_sum = 0u64;
    let mut cur = start;
    loop {
        if visited.contains(&cur) || frozen.contains(&cur) || seen.contains(&cur) {
            break;
        }
        seen.insert(cur);
        walk.push(cur);
        let mut best: Option<(f32, i32, usize)> = None;
        if let Some(cands) = adj.get(&cur) {
            for &cand in cands {
                if visited.contains(&cand.1) || frozen.contains(&cand.1) || seen.contains(&cand.1)
                {
                    continue;
                }
                if better_cand(cand, best) {
                    best = Some(cand);
                }
            }
        }
        match best {
            None => break,
            Some((_, w, e)) => {
                sup_sum += sup[e] as u64;
                cur = w;
            }
        }
    }
    (walk, sup_sum)
}

/// 按 ∝(gsc[e]+1e-9) 无放回采 k 个（f64 累计, 与 python probs 的 float64 一致）。
/// 逐次抽取（numpy choice 的 unique-ify 算法分布略异, 同为 PPS 无放回近似,
/// 验收口径允许 ±0.2pp 噪声）。
fn sample_starts(alive: &[usize], gsc: &[f32], k: usize, rng: &mut XorShift64) -> Vec<usize> {
    let mut pool: Vec<usize> = alive.to_vec();
    let mut out = Vec::with_capacity(k);
    for _ in 0..k {
        let total: f64 = pool.iter().map(|&e| gsc[e] as f64 + 1e-9).sum();
        let mut u = rng.next_f64() * total;
        let mut pos = pool.len() - 1;
        for (i, &e) in pool.iter().enumerate() {
            u -= gsc[e] as f64 + 1e-9;
            if u < 0.0 {
                pos = i;
                break;
            }
        }
        out.push(pool[pos]);
        pool.swap_remove(pos);
    }
    out
}

/// 单组件解码 → contig 序列集（语义见模块文档; cid 用于查支持表与混种子）。
pub fn decode_component(
    g: &GtGraph,
    prob: &[f32],
    support: &EdgeSupport,
    cid: u64,
    params: &GtasmParams,
) -> Vec<String> {
    let ne = g.n_edges();
    if ne == 0 {
        return Vec::new();
    }
    debug_assert_eq!(prob.len(), ne);

    // gsc_e = prob_e × (1 + log1p(sup_e))（全 f32, 与 python numpy float32 一致）
    let sup: Vec<u32> = (0..ne)
        .map(|ei| support.get(cid, g.vids[g.src[ei]], g.vids[g.dst[ei]]))
        .collect();
    let gsc: Vec<f32> = (0..ne)
        .map(|ei| prob[ei] * (1.0f32 + (sup[ei] as f32).ln_1p()))
        .collect();

    // 邻接: vid → [(gsc, 对端 vid, eid)]（eid 序插入; python defaultdict 同序）
    let mut succ: FxHashMap<i32, Vec<(f32, i32, usize)>> = FxHashMap::default();
    let mut pred: FxHashMap<i32, Vec<(f32, i32, usize)>> = FxHashMap::default();
    for (ei, (&s, &d)) in g.src.iter().zip(&g.dst).enumerate() {
        succ.entry(g.vids[s]).or_default().push((gsc[ei], g.vids[d], ei));
        pred.entry(g.vids[d]).or_default().push((gsc[ei], g.vids[s], ei));
    }

    // 逐 comp 确定性种子（splitmix64 混合）
    let mut seed_mix = params.seed ^ cid.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    seed_mix = (seed_mix ^ (seed_mix >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    seed_mix = (seed_mix ^ (seed_mix >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    seed_mix ^= seed_mix >> 31;
    let mut rng = XorShift64::new(seed_mix);

    let mut alive: Vec<usize> = (0..ne).collect();
    let mut visited: FxHashSet<i32> = FxHashSet::default();
    let mut contigs: Vec<String> = Vec::new();

    while !alive.is_empty() {
        let k = NB_STARTS.min(alive.len());
        let picks = sample_starts(&alive, &gsc, k, &mut rng);

        let mut best_walk: Option<Vec<i32>> = None;
        let mut best_key: u64 = 0;
        for e0 in picks {
            let s0 = g.vids[g.src[e0]];
            let d0 = g.vids[g.dst[e0]];

            // 前向: 从 d0 沿出边贪心; 反向: 从 s0 沿入边贪心（排除前向已走过节点）
            let (walk_f, sup_f) =
                greedy_extend(&succ, d0, &visited, &FxHashSet::default(), &sup);
            let seen_f: FxHashSet<i32> = walk_f.iter().copied().collect();
            let (mut walk_b, sup_b) = greedy_extend(&pred, s0, &visited, &seen_f, &sup);

            walk_b.reverse();
            walk_b.extend_from_slice(&walk_f);
            let walk = walk_b;
            if walk.len() < 2 {
                continue;
            }
            let key = sup_f + sup_b; // best_mode="sup_sum"
            if best_walk.is_none() || key > best_key {
                best_walk = Some(walk);
                best_key = key;
            }
        }

        let Some(walk) = best_walk else {
            break; // 无合法 walk → 终止
        };

        // 拼序列: 首个有名节点全 kmer + 后续节点末字符; 无名节点跳过
        let mut seq = String::new();
        for x in &walk {
            if let Some(kk) = g.names.get(x) {
                if seq.is_empty() {
                    seq.push_str(kk);
                } else {
                    seq.push(kk.as_bytes()[kk.len() - 1] as char);
                }
            }
        }
        if seq.len() >= params.min_len {
            contigs.push(seq);
        }
        for x in &walk {
            visited.insert(*x);
        }
        alive.retain(|&e| {
            !visited.contains(&g.vids[g.src[e]]) || !visited.contains(&g.vids[g.dst[e]])
        });
    }
    contigs
}

/// 单组件 gtasm 全链（池调用入口）: graph.out 文本 → allProbPaths FASTA 文本。
///
/// FASTA 头格式与 classic `print_final_paths` 同族（`>{name}_g{j}_i1 len=L
/// path=[gtasm]`）——harvest 的长度过滤/去重与 gene_trans_map 的
/// `c\d+_g\d+_i\d+` 解析直接兼容; 每条 contig 独立成 gene（walk 相互独立）。
pub fn run_component_gtasm(
    graph_text: &str,
    cid: u64,
    name: &str,
    min_len: usize,
    model: &GtModel,
    support: &EdgeSupport,
) -> Result<String> {
    let g = parse_component_graph(graph_text)?;
    if g.n_edges() == 0 {
        return Ok(String::new());
    }
    let (x, e) = build_gt_features(&g);
    let prob = model.forward(&x, &e, &g.src, &g.dst);
    let contigs = decode_component(
        &g,
        &prob,
        support,
        cid,
        &GtasmParams {
            min_len,
            seed: GTASM_DEFAULT_SEED,
        },
    );
    let mut out = String::new();
    for (j, seq) in contigs.iter().enumerate() {
        out.push_str(&format!(">{name}_g{}_i1 len={} path=[gtasm]\n", j + 1, seq.len()));
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

    /// 拼 graph.out 文本行。
    fn row(to: i32, from: i32, cnt: i64, kmer: &str) -> String {
        format!("{to}\t{from}\t{cnt}\t{kmer}\t1\n")
    }

    // -----------------------------------------------------------------------
    // parse_component_graph
    // -----------------------------------------------------------------------

    #[test]
    fn parse_load_comps_semantics() {
        // 含 root 行（from=-1）、重复 to（cnt 取 max、names 后写覆盖）、
        // <5 列行跳过、"Component" 头跳过
        let text = format!(
            "Component 7\n{}{}{}{}{}",
            row(10, -1, 3, "AAAAAAAAAAAAAAAAAAAAAAAA"),
            row(11, 10, 5, "AAAAAAAAAAAAAAAAAAAAAAAC"),
            row(11, 10, 2, "AAAAAAAAAAAAAAAAAAAAAAAC"), // 重复: cnt max(5,2)=5
            "garbage line without tabs\n",
            row(12, 11, 7, "AAAAAAAAAAAAAAAAAAAAAAAG"),
        );
        let g = parse_component_graph(&text).unwrap();
        // vids = {-1, 10, 11, 12}（-1 因 out 键入选）
        assert_eq!(g.vids, vec![-1, 10, 11, 12]);
        assert_eq!(g.cnt[&11], 5);
        assert_eq!(g.names[&10], "AAAAAAAAAAAAAAAAAAAAAAAA");
        assert!(!g.names.contains_key(&-1)); // -1 无名
        // 边序 = out 插入序: -1→10 先, 然后 10→11 ×2, 然后 11→12
        let edges: Vec<(i32, i32)> = g
            .src
            .iter()
            .zip(&g.dst)
            .map(|(&s, &d)| (g.vids[s], g.vids[d]))
            .collect();
        assert_eq!(
            edges,
            vec![(-1, 10), (10, 11), (10, 11), (11, 12)],
            "平行边保留且有序"
        );
    }

    #[test]
    fn features_hand_computed() {
        // 链: -1→10→11; cnt 10→3, 11→5
        let text = format!(
            "Component 0\n{}{}",
            row(10, -1, 3, "AAAAAAAAAAAAAAAAAAAAAAAA"),
            row(11, 10, 5, "AAAAAAAAAAAAAAAAAAAAAAAC"),
        );
        let g = parse_component_graph(&text).unwrap();
        let (x, e) = build_gt_features(&g);
        let log1p64 = |c: i64| (c as f64).ln_1p() as f32;
        // vids = [-1, 10, 11]; in_deg: -1:0, 10:1, 11:1; out_deg: -1:1, 10:1, 11:0
        assert_eq!(x[0], 0.0); // cnt[-1] 缺省 0 → log1p(0)=0
        assert_eq!(x[3], log1p64(3));
        assert_eq!(x[6], log1p64(5));
        assert_eq!(x[1], 0.0); // in_deg(-1)=0
        assert_eq!(x[4], 1.0f32.ln_1p());
        assert_eq!(x[8], 0.0); // out_deg(11)=0（x[6..9] = 节点 11）
        // 边 0: -1→10: [logc(-1)=0, logc(10), 0-logc(10), log1p(cnt[10])]
        assert_eq!(e[0], 0.0);
        assert_eq!(e[1], log1p64(3));
        assert_eq!(e[2], 0.0 - log1p64(3));
        assert_eq!(e[3], log1p64(3));
        // 边 1: 10→11
        assert_eq!(e[4], log1p64(3));
        assert_eq!(e[5], log1p64(5));
        assert_eq!(e[7], log1p64(5));
    }

    // -----------------------------------------------------------------------
    // EdgeSupport 存取往返
    // -----------------------------------------------------------------------

    #[test]
    fn edge_support_save_load_roundtrip() {
        let mut s = EdgeSupport::new();
        s.add(7, 10, 11);
        s.add(7, 10, 11);
        s.add(7, -1, 10); // -1 哨兵打包
        s.add(3, 5, 6);
        assert_eq!(s.get(7, 10, 11), 2);
        assert_eq!(s.get(7, -1, 10), 1);
        assert_eq!(s.get(7, 11, 10), 0); // 方向敏感
        assert_eq!(s.get(99, 1, 2), 0);
        let dir = std::env::temp_dir().join(format!("gts_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("sup.bin");
        s.save(&p).unwrap();
        let s2 = EdgeSupport::load(&p).unwrap();
        assert_eq!(s2.get(7, 10, 11), 2);
        assert_eq!(s2.get(7, -1, 10), 1);
        assert_eq!(s2.get(3, 5, 6), 1);
        assert_eq!(s2.len(), 3);
        let _ = std::fs::remove_dir_all(&dir);
        // 坏文件
        assert!(EdgeSupport::load(Path::new("/nonexistent/x.bin")).is_err());
    }

    // -----------------------------------------------------------------------
    // decode_component
    // -----------------------------------------------------------------------

    /// 手算小图: 链 -1→10→11→12（24-mer 各差一字母）, prob 全 0.5, 无支持。
    /// alive=3 边 ≤ NB_STARTS → 全采; 期望单条 walk 覆盖全链。
    fn chain_graph() -> GtGraph {
        let text = format!(
            "Component 0\n{}{}{}",
            row(10, -1, 3, "AAAAAAAAAAAAAAAAAAAAAAAA"),
            row(11, 10, 5, "AAAAAAAAAAAAAAAAAAAAAAAC"),
            row(12, 11, 7, "AAAAAAAAAAAAAAAAAAAAAAAG"),
        );
        parse_component_graph(&text).unwrap()
    }

    #[test]
    fn decode_chain_single_walk() {
        let g = chain_graph();
        let prob = vec![0.5f32; g.n_edges()];
        let sup = EdgeSupport::new();
        let contigs = decode_component(
            &g,
            &prob,
            &sup,
            0,
            &GtasmParams {
                min_len: 20,
                seed: 42,
            },
        );
        assert_eq!(contigs.len(), 1);
        // -1 无名跳过; 首名节点 10 全 24-mer + 11/12 各 1 字符 = 26bp
        assert_eq!(contigs[0], "AAAAAAAAAAAAAAAAAAAAAAAACG");
    }

    #[test]
    fn decode_min_len_filters() {
        let g = chain_graph();
        let prob = vec![0.5f32; g.n_edges()];
        let sup = EdgeSupport::new();
        let contigs = decode_component(
            &g,
            &prob,
            &sup,
            0,
            &GtasmParams {
                min_len: 200, // 26bp < 200 → 无输出
                seed: 42,
            },
        );
        assert!(contigs.is_empty());
    }

    #[test]
    fn greedy_extend_tuple_tie_breaks() {
        // vid 级平局: (0.5, 11, 0) vs (0.5, 12, 1) → vid 大者 12 胜
        let mut succ: FxHashMap<i32, Vec<(f32, i32, usize)>> = FxHashMap::default();
        succ.insert(10, vec![(0.5, 11, 0), (0.5, 12, 1)]);
        let sup = [0u32; 4];
        let (walk, _) = greedy_extend(&succ, 10, &FxHashSet::default(), &FxHashSet::default(), &sup);
        assert_eq!(walk, vec![10, 12]);

        // eid 级平局（平行边）: (0.5, 11, 0) vs (0.5, 11, 3) → eid 大者 3 胜,
        // sup_sum 计入 eid 3 的支持
        let mut succ2: FxHashMap<i32, Vec<(f32, i32, usize)>> = FxHashMap::default();
        succ2.insert(10, vec![(0.5, 11, 0), (0.5, 11, 3)]);
        let sup2 = [7u32, 0, 0, 9];
        let (walk2, sup_sum2) =
            greedy_extend(&succ2, 10, &FxHashSet::default(), &FxHashSet::default(), &sup2);
        assert_eq!(walk2, vec![10, 11]);
        assert_eq!(sup_sum2, 9);

        // frozen（前向已走）节点不可复走
        let frozen: FxHashSet<i32> = [12].into_iter().collect();
        let (walk3, _) = greedy_extend(&succ, 10, &FxHashSet::default(), &frozen, &sup);
        assert_eq!(walk3, vec![10, 11]);
    }

    #[test]
    fn decode_competition_by_sup_sum() {
        // 两条互斥路径竞争: 起点边不同, 支持总和决定胜负。
        // 图: -1→10→11 (长臂 A) 与 -1→20→21→22 (长臂 B)。
        // 支持: A 臂延伸边 sup=5, B 臂 sup=1 → A 赢（sup_sum 5 > 1+1=2? A 臂
        // 起点边不计 sup: 从 -1→10 起, 前向延伸 10→11 计 5; B 从 -1→20 起,
        // 延伸 20→21, 21→22 各 1 → 2。A 5 > B 2 → A 为首 walk）。
        let text = format!(
            "Component 0\n{}{}{}{}{}",
            row(10, -1, 3, "AAAAAAAAAAAAAAAAAAAAAAAA"),
            row(11, 10, 5, "AAAAAAAAAAAAAAAAAAAAAAAC"),
            row(20, -1, 3, "CCCCCCCCCCCCCCCCCCCCCCCC"),
            row(21, 20, 5, "CCCCCCCCCCCCCCCCCCCCCCCA"),
            row(22, 21, 5, "CCCCCCCCCCCCCCCCCCCCCCAG"),
        );
        let g = parse_component_graph(&text).unwrap();
        let prob = vec![0.5f32; g.n_edges()];
        let mut sup = EdgeSupport::new();
        sup.add(0, 10, 11);
        for _ in 0..4 {
            sup.add(0, 10, 11);
        }
        sup.add(0, 20, 21);
        sup.add(0, 21, 22);
        let contigs = decode_component(
            &g,
            &prob,
            &sup,
            0,
            &GtasmParams {
                min_len: 20,
                seed: 7,
            },
        );
        // 两条臂无共享节点 → visited 不交叉, 两轮各出一条
        assert_eq!(contigs.len(), 2);
        // A 臂先出（竞争胜）——顺序即胜出序; B 臂 = 20 的全 24-mer（24×C）
        // + 21 末字符 A + 22 末字符 G
        assert_eq!(contigs[0], "AAAAAAAAAAAAAAAAAAAAAAAAC");
        assert_eq!(contigs[1], "CCCCCCCCCCCCCCCCCCCCCCCCAG");
    }

    #[test]
    fn decode_empty_graph() {
        let g = parse_component_graph("Component 0\n").unwrap();
        assert_eq!(g.n_edges(), 0);
        let contigs = decode_component(
            &g,
            &[],
            &EdgeSupport::new(),
            0,
            &GtasmParams {
                min_len: 200,
                seed: 1,
            },
        );
        assert!(contigs.is_empty());
    }

    #[test]
    fn xorshift_deterministic() {
        let mut a = XorShift64::new(42);
        let mut b = XorShift64::new(42);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
        let mut c = XorShift64::new(43);
        assert_ne!(a.next_u64(), c.next_u64());
    }

    #[test]
    fn run_component_gtasm_fasta_format() {
        // 用零模型（全 0.5 概率）跑全链, 验证 FASTA 头兼容 harvest/gene_trans_map
        let m = crate::gtmodel::zero_model_for_test();
        // 手工把链拉长到 ≥200bp: 直接构造长链文本
        let mut text = String::from("Component 0\n");
        let kmers: Vec<String> = (0..220)
            .map(|i| {
                let base = [b'A', b'C', b'G', b'T'][(i % 4) as usize];
                let mut k = vec![b'A'; 24];
                k[23] = base;
                String::from_utf8(k).unwrap()
            })
            .collect();
        for i in 0..220i32 {
            text.push_str(&row(100 + i, 99 + i, 5, &kmers[i as usize]));
        }
        let fasta = run_component_gtasm(&text, 0, "c0", 200, &m, &EdgeSupport::new()).unwrap();
        assert!(fasta.starts_with(">c0_g1_i1 len="), "{fasta}");
        assert!(fasta.contains(" path=[gtasm]\n"));
        let seq: String = fasta.lines().filter(|l| !l.starts_with('>')).collect();
        assert!(seq.len() >= 200);
    }
}
