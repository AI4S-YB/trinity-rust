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
//! **visited 穿越（Rust 版扩展, [`GtasmParams::allow_visited_walk`]; CLI 默认开,
//! `--no-gtasm-visited-walk` 关）**: 原型把 visited 节点硬排除——强势 isoform 的
//! walk 先把共享外显子标 visited, 弱势 isoform 的 walk 到边界即截断, 完整路径被
//! 拆散吞没（dev1 实测 68 条引擎差距 principal 中 16 条 ≥95% 路径节点已被其他
//! gtasm 产出覆盖）。开启后:
//! - 候选两级过滤: 优先 `not in seen/frozen 且未 visited`（原型语义）; 一级为
//!   空再取 `not in seen/frozen` 的已 visited 候选（同 (gsc, vid, eid) 元组
//!   取 max）; walk 内 seen/frozen 防环逻辑不变;
//! - 重复产出防护: 组件内已产出序列的 64 位指纹集, 完全相同才丢弃（穿越可能
//!   重新发现与已产出完全相同的 walk; 共享部分节点的不同 isoform 照常产出）;
//! - 终止语义不变: 每轮胜出 walk 必含其起点边两端 → 该边必退出 alive; visited
//!   只增不减, 至多 ne 轮, 穿越 walk 不会使 alive 复活。
//!
//! **简单 bulge 双臂产出（Rust 版扩展, [`GtasmParams::bulge_dual`]; CLI
//! `--gtasm-bulge-dual`, 默认关）**: 前向延伸在节点 cur 遇到恰好 2 个一级
//! （未 visited、未 seen/frozen）候选边且指向不同 vid、并构成简单 bulge
//! （两支内部出度恒 1、≤ [`BULGE_MAX_CHAIN`] 节点内重汇合到同一节点 C 且
//! indeg(C)==2——rnaspades `bulge_merge` 语义）时, 不再按 (gsc,vid,eid)
//! 二选一, 而是分叉: 双臂各走一条, 重汇合后共享后缀继续贪心（串联 bulge
//! 可再分叉）。机制语义:
//! - 竞争: 每起点的全部 fork walk 以各自 sup_sum 的最大值为该起点 key 参与
//!   胜出竞争（严格大于、首达者胜——与单 walk 语义一致）; **胜出起点的全部
//!   fork 都参与产出**（各自独立拼序列/长度阈值/指纹去重）;
//! - visited: 并入胜出起点全部 fork 的节点并集——双臂同轮"办结", 不会下轮
//!   再各走一遍; 起点边两端仍必入 visited（fork 共享 back 前缀或同起于
//!   d0）, 逐轮 alive 单调退出, 终止性不变;
//! - 防环: 每 fork 独立 seen 集, walk 内节点不复走; fork 只发生在一级
//!   （未 visited）候选上, 已 visited 臂留给穿越机制;
//! - 安全阀: 单起点分叉上限 [`MAX_FORKS_PER_START`]（串联 bulge 2^n 截断,
//!   超限退回 argmax 单选）; 组件产出上限 [`MAX_COMPONENT_CONTIGS`]
//!   （rnaspades max_paths 同值）, 达到即停止解码;
//! - 去重: bulge_dual 下指纹去重随动开启（两臂 kmer 名碰巧拼出相同全序列
//!   时不重复产出）, 与 visited 穿越的去重共用同一指纹集;
//! - 非 bulge 双歧义 / ≥3 歧义维持 argmax 二选一现状; 反向延伸不分叉
//!   （分叉点是出度歧义, 上游汇合歧义由采样起点覆盖）。
//!
//! **read pair 分支证据（Rust 版扩展, [`GtasmParams::pair_boost`]; CLI
//! `--gtasm-pair-support`, 默认关）**: 单 read 76bp 够不着 300-400bp 的分支
//! 间距（look-back 绑定已证无效）, pair（insert ~300-400bp）够得着。支持
//! pass（trinity-cli `gt_support`, 近似静态版）对每对 mates: 一条的节点
//! 路径含分支出边 (v,w)（v 出度≥2）且另一条命中 v 的上游窗口（反向 BFS
//! ≤15 边）→ [`PairSupport`] (cid,v,w) += 1。解码消费为乘性变体:
//! gsc_e = prob_e × (1+ln(1+sup_e)) × (1+ln(1+pair_e))——与 sup 同构的
//! 静态因子, 起点采样与候选 argmax 同吃; 非分支出边 pair=0 → 因子 1,
//! pair_boost 关时逐点 = 现状。pair 与 sup 是不同证据, 表/文件/字段分离。
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

/// bulge 支链长度上限（节点数; 语义同 rnaspades `bulge_max_chain`, SNP bubble
/// 支长 ≈ k, 50 富余）。
pub const BULGE_MAX_CHAIN: usize = 50;

/// 单起点前向分叉安全阀（bulge_dual 开; 串联 bulge 每多一个路径数 ×2, 超限的
/// 分叉点退回 argmax 单选——参考 rnaspades max_paths 的防爆语义, 按起点计）。
pub const MAX_FORKS_PER_START: usize = 64;

/// 每组件产出安全阀（rnaspades `max_paths` 同值）: 达到即停止解码。
pub const MAX_COMPONENT_CONTIGS: usize = 4096;

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
// read pair 分支支持表（gt_support pair pass 的产物语义）
// ---------------------------------------------------------------------------

/// read pair 分支支持表: (cid, pack(v,w)) → 跨过分支点 v→w 的 read pair 数
/// （一条 mate 的节点路径含分支出边 (v,w), 另一条命中 v 的上游窗口——
/// 窗口 = 自 v 反向 BFS ≤15 边的节点集, 见 trinity-cli `gt_support`）。
///
/// **与 [`EdgeSupport`]（单 read 跨边计数）是不同证据**: 字段分离、文件分离
/// （并列 sidecar, magic "GTP1"; 记录布局与 GTS1 相同——(cid,u,v,cnt) 同构,
/// 仅 magic 区分, 防止误读混用）, 解码期走 [`GtasmParams::pair_boost`] 独立
/// 消费, 绝不混入 gsc 的 sup 项。
#[derive(Debug, Default)]
pub struct PairSupport {
    map: FxHashMap<(u64, u64), u32>,
}

impl PairSupport {
    pub fn new() -> Self {
        Self::default()
    }

    /// 计数 +1（pair 支持 pass 用; 调用方已对每对每边去重）。
    pub fn add(&mut self, cid: u64, v: i32, w: i32) {
        *self.map.entry((cid, pack_uv(v, w))).or_insert(0) += 1;
    }

    /// 查询（无记录 = 0, 语义同 [`EdgeSupport::get`]）。
    pub fn get(&self, cid: u64, v: i32, w: i32) -> u32 {
        self.map
            .get(&(cid, pack_uv(v, w)))
            .copied()
            .unwrap_or(0)
    }

    /// 并入另一张表（多线程分片合并用）。
    pub fn merge(&mut self, other: &PairSupport) {
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

    /// 落盘（checkpoint 语义同 [`EdgeSupport::save`]）。
    /// 格式（小端）: magic "GTP1" | u64 n | n × (u64 cid, i32 v, i32 w, u32 cnt),
    /// 记录按 (cid, v, w) 升序（确定性）。布局与 GTS1 相同, 仅 magic 区分。
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
        buf.extend_from_slice(b"GTP1");
        buf.extend_from_slice(&(recs.len() as u64).to_le_bytes());
        for (cid, v, w, c) in recs {
            buf.extend_from_slice(&cid.to_le_bytes());
            buf.extend_from_slice(&v.to_le_bytes());
            buf.extend_from_slice(&w.to_le_bytes());
            buf.extend_from_slice(&c.to_le_bytes());
        }
        std::fs::write(path, buf)?;
        Ok(())
    }

    /// 从落盘文件加载（[`PairSupport::save`] 的逆）。
    pub fn load(path: &Path) -> std::result::Result<PairSupport, CommonError> {
        let data = std::fs::read(path)?;
        if data.len() < 8 || &data[..4] != b"GTP1" {
            return Err(CommonError::Parse(format!(
                "{}: 不是 GTP1 pair 支持文件",
                path.display()
            )));
        }
        let n = u64::from_le_bytes(data[4..12].try_into().unwrap()) as usize;
        let rec_len = 20usize;
        if data.len() != 12 + n * rec_len {
            return Err(CommonError::Parse(format!(
                "{}: GTP1 长度不符（头 {} 字节 + {n}×{rec_len}, 实际 {}）",
                path.display(),
                12,
                data.len()
            )));
        }
        let mut map = FxHashMap::with_capacity_and_hasher(n, Default::default());
        for i in 0..n {
            let b = &data[12 + i * rec_len..12 + (i + 1) * rec_len];
            let cid = u64::from_le_bytes(b[0..8].try_into().unwrap());
            let v = i32::from_le_bytes(b[8..12].try_into().unwrap());
            let w = i32::from_le_bytes(b[12..16].try_into().unwrap());
            let c = u32::from_le_bytes(b[16..20].try_into().unwrap());
            map.insert((cid, pack_uv(v, w)), c);
        }
        Ok(PairSupport { map })
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
    /// visited 穿越（Rust 版扩展, 机制见模块文档）: true = walk 可穿越已
    /// visited 节点（两级候选过滤 + 完全重复产出丢弃）; false = 原型语义
    /// （visited 硬排除, 逐字节复刻 python）。
    pub allow_visited_walk: bool,
    /// 简单 bulge 双臂产出（Rust 版扩展, 机制见模块文档）: true = 前向延伸
    /// 遇简单 bulge 分叉点双臂各走一条、胜出起点的全部 fork 参与产出;
    /// false = 现状（argmax 二选一）。
    pub bulge_dual: bool,
    /// read pair 分支证据乘性消费（`--gtasm-pair-support`; 近似版, 静态
    /// gsc 层面）: true = gsc_e 追加 ×(1+log1p(pair_e)) 因子——pair_e 是
    /// [`PairSupport`] 的分支出边 pair 计数（非分支出边恒 0 → 因子 1,
    /// 不改变行为）; 起点采样与候选 argmax 同吃。false = 不消费（现状）。
    pub pair_boost: bool,
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
///
/// `allow_visited` = visited 穿越（[`GtasmParams::allow_visited_walk`]）:
/// 开时 visited 不再截断 walk——入列检查跳过 visited, 候选分两级, 优先未
/// visited（原型语义）, 一级为空再取已 visited 候选（同元组序取 max）;
/// seen/frozen 两级一致排除（walk 内防环不变）。关时逐字节 = 原型。
fn greedy_extend(
    adj: &FxHashMap<i32, Vec<(f32, i32, usize)>>,
    start: i32,
    visited: &FxHashSet<i32>,
    frozen: &FxHashSet<i32>,
    sup: &[u32],
    allow_visited: bool,
) -> (Vec<i32>, u64) {
    let mut walk = Vec::new();
    let mut seen: FxHashSet<i32> = FxHashSet::default();
    let mut sup_sum = 0u64;
    let mut cur = start;
    loop {
        if frozen.contains(&cur) || seen.contains(&cur) {
            break;
        }
        if !allow_visited && visited.contains(&cur) {
            break; // 原型语义: visited 硬截断（python `cur not in visited`）
        }
        seen.insert(cur);
        walk.push(cur);
        // 两级候选过滤: best = 未 visited 最优（原型）; best_cross = 已 visited
        // 最优（穿越开才收集）; 一级非空时一级优先。
        let mut best: Option<(f32, i32, usize)> = None;
        let mut best_cross: Option<(f32, i32, usize)> = None;
        if let Some(cands) = adj.get(&cur) {
            for &cand in cands {
                if frozen.contains(&cand.1) || seen.contains(&cand.1) {
                    continue;
                }
                if visited.contains(&cand.1) {
                    if allow_visited && better_cand(cand, best_cross) {
                        best_cross = Some(cand);
                    }
                } else if better_cand(cand, best) {
                    best = Some(cand);
                }
            }
        }
        match best.or(best_cross) {
            None => break,
            Some((_, w, e)) => {
                sup_sum += sup[e] as u64;
                cur = w;
            }
        }
    }
    (walk, sup_sum)
}

// ---------------------------------------------------------------------------
// 简单 bulge 双臂产出（--gtasm-bulge-dual）
// ---------------------------------------------------------------------------

/// bulge 判定用的结构邻接（节点级, 平行边去重; 语义同 rnaspades `build_adj`,
/// 在此独立实现以避免跨模块私有依赖）。
struct BulgeAdj {
    /// u → 去重后继（首现序）, 不含 dst<0 哨兵。
    succ: FxHashMap<i32, Vec<i32>>,
    /// 真实前驱数（不含 <0 哨兵 root 边）。
    indeg: FxHashMap<i32, u32>,
}

fn build_bulge_adj(g: &GtGraph) -> BulgeAdj {
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
    BulgeAdj { succ, indeg }
}

/// 简单 bulge 重汇合判定（rnaspades `bulge_merge` 语义）: 从两支各自沿
/// **唯一后继**前进（支内节点出度恒 1, 否则非简单 bulge）, 两支在
/// `max_chain` 节点内首次到达**同一**节点 C 且 C 的入度恰为 2（仅这两支
/// 汇合）→ `Some(C)`; 死端/支内再分叉/超上限 → `None`。C 是否已在当前
/// walk 上由调用方检查（排除集随分叉状态而变, 不在此固化）。
fn bulge_merge_point(adj: &BulgeAdj, b1: i32, b2: i32, max_chain: usize) -> Option<i32> {
    // 返回支链下游首个"汇合候选"节点（入度 != 1 者）; 支内死端/分叉/超限 → None
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
    if c1 == c2 && adj.indeg.get(&c1).copied().unwrap_or(0) == 2 {
        Some(c1)
    } else {
        None
    }
}

/// 前向延伸的一个活跃分叉（bulge_dual 用）。
struct ForkState {
    walk: Vec<i32>,
    seen: FxHashSet<i32>,
    sup_sum: u64,
    cur: i32,
}

/// 延伸共享上下文（邻接/visited/边支持/穿越开关）——收拢 dual 版参数用,
/// 字段语义与 [`greedy_extend`] 同名参数逐点一致。
struct ExtendCtx<'a> {
    adj: &'a FxHashMap<i32, Vec<(f32, i32, usize)>>,
    visited: &'a FxHashSet<i32>,
    sup: &'a [u32],
    allow_visited: bool,
}

/// 带 bulge 分叉的前向贪心延伸（[`greedy_extend`] 的双臂变体）。
///
/// 候选两级过滤与单 walk 版完全一致; 唯一差异: 当前节点的一级（未
/// visited、未 seen/frozen）候选恰好 2 条且指向不同 vid、并经
/// [`bulge_merge_point`] 判定为简单 bulge（汇合点不在本 walk 上）时,
/// 复制当前分叉走第二臂、本分叉走第一臂, 两臂独立延伸到终结。汇合点后
/// 两 fork 共享后缀（各自 seen 集不同, 下游贪心选择一致时才同构——若
/// 下游再遇 bulge 可各自再分叉, 串联 bubble 的多 isoform 由此产生）。
/// 返回全部终结 walk（≥1 条; 首条 = 主分叉）。
fn greedy_extend_dual(
    ctx: &ExtendCtx,
    start: i32,
    frozen: &FxHashSet<i32>,
    badj: &BulgeAdj,
    max_walks: usize,
) -> Vec<(Vec<i32>, u64)> {
    let (adj, visited, sup, allow_visited) = (ctx.adj, ctx.visited, ctx.sup, ctx.allow_visited);
    let mut states = vec![ForkState {
        walk: Vec::new(),
        seen: FxHashSet::default(),
        sup_sum: 0,
        cur: start,
    }];
    let mut done: Vec<(Vec<i32>, u64)> = Vec::new();
    while let Some(mut st) = states.pop() {
        // 沿单链延伸当前分叉, 直至终结或再分叉（与 greedy_extend 逐行同构）
        loop {
            if frozen.contains(&st.cur) || st.seen.contains(&st.cur) {
                break;
            }
            if !allow_visited && visited.contains(&st.cur) {
                break; // 原型语义: visited 硬截断
            }
            st.seen.insert(st.cur);
            st.walk.push(st.cur);
            let mut best: Option<(f32, i32, usize)> = None;
            let mut best_cross: Option<(f32, i32, usize)> = None;
            // 一级（未 visited）候选全量——bulge 判定用（平行边各自计入,
            // "恰 2 候选"按边数, 再要求 vid 不同, 与"分支点"语义一致）
            let mut first: Vec<(f32, i32, usize)> = Vec::new();
            if let Some(cands) = adj.get(&st.cur) {
                for &cand in cands {
                    if frozen.contains(&cand.1) || st.seen.contains(&cand.1) {
                        continue;
                    }
                    if visited.contains(&cand.1) {
                        if allow_visited && better_cand(cand, best_cross) {
                            best_cross = Some(cand);
                        }
                    } else {
                        if better_cand(cand, best) {
                            best = Some(cand);
                        }
                        first.push(cand);
                    }
                }
            }
            // bulge 分叉: 恰 2 个一级候选、vid 不同、构成简单 bulge、
            // 汇合点不在本 walk 上（防环）、分叉配额未超限
            if first.len() == 2
                && first[0].1 != first[1].1
                && done.len() + states.len() + 2 <= max_walks
            {
                let b1 = first[0];
                let b2 = first[1];
                if let Some(c) = bulge_merge_point(badj, b1.1, b2.1, BULGE_MAX_CHAIN) {
                    if !st.seen.contains(&c) && !frozen.contains(&c) {
                        // 复制走 b2（后写克隆）, 本分叉走 b1（邻接序稳定 = 确定性）
                        states.push(ForkState {
                            walk: st.walk.clone(),
                            seen: st.seen.clone(),
                            sup_sum: st.sup_sum + sup[b2.2] as u64,
                            cur: b2.1,
                        });
                        st.sup_sum += sup[b1.2] as u64;
                        st.cur = b1.1;
                        continue;
                    }
                }
            }
            match best.or(best_cross) {
                None => break,
                Some((_, w, e)) => {
                    st.sup_sum += sup[e] as u64;
                    st.cur = w;
                }
            }
        }
        done.push((st.walk, st.sup_sum));
    }
    done
}

/// 已产出序列指纹（重复产出防护; FxHasher 64 位——组件内 contig 至多数千条,
/// 碰撞概率 ~n²/2⁶⁵, 可忽略）。
fn seq_fingerprint(seq: &str) -> u64 {
    use std::hash::Hasher;
    let mut h = rustc_hash::FxHasher::default();
    h.write(seq.as_bytes());
    h.finish()
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
    pair: &PairSupport,
    cid: u64,
    params: &GtasmParams,
) -> Vec<String> {
    let ne = g.n_edges();
    if ne == 0 {
        return Vec::new();
    }
    debug_assert_eq!(prob.len(), ne);

    // gsc_e = prob_e × (1 + log1p(sup_e))（全 f32, 与 python numpy float32 一致）;
    // pair_boost 开时追加 ×(1 + log1p(pair_e))（read pair 分支证据, 见模块文档）
    let sup: Vec<u32> = (0..ne)
        .map(|ei| support.get(cid, g.vids[g.src[ei]], g.vids[g.dst[ei]]))
        .collect();
    let pair_cnt: Vec<u32> = if params.pair_boost {
        (0..ne)
            .map(|ei| pair.get(cid, g.vids[g.src[ei]], g.vids[g.dst[ei]]))
            .collect()
    } else {
        Vec::new()
    };
    let gsc: Vec<f32> = (0..ne)
        .map(|ei| {
            let base = prob[ei] * (1.0f32 + (sup[ei] as f32).ln_1p());
            if params.pair_boost {
                base * (1.0f32 + (pair_cnt[ei] as f32).ln_1p())
            } else {
                base
            }
        })
        .collect();

    // 邻接: vid → [(gsc, 对端 vid, eid)]（eid 序插入; python defaultdict 同序）
    let mut succ: FxHashMap<i32, Vec<(f32, i32, usize)>> = FxHashMap::default();
    let mut pred: FxHashMap<i32, Vec<(f32, i32, usize)>> = FxHashMap::default();
    for (ei, (&s, &d)) in g.src.iter().zip(&g.dst).enumerate() {
        succ.entry(g.vids[s]).or_default().push((gsc[ei], g.vids[d], ei));
        pred.entry(g.vids[d]).or_default().push((gsc[ei], g.vids[s], ei));
    }

    // bulge 结构邻接（bulge_dual 开才建; 前向分叉判定用）
    let badj = if params.bulge_dual {
        Some(build_bulge_adj(g))
    } else {
        None
    };

    // 逐 comp 确定性种子（splitmix64 混合）
    let mut seed_mix = params.seed ^ cid.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    seed_mix = (seed_mix ^ (seed_mix >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    seed_mix = (seed_mix ^ (seed_mix >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    seed_mix ^= seed_mix >> 31;
    let mut rng = XorShift64::new(seed_mix);

    let mut alive: Vec<usize> = (0..ne).collect();
    let mut visited: FxHashSet<i32> = FxHashSet::default();
    let mut contigs: Vec<String> = Vec::new();
    // 重复产出防护（visited 穿越开才启用）: 本组件已产出序列的指纹集。
    let mut emitted: FxHashSet<u64> = FxHashSet::default();

    while !alive.is_empty() {
        let k = NB_STARTS.min(alive.len());
        let picks = sample_starts(&alive, &gsc, k, &mut rng);

        // 胜出起点的全部 walk（bulge_dual 下可多条; 单 walk 模式恒 1 条）
        let mut best_walks: Option<Vec<Vec<i32>>> = None;
        let mut best_key: u64 = 0;
        for e0 in picks {
            let s0 = g.vids[g.src[e0]];
            let d0 = g.vids[g.dst[e0]];

            // 前向: 从 d0 沿出边贪心（bulge_dual 开时 bulge 分叉点双臂各走一条,
            // 返回全部 fork）; 反向: 从 s0 沿入边贪心（frozen = 全部前向 fork 的
            // 节点并集——任一 back+fork 组合都不交叉, 单 fork 时与原语义逐点一致）
            let fwd: Vec<(Vec<i32>, u64)> = match &badj {
                Some(b) => greedy_extend_dual(
                    &ExtendCtx {
                        adj: &succ,
                        visited: &visited,
                        sup: &sup,
                        allow_visited: params.allow_visited_walk,
                    },
                    d0,
                    &FxHashSet::default(),
                    b,
                    MAX_FORKS_PER_START,
                ),
                None => {
                    let (w, s) = greedy_extend(
                        &succ,
                        d0,
                        &visited,
                        &FxHashSet::default(),
                        &sup,
                        params.allow_visited_walk,
                    );
                    vec![(w, s)]
                }
            };
            let mut seen_f: FxHashSet<i32> = FxHashSet::default();
            for (w, _) in &fwd {
                seen_f.extend(w.iter().copied());
            }
            let (mut walk_b, sup_b) =
                greedy_extend(&pred, s0, &visited, &seen_f, &sup, params.allow_visited_walk);

            walk_b.reverse();
            // 组合 back+各 fork: 逐 fork 独立评估 sup_sum 与 ≥2 节点有效性
            let mut start_walks: Vec<(Vec<i32>, u64)> = Vec::with_capacity(fwd.len());
            for (wf, sf) in fwd {
                let mut w = walk_b.clone();
                w.extend_from_slice(&wf);
                if w.len() >= 2 {
                    start_walks.push((w, sf + sup_b)); // best_mode="sup_sum"
                }
            }
            if start_walks.is_empty() {
                continue;
            }
            // 起点 key = 其全部 fork 的 sup_sum 最大值; 首个达到最大者胜
            // （与原 `key > best_key` 严格大于一致）
            let start_key = start_walks.iter().map(|t| t.1).max().unwrap_or(0);
            if best_walks.is_none() || start_key > best_key {
                best_walks = Some(start_walks.into_iter().map(|t| t.0).collect());
                best_key = start_key;
            }
        }

        let Some(win_walks) = best_walks else {
            break; // 无合法 walk → 终止
        };

        // 胜出起点的全部 fork 参与产出（各自独立拼序列/长度阈值/去重）
        let mut hit_cap = false;
        for walk in &win_walks {
            // 拼序列: 首个有名节点全 kmer + 后续节点末字符; 无名节点跳过
            let mut seq = String::new();
            for x in walk {
                if let Some(kk) = g.names.get(x) {
                    if seq.is_empty() {
                        seq.push_str(kk);
                    } else {
                        seq.push(kk.as_bytes()[kk.len() - 1] as char);
                    }
                }
            }
            if seq.len() >= params.min_len {
                // 重复产出防护: 穿越/双臂重新发现的完全相同序列才丢弃（visited
                // 照常并入, 终止语义不受影响）; 不同 isoform 照常产出。两机制皆
                // 关时永不丢弃（原型语义, 短路不求指纹）。
                let dup = (params.allow_visited_walk || params.bulge_dual)
                    && !emitted.insert(seq_fingerprint(&seq));
                if !dup {
                    contigs.push(seq);
                }
            }
            for x in walk {
                visited.insert(*x);
            }
            if contigs.len() >= MAX_COMPONENT_CONTIGS {
                hit_cap = true; // 组件产出安全阀（防异常大组件产出爆炸）
                break;
            }
        }
        alive.retain(|&e| {
            !visited.contains(&g.vids[g.src[e]]) || !visited.contains(&g.vids[g.dst[e]])
        });
        if hit_cap {
            break;
        }
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
    model: &GtModel,
    support: &EdgeSupport,
    pair: &PairSupport,
    params: &GtasmParams,
) -> Result<String> {
    let g = parse_component_graph(graph_text)?;
    if g.n_edges() == 0 {
        return Ok(String::new());
    }
    let (x, e) = build_gt_features(&g);
    let prob = model.forward(&x, &e, &g.src, &g.dst);
    let contigs = decode_component(&g, &prob, support, pair, cid, params);
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
            &PairSupport::new(),
            0,
            &GtasmParams {
                min_len: 20,
                seed: 42,
                allow_visited_walk: false,
                bulge_dual: false,
                pair_boost: false,
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
            &PairSupport::new(),
            0,
            &GtasmParams {
                min_len: 200, // 26bp < 200 → 无输出
                seed: 42,
                allow_visited_walk: false,
                bulge_dual: false,
                pair_boost: false,
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
        let (walk, _) = greedy_extend(&succ, 10, &FxHashSet::default(), &FxHashSet::default(), &sup, false);
        assert_eq!(walk, vec![10, 12]);

        // eid 级平局（平行边）: (0.5, 11, 0) vs (0.5, 11, 3) → eid 大者 3 胜,
        // sup_sum 计入 eid 3 的支持
        let mut succ2: FxHashMap<i32, Vec<(f32, i32, usize)>> = FxHashMap::default();
        succ2.insert(10, vec![(0.5, 11, 0), (0.5, 11, 3)]);
        let sup2 = [7u32, 0, 0, 9];
        let (walk2, sup_sum2) =
            greedy_extend(&succ2, 10, &FxHashSet::default(), &FxHashSet::default(), &sup2, false);
        assert_eq!(walk2, vec![10, 11]);
        assert_eq!(sup_sum2, 9);

        // frozen（前向已走）节点不可复走
        let frozen: FxHashSet<i32> = [12].into_iter().collect();
        let (walk3, _) = greedy_extend(&succ, 10, &FxHashSet::default(), &frozen, &sup, false);
        assert_eq!(walk3, vec![10, 11]);
    }

    #[test]
    fn greedy_extend_visited_crossing_two_levels() {
        // 两级过滤: 有未 visited 候选时优先之（一级）; 一级为空才穿越（二级）。
        let mut succ: FxHashMap<i32, Vec<(f32, i32, usize)>> = FxHashMap::default();
        succ.insert(10, vec![(0.9, 11, 0), (0.1, 12, 1)]);
        let sup = [0u32; 4];
        let visited: FxHashSet<i32> = [11].into_iter().collect();
        // 11 visited 但 gsc 更高: 一级（12 未 visited）优先 → 走 12
        let (walk, _) = greedy_extend(&succ, 10, &visited, &FxHashSet::default(), &sup, true);
        assert_eq!(walk, vec![10, 12]);
        // 唯一候选已 visited: 穿越开 → 走 11; 关 → 迈出前截断（walk 止于起点）
        let mut succ2: FxHashMap<i32, Vec<(f32, i32, usize)>> = FxHashMap::default();
        succ2.insert(10, vec![(0.5, 11, 0)]);
        let (walk2, _) = greedy_extend(&succ2, 10, &visited, &FxHashSet::default(), &sup, true);
        assert_eq!(walk2, vec![10, 11]);
        let (walk3, _) = greedy_extend(&succ2, 10, &visited, &FxHashSet::default(), &sup, false);
        assert_eq!(walk3, vec![10]);
        // 起点自身 visited: 穿越开 → 照常入列延伸; 关 → 立即截断
        let visited2: FxHashSet<i32> = [10].into_iter().collect();
        let (walk4, _) = greedy_extend(&succ2, 10, &visited2, &FxHashSet::default(), &sup, true);
        assert_eq!(walk4, vec![10, 11]);
        let (walk5, _) = greedy_extend(&succ2, 10, &visited2, &FxHashSet::default(), &sup, false);
        assert!(walk5.is_empty());
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
            &PairSupport::new(),
            0,
            &GtasmParams {
                min_len: 20,
                seed: 7,
                allow_visited_walk: false,
                bulge_dual: false,
                pair_boost: false,
            },
        );
        // 两条臂无共享节点 → visited 不交叉, 两轮各出一条
        assert_eq!(contigs.len(), 2);
        // A 臂先出（竞争胜）——顺序即胜出序; B 臂 = 20 的全 24-mer（24×C）
        // + 21 末字符 A + 22 末字符 G
        assert_eq!(contigs[0], "AAAAAAAAAAAAAAAAAAAAAAAAC");
        assert_eq!(contigs[1], "CCCCCCCCCCCCCCCCCCCCCCCCAG");
    }

    // -----------------------------------------------------------------------
    // visited 穿越（--gtasm-visited-walk）
    // -----------------------------------------------------------------------

    /// 泡形双 isoform 图: 共享入口 S=10 / 汇 E=30 / 尾 T=50, 两臂 20-21（强,
    /// sup=5）与 40-41（弱, sup=1）。强臂 A = [-1,10,20,21,30,50], 弱臂
    /// B = [-1,10,40,41,30,50]——两条 walk 共享 {-1,10,30,50}。
    fn bubble_graph() -> (GtGraph, EdgeSupport) {
        let text = format!(
            "Component 0\n{}{}{}{}{}{}{}{}",
            row(10, -1, 3, "AAAAAAAAAAAAAAAAAAAAAAAA"),
            row(20, 10, 5, "AAAAAAAAAAAAAAAAAAAAAAAC"),
            row(21, 20, 5, "AAAAAAAAAAAAAAAAAAAAAAAG"),
            row(30, 21, 5, "AAAAAAAAAAAAAAAAAAAAAAAT"),
            row(50, 30, 3, "CCCCCCCCCCCCCCCCCCCCCCCG"),
            row(40, 10, 1, "CCCCCCCCCCCCCCCCCCCCCCCC"),
            row(41, 40, 1, "CCCCCCCCCCCCCCCCCCCCCCCA"),
            row(30, 41, 1, "AAAAAAAAAAAAAAAAAAAAAAAT"), // 汇的弱臂入口
        );
        let g = parse_component_graph(&text).unwrap();
        let mut sup = EdgeSupport::new();
        for _ in 0..5 {
            sup.add(0, 10, 20);
            sup.add(0, 20, 21);
            sup.add(0, 21, 30);
        }
        sup.add(0, 10, 40);
        sup.add(0, 40, 41);
        sup.add(0, 41, 30);
        (g, sup)
    }

    /// 强/弱臂的期望序列（手算: -1 无名跳过, 首名节点全 24-mer + 后续末字符）。
    const ISO_A: &str = "AAAAAAAAAAAAAAAAAAAAAAAACGTG"; // [-1,10,20,21,30,50]
    const ISO_B: &str = "AAAAAAAAAAAAAAAAAAAAAAAACATG"; // [-1,10,40,41,30,50]

    #[test]
    fn decode_shared_exon_both_isoforms() {
        // (a) 共享节点的两 isoform: 穿越开 → 强弱两臂都完整产出。
        // 竞争手算（prob 全 0.5）: 强臂延伸 sup=5×3, 弱臂 sup=1×3——
        // 轮 1 全采（8 边 ≤ NB_STARTS）,  sup_sum 最大 15 = 强臂全程, A 先出;
        // 轮 2 alive = 弱臂三边, walk 从 40/41 起穿越 visited 的 30/50/10/-1,
        // 完整 B 产出（序列 ≠ A, 去重不误伤）。
        let (g, sup) = bubble_graph();
        let prob = vec![0.5f32; g.n_edges()];
        let contigs = decode_component(
            &g,
            &prob,
            &sup,
            &PairSupport::new(),
            0,
            &GtasmParams {
                min_len: 20,
                seed: 42,
                allow_visited_walk: true,
                bulge_dual: false,
                pair_boost: false,
            },
        );
        assert_eq!(contigs, vec![ISO_A.to_string(), ISO_B.to_string()]);

        // 对照（原型语义）: 轮 2 walk 在 visited 边界截断, 弱臂只剩 [40,41]
        // 的 25bp 残段（40 全 24-mer + 41 末字符）——复现"完整路径被拆散吞没"
        // 的缺陷。
        let contigs_off = decode_component(
            &g,
            &prob,
            &sup,
            &PairSupport::new(),
            0,
            &GtasmParams {
                min_len: 20,
                seed: 42,
                allow_visited_walk: false,
                bulge_dual: false,
                pair_boost: false,
            },
        );
        assert_eq!(contigs_off.len(), 2);
        assert_eq!(contigs_off[0], ISO_A);
        assert_eq!(contigs_off[1], "CCCCCCCCCCCCCCCCCCCCCCCCA"); // 截断残段 ≠ ISO_B
    }

    #[test]
    fn decode_visited_walk_dedups_identical() {
        // (b) 去重: 同一 walk 经穿越重复发现不重复产出。
        // 链 -1→10→11→12 加一条无名平行入口 99→10（99 仅作 from, 无 kmer 名）。
        // 轮 1 从任一入口边走出全程 walk（[-1,...] 或 [99,...], 序列相同 =
        // K1+C+G）; 轮 2 alive 只剩另一条入口边, 穿越 visited 重新发现同一
        // 序列 → 完全相同, 丢弃。最终恰 1 条（与轮 1 采样序无关, 两种镜像
        // 殊途同归）。
        let text = format!(
            "Component 0\n{}{}{}{}",
            row(10, -1, 3, "AAAAAAAAAAAAAAAAAAAAAAAA"),
            row(11, 10, 3, "AAAAAAAAAAAAAAAAAAAAAAAC"),
            row(12, 11, 3, "AAAAAAAAAAAAAAAAAAAAAAAG"),
            row(10, 99, 3, "AAAAAAAAAAAAAAAAAAAAAAAA"), // 99→10 平行入口
        );
        let g = parse_component_graph(&text).unwrap();
        let prob = vec![0.5f32; g.n_edges()];
        let contigs = decode_component(
            &g,
            &prob,
            &EdgeSupport::new(),
            &PairSupport::new(),
            0,
            &GtasmParams {
                min_len: 20,
                seed: 42,
                allow_visited_walk: true,
                bulge_dual: false,
                pair_boost: false,
            },
        );
        assert_eq!(contigs, vec!["AAAAAAAAAAAAAAAAAAAAAAAACG".to_string()]);
    }

    #[test]
    fn decode_visited_walk_terminates_on_cycle() {
        // (c) 终止: 环图（10→11→12→10 + 入口 -1→10）不挂死。
        // walk 内 seen 防环 → 单条 walk 有限; 每轮 visited 只增不减且胜出
        // walk 必含起点边两端 → 至多 ne 轮。
        let text = format!(
            "Component 0\n{}{}{}{}",
            row(10, -1, 3, "AAAAAAAAAAAAAAAAAAAAAAAA"),
            row(11, 10, 3, "AAAAAAAAAAAAAAAAAAAAAAAC"),
            row(12, 11, 3, "AAAAAAAAAAAAAAAAAAAAAAAG"),
            row(10, 12, 3, "AAAAAAAAAAAAAAAAAAAAAAAA"), // 12→10 闭环
        );
        let g = parse_component_graph(&text).unwrap();
        let prob = vec![0.5f32; g.n_edges()];
        let contigs = decode_component(
            &g,
            &prob,
            &EdgeSupport::new(),
            &PairSupport::new(),
            0,
            &GtasmParams {
                min_len: 20,
                seed: 42,
                allow_visited_walk: true,
                bulge_dual: false,
                pair_boost: false,
            },
        );
        // 终止且产出有限: 轮 1 胜出 walk 取决于采样序（三种旋转之一, 序列
        // 均为 26bp）; 若 -1 未覆盖, 轮 2 穿越再出一条（序列不同则保留,
        // 相同则去重）→ 1~2 条, 绝不挂死。
        assert!(!contigs.is_empty() && contigs.len() <= 2, "{contigs:?}");
        assert!(contigs.iter().all(|s| s.len() == 26), "{contigs:?}");
    }

    #[test]
    fn decode_empty_graph() {
        let g = parse_component_graph("Component 0\n").unwrap();
        assert_eq!(g.n_edges(), 0);
        let contigs = decode_component(
            &g,
            &[],
            &EdgeSupport::new(),
            &PairSupport::new(),
            0,
            &GtasmParams {
                min_len: 200,
                seed: 1,
                allow_visited_walk: false,
                bulge_dual: false,
                pair_boost: false,
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
        let fasta = run_component_gtasm(
            &text,
            0,
            "c0",
            &m,
            &EdgeSupport::new(),
            &PairSupport::new(),
            &GtasmParams {
                min_len: 200,
                seed: GTASM_DEFAULT_SEED,
                allow_visited_walk: true,
                bulge_dual: false,
                pair_boost: false,
            },
        )
        .unwrap();
        assert!(fasta.starts_with(">c0_g1_i1 len="), "{fasta}");
        assert!(fasta.contains(" path=[gtasm]\n"));
        let seq: String = fasta.lines().filter(|l| !l.starts_with('>')).collect();
        assert!(seq.len() >= 200);
    }

    // -----------------------------------------------------------------------
    // 简单 bulge 双臂产出（--gtasm-bulge-dual）
    // -----------------------------------------------------------------------

    /// 简单 bulge 图: 入口 -1→10, 分支点 10 → 臂A 20→21 / 臂B 40→41,
    /// 重汇合于 30（indeg 2）, 尾 30→50。臂A 强（sup=5/边）, 臂B 弱（sup=1）。
    /// 单 walk 模式从 10 必选 A 臂（gsc 全等时 sup 不参选, (gsc,vid,eid)
    /// 元组 vid 大者优先 = 40!——故此处用支持差让 gsc 区分: A 臂 gsc 更高）。
    fn bulge_graph() -> (GtGraph, EdgeSupport) {
        let text = format!(
            "Component 0\n{}{}{}{}{}{}{}{}",
            row(10, -1, 3, "AAAAAAAAAAAAAAAAAAAAAAAA"),
            row(20, 10, 5, "AAAAAAAAAAAAAAAAAAAAAAAC"),
            row(21, 20, 5, "AAAAAAAAAAAAAAAAAAAAAAAG"),
            row(30, 21, 5, "AAAAAAAAAAAAAAAAAAAAAAAT"),
            row(50, 30, 3, "CCCCCCCCCCCCCCCCCCCCCCCG"),
            row(40, 10, 1, "CCCCCCCCCCCCCCCCCCCCCCCC"),
            row(41, 40, 1, "CCCCCCCCCCCCCCCCCCCCCCCA"),
            row(30, 41, 1, "AAAAAAAAAAAAAAAAAAAAAAAT"), // 汇的弱臂入口
        );
        let g = parse_component_graph(&text).unwrap();
        let mut sup = EdgeSupport::new();
        // A 臂每边 sup=5 → gsc_A = 0.5×(1+ln6) > gsc_B = 0.5×(1+ln2)
        for _ in 0..5 {
            sup.add(0, 10, 20);
            sup.add(0, 20, 21);
            sup.add(0, 21, 30);
        }
        sup.add(0, 10, 40);
        sup.add(0, 40, 41);
        sup.add(0, 41, 30);
        (g, sup)
    }

    /// 双臂期望序列（手算: -1 无名跳过, 首名节点全 24-mer + 后续末字符）。
    /// walkA = [-1,10,20,21,30,50], walkB = [-1,10,40,41,30,50]。
    const DUAL_A: &str = "AAAAAAAAAAAAAAAAAAAAAAAACGTG";
    const DUAL_B: &str = "AAAAAAAAAAAAAAAAAAAAAAAACATG";

    #[test]
    fn bulge_dual_emits_both_arms() {
        // bulge_dual 开: 分支点 10 遇简单 bulge → 单轮内双臂都产出。
        // 竞争手算: 轮 1 全采（8 边 ≤ NB_STARTS）; 从 A 臂边出发的 walk
        // sup_sum = 5×3 = 15 最大 → 胜出; 其前向延伸在 10 分叉（恰 2 一级
        // 候选 20/40, bulge_merge_point: 20→21→30 与 40→41→30 支内出度恒 1,
        // 同汇 30 且 indeg(30)=2）→ fork 出 [10,20,21,30,50] 与 [10,40,41,30,50],
        // 加上 back 前缀 [-1] → 双臂同轮产出, visited 并集覆盖全图 → 一轮收工。
        let (g, sup) = bulge_graph();
        let prob = vec![0.5f32; g.n_edges()];
        let contigs = decode_component(
            &g,
            &prob,
            &sup,
            &PairSupport::new(),
            0,
            &GtasmParams {
                min_len: 20,
                seed: 42,
                allow_visited_walk: true,
                bulge_dual: true,
                pair_boost: false,
            },
        );
        assert_eq!(
            contigs,
            vec![DUAL_A.to_string(), DUAL_B.to_string()],
            "简单 bulge 应同轮产出双臂"
        );
    }

    #[test]
    fn bulge_dual_off_keeps_single_choice() {
        // 对照: bulge_dual 关 → 分支点 argmax 二选一（gsc 高的 A 臂）;
        // 弱臂 B 轮 2 经穿越才补出（与 decode_shared_exon_both_isoforms 同型,
        // 此处锁定 bulge_dual=false 不改现状）。
        let (g, sup) = bulge_graph();
        let prob = vec![0.5f32; g.n_edges()];
        let contigs = decode_component(
            &g,
            &prob,
            &sup,
            &PairSupport::new(),
            0,
            &GtasmParams {
                min_len: 20,
                seed: 42,
                allow_visited_walk: true,
                bulge_dual: false,
                pair_boost: false,
            },
        );
        assert_eq!(
            contigs,
            vec![DUAL_A.to_string(), DUAL_B.to_string()],
            "穿越开时两轮各出一条（顺序 = 竞争胜出序）"
        );
        // 穿越也关（原型语义）: 弱臂只剩残段——bulge_dual 的价值对照
        let contigs_proto = decode_component(
            &g,
            &prob,
            &sup,
            &PairSupport::new(),
            0,
            &GtasmParams {
                min_len: 20,
                seed: 42,
                allow_visited_walk: false,
                bulge_dual: false,
                pair_boost: false,
            },
        );
        assert_eq!(contigs_proto[1], "CCCCCCCCCCCCCCCCCCCCCCCCA");
    }

    #[test]
    fn bulge_dual_non_bulge_not_forked() {
        // 非 bulge 不误产: 两支分叉后不再汇合（臂B 终点 41 无出边连到 30,
        // 改为死端）→ 分支点 argmax 单选, 弱臂不走; 且轮 2 起点的弱臂 walk
        // 穿越汇合不到主干, 独立产出残段（不翻倍、不死循环）。
        let text = format!(
            "Component 0\n{}{}{}{}{}{}{}",
            row(10, -1, 3, "AAAAAAAAAAAAAAAAAAAAAAAA"),
            row(20, 10, 5, "AAAAAAAAAAAAAAAAAAAAAAAC"),
            row(21, 20, 5, "AAAAAAAAAAAAAAAAAAAAAAAG"),
            row(30, 21, 5, "AAAAAAAAAAAAAAAAAAAAAAAT"),
            row(50, 30, 3, "CCCCCCCCCCCCCCCCCCCCCCCG"),
            row(40, 10, 1, "CCCCCCCCCCCCCCCCCCCCCCCC"),
            row(41, 40, 1, "CCCCCCCCCCCCCCCCCCCCCCCA"), // 41 死端: 无 41→30
        );
        let g = parse_component_graph(&text).unwrap();
        let prob = vec![0.5f32; g.n_edges()];
        let mut sup = EdgeSupport::new();
        for _ in 0..5 {
            sup.add(0, 10, 20);
            sup.add(0, 20, 21);
            sup.add(0, 21, 30);
        }
        sup.add(0, 10, 40);
        sup.add(0, 40, 41);
        let contigs = decode_component(
            &g,
            &prob,
            &sup,
            &PairSupport::new(),
            0,
            &GtasmParams {
                min_len: 20,
                seed: 42,
                allow_visited_walk: true,
                bulge_dual: true,
                pair_boost: false,
            },
        );
        // 恰 2 条: A 全程 + B 臂残段（轮 2 起点 10→40/40→41, 前向 40→41 后
        // 死端; back 穿越 visited 的 10→-1 得 [-1,10,40,41] = K(10)+C+A 26bp）
        assert_eq!(contigs.len(), 2, "{contigs:?}");
        assert_eq!(contigs[0], DUAL_A);
        assert_eq!(contigs[1], "AAAAAAAAAAAAAAAAAAAAAAAACA");
    }

    #[test]
    fn bulge_dual_serial_bulge_terminates() {
        // 串联 bulge + 终止性: bulge1 (10→{20,40}→30) 下游接 bulge2
        // (50→{60,70}→80)。双臂 × 双臂 = 恰 4 fork, 不挂死不爆炸。
        // 每延伸边 sup=1 → 起点 (-1,10) key=6 严格最大（其余起点 ≤5,
        // 手算见下）→ 轮 1 恒由 (-1,10) 胜出, 与采样序无关（确定性）。
        let text = format!(
            "Component 0\n{}{}{}{}{}{}{}{}{}{}{}",
            row(10, -1, 3, "AAAAAAAAAAAAAAAAAAAAAAAA"),
            row(20, 10, 5, "AAAAAAAAAAAAAAAAAAAAAAAC"),
            row(30, 20, 5, "AAAAAAAAAAAAAAAAAAAAAAAG"),
            row(40, 10, 5, "AAAAAAAAAAAAAAAAAAAAAAAT"),
            row(30, 40, 5, "AAAAAAAAAAAAAAAAAAAAAAAG"), // bulge1 汇 30
            row(50, 30, 5, "CCCCCCCCCCCCCCCCCCCCCCCC"),
            row(60, 50, 5, "CCCCCCCCCCCCCCCCCCCCCCCA"),
            row(80, 60, 5, "CCCCCCCCCCCCCCCCCCCCCCCG"),
            row(70, 50, 5, "CCCCCCCCCCCCCCCCCCCCCCCT"),
            row(80, 70, 5, "CCCCCCCCCCCCCCCCCCCCCCCG"), // bulge2 汇 80
            row(90, 80, 5, "GGGGGGGGGGGGGGGGGGGGGGGG"),
        );
        let g = parse_component_graph(&text).unwrap();
        let prob = vec![0.5f32; g.n_edges()];
        let mut sup = EdgeSupport::new();
        for (u, v) in [
            (10, 20),
            (10, 40),
            (20, 30),
            (40, 30),
            (30, 50),
            (50, 60),
            (50, 70),
            (60, 80),
            (70, 80),
            (80, 90),
        ] {
            sup.add(0, u, v);
        }
        let contigs = decode_component(
            &g,
            &prob,
            &sup,
            &PairSupport::new(),
            0,
            &GtasmParams {
                min_len: 20,
                seed: 42,
                allow_visited_walk: true,
                bulge_dual: true,
                pair_boost: false,
            },
        );
        // 4 条全组合 isoform 同轮产出（visited 并集覆盖全图 → 一轮收工）;
        // 序列 = K(10) + [C|T](20/40) + G(30) + C(50) + [A|T](60/70) + G(80)
        // + G(90), 24+6=30bp 四组合
        assert_eq!(contigs.len(), 4, "{contigs:?}");
        let mut sorted = contigs.clone();
        sorted.sort();
        let expect: Vec<String> = ["CGCAGG", "CGCTGG", "TGCAGG", "TGCTGG"]
            .iter()
            .map(|suf| format!("AAAAAAAAAAAAAAAAAAAAAAAA{suf}"))
            .collect();
        assert_eq!(sorted, expect);
    }

    #[test]
    fn bulge_merge_point_hand() {
        // bulge_merge_point 单元语义: 重汇合/死端/支内再分叉/超限。
        let (g, _sup) = bulge_graph();
        let badj = build_bulge_adj(&g);
        // 正常汇合: (20,40) → 30
        assert_eq!(bulge_merge_point(&badj, 20, 40, 50), Some(30));
        // 超限: max_chain=1 → 支内第一步 21/41 入度皆 1, 超上限 → None
        assert_eq!(bulge_merge_point(&badj, 20, 40, 1), None);
        // 死端: 50 无后继 → None
        assert_eq!(bulge_merge_point(&badj, 50, 20, 50), None);
        // 支内再分叉: 10 出度 2 → follow(10) 立即 None
        assert_eq!(bulge_merge_point(&badj, 10, 20, 50), None);
    }

    // -----------------------------------------------------------------------
    // read pair 分支证据（--gtasm-pair-support, 乘性）
    // -----------------------------------------------------------------------

    #[test]
    fn pair_support_save_load_roundtrip() {
        let mut s = PairSupport::new();
        s.add(7, 10, 20);
        s.add(7, 10, 20);
        s.add(7, 10, 40);
        s.add(3, 12, 20);
        assert_eq!(s.get(7, 10, 20), 2);
        assert_eq!(s.get(7, 10, 40), 1);
        assert_eq!(s.get(7, 20, 10), 0); // 方向敏感
        assert_eq!(s.get(99, 1, 2), 0);
        let dir = std::env::temp_dir().join(format!("gtp_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("pair.bin");
        s.save(&p).unwrap();
        let s2 = PairSupport::load(&p).unwrap();
        assert_eq!(s2.get(7, 10, 20), 2);
        assert_eq!(s2.get(3, 12, 20), 1);
        assert_eq!(s2.len(), 3);
        let _ = std::fs::remove_dir_all(&dir);
        // 坏 magic / 不存在
        assert!(PairSupport::load(Path::new("/nonexistent/x.bin")).is_err());
        let mut bad = EdgeSupport::new();
        bad.add(1, 2, 3);
        let bp = dir.join("e.bin");
        std::fs::create_dir_all(&dir).unwrap();
        bad.save(&bp).unwrap();
        assert!(PairSupport::load(&bp).is_err(), "GTS1 不可误读为 GTP1");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn pair_boost_flips_branch_choice() {
        // 分支出边 pair 计数乘进 gsc → 分支点选择翻转。
        // 图: -1→10, 10→{20,40} 双臂死端; 两臂 sup 皆为 5（gsc 全等）
        // → 无 pair 时 (gsc,vid,eid) 元组 vid 大者 40 胜（轮 1 出臂40,
        // 轮 2 穿越出臂20）; pair(10→20)=3 乘性 boost 后 gsc(10→20) 严格
        // 更大 → 轮 1 出臂20。起点 (-1,10) 的 key=5 严格大于臂内起点
        // (key=0) → 轮 1 恒由 (-1,10) 胜出, 与采样序无关（确定性）。
        let text = format!(
            "Component 0\n{}{}{}",
            row(10, -1, 3, "AAAAAAAAAAAAAAAAAAAAAAAA"),
            row(20, 10, 5, "AAAAAAAAAAAAAAAAAAAAAAAC"),
            row(40, 10, 5, "AAAAAAAAAAAAAAAAAAAAAAAT"),
        );
        let g = parse_component_graph(&text).unwrap();
        let prob = vec![0.5f32; g.n_edges()];
        let mut sup = EdgeSupport::new();
        for _ in 0..5 {
            sup.add(0, 10, 20);
            sup.add(0, 10, 40);
        }
        let arm20 = "AAAAAAAAAAAAAAAAAAAAAAAAC"; // [-1,10,20]: K(10)+C
        let arm40 = "AAAAAAAAAAAAAAAAAAAAAAAAT"; // [-1,10,40]: K(10)+T
        // 无 pair（pair_boost 开但表空 = 不翻）: 臂40 先出
        let contigs = decode_component(
            &g,
            &prob,
            &sup,
            &PairSupport::new(),
            0,
            &GtasmParams {
                min_len: 20,
                seed: 42,
                allow_visited_walk: true,
                bulge_dual: false,
                pair_boost: true,
            },
        );
        assert_eq!(contigs, vec![arm40.to_string(), arm20.to_string()]);
        // pair(10→20)=3 → gsc(10→20) 严格大 → 臂20 先出
        let mut pair = PairSupport::new();
        for _ in 0..3 {
            pair.add(0, 10, 20);
        }
        let contigs2 = decode_component(
            &g,
            &prob,
            &sup,
            &pair,
            0,
            &GtasmParams {
                min_len: 20,
                seed: 42,
                allow_visited_walk: true,
                bulge_dual: false,
                pair_boost: true,
            },
        );
        assert_eq!(contigs2, vec![arm20.to_string(), arm40.to_string()]);
        // pair_boost 关 → 同一张 pair 表不生效（逐点 = 现状）
        let contigs3 = decode_component(
            &g,
            &prob,
            &sup,
            &pair,
            0,
            &GtasmParams {
                min_len: 20,
                seed: 42,
                allow_visited_walk: true,
                bulge_dual: false,
                pair_boost: false,
            },
        );
        assert_eq!(contigs3, vec![arm40.to_string(), arm20.to_string()]);
    }
}
