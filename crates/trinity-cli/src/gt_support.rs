//! gtasm 边支持 pass——`calc_edge_support2.py` 的产品化（逐语义复刻）。
//!
//! 位置: chrysalis 之后、butterfly 池之前（仅 `--path-engine gtasm` 启用）。
//! 输入: Component_bins 的 graphs.pack（建 k2v）+ **原始 reads**
//! （`--left/--right/--single`, 归一化前——归一化 reads 深度被稀释,
//! 原型已证无效）。产物: 每 comp (u,v)→count 的 [`EdgeSupport`] 表,
//! 落盘 `gtasm_edge_support.bin`（checkpoint 复跑直接加载, 不重扫 reads）。
//!
//! 与 python 的逐点对应:
//! - k2v: 全部 comp 的 24-mer → (cid, vid), 正反义双写（revcomp 同键）,
//!   先写先赢（python setdefault）; 实现用 (2bit 值, 长度) 作键（等价字符串键,
//!   防前导 G 跨长度碰撞）。
//! - fastq 流式: **line_no % 4 语义**（python 原样, 含多行序列记录的收集与
//!   文件尾 flush）; 每 read 大写化后逐窗 k=K（K = 图 kmer 长, 原型恒 24;
//!   python `range(len(seq)-K)`——**不含最后一窗**, 原型 quirk 保留）。
//! - 相邻命中同 comp 且 (prev_vid, vid) 是该 comp 的真实图边 → count += 1;
//!   未命中重置 prev。python 会对非边转移也计数, 但解码期只查图边——
//!   计数时按边过滤对结果等价（内存小一个量级）。
//! - python 的 revcomp 回退查询在此省略: k2v 已正反义双写, 回退永不新增命中。
//!
//! ## read pair 分支支持 pass（`--gtasm-pair-support` 配套, 近似静态版）
//!
//! 单 read 76bp 够不着 300-400bp 的分支间距; pair（insert ~300-400bp）是
//! 够得着的证据源。规则（逐规格）: R1/R2 按 read 名配对（同序流式
//! lockstep）, 两条 mate 各自线程化（沿用同一 k2v, K=图 kmer 长）; 对每个
//! 分支出边 e=(v,w)（v 出度≥2, 按去重后继计）, 若一条 mate 的节点路径
//! 包含 e（相邻命中、真实图边）且另一条命中 "v 的上游窗口节点集"（自 v
//! 反向 BFS ≤ [`PAIR_WINDOW_EDGES`] 边的节点集合, 不含 v 自身）, 则
//! pair_support[(cid,v,w)] += 1。每对每边至多计 1（双取向各评一次,
//! 同对同边去重）。产物落盘 `gtasm_pair_support.bin`（GTP1 sidecar——
//! 与 GTS1 记录布局同构、magic 区分, 与单 read 边支持绝不混用同一字段）。
//!
//! 与边支持 pass 的差异（无 python 原型约束处, 从简）: fastq 按严格 4 行
//! 记录（不要 line%4 quirk）; kmer 窗口全扫（不要末窗 quirk）。

use std::fs;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use rustc_hash::{FxHashMap, FxHashSet};

use trinity_butterfly::{parse_component_graph, EdgeSupport, PairSupport};
use trinity_chrysalis::pack::{read_at, PackIndex, GRAPHS_PACK, KIND_GRAPH_OUT};
use trinity_common::error::CommonError;
use trinity_common::fasta::FastaReader;
use trinity_common::io_util::open_maybe_gz;
use trinity_common::kmer::{kmer_to_intval, revcomp_val, KmerId};

use crate::args::SeqType;

/// 支持 pass 统计（日志用）。
#[derive(Debug, Default, Clone, Copy)]
pub struct SupportStats {
    pub reads: u64,
    /// 有支持转移计数（python st["pairs"]——此处只含真实图边上的转移）。
    pub pairs: u64,
    pub edges_with_support: usize,
}

/// k2v 键: (kmer 2bit 值, kmer 长度); 值: (cid, vid)。
/// 键值类型 = KmerId（u128, 2026-08-30 扩容随动）: 本 pass 的 K 取自
/// chrysalis 图（原型恒 24）, 不随 --multi-k 变化（多 k 轮次各自建图,
/// gtasm 仅见最终轮图; multi-k×gtasm 默认被参数层回退为单 k——gtasm 只在
/// k=25 图上验证过, 显式 --multi-k 可强制组合, 见 args::resolve_multi_k）;
/// KmerId 域内任意 k ≤ 64 的图均可正确编码, 对 k>32 无隐性截断依赖。
type K2V = FxHashMap<(KmerId, u32), (u64, i32)>;

/// 每 comp 的真实图边集（pack_uv 键）——计数过滤用（见模块文档）。
type EdgeSets = FxHashMap<u64, FxHashSet<u64>>;

/// 分支点上游窗口的 BFS 深度上限（边数; pair pass 规格: 自 v 反向 BFS ≤15 边）。
pub const PAIR_WINDOW_EDGES: usize = 15;

/// 单窗口节点数安全阀（防汇合密集区内存爆炸; 超限截断为下近似, 实测远不到）。
const PAIR_WINDOW_MAX_NODES: usize = 1024;

/// 分支点上游窗口表: (cid, 分支节点 v) → 反向 BFS ≤ [`PAIR_WINDOW_EDGES`] 边
/// 的节点集（不含 v 自身; <0 哨兵不入——k2v 无其名, mate 命中不可能）。
/// 键只在 v 的去重后继出度 ≥ 2（真分支点）时存在——计数期据此判分支出边。
type BranchWindows = FxHashMap<(u64, i32), FxHashSet<i32>>;

/// 跑整个支持 pass: 建 k2v → 扫 reads → 落盘。返回统计。
pub fn run_support_pass(
    listing: &[(u64, PathBuf)],
    read_files: &[PathBuf],
    seq_type: SeqType,
    out_bin: &Path,
) -> Result<SupportStats, CommonError> {
    // ---- 1. 逐 comp 读 graph.out 记录, 建 k2v 与边集 ----
    let (k2v, edge_sets, k, _windows) = build_k2v(listing)?;
    eprintln!("-gtasm support: k2v {} 键 (K={k}), {} comps 边集", k2v.len(), edge_sets.len());

    // ---- 2. 扫原始 reads（逐文件并行, 末尾合并）----
    let mut stats = SupportStats::default();
    let merged = match read_files {
        [] => EdgeSupport::new(),
        [single] => {
            let (s, st) = scan_reads_file(single, seq_type, k, &k2v, &edge_sets)?;
            stats.reads += st.reads;
            stats.pairs += st.pairs;
            s
        }
        files => {
            let mut results: Vec<Result<(EdgeSupport, SupportStats), CommonError>> = Vec::new();
            std::thread::scope(|scope| {
                let handles: Vec<_> = files
                    .iter()
                    .map(|f| scope.spawn(|| scan_reads_file(f, seq_type, k, &k2v, &edge_sets)))
                    .collect();
                for h in handles {
                    results.push(h.join().unwrap_or_else(|_| {
                        Err(CommonError::Parse("gtasm support 扫描线程 panic".to_string()))
                    }));
                }
            });
            let mut merged = EdgeSupport::new();
            for r in results {
                let (s, st) = r?;
                stats.reads += st.reads;
                stats.pairs += st.pairs;
                merged.merge(&s);
            }
            merged
        }
    };
    stats.edges_with_support = merged.len();
    merged.save(out_bin)?;
    Ok(stats)
}

/// 从 graphs.pack 读全部 comp 的 graph.out → (k2v, 边集, K, 分支窗口)。
/// K = 首个 kmer 的长度（原型恒 24; 图内 kmer 长度一致由 butterfly 解析保证）。
/// 分支窗口只被 pair pass 消费, 但与图解析同源, 一并构建以免重读 pack。
///
/// 并行（性能）: listing 按连续同 bin 分组, 组内开一次 pack/idx 句柄,
/// rayon 组建局部表; 合并严格按组序（= listing 升序）→ k2v"先写先赢"
/// 语义与串行版逐点等价。
fn build_k2v(
    listing: &[(u64, PathBuf)],
) -> Result<(K2V, EdgeSets, usize, BranchWindows), CommonError> {
    use rayon::prelude::*;
    /// build_k2v 并行组的局部产物类型（clippy type_complexity 拆分）。
    type K2vPart = Result<(K2V, EdgeSets, BranchWindows, Option<usize>), CommonError>;
    // 连续同 dir 分组（listing 按 id 升序 → 同 Cbin 连续）
    let mut groups: Vec<&[(u64, PathBuf)]> = Vec::new();
    let mut start = 0;
    for i in 1..=listing.len() {
        let boundary = i == listing.len()
            || listing[i].1.parent() != listing[start].1.parent();
        if boundary {
            groups.push(&listing[start..i]);
            start = i;
        }
    }
    let parts: Vec<K2vPart> = groups
            .par_iter()
            .map(|grp| -> Result<(K2V, EdgeSets, BranchWindows, Option<usize>), CommonError> {
                let dir = grp[0]
                    .1
                    .parent()
                    .map(|d| d.to_path_buf())
                    .unwrap_or_default();
                let idx = PackIndex::load(&dir.join(format!("{GRAPHS_PACK}.idx")))?;
                let file = fs::File::open(dir.join(format!("{GRAPHS_PACK}.pack")))?;
                let mut k2v: K2V = FxHashMap::default();
                let mut edge_sets: EdgeSets = FxHashMap::default();
                let mut windows: BranchWindows = FxHashMap::default();
                let mut k_len: Option<usize> = None;
                for (id, _base) in *grp {
                    let (offset, len) = idx.find(*id, KIND_GRAPH_OUT).ok_or_else(|| {
                        CommonError::Parse(format!(
                            "graphs pack {} 缺记录 (id={id}, kind=GRAPH_OUT)",
                            dir.display()
                        ))
                    })?;
                    let bytes = read_at(&file, offset, len)?;
                    let text = String::from_utf8(bytes).map_err(|e| {
                        CommonError::Parse(format!("comp {id} graph.out 非 UTF-8: {e}"))
                    })?;
                    let g = parse_component_graph(&text).map_err(|e| {
                        CommonError::Parse(format!("comp {id} graph.out 解析失败: {e}"))
                    })?;
                    // k2v: 全部 names（含无边 comp——python 同; 先写先赢）
                    for (&vid, kmer) in &g.names {
                        let k = kmer.len();
                        let val = kmer_to_intval(kmer.as_bytes()).map_err(|e| {
                            CommonError::Parse(format!("comp {id} kmer 编码失败: {e}"))
                        })?;
                        k2v.entry((val, k as u32)).or_insert((*id, vid));
                        let rc = revcomp_val(val, k);
                        k2v.entry((rc, k as u32)).or_insert((*id, vid));
                        k_len = k_len.or(Some(k));
                    }
                    // 真实图边集
                    let set = edge_sets.entry(*id).or_default();
                    for ei in 0..g.n_edges() {
                        set.insert(trinity_butterfly::gtdecode::pack_uv(
                            g.vids[g.src[ei]],
                            g.vids[g.dst[ei]],
                        ));
                    }
                    // 分支点上游窗口（pair pass 用）: 去重后继 ≥2 的节点 v →
                    // 反向 BFS ≤ PAIR_WINDOW_EDGES 边的节点集
                    let mut succ: FxHashMap<i32, Vec<i32>> = FxHashMap::default();
                    let mut pred: FxHashMap<i32, Vec<i32>> = FxHashMap::default();
                    for ei in 0..g.n_edges() {
                        let (s, d) = (g.vids[g.src[ei]], g.vids[g.dst[ei]]);
                        if d < 0 {
                            continue; // 防御: 哨兵不作路径节点
                        }
                        let e = succ.entry(s).or_default();
                        if !e.contains(&d) {
                            e.push(d);
                        }
                        if s >= 0 {
                            let p = pred.entry(d).or_default();
                            if !p.contains(&s) {
                                p.push(s);
                            }
                        }
                    }
                    for (v, nexts) in &succ {
                        if nexts.len() >= 2 {
                            windows.insert(
                                (*id, *v),
                                upstream_window(
                                    &pred,
                                    *v,
                                    PAIR_WINDOW_EDGES,
                                    PAIR_WINDOW_MAX_NODES,
                                ),
                            );
                        }
                    }
                }
                Ok((k2v, edge_sets, windows, k_len))
            })
            .collect();
    // 按组序（= listing 升序）合并 → k2v 先写先赢
    let mut k2v: K2V = FxHashMap::default();
    let mut edge_sets: EdgeSets = FxHashMap::default();
    let mut windows: BranchWindows = FxHashMap::default();
    let mut k_len: Option<usize> = None;
    for part in parts {
        let (lk, le, lw, llen) = part?;
        for (key, val) in lk {
            k2v.entry(key).or_insert(val);
        }
        edge_sets.extend(le);
        windows.extend(lw);
        k_len = k_len.or(llen);
    }
    Ok((k2v, edge_sets, k_len.unwrap_or(24), windows))
}

/// 自 v 反向 BFS ≤ max_depth 边的上游节点集（不含 v 自身; seen 防环——
/// 有环分量也不会重复入集; 节点数超 cap 截断, 截断为下近似）。
fn upstream_window(
    pred: &FxHashMap<i32, Vec<i32>>,
    v: i32,
    max_depth: usize,
    cap: usize,
) -> FxHashSet<i32> {
    let mut seen: FxHashSet<i32> = [v].into_iter().collect();
    let mut out: FxHashSet<i32> = FxHashSet::default();
    let mut frontier = vec![v];
    for _ in 0..max_depth {
        let mut next = Vec::new();
        for &x in &frontier {
            if let Some(ps) = pred.get(&x) {
                for &p in ps {
                    if seen.insert(p) {
                        out.insert(p);
                        next.push(p);
                        if out.len() >= cap {
                            return out;
                        }
                    }
                }
            }
        }
        if next.is_empty() {
            break;
        }
        frontier = next;
    }
    out
}

/// 单序列窗口扫描（python `process` 语义, 含"不含最后一窗"quirk）。
/// 滚动 2-bit 编码（性能）: 与逐窗 `kmer_to_intval` 逐位等价——
/// 编码表 G=0,A=1,T=2,C=3（base_to_int 同表, 小写同受）, 窗口含非 GATC 时
/// 原版 Err→未命中, 此处 run<k 不查表; 两条路径的 prev 重置行为一致。
fn process_seq(
    seq: &[u8],
    k: usize,
    k2v: &K2V,
    edge_sets: &EdgeSets,
    sup: &mut EdgeSupport,
    pairs: &mut u64,
) {
    if seq.len() <= k {
        return;
    }
    let k32 = k as u32;
    let mask: KmerId = KmerId::low_bits(2 * k32);
    let mut val: KmerId = KmerId::ZERO;
    let mut run = 0usize; // 连续 GATC 长度（>=k 时当前窗有效）
    let mut prev: Option<(u64, i32)> = None;
    // python: for i in range(len(seq) - K)——窗起点 i ∈ [0, len-k-1];
    // 滚动指针 j = i+k-1（窗末字符）, j ∈ [k-1, len-2], 末字符不进循环
    for (j, &c) in seq[..seq.len() - 1].iter().enumerate() {
        match trinity_common::kmer::base_to_int(c) {
            Some(v) => {
                val = ((val << 2) | KmerId::from(v)) & mask;
                run += 1;
            }
            None => {
                val = KmerId::ZERO;
                run = 0;
            }
        }
        if j + 1 < k {
            continue;
        }
        let hit = if run >= k {
            k2v.get(&(val, k32))
        } else {
            None
        };
        match hit {
            None => prev = None,
            Some(&(cid, vid)) => {
                if let Some((pcid, pvid)) = prev {
                    if pcid == cid {
                        if let Some(es) = edge_sets.get(&cid) {
                            if es.contains(&trinity_butterfly::gtdecode::pack_uv(pvid, vid)) {
                                sup.add(cid, pvid, vid);
                                *pairs += 1;
                            }
                        }
                    }
                }
                prev = Some((cid, vid));
            }
        }
    }
}

/// 扫单个 reads 文件（fq: line%4 语义; fa: FASTA 记录）。
/// 流水线并行（性能）: 读取线程按原语义逐条产出完整序列（含 fq 多行收集/
/// 尾 flush quirk）打包入 channel; 消费端 rayon 批内并行 process_seq,
/// EdgeSupport::merge 计数可交换可结合, 合并序无关。
fn scan_reads_file(
    path: &Path,
    seq_type: SeqType,
    k: usize,
    k2v: &K2V,
    edge_sets: &EdgeSets,
) -> Result<(EdgeSupport, SupportStats), CommonError> {
    const BATCH: usize = 8192;
    let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<Vec<u8>>>(8);
    let path_owned = path.to_path_buf();
    let reader_handle = std::thread::spawn(move || -> Result<u64, CommonError> {
        let mut n_reads = 0u64;
        let mut batch: Vec<Vec<u8>> = Vec::with_capacity(BATCH);
        let mut reader = open_maybe_gz(&path_owned)?;
        match seq_type {
            SeqType::Fq => {
                // python calc_edge_support2 的 line_no%4 流式（含多行序列收集与
                // 文件尾 flush; 空行也计数行号, 不跳过）
                let mut line = String::new();
                let mut seq_buf: Vec<u8> = Vec::new();
                let mut line_no = 0u64;
                loop {
                    line.clear();
                    let n = reader.read_line(&mut line)?;
                    if n == 0 {
                        break;
                    }
                    match line_no % 4 {
                        0 => {
                            if !seq_buf.is_empty() {
                                batch.push(std::mem::take(&mut seq_buf));
                                if batch.len() >= BATCH
                                    && tx.send(std::mem::take(&mut batch)).is_err()
                                {
                                    return Ok(n_reads); // 消费端已挂
                                }
                            }
                            n_reads += 1;
                        }
                        1 => seq_buf.extend_from_slice(line.trim().as_bytes()),
                        _ => {}
                    }
                    line_no += 1;
                }
                if !seq_buf.is_empty() {
                    batch.push(seq_buf);
                }
            }
            SeqType::Fa => {
                let mut fr = FastaReader::new(reader);
                while let Some(rec) = fr.next_record()? {
                    n_reads += 1;
                    batch.push(rec.sequence.into_bytes());
                    if batch.len() >= BATCH && tx.send(std::mem::take(&mut batch)).is_err() {
                        return Ok(n_reads);
                    }
                }
            }
        }
        if !batch.is_empty() {
            let _ = tx.send(batch);
        }
        Ok(n_reads)
    });

    let mut sup = EdgeSupport::new();
    let mut pairs = 0u64;
    use rayon::prelude::*;
    while let Ok(batch) = rx.recv() {
        let (part_sup, part_pairs) = batch
            .par_iter()
            .fold(
                || (EdgeSupport::new(), 0u64),
                |mut acc, seq| {
                    process_seq(seq, k, k2v, edge_sets, &mut acc.0, &mut acc.1);
                    acc
                },
            )
            .reduce(
                || (EdgeSupport::new(), 0u64),
                |mut a, b| {
                    a.0.merge(&b.0);
                    a.1 += b.1;
                    a
                },
            );
        sup.merge(&part_sup);
        pairs += part_pairs;
    }
    let reads = reader_handle
        .join()
        .unwrap_or_else(|_| Err(CommonError::Parse("gtasm support 读取线程 panic".into())))?;
    Ok((
        sup,
        SupportStats {
            reads,
            pairs,
            edges_with_support: 0, // 由调用方以 merged.len() 回填
        },
    ))
}

// ---------------------------------------------------------------------------
// read pair 分支支持 pass（机制见模块文档）
// ---------------------------------------------------------------------------

/// pair 支持 pass 统计（日志用）。
#[derive(Debug, Default, Clone, Copy)]
pub struct PairSupportStats {
    /// 输入 read 对数。
    pub pairs: u64,
    /// 跨分支证据计数（每对每边去重后的总增量）。
    pub branch_hits: u64,
    /// 有 pair 支持的分支出边数（表大小）。
    pub edges_with_support: usize,
    /// lockstep 名不匹配的对数（>0 = 两文件可能错序, 警告用）。
    pub name_mismatches: u64,
}

/// 跑整个 pair 支持 pass: 建 k2v/边集/分支窗口 → 双端文件逐对同序流式
/// （逐文件对并行, 末尾合并）→ 落盘 GTP1。返回统计。
pub fn run_pair_support_pass(
    listing: &[(u64, PathBuf)],
    left_files: &[PathBuf],
    right_files: &[PathBuf],
    seq_type: SeqType,
    out_bin: &Path,
) -> Result<PairSupportStats, CommonError> {
    if left_files.len() != right_files.len() {
        return Err(CommonError::Parse(format!(
            "gtasm pair support: --left（{} 个）与 --right（{} 个）文件数不齐",
            left_files.len(),
            right_files.len()
        )));
    }
    let (k2v, edge_sets, k, windows) = build_k2v(listing)?;
    eprintln!(
        "-gtasm pair support: k2v {} 键 (K={k}), {} 分支点窗口",
        k2v.len(),
        windows.len()
    );

    let mut stats = PairSupportStats::default();
    let mut results: Vec<Result<(PairSupport, PairSupportStats), CommonError>> = Vec::new();
    std::thread::scope(|scope| {
        let handles: Vec<_> = left_files
            .iter()
            .zip(right_files.iter())
            .map(|(l, r)| {
                scope.spawn(|| scan_mate_pair_files(l, r, seq_type, k, &k2v, &edge_sets, &windows))
            })
            .collect();
        for h in handles {
            results.push(h.join().unwrap_or_else(|_| {
                Err(CommonError::Parse(
                    "gtasm pair support 扫描线程 panic".to_string(),
                ))
            }));
        }
    });
    let mut merged = PairSupport::new();
    for r in results {
        let (s, st) = r?;
        stats.pairs += st.pairs;
        stats.branch_hits += st.branch_hits;
        stats.name_mismatches += st.name_mismatches;
        merged.merge(&s);
    }
    stats.edges_with_support = merged.len();
    merged.save(out_bin)?;
    Ok(stats)
}

/// 单序列线程化（pair pass 用）: 逐窗 k2v 查询, 返回有序命中列
/// (窗口位置, cid, vid); 未命中的窗天然造成位置断档（等价 process_seq 的
/// prev 重置语义）。**全窗扫描**（pair pass 无 python 原型约束, 不保留
/// 末窗 quirk）。
fn thread_read(seq: &[u8], k: usize, k2v: &K2V) -> Vec<(u32, u64, i32)> {
    let mut hits = Vec::new();
    if seq.len() < k {
        return hits;
    }
    let k32 = k as u32;
    for i in 0..=(seq.len() - k) {
        if let Some(&(cid, vid)) = kmer_to_intval(&seq[i..i + k])
            .ok()
            .and_then(|val| k2v.get(&(val, k32)))
        {
            hits.push((i as u32, cid, vid));
        }
    }
    hits
}

/// 单对 mates 的分支证据计数（规格见模块文档）: 双取向各评一次——取向一
/// X 的节点路径含分支出边 (v,w)（位置相邻的同 comp 命中、真实图边、v 是
/// 分支点）且 Y 命中 v 的上游窗口 → 计入; 每对每边（双取向合并）至多 1。
/// 返回本对的增量计数（统计用）。
fn process_pair(
    xa: &[(u32, u64, i32)],
    xb: &[(u32, u64, i32)],
    edge_sets: &EdgeSets,
    windows: &BranchWindows,
    out: &mut PairSupport,
) -> u64 {
    let mut counted: FxHashSet<(u64, i32, i32)> = FxHashSet::default();
    for (x, y) in [(xa, xb), (xb, xa)] {
        for w in x.windows(2) {
            let (p0, c0, v0) = w[0];
            let (p1, c1, v1) = w[1];
            if c0 != c1 || p1 != p0 + 1 {
                continue; // 非同 comp 或位置断档（中间有未命中窗）→ 非相邻转移
            }
            let Some(win) = windows.get(&(c0, v0)) else {
                continue; // v0 非分支点
            };
            let Some(es) = edge_sets.get(&c0) else {
                continue;
            };
            if !es.contains(&trinity_butterfly::gtdecode::pack_uv(v0, v1)) {
                continue; // 非真实图边
            }
            if y.iter().any(|&(_, c, vid)| c == c0 && win.contains(&vid)) {
                counted.insert((c0, v0, v1));
            }
        }
    }
    let n = counted.len() as u64;
    for &(c, v, w) in &counted {
        out.add(c, v, w);
    }
    n
}

/// fastq 严格 4 行记录读取（pair pass 从简, 无 line%4 quirk）:
/// Ok(Some(四行)) = 一条记录; Ok(None) = 干净 EOF; Err = 截断。
pub(crate) fn read_fq_record(reader: &mut impl BufRead, path: &Path) -> Result<Option<[String; 4]>, CommonError> {
    let mut rec: [String; 4] = Default::default();
    let mut got = 0usize;
    for slot in rec.iter_mut() {
        slot.clear();
        if reader.read_line(slot)? == 0 {
            break;
        }
        got += 1;
    }
    match got {
        0 => Ok(None),
        4 => Ok(Some(rec)),
        _ => Err(CommonError::Parse(format!(
            "{}: fastq 记录截断（{got}/4 行）",
            path.display()
        ))),
    }
}

/// read 名归一化（配对校验用）: 首个空白前截断, 再去尾部 /1、/2。
pub(crate) fn norm_read_name(header: &str) -> &str {
    let t = header.split_whitespace().next().unwrap_or("");
    t.strip_suffix("/1").or_else(|| t.strip_suffix("/2")).unwrap_or(t)
}

/// 扫一对 mates 文件（fq 严格 4 行 lockstep / fa FASTA lockstep, 同序配对;
 /// 首条记录名不匹配计数, 末尾统计告警）。返回 (pair 支持表, 统计)。
fn scan_mate_pair_files(
    left: &Path,
    right: &Path,
    seq_type: SeqType,
    k: usize,
    k2v: &K2V,
    edge_sets: &EdgeSets,
    windows: &BranchWindows,
) -> Result<(PairSupport, PairSupportStats), CommonError> {
    let mut sup = PairSupport::new();
    let mut stats = PairSupportStats::default();
    let mut ra = open_maybe_gz(left)?;
    let mut rb = open_maybe_gz(right)?;
    match seq_type {
        SeqType::Fq => loop {
            let rec_a = read_fq_record(&mut ra, left)?;
            let rec_b = read_fq_record(&mut rb, right)?;
            match (rec_a, rec_b) {
                (Some(a), Some(b)) => {
                    if norm_read_name(&a[0]) != norm_read_name(&b[0]) {
                        stats.name_mismatches += 1;
                    }
                    stats.pairs += 1;
                    let xa = thread_read(a[1].trim().as_bytes(), k, k2v);
                    let xb = thread_read(b[1].trim().as_bytes(), k, k2v);
                    stats.branch_hits += process_pair(&xa, &xb, edge_sets, windows, &mut sup);
                }
                (None, None) => break,
                _ => {
                    return Err(CommonError::Parse(format!(
                        "gtasm pair support: {} 与 {} 记录数不齐（截断/错配）",
                        left.display(),
                        right.display()
                    )))
                }
            }
        },
        SeqType::Fa => {
            let mut fa = FastaReader::new(ra);
            let mut fb = FastaReader::new(rb);
            loop {
                let rec_a = fa.next_record()?;
                let rec_b = fb.next_record()?;
                match (rec_a, rec_b) {
                    (Some(a), Some(b)) => {
                        if norm_read_name(&a.header) != norm_read_name(&b.header) {
                            stats.name_mismatches += 1;
                        }
                        stats.pairs += 1;
                        let xa = thread_read(a.sequence.as_bytes(), k, k2v);
                        let xb = thread_read(b.sequence.as_bytes(), k, k2v);
                        stats.branch_hits += process_pair(&xa, &xb, edge_sets, windows, &mut sup);
                    }
                    (None, None) => break,
                    _ => {
                        return Err(CommonError::Parse(format!(
                            "gtasm pair support: {} 与 {} 记录数不齐（截断/错配）",
                            left.display(),
                            right.display()
                        )))
                    }
                }
            }
        }
    }
    Ok((sup, stats))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_graph(dir: &Path, id: u64, text: &str) {
        let mut w = trinity_chrysalis::pack::PackWriter::create(dir, GRAPHS_PACK).unwrap();
        w.append(id, KIND_GRAPH_OUT, text.as_bytes()).unwrap();
        w.finish().unwrap();
    }

    /// 两 comp 小图 + 手工 reads, 验证 k2v/窗扫/边过滤全链。
    #[test]
    fn support_pass_end_to_end_tiny() {
        let dir = std::env::temp_dir().join(format!("gt_sup_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let cbin = dir.join("Cbin0");
        std::fs::create_dir_all(&cbin).unwrap();
        // comp 0: 10→11（K=4 假装——K 由图推断, 测试用小 kmer）
        // kmer: 10:AAAA, 11:AAAC
        write_graph(
            &cbin,
            0,
            "Component 0\n10\t-1\t3\tAAAA\t1\n11\t10\t5\tAAAC\t1\n",
        );
        let listing = vec![(0u64, cbin.join("c0"))];
        let (k2v, edge_sets, k, _windows) = build_k2v(&listing).unwrap();
        assert_eq!(k, 4);
        // k2v 含正反义: AAAA 与 TTTT
        assert!(k2v.contains_key(&(kmer_to_intval(b"AAAA").unwrap(), 4)));
        assert!(k2v.contains_key(&(kmer_to_intval(b"TTTT").unwrap(), 4)));
        assert!(edge_sets[&0].contains(&trinity_butterfly::gtdecode::pack_uv(10, 11)));
        assert!(edge_sets[&0].contains(&trinity_butterfly::gtdecode::pack_uv(-1, 10)));

        // reads: "GAAAACT"（K=4: 窗 [0..3): GAAA, AAAA, AAAC——最后一窗
        // [3..7)=AACT 不查, quirk）→ 命中 AAAA(10)@i=1, AAAC(11)@i=2 → (10,11)+1
        let fq = dir.join("r.fq");
        std::fs::write(&fq, "@r1\nGAAAACT\n+\nIIIIIII\n").unwrap();
        let out_bin = dir.join("sup.bin");
        let st = run_support_pass(&listing, &[fq], SeqType::Fq, &out_bin).unwrap();
        assert_eq!(st.reads, 1);
        assert_eq!(st.pairs, 1);
        let sup = EdgeSupport::load(&out_bin).unwrap();
        assert_eq!(sup.get(0, 10, 11), 1);
        assert_eq!(sup.get(0, 11, 10), 0);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// 未命中重置 prev; 非边转移不计; 末窗不查（quirk 回归锁）。
    #[test]
    fn process_seq_semantics() {
        let mut k2v: K2V = FxHashMap::default();
        let aaaa = kmer_to_intval(b"AAAA").unwrap();
        let aaac = kmer_to_intval(b"AAAC").unwrap();
        k2v.insert((aaaa, 4), (0, 10));
        k2v.insert((aaac, 4), (0, 11));
        let mut edge_sets: EdgeSets = FxHashMap::default();
        edge_sets.insert(0u64, FxHashSet::default());
        // 边集为空（无任何边）→ 相邻命中也不计
        let mut sup = EdgeSupport::new();
        let mut pairs = 0u64;
        process_seq(b"GAAAACT", 4, &k2v, &edge_sets, &mut sup, &mut pairs);
        assert_eq!(pairs, 0);

        // 加上边 (10,11) → 计 1
        edge_sets
            .get_mut(&0)
            .unwrap()
            .insert(trinity_butterfly::gtdecode::pack_uv(10, 11));
        process_seq(b"GAAAACT", 4, &k2v, &edge_sets, &mut sup, &mut pairs);
        assert_eq!(pairs, 1);
        assert_eq!(sup.get(0, 10, 11), 1);

        // 含 N 的窗未命中 → prev 重置: "AAAANAAAC"——AAAA@0, 然后 N 窗未命中
        // ×2, AAAC@5 命中但 prev 已空 → 不计
        process_seq(b"AAAANAAAC", 4, &k2v, &edge_sets, &mut sup, &mut pairs);
        assert_eq!(pairs, 1);

        // 末窗 quirk: "AAAAC"（len 5, K 4）→ 只查 i=0: AAAA; i=1 的 AAAC 是
        // 最后一窗不查 → 无转移
        process_seq(b"AAAAC", 4, &k2v, &edge_sets, &mut sup, &mut pairs);
        assert_eq!(pairs, 1);
    }

    /// fastq line%4 语义: 多行序列收集 + 空行占行号（python 原样）。
    #[test]
    fn fastq_line_mod4_semantics() {
        let dir = std::env::temp_dir().join(format!("gt_sup_fq_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        // 两条记录, 第二条序列分两行（python 只收集 m==1 的行——第二行序列实际
        // 落在 m==2 位被丢弃: 这是 python 的（缺陷）语义, 逐样保留）
        let fq = dir.join("m.fq");
        std::fs::write(&fq, "@r1\nAAAA\n+\nIIII\n@r2\nAAAC\nGGGG\nIIII\n").unwrap();
        let mut k2v: K2V = FxHashMap::default();
        k2v.insert((kmer_to_intval(b"AAAA").unwrap(), 4), (0, 10));
        let edge_sets: EdgeSets = FxHashMap::default();
        let (sup, st) = scan_reads_file(&fq, SeqType::Fq, 4, &k2v, &edge_sets).unwrap();
        // r1: "AAAA"（len 4 → 无窗, len<=K）; r2: 只收 m==1 的 "AAAC"
        assert_eq!(st.reads, 2);
        assert_eq!(st.pairs, 0);
        drop(sup);
        let _ = std::fs::remove_dir_all(&dir);
    }

    // -----------------------------------------------------------------------
    // read pair 分支支持 pass
    // -----------------------------------------------------------------------

    /// BFS 窗口边界: 链 1→2→...→17→100（100 为分支点）。深度恰 15 边的节点
    /// 在窗口内, 16 边不在; v 自身不在; 深度上限截断生效。
    #[test]
    fn upstream_window_depth_boundary() {
        let mut pred: FxHashMap<i32, Vec<i32>> = FxHashMap::default();
        for i in 2..=17 {
            pred.insert(i, vec![i - 1]);
        }
        pred.insert(100, vec![17]);
        // 100←17(1 边)←16(2)...←3(15)←2(16)←1(17): 窗口 = {3..=17}
        let w = upstream_window(&pred, 100, PAIR_WINDOW_EDGES, PAIR_WINDOW_MAX_NODES);
        assert_eq!(w.len(), 15);
        assert!(w.contains(&3), "恰 15 边 → 在窗口");
        assert!(!w.contains(&2), "16 边 → 不在窗口");
        assert!(!w.contains(&100), "v 自身不在窗口");
        let w3 = upstream_window(&pred, 100, 3, PAIR_WINDOW_MAX_NODES);
        assert_eq!(w3, [17, 16, 15].into_iter().collect());
        // 环不挂死: 100↔17 互指
        let mut pred2: FxHashMap<i32, Vec<i32>> = FxHashMap::default();
        pred2.insert(100, vec![17]);
        pred2.insert(17, vec![100]);
        let wc = upstream_window(&pred2, 100, 15, 1024);
        assert_eq!(wc, [17].into_iter().collect());
    }

    /// pair pass 端到端小例（K=4 假装）: 分支点 12（出边 12→20 / 12→30）,
    /// 上游窗口 {10, 11}。mate 分落上游与候选边 = 计; 同侧（都在上游,
    /// 或另一条只命中分支点 v 自身）= 不计。
    #[test]
    fn pair_pass_counts_crossing_pair() {
        let dir = std::env::temp_dir().join(format!("gt_pair_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let cbin = dir.join("Cbin0");
        std::fs::create_dir_all(&cbin).unwrap();
        // comp 0: -1→10→11→12, 12→{20→21, 30→31}（kmer 按 (K-1) 重叠构造:
        // AACC→{ACCC, ACCG}, ACCC→CCCC, ACCG→CCGG）
        write_graph(
            &cbin,
            0,
            "Component 0\n10\t-1\t3\tAAAA\t1\n11\t10\t3\tAAAC\t1\n12\t11\t3\tAACC\t1\n\
             20\t12\t3\tACCC\t1\n21\t20\t3\tCCCC\t1\n30\t12\t3\tACCG\t1\n31\t30\t3\tCCGG\t1\n",
        );
        let listing = vec![(0u64, cbin.join("c0"))];
        // 4 对:
        // r1: L 上游链 (10,11,12) × R 含分支出边 (12,20) → 计 (0,12,20)
        // r2: 双侧都在上游 → 不计
        // r3: L 含 (12,20) × R 只命中分支点 12 自身（不在窗口） → 不计
        // r4: L 含 (12,30) × R 上游链 → 计 (0,12,30)
        let left = dir.join("l.fq");
        let right = dir.join("r.fq");
        std::fs::write(
            &left,
            "@r1/1\nAAAACC\n+\nIIIIII\n@r2/1\nAAACC\n+\nIIIII\n@r3/1\nAACCC\n+\nIIIII\n@r4/1\nAACCG\n+\nIIIII\n",
        )
        .unwrap();
        std::fs::write(
            &right,
            "@r1/2\nAACCC\n+\nIIIII\n@r2/2\nAAAACC\n+\nIIIIII\n@r3/2\nAACC\n+\nIIII\n@r4/2\nAAAACC\n+\nIIIIII\n",
        )
        .unwrap();
        let out_bin = dir.join("pair.bin");
        let st =
            run_pair_support_pass(&listing, &[left], &[right], SeqType::Fq, &out_bin).unwrap();
        assert_eq!(st.pairs, 4);
        assert_eq!(st.branch_hits, 2, "{st:?}");
        assert_eq!(st.name_mismatches, 0);
        assert_eq!(st.edges_with_support, 2);
        let sup = PairSupport::load(&out_bin).unwrap();
        assert_eq!(sup.get(0, 12, 20), 1);
        assert_eq!(sup.get(0, 12, 30), 1);
        assert_eq!(sup.get(0, 11, 12), 0); // 11 非分支点
        let _ = std::fs::remove_dir_all(&dir);
    }
}
