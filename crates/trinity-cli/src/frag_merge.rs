//! gtasm 碎片 pair-guided 合并（`--gtasm-frag-merge` 配套, 分支保真计划 (c)）。
//!
//! 位置: butterfly 组件池之后、harvest 之前（仅 `--path-engine gtasm` +
//! 双端 reads 启用）。残余引擎差距 = 中段碎片化: 完整路径在图内但 walk
//! 分支点选错臂, 正解被拆成 2+ 个 contig（参照坐标上重叠/紧邻）。本 pass
//! **只在同组件内、只在有 pair 证据时**做保守拼接, 不填序列（N 或重叠区
//! 直接连）——与已证伪的 SPAdes 式 gap closer（gap_close.rs 全图 tip 缝合）
//! 不同, 此处证据形状是 **contig 末端 × contig 末端的 mate 链接**, 且每
//! contig 至多参与一次合并。
//!
//! 机制（逐规格）:
//! 1. 收集本组件产出 contig 的末端序列（3'/5' 端各取 ≤ [`END_LEN`]bp）,
//!    建末端 k2v: kmer（K = 图 kmer 长, 自 graphs.pack 推断）→
//!    [(组件序, contig 序, 端)] token 列表（一 kmer 可属多个末端——碎片在
//!    参照坐标上重叠时共享 kmer, 不做先写先赢）, 正反义双写（revcomp 同键）。
//! 2. pair 链接: 原始 reads（--left/--right, lockstep 同序, 沿用
//!    gt_support 的严格 4 行 fq / FASTA 配对读取）全窗线程化; mate A 命中
//!    contig X 的 3' 末端且 mate B 命中同组件 contig Y 的 5' 末端 →
//!    link(X,Y) += 1（双取向各评一次, 同对同链接去重; 自链接除外）。
//! 3. 合并: 组件内按 link 权重降序（同重按 (X,Y) 字典序, 确定性）, 权重
//!    ≥ [`LINK_MIN`] 且两 contig 均未参与过合并 → 拼接: 先查末端重叠
//!    （X 后缀 × Y 前缀, ≥ [`OLC_MIN_OV`]bp 且一致度 ≥90% 则 OLC 合,
//!    取最长）; 无重叠则以 [`N_GAP`] 个 N 连接。产物 < [`MERGED_MIN_LEN`]
//!    拒收（两原件保留, 不消耗名额）。
//! 4. 输出契约与 harvest 完全一致（min_len 过滤 + 全序列去重 + 60 列折行;
//!    缺失 allProbPaths → Err）, 合并产物 header = X 原 header +
//!    ` fragmerge:+<Y acc>`（首 token 不变, gene_trans_map 解析不受响）。
//!
//! 安全阀: 单 mate 末端命中集 > [`MAX_END_HITS`] 的对放弃计数（低复杂度/
//! 重复序列 read 的命中集爆炸 → 二次方笛卡尔积, 且此类证据本无特异性）。

use std::fs;
use std::path::{Path, PathBuf};

use rustc_hash::{FxHashMap, FxHashSet};

use trinity_butterfly::parse_component_graph;
use trinity_chrysalis::pack::{read_at, PackIndex, GRAPHS_PACK, KIND_ALL_PROB_PATHS, KIND_GRAPH_OUT};
use trinity_common::error::CommonError;
use trinity_common::fasta::FastaReader;
use trinity_common::io_util::open_maybe_gz;
use trinity_common::kmer::{kmer_to_intval, revcomp_val, KmerId};

use crate::args::SeqType;
use crate::gt_support::{norm_read_name, read_fq_record};

/// contig 末端截取长度（bp）: 3'/5' 端各取 ≤150bp 建末端 k2v。
pub const END_LEN: usize = 150;
/// link 计数阈值: ≥ 该值才尝试合并（规格写死 2）。
pub const LINK_MIN: u32 = 2;
/// 无末端重叠时的 N 连接长度。
pub const N_GAP: usize = 50;
/// OLC 最小重叠（bp）。
pub const OLC_MIN_OV: usize = 40;
/// OLC 重叠扫描上限（bp）: 防长 contig 两两全扫; 末端区 150bp 之外的巨重叠
/// 罕见, 截断为下近似。
const OLC_MAX_OV: usize = 1000;
/// 合并产物最小长度（< 则拒收, 两原件保留且不消耗合并名额）。
pub const MERGED_MIN_LEN: usize = 200;
/// 单 mate 末端命中集安全阀（见模块文档）。
const MAX_END_HITS: usize = 64;

/// 端标记: 5' 端。
const END5: u64 = 0;
/// 端标记: 3' 端。
const END3: u64 = 1;

/// 末端 k2v 的 token 打包: comp_ord(高 31 位) | local_idx(32 位) | end(1 位)。
/// （组件数与单组件 contig 数均远小于位宽上限。）
fn pack_token(comp_ord: usize, local_idx: usize, end: u64) -> u64 {
    ((comp_ord as u64) << 33) | ((local_idx as u64) << 1) | end
}

/// token 解包 → (comp_ord, local_idx, end)。
fn unpack_token(t: u64) -> (u32, u32, u64) {
    ((t >> 33) as u32, ((t >> 1) & 0xffff_ffff) as u32, t & 1)
}

/// 末端 k2v: (kmer 2bit 值, kmer 长度) → token 列表（去重后有序）。
type EndK2V = FxHashMap<(KmerId, u32), Vec<u64>>;

/// 链接表: (comp_ord, x_idx, y_idx) → pair 计数（x 的 3' → y 的 5'）。
type LinkMap = FxHashMap<(u32, u32, u32), u32>;

/// 碎片合并统计（日志 + .stage_timings 统计行用）。
#[derive(Debug, Default, Clone, Copy)]
pub struct FragMergeStats {
    /// 输入 read 对数。
    pub pairs: u64,
    /// lockstep 名不匹配的对数（>0 = 两文件可能错序, 警告用）。
    pub name_mismatches: u64,
    /// 有 ≥1 对支持的候选链接数（去重后 (X,Y) 对数）。
    pub links_total: usize,
    /// 达阈值（≥ [`LINK_MIN`]）的链接数。
    pub links_kept: usize,
    /// OLC 重叠合并数。
    pub merged_olc: usize,
    /// N 连接合并数。
    pub merged_n: usize,
    /// 产物 < [`MERGED_MIN_LEN`] 被拒的合并尝试数。
    pub rejected_short: usize,
    /// 输入 contig 总数。
    pub contigs_in: usize,
    /// 合并后 contig 总数（min_len/去重前）。
    pub contigs_out: usize,
}

/// 合并版 harvest 的返回: 收录转录本数 + 统计。
pub struct FragMergeOutcome {
    pub n_transcripts: usize,
    pub stats: FragMergeStats,
}

/// 合并版 harvest: 读全部组件 contigs → 末端 k2v → 原始 reads pair 链接 →
/// 组件内保守拼接 → 按 harvest 契约落 `final_path`。
pub fn harvest_with_frag_merge(
    listing: &[(u64, PathBuf)],
    min_len: usize,
    final_path: &Path,
    left_files: &[PathBuf],
    right_files: &[PathBuf],
    seq_type: SeqType,
) -> Result<FragMergeOutcome, CommonError> {
    if left_files.len() != right_files.len() {
        return Err(CommonError::Parse(format!(
            "gtasm frag merge: --left（{} 个）与 --right（{} 个）文件数不齐",
            left_files.len(),
            right_files.len()
        )));
    }
    // ---- 1. 读全部组件的 allProbPaths 记录（缺失契约同 harvest: 收集后 Err;
    // 但在昂贵的 reads 扫描之前即报, 不白扫）----
    let comps = read_all_components(listing)?;
    let k = infer_graph_k(listing)?;
    let k2v = build_end_k2v(&comps, k);
    eprintln!(
        "-gtasm frag merge: {} comps {} contigs, 末端 k2v {} 键 (K={k})",
        comps.len(),
        comps.iter().map(Vec::len).sum::<usize>(),
        k2v.len()
    );

    // ---- 2. 扫原始 reads（逐文件对并行, 末尾合并链接表）----
    let (links, scan_stats) = scan_pair_files(left_files, right_files, seq_type, k, &k2v)?;

    // ---- 3. 组件内合并 + harvest 语义输出 ----
    let (n_transcripts, mut stats) = emit_merged(&comps, &links, min_len, final_path)?;
    stats.pairs = scan_stats.pairs;
    stats.name_mismatches = scan_stats.name_mismatches;
    stats.links_total = links.len();
    Ok(FragMergeOutcome { n_transcripts, stats })
}

/// 逐组件读 allProbPaths 记录入内存（缺失 → 同 harvest 契约: 处理完其余后
/// Err, 产物不落盘）。
fn read_all_components(listing: &[(u64, PathBuf)]) -> Result<Vec<Vec<(String, String)>>, CommonError> {
    let mut comps = Vec::with_capacity(listing.len());
    let mut missing: Vec<String> = Vec::new();
    // listing 按 id 升序 → 同 Cbin 连续; 当前 bin 的索引与句柄各开一次。
    let mut cur_bin: Option<(PathBuf, PackIndex, fs::File)> = None;
    for (id, base) in listing {
        let dir = base.parent().map(|d| d.to_path_buf()).unwrap_or_default();
        if cur_bin.as_ref().map(|(d, _, _)| d != &dir).unwrap_or(true) {
            cur_bin = crate::harvest::load_bin(&dir).ok();
        }
        let logical = format!("{}.graph.allProbPaths.fasta", base.display());
        let text = cur_bin
            .as_ref()
            .and_then(|(_, idx, file)| idx.find(*id, KIND_ALL_PROB_PATHS).map(|(o, l)| (file, o, l)))
            .and_then(|(file, o, l)| read_at(file, o, l).ok())
            .and_then(|bytes| String::from_utf8(bytes).ok());
        match text {
            Some(t) => comps.push(crate::harvest::records(&t)),
            None => {
                eprintln!("Error, no fasta file reported as: {logical}");
                missing.push(logical);
                comps.push(Vec::new());
            }
        }
    }
    if !missing.is_empty() {
        return Err(CommonError::Parse(format!(
            "Error, {} butterfly allProbPaths fasta file(s) missing (first: {})",
            missing.len(),
            missing[0]
        )));
    }
    Ok(comps)
}

/// 自 graphs.pack 推断图 kmer 长 K（首组件首 kmer 长度; 无 kmer 时回退 24,
/// 与 gt_support::build_k2v 的缺省一致）。
fn infer_graph_k(listing: &[(u64, PathBuf)]) -> Result<usize, CommonError> {
    let Some((id, base)) = listing.first() else {
        return Ok(24);
    };
    let dir = base.parent().map(|d| d.to_path_buf()).unwrap_or_default();
    let idx = PackIndex::load(&dir.join(format!("{GRAPHS_PACK}.idx")))?;
    let file = fs::File::open(dir.join(format!("{GRAPHS_PACK}.pack")))?;
    let (offset, len) = idx.find(*id, KIND_GRAPH_OUT).ok_or_else(|| {
        CommonError::Parse(format!(
            "graphs pack {} 缺记录 (id={id}, kind=GRAPH_OUT)",
            dir.display()
        ))
    })?;
    let bytes = read_at(&file, offset, len)?;
    let text = String::from_utf8(bytes)
        .map_err(|e| CommonError::Parse(format!("comp {id} graph.out 非 UTF-8: {e}")))?;
    let g = parse_component_graph(&text)
        .map_err(|e| CommonError::Parse(format!("comp {id} graph.out 解析失败: {e}")))?;
    Ok(g.names.values().next().map(|s| s.len()).unwrap_or(24))
}

/// 建末端 k2v: 每 contig 的 5'/3' 端各取 ≤ [`END_LEN`]bp, 全窗 kmer（正反义
/// 双写）→ token; 值列表构建后排序去重（同 kmer 重复命中同末端只留一份）。
fn build_end_k2v(comps: &[Vec<(String, String)>], k: usize) -> EndK2V {
    let mut k2v: EndK2V = FxHashMap::default();
    let k32 = k as u32;
    for (comp_ord, contigs) in comps.iter().enumerate() {
        for (local_idx, (_, seq)) in contigs.iter().enumerate() {
            let b = seq.as_bytes();
            if b.len() < k {
                continue; // 短于 K 的 contig 无窗, 不可链接
            }
            let end_len = b.len().min(END_LEN);
            for (region, end) in [(&b[..end_len], END5), (&b[b.len() - end_len..], END3)] {
                let token = pack_token(comp_ord, local_idx, end);
                for i in 0..=(region.len() - k) {
                    if let Ok(val) = kmer_to_intval(&region[i..i + k]) {
                        k2v.entry((val, k32)).or_default().push(token);
                        let rc = revcomp_val(val, k);
                        k2v.entry((rc, k32)).or_default().push(token);
                    }
                }
            }
        }
    }
    for v in k2v.values_mut() {
        v.sort_unstable();
        v.dedup();
    }
    k2v
}

/// 单 mate 线程化: 全窗 k2v 查询, 命中的 token 集（位置信息对本 pass 无用,
/// 任意命中即视为"命中该末端"）。
fn thread_mate(seq: &[u8], k: usize, k2v: &EndK2V) -> FxHashSet<u64> {
    let mut out = FxHashSet::default();
    if seq.len() < k {
        return out;
    }
    let k32 = k as u32;
    for i in 0..=(seq.len() - k) {
        if let Some(tokens) = kmer_to_intval(&seq[i..i + k])
            .ok()
            .and_then(|val| k2v.get(&(val, k32)))
        {
            out.extend(tokens.iter().copied());
        }
    }
    out
}

/// 单对 mates 的链接计数: SA × SB 双取向笛卡尔积——3' 端 × 5' 端、同组件、
/// 非自链接 → 候选 (comp, x, y); 每对每链接至多计 1（调用方去重后入表）。
/// 命中集超 [`MAX_END_HITS`] 的对整个放弃（安全阀, 见模块文档）。
fn process_pair_tokens(
    sa: &FxHashSet<u64>,
    sb: &FxHashSet<u64>,
    counted: &mut FxHashSet<(u32, u32, u32)>,
) {
    if sa.len() > MAX_END_HITS || sb.len() > MAX_END_HITS {
        return;
    }
    for &ta in sa {
        let (ca, ia, ea) = unpack_token(ta);
        for &tb in sb {
            let (cb, ib, eb) = unpack_token(tb);
            if ca != cb || ia == ib {
                continue; // 跨组件不链（合并只同组件内）; 自链接除外
            }
            match (ea, eb) {
                // A 命中 X 的 3' 端 × B 命中 Y 的 5' 端 → link(X, Y)
                (END3, END5) => {
                    counted.insert((ca, ia, ib));
                }
                // 对称取向: A 命中 Y 的 5' 端 × B 命中 X 的 3' 端 → link(X, Y)
                (END5, END3) => {
                    counted.insert((ca, ib, ia));
                }
                _ => {} // 同端组合不构链
            }
        }
    }
}

/// 扫全部 mates 文件对（逐文件对并行, 末尾合并）。返回 (链接表, 扫描统计)。
fn scan_pair_files(
    left_files: &[PathBuf],
    right_files: &[PathBuf],
    seq_type: SeqType,
    k: usize,
    k2v: &EndK2V,
) -> Result<(LinkMap, FragMergeStats), CommonError> {
    let mut results: Vec<Result<(LinkMap, FragMergeStats), CommonError>> = Vec::new();
    std::thread::scope(|scope| {
        let handles: Vec<_> = left_files
            .iter()
            .zip(right_files.iter())
            .map(|(l, r)| scope.spawn(|| scan_one_file_pair(l, r, seq_type, k, k2v)))
            .collect();
        for h in handles {
            results.push(h.join().unwrap_or_else(|_| {
                Err(CommonError::Parse("gtasm frag merge 扫描线程 panic".to_string()))
            }));
        }
    });
    let mut links: LinkMap = FxHashMap::default();
    let mut stats = FragMergeStats::default();
    for r in results {
        let (m, st) = r?;
        stats.pairs += st.pairs;
        stats.name_mismatches += st.name_mismatches;
        for (&key, &c) in &m {
            *links.entry(key).or_default() += c;
        }
    }
    Ok((links, stats))
}

/// 扫一对 mates 文件（fq 严格 4 行 lockstep / fa FASTA lockstep, 同序配对;
/// 名不匹配计数, 末尾统计告警）。返回 (本文件对的链接表, 统计)。
fn scan_one_file_pair(
    left: &Path,
    right: &Path,
    seq_type: SeqType,
    k: usize,
    k2v: &EndK2V,
) -> Result<(LinkMap, FragMergeStats), CommonError> {
    let mut links: LinkMap = FxHashMap::default();
    let mut stats = FragMergeStats::default();
    let mut counted: FxHashSet<(u32, u32, u32)> = FxHashSet::default();
    let mut ra = open_maybe_gz(left)?;
    let mut rb = open_maybe_gz(right)?;
    let mut tally = |sa: &FxHashSet<u64>, sb: &FxHashSet<u64>| {
        counted.clear();
        process_pair_tokens(sa, sb, &mut counted);
        for &key in &counted {
            *links.entry(key).or_default() += 1;
        }
    };
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
                    let sa = thread_mate(a[1].trim().as_bytes(), k, k2v);
                    let sb = thread_mate(b[1].trim().as_bytes(), k, k2v);
                    tally(&sa, &sb);
                }
                (None, None) => break,
                _ => {
                    return Err(CommonError::Parse(format!(
                        "gtasm frag merge: {} 与 {} 记录数不齐（截断/错配）",
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
                        let sa = thread_mate(a.sequence.as_bytes(), k, k2v);
                        let sb = thread_mate(b.sequence.as_bytes(), k, k2v);
                        tally(&sa, &sb);
                    }
                    (None, None) => break,
                    _ => {
                        return Err(CommonError::Parse(format!(
                            "gtasm frag merge: {} 与 {} 记录数不齐（截断/错配）",
                            left.display(),
                            right.display()
                        )))
                    }
                }
            }
        }
    }
    Ok((links, stats))
}

/// X 后缀 × Y 前缀的最长合格重叠: ≥ [`OLC_MIN_OV`]bp 且一致度 ≥90%
/// （mm*10 ≤ ov, 逐位早退）; 扫描上限 [`OLC_MAX_OV`]bp。返回重叠长。
fn find_overlap(x: &[u8], y: &[u8]) -> Option<usize> {
    let max_ov = x.len().min(y.len()).min(OLC_MAX_OV);
    if max_ov < OLC_MIN_OV {
        return None;
    }
    let mut ov = max_ov;
    while ov >= OLC_MIN_OV {
        let xs = &x[x.len() - ov..];
        let ys = &y[..ov];
        let mut mm = 0usize;
        let mut ok = true;
        for (&a, &b) in xs.iter().zip(ys.iter()) {
            if a != b {
                mm += 1;
                if mm * 10 > ov {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            return Some(ov);
        }
        ov -= 1;
    }
    None
}

/// 拼接两 contig: 有合格末端重叠 → OLC 合（Y 去重叠前缀）; 否则 N_GAP 个
/// N 连接。返回 (产物, 是否 OLC)。
fn join_contigs(x: &str, y: &str) -> (String, bool) {
    match find_overlap(x.as_bytes(), y.as_bytes()) {
        Some(ov) => {
            let mut s = String::with_capacity(x.len() + y.len() - ov);
            s.push_str(x);
            s.push_str(&y[ov..]);
            (s, true)
        }
        None => {
            let mut s = String::with_capacity(x.len() + N_GAP + y.len());
            s.push_str(x);
            s.push_str(&"N".repeat(N_GAP));
            s.push_str(y);
            (s, false)
        }
    }
}

/// 组件内合并 + harvest 语义输出（min_len 过滤 / 全序列去重 / 60 列折行）。
/// 返回 (收录转录本数, 统计)。
fn emit_merged(
    comps: &[Vec<(String, String)>],
    links: &LinkMap,
    min_len: usize,
    final_path: &Path,
) -> Result<(usize, FragMergeStats), CommonError> {
    let mut stats = FragMergeStats::default();
    // 组件序 → 达阈值链接（降序, 同重按 (x,y) 字典序, 确定性）
    let mut by_comp: FxHashMap<u32, Vec<(u32, u32, u32)>> = FxHashMap::default();
    for (&(c, x, y), &cnt) in links {
        if cnt >= LINK_MIN {
            by_comp.entry(c).or_default().push((x, y, cnt));
            stats.links_kept += 1;
        }
    }
    for v in by_comp.values_mut() {
        v.sort_by(|a, b| b.2.cmp(&a.2).then(a.0.cmp(&b.0)).then(a.1.cmp(&b.1)));
    }

    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut out = String::new();
    let mut n = 0usize;
    for (comp_ord, contigs) in comps.iter().enumerate() {
        stats.contigs_in += contigs.len();
        let empty = Vec::new();
        let comp_links = by_comp.get(&(comp_ord as u32)).unwrap_or(&empty);
        // 每 contig 至多参与一次合并（先到先得）; 产物挂在 X 槽位。
        let mut used = vec![false; contigs.len()];
        let mut products: Vec<Option<(String, String)>> = vec![None; contigs.len()];
        for &(x, y, _cnt) in comp_links {
            let (xi, yi) = (x as usize, y as usize);
            if used[xi] || used[yi] {
                continue;
            }
            let (merged, is_olc) = join_contigs(&contigs[xi].1, &contigs[yi].1);
            if merged.len() < MERGED_MIN_LEN {
                stats.rejected_short += 1;
                continue; // 不消耗名额, 两原件仍可参与后续链接
            }
            let y_acc = contigs[yi].0.split_whitespace().next().unwrap_or("");
            products[xi] = Some((format!("{} fragmerge:+{y_acc}", contigs[xi].0), merged));
            used[xi] = true;
            used[yi] = true;
            if is_olc {
                stats.merged_olc += 1;
            } else {
                stats.merged_n += 1;
            }
        }
        // 输出: X 槽位出产物（Y 被吞）, 未参与者原样出——同 harvest 的
        // min_len 过滤 + 全序列去重契约。
        for (i, (header, seq)) in contigs.iter().enumerate() {
            let (h, s) = match &products[i] {
                Some((ph, ps)) => (ph.as_str(), ps.as_str()),
                None if !used[i] => (header.as_str(), seq.as_str()),
                None => continue, // 作为 Y 被吞
            };
            if s.len() < min_len {
                continue;
            }
            if !seen.insert(s.to_string()) {
                eprintln!("-duplicate sequence detected, excluding it. (component {comp_ord})");
                continue;
            }
            out.push('>');
            out.push_str(h);
            out.push('\n');
            for chunk in s.as_bytes().chunks(60) {
                out.push_str(std::str::from_utf8(chunk).unwrap_or(""));
                out.push('\n');
            }
            n += 1;
        }
        stats.contigs_out += products.iter().filter(|p| p.is_some()).count()
            + used.iter().filter(|&&u| !u).count();
    }
    // harvest 无 checkpoint, 每次重写——本 pass 同（产物确定性）。
    fs::write(final_path, &out)?;
    Ok((n, stats))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 造 token 集查询: (comp, idx, end) → token。
    fn tok(comp: usize, idx: usize, end: u64) -> u64 {
        pack_token(comp, idx, end)
    }

    // -----------------------------------------------------------------------
    // pair 链接计数
    // -----------------------------------------------------------------------

    /// 末端 k2v 构建 + 线程化 + 链接计数全链（K=4 假装）:
    /// comp 0 两条 contig——X 的 3' 区独有 kmer "GGGG", Y 的 5' 区独有
    /// kmer "TTTT"; mate A 只含 GGGG 系窗, mate B 只含 TTTT 系窗 → link(X,Y)。
    #[test]
    fn pair_link_counting() {
        // X: 5' 区 AAAA..., 3' 区 ...GGGG; Y: 5' 区 TTTT..., 3' 区 ...CCCC
        // （K=4; 序列 >2*END_LEN 不现实, 小例全序列即两端——端内 kmer 用
        // 组成区分: X=AAAAGGGG, Y=TTTTCCCC, Z=ACACACAC）
        let comps = vec![vec![
            ("X".to_string(), "AAAAGGGG".to_string()),
            ("Y".to_string(), "TTTTCCCC".to_string()),
            ("Z".to_string(), "ACACACAC".to_string()),
        ]];
        let k2v = build_end_k2v(&comps, 4);
        // 两端同序列（长度<150, 全序列同为 5' 与 3' 端）: GGGG → X 的两端 token
        let g = kmer_to_intval(b"GGGG").unwrap();
        let toks = &k2v[&(g, 4)];
        assert!(toks.contains(&tok(0, 0, END5)));
        assert!(toks.contains(&tok(0, 0, END3)));
        // revcomp 双写: CCCC 同键（GGGG 的 rc = CCCC）→ 同时含 Y 的两端 token
        assert!(toks.contains(&tok(0, 1, END5)));
        assert!(toks.contains(&tok(0, 1, END3)));

        // mate A = "AAGGGG"（窗: AAGG, AGGG, GGGG → 命中的只有 GGGG → X 两端）
        // mate B = "AATTTT"（窗 TTTT → Y 两端; rc AAAA → X 两端也命中!）
        // → 用不含 rc 混淆的 read: B = "TTTTA"（窗: TTTT → Y; rc(TTTT)=AAAA → X）
        // 为避开 rc 互踩, 改用回文外设计: X 3' 独有 "GGGG"（rc CCCC 在 Y 3'）,
        // 链接方向由 5'/3' token 区分, rc 命中只会加 token, 方向断言按集合评。
        let sa = thread_mate(b"CCGGGG", 4, &k2v); // 窗 CCGG, CGGG, GGGG
        let sb = thread_mate(b"TTTTGG", 4, &k2v); // 窗 TTTT, TTGG
        let mut counted = FxHashSet::default();
        process_pair_tokens(&sa, &sb, &mut counted);
        // sa: GGGG → X5,X3,Y5(rc CCCC→GGGG 键),Y3; CCGG rc=CCGG, CGGG rc=CCCG 不在表
        // sb: TTTT → Y5,Y3, rc AAAA → X5,X3
        // 3'×5' 组合: X3×Y5 → (0,0,1); Y3×X5 → (0,1,0); X3×X5 自链排除; Y3×Y5 自链排除
        assert!(counted.contains(&(0, 0, 1)));
        assert!(counted.contains(&(0, 1, 0)));
        assert!(!counted.iter().any(|&(_, x, y)| x == y), "自链接不得入表");
        assert_eq!(counted.len(), 2);

        // Z 的 kmer 与 X/Y 无交集 → 无链接
        let sz = thread_mate(b"ACACAC", 4, &k2v);
        let mut c2 = FxHashSet::default();
        process_pair_tokens(&sz, &sz, &mut c2);
        assert!(c2.is_empty(), "Z 只命中自身, 自链接排除 → 空");
    }

    /// 命中集安全阀: 超 MAX_END_HITS 的对整个放弃。
    #[test]
    fn pair_link_safety_cap() {
        let big: FxHashSet<u64> = (0..(MAX_END_HITS as u64 + 1))
            .map(|i| tok(0, i as usize, END3))
            .collect();
        let one: FxHashSet<u64> = [tok(0, 999, END5)].into_iter().collect();
        let mut counted = FxHashSet::default();
        process_pair_tokens(&big, &one, &mut counted);
        assert!(counted.is_empty());
    }

    // -----------------------------------------------------------------------
    // OLC 重叠合并 / N 连接
    // -----------------------------------------------------------------------

    /// OLC: 60bp 重叠 1 处错配（59/60 ≥ 90%）→ 无 N 合并; 取最长合格重叠。
    #[test]
    fn olc_overlap_merge() {
        // X = 前缀 100A + 60bp 重叠区; Y = 重叠区(1 错配) + 后缀 80C
        let ov_x: String = "ACGT".repeat(15); // 60bp
        let mut ov_y = ov_x.clone().into_bytes();
        ov_y[30] = b'T'; // 1 处错配（原位是 G）
        let ov_y = String::from_utf8(ov_y).unwrap();
        let x = format!("{}{}", "A".repeat(100), ov_x);
        let y = format!("{}{}", ov_y, "C".repeat(80));
        let (merged, is_olc) = join_contigs(&x, &y);
        assert!(is_olc);
        assert_eq!(merged.len(), 100 + 60 + 80);
        assert!(merged.starts_with(&"A".repeat(100)));
        assert!(merged.ends_with(&"C".repeat(80)));
        assert!(!merged.contains('N'));

        // 错配超 10% → 不 OLC: 60bp 中 10 处错配（50/60 < 90%）→ 短重叠也不合格
        let mut ov_bad = ov_x.clone().into_bytes();
        for i in 0..10 {
            ov_bad[i * 6] = b'T';
        }
        let y_bad = format!("{}{}", String::from_utf8(ov_bad).unwrap(), "C".repeat(80));
        let (m2, olc2) = join_contigs(&x, &y_bad);
        assert!(!olc2, "10/60 错配且更短重叠也错配密集 → 落 N 连接");
        assert_eq!(m2.len(), x.len() + N_GAP + y_bad.len());
    }

    /// 无重叠 → N 连接（恰 N_GAP 个 N, 两侧原序列）。
    #[test]
    fn n_join_when_no_overlap() {
        let x = "A".repeat(300);
        let y = "C".repeat(250);
        let (merged, is_olc) = join_contigs(&x, &y);
        assert!(!is_olc);
        assert_eq!(merged.len(), 300 + N_GAP + 250);
        assert_eq!(&merged[300..300 + N_GAP], "N".repeat(N_GAP));
        assert!(merged.starts_with(&x));
        assert!(merged.ends_with(&y));
    }

    /// 重叠 < OLC_MIN_OV → 不 OLC（39bp 完全匹配也不够; 共享段用非周期
    /// 序列, 防 G 串错位滑出 38/40 的假合格重叠）。
    #[test]
    fn olc_min_overlap_enforced() {
        let ov = "AGCTTAGGCATGCAATCGGTACGATTCGAACTTGGCCAT"; // 39bp 非周期
        assert_eq!(ov.len(), 39);
        let x = format!("{}{}", "A".repeat(200), ov);
        let y = format!("{}{}", ov, "C".repeat(200));
        let (_, is_olc) = join_contigs(&x, &y);
        assert!(!is_olc);
    }

    // -----------------------------------------------------------------------
    // 每 contig 一次约束 / 阈值 / 短产物拒收
    // -----------------------------------------------------------------------

    fn comp3() -> Vec<Vec<(String, String)>> {
        // 三条互不重叠的 300bp contig（纯碱基区分, 保证无 ≥40bp 重叠 → N 连接）
        vec![vec![
            ("c0_g0_i0 len=300".to_string(), "A".repeat(300)),
            ("c0_g0_i1 len=300".to_string(), "C".repeat(300)),
            ("c0_g0_i2 len=300".to_string(), "G".repeat(300)),
        ]]
    }

    #[test]
    fn merge_each_contig_once() {
        let comps = comp3();
        // links: (0,1)=5, (1,2)=4, (0,2)=3 —— 全达阈值; 权重降序先 (0,1),
        // 此后 1 已用 → (1,2) 跳过; 0 已用 → (0,2) 跳过。只合并一次。
        let mut links: LinkMap = FxHashMap::default();
        links.insert((0, 0, 1), 5);
        links.insert((0, 1, 2), 4);
        links.insert((0, 0, 2), 3);
        let dir = std::env::temp_dir().join(format!("frag_out_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let out_fa = dir.join("t.fa");
        let (n, stats) = emit_merged(&comps, &links, 100, &out_fa).unwrap();
        assert_eq!(stats.links_kept, 3);
        assert_eq!(stats.merged_n, 1);
        assert_eq!(stats.merged_olc, 0);
        assert_eq!(stats.contigs_in, 3);
        assert_eq!(stats.contigs_out, 2);
        assert_eq!(n, 2); // 产物 650 + 幸存者 300
        let text = std::fs::read_to_string(&out_fa).unwrap();
        assert!(text.contains("fragmerge:+c0_g0_i1"), "产物 header 标来源");
        assert!(text.contains(&"N".repeat(N_GAP)));
        assert!(text.contains(">c0_g0_i2"), "未参与者原样出");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn merge_threshold_and_short_reject() {
        // 阈值: (0,1)=1 不达 LINK_MIN → 不合并
        let comps = comp3();
        let mut links: LinkMap = FxHashMap::default();
        links.insert((0, 0, 1), 1);
        let dir = std::env::temp_dir().join(format!("frag_thr_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let out_fa = dir.join("t.fa");
        let (n, stats) = emit_merged(&comps, &links, 100, &out_fa).unwrap();
        assert_eq!(stats.links_kept, 0);
        assert_eq!(stats.merged_n, 0);
        assert_eq!(n, 3);

        // 短产物拒收: 两条 60bp + 50N = 170 < 200 → 不合并, 原件保留,
        // 名额不消耗（后续 (0,2) 链接仍可成）。min_len=50 让短原件也能出。
        let comps2 = vec![vec![
            ("s0 len=60".to_string(), "A".repeat(60)),
            ("s1 len=60".to_string(), "C".repeat(60)),
            ("s2 len=300".to_string(), "G".repeat(300)),
        ]];
        let mut links2: LinkMap = FxHashMap::default();
        links2.insert((0, 0, 1), 9); // 产物 170 < 200 → 拒
        links2.insert((0, 0, 2), 2); // 0 未消耗 → 可成（产物 410）
        let out_fa2 = dir.join("t2.fa");
        let (_n2, stats2) = emit_merged(&comps2, &links2, 50, &out_fa2).unwrap();
        assert_eq!(stats2.rejected_short, 1);
        assert_eq!(stats2.merged_n, 1, "0 未被短拒消耗 → (0,2) 成");
        let text2 = std::fs::read_to_string(&out_fa2).unwrap();
        assert!(text2.contains(">s1"), "短拒原件保留");
        let _ = std::fs::remove_dir_all(&dir);
    }

    // -----------------------------------------------------------------------
    // 端到端小例（pack + fq 全链）
    // -----------------------------------------------------------------------

    /// 写 graphs.pack（只为 K 推断）与 butterfly.pack（allProbPaths 记录）。
    fn write_packs(dir: &Path, graph_text: &str, recs: &[(u64, &str)]) {
        let mut w = trinity_chrysalis::pack::PackWriter::create(dir, GRAPHS_PACK).unwrap();
        w.append(recs[0].0, KIND_GRAPH_OUT, graph_text.as_bytes()).unwrap();
        w.finish().unwrap();
        let mut w =
            trinity_chrysalis::pack::PackWriter::create(dir, trinity_chrysalis::pack::BUTTERFLY_PACK)
                .unwrap();
        for (id, text) in recs {
            w.append(*id, KIND_ALL_PROB_PATHS, text.as_bytes()).unwrap();
        }
        w.finish().unwrap();
    }

    /// comp 0 两条 200bp contig（无 ≥40bp 重叠）, 2 对 reads: L 落 X 3' 区,
    /// R 落 Y 5' 区 → link(X,Y)=2 达阈值 → N 连接产物 450bp 替换原件。
    #[test]
    fn end_to_end_tiny() {
        let dir = std::env::temp_dir().join(format!("frag_e2e_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let cbin = dir.join("Cbin0");
        std::fs::create_dir_all(&cbin).unwrap();
        // 真实 contig 序列: X = 200bp（3' 区带特异 8-mer GCGCGCGC）,
        // Y = 200bp（5' 区带特异 8-mer TATATATA）, 主体用互不重叠的低同源填充。
        // （周期填充序列使 rc 交叉命中不可避免——两方向链接都会达阈值,
        // 但合并先到先得只成一次, 断言据此放宽 links_kept ≥ 1。）
        let x = format!("{}{}", "AC".repeat(80), "GCGCGCGC".repeat(2) + &"AC".repeat(12));
        let y = format!("{}{}", "TATATATA".repeat(2) + &"GT".repeat(24), "GT".repeat(68));
        assert_eq!(x.len(), 200);
        assert_eq!(y.len(), 200);
        let allprob = format!(">c0_g0_i0 len=200\n{x}\n>c0_g0_i1 len=200\n{y}\n");
        // K=8 由图推断（kmer 内容只需长度合法）
        write_packs(
            &cbin,
            "Component 0\n10\t-1\t3\tAAAAAAAA\t1\n11\t10\t5\tAAAAAAAC\t1\n",
            &[(0, &allprob)],
        );
        let listing = vec![(0u64, cbin.join("c0"))];
        // reads: L = X[160..192]（X 3' 端 50..200 内, 不含 5' 端 0..150 → 只 X3）
        // R = Y[0..32]（Y 5' 端 0..150 内, 不含 3' 端 50..200 → 只 Y5）
        let l_seq = &x[160..192];
        let r_seq = &y[0..32];
        let left = dir.join("l.fq");
        let right = dir.join("r.fq");
        std::fs::write(
            &left,
            format!("@p1/1\n{l_seq}\n+\n{}\n@p2/1\n{l_seq}\n+\n{}\n", "I".repeat(32), "I".repeat(32)),
        )
        .unwrap();
        std::fs::write(
            &right,
            format!("@p1/2\n{r_seq}\n+\n{}\n@p2/2\n{r_seq}\n+\n{}\n", "I".repeat(32), "I".repeat(32)),
        )
        .unwrap();
        let out_fa = dir.join("out.fa");
        let oc = harvest_with_frag_merge(
            &listing,
            100,
            &out_fa,
            &[left],
            &[right],
            SeqType::Fq,
        )
        .unwrap();
        assert_eq!(oc.stats.pairs, 2);
        assert_eq!(oc.stats.name_mismatches, 0);
        assert!(oc.stats.links_kept >= 1, "{:?}", oc.stats);
        assert_eq!(oc.stats.merged_n, 1, "{:?}", oc.stats);
        assert_eq!(oc.stats.merged_olc, 0);
        assert_eq!(oc.n_transcripts, 1);
        let text = std::fs::read_to_string(&out_fa).unwrap();
        assert!(text.contains(">c0_g0_i0 len=200 fragmerge:+c0_g0_i1"));
        let merged_seq: String = text.lines().filter(|l| !l.starts_with('>')).collect();
        assert_eq!(merged_seq.len(), 200 + N_GAP + 200);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
