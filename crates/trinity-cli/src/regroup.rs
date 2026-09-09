//! 输出层严格家族拆分 pass（2026-09-07, python 原型 regroup_proto.py 内部化）。
//!
//! 机制背景: QM4 蛋白层诊断（docs/opt-log.md「严格家族拆分路线」）——基因内
//! "最长转录本"常被脏延伸/旁系混入的转录本抢位（换代表全量分类
//! better 2207 / worse 0, 拆分只赚不亏）。原型证明只用**基因内两两关系**即可
//! （跨基因归并收益为 0）, 故本 pass 是 split-only: 只拆旧基因, 不并跨基因。
//!
//! 判据（不依赖任何参照, 纯序列证据）: 基因内每对转录本算三个特征——
//! shared（比对共享碱基 = match + mism）、mism 密度（块内替换密度, **块状
//! indel/外显子跳越不进分母**）、双侧覆盖（比对区间并集 / 全长）。同位点对
//! 密度 p95=0.28% vs 异位点对 p50=4.2%, 密度是强判别特征:
//! - iso 边: 双侧 cov ≥ C 且 shared ≥ 200 且 dens < 0.005（同一 isoform 家族）;
//! - para 边（cannot-link）: shared ≥ 200 且 dens ≥ 0.02（旁系同源, 禁止同簇）;
//! - 其余: 无边（保守）。
//!
//! 内部比对（替代原型的 minimap2 all-vs-all）: 锚门控（共享 distinct 25-mer
//! <8 直接跳过, 索引侧出现 >4 次的 k-mer 丢弃防重复区干扰）+ 多对角带状局部
//! 比对（锚直方图取前 ≤5 条对角, 每条 ≥3 锚; 带宽 ±64, match +1 / mismatch
//! -1 / 线性 gap -1, 只统计两侧皆非 gap 的列）。比对范围钳制在带锚的窗口
//! 值域——锚是精确 25-mer 命中, 既定对角也定可信范围; 不钳制时局部比对会以
//! ~25% 命中的代价挤过相邻外显子的错框区把 dens 拖脏。对角差 ≤2×带宽 的
//! 选中带并链为一次 DP（旁系对小 indel 的累积漂移会把共享块挪出单管, 分带
//! 各自对齐会漏掉带间无锚的分歧段）。多带按锚数降序处理, 逐列占用位图精确
//! 去重（两侧皆非 gap 且未被先前带占用才计数; 区间两侧均覆盖 ≥80% 的带
//! 整条跳过）。较短侧正向与反向互补各跑一遍、双方向聚合（shared/mism 求和、
//! 区间并集）——对齐原型 all-vs-all 特征不分局的口径, 倒位重复家族单方向
//! 只有半边覆盖。
//!
//! 聚类: iso 边按 shared 降序贪心 union-find; 合并前检查两簇成员间是否存在
//! para 边（cannot-link, 拒绝则跳过并计数）。家族 = 簇 + 无边单例（每条转录本
//! 恰属于一个家族）, 编号 rg{:06} 按（成员最大长度降序, 成员名升序）——与原型
//! 一致保证可比。每家族输出最长代表（长度降序、id 升序 tiebreak）。
//!
//! reads 支持通道（Phase B, `--regroup-reads`）: 分组规则完全不变, 只改代表
//! 选择——argmax(支持度[<3 视为无证据], 长度, id)。支持度 = 原始 reads 的
//! k-mer 投票: 转录本 25-mer 索引（每 k-mer 出现在 >8 条转录本则丢弃, 重复
//! 序列防污染）, read 双链合计命中 ≥ 1/4 自身 k-mer 才投票给命中最多的
//! 转录本; 命中平手的候选以 read 内容哈希伪随机分摊（同 read 确定性）——
//! 包含关系是组装产物的常态, 集中式平手裁决会把家族票全灌给容器者/最前者,
//! 原型 minimap2 的逐 read 主命中天然分摊, 哈希分摊是对它的确定性近似。
//! QM4 mem9 实测: Rp 37.3 → 39.7 / P 49.1 → 53.0 / F1 39.2 → 43.3
//! （原型 minimap2 支持度: 40.4/53.5/44.3）; 换代表全量分类
//! better 4933 / worse 108（45.7:1, 原型 15.5:1）——长短颠倒的净收益为正,
//! 用真实读深压"长但无人读"的脏延伸, 代价是少量高噪声短者抢位。
//!
//! 产物（Trinity.fasta 同名侧车, 主 fasta/gene_trans_map 不动）:
//! `<name>.principal.fa`（>rg000001_c5_g1_i5, 60 列折行）与
//! `<name>.family_map`（rg000001\tc5_g1_i1, 每转录本一行）。
//!
//! principal 数量截断（`--principal-min-len N` / `--principal-min-supp M`,
//! 管线**转正默认 500/10**, `RegroupParams` 原语默认 0 = 关, 显式传 0/0 关回
//! 旧行为）: 门控只过滤 principal.fa——代表长度 / read 票数不达标的家族
//! 不落 principal（family_map 照旧全量, 分组与代表选择不变, worse=0 由构造
//! 保证）。QM4 画像依据（opt-log「principal 数量口径战役」）: 106,065 principal
//! 里 86% 是短货（真基因片段 + 无命中孤岛组件, ST 盲区真内容风险 ≤1%）;
//! 转正档 len500+supp10 组合保留 44,522, P 46.9→54+, F1 净升（丢代表级
//! ~5%）。票数门与代表选择同源（一 read 一票, **无 mapq 语义**——多映射区
//! 的真货不被冤枉, 与 minimap2 mapq≥10 口径的关键差异）; supp>0 隐含
//! --regroup-reads（无 reads 时 supp=0 会把 principal 剔空, 管线层防住——
//! 转正后裸 --regroup 也因此隐含 reads 通道）。example 探针默认 0/0
//! （探针无 reads 时 supp 门会剔空, 显式传参复现主线档位）。

use std::path::{Path, PathBuf};

use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use trinity_common::error::CommonError;

use crate::args::TrinityArgs;

#[derive(Debug, Clone)]
pub struct RegroupParams {
    /// 锚 k-mer 长（25 = 全管线一致窗口; 更短则随机命中泛滥）
    pub kmer_len: usize,
    /// 参与评估的最低序列长（两侧任一 < 此值跳过; 200 与 min_contig_length 对齐）
    pub min_pair_len: usize,
    /// 锚门控: 共享 distinct k-mer 数下限
    pub min_shared_kmers: usize,
    /// 索引侧 k-mer 出现次数 > 此值 → 丢弃（重复区防干扰）
    pub max_kmer_occ: u32,
    /// 最多参与比对的对角带数
    pub max_bands: usize,
    /// 单带最少锚数
    pub band_min_anchors: u32,
    /// 带宽（对角 ±）
    pub bandwidth: usize,
    /// 多带去重: 对齐区间任一侧被先前带覆盖 ≥ 此比例 → 丢弃
    pub band_overlap_drop: f64,
    /// iso 边双侧覆盖下限（--regroup-cov 暴露）
    pub cov: f32,
    /// iso 边 mismatch 密度上限
    pub iso_mism: f32,
    /// para 边 mismatch 密度下限
    pub para_mism: f32,
    /// iso/para 边共享碱基下限
    pub min_shared: usize,
    /// reads 通道: 支持度索引每 k-mer 允许的最大转录本数（超过即丢弃）
    pub max_kmer_tx: usize,
    /// reads 通道: 支持度索引 k-mer 长（独立于锚 k 的容错旋钮; 实测 19 与
    /// 25 等价——1% 错误率下 19-mer 只多回收 ~5% 窗口, 保持 25 与管线一致）
    pub support_kmer_len: usize,
    /// reads 通道: 代表选择的最低支持度（低于此视为无证据, 退回长度规则）
    pub min_support: u64,
    /// principal 输出截断: 代表长度下限（0 = 关）。不达标的家族不落
    /// principal.fa（family_map 照旧全量; 分组与代表选择不变）
    pub min_principal_len: usize,
    /// principal 输出截断: 代表 read 票数下限（0 = 关）。票数与代表选择
    /// 同源（一 read 一票, 无 mapq 语义）; 无 reads 时 supp=0, >0 会剔空
    pub min_principal_supp: u64,
}

impl Default for RegroupParams {
    fn default() -> Self {
        RegroupParams {
            kmer_len: 25,
            min_pair_len: 200,
            min_shared_kmers: 8,
            max_kmer_occ: 4,
            max_bands: 5,
            band_min_anchors: 3,
            bandwidth: 64,
            band_overlap_drop: 0.80,
            cov: 0.70,
            iso_mism: 0.005,
            para_mism: 0.02,
            min_shared: 200,
            max_kmer_tx: 8,
            support_kmer_len: 25,
            min_support: 3,
            min_principal_len: 0,
            min_principal_supp: 0,
        }
    }
}

#[derive(Debug, Default)]
pub struct RegroupStats {
    pub transcripts: usize,
    pub genes: usize,
    /// 基因内两侧均 ≥ min_pair_len 的候选对数
    pub pairs_total: usize,
    /// 过锚门控并完成比对的对数
    pub pairs_evaluated: usize,
    pub iso_edges: usize,
    pub para_edges: usize,
    pub blocked_merges: usize,
    pub families: usize,
    pub multi_member: usize,
    /// 成员落入 >1 个新家族的旧基因数
    pub genes_split: usize,
    /// reads 通道: 实际投出票的 read 数（0 = 通道未启用）
    pub reads_voted: u64,
    /// reads 通道: 支持度 ≥ min_support 的转录本数
    pub supp_tx: usize,
    /// principal 截断: 因长度/票数门控未落 principal.fa 的家族数
    pub principal_dropped: usize,
}

/// 统计行（-regroup: ...; 仿 -stitch 风格, maybe_regroup 与 example 共用）。
pub fn stats_line(st: &RegroupStats) -> String {
    let mut s = format!(
        "-regroup: tx {} gene {} → family {}（多成员 {}, 旧 gene 拆分 {}）; 评估对 {}（候选 {}）, iso 边 {}, para 边 {}, cannot-link 拒并 {}",
        st.transcripts,
        st.genes,
        st.families,
        st.multi_member,
        st.genes_split,
        st.pairs_evaluated,
        st.pairs_total,
        st.iso_edges,
        st.para_edges,
        st.blocked_merges
    );
    if st.reads_voted > 0 {
        s.push_str(&format!(
            "; reads 通道: 票 {} reads, {} 条转录本有支持",
            st.reads_voted, st.supp_tx
        ));
    }
    if st.principal_dropped > 0 {
        s.push_str(&format!(
            "; principal 截断: 保留 {}/{}（dropped {}）",
            st.families - st.principal_dropped,
            st.families,
            st.principal_dropped
        ));
    }
    s
}

/// 同名侧车路径: 剥最后一个扩展名后接后缀
/// （mem9.Trinity.fasta → mem9.Trinity.family_map）。
pub fn sidecar_path(fa: &Path, suffix: &str) -> PathBuf {
    let base = fa
        .file_stem()
        .and_then(|s| s.to_str())
        .filter(|s| !s.is_empty())
        .or_else(|| fa.file_name().and_then(|s| s.to_str()))
        .unwrap_or("trinity");
    fa.with_file_name(format!("{base}.{suffix}"))
}

/// orchestrate 收尾接线（--regroup 开启时）: 在 final_fa 旁落 principal/
/// family_map 侧车, 统计行 eprintln + 追加 .stage_timings（仿 stitch）。
/// 旗标关时**零副作用**（不读产物、不写任何文件、不写 .stage_timings）。
/// --regroup-reads 开启时把原始 reads（--left/--right/--single, 归一化前的
/// 原文件路径）一并传入走支持度通道。
pub fn maybe_regroup(
    final_fa: &Path,
    args: &TrinityArgs,
    stage_dir: &Path,
) -> Result<(), CommonError> {
    if !args.regroup {
        return Ok(());
    }
    let principal = sidecar_path(final_fa, "principal.fa");
    let map = sidecar_path(final_fa, "family_map");
    let params = RegroupParams {
        cov: args.regroup_cov,
        min_principal_len: args.regroup_principal_min_len,
        min_principal_supp: args.regroup_principal_min_supp,
        ..RegroupParams::default()
    };
    // reads 快照（沿 stitch 的做法: 输出层扫**原始** reads, 不碰归一化产物）
    let reads: Vec<PathBuf> = if args.regroup_reads {
        args.left
            .iter()
            .chain(&args.right)
            .chain(&args.single)
            .cloned()
            .collect()
    } else {
        Vec::new()
    };
    let st = trinity_common::timer::timed(stage_dir, "regroup", || {
        regroup_fasta(final_fa, &principal, &map, &reads, &params)
    })?;
    eprintln!("{}", stats_line(&st));
    if let Ok(mut fh) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(stage_dir.join(".stage_timings"))
    {
        use std::io::Write;
        let _ = writeln!(
            fh,
            "regroup_detail\tpairs={} iso={} para={} blocked={} genes={} families={} multi={} split={} reads={} supp_tx={} plen={} psupp={} pdropped={}",
            st.pairs_evaluated,
            st.iso_edges,
            st.para_edges,
            st.blocked_merges,
            st.genes,
            st.families,
            st.multi_member,
            st.genes_split,
            st.reads_voted,
            st.supp_tx,
            params.min_principal_len,
            params.min_principal_supp,
            st.principal_dropped
        );
    }
    Ok(())
}

/// 对 `in_fa` 执行严格家族拆分, 落盘 principal fasta 与 family map。
/// `reads` 非空时启用支持度通道（代表 = argmax(支持度, 长度, id)）;
/// 为空时行为与 Phase 1 完全一致（代表 = 最长）。
pub fn regroup_fasta(
    in_fa: &Path,
    out_principal: &Path,
    out_family_map: &Path,
    reads: &[PathBuf],
    params: &RegroupParams,
) -> Result<RegroupStats, CommonError> {
    let raw = std::fs::read(in_fa)?;
    let recs = parse_fasta(&raw)?;
    let mut stats = RegroupStats { transcripts: recs.len(), ..Default::default() };
    if recs.is_empty() {
        std::fs::write(out_principal, "")?;
        std::fs::write(out_family_map, "")?;
        return Ok(stats);
    }

    // ---- 0. reads 支持度通道（可选; 分组不依赖它）----
    let (supp, voted) = if reads.is_empty() {
        (FxHashMap::default(), 0u64)
    } else {
        count_support(&recs, reads, params)?
    };
    stats.reads_voted = voted;
    stats.supp_tx = supp.values().filter(|&&s| s >= params.min_support).count();
    // ---- 1. 基因分组（主规则 cX_gY_iZ; 兜底"最后一个 _i<数字> 前缀"）----
    let mut gene_of: Vec<String> = Vec::with_capacity(recs.len());
    let mut gene_ids: FxHashMap<String, u32> = FxHashMap::default();
    let mut genes: Vec<Vec<u32>> = Vec::new();
    for (i, (name, _)) in recs.iter().enumerate() {
        let g = parse_gene(name);
        let gi = *gene_ids.entry(g.clone()).or_insert_with(|| {
            genes.push(Vec::new());
            genes.len() as u32 - 1
        });
        genes[gi as usize].push(i as u32);
        gene_of.push(g);
    }
    stats.genes = genes.len();

    // ---- 2. 基因内两两特征 → 边（rayon 按基因外层并行）----
    struct Edge {
        a: u32,
        b: u32,
        shared: u32,
        iso: bool,
    }
    let per_gene: Vec<(usize, usize, Vec<Edge>)> = genes
        .par_iter()
        .map_init(Scratch::new, |sc, members| {
            let mut pairs_total = 0usize;
            let mut pairs_eval = 0usize;
            let mut edges: Vec<Edge> = Vec::new();
            // 长度降序（平局按输入序）: 逐一以"较长侧"建索引, 更短侧扫描
            let mut order = members.clone();
            order.sort_unstable_by(|&x, &y| {
                recs[y as usize].1.len().cmp(&recs[x as usize].1.len()).then(x.cmp(&y))
            });
            for (ii, &m) in order.iter().enumerate() {
                let m_rec = &recs[m as usize];
                if m_rec.1.len() < params.min_pair_len {
                    break; // 后续只会更短
                }
                let index = build_index(&m_rec.1, params.kmer_len, params.max_kmer_occ);
                for &s in &order[ii + 1..] {
                    let s_seq = &recs[s as usize].1;
                    if s_seq.len() < params.min_pair_len {
                        break;
                    }
                    pairs_total += 1;
                    let Some(f) = eval_pair(&m_rec.1, &index, s_seq, params, sc) else {
                        continue;
                    };
                    pairs_eval += 1;
                    let dens =
                        if f.shared > 0 { f.mism as f32 / f.shared as f32 } else { 1.0 };
                    // 端点按名字规范化（与原型 (a<b) 键一致, 保证边序可比）
                    let (a, b) = if m_rec.0 < recs[s as usize].0 { (m, s) } else { (s, m) };
                    if f.shared >= params.min_shared as u32 {
                        if f.cov_long.min(f.cov_short) >= params.cov && dens < params.iso_mism {
                            edges.push(Edge { a, b, shared: f.shared, iso: true });
                        } else if dens >= params.para_mism {
                            edges.push(Edge { a, b, shared: f.shared, iso: false });
                        }
                    }
                }
            }
            (pairs_total, pairs_eval, edges)
        })
        .collect();

    let mut iso_edges: Vec<&Edge> = Vec::new();
    let mut para_edges: Vec<&Edge> = Vec::new();
    for (pt, pe, edges) in &per_gene {
        stats.pairs_total += pt;
        stats.pairs_evaluated += pe;
        for e in edges {
            if e.iso {
                iso_edges.push(e);
            } else {
                para_edges.push(e);
            }
        }
    }
    stats.iso_edges = iso_edges.len();
    stats.para_edges = para_edges.len();

    // ---- 3. 聚类: iso 按 shared 降序贪心合并 + para cannot-link（同原型）----
    let n = recs.len();
    let mut uf = UnionFind::new(n);
    let mut para_adj: Vec<Vec<u32>> = vec![Vec::new(); n];
    for e in &para_edges {
        para_adj[e.a as usize].push(e.b);
        para_adj[e.b as usize].push(e.a);
    }
    // 序: (shared 降序, 名字降序) —— 原型 reverse=True 的 (s,a,b) 元组序
    iso_edges.sort_by(|x, y| {
        y.shared.cmp(&x.shared)
            .then_with(|| recs[y.a as usize].0.cmp(&recs[x.a as usize].0))
            .then_with(|| recs[y.b as usize].0.cmp(&recs[x.b as usize].0))
    });
    let mut members: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    for e in &iso_edges {
        let (ra, rb) = (uf.find(e.a), uf.find(e.b));
        if ra == rb {
            continue;
        }
        // cannot-link: 扫较小簇成员的 para 邻接, 对端落在大簇 → 拒并
        let (small_r, big_r) = if uf.size[ra as usize] <= uf.size[rb as usize] {
            (ra, rb)
        } else {
            (rb, ra)
        };
        let small: Vec<u32> = members
            .get(&small_r)
            .cloned()
            .unwrap_or_else(|| vec![small_r]);
        let big: FxHashSet<u32> = members
            .get(&big_r)
            .map(|v| v.iter().copied().collect())
            .unwrap_or_else(|| {
                let mut s = FxHashSet::default();
                s.insert(big_r);
                s
            });
        let hit = small
            .iter()
            .any(|&x| para_adj[x as usize].iter().any(|&y| big.contains(&y)));
        if hit {
            stats.blocked_merges += 1;
            continue;
        }
        uf.union(e.a, e.b);
        let root = uf.find(e.a);
        let mut merged = members.remove(&small_r).unwrap_or_else(|| vec![small_r]);
        merged.extend(members.remove(&big_r).unwrap_or_else(|| vec![big_r]));
        members.insert(root, merged);
    }

    // ---- 4. 家族 = 簇 + 无边单例; 排序编号; 代表 = 最长（id 升序 tiebreak）----
    let mut by_root: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    for t in 0..n as u32 {
        by_root.entry(uf.find(t)).or_default().push(t);
    }
    struct Family {
        members: Vec<u32>, // 名字升序（map 输出序）
        rep: u32,
        max_len: usize,
        min_name: String,
    }
    let mut fams: Vec<Family> = by_root
        .into_values()
        .map(|mut mem| {
            mem.sort_unstable_by(|&x, &y| recs[x as usize].0.cmp(&recs[y as usize].0));
            // 代表 = argmax(支持度[< min_support 视为 -1], 长度, id 升序):
            // reads 通道用真实读深压"长但无人读"的脏延伸; 无支持证据时
            // 全部并列 -1, 退化为 Phase 1 的最长规则。
            let supp_key = |t: u32| -> i64 {
                match supp.get(&t).copied().unwrap_or(0) {
                    s if s >= params.min_support => s as i64,
                    _ => -1,
                }
            };
            let rep = *mem
                .iter()
                .max_by(|&&x, &&y| {
                    supp_key(x)
                        .cmp(&supp_key(y))
                        .then_with(|| {
                            recs[x as usize]
                                .1
                                .len()
                                .cmp(&recs[y as usize].1.len())
                        })
                        .then_with(|| recs[y as usize].0.cmp(&recs[x as usize].0))
                })
                .unwrap();
            Family {
                max_len: mem.iter().map(|&t| recs[t as usize].1.len()).max().unwrap(),
                min_name: recs[mem[0] as usize].0.clone(),
                rep,
                members: mem,
            }
        })
        .collect();
    fams.sort_unstable_by(|x, y| {
        y.max_len.cmp(&x.max_len).then_with(|| x.min_name.cmp(&y.min_name))
    });
    stats.families = fams.len();
    stats.multi_member = fams.iter().filter(|f| f.members.len() > 1).count();

    // ---- 5. 落盘 ----
    let mut map_buf = String::new();
    let mut fa_buf = String::new();
    for (i, fam) in fams.iter().enumerate() {
        let gid = format!("rg{:06}", i + 1);
        for &t in &fam.members {
            map_buf.push_str(&gid);
            map_buf.push('\t');
            map_buf.push_str(&recs[t as usize].0);
            map_buf.push('\n');
        }
        // principal 截断门（--principal-min-len / --principal-min-supp, 默认 0=关）:
        // 只过滤 principal.fa 的输出行——分组与代表选择不变（worse=0 由构造保证）,
        // family_map 保持全量。长度门比较代表序列长; 票数门与代表选择同源
        // （无 reads 时 supp=0, 管线层由蕴含链防住剔空, 见 args.rs）。
        let rep_supp = supp.get(&fam.rep).copied().unwrap_or(0);
        if recs[fam.rep as usize].1.len() < params.min_principal_len
            || rep_supp < params.min_principal_supp
        {
            stats.principal_dropped += 1;
            continue;
        }
        let rep = &recs[fam.rep as usize];
        fa_buf.push('>');
        fa_buf.push_str(&gid);
        fa_buf.push('_');
        fa_buf.push_str(&rep.0);
        fa_buf.push('\n');
        for chunk in rep.1.chunks(60) {
            fa_buf.push_str(std::str::from_utf8(chunk).map_err(|_| {
                CommonError::Parse("regroup: 非 UTF-8 序列".to_string())
            })?);
            fa_buf.push('\n');
        }
    }
    std::fs::write(out_principal, fa_buf)?;
    std::fs::write(out_family_map, map_buf)?;

    // 旧 gene 拆分计数（成员落入 >1 家族）
    for members in &genes {
        let mut roots: Vec<u32> = members.iter().map(|&t| uf.find(t)).collect();
        roots.sort_unstable();
        roots.dedup();
        if roots.len() > 1 {
            stats.genes_split += 1;
        }
    }
    Ok(stats)
}

// ---------------------------------------------------------------------------
// 基因名解析
// ---------------------------------------------------------------------------

/// 主规则 cX_gY_iZ（orchestrate::write_gene_trans_map 同一正则）; 兜底
/// "最后一个 `_i<数字>`"前缀（多 k 等非规范命名的分组兜底）; 再不匹配以
/// 整名自成一基因（单成员家族）。
fn parse_gene(name: &str) -> String {
    if let Some((gene, _)) = crate::orchestrate::split_gene_trans(name) {
        return gene;
    }
    if let Some(i) = name.rfind("_i") {
        let tail = &name[i + 2..];
        if !tail.is_empty() && tail.bytes().all(|b| b.is_ascii_digit()) && i > 0 {
            return name[..i].to_string();
        }
    }
    name.to_string()
}

// ---------------------------------------------------------------------------
// 特征: 锚门控 + 多对角带状比对
// ---------------------------------------------------------------------------

/// 一对转录本的聚合特征（取 shared 更大的方向）。
#[derive(Debug)]
struct PairFeat {
    shared: u32,
    mism: u32,
    cov_long: f32,
    cov_short: f32,
}

/// 线程内复用的比对 scratch（dp 得分带 / 锚直方图 / 门控键集 / 逐列占用位图）。
struct Scratch {
    dp: Vec<i16>,
    hist: FxHashMap<i64, BandHist>,
    keys: FxHashSet<u64>,
    claim_a: Vec<u64>,
    claim_b: Vec<u64>,
}

impl Scratch {
    fn new() -> Self {
        Scratch {
            dp: Vec::new(),
            hist: FxHashMap::default(),
            keys: FxHashSet::default(),
            claim_a: Vec::new(),
            claim_b: Vec::new(),
        }
    }
}

/// 较长序列的 k-mer 索引: key → 出现位置（出现 > max_occ 次的键丢弃,
/// 防重复序列把锚堆到同一直方图桶）。含 N 的窗口不入键。
fn build_index(seq: &[u8], k: usize, max_occ: u32) -> FxHashMap<u64, Vec<u32>> {
    let mut idx: FxHashMap<u64, Vec<u32>> = FxHashMap::default();
    if seq.len() < k {
        return idx;
    }
    let mask = (1u64 << (2 * k as u32)) - 1;
    let mut valid = 0usize;
    let mut v: u64 = 0;
    for (i, &c) in seq.iter().enumerate() {
        let x = match c {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                valid = 0;
                v = 0;
                continue;
            }
        };
        v = ((v << 2) | x as u64) & mask;
        valid += 1;
        if valid >= k {
            idx.entry(v).or_default().push((i + 1 - k) as u32);
        }
    }
    idx.retain(|_, pos| pos.len() as u32 <= max_occ);
    idx
}

/// 单方向聚合: 带比对结果累加（shared/mism）+ 两侧对齐区间并集。
struct DirAgg {
    shared: u32,
    mism: u32,
    iv_long: Vec<(u32, u32)>,
    iv_short: Vec<(u32, u32)>,
}

/// 一对转录本双向评估: 较短侧正向与反向互补各跑一遍（锚门控 → 锚直方图
/// → 前 ≤max_bands 条对角带状局部比对, 多带去重）, 双方向**聚合**（shared/
/// mism 求和、区间并集）——对齐原型的跨链口径（其 all-vs-all 特征不分局,
/// 正反链证据同入一对比对; 倒位重复家族单方向只有一半覆盖）。
/// 双方向均未过门控 → None。
fn eval_pair(
    long: &[u8],
    index: &FxHashMap<u64, Vec<u32>>,
    short: &[u8],
    p: &RegroupParams,
    sc: &mut Scratch,
) -> Option<PairFeat> {
    let fwd = run_direction(long, index, short, p, sc);
    let mut rc_seq: Vec<u8> = Vec::with_capacity(short.len());
    rc_seq.extend(short.iter().rev().map(|&c| match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }));
    let rcagg = run_direction(long, index, &rc_seq, p, sc);
    let (mut shared, mut mism, mut iv_long, mut iv_short) =
        (0u32, 0u32, Vec::new(), Vec::new());
    for dir in [fwd, rcagg].into_iter().flatten() {
        shared += dir.shared;
        mism += dir.mism;
        iv_long.extend(dir.iv_long);
        iv_short.extend(dir.iv_short);
    }
    if shared == 0 && iv_long.is_empty() {
        return None;
    }
    let cov_long = union_len(&iv_long) as f32 / long.len() as f32;
    let cov_short = union_len(&iv_short) as f32 / short.len() as f32;
    Some(PairFeat { shared, mism, cov_long, cov_short })
}

/// 一个方向的完整流程: 扫描序列窗口 → 锚直方图 + 门控键集 → 取带 → 比对。
/// 门控不过（distinct 共享 k-mer < min_shared_kmers）→ None。
fn run_direction(
    long: &[u8],
    index: &FxHashMap<u64, Vec<u32>>,
    scan: &[u8],
    p: &RegroupParams,
    sc: &mut Scratch,
) -> Option<DirAgg> {
    let Scratch { hist, keys, dp, claim_a, claim_b, .. } = sc;
    let k = p.kmer_len;
    if scan.len() < k {
        return None;
    }
    hist.clear();
    keys.clear();
    let mask = (1u64 << (2 * k as u32)) - 1;
    let mut valid = 0usize;
    let mut v: u64 = 0;
    for (i, &c) in scan.iter().enumerate() {
        let x = match c {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                valid = 0;
                v = 0;
                continue;
            }
        };
        v = ((v << 2) | x as u64) & mask;
        valid += 1;
        if valid < k {
            continue;
        }
        let j = (i + 1 - k) as i64;
        if let Some(pos) = index.get(&v) {
            keys.insert(v);
            for &pi in pos {
                let e = hist.entry(pi as i64 - j).or_insert(BandHist {
                    anchors: 0,
                    i_lo: pi,
                    i_hi: pi,
                });
                e.anchors += 1;
                e.i_lo = e.i_lo.min(pi);
                e.i_hi = e.i_hi.max(pi);
            }
        }
    }
    if keys.len() < p.min_shared_kmers {
        return None;
    }
    // 前 ≤max_bands 条对角（锚数降序, 平局对角升序——确定性）
    let mut bands: Vec<(u32, i64, u32, u32)> = hist
        .iter()
        .filter(|(_, h)| h.anchors >= p.band_min_anchors)
        .map(|(&d, h)| (h.anchors, d, h.i_lo, h.i_hi))
        .collect();
    bands.sort_unstable_by(|x, y| y.0.cmp(&x.0).then_with(|| x.1.cmp(&y.1)));
    bands.truncate(p.max_bands);
    // 带内锚跨扩展: 管内（|d'-d| ≤ 带宽）其他对角的锚同属该带证据, 跨度并入。
    // 小 indel（≤ 带宽）把共享块挪到邻近对角时, 两簇锚之间的分歧连接段
    // （无锚）才不会被两侧带各自剪掉——mismatch 照记, 密度口径与原型一致。
    let w = p.bandwidth as i64;
    let bands: Vec<(u32, i64, u32, u32)> = bands
        .into_iter()
        .map(|(n, d, _, _)| {
            let (mut lo, mut hi) = (u32::MAX, 0u32);
            for (&d2, h) in hist.iter() {
                if (d2 - d).abs() <= w && h.anchors >= p.band_min_anchors {
                    lo = lo.min(h.i_lo);
                    hi = hi.max(h.i_hi);
                }
            }
            (n, d, lo, hi)
        })
        .collect();
    // 链式并带: 对角差 ≤ 2×带宽 的选中带并入同一次 DP（管加宽为
    // [d_min-w, d_max+w]）。旁系对的小 indel 累积漂移把共享块挪出单管,
    // 分带各自对齐会漏掉带间无锚的分歧段（密度被系统性洗白）; 并链后
    // 一条带内 DP 连续走完, mismatch 全额入账。跨度仍受锚跨钳制, 端部
    // slop 不因管加宽而回潮。
    let chain_gap = 2 * w;
    let chain_cap = 4 * w;
    let mut by_d = bands.clone();
    by_d.sort_unstable_by_key(|b| b.1);
    let mut chains: Vec<(u32, i64, i64, u32, u32)> = Vec::new(); // (锚数, d_lo, d_hi, i_lo, i_hi)
    for (n, d, i_lo, i_hi) in by_d {
        match chains.last_mut() {
            Some(c) if d - c.2 <= chain_gap && d - c.2 + (c.2 - c.1) <= chain_cap => {
                c.0 += n;
                c.2 = d;
                c.3 = c.3.min(i_lo);
                c.4 = c.4.max(i_hi);
            }
            _ => chains.push((n, d, d, i_lo, i_hi)),
        }
    }
    // 处理序: 链锚数降序, 平局链首对角升序（确定性）
    chains.sort_unstable_by(|x, y| y.0.cmp(&x.0).then_with(|| x.1.cmp(&y.1)));

    let mut agg = DirAgg {
        shared: 0,
        mism: 0,
        iv_long: Vec::new(),
        iv_short: Vec::new(),
    };
    // 逐列占用位图: 锚数多的链先认领自己的列, 后续链只统计两侧均未占用
    // 的新列——重复拷贝产生的交叉带（区间重叠 79% 之类）不再重复计数,
    // 方向间 shared 比较才有意义。
    let words_a = long.len().div_ceil(64);
    let words_b = scan.len().div_ceil(64);
    claim_a.clear();
    claim_a.resize(words_a, 0);
    claim_b.clear();
    claim_b.resize(words_b, 0);
    for &(_, d_lo, d_hi, i_lo, i_hi) in &chains {
        let d_mid = (d_lo + d_hi) / 2;
        let w_eff = (w + (d_hi - d_lo + 1) / 2) as usize;
        let Some(al) = banded_sw(
            long,
            scan,
            k,
            BandSel { d: d_mid, w: w_eff, i_lo, i_hi },
            dp,
            claim_a,
            claim_b,
        ) else {
            continue;
        };
        let (la, lb) = (al.a_hi - al.a_lo, al.b_hi - al.b_lo);
        if la > 0 && lb > 0 {
            // 多带去重（快速路径）: 两侧对齐区间均被先前带并集覆盖 ≥80% → 跳过
            let ca = covered_len(&agg.iv_long, al.a_lo, al.a_hi);
            let cb = covered_len(&agg.iv_short, al.b_lo, al.b_hi);
            if ca as f64 >= p.band_overlap_drop * la as f64
                && cb as f64 >= p.band_overlap_drop * lb as f64
            {
                continue;
            }
        }
        agg.shared += al.matches + al.mismatches;
        agg.mism += al.mismatches;
        agg.iv_long.push((al.a_lo, al.a_hi));
        agg.iv_short.push((al.b_lo, al.b_hi));
    }
    Some(agg)
}

/// 一条带的局部比对结果: 非 gap 列统计 + 两侧对齐区间（半开区间）。
#[derive(Debug)]
struct BandAlign {
    matches: u32,
    mismatches: u32,
    a_lo: u32,
    a_hi: u32,
    b_lo: u32,
    b_hi: u32,
}

fn bit_test(bits: &[u64], i: usize) -> bool {
    bits[i / 64] & (1u64 << (i % 64)) != 0
}

fn bit_set(bits: &mut [u64], i: usize) {
    bits[i / 64] |= 1u64 << (i % 64);
}

/// 一条对角的锚直方图项: 锚数 + 锚在较长侧的窗口起值域（比对区间钳制用）。
struct BandHist {
    anchors: u32,
    i_lo: u32,
    i_hi: u32,
}

/// 带状 DP 的几何（对角/带宽/行步长/窗界）与格子读取。
struct BandGeo {
    d: i64,
    w: usize,
    stride: usize,
    lb: usize,
    x_lo: usize,
    x_hi: usize,
}

impl BandGeo {
    /// 行 x 的带内 y 域 [y0, y1]（clamp 到 [0, lb]; 空行 = 带不与序列相交）。
    fn row_range(&self, x: usize) -> (usize, usize) {
        let lo = (x as i64 - self.d - self.w as i64).max(0);
        let hi = (x as i64 - self.d + self.w as i64).min(self.lb as i64);
        if lo > hi {
            (1, 0) // 空
        } else {
            (lo as usize, hi as usize)
        }
    }

    /// 读 H(x, y): y=0 是边界（恒 0, 不落盘）; 带外/窗外（未计算）→ None。
    /// 窗外 x 的 None 在 DP 侧 unwrap_or(0) = 窗缘自由重启（局部比对语义）。
    fn cell(&self, dp: &[i16], x: usize, y: usize) -> Option<i16> {
        if y == 0 || x == 0 {
            return Some(0);
        }
        if x < self.x_lo || x > self.x_hi {
            return None;
        }
        let (y0, y1) = self.row_range(x);
        if y < y0 || y > y1 {
            return None;
        }
        Some(dp[x * self.stride + (y - y0)])
    }
}

/// 对角 d（锚位置差 i-j）上的带状局部比对（Smith-Waterman, |x-y-d| ≤ w;
/// 一次带比对的选择参数（对角/带宽/锚跨窗, 由链结构导出）。
struct BandSel {
    d: i64,
    w: usize,
    i_lo: u32,
    i_hi: u32,
}

/// match +1 / mismatch -1 / 线性 gap -1）。比对范围钳制在该带锚的窗口值域
/// [i_lo, i_hi+k)——锚是精确 25-mer 命中, 既定对角也定可信范围; 不钳制时
/// 局部比对会以 ~25% 命中的代价挤过相邻外显子的错框区, 把 dens 拖脏。
/// 只统计两侧皆非 gap 的列（matches / mismatches）, 记录对齐区间; 块状
/// indel 走带内 gap 列、不进 shared 与密度分母。最优得分 0 → None。
fn banded_sw(
    a: &[u8],
    b: &[u8],
    k: usize,
    sel: BandSel,
    dp: &mut Vec<i16>,
    claim_a: &mut [u64],
    claim_b: &mut [u64],
) -> Option<BandAlign> {
    let BandSel { d, w, i_lo, i_hi } = sel;
    let (la, lb) = (a.len(), b.len());
    let stride = 2 * w + 1;
    let need = (la + 1) * stride;
    if dp.len() < need {
        dp.resize(need, 0);
    }
    let dp = &mut dp[..need];
    // 前缀坐标 x 消费 a[x-1]: 锚窗口 [i, i+k) → x ∈ [i_lo+1, i_hi+k]
    let geo = BandGeo {
        d,
        w,
        stride,
        lb,
        x_lo: (i_lo as usize + 1).max(1),
        x_hi: (i_hi as usize + k).min(la),
    };
    if geo.x_lo > geo.x_hi {
        return None;
    }
    let mut best = 0i16;
    let mut best_x = 0usize;
    let mut best_y = 0usize;
    for x in geo.x_lo..=geo.x_hi {
        let (y0, y1) = geo.row_range(x);
        let mut y = y0.max(1);
        while y <= y1 {
            let s: i16 = if a[x - 1] == b[y - 1] { 1 } else { -1 };
            let mut h: i32 = geo.cell(dp, x - 1, y - 1).unwrap_or(0) as i32 + s as i32;
            if let Some(up) = geo.cell(dp, x - 1, y) {
                h = h.max(up as i32 - 1);
            }
            if let Some(left) = geo.cell(dp, x, y - 1) {
                h = h.max(left as i32 - 1);
            }
            let h = h.max(0) as i16;
            dp[x * stride + (y - y0)] = h;
            if h > best {
                best = h;
                best_x = x;
                best_y = y;
            }
            y += 1;
        }
    }
    if best <= 0 {
        return None;
    }
    // 回溯（diag > up > left 平局序——确定性）: 只数非 gap 列; 已被先前带
    // 占用的列（任一侧）不再计数——多带去重的精确形态。
    let (mut x, mut y) = (best_x, best_y);
    let (mut a_lo, mut b_lo) = (best_x as u32, best_y as u32);
    let (mut matches, mut mismatches) = (0u32, 0u32);
    while geo.cell(dp, x, y).is_some_and(|h| h > 0) {
        let h = geo.cell(dp, x, y).unwrap_or(0);
        let s: i16 = if a[x - 1] == b[y - 1] { 1 } else { -1 };
        let diag = geo.cell(dp, x - 1, y - 1).unwrap_or(0);
        if diag as i32 + s as i32 == h as i32 {
            let (pa, pb) = (x - 1, y - 1);
            if !bit_test(claim_a, pa) && !bit_test(claim_b, pb) {
                bit_set(claim_a, pa);
                bit_set(claim_b, pb);
                if s == 1 {
                    matches += 1;
                } else {
                    mismatches += 1;
                }
            }
            x -= 1;
            y -= 1;
            a_lo = x as u32;
            b_lo = y as u32;
        } else if geo.cell(dp, x - 1, y).is_some_and(|v| v as i32 - 1 == h as i32) {
            x -= 1;
            a_lo = x as u32;
        } else {
            y -= 1;
            b_lo = y as u32;
        }
    }
    if matches + mismatches == 0 {
        return None;
    }
    Some(BandAlign {
        matches,
        mismatches,
        a_lo,
        a_hi: best_x as u32,
        b_lo,
        b_hi: best_y as u32,
    })
}

/// 区间并集总长。
fn union_len(ivs: &[(u32, u32)]) -> usize {
    covered_len(ivs, 0, u32::MAX) as usize
}

/// [lo, hi) 被区间并集覆盖的长度（≤5 个区间, 排序合并后求交）。
fn covered_len(ivs: &[(u32, u32)], lo: u32, hi: u32) -> u32 {
    let mut sorted = ivs.to_vec();
    sorted.sort_unstable();
    let mut tot = 0u32;
    let mut cur: Option<(u32, u32)> = None;
    for &(s, e) in &sorted {
        cur = Some(match cur {
            None => (s, e),
            Some((cs, ce)) => {
                if s > ce {
                    tot += intersect(cs, ce, lo, hi);
                    (s, e)
                } else {
                    (cs, ce.max(e))
                }
            }
        });
    }
    if let Some((cs, ce)) = cur {
        tot += intersect(cs, ce, lo, hi);
    }
    tot
}

fn intersect(cs: u32, ce: u32, lo: u32, hi: u32) -> u32 {
    let s = cs.max(lo);
    let e = ce.min(hi);
    e.saturating_sub(s)
}

struct UnionFind {
    parent: Vec<u32>,
    size: Vec<u32>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind { parent: (0..n as u32).collect(), size: vec![1; n] }
    }
    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            self.parent[x as usize] = self.parent[self.parent[x as usize] as usize];
            x = self.parent[x as usize];
        }
        x
    }
    fn union(&mut self, a: u32, b: u32) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        let (big, small) = if self.size[ra as usize] >= self.size[rb as usize] {
            (ra, rb)
        } else {
            (rb, ra)
        };
        self.parent[small as usize] = big;
        self.size[big as usize] += self.size[small as usize];
    }
}

// ---------------------------------------------------------------------------
// reads 支持度通道（Phase B, --regroup-reads）
// ---------------------------------------------------------------------------

/// 转录本 k-mer → 转录本列表索引（CSR: hash 定位 + 有序段表）。
/// 每 k-mer 出现在 > max_kmer_tx 条转录本的条目在构建时丢弃——共享重复
/// 段的票会同时泼给整族, 污染支持度。
struct KmerTxIndex {
    pos: FxHashMap<u64, u32>,
    /// 第 i 个键的转录本段 = vals[starts[i]..starts[i+1]]（starts 多一个哨兵）
    starts: Vec<u32>,
    vals: Vec<u32>,
}

impl KmerTxIndex {
    fn get(&self, key: u64) -> &[u32] {
        match self.pos.get(&key) {
            Some(&i) => {
                &self.vals[self.starts[i as usize] as usize..self.starts[i as usize + 1] as usize]
            }
            None => &[],
        }
    }
}

/// 构建支持度索引: 每条转录本先取 distinct k-mer（重复窗口去重, 防转录本
/// 内重复段自放大）, 全局排序分组后落 CSR; 组内转录本数 > max_kmer_tx 的键
/// 整组丢弃。rayon 分担逐转录本窗口提取与并行排序。
fn build_kmer_tx_index(recs: &[(String, Vec<u8>)], k: usize, max_kmer_tx: usize) -> KmerTxIndex {
    use rayon::slice::ParallelSliceMut;
    let kmer_tx: Vec<Vec<u64>> = recs
        .par_iter()
        .map(|(_, seq)| {
            let mut v: Vec<u64> = Vec::with_capacity(seq.len().saturating_sub(k - 1));
            if seq.len() >= k {
                let mask = (1u64 << (2 * k as u32)) - 1;
                let mut valid = 0usize;
                let mut cur: u64 = 0;
                for &c in seq {
                    let x = match c {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        b'T' => 3,
                        _ => {
                            valid = 0;
                            cur = 0;
                            continue;
                        }
                    };
                    cur = ((cur << 2) | x as u64) & mask;
                    valid += 1;
                    if valid >= k {
                        v.push(cur);
                    }
                }
            }
            v.sort_unstable();
            v.dedup();
            v
        })
        .collect();
    let total: usize = kmer_tx.iter().map(|v| v.len()).sum();
    let mut pairs: Vec<(u64, u32)> = Vec::with_capacity(total);
    for (tx, kmers) in kmer_tx.iter().enumerate() {
        for &key in kmers {
            pairs.push((key, tx as u32));
        }
    }
    drop(kmer_tx);
    pairs.par_sort_unstable();

    // 分组落 CSR（>max_kmer_tx 的组丢弃）
    let mut pos: FxHashMap<u64, u32> = FxHashMap::default();
    let mut starts: Vec<u32> = Vec::new();
    let mut vals: Vec<u32> = Vec::new();
    let mut i = 0usize;
    while i < pairs.len() {
        let key = pairs[i].0;
        let mut j = i + 1;
        while j < pairs.len() && pairs[j].0 == key {
            j += 1;
        }
        if j - i <= max_kmer_tx {
            pos.insert(key, starts.len() as u32);
            starts.push(vals.len() as u32);
            vals.extend(pairs[i..j].iter().map(|&(_, tx)| tx));
        }
        i = j;
    }
    starts.push(vals.len() as u32);
    KmerTxIndex { pos, starts, vals }
}

/// reads 文件的序列字节区间收集（fq 按首行 '@' / fa 按首行 '>' 自动识别;
/// fq 取 4 行组第 2 行, fa 逐非头行取——多行 fa 退化为逐行 read, 输入为
/// 单行 fa 时无差）。gzip 由 open_maybe_gz 魔数嗅探透明解压。
fn collect_read_ranges(data: &[u8]) -> Vec<(usize, usize)> {
    let mut is_fq: Option<bool> = None;
    let mut seqs: Vec<(usize, usize)> = Vec::new();
    let mut pos = 0usize;
    let mut line_no = 0usize;
    while pos < data.len() {
        let end = data[pos..]
            .iter()
            .position(|&b| b == b'\n')
            .map(|i| pos + i)
            .unwrap_or(data.len());
        let content_end = if end > pos && data[end - 1] == b'\r' { end - 1 } else { end };
        let is_blank = content_end == pos;
        let first = if !is_blank { data[pos] } else { 0 };
        if is_fq.is_none() && !is_blank {
            is_fq = Some(first == b'@');
        }
        match is_fq {
            Some(true) if line_no % 4 == 1 && !is_blank => seqs.push((pos, content_end)),
            Some(false) if !is_blank && first != b'>' => seqs.push((pos, content_end)),
            _ => {}
        }
        line_no += 1;
        pos = end + 1;
    }
    seqs
}

/// 单条 read 投票: 正链 + 反向互补双方向合计, 每个有效窗口对其命中的全部
/// 转录本记 1; 命中数 ≥ 自身窗口数的一半才成为候选, 命中最多者胜, 平手
/// 候选按 read 内容哈希分摊（见模块文档: 集中式裁决在包含关系家族会把票
/// 全灌给单个容器者, 与原型的逐 read 主命中分布相悖）。
fn vote_one_read(
    seq: &[u8],
    k: usize,
    index: &KmerTxIndex,
    hits: &mut FxHashMap<u32, u32>,
    rc_buf: &mut Vec<u8>,
) -> Option<u32> {
    if seq.len() < k {
        return None;
    }
    hits.clear();
    let mask = (1u64 << (2 * k as u32)) - 1;
    for flip in [false, true] {
        let strand: &[u8] = if flip {
            rc_buf.clear();
            rc_buf.extend(seq.iter().rev().map(|&c| match c {
                b'A' | b'a' => b'T',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                b'T' | b't' => b'A',
                _ => b'N',
            }));
            rc_buf
        } else {
            seq
        };
        let mut valid = 0usize;
        let mut cur: u64 = 0;
        for &c in strand {
            let x = match c {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => {
                    valid = 0;
                    cur = 0;
                    continue;
                }
            };
            cur = ((cur << 2) | x as u64) & mask;
            valid += 1;
            if valid >= k {
                for &tx in index.get(cur) {
                    *hits.entry(tx).or_insert(0) += 1;
                }
            }
        }
    }
    // 阈值 1/4: 横跨家族分歧区的 read 票数会分裂到多条转录本, 半数门槛会
    // 整票丢弃（minimap2 主命中容错语义下这些 read 仍投票）; 1/4 兼顾防误投。
    let threshold = ((seq.len() - k + 1).div_ceil(4)) as u32;
    // 平手集合内以 read 内容哈希定胜者——同 read 确定性, 不同 read 伪随机
    // 分摊（≈按比例分配选票）。三轮实验定案: 任何"集中到优选候选"的平票
    // 策略（最小 id 31.2 / 带内最长 34.0-34.2）都会饿死其余候选——分摊是
    // 简单策略最优; 与 minimap2 支持度（37.08）的剩余差距来自错配容忍打分
    // 本身, 需比对级归属（另立工程）。
    let top = hits.values().filter(|&&n| n >= threshold).max().copied();
    let top = top?;
    let mut tie: Vec<u32> = hits
        .iter()
        .filter(|&(_, &n)| n == top)
        .map(|(&tx, _)| tx)
        .collect();
    match tie.len() {
        0 => unreachable!("top 存在则 tie 非空"),
        1 => Some(tie[0]),
        _ => {
            tie.sort_unstable();
            let mut h: u64 = 0xcbf29ce484222325;
            for &b in seq {
                h = (h ^ b as u64).wrapping_mul(0x100000001b3);
            }
            Some(tie[(h % tie.len() as u64) as usize])
        }
    }
}

/// 支持度计数: 逐文件整读入内存 → rayon 分块投票 → 局部计数归并。
/// 返回（转录本下标 → 票数, 投票 read 总数）。
fn count_support(
    recs: &[(String, Vec<u8>)],
    reads: &[PathBuf],
    params: &RegroupParams,
) -> Result<(FxHashMap<u32, u64>, u64), CommonError> {
    let index = build_kmer_tx_index(recs, params.support_kmer_len, params.max_kmer_tx);
    let mut supp: FxHashMap<u32, u64> = FxHashMap::default();
    let mut voted = 0u64;
    for path in reads {
        let mut reader = trinity_common::io_util::open_maybe_gz(path)?;
        let mut data = Vec::new();
        std::io::Read::read_to_end(&mut reader, &mut data)?;
        let seqs = collect_read_ranges(&data);
        let partials: Vec<(FxHashMap<u32, u64>, u64)> = seqs
            .par_chunks(4096)
            .map(|chunk| {
                let mut local: FxHashMap<u32, u64> = FxHashMap::default();
                let mut hits: FxHashMap<u32, u32> = FxHashMap::default();
                let mut rc_buf: Vec<u8> = Vec::new();
                let mut n = 0u64;
                for &(s, e) in chunk {
                    if let Some(tx) = vote_one_read(
                        &data[s..e],
                        params.support_kmer_len,
                        &index,
                        &mut hits,
                        &mut rc_buf,
                    ) {
                        *local.entry(tx).or_insert(0) += 1;
                        n += 1;
                    }
                }
                (local, n)
            })
            .collect();
        for (local, n) in partials {
            for (tx, c) in local {
                *supp.entry(tx).or_insert(0) += c;
            }
            voted += n;
        }
    }
    Ok((supp, voted))
}

/// fasta 字节解析为 (name, seq)（名字取首 token; 序列大写化; 空序列跳过）。
fn parse_fasta(data: &[u8]) -> Result<Vec<(String, Vec<u8>)>, CommonError> {
    let mut out = Vec::new();
    let mut name: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();
    for line in data.split(|&b| b == b'\n') {
        let line = if line.last() == Some(&b'\r') { &line[..line.len() - 1] } else { line };
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if let Some(n) = name.take() {
                if !seq.is_empty() {
                    out.push((n, std::mem::take(&mut seq)));
                }
            }
            name = Some(String::from_utf8_lossy(&line[1..])
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string());
        } else {
            seq.extend(line.iter().map(|&b| b.to_ascii_uppercase()));
        }
    }
    if let Some(n) = name {
        if !seq.is_empty() {
            out.push((n, seq));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmpdir(tag: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("regroup_{tag}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    /// xorshift31 伪随机 ACGT 序列（不同 seed 序列间无偶然 25-mer 共享）。
    fn rand_seq(seed: u32, n: usize) -> Vec<u8> {
        let mut x = seed | 1;
        let mut s = Vec::with_capacity(n);
        for _ in 0..n {
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;
            s.push(b"ACGT"[(x % 4) as usize]);
        }
        s
    }

    /// 每 step 位换成一个**确定不同**的碱基（弥散分歧; step=30 → 密度恰 1/30
    /// ≈ 3.3%, 不靠随机换到同碱基）。
    fn mutate(seq: &[u8], step: usize) -> Vec<u8> {
        let mut out = seq.to_vec();
        for i in (step..seq.len()).step_by(step) {
            out[i] = match seq[i] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
        out
    }

    fn revcomp(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&c| match c {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect()
    }

    fn write_fa(dir: &Path, recs: &[(&str, Vec<u8>)]) -> std::path::PathBuf {
        let fa = dir.join("in.fasta");
        let mut txt = String::new();
        for (n, s) in recs {
            txt.push_str(&format!(">{n}\n{}\n", String::from_utf8_lossy(s)));
        }
        std::fs::write(&fa, txt).unwrap();
        fa
    }

    fn run(dir: &Path, recs: &[(&str, Vec<u8>)], cov: f32) -> RegroupStats {
        let fa = write_fa(dir, recs);
        let p = RegroupParams { cov, ..RegroupParams::default() };
        regroup_fasta(&fa, &dir.join("p.fa"), &dir.join("f.map"), &[], &p).unwrap()
    }

    /// 带 reads 支持度通道的运行: reads 文件路径相对 dir。
    fn run_reads(
        dir: &Path,
        recs: &[(&str, Vec<u8>)],
        cov: f32,
        reads: &[std::path::PathBuf],
    ) -> RegroupStats {
        let fa = write_fa(dir, recs);
        let p = RegroupParams { cov, ..RegroupParams::default() };
        regroup_fasta(&fa, &dir.join("p.fa"), &dir.join("f.map"), reads, &p).unwrap()
    }

    /// 写 fq（明文）。
    fn write_fq(dir: &Path, name: &str, seqs: &[Vec<u8>]) -> std::path::PathBuf {
        let p = dir.join(name);
        let mut txt = String::new();
        for (i, s) in seqs.iter().enumerate() {
            txt.push_str(&format!(
                "@r{i}\n{}\n+\n{}\n",
                String::from_utf8_lossy(s),
                "I".repeat(s.len())
            ));
        }
        std::fs::write(&p, txt).unwrap();
        p
    }

    /// 写 gzip fq（flate2 dev-dep, 与 trinity-common io_util 测试同款）。
    fn write_fq_gz(dir: &Path, name: &str, seqs: &[Vec<u8>]) -> std::path::PathBuf {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;
        let p = dir.join(name);
        let mut enc = GzEncoder::new(std::fs::File::create(&p).unwrap(), Compression::default());
        for (i, s) in seqs.iter().enumerate() {
            writeln!(
                enc,
                "@r{i}\n{}\n+\n{}",
                String::from_utf8_lossy(s),
                "I".repeat(s.len())
            )
            .unwrap();
        }
        enc.finish().unwrap();
        p
    }

    fn read_map(dir: &Path) -> Vec<(String, String)> {
        std::fs::read_to_string(dir.join("f.map"))
            .unwrap()
            .lines()
            .map(|l| {
                let mut it = l.split('\t');
                (
                    it.next().unwrap().to_string(),
                    it.next().unwrap().to_string(),
                )
            })
            .collect()
    }

    fn read_principal(dir: &Path) -> Vec<String> {
        std::fs::read_to_string(dir.join("p.fa"))
            .unwrap()
            .lines()
            .filter(|l| l.starts_with('>'))
            .map(|l| l[1..].to_string())
            .collect()
    }

    /// 直接评估一对序列（较先者作索引侧）, 断言用。
    fn feat_of(x: &[u8], y: &[u8]) -> Option<PairFeat> {
        let p = RegroupParams::default();
        let (long, short) = if x.len() >= y.len() { (x, y) } else { (y, x) };
        let index = build_index(long, p.kmer_len, p.max_kmer_occ);
        let mut sc = Scratch::new();
        eval_pair(long, &index, short, &p, &mut sc)
    }

    /// ① 同基因 isoform: 共享块完全一致 + 中段跳一外显子（大 indel）→
    /// iso 边, dens=0; gap 不进分母。
    #[test]
    fn isoform_pair_with_exon_skip_is_iso() {
        let e1 = rand_seq(11, 400);
        let e2 = rand_seq(23, 150);
        let e3 = rand_seq(37, 400);
        let i1 = [e1.clone(), e2.clone(), e3.clone()].concat(); // 950
        let i2 = [e1, e3].concat(); // 800（跳 e2）
        let f = feat_of(&i1, &i2).expect("锚门控应过");
        assert!(f.shared >= 700, "shared 应 ~800: {f:?}");
        assert_eq!(f.mism, 0, "共享块完全一致 dens=0");
        let dens = f.mism as f32 / f.shared as f32;
        assert!(dens < 0.005 && f.cov_long.min(f.cov_short) >= 0.7);
        // 簇层面: 同一家族, 最长者代表
        let d = tmpdir("iso");
        let st = run(&d, &[("c0_g1_i1", i1.clone()), ("c0_g1_i2", i2.clone())], 0.7);
        assert_eq!(st.families, 1);
        assert_eq!(st.iso_edges, 1);
        assert_eq!(st.multi_member, 1);
        assert_eq!(read_principal(&d), vec!["rg000001_c0_g1_i1".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    /// ② 共享块 ≥200bp 且 ~3% 弥散替换 → para 边, 不并簇。
    #[test]
    fn divergent_pair_is_para_and_split() {
        let x = rand_seq(101, 500);
        let y = mutate(&x, 30); // 每 30 位一换 ≈ 3.3%
        let f = feat_of(&x, &y).expect("锚门控应过");
        let dens = f.mism as f32 / f.shared as f32;
        assert!(f.shared >= 200, "shared: {f:?}");
        assert!(dens >= 0.02, "弥散分歧应判 para: dens={dens}");
        let d = tmpdir("para");
        let st = run(&d, &[("c0_g1_i1", x), ("c0_g1_i2", y)], 0.7);
        assert_eq!(st.families, 2, "para 不得并簇: {st:?}");
        assert_eq!(st.para_edges, 1);
        assert_eq!(st.iso_edges, 0);
        assert_eq!(st.genes_split, 1);
        let _ = std::fs::remove_dir_all(&d);
    }

    /// ③ 脏长 + 干净对: 干净俩并簇（代表取干净家族最长）, 脏者独立成家。
    #[test]
    fn dirty_long_split_from_clean_family() {
        let clean1 = rand_seq(211, 900);
        let clean2 = [clean1[..700].to_vec(), rand_seq(223, 100)].concat(); // 800
        let dirty = mutate(&[clean1[..900].to_vec(), rand_seq(227, 100).to_vec()].concat(), 30);
        let d = tmpdir("dirty");
        let st = run(
            &d,
            &[
                ("c0_g1_i1", clean1.clone()),
                ("c0_g1_i2", clean2),
                ("c0_g1_i3", dirty),
            ],
            0.7,
        );
        assert_eq!(st.families, 2, "脏者应被拆出: {st:?}");
        let map = read_map(&d);
        let fam_of = |t: &str| {
            map.iter().find(|(_, n)| n == t).unwrap().0.clone()
        };
        assert_eq!(fam_of("c0_g1_i1"), fam_of("c0_g1_i2"), "干净对同簇");
        assert_ne!(fam_of("c0_g1_i1"), fam_of("c0_g1_i3"), "脏者独立");
        let heads = read_principal(&d);
        assert_eq!(heads.len(), 2);
        // 干净家族代表 = 最长的干净转录本 i1（脏家族 i3 按长度先编号）
        assert!(
            heads.iter().any(|h| h.ends_with("_c0_g1_i1")),
            "干净家族代表应为 i1: {heads:?}"
        );
        assert!(
            heads.iter().any(|h| h.ends_with("_c0_g1_i3")),
            "脏家族代表应为 i3: {heads:?}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    /// ④ 无共享（锚门控拦下）→ 各自成家族。
    #[test]
    fn unrelated_pair_stays_separate() {
        let d = tmpdir("unrel");
        let st = run(
            &d,
            &[("c0_g1_i1", rand_seq(301, 600)), ("c0_g1_i2", rand_seq(401, 600))],
            0.7,
        );
        assert_eq!(st.families, 2);
        assert_eq!(st.pairs_evaluated, 0, "无共享不应评估: {st:?}");
        assert_eq!(st.iso_edges + st.para_edges, 0);
        let _ = std::fs::remove_dir_all(&d);
    }

    /// ⑤ 反向互补方向能识别: B = rc(A) → iso 同簇。
    #[test]
    fn reverse_complement_direction_detected() {
        let a = rand_seq(501, 600);
        let b = revcomp(&a);
        let d = tmpdir("rc");
        let st = run(&d, &[("c0_g1_i1", a), ("c0_g1_i2", b)], 0.7);
        assert_eq!(st.families, 1, "rc 方向应识别为 iso: {st:?}");
        assert_eq!(st.iso_edges, 1);
        // 等长 tiebreak: 代表取 id 升序
        assert_eq!(read_principal(&d), vec!["rg000001_c0_g1_i1".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    /// ⑥ cannot-link: A-B iso、B-C iso、A-C para → 不全并成一簇。
    /// B 同时与 A、C 各共享一块（A=s1+s3, B=s1+s2, C=s2+s3'）, 双块不相交
    /// 时 B 的双侧 cov 上限 0.5 → 本测试用 cov 0.3（机制测试, 非阈值测试）。
    #[test]
    fn cannot_link_blocks_transitive_merge() {
        let s1 = rand_seq(601, 1200); // A-B 干净共享
        let s2 = rand_seq(607, 1200); // B-C 干净共享
        let s3 = rand_seq(613, 400); // A-C 弥散共享
        let s3m = mutate(&s3, 30);
        let a = [s1.clone(), s3].concat(); // 1600
        let b = [s1, s2.clone()].concat(); // 2400
        let c = [s2, s3m].concat(); // 1600
        let d = tmpdir("cl");
        let st = run(
            &d,
            &[("c0_g1_i1", a), ("c0_g1_i2", b), ("c0_g1_i3", c)],
            0.3,
        );
        assert_eq!(st.iso_edges, 2, "A-B 与 B-C 均应 iso: {st:?}");
        assert_eq!(st.para_edges, 1, "A-C 应 para: {st:?}");
        assert_eq!(st.blocked_merges, 1, "A-C para 应拒并 B-C∪A-B 传递合并: {st:?}");
        assert_eq!(st.families, 2, "不得全并成一簇: {st:?}");
        let map = read_map(&d);
        let fam_of = |t: &str| map.iter().find(|(_, n)| n == t).unwrap().0.clone();
        // 序: B-C 先并（shared 平局按名字降序）, A-B 随后被 A-C para 拒并
        // → {B, C} 同簇, A 独立（传递合并被阻断）
        assert_eq!(fam_of("c0_g1_i2"), fam_of("c0_g1_i3"), "B-C 同簇");
        assert_ne!(fam_of("c0_g1_i1"), fam_of("c0_g1_i2"), "A 独立");
        let _ = std::fs::remove_dir_all(&d);
    }

    /// ⑦ family map 每转录本恰一次; 家族数 = principal 数。
    #[test]
    fn map_covers_each_transcript_once() {
        let g1a = rand_seq(701, 900);
        let g1b = [g1a[..700].to_vec(), rand_seq(709, 150)].concat();
        let recs = vec![
            ("c0_g1_i1", g1a),
            ("c0_g1_i2", g1b),
            ("c1_g1_i1", rand_seq(719, 500)),
            ("c2_g1_i1", rand_seq(727, 600)),
            ("c2_g1_i2", rand_seq(733, 600)),
        ];
        let d = tmpdir("map");
        let st = run(&d, &recs, 0.7);
        assert_eq!(st.transcripts, 5);
        let map = read_map(&d);
        assert_eq!(map.len(), 5, "每转录本恰一行");
        let mut names: Vec<&str> = map.iter().map(|(_, n)| n.as_str()).collect();
        names.sort_unstable();
        names.dedup();
        assert_eq!(names.len(), 5, "转录本不重不漏");
        let fams: std::collections::BTreeSet<&str> =
            map.iter().map(|(f, _)| f.as_str()).collect();
        assert_eq!(fams.len(), read_principal(&d).len(), "家族数 = principal 数");
        let _ = std::fs::remove_dir_all(&d);
    }

    /// ⑧ 确定性: 同输入两次运行输出逐字节一致。
    #[test]
    fn deterministic_output() {
        let x = rand_seq(809, 800);
        let y = mutate(&x, 45);
        let recs = vec![
            ("c0_g1_i1", x),
            ("c0_g1_i2", y),
            ("c3_g1_i1", rand_seq(811, 400)),
        ];
        let d1 = tmpdir("det1");
        let d2 = tmpdir("det2");
        run(&d1, &recs, 0.7);
        run(&d2, &recs, 0.7);
        assert_eq!(
            std::fs::read(d1.join("p.fa")).unwrap(),
            std::fs::read(d2.join("p.fa")).unwrap()
        );
        assert_eq!(
            std::fs::read(d1.join("f.map")).unwrap(),
            std::fs::read(d2.join("f.map")).unwrap()
        );
        let _ = std::fs::remove_dir_all(&d1);
        let _ = std::fs::remove_dir_all(&d2);
    }

    /// 非 cX_gY_iZ 命名兜底: 最后一个 `_i<数字>` 前缀分组。
    #[test]
    fn non_canonical_names_group_by_i_suffix() {
        let x = rand_seq(821, 700);
        let y = [x[..500].to_vec(), rand_seq(823, 100)].concat();
        let d = tmpdir("fallback");
        let st = run(&d, &[("weird_name_i1", x), ("weird_name_i2", y)], 0.7);
        assert_eq!(st.iso_edges, 1, "兜底分组下应识别 iso: {st:?}");
        assert_eq!(st.families, 1);
        assert_eq!(st.genes, 1);
        let _ = std::fs::remove_dir_all(&d);
    }

    /// 侧车命名: 剥最后扩展名。
    #[test]
    fn sidecar_paths_strip_extension() {
        assert_eq!(
            sidecar_path(Path::new("/a/b/mem9.Trinity.fasta"), "family_map"),
            PathBuf::from("/a/b/mem9.Trinity.family_map")
        );
        assert_eq!(
            sidecar_path(Path::new("/a/b/mem9.Trinity.fasta"), "principal.fa"),
            PathBuf::from("/a/b/mem9.Trinity.principal.fa")
        );
        assert_eq!(
            sidecar_path(Path::new("plain"), "family_map"),
            PathBuf::from("plain.family_map")
        );
    }

    // ---- Phase B: reads 支持度通道 ----

    /// 同基因一长一短（同簇）, reads 多数来自短者独有区 → 支持度压长度,
    /// 代表 = 短者。
    #[test]
    fn support_beats_length_for_representative() {
        let long = rand_seq(901, 900);
        let short = [long[..700].to_vec(), rand_seq(903, 100)].concat(); // 800
        let d = tmpdir("supp");
        // 5 条 read 横跨短者共享尾/独有尾: 独有窗口 ≥ 半数, 只投短者
        let reads: Vec<Vec<u8>> = (0..5).map(|i| short[650 + i..800].to_vec()).collect();
        let fq = write_fq(&d, "r.fq", &reads);
        let st = run_reads(
            &d,
            &[("c0_g1_i1", long.clone()), ("c0_g1_i2", short.clone())],
            0.7,
            &[fq],
        );
        assert_eq!(st.reads_voted, 5, "5 条 read 都应投票: {st:?}");
        assert_eq!(st.supp_tx, 1);
        let heads = read_principal(&d);
        assert_eq!(heads.len(), 1);
        assert!(
            heads[0].ends_with("_c0_g1_i2"),
            "支持度应压长度选短者: {heads:?}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    /// 支持度全部 < min_support（3）→ 退化为 Phase 1 最长规则。
    #[test]
    fn below_min_support_degrades_to_longest() {
        let long = rand_seq(911, 900);
        let short = [long[..700].to_vec(), rand_seq(913, 100)].concat();
        let d = tmpdir("supp2");
        let reads: Vec<Vec<u8>> = (0..2).map(|i| short[650 + i..800].to_vec()).collect();
        let fq = write_fq(&d, "r.fq", &reads);
        let st = run_reads(
            &d,
            &[("c0_g1_i1", long), ("c0_g1_i2", short)],
            0.7,
            &[fq],
        );
        assert_eq!(st.supp_tx, 0, "2 票 < 3 不算有支持: {st:?}");
        let heads = read_principal(&d);
        assert!(heads[0].ends_with("_c0_g1_i1"), "无证据退化为最长: {heads:?}");
        let _ = std::fs::remove_dir_all(&d);
    }

    /// 无 reads 路径 = Phase 1 等价锁: 空 reads 与不传 reads 输出逐字节一致。
    #[test]
    fn no_reads_path_matches_phase1() {
        let x = rand_seq(921, 800);
        let y = [x[..500].to_vec(), rand_seq(923, 100)].concat();
        let recs = vec![("c0_g1_i1", x), ("c0_g1_i2", y)];
        let d1 = tmpdir("noread1");
        let d2 = tmpdir("noread2");
        run(&d1, &recs, 0.7);
        run_reads(&d2, &recs, 0.7, &[]);
        assert_eq!(
            std::fs::read(d1.join("p.fa")).unwrap(),
            std::fs::read(d2.join("p.fa")).unwrap()
        );
        assert_eq!(
            std::fs::read(d1.join("f.map")).unwrap(),
            std::fs::read(d2.join("f.map")).unwrap()
        );
        let _ = std::fs::remove_dir_all(&d1);
        let _ = std::fs::remove_dir_all(&d2);
    }

    /// gz 与明文 fq 都能读; 反向互补 read 也能命中（rc 方向投票）。
    #[test]
    fn gz_and_plain_fq_with_rc_reads() {
        let long = rand_seq(931, 900);
        let short = [long[..700].to_vec(), rand_seq(933, 100)].concat();
        let d = tmpdir("gzrc");
        let fwd: Vec<Vec<u8>> = (0..3).map(|i| short[650 + i..800].to_vec()).collect();
        let rc: Vec<Vec<u8>> =
            (0..3).map(|i| revcomp(&short[650 + i..800])).collect();
        let plain = write_fq(&d, "plain.fq", &fwd);
        let gz = write_fq_gz(&d, "rc.fq.gz", &rc);
        let st = run_reads(
            &d,
            &[("c0_g1_i1", long), ("c0_g1_i2", short)],
            0.7,
            &[plain, gz],
        );
        assert_eq!(st.reads_voted, 6, "明文 + gz、正反链都应投票: {st:?}");
        let heads = read_principal(&d);
        assert!(heads[0].ends_with("_c0_g1_i2"), "rc read 命中短者: {heads:?}");
        let _ = std::fs::remove_dir_all(&d);
    }

    /// 索引重复 k-mer 剔除: 9 条转录本共享一段（> max_kmer_tx=8）→ 该段
    /// k-mer 全部不进索引, 落在共享段的 read 弃权（不泼票给全族）;
    /// 落在独有区的 read 正常投票。
    #[test]
    fn repeat_kmer_dropped_from_support_index() {
        let shared = rand_seq(941, 400);
        let names: Vec<String> = (0..9u32).map(|i| format!("c{i}_g1_i1")).collect();
        let uniques: Vec<Vec<u8>> = (0..9u32).map(|i| rand_seq(1001 + i, 400)).collect();
        let recs: Vec<(&str, Vec<u8>)> = names
            .iter()
            .zip(&uniques)
            .map(|(n, u)| (n.as_str(), [shared.clone(), u.clone()].concat()))
            .collect();
        let d = tmpdir("rep");
        // 3 条落在共享段（应弃权）+ 3 条落在 tx3 独有段（应投票 tx3）
        let shared_reads: Vec<Vec<u8>> =
            (0..3).map(|i| shared[100 + i..250 + i].to_vec()).collect();
        let uniq_reads: Vec<Vec<u8>> =
            (0..3).map(|i| uniques[3][100 + i..250 + i].to_vec()).collect();
        let mut all = shared_reads;
        all.extend(uniq_reads);
        let fq = write_fq(&d, "r.fq", &all);
        let st = run_reads(&d, &recs, 0.7, &[fq]);
        assert_eq!(st.reads_voted, 3, "共享段 read 应被 >8 剔除规则弃权: {st:?}");
        assert_eq!(st.supp_tx, 1);
        let _ = std::fs::remove_dir_all(&d);
    }

    /// principal 长度门: 只过滤 principal.fa（family_map 全量）, 边界 == 过。
    #[test]
    fn principal_min_len_gate_filters_output_only() {
        let d = tmpdir("plen");
        let fa = write_fa(
            &d,
            &[
                ("c0_g1_i1", rand_seq(951, 300)), // == 门槛, 过
                ("c1_g1_i1", rand_seq(953, 299)), // < 门槛, 剔
            ],
        );
        let p = RegroupParams {
            cov: 0.7,
            min_principal_len: 300,
            ..RegroupParams::default()
        };
        let st = regroup_fasta(&fa, &d.join("p.fa"), &d.join("f.map"), &[], &p).unwrap();
        assert_eq!(st.families, 2);
        assert_eq!(st.principal_dropped, 1);
        assert_eq!(read_principal(&d), vec!["rg000001_c0_g1_i1".to_string()]);
        assert_eq!(read_map(&d).len(), 2, "family_map 保持全量");
        let _ = std::fs::remove_dir_all(&d);
    }

    /// principal 票数门: 与代表选择同源的原始票数; 无票家族被剔。
    #[test]
    fn principal_min_supp_gate_uses_raw_votes() {
        let d = tmpdir("psupp");
        let voted = rand_seq(961, 800);
        let silent = rand_seq(963, 900);
        // 3 条 read 落在 voted 中段（150bp ≥ k, 命中 ≥ 1/4 窗口）→ 3 票
        let reads: Vec<Vec<u8>> = (0..3).map(|i| voted[300 + i..450 + i].to_vec()).collect();
        let fq = write_fq(&d, "r.fq", &reads);
        let fa = write_fa(&d, &[("c0_g1_i1", silent), ("c1_g1_i1", voted)]);
        let p = RegroupParams {
            cov: 0.7,
            min_principal_supp: 3,
            ..RegroupParams::default()
        };
        let st =
            regroup_fasta(&fa, &d.join("p.fa"), &d.join("f.map"), &[fq], &p).unwrap();
        assert_eq!(st.reads_voted, 3);
        assert_eq!(st.principal_dropped, 1, "silent 家族 0 票被剔: {st:?}");
        assert_eq!(read_principal(&d).len(), 1);
        assert!(read_principal(&d)[0].ends_with("_c1_g1_i1"));
        assert_eq!(read_map(&d).len(), 2, "family_map 保持全量");
        let _ = std::fs::remove_dir_all(&d);
    }

    /// 组合门 AND 语义: 长度门与票数门须同时达标。
    #[test]
    fn principal_combined_gates_are_conjunctive() {
        let d = tmpdir("pboth");
        let voted_short = rand_seq(971, 250); // 有票但 < 300
        let silent_long = rand_seq(973, 900); // ≥ 300 但 0 票
        let voted_long = rand_seq(975, 400); // 双过 → 唯一保留
        let reads: Vec<Vec<u8>> = (0..3)
            .map(|i| voted_short[50 + i..200 + i].to_vec())
            .chain((0..3).map(|i| voted_long[100 + i..250 + i].to_vec()))
            .collect();
        let fq = write_fq(&d, "r.fq", &reads);
        let fa = write_fa(
            &d,
            &[
                ("c0_g1_i1", silent_long),
                ("c1_g1_i1", voted_short),
                ("c2_g1_i1", voted_long),
            ],
        );
        let p = RegroupParams {
            cov: 0.7,
            min_principal_len: 300,
            min_principal_supp: 3,
            ..RegroupParams::default()
        };
        let st =
            regroup_fasta(&fa, &d.join("p.fa"), &d.join("f.map"), &[fq], &p).unwrap();
        assert_eq!(st.families, 3);
        assert_eq!(st.principal_dropped, 2, "长度/票数任一不达标即剔: {st:?}");
        let heads = read_principal(&d);
        assert_eq!(heads.len(), 1);
        assert!(heads[0].ends_with("_c2_g1_i1"), "双过者保留: {heads:?}");
        let _ = std::fs::remove_dir_all(&d);
    }

    /// 等价锁: 门控默认 0/0 → principal 数 == 家族数, dropped == 0
    /// （默认链与旧行为一致）。
    #[test]
    fn principal_gates_default_zero_keep_all() {
        let long = rand_seq(981, 900);
        let short = [long[..700].to_vec(), rand_seq(983, 100)].concat();
        let d = tmpdir("pdef");
        let st = run(&d, &[("c0_g1_i1", long), ("c0_g1_i2", short)], 0.7);
        assert_eq!(st.principal_dropped, 0);
        assert_eq!(read_principal(&d).len(), st.families);
        let _ = std::fs::remove_dir_all(&d);
    }
}
