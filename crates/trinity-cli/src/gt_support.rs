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

use std::fs;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use rustc_hash::{FxHashMap, FxHashSet};

use trinity_butterfly::{parse_component_graph, EdgeSupport};
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

/// 跑整个支持 pass: 建 k2v → 扫 reads → 落盘。返回统计。
pub fn run_support_pass(
    listing: &[(u64, PathBuf)],
    read_files: &[PathBuf],
    seq_type: SeqType,
    out_bin: &Path,
) -> Result<SupportStats, CommonError> {
    // ---- 1. 逐 comp 读 graph.out 记录, 建 k2v 与边集 ----
    let (k2v, edge_sets, k) = build_k2v(listing)?;
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

/// 从 graphs.pack 读全部 comp 的 graph.out → (k2v, 边集, K)。
/// K = 首个 kmer 的长度（原型恒 24; 图内 kmer 长度一致由 butterfly 解析保证）。
fn build_k2v(listing: &[(u64, PathBuf)]) -> Result<(K2V, EdgeSets, usize), CommonError> {
    let mut k2v: K2V = FxHashMap::default();
    let mut edge_sets: EdgeSets = FxHashMap::default();
    let mut k_len: Option<usize> = None;
    // listing 按 id 升序 → 同 Cbin 连续; 每 bin 的索引与句柄各开一次
    let mut cur_bin: Option<(PathBuf, PackIndex, fs::File)> = None;
    for (id, base) in listing {
        let dir = base.parent().map(|d| d.to_path_buf()).unwrap_or_default();
        if cur_bin
            .as_ref()
            .map(|(d, _, _)| d != &dir)
            .unwrap_or(true)
        {
            let idx = PackIndex::load(&dir.join(format!("{GRAPHS_PACK}.idx")))?;
            let file = fs::File::open(dir.join(format!("{GRAPHS_PACK}.pack")))?;
            cur_bin = Some((dir.clone(), idx, file));
        }
        let (_, idx, file) = cur_bin.as_ref().unwrap();
        let (offset, len) = idx.find(*id, KIND_GRAPH_OUT).ok_or_else(|| {
            CommonError::Parse(format!(
                "graphs pack {} 缺记录 (id={id}, kind=GRAPH_OUT)",
                dir.display()
            ))
        })?;
        let bytes = read_at(file, offset, len)?;
        let text = String::from_utf8(bytes)
            .map_err(|e| CommonError::Parse(format!("comp {id} graph.out 非 UTF-8: {e}")))?;
        let g = parse_component_graph(&text)
            .map_err(|e| CommonError::Parse(format!("comp {id} graph.out 解析失败: {e}")))?;
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
    }
    Ok((k2v, edge_sets, k_len.unwrap_or(24)))
}

/// 单序列窗口扫描（python `process` 语义, 含"不含最后一窗"quirk）。
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
    let mut prev: Option<(u64, i32)> = None;
    // python: for i in range(len(seq) - K)——最后一个合法窗 [len-K, len) 不查
    for i in 0..(seq.len() - k) {
        let hit = kmer_to_intval(&seq[i..i + k])
            .ok()
            .and_then(|val| k2v.get(&(val, k32)));
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
fn scan_reads_file(
    path: &Path,
    seq_type: SeqType,
    k: usize,
    k2v: &K2V,
    edge_sets: &EdgeSets,
) -> Result<(EdgeSupport, SupportStats), CommonError> {
    let mut sup = EdgeSupport::new();
    let mut stats = SupportStats::default();
    let mut reader = open_maybe_gz(path)?;
    match seq_type {
        SeqType::Fq => {
            // python calc_edge_support2 的 line_no%4 流式（含多行序列收集与
            // 文件尾 flush; 空行也计数行号, 不跳过）
            let mut line = String::new();
            let mut seq_buf = String::new();
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
                            process_seq(
                                seq_buf.as_bytes(),
                                k,
                                k2v,
                                edge_sets,
                                &mut sup,
                                &mut stats.pairs,
                            );
                            seq_buf.clear();
                        }
                        stats.reads += 1;
                    }
                    1 => seq_buf.push_str(line.trim()),
                    _ => {}
                }
                line_no += 1;
            }
            if !seq_buf.is_empty() {
                process_seq(
                    seq_buf.as_bytes(),
                    k,
                    k2v,
                    edge_sets,
                    &mut sup,
                    &mut stats.pairs,
                );
            }
        }
        SeqType::Fa => {
            let mut fr = FastaReader::new(reader);
            while let Some(rec) = fr.next_record()? {
                stats.reads += 1;
                process_seq(
                    rec.sequence.as_bytes(),
                    k,
                    k2v,
                    edge_sets,
                    &mut sup,
                    &mut stats.pairs,
                );
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
        let (k2v, edge_sets, k) = build_k2v(&listing).unwrap();
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
}
