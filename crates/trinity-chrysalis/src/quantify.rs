//! QuantifyGraph（Chrysalis/analysis/QuantifyGraph.cc）移植——P3-T6。
//!
//! 镜像语义：
//! - **读索引**（`add_kmers`，:225-235）：reads FASTA 以 shortName 读入
//!   （name = 首 token 含 '>'，序列整体大写），每条 read j 的每个 25-mer
//!   位置 i 入 `KmerEntry{read, pos}`——**含 N 的 k-mer 照入，无熵过滤**；
//!   排序键 = `KmerEntryCompare`（:186-223）对 read 原始字节从 pos 起逐
//!   字节比较 25 字符（跨 read 的字典序）。
//! - **第一遍扫图**（:367-380）：`first[node] = 第 4 列 kmer 串首字符`。
//! - **第二遍**（:386-473）：token <4 的行原样透传两个输出并冲刷上一组件
//!   的 SortPrint；否则 prevNode>=0 时拼 25-mer `sub = first[prevNode]+24mer`
//!   做二分收集（`BasesToNumberCountPlus`，:237-287），非 strand 模式下
//!   revcomp 再查一次并对新条目做 ori=-1 坐标翻转；第 3 列改写 n1+n2。
//!
//! **熵怪癖（有意复刻，勿"修复"）**：原版
//! `strncpy(kmerseq, &d[1], kmer_length)` 对 25 长的 sub 从下标 1 起拷 25
//! 字节——越界一字节。strncpy 语义（源不足 n 补 '\0'）下结果等价于
//! `sub[1..25]` 的 24 字符 + 1 个 NUL，共 25 字节进 `compute_entropy`：
//! **NUL 不进分子但计入分母**，故 24 个计数字符、分母 25。低复杂度
//! （`< 1.0`，QuantifyGraph.cc:28 的本地常量，**非** GraphFromFasta 的 1.3）
//! 的边直接跳过（n1/n2 保持原值）。
//!
//! **SortPrint 怪癖**（:45-146）：ids 按 (ori, id, start) 升序（ori=-1 组
//! 在前，KmerTable.h:94-109）；按 (id, ori) 分组，仅当 `lastStart > 组首
//! start`（**严格 >**，单 kmer 位置的 read 丢弃）时输出；输出行在
//! node2 与 seq 之间是**双 tab**——line 以 `name\tpos1\tnode1\t` 起拼，
//! 冲刷时追加 `pos2\tnode2\t` 后 `fprintf("%s\t")` 又补一个 tab
//! （`name \t pos1 \t node1 \t pos2 \t node2 \t \t seq \t +|-`），已与原版
//! 二进制产物逐字节对拍确认。

use std::collections::HashMap;

use rayon::prelude::*;
use trinity_common::error::CommonError;

use crate::dna_vector::{compute_entropy, read_fasta_bytes, revcomp, DnaSeq};

/// QuantifyGraph.cc:28 `static float MIN_KMER_ENTROPY = 1.0;`（本地于本文件，
/// 与 GraphFromFasta 的 1.3 不同）。
pub const QG_MIN_KMER_ENTROPY: f32 = 1.0;

/// CLI `-k`（默认 24），内部 k-mer 长 = k+1 = 25；`-strand`；`-max_reads`
/// （原版默认 -1 = 不限，Trinity 流水线实参 200000；>0 截取前 N 条 read）。
#[derive(Debug, Clone)]
pub struct QgParams {
    pub k: usize,
    pub strand: bool,
    pub max_reads: i64,
}

impl Default for QgParams {
    fn default() -> Self {
        QgParams {
            k: 24,
            strand: false,
            max_reads: 200_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QuantifyOutput {
    pub graph_out: String,
    pub reads_out: String,
}

/// QuantifyGraph.cc:163-184 `ReadsExt`：只在**最后 6 个字符**内找 '.'
/// （原版 `if (n - i > 6) break;` 先于判 '.'），从该 '.' 截断换 `.reads`；
/// 无则原样追加 `.reads`。
pub fn reads_ext_filename(path: &str) -> String {
    let b = path.as_bytes();
    let start = b.len().saturating_sub(6);
    for i in (start..b.len()).rev() {
        if b[i] == b'.' {
            let mut out = path[..i].to_string();
            out.push_str(".reads");
            return out;
        }
    }
    format!("{path}.reads")
}

/// **性能重构（2026-09-04 压时间战役, 产物逐字节不变）**: 原实现
/// `add_kmers` 对全部 25-mer 做间接寻址比较排序（`kmer_bytes` 双层间接 +
/// 25 字节 memcmp, O(N log N) 次 cache miss）, quantify 阶段 39→311s 的
/// 主热点。改为**两群物化键**:
/// - 主群: 全 ACGT 的 k-mer 按保序 2bit 映射（A=0,C=1,G=2,T=3 —— 与 ASCII
///   序 A<C<G<T 一致）打包 u64（kk ≤ 32 时恒可装下 2·32=64 bit）;
/// - 次群: 其余（含 N/其他字符）按原字节比较路径, 天然稀少。
///
/// **等价性**: 查询是精确匹配——sub 全 ACGT 时只可能与主群条目相等（次群
/// 至少一个字节非 ACGT）, 反之亦然, 两群不交叉, 无跨群序问题。同 kmer 组内
/// 收集序 = (read, pos) 升序: 原版稳定排序保持生成序（read 主序内 pos 升序）,
/// 新实现 `sort_unstable_by_key((key, read, pos))` 显式保证同键组内
/// (read, pos) 升序 —— 与生成序一致, ids push 顺序不变, 输出逐字节相同。
#[derive(Debug, Default)]
struct KmerTable {
    /// 主群（kk ≤ 32 且 k-mer 全 ACGT; kk > 32 时恒空, 全走次群）。
    main: Vec<MainEntry>,
    /// 次群（key = 原始字节, 稀少: Q30 数据 N 率 ~0.01%）。
    alt: Vec<(Vec<u8>, u32, u32)>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct MainEntry {
    key: u64,
    read: u32,
    pos: u32,
}

/// 保序 2bit 打包: 全 ACGT 的定长字节串 → u64（None = 含其他字符或 kk > 32）。
/// A=0,C=1,G=2,T=3 与 ASCII 序一致 → u64 字典序 == 字节字典序。
#[inline]
fn pack_acgt(sub: &[u8]) -> Option<u64> {
    if sub.len() > 32 {
        return None;
    }
    let mut k = 0u64;
    for &c in sub {
        let v = match c {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        k = (k << 2) | v as u64;
    }
    Some(k)
}

/// KmerTable.h:63-135 `IDS`：ori=-1 表示 revcomp 命中（坐标已翻转）。
/// 字段紧凑化（压时间战役 2026-09-04, 值域论证）: id ≤ max_reads 截断
/// （u32 充裕）; start < read 长（i32）; edge = 组件内节点号（组件图行数
/// 远达不到 i32 域）; ori ∈ {1,-1}（i8）。40B → 16B, 亿级命中条目的
/// 排序比较与内存流量同步降。排序仍为稳定 sort_by_key((ori,id,start))
/// ——tie（同 (ori,id,start) 不同 edge）保持收集序, 与原版一致。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Ids {
    id: u32,
    start: i32,
    edge: i32,
    ori: i8,
}

/// QuantifyGraph 主流程（两遍扫图）。输入 `graph_tmp` / `reads_tmp`（分区器
/// 产物，reads_tmp 为 `>name pct\nseq\n` 的 FASTA）。输出名镜像 Trinity:2352
/// 实参形态：`-o cN.graph.out` → `graph_out = 去 .tmp 后缀 + ".graph.out"`，
/// `reads_out = reads_ext_filename(graph_out)`（→ `cN.graph.reads`）。
/// 返回两条输出路径。不删除输入（原版默认 unlink 输入，库层不复刻）。
pub fn quantify_graph(
    graph_tmp: &str,
    reads_tmp: &str,
    p: &QgParams,
) -> Result<QuantifyOutput, CommonError> {
    let graph = std::fs::read_to_string(graph_tmp)?;
    let reads = std::fs::read(reads_tmp)?;
    let (g, r) = quantify_graph_mem(&graph, &reads, p)?;
    let graph_out = match graph_tmp.strip_suffix(".tmp") {
        Some(stem) => format!("{stem}.graph.out"),
        None => format!("{graph_tmp}.graph.out"),
    };
    let reads_out = reads_ext_filename(&graph_out);
    std::fs::write(&graph_out, &g)?;
    std::fs::write(&reads_out, &r)?;
    Ok(QuantifyOutput {
        graph_out,
        reads_out,
    })
}

/// 内存版 QuantifyGraph（IO 优化: 管线层直写打包文件, 跳过 .graph.out /
/// .graph.reads 落盘）。输入输出与 [`quantify_graph`] 逐字节一致——
/// `read_fasta_short_names(path)` ≡ `read_fasta_bytes(bytes, true)`。
/// 返回 (graph_out 文本, reads 文本)。
pub fn quantify_graph_mem(
    graph_text: &str,
    reads_data: &[u8],
    p: &QgParams,
) -> Result<(String, String), CommonError> {
    let kk = p.k + 1; // 内部 25
    let t_parse = std::time::Instant::now();
    let mut seqs = read_fasta_bytes(reads_data, true);
    if p.max_reads > 0 {
        seqs.truncate(p.max_reads as usize);
    }

    let t_add = std::time::Instant::now();
    let kmers = add_kmers(&seqs, kk);
    let t_join = std::time::Instant::now();
    if qg_profile_enabled() {
        QG_STATS.parse_nanos.fetch_add(
            t_parse.elapsed().as_nanos() as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
        QG_STATS.add_nanos.fetch_add(
            t_add.elapsed().as_nanos() as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
    }

    let graph = graph_text;

    let mut first: HashMap<i64, u8> = HashMap::new();
    for line in graph.lines() {
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() >= 4 {
            let node: i64 = toks[0].parse().unwrap_or(0);
            first.insert(node, toks[3].as_bytes()[0]);
        }
    }

    let (mut g, mut r) = (String::new(), String::new());
    let mut ids: Vec<Ids> = Vec::new();
    for line in graph.lines() {
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() < 4 {
            g.push_str(line);
            g.push('\n');
            if !ids.is_empty() {
                sort_print(&mut r, &mut ids, &seqs);
            }
            r.push_str(line);
            r.push('\n');
            ids.clear();
            continue;
        }
        let prev_node: i64 = toks[1].parse().unwrap_or(0);
        let mut n1: i64 = 0;
        let mut n2: i64 = 0;
        if prev_node >= 0 {
            let kmer = toks[3].as_bytes();
            // sub = first[prevNode] + kmer（25 字节；原版 DNAVector 下标 0 起拼）
            let lead = match first.get(&prev_node) {
                Some(&c) => c,
                None => {
                    println!("ERROR!! first[prevNode] where prevNode = {prev_node} unset");
                    b'N'
                }
            };
            let mut sub = Vec::with_capacity(kk);
            sub.push(lead);
            sub.extend_from_slice(kmer);

            bases_to_number_count_plus(&kmers, &seqs, &mut ids, &mut n1, &sub, prev_node, kk);

            if !p.strand {
                let rc = revcomp(&sub);
                let from = ids.len();
                bases_to_number_count_plus(&kmers, &seqs, &mut ids, &mut n2, &rc, prev_node, kk);
                if n1 + n2 < i32::MAX as i64 {
                    for e in &mut ids[from..] {
                        e.ori = -1;
                        let len = seqs[e.id as usize].seq.len() as i32;
                        let pos = e.start + 1;
                        // 原版 SetStart(len - pos - k + 1) = len - start - 25
                        e.start = len - pos - kk as i32 + 1;
                    }
                } else {
                    println!("WARNING: k-mer overflow, n={} . Discarding.", n1 + n2);
                    // 原版只清零计数，已入 ids 的条目保留（ori 仍 1）——照抄
                    n1 = 0;
                    n2 = 0;
                }
            }
        }

        for (i, t) in toks.iter().enumerate() {
            if i > 0 {
                g.push('\t');
            }
            if i == 2 {
                g.push_str(&format!("{}", (n1 + n2) as i32));
            } else {
                g.push_str(t);
            }
        }
        g.push('\n');
    }
    if !ids.is_empty() {
        sort_print(&mut r, &mut ids, &seqs);
    }
    if qg_profile_enabled() {
        QG_STATS.join_nanos.fetch_add(
            t_join.elapsed().as_nanos() as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
    }

    Ok((g, r))
}

/// QG_PROFILE=1 时的分段累计（压时间战役 2026-09-04, 不对外承诺格式）:
/// parse = read_fasta_bytes + 截断, add = add_kmers（生成+排序）,
/// join = 其余（两遍扫图 + 二分收集 + sort_print）。
pub struct QgStats {
    pub parse_nanos: std::sync::atomic::AtomicU64,
    pub add_nanos: std::sync::atomic::AtomicU64,
    pub join_nanos: std::sync::atomic::AtomicU64,
}

pub static QG_STATS: QgStats = QgStats {
    parse_nanos: std::sync::atomic::AtomicU64::new(0),
    add_nanos: std::sync::atomic::AtomicU64::new(0),
    join_nanos: std::sync::atomic::AtomicU64::new(0),
};

pub fn qg_profile_enabled() -> bool {
    static P: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *P.get_or_init(|| std::env::var_os("QG_PROFILE").is_some())
}

/// `add_kmers`（:225-235）: 全部 k-mer 入两群表（见 [`KmerTable`]）,
/// 各群按 (键, read, pos) 排序——同 kmer 组内 (read, pos) 升序 = 生成序,
/// 与原单一数组 + 稳定排序的收集序逐条目一致。
fn add_kmers(seqs: &[DnaSeq], kk: usize) -> KmerTable {
    let mut t = KmerTable::default();
    for (j, s) in seqs.iter().enumerate() {
        if s.seq.len() >= kk {
            for i in 0..=s.seq.len() - kk {
                let sub = &s.seq[i..i + kk];
                match pack_acgt(sub) {
                    Some(key) => t.main.push(MainEntry {
                        key,
                        read: j as u32,
                        pos: i as u32,
                    }),
                    None => t.alt.push((sub.to_vec(), j as u32, i as u32)),
                }
            }
        }
    }
    // par_sort: 巨型组件（200k read ≈ 25M k-mer）排序是大头; 组件间已在
    // rayon 池并行, 长尾大组件经 work-stealing 吃到空闲核。(key, read, pos)
    // 为全序无重复键 → 并行排序结果确定, 与串行 sort_unstable 逐字节一致。
    t.main.par_sort_unstable_by_key(|e| (e.key, e.read, e.pos));
    t.alt
        .sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    t
}

/// `BasesToNumberCountPlus`（:237-287）：
/// 1. 熵怪癖：对 `sub[1..25]`（24 字符）补一个 NUL 成 25 字节算熵，
///    `< 1.0` 直接返回（count 不动）；
/// 2. lower_bound（键 = sub 完整 25 字符）起线性收集相等条目入 ids，
///    count = 命中数；查不到则 count 不动（原版 return -1 前不赋值）。
///    主/次群按 [`pack_acgt`] 归位（精确匹配下两群不交叉）。
fn bases_to_number_count_plus(
    table: &KmerTable,
    _seqs: &[DnaSeq],
    ids: &mut Vec<Ids>,
    count: &mut i64,
    sub: &[u8],
    edge: i64,
    kk: usize,
) {
    // strncpy(&d[1], kmer_length) 越界一字节的等价化：24 字符 + NUL，分母 25
    let mut ent = sub[1..kk].to_vec();
    ent.push(0);
    if compute_entropy(&ent) < QG_MIN_KMER_ENTROPY {
        return;
    }

    let mut n = 0i64;
    if let Some(sk) = pack_acgt(sub) {
        let m = &table.main;
        let mut i = m.partition_point(|e| e.key < sk);
        while i < m.len() && m[i].key == sk {
            ids.push(Ids {
                id: m[i].read,
                start: m[i].pos as i32,
                edge: edge as i32,
                ori: 1,
            });
            n += 1;
            i += 1;
        }
    } else {
        let a = &table.alt;
        let mut i = a.partition_point(|e| e.0.as_slice() < sub);
        while i < a.len() && a[i].0.as_slice() == sub {
            ids.push(Ids {
                id: a[i].1,
                start: a[i].2 as i32,
                edge: edge as i32,
                ori: 1,
            });
            n += 1;
            i += 1;
        }
    }
    // 命中（原版此路径 n ≥ 1）才覆盖 count; 未命中不动（原版 return -1 语义）
    if n > 0 {
        *count = n;
    }
}

/// `SortPrint`（:45-146）。输出直接追加到 `out`（graph/read 两个文件共享的
/// 文本流由调用方拆分——见 quantify_graph）。
fn sort_print(out: &mut String, ids: &mut [Ids], seqs: &[DnaSeq]) {
    // KmerTable.h IDS::operator< ：(ori, id, start) 升序，ori=-1 组在前
    ids.sort_by_key(|e| (e.ori, e.id, e.start));

    let mut last_id: i64 = -1;
    let mut id: i64 = -1;
    let mut start: i64 = -1;
    let mut edge: i64 = -1;
    let mut last_start: i64 = -1;
    let mut last_edge: i64 = -1;
    let mut ori: i32 = 1;
    let mut line = String::new();
    let mut last_start_temp: i64 = -1;
    let mut last_ori: i32 = 1;

    let flush = |out: &mut String, line: &str, cond: bool, seq: &[u8], rev: bool, last_ori: i32| {
        if !cond {
            return;
        }
        out.push_str(line);
        out.push('\t'); // fprintf("%s\t")：line 已以 '\t' 结尾 → 双 tab
        let d = if rev { revcomp(seq) } else { seq.to_vec() };
        out.push_str(&String::from_utf8_lossy(&d));
        out.push_str(if last_ori == -1 { "\t-\n" } else { "\t+\n" });
    };

    for e in ids.iter() {
        id = e.id as i64;
        ori = e.ori as i32;
        start = e.start as i64;
        edge = e.edge as i64;
        if id != last_id || ori != last_ori {
            if last_id != -1 {
                line.push_str(&format!("{last_start}\t{last_edge}\t"));
                flush(
                    out,
                    &line,
                    last_start > last_start_temp,
                    &seqs[last_id as usize].seq,
                    last_ori == -1,
                    last_ori,
                );
            }
            line = format!("{}\t{start}\t{edge}\t", seqs[e.id as usize].name);
            last_start_temp = start;
        }
        last_id = id;
        last_start = start;
        last_edge = edge;
        last_ori = ori;
    }

    if id != -1 {
        line.push_str(&format!("{start}\t{edge}\t"));
        flush(
            out,
            &line,
            last_start > last_start_temp,
            &seqs[id as usize].seq,
            ori == -1,
            last_ori,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(name: &str, s: &str) -> DnaSeq {
        DnaSeq {
            name: name.into(),
            seq: s.as_bytes().to_vec(),
        }
    }

    #[test]
    fn reads_ext_variants() {
        // '.' 在最后 6 字符内 → 截断换 .reads
        assert_eq!(reads_ext_filename("c0.graph.out"), "c0.graph.reads");
        assert_eq!(reads_ext_filename("qg_out_0.graph"), "qg_out_0.reads");
        // 最后 6 字符内无 '.' → 追加
        assert_eq!(reads_ext_filename("graph_out"), "graph_out.reads");
        // 恰在第 7 个位置（越界一侧）的 '.' 不算："abcdefg.h" 的 '.' 距尾 2
        assert_eq!(reads_ext_filename("abcdefg.h"), "abcdefg.reads");
        // 距尾 3 的 '.'（"abcdef.gh"）在最后 6 字符内 → 截断
        assert_eq!(reads_ext_filename("abcdef.gh"), "abcdef.reads");
        // 距尾 7 的 '.' 检查不到："abcdefg.hij"（'.' 是第 8 个字符）
        assert_eq!(reads_ext_filename("abcdef.ghijkl"), "abcdef.ghijkl.reads");
    }

    #[test]
    fn add_kmers_two_group_ordering_equivalent_to_lexicographic() {
        // 原语义: 全部 k-mer 按原始字节字典序稳定排序（N 的 ASCII 在 G/T 之间）。
        // 新实现: 主群（全 ACGT, 2bit u64 保序）+ 次群（含 N 按字节）,
        // 各群内 (键, read, pos) 升序。两群精确匹配不交叉 → 收集序等价。
        let seqs = vec![
            seq(">a", "ACGTN"), // k=4（kk=5）: ACGTN、CGTN…（ACGTN 含 N → 次群）
            seq(">b", "ACGTA"),
            seq(">c", "ACGTG"), // 主群序: ACGTA < ACGTG; 次群: ACGTN
        ];
        let t = add_kmers(&seqs, 5);
        // 主群键序 == 字节序（A<C<G<T 保序打包）
        let main_keys: Vec<u64> = t.main.iter().map(|e| e.key).collect();
        let mut sorted = main_keys.clone();
        sorted.sort_unstable();
        assert_eq!(main_keys, sorted, "主群须按 key 升序");
        assert!(t.main.iter().all(|e| e.pos < 5));
        // 次群字节序
        let alt_keys: Vec<&[u8]> = t.alt.iter().map(|e| e.0.as_slice()).collect();
        let mut sorted = alt_keys.clone();
        sorted.sort();
        assert_eq!(alt_keys, sorted, "次群须按字节升序");
        // 含 N 的 5-mer 全在次群（ACGTN、CGTN…）
        assert!(t.alt.iter().all(|e| e.0.contains(&b'N')));
        assert!(!t.main.is_empty());

        // 同 kmer 组内 (read, pos) 升序 = 生成序（原稳定排序的收集序）
        let seqs = vec![seq(">x", "ACGTAACGTA"), seq(">y", "ACGTAACGTA")];
        let t = add_kmers(&seqs, 5);
        let hits: Vec<(u32, u32)> = t
            .main
            .iter()
            .filter(|e| e.key == pack_acgt(b"ACGTA").unwrap())
            .map(|e| (e.read, e.pos))
            .collect();
        assert_eq!(
            hits,
            vec![(0, 0), (0, 5), (1, 0), (1, 5)],
            "同 kmer 组内 (read,pos) 须为生成序"
        );
    }

    // 手推拼查询 + ori=-1 坐标换算的端到端小向量：
    // first[10]='G'（节点 10 的 kmer 以 G 开头），边 11←10 的 kmer 列
    // "ACGTT…"，sub = "GACGTT…"。read 长度 30，正链命中位置 2；revcomp 命中
    // 位置 5 → 翻转后 start = len - start - 25 = 30-5-25 = 0。
    fn write_inputs(dir: &std::path::Path) {
        let reads = ">r1 100%\nTTGACGTTACGTTACGTTACGTTACGTACGTA\n";
        std::fs::write(dir.join("c0.reads.tmp"), reads).unwrap();
        let graph = "Component 0\n10\t-1\t1\tGCGTTACGTTACGTTACGTTACGT\n11\t10\t7\tACGTTACGTTACGTTACGTTACGT\n";
        std::fs::write(dir.join("c0.graph.tmp"), graph).unwrap();
    }

    #[test]
    fn quantify_join_flip_and_output() {
        let dir = std::env::temp_dir().join("p3_t6_unit_qg");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        write_inputs(&dir);
        let out = quantify_graph(
            dir.join("c0.graph.tmp").to_str().unwrap(),
            dir.join("c0.reads.tmp").to_str().unwrap(),
            &QgParams {
                k: 24,
                strand: false,
                max_reads: -1,
            },
        )
        .unwrap();
        let g = std::fs::read_to_string(&out.graph_out).unwrap();
        // sub = "G"+"ACGTT…"*4+"ACG"？"GACGTTACGTTACGTTACGTTACG"（低复杂度？
        // 熵 = ACGT 各 ~6/25 → 远大于 1，不跳过）。read 内 "…TT|GACGTT…|ACGTA"
        // 命中 1 次；revcomp 不含 → 计数 1。
        let lines: Vec<&str> = g.lines().collect();
        assert_eq!(lines[0], "Component 0");
        assert_eq!(lines[1], "10\t-1\t0\tGCGTTACGTTACGTTACGTTACGT");
        let f: Vec<&str> = lines[2].split('\t').collect();
        assert_eq!(f[0], "11");
        assert_eq!(f[1], "10");
        assert_eq!(f[3], "ACGTTACGTTACGTTACGTTACGT");
        assert_eq!(f[2], "1", "唯一正链命中（revcomp 不在 read 中）");
        let r = std::fs::read_to_string(&out.reads_out).unwrap();
        // 只有 1 个命中位置 → 单 kmer 位置丢弃（严格 >），无记录行
        let recs: Vec<&str> = r.lines().filter(|l| l.starts_with('>')).collect();
        assert!(recs.is_empty(), "单一命中位置的 read 必须被丢弃：{r}");
    }

    #[test]
    fn sort_print_grouping_and_double_tab() {
        // 直接构造 ids：read 0（len 60）在 edge 7 有两个 ori=1 位置 3、9，
        // 两个 ori=-1 位置（翻转后 10、20），read 1 单位置（必须丢弃）
        let seqs = vec![seq(">r0", &"AC".repeat(30)), seq(">r1", &"GT".repeat(30))];
        let mut ids = vec![
            Ids {
                id: 0,
                start: 3,
                edge: 7,
                ori: 1,
            },
            Ids {
                id: 0,
                start: 9,
                edge: 7,
                ori: 1,
            },
            Ids {
                id: 1,
                start: 5,
                edge: 8,
                ori: 1,
            }, // 单位置 → 丢弃
            Ids {
                id: 0,
                start: 20,
                edge: 7,
                ori: -1,
            },
            Ids {
                id: 0,
                start: 10,
                edge: 7,
                ori: -1,
            },
        ];
        let mut out = String::new();
        sort_print(&mut out, &mut ids, &seqs);
        let lines: Vec<&str> = out.lines().collect();
        // ori=-1 组在前；revcomp 输出 + '-'；pos1=组首(10)、pos2=组末(20)
        assert_eq!(lines.len(), 2);
        let f: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(f[0], ">r0");
        assert_eq!(f[1], "10");
        assert_eq!(f[2], "7");
        assert_eq!(f[3], "20");
        assert_eq!(f[4], "7");
        assert_eq!(f[5], ""); // 双 tab 怪癖
        assert_eq!(f[6], "GT".repeat(30)); // revcomp(AC*30)
        assert_eq!(f[7], "-");
        let f2: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!((f2[1], f2[3], f2[7]), ("3", "9", "+"));
        assert_eq!(f2[6], "AC".repeat(30));
    }

    #[test]
    fn entropy_quirk_skips_low_complexity_edge() {
        // 24 字符 + NUL、分母 25 的怪癖：poly-T 边熵 < 1.0 被跳过，
        // 即使 read 中存在完整 25-mer 命中
        let dir = std::env::temp_dir().join("p3_t6_unit_ent");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(
            dir.join("c0.reads.tmp"),
            ">r1 100%\nGTTTTTTTTTTTTTTTTTTTTTTTTAAAAA\n",
        )
        .unwrap();
        let graph = format!(
            "Component 0\n10\t-1\t1\tG{}\n11\t10\t7\t{}\n",
            "T".repeat(23),
            "T".repeat(24)
        );
        std::fs::write(dir.join("c0.graph.tmp"), graph).unwrap();
        let out = quantify_graph(
            dir.join("c0.graph.tmp").to_str().unwrap(),
            dir.join("c0.reads.tmp").to_str().unwrap(),
            &QgParams {
                k: 24,
                strand: false,
                max_reads: -1,
            },
        )
        .unwrap();
        let g = std::fs::read_to_string(&out.graph_out).unwrap();
        // sub = "G"+24T；熵输入 = 本行 kmer 列（24T）+ NUL（分母 25）= 0 → 跳过，
        // 即使 read "G"+24T+"AAAAA" 含完整 25-mer 命中，计数仍 0
        assert!(g.contains(&format!("11\t10\t0\t{}", "T".repeat(24))), "{g}");
    }

    #[test]
    fn ori_flip_coordinate_math() {
        // len=50，start=7 → pos=8 → new = 50-8-25+1 = 18
        let seqs = [seq(">x", &"A".repeat(50))];
        let mut ids = [Ids {
            id: 0,
            start: 7,
            edge: 1,
            ori: 1,
        }];
        for e in ids.iter_mut() {
            e.ori = -1;
            let len = seqs[e.id as usize].seq.len() as i32;
            let pos = e.start + 1;
            e.start = len - pos - 25 + 1;
        }
        assert_eq!(ids[0].start, 18);
    }
}
