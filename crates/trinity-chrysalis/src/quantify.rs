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

use trinity_common::error::CommonError;

use crate::dna_vector::{compute_entropy, read_fasta_short_names, revcomp, DnaSeq};

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

/// KmerTable.h:31-59 `KmerEntry`（只留 Index/Pos）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct KmerEntry {
    read: usize,
    pos: usize,
}

/// KmerTable.h:63-135 `IDS`：ori=-1 表示 revcomp 命中（坐标已翻转）。
#[derive(Debug, Clone, PartialEq, Eq)]
struct Ids {
    id: usize,
    start: i64,
    edge: i64,
    ori: i32,
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
    let kk = p.k + 1; // 内部 25
    let mut seqs = read_fasta_short_names(reads_tmp)?;
    if p.max_reads > 0 {
        seqs.truncate(p.max_reads as usize);
    }

    let kmers = add_kmers(&seqs, kk);

    let graph = std::fs::read_to_string(graph_tmp)?;
    let graph_out = match graph_tmp.strip_suffix(".tmp") {
        Some(stem) => format!("{stem}.graph.out"),
        None => format!("{graph_tmp}.graph.out"),
    };
    let reads_out = reads_ext_filename(&graph_out);

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
                        let len = seqs[e.id].seq.len() as i64;
                        let pos = e.start + 1;
                        // 原版 SetStart(len - pos - k + 1) = len - start - 25
                        e.start = len - pos - kk as i64 + 1;
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

    std::fs::write(&graph_out, &g)?;
    std::fs::write(&reads_out, &r)?;
    Ok(QuantifyOutput {
        graph_out,
        reads_out,
    })
}

/// `add_kmers`（:225-235）：全部 25-mer 入表后按 25 字节字典序排序。
fn add_kmers(seqs: &[DnaSeq], kk: usize) -> Vec<KmerEntry> {
    let mut v: Vec<KmerEntry> = Vec::new();
    for (j, s) in seqs.iter().enumerate() {
        if s.seq.len() >= kk {
            for i in 0..=s.seq.len() - kk {
                v.push(KmerEntry { read: j, pos: i });
            }
        }
    }
    v.sort_by(|a, b| kmer_bytes(seqs, a, kk).cmp(kmer_bytes(seqs, b, kk)));
    v
}

#[inline]
fn kmer_bytes<'a>(seqs: &'a [DnaSeq], e: &KmerEntry, kk: usize) -> &'a [u8] {
    &seqs[e.read].seq[e.pos..e.pos + kk]
}

/// `BasesToNumberCountPlus`（:237-287）：
/// 1. 熵怪癖：对 `sub[1..25]`（24 字符）补一个 NUL 成 25 字节算熵，
///    `< 1.0` 直接返回（count 不动）；
/// 2. lower_bound（键 = sub 完整 25 字符）起线性收集相等条目入 ids，
///    count = 命中数；查不到则 count 不动（原版 return -1 前不赋值）。
fn bases_to_number_count_plus(
    kmers: &[KmerEntry],
    seqs: &[DnaSeq],
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

    let idx = kmers.partition_point(|e| kmer_bytes(seqs, e, kk) < sub);
    if idx == kmers.len() || sub < kmer_bytes(seqs, &kmers[idx], kk) {
        return; // 未命中
    }
    let mut n = 0i64;
    for e in &kmers[idx..] {
        if kmer_bytes(seqs, e, kk) != sub {
            break;
        }
        ids.push(Ids {
            id: e.read,
            start: e.pos as i64,
            edge,
            ori: 1,
        });
        n += 1;
    }
    *count = n;
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
        ori = e.ori;
        start = e.start;
        edge = e.edge;
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
            line = format!("{}\t{start}\t{edge}\t", seqs[e.id].name);
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
    fn add_kmers_sorts_by_25_byte_lexicographic_including_n() {
        // N 的 ASCII ('N'=0x4E) > 'GATC'，混 N 的 k-mer 按原始字节序排
        let seqs = vec![
            seq(">a", "TTACG"),
            seq(">b", "AAAN"), // 短于 25 的 read 不入表
        ];
        let kmers = add_kmers(&seqs, 25);
        // "TTACG" 25-mer 不存在（len 5 < 25）→ 空
        assert!(kmers.is_empty());

        let seqs = vec![
            seq(">a", "ACGTN"), // k=4（kk=5 测试）：ACGTN、CGTN…
            seq(">b", "ACGTA"),
        ];
        let kmers = add_kmers(&seqs, 5);
        let keys: Vec<&[u8]> = kmers.iter().map(|e| kmer_bytes(&seqs, e, 5)).collect();
        assert_eq!(
            keys,
            vec![
                &b"ACGTA"[..], // read b（'A' < 'N'）
                &b"ACGTN"[..], // read a
            ]
        );
        assert_eq!(kmers[0].read, 1);
        assert_eq!(kmers[1].read, 0);
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
            let len = seqs[e.id].seq.len() as i64;
            let pos = e.start + 1;
            e.start = len - pos - 25 + 1;
        }
        assert_eq!(ids[0].start, 18);
    }
}
