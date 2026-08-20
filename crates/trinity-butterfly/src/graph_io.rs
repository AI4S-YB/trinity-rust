//! graph 文件 IO：镜像 `TransAssembly_allProbPaths.java` 的
//! preProcessGraphFile (L12694) / buildNewGraphUseKmers (L12742) /
//! getReadStarts+readAndMapSingleRead 解析部分 (L11770) / writeDotFile (L12837)。

use std::fs;
use std::path::Path;

use rustc_hash::FxHashMap;

use crate::context::BflyContext;
use crate::graph::{DiGraph, SeqVertex, SimpleEdge, T_VERTEX_ID, VERTEX_ROOT_ID};
use crate::{BflyError, Result};

/// Java `Math.round(double)` 语义：floor(x + 0.5)。
pub fn java_math_round(x: f64) -> i64 {
    (x + 0.5).floor() as i64
}

// ---------------------------------------------------------------------------
// preProcessGraphFile
// ---------------------------------------------------------------------------

/// 统计每个节点的 in/out 流量并记录 kmer（首行组件头丢弃）。
/// 每行 5 列：`to from supp kmer flag`。
pub fn pre_process_graph_file(
    text: &str,
) -> (
    FxHashMap<i32, f64>,
    FxHashMap<i32, f64>,
    FxHashMap<i32, String>,
) {
    let mut out_flow: FxHashMap<i32, f64> = FxHashMap::default();
    let mut in_flow: FxHashMap<i32, f64> = FxHashMap::default();
    let mut kmers: FxHashMap<i32, String> = FxHashMap::default();

    for (i, line) in text.lines().enumerate() {
        if i == 0 || line.is_empty() {
            continue; // header of component
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            continue;
        }
        let to: i32 = fields[0]
            .parse()
            .unwrap_or_else(|_| panic!("bad to: {line}"));
        let from: i32 = fields[1]
            .parse()
            .unwrap_or_else(|_| panic!("bad from: {line}"));
        let supp: f64 = fields[2]
            .parse()
            .unwrap_or_else(|_| panic!("bad supp: {line}"));
        *out_flow.entry(from).or_insert(0.0) += supp;
        *in_flow.entry(to).or_insert(0.0) += supp;
        kmers.insert(to, fields[3].to_string());
    }
    (in_flow, out_flow, kmers)
}

// ---------------------------------------------------------------------------
// buildNewGraphUseKmers
// ---------------------------------------------------------------------------

/// 建图产物。
pub struct BuildResult {
    pub graph: DiGraph,
    pub ctx: BflyContext,
    /// rootIDs：from<0 的行对应的 to 节点（每碱基权重 = supp）。
    pub root_ids: Vec<i32>,
}

/// `getSeqVertex` 语义：nodeTracker 命中且节点仍在图中。
fn get_seq_vertex(graph: &DiGraph, tracker: &FxHashMap<i32, i32>, id: i32) -> bool {
    tracker.contains_key(&id) && graph.contains_vertex(id)
}

/// 从 graph.out 文本建 de Bruijn 图（保留每个 kmer 首字母的版本）。
///
/// - supp < 0 的行跳过（INITIAL_EDGE_ABS_THR = 0）
/// - LAST_ID = 所有行 max(from, to)
/// - KMER_SIZE 由首个 kmer 推断；不一致则 Err
/// - root（from < 0）：toV = SeqVertex(to, kmer, supp)（每碱基权重 supp），无入边
/// - 非 root：fromV/toV 为普通 SeqVertex（无权重），加 SimpleEdge(supp)
pub fn build_new_graph_use_kmers(text: &str) -> Result<BuildResult> {
    let (_, _, kmers) = pre_process_graph_file(text);

    let mut graph = DiGraph::new();
    let mut ctx = BflyContext::new();
    let mut root_ids: Vec<i32> = Vec::new();

    for (i, line) in text.lines().enumerate() {
        if i == 0 || line.is_empty() {
            continue; // header of component
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            return Err(BflyError::GraphFileFormat(format!(
                "line {}: {:?}",
                i + 1,
                line
            )));
        }
        let to: i32 = fields[0].parse().map_err(|_| {
            BflyError::GraphFileFormat(format!("bad to at line {}: {}", i + 1, line))
        })?;
        let from: i32 = fields[1].parse().map_err(|_| {
            BflyError::GraphFileFormat(format!("bad from at line {}: {}", i + 1, line))
        })?;
        let supp: f64 = fields[2].parse().map_err(|_| {
            BflyError::GraphFileFormat(format!("bad supp at line {}: {}", i + 1, line))
        })?;
        if supp < 0.0 {
            continue; // INITIAL_EDGE_ABS_THR
        }

        if from > ctx.last_id {
            ctx.last_id = from;
        }
        if to > ctx.last_id {
            ctx.last_id = to;
        }

        let kmer = fields[3];
        if ctx.kmer_size == 0 {
            ctx.kmer_size = kmer.len();
        } else if ctx.kmer_size != kmer.len() {
            return Err(BflyError::GraphFileFormat(format!(
                "Error, discrepancy among kmer lengths.  Stored: {}, found: {}\n{}",
                ctx.kmer_size,
                kmer.len(),
                line
            )));
        }

        // fromV = getSeqVertex(graph, from)
        let mut from_v_exists = get_seq_vertex(&graph, &ctx.node_tracker, from);
        if !from_v_exists && from >= 0 {
            let name = kmers.get(&from).cloned().unwrap_or_default();
            let v = SeqVertex::new(from, name);
            ctx.node_tracker.insert(from, 1);
            graph.add_vertex(v);
            from_v_exists = true;
        }

        let is_root = from < 0 || !from_v_exists;

        // important to look up toV after possibly creating fromV (fromV == toV bugfix)
        let to_v_exists = get_seq_vertex(&graph, &ctx.node_tracker, to);

        if is_root {
            if !to_v_exists {
                // root 节点：每碱基权重 supp
                let v = SeqVertex::with_per_base_weight(to, kmer, supp);
                ctx.node_tracker.insert(to, 1);
                graph.add_vertex(v);
                root_ids.push(to);
            }
        } else {
            if !to_v_exists {
                let v = SeqVertex::new(to, kmer);
                ctx.node_tracker.insert(to, 1);
                graph.add_vertex(v);
            }
            let e = SimpleEdge::new(supp, from, to);
            graph.add_edge(from, to, e);
        }
    }

    // Java main L739：图建完后 LAST_REAL_ID = LAST_ID（此后新建的 id 均为"合成"节点）
    ctx.last_real_id = ctx.last_id;

    Ok(BuildResult {
        graph,
        ctx,
        root_ids,
    })
}

// ---------------------------------------------------------------------------
// reads 解析（getReadStarts / readAndMapSingleRead 的字段解析部分；穿线在后续任务）
// ---------------------------------------------------------------------------

/// read 解析错误。
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum ReadParseError {
    #[error("read 行字段不足（需 ≥7，含空 field5）: {0}")]
    TooFewFields(String),
    #[error("read 行数字字段非法: {0}")]
    BadNumber(String),
}

/// graph.reads 行的原始解析结果（未穿线）。
#[derive(Debug, Clone, PartialEq)]
pub struct RawRead {
    pub name: String,
    pub seq: String,
    pub start_in_read: i64,
    /// fields[3] + KMER_SIZE（FIXME off-by-one，镜像 Java）。
    pub end_in_read: i64,
    pub from_orig_v: i32,
}

fn strip_read_name(name: &str, treat_pairs_as_single: bool) -> String {
    let mut name = name.strip_prefix('>').unwrap_or(name);
    if !treat_pairs_as_single {
        for suffix in ["/1", "/2", "\\1", "\\2", ":1", ":2"] {
            if let Some(stripped) = name.strip_suffix(suffix) {
                name = stripped;
                break;
            }
        }
    }
    name.to_string()
}

/// 解析 graph.reads（首行组件头丢弃）。每行 8 列：
/// `>name start fromOrigV end endFull <empty> seq strand`（seq = fields[6]）。
///
/// **只解析不穿线**（穿线 findPathInGraph 属后续任务）。
pub fn parse_graph_reads(
    text: &str,
    kmer_size: usize,
    treat_pairs_as_single: bool,
) -> std::result::Result<Vec<RawRead>, ReadParseError> {
    let mut reads = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if i == 0 || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 7 {
            return Err(ReadParseError::TooFewFields(line.to_string()));
        }
        let name = strip_read_name(fields[0], treat_pairs_as_single);
        let start_in_read: i64 = fields[1]
            .parse()
            .map_err(|_| ReadParseError::BadNumber(line.to_string()))?;
        let from_orig_v: i32 = fields[2]
            .parse()
            .map_err(|_| ReadParseError::BadNumber(line.to_string()))?;
        let end_field: i64 = fields[3]
            .parse()
            .map_err(|_| ReadParseError::BadNumber(line.to_string()))?;
        // Java FIXME: endInRead = fields[3] + KMER_SIZE（chrysalis 恒定 off-by-one）
        let end_in_read = end_field + kmer_size as i64;
        let seq = fields[6].to_string();
        reads.push(RawRead {
            name,
            seq,
            start_in_read,
            end_in_read,
            from_orig_v,
        });
    }
    Ok(reads)
}

// ---------------------------------------------------------------------------
// writeDotFile
// ---------------------------------------------------------------------------

/// 生成 DOT 文本（镜像 writeDotFile L12837-12927）。
///
/// - 头 `digraph G {`
/// - 节点行 `{id} [label="{shortSeqWID}[L:{adjLen}][T:{discTime}]"]`，
///   print_full_seq 时用 `{longSeqWID}[L:{adjLen}]`
/// - weightAvg > 25 → ` ,style=bold,color="#AF0000"`
/// - 边行 `{u}->{v}[label={w}]`；w > 20 → `[style=bold,label={w},color="#AF0000"]`
/// - ROOT(-1)/T(-2) 的节点行与边行均不输出
pub fn write_dot_string(graph: &DiGraph, kmer_size: usize, print_full_seq: bool) -> String {
    let mut out = String::new();
    out.push_str("digraph G {\n");

    for &id in graph.vertex_ids() {
        let vertex = graph.get_vertex(id).expect("vertex in vertex_ids");
        let has_preds = graph.in_degree(id) > 0;
        let adj = vertex.get_name_kmer_adj(kmer_size, has_preds);

        let label = if print_full_seq {
            format!(
                "{}[L:{}]",
                vertex.get_long_seq_wid(kmer_size, has_preds),
                adj.len()
            )
        } else {
            format!(
                "{}[L:{}][T:{}]",
                vertex.get_short_seq_wid(kmer_size, has_preds),
                adj.len(),
                vertex.dfs_discovery_time
            )
        };

        let mut ver_desc = format!("{} [label=\"{}\"", id, label);
        if vertex.get_weight_avg() > 25 {
            ver_desc.push_str(" ,style=bold,color=\"#AF0000\"");
        }
        ver_desc.push(']');

        if id != T_VERTEX_ID && id != VERTEX_ROOT_ID {
            out.push_str(&ver_desc);
            out.push('\n');
        }

        for &succ in graph.get_successors(id) {
            let Some(edge) = graph.find_edge(id, succ) else {
                continue;
            };
            let weight = java_math_round(edge.weight);
            let edge_style = if weight > 20 {
                format!("[style=bold,label={weight},color=\"#AF0000\"]")
            } else {
                format!("[label={weight}]")
            };
            if succ != T_VERTEX_ID && id != VERTEX_ROOT_ID {
                out.push_str(&format!("{id}->{succ}{edge_style}\n"));
            }
        }
    }

    out.push_str("}\n");
    out
}

/// 写 DOT 文件到 path。
pub fn write_dot_file(
    graph: &DiGraph,
    kmer_size: usize,
    path: &Path,
    print_full_seq: bool,
) -> std::io::Result<()> {
    fs::write(path, write_dot_string(graph, kmer_size, print_full_seq))
}

#[cfg(test)]
mod tests {
    use super::*;

    const SMALL: &str = "Component 0\n\
        2341\t-1\t0\tTTTTCTTGCAATACACAAAAGTTT\t1\n\
        2340\t2341\t1\tTTTCTTGCAATACACAAAAGTTTA\t1\n\
        2339\t2340\t2\tTTCTTGCAATACACAAAAGTTTAT\t1\n";

    #[test]
    fn pre_process_counts_flows_and_kmers() {
        let (in_flow, out_flow, kmers) = pre_process_graph_file(SMALL);
        assert_eq!(out_flow[&-1], 0.0);
        assert_eq!(out_flow[&2341], 1.0);
        assert_eq!(out_flow[&2340], 2.0);
        assert_eq!(in_flow[&2341], 0.0);
        assert_eq!(in_flow[&2340], 1.0);
        assert_eq!(in_flow[&2339], 2.0);
        assert_eq!(kmers[&2341], "TTTTCTTGCAATACACAAAAGTTT");
        assert!(!kmers.contains_key(&-1));
    }

    #[test]
    fn build_small_graph_semantics() {
        let br = build_new_graph_use_kmers(SMALL).unwrap();
        let g = &br.graph;
        // KMER_SIZE = 24（首 kmer 长度）
        assert_eq!(br.ctx.kmer_size, 24);
        // LAST_ID = max(from,to) = 2341
        assert_eq!(br.ctx.last_id, 2341);
        // 节点：2341(root), 2340, 2339；root 无入边
        assert_eq!(g.vertex_count(), 3);
        assert_eq!(br.root_ids, vec![2341]);
        assert_eq!(g.in_degree(2341), 0);
        // root 权重：每碱基 supp=0 → weightAvg = round(0)=0，weights len 24
        let root = g.get_vertex(2341).unwrap();
        assert_eq!(root.weights.len(), 24);
        assert!(root.weights.iter().all(|&w| w == 0.0));
        // 非 root 无权重 → weightAvg = -1
        assert_eq!(g.get_vertex(2340).unwrap().get_weight_avg(), -1);
        // 边权
        assert_eq!(g.find_edge(2341, 2340).unwrap().weight, 1.0);
        assert_eq!(g.find_edge(2340, 2339).unwrap().weight, 2.0);
        assert_eq!(g.edge_count(), 2);
        // node_tracker 覆盖所有建过的节点
        for id in [2341, 2340, 2339] {
            assert!(br.ctx.node_tracker.contains_key(&id));
        }
    }

    #[test]
    fn build_root_weight_per_base() {
        let text = "Component 0\n5\t-1\t7\tACGTACGTAC\t1\n";
        let br = build_new_graph_use_kmers(text).unwrap();
        let root = br.graph.get_vertex(5).unwrap();
        assert_eq!(root.weights, vec![7.0; 10]);
        assert_eq!(root.get_weight_avg(), 7);
    }

    #[test]
    fn build_skips_negative_supp_and_kmer_mismatch_err() {
        let text = "Component 0\n5\t-1\t7\tACGTACGTAC\t1\n6\t5\t-1\tCGTACGTACA\t1\n";
        let br = build_new_graph_use_kmers(text).unwrap();
        assert_eq!(br.graph.edge_count(), 0);
        assert_eq!(br.ctx.last_id, 5);

        let mismatch = "Component 0\n5\t-1\t7\tACGTACGTAC\t1\n6\t5\t3\tCGTACGTACAT\t1\n";
        assert!(build_new_graph_use_kmers(mismatch).is_err());
    }

    #[test]
    fn dot_output_small() {
        let br = build_new_graph_use_kmers(SMALL).unwrap();
        let dot = write_dot_string(&br.graph, br.ctx.kmer_size, false);
        let lines: Vec<&str> = dot.lines().collect();
        assert_eq!(lines[0], "digraph G {");
        assert_eq!(*lines.last().unwrap(), "}");
        // root 节点：L:24（全名，无前驱）
        assert!(dot.contains("2341 [label=\"TTTTCTTGCAATACACAAAAGTTT:W0(V2341_D-1)[L:24][T:-1]\"]"));
        // 非 root：kmerAdj 单字母，W-1
        assert!(dot.contains("2340 [label=\"A:W-1(V2340_D-1)[L:1][T:-1]\"]"));
        assert!(dot.contains("2339 [label=\"T:W-1(V2339_D-1)[L:1][T:-1]\"]"));
        assert!(dot.contains("2341->2340[label=1]"));
        assert!(dot.contains("2340->2339[label=2]"));
        // print_full_seq 变体
        let dot_full = write_dot_string(&br.graph, br.ctx.kmer_size, true);
        assert!(dot_full.contains("2340 [label=\"A(V2340)[L:1]\"]"));
    }

    #[test]
    fn dot_bold_thresholds() {
        let mut g = DiGraph::new();
        g.add_vertex(SeqVertex::with_per_base_weight(1, "ACGT", 30.0)); // avg 30 > 25
        g.add_vertex(SeqVertex::new(2, "CGTA"));
        g.add_edge(1, 2, SimpleEdge::new(21.0, 1, 2)); // round 21 > 20
        g.add_vertex(SeqVertex::new(3, "GTAC"));
        g.add_edge(1, 3, SimpleEdge::new(20.4, 1, 3)); // round 20 → 普通
        let dot = write_dot_string(&g, 4, false);
        assert!(
            dot.contains("1 [label=\"ACGT:W30(V1_D-1)[L:4][T:-1]\" ,style=bold,color=\"#AF0000\"]")
        );
        assert!(dot.contains("1->2[style=bold,label=21,color=\"#AF0000\"]"));
        assert!(dot.contains("1->3[label=20]"));
    }

    #[test]
    fn dot_root_and_t_suppressed() {
        let mut g = DiGraph::new();
        g.add_vertex(SeqVertex::new(VERTEX_ROOT_ID, "S"));
        g.add_vertex(SeqVertex::with_per_base_weight(1, "ACGT", 5.0));
        g.add_vertex(SeqVertex::new(T_VERTEX_ID, "E"));
        g.add_edge(VERTEX_ROOT_ID, 1, SimpleEdge::new(5.0, VERTEX_ROOT_ID, 1));
        g.add_edge(1, T_VERTEX_ID, SimpleEdge::new(5.0, 1, T_VERTEX_ID));
        let dot = write_dot_string(&g, 4, false);
        assert_eq!(dot.lines().count(), 3); // 头 + 节点1 + }
                                            // 节点 1 有前驱（ROOT），故只输出 kmerAdj 尾部
        assert!(dot.contains("1 [label=\"T:W5(V1_D-1)[L:1][T:-1]\"]"));
        // ROOT/T 的节点行与边行都被抑制（剩余行只能是头/节点1/}）
        for l in dot.lines() {
            assert!(
                !l.starts_with("-1 ") && !l.starts_with("-2 "),
                "unexpected line: {l}"
            );
            assert!(
                !l.contains("->-1")
                    && !l.contains("->-2")
                    && !l.starts_with("-1->")
                    && !l.starts_with("-2->")
            );
        }
    }

    #[test]
    fn parse_reads_off_by_one_and_name_suffix() {
        let text = "Component 0\n\
            >readA/2\t11\t101393\t36\t101418\t\tGAAAGACTGTCACCCTTGAGGTGGAGTCCTCTGAC\t-\n\
            >readB:1\t0\t5\t10\t20\t\tACGTACGTACGTACGT\t-\n";
        let reads = parse_graph_reads(text, 25, false).unwrap();
        assert_eq!(reads.len(), 2);
        let a = &reads[0];
        assert_eq!(a.name, "readA"); // /2 去掉
        assert_eq!(a.start_in_read, 11);
        assert_eq!(a.end_in_read, 36 + 25); // off-by-one: f3 + K
        assert_eq!(a.from_orig_v, 101393);
        assert_eq!(a.seq, "GAAAGACTGTCACCCTTGAGGTGGAGTCCTCTGAC");
        let b = &reads[1];
        assert_eq!(b.name, "readB"); // :1 去掉
        assert_eq!(b.end_in_read, 10 + 25);
        // TREAT_PAIRS_AS_SINGLE = true → 后缀保留
        let reads = parse_graph_reads(text, 25, true).unwrap();
        assert_eq!(reads[0].name, "readA/2");
        assert_eq!(reads[1].name, "readB:1");
    }

    #[test]
    fn parse_reads_errors() {
        let bad = "Component 0\n>short\t1\t2\t3\n";
        assert!(matches!(
            parse_graph_reads(bad, 25, false),
            Err(ReadParseError::TooFewFields(_))
        ));
    }
}
