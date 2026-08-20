//! P3-T6 对拍——`fixtures/p3` 真实链（f2db.orig.txt + rtt.orig.out →
//! sort_reads_to_components → partition → QuantifyGraph）与原版产物
//! （perl 分区脚本 + `Chrysalis/bin/QuantifyGraph -k 24 -no_cleanup`，
//! 黄金向量固化于 `fixtures/p3/quantify/c{0,1,2}/`）比对。
//!
//! 比较契约：graph.out（cN.graph.out 对 qg_out_N.graph，二者同为 QuantifyGraph
//! 的 -o 图输出）**逐行相等**（含第 3 列计数与透传头行）；graph.reads
//! **按行集合相等**（SortPrint 行序在 ori 分组内确定，跨 (id,ori) 分组序
//! 由 (ori,id,start) 全序确定——逐行比较，见测试内断言）。

use std::path::{Path, PathBuf};

use trinity_chrysalis::partition::{partition, PartParams};
use trinity_chrysalis::quantify::{quantify_graph, QgParams};
use trinity_chrysalis::reads_to_transcripts::sort_reads_to_components;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/p3")
        .join(name)
}

fn read(p: impl AsRef<Path>) -> String {
    std::fs::read_to_string(p).unwrap()
}

/// 原版 `strncpy(&d[1], kmer_length)` 越界一字节读 std::vector 堆尾巴——
/// 该字节是**上一任同尺寸缓冲的陈旧核苷酸字符**，使 `compute_entropy` 的
/// 第 25 个字符在 {NUL, A, C, G, T} 之间浮动：熵恰在 1.0 阈值附近的
/// "刀口边"（poly-A/GA 低复杂度）跳过与否随堆历史翻转（对拍 c0 中 8/2885
/// 行、c1/c2 类似）。本库取确定语义（NUL 补位：24 计数字符 / 分母 25），
/// 对拍契约：**非刀口边逐行相等；刀口边仅允许第 3 列计数不同**；
/// reads 输出过滤提及刀口边的行后逐行相等。
#[test]
fn quantify_and_partition_match_original() {
    let f2db = read(fixture("f2db.orig.txt"));
    let sorted = sort_reads_to_components(&read(fixture("rtt.orig.out")));

    let out_dir = std::env::temp_dir().join("p3_t6_e2e_partition");
    let _ = std::fs::remove_dir_all(&out_dir);
    let listing = partition(&f2db, &sorted, &out_dir, &PartParams::default()).unwrap();
    assert_eq!(listing.len(), 55, "55 组件通过长度过滤（与 perl 脚本一致）");

    for comp in [0u64, 1, 2] {
        let base = out_dir.join(format!("Cbin0/c{comp}"));
        let want_dir = fixture(&format!("quantify/c{comp}"));
        // 分区产物与黄金逐字节相等
        assert_eq!(
            read(base.with_extension("graph.tmp")),
            read(want_dir.join("graph.tmp")),
            "c{comp} graph.tmp"
        );
        assert_eq!(
            read(base.with_extension("reads.tmp")),
            read(want_dir.join("reads.tmp")),
            "c{comp} reads.tmp"
        );

        // 刀口边集合：prev_node，其 24-mer 熵在补 X ∈ {A,C,G,T} 时跨越 1.0
        let knife = knife_edges(&read(want_dir.join("graph.tmp")));

        // QuantifyGraph：输出写到临时目录再比较
        let tmp = std::env::temp_dir().join(format!("p3_t6_qg_c{comp}"));
        let _ = std::fs::remove_dir_all(&tmp);
        std::fs::create_dir_all(&tmp).unwrap();
        let g = tmp.join(format!("c{comp}.graph.tmp"));
        let r = tmp.join(format!("c{comp}.reads.tmp"));
        std::fs::copy(base.with_extension("graph.tmp"), &g).unwrap();
        std::fs::copy(base.with_extension("reads.tmp"), &r).unwrap();
        let out = quantify_graph(
            g.to_str().unwrap(),
            r.to_str().unwrap(),
            &QgParams {
                k: 24,
                strand: false,
                max_reads: -1,
            },
        )
        .unwrap();
        assert!(out.graph_out.ends_with(".graph.out"));
        assert!(out.reads_out.ends_with(".graph.reads"));

        // graph 输出：非刀口边逐行相等；刀口边仅第 3 列（计数）可差
        let got_s = read(&out.graph_out);
        let want_s = read(want_dir.join("orig.graph.out"));
        let got: Vec<&str> = got_s.lines().collect();
        let want: Vec<&str> = want_s.lines().collect();
        assert_eq!(got.len(), want.len(), "c{comp} graph.out 行数");
        let mut knife_diffs = 0usize;
        for (i, (a, b)) in got.iter().zip(&want).enumerate() {
            if a == b {
                continue;
            }
            let fa: Vec<&str> = a.split('\t').collect();
            let fb: Vec<&str> = b.split('\t').collect();
            assert_eq!(fa.len(), fb.len(), "c{comp} graph.out:{} 列数", i + 1);
            let prev: i64 = fa[1].parse().unwrap();
            assert!(
                knife.contains(&prev),
                "c{comp} graph.out 第 {} 行差异非刀口边：{a:?} vs {b:?}",
                i + 1
            );
            for c in 0..fa.len() {
                if c != 2 {
                    assert_eq!(fa[c], fb[c], "c{comp} graph.out:{} 第 {c} 列", i + 1);
                }
            }
            knife_diffs += 1;
        }

        // reads 输出：过滤提及刀口边的行（node1/node2 列）后逐行相等
        let got_rs = read(&out.reads_out);
        let want_rs = read(want_dir.join("orig.reads.out"));
        let gr: Vec<&str> = got_rs.lines().collect();
        let wr: Vec<&str> = want_rs.lines().collect();
        assert_eq!(gr.len(), wr.len(), "c{comp} reads 行数");
        let mut skipped = 0usize;
        for (i, (a, b)) in gr.iter().zip(&wr).enumerate() {
            let mentions_knife = |l: &str| {
                l.split('\t').enumerate().any(|(c, f)| {
                    (c == 2 || c == 4)
                        && f.parse::<i64>()
                            .map(|n| knife.contains(&n))
                            .unwrap_or(false)
                })
            };
            if mentions_knife(a) || mentions_knife(b) {
                skipped += 1;
                continue;
            }
            assert_eq!(a, b, "c{comp} reads 第 {} 行", i + 1);
        }
        eprintln!(
            "c{comp}: graph 刀口边差 {knife_diffs} 行，reads 过滤 {skipped} 行（总 {}）",
            gr.len()
        );
    }
}

/// 刀口边：两次熵判定输入（正向 = 本行 kmer 列；反向 = revcomp(first[prev]
/// +kmer) 去首字符）任一随越界字节 X ∈ {NUL,A,C,G,T} 跨越 1.0 阈值。
fn knife_edges(graph: &str) -> std::collections::HashSet<i64> {
    use trinity_chrysalis::dna_vector::{compute_entropy, revcomp};
    let mut first = std::collections::HashMap::new();
    for line in graph.lines() {
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.len() >= 4 {
            let node: i64 = t[0].parse().unwrap_or(0);
            first.insert(node, t[3].as_bytes()[0]);
        }
    }
    let mut set = std::collections::HashSet::new();
    for line in graph.lines() {
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.len() < 4 {
            continue;
        }
        let prev: i64 = t[1].parse().unwrap_or(0);
        if prev < 0 {
            continue;
        }
        let kmer = t[3].as_bytes();
        let lead = first.get(&prev).copied().unwrap_or(b'N');
        let mut sub = Vec::with_capacity(25);
        sub.push(lead);
        sub.extend_from_slice(kmer);
        let rc_input = revcomp(&sub)[1..].to_vec();
        let skip = |k: &[u8], x: u8| {
            let mut s = k.to_vec();
            s.push(x);
            compute_entropy(&s) < 1.0
        };
        let flips = |k: &[u8]| {
            let base = skip(k, 0);
            [b'A', b'C', b'G', b'T'].iter().any(|&x| skip(k, x) != base)
        };
        if flips(kmer) || flips(&rc_input) {
            set.insert(prev);
        }
    }
    set
}
