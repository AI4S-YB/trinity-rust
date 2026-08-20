//! c0 组件全链路径搜索对拍：穿线 → SuffStats → POG/zipping（T7）→
//! reorganizeReadPairings → 分量 → triplet → getAllProbablePaths → fasta。
//!
//! 黄金：`fixtures/p4/c0/allprobpaths.orig.fasta`（Butterfly.jar 完整输出）：
//! ```text
//! java -jar Butterfly.jar -N 4342 -L 200 -F 10000 -R 2 -C c0.graph
//! ```
//! jar 输出含 T9 后处理（assignCompatibleReadsToPaths / cd-hit 类去冗余 /
//! EM reduce / 基因分组）；本任务的管线止于 remove_identical_subseqs +
//! remove_short_seqs 直出。因此比较策略：
//! 1. **序列多重集**必须一致（核心目标）；
//! 2. 数量可能因 T9 后处理（未分配 read 的路径删除等）有差异——本测试断言
//!    我们的序列集 ⊆ 黄金序列集（T9 只删不增），T9 完成后回归为全等。

use std::fs;
use std::path::Path;

use trinity_butterfly::{run_component, ComponentParams};

fn repo_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
}

/// fasta → (header, seq) 列表。
fn parse_fasta(text: &str) -> Vec<(String, String)> {
    let mut res = Vec::new();
    let mut header = String::new();
    let mut seq = String::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if !header.is_empty() {
                res.push((std::mem::take(&mut header), std::mem::take(&mut seq)));
            }
            header = h.to_string();
        } else {
            seq.push_str(line.trim());
        }
    }
    if !header.is_empty() {
        res.push((header, seq));
    }
    res
}

/// c0 全链：走库层 run_component（CLI/测试同入口，防漂移）。
/// 黄金命令：`java -jar Butterfly.jar -N 4342 -L 200 -F 10000 -R 2 -C c0.graph`
/// （jar 默认 EM 形态 = no_em_reduce:false；中间检查点对拍见
/// threading_c0 / pair_stats_c0 / pog_c0）。
fn run_c0_chain() -> Vec<(String, String)> {
    let graph_text =
        fs::read_to_string(repo_root().join("fixtures/p3/quantify/c0/orig.graph.out")).unwrap();
    let reads_text =
        fs::read_to_string(repo_root().join("fixtures/p3/quantify/c0/orig.reads.out")).unwrap();
    let result = run_component(
        &graph_text,
        &reads_text,
        &ComponentParams {
            n: 4342,
            name: "c0".to_string(),
            ..ComponentParams::default()
        },
    )
    .unwrap();
    parse_fasta(&result.all_prob_paths_fasta)
}

/// T9 全量对拍：header + 序列全行多重集与黄金 allprobpaths.orig.fasta
/// 完全一致（含 `_g1_i{1..5}` 命名 / len= / path=[...] MISO 坐标 +
/// 尾部完整路径 List 串）。
#[test]
fn c0_all_prob_paths_full_lines_match_java() {
    let ours = run_c0_chain();
    let golden_text =
        fs::read_to_string(repo_root().join("fixtures/p4/c0/allprobpaths.orig.fasta")).unwrap();
    let golden = parse_fasta(&golden_text);
    assert_eq!(golden.len(), 5, "c0 黄金应有 5 条转录本");

    // ---- 全行（header + 序列）多重集相等 ----
    let norm =
        |recs: &Vec<(String, String)>| -> std::collections::BTreeMap<(String, String), usize> {
            let mut m = std::collections::BTreeMap::new();
            for (h, s) in recs {
                *m.entry((h.trim().to_string(), s.clone())).or_default() += 1;
            }
            m
        };
    let (o, g) = (norm(&ours), norm(&golden));
    assert_eq!(o, g, "c0 全行多重集必须与黄金完全一致（T9 闭合）");
}

/// c0 输出按 Java HashMap 桶序逐条与黄金顺序一致（i1..i5 的 len 序：
/// 1034/1287/1078/1103/914——非长度序，来自 HashMap 迭代序；EM 路径下
/// 最终 map 为 putAll 构建，表容量 tableSizeFor(5/0.75+1)=8，i1/i2 同桶）。
#[test]
fn c0_isoform_order_and_naming_match_java() {
    let ours = run_c0_chain();
    let golden_text =
        fs::read_to_string(repo_root().join("fixtures/p4/c0/allprobpaths.orig.fasta")).unwrap();
    let golden = parse_fasta(&golden_text);
    assert_eq!(ours.len(), golden.len());
    // 逐条（顺序敏感）：header 与序列都应一致——顺序即 _g1_i1..i5 编号序
    for (i, ((oh, os), (gh, gs))) in ours.iter().zip(golden.iter()).enumerate() {
        assert_eq!(
            oh,
            gh,
            "第 {i} 条 header 应逐字段一致（含 _g1_i{} 编号）",
            i + 1
        );
        assert_eq!(os, gs, "第 {i} 条序列应一致");
        assert!(
            gh.contains(&format!("_g1_i{} ", i + 1)) || gh.contains(&format!("_g1_i{}]", i + 1))
        );
    }
}

/// 1136 变体删除路径复现：T8 直出的 6 条中，len 1136 变体
/// （[2341 → 22bp 变体节点 → 1599]）经 reduce_cdhit_like 因与
/// len 1287 的 [2341, 4651, 1599] 过于相似（低支持）被删。
/// 间接验证：最终集合不含 1136 长度序列。
#[test]
fn c0_len1136_variant_removed() {
    let ours = run_c0_chain();
    assert_eq!(ours.len(), 5);
    assert!(
        !ours.iter().any(|(_, s)| s.len() == 1136),
        "len 1136 变体应被 cd-hit 式去冗余删除"
    );
}
