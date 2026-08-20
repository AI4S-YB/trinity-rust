//! P3-T3 对拍——`fixtures/p3` 真实链（gff.welds.orig.txt → sort →
//! BubbleUpClustering → CreateIwormFastaBundle）与原版二进制产物比对。
//!
//! 比较契约：
//! - **COMPONENT 块多重集**：按成员集合分组的块（块内成员序 + #POOL_INFO +
//!   序列折行）多重集相等——component 编号非契约（原版池遍历序受
//!   weld 图排序稳定性影响，本库 `sort_weld_graph` 为稳定排序）；
//! - **bundle 逐字节**：`>s_<no>` 编号来自 COMPONENT 块编号，原版产物
//!   确定性 → 逐字节相等意味着上一步的编号在本链上实际也一致。

use std::collections::BTreeMap;
use std::path::Path;

use trinity_chrysalis::bubble_up::{bubble_up_clustering, BubbleParams};
use trinity_chrysalis::bundle::create_iworm_fasta_bundle;
use trinity_chrysalis::dna_vector::read_fasta;
use trinity_chrysalis::graph_from_fasta::sort_weld_graph;

fn fixture(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/p3")
        .join(name)
}

/// 把 COMPONENT 输出切成块：键 = 成员 iworm 下标序列（含重复），
/// 值 = 整块原文（去掉 component 号本身——`COMPONENT <id>` 行、
/// `>Component_<id>` 前缀、`#POOL_INFO <id>` 中的编号均替换为占位）。
fn component_blocks(text: &str) -> BTreeMap<String, Vec<String>> {
    let mut blocks: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut cur_members: Vec<String> = Vec::new();
    let mut cur_body: Vec<String> = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split_whitespace().collect();
        if !f.is_empty() && f[0] == "COMPONENT" {
            cur_members.clear();
            cur_body.clear();
        } else if !f.is_empty() && f[0].starts_with(">Component") {
            // >Component_<id> <size> <z> [iworm>name]
            cur_members.push(f[2].to_string());
            cur_body.push(format!("MEMBER {} {}", f[1], f[3]));
        } else if !f.is_empty() && f[0] == "#POOL_INFO" {
            cur_body.push(format!(
                "#POOL_INFO {}",
                f.get(2..).map(|s| s.join(" ")).unwrap_or_default()
            ));
        } else if !f.is_empty() && f[0] == "END" {
            let key = cur_members.join(",");
            blocks.entry(key).or_default().push(cur_body.join("\n"));
        } else {
            cur_body.push(line.to_string()); // 序列折行
        }
    }
    blocks
}

#[test]
fn bubble_up_matches_original_block_multiset() {
    let iworm = read_fasta(fixture("gff.iworm.fa")).unwrap();
    let sorted = sort_weld_graph(&std::fs::read_to_string(fixture("gff.welds.orig.txt")).unwrap());
    let params = BubbleParams {
        min_contig_length: 200,
        max_cluster_size: 25,
        debug_weld_all: false,
    };
    let ours = bubble_up_clustering(&iworm, &sorted, &params).unwrap();
    let orig = std::fs::read_to_string(fixture("bubble.orig.out")).unwrap();
    let a = component_blocks(&ours);
    let b = component_blocks(&orig);
    // 每个成员集合的块数（多重集）与块体逐一比较
    let keys: std::collections::BTreeSet<String> = a.keys().chain(b.keys()).cloned().collect();
    for k in keys {
        let va = a.get(&k).map(|v| v.as_slice()).unwrap_or(&[]);
        let vb = b.get(&k).map(|v| v.as_slice()).unwrap_or(&[]);
        assert_eq!(va.len(), vb.len(), "成员集合 [{k}] 块数不一致");
        for (x, y) in va.iter().zip(vb.iter()) {
            assert_eq!(x, y, "成员集合 [{k}] 块体不一致");
        }
    }
}

#[test]
fn bundle_matches_original_byte_for_byte() {
    let bubble = {
        let iworm = read_fasta(fixture("gff.iworm.fa")).unwrap();
        let sorted =
            sort_weld_graph(&std::fs::read_to_string(fixture("gff.welds.orig.txt")).unwrap());
        let params = BubbleParams {
            min_contig_length: 200,
            max_cluster_size: 25,
            debug_weld_all: false,
        };
        bubble_up_clustering(&iworm, &sorted, &params).unwrap()
    };
    let ours = create_iworm_fasta_bundle(&bubble, 200).unwrap();
    let orig = std::fs::read_to_string(fixture("bundle.orig.fa")).unwrap();
    assert_eq!(ours, orig);
}
