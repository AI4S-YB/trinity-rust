# M3 Read-Support 加权流分解 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在 PogFlow 的流分解中注入 PairPath read 支持——边权从纯 POG 权重升级为「POG 权重 × read 路径覆盖」，后过滤零支持路径。目标 dev1 R ≥ 9%（当前 8.1%）且 P ≥ 84%（保持优势）。

**Architecture:** POG 图构建后，用 `combined_read_hash` 的 PairPath 数据统计每条边被多少 read 路径覆盖，将其作为额外权重注入 FlowGraph。剥离后对每条产出路径计数 PairPath 兼容数，低于 `min_read_support` 的路径被过滤。经典引擎的 read-reinforcement 语义以简化形式复现。

**Tech Stack:** 既有 trinity-butterfly + trinity-flow；无新依赖。

---

## 核心设计

```
POG 图（现有 seqvertex_graph）
  + combined_read_hash（现有 PairPath → count）
  → ① 边 read 覆盖统计: 遍历所有 PairPath 的 path1/path2,
       对路径中每对相邻节点 (u,v) 在 edge_read_count[(u,v)] += count
  → ② 权重融合: edge_weight = max(pog_weight, edge_read_count)
       （read 覆盖高的边不被 POG 低权截断——治零权/低权边问题）
  → ③ 流分解（现有 peel_component）
  → ④ 路径 read 支持过滤: 对每条产出路径,
       count = Σ PairPath.is_compatible_and_contained_by_single_path(path) 的 count
       count < min_read_support → 过滤
  → ⑤ 输出（现有 dedup + FASTA）
```

**为什么 max 而非 sum**: POG 权重和 read 计数的量纲不同（POG 权重来自路径重叠，read 计数是 read 条数）。取 max 保证两者中较高者生效——read 覆盖高时不被 POG 低权截断，POG 权重高时保持稳定。

## 文件结构

```
crates/trinity-butterfly/src/component.rs  # PogFlow 分支扩展
crates/trinity-butterfly/src/pair_paths.rs  # PairPath::edge_coverage 公共方法
docs/opt-log.md
```

---

### Task 1: PairPath 边覆盖统计 + 路径兼容计数（TDD）

**Files:** Modify `crates/trinity-butterfly/src/pair_paths.rs`

- [ ] **Step 1: 失败测试**

```rust
    #[test]
    fn edge_coverage_counts() {
        // PairPath path1=[1,2,3], path2=[4,5] → 边覆盖: (1,2)+1, (2,3)+1, (4,5)+1
        let mut s = SuffStats::default();
        let pp = PairPath::new(vec![1,2,3]);
        s.combined_read_hash.entry(1).or_default().insert(pp.clone(), 5);
        let cov = edge_read_coverage(&s);
        assert_eq!(cov.get(&(1,2)), Some(&5));
        assert_eq!(cov.get(&(2,3)), Some(&5));
        assert!(!cov.contains_key(&(1,3)));  // 不相邻的边不计数
    }

    #[test]
    fn path_read_support_counts_compatible() {
        // 路径 [1,2,3,4] 包含 PairPath path1=[1,2,3] → support=5
        let mut s = SuffStats::default();
        s.combined_read_hash.entry(1).or_default().insert(PairPath::new(vec![1,2,3]), 5);
        s.combined_read_hash.entry(1).or_default().insert(PairPath::new(vec![9,10]), 3);
        let support = path_read_support(&s, &[1,2,3,4]);
        assert_eq!(support, 5);  // 只有 [1,2,3] 被包含
    }
```

- [ ] **Step 2: 实现**

```rust
/// 统计所有 PairPath 的路径中，每对相邻节点 (u,v) 被多少 read 覆盖。
pub fn edge_read_coverage(s: &SuffStats) -> FxHashMap<(i32, i32), i64> {
    let mut cov: FxHashMap<(i32, i32), i64> = FxHashMap::default();
    for inner in s.combined_read_hash.values() {
        for (pp, &count) in inner {
            for path in [&pp.path1, &pp.path2] {
                for w in path.windows(2) {
                    *cov.entry((w[0], w[1])).or_insert(0) += count;
                }
            }
        }
    }
    cov
}

/// 统计一条产出路径被多少 PairPath 兼容包含。
pub fn path_read_support(s: &SuffStats, path: &[i32]) -> i64 {
    let mut total = 0;
    for inner in s.combined_read_hash.values() {
        for (pp, &count) in inner {
            if pp.is_compatible_and_contained_by_single_path(path) {
                total += count;
            }
        }
    }
    total
}
```

- [ ] **Step 3: 绿 + 提交**

---

### Task 2: PogFlow 分支注入 read 支持（TDD + 验证）

**Files:** Modify `crates/trinity-butterfly/src/component.rs`

- [ ] **Step 1: digraph_to_flow 加 read 覆盖参数**

```rust
fn digraph_to_flow(
    g: &DiGraph,
    read_cov: &FxHashMap<(i32, i32), i64>,
) -> trinity_flow::graph::FlowGraph {
    let mut fg = trinity_flow::graph::FlowGraph::default();
    for &vid in g.vertex_ids() {
        if let Some(v) = g.get_vertex(vid) {
            fg.nodes.insert(vid as i64, trinity_flow::graph::FlowNode {
                kmer: v.name.as_bytes().to_vec(),
            });
        }
    }
    for &uid in g.vertex_ids() {
        for &vid in g.get_successors(uid) {
            if let Some(e) = g.find_edge(uid, vid) {
                let pog_w = e.weight.max(0.0).round() as u64;
                let read_w = read_cov.get(&(uid, vid)).copied().unwrap_or(0) as u64;
                // max 融合: read 覆盖高时不受 POG 低权截断
                let w = pog_w.max(read_w).max(1); // 至少 1（M2.1 终审建议）
                fg.edges.insert((uid as i64, vid as i64), w);
            }
        }
    }
    fg.component_id = 0;
    fg
}
```

- [ ] **Step 2: PogFlow 分支改造**

```rust
    if params.path_engine == PathEngine::PogFlow {
        // ① PairPath read 边覆盖（M3: read 支持注入流分解）
        let read_cov = edge_read_coverage(&suff);
        // ② POG 图 + read 覆盖融合权重
        let flow_g = digraph_to_flow(&result.seqvertex_graph, &read_cov);
        let kmers: Vec<(i64, Vec<u8>)> = flow_g.nodes.iter()
            .map(|(&id, n)| (id, n.kmer.clone())).collect();
        let mut all: Vec<trinity_flow::paths::PathOutput> = Vec::new();
        for comp in flow_g.weak_components() {
            // 单节点分量直接输出（M2.1 已有）
            if comp.len() == 1 {
                let seq = flow_g.nodes[&comp[0]].kmer.clone();
                if seq.len() >= params.min_len {
                    let support = path_read_support(&suff, &comp);
                    if support >= params.min_read_support {
                        all.push(trinity_flow::paths::PathOutput {
                            seq, support: support as u64, path: comp.to_vec(),
                        });
                    }
                }
                continue;
            }
            let paths = trinity_flow::peel::peel_component(&flow_g, &comp,
                &trinity_flow::peel::PeelParams { min_support: 1, min_path_nodes: 2 });
            for path in paths {
                let seq = trinity_flow::paths::path_to_seq(&path, &kmers, ctx.kmer_size);
                if seq.len() < params.min_len { continue; }
                // ③ read 支持过滤（M3 核心变更）
                let support = path_read_support(&suff, &path);
                if support < params.min_read_support { continue; }
                all.push(trinity_flow::paths::PathOutput {
                    seq, support: support as u64, path,
                });
            }
        }
        // dedup + FASTA（现有）
        let deduped = trinity_flow::paths::dedup_similar(all, params.flow_dedup_threshold);
        let fasta = trinity_flow::paths::to_fasta(&deduped, &params.name, 0, 60);
        return Ok(ComponentResult { all_prob_paths_fasta: fasta, num_paths: deduped.len() });
    }
```

**注意**: `suff` 变量已在 `run_component_inner` 中可用（`let suff = get_suff_stats_w_pairs(&graph, &ordered);`）。`params.min_read_support` 默认 2（与 classic 一致）。

- [ ] **Step 3: c0 fixture 测试仍通过**（read 支持注入后产出可能变化但应非空）

- [ ] **Step 4: 提交**

---

### Task 3: dev1 验收

```bash
A=/storage/home/senior007/test/athal/bench
cargo build --release -p trinity-cli && cargo build --release -p xtask
cargo run -q --release -p xtask -- bench --tag m3-dev1 \
  --reads $A/reads/dev1_1.fq.gz,$A/reads/dev1_2.fq.gz --ref $A/refs/dev1.fa \
  -- --path-engine pogflow --max_memory 8G
```

对照表:
| 引擎 | 期望 R | 期望 P | 期望 F1 |
|---|---|---|---|
| classic | 8.2% | 80.5% | 14.9% |
| pogflow (M2.1) | 8.1% | 86.0% | 14.8% |
| **m3** | **≥9%** | **≥84%** | **≥15.5%** |

Gate: m3 F1 > classic 14.9% 且 P ≥ 84% → 跑 dev2/3/holdout
      否则记录分析（read 支持过滤是否过严/过松）
