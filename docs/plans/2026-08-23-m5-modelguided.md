# M5 模型引导组装 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subunit-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]` syntax for tracking.

**Goal:** `--model-guided` 模式：reads 招募到模型转录本的基因分组 → 每组定向重组装 → 合并。dev1 F1 ≥ 15.5%。

> **状态（2026-08-24 回填）:** 全部完成，Gate PASS。dev1 F1 15.52%（≥15.5 目标，历史最佳），dev2 14.4% / dev3 14.1% 均超 classic，holdout 3.9% 持平。详见 docs/opt-log.md M5 段。

**Architecture:** orchestrate 加 model-guided 分支：归一化+prep 共享 → [自引导 or 外部模型] → minimap2 招募 → 基因分组 → 每组 mini-assembly（复用 run_mini 但输入为组内 reads）→ 合并 harvest。

**Tech Stack:** 既有组件 + minimap2（外部进程，读 PAF）。

---

## 文件结构

```
crates/trinity-cli/src/{args,orchestrate,model_guided}.rs
docs/opt-log.md
```

---

### Task 1: model_guided.rs 核心模块

**Files:** Create `crates/trinity-cli/src/model_guided.rs`

- [x] **Step 1: 基因聚类**

```rust
/// 模型转录本按 ≥50% 序列重叠聚类成"基因"（弱连通分量）。
/// 输入: 模型 fasta 字节; 输出: Vec<Vec<usize>>（每组是 transcript 索引列表）。
pub fn cluster_transcripts_to_genes(model_fa: &[u8]) -> Result<Vec<Vec<usize>>, CommonError>
```
实现：minimap2 self-alignment（`-x asm5 --secondary=no`），PAF 中 identity ≥80% 且 target 覆盖 ≥50% 的 transcript 对建边 → 并查集弱连通分量。

- [x] **Step 2: reads 招募**

```rust
/// reads 比对到模型 → 按 best hit 分配到基因。
/// 返回: (组→reads 字节, orphan reads 字节)。
pub fn recruit_reads(
    reads_fa: &[u8],           // both.fa 或 single.fa
    model_fa: &[u8],
    gene_groups: &[Vec<usize>], // 基因聚类结果
) -> Result<(Vec<Vec<u8>>, Vec<u8>), CommonError>
```
实现：minimap2 `-x sr -t 8` reads vs model → PAF → 每 read 的 best hit（最高 identity×coverage）→ 找到 transcript → 找到其基因 → 加入组。PE 对：read 名去 /1//2 后缀匹配，两端命中取并集（任一端命中即入组）。

- [x] **Step 3: 组装合并**

```rust
/// 每组跑 mini-assembly（inchworm→chrysalis→pogflow）→ 合并 + 去冗余。
pub fn assemble_groups(
    groups: &[Vec<u8>],        // 每组的 reads fasta 字节
    orphan: &[u8],
    params: &MiniParams,
) -> Result<Vec<u8>, CommonError>
```
实现：每组 `run_mini`（recursive.rs 已有），结果拼合 → 全局去冗余（min_len + 全序列 HashSet）。

- [x] **Step 4: 测试**（合成集手推; 实际落了 3 个单测: cluster_and_recruit_synthetic / pe_union_uses_either_mate / orphan_reads_assembled）

```rust
#[test]
fn cluster_and_recruit_synthetic() {
    // 3 条转录本: 2 条属于同一基因（60% 重叠）, 1 条独立
    // 10 对 reads 平铺到各转录本
    // 验证: 基因聚为 2 组; reads 正确分配; orphan 组为空
}

#[test]
fn orphan_reads_assembled() {
    // 2 条模型转录本 + 1 对 reads 来自第三条"新"转录本
    // 验证: orphan 组产出该 reads 的组装
}
```

---

### Task 2: orchestrate 接线

**Files:** Modify `crates/trinity-cli/src/{args,orchestrate}.rs`

- [x] **Step 1: args.rs** `pub model_guided: bool`（`--model-guided`）+ `pub model: Option<PathBuf>`（`--model <fa>`）

- [x] **Step 2: orchestrate.rs** 在归一化+prep 后（both.fa 就绪后）分支：

```rust
if args.model_guided {
    return run_model_guided(a, args, &outdir, &both_fa_bytes, &ss, ...);
}
```

`run_model_guided` 流程：
1. 模型获取：`args.model` 存在则读文件；否则先跑完整 PogFlow 流程（复用 classic 段代码但 `--path-engine pogflow`）→ 取 Trinity.fasta 作模型
2. `cluster_transcripts_to_genes` → `recruit_reads` → `assemble_groups`
3. 输出命名/gene_trans_map 与现有格式一致
4. checkpoint `.model_guided.ok`

- [x] **Step 3: e2e 测试**（合成集; 当初以手工 e2e 验证并入档 opt-log——3 转录本 → 2 基因组, 150/150 reads 入组, 产出 3 条转录本 + gene_trans_map。**2026-08-24 补齐自动化用例 tests/model_guided_e2e.rs**, 对齐 M4 的 super_reads_e2e.rs 惯例）

---

### Task 3: dev1 验收 + 入档

```bash
A=/storage/home/senior007/test/athal/bench
cargo run -q --release -p xtask -- bench --tag m5-dev1 \
  --reads $A/reads/dev1_1.fq.gz,$A/reads/dev1_2.fq.gz --ref $A/refs/dev1.fa \
  -- --path-engine pogflow --model-guided --max_memory 16G
```
Gate: F1 ≥ 15.5% → dev2/3；否则记录分析。

> **结果: PASS** — dev1 F1 **15.52%**（R=8.6 P=83.6, N50=1735 n=8165）→ 跑齐 dev2（14.4%）/ dev3（14.1%）/ holdout（3.9%, 持平 classic）。全部入档 opt-log M5 段。
