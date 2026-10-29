# M4 K19→K25 超级 read 递进 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `--super_reads` 旗标：K19 预组装的 ≥200bp contig 作为超级 read 加入 K25 Chrysalis reads 池（参与 RTT 投票 + QuantifyGraph 边权），dev1 F1 ≥ 15.4%（M3 基线 15.1%）。

> **状态（2026-08-24 回填）:** 实现全部完成并提交（da5b19d / 3619bb8），**Gate FAIL**——dev1 F1 14.8% < 15.4% 目标且低于 M3 基线 15.1%，按计划记录分析、未跑 dev2/3，`--super_reads` 保留为默认 off 的 opt-in flag。详见 docs/opt-log.md M4 段。

**Architecture:** orchestrate 归一化后插入 K19 预组装段（k-mer 计数 + inchworm + min200 过滤），产出 `super_reads.fa`；Chrysalis 阶段将 both.fa 与 super_reads.fa 合并喂入。Butterfly/PogFlow 无改动。

**Tech Stack:** 既有全部组件；无新 crate。

---

## 文件结构

```
crates/trinity-cli/src/{args,orchestrate}.rs  # --super_reads 参数 + 预组装段 + reads 合并
docs/opt-log.md                                # 验收数字
```

---

### Task 1: `--super_reads` 参数 + K19 预组装段

**Files:** Modify `crates/trinity-cli/src/{args,orchestrate}.rs`

- [x] **Step 1: args.rs** `pub super_reads: bool`（`--super_reads` 解析，默认 false）

- [x] **Step 2: orchestrate.rs 归一化后、K25 计数前插入**（实现偏差见下: 过滤函数换 `filter_min_len`; 断点命中时从盘补读 super_reads.fa）

```rust
// ---- 4a. K19 预组装（M4 超级 read 递进）----
let mut super_reads_data: Vec<u8> = Vec::new();
if args.super_reads {
    let sr_ok = outdir.join(".k19_preasm.ok");
    let sr_fa = outdir.join("super_reads.fa");
    run_with_checkpoint(&sr_ok, "K19 pre-assembly (super reads)", || {
        if checkpoint_exists(&sr_fa) {
            super_reads_data = fs::read(&sr_fa)?;
            return Ok(());
        }
        // K=19 k-mer 计数（归一化后的 both.fa/single.fa）
        let k19 = 19usize;
        let counts19 = KmerCountTable::count_fasta_data(&reads_data, k19, !ss.is_none());
        // K=19 inchworm（排序种子单线程——预组装无需并行）
        let mut counter19 = populate_from_kmers(&dump_bytes(&counts19, k19, 1), k19, !ss.is_none())?.0;
        prune_some_kmers(&mut counter19, args.min_kmer_count.max(1), 1.0, false, 0.005);
        let mut buf = Vec::new();
        trinity_inchworm::irke::compute_sequence_assemblies(
            &mut counter19, &IrkeParams::default(),
            &AssemblyParams { min_assembly_length: k19, ..Default::default() },
            &Monitor::new(0), true, &mut buf)?;
        // ≥200bp 过滤（spec: super_reads 定义）
        let filtered = filter_iworm_min_len_or_cov(&buf, 200, 0)?;
        // header 改写: >a<N>;cov → >sr_<N>（避免与 /1//2 配对冲突）
        let sr = rewrite_super_read_headers(&filtered)?;
        fs::write(&sr_fa, &sr)?;
        super_reads_data = sr;
        Ok(())
    })?;
}
```

- [x] **Step 3: `rewrite_super_read_headers` 辅助函数**（orchestrate.rs; 实际落成 `pub`, 计划骨架为私有）

```rust
/// >aN;cov total_counts:... → >sr_N（M4: 超级 read 无 PE 后缀, 防 Butterfly 配对干扰）
fn rewrite_super_read_headers(data: &[u8]) -> Result<Vec<u8>, CommonError> {
    let text = std::str::from_utf8(data)
        .map_err(|e| CommonError::Parse(format!("super reads UTF-8: {e}")))?;
    let mut out = Vec::with_capacity(data.len());
    for line in text.lines() {
        if line.starts_with('>') {
            // >a123;45 → >sr_123
            let acc = line[1..].split(';').next().unwrap_or("");
            let num = acc.strip_prefix('a').unwrap_or(acc);
            out.extend_from_slice(format!(">sr_{num}\n").as_bytes());
        } else {
            out.extend_from_slice(line.as_bytes());
            out.push(b'\n');
        }
    }
    Ok(out)
}
```

- [x] **Step 4: Chrysalis 阶段合并 reads**（orchestrate.rs `chrysalis_reads` 段）

在 `run_chrysalis_pipeline` 调用处，将 `super_reads_data` 追加到 `reads_data`：

```rust
let mut chrysalis_reads = reads_data.clone();
if args.super_reads && !super_reads_data.is_empty() {
    chrysalis_reads.extend_from_slice(b"\n");
    chrysalis_reads.extend_from_slice(&super_reads_data);
}
// 用 chrysalis_reads 替代原来的 reads_data 喂 chrysalis
```

- [x] **Step 5: 测试**（tests/super_reads_e2e.rs: on/off + classic 不受影响）

合成集集成测试（复用 recursive_e2e 手法）：
- `--super_reads` on → `super_reads.fa` 存在且 header 含 `>sr_`；`.k19_preasm.ok` 存在
- 不开 super_reads → 无 super_reads.fa
- classic 路径不受影响

- [x] **Step 6: 提交**（da5b19d; 过滤 bug 修复 3619bb8 + fmt c31b647）

```bash
git add crates/trinity-cli
git commit -m "feat(cli): --super_reads K19 预组装超级 read 递进（M4）"
```

---

### Task 2: dev1 验收

```bash
A=/storage/home/senior007/test/athal/bench
cargo build --release -p trinity-cli && cargo build --release -p xtask
cargo run -q --release -p xtask -- bench --tag m4-dev1 \
  --reads $A/reads/dev1_1.fq.gz,$A/reads/dev1_2.fq.gz --ref $A/refs/dev1.fa \
  -- --path-engine pogflow --super_reads --max_memory 8G
```

Gate: F1 ≥ 15.4% → dev2/3 跑齐；否则记录分析（超级 read 数量/效果）

> **结果: FAIL** — dev1 F1 **14.8%**（R=8.1 P=84.1, N50=1727 n=7874）< 15.4% 目标, 亦低于 M3 基线 15.1%。按计划记录分析、不跑 dev2/3。首跑踩坑: 计划骨架的 `filter_iworm_min_len_or_cov(buf, 200, 0)` 因 min_cov≥0 恒真放行全部 112,666 条 contig——改用专用 `filter_min_len`（3619bb8）后 F1 仍 14.8%, 判定超级 read 本身在 dev1 高深度 PE 下无净增益（分析见 opt-log M4 段）。

### Task 3: 验收入档 + 并 main

opt-log 对照表（M3 vs M4）+ holdout + 最终审查

> **结果: 完成** — opt-log M4 段入档（7acb2b2 验收 + f2c1605 负结果分析）。工作直接提交在 main 上, "并 main" 无需额外合并。holdout 未跑（Gate FAIL 后按计划止步 dev1）。
