# P4 trinity-butterfly — Trinity Rust 复刻实施计划

> **状态（2026-08-18）: 11/11 完成。** 验证门全过：build OK / 481 测试 0 失败 /
> clippy -D warnings 干净 / fmt OK / xcheck-butterfly 3 PASS + 3 PASS-WARN（6 检查点）/
> xcheck-chrysalis 7/7 / xcheck-inchworm 4/4 / xcheck-kmer 3/3。
> **jar 裁定**（docs/setup.md）：oracle = 发布版 `Butterfly/Butterfly.jar`（md5 312f…）；
> 内层 jar tarball 原始 md5 3793…（计划早期"tarball 内无此文件"假设有误，已核实），
> 工作树内层 jar 为本地构建（794d…）覆盖了它，仅作对照。
> **c2 已知差异记录**: em/noem 各多 1 条 [2182] 短异构本——jar 的路径搜索/过滤
> 丢弃之，黄金序列全覆盖（PASS-WARN）；c0/em 仅 header/输出顺序差（Java HashMap
> 迭代序）。基准（docs/benchmarks.md）：c0 ~40×/RSS 1/27、c1/c2 ~2.5×/RSS 持平。

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]` syntax for tracking.

**Goal:** 实现 trinity-butterfly crate：TransAssembly_allProbPaths（15,920 行 Java 上帝类）+ 内嵌 jaligner 的完整移植，按 9 层依赖链推进；每层用原版 jar 的 8 个 DOT 检查点与组件级 allProbPaths 输出对拍，最终在 fixtures/p3 的 Chrysalis 产物上端到端复现 Butterfly 的转录本输出。

**Architecture:** 镜像式移植。图结构自研有向图（JUNG DirectedSparseGraph 语义：禁平行边、O(1) findEdge）；全部静态全局收进 `BflyContext`；f32 打分贯穿比对器；read 穿线递归保留（组件线程大栈）。Oracle：原版 Butterfly.jar（DOT 检查点 + 样本组件 allProbPaths）+ 小 Java driver 直调 jar 内 jaligner 类产黄金向量。

**Tech Stack:** 既有 workspace。Java oracle：`/public/home/senior007/miniconda3/envs/trinity/bin/java` + Butterfly.jar；javac 探测（jaligner 黄金 driver 需要；无 JDK 则手造小用例 + jar 行为级验证）。

---

## 移植契约（精读报告浓缩；各任务派发内嵌完整规格）

### 全局状态收编（BflyContext）
LAST_ID/LAST_REAL_ID、EDGE_THR=0.02/FLOW_THR=0.02、MAX_MM_ALLOWED（**每 read 重写的静态——作为穿线上下文参数传递**）、MAX_READ_SEQ_DIVERGENCE=0.05、MAX_READ_LOCAL_SEQ_DIVERGENCE=0.1、EXTREME_EDGE_FLOW_FACTOR=200、READ_END_PATH_TRIM_LENGTH、USE_DP_READ_TO_VERTEX_ALIGN、nodeTracker/origIDnodeTracker、原 KMER_SIZE

### 逐字保留的风险点（精读报告清单）
1. graph.reads 解析 `endInRead = fields[3] + KMER_SIZE` 的 off-by-one（FIXME 保留）
2. removeLightFlowEdges 严格 `<` vs In/Out 版 `<=`；In/Out 删除中不重算 total
3. removeSingleNtBubbles 平局保 v2；addToPrevIDs 的 `id >= lastRealID` 疑似 bug 方向**保留**
4. NW banded：行末 vDiagonal=0 重置、首行宽松初始化、traceback 提前终止、**float** 平局 tie-break（Diag>UP>Left）
5. updatePathRecursively 后继选择 `<=`（平局取迭代序**最后一条**）；memo 命中返回**深拷贝**；tied 分支直接改子对象
6. ZipperAlignment 右锚循环的 i 下标 quirk
7. compactLinearPaths 迭代中改图（快照+while 语义等价）
8. My_DFS 双向递归 + min-depth 平移 + down-only 两阶段（最终深度来自第二阶段）
9. concatVertex 的 weights/prevID 合并序；getNameKmerAdj 去 K-1 前缀
10. anchors 再锚定滑窗从 i=2 起

### DOT 检查点对拍（验证基础设施）
A 建图 → B 剪轻边 → C 压缩+DFS → H SNP 坍缩 → D 小组件清除 →（POG 系）→ I zipping →（repeat unroll rN）→ Z 最终
节点集+边集（含权重）语义比较；jar 需 `--generate_intermediate_dot_files`

### 组件级端到端
fixtures/p3/quantify 的 c0/c1/c2（.graph.out/.graph.reads 已有）+ 原版 jar 跑出的 allProbPaths 固化为黄金；Butterfly 自带 sample_data/c1.graph 亦为起点 fixture。

## 文件结构

```
crates/trinity-butterfly/
├── Cargo.toml
└── src/
    ├── lib.rs
    ├── context.rs        # T1: BflyContext（全局状态收编）
    ├── graph.rs          # T1: DiGraph（JUNG 语义）+ SeqVertex + SimpleEdge
    ├── graph_io.rs       # T1: preProcess/buildNew + reads 解析 + DOT 写出
    ├── pruning.rs        # T2: 剪枝链 + compactLinear + SNP bubbles + 小组件
    ├── dfs.rs            # T3: My_DFS 移植
    ├── align.rs          # T4: NW-Gotoh banded(f32)/zipper/AlignmentStats
    ├── read_threading.rs # T5: updatePathRecursively 全链
    ├── pair_paths.rs     # T6: PairPath/Read 兼容性族 + getSuffStats_wPairs
    ├── pog.rs            # T7: POG 构建/破环/转 DAG/zipping
    ├── paths.rs          # T8: getAllProbablePaths + triplet
    ├── postprocess.rs    # T9: cdhit 式去冗余/基因分组/输出
    └── bin/butterfly.rs  # T10: CLI（原版参数面）
fixtures/p4/              # DOT 检查点黄金 + jaligner 黄金 + 组件输出黄金
xtask/src/main.rs         # T10: xcheck-butterfly
docs/porting-map.md       # T11
```

### 任务划分（9 层依赖链 + 3 个支撑任务）

- **T1** 图基础：DiGraph/SeqVertex/SimpleEdge/BflyContext/graph IO/DOT 写出 → 检查点 A 对拍
- **T2** 剪枝链（fix/removeLight×3/compact/SNP/小组件）→ 检查点 B/C/H/D 对拍
- **T3** My_DFS（visitVertex2 双向 + down-only 两阶段）
- **T4** 比对器（NW banded f32/zipper/AlignmentStats）+ jaligner 黄金（Java driver 或手造）
- **T5** read 穿线（最重：600 行递归 + memo + 三级比对切换）→ graph.reads 消费 + read paths 抽样对拍（jar -V 输出?）
- **T6** PairPath 族 + 配对组装
- **T7** POG（overlap layout/破环/zipping）→ POG/I 系检查点
- **T8** getAllProbablePaths + triplet → allProbPaths 雏形
- **T9** 后处理（去冗余/基因分组/_g_i 命名/输出格式）
- **T10** CLI + xcheck-butterfly（DOT 检查点全链 + 组件端到端 + eval 统计）
- **T11** 验证门 + porting-map/benchmarks + 最终审查 + 并回 main

### 对拍判定
- 早期层：DOT 节点/边集（含 round(weight)）
- read 穿线：无直接可观测输出（中间态）——通过 T8 输出与 POG 检查点间接验证 + jar 的 -V debug 输出抽样比对
- 终局：组件 allProbPaths.fasta 的序列多重集（header 的 path=[...] 可能因 tie-break 差异——按 spec 允许集合级比较）+ 基因/异构体计数统计
- 验收（spec §11 P4）：c1.graph 检查点全过 → sample_data（fixtures/p3）全组件 allProbPaths 多重集高一致

### 已知不可复刻项（预先声明）
- Java HashMap/JUNG 迭代序相关的 tie-break（read 穿线后继平局取末、POG 边序）——输出集合稳定、个别路径选择可能不同；按序列多重集评估
- 线程模型：原版进程级并行（ParaFly 每组件一 JVM）；我们组件线程池内单线程 + 大栈（512MB 可配）
