# P3 trinity-chrysalis — Trinity Rust 复刻实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **状态（2026-08-17，P3 验证门通过）**: **8/8 任务全部完成**（Task 0–8，分支 p3-chrysalis）。
> 验证门：build ok、test 340/340、clippy -D warnings 通过、fmt 通过、xcheck-chrysalis 7/7、
> xcheck-inchworm 回归 4/4、xcheck-kmer 回归 3/3。基准与对拍数字见 `docs/benchmarks.md`
> （2026-08-17 chrysalis 全链路条：-t 1 快 ~1.2× / RSS ~1/3.2，-t 4 慢于原版 OMP——
> ReadsToTranscripts 并行度为 P4 优化项），移植映射见 `docs/porting-map.md` trinity-chrysalis
> 节。本计划以 `### Task N` 标题组织（无 `- [ ]` 步骤 checkbox，故无勾选动作——
> 以本状态行代之）。

**Goal:** 实现 trinity-chrysalis crate：六个子阶段的完整移植（GraphFromFasta 焊接聚类 / BubbleUpClustering / CreateIwormFastaBundle / ReadsToTranscripts / FastaToDeBruijn+DeBruijnGraph / QuantifyGraph + 分区器），CLI 与原版二进制同名同参；每阶段输出与原版对拍，全链路（inchworm.fa → c*.graph.out/.reads）可跑通且中间文件格式逐字段兼容。

**Architecture:** 镜像 `Chrysalis/analysis/*` + `Inchworm/src/{FastaToDeBruijn,DeBruijnGraph}`（后者物理在 Inchworm，按 spec 放 trinity-inchworm crate 的 debruijn.rs，由 chrysalis 编排调用）。排序责任在管线（原版 shell sort）——本阶段提供库内等价排序（weld 图 total 降序；readsToComponents 三键）。逐位复刻原版怪癖（清单见下）。

**Tech Stack:** 既有 workspace。Oracle：已编译原版 Chrysalis/bin 六个二进制 + Inchworm/bin/FastaToDeBruijn。

---

## 移植契约（两份精读报告浓缩；完整规格在各任务派发时内嵌）

### 管线数据流与排序责任
```
inchworm.fa.min100 → [GraphFromFasta] → welds_graph.txt
  → sort -k9,9gr（total 降序） → [BubbleUpClustering] → GraphFromIwormFasta.out（COMPONENT 块）
  → [CreateIwormFastaBundle -min 200] → bundled_iworm_contigs.fasta（s_N + 'X' 连接，无折行）
  → [ReadsToTranscripts -p 50 -min_kmer_entropy 1.5] → readsToComponents.out
  → sort -k1,1n -k3,3nr -k2,2 → .sort
  → [FastaToDeBruijn --graph_per_record -K 24] → .deBruijn（5 列图，节点 id 1-based）
  → [partition -N 1000 -L 200] → Cbin_*/c*.graph.tmp + c*.reads.tmp + component_base_listing.txt
  → [QuantifyGraph -k 24 -max_reads 200000] → c*.graph.out + c*.graph.reads（Butterfly 输入）
```
注意：`-p` 主线钳制到 50（Trinity:1044）；`-L` 主线 200。

### 必须逐位复刻的怪癖（精读报告"易踩坑"清单）
1. **ReadsToTranscripts 最长 run 的 off-by-one/two**：`Sort(comp)` 后组尾判定 `run=m-1`（非末组）、`run=m-2`（末组）；大小 1 的组永不成为 best
2. **BubbleUpClustering 合并条件 `A+B+2 <= MAX_CLUSTER_SIZE`**（合并前 size 和再 +2）；单端加入是"加入前判定 `<`"；`while(!eof())` 产生 `0->0` 伪边（需决策复刻或修复并文档化）
3. **QuantifyGraph 熵检查 `strncpy(&d[1],25)` 越界一字节**（24 真实碱基+1 任意字节，分母 25）——决策：按"24 字符+分母 25"复刻并注释
4. **SortPrint 输出条件 `lastStart > lastStartTemp`**（组内末 start 严格大于首 start——单 kmer 位置命中的 read 被丢弃）；ori=-1 组排前；pos2/node2=组内最后命中
5. **PrintSeq 80 整数倍长度多一空行**；`#POOL_INFO` 各成员带尾随空格
6. DeBruijn：'-' 定向 id = id + N（非负数！）；-1=终端哨兵；weight/flag 初值恒 1
7. minCov 截断 `(int)`；IsShadow 的整除链（`len/25 - 1` 整除、`nn < n/5`）；per_id `> 97` 严格；`a/b > 0.05` 严格；`d.len/10 > dd.len` 整除
8. WeldableKmer 越界 `stopB >= b.len()`（不是 >）
9. Phase 1 IsSimple 无条件跳过 vs Phase 2 `IsSimple && 总命中>1` 才跳过
10. NonRedKmerTable 计数前 read toupper；DNAStringStreamFast 多行序列无分隔拼接
11. ReadsToTranscripts 的 `pct = (int)((float)max/num_kmer_pos*100 + 0.5)`，分母只算正向 k-mer 位置数
12. component_id 链：BubbleUpClustering 输出序号（过滤留空洞，不重编号）= s_N 的 N = readsToComponents 第 1 列

### 非确定性点（输出集合稳定、行序可变，注释声明）
- GraphFromFasta weld 图行序（下游会重排）——我们选确定性序
- FastaToDeBruijn Component 块序（omp）——我们顺序输出
- ReadsToTranscripts 共享 25-mer 的 bundle 标签竞争——按单线程"后写者赢"语义
- priority_queue 平局弹出序——行序不影响行集合

## 文件结构

```
crates/trinity-common/src/cli.rs        # ★ T0: 从 kmer/inchworm 两份手写解析器抽取的共享模块
crates/trinity-inchworm/src/debruijn.rs # T5: DeBruijnGraph + toChrysalisFormat（原版物理位置）
crates/trinity-chrysalis/
├── Cargo.toml                          # trinity-common, trinity-inchworm, rayon, rustc-hash
└── src/
    ├── lib.rs
    ├── kmer_align.rs                   # T1: KmerAlignCore（12-mer CSR 倒排索引）
    ├── nonred_table.rs                 # T1: NonRedKmerTable（排序数组+二分，双语义）
    ├── dna_vector.rs                   # T1: vec 读入（header '_' 连接/大写）、DNAStringStreamFast 等价流
    ├── graph_from_fasta.rs             # T2: Welder/Phase1/计数/Phase2/report + sort_welds
    ├── bubble_up.rs                    # T3: 贪心聚类 + COMPONENT 输出
    ├── bundle.rs                       # T3: COMPONENT→bundle
    ├── reads_to_transcripts.rs         # T4: k-mer 投票
    ├── quantify.rs                     # T6: QuantifyGraph
    ├── partition.rs                    # T6: Cbin 分区 + listing
    └── bin/trinity-chrysalis.rs        # T7: 六子命令 CLI
xtask/src/main.rs                       # T7: xcheck-chrysalis
fixtures/p3/                            # oracle 产物
docs/porting-map.md / benchmarks.md     # T8
```

### 任务划分（每任务一个实现者 + spec/质量双审）

- **T0** CLI 解析器抽取（trinity-common/src/cli.rs；kmer/inchworm 两个 bin 改用；纯机械）
- **T1** 基础层：dna_vector（读入/流式/revcomp/per_id/encapsulates/IsSimple/熵）+ KmerAlignCore（CSR）+ NonRedKmerTable
- **T2** GraphFromFasta（最重：Phase1 候选 weldmer → 计数 → Phase2 验证建边 → report + `sort_weld_graph` 库函数）
- **T3** BubbleUpClustering + CreateIwormFastaBundle（COMPONENT 读写闭环）
- **T4** ReadsToTranscripts（含最长 run 怪癖 + `sort_reads_to_components` 库函数）
- **T5** FastaToDeBruijn + DeBruijnGraph（inchworm crate 的 debruijn.rs；canonical 方向翻转/负向 id 偏移/优先队列遍历）
- **T6** QuantifyGraph（25-mer 拼查询/二分收集/SortPrint）+ partition
- **T7** CLI bin 六子命令 + xcheck-chrysalis（每阶段 oracle 对拍 + 全链路端到端 + 喂 Butterfly jar 冒烟）
- **T8** P3 验证门（test/clippy/fmt/xcheck 4+3 全绿 + 链路跑通）+ porting-map/benchmarks + 最终审查 + 并回 main

### 每阶段对拍判定
- GraphFromFasta：weld 边集合（(A,B) 无序对多重集 + weld/scaff 计数）相等；行序不管（下游重排）
- BubbleUp：COMPONENT 块多重集（按成员集合分组）+ 组件数
- Bundle：逐字节（确定性输出）
- ReadsToTranscripts：(comp, read) 集合 + pct 值分布；行序按我们的 sort 键自排
- FastaToDeBruijn：逐 Component 块的边集合（5 列去 id 偏移影响后的多重集）——块序不管
- QuantifyGraph：graph.out 逐行（除第 3 列）+ graph.reads 按行集合（SortPrint 顺序敏感处按集合）
- 端到端：我们的 inchworm.fa（P2 产物）→ 全链路 → c*.graph.out 喂原版 Butterfly jar 跑通一个组件

### 验收（spec §11 P3）
每子阶段输出 vs 原版逐个比对通过；xcheck-chrysalis 一键全绿；链路端到端 + Butterfly 冒烟；基准入档。
