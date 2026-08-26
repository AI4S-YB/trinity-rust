# Trinity RNA-Seq 组装器纯 Rust 复刻 — 设计文档

日期：2026-08-15
状态：已与用户逐节确认
原版参考：trinityrnaseq v2.15.2（源码位于 `/storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2/`，含 FULL 源码与测试数据）

## 1. 背景与目标

Trinity 是 RNA-seq de novo 转录组组装的标准工具，由三阶段流水线构成：**Inchworm**（贪心 k-mer 延伸组装初生 contig）→ **Chrysalis**（contig 聚类成组件并建 de Bruijn 图、reads 回配）→ **Butterfly**（在组件图上用 read 支持度恢复全长异构体）。原版实现为 Perl 主控 + C++（Inchworm/Chrysalis）+ Java（Butterfly），并依赖 jellyfish、bowtie2、samtools、ParaFly、seqtk、salmon 等外部工具。

本项目的目标：**用纯 Rust 复刻 Trinity 组装主线，达到生产可用**（正确性 > 健壮性 > 速度），交付单二进制 `trinity`，不依赖任何外部工具与 C 库。

原版核心代码规模（已实测）：

| 模块 | 语言 | 行数 | 形态 |
|---|---|---|---|
| Inchworm | C++ | 7,231 | 单二进制（IRKE）+ 辅助工具（FastaToDeBruijn 等） |
| Chrysalis | C++ | 4,162 | 6 个独立二进制 |
| Butterfly | Java | 28,049 | 单 jar（含 15,920 行上帝类 TransAssembly_allProbPaths + 内嵌 jaligner） |
| Trinity 主控 | Perl | 4,278 | 编排脚本 |

预估 Rust 规模：1.8 万–2.5 万行（含测试）。

## 2. 范围

**纳入**：
- 组装主线全流程：读预处理（SS 方向规整）→ in silico 归一化（DigiNorm）→ k-mer 计数（jellyfish 替代）→ Inchworm → Chrysalis（6 子阶段）→ Butterfly → 汇总输出 `Trinity.fasta` + gene_trans_map；
- 输入模式：PE（`--left/--right`）与 SE（`--single`）；链特异（`--SS_lib_type` F/R/FR/RF）与默认非链特异（DS，即 canonical k-mer 模式）；
- FASTA 与 FASTQ（含 gzip）输入；
- 断点续跑（复刻 `.ok` checkpoint 语义）。

**排除**（有意裁剪，非主线或有等价外部依赖）：
- 基因组指导模式（--genome_guided）、DNA 组装模式；
- PasaFly / CuffFly / cuff_no_extend 等 Butterfly 变体模式（仅实现默认路径 + `--all_possible_paths` 等四档路径扩展模式）；
- bowtie2 PE scaffolding（Chrysalis Stage 0；需要 bowtie2 级别比对器，工程量不匹配）；
- salmon 表达量过滤（原版 2.15.2 默认开启，但 salmon 为独立 C++ 项目；**本复刻默认等价于 `--no_salmon` 行为**——这是与原版默认行为的唯一已知功能性偏差）；
- Jaccard clip、长读（LR）支持、EM 表达量估计（原版主流程默认 `--NO_EM_REDUCE`，EM 代码保留为可选后置项）；
- grid/HPC 调度集成。

## 3. 兼容性标准（验收定义）

**中间文件格式逐字段兼容**：各阶段中间文件（`jellyfish.kmers.K.fa`、`inchworm.kmers.fa`、`GraphFromIwormFasta.out`、`bundled_iworm_contigs.fasta`、`.deBruijn`、`readsToComponents.out`、`c*.graph.out`、`c*.graph.reads`、`*.allProbPaths.fasta` 等）与原版格式完全一致，可与原版模块**双向互喂**交叉验证。

**不追求**逐字节相同的最终输出——原版自身受哈希迭代序、OpenMP 竞争、`rand()` 平局打破影响不可复现。最终结果验收用序列多重集一致性量化（见 §11）。

## 4. 总体架构

### 4.1 Cargo workspace 布局

```
trinity-rust/  (Cargo workspace)
├── crates/
│   ├── trinity-common    # 共享基础：2bit k-mer 编码、revcomp、熵、FASTA/FASTQ IO
│   ├── trinity-kmer      # jellyfish 替代：并行 k-mer 计数 + dump + DigiNorm
│   ├── trinity-inchworm  # IRKE 移植 + FastaToDeBruijn
│   ├── trinity-chrysalis # 6 子阶段
│   ├── trinity-butterfly # TransAssembly_allProbPaths 移植
│   └── trinity-cli       # `trinity` 主控（替代 Perl 脚本）
├── xtask/                # 交叉验证/评估工具（xcheck、gen-fixtures、eval）
├── fixtures/             # 黄金向量与小型测试数据
└── docs/
```

每个阶段 crate 同时提供 **library + 可独立运行的 bin**：库供主控进程内调用（无子进程开销），bin 供与原版互喂交叉验证。

### 4.2 流水线阶段映射

| # | 阶段 | 原版实现 | Rust 实现 | 关键输出（格式兼容） |
|---|---|---|---|---|
| 0 | 预处理（FQ→FA、SS 规整、PE 合并） | Perl prep_seqs | trinity-cli | `both.fa` / `single.fa` |
| 0.5 | DigiNorm 归一化 | jellyfish + Perl + seqtk | trinity-kmer | `*.normalized.fa` |
| 1 | k-mer 计数 | jellyfish count/dump | trinity-kmer | `jellyfish.kmers.K.fa`（header=计数） |
| 2 | Inchworm 组装 | inchworm 二进制 | trinity-inchworm | `inchworm.kmers.fa`（`>aN;cov` header） |
| 3 | Chrysalis | 6 C++ 二进制 + sort + Perl 分区 | trinity-chrysalis | `Cbin*/c*.graph.out` + `.graph.reads` |
| 4 | Butterfly | ParaFly + Java jar | trinity-butterfly（组件线程池） | `*.allProbPaths.fasta` |
| 5 | 汇总 | print_butterfly_assemblies.pl | trinity-cli | `Trinity.fasta` + `.gene_trans_map` |

### 4.3 关键技术决策

- **并行模型**：rayon 替代 OpenMP（k-mer 计数、read 解析、inchworm PARALLEL_IWORM）；Butterfly 组件间并行用线程池替代 ParaFly，**每组件单线程**（镜像原版"一组件一 JVM"语义，规避其全局可变状态）；
- **进程内编排**：主控直接调用各阶段库 API，各 bin 仅用于交叉验证与手工调试；
- **断点续跑**：`.ok` checkpoint 文件与原版同名同位；
- **CLI**：参数与原版同名同义（`--seqType --left --right --single --SS_lib_type --CPU --max_memory --KMER_SIZE --min_kmer_count --min_contig_length --output --no_normalize_reads ...`），未知参数报错；
- **依赖策略**：纯 Rust 依赖树（flate2 用 rust backend = miniz_oxide；无任何 C 绑定），发行单二进制；
- **MSRV**：构建时锁定的当前稳定版。

## 5. 各 crate 内部设计

**镜像原则**：每个 Rust 模块/函数对应原版文件/方法，注释标注原版行号（如 `// IRKE.cpp:933 inchworm_step`）。

### 5.1 trinity-common（~1.5k 行）
- `kmer.rs`：`KmerId = u64`；编码 **G=0, A=1, T=2, C=3**（互补 = 按位取反）；`revcomp_val`（`~kmer` + 2-bit 组反转）；`canonical = max(k, rc(k))`；香农熵（log2）。直译 `Inchworm/src/sequenceUtil.cpp`。
- `fasta.rs`：流式 reader（并行分块）；writer（60/80 列折行，对齐 `add_fasta_seq_line_breaks`）。
- `fastq.rs`：FASTQ/gz 读取。
- 单测：用原版行为黄金向量断言。

### 5.2 trinity-kmer（~2k 行）
- `counter.rs`：分片开放寻址哈希表（u64 kmer + u32 count ≈ 12B/kmer；按 kmer 哈希分 N 片、片级独立锁）；rayon 按片并行；canonical 模式入口规范化（镜像 `jellyfish --canonical`）；内存护栏复刻原版 `(jellyfish_ram − file_size)/7` 估算逻辑，超限清晰报错。
- `dump.rs`：jellyfish dump 兼容输出（`>count\nkmer` FASTA，`-L` 下限过滤）。
- `diginorm.rs`：每 read 覆盖统计（均值/中位数/标准差）→ nbkc 选择规则：`med_cov ≥ min_kmer_cov(=2)`、`mean_cov > 0`、`CV ≤ max_CV`，接受概率 `max_cov/med_cov`（默认 max_cov=50）。注：概率接受用 Rust rand（原版 Perl `srand(12345)` 无法逐值复刻，属协议允许的固有不确定性；保留统计等价）。

### 5.3 trinity-inchworm（~3k 行）
- `kmer_counter.rs`：惰性删除（count 置 0）；`prune_some_kmers`（min_count / min_entropy / 错误 kmer 比率 0.005）；DS 模式所有入口 canonical 化。
- `irke.rs`：种子列表（count 降序；PARALLEL 模式哈希迭代序）；`compute_sequence_assemblies`：rayon 并行 + 每线程私有缓冲 + `generateHash` 去重（镜像 PARALLEL_IWORM 两阶段种子策略）。
- `greedy.rs`：`inchworm_step` 候选位运算生成（前向 `(seed << (33−K)·2) >> (32−K)·2 | i`；后向 `(i << (K·2−2)) | (seed >> 2)`；候选序 G,A,T,C）；平局递归加深（MAX_RECURSION 1 → 硬上限 50）→ 随机打破（**固定种子 RNG，保证单线程完全可复现**）。
- `debruijn.rs`：FastaToDeBruijn + DeBruijnGraph（24-mer 节点、4-bit prev/next 掩码、`toChrysalisFormat` 优先队列遍历、DS 模式负向节点 id 偏移）。
- bins：`inchworm`、`FastaToDeBruijn`（名字对齐原版）。
- 关键参数默认值：K=25、MIN_SEED_ENTROPY=1.5、MIN_SEED_COVERAGE=2、MIN_ASSEMBLY_COVERAGE=2、`exceeds_min_connectivity` 的 `<1e5` 短路原样保留。

### 5.4 trinity-chrysalis（~4k 行，6 子命令）
```
bin 子命令: graph-from-fasta | bubble-up-clustering |
            create-iworm-fasta-bundle | reads-to-transcripts | quantify-graph
（partition 逻辑内置于库, 由主控调用）
```
- `kmer_align.rs`：KmerAlignCore 移植——2×12-mer（4^12 桶）索引求交得 24-mer 精确匹配。
- `nonred_table.rs`：NonRedKmerTable——排序字符串数组 + 二分；双语义（weldmer 计数 / k-mer→组件索引）。
- `graph_from_fasta.rs`：两遍扫描。Phase1 找候选 weldmer（48-mer = 匹配 24-mer 两侧各延伸 12bp 跨 contig 拼接；跳过低熵 k-mer，熵阈 1.3）；计数阶段用 reads 回验；Phase2 动态阈值 `minCov = max(min_glue=2, ceil(max_covA,covB × glue_factor=0.05))` → shadow/包含（per_id>97）剔除 → 覆盖比检查（min_iso_ratio=0.05）→ 双向连边。输出 weld 图文本。
- `bubble_up.rs`：支持度降序贪心单链接聚类（外部 `sort -k9,9gr` 语义内建），簇容量 ≤ 25，总长 < min_contig_length 的簇丢弃，未聚类 contig 各自成簇；输出 COMPONENT 块。
- `bundle.rs`：COMPONENT 解析 → 'X' 连接的 bundle FASTA（`>s_<comp>` header 带 cov）。
- `reads_to_transcripts.rs`：每 read 全部 25-mer（熵 ≥1.5）查表 → 排序取最长 run（多数投票）→ `pct = max/num × 100 ≥ -p` 才分配；输出 `comp_id\tread_name\tNN%\tseq`。
- `quantify.rs`：QuantifyGraph 移植——read 25-mer（K+1）排序数组二分命中收集（含反向命中坐标换算）；graph 第 3 列改写为跨边 read 计数；输出 `graph.reads`（`name\tpos1\tnode1\tpos2\tnode2\tseq\t+|-`）。
- 分区器：Cbin 目录（每 1000 组件）、reads 分发、`component_base_listing.txt`。

### 5.5 trinity-butterfly（~8k 行，最难）

上帝类按方法组拆分：
```
context.rs          显式上下文,收编全部静态全局: nodeTracker / origIDnodeTracker /
                    NUM_MISMATCHES_HASH / MAX_MM_ALLOWED / LAST_ID
graph.rs            SeqVertex(_prevVerticesID 合并历史/逐碱基权重/节点深度)、SimpleEdge、
                    自研有向图(O(1) findEdge、禁平行边——镜像 JUNG DirectedSparseGraph 语义;
                    虚拟节点 ROOT=-1 / T_VERTEX=-2)
graph_io.rs         graph.out(5列) / graph.reads 解析
pruning.rs          fixExtremelyHighSingleEdges → removeLightEdges(姐妹边/流量阈值 0.02)
                    → compactLinearPaths → removeSingleNtBubbles(SNP 气泡)
dfs.rs              My_DFS(四色标记/拓扑序/finish time/节点深度)
align.rs            NW-Gotoh / SW-Gotoh / banded-NW 移植(f32 打分,对齐 jaligner;
                    match=4 mismatch=−5 gap_open=10 gap_extend=1)、ZipperAlignment、
                    AlignmentStats、比对缓存("path1;path2" 字典序键)
read_threading.rs   updatePathRecursively: 递归+记忆化(键 "{vertexID}_{locInNode}_{locInSeq}")
                    + zipper → 100bp 短 NW 预检 → banded NW 三级切换
                    + 端部 gap 归属修正 + 平局保留第一条
pair_paths.rs       PairPath 兼容性判定族(isCompatible / isCompatibleAndContainedBySinglePath /
                    containsSubPath / node_is_contained_or_possibly_in_gap)
pog.rs              read-path 重叠图(POG)构建 → 破环 → 转回 SeqVertex DAG → PairPath 重映射
paths.rs            getAllProbablePaths: 深度优先队列 + 父节点延迟 + triplet/扩展 triplet 锁定
                    + 四档 read 支持判定(original > compatible(默认) > lenient > all_possible)
                    + 裸边保底播种 + MAX_NUM_PATHS_PER_NODE(100/25) 上限
postprocess.rs      CD-HIT 式去冗余(twoPathsAreTooSimilar: findLastSharedNode 三段分解 +
                    图首/图尾/内部三种 gap 归属规则) → 基因分组(_g<i>_i<j>, 30% 长度重叠)
                    → EM 降维(可选,默认关闭,对齐 --NO_EM_REDUCE)
```
- 路径序列拼接语义：首节点全序列，后续节点去前 K−1 碱基（`getNameKmerAdj`）；orig-id 回溯（`_prevVerticesID` 嵌套结构）精确移植。
- 路径强化回看距离 = `PATH_REINFORCEMENT_DISTANCE_PERCENT(25%) × MAX_PAIR_DISTANCE`；`MIN_READ_SUPPORT_THR` 默认 1。
- 输出：`>{comp}_g{gene}_i{iso} len={len} path=[...]`，60 列折行。

### 5.6 trinity-cli（~1.5k 行）
- 参数解析（原版同名参数）；prep_seqs（SS 时 R 端 revcomp 规整、FQ→FA）；
- 调度 DigiNorm → k-mer → inchworm → chrysalis → butterfly（组件命令生成 + 线程池）→ 汇总；
- `.ok` checkpoint 生命周期管理与断点续跑；
- stderr 进度输出贴近原版格式；`-V` 透传给 Butterfly；
- 最终输出 `Trinity.fasta` + `Trinity.fasta.gene_trans_map`。

## 6. Butterfly 正确性保障（最大风险点）

**分层移植 + 图状态快照对比**。原版 `--generate_intermediate_dot_files` 提供 8 个检查点：`_deBruijn.A → _pruned.B → _compacted.C → _SNPs.H → _POG → _POG.cyclesRemoved.rN → _vertexDAG.I → _final.Z`。Rust 同样实现各检查点导出，逐阶段 diff 图状态（节点/边集合语义），把 15k 行黑盒验证拆成 8 个可独立对比的中间态。

移植顺序按依赖链：`graph_io → pruning → dfs → read_threading → pair_paths → pog → paths → postprocess`，每层先过原版自带 `sample_data/c1.graph` 检查点，再上大数据。

四个已识别陷阱与对策：

| 陷阱 | 对策 |
|---|---|
| 全局可变状态（含 `MAX_MM_ALLOWED` 每条 read 隐式改写的重入语义） | 全部收进 `BflyContext` 显式传递；隐式语义原样保留并注释 |
| jaligner 用 float(f32) 打分，路径等价判定依赖浮点精确相等 | Rust 同用 `f32`；缓存键字典序语义原样移植 |
| `updatePathRecursively` 600 行递归（Java 靠 `-Xss`） | 保留递归，组件线程 `stack_size` 可配（默认 256MB）；实测溢出则转显式栈（唯一允许的结构性偏离，需注释说明） |
| `getPrevCalcNumMismatches` 三段分解 + 三种 gap 归属规则 | 逐分支移植 + 每分支最小复现用例，测试密度加倍 |

## 7. 确定性策略

- 单线程模式完全可复现（平局打破固定种子 RNG）——比原版更严格，利于回归；
- 多线程模式镜像原版语义（并行冗余用同样去重机制收敛）；
- 组件级并行天然安全（组件独立），是 Butterfly 主要加速来源。

## 8. 错误处理与健壮性

- 库层：每 crate `thiserror` 错误类型，`Result` 贯穿；CLI 层：`anyhow` + 上下文链，错误消息贴近原版（`Error, ...` 前缀）；
- 输入校验：FASTQ 逐条校验（长度一致、合法字符），报错带 read 名与行号；
- 内存护栏：k-mer 表按 `--max_memory` 预估（复刻 `(jellyfish_ram − file_size)/7`），超限清晰报错不 OOM；
- 组件级故障隔离：Butterfly 线程池单组件 panic `catch_unwind` 捕获，标记失败、其余继续，结束非零退出并列出失败组件（`.ok` checkpoint 支持重跑续传）；
- 日志：stderr 贴近原版格式。

## 9. 性能（正确性达标后）

- FASTQ：memmap + 字节扫描零拷贝；
- k-mer 计数：分片开放寻址（~12B/kmer）、rayon 按片并行无锁竞争；
- Butterfly：组件线程池 + 比对缓存 + 路径包含缓存；
- 基准流程：`sample_data`（小）→ `trinity_ext_sample_data`（真实规模），对比原版 wall-clock + 峰值 RSS（`/usr/bin/time -v`），结果入 `docs/benchmarks.md`。性能指标为跟踪项，非硬验收条件。

## 10. 验证体系（三层）

**第 1 层 单元黄金向量**：从原版行为提取（编码/revcomp/熵/平局/比对打分样例），`xtask gen-fixtures` 用原版二进制批量生成固化到 `fixtures/`，单测断言。

**第 2 层 交叉验证（`xtask xcheck`）**：
```
xcheck inchworm | chrysalis | butterfly
```
- 同输入双跑（原版 vs Rust）→ 格式感知比较器：`graph.out` 逐行精确；inchworm contigs 多重集；`readsToComponents` 按 (comp, read) 集合；DOT 图按节点/边集合；
- 双向互喂：Rust 输出 → 原版下一阶段二进制；原版输出 → Rust 下一阶段；
- Butterfly 用 8 检查点 DOT diff；
- 前提：本机编译原版 Trinity（`docs/setup.md` 记录步骤）。

**第 3 层 端到端**：`sample_data`（SP + PE）→ `trinity_ext_sample_data`（真实规模）；`xtask eval` 报告：序列多重集比对（精确匹配 / ≥99% 一致性聚类 / 仅原版有 / 仅 Rust 有）+ 基因数统计对照。差异量化报告而非硬失败。

## 11. 实施阶段（每阶段有验证门）

| 阶段 | 内容 | 验证门 |
|---|---|---|
| P0 | workspace 脚手架 + trinity-common + 编译原版 Trinity | 黄金向量单测全绿；原版二进制可用 |
| P1 | trinity-kmer（计数/dump/DigiNorm） | dump 可喂原版 inchworm；DigiNorm 保留率与原版同数量级（±2%） |
| P2 | trinity-inchworm | 同输入 contig 多重集高一致（单线程模式）；输出可喂原版 Chrysalis 跑通 |
| P3 | trinity-chrysalis（6 子阶段） | 每子阶段输出 vs 原版逐个比对通过 |
| P4 | trinity-butterfly（按依赖链分层） | c1.graph 检查点全过 → sample_data 全组件通过 |
| P5 | trinity-cli + 端到端 + 断点续跑 + 性能 + 文档 | 三层验证全绿；基准入档 |

工作量重心：P4 ≈ 40%，P3 ≈ 25%。

## 12. 风险与对策

| 风险 | 对策 |
|---|---|
| POG/破环逻辑复杂 | 检查点 diff 把问题隔离在单层 |
| 深数据集 k-mer 内存 | `--max_memory` 护栏 + 内存对照文档 |
| 两实现固有不确定性差异 | 比较全部用多重集/集合语义；评估量化而非硬失败 |
| 单阶段卡死 | 各阶段独立 bin，可手工喂数据调试 |
| 递归栈溢出 | 组件线程大栈可配；兜底转显式栈 |

## 13. 文档

`README.md`（用法）、`docs/porting-map.md`（Rust↔原版模块/函数/行号映射表，持续维护）、`docs/setup.md`（含原版编译步骤）、`docs/benchmarks.md`。

## 附录 A：原版关键事实速查（复刻必须遵守的常量）

- k-mer 编码：G=0, A=1, T=2, C=3；互补 = 按位取反；K 上限 32；主线 K=25（Chrysalis 传递 K−1=24 为节点长，25 为 read 匹配单位；weldmer kk=48=2×24）
- jellyfish：DS（无 SS_lib_type）加 `--canonical`；dump 格式 header 为裸计数；`dump -L $min_kmer_cov`
- Inchworm 默认：K=25、min_kmer_count=1、MIN_SEED_ENTROPY=1.5、MIN_SEED_COVERAGE=2、MIN_ASSEMBLY_COVERAGE=2、min_ratio_non_error=0.005、MAX_RECURSION=1、硬停 50；输出 header `>a<N>;<avg_cov> total_counts: <tc> Seed: <seed_kmer_count> K: <K> length: <len>`，60 列折行；主流程默认 `--no_prune_error_kmers`
- GraphFromFasta：min_glue=2、glue_factor=0.05、min_iso_ratio=0.05、熵阈 1.3、per_id 包含阈 97
- BubbleUpClustering：max_cluster_size=25、min_contig_length=200（bundle min=200）
- ReadsToTranscripts：k=25、min_kmer_entropy=1.5、pct 阈 -p=50（参数默认 10，主流程钳制到 50）、`-max_mem_reads` 默认 50M
- QuantifyGraph：k=24（用 K+1=25 匹配）、max_reads=200000
- Butterfly（Trinity 主流程调用）：`-N 100000 -F pair_dist`、`--path_reinforcement_distance`、`--NO_EM_REDUCE`；EDGE_THR=FLOW_THR=0.02；比对 4/−5/10/1；MAX_READ_SEQ_DIVERGENCE=0.05；max_paths_per_node 100/25；min_per_id_same_path=98%
- 归一化：K=25、jellyfish `dump -L 2`（MIN_KMER_COV_CONST=2 硬编码，是**计数表**下限）；nbkc 选择用 `min_cov = min_kmer_cov`（Trinity 默认 1）、`max_cov=200`（Trinity L214）、`max_CV=10000`（事实关闭）；接受概率 `max_cov/med_cov`，过滤 `med_cov<min_cov`、`mean≤0`、`CV>max_CV`；PE 默认 together 模式（`--pairs_together --PARALLEL_STATS`），stats 合并取左右 `sprintf("%.1f",(l+r)/2)`；stats 由 `Inchworm/bin/fastaToKmerCoverageStats` 生成（median 偶数个为**整数截断除法**、mean/stdev 全程 f32、stdev 除 n-1，n=1→NaN 可通过 CV 过滤、缺失/含 N kmer 一律按 1）
- 输出目录布局与 `.ok` checkpoint 命名与原版一致（`chrysalis/Component_bins/Cbin<i>/c<id>.graph.out` 等）

---

## 实施总结（2026-08-19，P0-P5 全部完成）

| 阶段 | 内容 | 验证等级 |
|---|---|---|
| P0 | workspace 脚手架 + trinity-common + 原版 oracle 编译 | 黄金向量单测全绿 |
| P1 | trinity-kmer（jellyfish 计数/dump/DigiNorm） | dump/coverage-stats/diginorm 端到端**逐字节**一致（xcheck-kmer 3/3） |
| P2 | trinity-inchworm | 同种子序重放（--monitor 2 抓取）三组 stdout **BYTE-MATCH**；PARALLEL/--reads 多重集对拍 4/4 |
| P3 | trinity-chrysalis 六子阶段 | 全链对拍：bundle/partition 逐字节、GFF/RTT/F2DB 多重集、quantify 刀口边白名单（7/7） |
| P4 | trinity-butterfly | 组件级 c0/c1/c2 × em/noem 对拍发布版 jar：3 PASS + 3 PASS-WARN（黄金序列全覆盖，仅多 1 条保守方向差异） |
| P5 | trinity-cli 主控 + 端到端 + 断点续跑 | 三层验证全绿；xcheck-trinity 双向覆盖率 57.4%/70.7%（阈值 50%）PASS；断点续跑演示入档 benchmarks.md |

- **已知差异带**（均为 tie-break/迭代序类，非移植缺口）：inchworm 并行种子平局序（原版固有竞态）；QuantifyGraph 刀口边 strncpy 越界读堆尾巴（真 UB，本库取 NUL 补位确定语义）；Butterfly c2 裸边播种迭代序（JUNG HashSet 后继序）→ 各多 1 条保守方向短异构体；全管线 header 内图节点编号随并行布局浮动（序列多重集稳定）。
- **backlog**：分片开放寻址 / 零拷贝流式 / diginorm 管道化等性能项已按 T4 实测裁定暂缓，依据与 revisit 触发条件见 `docs/backlog.md`。
- 验证资产与逐函数行号映射见 `docs/porting-map.md`；基准见 `docs/benchmarks.md`；端到端判定校准见 `docs/xcheck-trinity-report.md`。
