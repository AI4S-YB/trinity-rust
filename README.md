# trinity-rust

**[English](README_English.md)** | 中文

Trinity RNA-Seq de novo assembler（原版 [Trinity v2.15.2](https://github.com/trinityrnaseq/trinityrnaseq)）主装配管线的 Rust 移植——从 in silico read normalization 到 Butterfly 转录本输出，端到端单二进制、无 Perl/Java/C++ 外部依赖（jellyfish / seqtk / ParaFly / Butterfly.jar 均不再需要）。

## 用法

```bash
cargo build --release
target/release/trinity-cli \
  --seqType fq --left reads.left.fq.gz --right reads.right.fq.gz \
  --CPU 8 --max_memory 2G --output out
# 产出: out.Trinity.fasta + out.Trinity.fasta.gene_trans_map
```

常用参数（与原版 `Trinity` 同名）：

| 参数 | 说明 |
|---|---|
| `--seqType fq\|fa` | 输入类型（`.gz` 自动解压） |
| `--left/--right` | PE reads（逗号分隔多文件）；`--single` 为 SE |
| `--SS_lib_type F\|R\|RF\|FR` | 链特异性库（RF/FR = PE） |
| `--CPU/--max_memory` | 线程数 / 内存上限护栏 |
| `--KMER_SIZE`（别名 `__KMER_SIZE`，默认 25） | inchworm k-mer 大小 |
| `--min_kmer_cov`（别名 `--min_kmer_count`） | k-mer 计数下限 |
| `--normalize_max_read_cov` / `--no_normalize_reads` | 归一化上限 / 跳过归一化 |
| `--bfly_stack_mb` | Butterfly 线程栈（本实现扩展，默认对齐原版 JVM 栈行为） |
| `--no_cleanup` | 保留中间产物（both.fa / left.fa / right.fa、归一化 norm.fq 等；位置在 scratch 目录） |
| `--scratch <dir>` | 中间产物写到 `<scratch>/<output 目录名>/`（显式指定时不做空间检查） |

断点续跑：各阶段 `.ok` checkpoint 文件，重跑同 `--output` 自动跳过已完成阶段。

**中间产物的默认位置（本实现的既定默认，与原版不同）**：不传 `--scratch` 时，全部中间产物（归一化目录 / both.fa / mer_counts / inchworm / chrysalis / butterfly）写到本地临时目录 `${TMPDIR:-/tmp}/trinity-rust-<输出目录路径哈希>/`，`--output` 目录只保留最终产物（`<output>.Trinity.fasta` + gene_trans_map 照常落在其旁）。路径按输出目录哈希确定，同一 `--output` 重跑命中同一 scratch 目录，续跑语义不变。启动时做 statvfs 空间检查（粗估需求 = 压缩输入总大小 × 6，下限 4GB）；TMPDIR 不可写或空间不足时打印警告并回落旧行为（中间产物直接放 `--output`）。scratch 视为易失存储：删掉它等于从头重算；运行成功后该目录不自动清理，可手动删除。

`--recursive` 模式（分区重组装）也可用。

## 架构

五个 crate + 编排 CLI，数据流自上而下镜像原版：

```
                 trinity-cli（编排：参数面 / checkpoint / 逐阶段调用 / 汇总）
                                    |
   归一化（diginorm, K25 maxC200） → both.fa
                                    |
   trinity-kmer        jellyfish count/dump 等价（计数 + -L 过滤 + 多文件合并）
                                    |
   trinity-inchworm    线性 contig 构建（k-mer 种子延伸, PARALLEL 模式）
                                    |
   trinity-chrysalis   GraphFromFasta → BubbleUpClustering → ReadsToTranscripts
                       → FastaToDeBruijn → QuantifyGraph → 组件分区
                                    |
   trinity-butterfly   每组件图搜索 + EM + 路径输出（allProbPaths）
                                    |
   harvest             <out>.Trinity.fasta + gene_trans_map
                                    |
   trinity-common      2-bit k-mer 编码 / FASTA·FASTQ 读 / sdbm seq_hash /
                       drand48 / seqtk 读名改写 等底层原语
```

## 与原版的差异（既定清单）

**未移植（明确不支持）**：salmon 表达定量、bowtie（`--no_bowtie` 语义恒成立）、jaccard 剪枝读聚类、长读（`--long_read`）、DNA 模式（`--genome_guided` 之外的 DNA 组装）。以上入口直接报错或忽略并警告。

**参数别名**：`__KMER_SIZE` → `--KMER_SIZE`、`--min_kmer_cov` → `--min_kmer_count`（两个名字都收，与原版主程序名兼容）。

**扩展**：`--bfly_stack_mb` 显式控制 butterfly 工作线程栈上限（原版由 JVM `-Xss` 隐式决定）。

**中间产物布局（性能优化，既定差异）**：`chrysalis/Component_bins/Cbin*/` 下不再逐组件写 `cN.graph.out` / `cN.graph.reads` / `cN.graph.allProbPaths.fasta` 小文件（2.8 万组件 ≈ 8.5 万个文件，NFS 上元数据往返极贵），改为每个 Cbin 两个打包文件 `graphs.pack`/`graphs.idx`（QuantifyGraph 产物）与 `butterfly.pack`/`butterfly.idx`（Butterfly 产物）；`.idx` 为定长 32B 记录（component id / kind / offset / len，小端）。归一化的 `left.norm.*` 便捷名是长名文件的硬链接（原版为符号链接），不再物理双写。

**归一化产物格式（fq 输入时）**：fq 输入直接产出归一化 **fa**（`left.norm.fa` / `{名}.normalized_K…​.fa`），跳过 `norm.fq` 的全量物化与下游 fq→fa 往返——归一化内部本就要做 fq→fa，选中记录直接从内存 fa 抽取，与"写 fq 再转 fa"逐字节等价（有测试锁定）。原版会把归一化 fq 留给用户；本实现**默认不再保留 norm.fq**，需要时用 `--no_cleanup`（此时 fq 与 fa 都写）。

**已知 tie-break 差异（既定差异带）**：

- inchworm 并行种子平局序与原版 `--PARALLEL_IWORM` 一样非确定——同实现两次跑的转录本多重集也会漂移（实测覆盖率带 78~94%）；跨实现端到端双向覆盖率实测 57~71%（见 `docs/xcheck-trinity-report.md` 阈值校准说明）。差异主体是端部截断的同序列（100% 一致的包含关系），非移植错位。
- 逐字节比对不适用于全管线输出；等价性按"序列多重集（含 revcomp 归一）+ ≥99% 聚类"判定（`cargo xtask eval-trinity`）。
- 输出遍历序（哈希序）不保证与原版一致，等价性均按多重集判定。

## 性能摘要

详见 `docs/benchmarks.md`（分阶段对拍）与 `docs/opt-log.md`（端到端优化记录）。

**端到端（真实数据：拟南芥 DRR095137 前缀子集，1000 万 PE reads，~1.1 GB gz）：**

| 版本 | 8 线程 | 64 线程 | 说明 |
|---|---|---|---|
| v0.0.1（直移植版） | ~54 min | 未测 | 移植基线，指标锚定原版 |
| v0.3.0 | ~17 min | 12:29 | 算法级优化（见下） |
| **v0.4.0** | **11:20（-79%）** | **7:37** | IO/系统开销优化 |

端到端小集（原版 sample_data 全量 30575 PE reads，8 线程）：本实现 ~9 s vs 原版（含 Perl/JVM）~26 s。

**分阶段对比原版 Trinity v2.15.2（同等线程，3 次复测取代表值）：**

| 阶段 | 原版 wall / 峰值 RSS | 本实现 wall / 峰值 RSS | 加速 | 内存 |
|---|---|---|---|---|
| k-mer 计数+dump（100k reads，4 线程） | jellyfish 0.95 s / 626 MB | 0.30 s / 86 MB | ~3.2× | ~1/7 |
| inchworm（28.8 万 k-mer，单线程） | 0.72 s / 26.1 MB | 0.37 s / 22.5 MB | ~2.0× | -14% |
| chrysalis 全链（55 组件，单线程） | 11.6 s / 466 MB | 9.9 s / 147 MB | ~1.2× | ~1/3.2 |
| butterfly 小组件（c0, N=4342） | 1.15 s / 244 MB | 0.03 s / 9.2 MB | ~40× | ~1/27 |
| butterfly 中组件（c1/c2, N≈1 万） | 1.7 s / ~420 MB | 0.6 s / ~380 MB | ~2.5× | 持平 |

注：chrysalis 的 ReadsToTranscripts 阶段在移植初期多线程未生效（-t 4 一度慢于原版），v0.3.0 已重做其查询结构（16.3×）。butterfly 小组件的悬殊差距主要来自原版 JVM 启动与基线开销。

**线程扩展性（v0.4.0，真实数据，`--max_memory 64G`）：**

| 线程 | wall | 峰值 RSS | 加速比 |
|---|---|---|---|
| 8 | 11:20 | 22.7 GB | 1.00× |
| 16 | 9:12 | 28.5 GB | 1.23× |
| 32 | 8:00 | 33.7 GB | 1.42× |
| 64 | 7:37 | 37.2 GB | 1.49× |

32 线程后趋于饱和（32→64 仅 1.05×）；峰值内存随线程数上升（k-mer 计数器并发累加器效应）。本数据规模推荐 16-32 线程。

**速度优化技术分解：**

- v0.3.0（算法级）：k-mer dump 落盘消除（二进制直桥，省 ~8 min）、ReadsToTranscripts 打包哈希索引 + 滑窗熵（16.3×）、归一化并行覆盖统计 + 滚动 revcomp（3.6×）、组件分区批量写出（5.5×，原实现每行裸写 ~45M 次系统调用）。
- v0.4.0（IO/系统开销）：归一化双写消重（-2.6 GB 写）、中间产物默认落本地临时目录（`--scratch` 默认化，输出目录零中间文件）、组件数据全量打包（消灭 ~14 万 per-component 小文件，butterfly 阶段 136 s → 31 s）、归一化直出 fa、mimalloc 全局分配器（sys time -85%，页面错误 -96%）。

**IO 开销治理（v0.4.0，真实数据 64 线程实测前后对比）：**

| 指标 | v0.3.0 | v0.4.0 | 变化 |
|---|---|---|---|
| 中间产物写入量 | ~13.4 GB | ~9.6 GB | -29% |
| per-component 小文件数 | ~14 万个 | 每 Cbin 4-6 个包文件 | ≈ 消除 |
| sys time 占比 | 2056 s（31%） | 561 s（13%） | -73% |
| minor page fault | 1.97 亿 | 3540 万 | -82% |
| butterfly 阶段 wall | 136 s | 31 s | -77% |
| wall（默认即 scratch） | 12:29 | 7:23 | -41% |

背景：组装的输出目录常挂在网络盘（NFS）上，原版式的"每组件一套小文件 + 阶段间全量文件倒手"在 NFS 上每个 create/write/close 都是网络元数据往返。v0.4.0 的治理路径是三条：中间产物整体移出网络盘（`--scratch` 默认化，见上文"中间产物的默认位置"）、小文件打包成顺序读写的 pack + 定长偏移索引、以及同内容文件硬链接去重。峰值 RSS 不受影响（mimalloc 前后持平）。

所有优化均带 oracle/字节级等价验证，组装指标与移植基线一致（并行平局序不确定带见"差异"节）。完整记录见 `docs/benchmarks.md` 与 `docs/opt-log.md`。

## 验证

- `cargo test --workspace`：543 项单测/组件对拍全绿。
- 三层交叉验证（`cargo xtask`）：
  - `xcheck-kmer / xcheck-inchworm / xcheck-chrysalis / xcheck-butterfly`——第 1、2 层（单元/管线级，对拍原版二进制或黄金）；
  - `xcheck-trinity [--full]`——第 3 层端到端：截断（默认 50000 PE reads）/全量双侧全管线对拍 + eval 统计 + both.fa 互喂抽查 + SS(RF) 合成小集；判定阈值与校准依据见 `docs/xcheck-trinity-report.md`。
  - `eval-trinity <ours> <orig>`——单独出 eval 报告（精确匹配 / 99% 聚类 / 双向覆盖率 / gene 数）。

原版工具链环境（对拍用）：`TRINITY_SRC`（原版源码树）、jellyfish/java 在 conda env `trinity`，详见 `docs/setup.md`。

## 开发文档

- `docs/porting-map.md`——Rust ↔ 原版源码逐函数映射表（含行号与已证差异）
- `docs/benchmarks.md`——基准记录
- `docs/backlog.md`——积压与裁定记录
- `docs/setup.md`——环境与原版工具链准备
- `docs/plans/` / `docs/specs/`——阶段实施计划与设计文档（P0–P5 移植、M0/M1 基准与递归分区、M7 速度优化）
