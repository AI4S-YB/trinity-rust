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
| `--KMER_SIZE`（别名 `__KMER_SIZE`，默认 25） | inchworm k-mer 大小；单独给出（无 `--multi-k`/`--no-multi-k`）= 隐式单 k 模式并打印提示警告 |
| `--multi-k auto\|k1,k2,...`（扩展，**默认 auto**） | rnaSPAdes 式多 k 迭代 + trusted contig 注入；`auto` 按读长公式选 k（RL=76 → [29,37]），显式列表须严格递升且 k∈[15,63] |
| `--no-multi-k`（扩展） | 回退 classic 单 k（k = `--KMER_SIZE` 值）；与显式 `--multi-k` 互斥 |
| `--min_kmer_cov`（别名 `--min_kmer_count`） | k-mer 计数下限 |
| `--normalize_max_read_cov` / `--no_normalize_reads` | 归一化上限 / 跳过归一化 |
| `--bfly_stack_mb` | Butterfly 线程栈（本实现扩展，默认对齐原版 JVM 栈行为） |
| `--path-engine classic\|gtasm\|rnaspades`（扩展，默认 classic） | butterfly 组件的路径引擎；`gtasm` = GT 边打分（GTW1 权重）× GTasm 式解码 × 原始 reads 边支持（只在 k=25 图上验证过：未显式 `--multi-k` 时强制单 k 并警告，显式 `--multi-k` 可强制组合） |
| `--gt-weights <path>`（扩展） | gtasm 引擎的 GTW1 模型权重（`--path-engine gtasm` 时必填） |
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

**扩展**：`--bfly_stack_mb` 显式控制 butterfly 工作线程栈上限（原版由 JVM `-Xss` 隐式决定）。`--path-engine gtasm` + `--gt-weights` 启用 GT 打分解图引擎（python GTasm 原型产品化；chrysalis 之后多一个原始 reads 边支持 pass，产物 `gtasm_edge_support.bin` 随 scratch 落盘并受 `.gtasm_support.ok` checkpoint 保护；gtasm 与 `--recursive` 互斥）。`--multi-k auto|k1,k2,...` 启用 rnaSPAdes 式多 k 迭代 + trusted contig 注入（中间轮 `iter_k<k>/` 只构图、最终轮完整管线，三数据集交叉验证为净增益）——**自 2026-08-30 起为默认行为（缺省 = `auto`）**；`--no-multi-k` 回退 classic 单 k；单独 `--KMER_SIZE` = 隐式单 k + 提示警告；显式 `--multi-k` 与 `--KMER_SIZE`/`--no-multi-k`/`--recursive` 互斥。**gtasm 特例**：`--path-engine gtasm` 未显式给 `--multi-k` 时强制单 k 并警告（gtasm 只在 k=25 图上验证过），显式给则放行并警告组合未验证；`--path-engine rnaspades` 为 rnaSPAdes MultiExtender 规则复刻引擎（高精度低召回，默认多 k 下可正常组合）。管线启动时打印生效 k 序列（`multi-k: auto → [29, 37]` / `single-k: 25`）便于确认模式。

**中间产物布局（性能优化，既定差异）**：`chrysalis/Component_bins/Cbin*/` 下不再逐组件写 `cN.graph.out` / `cN.graph.reads` / `cN.graph.allProbPaths.fasta` 小文件（2.8 万组件 ≈ 8.5 万个文件，NFS 上元数据往返极贵），改为每个 Cbin 两个打包文件 `graphs.pack`/`graphs.idx`（QuantifyGraph 产物）与 `butterfly.pack`/`butterfly.idx`（Butterfly 产物）；`.idx` 为定长 32B 记录（component id / kind / offset / len，小端）。归一化的 `left.norm.*` 便捷名是长名文件的硬链接（原版为符号链接），不再物理双写。

**归一化产物格式（fq 输入时）**：fq 输入直接产出归一化 **fa**（`left.norm.fa` / `{名}.normalized_K…​.fa`），跳过 `norm.fq` 的全量物化与下游 fq→fa 往返——归一化内部本就要做 fq→fa，选中记录直接从内存 fa 抽取，与"写 fq 再转 fa"逐字节等价（有测试锁定）。原版会把归一化 fq 留给用户；本实现**默认不再保留 norm.fq**，需要时用 `--no_cleanup`（此时 fq 与 fa 都写）。

**已知 tie-break 差异（既定差异带）**：

- inchworm 并行种子平局序与原版 `--PARALLEL_IWORM` 一样非确定——同实现两次跑的转录本多重集也会漂移（实测覆盖率带 78~94%）；跨实现端到端双向覆盖率实测 57~71%（见 `docs/xcheck-trinity-report.md` 阈值校准说明）。差异主体是端部截断的同序列（100% 一致的包含关系），非移植错位。
- 逐字节比对不适用于全管线输出；等价性按"序列多重集（含 revcomp 归一）+ ≥99% 聚类"判定（`cargo xtask eval-trinity`）。
- 输出遍历序（哈希序）不保证与原版一致，等价性均按多重集判定。

## 性能摘要

详见 `docs/benchmarks.md`（分阶段对拍）与 `docs/opt-log.md`（端到端优化记录）。

**相较原版 Trinity v2.15.2（C++/Perl/Java 工具链）：**

- k-mer 计数+dump：比 jellyfish count+dump 快 ~3.2×、峰值内存 ~1/7，输出多重集逐字节等价（sample_data PE 100k reads，4 线程）。
- inchworm：约 2× 快，峰值内存低 ~14%（单线程）。
- chrysalis 全链：峰值内存 ~1/3.2（147 MB vs 466 MB）。
- butterfly：小组件快近一个量级、峰值内存低 ~27×（无 JVM 启动与基线开销）。
- 端到端 sample_data 全量（30575 PE reads，8 线程）：本实现 ~9 s vs 原版（含 Perl/JVM）~26 s。

**相较本项目的 v0.0.1 移植版（当前版本的端到端速度优化）：**

- 拟南芥真实数据（PE，`--CPU 8`）端到端 wall time **54 min → 17 min（-69%）**，组装指标与基线一致。
- 四项优化：k-mer dump 落盘消除（二进制直桥，省 ~8 min）、ReadsToTranscripts 查询结构 16.3×、归一化 3.6×、组件分区写出 5.5×；全部带 oracle/字节级等价验证。

## 验证

- `cargo test --workspace`：553 项单测/组件对拍全绿。
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
- `docs/superpowers/plans/` / `docs/superpowers/specs/`——阶段实施计划与设计文档（P0–P5 移植、M0/M1 基准与递归分区、M7 速度优化）
