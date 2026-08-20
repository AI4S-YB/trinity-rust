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
| `--no_cleanup` | 保留中间产物（both.fa / left.fa / right.fa 等） |

断点续跑：各阶段 `.ok` checkpoint 文件，重跑同 `--output` 自动跳过已完成阶段。

## 架构

六个 crate + 编排 CLI，数据流自上而下镜像原版：

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

**已知 tie-break 差异（既定差异带）**：

- inchworm 并行种子平局序与原版 `--PARALLEL_IWORM` 一样非确定——同实现两次跑的转录本多重集也会漂移（实测覆盖率带 78~94%）；跨实现端到端双向覆盖率实测 57~71%（见 `docs/xcheck-trinity-report.md` 阈值校准说明）。差异主体是端部截断的同序列（100% 一致的包含关系），非移植错位。
- 逐字节比对不适用于全管线输出；等价性按"序列多重集（含 revcomp 归一）+ ≥99% 聚类"判定（`cargo xtask eval-trinity`）。
- 输出遍历序（哈希序）不保证与原版一致，等价性均按多重集判定。

## 性能摘要

详见 `docs/benchmarks.md`。代表数据（sample_data PE 100k reads，4 线程）：

- k-mer 计数+dump：比 jellyfish count+dump 快 ~3.2×、峰值内存 ~1/7，输出多重集逐字节等价。
- 端到端 sample_data 全量（30575 PE reads，8 线程）：本实现 ~9 s vs 原版（含 Perl/JVM）~26 s。

## 验证

- `cargo test --workspace`：509 项单测/组件对拍全绿。
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
- 阶段 spec / plans（P1–P5 工作集）不在仓库内，随会话归档；结论均已落上述文档。
