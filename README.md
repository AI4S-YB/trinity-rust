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
#        + out.Trinity.principal.fa + out.Trinity.family_map（去冗余交付链, 默认开）
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
| `--no-gtasm-visited-walk`（扩展） | 关闭 gtasm 解码的 visited 穿越（默认开：walk 可穿越已 visited 节点救回共享外显子的弱势 isoform，完全重复产出自动丢弃）；给此旗标回退 python 原型语义 |
| `--no_cleanup` | 保留中间产物（both.fa / left.fa / right.fa、归一化 norm.fq 等；位置在 scratch 目录） |
| `--no-regroup`（扩展） | 关闭输出层去冗余交付链（**默认开**：严格家族拆分 + reads 支持度代表 + principal 截断 len≥500/supp≥10，落 `*.principal.fa`/`*.family_map` 侧车，主 fasta 不动；`--principal-min-len/--principal-min-supp` 调档，显式 0/0 关截断门） |
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

**输出层去冗余交付链（regroup，2026-10-29 起默认开）**：Trinity.fasta 落盘后做基因内严格家族拆分（纯序列证据区分旁系同源/isoform：锚门控 + 多对角带状比对，iso 边贪心合并 + para 边 cannot-link，split-only 不跨基因并），reads 支持度选代表（原始 reads k-mer 投票压"长但无人读"的脏延伸），principal 截断（默认代表长 ≥500 且票数 ≥10）——`<out>.Trinity.principal.fa`（每家族一条代表）+ `<out>.Trinity.family_map`（全量转录本→家族映射）随主产物落盘，主 fasta 与 gene_trans_map 不动。QM4 BUSCO 对照：全集 D 41.0→9.5、转正档 44,522 条 C 81.3 近乎无损（详见 docs/opt-log.md 末节）。`--no-regroup` 整体回退（零副作用）；`--principal-min-len 0 --principal-min-supp 0` 只关截断门（分组+代表照旧）；`--regroup-cov` 调拆分阈值。

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

**默认开启的 QM4(150bp) 验证项（2026-09-04 转正）**：

- `--rtt-multi-kmer`（**默认开**，`--no-rtt-multi-kmer` 回退）：ReadsToTranscripts 25-mer 索引多值化——共享 kmer 记全部组件，救跨 weld 口 read 被其他组件抢占（单值后写者赢）导致的组件断口无桥。QM4（150bp 真实）：R +0.3、缺失转录本召回 +99、230 个原零 read 组件进图；dev1/dev3（76bp 模拟）：R/P 零回退。
- `--weld-min-glue`（**默认 auto**）：读长 ≥150 → 1（mg1，QM4 组合验证 R+0.6/F1+0.2），短读 → 2（原版阈值，dev1 76bp 上单 read 证据 P−1.7 系噪声）；显式值覆盖。分界常数 `multik::WELD_MG1_READ_LEN_MIN`。
- quantify 阶段随 rtt-multi-kmer 放大的耗时已由物化键排序优化压回（间接寻址 memcmp → ACGT 2bit u64 主群 + 含 N 次群，微基准 3.5×，产物逐字节一致）。
- diginorm 计数表 k≤32 键层再降一档（u128→u64，键带宽/哈希域减半）：QM4 全量微基准 602→562s（热段 −15%），产物逐字节一致。
- inchworm 轮 kmer_count 分层窄键计数（k≤32 u64 / k≤64 u128 + 终点 widen 回 KmerId，消费方接口不变）：微基准 k=49 1.38×，QM4 端到端 34:07→**29:48**（kmer_count 两轮 -241s，峰值 RSS 118→88GB），计数表亿键逐键对拍等价。
- diginorm 计数改分片共享表（Mutex 片表 + 终点精确容量折叠）：diginorm 同负载 677→516s（-24%），产物逐字节一致。当时"锁护航"裁定轮计数维持局部表+reduce——后被内存战役证实系分片缓冲病污染的误判，轮计数已统一切分片（见下方内存战役节）。
- mimalloc 全局分配器补 `PURGE_DELAY=0`（程序内默认）：分片归一化前置 + 轮计数大表翻倍在默认滞留下全量 QM4 必现 OOM（45s 内 35→455GB），置 0 后全程通过（峰值 85GB）；mimalloc 相对 glibc 切片口径快 22%。

**内存战役（2026-09-05/06 两期，QM4 峰值 83→21.1GB，目标 32→达成后压至 ~20G 档）**：

- 第一期六刀（83→28.7GB）：分片冲刷字节定界+`with_min_len`（缓冲病根修）；轮计数统一切分片（旧"锁护航回退"裁定系缓冲病污染，作废——count_bench k49 249→**87s**）；gff 改 omp atomic 共享计数 + crossover 即释 + Cow 借用 + 共享 rc_seqs（102.8→24.1GB，weld 图逐字节一致）；种子排序 `sort_unstable`；种子收集精确容量；DashMap 重建提前到空窗+种子只取键。
- 第二期六刀（28.7→21.1GB）：diginorm 侧计数串行 + 表/行集分阶段释放 + emit 流式化（span 元数据 + 原位切片流写）；**轮链窄键化**——33≤k≤64 全程 u128（`count_fasta_data_layered128` 直通 + `save_counts_bin128` 流式落盘 + prune 泛型化 + `SyncKmerCounter128`，KmerId 宽表不再物化）；diginorm read+prep 按侧串行（raw/fa 不四份并存）；gff phase1 平铺竞技场 + 窗口直指竞技场（候选 11.7→4.5GB，免二次平铺）；diginorm 折叠带 min_kmer_cov 过滤（单例错误 k-mer 过半不进表）。
- 等价性：每刀带等价锁（窄目录查询/置零逐项一致、窄链 prune 多重集一致、flat/Vec 模板表逐条一致、折叠过滤=retain_min、sharded=locals 全键对拍），产物 n/N50 全程带内（117.7-118.2K，平局非确定带宽）。
- 附带收益：diginorm 432→~290s、gff 175→~170s；`save_counts_bin` 流式。端到端 34:57→25-28min（重负载窗）。685 项测试全绿。
- `--bfly-no-path-merging` / `--bfly-empty-retry` / `--weld-reads <fa>` / `--weld-glue-factor` / `--cluster-max-size`：机制探索旗标（详见 `docs/opt-log.md`），未达转正线，保留供复现实验。

## 验证

- `cargo test --workspace`：685 项单测/组件对拍全绿。
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
