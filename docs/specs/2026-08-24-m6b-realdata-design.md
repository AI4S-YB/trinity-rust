# M6b 真实数据基准 — 设计文档

日期：2026-08-24
状态：已执行完毕（2026-08-25）——Gate G1/G2/G3 全 PASS, 详见 opt-log「M6b 真实数据基准」段

## 核心

拟南芥真实 RNA-seq（抽样 10M PE 对）上全套引擎重验：classic / PogFlow /
M5 / M6a 四配置各跑一次 xtask bench，验证 opt-log 结论体系（引擎排序、
M5 收益、M6a 中性）能否泛化到真实数据。xtask bench 沉淀**蛋白级全长回收**
（P_full/P70）为常设指标。M6a 的引导库换 *A. lyrata* 近缘种蛋白组——
真实场景、与被评物种零重叠（天然公平，无循环）。

## 背景与动机

M6a 探路四轮收敛为中性偏负（F1 15.1~15.3 ≈ 基线，opt-log 已入档）;
M6b 原定方向（域语法 + CDS 完整性度量）立项前探索证伪——dev1 基线蛋白
回收已 68.5%，70-95% 带缺口 96% 的序列在原始 reads 中不存在（数据上限，
非机制问题）。故 M6b 转向：换数据而非换机制。

## 已确认决策

| 决策点 | 选择 |
|---|---|
| 数据集 | 拟南芥真实 RNA-seq（ENA，PE Illumina，≥20M 对，选样后记录 accession） |
| 规模 | 确定性抽样 10M PE 对（固定 seed） |
| 范围 | 全套引擎重验（classic / PogFlow / M5 / M6a 四配置） |
| Gate | 泛化验证型（排序复现，不做硬提升） |
| 实现路径 | A：xtask bench 增强蛋白回收指标 + 复用现有框架，主工程零改动 |

## 数据准备（一次性脚本，落仓库 `scripts/prepare_realdata.sh`）

产物落 `/storage/home/senior007/test/athal/bench/real/`，幂等（存在即跳过）：

| 产物 | 来源 |
|---|---|
| `R1.fq.gz` / `R2.fq.gz`（各 10M） | ENA 选定数据集，确定性抽样（seed 固定） |
| `atha_cdna.fa` | Ensembl Plants（TAIR10 cDNA）——转录本参考 |
| `atha_pep.fa` | Ensembl Plants（TAIR10 pep.all）——蛋白参考 |
| `alyr_pep.fa` | Ensembl Plants（*A. lyrata* pep.all）——M6a 引导库 |

选样标准：paired-end、Illumina、≥20M 对、公开发表、metadata 无 adapter/rRNA
污染警示；经 ENA portal API 查询选定，accession 与链接入档 opt-log。

## xtask bench 蛋白回收指标

- 新旗标 `--prot-ref <fa>`；缺省时将 `--ref` 转录本最长 ORF 翻译生成
  （翻译复用 `trinity-butterfly::protein_extend::translate_codon`，xtask
  为 workspace 成员可直接依赖）
- 组装产物 DIAMOND blastx vs 蛋白参考（env `TRINITY_DIAMOND`，默认
  diamond2go 路径；makedb 按参考内容缓存）
- 指标：同一 (query, subject) 对的 HSP 求 subject 并集跨度 → coverage =
  并集 / slen，identity = 该对最大 pident; 每参考蛋白取所有 query 中
  (coverage × identity) 最高者;
  `P_full` = ≥95%cov & ≥95%id 比例，`P70` = ≥70%cov & ≥90%id 比例
- 输出进 stdout 汇总行与 metrics.json（新增 `p_full` / `p70` 字段）

## 运行矩阵

四配置统一 `--max_memory 16G --CPU 8`：

| tag | 参数 |
|---|---|
| `real-classic` | （默认引擎） |
| `real-pogflow` | `--path-engine pogflow` |
| `real-m5` | `--path-engine pogflow --model-guided` |
| `real-m6a` | `--path-engine pogflow --protein-guided --protein-db real/alyr_pep.fa` |

## Gate（泛化验证型）

- **G1**: `real-pogflow` F1 ≥ `real-classic`
- **G2**: `real-m5` F1 ≥ `real-classic`
- **G3**: `real-m6a` F1 ≥ `real-classic` − 0.5pp
- 副报告（不设门）：各配置 P_full/P70 对比；m6a 延长统计
  （extended/accepted/+nt）与 same-ref/diff-ref 覆盖度诊断重演——
  与 dev 探路结论背离时提供第一手证据

## 风险与对策

| 风险 | 对策 |
|---|---|
| adapter/rRNA 污染 | 选样阶段按 metadata 过滤；不自建清洗（YAGNI） |
| 组装时长（~35 分钟/配置 × 4 ≈ 2.5h） | bench 缓存/断点惯例 |
| 真实数据 R_full 绝对值低（isoform 复杂度） | 预期内；Gate 只看配置间排序 |
| ENA 下载中断 | 准备脚本幂等 |
| Ensembl release 路径变动 | 脚本内 URL 常量集中一处，失败即报错 |

## 与 M6a 的关系

M6a 已收档（中性偏负，`--protein-guided` 留 opt-in flag）。M6b 不改主工程
机制——四配置中的 real-m6a 用现状代码跑真实数据，若副报告显示真实数据上
延长统计与 dev 背离（如 accepted 率显著不同），机制迭代另立 M6c 评估。
