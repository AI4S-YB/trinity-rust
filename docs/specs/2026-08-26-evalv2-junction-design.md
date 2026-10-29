# 评估 v2 — junction 级双维度（StringTie 参照）设计文档

日期：2026-08-26
状态：已执行完毕（2026-08-26）——四配置+dev 双指标全景入档 opt-log, 区分力/稳定性验收 PASS

## 核心

以 **StringTie 组装 GTF**（样本自身转录本结构模型）为参照，用**内含子链
（junction chain）精确比对**评估 de novo 组装的两个正交维度：
**完整度**（样本感知基因级召回 Jgene）与**准确度**（junction 级精确率
Jprec）。旧指标（R/P/F1/P_full）全部保留，历史连续。

## 背景与动机

现行 F1 体系的实测盲点（2026-08-25 讨论入档）：R_full 实为注释符合度
（单样本 vs 全注释, 真实数据 3.5% 天花板）; P 碱基级对结构无感（嵌合/旁系
拿部分分）; F1 混合两种粒度。用户裁定按"准不准 + 装得完整"两维度重构,
准确度取 junction 级（重）, 参照用 StringTie 组装结果而非基因组注释——
StringTie 模型天然只含本样本表达转录本, 样本感知分母自动解决。

## 已确认决策

| 决策点 | 选择 |
|---|---|
| 完整度分母 | 样本感知基因级（StringTie 基因, 表达集天然） |
| 准确度准则 | junction 链级（内含子链与参考模型精确一致） |
| 参照来源 | StringTie 组装 GTF（dev 现成 / real 一次性新跑） |
| 实现路径 | A: 自研判定规则（minimap2 splice + GTF 解析, 工具全现成） |
| 旧指标 | 全部保留, junction 指标 opt-in（--junc-ref/--genome） |

## 数据资产（已核实）

- dev 参照: `athal/stringtie.gtf` + `stringtie.fa`（STRG.* 参考之源）
- 基因组: `athal/athal5.fa` + HISAT2 索引 `athal5.*.ht2`
- 工具: hisat2（`rnaseq` 环境）/ stringtie（`stringtie-test` 环境）/
  minimap2（cactus 环境, splice 模式）
- dev 管线原始脚本: `athal/run_orig.sh` / `hd_pipeline.sh`（StringTie
  调用参数的镜像来源）
- real reads: `bench/real/{R1,R2}.fq.gz`（10M 对, DRR095137）

## 指标定义

### junction 链提取

- 参考侧: GTF exon 行按 transcript 分组、坐标排序 → 内含子链
  [(start,end)...]（相邻外显子间隔区间）; 按 gene_id 分组（STRG.G 聚类）
- 组装侧: `minimap2 -x splice --secondary=no <genome> <assembly.fa>` →
  PAF cg 标签 `N` 操作 = 内含子 → 每组装本内含子链 + 链方向

### 完整度（Jgene / Jtx）

- 基因找回: 该基因 ≥1 条参考转录本的内含子链被某组装本**精确复现**
  （内含子集合一致且同链）; 单外显子参考转录本回退（被某组装本覆盖
  ≥70% 且 id ≥95%, asm5 语义）
- **Jgene = 找回基因数 / GTF 基因总数**; 副报 Jtx = 精确复现的参考
  转录本数 / 参考转录本总数

### 准确度（Jprec / novelJ / single）

- 多外显子组装本: 内含子链与某参考转录本完全一致 → junction-correct;
  含任何未知内含子 → 不正确
- **Jprec = junction-correct / 多外显子组装本总数**（单外显子本不进分母）
- novelJ = 含 ≥1 novel 内含子的组装本比例（单列, 不隐含对错）
- single = 单外显子组装本比例（单列）

## xtask bench 落点

- 新旗标: `--junc-ref <gtf>` + `--genome <fa>`（两者需同时提供）
- 提供时计算并输出: stdout 行 `Jgene=..% Jtx=..% Jprec=..% novelJ=..%
  single=..%`; metrics.json 增 `j_gene/j_tx/j_prec/novel_j/single_exon`
- 不提供时行为与现状完全一致; 缓存命中的跑同样计算（回溯重评零组装成本）

## real 参照构建（一次性 `scripts/prepare_stringtie_real.sh`）

```
real/{R1,R2}.fq.gz → hisat2 -x athal5 -p 8 → samtools sort → real/hisat2.bam
  → stringtie（参数镜像 hd_pipeline.sh 原调用）→ real/real_stringtie.gtf
```

幂等; GTF 基因/转录本数入档。hisat2 索引与 athal5.fa 匹配性先以单条测试
read 比对验证。

## 验收（评估体系自身的 Gate）

1. **区分力**: junction 指标复现已知排序——real 上 M5 > classic（F1 体系
   同向）; dev1 上 pogflow ≥ classic
2. **稳定性**: 同输入重放两次指标逐位一致（纯解析+集合比较）
3. **单测**: 判定规则全覆盖（链提取/精确复现/novel/基因分组/单外显子
   回退/反链/STRG id 解析, 合成用手算用例）
4. 入档 opt-log「评估 v2」段: 新旧指标同表对照四套 real 配置 + 关键 dev
   配置, 背离点注明

## 风险

| 风险 | 对策 |
|---|---|
| StringTie 参数口径漂移 | 实现第一步读 run_orig.sh/hd_pipeline.sh 原样镜像 |
| minimap2 splice 错配 | --secondary=no + 链方向校验; dev 轨道 sanity（Jgene 应显著高于 3.5% 量级） |
| hisat2 索引/fa 不匹配 | 构建脚本先单 read 验证 |
| novel junction 误伤真新剪接 | v1 定位明确: Jprec 度量与 StringTie 模型的一致度, novel 单列不判对错 |
