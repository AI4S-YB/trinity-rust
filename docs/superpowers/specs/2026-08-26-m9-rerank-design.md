# M9 候选重排（candidate reranking）— 设计文档

日期：2026-08-26
状态：已与用户确认（Transformer 讨论的收敛形态：DeepVariant 范式——引擎产候选、模型做全长级判别）

## 核心

PogFlow 放宽过滤多产候选 → 每条候选提取**全长级特征**（序列 + 覆盖度形状 + read 支持）
→ 模型判别"是否真转录本"。监督 = 参考转录本锚定（minimap2, eval 同口径）。
第一阶段只验证信号：数据集 + GBM 基线；有序列结构信号再上 Transformer（3090 可用）。

## 与 M8 探针的关系

M8 探针（转换级）证明：边/转换级的判别信号 ≈ 覆盖度，序列浅特征无增量。
M9 换了问题层级：**全长转录本级**——ORF/密码子结构、嵌合体双 ORF、重复区、
末端形态、覆盖度形状等信号只在全长层级存在，是 M8 探针原理上摸不到的。

## 数据集构建（xtask export-candidates）

- 候选来源：xtask 直接以库调用 `trinity_butterfly::run_component`
  （PogFlow 引擎, min_read_support=1 放宽, min_len=200, flow_cov=None），
  离线消费 bench 缓存的 `c*.graph.out/.graph.reads`——**管线零改动、不重跑组装**；
- 特征列（全部推理时可得）：len、support（read 支持计数）、GC、熵、
  最长 homopolymer、最长 ORF（aa 长度/占比/起始 M 有无, 3 框架正向,
  bench.rs translate_longest_orf 同语义）、覆盖度形状（沿序列 25-mer 查计数表:
  均值/中位/最小/cov≥2 窗口占比/变异系数, 计数表复用 export-edges 的装载）;
- 标签列（oracle, 仅训练/切片）: minimap2 候选→参考, best hit 的 cov/ident,
  正例 = cov≥0.95 且 ident≥0.95（eval R_full 同口径）;
- 输出: candidates.<tag>.tsv, 每候选一行。

## 评估

- GBM（HistGradientBoosting）: dev1+dev2 训 / dev3 测 / holdout 终测;
  特征重要性排列——重点看 ORF/覆盖度形状是否真有增量;
- **操作口径**（直接模拟引擎换过滤器）: 按模型分过滤候选 → 对参考算 R_full/P/F1,
  与现状过滤器（support≥2）的同口径对照——无需重跑组装;
- 判定: 模型过滤 F1 明显超 support≥2 基线（>0.5pp）且有序列特征重要性
  → 进 Transformer 阶段; GBM 已打平 → 结论"全长信号也可被浅特征吸收", 收档。

## 风险（沿用既有结论, 如实列）

1. 标签噪声: StringTie 参考漏报的真实转录本 → 假阴性（M8 已量: 高覆盖"负例"
   63% 属此类）——训练参考用 StringTie GTF 对应 fa, 必要时换 Ensembl 复评;
2. 模拟/真实域差: 训练在 dev（模拟 reads）, real 轨道留作域外验证;
3. 天花板: 若残差主要是参照不全, 完美过滤器在 F1 上也只显小增益——
   真正的奖品可能是 R（找回 support 过滤误杀的低丰度真转录本）。
