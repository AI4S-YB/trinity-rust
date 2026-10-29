# M6 蛋白引导组装 — 设计文档

日期：2026-08-24
状态：已与用户确认

## 核心

已知蛋白（UniProt 风格库）作模板：第一遍正常枚举的初步转录本整条翻译后
DIAMOND 比对选模板 → 算 N/M/C 三段覆盖缺口 → 第二遍从路径端点在组件图上
**引导式延长**（穿越低 read 支持区，上限为模板剩余跨度）→ 延长产物重新比对
验收，分不涨则保留原样。Pfam 域架构语法（域顺序 + 位置偏好）作为延长打分的
合法性项，probe 阶段可缺省降级为纯覆盖轮廓。

单组件级别机制上不可能比基线差（延长失败 = 保留初步转录本）。

## 已确认决策

| 决策点 | 选择 |
|---|---|
| 插入点 | 枚举内引导延长（非事后重打分） |
| 域建模粒度 | 真域架构语法（Pfam），probe 期可降级 |
| 引导库 | 现库（含拟南芥）迭代 + 近缘种（A. lyrata）公平验收 |
| Gate | 先探后定：dev1 现库 probe 看信号再定数字 |
| 实现路径 | B：两遍模板对照 + 端点引导延长 |

## 流程

```
PogFlow 组件级初步转录本（现状枚举，零改动）
  → 6-frame 最长 ORF 翻译（新写小型翻译器，无外部依赖）
  → 全组件翻译产物合并 → 一次 DIAMOND blastx（批量单调用）
  → 每组件最佳模板（coverage × identity 最高）
  → 缺口轮廓: 模板 N/M/C 三段覆盖 vs 初步转录本覆盖
  → 显著缺口组件: 重载组件图（graph_io）→ 从路径端点引导延长
      候选边得分 = 蛋白覆盖增益 × read 支持（门槛放宽至 ≥1 对，默认门槛更高）× 域合法性
      终止: 模板两端已覆盖 / 无正增益扩展 / 达模板剩余跨度上限
  → 延长产物重新 DIAMOND: 综合分须改善，否则保留原转录本
  → harvest（不变）
```

## 模块

- `crates/trinity-cli/src/protein_guided.rs` — 编排：翻译调度、DIAMOND
  外部进程（env `TRINITY_DIAMOND` 覆盖，默认 diamond2go 环境路径）、模板
  选择、缺口轮廓、验收、断点 `.protein_guided.ok`
- `crates/trinity-butterfly/src/protein_extend.rs` — 组件图端点的支持感知
  延长走图（复用 `pair_paths.rs` 内部 API）
- `xtask pfam-prep` — 离线域语法表生成（不进运行时热路径）

## 域语法表（xtask pfam-prep，版本化 tsv）

- 蛋白侧注释（protein.fa 为裸 UniProt AC，无现成注释）二选一：
  (a) 下载 EBI 预计算注释按 AC 过滤；(b) Pfam-A.hmm + busco 环境 HMMER
  对 43,823 条蛋白 hmmscan（一次性，8 线程估 1–3 小时）
- 表内容：① 每域相对位置统计 → N 偏好/中间/C 偏好分类；② 域邻接条件
  概率 P(下一域 | 当前域)；③ 每蛋白架构 AC → [(pfam_acc, start, end)]
- 运行时：延长产物隐含架构按「顺序合法性 + 模板域集覆盖 + N/M/C 轮廓
  匹配」打分；表查找顺序: `--domain-table` > env `TRINITY_DOMAIN_TABLE` >
  均无则降级为纯覆盖轮廓 + 警告

## CLI

`--protein-guided`（触发，默认 off）；`--protein-db <fa>`（蛋白库，必填
当旗标开）；`--domain-table <tsv>`（可选，缺省走 env `TRINITY_DOMAIN_TABLE`）

## 错误处理

| 情形 | 行为 |
|---|---|
| DIAMOND 缺失/失败 | 硬错误（opt-in 旗标，同 minimap2 惯例） |
| `--protein-guided` 开而 `--protein-db` 缺/文件不存在 | 硬错误（args 校验） |
| 组件翻译无蛋白命中 | 原样通过（预期多数） |
| 延长验收不通过 | 保留初步转录本 |
| 域表缺失 | 降级覆盖轮廓模式 + 警告 |
| 组件图重载失败 | 跳过该组件延长，警告不炸 run |

## 与现有特性关系

- classic / PogFlow / recursive 现有路径零改动（旗标 off 时无新代码路径）
- M5 `--model-guided` 的 mini 组装内**不**嵌套蛋白引导（正交，组合留
  M6c）

## 测试

- 单测：翻译器（正/反链、无 ORF、边界）；缺口轮廓；域顺序打分；
  合成图延长走图（带分支小图，低支持边在模板缺口指示下被穿越）
- e2e `tests/protein_guided_e2e.rs`：合成基因中间窗口仅 1 对 reads、
  其余平铺 → vanilla 在低覆盖区截断；蛋白库含真蛋白 → 断言旗标开恢复
  全长、旗标关零变化、DIAMOND 缺失 skip（仿 model_guided_e2e 惯例）
- probe：dev1 现库 `m6-probe` → R/P/F1 + 延长统计（延长组件数、平均
  延长长度、验收率）入 opt-log，作为 Gate 定标依据

## 里程碑

1. **M6a 探路**：翻译器 + DIAMOND 接线 + 模板选择 + 覆盖轮廓延长 +
   合成 e2e + dev1 现库 probe → 定 Gate
2. **M6b 验收**：Pfam 域语法表 + 近缘种库（A. lyrata 下载）+ 四集验收
3. **M6c 可选**：与 M5 组合

## 风险

- Pfam-A.hmm / 近缘种蛋白组 / EBI 注释均需网络下载，本机网络可用性
  未验证——M6b 前须确认；不可达则域语法降级、验收改用现库剔除方案
  （裸 AC 需先下载 idmapping 做物种过滤）
- 延长穿越低支持区可能引入错误序列掉 P——验收门槛 + 模板跨度上限 +
  probe 统计三重控制
