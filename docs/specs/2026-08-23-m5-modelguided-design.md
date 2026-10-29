# M5 模型引导组装 — 设计文档

日期：2026-08-23
状态：已与用户确认

## 核心

reads 招募到模型转录本 → 基因级分组 → 每组定向重组装 → 合并输出。
模型可来自第一遍自组装（自引导）或外部参照（可选）。

## 流程

```
reads → [外部模型 or 第一遍 PogFlow] → 模型转录本
  → 基因聚类（≥50% 重叠 → 弱连通分量）
  → minimap2 比对 reads → 模型（identity ≥80%, best hit）
  → PE 对取两端命中并集 → 分配到基因
  → 每基因组: 迷你组装（inchworm→chrysalis→pogflow）
  → orphan 组: 同流程
  → 合并 + 去冗余 → Trinity.fasta
```

## CLI

`--model-guided`（触发）；`--model <fa>`（可选外部参照，无则自引导第一遍）

## 验收

dev1 F1 ≥ 15.5%（M3 基线 15.1% + ≥0.4pp——分组组装的简化收益）
