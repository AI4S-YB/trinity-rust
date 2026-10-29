# M4 K19→K25 超级 read 递进 — 设计文档

日期：2026-08-22
状态：已与用户确认

## 背景

M3 PogFlow 达成 dev1 F1 15.1%（首次超越 classic 14.9%）。诊断显示 K 间 73% k-mer 重合但 15% 独有——多 K 的收益来自小 K 的连通性信息（更多 k-mer 共享→更长的初始 contig）注入大 K 的 read 证据链。

## 设计

### 流程变更

```
现有: 归一化 → K25 计数 → K25 Inchworm → Chrysalis(reads=原始) → PogFlow
M4:   归一化 → K19 计数 → K19 Inchworm → ≥200bp contig = 超级 reads
              → K25 计数 → K25 Inchworm
              → Chrysalis(reads=原始+超级) → PogFlow
```

### 超级 read 机制

- **RTT**：超级 read 的 25-mer 参与 k-mer 投票（帮助低覆盖转录本分配到正确组件）
- **QuantifyGraph**：25-mer 增加边权重（更多证据→流分解权重更准确）
- **Butterfly 穿线**：长 read 跨越更多节点（路径连续性更好）

### 实现要点

- `--super_reads` 旗标（默认 off）；触发 K19 预组装（K 固定 19，不暴露给用户）
- 超级 read header: `>sr_<N>` 格式（无 /1//2 后缀→Butterfly PE 配对不受干扰）
- 追加到 both.fa（或单独文件后在 RTT/Quantify 前合并）
- checkpoint: `.k19_preasm.ok`
- 运行时 ~1.5×（K19 inchworm 快，归一化共享）

## 验收

dev1 F1 ≥ 15.4%（M3 基线 15.1% + 至少 +0.3pp）；dev2/3 同向。
