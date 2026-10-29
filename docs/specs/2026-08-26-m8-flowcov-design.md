# M8 覆盖度加权流分解 — 设计文档

日期：2026-08-26
状态：已与用户确认（探针结论驱动，见 opt-log M8 段）

## 核心

PogFlow 流分解的边权兜底从常数 1 改为 **25-mer 覆盖度（封顶融合）**：
read 证据优先，无 read 证据的边用 `min(cov25, CAP)` 兜底——
流分解在分叉处不再"盲选"（现状：POG 边权/PairPath 覆盖均为 0 时一律 floor 1，
DP 只剩路径长度可分）。

一句话：**read 证据管不了的边，覆盖度管。**

## 依据（M8 探针，2026-08-26 已入档 opt-log）

- 决策点出边：正例 e2_cov 中位 67 vs 负例中位 8（8× 间隔）；
- GBM ≈ 单特征覆盖度规则（dev3 决策级 76.9% vs 76.6%），其余 17 个特征
  合计重要性 <0.02——DL/GNN 路线已否决，信号收敛为覆盖度；
- 纯配对切片已被覆盖度排序解决 95.9%；混合切片"错误"63% 为 StringTie
  漏报的高覆盖真实分支（标签噪声，非特征问题）。

## 机制

融合点：`crates/trinity-butterfly/src/component.rs::digraph_to_flow`

```
现状:  w = max(pog_w, read_w, 1)            // 无证据边恒 1
改后:  w = max(pog_w, read_w)
       if w == 0 { w = max(1, min(cov25, CAP)) }   // 覆盖度只补 floor，不压 read 证据
```

- **junction 25-mer**（POG 边 (u→v)，节点为 zipping 后变长序列）：
  取 `u.name` 末 24 碱基 + `v` 经最长后缀-前缀重叠后的首个延伸碱基
  （与 `trinity_flow::paths::path_to_seq` 同一重叠函数）；
  `u.name` 不足 24 时向前借 v 的延伸补足 25；v 无延伸（包含关系）或拼接
  不足 25 → cov=0（等价现状 floor 1，不产生行为差异）。
- **计数表**：`outdir/mer_counts.25.asm.rs.bin`（butterfly 阶段尚在——收尾
  清理在 harvest 之后才删）；缺失 → 回退 `both.fa` 重数（同在收尾才删）；
  都没有 → **优雅降级**：禁用融合、stderr 警告一次，行为与现状逐边一致。
  SS 数据集用对应 `mer_counts.25.ss.rs.bin`（k/ds 校验 load_counts_bin 内建）。
- **CAP**：唯一待标定旋钮。read support 尺度小（min_read_support=2）、
  覆盖度中位 67，CAP 限定覆盖度只影响"无 read 证据的边之间"的相对排序，
  不颠覆 read 证据主导权。候选扫描 {2, 5, 10, 50}，dev1sub 排序 + dev1 定稿。
- 实现形态：计数表在 butterfly_pool 装载一次，`Arc<CountMap>` 跨 rayon
  线程共享（只读）；组件级零拷贝。

## CLI

`--flow-cov-fusion`（opt-in，默认关——classic/PogFlow 现状路径零改动）；
`--flow-cov-cap N`（默认 10，扫描后定稿）。

## 验收

- **Gate**：dev1 `--path-engine pogflow --flow-cov-fusion` F1 > 15.1%
  （M3 基线）且 P ≥ 84%（M3 的 84.5% 口径）；
- 泛化：dev2/dev3/holdout 同配置复跑，F1 不低于 M3 对应集 −0.2pp；
- **R/P/n 分解审视**（循环性警示）：若 P 降 R 升且 n 增，需判定是"找回
  StringTie 漏报的真实分支"（可用 Ensembl r60 cdna 参照复评 P 佐证）还是
  真精度损失——探针已证 63% 的"参照负例高覆盖边"属前者；
- 回归：`--flow-cov-fusion` 关闭时全部产物逐字节同现状（e2e 断言）。

## 风险

1. CAP 过大 → 覆盖度压过 read 证据，退回"覆盖度独裁"（探针显示其决策级
   上限 ~77%，不如 read 证据+覆盖度组合）；扫描兜底。
2. POG 拼接产生 raw reads 中不存在的 junction 25-mer（cov=0）——不救也不伤，
   与现状一致。
3. checkpoint 断点续跑场景 bin/both.fa 已删 → 降级路径必须逐边等价（测试锁定）。
