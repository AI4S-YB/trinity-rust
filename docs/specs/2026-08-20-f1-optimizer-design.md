# trinity-rust 质量优化（F1 → 逼近 StringTie）— 设计文档

日期：2026-08-20
状态：已与用户逐节确认（§1-§4）
前置：P0-P5 复刻项目收官（v0.1.0）。本项目**放弃镜像保真目标**，转向组装质量优化。

## 1. 背景与诊断

高深度拟南芥 chr5 基准（`docs/athaliana-benchmark-hd.md`，~70×，StringTie 参照 8,712 条）：

| | StringTie | trinity-rust | 原版 Trinity |
|---|---|---|---|
| 转录本 | 8,712 | 7,831 | 8,847 |
| N50 / 平均长 | 2,248 / 1,850 | 1,833 / **1,201** | 1,952 / 1,308 |
| 全长恢复（≥95%） | — | 8.3% | 9.1% |
| 碱基精确度 | — | 80.5% | 76.9% |

参照转录本命运分解（我们）：**全长 15.7% / 部分(50-95%) 49.5% / 低覆盖 16.6% / 无命中 18.2%**。

**核心病灶**：不是基因检测不到，而是**碎片化与边界截断**（平均短 650bp；最大损失桶
是"部分覆盖"）。次生病灶：18% 无命中（低表达基因覆盖不足或错误剔除）。

## 2. 目标与约束（用户确认）

- **纯 de novo**：组装过程不接触参考基因组/注释；StringTie 仅作事后评分参照
- **主指标 F1** = 2·R·P/(R+P)，其中 R = 全长恢复率（参照被 ≥95%cov&≥95%id 重建比例）、
  P = 碱基精确度（de novo 碱基 ≥95% id 比上参照的比例）。评估口径固定为现有
  minimap2 asm5 流水线（`eval_vs_ref.py` 逻辑），全程可比
- **改动自由度：激进**——核心模块可替换；但保留 classic 路径作回退
- **评估数据**：开发/保留双集。dev-1=chr5 高深度（现有）；dev-2/dev-3=chr1/chr3 高深度
  （从现有 hd.bam 提取，零下载）；holdout=chr5 低深度（另一样本另研究，仅里程碑时碰）

## 3. 分层路线（用户选定）

### Phase 0：基准设施（M0，~0.5 天）

- `cargo xtask bench --reads L,R --ref ref.fa --tag 名 [-- trinity-cli 参数...]`：
  组装（缓存：reads+参数内容哈希）→ minimap2 → R/P/F1 + 参考列（≥70%、N50、条数）
- StringTie 参照生成脚本化（每染色体 GTF 重命名 → stringtie → gffread 去重），
  产物入 `athal/bench/`
- 实验日志 `docs/opt-log.md`：日期/参数/三集 F1/结论逐次追加
- **验收**：dev 三集 + holdout 的基线 F1 表落档

### Phase 1：已知杠杆（M1，~1-2 天）

1. **递归分区模式** `--recursive`：chrysalis ReadsToTranscripts 分区 → 逐分区进程内迷你
   流程（inchworm `--reads` 装载 + chrysalis 全链 + butterfly `run_component`，均已具备）
   → 输出拼接；`.ok` 断点沿用。预期 R +~10%（原版默认架构即此）
2. **参数扫描**：一档网格 diginorm max_cov {off,50,100,200} × K {19,23,25,29}；
   二档邻域细化（path_reinforcement_distance、min_contig_length、min_kmer_count）
- **防过拟合**：胜出参数须 dev-2/3 复现方向；**验收门**：dev-1 R ≥ 13% 且 F1 同向
- **保真回归策略**：新行为全部走旗标；`--classic` 路径不动；既有测试与 4 个 xcheck
  继续守护镜像路径（必要时显式 `--classic` 调用）

### 里程碑与验收

- M0：harness + 基线表
- M1：递归+调参收官，dev-1 F1 显著超 classic 基线（R ≥ 13%）

## 4. 评估协议与工程组织

- 调参/开发只看 dev-1/2/3；holdout 仅里程碑时碰
- 合并标准：dev-1 F1 升且 dev-2/3 倒退 ≤1pp（噪声带）
- 每实验 `opt/*` 分支，达标才并 main；提交带实验号与 F1 变化
- 既有镜像资产（509 测试、4 xcheck、fixtures）继续生效，守护 classic 路径

## 5. 风险

| 风险 | 对策 |
|---|---|
| 参数在 dev 上过拟合 | dev-2/3 复现规则 + holdout 里程碑制 |
| 递归模式内存/耗时膨胀 | 分区天然限制单分量规模；沿用组件池并行 |
| F1 口径被质疑（P 是碱基级、R 是转录本级） | 口径在 M0 基线表一次性定档，全程不变；参考列齐全可复查 |
