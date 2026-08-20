# 架构发现：v2.15.2 默认流程是两阶段递归式（2026-08-19 真实数据实验期间发现）

## 事实（Trinity 主脚本源码）

`run_chrysalis`（Trinity:2265-2273）中：

```perl
if (! $TRINITY_COMPLETE_FLAG) {
    ## Trinity Phase-1 only!
    &run_recursive_trinity($sorted_reads_to_components_file);
    return ("RECURSIVE_TRINITY_COMPLETE");
}
## Trinity phase-2: full assembly below （仅 --trinity_complete 时到达）
```

即 v2.15.2 **默认**（无 --trinity_complete）流程为：

1. **Phase-1（全局）**：归一化 → jellyfish+Inchworm → Chrysalis **前半**（GraphFromFasta 聚类 →
   BubbleUp → bundle → ReadsToTranscripts → reads 按组件分区落盘）——**不跑** FastaToDeBruijn /
   QuantifyGraph / Butterfly（butterfly_commands 路径只在 phase-2 内生成）
2. **递归阶段**：每个 read 分区（数以千计）各自运行一次完整**迷你 Trinity**
   （`Trinity --trinity_complete --full_cleanup`）：内部启用 TRINITY_COMPLETE_FLAG 语义
   （L1028-1037：inchworm_cpu=1、NO_PARALLEL_IWORM、不归一化、**FORCE_INCHWORM_KMER_METHOD
   =inchworm 内置 kmer 目录**——即 `--reads` 装载模式），跑自己的 prep→inchworm→chrysalis
   phase-2（deBruijn/Quantify）→butterfly
3. 最终 Trinity.fasta = 各分区输出拼接

## 对本项目的影响

- 我们移植的"经典单遍管线"（全局 chrysalis phase-2 + butterfly）对应**递归实例内部**的
  phase-2 路径，而非外层默认编排。P5 端到端对比（含 rep10 十次实验）实际是**跨架构对比**
  （我们=经典式，原版=递归式）——这是 ~50% 覆盖带的首要解释，此前归因"哈希迭代序系统性
  偏好"仅是次要因素。rep10 报告的判定需要此修正。
- 已具备的组件覆盖了递归模式所需的**全部内部阶段**（inchworm --reads 装载、chrysalis 全链、
  butterfly 组件驱动、prep）。缺口仅是**外层编排**：分区写盘（run_recursive_trinity 的
  read_partitions 布局）+ 逐分区迷你流程调度 + 输出拼接（~1-2 天工程量）。
- P1-P4 的阶段级对拍结论**不受影响**（同输入对同阶段，字节级证据仍有效）。

## 待办

- [ ] P7（若立项）：实现递归编排层（`--recursive` 模式或直接作为默认），对齐 v2.15.2 默认行为
- [ ] rep10-report.md 的"判定"段补引本文档的修正
