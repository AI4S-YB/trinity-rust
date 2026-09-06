# M7 de novo 组装提速 — 设计文档

日期：2026-08-25
状态：已与用户确认

## 核心

真实数据基准（10M PE，--CPU 8）的 classic 全流程从 **44 分钟降到 ≤15 分钟
（-66%）**，质量指标不退。路径：阶段计时探针先行 → 消灭 dump 中转 →
F2DB 串行段并行化 → 中期 checkpoint → 按探针数据定向取用优化池。

## 基线画像（2026-08-25 实测，real-classic，三套运行一致）

| 阶段 | 耗时 | 占比 |
|---|---|---|
| 归一化 + prep | 1.0 min | 2% |
| K25 计数 | 1.2 min | 3% |
| **kmer dump 落盘** | **8.0 min** | 18% |
| inchworm 本体 | 2.0 min | 5% |
| GraphFromFasta+Bundle | 1.6 min | 4% |
| **RTT + F2DB + partition** | **24.3 min** | 55% |
| QuantifyGraph | 1.5 min | 3% |
| butterfly + harvest | ~3 min | 7% |
| 蛋白 stage | 0.7 min | 2% |

已知结构问题：F2DB（`graph_per_record`）移植时被串行化（原版 OpenMP 并行）;
RTT 已并行但仅 8 线程; dump 是 计数表→落盘→inchworm 读回 的三重浪费。

## 已确认决策

| 决策点 | 选择 |
|---|---|
| 确定性红线 | **指标等价即可**（R/P/F1/P_full 不退; 能免费保住逐字节一致处顺手保住） |
| CPU 政策 | **固定 --CPU 8**（与全部历史数字可比; 机器 128 核不动用） |
| Gate | real-classic ≤ **15 分钟**（-66%） |
| 实现路径 | A：探针先行 + 三层递进（消浪费→并行化→定向算法） |

## 组件设计

### ① 阶段计时探针

- `orchestrate.rs` + `trinity-chrysalis/src/pipeline.rs` 每个子阶段
  `Instant` 计时; stdout 的 `-running X...` 行尾带耗时
- 落盘 `outdir/.stage_timings`（追加行 `阶段名\t秒`）
- xtask bench 结束打印 top-5 耗时段（从该文件读）

### ② 消 dump 中转（两档, 探针定档）

- 第一档: dump 字节留内存直传 `populate_from_kmers`（不落盘）
- 第二档（若 dump 文本生成本身是大头）: 新增 `populate_from_table`
  直接从计数 HashMap 构建, 绕过文本格式
- 断点语义: `.jellyfish_count.ok`/`.jellyfish_dump.ok` 合并为内存段
  （重跑重算, 不保留中间文件断点）

### ③ F2DB 并行化

`graph_per_record` 逐 bundle 记录独立 → rayon `par_iter`（`p.threads` 线程）
收集后**按输入序拼接** → 输出逐字节不变。

### ④ 中期 checkpoint（先探后定内置）

①②③ 完成后冷跑 `real-classic` 量 wall:
- 已 ≤15 min → 直接验收
- 缺口 <5 min → 优化池定向补刀
- **缺口 >8 min → 带探针数据回谈 Gate**（止损, 不硬凑）

### ⑤ 定向优化池（按探针取用, 不预做）

| 候选 | 针对场景 | 风险 |
|---|---|---|
| pipeline 双份 reads 拷贝消除（`reads_text: String` + `reads_seqs: Vec` 并存） | RTT 段 | 低 |
| RTT kmer 查询结构优化（FxHashMap → 排序数组二分/完美哈希） | RTT 主导 | 中 |
| partition 写盘合并（6388 小文件批量 buffered write） | partition 主导 | 低 |
| inchworm 种子排序并行化 | inchworm 段（~2 min, 池底） | 中 |

## 验证矩阵与 Gate

- **Gate**: real-classic（10M PE, --CPU 8, 冷跑）wall ≤ 15 min
- **护栏**:
  - real-classic 与 real-pogflow 的 R/P/F1/P_full 对比基线 ±0.1pp
    （基线: classic 3.5/82.8/6.7/40.1/52.6; pogflow 3.5/84.2/6.7/40.1/52.5）
  - dev1 classic 复跑指标持平
  - workspace 测试全绿
- 速度与指标入档 opt-log「M7 性能」段

## 风险

- -66% 为激进目标: checkpoint 缺口 >8 min 时回谈, 不牺牲指标护栏凑数字
- RTT/F2DB 块内比例未知: 探针先行即为此设计; 若 RTT 占 20+ min 且算法优化
  不达, Gate 可能需回调——届时以数据回谈
- F2DB 并行输出序拼接若引入隐性依赖（记录间共享状态）: 记录独立已由
  原版 OpenMP 并行佐证, 实现时以逐字节对拍兜底
