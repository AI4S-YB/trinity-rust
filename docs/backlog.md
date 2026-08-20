# 性能积压裁定（P5-T4）

P1 审查 I1 的三项优化候选——**分片开放寻址表 / 零拷贝流式记录迭代 / diginorm 管道化**——
基于 T4 实测裁定为**暂缓**（记入积压，P5 不做）。本文记录依据与 revisit 触发条件。

## 实测依据（2026-08-18，128 核登录节点）

### 1. 全管线峰值 RSS vs --max_memory（T4 基准）

数据：sample_data/test_Trinity_Assembly 全量 PE（reads + reads2 各侧合并，
fq.gz 合计 ~5.0 GB 解压前 23.4 MB gz / 118 MB fq），diginorm 后 both.fa ≈ 21.7 MB。

| 配置 | wall | 峰值 RSS | --max_memory | RSS / max_memory |
|---|---|---|---|---|
| `--CPU 8 --max_memory 2G` | 17.4 s | **868 MB** | 2147 MB | **0.40×** |
| `--CPU 1 --max_memory 2G`（同数据） | ~35 s | （低于左值） | 2147 MB | <0.40× |

产出 92 条转录本；`--CPU 1` 与 `--CPU 8 --inchworm_cpu 1` 最终
`*.Trinity.fasta` **逐字节一致**（GFF/RTT/quantify 的 -t 确定性验证；
inchworm 自身 PARALLEL 轮的种子平局序差异为原版固有，与本裁定无关）。

### 2. 分阶段已知峰值（docs/benchmarks.md P1 审查补记）

- diginorm 全流程 88 MB PE fq 输入（~7.9M distinct 25-mer）：128 线程 2.32 GB、
  4 线程 1.04 GB——峰值随 rayon 并发段数放大是**结构性**的（fold 累加器趋近全表），
  但绝对值在 max_memory 约束内可控（--CPU 已接入全局池，用户可压线程数）。

## 裁定

- **RSS 可控**：全管线峰值 0.40× max_memory（<2× 阈值），绝对值 <1 GB @ 全量小数据
  → 三项优化全部记入积压，P5 不做。
- 已落地的替代缓解（T4）：
  1. `--max_memory` 硬护栏（diginorm/计数, `orchestrate::memory_guard_error`）
     + 软警告（chrysalis/butterfly）——超限时报错而非 OOM;
  2. `estimate_hash_size` → 计数 HashMap 每线程预 reserve（`count_fasta_data_with_capacity`,
     上限 4M 槽/线程），抑制增长期 rehash 峰值;
  3. rayon 全局池按 `--CPU` 一次设定 + 各并行入口显式 scoped 池——`-t 1` 真单线程。

## Revisit 触发条件（任一命中即重启三项之一）

| 条件 | 首选项 |
|---|---|
| 单文件输入（fq/fa）> 2 GB，或 both.fa > 4 GB | 零拷贝流式记录迭代（trinity-common） |
| 全管线峰值 RSS > 2× --max_memory（护栏放行后仍被 OOM kill） | 分片开放寻址表（按 kmer 哈希分片, 内存=1×表） |
| diginorm 阶段 RSS > 8× 输入字节数且 --CPU 无法压下 | diginorm 管道化（raw/fa/both/stats 多副本流式化） |
| dense 数据（distinct k-mer ≈ 输入碱基数）下计数 RSS > 12× 表大小 | 分片开放寻址表 |

## 其他遗留（P5 范围外，顺带记录）

- jellyfish histo 及其 `.ok` 未复刻（无下游消费）;
- chrysalis 六子阶段用单一 `.quantify_graph.ok` 包裹（粒度差异见 orchestrate 模块文档）;
- inchworm PARALLEL 模式的种子平局序与原版同为非确定（跨线程数不可复现）。
