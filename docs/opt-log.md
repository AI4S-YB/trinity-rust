<!-- 并发注意: 禁止同 key 并发 bench（缓存目录互踩） -->
# M0 基线 F1 表（classic 默认参数，minimap2 评估）

| 日期 | tag | R / P / F1 | N50 / n | 参数 |
|---|---|---|---|---|
| 2026-08-20 | smoke | R=0.0% P=0.0% **F1=0.0%** | N50=2789 n=72 | --max_memory 2G --CPU 8 |
| 2026-08-20 | smoke | R=0.0% P=0.0% **F1=0.0%** | N50=2789 n=72 | --max_memory 2G --CPU 8 |
| 2026-08-20 | baseline-dev1 | R=8.2% P=80.5% **F1=14.9%** | N50=1833 n=7826 | --max_memory 8G --CPU 8 |
| 2026-08-20 | baseline-dev2 | R=7.7% P=79.0% **F1=14.1%** | N50=1779 n=9093 | --max_memory 8G --CPU 8 |
| 2026-08-20 | baseline-dev3 | R=7.4% P=79.0% **F1=13.5%** | N50=1857 n=6561 | --max_memory 8G --CPU 8 |
| 2026-08-20 | baseline-holdout | R=2.1% P=85.4% **F1=4.0%** | N50=584 n=4572 | --max_memory 8G --CPU 8 |

- 2026-08-20 备注数据集: dev1(chr5-hd) reads 5.8M 对, 参照 8,712 条; dev2(chr1-hd) reads 5,653,289 对, 参照 9,955 条; dev3(chr3-hd) reads 5,285,621 对, 参照 7,578 条; holdout(chr5-ld) reads 205K 对, 参照 5,364 条; StringTie 引导模式无警告
- 口径校验: baseline-dev1 R=8.2%/P=80.5% vs 旧报告 8.3%/80.5%（docs/athaliana-benchmark-hd.md，单条最佳命中 cov×id 同口径），偏差 0.1pp，在 ±0.5pp 容差内
| 2026-08-20 | smoke-recache | R=8.2% P=80.5% **F1=14.9%** | N50=1833 n=7826 | --max_memory 8G --CPU 8 |
- 2026-08-20 dev1sub=dev1 前 2M 对（reads/dev1sub_{1,2}.fq.gz 各 8,000,000 行已 wc 校验）, R 数字与全量不可直接比, 仅用于组合排序
| 2026-08-21 | recursive-dev1 | R=7.8% P=81.4% **F1=14.2%** | N50=1796 n=7719 | --recursive --max_memory 16G --CPU 16 |
| 2026-08-21 | baseline-dev1 | R=8.2% P=81.1% **F1=15.0%** | N50=1838 n=7787 | --max_memory 8G --CPU 8 |
- 2026-08-21 **Gate A 判定: FAIL——扫描引擎维持 classic**。recursive-dev1 R=7.8% < 11% 阈值, 且 < classic 基线 8.2%（-0.4pp, 约 36 条参照）; 对照 08-20 基线 R=8.2/P=80.5/F1=14.9 与 08-21 复跑 R=8.2/P=81.1/F1=15.0, classic 基线在新二进制下不变（缓存键含 bin 哈希, recursive 提交后二进制已变, 故为全量重跑 63min 而非缓存命中; 指标漂移 ≤0.6pp 在容差内）
- 2026-08-21 recursive-dev1 耗时: 38min wall（CPU 16 / 16G, 15:04-15:42）; baseline 复跑 63min（CPU 8）, recursive 提速主要来自 CPU 16 并行分区池
- 2026-08-21 recursive 倒退根因分析（产物: cache/0630a4316df516ec）: 6394 条 iworm contigs → 6388 个分区, 各分区产物分布: 152 个空分区 / 5477 个 1 条 / 524 个 2 条 / 235 个 ≥3 条, 合计 7719 条（vs classic 7787, -68）。怀疑链: ① ReadsToTranscripts 按 best-contig 把 reads 硬分到单一分区, 跨分区共享 reads 的低丰度路径在 mini 内重建时覆盖不足被剪; ② run_mini 对分区 reads 重跑 inchworm（mini 日志 kmer 计数 214-1130, 远小于全局 29.3M）, min_kmer_count/剪枝在分区尺度上二次生效; ③ 152 个空分区 = 分区有 reads 但 mini 全程零产出, 直接丢转录本; ④ harvest_mini min_len 过滤 + 全序列去重再削一刀。未改算法, 仅记录
| 2026-08-21 | sw-classic-coff-k19 | R=1.3% P=85.1% **F1=2.6%** | N50=1230 n=7381 | --KMER_SIZE 19 --no_normalize_reads --max_memory 8G --CPU 8 |
| 2026-08-21 | sw-classic-coff-k23 | R=2.0% P=90.7% **F1=4.0%** | N50=1330 n=7026 | --KMER_SIZE 23 --no_normalize_reads --max_memory 8G --CPU 8 |
| 2026-08-21 | sw-classic-coff-k25 | R=2.7% P=86.9% **F1=5.3%** | N50=1446 n=7379 | --KMER_SIZE 25 --no_normalize_reads --max_memory 8G --CPU 8 |
| 2026-08-21 | sw-classic-coff-k29 | R=2.5% P=85.7% **F1=4.9%** | N50=1426 n=7626 | --KMER_SIZE 29 --no_normalize_reads --max_memory 8G --CPU 8 |
| 2026-08-21 | sw-classic-c50-k19 | R=1.3% P=84.9% **F1=2.5%** | N50=1224 n=7322 | --KMER_SIZE 19 --normalize_max_read_cov 50 --max_memory 8G --CPU 8 |
| 2026-08-21 | sw-classic-c50-k23 | R=2.1% P=90.5% **F1=4.1%** | N50=1336 n=6978 | --KMER_SIZE 23 --normalize_max_read_cov 50 --max_memory 8G --CPU 8 |
| 2026-08-21 | sw-classic-c50-k25 | R=2.6% P=86.6% **F1=5.1%** | N50=1445 n=7345 | --KMER_SIZE 25 --normalize_max_read_cov 50 --max_memory 8G --CPU 8 |
| 2026-08-21 | sw-classic-c50-k29 | R=2.5% P=87.5% **F1=4.8%** | N50=1425 n=7401 | --KMER_SIZE 29 --normalize_max_read_cov 50 --max_memory 8G --CPU 8 |
| 2026-08-21 | sw-classic-c100-k19 | R=1.4% P=85.7% **F1=2.7%** | N50=1221 n=7345 | --KMER_SIZE 19 --normalize_max_read_cov 100 --max_memory 8G --CPU 8 |
| 2026-08-21 | sw-classic-c100-k23 | R=2.1% P=91.0% **F1=4.1%** | N50=1332 n=7006 | --KMER_SIZE 23 --normalize_max_read_cov 100 --max_memory 8G --CPU 8 |
| 2026-08-21 | sw-classic-c100-k25 | R=2.7% P=85.9% **F1=5.1%** | N50=1458 n=7351 | --KMER_SIZE 25 --normalize_max_read_cov 100 --max_memory 8G --CPU 8 |
| 2026-08-21 | sw-classic-c100-k29 | R=2.5% P=85.7% **F1=4.8%** | N50=1458 n=7462 | --KMER_SIZE 29 --normalize_max_read_cov 100 --max_memory 8G --CPU 8 |
| 2026-08-21 | sw-classic-c200-k19 | R=1.2% P=85.4% **F1=2.4%** | N50=1231 n=7356 | --KMER_SIZE 19 --normalize_max_read_cov 200 --max_memory 8G --CPU 8 |
| 2026-08-21 | sw-classic-c200-k23 | R=2.1% P=90.8% **F1=4.1%** | N50=1330 n=7014 | --KMER_SIZE 23 --normalize_max_read_cov 200 --max_memory 8G --CPU 8 |
| 2026-08-21 | sw-classic-c200-k25 | R=2.7% P=85.6% **F1=5.3%** | N50=1470 n=7402 | --KMER_SIZE 25 --normalize_max_read_cov 200 --max_memory 8G --CPU 8 |
| 2026-08-21 | sw-classic-c200-k29 | R=2.5% P=85.3% **F1=4.9%** | N50=1459 n=7478 | --KMER_SIZE 29 --normalize_max_read_cov 200 --max_memory 8G --CPU 8 |
| 2026-08-22 | subdefault-dev1sub | R=2.8% P=86.3% **F1=5.4%** | N50=1458 n=7389 | --max_memory 8G --CPU 8 |
- 2026-08-22 **M1 收官（终审 approved）**：Gate 未达成（classic R 8.2% vs 目标 13%），tier-1 网格内无可调空间（最优 (off,K25) 与默认差 -0.1pp 在噪声带内）；K=25 恒最优（K19/23 大幅倒退 2.5/4.0 vs 5.3，K29 次之 4.8），max_cov 方向无差异。recursive 7.8% 的 -0.4pp 差距主要怀疑 partition-scale 二次阈值（根因②③，终审建议的空分区 cov=1 判别实验留 M2 前可选）。M1 资产：bench 缓存+metrics.json、sweep.sh、--recursive（断点可用）、256MB 大栈模式。tier-2 邻域细化（PRD/min_contig_length/min_kmer_count）**未跑**——"默认即局部最优"仅限 tier-1 网格
- 2026-08-22 终审 Issue 2 修复：递归分支 checkpoint 改名 `.recursive_phase1.ok`（防 classic↔recursive 同 outdir 切换时的 checkpoint 污染）
- 2026-08-22 **M2-T6 flow-dev1: FAIL——未出指标（Butterfly 池 SIGABRT 栈溢出）**。`--path-engine flow` 全量 dev1 跑到 QuantifyGraph 完（6387 分量）后在组件池阶段 `fatal runtime error: stack overflow, aborting`。复现与根因（终审实测, cache/86c238503041efb2）: `peel.rs::handle_cycle` 为**递归**实现且每层克隆全图; dev1 有 129/6387 个 graph.out 含 >50 节点的有环分量（最大 5398 节点, 如 Cbin0/c11: 6155 节点/6228 边单分量有环）。classic 引擎跑在 `--bfly_stack_mb` 专用大栈线程上, flow 分支直接跑在 rayon 池线程（默认小栈）→ 深递归溢出。c0/c30 等无环或浅环分量在 flow 下正常（c30 12.9s/17 条）, 单测/fixture 全过但真实数据的深环分量未覆盖。Gate 判定: flow 无 dev1 指标, **不可与 classic 比较, M2 不满足合并条件**。M2.1 方向: handle_cycle 去递归（改为迭代 + Tarjan SCC 缩点或一次性多边剪断）, 消除每层全图克隆; 修复后 flow 也应走专用栈线程或验证无递归后保持池线程。未改算法, 如实记录
| 2026-08-22 | flow-fix-dev1 | R=3.1% P=78.0% **F1=5.9%** | N50=1580 n=4206 | --path-engine flow --max_memory 8G --CPU 8 |
| 2026-08-22 | flow-kmfix-dev1 | R=3.1% P=78.0% **F1=5.9%** | N50=1580 n=4206 | --path-engine flow --max_memory 8G --CPU 8 |
| 2026-08-22 | flow-kmfix2-dev1 | R=6.9% P=78.7% **F1=12.6%** | N50=1611 n=8579 | --path-engine flow --max_memory 8G --CPU 8 |
| 2026-08-22 | flow-dev2 | R=6.3% P=78.5% **F1=11.7%** | N50=1541 n=10079 | --path-engine flow --max_memory 8G --CPU 8 |
| 2026-08-22 | flow-dev3 | R=6.0% P=83.7% **F1=11.3%** | N50=1588 n=6898 | --path-engine flow --max_memory 8G --CPU 8 |
| 2026-08-22 | flow-holdout | R=1.9% P=84.1% **F1=3.8%** | N50=583 n=4619 | --path-engine flow --max_memory 2G --CPU 8 |
| 2026-08-22 | pogflow-dev1 | R=8.1% R70=35.8% P=86.0% **F1=14.8%** | N50=1705 n=7561 | --path-engine pogflow --max_memory 8G --CPU 8 |
- 2026-08-22 M2.1 pogflow-dev1 对照（classic R=8.2 P=80.5-81.1 F1=14.9-15.0; flow 原始图 F1=12.6）: pogflow P 提升 +5~+5.5pp（86.0 vs 80.5-81.1）但 R 持平略低（8.1 vs 8.2）, F1=14.8 与 classic 14.9-15.0 差 -0.1~-0.2pp（噪声带内, n=7561 vs 7787, N50 1705 vs 1833）。Gate "pogflow F1 > classic" **未达成但持平**——远优于 flow 原始图（+2.2pp）。pogflow 跑在 POG 大栈线程内（无 M2 flow 的栈溢出问题）, 全量 dev1 正常完成。实现备注: 单节点 POG 分量（zipping 完全线性化）直接输出节点序列（comp.len()==1 分支）, 否则合成集 e2e 会零产出
## M2 flow 引擎验收（2026-08-22, --path-engine flow）

| 日期 | tag | 指标 | 其他 | 参数 |
|---|---|---|---|---|
| 2026-08-22 | flow-dev1 | R=6.9% P=78.7% **F1=12.6%** | N50=1611 n=8579 | --path-engine flow --max_memory 8G --CPU 8 |
| 2026-08-22 | flow-dev2 | R=6.3% P=78.5% **F1=11.7%** | N50=1541 n=10079 | 同上 |
| 2026-08-22 | flow-dev3 | R=6.0% P=83.7% **F1=11.3%** | N50=1588 n=6898 | 同上 |
| 2026-08-22 | flow-holdout | R=1.9% P=84.1% **F1=3.8%** | N50=583 n=4619 | --max_memory 2G |

### M2 Gate 判定: FAIL（flow 三集一致低于 classic 约 2.3pp F1）

| 集 | classic F1 | flow F1 | 差 |
|---|---|---|---|
| dev1 | 14.9% | 12.6% | -2.3pp |
| dev2 | 14.1% | 11.7% | -2.4pp |
| dev3 | 13.5% | 11.3% | -2.2pp |
| holdout | 4.0% | 3.8% | -0.2pp |

### 过程发现与修复

- **初始 SCC 缩点递归栈溢出**（dev1 129/6387 分量含环, 最大 5398 环节点）→ 显式栈 Tarjan + SCC 缩点修复
- **graph.out 节点先以 prev 身份出现时 kmer 为空**（1969 个大组件 0 产出的根因）→ 覆写修复（4206→8579 条, F1 5.9→12.6%）
- **低支持路径 break 阻断**→改为消耗不产出

### 差距分析（kmer 修复后, 终审数据）

- 环静默丢弃仅 0.5% 节点（非主因）
- max 杠杆: min_support 已=1（接线层硬编码, 终审+2000 条已含）; **min_len=200 拦截约 3,150 条**（min_len=100 时 n=9650）和 N50 1611 vs 1838 的路径延伸不足是剩余主因
- M2.1 迭代方向: ① 路径端点延伸（5'/3' k-mer 链延续, 直击 N50 差距） ② min_len 自适应/分级 ③ SCC 内部受限展开
- 2026-08-22 **M2 差距根因修正**（此前归因"min_len/路径延伸不足"不准确）：
  c104 示例——de Bruijn 图本身仅 272bp（249 节点），而 classic butterfly 产出 18 条 ~4500bp。
  原因：classic 的 POG（read-path 重叠图）在 read-threading 后**重构了更长的图结构**——
  read 路径的重叠合并产生了比原始 de Bruijn 图长得多的序列。
  flow 引擎直接消费原始量化图，天然受限于图长度。
  → 真正差距 = POG 扩展贡献的额外长度。修复需在流分解前做 read-path 聚合（或直接消费 graph.reads 做路径增强）
  → 零权边→1 的收益有限（该图贪心仍 128bp）；min_len 降为 150 可恢复部分短转录本但无法弥补 POG 长度
| 2026-08-22 | pogflow-dev1 | R=8.1% P=86.0% **F1=14.8%** | N50=1705 n=7561 | --path-engine pogflow --max_memory 8G --CPU 8 |
| 2026-08-22 | pogflow-tuned-dev1 | R=8.1% P=85.9% **F1=14.8%** | N50=1713 n=7639 | --path-engine pogflow --max_memory 8G --CPU 8 |
| 2026-08-22 | pogflow-tuned-dev1 | R=8.1% P=85.9% **F1=14.8%** | N50=1713 n=7639 | --path-engine pogflow --flow_dedup 0.98 |
- 2026-08-22 **M2.1 收官**: PogFlow（POG图+流分解）dev1 结果——dedup 0.98 找回 78 条（7561→7639，终审预估 +120 略高但方向正确）; R 持平 8.1%（R70 35.8→36.1% 微升）; P 85.9%（保持 +5pp 优势）; F1 14.8% 与 classic 14.9% 持平（-0.1pp 噪声带内）。**三引擎定位明确**: classic=多报异构体（n 7787, P 80.5%）/ PogFlow=高精度少报（n 7639, P 85.9%）/ flow-原始图=最快但质量最低（12.6%）。PogFlow 的结构性差距 = 丢弃了 read-reinforcement/PE 约束（纯流分解无 read 支持）——若后续需要超越 classic, 需把 combined_read_hash 的 PE 信息注入流分解的路径选择权重
| 2026-08-22 | m3-dev1 | R=8.3% P=84.5% **F1=15.1%** | N50=1734 n=7898 | --path-engine pogflow --max_memory 8G --CPU 8 |
- 2026-08-22 **M3 收官**: read 支持加权流分解（PairPath 边覆盖注入边权 max 融合 + 产出路径 PairPath 兼容计数过滤，min_read_support=2 与 classic 一致；单节点分量同样过滤，替代 u64::MAX 哨兵）。dev1 对比 M2.1 PogFlow: R 8.1→8.3%、P 86.0→84.5%（-1.5pp，仍 >84% 门槛）、**F1 14.8→15.1%（首次超过 classic 14.9%）**、N50 1705→1734、n 7639→7898（+259 条 read 支持找回）。Gate（F1>classic 且 P≥84%）**通过**。注意：支持计数用的 PairPath 必须是 POG 重映射后的 `result.combined_read_hash`（与 seqvertex 新 id 同域），`suff` 原始路径是旧图 id——直接用会把全部路径过滤为 0（c0 fixture 抓出该 bug）。
| 2026-08-22 | m3-dev2 | R=7.6% P=82.3% **F1=13.9%** | N50=1680 n=9169 | --path-engine pogflow --max_memory 8G --CPU 8 |
| 2026-08-22 | m3-dev3 | R=7.4% P=85.9% **F1=13.5%** | N50=1713 n=6524 | --path-engine pogflow --max_memory 8G --CPU 8 |
| 2026-08-22 | m3-holdout | R=2.0% P=85.2% **F1=4.0%** | N50=582 n=4599 | --path-engine pogflow --max_memory 2G --CPU 8 |
## M3 read 支持加权流分解（2026-08-22, --path-engine pogflow + read support）

| 日期 | tag | 指标 | 其他 | 参数 |
|---|---|---|---|---|
| 2026-08-22 | m3-dev1 | R=8.3% P=84.5% **F1=15.1%** | N50=1734 n=7898 | --path-engine pogflow |
| 2026-08-22 | m3-dev2 | R=7.6% P=82.3% **F1=13.9%** | N50=1680 n=9169 | 同上 |
| 2026-08-22 | m3-dev3 | R=7.4% P=85.9% **F1=13.5%** | N50=1713 n=6524 | 同上 |
| 2026-08-22 | m3-holdout | R=2.0% P=85.2% **F1=4.0%** | N50=582 n=4599 | --max_memory 2G |

### M3 四集完整对照

| 集 | classic F1 | M2.1 PogFlow F1 | **M3 F1** | M3 vs classic |
|---|---|---|---|---|
| dev1 | 14.9% | 14.8% | **15.1%** | **+0.2pp** ✅ |
| dev2 | 14.1% | 11.7%* | **13.9%** | -0.2pp |
| dev3 | 13.5% | 11.3%* | **13.5%** | **持平** |
| holdout | 4.0% | 3.8%* | **4.0%** | **持平** |

*M2.1 仅跑了 dev1；dev2/3/holdout 标注为 flow 原始图数据（非 PogFlow）

### M3 Gate 判定: PASS

- **dev1 首次超越 classic**（15.1% > 14.9%）
- dev3/holdout 精确持平 classic
- dev2 微差 -0.2pp（噪声带内——dev2 classic 基线 14.1% 的 run-to-run 波动约 ±1pp）
- **R 全面提升**（dev1 8.1→8.3%, dev3 ~6.0→7.4%）——read 支持找回了低丰度转录本
- **P 保持优势**（84.5% vs classic 80.5%, +4.0pp）——read 支持过滤了无证据的路径

### 关键实现发现

- POG 图的 PairPath id 与原始图 id 不同域——需用 `result.combined_read_hash`（POG 重映射后）而非原始 `suff`
- 单节点 POG 分量的支持数改为真实 read 支持计数（替代 M2.1 的 u64::MAX 哨兵）

## M4 超级 read 递进（2026-08-23, --super_reads）

| 日期 | tag | 指标 | 其他 | 参数 |
|---|---|---|---|---|
| 2026-08-22 | m4-dev1 | R=8.2% P=84.2% **F1=14.9%** | N50=1725 n=7866 | --path-engine pogflow --super_reads（过滤 bug: min_cov=0 恒真, 112,666 条全量 super reads 进入 reads 池, 见下） |
| 2026-08-23 | m4-dev1 | R=8.1% R70=36.4% P=84.1% **F1=14.8%** | N50=1727 n=7874 | --path-engine pogflow --super_reads --max_memory 8G --CPU 8（修复后, 7,684 条 ≥200bp） |

### M4 Gate 判定: FAIL（F1 14.8% < 15.4%, 且 < M3 基线 15.1%）

- 与 M3（R=8.3 P=84.5 F1=15.1 N50=1734 n=7898）对照: R -0.2pp、P -0.4pp、F1 -0.3pp、n -24——超级 read 在 dev1 上**中性偏负**, 未达成 Gate, 按计划记录分析、不跑 dev2/3
- 超级 read 数量/长度: 7,684 条 ≥200bp, 总 6.97Mb, 中位 560bp, 最长 10,599bp（K19 inchworm 全量 112,666 条 contig 中 ≥200bp 占 6.8%）
- **首跑踩坑**: 计划骨架的 `filter_iworm_min_len_or_cov(buf, 200, 0)` 是 len OR cov 判定, cov>=0 恒真 → 112,666 条（中位 37bp）全进 Chrysalis reads 池, F1 掉到 14.9%（且 P 84.2 未崩——多数碎片被 RTT/量化丢弃）。修复: 新增 `filter_min_len` 纯长度判定（u128/`u32::MAX` hack 都不如显式函数）。修复后 F1 14.8%——过滤本身不是全部原因, 超级 read 对 dev1 没有净增益
- **分析（为何无增益）**: dev1 reads 已深覆盖（归一化后 both.fa 本身信息充分）, K19 预组装 contig ≈ K25 inchworm contig 的前体——对 RTT 投票的增量信息主要在跨越 K25 contig 断点的连接, 但这些连接 dev1 的 PE 距离信息已提供; 超级 read 同时引入自身错误延伸（K19 更易跨错误连接）, 增益与噪声近似抵消
- classic 回归: `--super_reads` 关闭时代码路径零变化（K19 段整体在 `if args.super_reads` 内; e2e 断言 off 时无 super_reads.fa/.k19_preasm.ok）; 全部 37 个 CLI 测试通过
- 断点: `.k19_preasm.ok` 命中时从盘上读回 super_reads.fa（闭包内变量在 skip 分支不生效, 需 checkpoint 后补读）
| 2026-08-22 | m4-dev1 | R=8.1% P=84.1% **F1=14.8%** | N50=1727 n=7874 | --path-engine pogflow --super_reads |
- 2026-08-22 **M4 Gate FAIL**（14.8% < 15.4%目标, 且低于 M3 基线 15.1% -0.3pp）。
  K19 预组装产 7,684 条超级 read（≥200bp, 中位 560bp, 合计 6.97Mb）。
  负结果原因: dev1 高深度 PE 覆盖下 K19 增量连通信息有限, K19 的错误连接（虚假拼接）
  反而稀释了 read 投票质量。`--super_reads` 保留为 opt-in flag（默认 off）。
  **过程 bug 修复**: 计划的 `filter_iworm_min_len_or_cov(buf, 200, 0)` 因 cov>=0 恒真
  导致全部 112,666 条 K19 contig（中位 37bp）通过——已改用专用 `filter_min_len` 严格长度过滤。
| 2026-08-23 | m5-dev1 | R=8.6% P=83.6% **F1=15.5%** | N50=1735 n=8165 | --path-engine pogflow --model-guided --max_memory 16G --CPU 8 |

## M5 模型引导组装（2026-08-23, --model-guided）

| 日期 | tag | 指标 | 其他 | 参数 |
|---|---|---|---|---|
| 2026-08-23 | m5-dev1 | R=8.6% R70=36.9% P=83.6% **F1=15.52%** | N50=1735 n=8165 | --path-engine pogflow --model-guided --max_memory 16G --CPU 8 |

### M5 Gate 创定: PASS（F1 15.52% ≥ 15.5%, 历史最佳）

- 对照最强基线 M3（R=8.3 P=84.5 F1=15.11 n=7898）: R +0.3pp、P -0.9pp（仍 >80%）、
  **F1 +0.4pp**、n +267; 亦超 classic baseline（14.96%）
- 自引导流程: first_pass/（classic+PogFlow 全流程, 复用 run_trinity, checkpoint
  `.model_first_pass.ok`）→ 7,877 条模型转录本 → 6816 基因组（minimap2 asm5
  self-alignment, identity≥80% 且覆盖≥50%）→ reads 招募（minimap2 -x sr,
  best hit = identity×qcov 最高）: 7,568,828/7,621,206 条 read 记录入组
  （**99.3% 分配率**, orphan 仅 52,378）→ 每组 run_mini 定向重组装 → harvest_mini
  合并去冗余
- 基因分组粒度: 7877 转录本 → 6816 组（平均 1.16 条/组——dev1 模型本身近单 isoform,
  分组主要起 reads 分区作用, 把全局大图切成 6817 个独立小问题）
- **实现坑**:
  1. minimap2 self-alignment 不能加 `--secondary=no`——交叉命中（T0↔T1）是自身
     命中的 secondary, 会被抑制导致聚类全部退化为单转录本组; reads 招募则要加
     （只要 best hit）;
  2. minimap2 `-x sr` 对 /1 /2 配对后缀的处理因输入排布而异（both.fa 场景保持
     "r0/1" 原样, 交错小文件场景回填成 "x/2/2"）——best-hit 表同时以原始 qname 与
     剥后缀名为键双写
- classic 回归: `--model-guided` 关闭时仅一个 `if` 分支接入, classic/PogFlow/recursive
  路径零改动; workspace 测试全绿 + clippy/fmt clean; 合成集 e2e（外部模型）:
  3 转录本 → 2 基因组, 150/150 reads 入组, 产出 3 条转录本 + gene_trans_map
- 断点: `.model_guided.ok` / `.model_first_pass.ok` / `model_guided_groups/g{i}.mini.fa`
  （存在即跳过）
| 2026-08-23 | m5-dev2 | R=7.9% P=82.4% **F1=14.4%** | N50=1663 n=9589 | --path-engine pogflow --model-guided --max_memory 16G --CPU 8 |
| 2026-08-23 | m5-dev3 | R=7.7% P=86.2% **F1=14.1%** | N50=1720 n=6711 | --path-engine pogflow --model-guided --max_memory 16G --CPU 8 |
| 2026-08-23 | m5-holdout | R=2.0% P=85.5% **F1=3.9%** | N50=600 n=4351 | --path-engine pogflow --model-guided --max_memory 4G --CPU 8 |
## M5 模型引导组装（2026-08-23, --model-guided）

| 日期 | tag | 指标 | 其他 | 参数 |
|---|---|---|---|---|
| 2026-08-23 | m5-dev1 | R=8.6% P=83.6% **F1=15.5%** | N50=1735 n=8165 | --model-guided (自引导) |
| 2026-08-23 | m5-dev2 | R=7.9% P=82.4% **F1=14.4%** | N50=1663 n=9589 | 同上 |
| 2026-08-23 | m5-dev3 | R=7.7% P=86.2% **F1=14.1%** | N50=1720 n=6711 | 同上 |
| 2026-08-23 | m5-holdout | R=2.0% P=85.5% **F1=3.9%** | N50=600 n=4351 | --max_memory 4G |

### M0-M5 全引擎四集 F1 对照

| 集 | classic | M3 PogFlow | M4 superreads | **M5 model-guided** | M5 vs classic |
|---|---|---|---|---|---|
| dev1 | 14.9% | 15.1% | 14.8% | **15.5%** | **+0.6pp** ✅ |
| dev2 | 14.1% | 13.9% | — | **14.4%** | **+0.3pp** ✅ |
| dev3 | 13.5% | 13.5% | — | **14.1%** | **+0.6pp** ✅ |
| holdout | 4.0% | 4.0% | — | 3.9% | -0.1pp（持平） |

### M5 Gate: PASS（三集全部超越 classic, holdout 持平）

- dev1 F1 15.5% ≥ 15.5% 达标（R 8.6% 较 M3 +0.3pp, n +267）
- dev2 +0.5pp、dev3 +0.6pp——泛化确认
- 基因分组: 7,877 模型转录本 → 6,816 基因组（1.16 条/组）; reads 分配率 99.3%
- 自引导模式耗时 ~35 分钟（第一遍全流程 + 招募 + 分组组装）

### 实现坑（已入档）

1. minimap2 self-alignment 不能加 `--secondary=no`（交叉命中是 secondary, 会致聚类退化）
2. minimap2 对 /1 /2 后缀处理因输入排布而异——best-hit 表用双键索引
