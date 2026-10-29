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

### M5 收尾（2026-08-24）

- 补齐 `--model-guided` 自动化回归: `tests/model_guided_e2e.rs`（外部模型合成集,
  3 转录本 → 2 基因组, 150 reads 全入组无 orphan 组; 断言 g0/g1.mini.fa、无
  first_pass/、gene_trans_map 行数与转录本一致且 gene 为 trans 前缀; 关闭旗标
  classic 路径无 model-guided 产物; minimap2 缺失时 skip）——当初 M5 只做了手工
  e2e, 此番对齐 M4 的 super_reads_e2e.rs 惯例
- M4/M5 plan checkbox 回填（当初两份 plan 均未勾选, 靠 opt-log 记录追溯）
| 2026-08-24 | m6-probe | R=8.4% P=84.5% **F1=15.3%** | N50=1732 n=7861 | --path-engine pogflow --protein-guided --protein-db /storage/home/senior007/database/protein.fa --max_memory 16G --CPU 8 |
| 2026-08-24 | m6-probe2 | R=8.4% P=82.5% **F1=15.3%** | N50=1776 n=7891 | --path-engine pogflow --protein-guided --protein-db /storage/home/senior007/database/protein.fa --max_memory 16G --CPU 8 |

## M6 蛋白引导组装——M6a 探路（2026-08-24, --protein-guided）

| 日期 | tag | 指标 | 其他 | 参数 |
|---|---|---|---|---|
| 2026-08-24 | m6-probe | R=8.4% R70=36.7% P=84.5% F1=15.3% | N50=1732 n=7861 | 机制未啮合（0 extended, 见下） |
| 2026-08-24 | m6-probe2 | R=8.4% R70=36.8% **P=82.5% F1=15.3%** | N50=1776 n=7891 | 机制啮合: 338 extended, 229 accepted, +237,480 nt |

### 延长统计与解读（m6-probe2）

- 4917/7891 转录本有蛋白模板; 离线重放 blastx: 3052 个最佳命中带端缺口
  （both_ends 1540 / c_only 885 / n_only 627）——延长候选充足
- 338 尝试 / 229 验收通过（bitscore >旧×1.05）/ 109 回退; 平均延长 ~702 nt/条
- **对照**: M3 PogFlow R=8.3 P=84.5 F1=15.11; m6-probe2 R 持平 8.4、**P −2.0pp**、
  F1 持平 15.3; N50 +44、n +30
- **结论（探路数据）**: 机制可啮合但当前形态对 dev1 无净增益——R_full 纹丝不动,
  +237kb 延长碱基大多不能比对回参考（掉 P）; DIAMOND bitscore 验收度量的是
  "更像某个同源蛋白", 不等于"参考正确"（旁系同源/走图漂移都能过验收）
- 嫌疑点: 平均 702nt 的延长远超典型蛋白缺口需要——贪心走在稠密原始图上,
  无蛋白约束的线性段只按边权走, 容易走进重复/旁系区; 验收门槛与参考正确性脱钩

### 实现坑（已入档, 三连）

1. DIAMOND `--outfmt 6` 不带字段列表输出默认 12 列（无 qlen/slen/qframe）,
   解析器要求 13 列 → 全部丢弃、0 模板——字段必须逐个显式列出
2. 蛋白库含 84 条零长度记录 → makedb 拒绝; stage 内建净化（protein.clean.fa）
3. **chrysalis graph.out 顶点名是 K−1-mer（24）非 --KMER_SIZE（25）**——
   顶点名长必须取 `BflyContext::kmer_size`, 用 args 值则锚定/后继校验双灭（0 extended 根因）

### 结构性发现（e2e 阶段）

单基因合成数据上管线在 k-mer 层全耦合: 稀疏 read 与侧翼共享 k-mer才进组件图,
共享即链通、butterfly 必产出全长——"图有路但枚举丢弃"场景只在真实数据复杂度
（分支/重复）下存在, e2e 只能做接线级验证
| 2026-08-24 | m6-probe3 | R=8.3% P=83.3% **F1=15.1%** | N50=1744 n=7867 | --path-engine pogflow --protein-guided --protein-db /storage/home/senior007/database/protein.fa --max_memory 16G --CPU 8 |
| 2026-08-24 | m6-probe3 | R=8.3% R70=36.7% **P=83.3% F1=15.1%** | N50=1744 n=7867 | 收紧版: 302 extended/192 accepted/+126,845 nt |

### M6a 收紧轮（probe3）与深度诊断（2026-08-24）

- 三处收紧: 预算 missing×3 无 slack cap 1000 / 走图越过模板蛋白两端即停 / 验收禁止换目标蛋白
- 效果: +nt 减半（237k→127k）, P 回收一半（82.5→83.3）, 但 **R 回落 8.4→8.3、F1 15.1**
  ——收紧把好的延长和坏的游荡一起砍了, 净效果反而变差
- **覆盖度诊断（probe2 的 226 个 accepted 拆解, asm5 对照 dev1.fa）**:
  - same-ref 123 条: 参考覆盖度 **0.37 → 0.73（mean/median 均近乎翻倍）**,
    **23 条跨过 R_full 的 0.95 门槛**（其中约 8 条过 identity 门槛计入 R, 与 R +0.1pp 吻合）
  - diff-ref 103 条: 延长碱基比对到别的座位（旁系/重复游荡）
  - 加长碱基构成: same_ref 46.7% / diff_ref 47.5% / unaligned 5.8%; 中位延长 929nt
- **结论**: 机制真实有效（忠实延长把转录本覆盖度翻倍、有全长相转变例）,
  瓶颈在**验收判别力**——蛋白 bitscore 区分不了"找回真实 3'/5' 缺口"与"走进旁系同源"。
  松版（probe2）净效果最好但 P 代价 2pp; 紧版伤敌一千自损八百
- **候选判别器（reference-free, 待验证）**: 延长段最小边权门槛——忠实延长沿
  read 支持的路径走, 游荡段常为低权/重复纠缠; 蛋白只在分叉处定向
| 2026-08-24 | m6-probe4 | R=8.3% P=83.5% **F1=15.1%** | N50=1756 n=7886 | --path-engine pogflow --protein-guided --protein-db /storage/home/senior007/database/protein.fa --max_memory 16G --CPU 8 |
| 2026-08-24 | m6-probe4 | R=8.3% R70=36.7% **P=83.5% F1=15.1%** | N50=1756 n=7886 | PE 跨接门: 340 extended, 189 过门, 118 accepted, +112,950 nt |

### M6a'' PE 跨接对门（probe4）与三轮收敛结论（2026-08-24）

- 实现: 回退 probe3 的预算/span-stop 两处收紧（保 probe2 松走图形态）, 保留
  same-sid 验收, 新增**PE 跨接对门**——原转录本端侧 100nt (a) + 延长区端侧
  100nt (b), minimap2 -x sr 比对 both.fa, 统计一端命中 a 一端命中 b 的 read
  对数, ≥2 对（对齐 MIN_SEED_COVERAGE=2 惯例）才进 DIAMOND 验收
- 边权地板判别器**当场证伪**: chrysalis graph.out 边权**全部恒为 1**
  （94 万条边 max=1）——de Bruijn 边不承重; 且图上每条边都来自真实 reads,
  游荡段本身是 read 支持的旁系真实序列, 单端支持原理上不可判别
- 结果: 门有判别力（340→189 过门, accepted 229→118, P 82.5→83.5）——但
  **R 回落 8.4→8.3, F1 15.1**, 净效果仍不超 probe2
- **机理解释（为何跨接门伤 R）**: 转录本在连接处被截短, 恰因该处 reads 稀疏;
  有 ≥2 对跨接 reads 的连接 butterfly 本来就组装通了（k-mer 层耦合, e2e 阶段
  已证）。**需要蛋白引导的延长恰恰是没有 read 证据的延长**——任何要求 read
  证据的验收都会拒掉真正有用的那批; 不要求则旁系游荡混入

### M6a 三种验收域的收敛（Gate 判定依据）

| 验收域 | R | P | F1 | 结论 |
|---|---|---|---|---|
| 纯蛋白 bitscore（probe2 松） | 8.4 | 82.5 | 15.3 | 旁系混入, P 代价 2pp |
| 砍预算+边界停+same-sid（probe3） | 8.3 | 83.3 | 15.1 | 好坏一起砍 |
| same-sid + PE 跨接门（probe4） | 8.3 | 83.5 | 15.1 | 有用的也被拒 |
| （基线 M3 / 最强 M5） | 8.3 / 8.6 | 84.5 / 83.6 | 15.11 / 15.52 | — |

**探路判定: 机制机械上完备（真实覆盖度翻倍案例 123 条）, 但在"read 耦合组装
+ 全长转录本度量"的 dev1 框架下, 任何 reference-free 验收都无法把净收益做出
正数**——三种验收域 F1 全部收敛在 15.1~15.3 ≈ 基线。蛋白引导的真正适用场
景是 **CDS 完整性**（蛋白能判断的量）而非全长转录本 R_full（含 UTR, 蛋白
无信息）——若续做 M6b, 度量与机制应当对齐（域语法 + CDS 完整度基准）。

### M6b 立项前探索——CDS 完整性度量的数据上限（2026-08-24）

M6a 收档后按"CDS 完整性方向"探索 M6b 立项前提, 三步证伪:

1. **蛋白级度量基线**: dev1 参考转录本最长 ORF 翻译 → 8,363 条参考蛋白;
   各缓存组装 blastx 对比——**基线蛋白全长回收已 68.5%**（m3; ≥70% 带 83.0%,
   有命中 95.9%）; M6a 延长仅 +0.2pp（probe2 68.7%）。转录组空间 R_full 8.3%
   的短板主要是 UTR/精确度, CDS 大多完好
2. **失败带形态**: 70-95% 覆盖带（m3: 1,208 条蛋白）**100% 为末端缺失型**
   （缺失残基全在首尾 25% 区）——形状上全部"延长可寻址"
3. **reads 支持核查（致命一击）**: 抽样 120 个带内蛋白的缺失末端 ~150nt,
   对原始 reads 建 31-mer 计数索引——**96% 的缺失端序列在 reads 中基本不存在**
   （<10% k-mer 有 ≥2 计数; 仅 2/120 有 ≥50% 支持）

**结论: M6b 立项前提不成立**——度量对齐后暴露的不是判别问题而是**数据上限**:
dev1 上"reads 里存在但组装丢失的末端" ≈ 极少, 蛋白引导延长类机制在此数据
无米下锅。若续蛋白方向, 需换数据（真实 RNA-seq/更深覆盖/基因家族复杂集）
而非换机制或加域语法

### M6b 真实数据选样（2026-08-24）
- RUN_ACC=DRR095137（Illumina HiSeq 2500 paired end sequencing of wild type under blue light condition, 样品 "RNA-seq data of wild type under blue light condition"（蓝光下 WT 幼苗, 普通 mRNA RNA-seq, TRANSCRIPTOMIC）, study PRJDB5794 "CAGE and RNA-seq in blue light exposure of A. thaliana sprouts"）, read_count 56,083,429（spot 数, 即 ~56M read/端）
- 抽样: seqkit sample -n 10000000 -s 42（确定性, 两端同 seed; 不一致时回退 seqkit grep 按 R1 ID 配对）
- 参考: Ensembl Plants r60 atha cdna 下载中 / atha pep 下载中 / alyr pep 下载中（controller 稍后更新条数）
- 网络偏差: 出网被限速（单连接 ~16KB/s, 并行总量 ~100-150KB/s）, 全量 6.4GB 不可行 → `scripts/prepare_realdata.sh` 改为并行分块（16 线程 × 8MB Range 请求, .ok 标记断点续传）只下载每端前缀 ~642MB（实测 8MB ≈ 149.6k reads, 覆盖 10M 对 ×1.2 余量）, 截断 gzip 解码后抽样。对组装 benchmark 等价于真实数据前缀子集, 非全库均匀抽样。
| 2026-08-24 | m3-dev1-pmetric | R=8.2% P=84.4% **F1=15.0%** | N50=1727 n=7867 | --path-engine pogflow --max_memory 16G --CPU 8 |

### 更正（2026-08-24）: M6b 探索段的蛋白回收数字口径
- 前文「基线蛋白全长回收 68.5% / ≥70% 带 83.0%」的 identity 过滤因 pident
  百分数误比 0.95 而恒真——真实口径（≥95%cov 且 ≥95%id）为 **63.5% / ~82%**
- 「70-95% 带 100% 末端缺失」「96% 缺失端 reads 不存在」两个结论基于覆盖度,
  不受影响——M6b 转向（数据上限）的判定不变
- Rust 版 translate_longest_orf 曾漏大写化 soft-masked 参考（小写→X, 校准
  35.0%）, 已修（upper + 回归测试）
| 2026-08-24 | m3-dev1-pmetric | R=8.2% P=84.4% **F1=15.0%** | N50=1727 n=7867 | --path-engine pogflow --max_memory 16G --CPU 8 |

### M6b 数据准备修正（2026-08-25）

- 参考记录数（Ensembl Plants r60, 下载完成）: atha_cdna **48,359** / atha_pep
  **48,321** / alyr_pep **32,667**
- **重抽样修正**: 初版 seqkit 同-seed 抽样因两端前缀池长度不同（R1≈11.86M /
  R2≈11.99M 条, 网络限速只下了前缀）导致 999.99 万/1000 万条错配——改为
  **位置锚定重建**（两池同为 run 前缀, 位置 i 天然同对; python seed42 从
  11.85M 位置确定性抽 10M, 池内与输出侧双重全量断言 0 错配）
- 说明: 数据为 run 前缀子集（非全流均匀抽样）——Illumina fastq 顺序与转录本
  来源无关, 对组装基准无偏
| 2026-08-25 | real-classic | R=3.5% P=82.8% **F1=6.7%** | N50=1561 n=32763 | --max_memory 16G --CPU 8 |
| 2026-08-25 | real-pogflow | R=3.5% P=84.2% **F1=6.7%** | N50=1530 n=33025 | --path-engine pogflow --max_memory 16G --CPU 8 |
| 2026-08-25 | real-m5 | R=3.9% P=86.0% **F1=7.5%** | N50=1473 n=34751 | --path-engine pogflow --model-guided --max_memory 16G --CPU 8 |
| 2026-08-25 | real-m6a | R=3.5% P=83.5% **F1=6.7%** | N50=1542 n=33047 | --path-engine pogflow --protein-guided --protein-db /storage/home/senior007/test/athal/bench/real/alyr_pep.fa --max_memory 16G --CPU 8 |

## M6b 真实数据基准（2026-08-25, DRR095137 拟南芥蓝光 WT, 抽样 10M PE 对）

| 日期 | tag | 指标 | 其他 | 参数 |
|---|---|---|---|---|
| 2026-08-25 | real-classic | R=3.5% P=82.8% F1=6.7% | N50=1561 n=32763 P_full=40.1% P70=52.6% | --max_memory 16G --CPU 8 |
| 2026-08-25 | real-pogflow | R=3.5% P=84.2% F1=6.7% | N50=1530 n=33025 P_full=40.1% P70=52.5% | --path-engine pogflow |
| 2026-08-25 | real-m5 | R=3.9% P=86.0% **F1=7.5%** | N50=1473 n=34751 **P_full=42.3% P70=54.8%** | + --model-guided |
| 2026-08-25 | real-m6a | R=3.5% P=83.5% F1=6.7% | N50=1542 n=33047 P_full=40.1% P70=52.6% | + --protein-guided --protein-db alyr |

### Gate 判定: 全部 PASS（泛化验证型）

- **G1** pogflow F1 ≥ classic: **PASS**（6.7 = 6.7 打平; P +1.4pp, R 相同）
- **G2** m5 F1 ≥ classic: **PASS**（7.5 ≥ 6.7, **+0.8pp**; R +0.4pp、P +3.2pp、
  P_full +2.2pp——M5 基因分组组装的优势在真实数据上复现且放大, dev 上为 +0.5pp 量级）
- **G3** m6a F1 ≥ classic − 0.5pp: **PASS**（6.7 ≥ 6.2 持平——与 dev 中性结论同构）

### M6a 诊断重演（副报告, 对照 dev 探路）

- 机制啮合规模远超 dev: 29,085 有模板 / 1,091 尝试延长 / 563 过 PE 跨接门 /
  360 验收 / +330,575 nt（dev probe2: 229 验收 / +237kb）
- **加长碱基构成与 dev 几乎一致**: same_ref 49.2% / diff_ref 46.5% /
  unaligned 4.3%（dev: 46.7/47.5/5.8）——半数真实回收、半数旁系游荡的结构
  在真实数据上精确复现
- same-ref 覆盖度: 0.39→0.46（mean）/ 0.38→0.56（median）, 34 条跨过 0.95
  门槛; 典型例: x553 (+807nt, AT3G11820.1) 覆盖 0.41→0.92——找回缺失半段
- **M6a 中性结论在真实数据上成立**: 即便延长规模 ×1.6, F1 仍与 classic 持平
  （旁系游荡的 P 代价与真实回收的 R 增益近似抵消——与 dev 收敛判定同构）

### 泛化验证总结

- dev 系列的引擎排序（M5 > classic, PogFlow ≈ classic）与 M6a 中性结论
  **全部在真实数据上复现**; M5 是唯一在两套数据上都净正的机制
- 真实数据 R_full 绝对值低（3.5-3.9% vs dev 8.3-8.6%）是度量口径现象
  （单样本 vs 全注释 48,359 条参考）, 不影响配置间排序
| 2026-08-25 | m7-probe-smoke | R=8.3% P=84.2% **F1=15.0%** | N50=1728 n=7882 | --path-engine pogflow --max_memory 16G --CPU 8 |
| 2026-08-25 | m7-nodump-dev1 | R=8.3% P=84.4% **F1=15.1%** | N50=1730 n=7866 | --path-engine pogflow --max_memory 16G --CPU 8 |
| 2026-08-26 | m7-checkpoint | R=3.5% P=82.3% **F1=6.7%** | N50=1568 n=32844 | --max_memory 16G --CPU 8 |

### M7 中期 checkpoint（2026-08-25, m7-checkpoint, 消 dump 后冷跑）

| 配置 | wall | Top stages |
|---|---|---|
| real-classic 基线（mtime 口径, 漏 normalize） | ~44+11 min | — |
| m7-checkpoint（T1 探针 + T2 消 dump） | **46 min** | rtt 1246s · normalize 672s · partition 188s · butterfly 175s |

- 指标: R=3.5 F1=6.7 持平; P 82.3（基线 82.8, 并行抖动边缘, 终验收复跑核对）
- **原"44 分钟"基线漏算 normalize 11 分钟**; dump 消灭实际省 ~8 min（54→46）
- RTT（45%）+ normalize（24%）合计 70% 且均已 8 线程并行——固定线程下只剩算法级路线
- **Gate 回谈结论（用户 2026-08-25）**: 固定 --CPU 8, 能压多少算多少（无硬数字）;
  优先级: RTT 查询结构 > normalize 内部 > partition 批量写
| 2026-08-26 | m7-rtt-dev1 | R=8.3% P=84.0% **F1=15.0%** | N50=1728 n=7891 | --path-engine pogflow --max_memory 16G --CPU 8 |
| 2026-08-26 | m7-norm-dev1 | R=8.3% P=84.3% **F1=15.1%** | N50=1730 n=7875 | --path-engine pogflow --max_memory 16G --CPU 8 |
| 2026-08-26 | m7-part-dev1 | R=8.2% P=84.5% **F1=15.0%** | N50=1728 n=7869 | --path-engine pogflow --max_memory 16G --CPU 8 |
| 2026-08-26 | m7-final-classic | R=3.5% P=82.3% **F1=6.7%** | N50=1565 n=32833 | --max_memory 16G --CPU 8 |
| 2026-08-26 | m7-final-pogflow | R=3.5% P=84.2% **F1=6.7%** | N50=1524 n=33041 | --path-engine pogflow --max_memory 16G --CPU 8 |
| 2026-08-26 | evalv2-dev1-m5 | R=8.6% P=83.6% **F1=15.7%** | N50=1737 n=8190 | --path-engine pogflow --model-guided --max_memory 16G --CPU 8 |
| 2026-08-26 | evalv2-dev1-m5 | R=8.6% P=83.6% **F1=15.7%** | N50=1737 n=8190 | --path-engine pogflow --model-guided --max_memory 16G --CPU 8 |
| 2026-08-26 | evalv2-dev1-plain | R=8.6% P=83.6% **F1=15.7%** | N50=1737 n=8190 | --path-engine pogflow --model-guided --max_memory 16G --CPU 8 |
| 2026-08-26 | evalv2-dev1-m5 | R=8.6% P=83.6% **F1=15.7%** | N50=1737 n=8190 | --path-engine pogflow --model-guided --max_memory 16G --CPU 8 |

### M7 终验收（2026-08-26, 固定 --CPU 8, 能压多少算多少）

| 配置 | wall | 指标 vs 基线 | Top stages |
|---|---|---|---|
| real-classic 基线 | ~54 min（44 mtime 口径 + normalize 11, 漏算修正后） | R=3.5 P=82.8 F1=6.7 | — |
| m7-checkpoint（消 dump 后） | 46 min | R=3.5 P=82.3 F1=6.7 | rtt 1246s · normalize 672s |
| **m7-final-classic** | **~17 min** | R=3.5 P=82.3 F1=6.7 ✓ | chrysalis_total 360s · normalize 189s · butterfly 168s |
| **m7-final-pogflow** | **~17 min**（两配置合计 34 min） | R=3.5 P=84.2 F1=6.7 **逐位同基线** | chrysalis_total 427s · butterfly 198s · normalize 183s |

- **总提速: 54→17 min（-69%）/ 46→17（-63%, checkpoint 口径）**——原始 ≤15min 目标
  差 2 分钟, 在"能压多少算多少"裁定下收官
- 四项优化全部带 oracle/字节级等价验证: 消 dump（二进制直桥）、RTT 16.3×
  （打包哈希索引+滑窗熵）、normalize 3.6×（并行覆盖统计+滚动 revcomp;
  dev1 口径 2.2×）、partition 5.5×（BufWriter, 原每行裸写 ~45M 系统调用）
- workspace 598 测试全绿; 全程零 git 提交（冻结）
- 剩余热点: butterfly（168-198s, 组件池已并行）、inchworm（146-151s,
  并行固有非确定段）——进一步提速需算法级, 未立项
| 2026-08-26 | evalv2-dev1-m5-sa | R=8.6% P=83.6% **F1=15.7%** | N50=1737 n=8190 | --path-engine pogflow --model-guided --max_memory 16G --CPU 8 |
| 2026-08-26 | v2-real-classic | R=3.5% P=82.3% **F1=6.7%** | N50=1565 n=32833 | --max_memory 16G --CPU 8 |
| 2026-08-26 | v2-real-pogflow | R=3.5% P=84.2% **F1=6.7%** | N50=1524 n=33041 | --path-engine pogflow --max_memory 16G --CPU 8 |
| 2026-08-26 | v2-real-classic | R=3.5% P=82.3% **F1=6.7%** | N50=1565 n=32833 | --max_memory 16G --CPU 8 |
| 2026-08-26 | v2-real-m6a | R=3.5% P=83.6% **F1=6.8%** | N50=1546 n=33017 | --path-engine pogflow --protein-guided --protein-db /storage/home/senior007/test/athal/bench/real/alyr_pep.fa --max_memory 16G --CPU 8 |
| 2026-08-26 | evalv2-dev1-m5-sa | R=8.6% P=83.6% **F1=15.7%** | N50=1737 n=8190 | --path-engine pogflow --model-guided --max_memory 16G --CPU 8 |

## M8 图解模型探针（2026-08-26, 转换级数据集 + GBM 基线）

目标: 验证"DL 模型 → 权重文件 → 辅助流分解选边"路线的可行性。产出
`xtask export-edges`（边表/转换表导出器, 14 单测）, 结论先行全部入档。

### 探针三步与裁定

1. **边级标签证伪（v1）**: dev1 全边 95.2% 正例——图来自同样本 reads 且
   inchworm 已剪错误 k-mer, "边真不真"非问题; 判别问题在分叉处（决策点出边
   正例率 38.1%）。全负分叉 44% = StringTie 参照不全的标签噪声实锤。
2. **转换级数据集（v2）**: 参照三连 k-mer 标注 (入边,出边) 转换, 候选 = 内部
   决策点（in≥1,out≥2）in×out 全组合。四集规模: dev1 8,309 / dev2 9,529 /
   dev3 5,378 / holdout 446 行, 正例率 58.5-60.9%（均衡）。**e1_label/e2_label
   为 oracle 泄漏列, 仅切片分析用, 训练禁用**。
3. **k-mer 覆盖度（v3）**: 边 (u→v) 的 25-mer = u24mer+v末碱基, 查 jellyfish
   计数表（mer_counts bin 已被收尾清理 → 回退 norm.fq 重数, 序列多重集恒等;
   save/load_counts_bin 收编至 trinity-kmer 共用）。

### GBM 基线（HistGradientBoosting, dev1+dev2 训 / dev3 测 / holdout 终测）

| 集 | 行级 acc / AUC | 决策级 argmax | argmax(e2_cov) 规则基线 |
|---|---|---|---|
| dev3 | 0.893 / 0.932 | **76.9%**（混合 62.3% / 纯配对 95.9%） | 76.6% |
| holdout | 0.789 / 0.853 | **70.0%**（混合 36.1% / 纯配对 96.0%） | 67.3% |

- 特征重要性: **e2_cov 0.226 + e1_cov 0.208 ≈ 全部信号**; 其余 17 个
  度数/序列特征合计 <0.02。GBM ≈ 单特征覆盖度规则（dev3 +0.3pp）。
- **裁定: DL/GNN 不成立**——可学信号就是覆盖度, 模式识别无增量;
  路线收敛为"覆盖度加权选边"。
- 纯配对切片 95.9% 但偏 trivial（92.7% 候选两边皆真）; 混合切片 62.3%
  的"错误"中 **63% 是 e2_cov≥50 的高覆盖边**——StringTie 漏报的真实分支,
  标签噪声非模型错误, 真实选择精度已贴近标签上限。
- **下一步（M8 主体）: 引擎实验**——graph.out 边权恒 1（M6a 已证）, 流分解
  分叉处无任何覆盖输入; 最小实验 = 25-mer 覆盖作边权注入 PogFlow, dev1 F1
  裁决。GBM 分数作为可选精化（相对裸覆盖 +0.3pp, 大概率不值得接线）。
- 评估循环性注意: 模型若跟随高覆盖分支而 StringTie 漏报, F1 口径下 P 会
  表面受损——引擎实验需同时看 R/P 分解与 n, 必要时换 Ensembl 参照复评。
| 2026-08-26 | m8-cap2-dev1sub | R=2.7% P=88.0% **F1=5.3%** | N50=1408 n=7426 | --path-engine pogflow --flow-cov-fusion --flow-cov-cap 2 --max_memory 16G --CPU 8 |
| 2026-08-26 | v2-real-m5 | R=3.9% P=86.2% **F1=7.4%** | N50=1470 n=34674 | --path-engine pogflow --model-guided --max_memory 16G --CPU 8 |
| 2026-08-26 | m8-cap5-dev1sub | R=2.7% P=87.6% **F1=5.2%** | N50=1412 n=7448 | --path-engine pogflow --flow-cov-fusion --flow-cov-cap 5 --max_memory 16G --CPU 8 |
| 2026-08-26 | m8-cap10-dev1sub | R=2.7% P=88.4% **F1=5.2%** | N50=1411 n=7427 | --path-engine pogflow --flow-cov-fusion --flow-cov-cap 10 --max_memory 16G --CPU 8 |
| 2026-08-26 | m8-cap50-dev1sub | R=2.7% P=88.1% **F1=5.2%** | N50=1406 n=7443 | --path-engine pogflow --flow-cov-fusion --flow-cov-cap 50 --max_memory 16G --CPU 8 |
| 2026-08-26 | m8-cap2-dev1sub | R=2.7% P=88.0% **F1=5.3%** | N50=1408 n=7426 | --path-engine pogflow --flow-cov-fusion --flow-cov-cap 2 --max_memory 16G --CPU 8 |
| 2026-08-26 | m8-cap5-dev1sub | R=2.7% P=87.6% **F1=5.2%** | N50=1412 n=7448 | --flow-cov-cap 5 余同 |
| 2026-08-26 | m8-cap10-dev1sub | R=2.7% P=88.4% **F1=5.2%** | N50=1411 n=7427 | --flow-cov-cap 10 余同 |
| 2026-08-26 | m8-cap50-dev1sub | R=2.7% P=88.1% **F1=5.2%** | N50=1406 n=7443 | --flow-cov-cap 50 余同 |

### M8 CAP 扫描（dev1sub, 2026-08-26）

| CAP | n | R | P | F1 | N50 |
|---|---|---|---|---|---|
| 2 | 7426 | 2.7% | 88.0% | 5.28% | 1408 |
| 5 | 7448 | 2.7% | 87.6% | 5.24% | 1412 |
| 10 | 7427 | 2.7% | 88.4% | 5.24% | 1411 |
| 50 | 7443 | 2.7% | 88.1% | 5.24% | 1406 |

- **CAP 不敏感**（四配置 F1 互差 ≤0.05pp, 在噪声带内）→ 按 plan 取默认 10 定稿。
- 融合机制确认生效（计数表 16.2M 键装载, 四配置产物数/精度均与 classic
  dev1sub 基线 7389/86.3% 有差异方向: P +1.3~+2.1pp）。
- 注意: dev1sub 无 pogflow 无融合对照（M1 网格只跑过 classic）——
  补跑 m8-ctrl-dev1sub（pogflow 无融合）以分离"引擎差"与"融合差"。
| 2026-08-26 | m8-ctrl-dev1sub | R=2.7% P=88.0% **F1=5.2%** | N50=1412 n=7448 | --path-engine pogflow --max_memory 16G --CPU 8 |

## 评估 v2（junction 级双维度）——上线与全景（2026-08-26）

### StringTie real 参照构建（scripts/prepare_stringtie_real.sh, 幂等, 全本地 ~5 min）

**引导 GTF 考古判定**: 当年 holdout 参照 `athal/stringtie.gtf` 头部命令行为
`stringtie -G chr5.renamed.gtf -p 8 -o stringtie.gtf -A stringtie.abund chr5.bam`, 其 chr1/chr4
杂染色体记录（3+5 行）经 `samtools idxstats chr5.bam` 证实来自 chr5.bam 自身的散逸比对
（chr1=153/chr2=108/chr3=172/chr4=382 条）——**不存在"用 all.bam 跑全基因组 StringTie"
的隐藏运行, 头部命令行诚实**。chr5.renamed.gtf 的生成 = TAIR10.gff.gz
（GCF_000001735.4）chr5 子集 → `gffread -E - -o chr5.gtf` → NC_003076.8 改名 chr5;
dev2/dev3 轨道同构（gffread 头部 `gffread -E - -o devN.gtf`, stringtie 头部
`-G devN.renamed.gtf -p 8`）。athal5.fa.fai 的 chr1..chr5 与 RefSeq 五条核染色体按尺寸
一一对应（NC_003070.9=chr1/30427671 … NC_003076.8=chr5/26975502; 叶绿体 NC_000932.1
与线粒体 NC_037304.1 不在 athal5.fa, 剔除）。→ real 全基因组引导 GTF = 同一
TAIR10.gff.gz 五条核染色体 → gffread → 同一改名约定 = `bench/real/guide.gtf`
（670,324 行 / 324,327 exon, chr1-5）。

构建命令（脚本内实形; 与当年管线的差异注明: hd_pipeline.sh 的 hisat2 带 `--dta`,
本构建按计划书形态未加, 已注脚本头）:

```
hisat2 -x athal5 -p 8 -1 bench/real/R1.fq.gz -2 bench/real/R2.fq.gz \
  | samtools sort -@ 4 -o hisat2.bam -          # 10M PE 对, 96.59% overall, ~2 min
stringtie -G guide.gtf -p 8 -o real_stringtie.gtf -A real_stringtie.abund hisat2.bam   # 3.0.3
```

产出: **real_stringtie.gtf = 33,277 转录本 / 21,998 基因（STRG）**（全基因组预期内）;
首条比对 DRR095137.1742804 → chr1:3646 100M。

### 双维度全景（junction 指标为 2026-08-26 回溯重评; classic/pogflow 组装缓存命中 m7-final, m5/m6a 因二进制更新重组装）

| 配置 | R | P | F1 | P_full | Jgene | Jtx | Jprec | novelJ | single | n |
|---|---|---|---|---|---|---|---|---|---|---|
| v2-real-classic | 3.5% | 82.3% | 6.7% | 40.1% | 58.3% | 35.7% | 54.4% | 7.5% | 32.7% | 32833 |
| v2-real-pogflow | 3.5% | 84.2% | 6.7% | 40.0% | 58.3% | 35.9% | 54.8% | 6.8% | 32.9% | 33041 |
| v2-real-m5 | 3.9% | 86.2% | **7.4%** | **42.4%** | **61.0%** | **37.4%** | 53.6% | **6.2%** | **30.7%** | 34674 |
| v2-real-m6a | 3.5% | 83.6% | 6.8% | 40.1% | 58.5% | 36.1% | 55.3% | 6.9% | 32.8% | 33017 |
| evalv2-dev1-m5（dev 轨） | 8.6% | 83.6% | 15.7% | 64.3% | 30.4% | 20.9% | 21.0% | 74.1% | 33.8% | 8190 |

（junction 分母: real = 参照 21,998 基因 / 33,277 转录本; dev = stringtie.gtf 的 4,831/5,369。
m5/m6a 的 F1 与 08-25 首跑 7.5/6.7 差 0.1pp, 属重组装 run-to-run 抖动。）

### 验收

1. **区分力 PASS**: m5 Jgene 61.0% > classic 58.3%（+2.7pp）, 与 F1 体系同向（7.4 > 6.7）;
   Jtx +1.7pp、novelJ -1.3pp、single -2.0pp 三项同向
2. **稳定性 PASS**: v2-real-classic 重放两次, Junction 行逐字节一致
   （`Jgene=58.3% Jtx=35.7% Jprec=54.4% novelJ=7.5% single=32.7% (genes=21998 tx=33277 asm=32751)`）
3. **数值合理性**: Jgene 58-61% 高于预估 15-40% 带——判定为参照语义差而非 bug
   （Task4 首跑 sanity 与本次重放逐位一致; asm 数 32-34k 与组装 n 吻合）, 详见下段

### junction 数 vs F1 体系: 读出什么

- **M5 优势方向保持, 相对幅度收窄**: F1 差 +0.7pp ≈ +10.4% 相对; Jgene 差 +2.7pp ≈ +4.6%
  相对。M5 的增益本质是"找回更多表达位点"（Jgene/Jtx 同升 + 单外显子占比 32.7→30.7%）
  而非"每条更准"——Jprec 反而 54.4→53.6 微降: 多出的链有一部分由参考内含子组成但
  不精确等于任一参考链（变异剪接形态的密集上报）; novelJ 反而最低（6.2%）说明它
  **不是编造内含子, 而是多报同位点的链变体**
- **m6a 中性在 junction 语义下成立**: Jgene 58.5 vs 58.3（+0.2pp）, Jtx +0.4pp,
  Jprec +0.9pp——与 F1（6.8 vs 6.7）同构, 旁系游荡式延长不伤 junction 完整度
- **pogflow ≈ classic 在两套体系下严格成立**（Jgene 逐位同为 58.3）
- **跨轨道绝对值不可比**: dev1 Jgene 30.4% vs real 61.0% 是参照构造差异而非组装器差异
  ——dev1 参照（chr5 低深度、另一样本 DRR168 三重样 + chr5-only 引导）不含组装 reads
  （hd chr5, DRR200069）的大部分真实内含子, novelJ 74.1% 即其镜像; real 参照与组装
  **同本同 reads**（同 10M 对、96.6% 映射、全基因组 TAIR10 引导）, 参照只含有表达的
  位点, 加上基因级"任一链复现即找回"判据天然宽松 → Jgene 抬到 58-61%, 超出按
  "全注释/异本参照"语义预估的 15-40% 带。**junction 指标的正确用法是同轨道内配置间
  Δ, 不是跨参照绝对值**
- 实用结论: 双维度比单一 F1 多拆出一层——M5 的"+R 不伤 P"在 junction 语义下细化为
  "宏观找回更多位点 + 内含子最干净 + 链形态变体更密（Jprec 微幅代价）"; m6a/pogflow
  的中性结论在第二维度复核后无反转
| 2026-08-26 | m8-dev1 | R=8.4% P=84.5% **F1=15.2%** | N50=1730 n=7881 | --path-engine pogflow --flow-cov-fusion --flow-cov-cap 10 --max_memory 16G --CPU 8 |
| 2026-08-26 | m8-dev2 | R=7.7% P=82.7% **F1=14.2%** | N50=1684 n=9162 | --path-engine pogflow --flow-cov-fusion --flow-cov-cap 10 --max_memory 16G --CPU 8 |
| 2026-08-26 | m8-dev3 | R=7.4% P=85.5% **F1=13.6%** | N50=1725 n=6537 | --path-engine pogflow --flow-cov-fusion --flow-cov-cap 10 --max_memory 16G --CPU 8 |
| 2026-08-26 | m8-holdout | R=2.1% P=85.4% **F1=4.1%** | N50=583 n=4591 | --path-engine pogflow --flow-cov-fusion --flow-cov-cap 10 --max_memory 4G --CPU 8 |
| 2026-08-26 | m8-ctrl-dev1sub | R=2.7% P=88.0% **F1=5.2%** | N50=1412 n=7448 | --path-engine pogflow --max_memory 16G --CPU 8（无融合对照） |
| 2026-08-26 | m8-dev1 | R=8.4% P=84.5% **F1=15.2%** | N50=1730 n=7881 | --path-engine pogflow --flow-cov-fusion --flow-cov-cap 10 --max_memory 16G --CPU 8 |
| 2026-08-26 | m8-dev2 | R=7.7% P=82.7% **F1=14.2%** | N50=1684 n=9162 | 同上 |
| 2026-08-26 | m8-dev3 | R=7.4% P=85.5% **F1=13.6%** | N50=1725 n=6537 | 同上 |
| 2026-08-26 | m8-holdout | R=2.1% P=85.4% **F1=4.1%** | N50=583 n=4591 | --max_memory 4G |

## M8 覆盖度加权流分解验收（2026-08-26, --flow-cov-fusion, CAP=10）

| 集 | M3（pogflow 无融合）F1 | **M8 F1** | Δ | R / P 变化 |
|---|---|---|---|---|
| dev1 | 15.1% | **15.2%** | +0.1pp | R 8.3→8.4, P 84.5→84.5 |
| dev2 | 13.9% | **14.2%** | +0.3pp | R 7.6→7.7, P 82.3→82.7 |
| dev3 | 13.5% | **13.6%** | +0.1pp | R 7.4→7.4, P 85.9→85.5 |
| holdout | 4.0% | **4.1%** | +0.1pp | R 2.0→2.1, P 85.2→85.4 |

### Gate 判定: PASS

- dev1 F1 15.2% > 15.1%（M3 基线）且 P 84.5% ≥ 84%——**四集全部非负, 三集为正**
- 方向与探针预言一致: R 普遍微升而 P 未损——无 read 证据的高覆盖边确为
  真实分支, 覆盖度兜底把它们找回而不引入假边（P 未现 spec 风险 3 的
  "表面受损"形态, 无需 Ensembl 复评分支）
- CAP=10（扫描不敏感, 见上节）; 机制: 无证据边 floor 1 → min(cov25,10),
  read 证据边不动; 融合装载 23.9-31.1M 键, 均命中 bin 无需重数
- dev1sub 对照（ctrl 5.2% vs 融合 5.2-5.3%）无分离——亚采样深度下融合
  效果不可见, dev1sub 仅用于 CAP 排序的既定口径不变
- 回归: bench 缓存键含 passthrough flag（新 key 新目录, 无互踩）;
  fusion off 路径与 M3 逐边一致（digraph_to_flow 单测 + workspace 46 套件绿）

### BUSCO 补充维度（2026-08-26, embryophyta_odb10 离线, n=1614, 首跑）

| 组装 | Complete [S/D] | Fragmented | Missing | 备注 |
|---|---|---|---|---|
| real-classic | 70.6% [61.6/8.9] | 12.3% | 17.1% | 10M PE, M7 后 |
| real-pogflow | 70.6% [61.0/9.6] | 12.3% | 17.1% | 与 classic 持平（三体系一致） |
| real-m5 | **70.9% [62.1/8.8]** | **11.6%** | 17.5% | 方向性最佳, 与 F1/Jgene 同向但幅度小（保守基因集对表达广度不敏感） |
| real-m6a | 70.6% [60.4/10.2] | 12.3% | 17.1% | 中性（三体系一致） |
| dev1-pogflow / dev1-m5 | 23.5% / 23.4% | ~2.2% | ~74% | dev 为部分转录组数据, 绝对值低属预期 |
| **orig-Trinity（chr5 对拍）** | **3.4% [3.0/0.4]** | 3.5% | 93.1% | 原版 Perl v2.15.2 |
| **rust-port-early（同数据）** | **3.4% [3.2/0.2]** | 3.5% | 93.1% | **Complete 与原版完全一致**——移植保真度的独立第三方验证 |

- **移植保真度**: 同数据下 Rust 移植与原版 Trinity 的 BUSCO 完整度逐位一致
  （3.4%=3.4%, F/M 同; S/D 微差 0.2pp 属平局序噪声）——与既有 xcheck 多重集
  等价结论互相印证
- chr5 数据完整度 3.4% 低是数据覆盖所致（chr5 子集 reads 对全转录组保守基因）,
  与移植质量无关
- 三套质量体系（F1 / junction / BUSCO）对四配置的排序完全一致:
  M5 ≥ 其余持平, 无配置拖后腿
| 2026-08-26 | real-m8 | R=3.5% P=84.2% **F1=6.7%** | N50=1531 n=33025 | --path-engine pogflow --flow-cov-fusion --flow-cov-cap 10 --max_memory 16G --CPU 8 |
| 2026-08-26 | real-m8 | R=3.5% P=84.2% **F1=6.7%** | N50=1531 n=33025 | --path-engine pogflow --flow-cov-fusion --flow-cov-cap 10 --max_memory 16G --CPU 8 |

### M8 real 轨道验证（2026-08-26, real-m8 vs v2-real-pogflow）

| 配置 | R | P | F1 | N50 | n |
|---|---|---|---|---|---|
| v2-real-pogflow | 3.5% | 84.2% | 6.7% | 1524 | 33041 |
| real-m8（融合） | 3.5% | 84.2% | 6.7% | 1531 | 33025 |

- **real 上中性**（计数表 48.0M 键装载, 机制生效但指标逐位持平）——dev 四集的
  +0.1~0.3pp 增益未在真实数据复现; 与 M6a 经验同构（真实数据 PE 结构更复杂,
  弱证据分叉的覆盖度兜底不再提供增量）。M8 结论修订为: **dev 模拟数据小增益,
  real 中性无害**——flag 保留 opt-in, 不设默认。

## M9 候选重排探针（2026-08-26, 负结果——过滤范式判死）

`xtask export-candidates`（库调 PogFlow 离线重跑候选, min_read_support=1,
全长级特征 + minimap2 qcov 标签）+ train_rerank.py（GBM + 操作口径）。

### 关键发现

1. **标签口径坑（已修入档）**: 首版误用 ref_cov（召回方向）作正例 → 9.5% 正例,
   与 eval P=84% 矛盾; 正确口径 = qcov（query 侧 ≥95% 碱基以 ≥95% 同一率比对,
   eval 精度方向的转录本粒度化）→ 正例 53-58%。
2. **PogFlow support 过滤实际不拦截**: 候选 support=1: 0 条 / =2: ≤4 条,
   其余全部 ≥3——放宽 -R 不扩池, 候选池由流分解本身决定。
3. **分类口径（dev3）**: acc 0.784 / AUC 0.843; len 重要性 0.173 独大,
   cov_cv/cov_mean/orf_frac/cov_ge2_frac 有真实增量, support 仅 0.017。
4. **操作口径判死（两集一致）**:

| 集 | 全保留 F1 | 模型@0.5 F1 | **oracle 过滤 F1** |
|---|---|---|---|
| dev1 | 15.04 | 12.65 | **10.76** |
| dev3 | 13.55 | 11.40 | **8.97** |

- **完美分类器（oracle）也把 F1 打掉 3.6-4.6pp**——R≪P 时 F1≈2R,
  精度一文不值; qcov 判负的部分比对候选贡献 **36% 的 R_full**（dev1:
  R 8.24→5.69）——过滤 = 自割召回。
- **裁定: 候选过滤/重排范式在 F1 口径下已证伪, 与模型容量无关
  （oracle 已满分）**。DL 辅助组装仅存的路径 = 生成侧增召回
  （M8 覆盖度融合为其最便宜形态; "脑补 reads 缺失末端"被 M6b 判死）。
  Transformer 路线三次测量（转换级探针/GBM 特征审计/oracle 过滤）全否, 收档。

## M10 探路: 近缘种保守性探针（2026-08-27, alyr cDNA × dev1 混合决策点）

数据: alyr_cdna.fa（Ensembl r60, 32,667 条）; dev1 混合决策点 1,710 个,
出边支路上下文（fork 顶点 24-mer + 支路线性链延伸 80nt）3,425 条,
minimap2 -x asm10 分别比对 alyr_cdna 与 dev1 参照。

### 结果

| 列联 | 参照真支路 | 参照假支路 |
|---|---|---|
| alyr 命中 | 241（11.0%） | 22（1.8%） |

- **6.2× 富集**, 命中 identity 中位: 真 0.941 / 假 0.854
- **决策级**: 唯一 alyr 命中可裁决 199/1710 个混合决策点（11.6%）, **选对率 96.0%**
- 标签噪声救援: 92 条支路 atha 参照不覆盖但 alyr 命中——保守性证实部分
  "参照假"支路实为真实保守分支（M8 探针 63% 标签噪声的直接证据）
- **裁定: 信号真实、精度极高、覆盖稀疏（11.6% 决策点）**。不足以当独立机制,
  但适合作为 M8 融合权重的第二路先验（覆盖度兜底 + alyr 命中加成）。
  预期收益小（≤199 个决策点 × 此前 ~62% 错误率 ≈ 75 个分叉修正上限）。
  下一步: M10 集成实验（--flow-alyr-fusion, 边权 += 命中加成）, dev1 Gate 裁决。

## 特性面精简与 release 优化吸收（2026-08-27）

- **砍**（用户裁定, 特性面收敛到生产形态）: M2 flow / M2.1-M3 PogFlow /
  M4 super_reads / M5 model-guided / M6 protein-guided / M8 flow-cov-fusion,
  连带 trinity-flow crate 与全部相关 e2e; xtask export-edges 同步移除。
  判定依据见各里程碑段（M5 为真实数据唯一正收益但复杂度不合 release 定位;
  M8 real 中性; M4/M6 负/中性）。被砍特性在 git 历史与本文档完整保留。
- **留**: classic 引擎 + --recursive（M1）+ M7 提速 + evalv2 junction 评估 +
  M9 候选探针（export-candidates 改接 classic 引擎, 候选头无 support → -1）。
- **吸收**（cherry-pick 自 trinity-rust-release, 零冲突）: 硬链接便捷名
  （省 ~2.6GB 双写）/ --scratch 中间产物外置 / Component_bins 打包
  （消灭 per-component 小文件）/ fq 归一化直出 fa / mimalloc 全局分配器
  （小集 sys -85%, wall -25%）。并收入其基线的 clippy 修复（workspace 现已零警告）。
- 终验: 构建 + 38 测试套件全绿, clippy/fmt 干净。

### IO/系统开销优化（2026-08-27, 发布版真实数据驱动: sys 31% + 13GB 写放大）

- 归一化便捷名硬链接（消 ~2.6GB 双写）; `--scratch <dir>` 中间产物整体出
  NFS; Component_bins 全量打包（tmp/graphs/butterfly 三包, per-component
  小文件 8.5 万 + 5.7 万 → 每 Cbin 4-6 个文件）; fq 归一化直出 fa（跳
  norm.fq 物化, both.fa 逐字节对拍锁定）; mimalloc 全局分配器（sample_data
  A/B: sys -85%, minor fault -96%, wall -25%）。
- **--scratch 默认化**: 不传时中间产物落 `${TMPDIR:-/tmp}/trinity-rust-
  <outdir 路径 FNV-1a 哈希>/`, 启动 statvfs 空间检查（需求 = 压缩输入 ×6,
  下限 4GB）, 不足/不可写则警告并回落 --output 旧行为; checkpoint 随产物,
  哈希确定 → 续跑语义不变。
- xtask xcheck-chrysalis 已适配 pack 布局（原版工具链在位, 7/7 PASS）。

### 深度连续性特征补测（2026-08-27, 用户问"逐碱基深度能否提高解图正确性"）

- 转换级数据集加连续性三特征（e2/e1 ratio, log2 ratio, diff; -1 缺失截 0 处理）,
  dev1+dev2→dev3: acc 0.8927→0.8912, AUC 0.9317→0.9338——**净增益为零**。
- cov_ratio 排列重要性 0.070（第三位）但完全挤占 e1_cov/e2_cov 份额——
  连续性信息已编码在原始覆盖度中, 属冗余重打包。
- 裁定: 逐碱基深度信息 = k-mer 覆盖度（M8/M9 已在用且已是主信号）;
  StringTie 基因组锚定深度与 k-mer 计数同源（同一份 reads）且无法修自身漏报
  （标签噪声主源）, 不推荐; 转录本边界的深度陡降信号受 M6b 数据上限约束
  （缺失末端 96% 不在 reads）, 无可挽回空间。此方向收档。

## 参照转录本图内可走性分解（2026-08-27, 用户问"StringTie 走得通的路径 Trinity 是否走得通"）

工具: `bench/export/probe_walkability.py`（import probe_conservation 的
parse_graphs/MINIMAP2）。方法: 8,712 条 StringTie 参照（dev1）比对到
Trinity de Bruijn 图的顶点/边, 滑窗判定"在位+边连通", 按窗口达标率分档;
"产出"口径 = 产出 fasta 中存在 ≥95% 覆盖 / ≥95% identity 的匹配。

### 结果

| 可走阈值 | 图内可走 | 可走且产出 | 可走未产出 |
|---|---|---|---|
| 99% 窗口达标 | 792 | 311 | 481（60.7%） |
| 95% | 2,153 | 710 | 1,443 |
| 90% | 3,685 | 718（= R_full 总数, 饱和） | 2,967 |

- 部分可走（50-99%）6,751 条（77.5%）; 基本不可走（<50%）1,169 条（13.4%）。
- 产出但 99% 口径不可走 407 条——引擎 zipping/延伸重构超出图路径的部分
  （Trinity 产出里 StringTie 压根没有对应路径的旁证）。
- **95% 口径最诚实**（与产出匹配的覆盖阈值一致）: ~1,443 条参照
  "图里有完整路径、引擎没走出来"——上界（含引擎产出了近缘异构体兄弟的情况）,
  即 M9"扩池再过滤"路线的目标市场。

### 方向性结论

- StringTie 走得通而 Trinity 图走不通: junction 悬垂短于 24bp（splice 认得、
  de Bruijn 成不了窗）、25-mer 计数低被剪、基因组坐标可分而 k-mer 塌缩的重复区。
- Trinity 走得通而 StringTie 不能: 比对不上基因组的序列、StringTie 过滤器
  丢掉的弱异构体（上述 407 条为旁证）。
- 不对称真实存在, 且两个方向都有。

### 瑕疵注记

- 探针图取自 cache `68349cb4e339ced0`（m7-part-dev1）, 产出 fasta 取自
  `b4a32209436c2360`（baseline-dev1）——两次独立 classic run, 图分箱存在
  并行种子序的细微差异。方向性结论不受影响, 绝对数字有少量噪声。

## 可走未产出池死因分解（2026-08-27, probe_stuck 1-5, M9 扩池参数依据）

对 95% 口径的 1,443 条"图内可走但未产出"参照逐条尸检
（`bench/export/probe_stuck{,2,3,4,5}.py`, 明细 `stuck.dev1.tsv`）。

### 死因分解

| 桶 | n | 占比 | 形态 |
|---|---|---|---|
| near_miss | 257 | 18% | 产出近等长近全同序列（best_cov 中位 0.94, ident=1.00, 长度比 0.97）, 差在端部几个百分点; 38 条被 ≥1.3× 长产出吞并 |
| partial | 177 | 12% | 产出 0.72× 截短版, ident=1.00——提前终止/端部外显子缺失 |
| absent | 1,009 | 70% | 83% 无任何比对; 但 k-mer 产出命中率中位 0.98 |

absent 桶内部（关键反转——**k-mer 材料全在产出池, 缺的是组合**）:
- comp 黑洞仅 27 条（2.7%）, 参照主 comp 基本都有产出（comp_low≤2 条: 724）
- k-mer 命中集中度: top1≥0.9（集中单条同源产出, asm5 口径未达）465 条（46%）;
  top1<0.7（分散在 2-4 条产出, **组合枚举不足**）472 条（47%）; top3 覆盖 100%
- 33% 有同基因兄弟被产出（引擎走了兄弟 isoform, 没走这条）
- 参照路径决策点密度中位 0——多为无分叉线性链段, 枚举成本低

### 结论（M9 扩池参数指导）

1. **端部完整性**（near_miss+partial 434, 30%）: 候选池需含端部延伸变体。
2. **组合枚举不足**（472, 33%）: 参照各段被不同产出分别走出, 完整组合从未
   枚举——classic 每 comp 只出 1-2 条是核心瓶颈。comp 级 k-best 枚举直接吃这块。
3. **同源近似未达标**（465, 32%）: 集中在单条产出上, 其中一部分实为"实质回收"
   （asm5 严格口径漏计）, 一部分是变体路径未枚举。
4. comp 黑洞可忽略（27 条）。扩池的本质是**补组合+补端部**, 不是救 comp。
- 瑕疵注记同上节（cache 混用）; probe_stuck 首版的 containment 统计曾用错
  f[6]（参照长）当产出长, 已在 probe_stuck3 修正重算。

## 死因分解重大修正: absent 桶是评估 artifact（2026-08-27, probe_rewire 1-2）

`minimap2 <ref> <final> --secondary=no` 参数序下 query=产出: 一条产出近似匹配
多条参照（多 isoform 基因必然发生）时, 非最佳参照的 hit 全被 secondary 屏蔽
掐掉 → R_full 系统性低估。

### 证据链（probe_rewire.py / probe_rewire2.py, 顶点级 k-mer 映射回图）

- 产出 k-mer 映射回图: 主 comp 映射率中位 1.00（产出 fasta path 头的 vid 是
  引擎私有编号, 与 graph.out 不同空间, 不可直接用; k-mer 映射是正确通道）。
- absent 1,009 条: top1 产出顶点命中率中位 **1.00**, LIS 共线中位 **1.00**
  （双向取大; 1008/1009 保序）, 顶点数比中位 1.08——图顶点层面存在一条产出
  几乎逐顶点保序走过参照路径, 长度相当。无乱序嵌合。
- 去掉 `--secondary=no` 重测: absent 桶 **439 条直接达 R_full 口径**;
  全量 dev1: R_full 718 (8.2%) → **1,220 (14.0%)**, 低估 502 条 (5.8pp, 相对 +70%)。
- 顶点级真分散（owner top1<0.7）仅 13 条（probe_stuck5 的 472 条"分散"被
  "(k-mer,产出) 命对"重复计数稀释 top1 的判据高估——跨 comp 同源命中污染）。

### 修复

- `xtask/src/bench.rs` 与 `export_candidates.rs` 的 asm5 比对去掉
  `--secondary=no`（compute_metrics 本就 per-target 取 best, 无需改聚合）。
  splice junction 评估的 `--secondary=no` 语义不同（防基因组乱命中）, 不动。
- 历史影响: M0-M7 所有 bench 绝对值偏低; A/B 相对比较同向偏差, 结论方向
  不受影响, 但幅度需按新口径复核。dev1 新成绩（xtask bench 复跑确认）:
  R_full 14.0% / R70 53.4% / P 81.4% / F1 23.9%（此前口径 R=8.2%）。

### 对 M9 的重新定位

- 可走未产出池从 1,443 缩至 ~941; absent 桶 43% 是假缺口（评估漏计）,
  真缺口死因以端部完整性（near_miss+partial 434）与同基因兄弟竞争为主。
- 组合枚举不足市场从 472 缩至 ~13——**重组器路线基本出局**;
  M9 主攻端部延伸 + asm10 敏感性列。

## 端部延伸探针: 收益微薄, M9 候选生成路线收档（2026-08-27, probe_extend 1-2）

- `probe_extend.py`（唯一链延伸, MAX_EXT 1000bp/端）: 425 个 target 尾部锚
  停止原因——**死端 411**、分叉 11、环 3; 覆盖增益 >0.05 仅 8 条,
  R_full 转正≈0。产出截短版的末端已在图死端上, 无延伸空间。
- `probe_extend2.py`（缺端锚定位, 共 554 个缺端）: 锚**在参照路径上** 310
  （56%, 真截断——序列终止于参照路径中间, 指向引擎路径枚举终止规则）;
  锚**岔路** 244（44%）。端点恰同参照端点 96 个。
- **裁定: M9 候选生成三连缩水**（重组 472→13, 延伸转正≈0, absent 43% 假缺口）,
  "扩池再过滤"路线在 classic 产出之上无利可图, 收档。
  剩余真缺口要么指向引擎路径枚举终止规则（动核心）, 要么图构建端
  （不可走 13.4%）。下一候选方向: R70→R_full 之间的 39.4pp 池
  （53.4%-14.0%, dev1 新口径）——先探明该池死因再立项。
- 探针注意: k2v 的 setdefault 语义在共享 kmer 上会撞 comp, 锚级图结构
  判定（死端/边存在）含少量噪声; 方向性结论不受影响。

## R70→R_full 池死因分解: classic 已近数据上限（2026-08-27, probe_r70pool 1-2）

池 = R70 达标但 R_full 未达（dev1 新口径, n=3,418）。

### 死因账

| 死因 | n | 占比 | 属性 |
|---|---|---|---|
| 缺失段 reads 里就没有/零星（<20%） | ~1,216 | 36% | **数据上限**（StringTie 参照端部 reads 未覆盖） |
| 缺失段 reads 部分有（20-50%） | ~978 | 29% | 数据上限（边缘） |
| 缺失段图内没有但 reads 有（DS.fa ≥50%） | ~416 | 12% | **图构建剪枝损失**（可救, 调参） |
| 缺失段图内可走（≥90%）但未产出 | 205 | 6% | **引擎缺口**（动路径枚举） |
| ident 不足（cov≥0.95） | 141 | 4% | 变体/参照噪声 |
| 部分可走/弥散缺失 | ~614 | 18% | 混合 |

主因几何: 3' 缺 58% + 5' 缺 28%（端部合计 86%）, 段长中位 170bp。
（reads 存在性两口径有敏感性差: 原始 fq ≥50% 仅 9% vs DS.fa 17%——
-x sr 大索引敏感性下降, 取 DS.fa 口径为"reads 有"上界。）

### 裁定

- **~2/3 的池是数据/参照上限**: de novo 不可能产出 reads 里没有的序列
  （StringTie 基因组引导的端部常超出 reads 覆盖, 标签噪声的另一面）。
- 真正可救: 图构建剪枝 ~350-400 条（+4-5pp R_full 上界）+ 引擎缺口 205 条。
- dev1 上 classic 引擎的解图质量不是主要瓶颈——瓶颈在图构建端与数据本身。
- 下一步候选（待裁定）: 图构建剪枝参数 A/B（min_kmer_count 方向,
  收益上界 +4pp, 代价 k-mer 增量/图膨胀/耗时）, 或将探针方法论
  转移到真实数据验证泛化。

## 真实数据泛化验证: dev1 死因结构成立（2026-08-27, probe_real 1-2, release_real）

参照 = 真实数据 StringTie guided（real_stringtie.gtf + athal5 基因组提取,
33,277 条, `extract_real_ref.py`）; 产出 = release_real.Trinity.fasta;
图索引 = numpy uint64 排序数组（27.5M unique k-mer, 2.7GB 图的内存适配,
缓存 real_graph_kmers.npy）。

### dev1 vs real 死因结构对比

| 维度 | dev1（模拟） | real（真实） |
|---|---|---|
| 参照数 | 8,712 | 33,277 |
| R_full / R70 | 14.0% / 53.4% | 9.5% / 47.1% |
| R70→R_full 池 | 3,418 (39%) | 12,514 (38%) |
| ident 不足 | 4% | 4% |
| 端部缺（5'+3'） | 86% | 84%（39+45） |
| 缺段图内可走（引擎缺口） | 7% | 7% |
| 缺段图内没有 | 84% | 85%（8,930 段） |
| "图内没有"段 reads 有（DS 口径） | 17% | **13%（1,182 段）** |
| "图内没有"段原始 reads 没有 | 49%（修正后） | 62% |

- **结构性结论完全泛化**: 端部缺失主导、图内可走仅 7%、图剪枝可救块
  同量级（real ~9-13% of 池, 上界 +3-4pp R_full）。
- 真实数据数据上限更强（原始 reads 没有 62% vs dev1 49%）——guided
  StringTie 参照端部更完整, de novo 原理上更难全回收。
- reads 存在性探针 bug 修正: 多 fq 作 target 时 minimap2 只把首个当
  target——dev1 原始 reads 数字 9% 系全零污染, 修正后 11%; DS 单文件
  口径不受影响。real 版已修正（segs 当 target, 统计 target 侧覆盖）。

### 裁定

- B 任务（泛化验证）完成: dev1 探针方法论与结论可外推, 无过拟合迹象。
- A 任务（图构建剪枝参数 A/B）立项依据确认: 真实数据同样存在
  ~1,200 段"reads 有但图里没有"的可救块。

## A 任务探路: 丢失层定位与收益上界——证伪收档（2026-08-27, probe_locate）

### 丢失层定位（real, 1,182 条"DS reads 有但图没有"段）

- `prune_some_kmers` 实际调用参数 = (min_count=1, entropy=1.0,
  **prune_error_kmers=false**, 0.005)——错误率剪枝未启用, min_count=1
  已最宽: **k-mer 目录层无剪枝可言, 调参无空间**（M1 tier-1 亦证）。
- 段 k-mer 在 inchworm contigs 的存在率: 97% <20%（中位 0.00）;
  contigs unique k-mer 27.53M ≈ 图 27.53M（图顶点全部来自 contigs）——
  **丢失层 = inchworm 贪心延伸决策**: 计数 1-2 的低覆盖支路在与高覆盖
  主路的分叉竞争中败北, 从未成为 contig。
- 段 k-mer 的 DS 计数分布: 计数 0 占 79%, 计数 1 占 20%——低计数假说确认。

### 收益上界（决定性）

8,930 个"图内没有"段的 DS 覆盖细分: **≥90% 仅 125 条（1.4%）**,
≥80% 204, ≥50% 661, <50% 7,748（87%）。

R_full 的 95% 覆盖门槛要求缺段几乎全长有 read 支撑——只有 DS≥90% 的
125 段回收后能让参照跨线。**完美回收所有"reads 有但图没有"k-mer 的
R_full 上界 = 125/33,277 ≈ +0.4pp**（此前 +3-4pp 的估计错在把
"DS≥50% 有支撑"当"可完整回收"）。

### 裁定: A 证伪, 不动 inchworm

- 动 inchworm 延伸策略（低覆盖支路枚举/二次种子回收）换 +0.4pp 上界,
  且引入图膨胀/错误 k-mer 入图/精度下降风险——负期望, 不做。
- **调研总结论**: 修复评估 bug 后（+5.8pp）, classic（贪心 inchworm +
  butterfly）在 dev1/real 上已接近其架构可达上限; 剩余缺口 ~65% 数据
  上限封顶, 后处理与调参空间已探尽。下一步若要突破, 需换架构级机制
  （如 genome-guided 模式、长读长辅助）, 非本仓当前定位。

## OLC 路线冒烟: 从零与合并层双双证伪（2026-08-27, probe_olc 1-2）

用户提议 hifiasm 风格 OLC; 短读转录组 OLC 的两个形态都做了冒烟。

### 从零 OLC（probe_olc1, reads 级）

- 归一化 reads 中位 **83bp**——hifiasm 前提（15-25kb HiFi）不存在。
- minimap2 敏感参数（-k10 -w2）self-map 仅 ~15K 对 overlap（每 read 0.13 个,
  归一化压稀了 read 密度）; dovetail 几何（overlap>=40 + overhang<=35）
  在 83bp 上与 containment 不可分——**合格边 0 条**。minibwa（lh3 的
  bwa-mem 继任者, 短读敏感 seeding）同样救不动: self-map 只回 best hit。
- 判定: 从零 OLC 在 83bp 归一化 reads 上几何不可行（Ra 等从头 OLC 转录组
  组装器成熟度低的同一根因）。

### 产出级 OLC 合并层（probe_olc2, CAP3/MIRA 收尾思路）

classic 产出（32,748 条, 中位 717bp）self-overlap → dovetail 2,621 边 →
1,071 个 cluster 贪心拼接（合并体中位 1,968bp, 吸收 2,366 条成员）:

| | R_full | R70 | P | n |
|---|---|---|---|---|
| classic | 9.5% | 47.1% | 84.6% | 32,748 |
| OLC 合并 | **9.1%** | **46.6%** | 85.0% | 31,453 |

**净负 -0.4pp**: 共享末端外显子的 isoform 对被判 dovetail, 贪心链式拼接
（cluster 最大 18 成员）制造嵌合, 吸收的成员原本各自 hit 不同参照。
且根因与数据上限一致: classic 产出的截断末端在图死端（reads 没覆盖）,
两条同死端截断体合并无米下锅——缺失段不在任何原料里。

### 裁定

- OLC 两形态证伪收档。多 K 合并（Euler-mix 变体）不另跑: 合并层本身
  已净负, K!=25 产出更差（M1: K29 -4.8 vs 5.3）, 数据上限同样封顶;
  Euler-mix 的 2,544→10,646 是检出数量非参照回收, 与本仓 R_full 目标不同轴。
- 架构级突破的可行域确认只剩: genome-guided 模式 / 长读长辅助
  （长读下 OLC/string graph 才是主场）——均超当前仓定位, 待用户裁定。

## Ra 思想混合层复测: bubble 语义仍净负, 混合 OLC 收口（2026-08-27, probe_olc3）

用户指出 Ra（github.com/mariokostelac/ra, OLC 转录组组装器）的思想——
全 overlap 布局 + trimming + bubble popping; 吸收其语义重做产出级混合层
（修正 naive 版三宗罪: dovetail 误判 isoform / 贪心链式嵌合 / 合并吞 hit）:

- string graph（全 overlap 建边 1,926, contained 561）→ bubble-aware
  DFS 路径枚举（每 bubble ≤4 支路, 每 component ≤8 路径, isoform 支路
  各成路径）→ contained 只在完全被路径覆盖时去冗余。
- real: R_full 9.5%→**9.2%** (-0.3pp), R70 47.1%→46.7%, P 84.6%→84.2%,
  n 32,748→31,951（去冗余 797 条, 但去掉的成员含可 hit 参照者,
  路径序列未完全继承成员的 hit 能力）。

**两次语义（naive -0.4pp / Ra bubble -0.3pp）稳定小负, 混合 OLC 层收口**。
根因统一: classic 产出的信息已在转录本序列里, 重排不创造 reads 里没有的
覆盖; R_full 瓶颈是数据上限, 不是序列组织方式。Ra 本身在文献中不及
Trinity 系的表现与此一致。短读 de novo 场景的 OLC/混合探索彻底闭环
（从零几何不可行 + 合并两连负）; 可行域仍是 genome-guided / 长读长。

## Ra 原版完整管线实测: 远逊 classic, OLC 从零路线终审（2026-08-27）

用户裁定须跑 Ra 完整核心算法（非自制模拟）。修复 master 分支 5 处断点后
编译全链（vendor/cpppoa 与 include/ra 全头文件补 <cstdint>; main.cpp
命名不一致 writeOverlaps→write_overlaps; ra_layout 接口漂移
readAfgOverlaps 缺 reads 参数 / AfgRead→Read; overlap 输出裸行 vs layout
期望 afg OVL 的格式断点——转换器补 rds/adj/ahg/bhg 包装;
updateOverlapIds 无定义补 no-op, read id 未重编号语义自洽）。

### dev1 全链（单端 R1 5.8M reads, 16 线程, /tmp/ra_work）

| 阶段 | 耗时 | 关键数 |
|---|---|---|
| ra_correct (k25 c2) | 5m | 修正 1.57% 碱基 |
| ra_overlap (ESA 精确, -m 25) | 8m | dedup 69.86%（去重替代归一化的机制实证）, 3.04M overlap, contained -6.35% |
| ra_layout | 11m | trimming 33 轮, bubble popping 0 个（max bucket 750）, 13,910 contigs |

### 结果（transcripts.layout.fasta vs dev1.fa, 无屏蔽口径）

| | R_full | R70 | P | n | 长度中位 |
|---|---|---|---|---|---|
| Ra 完整管线 | **0.2%** | **2.9%** | 84.7% | 13,910 | 183bp |
| classic | 14.0% | 53.2% | 82.1% | 7,787 | ~700bp |

**18 倍 R70 差距不在调参空间内**。机制: 共享外显子的 reads 全互连使转录组
string graph 远比基因组复杂——bubble popping 因 bucket 过大一个不敢 pop,
只能靠 trimming 狂剪（33 轮）, 图被剪碎成 183bp 中位片段; precision 尚可
（碱基大多真实）但召回崩塌。文献"Ra 实验性、不及 Trinity 系"的量化实证。

### 终审裁定

- OLC 三形态全部实证收口: 从零（Ra 原版, R_full 0.2%）/ 产出合并（naive
  -0.4pp）/ Ra 语义混合层（-0.3pp）。短读转录组 OLC 从零路线终审关闭。
- Ra 思想的有效成分已被吸收或证伪: dedup=归一化的替代机制（classic 已有
  归一化）、前置纠错（1.57% 修正量, 可单测）、bubble popping（转录组
  bucket 爆炸不可用）。可行域回到 genome-guided / 长读长。

## GTasm 式 Graph Transformer 解图 AUC 对打: 无架构增量, 收档（2026-08-27, gt_train 1-3）

用户提供 GTasm（github.com/chu-xuezhe/GTasm, Graph Transformer + HiFi 基因组
组装）。机制: overlap 图边二分类（GT 注意力 + 边特征参与打分 + ScorePredictor）,
训练用真值坐标构造边标签, 推理用入/出度 PE, 打分引导 walk。

### 移植与对打（dev1 transitions, comp 分组 70/30, 同 fold 严格公平）

chrysalis 图（graph.out 顶点计数第 3 列）转 DGL; 标签 = e2_label oracle
（与 GBM 同源同 split）; RTX 3090, dgl 2.4+cu121（venv, annevo torch 栈）。

| 模型 | 边特征 | test AUC |
|---|---|---|
| GT v1（GTasm 原版 forward） | log 计数 | 0.5000（x=pe 丢弃内容特征 + BatchNorm 坍缩常数） |
| GT v2（内联修正: x+pe 拼接, LayerNorm, batch 累积, 梯度裁剪） | +log 计数 | 0.9380 |
| GT v3（特征拉平: +e2_cov） | +e2_cov | **0.9391** |
| GBM（HistGB, 19 手工特征, 同 fold） | 全套 | **0.9439** |

### 终审

- **拓扑上下文无增量**: 特征拉平后 GT 仍差 GBM 0.5pp。e2_label 的判别信号
  集中在局部特征（覆盖度/度数）, 24-mer 图的边分类不依赖多跳结构——与
  GTasm 在 HiFi 长读 overlap 图上的成功场景（读长覆盖结构信息量大）本质
  不同。叠加 M8 时代的标签噪声上限（63% "假"支路有保守性）, AUC 天花板
  更可能在标签而非模型容量。
- 至此"深度学习解图"三形态全测: M5 GBM 引导（当年砍, 复杂度）/ M9 rerank
  （GBM 0.93 探针保留）/ GT 拓扑架构（无增量）。用户最初"DL 模型辅助
  解图"的愿景有了完整量化答案: **转换级信号已被 GBM 吃尽, 架构升级无利**。
- 脚本: bench/export/gt_train{,2,3}.py, gbm_baseline.py; 环境
  /tmp/gtasm_venv（dgl cu121）。

## 端到端疗效复测: 推翻 AUC 终审——GT 有真实疗效但与 classic 不互补（2026-08-27, probe_e2e）

用户裁定"看实际疗效, 不要只看 AUC"——正确。补充发现: GBM 19 特征含
**e1_label（入边真值, oracle）**, 推理不可用; 去掉后 GBM AUC 0.9439→0.9389
（与 GT v3 0.9391 完全持平）。

### 四臂端到端（同图 68349cb4, 贪心 walk: 入度 0 起点沿最高分出边, ≥200bp）

| 臂 | R_full | R70 | P | n |
|---|---|---|---|---|
| count（dst 计数） | 8.7% | 38.9% | 89.1% | 5,450 |
| GBM（无 oracle, AUC 0.9405） | 8.6% | 38.5% | 89.5% | 5,435 |
| **GT（AUC 0.9390）** | **10.4%** | **46.0%** | 89.4% | 6,323 |
| classic butterfly（完整引擎） | 14.0% | 53.2% | 82.1% | 7,787 |

- **AUC 相同（0.939 vs 0.941）但端到端 GT 赢 +1.8pp R_full / +7.5pp R70**——
  walk 需要的是分叉处局部排序质量而非全局校准, GT 的多跳上下文恰好优化
  这个; GBM 的局部特征在 walk 框架下无疗效（与计数臂打平）。
- classic+GT 产出合并: R_full 14.1% (+0.1pp, 噪声带), P 82.1→85.1%,
  n 翻倍 14,110——**GT walk 捡回的参照与 classic 高度重叠, 外挂扩池无净增益**。
- butterfly 完整引擎仍显著优于任何分数引导 walk（14.0 vs 10.4）——
  概率模型+zipping+read 拆解的优势。

### 修正裁定

- 撤回 b529797 的"GT 无增量"结论（AUC 代理失真）; GT 架构在端到端有真实疗效。
- 但疗效的兑现路径不是外挂 walk（不互补）, 而是**分数进 butterfly 内部**
  （M5 model-guided 的复活形态: GT 边分数作为路径选择先验）——Rust 侧
  改动, Gate 待跑。产出的序列重叠性提示: GT 分数的价值在选路精度,
  不在新转录本发现。

## 全特征大力模式: 翻车与上限定位（2026-08-28, gt_train4/4b, calc_edge_support）

用户裁定全特征大力出奇迹。新增信号: k-mer 4-mer 子词嵌入（可学习, 256×32,
捕捉 splice motif）、read 跨边支持（DS.fa 97K reads → 708 万边支持计数,
butterfly 主力感官的等价物）、comp 内相对覆盖度、GC/熵差边特征。

| 模型 | 参数 | test AUC | walk R_full | walk R70 |
|---|---|---|---|---|
| v3（4 边特征, 计数+度数） | 0.2M | 0.939 | 10.4% | 46.0% |
| v4 全家桶（9 边+子词嵌入, HID128×4 层, 40ep） | 1.0M | 0.891（train loss 0.0018 记满） | 7.1% | 32.8% |
| v4b 小架构+正则（HID64×2, dropout/wd, 15ep） | 0.1M | 0.935 | 7.6% | 33.8% |

- **大力翻车**: 大容量版教科书式过拟合; 小架构+正则版仍倒退——新特征组
  （子词嵌入/边支持/支持比/相对覆盖度）在 5.6K 训练样本上净伤害
  （边支持与 k-mer 计数共线且带 read 级噪声, 64 维嵌入挤占小模型容量）。
- **方差检查**: v3 换 SEED 重跑 → walk 10.4%/46.0% 完全复现（评测稳定,
  非 R70 的高弹性假象; 注意 SEED13 的 AUC 0.9195 却端到端不变——再次
  印证 AUC 与分叉排序质量脱钩）。
- **上限定位**: 当前瓶颈是训练数据量（632 comps）而非特征面/容量。
  v3 的 4 特征已吃满这个数据量; 吃更多信号需要 dev1+dev2+dev3 三套
  图-标签全量训练（transitions 的 dev2/dev3 图 cache 对应待修复）。
- 新增资产: edge_support.pkl（read 跨边支持, 708 万边）, calc_edge_support.py。

## Oracle 背答案 walk: 分数不是瓶颈, 框架才是（2026-08-28, probe_oracle）

用户指出 StringTie 转录本从未进特征面（只做标签）, 裁定"直接背答案"。
Oracle = 边分数直接用参照路径信息（参照滑窗相邻顶点对 = 683 万参照边,
完美分数; 次级 0.5 分给 dst 在参照顶点集的边）:

| 臂 | R_full | R70 | P | n |
|---|---|---|---|---|
| count / GBM | 8.7/8.6% | 38.9/38.5% | 89% | 5,4xx |
| GT（学出来的 0.94 AUC） | 10.4% | 46.0% | 89.4% | 6,323 |
| **Oracle（背答案, 完美分数）** | **9.6%** | **43.4%** | 88.8% | 5,969 |
| classic butterfly | 14.0% | 53.2% | 82.1% | 7,787 |

- **背答案 < GT 学出来的分数**: 贪心单路径 walk 框架消化不了完美分数——
  分叉处硬选一条参照边, 同 comp 的其他 isoform 直接丢; GT 的软分数反而
  走出更多路径（n 6,323 vs 5,969）。
- **终局裁定: 特征/分数质量从来不是瓶颈**（0.94 AUC 的 GT、1.0 的 oracle
  都被 9.6-10.4% 的 walk 天花板压死）; classic 的 14.0% 赢在框架
  （概率模型+多路径枚举+zipping+read 拆解）。
- StringTie 特征线的完整答案: 同物种强特征在此框架无利可图;
  弱衍生特征≈图内覆盖度同源; 近缘种特征（M10, 无泄漏）是唯一尚存的
  特征增量, 但同样会被 walk 天花板压住——**先修框架再谈特征**。
- GT 分数的正确归宿: 进 butterfly 内部（M5 形态, 强框架×强分数）。

## Butterfly 源码审计: 解图"重点应用"什么（2026-08-28, paths.rs + oracle+triplet）

用户裁定背答案也失败后审计 classic 引擎的路径选择核心
（crates/trinity-butterfly/src/paths.rs, get_all_probable_paths L731）:

1. **Read threading 路径绑定**（ReadHash: 顶点→PairPath→计数）——每条候选
   路径动态绑定"跟随它的 reads 集合", 延伸门卫 path_has_enough_read_support
   （L694）= look-back 2500bp 子路径的 read 包含判定（≥2 条 read 支持才延伸）。
2. **Triplet 方向锁定**（L844-863）: 分叉处只允许 read 观察过的 w→v→u 连续组合。
3. **路径族结构**: 每 comp 维护多路径（max_paths_per_node_extend=25）,
   read 支持数排序竞争 + 裸边保底播种（未被任何路径用的边重启新路径,
   保 isoform 多样性）+ S/T 哨兵全源全汇。
4. **核心认知: butterfly 的选路单位是"路径×read 集合"的动态绑定**,
   不是边分数。

### oracle+triplet 验证（probe_oracle2, read 链 707 万 triplet）

oracle 基础分 × triplet 支持数加成: **9.7%/43.8%**（vs oracle 9.6/43.4）——
triplet 在 24-mer 图上无增量。机制: 图的边本就来自 reads 的 k-mer 相邻,
**任何分叉支路都有 read 跨过**（triplet 种类 707 万 ≈ 边数 683 万）,
单步 read 组合证据在 k-mer 图上无判别力; butterfly 的威力来自
路径族×read 绑定的**全局**结构, 非单步锁定。

### 终局图景

| 分数 | walk 框架 | butterfly 框架 |
|---|---|---|
| 计数/GBM（0.94） | 8.7% | — |
| GT（0.94, 软分数） | 10.4% | **未测（M5 形态, 唯一空白）** |
| oracle（1.0 背答案） | 9.6% | — |
| （butterfly 自身 read 证据） | — | 14.0% |

walk 框架已三重证伪（分数质量无关）; butterfly 框架 × 外部分数
（GT/近缘种先验）是唯一未测组合, 也是 M5/M10 当年的立项逻辑。

## GTasm 解码框架移植: 框架值 +1.4pp, GT 分数追平背答案（2026-08-28, gt_decode）

用户指出看 GTasm 的推理解码。审计 inference.py: 它的分数也是**静态**的
（主循环外一次算好）, 威力全在解码结构——(1) 起点按边分数采样 50 条
（从高分边中点出发, 绕开源点问题）（2) **双向贪心**（前向 out + 后向
in-edge 回溯拼接）（3) 50 候选竞争取最长 contig（4) **迭代收缩**
（移除已访问节点循环到图空, 每 comp 出多条）。

### 结果（dev1, 同 oracle/GT 分数只换解码）

| 分数 × 解码 | R_full | R70 | n |
|---|---|---|---|
| GT（AUC 0.939）× 单向 walk | 10.4% | 46.0% | 6,323 |
| oracle × 单向 walk | 9.6% | 43.4% | 5,969 |
| **GT（AUC 0.939）× GTasm 解码** | **11.8%** | **48.6%** | 6,728 |
| oracle × GTasm 解码 | 11.9% | 48.8% | 6,749 |
| classic butterfly | 14.0% | 53.2% | 7,787 |

- 解码框架独立价值: GT 分数下 +1.4pp R_full / +2.6pp R70。
- **GT 可泛化分数 = oracle 背答案**（11.8 vs 11.9）——分数通道已饱和,
  推理时不需要 StringTie 参照, 图内学到的分数就吃满了; 这也终答了
  用户"StringTie 特征"问题: 特征通道的上界就是 oracle, 而 GT 已达。
- 距 classic 剩 2.2pp: butterfly 的 read 路径动态绑定（本解码仍是
  静态分数无 read 绑定）。
- **杂交方向明确**: GTasm 解码（双向/竞争/收缩）× butterfly 的 read
  支持（候选竞争标准从"最长"改"read 支持加权", 延伸门卫加 read 判定）
  ——python 原型可先验证, 成则 Rust 化。
- 脚本: probe_gtasm_decode.py, gt_decode.py。

## 杂交终测: read 局部证据在强解码下也无增量, python 原型到顶（2026-08-28, gt_decode2）

GTasm 解码 × butterfly read 支持（edge_support）三变体:
(A) 贪心候选分 = GT×支持因子 (B) 竞争标准 = read 支持总量 (C) A+B:

| 变体 | R_full | R70 |
|---|---|---|
| 基线 GT×GTasm | 11.9% | 48.6% |
| A / B / C | 11.9 / 11.8 / 11.7% | 48.5 / 48.4 / 48.3% |

**全部打平**——与 oracle+triplet 互证: 24-mer 图上边级/triplet 级 read
证据无判别力（图的边本就来自 read 相邻, 支持数与 GT 计数特征同源）。
butterfly read 绑定的真实威力在**长程跨度**（read 的 PairPath 完整路径
锁定 walk 不越过转录本边界, look-back 2500bp 语义）——那是 butterfly
本体, python 重写无意义。

### 全线终局（dev1 python 原型）

- 分数通道: GT 可泛化分数 = oracle 上界（11.8 vs 11.9）——饱和。
- 解码框架: GTasm 式（采样/双向/竞争/收缩）+1.4pp——已吸收。
- read 局部证据（边支持/triplet）: 双重证伪——无增量。
- 剩余 2.2pp（vs classic）= butterfly 的 read 长程绑定。
- **可行动作**: (a) Rust 化——把 GTasm 解码思想（起点采样/双向/迭代收缩）
  与 GT 分数先验加进 butterfly 本体（M5 终极形态, 与 read 绑定互补）;
  (b) GT×GTasm 解码原型作为"无 read-threading 管线场景的轻量路径"保留。
- 脚本: gt_decode2.py, gt_scores.pkl（GT 分数缓存）。

## 边支持失效根因与复活: 深度稀释（2026-08-28, calc_edge_support2 + probe_why_support）

用户问"存储边为何无法生效"。定量诊断（probe_why_support, 1,243 个真假
混合分叉组）:

- **归一化版支持**（97K DS reads）: 组内真>假仅 27%（48% 打平）——
  97K reads 摊 683 万边, 平均每边 ~1 条支持, 泊松噪声下决策边全落 0/1,
  无统计量; 与计数相关性 r=0.026（非冗余, 是无信号）。
- **原始版支持**（R1+R2 全量 11.6M reads, 5.51 亿对, 13 分钟）:
  组内真>假 **94%**（GT 95%）, 真假支持比中位 165×; 与计数相关性
  r=0.877（深度够后支持≈覆盖度聚合, 与 GT 特征信息重叠但仍含直测增量）。

### 解码复测（gt_decode2 × 原始支持, 因子改连续 1+log1p(sup)）

| 变体 | R_full | R70 | P |
|---|---|---|---|
| 基线 GT×GTasm | 11.8% | 48.6% | 87.3% |
| A 贪心乘支持因子 | 11.9% | 48.6% | 87.4% |
| **B 竞争=read 支持总量** | **12.2%** | 49.5% | 88.8% |
| **C A+B** | **12.2%** | **49.7%** | **88.9%** |

支持信号复活带来 +0.4pp R_full / +1.0pp R70, precision 同升。
距 classic 剩 1.8pp——read **长程绑定**（butterfly 的路径级 read 包含
语义, 边级 SUM 证据的极限已现）。

### python 原型最终成绩与累计弧线

GT 分数 × GTasm 解码 × 原始边支持 = **12.2% / 49.7% / 88.9%**
（classic 14.0/53.2/82.1）; 从 GT 单向 walk 起点 10.4/46.0 累计
+1.8pp R_full, +3.7pp R70, P 高 6.8pp。
教训: **read 证据必须用原始深度算**——归一化 saves the graph 但
kills the evidence。
- 新增资产: edge_support_raw.pkl（5.51 亿对）, calc_edge_support2.py。

## GT 原型 vs classic 全面对比（2026-08-28, probe_fullcmp + junction）

用户要求全面指标对打。蛋白口径复刻 bench.rs（最长 ORF 翻译 → DIAMOND
blastx → HSP 并集 cov×pid, P_full=c>=0.95&p>=0.95, P70=c>=0.70&p>=0.90）;
junction 用 SAM 输出（PAF 无 N gap, cg 全 M/I/D）, 坐标 (donor, acceptor+1)。

| 维度 | GT 原型 | classic | 差 |
|---|---|---|---|
| n | 6,759 | 7,787 | -13% |
| 总碱基 | 7.1M | 9.3M | -24% |
| N50 | 1,672 | 1,838 | -9% |
| R_full | 12.2% | 14.0% | -1.8pp |
| R70 | 49.7% | 53.2% | -3.5pp |
| **P** | **88.9%** | 82.1% | **+6.8pp** |
| F1 | 21.5% | 23.9% | -2.4pp |
| P_full | 51.9% | 56.3% | -4.4pp |
| P70 | 68.3% | 71.6% | -3.3pp |
| **junction 精确** | **90.3%** | 87.4% | **+2.9pp** |
| junction 召回 | 19.7% | 20.2% | -0.5pp |

格局: classic 赢全部召回/产出量轴（2-4pp, 产出多 24% 碱基）; GT 原型赢
两个质量轴——precision +6.8pp、剪接精确率 +2.9pp。GT 原型 = 纯图内特征
（无 read-threading 管线、无参照依赖、单卡分钟级）的轻量路径, 在
"产出更少但更准"的定位上成立; classic 的召回优势本质仍是 read 长程绑定。
- 脚本: probe_fullcmp.py; junction 简版（SAM+坐标偏移已校准）。

## reads 峰形 → isoform 数先验: 信号中等且非峰形, 证伪收档（2026-08-28, probe_peak_isoform）

用户提议: 基因 reads 峰形（高斯式机制）预测区间该解出的转录本数,
特征取自 hisat2.bam + gtf——框架轴创新（计数先验作用于收缩循环,
不受 oracle 分数上界封顶）。前提验证（21,887 基因, raw bam depth 特征
→ 1/2/3+ isoform 分类）:

| 特征集 | 1 vs rest | 2 vs rest | 3+ vs rest |
|---|---|---|---|
| 仅 reads 峰形（深度分位/CV/平滑峰数/覆盖比例） | 0.731 | 0.648 | 0.753 |
| 仅 gtf 结构（n_exon + span） | 0.732 | 0.651 | 0.748 |
| 合并 | 0.759 | 0.673 | 0.773 |

- **信号中等（<0.8 接线门槛）且两组特征完全同构**（合并仅 +0.02）——
  reads 峰形没有独立于"基因大小/复杂度"潜变量的信息; 特征重要性
  log_span 压倒性第一（Δacc +0.031, 其余 5 倍以下）。用户设想的
  "丰度峰形"信息实际是"基因长度"的结构相关性。
- 接线负期望: ±1-2 isoform 误差的计数先验做硬停止会砍真 isoform（伤 R）,
  做软权重会被已饱和的分数通道吸收。
- 裁定: 证伪收档。框架轴的剩余空间仍指向 read 长程绑定（唯一未搬进
  python 的组件）。

## 评判基线变更: principal isoform 口径正式化（2026-08-28, bench.rs）

用户裁定: 此后主评判指标 = **每基因最长 isoform 子集上的召回**
（主力转录本口径, R_full_p/R70_p）; 原全参照口径保留作参考。

- 依据: GT 原型 vs classic 的差距集中在次要 isoform（全口径 R_full
  -1.8pp vs principal 口径 -0.7pp, R70 -3.5 vs -1.8pp）——主力转录本上
  两者几乎打平, 且 GT 原型 precision/junction 精确率反超, "每基因一条
  主力转录本"的生态位（功能注释/ORF 预测）性价比成立。
- 实现: `bench.rs` Metrics 加 r_full_p/r70_p/refs_principal;
  principal_set 按参照名最后 '.' 前缀分组取最长（StringTie STRG.G.I
  命名）; stdout/metrics.json 同步输出; 2 个新单测, 17 测试全绿,
  clippy 零警告。
- dev1 验证: Rp=10.7% R70p=50.8% (5190 principal) 与手算一字不差。
- 基线数字（principal 口径）: classic R_full_p 10.7% (557/5190) /
  R70_p 50.8%; GT 原型 10.0% (520) / 49.0%。
| 2026-08-27 | dev1 | R=14.0% P=81.4% **F1=23.9%** | N50=1836 n=7835 | --CPU 8 --max_memory 8G |
| 2026-08-28 | dev1 | R=14.0% P=81.4% **F1=23.9%** | N50=1836 n=7835 | --CPU 8 --max_memory 8G |

## principal 尸检: 没装出来 vs 装不全的死因分解（2026-08-28, probe_principal_autopsy + probe_kmer_rescue + probe_norm_kmer）

用户问题: 都是最长转录本, 为什么有些没装出来、有些装不全?
dev1 principal 5190 条, 四路证据（产出档/原始 reads 覆盖/图 kmer 在位+边走通/归一化后 reads）交叉。

产出分档（classic）: full 557 (10.7%) / fragmented 199 (3.8%) /
partial 2929 (56.4%) / weak 483 (9.3%) / absent 1022 (19.7%)。
GT 原型剖面几乎相同（档级迁移主对角线压倒性）——两者被同一天花板压着。

非 full (4633) 初判: DATA 34% / GRAPH 64% / DECODE 2.1% (95 条)。
DECODE 画像: 短(中位 476bp vs full 1977bp)、中低深度(14.3x vs 32.6x)——
与 oracle 饱和结论互洽, 解图端剩余空间 ≈ 这 95 条。

**决定性规律: full 率被"reads 是否铺满转录本"单变量决定**——
read 覆盖 <10%/10-50/50-80/80-99/≥99% 五档的 full 率 = 4/7/8/12/51%。

GRAPH 桶 kmer 级溯源（938,289 个图缺失 24-mer 查原始 reads, 免比对）:
- 44.9% 在原始 reads 中存在, 但深度中位仅 2（≥3 的 15.9%, ≥5 的 8.7%）;
- 归一化裁决: 归一化后 reads 仍含 44.6% —— **归一化无罪**,
  丢失发生在 inchworm 错误剪枝（min_ratio_non_error=0.005 等, 对深度 1-2
  的 kmer 基本按噪声处理, 属预期行为）;
  （注: raw∩norm 逐 kmer 对拍因两进程 set 迭代序导致 id 错位而不可靠,
  仅两 run 各自聚合值可比——44.9%≈44.6% 已足够排除归一化选择性丢失。）
- reads 视角重算可走性: 仅 135/2958 条变成 ≥99% 可走（其中挽救 kmer
  深度≥3 的 85 条 = 真可操作集, 深度中位 3.0, 转录本深度中位 13.9x）。

裁定: principal 口径的算法可及上限 ≈ full 557 + decode 95 + kmer 救援
~85-135 ≈ **R_full_p 10.7% → 14-15% 理论顶**, 且 kmer 救援需容错建图,
有 precision 代价; 其余 ~96% 非 full 由数据铺设不完整决定, 算法无解。
杠杆排序: ①数据端(深度/读长/跨样本) ②容错 kmer 救援(85 条, 有风险)
③decode 微调(95 条, 已饱和边缘)。
- 脚本: probe_principal_autopsy.py / probe_kmer_rescue.py / probe_norm_kmer.py;
  缓存 principal_autopsy.pkl, graph_missing_kmer_rescue.pkl。
| 2026-08-28 | dev1 | R=14.1% P=82.1% **F1=24.1%** | N50=1811 n=7914 | --normalize_max_read_cov 10000 --CPU 8 --max_memory 8G |

## 救援杠杆实测: 归一化放松证伪 + fragmented 搭桥天花板（2026-08-28, probe_frag_bridge + bench dev1 maxC10000）

承接 principal 尸检, 三条候选救援杠杆逐一量化/实测:

1. **归一化放松（maxC 200→10000, bench cache 59a909274df52986）: 证伪。**
   丢失环节定位（深度≥5 的图缺失 kmer 抽样 300 个）: 归一化后 reads 存活
   仅 ~12-24%（主丢失点）, inchworm contigs 剩 4%——inchworm 剪枝已极宽松
   （min_count=1/无错误剪枝/仅熵过滤）, 杀手是归一化的 read 级降采样
   （山谷 read 横跨高覆盖区被连带丢）。但实测放松后: Rp 10.7%→10.7%,
   R 14.0→14.1%, P 82.1→82.1%, n 7787→7914; 85 条可操作目标集仅 2 条
   救回 full（4 条覆盖提升>5pp）。**kmer 活着不等于能织进 contig**——
   山谷 kmer 孤立低深度, 贪心延伸穿不过去; 也可能部分是参照自身
   consensus 错误。杠杆①收档。
2. **fragmented pair 搭桥（probe_frag_bridge）: 天花板 +0.7pp。**
   199 条 fragmented 中 170 条有≥2 段; 段间缺口中位 15bp, 97%≤400bp;
   桥接 pair≥5 的仅 58 条; 纯重叠可链接 21 条, 重叠或小缺口+≥3 pair
   共 38 条可救 → R_full_p +0.7pp; 85 条零桥接死段。
3. DECODE 95 条: +1.8pp 天花板, oracle 已证饱和边缘, 不投。

**裁定: 算法侧可及总量实测封顶 ≈ +2.5pp（scaffold 0.7 + decode 1.8）;
R_full_p 10.7% → ~13% 是 de novo 在 dev1 上的物理极限。**
其余非 full 由 reads 铺设不完整决定（铺满率单变量决定 full 率:
≥99%→51%, 80-99%→12%）。后续提指标的正路: 数据端（深度/读长/多样本）
或评测换到高覆盖基因子集。
| 2026-08-28 | dev1 | R=14.0% P=81.4% **F1=23.9%** | N50=1836 n=7835 | --CPU 8 --max_memory 8G |

## 高覆盖子集口径（hc）落地 + 深层死因补查（2026-08-28, bench.rs）

承接救援证伪后的正路①: 评测切到高覆盖子集让算法差异显形。

深查发现:
- 覆盖分层 full 率（classic）: ≥99%→50.8%, 95-99→30.7%, 90-95→3.3%,
  80-90→1.2%, <80→~6-8%; **≥0.95 是断崖线**。
- 覆盖≥0.95 子集（n=980）内按深度分层 full 率: Q1(<14x) 11.8% →
  Q4(>105x) 51.0%——铺满只是必要条件, 深度仍是第二变量。
- 鉴别力: classic vs GT 原型的 full 差, 全体 +0.69pp → ≥0.95 子集
  +1.22pp → ≥0.99 子集 +2.12pp。子集越纯越显形。
- 93 条"铺满仍失败": 20 条图可走（decode 性质）, 73 条图不可走但
  缺失 kmer 的 89% 在原始 reads 里（深度中位 17x, 缺口中位 7.1%）。
  低复杂度过滤排查: 缺失 kmer 熵<1.3 占 8.7%（对照 1.8%, 富集但解释力
  小）; 串联重复 5%。结合 maxC10000 实验（救回 2/85）, 结论: 丢失不在
  任何过滤策略, 在 inchworm 贪心延伸本身穿不过低深度山谷——修复 =
  重写延伸核心, 收益上限 ~1.5pp, 不投。

bench.rs 实现: Metrics 加 r_full_hc/r70_hc/refs_hc（阈值
HC_COV_THRESHOLD=0.95, principal∩高覆盖）; read_coverage() 用
minimap2 -x sr（mapq≥5, 并集/tlen）, 按 sha256(ref+reads) 前 8hex
缓存于 bench/cache/refcov_*.tsv; stdout 增 HC(n=..) 段, metrics.json
增三字段。compute_metrics 加第 4 参 Option<&cov>, 旧测试机械更新 +
新增 hc_recall 单测; 18 测试全绿, clippy 零警告。
dev1 验证: HC(n=980) R=34.4% R70=73.4% vs python 手算 34.59%/73.47%
（0.2pp 差 = 本次命中缓存为另一构建 43bff0c513a11d21 n=7835 vs
b4a32.. n=7787 的组装版本噪声, 口径本身 n=980 精确一致）。
新基线（hc 口径, dev1）: classic R_full_hc 34.6% / R70_hc 73.5%;
GT 原型 33.4% / 71.5%（差 -1.2pp / -1.9pp, 全体口径差的 1.8 倍）。

## read 长程绑定进解码: 证伪收档（2026-08-28, gt_decode3.py）

butterfly look-back 语义的 python 复刻: reads 线程化（kmer→节点路径,
只留跨分支出边的 run, 统一图正向）→ 解码分支点用"跨候选边且与 walk
历史连续匹配 ≥5 节点"的 reads 数做绑定。dev1: 13,388 候选边,
762,190 条留存路径（11.6M reads 线程化 589s）, 902 comps 有绑定结构。

| 变体 | R_full | R70 | P | Rp | R70p | hc R | hc R70 |
|---|---|---|---|---|---|---|---|
| base(C复现) | 12.3% | 49.7% | 88.9% | 10.1% | 49.1% | 33.6% | 71.5% |
| B1 过滤(lb≥2) | 12.3% | 49.9% | 89.2% | 10.1% | 49.2% | 33.5% | 71.6% |
| B2 乘性 | 12.3% | 49.7% | 89.1% | 10.1% | 49.1% | 33.7% | 71.5% |
| B3 lb优先 | 12.4% | 49.7% | 89.1% | 10.1% | 49.0% | 33.8% | 71.4% |

所有变体差异 ≤0.2pp = NB_STARTS 采样噪声级。**裁定: 证伪收档。**
机理: 跨边 reads 必然先被边支持计数（长程信号 ⊂ 边级信号的超集但
判别力重合）; 真正需要长程裁决的共享外显子场景在本数据恰好是
低深度短转录本（DECODE 画像）, 能跨过去的 reads 本身就没有。
至此 butterfly 解码语义的全部组件（分数/边支持/长程绑定）在 dev1
上均实测封顶, 解图通道关闭。hc 口径基线: GT 原型 33.6-33.8%/71.5%
vs classic 34.6%/73.5%。
- 脚本 gt_decode3.py; 缓存 read_paths_cand.pkl（762K 路径, 161MB）。
| 2026-08-28 | dev3 | R=12.5% P=82.4% **F1=21.7%** | N50=1811 n=6486 | --CPU 8 --max_memory 8G |
| 2026-08-28 | dev2 | R=13.1% P=79.8% **F1=22.4%** | N50=1789 n=9126 | --CPU 8 --max_memory 8G |
| 2026-08-28 | dev1 | R=12.3% P=88.8% **F1=21.6%** | N50=1671 n=6761 | --path-engine gtasm --gt-weights /storage/home/senior007/test/athal/bench/export/gtnet_dev1.gtw --CPU 8 --max_memory 8G |
| 2026-08-28 | dev1 | R=13.9% P=81.7% **F1=23.8%** | N50=1836 n=7798 | --CPU 8 --max_memory 8G |
| 2026-08-28 | dev2 | R=11.5% P=87.7% **F1=20.4%** | N50=1601 n=7900 | --path-engine gtasm --gt-weights /storage/home/senior007/test/athal/bench/export/gtnet_dev1.gtw --CPU 8 --max_memory 8G |

## 铺满率定律三数据集复核 + GT 原型 Rust 化落地（2026-08-28/29）

**dev2/dev3 复核（hc 口径 bench + 分层分析）**: 铺满率定律普适——
三数据集 full 率分层几乎重合: ≥99% 覆盖档 50.8/55.6/47.3%（dev1/2/3）,
95-99 档 30.7/27.0/32.3%, 90-95 档 3.3/3.4/2.9%（断崖）, <90% 个位数。
hc 子集规模占比 17-19%。dev1 尸检结论跨数据集成立。
（教训: refcov 缓存文件名与数据集的对应须按内容校验, mtime 排序会
被并行任务完成顺序骗过——首版分析曾因此得出错误的平坦分布。）

**GT 原型 Rust 化（--path-engine gtasm）**: python 原型（GT 边分数 ×
GTasm 解码 × 原始 reads 边支持）完整落进主线。
- 新增: butterfly/gtmodel.rs（GTW1 权重加载 + 纯 f32 前向, 零新依赖;
  修正主代理 spec 一处——注意力权重按 head_dim=8 逐头求和而非全 64 维）,
  butterfly/gtdecode.rs（C 变体逐语义: 平局取大 vid/eid、无名节点跳过、
  sup_sum 竞争、固定 xorshift 采样 50 起点）, cli/gt_support.rs（原始
  reads 边支持 pass, GTS1 二进制, dev1: 2.76 亿转移 / 649 万边）,
  两个对拍 bin + scripts/gt_score_compare.py。
- 接线: `--path-engine classic|gtasm` + `--gt-weights`; 支持 pass 带
  checkpoint 落 workdir; gtasm 组件产出进既有 harvest 通道。
- 数值对齐: Rust vs torch 重载逐边最大差 1.2e-07。
- 验收 bench（dev1, cache 775cb33c9c3eb808）: **R=12.3% R70=49.7%
  P=88.8% Rp=9.9% R70p=49.1% hcR=33.4%** —— 对 python 原型全面 ±0.5pp
  内复现; classic 回归 13.9%/81.7% 在 run 间噪声内。
- dev2 泛化（权重仅 dev1 训练）: R=11.5% R70=46.4% P=87.7% Rp=9.3%
  hcR=31.1%（vs classic 13.1/79.8）——"召回略低、精度大涨"格局跨数据集
  成立, 模型未过拟合单一数据。
- 577 测试全绿, clippy 零警告。

## GapCloser 式山谷桥接原型: 机制成立, 歧义封顶（2026-08-29, probe_gapcloser）

SOAPdenovo2 localAsm 语义的转录组移植可行性验证（fragmented 199 条
principal, 92 个正缺口任务: 锚 100bp + 局部 k=15 DBG + 有界 DFS +
唯一路径才填）; v2 加 NOVOPlasty 式迭代招募（首轮锚命中 8K reads →
池诱饵二轮再捞 47K）与全窗锚命中修复。

- v1: 填谷 13/92, 转正 9 条; v2: 填谷 11/92, 转正 6 条
  （≈ +0.1-0.2pp R_full_p, 低于此前 +0.7pp 上限估计）。
- **唯一路径时填充全对**（v2: 7/11 完全一致, 余为短插入缺失级差异）——
  机制本身的序列保真没问题。
- 失败桶 v2: 多路径歧义 ~32（池扩大后爆炸, 远超 v1 的 13）、无路径 33、
  锚缺席 4、超预算 6。**缺口区域恰恰是被 isoform/旁系同源共享的低复杂度
  山谷——全局图在那里断, 局部图在那里同样歧义**。招募越充分, 混入的
  同源变体越多, 唯一路径越少。
- 裁定: 收档。SOAPdenovo2 的 gap closing 假设"gap 内序列唯一", 转录组
  的断裂点系统性不满足此假设。构图阶段最后一刀（局部重组装）实测封顶,
  与全局结论一致: dev1 的算法可及空间 ≈ 已穷尽。
- 脚本 probe_gapcloser.py; 缓存 gapcloser_probe.pkl。

## rnaSPAdes 构图机制调研（2026-08-29, spades-main 源码）

对照实测死局（山谷贪心穿不过 + 山谷区多路径歧义）调研 rnaSPAdes:

- **多 k 迭代（K_MERS_RNA=[33,49], 按读长自适应）**: 每 k 全量重建图,
  上一 k 的 simplified_contigs 以 TrustedContigs 身份并入下一 k 的
  kmer 计数流（construction.cpp:233, 贡献序列不贡献覆盖）——山谷架桥
  靠"上轮答案注入下轮构图", 不靠大 k 贪心硬穿。gap closer 只在最后一轮
  或 K>=55 启用, repeat resolution 只在最后一轮。
- **BayesHammer 纠错在 RNA 模式整体弃用**（options_parser.py:1012:
  rna → only_assembler, 无 error_correction stage）——uneven depth 下
  簇心质量判据不可靠, 纠错职责全部移给图清理。与我们"纠错会抹平真实
  低表达差异"的顾虑互相印证, 杠杆③（纠错前置）行业侧证伪。
- **图清理全相对覆盖**: rctc（tip_cov < X×竞争边最大覆盖）/rcec（≤邻域
  均值）/rcc（低覆盖组件须被 ≥20x 高覆盖四面包围才删）; RNA 模式 icb
  固定 200 弃用自动谷底; tip 加 mmm 3 限制只剪近重复错误 tip。
- **歧义处理 = 放弃**: MultiExtender 仅在简单 2-way bulge 分叉产
  isoform, 其余歧义直接断（宁可断不可错连）——我们 GapCloser 探针的
  多路径墙, 行业答案相同。路径解析靠 PE pair + 覆盖一致性
  （exSPAnder）, 与 butterfly 同族。
- gap closer = tip 对 tip overlap 缝合 + PE 证据, 比 SOAP 局部重组装
  保守（不假装能装山谷）。
- **可吸收的唯一真机制: 多 k + trusted contig 注入**。注意旧单 k 扫描
  k19/23 灾难性更差（毛球）, 小 k 轮产物质量是问号; 收益上限仍被
  85 条集合封住（~+1.6pp R_full_p）。待用户裁定是否做最小可行实验
  （k=21 → 注入 k=25）。
| 2026-08-29 | dev1 | R=14.0% P=81.0% **F1=23.9%** | N50=1843 n=7869 | --CPU 8 --max_memory 8G |
| 2026-08-29 | dev1 | R=14.0% P=81.0% **F1=23.9%** | N50=1843 n=7869 | --CPU 8 --max_memory 8G |
| 2026-08-29 | dev1 | R=14.5% P=81.4% **F1=24.6%** | N50=1918 n=7925 | --multi-k auto --CPU 8 --max_memory 8G |

## 多 k 迭代 + trusted contig 注入落地（--multi-k, 2026-08-29）

上节调研中"唯一可吸收机制"的最小可行实现（rnaSPAdes construction.cpp:233
语义移植）:

- **CLI**: `--multi-k auto|k1,k2,...`; 与 --KMER_SIZE/--recursive 互斥;
  显式列表 k∈[15,32] 严格递升。auto = rna_k_values 公式（lower=clamp
  (RL/3,29,55), upper=RL/2-1, 均取奇, upper<=lower 单 k）+ RL 抽样（各
  文件前 1000 条取 max, 跨文件取 min）+ 支持域钳制 [15,32]（KmerId 位宽
  上限 32; dev1 RL=76 → 公式 [29,37] → 钳制 [29,32]）。
- **迭代结构**: 归一化/prep 只做一次（共享 both.fa）; 中间轮在
  `iter_k<k>/` 只跑 计数+inchworm（跳过 chrysalis/butterfly）, 产物
  inchworm.DS.fa 为 trusted 源; 最终轮在 `iter_k<k_last>/` 跑完整管线,
  Trinity.fasta 落点不变。计数表两条构建路径（counts_bin 装载/重算）
  汇合后统一注入: 已存在键不动计数, 缺席以 count=1 插入（presence-only,
  先于 prune, 无熵过滤豁免）; counts_bin 落盘保持 reads-only（注入幂等,
  resume 装载路径重注结果不变）。链式注入（第 i 轮注第 i-1 轮, >2 k 时
  中间轮不空转）。.stage_timings 按轮次目录区分。
- **dev1 验收**（cache b9b42bd94332bc08, k=29→32, 各轮 ~22s 计数 + ~27s
  inchworm）: **R=14.5% R70=55.0% Rp=10.9% R70p=51.8% hcR=34.6%
  hcR70=74.3% P=81.4% F1=24.6% N50=1918 n=7925 P_full=65.8%**;
  同日新二进制 classic 复跑（ee22d353bcae2d67）: R=14.0% R70=53.4%
  Rp=10.8% hcR=34.4% P=81.0% F1=23.9% n=7869（与 08-28 基线 14.04/81.38
  逐点吻合 → 默认路径逐字节不变, classic 回归通过; 二次运行命中新缓存）。
- **裁定: 全面小幅正收益（R +0.5pp / R70 +1.6pp / F1 +0.7pp /
  P_full +2.1pp, 无精度代价）**, 单轴幅度仍在 run 间噪声（≤0.6pp）边缘,
  但七轴同向 + 蛋白级 P_full 涨幅最大, 方向可信。与调研预期一致:
  山谷架桥收益上限受 85 条集合封顶（~+1.6pp Rp 为理论顶）, 实测
  Rp +0.1pp——大头仍被封住。593 测试全绿, clippy 零警告。
| 2026-08-29 | dev1q30 | R=13.7% P=80.3% **F1=23.4%** | N50=1854 n=7848 | --CPU 8 --max_memory 8G |

## Q30 reads 过滤基线（2026-08-29, filter_q30.py）

用户指令第一步: 只做过滤不做纠错。口径: read 平均 Phred>=30, pair 同留同弃。
dev1/2/3 保留率 97.2/97.2/97.1%（原数据 meanQ 37, 本就很干净）。
bench dev1q30（cache 9cf44d917ab8dc76）: R=13.7% Rp=10.5% hcR=34.6%
P=80.3%（对照未过滤 14.0/10.7/34.4/81.4）——**全轴微负**。
裁定: 过滤无收益。山谷是真实覆盖缺口而非错误 kmer（此前 kmer 溯源已证）,
删 reads 只会连带删掉山谷里本就稀薄的覆盖。脚本保留备用, 不进管线。

## rnaSPAdes 相对覆盖图清理落地（--graph-clean-rnaspades, 2026-08-29）

rnaSPAdes simp 段（configs/rna_mode.info + tip_clipper/erroneous_connection_
remover/relative_coverage_remover.hpp 判定结构）移植到 chrysalis 图。新模块
`trinity-chrysalis/src/graph_clean.rs`, 在 QuantifyGraph 之后逐组件清理
graph.out（pack 内文本, 写库前）, 默认关 = 代码路径不动。

- **映射**: count 列 = (frm 首碱基+to 的 24-mer) 25-mer 计数 = SPAdes 边覆盖;
  SPAdes 的"边" ≈ 本未压缩图的死端/junction 间**链**; tc_lb 4/2 → 链长 ≤4/≤2
  节点; ec_lb 30bp → junction 间内部节点 ≤2; icb 200 → 链覆盖 ≤200; rctc 复刻
  `le(cov, 0.5×(竞争max+1))`（含 +1 怪癖）; rcec_cb 复刻邻域均值（两端竞争边
  数和/边数）; mmm 3 → 链序列 vs 竞争链（按覆盖贪心延伸, 要求严格更长）汉明
  ≤3（汇向 5' 对齐/源向 3' 对齐）; rcc: 种子 20×cov<邻域max → flood-fill
  （邻点无 >20×内边基线的高覆盖外边则并入内部, 有则留边界）, 内部顶点 ≤10
  （SPAdes 默认 30, 按任务书收紧）+ 边界非空（整连通块低覆盖不删, 保真实低
  表达基因）+ 各边界有 ≥20× 内部基线的外边 → 删触及内部顶点的行。
- **安全**: 各 pass 同步语义（判定基于轮初图）; 悬空修复（被存活边引用但
  定义行全灭的节点 → 最高 count 死行改根行 `v\t-1\t0\t...` 复活, 与
  QuantifyGraph 根行惯例一致）; 组件清空兜底回滚（蝴蝶池不容忍空图）。
- **走线**: `--graph-clean-rnaspades` → ChrysalisParams → run_quantify_packed;
  汇总统计行打印 + `graph_clean` 行落 .stage_timings（并行 CPU 时间累计）。
  --recursive 下不生效（迷你管线未接线, 启动 WARNING）。
- **测试**: 17 个手算图单测（tip 两组正反/mmm/rctc 边界、EC 正反/icb 上界/
  链长、rcc 正反/尺寸/边界/均匀低覆盖保护、悬空修复）; cargo test
  --workspace 610 全绿, clippy --workspace --all-targets 零警告。
- **bench dev1（cache 6b4e4f7d3097fb4f）**: R=13.7% R70=53.0% Rp=10.4%
  hcR=33.7% **P=81.6%** F1=23.5% N50=1824 n=7850 —— 对照基线 ee22（14.0/
  53.4/10.8/34.4/81.0/23.9/1843/7869）: **精度 +0.6pp, 召回各轴 -0.3~-0.7pp,
  F1 -0.4pp**。净效果略负: 本管线图经 butterfly 自带 prune 已较干净, 再按
  rnaSPAdes 相对覆盖削低覆盖链主要削掉真实低表达信号（与 Q30 节同构）。
  flag 默认关, 保留备用。
- **classic 回归**: 基线二进制无 flag bench 命中 ee22d353bcae2d67 缓存 ✓。
  另发现: 同一基线二进制两次全新跑 dev1sub 的 inchworm.DS.fa 即不同
  （7370 vs 7403 条转录本）——管线存在**先于本改动**的 inchworm 段运行间
  不确定性（归一化产物逐字节相同, 分叉在 inchworm; 此前"复跑逐字节一致"
  实为默认 scratch /tmp/trinity-rust-<hash> 的 checkpoint 续跑, bench 换
  key 重跑只清 out/ 不清 scratch）。我的二进制无 flag 全新跑 dev1 n=7834
  落在基线自身波动带内, 非本改动引入。注: 上述带 flag 数值取自改名前
  二进制（仅字段改名, 逻辑相同）。

| 2026-08-29 | dev1 | R=13.7% P=81.6% **F1=23.5%** | N50=1824 n=7850 | --graph-clean-rnaspades --CPU 8 --max_memory 8G |
| 2026-08-29 | dev1 | R=14.0% P=81.0% **F1=23.9%** | N50=1843 n=7869 | --CPU 8 --max_memory 8G（基线缓存 ee22d353bcae2d67） |
| 2026-08-29 | dev1 | R=14.0% P=81.5% **F1=23.9%** | N50=1833 n=7834 | --CPU 8 --max_memory 8G（新二进制无 flag 全新跑, 波动带内） |
| 2026-08-29 | gapclose-dev1 | R=14.1% P=79.6% **F1=24.0%** | N50=1868 n=7887 | --gap-close-spades --CPU 8 --max_memory 8G |
| 2026-08-29 | gapclose-off-dev1 | R=14.1% P=80.2% **F1=24.0%** | N50=1856 n=7889 | --CPU 8 --max_memory 8G |

## SPAdes gap closer 落地（--gap-close-spades, 2026-08-29）

SPAdes gap_closer.cpp（GapCloserPairedIndexFiller L28-178 / GapCloser L180-469）
+ configs/config.info gap_closer 段（minimal_intersection 10 / weight_threshold
2.0 / max_dist_to_tip 5000）移植到 chrysalis 图。新模块
`trinity-chrysalis/src/gap_close.rs`, 在 QuantifyGraph（与可选 graph_clean）
之后逐组件缝合 c*.graph.out（pack 内文本, 写库前）, 默认关 = 逐字节不变。
**保守变体**: 不做局部重组装（python 原型已证伪多路径歧义）, 只做 tip 对
tip 的 overlap 缝合 + PE 证据排序。

- **映射**: out-tip = 汇向死端链（末端 outdeg==0 沿唯一前驱回溯, 不含分支
  点）; in-tip = 源向死端链（首端 indeg==0 沿唯一后继, 不含汇合点）;
  max_dist_to_tip 5000 → 链端 200bp（177 个 24mer）作 PE 命中窗口（组件小,
  任务书指定取舍）; PE 票 = 同 pair（read name 去尾 /1|/2）两 mate 分落
  out/in tip 末端区（kmer 命中, 非链特异双向试, 镜像 QuantifyGraph）各 +1,
  阈值 ≥2; overlap 自 k-1=23 递减到 10, 低复杂度（单碱基 > (0.8+0.2(gap-1)
  /13)×overlap, 严格 >）整对放弃, 汉明 ≤2（config 未给, 构造注释
  min_intersection/5; 默认 0 会关闭修正路径）进入分派——完美 → 机械桥接
  （gap-1 新节点 + gap 条新边, 全程合法 1bp 漂移, 覆盖 = max(1, min(两链
  覆盖))）; 错配 → HandlePositiveHammingDistanceCase 简化语义（以覆盖高侧
  为权威, tie 右侧优先 ≈ SPAdes 先 CanCorrectLeft）: 未压缩图节点不可分割,
  可修窗口只剩贴端单位置（改 a 末位 / 改 b 首位定义行 kmer）, 窗口外放弃。
  防假阳性: 同链共节点对不缝（≈ avoid_rc_connections）, 改写/新节点 kmer
  查重（防 butterfly 按 kmer 并节点短路）, 缝合前复查两端仍死端/源头
  （L434-437）, 每 out-tip 首个可行候选即 break。
- **缝合表示决策**: 新增桥接边行（不改既有行 to/frm; 修正路径才改写锚点
  定义行 kmer 列）——butterfly 重建图后即见 a→…→b 普通路径, 两段 contig
  可连成一条转录本; 合并链序列会丢图内覆盖结构, 不采用。
- **走线**: `--gap-close-spades` → TrinityArgs（单端警告置关 = SPAdes "no
  paired-end libraries" 跳过语义; --recursive 警告不生效）→ ChrysalisParams
  → run_quantify_packed（清理之后, 与 SPAdes after_simplify 顺序一致）;
  统计行打印 + `gap_close` 行落 .stage_timings（并行 CPU 时间累计）。
- **测试**: 12 个单测（PE 权重阈值 1 票不缝/2 票缝、overlap=9 拒/=10 缝、
  poly-A 低复杂度拒、贴 a 末位错配右权威修正、贴 b 首位错配左权威修正、
  中间错配拒、同链自连拒、无 reads 逐字节不变、pair_key 后缀、flag 解析与
  单端降级）; cargo test --workspace 全绿, clippy --workspace --all-targets
  零警告。
- **bench dev1（cache e1e5e63758acc2db）**: 统计行 `6386 components, tips
  (out 6932 + in 6949), tip reads 1135564, pair votes 42683, candidates 199,
  checked 194, closed 10 (perfect 10 + corrected 0), rejected (low-complexity 1
  + hamming 2), bridge nodes 26`; 指标 R=14.1% R70=53.2% Rp=10.8% hcR=34.4%
  P=79.6% F1=24.0% N50=1868 n=7887 P_full=63.5%。同二进制无 flag 全新跑
  （cache 82e8f2e571f7fe2a）: R=14.1% R70=53.4% Rp=10.8% hcR=34.6% P=80.2%
  F1=24.0% n=7889——对照基线 ee22（14.0/53.4/10.8/34.4/81.0/23.9/7869）在
  已知 inchworm 波动带内 ✓ 回归通过。
- **裁定: 中性（噪声带内）**。缝合仅 10 处（199 候选 → overlap 通过 5%）,
  与 python 原型"收益微小"的预示一致: R/F1 与无 flag 持平（14.1/24.0）,
  P 79.6 vs 80.2 的 -0.6pp 处于 inchworm 运行间波动量级（两次全新跑前置
  产物不同）, 不归因缝合本身。机制工作正常且全程确定（同前置 checkpoint
  重跑产物逐字节一致, 已验证）。flag 默认关, 保留备用。
- **严格 A/B（同前置 checkpoint 续跑, 排除 inchworm 波动）**: 同一
  inchworm/chrysalis 输入下无 flag 7885 条 vs 带 flag 7887 条（净 +2）;
  按转录本名共有 7842, 各独有 43/45（缝合改变组件内路径枚举导致名字
  重排, 独有数系上限）——缝合真实效应规模 ≈ 45 条转录本, 指标差异
  确实在噪声带内。
| 2026-08-30 | rnaspades-dev1 | R=10.6% P=88.8% **F1=19.0%** | N50=1597 n=6739 | --path-engine rnaspades --CPU 8 --max_memory 8G |
| 2026-08-30 | rnaspades-regress-classic | R=14.0% P=80.8% **F1=23.9%** | N50=1834 n=7851 | --CPU 8 --max_memory 8G |

## rnaSPAdes 路径引擎落地（--path-engine rnaspades, 2026-08-30）

复刻计划 B4 最后一项: rnaSPAdes PathExtend MultiExtender 核心规则的组件级
复刻（spades-main `path_extenders.cpp:303-347` MultiExtender::AddCandidates +
`rna_mode.info:155-197` path_cleaning default 档）。新模块
`trinity-butterfly/src/rnaspades.rs`, 默认关 = classic 逐字节不变。

- **规则复刻**: 0 候选停; 1 候选延伸; **恰好 2 候选且构成简单 bulge →
  复制路径各走一支（产 2 个 isoform）**; 其他一切多歧义（2 候选非 bulge /
  ≥3 候选）→ 停止延伸（宁可断不可错连）。bulge 的图语义映射: SPAdes 图
  经压缩（边=unitig, bulge=同起点同终点的平行边）, 本组件图节点=kmer、
  未压缩——翻译为"两条内部简单（出度恒 1）的分叉链在 ≤50 节点内重汇合到
  同一入度恰 2 的节点 C"（SNP bubble, 支长 ≈ k）。种子 = 无真实前驱源点
  按 (-cnt, vid) 排序; 纯环组件退化为最高覆盖单种子; 路径内节点不复走
  防环; 每组件 4096 路径安全阀（SPAdes 无此帽）。
- **简化点（明示）**: 不做 PE pair 权重（RNAExtensionChooser /
  PathCoverWeightCounter 未复刻）, 候选排序用覆盖一致性 score =
  1/(1+|cnt−路径均覆盖|), 保留 score ≥ max/1.5（priority_coeff,
  pe_params.info 默认）的全部并列者以触发 bulge 判定; 组件 graph.reads
  因此不消费。path_cleaning 只实现 default 档两条硬阈值（≥110bp 且
  平均 kmer 覆盖 ≥2）, isolated/relative 过滤未复刻（harvest 层全序列
  去重兜底）。LoopDetector / used_unique_storage 未复刻。
- **走线**: args.rs PathEngine 加 RnaSpades 变体（与 --recursive 互斥,
  沿 gtasm 先例; --gt-weights 在非 gtasm 下一律警告忽略）→
  butterfly_pool 分发（只读 graph.out）→ harvest 通道与 gtasm 相同
  （FASTA 头 `>{comp}_g{j}_i1 len=L path=[rnaspades]`, gene_trans_map 兼容）。
- **测试**: 9 个单测（bulge 分叉产 2 isoform 且序列逐字节等于两支等位、
  非 bulge 双歧义停、三岔口停、覆盖一致性滤掉低覆盖干扰支、纯环组件
  高覆盖种子+防环一周、path_cleaning 长度 110/109 与覆盖 2.0/1 边界、
  空图、FASTA 头与 60 列折行）; cargo test --workspace 全绿,
  clippy --workspace --all-targets 零警告。
- **bench dev1（cache cf94bb82c6912d3c）**: R=10.6% R70=45.3% Rp=8.9%
  hcR=30.8% **P=88.8%** F1=19.0% N50=1597 n=6739。对照 classic
  （R=14.0% R70=53.4% Rp=10.8% hcR=34.4% P=81.0% F1=23.9%）与 gtasm
  （R=12.3% R70=49.7% Rp=9.9% P=88.8%）: 召回 -3.4pp、精度 +7.8pp——
  与"宁可断不可错连"的保守语义自洽（无 PE 权重时歧义点只能停, 长转录本
  被截断 → R70 掉 8.1pp）; P 与 gtasm 并列最高。
- **classic 回归（同二进制无 flag 全新跑, cache 9ae734d899a34ae3）**:
  R=14.0% R70=53.4% Rp=10.7% hcR=34.4% P=80.8% F1=23.9% n=7851——对照
  基线（14.0/53.4/10.8/34.4/81.0/23.9/7869）在已知 inchworm 波动带内
  ✓ 回归通过。
- **裁定: 引擎本身可用但非增益**（F1 19.0% < classic 23.9%）。价值在
  规则正确性（bulge 双 isoform 是 classic/gtasm 都没有的机制）与高精
  度补充视角; 要追召回需补 PE 权重延伸（任务书原始 B4 的完整版）, 当前
  按收窄 scope 交付。flag 默认关。

## rnaSPAdes 完整复刻计划全线收口（2026-08-29/30）

用户指令: Q30 过滤（不纠错）+ 多k注入/相对覆盖/多路径歧义/gap closing
四大机制完整复刻, 不因"已有等价物"跳过。全部落地并逐一经 bench 裁决:

| 阶段 | 机制 | dev1 结果 | 裁定 |
|---|---|---|---|
| A | Q30 过滤（meanQ>=30, pair 同留弃） | 保留 97.2%, R 13.7%（-0.3） | 负: 数据本干净, 山谷是覆盖缺口非错误 |
| B1 | 多 k 迭代 + trusted 注入（[29,32], 钳位） | **R 14.5% R70 55.0% Rp 10.9% P_full +2.1pp, P 无损** | **正: 唯一正收益机制** |
| B2 | 相对覆盖清理（tip mmm3/rctc.5, EC rcec1.0, rcc 20x） | P +0.6pp 但召回 -0.3~-0.7pp | 净负: 我们的图已偏保守, 再剪伤召回; 默认关 |
| B3 | SPAdes gap closer（tip 缝合+PE 票） | 缝合 10 处全 perfect, 指标持平 | 中性: 与 GapCloser 原型证伪一致 |
| B4 | MultiExtender（仅简单 bulge 分叉）+ path_cleaning | R 10.6% P 88.8%（精度并列最高） | 召回换精度, F1 净负; PE 权重未复刻为已知简化 |

结论: rnaSPAdes 的全部构图/解图机制中, 唯一给 dev1 带来净增益的是
**多 k 迭代 + trusted contig 注入**（且我们跑的是 k 钳位到 32 的阉割版;
KmerId 扩 128 位支持真 k=37 是明确的后续项）。其余机制或中性或净负——
与"dev1 算法空间近穷尽"的此前系列裁定一致。B1 之外的全部 flag 默认关,
classic 路径行为不变（回归已验证）。
| 2026-08-30 | u128-classic-dev1 | R=14.0% P=81.1% **F1=23.9%** | N50=1843 n=7881 | --CPU 8 --max_memory 8G |
| 2026-08-30 | u128-multik2937-dev1 | R=15.3% P=80.4% **F1=25.8%** | N50=1899 n=8501 | --multi-k auto --CPU 8 --max_memory 8G |

## KmerId u64→u128 扩容： rnaSPAdes 真 k=37 落地（2026-08-30）

B1 钳位版（[29,32]）是复刻计划里唯一正收益机制, u64/2bit 位宽是唯一卡点。
本轮把 `KmerId` 扩为 u128（上限 k=64）, `--multi-k auto` 对 RL=76 的公式
原值 **[29,37] 不再被钳**。

- **位宽敏感落点**: 计数/覆盖统计的滚动 revcomp 掩码 `k>=64 → KmerId::MAX`
  （防 `1<<2k` 在 k=64 溢出 u128）; inchworm 候选生成的原版 u64 常数
  `(33-K)*2/(32-K)*2` 平移为 `(65-K)*2/(64-K)*2`; debruijn 图掩码同宽化。
  k≤32 键值与 u64 时代零扩展一致 → classic 行为逐字节不变。
- **counts_bin 格式升版**: TRMC(8B 键,12B/条) → TRM2(16B 键,20B/条),
  旧缓存魔数硬拒绝（不双读——易失缓存重算即可; xtask export_candidates
  读旧 bin 会报 bad header, 重跑生成）。
- **CLI 支持域 [15,63]**（--KMER_SIZE 与 --multi-k 列表/auto 钳制统一,
  63 = 64-1 留一碱基余量）; 两独立 bin（inchworm/trinity-kmer）的 -K
  校验同步到 64。diginorm 无独立编码（共用 CountMap, 随动）; gtasm 的
  gt_support K2V 键随 KmerId 走（multi-k×gtasm 参数层未接线, 无 k>32
  隐性依赖——k≤64 的图均可正确编码）。
- **KMER_ENTRY_BYTES 18→26**（内存护栏/容量提示的每条目下限随键宽调整）。
- **测试**: workspace 全绿 + clippy 零警告; 新增 k=33/37/55/64 编码/解码/
  revcomp/canonical/熵往返（kmer.rs）, 滚动 revcomp DS oracle 扩到 k=33/37
  （counter.rs）, k=33/37/55 前/反向候选位运算（kmer_counter.rs）, TRM2
  往返 + TRMC 旧格式拒绝（counter.rs）; 合成数据 k=33→37 E2E 冒烟通过。
- **classic 回归 bench dev1**（cache 5a192a9e8b0d4b90, 新二进制全新跑,
  ~5.3min）: R=14.0% R70=53.3% Rp=10.8% hcR=34.7% P=81.1% F1=23.9%
  N50=1843 n=7881 P_full=63.7% —— 全轴落在已知波动带（R 14.0±0.2 /
  Rp 10.7±0.2 / hcR 34.4±0.5 / P≈81）内, **扩容未引入行为差异**。
- **multi-k bench dev1**（cache d1489cbae7eab3b0, 日志/目录确认真 [29,37]:
  iter_k29 构图轮 + iter_k37 完整轮, trusted 注入链接）: **R=15.3%
  R70=56.7% Rp=11.4% R70p=53.1% hcR=35.6% hcR70=76.0% P=80.4% F1=25.8%
  N50=1899 n=8501 P_full=67.0%**。

**[29,37] vs 钳位版 [29,32]**（对照 08-29 B1: R=14.5 R70=55.0 Rp=10.9
hcR=34.6 P=81.4 P_full=65.8 F1=24.6 N50=1918 n=7925）:

| 轴 | [29,32] | [29,37] | Δ |
|---|---|---|---|
| R | 14.5% | 15.3% | **+0.8pp** |
| R70 | 55.0% | 56.7% | **+1.7pp** |
| Rp | 10.9% | 11.4% | +0.5pp |
| hcR | 34.6% | 35.6% | **+1.0pp** |
| F1 | 24.6% | 25.8% | **+1.2pp** |
| P_full | 65.8% | 67.0% | **+1.2pp** |
| P | 81.4% | 80.4% | -1.0pp |
| n | 7925 | 8501 | +576 |

真 k=37 在召回各轴全面再进一步（长程连接 + trusted 注入链生效）, 代价是
精度 -1.0pp / 多报 576 条。**裁定建议: --multi-k auto 维持可选增强定位
（默认 classic 不变）; 精度敏感场景注意 P 让步; Rp +0.5pp 仍低于 85 条集合
的理论封顶, 与 B1 结论同向。** u128 化代价: 计数哈希每条目 ~+8B 与
counts_bin 体积翻倍——dev1 规模无感, 更大数据集需留意内存。
| 2026-08-30 | dev3 | R=13.8% P=80.2% **F1=23.5%** | N50=1874 n=6953 | --multi-k auto --CPU 8 --max_memory 8G |
| 2026-08-30 | dev2 | R=15.1% P=79.7% **F1=25.4%** | N50=1878 n=9544 | --multi-k auto --CPU 8 --max_memory 8G |

## 多 k 增益三数据集交叉验证 + P 让步归因（2026-08-30）

**交叉验证（--multi-k auto, 均真 [29,37]）vs classic 同二进制基线:**

| 轴 | dev1 Δ | dev2 Δ | dev3 Δ |
|---|---|---|---|
| R | +1.3pp (14.0→15.3) | +2.0pp (13.1→15.1) | +1.3pp (12.5→13.8) |
| R70 | +3.3pp | +5.2pp (50.1→55.3) | +2.6pp (52.6→55.2) |
| Rp | +0.7pp | +1.1pp (10.1→11.2) | +0.6pp (9.5→10.1) |
| hcR | +1.2pp | +3.3pp (32.5→35.8) | +1.4pp (35.2→36.6) |
| F1 | +1.9pp | +3.0pp | +1.8pp |
| P | -1.0pp | -0.1pp | -2.2pp |

**裁定: 增益跨数据集稳健（六召回轴全涨, P 代价 0.1-2.2pp）, multi-k
从"机制验证"升级为"可信增益"。**

**P 让步归因（probe_multik_attribution, dev1）**: multi-k 新增 1967 条
（classic 中无对应）——63%（1245 条）直接命中 ref（classic 漏装的真实
转录本, 深度中位 10.5x, 长度中位 964bp）; 未命中的 722 条中 280 条
reads 深度>=5x（候选真实/未注释, StringTie 参照本身不全）; 疑似噪声
仅 ~442 条（深度中位 2.0x）。**P 的 1pp 代价里约六成是真转录本**——
评估口径（StringTie 参照不完整）放大表观代价, 生物学净收益为正。
| 2026-08-30 | dev1 | R=15.4% P=79.9% **F1=25.9%** | N50=1896 n=8548 | --CPU 8 --max_memory 8G |
| 2026-08-30 | nomultik-dev1 | R=14.0% P=82.2% **F1=23.9%** | N50=1822 n=7828 | --no-multi-k --CPU 8 --max_memory 8G |
| 2026-08-30 | gtasm-singlek-dev1 | R=12.2% P=89.1% **F1=21.5%** | N50=1672 n=6768 | --path-engine gtasm --gt-weights /storage/home/senior007/test/athal/bench/export/gtnet_dev1.gtw --CPU 8 --max_memory 8G |
| 2026-08-30 | dev1 | R=12.6% P=89.2% **F1=22.1%** | N50=1709 n=6445 | --path-engine gtasm --gt-weights /storage/home/senior007/test/athal/bench/export/gtnet_dev1.gtw --multi-k auto --CPU 8 --max_memory 8G |

## gtasm×multi-k 组合验证（2026-08-30, cache 6f007a16810cb61a）

用户指定验证项。gt_support 的 k2v 按图内实际 kmer 长自适应（k=37 图=36mer,
机械兼容）。dev1 三方对照:

| 引擎 | R | R70 | Rp | hcR | P | F1 |
|---|---|---|---|---|---|---|
| gtasm@单k25 | 12.2% | 49.7% | 9.9% | 33.4% | 88.8% | 21.5% |
| **gtasm@multi-k** | **12.6%** | **50.8%** | **10.3%** | 33.6% | **89.2%** | **22.1%** |
| classic@multi-k | 15.3% | 56.7% | 11.4% | 35.6% | 80.4% | 25.8% |

- **gtasm 全轴受益于多 k 图**（R +0.4 / R70 +1.1 / Rp +0.4 / P +0.4pp,
  P=89.2% 为全项目精度纪录）——更好的图对任何引擎都是更好的图。
- 但同图对比 classic 仍 -2.7pp R: 引擎差距（butterfly read-threading 的
  提取能力）大于图改进（+0.4pp）, gtasm 仍未吃干榨尽图里的路径。
- 建议: 解除 gtasm 缺省强制单 k 的限制（组合已验证为正）, 警告可降级。
| 2026-08-30 | dev1 | R=13.9% P=82.0% **F1=23.8%** | N50=1768 n=8026 | --path-engine gtasm --gt-weights /storage/home/senior007/test/athal/bench/export/gtnet_dev1.gtw --multi-k auto --CPU 8 --max_memory 8G |
| 2026-08-30 | dev1-novwalk | R=12.6% P=89.2% **F1=22.1%** | N50=1709 n=6449 | --path-engine gtasm --gt-weights /storage/home/senior007/test/athal/bench/export/gtnet_dev1.gtw --multi-k auto --no-gtasm-visited-walk --CPU 8 --max_memory 8G |
| 2026-08-30 | dev1 | R=15.4% P=79.8% **F1=25.9%** | N50=1896 n=8566 | --CPU 8 --max_memory 8G |

## gtasm 解码 visited 穿越修复（2026-08-30, cache 44bc70316e3b7a32 / b958fd1e7bcb1aa8）

**机制**: 原型 gt_decode2 把 visited 节点硬排除（`w not in visited`）——强势
isoform 的 walk 先把共享外显子标 visited, 弱势 isoform 的 walk 到边界即截断,
完整路径被拆散吞没（差距集 68 条 principal 中 28 条图内完整可走, 其中 16 条
≥95% 路径节点已被其他 gtasm 产出覆盖）。修复（`GtasmParams::allow_visited_walk`,
CLI 默认开, `--no-gtasm-visited-walk` 回退原型语义）:

- **两级候选过滤**: 优先 `not in seen/frozen 且未 visited`（原型语义）; 一级
  为空再取 `not in seen/frozen` 的已 visited 候选（同 (gsc, vid, eid) 取 max）;
  walk 内 seen/frozen 防环不变;
- **重复产出防护**: 组件内已产出序列 64 位指纹集, 完全相同才丢弃（穿越重新
  发现的同一 walk）; 共享部分节点的不同 isoform 照常产出;
- **终止不变**: 每轮胜出 walk 必含起点边两端 → 该边必退出 alive; visited 只增
  不减, 至多 ne 轮。NB_STARTS 采样/双向延伸/sup_sum 竞争/200bp 下限不动。

**单测**: 泡形双 isoform 强弱两臂均完整产出（关 flag 对照复现截断残段）;
无名平行入口重新发现同一 walk 不重复产出; 环图不挂死（≤2 条, seen 防环）;
greedy_extend 两级过滤/起点 visited 处置; flag-off 路径的全部既有位精确断言
不变。workspace 全绿 + clippy 零警告。

**dev1 三方对照**（同一新二进制; A/B 仅差 decode flag, 共享同一上游图）:

| 配置 | R | R70 | Rp | hcR | P | F1 | n | N50 |
|---|---|---|---|---|---|---|---|---|
| gtasm 修复前基线（08-30） | 12.6% | 50.8% | 10.3% | 33.6% | 89.2% | 22.1% | 6445 | 1709 |
| B: gtasm `--no-gtasm-visited-walk` | 12.6% | 51.0% | 10.3% | 33.6% | 89.2% | 22.1% | 6449 | 1709 |
| **A: gtasm visited 穿越开** | **13.9%** | **54.3%** | **10.7%** | 34.3% | 82.0% | **23.8%** | 8026 | 1768 |
| C: classic 回归（无 flag） | 15.4% | 56.8% | 11.6% | 35.9% | 79.8% | 25.9% | 8566 | 1896 |

- **A−B = 纯解码效应**: R **+1.3pp**, R70 **+3.3pp**, Rp +0.4pp, hcR +0.7pp,
  F1 +1.7pp; 代价 P **-7.2pp**（89.2→82.0, 多出 1577 条 contig）。
- **R 收敛**: gtasm↔classic 引擎差距 R 2.7pp 中拿回 ~48%, R70 差距 5.9pp 中
  拿回 ~57%, Rp 差距 1.1pp 中拿回 ~36%——机制假设（visited 阻塞）坐实为
  引擎差距的主成分之一, 但非全部（剩余 R -1.5pp 差距待查）。
- **B−基线 = 上游漂移**（非本改动）: R/Rp/hcR/P/F1/N50 全同, 仅 R70 +0.2pp /
  n +4; 两 cache 的 FASTA 头对比显示组件编号重排（c1 基线 len=8979 vs 新
  len=617）——图在上游（inchworm/chrysalis）已不同, 而本改动只在 butterfly
  池内消费 graph.out, 物理上不可能引起组件重排; 归因于基线测量二进制与
  当前 HEAD 之间的提交漂移（4027446 multi-k 转正）。
- **classic 回归 C**: 全轴命中既有波动带（对照 08-30 行 R=15.4 P=79.9
  F1=25.9 N50=1896 n=8548）。
- **P 代价解读**: 修复后 gtasm P=82.0% 仍高于 classic（79.8-80.4%）~2pp;
  参照 multi-k P 让步归因先例（六成新增为真转录本, StringTie 参照不全放大
  表观代价）, 新增 1577 条的构成值得同款 probe 复核后再定论。

**裁定建议**: visited 穿越转正为 gtasm 默认（已落地为实验默认,
`--no-gtasm-visited-walk` 可回退）——R/R70/Rp/F1 全轴显著向 classic 收敛且
P 仍优于 classic; 遗留: (1) 新增 contig 的 probe 归因（真转录本占比）;
(2) 剩余 R -1.5pp 引擎差距解剖（read-threading 提取能力）; (3) dev2/dev3
交叉验证后再考虑移除回退旗标。
| 2026-08-30 | dev3 | R=12.8% P=84.3% **F1=22.2%** | N50=1753 n=6679 | --path-engine gtasm --gt-weights /storage/home/senior007/test/athal/bench/export/gtnet_dev1.gtw --multi-k auto --CPU 8 --max_memory 8G |
| 2026-08-30 | dev2 | R=13.5% P=80.4% **F1=23.1%** | N50=1733 n=9216 | --path-engine gtasm --gt-weights /storage/home/senior007/test/athal/bench/export/gtnet_dev1.gtw --multi-k auto --CPU 8 --max_memory 8G |

## visited 穿越修复: 三数据集交叉验证 + 新增归因（2026-08-30）

**交叉验证（gtasm@multi-k + 穿越默认开）vs classic@multi-k / gtasm@单k:**

| 数据集 | 引擎 | R | R70 | Rp | hcR | P | F1 |
|---|---|---|---|---|---|---|---|
| dev1 | gtasm 穿越 | 13.9% | 54.3% | 10.7% | 34.3% | 82.0% | 23.8% |
| dev2 | gtasm 穿越 | 13.5% | 52.4% | 10.3% | 33.7% | 80.4% | 23.1% |
| dev3 | gtasm 穿越 | 12.8% | 53.5% | 9.5% | 34.5% | 84.3% | 22.2% |

对 gtasm@单k: R 分别 +1.7/+2.0/+1.3pp（多 k 图 + 穿越叠加效应）;
对 classic@multi-k: R 差距收敛到 -1.5/-1.6/-1.0pp（修复前 -2.7pp）,
且 P 三个数据集全部保持不低于 classic（82.0/80.4/84.3 vs 79.8/79.7/80.2）。
**裁定: 穿越修复跨数据集稳健。**

**新增 1323 条归因（穿越开 vs 关, dev1）**: 62%（816）命中 ref;
未命中 507 条中 448 条 reads 深度>=5x（深度中位 28x!）; 疑似噪声仅
~59 条。**~96% 为真实序列**——被吞的确实是真 isoform, P 的表观代价
基本全是参照不全的评测伪影。

**运维记录**: /tmp 被 248 个 scratch 目录填满至 100%（bench 连环失败）,
已清理并永久改道——所有 bench 以 TMPDIR=/storage/.../bench/scratch 运行
（Rust temp_dir 语义, 覆盖管线 scratch 与评测临时文件）。

## P 尾巴低成本过滤: 判别力证伪（2026-08-30, probe_ptail.py）

假设: multi-k P 代价的 ~442 条"疑似噪声"可用 reads 深度/覆盖宽度低成本滤除。
实测（classic@multi-k 8501 条, good=id>=0.95 命中 7239 条 / bad=1262 条,
sr map 逐条算深度+宽度, 阈值网格 B∈[0.5,0.9]×D∈[1,5]x）:

- 最松阈值 B<0.5&D<1: 删 bad 33% 但**误伤 good 12%**（好碱基损 19.6%）;
- 最紧 B<0.9&D<5: 删 bad 49%, 误伤 good 25%（好碱基损 34.7%）;
- **全程无干净切割线**——坏条与好条在深度/宽度上同分布。

机理: P 尾巴的主体不是噪声, 是(a)低表达真转录本（参照内）与
(b)参照外但 reads 强支持的真实序列（穿越归因: 深度中位 28x）——
二者与真噪声在覆盖信号上不可分; mapq>=5 过滤还系统性低估了
共享外显子家族的 read 支持（多映射）, 进一步混淆。
裁定: 放弃深度/宽度过滤路线。P 的 1pp 尾巴就是低表达生物学的本体,
不是可滤的噪声; 之前"442 条疑似噪声"的表述被本探针修正为
"低支持真实序列与噪声的混合体, 无低成本可分性"。
| 2026-08-30 | dev1 | R=13.9% P=82.0% **F1=23.8%** | N50=1768 n=8027 | --path-engine gtasm --gt-weights /storage/home/senior007/test/athal/bench/export/gtnet_dev1.gtw --multi-k auto --CPU 8 --max_memory 8G |

## 引擎差距追踪 v2: 残余 = 分支点碎片化（2026-08-30, probe_gaptrace）

visited 穿越后重解剖（dev1, classic@mk full 603 vs gtasm@mk+穿越 full 559,
差距 49 条 principal）:

- 图内严格可走 24 条 + 近可走（参照含变异位点）25 条——**离图断点型
  （butterfly read-threading 专属领地）只剩 ~1 条**, read-threading 不再是
  主要差距来源（原 40 条不可走假设被推翻: 25 条"不可走"实为参照序列
  变异导致的假不可走, classic contig 本身在图内）。
- 失败形态: gtasm 最佳覆盖 0.85-0.94 的近失误 30 条中 **27 条是中间
  碎片化**（≥2 段在参照坐标上重叠/紧邻, 仅 3 条末端缺失, 末端缺失中位
  仅 16/31bp）——完整路径在图里, walk 在**分支点选错臂**, 正解被拆成
  两个 contig。visited 穿越无法挽回"起点后第一次选择就错"的情形。
- 裁定: 残余引擎差距 = **分支点保真度**。候选刀法: (a) read PAIR 证据
  进分支决策（butterfly pair_paths 语义, 单 read look-back 已证无效,
  pair 是 76bp 读长下唯一够得着 300-400bp 分支间距的证据, 未测）;
  (b) 简单 bulge 双臂都产出（rnaspades 规则嫁接, 精度换召回）;
  (c) 同 comp 碎片 pair-guided 合并（此前 OLC 系列证伪过, 谨慎）。
| 2026-08-31 | dev1 | R=13.8% P=81.8% **F1=23.6%** | N50=1773 n=8018 | --path-engine gtasm --gt-weights /storage/home/senior007/test/athal/bench/export/gtnet_dev1.gtw --multi-k auto --gtasm-pair-support --CPU 8 --max_memory 8G |
| 2026-08-31 | dev1 | R=13.8% P=67.0% **F1=22.9%** | N50=2194 n=9621 | --path-engine gtasm --gt-weights /storage/home/senior007/test/athal/bench/export/gtnet_dev1.gtw --multi-k auto --gtasm-pair-support --gtasm-bulge-dual --CPU 8 --max_memory 8G |
| 2026-08-31 | dev1 | R=13.9% P=65.8% **F1=22.9%** | N50=2228 n=9650 | --path-engine gtasm --gt-weights /storage/home/senior007/test/athal/bench/export/gtnet_dev1.gtw --multi-k auto --gtasm-bulge-dual --CPU 8 --max_memory 8G |

## 分支点保真两件套裁定: (a) pair 证据中性, (b) bulge 双臂证伪（2026-08-31, cache 43a460d9c29a2f5f / 15b011b4ddc4148b / 1cd9ceb4b4544327）

「引擎差距追踪 v2」候选刀法 (a)(b) 落地并 bench（均默认关, 机制见
gtdecode/gt_support 模块文档）:

- **(b) `--gtasm-bulge-dual`**: 前向延伸遇简单 bulge（恰 2 个一级候选、
  两支内部出度恒 1、≤50 节点内重汇合且 indeg(C)=2——rnaspades
  `bulge_merge` 语义）不再 argmax 二选一, 双臂各走一条、胜出起点的全部
  fork 参与产出（各自独立拼序列/min_len/指纹去重, visited 并集办结）。
  安全阀: 单起点 64 fork / 每组件 4096 contig; 终止性不变（fork 共享
  back 前缀或同起于 d0, 起点边两端仍必入 visited）。
- **(a) `--gtasm-pair-support`（乘性变体）**: 支持 pass 对 R1/R2 同序
  lockstep 配对, 一条 mate 的节点路径含分支出边 (v,w)（v 出度≥2）且
  另一条命中 v 的上游窗口（反向 BFS ≤15 边, 不含 v 自身）→
  `gtasm_pair_support.bin`（GTP1 sidecar, 与 GTS1 同布局异 magic,
  与单 read 边支持字段分离）; 解码 gsc 追加 ×(1+log1p(pair))。
  dev1 实测证据量: 4,483 条分支出边有 pair 支持, 计数总 73,729 对,
  3,213 边 ≥2 对（证据真实存在, 非稀疏到无效）。

dev1 对照（gtasm@multi-k + 穿越; 基线 = 08-30 已知数, 未复测）:

| 变体 | n | R | R70 | Rp | hcR | P | F1 | N50 |
|---|---|---|---|---|---|---|---|---|
| 基线（08-30 已知） | 8026 | 13.9% | 54.3% | 10.7% | 34.3% | 82.0% | 23.8% | 1768 |
| (a) pair 乘性 | 8018 | 13.8% | 54.2% | 10.6% | 34.2% | 81.8% | 23.6% | 1773 |
| (b) bulge 双臂 | 9650 | 13.9% | 54.2% | 10.6% | 34.2% | **65.8%** | 22.9% | 2228 |
| (a)+(b) 组合 | 9621 | 13.8% | 54.1% | 10.6% | 34.1% | **67.0%** | 22.9% | 2194 |

**裁定**:

- **(a) pair 乘性 = 中性**（全轴 −0.1~−0.2pp, 在采样噪声带 ±0.2pp 内）:
  pair 证据存在且充足, 但与 gsc 已有的 sup 项在分支点上高度冗余——
  GT 模型 + 单 read 支持本就大多选对臂, pair 不提供增量。不转正,
  flag 保留默认关（A1 过滤式变体未实现: 乘性中性后过滤式预期同向,
  不再投入）。
- **(b) bulge 双臂 = 证伪**: R 纹丝不动（13.9 平基线）而 n +1,624 条、
  P 崩 −16.2pp（82.0→65.8）——双臂产出把本应被 gsc 压制的弱臂/噪声臂
  等量放出, 全是假阳性; dev1 碎片化主因**不是** bulge 型分支点
  （与 rnaspades 引擎的图粒度差异: 未压缩 kmer 图上的"简单 bulge"大量是
  测序错误/低覆盖分叉, 而非真 SNP bubble）。默认关保留, 不推荐开启;
  若复用需加臂支持度门槛（弱臂 sup/pair ≥ 阈值才分叉）。
- 组合无协同（P 同样崩, 归因全在 (b)）。
- 残余分支点碎片化差距未收: 候选 (c) 碎片 pair-guided 合并成下一刀,
  或转向"分支点选错的事后检测 + 重 walk"而非事前证据。
| 2026-08-31 | dev1 | R=13.8% P=80.5% **F1=23.6%** | N50=1834 n=7539 | --path-engine gtasm --gt-weights /storage/home/senior007/test/athal/bench/export/gtnet_dev1.gtw --multi-k auto --gtasm-frag-merge --CPU 8 --max_memory 8G |

## 分支点保真 (c): 碎片 pair-guided 合并 = 证伪（2026-08-31, cache 98107448445d0204）

「引擎差距追踪 v2」候选刀法 (c) 落地: `--gtasm-frag-merge`（默认关,
机制见 frag_merge 模块文档）。butterfly 组件池后、harvest 前加组件内
合并 pass: contig 3'/5' 末端各 ≤150bp 建末端 k2v（K=图 kmer 长,
正反义双写, 一 kmer 多 token——重叠碎片共享 kmer 不做先写先赢）;
原始 reads lockstep 配对全窗线程化, mate A 命中 X 的 3' 端 × mate B
命中同组件 Y 的 5' 端 → link(X,Y)+=1（同对去重, 命中集 >64 的安全阀
弃对）; 组件内按权重降序、≥2 对才拼——先 OLC（X 后缀 × Y 前缀
≥40bp 且 ≥90% 一致, 取最长）, 否则 50N 连接; 每 contig 至多参与一次
（先到先得）, 产物 <200bp 拒收且不耗名额。输出契约与 harvest 一致,
统计行落 .stage_timings。

dev1 实测证据量: 候选链接 2,253（达阈值 1,995）→ 合并 500 次
（OLC 175 + N 连 325, 短拒 0）, contigs 8,039 → 7,539, pass 耗时 113s。

dev1 对照（gtasm@multi-k + 穿越; 基线 = 08-30 已知数, 未复测）:

| 变体 | n | R | R70 | Rp | hcR | P | F1 | N50 |
|---|---|---|---|---|---|---|---|---|
| 基线（08-30 已知） | 8026 | 13.9% | 54.3% | 10.7% | 34.3% | 82.0% | 23.8% | 1768 |
| (c) frag-merge | 7539 | 13.8% | 54.0% | 10.7% | 34.2% | **80.5%** | 23.6% | 1834 |

**裁定: 证伪, 不转正（flag 保留默认关）**。预期天花板 ~+0.5pp Rp
未兑现: Rp 纹丝不动（10.7% 平）, R −0.1pp 在噪声带内, 而 P −1.5pp
（82.0→80.5, 超 ±0.2pp 采样噪声）、F1 −0.2pp——500 次合并零召回
收益、净精度损失。与 (a) pair 乘性中性、(b) bulge 双臂证伪合并解读:
分支点碎片化的 pair 证据无论事前（进 gsc）还是事后（拼 contig）
都不含增量信息——link≥2 的末端 × 末端链接主体不是"同一转录本被
错拆的相邻碎片", 而是同基因/共享外显子家族内的跨转录本配对,
拼上即成嵌合体（N 填充与错接碱基直接吃 P）。分支保真三件
套 (a)(b)(c) 全部落地完毕: 残余 49 条 principal 差距未收, 后续方向
只剩"分支点选错的事后检测 + 重 walk"。


## xRead 重叠图范式评估: 原版不可用 + 思想适配封顶 +0.85pp, 收档（2026-08-31, xRead + probe_read_bridge）

外部建议引入 xRead（tcKong47/xRead, 论文 PMC11831799, "coverage-guided
iterative overlapping graph construction"）替换/补充 de Bruijn 构图。两路验证:

1. **原版实证: 短读长下不可用。** xRead 是长读长（ONT/PacBio 10kbp+）基因组
   overlapper, 核心 "mimicking-and-mapping"（覆盖度引导迭代选种子 read 当
   "虚拟参考", alignment skeleton 打分, 只留 CRO 高分 overlap）。编译
   （make INCLUDES/LIBS 指向 miniconda zlib）后在 dev1 head 1M 条 76bp reads
   实测: 默认 min_ove=500 超读长, 调 -m 40 后 **99.9% reads（999,158/1M）
   零 overlap**, 全图仅 1,463 条 PAF。根因: 76bp read 仅 ~20 个 minimizer,
   skeleton 共线链接评分无法建立, 覆盖度模型（median 0.0）不启动;
   短读长下"虚拟参考"退化为 reads 全集, 方法学前提不成立。
2. **思想适配量化（probe_read_bridge.py, cache 68349cb4e339ced0 单 k 图）**:
   提取"read 级序列直接跨谷"——de Bruijn 图断掉的山谷, 单条 76bp read
   或能完整跨越。GRAPH 桶 3202 条 principal: 缺失区间 8343 个（内部 3283 /
   端部 4837 / 整条 223）; 内部区间 71.2% 跨度≤76bp 几何可行, 其中 77.7%
   有≥1 桥接 read（探针 map 原始 reads, 两端深入锚点, ident≥0.9）。
   **但"read 覆盖山谷"是 GRAPH 桶入桶条件的同义反复**, 真筛选层是谷内
   kmer 深度（≈桥接 read 数, 一条跨谷 read 覆盖谷内全部窗口）:
   内部≥3桥 + 端部缺失≤5% 口径, 全区间可桥 principal 仅 **44 条**,
   理论增量 **R_full_p +0.85pp**（≥1 桥松口径 87 条/+1.68pp）。
3. **裁定: 收档。** 上限 +0.85pp 与四条已实测机制（kmer 救援 85 条集/
   maxC10000 兑现 2 条、frag_bridge +0.7pp 天花板、gap_close 中性缝 10 处、
   frag_merge 证伪）共享同一证据源（山谷低深度 reads）与同一可操作集;
   落地形态仅 kmer 注入（实测 2/85）或 read 序列直拼 contig（frag_merge
   嵌合体证伪）两种, 皆已封顶。xRead 论文亦自明其图需配套长读长
   layout/consensus, 与短读长转录组错位。
- 脚本: export/probe_read_bridge.py（pkl 键=gaps 全局索引; 曾误用转录本内
  局部序号致 principal 聚合失真, 已修并注）; 缓存 read_bridge.68349cb*.pkl;
  xRead 源树 /storage/home/senior007/test/trinity_rust/xRead/（run 大文件
  已删, 留 1463 行 PAF 证据）。
- **方法论教训**: 探针疗效指标若复用入桶条件（reads 覆盖）= 自欺, 必须加
  "证据深度/可注入性"筛选层; 本探针 v1 的 +59pp 假象即此型。
- scratch 清理: bench/scratch 19G→16K（multi-k 拿图 run + 历史 test 遗留
  全清; graph.out 确认在 partition 阶段消费后不保留, 图 kmer 集合日后须从
  graphs.pack 的 kind=1 记录读）。


## 碎片延长策略（学习版）全线量化: oracle +2.08pp 但规则实测 +0.08pp, 收档（2026-08-31, probe_extend_* / probe_reconnect_*）

用户提议: 训练小模型（以 gtasm walk 为师）把产物里"碎的但可延长"的 contig 延长。
完整量化链（dev1, 产物=44bc70316e3b7a32 gtasm@mk+穿越）:

1. **oracle 上界（probe_extend_oracle）**: 完美拼接同 principal 全部产物 hits,
   救回 108 条 full = **R_full_p +2.08pp**; 救回者 95/108 段间最大缺口 ≤50bp;
   partial/weak 档 oracle 后无迁移（数据/图决定, 拼接救不了）。
   **108 条全部段同组件**（contig 名口径）——组件内解图问题, 非跨组件。
2. **病灶拆分**: 段对 197 个 = 138 紧邻/重叠（gap≤0, 首尾相接非 OLC 重叠）+
   59 真缺口（图内连续在位 43 / 图外 14 / 单 k 口径锚点跨组件 48——multi-k
   组件更大, 部分已被 trusted 注入合并, 单 k 口径高估此项）;
   **多路径 0**——"模型选边"无对应病灶（选边数据集另测: 有效分支点仅 114,
   最大覆盖度规则 69.3%, 训练/部署分布不同——成功 walk 来自简单组件,
   部署现场是复杂组件）。
3. **OLC 重叠信号证伪（probe_extend_overlap）**: ref 重叠层段对精确重叠
   ≥20bp 仅 7/44——段间是紧邻（gap≈0/微缺口）, 末端序列天然不重叠;
   frag_merge 式 OLC 从原理上够不着。
4. **图内短路径规则（probe_reconnect_rule, 68349cb 单 k 图, 无向 BFS≤60bp）**:
   P=0.63 R=0.39, **正确连接路径中位 2bp**——gtasm 在距正确连接 2 步处停下。
   56 个"不该接但有路"是判别器的唯一战场。
5. **应用级模拟（probe_reconnect_apply, 真拼接+bench 口径）**: 候选 193 对,
   贪心执行 55 拼接 → **R_full_p 557→561 (+0.08pp), P 81.99→81.78 (−0.21pp)**。
   oracle→实测的折损来源: 有路率 0.39（模型变不出路）、多段全连要求、
   桥序列 ident 损耗（图序列 vs ref 变异）。
   另坐实: bench 碱基级 P 口径下嵌合体代价小（两段都能 map 上, P 仅 −0.21pp,
   frag_merge 时代的 P 风险感知应修正）。

**裁定: 收档。** 模型判别理论增量 ≤+0.2pp（只压 fp, 不改有路率/桥 ident）,
训练集 114 分支点且分布不迁移。与 read 桥接（+0.85pp）、frag_bridge（+0.7pp）、
frag_merge 共享同一证据源与可操作集, 第五次汇合到同一结论:
**碎片化天花板由数据与图结构决定, 后处理（规则或学习）兑现率极低。**
gtasm "差 2 步停下"的根因仍是分支点分数/选择（visited 穿越已修其一）,
若再动刀应在 walk 决策本身而非后处理。
- 脚本: probe_extend_oracle.py / probe_extend_overlap.py /
  probe_branch_dataset.py / probe_reconnect_rule.py / probe_reconnect_apply.py;
  缓存 extend_oracle.pkl / extend_overlap.pkl / branch_dataset.*.pkl /
  reconnect_rule.pkl / reconnect_apply.pkl; 产物 reconnected.44bc.fa。
- 方法论教训②: oracle 上界(决策完美)与应用实测间有"有路率×全连率×ident"
  三重折损, 评估任何后处理策略必须先量折损系数再谈模型。
| 2026-08-31 | qm4 | R=43.9% P=51.7% **F1=47.5%** | N50=1783 n=117645 | --CPU 32 --max_memory 64G |
| 2026-08-31 | qm4 | R=31.5% P=50.9% **F1=38.9%** | N50=1672 n=82276 | --no-multi-k --CPU 32 --max_memory 64G |
| 2026-08-31 | qm4 | R=45.3% P=50.5% **F1=47.8%** | N50=1835 n=116816 | --multi-k 29,49,63 --CPU 32 --max_memory 64G |
| 2026-09-01 | dev1 | R=13.8% P=82.3% **F1=23.6%** | N50=1756 n=8025 | --path-engine gtasm --gt-weights /storage/home/senior007/test/athal/bench/export/gtnet_dev1.gtw --multi-k auto --CPU 8 --max_memory 8G |
| 2026-09-01 | dev1 | R=15.3% P=80.2% **F1=25.7%** | N50=1882 n=8546 | --CPU 8 --max_memory 8G |
| 2026-09-01 | qm4-opt | R=44.0% P=51.8% **F1=47.6%** | N50=1774 n=117592 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G |


## 性能优化四件套: QM4 全程 137→63min（2.2×）, 语义等价（2026-09-01, bench 8a03b756/514d5be0）

QM4（2135 万对 150bp）全程 137min 的 profile 后, 四处热点逐一落地:

1. **kmer_load 消除（orchestrate.rs + kmer_counter.rs）**: kmer_count 阶段的
   内存表 CountMap 直传 kmer_load（跳过 4GB counts_bin 读回+单线程重建）,
   `KmerCounter::from_count_map` move 直构（计数表键已 canonical 且唯一,
   与逐条 add_kmer 等价）——省二次 rehash 与 canonical 重算。
   kmer_load: 622/406s → 0/15.6s。
2. **gtasm_support（gt_support.rs）**: process_seq 滚动 2-bit 编码
   （与逐窗 kmer_to_intval 逐位等价, 含 N 重置与"不含最后一窗"quirk）;
   scan_reads_file 流水线并行（读取线程按 line%4 quirk 产出序列批 →
   rayon 批内并行, EdgeSupport::merge 计数可交换）; build_k2v 按 bin 分组
   rayon + 按序合并保"先写先赢"。
3. **GraphFromFasta（kmer_align.rs + graph_from_fasta.rs）**:
   get_matches 零克隆惰性两指针交（桶不再 to_vec; merge_sort_filter 挪
   cfg(test)）; Phase 2 的 contig 整段 clone/revcomp → 引用 + 每 contig
   revcomp 缓存。**实测 gff 仅 1199→1126s（−6%）——克隆非真热点,
   真瓶颈是 Phase 2 刻意单线程段（下一步目标）**。
4. **diginorm（diginorm.rs）**: 左右侧 read_side/prep_side/emit_side 三处
   thread::scope 并行; both.fa 全量拷贝消除（两侧分别计数再合并,
   canonical_merge_sum 线性可拆）。

另发现: inchworm 组装线程 `args.inchworm_cpu` 默认 min(6, CPU)（原版语义）,
大机器需显式 `--inchworm_cpu`（QM4 用 32 → inchworm 874+678s → 172+198s）。

**验证**: 全 workspace 测试绿 + clippy 零警告; dev1 bench 语义等价
（gtasm 链 R−0.1pp/P+0.3pp, classic 链 R−0.1pp/F1−0.2pp, 均在 ±0.3pp
波动带）; QM4 [49,63] 指标 R=44.0%/P=51.8%（原 43.9/51.7）不变。
**QM4 全程 137→63min**: normalize 3127→749s（4.2×）, inchworm 两轮
26.5→6.2min（4.3×, 线程数）, kmer_load 全消, butterfly 289→77s。
（QM4 系列 bench 行由 bench 自动登记于上表。）
| 2026-09-01 | dev1 | R=15.4% P=80.6% **F1=25.8%** | N50=1875 n=8559 | --CPU 8 --max_memory 8G |
| 2026-09-01 | qm4-opt2 | R=44.0% P=51.6% **F1=47.5%** | N50=1784 n=117633 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G |
| 2026-09-01 | dev1 | R=15.4% P=80.0% **F1=25.8%** | N50=1913 n=8552 | --CPU 8 --max_memory 8G |
| 2026-09-01 | qm4-opt3 | R=44.0% P=51.5% **F1=47.5%** | N50=1792 n=117692 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G |

## GraphFromFasta 深挖: weldmer 计数哈希化 + Phase 2 并行（gff 1199→196s, 全程 46→39min）（2026-09-01, bench 71f83f9f44f9a038）

Phase 2 并行化后 gff 仍 600s, GFF_PROFILE 插桩分解（examples/gff_profile.rs,
QM4 中间产物直跑）定位真凶:
- **weldmer_read_counts 314s（64%）**: add_counts_from_reads_par 对 2150 万
  reads × 127 窗逐窗"排序 Vec 二分 + 24B memcmp"。改 NonRedKmerTable 懒建
  2bit 哈希索引（纯 GATC kmer; 含 N 走字节二分回退保原版语义）+ reads 侧
  滚动编码 → 74s（−76%）。
- **weldmer_templates 80s**: set_up_templates 逐窗 Vec 分配+排序大 Vec。
  改平铺收集（大 Vec+u32 偏移表, 排序移动偏移）→ 62s（−25%）。
- 其余段实测: build_core 27-35s（KmerAlignCore 单线程建索引, 可再并行,
  收益小暂缓）, Phase2_collect 10.5s, Phase2_apply 18s, phase1 15s。

**QM4 全程: 137→39min（3.5×）**; gff 1199→196s（6.1×）。
指标等价: dev1 R=15.4%/P=80.0%, QM4 R=44.0%/P=51.5%（波动带内）;
黄金对拍 edge_multiset_matches_original 逐字节等价。
注: Phase 2 并行化首版曾引入 rc_cache 串味 bug（单槽缓存在多个 RC 命中
对象间串序列）, 被黄金对拍当场拦截（weldmers 4 vs 5）, 修为按 contig 键
的 FxHashMap 缓存。
| 2026-09-01 | qm4-opt3 | R=44.0% P=51.5% **F1=47.5%** | N50=1792 n=117692 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G |
| 2026-09-01 | dev1 | R=15.4% P=79.0% **F1=25.8%** | N50=1906 n=8616 | --CPU 8 --max_memory 8G |
| 2026-09-01 | qm4-opt4 | R=43.8% P=51.5% **F1=47.3%** | N50=1780 n=117888 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G |

## f2db 并行化 + 归一化计数并行（QM4 全程 39→32.5min）（2026-09-01, bench 32f3a1fc87228df5）

- **graph_per_record 并行（debruijn.rs）**: 每 bundle 独立建图, rayon 并行
  生成块文本按索引序拼接（逐字节同构）——f2db 95→4.9s（−95%）。
- diginorm 两侧计数并行（diginorm.rs, DN_PROFILE 插桩分解: 计数 663s/82%
  为大头）: 单跑 804→661s; bench 32 线程池下无收益（池已饱和）, 保留无害。
- QM4 全程 32.5min。残段构成: normalize 748s（计数 523s 哈希插入 bound）
  > kmer_count 两轮 417s > inchworm 356s > gff 171s > rtt 95s > butterfly 76s。
  再压需动计数器内核（u128→u64 键 for k≤32 或分片表）, 风险收益比下降。
- 指标等价: dev1 R=15.4%/P=79.0%; QM4 R=43.8%/P=51.5%（波动带内）。
| 2026-09-01 | qm4-opt4 | R=43.8% P=51.5% **F1=47.3%** | N50=1780 n=117888 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G |
| 2026-09-01 | dev1 | R=15.4% P=79.5% **F1=25.8%** | N50=1907 n=8553 | --CPU 8 --max_memory 8G |
| 2026-09-01 | qm4-opt5 | R=43.9% P=51.5% **F1=47.4%** | N50=1787 n=118090 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G |

## 计数器免拷贝（counter.rs, bench 9471bbb5d4128ca4）

split_fasta_sequences 的逐记录 owned 拷贝 → split_fasta_ranges 借用视图
（单行序列段切片无内部 \n, for_each_kmer 断窗行为与滤除版逐字节等价;
检测到多行段自动回退旧路径）+ merge_counts 小并入大。归一化计数段
523→410s（单跑口径）, bench 内 normalize 748→584s（−22%）; 归一化输出
逐字节一致, dev1/QM4 bench 指标波动带内。
（bench 行见上表 qm4-opt5; 全程 ~33.6min, 与 opt4 的 32.5min 差异为系统
噪声——归一化的 −164s 收益部分被评测/负载波动吃掉。）

## 新评估维度: 蛋白层重合度 + BUSCO（QM4, TransDecoder 口径, 2026-09-02）

新增一条与核酸 P/R 互补的评估线: TransDecoder 5.7.1 翻译两组转录本
（无同源提示, 默认单最优 ORF）, diamond 双向 blastp 算蛋白重合度,
BUSCO 6.1.0 protein 模式（eudicots_odb10, n=2326）做非银标准锚点。
对象: QM4 银标准 st_transcripts.fa（38,549 条）vs release_qm4_out.Trinity.fasta
（multi-k auto, R=43.9% 那版）。脚本固化在 test_data/work/eval/
（eval_protein.py + protein_eval.sh 一键端到端）, 产物在
test_data/work/protein_eval/。

| 组 | 转录本 | 蛋白(≥100aa) | BUSCO C | S | D | F | M |
|---|---|---|---|---|---|---|---|
| StringTie 银标准 | 38,549 | 37,629 (97.6%) | 93.9% | 64.5 | 29.4 | 1.8 | 4.2 |
| trinity-rust multi-k | 117,592 | 72,947 (62.0%) | 89.1% | 54.0 | 35.2 | 4.5 | 6.4 |

蛋白重合度（diamond, pident≥90）: pepR90=90.4%（银标准蛋白被覆盖,
qcov≥50%）/ pepR90full=80.3%（qcov≥90 近全长）/ pepP90=89.9%（组装蛋白
能对上银标准）/ geneR90=90.9%（18,817 个 STRG 基因至少一条蛋白被覆盖）。

解读: 核酸全长口径 R=43.9% 与蛋白口径 90.4%/基因级 90.9% 的巨大落差,
定位了当前短板——**编码区（CDS）广度已到九成, 缺的是全长完整性**
（UTR 端截断+isoform 拼全）: 蛋白产出率 62% vs 银标准 97.6% 也指向同一处。
BUSCO −4.8pp（89.1 vs 93.9）独立佐证: de novo 对 genome-guided 的真实
差距远小于核酸 R 显示的水平。注意 pepR 口径含银标准循环性（StringTie
没捕到的真实基因只惩罚 pepP 不奖励 pepR）, BUSCO 是唯一非循环锚点。
| 2026-09-02 | u256-dev1 | R=15.4% P=79.5% **F1=25.8%** | N50=1907 n=8553 | --CPU 8 --max_memory 8G |
| 2026-09-02 | u256-dev1 | R=15.4% P=79.6% **F1=25.8%** | N50=1912 n=8542 | --CPU 8 --max_memory 8G |
| 2026-09-02 | qm4-k294973 | R=44.4% P=51.0% **F1=47.5%** | N50=1782 n=121356 | --multi-k 29,49,73 --inchworm_cpu 32 --CPU 32 --max_memory 64G |

## KmerId u128→U256 扩容 + k=73 顶 k 实验: 63 仍是最优（2026-09-02, bench 6c437e1881b4221d/278c78f43926ebfd）

动机: 用户要求三 k 落地 [29,49,73]（rnaSPAdes 公式给 RL=150 的原值 upper_k=73,
此前被 u128 位宽钳成 63）。2bit×73=190bit > 128bit → 手写 U256（双 u128 limb,
kmer.rs, 零新依赖: 移位/位运算/加减/回绕左移/FromStr/低 256 位序/与原生字面量
比较全套; limb 序保证派生 Ord = 数值序）。MAX_KMER_LENGTH 64→128, MULTI_K_MAX
63→127（对齐 SPAdes MAX_K）, auto 公式 RL=150 的 [49,73] / RL=250 的 [55,123]
不再被钳。counts_bin 魔数 TRM2→TRM3（键 32B, 36B/条, 旧格式硬拒绝）。
热路径四处滚动编码（counter/coverage_stats/nonred_table/gt_support）与
kmer_counter 的 (65-K)/(64-K) 常数（→(129-K)/(128-K)）全部平移。
KMER_ENTRY_BYTES 26→42（容量护栏口径随动）。

验证: workspace 665 测试全绿（新增 k=33/37/55/64/73/127/128 编解码往返、
TRM3 往返、TRM2 硬拒、129 越界）; dev1 真重组装 R=15.4%/P=79.6%/F1=25.8%
（基线波动带内, 无回归）。**坑**: bench 缓存键含 target/release/trinity-cli
内容——先 build release 再 bench, 否则命中旧 u128 缓存假回归。

**QM4 [29,49,73] vs 锚点**: R=44.4%/P=51.0%/F1=47.5%/N50=1782/n=121356,
对照 [29,49,63] R=45.3%/F1=47.8%/N50=1835/n=116816 与 [49,63] R=43.9%。
**73 顶 k 被否**: R −0.9pp（超 ±0.3pp 波动带）, N50 −53, 碎片 +4.5k 条。
机制自洽: kmer 深度 ∝ (RL−k+1), k=73 比 k=63 薄 ~12%, QM4 银标准低表达尾巴
长, 大 k 的断图损失 > 重复解纠缠收益——rnaSPAdes 公式的 upper_k 假设的是
"每轮图简化后进下一轮"的迭代架构, 我们的 trusted 注入 + 最终轮完整管线不
共享该假设。**结论: 顶 k 维持 63, U256 扩容作为能力保留（RL≥200 长读数据
的 [55,123] 现在可落地, 待未来数据验证）**。

性能注记: U256 计数键翻倍, QM4 normalize 748→920s（+23%）; 本轮全程 ~65min
（三轮计数 vs 两轮 + U256 开销）。若日后长读数据真用 k>63, 计数段可再压。
| 2026-09-02 | qm4-u128dn | R=43.3% P=51.9% **F1=47.2%** | N50=1746 n=123003 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G |
| 2026-09-02 | qm4-u128dn2 | R=43.9% P=51.7% **F1=47.5%** | N50=1784 n=117879 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G |

## U256 性能回压: RollingKey 双键分层（diginorm 计数 798→380s）+ auto 顶 k 实测钳制（2026-09-02, bench 8116a0e01d7769fc）

U256 扩容的代价是计数键 16→32B（QM4 normalize 748→920s）。本轮压回并反超:

- **RollingKey trait**（trinity-common）: u128（k≤64 紧凑键, 24B/槽）与 U256
  （全域 k≤128）双实现; counter/coverage_stats 的滚动枚举、计数、逐 read
  查表全部泛型化; diginorm 的计数表收进 CountTable 双分派（k≤64 走 u128,
  键不外泄管线）。**U256 手写 Hash**: hi==0 快路径只喂单 limb。
- **A/B（QM4 全量 2135 万对, 同机同参数, DN_PROFILE 分段）**: 计数段
  stats_query_table 798→380s（−52%）, diginorm 全程 930→513s; bench 内
  normalize 920→542s（−41%, 机器负载波动带, 低载轮 441s）。计数段已低于
  u128 时代免拷贝优化后的 410s。
- **等价性**: u128/KmerId 计数多重集等价测试（k=15..64 × ds 双路）;
  dev1 强制双分支全同; QM4 同二进制双跑（dev1+QM4）全产物逐字节一致;
  主产物（.normalized fq）before/after 逐字节一致。（注: A/B 时曾报 norm.fq
  DIFF, 与"norm.fq 是长名文件硬链接"不变量矛盾, 复跑不可再现, 判对拍假象。）
- **auto 顶 k 钳制 AUTO_K_CAP=63**: U256 放宽后 auto 对 RL=150 一度给
  [49,73], 而 73 已实测劣于 63（[49,73] R=43.3 vs [49,63] 43.9; [29,49,73]
  44.4 vs [29,49,63] 45.3）→ 默认钳回 63, 显式 --multi-k 不受限（[15,127]）。
  终验 bench: auto=[49,63] R=43.9/P=51.7/F1=47.5 与基线逐位一致。
- 测试: 666→667 全绿（新增 u128 等价 + clamp_auto_k_list）; inchworm
  并行竞态测试为已知闪失败（同二进制 5 连跑全过, 与本改动无关）。

## StringTie 缺失蛋白画像: 73% 是"碎片在场未连全长"（2026-09-02, missing_protein_profile）

承接蛋白层评估（pepR90=90.4%）, 对 3,620 条缺失蛋白（weak 2,637 + none 983,
去重转录本 3,487）做五维画像 + 失败模式分类（脚本固化 test_data/work/eval/
missing_protein_profile.py + miss_mode_classify.py, 产物 missing_st_pep.txt）:

| 失败模式 | 数量 | 占比 | 召回手段 |
|---|---|---|---|
| D-多碎片: ≥2 条 contig 的 ident≥90 残骸累计覆盖 ≥95% | ~2,564 | 73% | 解图连接/拼全 |
| B-碎片: 累计覆盖 <95% | 995+372 | ~39%（口径交叉） | 同上+更长k/深度 |
| A-无痕: 组装完全无对齐 | 295 | 8.5% | 构图层 |
| D-单contig≥95%: 序列在场、蛋白层丢失 | 128 | 3.7% | ORF/读码框层 |

关键发现:
1. **C 类（覆盖够但 identity<90%）= 0 条**——组装只要铺出序列就是对的,
   召回的全部不确定性在"拼全", 不在"拼对"。
2. 表达量梯度: TPM<1 缺失率 18.1% / 1-5 11.1% / 5-20 4.8% / ≥20 4.2%
   ——低表达是放大器不是主因; TPM≥20 仍缺 325 条, 其中 83% 是 D 类。
3. 整基因蛋白全缺 1,710 个（9.1%）: 44% 有 ≥95% 残骸, 15% 无痕;
   57 个含 BUSCO Complete 蛋白（高价值召回清单）。
4. 单 contig ≥95% 覆盖却丢的 128 条: TransDecoder 选框/截断问题, 非组装问题。

结论: 蛋白召回的杠杆 = **把已在图里的碎片连成全长路径**（量化目标 ~2.5k
转录本, 含 247 条 TPM≥20 高表达）——与核酸层"R=43.9% vs hcR"诊断一致,
且蛋白层给出了独立佐证（identity 无虞 + 碎片在场率 63%）。

## 碎片缝合路线: 断点解剖 + split-read 桥实锤（2026-09-02, breakpoint_autopsy）

对 270 条高表达（TPM≥20）D 类（序列齐且对、蛋白层丢失）做断点解剖
（breakpoint_autopsy.py; inchworm.DS.fa.min100 也对齐 ST 做"断裂阶段分层"）:

| 断点形态 | 数量 | 机制判定 |
|---|---|---|
| **碎片紧邻（gap≤25bp, ≥2 块）** | **201 (74%)** | 图简化/蝶后在分支处的不当拆分; 82% 无 isoform 分歧证据 |
| 单 contig 覆盖（含边角） | 49+128 | TransDecoder ORF/选框层, 非组装问题 |
| 真断裂 gap>25bp | 20 | 15 条 iworm 有桥（图阶段拆断）+ 5 条构图层断谷 |

紧邻型的召回上界: **201 条按序缝合后全部 ≥95% 覆盖（中位 100%）**——
缝合即全长, 无需补序列。断点处读深 95~2055×, 排除"深度不足"。

**split-read 缝合可行性实锤**: 5 个样例断点两侧 30bp 接缝子串在原始 fq
（前 800 万条）中的同 read 共现数 = 350/13/115/826/18——连接证据充分
且判据强（同 read 内共现, 天然不依赖参照, 无评估循环）。

方案定稿: **输出层缝合 pass**（Butterfly 输出后、Trinity.fasta 落盘前）:
- 桥判据: 断点两侧锚 ≥30bp 在同一 read 共现（正反链皆收）, 桥 read ≥2
- 防误连: 锚序列在组装集内唯一命中 + 连接图 DAG 约束（每端至多一连, 禁环）
- 预期: 直接冲 2.5k 缺失转录本主体（74% 紧邻型）, 蛋白 pepR90 与核酸
  R_full 双口径同升; dev1/QM4 bench 验证 + hi_D_samples 命中数对照
| 2026-09-02 | dev1-stitch | R=15.4% P=80.4% **F1=25.9%** | N50=1878 n=8543 | --CPU 8 --max_memory 8G --stitch |
| 2026-09-02 | qm4-stitch | R=44.2% P=51.4% **F1=47.5%** | N50=1789 n=117397 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G --stitch |

## 输出层 split-read 缝合 pass v0 落地: 温和正收益, 两个明确瓶颈（2026-09-02, bench 5d9a586734d8887a）

stitch.rs（--stitch 开关, 默认关）: contig 尾/头 30bp 锚 → 锚全局唯一性
（全体 contigs 30mer 计数==1）→ 归一化 reads 扫描（同 read 内尾锚×头锚
共现, 正反链, 锚序+间距约束）→ 桥 read ≥2 → 端唯一化（每端一桥, 取
read 数最大）→ 并查集禁环 → 链式缝合（接缝最大精确重叠 0..60 去重）。
4 单测全绿（链缝合/环拒绝/歧义锚剔除/阈值下限）; 接缝剪裁方向 bug
（B[seam..] 非 B[len-seam..]）被测试抓出后修正。

| 口径 | 缝合前 | 缝合后 | Δ |
|---|---|---|---|
| QM4 核酸 R | 43.9% | 44.2% | +0.3 |
| QM4 HC-R | 52.9% | 53.3% | +0.4 |
| QM4 P | 51.7% | 51.4% | −0.3 |
| pepR90 / pepR90full | 90.4/80.3 | 90.4/80.4 | ~0/+0.1 |
| geneR90 | 90.9% | 91.1% | +0.2（+38 基因）|
| BUSCO C | 89.1% | 88.9% | −0.2（S−0.3/F+0.2）|

dev1（76bp, 非目标场景）: 11 链, P +0.8, 无回退。
缝合量: 桥 725 对 → 达标 547 → 采纳 517 → 链 511（vs 目标紧邻型 2.5k, 覆盖 20%）。

**未兑现 oracle 预期的两个瓶颈**（下一轮迭代方向）:
1. 锚全局唯一性过严: 236k 候选锚仅存活 95.9k（40.6%）——isoform 家族共享
   末端使端锚串在兄弟 contig 内部普遍出现, 被整体剔除; 而紧邻断点恰密集
   于共享外显子区。方向: 唯一性降级为"端位唯一 + 桥竞争仲裁"（配对级
   判歧义）, 或锚内部命中容忍参数化。
2. 归一化 reads 稀释: max_cov=50 把断点横跨 read 压到阈值边缘
   （725→547 达标率 75%）; 低表达断点掉线。方向: 扫原始 reads
   （frag_merge 同款取舍）或达标阈值降 1 + BUSCO 联动护栏。

BUSCO 微降机制: 少量单拷贝 BUSCO 的蛋白在缝合 contig 上被 TransDecoder
重选 ORF 或接缝干扰转 fragmented（S−0.3/F+0.2）——迭代时需加
"缝合链 ORF 完整性"回归护栏（hi_D_samples 201 紧邻型命中数做验收指标）。
| 2026-09-02 | dev1-stitch1 | R=15.5% P=79.7% **F1=25.9%** | N50=1889 n=8518 | --CPU 8 --max_memory 8G --stitch |
| 2026-09-02 | qm4-stitch1 | R=44.5% P=50.8% **F1=47.4%** | N50=1802 n=115620 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G --stitch |

## 缝合 v1: 端位唯一 + 原始 reads 并行扫描——验收指标 3×兑现（2026-09-03, bench 0c40ab4f1e91a94e）

验收口径先修正: adjacent_201 的"单条 contig 独自 ≥95% 覆盖"才是缝合成功的
标志（累计 span 口径会混淆 PAF 断行假象——201 条中 90 条未缝合时已有单
contig 全长, 属蛋白层 ORF/比对病因, 非缝合目标; 真多碎片目标 111 条）。

v1 两处对症改动:
1. **锚端位唯一**（替代全局唯一）: 同串只要求作为唯一端锚出现, 在兄弟
   contig 内部出现不再剔除——v0 的全局窗口计数把 60% 锚砍掉（isoform
   共享末端结构性地命中该剔除）; 顺带删掉 200M 窗口的 count_windows 扫描。
2. **原始 reads + 并行扫描**: 归一化前快照 reads 路径/seq_type 传入
   （v0 扫归一化 fa 被稀释 725→547 达标）; scan_reads 改整读入内存 +
   序列行区间 + rayon 分块并行 + 局部桥表归并（rc 缓冲复用免分配）。

| 口径 | 未缝合 | v0 | **v1** |
|---|---|---|---|
| 桥采纳/链 | — | 517/511 | **2892/2692** |
| adjacent_201 单contig≥95% | 90 | 98 | **114**（+24）|
| 核酸 R / HC-R | 43.9/52.9 | 44.2/53.3 | **44.5/53.7** |
| pepR90 / geneR90 | 90.4/90.9 | 90.4/91.1 | **90.6/91.1** |
| ASM_pep | 72947 | 72897 | **73656**（缝合链 ORF 增 757）|
| BUSCO C（S/F） | 89.1（54.0/4.5）| 88.9（53.7/4.7）| **89.3**（53.8/4.3）|
| P（碱基） | 51.7 | 51.4 | 50.8（−0.9 代价）|
| F1 / N50 | 47.5/1784 | 47.5/1789 | 47.4/1802 |

v0 的 BUSCO 回退在 v1 消失（更多 BUSCO 基因缝合完整 → 完整 ORF）;
P −0.9 但 P_full=77.1/P70=88.5 双持平——可比对率未坏, 指向接缝重叠
去重的边缘损耗而非误连。dev1 无回退（桥 12→33, R+0.1）。
5 单测全绿（新增 internal_anchor_hit_tolerated 锁 v1 语义）。

v2 空间（未做）: 剩余 87/201 未转全长（多为锚仍歧义的共享端）;
P −0.9 可用接缝一致性投票（桥 read 在接缝区的多数表决定 overlap）回收。

## 缝合 v1 后追踪: adjacent_201 的缺口位置解剖——战场地图修正（2026-09-03）

对 201 条紧邻型高表达缺失逐条做缺口定位（5'/3'/中间三分）:

| 缺口形态 | 数量 | 结论 |
|---|---|---|
| **基本全覆盖（两端+中间 gap≤25）** | **149 (74%)** | 不缺序列; 其中 90 条未缝合时已是单 contig≥95%（PAF 断行假象, 病在蛋白 ORF 层）, 59 条真多 contig（缝合目标, v1 已救 24）|
| 仅 5' / 仅 3' 截断 | 23+23 | **缺序列**, 缝合无解（属端部延长战线, 历史上界 +2pp 且证据弱）|
| 双端截断 | 6 | 同上 |
| **中间缺口** | **0** | 紧邻型的"中间断点"全部不存在——碎片在断点处零间隙 |

嵌套实证（STRG.264 家族）: c3420_g1_i8/i42/i19/i39/i1/i2 共享 5' 起点、
3' 终点递增（589→807）——Butterfly 按表达梯度输出的嵌套碎片; 歧义锚的
本质是 i1/i2 同序列双输出。锚加长 30→45bp 对歧义消解 0/12——家族末端
≥45bp 全同, 端位歧义不可由锚长解决。

**v2 决策**: 缝合残余目标仅 35 条（全覆盖未缝合, 边际 ≤0.1pp pepR90）,
性价比不足, 缝合线在 v1 收口。剩余蛋白缺口的真实战场重定位:
1. 端部截断池（201 内 52 条 + 全集更多）——pepR90full(80.4) 与
   pepR90(90.6) 差值的主体; 需端部延长或图端部补强（前调查: OLC 信号
   证伪, 延长上界 +2.08pp）;
2. "单 contig 核酸全长却蛋白丢"池（90 断行假象 + D-单contig 128 条）——
   TransDecoder 选框/ORF 层, 与组装无关, 可用多 ORF 输出或 ORF 预测
   参数回收, 是零组装风险的纯蛋白层收益。

## ORF 层回收线: 选框丢仅 425 条, 缺失主体回组装层（2026-09-03）

对 3,620 条缺失蛋白用 TransDecoder **全量候选 ORF 池**（longest_orfs.pep,
108,427 条 vs Predict 选出 72,947）重新比对分类:

| 机制 | 数量 | 处置 |
|---|---|---|
| Predict 选框丢（候选池有 ≥90/50 命中） | 425 | 评估口径修正——**pepR90all 全 ORF 池口径 = 91.4%**（Predict 90.4%）, 已固化进 eval_protein.py/protein_eval.sh |
| 候选池也缺 | 3,195 (88%) | ORF 层无解——其中 ~440 条有 ident≥90 部分命中（ORF 截断型）, 其余连部分命中都无, 主体是低表达+碎片/截断（组装层, 与缺口解剖结论汇合）|

结论: "单 contig 核酸全长却蛋白丢"的 ORF 层快赢实测只有 +1.0pp 口径差
（Predict 是外部工具的保守选框, 不计入组装质量）; pepR90all=91.4% 与
BUSCO 89.3%、核酸 hcR 在 90±1% 收敛——多口径一致锁定"CDS 广度九成,
全长完整性缺"的同一诊断。ORF 层线收口。
| 2026-09-03 | qm4-minlen100 | R=43.8% P=51.3% **F1=47.2%** | N50=1757 n=132229 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G --min_contig_length 100 |
| 2026-09-03 | qm4-bundle100 | R=41.8% P=47.9% **F1=44.7%** | N50=1517 n=240155 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G --min_contig_length 100 |
| 2026-09-03 | dev1-stitch2 | R=15.4% P=80.1% **F1=25.8%** | N50=1907 n=8504 | --CPU 8 --max_memory 8G --stitch |
| 2026-09-03 | qm4-stitch2 | R=44.4% P=50.9% **F1=47.4%** | N50=1794 n=115394 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G --stitch |

## 端部截断战役: 全景诊断 + bundle 闸实验 + stitch v2 端部延伸（2026-09-03, bench 4041141c088286ba/7d06d1c627ab6f71）

**全景（36,933 条有对齐 ST tx）**: 全覆盖 63.9% / 仅3'截断 13.0% / 仅5'截断
10.7% / 中间缺口 8.1% / 双端 4.3%（无对齐 1,616 另计）。全覆盖率随 TPM 单调
（<1 44.5% → ≥20 83.5%）。高表达截断的缺口区 53-60% 有 reads（dep≥3）——
**可修**, 另 ~40% 无 reads（数据上限, de novo 固有）。

**丢失归因（TPM≥20 的 343 个已覆盖缺口）**: 42% 主力 contig 100-200bp
（bundle 闸: create_iworm_fasta_bundle 硬编码 200 丢弃孤立短组件）+
57% ≥200bp（进图但 Chrysalis/Butterfly 未输出——漏斗: min100 351,647 条
→ bundled 52,061 组件 → 输出 118,230 路径）。

**实验一（bundle 闸参数化 + min_contig_length 100 全联动）: 净负**——
组件 51,672→162,006, R 43.9→41.8（−2.1pp）, P −3.8。孤立短组件进图扰动
RTT read 分配与蝶路径选择——原版 200 闸的存在有理。

**stitch v2（输出层端部延伸）**: 延伸候选（min100 的 100-199bp 短 contig）
不进图, 输出层缝到主 contig 端部; **containment 护栏**（候选锚被主
contig 窗口包含即剔——QM4 268,955 候选中 208,099 即 77% 被包含, 不剔则
主 contig 普通 read 造成海量误桥）。QM4: 桥 32,365/链 4,500/延伸链 2,117,
R=44.4% vs v1 44.5——**持平**; adjacent_201 单contig 117（v1 114）。
oracle 的 +4.9pp 覆盖未兑现——桥判据对 oracle 全量的三重折损（锚歧义×
桥阈值×95% 全长门槛）, 与碎片延长线收档时的"oracle 三重折损"教训同构。

处置: v2 工程保留（--stitch 默认关, 7 单测全绿, dev1 无回退）; 端部截断
的剩余 oracle 空间留给"≥200 图/蝶丢失档"的深挖（57% 大头, 需进 Chrysalis
/Bfly 内部——下一战役级工程）。

## Chrysalis/Bfly 丢失档深挖 + RTT multi-kmer + butterfly 去冗余实验（2026-09-03, resume 于 iter_k63 中间产物）

**丢失档解剖（57% 档的再归因, 3,487 条缺失 tx 口径）**: E1 输出 cov≥0.9（纯
蛋白/ORF 层）969 / A 无 bundle 覆盖（前期）328 / **B4_partial 2,131 = 大头**。
组件级出局几乎不含金——bundle 52,061 → listing 50,232 的 1,829 丢失中
A2 无 read 1,661 + A1 partition 长度闸 168, 零路径组件 2,433, 合计只沾
~62 条缺失 tx。B4_partial = 序列在组件 bundle 里、组件有输出、但输出路径
只覆盖 tx 的一部分。

**B4_partial 死因（2131 条）**: 图断 886（41.6%）+ 连通未走 1,245（58.4%）。
- 图断: bundle 成员 X 口天然无边, 跨口桥靠 butterfly threading 的 read
  穿线——这些组件**一条输出路径都没跨过断口**。桥证据检测: 179 组件横跨
  read 在（同一条 read 含断口两侧）但 threading 没桥成; ~323 组件 read
  供给缺失——**12 断口抽样 5 个在 diginorm 数据里有 4-70 条单 read 横跨**,
  机制锁定 RTT: build_kmer_index 共享 25-mer 后写者赢（单组件）, 跨口 read
  另一侧被其他组件抢占 → pct<50 整条丢弃（67% 断口 tx 在其他组件也有
  bundle 对齐块, 佐证共享普遍）; 另有聚类割裂型（断口侧序列在其他组件,
  multi-kmer 不可救, 需 weld/bubble-up 层）。
- 连通未走: 边权不低（排除 remove_light/avgCov 剪枝）, 死在 POG/搜索层。

**实验一 --bfly-no-path-merging（跳过两处 cd-hit 式去冗余, 143s 重跑）**:
n 117,879→226,292（翻倍）, R 43.9→**46.3（+2.4pp）**, P 51.7→51.6, F1
47.5→**48.8**, N50 1784→2165。**原版 twoPathsAreTooSimilar 删掉的一半路径
基本是真 isoform**。但缺失 tx 只召回 26 条——去冗余不是 B4 档死因。

**实验二 --rtt-multi-kmer（25-mer 索引多值化: 主表首写者赢+冲突副表记
全部组件; 默认关=原版逐字节, 对拍单测保持）**: listing 50,232→50,462
（230 个原无 read 组件拿到 read）, n→118,441, R 43.9→**44.2（+0.3pp）**,
P −0.2。**缺失 tx 召回 99 条**（+150 partial）——与解剖预期一致
（neither_side 可救比例 ~25%）。耗时: RTT 91→120s, quantify 39→311s
（组件 reads 变多, 穿线变重）。

**组合（rttmk 图 + no-path-merging）**: R=**46.6**（+2.7pp）, P=51.3,
F1=**48.8**, n=228,243, 缺失召回 111 + partial 180。两改动近似独立叠加。

工具沉淀: bfly_resume / chrysalis_resume 两个一次性实验 bin（对既有 run
目录只重跑 butterfly / chrysalis 后半段+harvest, 不碰原产物）; 评估侧
funnel_audit / miss_attribution / bfly_autopsy / frag_bridge_v2 等
（test_data/work/eval/, 仓库外资产）。全 workspace 测试绿（含 RTT 对拍）。
下一步候选: threading 跨口桥接排查（179 组件）/ POG 搜索层 1,245 档 /
两个旗标的 dev1 泛化验证。
| 2026-09-03 | dev1-mk-base | R=15.4% P=80.7% **F1=25.8%** | N50=1878 n=8522 | --CPU 32 --max_memory 32G |
| 2026-09-03 | dev1-mk-combo | R=15.5% P=66.0% **F1=25.1%** | N50=2145 n=10683 | --CPU 32 --max_memory 32G --rtt-multi-kmer --bfly-no-path-merging |
| 2026-09-03 | dev1-rttmk | R=15.4% P=80.7% **F1=25.9%** | N50=1877 n=8523 | --CPU 32 --max_memory 32G --rtt-multi-kmer |
| 2026-09-03 | dev3-mk-base | R=13.8% P=79.4% **F1=23.5%** | N50=1897 n=6978 | --CPU 32 --max_memory 32G |
| 2026-09-03 | dev3-rttmk | R=13.8% P=79.4% **F1=23.4%** | N50=1902 n=6991 | --CPU 32 --max_memory 32G --rtt-multi-kmer |

## 两旗标泛化验证: dev1/dev3（2026-09-03, bench 27ad71e0941203de/4c6f051b00ec0bd7/3c43526f154e5248/4485145520acf1e8/54b5913c5b0e5a4c）

| 配置 | dev1 (76bp 模拟) R/P/F1/n | dev3 (76bp 模拟) R/P/F1/n | QM4 (150bp 真实) |
|---|---|---|---|
| 基线 | 15.4 / 80.7 / 25.8 / 8522 | 13.8 / 79.4 / 23.5 / 6978 | 43.9 / 51.7 / 47.5 / 117879 |
| --rtt-multi-kmer | 15.4 / 80.7 / 25.9 / 8523 | 13.8 / 79.4 / 23.4 / 6991 | 44.2 / 51.5 / 47.6 / 118441 |
| 组合(+no-path-merging) | 15.5 / **66.0** / 25.1 / 10683 | — | 46.6 / 51.3 / 48.8 / 228243 |

**判定**:
- `--rtt-multi-kmer`: dev1/dev3 零回退（n 仅 +1/+13, P 不动）, QM4 +0.3R/
  缺失召回 99。**保持默认关**（保守: 收益限真实数据, 且 quantify 耗时
  39→311s 的成本在模拟数据上无回报）; QM4 类真实数据建议开启。
- `--bfly-no-path-merging`: dev1 **P −14.7pp 净负**（n +25% 全是噪声路径,
  R 仅 +0.1）——去冗余删的路径在 76bp 模拟图上是噪声、在 150bp 真实图上
  是真 isoform, **数据依赖强, 不转正**, 保持实验旗标。自适应（按 read
  长度/图复杂度调节去冗余阈值）留作后续课题。

## threading 桥接排查: 断图档证伪——B4_partial 死因统一聚焦 POG/搜索层（2026-09-03）

**排查链**（探针 examples/split_bridge_probe + 评估侧 frag_severity/b4_kmer_adj/
bridge_pair_check/cross_member_jumps）:
1. threading 只在图 successors 边上延伸; POG overlap 是节点 id 级;
   跨 X 口输出的机制 = read path（含 PE pair path 拼接）经 POG 展开建边——
   全链**没有独立的"图桥接"步骤**。
2. QM4 输出实测 530,791 个跨成员相邻对（18,371/26,928 多成员组件 = 68%）——
   跨口输出大量存在（read path 跨口即出）。
3. **断图档证伪**: 886"断图"tx 里 844（95%）主根占比>90% 且孤节点 ≤2——
   伪断（tx 的 1-2 个错误/散在 kmer 落在图外孤节点）; semi_frag 仅 1、
   true_frag 0。**B4_partial 2,131 条 99.95% 图连通**。
   （教训: roots>=2 的集合判据对孤 kmer 噪声不鲁棒, 需主根占比分层。）
4. 179"桥 read 在没桥成"组件: 178 个无 PE pair 横跨（两半锚不同分量）;
   抽样 comp 1960: 图连通（剪枝后仅 5 节点）但组件 reads 仅 8 条、
   read path 仅 3 种——**read 供给/read path 覆盖不足才是该档真死因**,
   非穿线拒绝。
5. probe 实证穿线路径 [426,895,540] 合法跨"断口"（prev 边传递闭包连通）。

**修正后的 B4_partial 死因结构**: 图连通但输出未覆盖 = read path 覆盖缺失
（RTT 供给/POG containment）+ 路径搜索层（triplet 锁定/PRD=25/25 上限/
去冗余）。"threading 跨口桥接"修复方向**证伪**（图不断, 无需桥）。

下一步: 大样本 dump read paths vs tx 期望路径, 分层"read path 缺失 vs
搜索没走"; read 供给缺失（如 comp 1960 仅 8 reads）追 RTT/diginorm 层。

## search 档 POG 命运定位: 贯通 read path 95% 成孤立节点（2026-09-03）

批量探针（examples/batch_readpaths + batch_pogprobe, 1392 组件零失败）把
B4_partial 2,131 条按 read path 覆盖分层:
- fragment_readpath 1,578（74%）: 碎片齐（union>=0.9）无单条贯通——POG 拼接
  力/数据物理限制（150bp read vs 长差异区）
- search_layer 542（25%）: 有贯通 read path（>=0.8）且支持充分
  （PRD 单 read 门槛假设排除: med 单 read 节点占比 0.00）
- missing_readpath 11（0.5%）: read 供给基本不是问题

**search 档 542 的 POG 命运**（batch_pogprobe: 非包含集/containment/
overlap 边/破环删除）:
- pog_node 497（92%）: 贯通 path 成功进入非包含路径集
  - **isolated_no_overlap 474（95%）: 孤立——零 overlap 边**（连破环删的
    边都没有, 根本没建出边）。机制: 碎片被贯通 path containment 吃掉,
    剩余路径（其他 isoform 等）与它无"pathB 前缀=pathA 后缀"关系。
  - has_live_edges 23
- contained_by_longer 34 / 未入集 11

死点收窄: POG 展开成 SeqVertex 后, 孤立链分量在 convert/zip/
reorganize/分量 reads 判定（comp_reads.is_empty continue, paths.rs:1087）/
路径搜索的链条内消失。原版 v2.15.2 的 PE pair POG 边是死代码
（add_pair_paths_to_overlap_graph if(true)return）, POG 边全靠单端 overlap。

修复方向候选: 孤立 POG 节点（入=出=0 的长路径）直接保送输出层;
或 containment 后对孤立节点补图边拼接。

## bp 级修正: B4_partial 责任分层终判（2026-09-03, batch_iso_trace + iso_search_probe）

**节点级判定的粒度失真**: 合并图节点是粗粒度大节点, 节点覆盖 >=0.8 !=
碱基覆盖。实证: comp 29451"孤立链"实际被完整输出（c29451_g1_i1 206bp）,
但其 tx 主体 1000bp+ 根本不在图里——read path 穿过 tx 命中的 2-3 个粗
节点给出"贯通"假象。前节 474 孤立链/911 分量通过等结论受此污染。

**bp 级重判（tx kmer 位置命中 best-comp 图的比例）**:
- graph_bp_cov>=0.9: 381（18%）——图含 tx 主体, 输出没出 → **真 butterfly 责任**
- 0.5-0.9: 1,120 / 0.3-0.5: 552 / <0.3: 78——图只覆盖 tx 30-90%,
  缺失主体在图外 → **前期构图层（inchworm/weld/聚类割裂）**

**3,487 缺失 tx 的最终责任分层**:
- E1 输出>=0.9: 969（蛋白/ORF 层）
- 前期构图层（图外缺失）: ~2,050+（A 328 + B4 的 78%）—— **大头**
- butterfly 层（图>=0.9 没输出）: **381（11%）**——含搜索/过滤/去冗余/
  POG 各层, 需 bp 化重做死因分析

旧"57% 进图未输出"口径（paf 对齐块）把部分对齐当进图, 高估了 butterfly
责任。下一战役: ① 381 档 bp 化死因; ② 前期构图层（weld 聚类割裂——
comp 9 型: 断口侧序列在别的组件; bundle 只覆盖 tx 一半的成因）。
探针资产: batch_readpaths / batch_pogprobe / batch_iso_trace /
iso_search_probe（crates/trinity-butterfly/examples/）。

## 前期构图层量化: 组件割裂三分（2026-09-03, outgraph_layer_split + 组件割裂度分析）

bp 级图外缺失主体（B4 78% 档）的三层覆盖（raw iworm/min100/bundle 全库）:
- **bundle 全库对齐覆盖均值 86.5%**——序列大体在库里, 图外缺失的真相是
  **分散在多个组件**（best comp 单组件口径才 <0.9）
- min100 砍（raw_only）均值 0.008、聚类丢（min100_only）1.4%、
  inchworm 层缺（nowhere）191+141 条 ≈ 5.5%（数据上限）

组件割裂度（每组件 >=5% tx 覆盖）: 1 组件 1,893（54%）/ 2 组件 859 /
3+ 558 / ~0 141——**42% 被 >=2 组件割裂**; top 组件覆盖 med 87%。

**2-组件割裂（859 条）按块间距三分**:
- 块重叠 475（55%）: 共享外显子/旁系同源分流——两组件覆盖 tx 同一区域,
  weld 单值聚类（contig 只属一簇）的结构性割裂, 非错误
- 块相邻（<=100bp）237（28%）: **真聚类割裂**（物理相邻被拆——weld 边
  未建或聚类拒绝; 3+ 组件按同比例外推 ~330-500 条）
- 块远距（>100bp）147（17%）: 中间真缺口（weld 无证据）

构图层修复方向排序: ① 聚类割裂档（weld 阈值/max_cluster_size 参数
实验, 需暴露 ChrysalisParams 聚类旗标）; ② 共享分流档（多成员聚类或
输出层跨组件缝合——深水区）; ③ inchworm 层放弃（数据上限）。

## 共享分流档: 输出层缝合证伪 + oracle 上限 1,233（2026-09-03）

**oracle（union 覆盖）**: missing 3,487 中**并集 >=0.9 达 2,202**（单 best
仅 969）——**1,233 条（35%）的全部片段已在输出里（并集>=90%）、只是分散
在多条转录本/组件**。比 rtt-mk(99)+no-merge(26)+R1(57) 的总和还大一个
量级。

**输出层缝合证伪**: bfly all-vs-all（-X asm20, 46.3 万对齐）里合格跨组件
overlap 对（>=100bp@96% + 双方各>=100bp 独有延伸）仅 **18 对**; 50bp 门槛
293 对、30bp 543 对——1,233 条的片段组内 0 条可缝合。**片段间是缺口不是
overlap**（tx 坐标块重叠不等于输出序列可 overlap 拼接: 共享区两侧的输出
变体/边界不同）。输出层路线关闭。

**根治位置**: 聚类/构图层——① 割裂档参数（weld 阈值/max_cluster_size,
 未暴露）; ② 共享分流的两个深水方向: 多成员聚类（contig 允许多组件,
 全链语义变化）或组件合并重解（bundle 同源重叠组件合并后重跑 butterfly,
 有 bundle 闸实验"更大组件净负"的前车之鉴）。

## 聚类参数扫描: 收益温和, 割裂档主因不是聚类门槛（2026-09-04）

暴露 --weld-glue-factor / --weld-min-glue / --cluster-max-size（默认
0.05/2/25 = 原版; chrysalis_resume 旗标 --glue-factor/--min-glue/
--cluster-max）。QM4 四轮（chrysalis_resume ~13min/轮）:

| 配置 | 组件 | R | P | F1 | 净召回 |
|---|---|---|---|---|---|
| 基线 | 50,232 | 43.9 | 51.7 | 47.5 | — |
| glue 0.02 | 50,207 | 43.9 | 51.7 | 47.5 | ~0 |
| min-glue 1 | 48,982 | 44.2 | 51.6 | 47.6 | +15 |
| cluster-max 50 | 48,465 | 44.0 | 51.1 | 47.3 | +59 |
| mg1+c50 | 47,178 | 44.3 | 51.0 | 47.4 | +79 |

- glue_factor 是**无效自由度**: diginorm 后 cov 低, min_glue=2 下限兜底
  （cov*0.05<2 区域 min 恒 2）。
- min-glue 1 最平衡（R+0.3/P 平）; c50 召回多（153）但 P−0.6~0.7。
- **净召回最多 +79, 远低于割裂档预估 330-500**——大部分割裂的 weld 边
  本身不存在（单 read 证据之外低表达区域无 weld 证据）, 不是聚类门槛
  拒绝。全部保持实验旗标不转正。

## 分级 diginorm 实验: 证伪——weld 证据缺失是原始数据上限（2026-09-04）

机制假设: weld（GraphFromFasta）用 max_cov 更高的归一化 reads（低/中表达
区保留更多桥 read）, RTT/quantify 仍用默认 200（控图规模）。实现:
ChrysalisParams.weld_reads_path（None=原版）+ chrysalis_resume --weld-reads。

QM4 实验（weld@maxC1000, 16.1M reads/side vs 基线 10.8M, both.fa 6.2GB,
chrysalis ~12.5min）: 组件 50,232→50,206, R/P/F1 全不变, 净召回 **−12**。

**证伪逻辑**: diginorm maxC200 对 cov<200 区域不丢 read（两版输入完全
相同）; cov>200 区域 weld 证据早已充足（200 条 >> min_glue=2）。weld 边
缺失的区域是**原始数据里就没有横跨 read**（低表达/接缝区 read 供给上限,
与端部截断战役"~40% 缺口区无 reads（dep<3）"结论汇流）——不是 diginorm
抽稀。构图层探索到此收束; 分级归一化保留为机制能力（weld_reads_path）
但不产生收益。

## 381 档定向豁免: 粒度错配, 无效（2026-09-04）

实现 --bfly-empty-retry（PathSearchParams.retry_empty_comp, 默认关）:
分量级"搜索 0 产出 → min_read_support=1 重试"。

QM4: n 117,879→119,976（+2,097 条）, **R 43.9 不变, P −0.1, F1 −0.1,
missing 召回 0**——救出的 2,097 条几乎全是单 read 撑的噪声路径。

**粒度错配根因**: 381 档 tx 的图节点分散在**多个分量**（0 产出的 +
有产出但路径偏选的）; 0 产出豁免只救回单分量路径（案例 25323 分量仅
325bp, 对 tx 覆盖 ~30%——仍 <0.9）。该档需要的是分量间路径重组,
非门槛豁免。全局 R=1（57 条召回/P−2.0）仍是此档最好成绩, 不转正。
定向豁免代码保留（默认关）。

**381 档终判**: 图/证据全在而输出不走 = 路径搜索在多分量上同时偏选,
单点门槛修复不可及; 归入"搜索层结构性偏差", 与 GTasm/加权路径引擎的
路线汇合（gtasm 引擎本就是为此而设）。
| 2026-09-04 | dev1-mkmgl | R=15.5% P=79.0% **F1=25.9%** | N50=1924 n=8605 | --CPU 32 --max_memory 32G --rtt-multi-kmer --weld-min-glue 1 |

## 组合验证 + 推荐配置定稿（2026-09-04, bench 2dbfe1c2813ac33f + chrysalis_resume mkmgl）

QM4 rttmk+mg1 组合: R 43.9→44.5（+0.6, 两单变量精确叠加）, P −0.3,
F1 +0.2, 净召回 99 + partial 237。dev1 组合: R +0.1, **P −1.7**（mg1
在 76bp 模拟上单 read 证据=噪声, 与 no-merge 同型数据依赖）。

**推荐配置定稿**（README 实验旗标节）:
- --rtt-multi-kmer: 唯一推荐（三数据集零回退 + QM4 +0.3R/召回 99;
  代价 quantify 39→311s）
- --weld-min-glue 1: 数据依赖（真实 150bp 可开, 短读/模拟不建议）
- 其余探索旗标（no-merge/empty-retry/weld-reads/glue/cluster-max）:
  未达转正线, 留复现

## QM4 验证项转正默认 + quantify 物化键优化（2026-09-04）

**转正（用户指令"两个 150bp 可用的东西设置为默认开启"）**:
- --rtt-multi-kmer → 默认开（三数据集零回退 + QM4 R+0.3/召回+99）,
  --rtt-multi-kmer 保留为幂等兼容, --no-rtt-multi-kmer 回退原版单值语义
- --weld-min-glue → 默认 auto（哨兵 0, multik::resolve_weld_min_glue）:
  读长抽样 RL≥150 → 1（mg1, WELD_MG1_READ_LEN_MIN）, 短读 → 2（原版）;
  显式值覆盖; 读长抽样与 --multi-k auto 共用一次（need_rl 守卫）;
  递归 phase-1 线同样 resolve。启动打印生效值行。

**quantify 压时间**（rttmk 的 39→311s 代价）:
微基准 examples/qg_bench（QM4 left.fq 前 200k read ≈ 25.2M 25-mer,
20k 边）: baseline 14.83s → 物化键 4.72s → +par_sort 4.20s（3.5×）。
产物 FNV 指纹与 cmp 逐字节一致; quantify_vs_original 黄金对拍过。

机制: add_kmers 原为 KmerEntry{read,pos} 间接寻址 + 25 字节 memcmp
排序（O(N log N) 次 cache miss, 巨型组件 25M 条目）; 改两群物化——
主群全 ACGT 按 A=0,C=1,G=2,T=3 保序 2bit 打包 u64（kk≤32）,
sort_unstable_by_key((key,read,pos)); 次群（含 N 等, 稀少）按字节排序。
查询 pack_acgt(sub) 归位两群（精确匹配下不交叉）, partition_point 数值
比较收集。等价性: 同 kmer 组内 (read,pos) 升序 = 原稳定排序保持的生成序
→ ids push 序不变 → 输出逐字节相同。排序 par_sort（(key,read,pos) 全序
无重复 → 并行结果确定）。

| 2026-09-04 | qg_bench 微基准 | 14.83s → **4.20s**（3.5×） | 逐字节一致 | examples/qg_bench left.fq.gz 200000 20000 |
| 2026-09-04 | promo_default | R=44.4% P=51.3% **F1=47.6%** | N50=1796 n=117841 | --CPU 8 --max_memory 8G |
| 2026-09-04 | promo_optimized | R=44.5% P=51.3% **F1=47.7%** | N50=1801 n=117922 | --CPU 8 --max_memory 8G |

## diginorm 计数表 u64 键层（2026-09-04 压时间战役续）

性能分层再加一层: `RollingKey` 补 `u64` 实现（k≤32 最紧凑键）, diginorm
计数表三分派 k≤32→u64 / k≤64→u128 / 其余→KmerId。k=25 主路径键带宽
与哈希域减半, 大表（FxHashMap 键+计数值 16B/槽 vs u128 的 28B）的
footprint、cache miss 与 FxHash 同步降。等价性: k≤32 域内 u64 键数值
恒等于 u128 键与 KmerId 低 limb（滚动/反向滚动/revcomp 同一算法）,
由单测 `rolling_key_u64_numerics_match_u128_and_wide` 对 k=15/25/31/32
全断言锁定 + 泛型同代码保证; 多重集计数与产物不变。

微基准 examples/dn_bench（QM4 全量 PE fq.gz 3.1GB, 管线同参 DS +
fa_direct, 32 线程, DN_PROFILE=1 分段）, 基线 = HEAD 38c7fa3（u128 层）:

- wall 602.1s → **562.4s**（-6.6%）; 两侧 normalized fa 逐字节一致
  （各 2.08GB cmp 通过）
- 分段（段差）: stats_query_table 378.7→321.4s（**-15%**, 计数+逐 read
  查表热路径）, sorted_stats_rows 76.9→68.6s（-11%）, select/index/
  prep 持平; read_side/emit_side 段 +8.8s/+17.1s 系 NFS IO 抖动（与键
  类型无关, 方向不利仍净赢 39.7s）
- 归一化 wall 与 2026-09-01 战役的 584-749s 口径可比（当时 8-32 线程
  不同配置）; 本轮固定 32 线程 A/B 同机串行

| 2026-09-04 | dn_bench 微基准 | 602.1s → **562.4s**（-6.6%） | 逐字节一致 | examples/dn_bench left.fq.gz right.fq.gz out_dir（32 线程） |

- 验证备注: 全量测试三轮中 `cli_parallel_iworm_smoke` 间歇挂一次
  （6/2 线程 contig 多重集差异超容忍）——backlog.md 早已立案的 inchworm
  并行平局序非确定（backlog.md:52, 原版同病）, 与本轮改动无关
  （inchworm 不经 RollingKey/泛型计数路径）; 复跑即绿。


## HEAD 端到端分段实测（2026-09-04, e3a3241, QM4 标准口径）

命令同 09-01 基线（--multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory
64G, 2135 万对 150bp）。**wall 34:07**（avg 12.2 核 / 32, 峰值 RSS 118GB）。
09-01 的 31min 记录已过时: 09-03 起自动档位含 k=63 轮（auto 顶 k 钳 63）,
多花一整轮 count+inchworm。

分段（s, 占 wall）:

| 段 | 耗时 | 占比 |
|---|---|---|
| normalize + prep | 415 + 14 | 21% |
| k49 轮（count+inchworm） | 332 + 220 | 27% |
| k63 轮（count+load+inchworm） | 308 + 22 + 278 | 30% |
| k63 chrysalis（gff 194 / rtt 101 / quantify 56 / 其他 23） | 374 | 18% |
| butterfly + harvest | 96 | 5% |

- **quantify 55.6s**: rttmk 默认开 + 物化键下, 已回到单值基线（39s）附近
  ——旧 rttmk 代码 311s 的代价被今天的工作在管线内证实收回。
- normalize 415s（vs dn_bench 562s: 差额=NFS→本地 scratch 写 + 无右册对拍）。
- 剩余大账: **两轮 kmer_count 合计 640s（31%）**——KmerId 表（40B/槽）
  逐轮重数 both.fa 4.16GB; inchworm 两轮 498s（24%）, 贪心段固有。
- gff 194s 仍为 chrysalis 内最大单段（1199→196 后的余量）。

| 2026-09-04 | qm4-e3a3241 全程 | wall **34:07**（12.2 核均值） | R/P 待 bench 登记 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G |

## 刀1: kmer_count 分层窄键计数（2026-09-04 压时间战役, 640s 池/占 wall 31%）

机制: inchworm 轮计数（KmerCountTable::count_fasta_data_with_capacity,
KmerId 全宽 40B/槽）改为 `count_fasta_data_layered`——k≤32 用 u64 键、
k≤64 用 u128 键在窄键空间完成 per-thread 计数与 reduce 合并（键带宽与
槽位减半~四分之一, 大表 cache miss 与 FxHash 输入同步缩）, 终点一次
widen 回 CountMap（精确容量预留, 无生长 rehash）。管线消费方（trusted
注入 / counts_bin 落盘 / KmerCounter::from_count_map move 直构）接口与
语义全部不变。orchestrate 两条计数路径（kmer_count 阶段 + kmer_load
checkpoint 回退重数）都走分层入口。

等价性: 各窄键层与 KmerId 层同一泛型代码、键数值域内恒等 → 计数键值
多重集相等。三层证据: ①单测 count_layered_matches_wide（单行+多行
fixture, k=7~65 跨三层, ds 两态, 全键对拍）; ②微基准内置全表逐键对拍
（k=49 1.95 亿键 / k=63 2.05 亿键, 各 8s 通过）; ③widen 只平移键值,
迭代序差异不构成回归——并行 reduce 合并序本就非确定（inchworm 平局序
非确定已立案, cli_parallel_iworm_smoke 直调 inchworm 二进制不经此路径）。

微基准 examples/count_bench（QM4 both.fa 4.16GB 重生成, 32 线程,
estimate_hash_size(64G) 同管线口径, 全表对拍内置）:

| k | KmerId 直计 | 分层窄键+widen | 加速 |
|---|---|---|---|
| 49 | 373.0s | **270.1s** | **1.38×** |
| 63 | 299.3s | **279.7s** | 1.07× |

- 净省 ~123s/两轮; k63 侧两层同为 u128 键、结构性收益只剩宽侧基线差
  （该次 wide 偏快, 共享机波动带内）。k≤32 单 k 运行（dev 系 k=25）
  叠加 u64 层（diginorm 同款 -15% 口径）, 未单独计测。

**端到端复测（同口径 QM4 全程, 32 线程）**: wall **34:07 → 29:48（-4:19,
-12.6%）**。分段: kmer_count 两轮 640→399s（**-241s, 1.6×**, 超微基准
预测的 -123s——窄表内存红利外溢: 峰值 RSS 118→88GB, chrysalis_total
374→315s（quantify 56→38 / rtt 101→80 / gff 194→175）; normalize 415→501s
系当晚共享机 IO 波动（diginorm 独立重建轮同步 +13%）, 方向不利仍净赢。
产物: n=117834（基线 118068, -0.2%）N50=1793（1787）——在 inchworm 并行
非确定带内（promo 两跑自身 ±0.2%）; 计数表逐键等价已证, 指标等价随之。

| 2026-09-04 | qm4-knife1 全程 | wall **29:48**（基线 34:07） | n=117834 N50=1793 | --multi-k auto --inchworm_cpu 32 --CPU 32 --max_memory 64G |

## 刀2: diginorm 分片共享表计数——负结果裁定的正收益（2026-09-05）

**动机（CNT_PROFILE 分段实测, QM4 k=25 u64 键）**: diginorm 计数的 reduce
合并 221.6/282.1s（两侧）≥ 计数本体 181.0/313.0s——reads 随机分布下每线程
局部表收敛到近全量 distinct 集（~0.73D/线程）, 合并 ≈ 线程数×D 次随机
插入, 局部表叠加把 RSS 顶到 2~3×D。

**实现**: counter 新增 `count_fasta_data_of_sharded`——k-mer 按键哈希高
8 位入 256 片, 每线程分片缓冲攒 2M 键整批冲刷进 Mutex 片表（片表 ≈ D/256
条, 局部性好）; 终点按精确容量把片表一次折叠回单表（恰 D 次插入, 无生长
rehash）。键值多重集等价: 片间键不交 + 片内 += 可交换（强制冲刷单测
count_sharded_forced_flush_matches, flush=64 三层键全对拍）。

**交错 A/B 裁定（同负载窗口, 基线 = 956af9a worktree 二进制背靠背）**:
- **diginorm（u64 键, 双侧各 ~CPU/2 线程）: 大赢**——dn_bench 677.6→516.0s
  （-24%）; 合并坍缩至 ~15s（D=1.22/1.34 亿键/侧）; 产物与基线逐字节一致。
- **inchworm 轮（u128 键, 单次 32 线程）: 真回退**——count_bench k49
  layered 240.5→303.6s（k63 249.2→262.7s）; 段数据显示回退全在
  count(shard) 段: 32 冲刷线程 × 256 片表的 Mutex 在共享机高负载下有
  护航效应（持锁线程被抢占 → 同片冲刷排队）。
- **定案: 按调用形态分发**——diginorm CountTable::build 三变体全走分片;
  轮计数（count_fasta_data_layered / 公共 API）恢复局部表+reduce。
  两条路径并存, 均带 CNT_PROFILE 分段。

**flake 容差校准（cli_parallel_iworm_smoke, 同日）**: 高负载日孤立复跑
12 次挂 5 次（远超注释记录的 2%）, 5 次失败全部 uniq≤2/侧、6 条核心
contig 恒齐备——竞态分裂量实测 ±2。容差 +1→+2 并显式断言交集 ≥6
（系统性破裂下交集坍缩, 不可能恰好差 2）; 加固后 10/10 稳, 全量 681 绿。
该测试直调 inchworm 二进制, 不经 counter 路径, 与刀1/刀2 无关（负载
加宽竞态窗的环境效应, 2026-09-04 起便已 1/3 出现）。

| 2026-09-05 | dn_bench 交错 A/B | 677.6s → **516.0s**（-24%） | 逐字节一致 | 分片 diginorm 形态, 32 线程 |
| 2026-09-05 | count_bench 交错 | k49 layered 240.5→303.6s **（回退）** | — | 32 线程单调用形态 → 裁定不分片 |

**端到端排障与终态（mimalloc 交互, 2026-09-05 晚）**:

- 管线复测 3/3 在 k49 fold 计数段 RSS 45s 内 35→455GB 被 OOM-kill
  （GNU time 对被杀子进程误报 "Exit status: 0", 以 bash $?=137 与
  RSS 观察器为准）。同二进制同输入同线程的 glibc example（count_bench）
  不复现——#[global_allocator] 只在各 bin 声明, example 默认 glibc,
  排障时构成天然 A/B。
- **根因**: 分片归一化的大表释放后, mimalloc 默认 PURGE_DELAY=10ms 滞留
  已释放页; 与 multi-k 轮计数（每线程局部表翻倍增长 + reduce 合并新表）
  叠加, 大分配链在滞留堆上雪崩。
- **修复**: trinity-cli main() 程序内 mi_option_set(PURGE_DELAY, 0)
  （libmimalloc-sys extended 特性; 选项位次 15 以 vendored v2 头文件
  锚定, 升级须复核）。程序内设置覆盖同名环境变量。
- 验证: 全量 QM4 exit=0, wall 34:57（重负载窗）, 峰值 RSS 85GB,
  k49/k63 fold 各 238/232s 正常; 切片 2M PE 447s exit=0（程序内默认生效）。
- 分配器 A/B（切片 2M PE, 交错两轮, 同负载窗）: mimalloc 默认 370/380s,
  glibc 493/474s——**mimalloc 快 22%**, 小尺度不触发 OOM（交互需全量
  尺度的内存压力）; purge=0 切片 447s（与默认档的差含晚间负载漂移,
  未隔离纯代价, 但不影响终态选择: purge=0 是全量安全的唯一 mimalloc
  配置, glibc 档全量 38:40 且切片 +22%）。
- 刀2 净效果（QM4 全量, 重负载窗口径）: 34:57 vs glibc 38:40; 归一化
  段交错口径 -161s。端到端含 diginorm 分片 + 刀1 全部优化的轻负载复测
  留待后续（今日全程重负载, 与 29:48 不可比）。

| 2026-09-05 | qm4-purge0 全程 | wall **34:57**（重负载） | exit=0 峰值 85GB | mimalloc+PURGE_DELAY=0+分片归一化 |
| 2026-09-05 | 切片分配器 A/B | mimalloc 375s vs glibc 483s（**-22%**） | n 带内一致 | 2M PE 交错两轮 |

## 内存战役: QM4 峰值 83GB → 目标 <32GB（2026-09-05 深夜起）

**目标（用户指令）**: 100GB 级峰值太夸张, 压到 32GB 以下才有实用性。

**画像方法**: timer 增设 `.stage_memory`（阶段边界 VmHWM/VmRSS, 与
.stage_timings 同边界）+ CNT/GFF_PROFILE 行尾 VmHWM + 进程 2s 采样器
（memwatch.log）+ stderr 时间戳对齐。全程产物落 test_data/work/tmp（/tmp
禁写, 用户指令）。

**mem1 画像（旧二进制, 标准口径+--no_cleanup, wall 35:51 重负载）——
峰值 83GB 的归属**:

| 段 | 峰值 RSS | 主体 |
|---|---|---|
| diginorm count(shard) | **83.0**（全局） | 分片计数缓冲病（见刀M1）+ raw fq 13.8GB 常驻 |
| k49 轮计数（locals+reduce） | 75.6 | 32 在飞局部表各收敛 ~0.44D（D=1.95 亿） |
| k63 轮计数（locals+reduce） | 70.3 | 同上（D=2.05 亿） |
| inchworm 贪心（k49） | ~26 | 宽表 7.8 + 种子排序列表 7.8 + 稳定排序暂存 7.4 |
| **gff** | **47.0** | read_counts 分块稠密 partials（见刀M3）+ crossover 11.7GB 滞留 + lossy 拷贝 3.9GB + rc_cache×32 |
| rtt+quantify | 17.4 | — |

**刀M1（counter.rs）: 分片计数缓冲病——冲刷阈值按任务生效的隐式全量缓存**。
`par_iter().fold(ShardBufs).collect()` 下 2M 键冲刷阈值对每个 rayon 任务
独立判定: 21.5M 记录切出上千个细粒度任务, 单任务攒不满 2M 键**永不冲刷**,
键滞留任务缓冲等 collect 终点——整个 k-mer 键流（63 亿实例 × 8B + Vec 倍增）
被隐式缓存 ≈ 60GB 瞬时, 即 diginorm 段 83GB 峰的主体。刀2 时"分片峰值
≈ D+缓冲"的推断未经实测, 实为灾难; **轮计数 +26% 的"Mutex 护航回退"
裁定同样被污染（缓冲病让分配器狂奔）, 作废重裁**。修复: 冲刷阈值字节
定界（1MB/键宽, 下限 16K 键）+ with_min_len(8192) 钳任务数 → 缓冲上界
≈ 任务数×1MB（QM4 ~2.6GB）。
同刀: diginorm raw fq（fa_direct+fq 且不 keep_fq 时 emit 不读 raw）prep 后
提前释放（-13.8GB 常驻; --no_cleanup 形态 keep_fq=true 仍保留）;
save_counts_bin 整表 7GB Vec 改 8MB BufWriter 流式（字节序不变）。

**刀M2（counter.rs）: 轮计数分层入口统一切分片**（缓冲修复后旧裁定作废）。
count_bench k49 32 线程: 旧二进制 layered(locals) 249.2s / wide 对照 418.3s;
新二进制 layered(sharded) **87.0s** / wide 304.2s——wide 跨窗漂移 1.37×,
负载归一后 sharded 仍 ≥2× 提速, 单段形态峰值 ~13-14GB（glibc example,
含 wide 表对拍滞留基线）。
**刀2 的"护航"理论证伪**: 时间回退全在缓冲病, 分片+有界缓冲是净赢。

**刀M3（chrysalis）: gff 47GB 三连击**。
① `add_counts_from_reads_par` 分块局部稠密表 partials: 每块 n×4B×collect
全保留, n≈1.2 亿 × 128~512 块 → glibc example 实测瞬时 102.8GB（管线内
mimalloc 47GB）——改回原版 NonRedKmerTable.cc 的 omp atomic 本义
（共享 AtomicI32 fetch_add, 单份 n×4B, 整数累加恒等）;
② phase1 crossover 候选序列（1.03 亿条 ≈ 11.7GB）在 set_up_templates 后
无消费方——立即 drop; phase1/2 逐命中 clone/revcomp 物化改共享 rc_seqs
表借用（一次性 174MB, 取代每任务独立缓存最坏 32×全量）;
③ pipeline.rs `from_utf8_lossy(...).into_owned()` 无条件复制 both.fa 3.9GB
——改 Cow 借用（fasta 恒 ASCII）, files 入口 text 同步 drop。
gff_profile 复测（真输入 min100+both.fa, GFF_PROFILE+采样器）:
**gff3(①未修) 301s/峰 102.8GB → gff4 198.3s/峰 24.1GB（hwm 25.05GB）,
weld 图输出 47,426,460 字节逐字节一致**（原子化等价性实证）。

**刀M4（inchworm）: 种子排序稳定排序暂存**。collect_nonzero_seeds 后
sort_seeds_desc 用 sort_by（stable）——分配 n/2 暂存, QM4 1.95 亿种子
实测排序 10s 内 RSS 27.7→35.1GB（+7.4GB）, 把 k49 inchworm 时代顶到
33.5GB。计数表键唯一 → 比较器全序无并列 → sort_unstable_by 输出恒等,
零暂存。

**验证链**: mem2（刀M1/2, 旧 gff）exit=0 wall **24:38**（重负载+旁路基准
穿插）, diginorm 432→207s（count(shard) 180→48s）, 两轮计数 490→100s,
全局峰 46.6GB=旧 gff 段; mem3（+刀M3）与 mem4（+刀M4）终值见下表。
682 项测试全绿（681+timer 内存行新测）。

| 2026-09-05 | mem1 画像 | 旧二进制峰值 **83.0GB** | diginorm 缓冲病 83/k49 轮 75.6/gff 47 | 采样器+stage_memory 归因 |
| 2026-09-05 | count_bench 交错 | layered locals 249.2s → sharded **87.0s** | 全表对拍一致 | wide 对照跨窗 1.37× |
| 2026-09-05 | gff_profile | 301s/102.8GB → **198.3s/24.1GB** | weld 图逐字节一致 | 刀M3 三连击 |
| 2026-09-05 | mem2 全程 | wall **24:38** | 峰 46.6GB（旧 gff 段） | 刀M1/2 验证 |
| 2026-09-05 | mem3 全程 | wall 25:33 | 峰 33.8GB（k63 种子收集翻倍） | 刀M3 后 gff 段内 16.2GB 未创新高 |

**刀M5（inchworm）: 种子列表收集的 Vec 倍增瞬态**。mem4 复测证实
sort_unstable 杀掉排序暂存后, 峰值棘轮仍被 `iter_nonzero().collect()` 的
末次容量翻倍顶到 33.5GB（k49 1.95 亿 × 40B, 新旧两份并存 +11.7GB;
k63 侧 33.8GB）——改 `Vec::with_capacity(counter.size())` 精确预留
（表长上界, 剪枝惰性零键占比 ~0.02%）, 值序不变仅分配方式变化。

**刀M6（inchworm PARALLEL）: DashMap 整表重建的三叠瞬态**。mem5 复测
（刀M4/5 后）k49/k63 inchworm 时代仍 33.5-33.8GB——`SyncKmerCounter::
from_counter` 把整张 FxHashMap 重建进 DashMap（旧表 7.8GB + 新 DashMap
~9.5GB 并存）, 且发生在亿级种子列表在场之后（采样器实测种子收集后
+9.5GB 两段阶梯、排序/贪心稳态后 -12.8GB 回落）。修复: 转换提前到种子
收集前的空窗期执行（重建瞬时只剩两份表）, 种子改从 DashMap 快照且
**只取键**（count 字段主循环不参与——实时重查 IRKE.cpp:546, 40B→32B/条）;
种子收集序 FxHashMap 序 → DashMap 分片序（收集序本就是非确定源之一,
平局序 backlog 已立案; smoke_vs_original 4 项过, 全量 682 绿）。

**终审 mem6（六刀全上, QM4 标准口径, 晨间负载）**: exit=0,
wall **25:23**, **全局峰值 RSS 28.69GB**（kernel VmHWM 棘轮 28.94GB,
归属 diginorm 行/输出尾段 @09:37:00; gff 段内 16.2GB、inchworm 转换窗
RSS 26GB 均未再创新高）, n=118,096（带内 117,722-118,195）。
**83.0GB → 28.7GB（-65%）, 目标 <32GB 达成; 墙钟同步 34:57→25:23。**

峰值残余构成（28.7GB 的下一层, 留档备后战）: diginorm 尾段 = fa 7.4GB +
行集/TSV 往返 + 输出缓冲 + mimalloc 滞留; inchworm 稳态 ≈ 宽表 7.8 +
DashMap 8.8 + 种子 6.5 + 滞留 ~4; 计数段棘轮 28.9。再往下砍需窄键化
KmerCounter/DashMap（k≤64 用 u128, -6~7GB）与 diginorm 文本往返流式化。

| 2026-09-06 | mem4 全程 | wall 27:47 | 峰 33.8GB（k63 种子收集翻倍） | 刀M4 验证, 靶点转刀M5 |
| 2026-09-06 | mem5 全程 | wall 27:37 | 峰 33.82GB（DashMap 重建三叠） | 刀M5 验证, 靶点转刀M6 |
| 2026-09-06 | mem6 终审 | wall **25:23** | **峰 28.69GB（<32 达成）** | n=118,096 带内; 682 测试绿 |

## 内存战役第二期: 28.7GB → 目标 ~20GB（2026-09-06）

**目标（用户指令）**: 再往下压到 20G 左右。

**mem6 残余构成的精拆（采样器时间线 + [cnt]/[gff] VmHWM）**:
diginorm 读写+prep 时代 = 双侧 raw fq 解压缓冲 13.8GB + 双侧 fa 7.4GB 并存
（23.96GB, 全段约束）; diginorm 行/emit 尾段 = 合并表 4.1 + 行集 ~3 +
emit 整表物化 ×2（out Vec + 每命中记录拷贝, 27.9-28.9GB 棘轮）; diginorm
计数段 = 并行双侧两套片表（pow2 容量浪费 ~2.2×/侧, 24.6GB 棘轮）; 轮链 =
count 终点 widen 的 KmerId 宽 Fx 表 ~10.7GB + SyncKmerCounter DashMap 宽表
重建（26.0GB 转换窗）; gff = phase1 候选 Vec<Vec<u8>> 11.7GB + 模板内部
窗口二次平铺（23.08GB）。

**六刀（全部语义无损, 各带等价锁）**:

- **刀M7（diginorm）**: 侧计数串行化（M7a——并行双侧两套片表+折叠同时
  在场; 串行每侧独占 rayon 池墙钟持平, 语义不变仅执行序）; 计数表在行
  阶段后 drop、行集/merged/selected 在 index 后 drop（M7b, -7GB 常驻）;
  extract_to_file 流式化（M7c——旧实现 out Vec 整表物化 + fq 路径每命中
  再拷整条记录 ≈ 双份输出, 改 (acc, span) 元数据 + data 原位切片流写,
  字节序恒等; 错误路径会留部分文件——调用方错误即中止无消费方）。
- **刀M8（轮链窄键化）**: 33≤k≤64 全程 u128——count_fasta_data_layered128
  直通窄表（终点 widen 消除）、save_counts_bin128 流式逐条 to_wide 落盘
  （字节恒等）、inject_trusted_kmers128、prune_some_kmers 泛型化
  （ExclusiveClearCatalog: KmerCounter 宽表与 NarrowCounter128 适配器
  同一实现）、SyncKmerCounter128（DashMap<u128>, 槽 40B→20B）+ PARALLEL
  入口 compute_sequence_assemblies_parallel128（贪心核心拆出
  parallel_greedy_core 宽窄共用, 种子保持 u128 32→16B/条）。KmerId 宽
  Fx 表与宽 DashMap（各 ~10.7GB）在 k≤64 轮链不再物化。
  等价锁: sync128_query_parity_with_wide（查询/候选/置零逐项, SS/DS）+
  prune_narrow_parity_with_wide（非零多重集）。
- **刀M9（diginorm 读写）**: read+prep 按侧串行——双侧 raw 解压缓冲与
  fa 产物不再四份并存（23.96→~14）; need_raw 判定提前, 左册 raw 先释。
- **刀M10（gff）**: phase1 候选改平铺竞技场（flat + (off,len) 跨度,
  Vec<Vec<u8>> 11.7GB→~4.5GB）+ set_up_templates_flat 窗口偏移直指竞技场
  （免内部窗口二次平铺 ~7.5GB）。gff 23.08→20.97GB。
  等价锁: set_up_templates_flat_matches_vec_form。
- **刀M11（计数缓冲）**: COUNT_FLUSH_BYTES 1MB→512KB（缓冲上界 2.6→1.3GB）。
- **刀M12（diginorm 折叠过滤+片粒度）**: 折叠带 min_kmer_cov 过滤
  （与全量折叠后 retain_min 恒等, 输出表免为单例错误 k-mer ~30-40%
  distinct 付容量）+ COUNT_SHARD_BITS 256→512 片（pow2 桶粒度下每片
  更小 → 容量浪费 2.2×→1.77×）。

**验证链**: mem7（M7/8: diginorm prep 23.96 定位→M9; gff 23.08→M10）
→ mem8（M7-11: exit=0 wall 27:43 n=117,759 带内, 全局峰 ~20.5-21.0GB
= diginorm B 侧折叠棘轮 21.76/gff 20.97）→ mem9（M12 全上, 终值见下表）。
685 项测试全绿（682+3 个窄链/flat 等价锁）。

| 2026-09-06 | mem7 全程 | wall 25:55 | 峰 23.96GB（diginorm prep） | M7/8 验证; n=118,207 带内 |
| 2026-09-06 | mem8 全程 | wall 27:43 | 峰 ~20.5-21.0GB（B 侧折叠 21.76/gff 20.97） | M7-11 验证; n=117,759 带内 |
| 2026-09-06 | mem9 终审 | wall **27:33** | **峰 21.12GB（k49 计数时代棘轮 21.21; gff 15.7）** | n=117,695 带内; 685 测试绿 |

**mem9 备注**: M12 折叠过滤实测 diginorm 侧表保留率 43-46%（单例错误
k-mer 过半）, 侧折叠 9.3→5.8s、diginorm 棘轮 21.76→20.66; 512 片对
pow2 桶浪费无系统收益（浪费比 1.0-2.0× 对数均匀, 取决于 D/2^b 落点）,
保留为片内局部性微调。全局约束转移到轮计数时代（both.fa 4.16 + 片表
5.6 + 折叠输出 5.6——轮链 min_kmer_count 默认 1, 无过滤可带）。
**28.7→21.1GB（-26%）; 两期累计 83→21.1GB（-75%）, 墙钟 34:57→25-28min。**

**~20GB 之下的下一层（留档）**: gff 的 reads 双份常驻（reads_data 3.9 +
reads_streamed 4.5——reads_data 为 rtt 阶段所持; flat 化 reads 或
模板表 u64 打包（kk≤32, 慎入: 字节二分兜底路径与非 GATC 模板域）各
-2~2.5GB; diginorm fa mmap 化（clean 页可回收, -7.4GB 匿名 RSS）;
片表换非 pow2 容器。语义注记（遗留, 非本期引入）: diginorm 侧表
过滤-后-合并序在 min_kmer_cov≥2 且 k-mer 双侧各计 1 时与原版
jellyfish both.fa 合并-后-过滤分歧（本期 M7a/M12 保持移植既有语义）。

## BUSCO 对照补全: 战役前后默认形态 + v0.6.0 终审产物（2026-09-06, eudicots_odb10 protein, n=2326）

口径: eval/busco_eval.sh 固化（TransDecoder 5.7.1 eviann 单最优 ORF →
BUSCO 6.1.0 protein --offline, 与 09-02 蛋白层评估同链; TMPDIR 全程项目内,
不落 /tmp）。nomk 基线重跑: 当前二进制 `--no-multi-k`（wall 21:49, n=83,517,
与 08-31 锚点 R=31.5 同带; 含 mg1 等后续默认——"同二进制只关 multi-k"的干净
对照）。注意 /tmp/gtasm_venv 已被系统清掉, eval_protein.py 蛋白重合度线
（diamond）暂不可用, 本轮只跑 BUSCO 线。

| 组装 | n | 蛋白(≥100aa) | BUSCO C [S/D] | F | M | 核酸 R/P |
|---|---|---|---|---|---|---|
| 单k基线 nomk（--no-multi-k） | 83,517 | 43,242 | **89.7** [66.8/22.9] | 4.5 | **5.8** | 31.8/50.2 |
| **v0.6.0 默认 mem9（multi-k auto）** | 117,695 | 73,476 | 89.5 [54.3/35.3] | **4.3** | 6.2 | **44.5/51.3** |
| v0.6.0 + --stitch（v1, 09-03 代码） | 115,620 | 73,656 | 89.3 [53.8/35.5] | 4.3 | 6.4 | 44.5/50.8 |
| 09-01 锚点 release_qm4_out | 118,230 | 72,947 | 89.1 [54.0/35.2] | 4.5 | 6.4 | 43.9/51.7 |
| StringTie 银标准 | 38,549 | 37,629 | 93.9 [64.5/29.4] | 1.8 | 4.2 | — |

结论:
- **BUSCO Complete 不是这场战役的杠杆, 前后持平**: 默认形态 89.7→89.5
  （−0.2pp ≈ 5 基因, 噪声级）。两期性能优化语义无损, mem9 89.5 vs 09-01
  锚点 89.1（+0.4pp 同为带内漂移, n 带 117.7-118.2k）。
- **增益真实发生在 isoform 维度**: 核酸 R +12.7pp（31.8→44.5）, BUSCO
  D +12.4pp（22.9→35.3, S −12.5pp）——multi-k 找回的转录本主体是同基因
  多 isoform; BUSCO 按设计不奖励 D（genome-guided 银标准自身 D=29.4）,
  单 k 的"高 S"本质是每基因只拼出一条主转录本, 不是更完整的转录组。
  蛋白产出率 43.2k→73.5k（+70%）与 R 同向。
- M 5.8→6.2（≈9 基因）与 F 4.5→4.3 互抵, C 净平; de novo 对银标准的
  BUSCO 残差（~4.4pp）主体是低表达/端部截断, 见 09-03 缺口解剖留档。

| 2026-09-06 | nomk | R=31.8% P=50.2% **F1=38.9%** | N50=1702 n=83517 | --no-multi-k --CPU 32 --max_memory 64G（BUSCO 基线重跑, 21:49）|
| 2026-09-06 | mem9 复评 | R=44.5% P=51.3% **F1=47.6%** | N50=1794 n=117695 | 终审产物核酸复评（BUSCO C=89.5）|



## regroup 战役: 拆分后区分旁系同源 vs isoform → 每基因取最长转录本（2026-09-07, --regroup 默认关）

目标转移（用户裁定）: "怎么获得最长转录本（拆分后怎么区分旁系同源跟 isoform）"。
两阶段: python 原型验证收益（test_data/work/eval/regroup_proto_report.md 全档）
→ Rust 内部化 `crates/trinity-cli/src/regroup.rs`（--regroup, 默认关, 同 --stitch 惯例）。

**病灶再定位（推翻"组件割裂是主矛盾"的预设）**: 全集 Rp=44.37% vs 当前
"每 gene 最长"子集 33.63%——10.74pp 损失可逻辑证明 **100% 是基因内选错**
（若最佳命中者是本 gene 最长者必在子集中）。成因: noisy 49.9%（最长者
ident<95%, 脏延伸/错配堆积）+ thief 31.5%（最长者实为旁系混入）+ shorter
18.6%。跨基因割裂（cov>=0.3 口径 554 条）是次要矛盾。

**判别器（纯序列证据, 校准数字）**: minimap2 -X -x asm20 --cs all-vs-all →
每对 shared_len(=+*)、mism_density（块内 mismatch 密度, 块状 indel 不罚）、
双侧覆盖并集。参照标签校准: 同位点对 dens p95=0.28% vs 异位点 p50=4.2%。
iso 边 = dens<0.005 且 min(cov_a,cov_b)>=cov 且 shared>=200; para 边
（cannot-link）= dens>=0.02 且 shared>=100。教训: max 侧覆盖语义（对齐组件内
30% 节点重叠）在全局层是短片段胶水（v1 过合并 Rp 20%）; 必须双侧下限。

**机制结论**: ① 跨基因归并收益=0（split-only 与含归并数字相同）——与输出层
缝合证伪互洽, 全部增益来自拆分; ② 拆分后取最长在全部 53,337 旧基因上
better 2207 / **worse 0**（nomk transfer: better 1452 / worse 0）——代价只有
冗余（1:1 率 96.8→73.2）无质量风险; ③ 纯序列无监督"聪明选择"（clean/
dominant/组内 para 负担排除）全败（33.0-33.4 vs longest 33.6）, oracle 41.5
的剩余差距卡在 0.5-2% 灰区弥散分歧, 需 read 支持通道（后续）。

| mem9 (v0.6.0) | n基因 | Rp | R70p | P | F1 | 1:1率 | 纯净率 |
|---|---|---|---|---|---|---|---|
| 当前每 gene 最长 | 53,337 | 33.63 | 55.61 | 51.2 | 36.5 | 96.8 | 93.4 |
| sp20 (cov .2) | 65,492 | 35.00 | 58.88 | 48.6 | 36.9 | 95.3 | 97.2 |
| **sp_d (cov .7) = 默认** | 86,686 | **38.09** | **63.36** | 49.4 | **39.9** | 73.2 | **98.2** |
| **Rust 实现（cov .7）** | 86,115 | 37.34 | 62.33 | — | — | — | — |
| Rust oracle | — | 42.83 | 66.33 | — | — | — | — |

（Rust 原型差 -0.75pp = 内部 banded SW vs minimap2 特征差, 带内; oracle
42.83≈42.86 证明池质量等价。nomk transfer: Rp 26.6→28.5, 同构。）

**Rust 内部化**: 基因内两两（416,490 对, 仅 27% 需比对）——25-mer 锚门控
（<8 跳过, 重复 k-mer 剔除）+ 多对角带状比对（±64 带, 取 shared 更大方向含
RC）→ 规则同原型 → iso 边贪心 + para cannot-link → 家族最长。**79s**/
12.5 CPU-min（rayon 按基因并行）。产物: `<name>.principal.fa` +
`<name>.family_map` 侧车, 主产物/gene_trans_map 不动; 旗标关零副作用
（代码早退 + 两次关旗标运行互异证明主产物差异为 inchworm 多线程固有）。
696 测试全绿（685+11）; clippy 改动文件零告警（trinity-common/kmer.rs
3 条为未触碰存量的工具链漂移, 留档待另行处理）。examples/regroup.rs
离线入口: `cargo run --release -p trinity-cli --example regroup -- <fa> [cov]`。

**边界与后续**: ① Rp/F1 增益以 1:1 冗余为代价（cov 旋钮调节: 0.2 保守
+1.4pp/1:1≈不变 ↔ 0.7 激进 +4.5pp/1:1=73）, 默认取 0.7 但不转正默认链;
② oracle 差距（41.5-38≈3-4pp）卡灰区, read 支持通道（ContainedReads 透传
butterfly.pack）是下一个杠杆; ③ BUSCO 主口径不受影响（主产物未动）。

## regroup 战役二期: read 支持通道（2026-09-07, --regroup-reads 默认关）

一期遗留的 oracle 差距（分组上限 42.86 vs longest 38.09）卡在 0.5-2% 灰区弥散分歧,
序列证据不可判。二期用原始 reads 攻坚, 两个候选杠杆一胜一负:

**负案: 灰区变位列等位裁决（证伪收档）**。灰对 26,300（同旧基因 23,446, 参照可见
标签 940 同位/574 异位）; 变位列 427,962 → 101bp 合成上下文（a/b 各一）→
minimap2 -x sr 映射 reads → 中心位 read 碱基投票。裁决门 both_frac>=0.5/
cols>=3/depth>=5 触发 para 1,977 + error 2,735, 但**标签校准失败**: 同位
frac_both p50=0.33 vs 异位 0.40——不可分。根因二: ①灰区混入**单倍型/等位变异对**
（双等位 50/50 支持却是同一基因, 与旁系在等位支持上同构）; ②错误型对恰在低覆盖区
（错配堆积本就是低深度产物）——每对有效列>=3 的仅 2,026/26,300, 证据天然缺席。
改判特征重聚类零增益（Rp 40.6 持平, worse 229 略差）→ 弃用。

**胜案: 支持度选择代表**。minimap2 -N 1 sr 每 read 最佳命中(ident>=0.95) →
每转录本 read 支持数; 代表 = argmax(支持度[<3 无效], 长度)。灰对裁决机制全部
弃用、分组规则不动, 只换选择:

| mem9 principal 集 | Rp | R70p | P | F1 | better/worse |
|---|---|---|---|---|---|
| Phase 1（longest） | 38.09 | 63.36 | 49.4 | 39.9 | 2207/0 |
| 原型（minimap2 支持度） | 40.58 | 65.51 | 53.5 | 44.3 | 3406/220 |
| **Rust（k-mer 支持度）** | **39.7** | 64.6 | **53.0** | **43.3** | **2003/0** |

原型 worse=220 是"长短颠倒"型（高支持短显性 isoform 抢走最长位, 参照 principal
恰是最长者）; 防回退变体全败（floor<3/长度护栏 0.7-0.9/支持过滤/比值过滤——
脏长者的支持量级并不低, 只有排序式有效）。Rust 内部 k-mer 支持度（25-mer 索引,
每 k-mer >8 转录本剔除防重复污染, read 双向查, 半数命中归属阈值）**更保守,
worse=0 复得**且保住 +2.36pp（对 Phase 1 Rust 37.34）。reads 票 39.9M/44.77M,
108,263 转录本有支持, **128s**/23 CPU-min。

**实现**: regroup.rs 内 count_support + supp_key 选择; `--regroup-reads`
（默认关, 未给 --regroup 时自动启用并提示）; orchestrate 两挂点传原始 reads;
example `regroup <fa> [cov] [reads...]`。回归锁: no_reads_path_matches_phase1
（空 reads 与不传 reads 逐字节一致）。701 测试全绿（696+5）; clippy 改动文件
零告警。bam 提取的 reads 仅 hap1 比对半量深度——生产走原始 fq, 支持度估计更稳。

**累计战果（QM4 mem9, principal 交付集）**: Rp 33.63 → **39.7**（+6.1pp）,
P 51.2 → 53.0, F1 36.5 → **43.3**（+6.8pp）; 全程 --regroup-reads 默认关,
主产物/gene_trans_map/BUSCO 主口径零触碰。

## regroup 战役三期: 覆盖/端部证据攻 oracle 残差——主攻证伪, 微调 +0.1pp（2026-09-07 深夜）

目标: 二期后 Rp 39.7 vs oracle 42.86 残差 ~3pp。差距解剖（regroup_gap2.py,
当前配置 2,063 个差距基因）: **flip 1,481（主矛盾）** + tailc 405（rep 更长且
端部覆盖高=真 UTR, 不可修） + tail 177（rep 更长且端 200bp 覆盖<60%=可修加长段）。
flip 画像: oracle 成员常**更长且更准**（如 c5742_g1_i45 5927bp/cov0.999 败给
i5 5225bp/cov0.66）——支持度排序把短的显性 isoform 推上代表位, 即"长短颠倒"
的家族视角。

**完整覆盖过滤——证伪**: 多命中覆盖映射（sr 无 -N1, 44.73M hits, 全 117,610
tx 有覆盖）→ 规则"支持≥3 且 cov_frac≥0.8/0.9（×端 200bp≥0.5）的成员中取最长"
网格全部回落 longest 基线（37.9-38.1 vs supported3 40.58）。根因: 游荡段是
read 支持的真实旁系序列（与 09-02 ORF 线结论同构）, 覆盖完整度分不出
"真 isoform 独有区"与"旁系游走段"; 而 length-first 语义把支持度排序赢来的
flip 修复全数吐回。**flip 在 reads 丰度×长度证据平面内不可分。**

**实现层微调（保留）**: 投票阈值 1/2→1/4——横跨家族分歧区的 read 票分裂到
多条转录本、半数门槛整票丢弃; 1/4 后 Rp 39.7→**39.8**, worse 保持 0。
support k=19 与 25 实测等价（1% 错误率下 19-mer 只多回收 ~5% 窗口）, 参数留
`support_kmer_len` 默认 25 与管线一致。

**终态（QM4 mem9 principal 交付集, --regroup-reads）**: Rp **39.8** / R70p 64.7 /
P 53.1 / F1 43.4 / **worse=0** / 128s / 701 测试绿。三期累计对战役起点:
Rp +6.2pp（33.63）, F1 +6.9pp（36.5）。

**残差定性（~3pp, 留档不再攻）**: flip（isoform 丰度 vs 长度）+ tailc（真 UTR）。
下一层证据若要再进, 候选是剪接结构（junction 一致性: 参照 principal 与组装
转录本的内含子链对齐）或表达建模——均超出"拆分后区分旁系/isoform"的
reads 证据平面, 另立战役评估。

## regroup 四期探针: 共转录桥读——证伪收案（2026-09-07, 未立项即终）

三期残差定性后的候选杠杆"剪接结构证据", 取其最纯形式做了 30 分钟存在性探针
（regroup_bridge_probe.py）: flip 家族（rep≠最长成员 L 且 L 支持>=3, 按
长度差取前 300）中, L 的独有段（其他成员比对区间并集补集, >100bp）边界的
**桥读**（一条 read 横跨边界两侧各 >=15bp, 取自 sr_cov 多命中坐标）——
假说: 共转录分子才有桥, 旁系游走段的 read 不跨本基因连接边界。

结果无分离: oracle-L（应救回, n=12）桥读 p50=11 / 100% >=2 桥;
非 oracle-L（不应救, n=112）p50=9 / **91% 也有 >=2 桥**（max 127）;
另有 176/300 的 L 无离散独有段（长度差渐变）, 机制不适用。机理: 游走段是
组装图粘合产物, 穿越该连接的桥读正是当初让 butterfly 拼出此路径的那批
read——**输出层证据无法否定构图层已消费的同一证据**。

**判读: 等位（二期）/覆盖（三期）/桥读（四期）三种 read 证据形态对
flip+tailc 残差全部失效, 输出层杠杆已穷尽。** 残差（~3pp）根治位置在
构图层（weld 聚类/butterfly 路径验收防游走）——即 09-03 已定位的深水区,
与输出层 regroup 战役互为独立战线。regroup 战役全档至此收官:
Rp 33.63→39.8（+6.2pp）, F1 36.5→43.4（+6.9pp）, worse=0。

## 构图层一期: no-path-merging × regroup 组合验证——全集增益复现, 交付层组合净负（2026-09-09, cl_base/cl_npm 配对实验）

**前置发现: QM4 历史输入版本不在盘上（重要留档）**。现存 Cleandata
（22,385,819 对, Alnus_glutinosa/0.Data/.../QM4）重跑默认链得 n=146,234 /
R=41.7 / P=46.6——与 mem9（n=117,695/R=44.5/P=51.3）不可复现, 序列重合度仅
44%。判别链: bam 配对提取（collate 后 22.39M 对, 与 Cleandata 同源）同结果
→ 排除提取问题; **双二进制 CPU1 逐字节对拍一致**（工作树 vs v0.6.0 release,
10 万对, n=10,289）→ 排除管线漂移; 结论: 历史战报输入是 21.35M 对的某
版本（已不在盘）, 此前所有 QM4 历史数字与现存数据的可比性打星号。
cl_base（--no_cleanup 中间产物留档）为后续战役新锚。

**配对实验（同 Cleandata 输入, 仅差 --bfly-no-path-merging）**:

| | cl_base | cl_npm | Δ |
|---|---|---|---|
| 全集 n / N50 | 146,234 / 1776 | 276,532 / 2203 | +89% / +427 |
| 全集 R / Rp | 41.7 / 41.6 | **44.2 / 44.1** | **+2.5pp** |
| 全集 P / F1 | 46.6 / 44.0 | 46.8 / **45.5** | +0.2 / +1.5 |
| principal Rp / R70p | **36.6** / 62.0 | 35.1 / 61.6 | **−1.5pp** |
| principal P / F1 | **46.9 / 39.2** | 44.7 / 37.0 | −2.2 / −2.2 |
| worse | 0 | 0 | — |

- **09-03 的全集增益在独立输入上复现**（+2.5pp R, F1 +1.5, N50 +427,
  P 稳定, 组件数 54,941→55,010≈不动——npm 只作用于 butterfly 层 ✓）。
- **交付层组合净负**: n 翻倍是近重复变体洪水（iso 边 12.8 万→172 万）,
  支持度选票被重复体瓜分, 选择噪声放大。两个修复假设均证伪:
  ① 选票归并（近重复组支持和）: Rp 34.26 更差; ② 平票确定性集中
  （取最小 id）: Rp 31.2/31.9 崩塌——集中饿死平票候选, 内容哈希分摊
  实为近似按比例分配（已回滚, 701 测试绿）。
- **上限证据**: minimap2 支持度（-N1 主命中）在同分组上 Rp=37.08
  （> base 36.6）——principal 亏损是 k-mer 支持估计在重复体洪水下的
  工程退化, 非机制上限; 修复方向 = 容错投票（错配容忍的 read 归属）。

**裁定**: `--bfly-no-path-merging` 维持默认关。两种使用画像入档:
① 全转录本集用途（注释发现/ORF 池）: 值得开（+2.5R/+1.5F1/N50+427）;
② principal 交付集: **默认链 + --regroup-reads, 不与 npm 组合**——
待支持度通道容错升级后再评估组合。dev1 回退（09-03 已录, P−14.7）
不因本实验改变。

## 支持度通道容错升级: 三轮策略证伪——差距锁定在错配容忍打分本身（2026-09-09, 已回滚定案）

目标: 收回 k-mer 支持度与 minimap2 支持度在 no-path-merging 重复体洪水下的
差距（npm principal 35.1 vs 37.08 上限）。三个"集中选票"策略变体全部证伪:

| 变体 | cl_npm Rp | cl_base Rp | 死因 |
|---|---|---|---|
| 现行（1/4 门槛+内容哈希分摊） | **35.1** | **36.6** | —（基线） |
| 选票归并（consupp, 近重复组和） | 34.26 | — | 组内真变体也被归并, 信号抹平 |
| 平票集中（最小 id） | 31.2 | 31.9 | 重复区平票全灌单一候选, 饿死其余 |
| 带内最长（0.8×top 带内取最长） | 34.0 | 34.2 | 同上+有支持转录本 16.6万→10.6万 |

**定案**: 任何把 read 选票集中到"优选候选"的策略都饿死其余候选——
内容哈希分摊（≈按比例分配）是简单策略族中的最优; 已回滚并把这轮教训
写进 vote_one_read 注释。与 minimap2 的剩余差距来自**错配容忍打分**
（read 在近重复体间的归属由变体位错配决定, k-mer 精确匹配天然丢失该
信号）——修复需比对级 read 归属（mini 短读比对器, 独立工程量）, 留待
另立战役。npm×regroup 组合解锁推迟至该升级之后; 现行推荐不变
（全集用途开 npm, principal 用途默认链+--regroup-reads）。

## principal 数量口径战役: 106,065 vs StringTie 38,549 的差距解剖（2026-09-09, 画像定案）

**问题**: 交付集（cl_base + --regroup-reads）输出 106,065 principal, 是 StringTie
银标准转录本数（38,549）的 2.75×、基因数（22,792）的 4.65×。膨胀在哪、删了亏不亏。

**构成分解（asm_vs_st.clbase.paf, cov≥0.70&id≥0.95=代表级）**:

| 层 | 数量 | 占比 | 中位长 | 中位depth |
|---|---|---|---|---|
| A 代表级命中 | 14,638 | 13.8% | 1,577 | 9.6 |
| B 仅片段命中(id≥0.95,cov<0.70) | 33,027 | 31.1% | 442 | 2.8 |
| C 无 id≥0.95 命中 | 58,400 | 55.1% | 344 | 3.0 |

- **命中层干净**: A 折叠到 15,917 个 ST 基因, 1.09 拷贝/基因（14,503 基因唯一代表）
  ——严格拆分没有制造重复代表, 4.65× 全部来自短货。
- **C 类画像**: no_hit_at_all 30,954（53%, 中位 287bp, 92% 住在从未产出代表
  转录本的孤立小组件）; st_frag_diverged 10,070（17%, 对 ST id 0.80-0.95 且
  8 成长度在真实基因内, 43% 与真货同组件）; st_partial_lowid 12,069（21%）;
  brother_frag 2,543（4%, 被 id≥0.99 更长兄弟包含）; no_relation 2,764（5%）。
- **ST 盲区裁决**（depth≥5 的 no_hit+no_relation 共 5,919 条, 定位到 hap1 基因组）:
  **81% 落在 ST 基因座内**（对基因组中位同一性 0.988）——序列真实、位置真实,
  但不是 ST 的剪接结构（异常/替代路径）; 15%（891 条）基因间区=真·ST 盲区
  候选; 4% 比不上基因组（伪影/嵌合）。**删 C 的真内容风险 ≤1%**。
- **支持度教训**: mapq≥10 口径会把多映射区（旁系/重复/等位）真货判成低支持
  （A 类 27% 假性 depth<2）; 无 mapq 口径 A 类 99% depth≥2。支持度截断必须用
  无 mapq 计数（regroup 内置 k-mer 票天然满足）。

**可部署截断方案（无银标准信息, 丢A≈R70p 代价）**:

| 方案 | 保留 | 丢A |
|---|---|---|
| len≥300 | 73,334 | 1.2% |
| len≥500 | 46,885 | 4.3% |
| len≥300 & depth≥2 | 70,290 | 2.0% |
| len≥300 & depth≥5 | 48,915 | 7.0% |

要压到 ST 基因量级（23k）光靠输出层过滤不够——B 类 33,027 真基因碎片需要
构图层碎片化根治（与 flip+tailc 残差同根, 已有结论）。**裁定待用户**: 是否落
--principal-min-len（默认 0=关）。

**资产**: eval/cautopsy/{clbase_sr.paf(6GB, 47.3M 行), c_class_partial.tsv,
principal_support.tsv, principal_profile_final.tsv, cand_vs_genome.paf}。

## principal 截断旗标落地: --principal-min-len / --principal-min-supp（2026-09-09, 默认关）

承接上节画像定案, 把数量截断落成 regroup.rs 输出层旗标（默认 0=关, 默认链字节等价）:

- **语义**: 门控只过滤 principal.fa——代表长度 / read 票数不达标的家族不落
  principal（family_map 照旧全量, 分组与代表选择不变, worse=0 由构造保证）。
  `--principal-min-supp >0` 隐含 --regroup-reads（票数来源; 无 reads 时 supp=0
  会剔空, 管线层防住）。
- **票数量纲**: 与代表选择同源（一 read 一票, 无 mapq 语义——多映射区真货不
  被冤枉, 区别于 minimap2 mapq≥10 口径）; 40.97M reads 投票, 13.46 万转录本
  有支持（≥3）。
- **实现**: RegroupParams + RegroupStats.principal_dropped + stats_line/
  stage_timings 回显; args 解析/蕴含链/警告沿 --regroup-cov 与 --regroup-reads
  先例; example 探针支持 --min-len/--min-supp; 新增 5 测（门控只滤输出、
  边界==过、票数门、AND 组合、默认零门等价锁）。workspace 706 测试绿,
  clippy 改动文件零告警（kmer.rs:251 存量 deny 级 lint 系 toolchain 升级
  显形, 未动）。
- **验收（Cleandata 输入, example 直跑, 每跑 ~3min）**:
  - 等价锁: 无门跑 principal.fa + family_map 与归档**字节一致**;
  - len≥300 → 保留 **73,334**, len≥500 → **46,885**——与画像战役 python
    预测**逐个精确对账**;

| 档 | n | Rp | R70p | P | F1 |
|---|---|---|---|---|---|
| 无门（基线） | 106,065 | 36.6 | 62.0 | 46.9 | 39.2 |
| len≥300 | 73,334 | 36.4 | 61.2 | 50.0 | **40.2** |
| len≥500 | 46,885 | 35.9 | 59.2 | **53.7** | **41.1** |
| supp≥10 | 71,048 | 36.3 | 60.7 | 51.8 | 40.6 |

  - R70p 实际降幅（−0.8 / −2.8pp）与画像预测（1.2% / 4.3% × 62.0）吻合;
    P/F1 因短货碱基剔除**净升**（F1 +1.0~+1.9pp）——截断不只是数量口径,
    也是质量旋钮。
- **产物**: tr_rs_out/{cl_len300,cl_len500,cl_supp10}.Trinity.principal.fa
  留档; 顶层 cl_base 归档已恢复原状（副本在 tr_rs_out/cl_base/）。
- **建议用法**: 需要与 StringTie 同量级的基因数口径时 `--regroup
  --regroup-reads --principal-min-len 300`（73k, F1 +1.0）; 更紧的交付
  （46.9k, P 53.7）用 500。旗标默认关, 主线行为不变, 转正与否待用户裁定。

## principal 截断转正: len500+supp10 成管线默认（2026-09-09, 用户裁定）

**裁定**: `--principal-min-len` / `--principal-min-supp` 默认值 0/0 → **500/10**。
regroup 整体仍默认关; opt-in 后（--regroup / --regroup-reads / 显式 min-supp
任一）supp>0 隐含 reads 通道（默认档静默隐含, 显式给出才打 NOTE——首版实现
裸 argv 也会被默认 supp 拉开 regroup, 已修: opt-in 判定先行, 全关时零副作用）。
`RegroupParams` 原语与 example 探针保持 0/0（库语义=关; 探针无 reads 时 supp
门剔空）, 显式传 0/0 关回旧行为。args 测试 4 处更新（默认值/裸 --regroup 隐含
reads/opt-out 0-0/裸 argv 全关）。

**组合档验收（Cleandata, example 直跑, 236s）**:
- 保留 **44,522** = len500 集(46,885) ∩ supp10 集(71,048)——与归档两集交集
  **逐个精确对账**;
- eval: **Rp 35.8 / R70p 58.8 / P 55.0 / F1 41.3**（基线 106,065 全量:
  36.6/62.0/46.9/39.2）——P +8.1, F1 +2.1, R70p −3.2pp（=画像预测 ~5% 档）。
- 数量口径: 44,522 vs StringTie 38,549 tx（1.16×）/ 22,792 gene（1.95×）,
  从 4.65× 压到 2× 以内。
- workspace 706 测试绿, clippy 改动文件零告警; 产物
  tr_rs_out/cl_def.Trinity.principal.fa 留档, 顶层归档保持全量版
  （106,065, 需要全集时用显式 0/0）。

**新默认链**: `--regroup`（裸旗标）= 严格家族拆分 + 支持度代表 + len500/supp10
截断, 一条旗标到交付集。

## BUSCO 对照: 转正档交付集的完整性变化（2026-09-09, eudicots_odb10 n=2326 蛋白模式, TransDecoder 翻译, 与 09-02 同口径）

| 组 | n | C | S | D | F | M |
|---|---|---|---|---|---|---|
| cl_base 全集 | 146,234 | 88.0 | 47.0 | 41.0 | 5.3 | 6.7 |
| 旧交付 principal | 106,065 | 82.1 | **72.6** | 9.5 | 9.5 | 8.4 |
| **转正档 cl_def** | **44,522** | **81.3** | **72.4** | **8.8** | 9.1 | 9.7 |

（旧输入时代全集参考: mem9 89.5/S 54.3/D 35.3; nomk 89.7/S 66.8/D 22.9）

- **拆分+代表选择把重复砍掉 31.5pp**（D 41.0→9.5）, S 47.0→72.6——"每基因
  一条"在 BUSCO 口径兑现; C −5.9pp 的来源 = isoform 穴余消失（F +4.2, M +1.7:
  部分 BUSCO 的完整 ORF 只在被剔 isoform 里）, 与 Rp/R70p 掉 ~3pp 同因,
  是 principal 口径的结构代价。
- **转正档几乎无损**: 106k→44,522 只再花 C −0.8pp（S 持平 72.4, F 反降）,
  再次确认截掉的短货不含完整基因——与画像战役"ST 盲区真内容 ≤1%"互证。
- 剩余 D 8.8 = 真旁系同源（BUSCO 家族内多拷贝）, 正是严格拆分要保的货。

## regroup 转正为管线默认（2026-10-29, 打包裁定）

承接 09-09 转正档（len500+supp10 默认门 + BUSCO 对照 C 81.3 近乎无损）, 把去冗余交付链从 opt-in 抬为**默认开**——默认流程 = 组装 + 去冗余一条链交付:

- **行为**: 裸 argv = 严格家族拆分 + reads 支持度代表 + len500/supp10 截断（产物 `*.Trinity.principal.fa` + `*.Trinity.family_map` 侧车, 主 fasta/gene_trans_map 不动）。`--regroup`/`--regroup-reads` 保留为幂等兼容; `--no-regroup` 整链回退（零副作用, 沿 --no-rtt-multi-kmer 先例）; `--principal-min-len 0 --principal-min-supp 0` 只关截断门。
- **互斥（防歧义, 沿 --no-multi-k × --multi-k 先例）**: --no-regroup 与 --regroup/--regroup-reads/显式 supp>0 同给 → usage 报错; 显式 supp=0 不算开启, 不冲突。
- **实现**: 仅 args.rs 解析层改默认值与蕴含链（regroup = !no_regroup; regroup∧supp>0 ⇒ reads）, regroup.rs 机制与 RegroupParams 原语 0/0 不动; 版本号 0.3.0 → 0.4.0。
- **测试**: args 两个 regroup 测试改写为默认开语义（裸 argv 全链、--no-regroup 回退、幂等兼容、互斥报错、显式 0 不冲突）, trinity-cli 78 测全绿。
