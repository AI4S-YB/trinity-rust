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
