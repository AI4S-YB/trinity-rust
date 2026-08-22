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
