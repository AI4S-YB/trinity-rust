# 基准记录

## 2026-08-16 · trinity-kmer count vs jellyfish count+dump（首条）

- **数据**: sample_data/test_Trinity_Assembly PE reads（reads.left/right + reads2.left/right 合并，
  按 DigiNorm both.fa = cat left right 的真实用法），取前 400k 行 = **100,000 条 reads**，
  fq→fa 后 `/tmp/bench.fa` = 11.3 MB。
- **参数**: `-m 25 --canonical`；线程对齐 4（jellyfish `-t 4`，Rust `RAYON_NUM_THREADS=4`）。
- **规模**: 去重 canonical 25-mer（计数 ≥1）= **518,514**；dump 输出 15.1 MB。
- **机器**: 128 核登录节点（本机）。

| 步骤 | wall（3 次取代表） | RSS |
|---|---|---|
| jellyfish count -s 100M -t 4 | 0.84 s（0.82/0.85/0.84） | 626 MB |
| jellyfish dump -L 1 | 0.11 s | 3 MB |
| **jellyfish count+dump 合计** | **~0.95 s** | **~626 MB** |
| Rust trinity-kmer count（含 dump 输出，单进程） | **0.30 s**（0.31/0.30/0.29） | **86 MB** |

注：spec 命令原样 `-s 2G` 时 jellyfish 为 8.82 s / 8.7 GB RSS——大头是 2G 槽位哈希预分配
（518k k-mer 用不到），对本文数据无代表性，故公平对比取尺寸合适的 `-s 100M`。

- **等价性**: `cmp <(sort b.jf.dump) <(sort b.rs.dump)` → **BENCH-EQUIV**（排序后逐字节相等；
  两方遍历序均为哈希序不定，等价性按多重集判定）。
- **结论**: 同等 4 线程下 Rust 计数+落 dump 一体化比 jellyfish count+dump 两步快 ~3.2×、
  峰值内存 ~1/7，输出多重集逐字节等价。

## P1 最终审查补记（2026-08-16）：计数器内存随线程数放大（P5 工作项依据）

稠密数据实测（200k PE reads/100bp，~7.9M distinct 25-mer，45.8MB both.fa）：

| 场景 | 峰值 RSS |
|---|---|
| diginorm 全流程（88MB PE fq，默认 128 线程） | 2.32 GB（~26× 输入） |
| 同上 RAYON_NUM_THREADS=4 | 1.04 GB（~12×） |
| 仅 count：4T / 16T / 64T / 128T | 824 MB / 1.25 GB / 1.94 GB / 2.22 GB |

**结论**：rayon fold 的每个并发累加器在稠密数据上趋近全表大小，峰值 ≈ 表 × 并发段数，随核数线性增长。P1 计划偏离注记的"~18B/kmer"被此实测证伪（该估算仅在固定 4 线程的基准下近似成立）。外推 10M PE reads（4.4GB fq）峰值 25–45+ GB。

**P5 正式工作项**（据此立项）：
1. 分片开放寻址表（按 kmer 哈希分片，内存 = 1× 表、与线程数无关）——spec §5.2 原设计
2. trinity-common 按记录块边界的零拷贝迭代器（替代整文件读入 + Vec<Vec<u8>> 拷贝）
3. diginorm 各阶段管道化（raw/fa/both/stats 行/排序副本多份并存 → 流式）
4. --JM/--max_memory 护栏接线（estimate_hash_size 目前无生产调用方）
5. 多文件逗号列表 / --left_list（ReadsInput 单路径 → 列表）

## 2026-08-16 · trinity-inchworm vs 原版 inchworm（P2 门）

- **数据**: sample_data/test_Trinity_Assembly `reads.left.fq.gz` 全量 fq→fa（3.4 MB）→
  jellyfish `-m 25 --canonical -s 100M` dump = **288,470 kmers**（8.4 MB dump 文件）。
- **参数**（主线形态，Trinity:2654-2712 同形）: `--kmers <dump> --run_inchworm -K 25
  --monitor 1 --DS --min_any_entropy 1.0 -L 25 --no_prune_error_kmers`；并行轮加
  `--num_threads 6 --PARALLEL_IWORM`（单线程轮 `--num_threads 1`）。
- **oracle**: `trinityrnaseq-v2.15.2/Inchworm/bin/inchworm`（本机已 make）。
- **方法**: 每配置 3 次独立复测，`/usr/bin/time -v` 取 wall 与峰值 RSS；
  机器同前（128 核登录节点）。二进制均为 release/profile 构建。

| 配置 | 原版 wall ×3 | 原版 RSS 峰值 | Rust wall ×3 | Rust RSS 峰值 |
|---|---|---|---|---|
| 单线程 | 0.72 / 0.73 / 0.52 s | 26.1 MB | **0.36 / 0.41 / 0.35 s** | **22.5 MB** |
| 6 线程 PARALLEL | 0.72 / 0.81 / 0.80 s | 22.1–25.0 MB（随竞态波动） | **0.36 / 0.37 / 0.41 s** | 24.2–25.2 MB |

- 单线程 ~**2.0×** 快、RSS 低 ~14%；6 线程 PARALLEL ~**2.1×** 快、RSS 持平。
  该输入规模（28 万 k-mer，亚秒级）下并行对**两侧都无收益**——线程池/竞态开销主导
  （原版还多 per-thread tmp 文件往返）；扩展数据下的并行收益待 P3+ 复测。
- **contig 数**: 单线程 824 vs 824（两侧同）；PARALLEL 双方 run-to-run 漂移
  （原版 838/836/833，我们 834/828/841）——弱一致清零竞态，by-design
  （各自两次运行的 rc 不变多重集自差也各有 26/838、22/834）。
- **等价性**: 单线程 rc 不变序列多重集差 **26/824（3.2%）**，在既定种子平局带内
  （T5/T7 轮 822 contig 实测 **3.6%**（kmers 模式）/ 4.4%（reads 模式）；本轮 dump
  为 jellyfish 重建、哈希序漂移致数字小幅波动）；以原版 `--monitor 2` 种子序**重放则
  逐字节一致**（smoke / sample_data kmers / reads 三组 BYTE-MATCH，T5/T7 结论）——
  残余差异全部归因种子平局排序序，非贪心核心。
- **结论**: 主线参数形态下 Rust inchworm 单线程与 6 线程 PARALLEL 均约 2× 于原版、
  峰值内存单线程低 14%（并行持平），正确性由重放逐字节 + 多重集平局带双锚定。

## 2026-08-17 · trinity-chrysalis 全链路 vs 原版六程序串联（P3 门）

- **数据**: fixtures/p3 链（sample_data 真实数据）：`gff.iworm.fa`（219 contig，
  原版 inchworm 产物过滤所得）+ `gff.reads.fa`（61,150 reads）→ 55 组件 / Cbin0。
- **双方管线**: 原版 = GraphFromFasta → `sort -k9,9gr` → BubbleUpClustering →
  CreateIwormFastaBundle → ReadsToTranscripts → `LC_ALL=C sort` → FastaToDeBruijn
  → partition perl → QuantifyGraph ×55（shell 串联，`/usr/bin/time -v` 计全程）；
  Rust = `trinity-chrysalis chrysalis-all -t N` 单进程全链（内含排序与逐组件 quantify）。
  参数同 Trinity:2180 实参形态。二进制均为 release 构建；每配置 3 次独立复测。

| 配置 | 原版 wall ×3 | 原版 RSS 峰 | Rust wall ×3 | Rust RSS 峰 |
|---|---|---|---|---|
| -t 1 全链 | 12.09 / 11.39 / 11.45 s | 466 MB | **9.31 / 9.98 / 10.35 s** | **147 MB** |
| -t 4 全链 | 4.40 / 4.69 / 4.29 s | 466 MB | 8.76 / 9.04 / 8.77 s | 147 MB |

- 分阶段（-t 4，大头三项；3 次代表值）：

| 阶段 | 原版 wall | Rust wall | 备注 |
|---|---|---|---|
| graph-from-fasta | 0.97–1.04 s | **0.49–0.56 s**（~1.9× 快） | |
| reads-to-transcripts | **2.52–2.98 s** | 6.26–6.57 s（~2.4× 慢） | Rust t4 与 t1 几乎持平——本版 rtt 并行度不足，P4 优化项 |
| quantify ×55 组件 | 4.19–4.24 s | **2.22–2.39 s**（~1.8× 快） | 逐组件串行，双方同形态 |

- **RSS**: 全链峰值 147 MB vs 466 MB（~1/3.2）——原版大头在 GraphFromFasta 的
  weldmer/对齐结构（其 RSS 不随 -t 变化）。
- **等价性**: xcheck-chrysalis 7/7（边多重集 24 边全等 / bundle 逐字节 / 分区逐字节
  55 组件 / quantify 刀口边白名单 c0 8 + c2 4 行 / chrysalis-all 自一致）。
- **结论**: 单线程全链 Rust 快 ~1.2×、峰值内存 ~1/3.2；-t 4 下原版 OMP 并行
  （GFF + ReadsToTranscripts）占优（~2.1× 于我们）——我们的 ReadsToTranscripts
  多线程未生效是唯一短板，已记 P4 优化项；正确性由七重对拍锚定。

## 2026-08-18 · trinity-butterfly vs 发布版 Butterfly.jar（P4 门）

- **数据**: fixtures/p3/quantify 的 c0/c1/c2（orig.graph.out + orig.reads.out，
  原版 Chrysalis 真实产物；reads 条数 N = 4342 / 9773 / 7669）。
- **形态**: EM 主线（`-N <reads> -L 200 -F 10000 -R 2 -C prefix`，jar 加 `-V 10`）。
  jar = `$TRINITY_SRC/Butterfly/Butterfly.jar`（发布版 312f…，裁定见 setup.md），
  `java -Xmx4g`；Rust = `target/release/butterfly`。每配置 3 次独立复测，
  `/usr/bin/time -v` 取 wall 与峰值 RSS；机器同前（128 核登录节点）。

| 组件 | jar wall ×3 | jar RSS 峰 | Rust wall ×3 | Rust RSS 峰 |
|---|---|---|---|---|
| c0 (N=4342) | 1.28 / 1.09 / 1.11 s | 239–249 MB | **0.03 / 0.03 / 0.03 s** | **9.2 MB** |
| c1 (N=9773) | 1.82 / 1.54 / 1.99 s | 358–493 MB | **0.55 / 0.55 / 0.55 s** | 355 MB |
| c2 (N=7669) | 1.71 / 1.79 / 1.65 s | 347–360 MB | **0.66 / 0.67 / 0.66 s** | 403 MB |

- **wall**: Rust 快 **~2.5×（c1/c2）至 ~40×（c0）**——c0 上 jar 的 JVM 启动 + 类加载
  （~1s）占绝对大头；Rust 侧无固定开销，小组件差距被放大。
- **RSS**: c0 上 Rust 9.2 MB（jar 的 ~1/27）；c1/c2 上两侧同量级（Rust 略高 ~10%，
  组件线程大栈 + 全量 read 装载），jar 的 JVM 基线 (~200MB) 使其下限更高。
- **等价性**: xcheck-butterfly 6 检查点 = 3 PASS + 3 PASS-WARN
  （c0/em 序列多重集全等仅顺序差；c2 各多 1 条 [2182] 短异构本——jar 侧丢弃，
  黄金序列全覆盖，见 porting-map 验证资产节）。
- **结论**: 主线参数形态下 Rust butterfly 对小组件快一个量级、对中组件 ~2.5×，
  峰值内存小组件低 ~27×、中组件持平；正确性由 8 层对拍 + 组件端到端锚定。

## 2026-08-18 · P5-T4：全管线 RSS 基准（护栏接线 + backlog 裁定依据）

- **数据**: sample_data/test_Trinity_Assembly 全量 PE（reads+reads2 各侧合并，
  fq.gz 4 文件合计 ~5.0 MB gz / 118 MB 解压 fq；diginorm 后 both.fa ≈ 21.7 MB）。
- **命令**: `trinity-cli --seqType fq --left ... --right ... --max_memory 2G --CPU 8 --output ...`
  （`/usr/bin/time -v` 取 wall 与峰值 RSS；release 构建）。

| 配置 | wall | 峰值 RSS | 产出 |
|---|---|---|---|
| `--CPU 8`（inchworm 6T PARALLEL） | **17.4 s** | **868 MB**（0.40× max_memory） | 92 transcripts |
| `--CPU 1` | ~35 s | <868 MB | 91 transcripts（inchworm 种子平局序差异, 原版固有） |

- **确定性**: `--CPU 1` vs `--CPU 8 --inchworm_cpu 1` 最终 `*.Trinity.fasta`
  **逐字节一致**——GFF Phase1 / RTT 查询 / quantify 组件池的显式 -t 池
  均与线程数无关。
- **RTT 查询并行**（单元基准, 2000 reads × 300nt × 4 bundles）: -t1 3.50 s →
  -t4 0.99 s（**3.5×**），输出逐字节一致。
- **裁定**: RSS 0.40× max_memory → 分片表/流式/diginorm 管道化记入
  `docs/backlog.md`（含 revisit 触发条件），P5 不做。

## 2026-08-19 · P5 收官：端到端 wall + 断点续跑演示

### 全管线 wall 对照（sample_data，fast 档 = 截断 50000 PE reads）

`cargo xtask xcheck-trinity` 同机双侧计时（release，--CPU 8 --max_memory 2G）：

| 侧 | wall |
|---|---|
| 原版 Trinity（--no_bowtie --no_salmon） | 32.3 s |
| trinity-cli | **24.3 s**（独立直跑冷缓存下同输入 9.9 s，见下行演示首跑） |

全量（T4 记录，同上表命令不变仅去截断）：`--CPU 8` **17.4 s** / 峰值 RSS 868 MB / 92 transcripts。

### 断点续跑演示（`/tmp/p5-resume-fast`，50000 PE reads，--CPU 8）

同一 workdir 三次运行（`.ok` checkpoint 与原版同名同位）：

1. **全量首跑**（空目录）→ `FULL_WALL=9.90 s`，RSS 798 MB，94 transcripts。
2. **删 chrysalis 树 + `.quantify_graph.ok`/`.butterfly.ok` + 最终 fa 后重跑**——
   stderr 三行跳过（归一化/jellyfish/inchworm 链）：
   ```
   ---- Trinity (rust) checkpoint found, skipping: .../insilico_read_normalization/normalization.ok
   ---- Trinity (rust) checkpoint found, skipping: .../.iworm.25.asm.ok
   ---- Trinity (rust) checkpoint found, skipping: .../iworm_renamed.25.asm.ok
   ```
   → `RESUME_WALL=5.71 s`；最终 `*.Trinity.fasta` **序列多重集与首跑完全相等**
   （header 内图节点编号为并行编号差异，属已知带）。
3. **什么都不删再跑**——5 个 checkpoint 全部跳过，只剩 harvest 汇总 →
   `NOSKIP_NEEDED_WALL=0.14 s`。

反例（语义边界）：若删除 `inchworm.fa` 但保留 `.iworm.25.asm.ok`，下游
chrysalis 报 `no Inchworm output is detected` 非零退出——与原版 Pipeliner
"产物缺失但 .ok 在"行为一致（prep 侧会重建，chrysalis 侧要求上游产物在）。
