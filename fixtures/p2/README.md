# fixtures/p2 — P2 trinity-inchworm 黄金向量

统一再生成入口：`cargo xtask gen-fixtures`（同时再生 P0/P1 的 kmer_golden）。
环境约定同 `fixtures/p1/README.md`（`TRINITY_SRC`、conda gxx 可用 `TRINITY_GXX` 覆盖）。

| 文件 | 来源 | 再生成 |
|---|---|---|
| `hash_golden_input.txt` | 手造 20 行序列：单碱基 GATC、GATC/CATG/ACGT 置换、小写/混合大小写、N 与含 N 序列、homopolymer、空串（1 行）、15/27/36-mer（36-mer 触发 32 位回绕） | 手写，勿再生成 |
| `hash_golden.tsv` | C harness（`xtask/fixtures-src/dump_hash_golden.cpp`，链原版 sequenceUtil.cpp）对上行输入的 `generateHash` 输出，TSV `seq\thash_u64` | `cargo xtask gen-fixtures` |
| `glibcrand_seed1.txt` | glibc `srand(1); rand()` 前 100 值（原版 inchworm 未 srandom → 运行时 rand() 即此序列） | `cargo xtask gen-fixtures` |
| `glibcrand_mod2.txt` | glibc `srand(1); rand()%2` 前 50 值（inchworm_step 真平局二选一路径） | `cargo xtask gen-fixtures` |
| `smoke.kmers.fa` | P1 `smoke.kmers.fa` 的原样拷贝（jellyfish K25 canonical dump of smoke.fa，1887 kmers）——xcheck-inchworm 的默认输入 | `cp ../p1/smoke.kmers.fa .` |
| `smoke.orig.fa` | 原版 `Inchworm/bin/inchworm --kmers ../p1/smoke.kmers.fa --run_inchworm -K 25 --monitor 1 --DS --num_threads 1` 的 stdout（6 contig） | 原版二进制重跑（命令见 `trinity-inchworm/tests/smoke_vs_original.rs` 头注） |
| `smoke.seed_order.orig.tsv` | 同上运行 `--monitor 2` 抓取的通过闸门种子序（`SEED kmer: X, count: N` → `X\tN`，12 行）——**贪心核心逐字节锁定的重放输入** | 同上 |

harness 说明：

- `dump_hash_golden.cpp` **不跳过空行**——输入中间的空行即空串样本（原版
  generateHash("") = 0）；
- `dump_glibcrand_golden.cpp` 是独立小程序（不链 Trinity 源），oracle 即本机
  glibc；参数 `raw`/`mod2` 两模式各自独立 `srand(1)`。

Rust 对拍：`trinity-common/tests/seq_hash_golden.rs`（20 行逐行 u64 相等）与
`trinity-inchworm/src/glibc_rand.rs` 内嵌黄金测试（100 + 50 值位级一致）。

inchworm 冒烟对拍（P2-T5，`trinity-inchworm/tests/smoke_vs_original.rs`）:

- **重放测试**（`replay_with_original_seed_order_is_byte_identical`）: 以
  `smoke.seed_order.orig.tsv` 的原版种子序重放 populate → 默认剪枝 →
  `compute_sequence_assemblies_from_seeds`，stdout 与 `smoke.orig.fa`
  **逐字节一致**——证明贪心延伸/tie 打破/glibc rand/清零/去重/格式化全链路
  位级复刻，与原版的残余分歧只能来自种子平局序。该实验在 sample_data 全量
  （jellyfish K25 canonical dump，288,470 kmers，822 contig）与 `--reads`
  模式（默认 prune_error_kmers 路径）同样 BYTE-MATCH（见测试头注）。
- **CLI 测试**: 默认种子序（count 降序 + kmer 值降序平局）跑 CLI，断言 contig
  数、header 多重集（去 aN）与 **rc 不变序列多重集**与原版一致（DS 语义下
  contig 与其 revcomp 同义——smoke fixture 上 rc 不变多重集完全相等，残余
  仅链方向;sample_data 上 rc 不变差 15-18/822，为分支点平局划分差异，已由
  重放实验归因）。

## xcheck-inchworm 四重对拍（live，不用上表黄金）

`cargo xtask xcheck-inchworm [--kmers <fa>] [--reads <fa>]`（默认
`fixtures/p2/smoke.kmers.fa` 与 `fixtures/p1/smoke.fa`;需 `TRINITY_SRC` 下原版
`Inchworm/bin/inchworm` 与 `Chrysalis/bin/GraphFromFasta` 已 make;临时产物在
`target/xcheck/`）:

1. **单线程对拍**: 两侧 `--kmers <f> --run_inchworm -K 25 --monitor 1 --DS
   --num_threads 1`——rc 不变序列多重集 + header（去 aN）多重集，smoke 期望全等;
2. **PARALLEL 对拍**: 两侧 `--num_threads 4 --PARALLEL_IWORM -L 25
   --no_prune_error_kmers`——rc 不变多重集（header 不比: PARALLEL 下同一 contig
   的 Seed 值随 chunk 划分漂移）;两侧同为多线程竞态，严格模式重试至多 10 次
   命中「竞态窗口外」的全等;
3. **--reads 模式对拍**: 两侧 `--reads <f> --run_inchworm -K 25 --monitor 1 --DS
   --num_threads 1`（默认 prune_error_kmers=true 路径），判定同 [1];
4. **喂原版 Chrysalis**（P2 门）: 我们的 [1] 输出 → `GraphFromFasta -i <ours.fa>
   -r <reads.fa> -k 24 -kk 48 -min_glue 2 -glue_factor 0.05 -min_iso_ratio 0.05
   -t 1`（参数即 Trinity 主脚本 `Trinity:2180` 的默认: -k K-1=24、-kk 2(K-1)=48），
   要求 exit 0 且 stdout（weld 图）非空;smoke 输入本就无 kk=48 重叠候选
   （原版输出同样 0 行），以原版输出同参对照判定。

`--kmers`/`--reads` 显式指定（大输入）时对应检查降为统计报告: rc 多重集差异率
<5% 为既定平局带（报告注明，不 FAIL）;≥5% 打警告。大输入轮（sample_data 全量，
288,470 kmers / 822 contig）重建:

```bash
JF=/public/home/senior007/miniconda3/envs/trinity/bin/jellyfish   # 可用 JELLYFISH 覆盖
SD=/storage/home/senior007/test/trinity_rust/trinityrnaseq-Trinity-v2.15.2/sample_data/test_Trinity_Assembly
mkdir -p /tmp/t5
zcat $SD/reads.left.fq.gz | awk 'NR%4==1{sub(/^@/,">");print} NR%4==2{print}' > /tmp/t5/reads.fa
$JF count -m 25 --canonical -s 100M -t 4 /tmp/t5/reads.fa   # 产出 mer_counts.jf
$JF dump -L 1 mer_counts.jf > /tmp/t5/kmers.fa
cargo xtask xcheck-inchworm --kmers /tmp/t5/kmers.fa --reads /tmp/t5/reads.fa
```

（重建的 dump 与原 /tmp/t5 逐字节不同——jellyfish dump 序随运行漂移，k-mer
多重集一致（288,470 条 `sort` 后全等），单线程对拍结果不受影响。）
