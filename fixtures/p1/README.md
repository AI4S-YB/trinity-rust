# fixtures/p1 — P1 trinity-kmer 对拍 fixtures

各 fixture 的来源与再生成命令。统一入口：`cargo xtask xcheck-kmer`（三项 oracle 对拍，
临时产物在 `target/xcheck/`，不落仓库）。

环境约定（再生成命令中引用）：

```bash
JF=/public/home/senior007/miniconda3/envs/trinity/bin/jellyfish   # 可用 JELLYFISH 覆盖
TRINITY_SRC=/storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2
STATS=$TRINITY_SRC/Inchworm/bin/fastaToKmerCoverageStats          # 需已 make
NORM=$TRINITY_SRC/util/insilico_read_normalization.pl
SEQTK=$TRINITY_SRC/trinity-plugins/seqtk-trinity/seqtk-trinity    # 需已 make
# perl 一律 /usr/bin/perl（含 DB_File；conda perl 缺该模块）
# 原版管线 PATH 必须让 tarball 编译版 seqtk-trinity 优先于 conda env 的（后者段错误）
```

## 输入 reads

| 文件 | 来源 | 再生成 |
|---|---|---|
| `smoke.fa` | trinityrnaseq `sample_data/test_Trinity_Assembly/__indiv_ex_sample_derived/ex01/ex01.reads.left.fq`（在 trinityrnaseq-Trinity-v2.15.2 解包内）前 50 条 reads 转 FASTA（名/序列原样，均 76bp） | `awk 'NR%4==1{sub(/^@/,">");print} NR%4==2{print}' <ex01.reads.left.fq \| head -100 > smoke.fa` |
| `edge.fa` | 手造边界样本：20bp 短 read、含 N、3 份相同序列（dup）、小写碱基、多行折行、恰好 25bp 的 homopolymer | 手写，勿再生成 |
| `diginorm/pe.l.fq`, `diginorm/pe.r.fq` | 手造 6 对 PE（pair1/2 完全同序列、pair3 含 N、pair4 20bp 短 read、pair5/6 唯一）；diginorm_e2e 测试头注释有手推名单 | 手写，勿再生成 |
| `diginorm/ss.pe.l.fq`, `diginorm/ss.pe.r.fq` | SS F 互补链（左 X / 右 revcomp(X)，dUTP 型）：ssA×2/侧、ssB 左1右2、ssC×3/侧 | 手写，勿再生成 |
| `diginorm/ss.single.fq` | 单端混链（同一文件含 X 与 revcomp(X)）：ssA 正/反各2、ssB 正1反2、ssC 正/反各3 | 手写，勿再生成 |
| `diginorm/pe.l.fa`, `diginorm/pe.r.fa` | 原版 prep 产物：`seqtk-trinity seq -A -R 1/-R 2` 对 pe.{l,r}.fq 的名规范化中间 fa | `$SEQTK seq -A -R 1 pe.l.fq > pe.l.fa`（`-R 2` 同理） |

## oracle 黄金产物（对拍基准，不重新生成除非 oracle 升级）

| 文件 | 来源 | 再生成 |
|---|---|---|
| `smoke.kmers.fa`, `edge.kmers.fa` | `$JF count -m 25 --canonical -s 100M <reads> && $JF dump -L 1 mer_counts.jf` | 见左（xcheck [1] 每次现算，不依赖此文件） |
| `smoke.stats.orig.tsv`, `edge.stats.orig.tsv` | `$STATS --reads <reads> --kmers <kmers.fa> --kmer_size 25 --num_threads 1`（DS 默认） | 见左 |
| `g6_golden.tsv` | C++ `%g`（6 位有效数字）格式化黄金：printf 直出 | `printf '%g\n' <值>` 系列，锁定 -0/NaN/指数形态 |
| `avg1f_golden.txt` | perl `sprintf("%.1f", (a+b)/2)` 黄金（nbkc merge 的 avg_1f 语义） | `perl -e 'printf "%.1f\n", (1.25+1.35)/2'` 等 |
| `perl_rand_12345.txt` | `perl -e 'srand(12345); print rand(1), "\n" for 1..20'`（perl ≥5.20 即 drand48，与 Drand48::new(12345) 位级一致） | 见左 |
| `seqtk_names.fa` | `$SEQTK seq -A -R 1 /tmp/p1_names.fq`（名规范化样本：/1 尾巴、_forward/_reverse、空白描述） | 见左（输入 fq 手写） |
| `diginorm/pe.{l,r}.stats.sort.golden.tsv` | 原版管线中间产物 `left/right.fa.K25.stats.sort` 切前 4 列（tid 列切除） | 原版管线 `--no_cleanup` 跑 pe fixture 后从 `tmp_normalized_reads/` 切列 |
| `diginorm/pairs.stats.golden.tsv` | 原版管线 `pairs.K25.stats` 的 acc+合成 3 列（cut -f1,10,11,12） | 同上 |
| `diginorm/ss.pe.left.norm.golden.fq`, `ss.pe.right.norm.golden.fq` | 原版管线 SS F PE 输出（`--SS_lib_type F --pairs_together --max_cov 200 --min_cov 3`） | 见下「原版 diginorm 管线」 |
| `diginorm/ss.single.norm.golden.fq` | 原版管线单端混链输出（同上参数，`--single`） | 同上 |

### 原版 diginorm 管线（golden 的生成方式）

```bash
export PATH="$TRINITY_SRC/trinity-plugins/seqtk-trinity:$(dirname $JF):$PATH"
/usr/bin/perl $NORM --seqType fq --JM 1G --max_cov 200 --min_cov 3 --CPU 2 \
  --output <outdir> --SS_lib_type F --pairs_together \
  --left fixtures/p1/diginorm/ss.pe.l.fq --right fixtures/p1/diginorm/ss.pe.r.fq
# <outdir>/left.norm.fq → ss.pe.left.norm.golden.fq（right 同理; 单端把 --left/--right 换 --single）
```

## xcheck-kmer 三项对拍（live，不用上表黄金）

1. dump 多重集：`$JF count -m 25 -s 100M [-canonical] && dump -L 1` vs
   `trinity-kmer count`——DS/canonical 双模式 `sort` 后逐字节；
2. coverage-stats：`$STATS ... --num_threads 1` vs `trinity-kmer coverage-stats`——
   tid 列切除后排序比较（smoke 用 [1] 的 jellyfish dump，edge 用 checked-in kmers）；
3. diginorm 端到端：`$NORM --pairs_together`（perl/PATH 见上）vs `trinity-kmer diginorm`
   ——三组 PE（DS maxC200 / SS-F maxC200 / DS maxC2 rand 路径）left/right.norm.fq 逐字节。
