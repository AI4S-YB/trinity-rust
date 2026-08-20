# 环境搭建

## 原版 Trinity（仅供交叉验证）

- 源码: /storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2 (v2.15.2 FULL tarball)
- 构建: `make inchworm_target chrysalis_target`（Butterfly.jar 预编译；两者均为 cmake 包装的 Makefile）
- 关键二进制:
  - Inchworm/bin/inchworm, Inchworm/bin/FastaToDeBruijn
  - Chrysalis/bin/{GraphFromFasta,BubbleUpClustering,CreateIwormFastaBundle,ReadsToTranscripts,QuantifyGraph}
  - Butterfly/Butterfly.jar
- 环境变量 `TRINITY_SRC` 可覆盖源码路径（默认上述绝对路径），xtask 使用
- 工具链版本:
  - cmake 3.28.3 (/usr/bin/cmake，系统自带)
  - g++ (conda-forge GCC) 13.4.0 — 系统只有 gcc 没有 g++，用 conda env `trinity-build`
    (`/public/home/senior007/miniconda3/envs/trinity-build/bin/x86_64-conda-linux-gnu-g++`)
  - java 17.0.8 (`module load apps/java/17.0.8`，仅运行 Butterfly 需要，编译不需要)
- 构建完成日期: 2026-08-15

## 重新构建（备忘）

本机（Ubuntu 24.04, HPC login 节点）没有系统 g++，也无 sudo，构建时需把 CC/CXX 指到 conda env：

```bash
export CC=/public/home/senior007/miniconda3/envs/trinity-build/bin/x86_64-conda-linux-gnu-gcc
export CXX=/public/home/senior007/miniconda3/envs/trinity-build/bin/x86_64-conda-linux-gnu-g++
cd /storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2
make inchworm_target chrysalis_target
```

env 创建命令（已执行，一次性）：

```bash
mamba create -y -n trinity-build -c conda-forge gxx_linux-64=13
```

注: conda g++ 链接时会把 RPATH 写进二进制（指向 env 的 libstdc++），产物在 env 存在的前提下可独立运行，无需激活 env。

## smoke 测试记录（inchworm）

输入（注意实际为 28bp，任务描述里写的 30bp 有误）：

```
>r1
ACGTACGTACGTACGTACGTACGTACGT
```

任务给的原样命令（`--monitor 1`，其余默认）**不输出任何 contig**：该 read 的 4 个 25-mer
count 均为 1，而默认 `min_seed_coverage=2`、`min_assembly_coverage=2`（见 Inchworm/src/IRKE.hpp
与 IRKE_run.cpp），种子全被 "insufficient coverage" 拒绝；进程 exit 0。

能产出 contig 的两种等价方式：

1. 放宽阈值（read 保持 1 条）：

```bash
Inchworm/bin/inchworm --reads tiny.fa --run_inchworm -K 25 --monitor 1 \
    --min_seed_coverage 1 --min_assembly_coverage 1
```

stdout：

```
>a1;1 total_counts: 4 Seed: 1 K: 25 length: 28
TACGTACGTACGTACGTACGTACGTACG
```

2. read 写两份（全部默认参数）：

```
>a1;2 total_counts: 8 Seed: 2 K: 25 length: 28
TACGTACGTACGTACGTACGTACGTACG
```

contig header 格式（后续黄金向量任务要对齐）：

```
>a<N>;<avg_cov> total_counts: <total_counts> Seed: <seed_kmer_count> K: <kmer_length> length: <seq_length>
```

- 序列 60 字符/行折行（`add_fasta_seq_line_breaks(seq, 60)`）
- 本例 contig 是 read 的循环移位（周期 4 序列的 de Bruijn 图是环，起始由被选中的种子决定，
  此处种子为 `TACGTACGTACGTACGTACGTACGT`），所以不是原 read 的前缀
- 进度信息全部走 stderr，stdout 只有 FASTA

## Butterfly.jar 裁定（P4 发现，2026-08-18 核实）

源码树里存在**两个语义不同的 Butterfly.jar**：

| jar | md5 | 来源 | 语义 |
|---|---|---|---|
| `Butterfly/Butterfly.jar`（发布版） | `312f253d3b2fcfe24d6e96c025b744ba` | 原始 tarball 自带，未改动（1,460,155 字节） | `TransAssembly_allProbPaths` 中 **combinePaths 两处调用**（`invokestatic` 偏移 108/487）均在执行；**新版 DFS_add_path_to_graph** |
| `Butterfly/Butterfly/Butterfly.jar`（源码树内层） | 工作树 `794d6f94daf036b6a4373f872e4cce28`（2026-08-17 本地构建）；tarball 原始版为 `379324fdb646a374e2ab89126ef6fa0d` | tarball **本身就有**这个内层 jar（与计划早期"tarball 内无此文件"的假设不符，已核实）；工作树版本是 P4 期间代理从源码树重新编译、**覆盖了 tarball 原始文件**的构建产物 | combinePaths **仅一处调用**（偏移 100，另一处调用在源码树中被注释）；内层 jar 与发布版在 combinePaths / DFS_add_path_to_graph 上均有行为差异 |

裁定：

- **发布版裁判 jar = `$TRINITY_SRC/Butterfly/Butterfly.jar`（312f…）**。所有黄金固化与
  `xcheck-butterfly` 对拍均以它为准——xtask 常量见 `xtask/src/main.rs::butterfly_jar()`
  （`trinity_src().join("Butterfly/Butterfly.jar")`）。
- 内层 jar（无论 tarball 原始 3793… 还是本地构建 794d…）**仅作源码树行为对照，不得再当 oracle**。
- 若需恢复 tarball 原始内层 jar：从 tarball 解出
  `trinityrnaseq-v2.15.2/Butterfly/Butterfly/Butterfly.jar` 覆盖即可（md5 应为 3793…）。

语义差异明细（发布版 vs 源码树/内层 jar，均经字节码与行为实证）：

1. **combinePaths**：`getSuffStats_wPairs` 中发布版 jar 实际执行 PairPath 合并
   （源码树该方法的相关调用被注释）；详见 `crates/trinity-butterfly/src/pair_paths.rs` 头注释的对拍记录。
2. **DFS_add_path_to_graph**：发布版为更新版实现（POG 回写图）；源码树为旧版；
   详见 `crates/trinity-butterfly/src/pog.rs` 头注释（c0 对拍据此调和）。
