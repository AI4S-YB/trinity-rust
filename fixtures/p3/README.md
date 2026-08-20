# fixtures/p3 — P3-T2 GraphFromFasta 对拍向量

sample_data（trinityrnaseq-Trinity-v2.15.2/sample_data/test_Trinity_Assembly）
真实数据的焊接聚类黄金（2026-08-17 生成，逐行命令见下——产物入库，原版
二进制不随测试运行）：

| 文件 | 来源 | 再生成 |
|---|---|---|
| `gff.reads.fa` | `zcat reads.left.fa.gz reads.right.fa.gz`（61,150 read） | 手跑 |
| `gff.iworm.fa` | 原版 `Inchworm/bin/inchworm --reads gff.reads.fa --run_inchworm -K 25 --DS --num_threads 4` 的 1287 contig 过 `filter_iworm_by_min_length_or_cov.pl iworm.fa 100 10`（219 条） | 原版二进制（inchworm 种子序非确定，**再生成会变**——本表为黄金向量，勿再生） |
| `gff.welds.orig.txt` | 原版 `Chrysalis/bin/GraphFromFasta -i gff.iworm.fa -r gff.reads.fa -min_contig_length 200 -min_glue 2 -glue_factor 0.05 -min_iso_ratio 0.05 -t 1 -k 24 -kk 48`（Trinity:2180 实参形态；非 SS 无 `-strand`；**-t 1** 消除 OMP toasted 竞态） | 原版二进制 |

对拍结论（`trinity-chrysalis/tests/gff_vs_original.rs`）：边多重集（(A,B) 有
向对 + weldmers/total/min_len 全字段）**完全相等**（24 边）；行序差异在预期
内（原版 report 用非稳定 sort + map 序，本库取确定性 (size, id, 插入序)，
行序非下游契约——BubbleUpClustering 消费前先过 `sort -k9,9gr`）。

P3-T3 追加（BubbleUpClustering / CreateIwormFastaBundle 对拍，同上 2026-08-17）：

| 文件 | 来源 | 再生成 |
|---|---|---|
| `bubble.orig.out` | `sort -k9,9gr gff.welds.orig.txt \| Chrysalis/bin/BubbleUpClustering -i gff.iworm.fa -weld_graph /dev/stdin -min_contig_length 200 -max_cluster_size 25` | 原版二进制（确定性，可再生成） |
| `bundle.orig.fa` | `Chrysalis/bin/CreateIwormFastaBundle -i bubble.orig.out -o bundle.orig.fa -min 200` | 原版二进制（确定性） |

对拍结论（`tests/bubble_vs_original.rs`）：COMPONENT 块多重集（按成员集合
分组，块内成员序/#POOL_INFO/折行一并比较）**完全相等**（55 组件）；
bundle 输出**逐字节相等**（110 行）。注意 `sort -k9,9gr` 无 `-s` 时对
key 相等的行做 last-resort 整行字节序比较——簇序（component 编号）依赖
该 tie-break，`sort_weld_graph` 已按此实现（否则本链 s_3/s_4 会互换）。

P3-T6 追加（QuantifyGraph + 分区器对拍，同上）：

| 文件 | 来源 | 再生成 |
|---|---|---|
| `quantify/c{0,1,2}/graph.tmp`、`reads.tmp` | `LC_ALL=C sort -k1,1n -k3,3nr -k2,2 rtt.orig.out`（**必须 C locale**——en_US.UTF-8 collation 会重排 key2）→ 原版 `partition_chrysalis_graphs_n_reads.pl --deBruijns f2db.orig.txt --componentReads sorted -N 1000 -L 200` | perl 脚本（确定性） |
| `quantify/c{0,1,2}/orig.graph.out`、`orig.reads.out` | 原版 `Chrysalis/bin/QuantifyGraph -i cN.reads.tmp -g cN.graph.tmp -o qg_out_N.graph -k 24`（读输出即 ReadsExt → `qg_out_N.reads`；注意无 `-no_cleanup` 时输入被删，先复制） | 原版二进制 |

对拍结论（`tests/quantify_vs_original.rs`）：分区产物（graph.tmp / reads.tmp）
与我们的 `partition`（先 `sort_reads_to_components`）**逐字节相等**（55 组件）；
QuantifyGraph graph.out **逐行相等，除"刀口边"仅第 3 列计数可差**（c0 8 行 /
c2 4 行 / c1 0 行），reads 输出过滤提及刀口边的行后逐行相等（c0 26/4343、
c2 5/7670）。刀口边 = 原版 `strncpy(&d[1], kmer_length)` 越界读
std::vector 堆尾巴一字节（上一任同尺寸缓冲的陈旧核苷酸字符），使熵判定的
第 25 字符在 {NUL,A,C,G,T} 间浮动——poly-A/GA 等熵 ≈ 1.0 的边跳过与否随堆
历史翻转（真 UB，不可精确复刻）；本库取确定语义（NUL 补位：24 计数字符 /
分母 25），非刀口边与原版完全一致。SortPrint 的 node2-seq 双 tab、单一
kmer 位置丢弃、ori=-1 组前置均已逐字节对拍锁定。
