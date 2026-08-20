# xcheck-trinity 报告（fast(50000)）

由 `cargo xtask xcheck-trinity` 生成。

# eval-trinity 报告

- 我们: `/storage/home/senior007/test/trinity_rust/trinity-rust/target/xcheck-trinity/out.Trinity.fasta`
- 原版: `/tmp/xcheck-trinity-orig/trinity_out.Trinity.fasta`

| 指标 | 我们 | 原版 |
|---|---|---|
| 转录本总数 | 89 | 80 |
| 总 bp | 143077 | 121309 |
| gene 数（gene_trans_map） | 59 | 60 |
| 双向覆盖率（精确+99% 聚类合并） | 62.9% | 70.0% |

- 精确匹配（rc 多重集交集）: **38** 条
- ≥99% 一致性聚类对: **18** 对（长度差 ≤10%、matches≥90% 较长链）
- 仅原版有: 42 条，长度 top10: [4351, 3969, 3918, 3739, 2817, 2812, 2792, 2758, 2746, 2683]
- 仅我们有: 51 条，长度 top10: [8314, 8254, 4262, 3966, 3915, 2789, 2783, 2762, 2730, 2712]

## both.fa 互喂抽查

- 归一化逐字节一致（记录多重集完全相等）

## SS(RF) 合成小集（5000 reads）

- 我们 4 条 / 原版 4 条；双向覆盖率 我们 75.0% / 原版 75.0%

## 判定

- 双向覆盖率 我们 62.9% / 原版 70.0%，阈值 50% → **PASS**
