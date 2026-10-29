# 评估 v2（junction 级双维度）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **提交冻结（2026-08-24 用户指示）:** 执行本计划期间**跳过所有 Commit 步骤**（本地也不提交），直到用户解除冻结。其余步骤照常。

**Goal:** xtask bench 新增 junction 级双维度指标（Jgene 完整度 / Jprec 准确度, StringTie GTF 参照, opt-in），real 轨道构建 StringTie 参照并回溯重评全部配置。

**Architecture:** 纯函数三件套（GTF→参考内含子链、PAF cg→组装本内含子链、判定规则）全部 TDD 手算用例锁定；bench 以 `--junc-ref`+`--genome` 旗标接入（不给旗标行为不变, 缓存命中同样计算）; real 参照用一次性幂等脚本镜像当年管线（HISAT2 → StringTie `-G` 引导模式）。

**Tech Stack:** xtask（Rust, 无新依赖）+ minimap2（已有 MINIMAP2 常量, `-x splice`）+ hisat2（rnaseq 环境）/ stringtie（stringtie-test 环境, 镜像当年 `-G` 调用）。

---

## 文件结构

```
xtask/src/junction.rs        # 新建: GTF 解析 / cg 提取 / 判定规则（纯函数 + 测试）
xtask/src/bench.rs           # 接线: --junc-ref/--genome + splice 调用 + 指标输出
xtask/src/main.rs            # +1 mod 声明（若 bench.rs 不是唯一使用处则并入 bench）
scripts/prepare_stringtie_real.sh  # 新建: real 轨道参照构建（幂等）
docs/opt-log.md              # 「评估 v2」段
```

---

### Task 1: GTF 解析 + 参考内含子链（纯函数, TDD）

**Files:** Create `xtask/src/junction.rs`; Modify `xtask/src/main.rs`（`mod junction;`）

- [x] **Step 1: 写失败测试（junction.rs 底部 tests）**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    /// 合成 GTF: 2 基因。
    /// STRG.1.1 = 3 外显子(chr1, +): [100,200],[301,400],[501,600]
    ///   → 内含子链 [(201,300),(401,500)]
    /// STRG.1.2 = 2 外显子(chr1, +): [100,200],[501,600] → 链 [(201,500)]
    /// STRG.2.1 = 单外显子(chr4, -): [1000,1500] → 空链
    /// exon 行故意乱序给 1.1（坐标排序后正确）。
    const GTF: &str = "# header\n\
chr1\tStringTie\ttranscript\t100\t600\t1000\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.1\";\n\
chr1\tStringTie\texon\t501\t600\t1000\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.1\"; exon_number \"3\";\n\
chr1\tStringTie\texon\t100\t200\t1000\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.1\"; exon_number \"1\";\n\
chr1\tStringTie\texon\t301\t400\t1000\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.1\"; exon_number \"2\";\n\
chr1\tStringTie\ttranscript\t100\t600\t1000\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.2\";\n\
chr1\tStringTie\texon\t100\t200\t1000\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.2\"; exon_number \"1\";\n\
chr1\tStringTie\texon\t501\t600\t1000\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.2\"; exon_number \"2\";\n\
chr4\tStringTie\ttranscript\t1000\t1500\t1000\t-\t.\tgene_id \"STRG.2\"; transcript_id \"STRG.2.1\";\n\
chr4\tStringTie\texon\t1000\t1500\t1000\t-\t.\tgene_id \"STRG.2\"; transcript_id \"STRG.2.1\"; exon_number \"1\";\n";

    #[test]
    fn parse_gtf_and_chains() {
        let model = parse_reference_gtf(GTF);
        assert_eq!(model.genes.len(), 2);
        let g1 = &model.genes["STRG.1"];
        assert_eq!(g1.len(), 2); // 两条转录本
        let tx11 = g1.iter().find(|t| t.tx == "STRG.1.1").unwrap();
        assert_eq!(tx11.chain, vec![(201, 300), (401, 500)]);
        assert_eq!(tx11.chrom, "chr1");
        assert_eq!(tx11.strand, '+');
        let tx12 = g1.iter().find(|t| t.tx == "STRG.1.2").unwrap();
        assert_eq!(tx12.chain, vec![(201, 500)]);
        let g2 = &model.genes["STRG.2"];
        assert!(g2[0].chain.is_empty(), "单外显子 → 空链");
        assert_eq!(g2[0].span(), (1000, 1500));
    }

    #[test]
    fn intron_set_dedupes_across_transcripts() {
        let model = parse_reference_gtf(GTF);
        let set = model.intron_set();
        assert!(set.contains(&(String::from("chr1"), '+', (201, 300))));
        assert!(set.contains(&(String::from("chr1"), '+', (401, 500))));
        assert!(set.contains(&(String::from("chr1"), '+', (201, 500))));
        assert_eq!(set.len(), 3);
    }
}
```

- [x] **Step 2: 跑测试确认失败**

Run: `cargo test -p xtask junction`
Expected: FAIL（模块不存在）

- [x] **Step 3: 实现**

```rust
//! 评估 v2（junction 级）: StringTie GTF 参考模型 + 组装本内含子链 + 判定规则。
//! 全部纯函数, 手算用例锁定。

use std::collections::{HashMap, HashSet};

/// 一条参考转录本的外显子模型（坐标排序后）。
#[derive(Debug, Clone)]
pub struct RefTx {
    pub tx: String,
    pub chrom: String,
    pub strand: char,
    pub exons: Vec<(u64, u64)>, // 升序, 1-based inclusive（GTF 原语义）
    pub chain: Vec<(u64, u64)>, // 内含子 = 相邻外显子间隔
}

impl RefTx {
    /// 基因组跨度（单外显子回退准则用）。
    pub fn span(&self) -> (u64, u64) {
        (self.exons.first().map(|e| e.0).unwrap_or(0),
         self.exons.last().map(|e| e.1).unwrap_or(0))
    }
}

/// 参考模型: gene_id → 转录本列表 + 全体内含子集合。
pub struct RefModel {
    pub genes: HashMap<String, Vec<RefTx>>,
}

impl RefModel {
    /// (chrom, strand, intron) 全集合（novel 判定用）。
    pub fn intron_set(&self) -> HashSet<(String, char, (u64, u64))> {
        let mut s = HashSet::new();
        for txs in self.genes.values() {
            for tx in txs {
                for &iv in &tx.chain {
                    s.insert((tx.chrom.clone(), tx.strand, iv));
                }
            }
        }
        s
    }

    /// 多外显子参考链集合（chrom, strand, 链）——Jtx/Jgene 判定用。
    pub fn chain_set(&self) -> HashSet<(String, char, Vec<(u64, u64)>)> {
        let mut s = HashSet::new();
        for txs in self.genes.values() {
            for tx in txs {
                if !tx.chain.is_empty() {
                    s.insert((tx.chrom.clone(), tx.strand, tx.chain.clone()));
                }
            }
        }
        s
    }
}

/// GTF 解析: 只取 exon 行, 按 transcript_id 聚合; gene_id 从属性列提取。
/// exon 坐标行内升序排序（GTF 天然有序, 乱序防御性排序）。
pub fn parse_reference_gtf(text: &str) -> RefModel {
    let mut txs: HashMap<String, RefTx> = HashMap::new();
    for line in text.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 || f[2] != "exon" {
            continue;
        }
        let (Ok(start), Ok(end)) = (f[3].parse::<u64>(), f[4].parse::<u64>()) else { continue };
        let gene = attr(f[8], "gene_id");
        let txid = attr(f[8], "transcript_id");
        let (Some(gene), Some(txid)) = (gene, txid) else { continue };
        let e = txs.entry(txid.clone()).or_insert_with(|| RefTx {
            tx: txid,
            chrom: f[0].to_string(),
            strand: f[6].chars().next().unwrap_or('.'),
            exons: Vec::new(),
            chain: Vec::new(),
        });
        e.exons.push((start, end));
    }
    let mut genes: HashMap<String, Vec<RefTx>> = HashMap::new();
    // 需要 gene 归属: 重新扫 transcript 行拿 gene_id↔tx 映射（exon 行也带 gene_id,
    // 但上面 entry 只按 txid——补: 解析时先记 tx→gene）
    //（实现注: 上面 attr 已取 gene——改为 entry 时一并存。此处示意, 实现按此逻辑完成。）
    for (txid, mut tx) in txs {
        tx.exons.sort_unstable();
        tx.chain = tx.exons.windows(2).map(|w| (w[0].1 + 1, w[1].0 - 1)).collect();
        genes.entry(tx_gene_map(&txid)).or_default().push(tx);
    }
    RefModel { genes }
}
```

（实现注: `tx_gene_map` 是上面推导的示意残留——**正确形态**是 `txs` 的 entry 用
`(txid, gene)` 二元组存储, 结束时按 gene 分组, 不需要二次查找。属性提取:

```rust
fn attr<'a>(field: &'a str, key: &str) -> Option<String> {
    let pat = format!("{key} \"");
    field.split(';').find_map(|seg| {
        let seg = seg.trim();
        seg.strip_prefix(&pat)
            .map(|rest| rest.trim_end_matches('"').to_string())
    })
}
```

）

- [x] **Step 4: 跑测试确认通过**

Run: `cargo test -p xtask junction`
Expected: PASS (2)

- [x] **Step 5: clippy + fmt + Commit（冻结跳过 commit）**

---

### Task 2: PAF cg → 组装本内含子链（纯函数, TDD）

**Files:** Modify `xtask/src/junction.rs`

- [x] **Step 1: 写失败测试（tests 模块追加）**

```rust
    /// cg 遍历: tstart=1000(0-based), cg=100M200N50M30I20M
    /// → 内含子 (1101, 1300)（1-based inclusive: 1000+100=1100 位后起, N 长 200）
    #[test]
    fn introns_from_cigar() {
        let chain = introns_from_cigar("100M200N50M30I20M", 1000);
        assert_eq!(chain, vec![(1101, 1300)]);
        let none = introns_from_cigar("300M", 1000);
        assert!(none.is_empty());
        // D 推进参考坐标但不产生内含子
        let d = introns_from_cigar("100M10D100M", 0);
        assert!(d.is_empty());
    }

    /// PAF 行解析: qid/链方向/tname/链提取。
    #[test]
    fn parse_splice_paf_lines() {
        let paf = "a1\t600\t0\t599\t+\tchr1\t100000\t900\t1499\t540\t560\t60\ttp:A:P\tcg:Z:300M200N300M\n\
                   a2\t300\t0\t299\t-\tchr4\t2000\t2000\t2299\t280\t300\t60\ttp:A:P\tcg:Z:300M\n";
        let asm = parse_splice_paf(paf);
        assert_eq!(asm.len(), 2);
        assert_eq!(asm[0].qid, "a1");
        assert_eq!(asm[0].chrom, "chr1");
        assert_eq!(asm[0].strand, '+');
        assert_eq!(asm[0].chain, vec![(1201, 1400)]);
        // ident = mch/aln = 540/560
        assert!((asm[0].ident - 540.0 / 560.0).abs() < 1e-9);
        assert!(asm[1].chain.is_empty(), "单外显子组装本");
    }
```

- [x] **Step 2: 跑测试确认失败**

Run: `cargo test -p xtask junction`
Expected: FAIL

- [x] **Step 3: 实现**

```rust
/// 一条组装本的主比对（splice PAF 主线）。
pub struct AsmTx {
    pub qid: String,
    pub chrom: String,
    pub strand: char,
    pub chain: Vec<(u64, u64)>, // 内含子（1-based inclusive, 基因组坐标）
    pub ident: f64,             // mch/aln（单外显子回退准则用）
    pub blocks: Vec<(u64, u64)>,// 比对块（基因组坐标, 覆盖度用）
}

/// cg 标签 → 内含子链。tstart 为 PAF 第 8 列（0-based）。
/// 内含子 1-based inclusive: [tstart+consumed+1, tstart+consumed+N]。
pub fn introns_from_cigar(cg: &str, tstart: u64) -> Vec<(u64, u64)> {
    let mut out = Vec::new();
    let bytes = cg.as_bytes();
    let mut i = 0;
    let mut ref_pos = tstart; // 0-based
    while i < bytes.len() {
        let num_start = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() { i += 1; }
        if i == num_start || i >= bytes.len() { break; }
        let n: u64 = cg[num_start..i].parse().unwrap_or(0);
        let op = bytes[i] as char;
        i += 1;
        match op {
            'N' => { out.push((ref_pos + 1, ref_pos + n)); ref_pos += n; }
            'M' | 'D' | '=' | 'X' => ref_pos += n,
            _ => {} // I/S/H 不推进参考
        }
    }
    out
}

/// splice PAF 全文 → 每查询一条主比对（tp:A:P 或首行）; cg 缺失的行跳过。
pub fn parse_splice_paf(text: &str) -> Vec<AsmTx> {
    let mut best: HashMap<String, AsmTx> = HashMap::new();
    let mut order: Vec<String> = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 { continue; }
        let cg = f.iter().find_map(|t| t.strip_prefix("cg:Z:"));
        let Some(cg) = cg else { continue };
        let qid = f[0].to_string();
        let tstart: u64 = f[7].parse().unwrap_or(0);
        let mch: u64 = f[9].parse().unwrap_or(0);
        let aln: u64 = f[10].parse().unwrap_or(1);
        let aln_len = f[10].parse::<u64>().unwrap_or(0); // 以 query 对齐长度论主次
        let better = match best.get(&qid) {
            Some(b) => aln_len > b.blocks.iter().map(|_| 0).count() as u64 + b.chain_len(),
            None => true,
        };
        //（实现注: 主次判定以对齐总长（query 侧 qe-qs）为准——直接用 f[3]-f[2];
        // 上面示意残留不采用, 按此实现）
        let qs: u64 = f[2].parse().unwrap_or(0);
        let qe: u64 = f[3].parse().unwrap_or(0);
        let qlen_hit = qe.saturating_sub(qs);
        let take = match best.get(&qid) {
            Some(b) => qlen_hit > b.qlen_hit,
            None => true,
        };
        if !take { continue; }
        if !best.contains_key(&qid) { order.push(qid.clone()); }
        let chain = introns_from_cigar(cg, tstart);
        // blocks: 由 cg 重建 M 块（基因组坐标）
        let mut blocks = Vec::new();
        let mut pos = tstart;
        let bytes = cg.as_bytes();
        let mut i = 0;
        while i < bytes.len() {
            let ns = i;
            while i < bytes.len() && bytes[i].is_ascii_digit() { i += 1; }
            if i == ns || i >= bytes.len() { break; }
            let n: u64 = cg[ns..i].parse().unwrap_or(0);
            match bytes[i] as char {
                'M' | '=' | 'X' => { blocks.push((pos + 1, pos + n)); pos += n; }
                'N' | 'D' => pos += n,
                _ => {}
            }
            i += 1;
        }
        best.insert(qid.clone(), AsmTx {
            qid,
            chrom: f[5].to_string(),
            strand: f[4].chars().next().unwrap_or('.'),
            chain,
            ident: mch as f64 / aln.max(1) as f64,
            blocks,
            qlen_hit,
        });
    }
    order.into_iter().filter_map(|q| best.remove(&q)).collect()
}
```

（`AsmTx` 需加 `qlen_hit: u64` 字段; `chain_len()` 残留删除。）

- [x] **Step 4: 跑测试确认通过**

Run: `cargo test -p xtask junction`
Expected: PASS (4)

- [x] **Step 5: clippy + fmt + Commit（冻结跳过）**

---

### Task 3: 判定规则（纯函数, TDD）

**Files:** Modify `xtask/src/junction.rs`

- [x] **Step 1: 写失败测试（追加）**

```rust
    /// 手算判定用例（复用 Task1 GTF 模型）:
    /// 组装本:
    ///   A1 = chr1,+ 链 [(201,300),(401,500)] → 精确复现 STRG.1.1 ✓
    ///   A2 = chr1,+ 链 [(201,999)] → novel 内含子 → 不正确
    ///   A3 = chr4,- 无链（单外显子, blocks 覆盖 STRG.2.1 的 100%）
    ///   A4 = chr1,+ 无链（单外显子, 不覆盖任何）
    /// 期望: Jprec = 1/2; Jgene = 2/2（STRG.1 经 A1 找回, STRG.2 经 A3 回退找回）;
    /// Jtx = 1/3; novelJ = 1/2; single = 2/4
    #[test]
    fn eval_junction_hand_case() {
        let model = parse_reference_gtf(GTF);
        let asm = vec![
            AsmTx { qid: "A1".into(), chrom: "chr1".into(), strand: '+',
                    chain: vec![(201, 300), (401, 500)], ident: 0.99,
                    blocks: vec![(100, 600)], qlen_hit: 500 },
            AsmTx { qid: "A2".into(), chrom: "chr1".into(), strand: '+',
                    chain: vec![(201, 999)], ident: 0.99,
                    blocks: vec![(100, 1000)], qlen_hit: 900 },
            AsmTx { qid: "A3".into(), chrom: "chr4".into(), strand: '-',
                    chain: vec![], ident: 0.99,
                    blocks: vec![(1000, 1500)], qlen_hit: 501 },
            AsmTx { qid: "A4".into(), chrom: "chr1".into(), strand: '+',
                    chain: vec![], ident: 0.99,
                    blocks: vec![(5000, 5200)], qlen_hit: 201 },
        ];
        let m = eval_junction(&asm, &model);
        assert!((m.j_prec - 0.5).abs() < 1e-9);
        assert!((m.j_gene - 1.0).abs() < 1e-9);
        assert!((m.j_tx - 1.0 / 3.0).abs() < 1e-9);
        assert!((m.novel_j - 0.5).abs() < 1e-9);
        assert!((m.single_exon - 0.5).abs() < 1e-9);
    }
```

- [x] **Step 2: 跑测试确认失败**

Run: `cargo test -p xtask junction`
Expected: FAIL

- [x] **Step 3: 实现**

```rust
/// junction 双维度指标。
pub struct JunctionMetrics {
    pub j_gene: f64,   // 找回基因 / 基因总数
    pub j_tx: f64,     // 精确复现参考转录本 / 参考转录本总数
    pub j_prec: f64,   // junction-correct 多外显子组装本 / 多外显子组装本总数
    pub novel_j: f64,  // 含 novel 内含子的多外显子组装本比例
    pub single_exon: f64, // 单外显子组装本比例
    pub n_gene: usize,
    pub n_tx: usize,
    pub n_asm: usize,
}

/// 判定（全部规则见 spec）:
/// - 正确: 组装本链 == 某参考链（chrom+strand+内含子集合完全一致）
/// - 基因找回: 任一多外显子参考链被复现, 或任一单外显子参考转录本被
///   覆盖 ≥70%（同 chrom 同链组装本 blocks 并集）且该组装本 ident ≥95%
/// - Jtx: 精确复现的参考转录本 / 总数（单外显子参考按回退找回计）
pub fn eval_junction(asm: &[AsmTx], model: &RefModel) -> JunctionMetrics {
    let chains = model.chain_set();
    let introns = model.intron_set();

    let multi: Vec<&AsmTx> = asm.iter().filter(|a| !a.chain.is_empty()).collect();
    let single = asm.len().saturating_sub(multi.len());
    let correct = multi.iter()
        .filter(|a| chains.contains(&(a.chrom.clone(), a.strand, a.chain.clone())))
        .count();
    let novel = multi.iter()
        .filter(|a| a.chain.iter().any(|iv| !introns.contains(&(a.chrom.clone(), a.strand, *iv))))
        .count();

    // 参考侧: 逐基因判定
    let mut genes_hit = 0usize;
    let mut tx_hit = 0usize;
    let n_tx_total: usize = model.genes.values().map(|v| v.len()).sum();
    for txs in model.genes.values() {
        let mut gene_hit = false;
        for tx in txs {
            let hit = if !tx.chain.is_empty() {
                asm.iter().any(|a| a.chain == tx.chain && a.chrom == tx.chrom && a.strand == tx.strand)
            } else {
                // 单外显子回退: 覆盖 ≥70% 且组装本 ident ≥0.95
                let (s, e) = tx.span();
                let span_len = e.saturating_sub(s) + 1;
                asm.iter().any(|a| {
                    a.ident >= 0.95 && a.chrom == tx.chrom && a.strand == tx.strand && {
                        let cov: u64 = a.blocks.iter()
                            .map(|&(bs, be)| be.min(e).saturating_sub(bs.max(s)).saturating_add(1))
                            .sum();
                        cov * 10 >= span_len * 7
                    }
                })
            };
            if hit { tx_hit += 1; gene_hit = true; }
        }
        if gene_hit { genes_hit += 1; }
    }

    let n = asm.len().max(1);
    JunctionMetrics {
        j_gene: genes_hit as f64 / model.genes.len().max(1) as f64,
        j_tx: tx_hit as f64 / n_tx_total.max(1) as f64,
        j_prec: correct as f64 / multi.len().max(1) as f64,
        novel_j: novel as f64 / multi.len().max(1) as f64,
        single_exon: single as f64 / n as f64,
        n_gene: model.genes.len(),
        n_tx: n_tx_total,
        n_asm: asm.len(),
    }
}
```

- [x] **Step 4: 跑测试确认通过**

Run: `cargo test -p xtask junction`
Expected: PASS (5)

- [x] **Step 5: clippy + fmt + Commit（冻结跳过）**

---

### Task 4: bench 接线（--junc-ref/--genome + splice 调用 + 输出）

**Files:** Modify `xtask/src/bench.rs`

- [x] **Step 1: 旗标解析（仿 --prot-ref 模式）**

```rust
    let mut junc_ref = String::new();
    let mut genome = String::new();
    // match 分支:
            "--junc-ref" => { junc_ref = argv.get(i + 1).cloned().unwrap_or_default(); i += 2; }
            "--genome" => { genome = argv.get(i + 1).cloned().unwrap_or_default(); i += 2; }
```

校验: 一方给另一方缺 → 用法报错（返回 2）。

- [x] **Step 2: 计算与输出（metrics.json 写入之前插入）**

```rust
    // 评估 v2: junction 双维度（opt-in: --junc-ref + --genome）
    if !junc_ref.is_empty() && !genome.is_empty() {
        let model = junction::parse_reference_gtf(
            &std::fs::read_to_string(&junc_ref).unwrap_or_else(|e| panic!("read {junc_ref}: {e}")));
        let paf = std::process::Command::new(MINIMAP2)
            .args(["-x", "splice", "-t", "8", "--secondary=no"])
            .arg(&genome).arg(&final_fa)
            .output()
            .expect("无法启动 minimap2");
        if !paf.status.success() {
            eprintln!("minimap2 splice 失败: {}", String::from_utf8_lossy(&paf.stderr));
            return 1;
        }
        let asm = junction::parse_splice_paf(&String::from_utf8_lossy(&paf.stdout));
        let m = junction::eval_junction(&asm, &model);
        println!("  Junction: Jgene={:.1}% Jtx={:.1}% Jprec={:.1}% novelJ={:.1}% single={:.1}% (genes={} tx={} asm={})",
            100.0 * m.j_gene, 100.0 * m.j_tx, 100.0 * m.j_prec, 100.0 * m.novel_j,
            100.0 * m.single_exon, m.n_gene, m.n_tx, m.n_asm);
        // metrics.json 行追加 ,"j_gene":{:.4},"j_tx":{:.4},"j_prec":{:.4},"novel_j":{:.4},"single_exon":{:.4}
        //（对应实参 m.j_gene 等——并把这些值先存到外层变量供两处使用）
    }
```

（metrics.json 的 format 串与实参列表同步扩展; stdout 行放在既有汇总行之后。bench.rs 内引用用 `use crate::junction;` 或 `crate::junction::` 前缀——junction 模块在 main.rs 声明为 `mod junction;`, 若 bench.rs 是独立编译单元则按 xtask 现有模块引用惯例调整。）

- [x] **Step 3: dev 缓存重放 sanity（缓存命中, 秒级）**

Run: `A=/storage/home/senior007/test/athal && cargo run -q --release -p xtask -- bench --tag evalv2-dev1-m5 --reads $A/bench/reads/dev1_1.fq.gz,$A/bench/reads/dev1_2.fq.gz --ref $A/bench/refs/dev1.fa -- --path-engine pogflow --model-guided --max_memory 16G --CPU 8 --junc-ref $A/stringtie.gtf --genome $A/athal5.fa`

（junc-ref/genome 旗标不进 passthrough——它们是 bench 的; 传参位置在 `--` 之前。若 bench 参数序要求在 -- 之前则如上。）

Expected: `[bench] 组装缓存命中` + `Junction:` 行; **Jgene 应 ≥30%**（dev 参考即 stringtie 子集, 高找回是 sanity 底线; 若 <10% 说明链提取/坐标系有 bug, 停下检查）。注意 dev 组装参考 fasta 是 stringtie.fa 子集而 GTF 是全量——Jgene 分母含未纳入 dev 子集的基因, 数字低于 100% 是预期的。

- [x] **Step 4: clippy + fmt + 全量测试 + Commit（冻结跳过）**

---

### Task 5: real 参照构建 + 回溯重评 + 验收入档

**Files:** Create `scripts/prepare_stringtie_real.sh`; Modify `docs/opt-log.md`

- [x] **Step 1: 读当年管线, 镜像调用**

Run: `cat /storage/home/senior007/test/athal/run_orig.sh`（hd_pipeline.sh 已知: `stringtie -G chr5.renamed.gtf -p 8 -o ... -A ...`）——记下全基因组版的实际参数（引导 GTF 的生成方式 + 线程数）。**若 run_orig.sh 与 hd 版不同, 以全基因组版为准并在脚本注释里注明来源行号。**

- [x] **Step 2: 写 `scripts/prepare_stringtie_real.sh`（幂等）**

```bash
#!/usr/bin/env bash
# 评估 v2 real 参照构建（幂等）。产物: bench/real/{hisat2.bam,real_stringtie.gtf}
set -euo pipefail
ATHAL=/storage/home/senior007/test/athal
R=$ATHAL/bench/real
H2=/public/home/senior007/miniconda3/envs/rnaseq/bin/hisat2
ST=/public/home/senior007/miniconda3/envs/stringtie-test/bin/stringtie
SAM=/public/home/senior007/miniconda3/envs/stringtie-test/bin/samtools
mkdir -p "$R"

# 引导 GTF: 镜像当年 -G（若当年用 chr5.renamed.gtf 的全基因组版不存在,
# 用 gffread 由 TAIR10.gff.gz 转换并改染色体名对齐 athal5.fa——先对比
# athal5.fa 与 TAIR10.fna.gz 的染色体命名再定）
GUIDE=$R/guide.gtf
if [ ! -f "$GUIDE" ]; then
  #（按 Step 1 的实际发现实现: 优先复用当年已有文件; 否则 gffread 转换+改名）
  :
fi

if [ ! -f "$R/hisat2.bam" ]; then
  $H2 -x $ATHAL/athal5 -p 8 -1 $R/R1.fq.gz -2 $R/R2.fq.gz 2> $R/hisat2.log \
    | $SAM sort -@ 4 -o $R/hisat2.bam -
fi
if [ ! -f "$R/real_stringtie.gtf" ]; then
  $ST -G "$GUIDE" -p 8 -o $R/real_stringtie.gtf -A $R/real_stringtie.abund $R/hisat2.bam
fi
echo "genes: $(grep -c $'\ttranscript\t' $R/real_stringtie.gtf) transcripts"
```

（`:` 占位的 GUIDE 段必须按 Step 1 的实际发现填成可执行代码——这是本 Task 的核心判断点, 不许留空。）

- [x] **Step 3: 执行（后台, HISAT2 ~10-20 min + StringTie ~5-10 min）**

Run: `chmod +x scripts/prepare_stringtie_real.sh && ./scripts/prepare_stringtie_real.sh`
Expected: 转录本数量级 ~3-5 万（与 dev stringtie.gtf 同量级）; 单 read sanity: `$SAM view hisat2.bam | head -1` 有比对。

- [x] **Step 4: 回溯重评（全部缓存命中, 秒级/次）**

```bash
A=/storage/home/senior007/test/athal/bench/real
R=$A/R1.fq.gz,$A/R2.fq.gz; REF=$A/atha_cdna.fa; JR=$A/real_stringtie.gtf; GM=$ATHAL/athal5.fa
cd /storage/home/senior007/test/trinity_rust/trinity-rust
cargo run -q --release -p xtask -- bench --tag v2-real-classic  --reads $R --ref $REF --junc-ref $JR --genome $GM -- --max_memory 16G --CPU 8
cargo run -q --release -p xtask -- bench --tag v2-real-pogflow --reads $R --ref $REF --junc-ref $JR --genome $GM -- --path-engine pogflow --max_memory 16G --CPU 8
cargo run -q --release -p xtask -- bench --tag v2-real-m5      --reads $R --ref $REF --junc-ref $JR --genome $GM -- --path-engine pogflow --model-guided --max_memory 16G --CPU 8
cargo run -q --release -p xtask -- bench --tag v2-real-m6a     --reads $R --ref $REF --junc-ref $JR --genome $GM -- --path-engine pogflow --protein-guided --protein-db $A/alyr_pep.fa --max_memory 16G --CPU 8
```

注意: bench 缓存 key 只含 passthrough+reads+binary → 四 tag 的组装应缓存命中（若 miss 则照常跑, 记录 wall）。

- [x] **Step 5: 验收检查 + 入档**

1. **区分力**: v2-real-m5 的 Jgene > v2-real-classic 的 Jgene（F1 体系同向）
2. **稳定性**: 任一 tag 重放第二次, Junction 行逐字一致
3. opt-log「评估 v2」段: 新旧指标同表（R/P/F1/P_full + Jgene/Jtx/Jprec/novelJ/single）四配置 + evalv2-dev1-m5 一行; 背离点注明; StringTie 参照构建命令与产出数入档
4. 向用户汇报双维度全景

- [x] **Step 6: Commit（冻结跳过）**

---


> **执行记录（2026-08-26）:** 三处证据驱动的修正——① minimap2 splice 需 `-c` 才出 cg 标签
> （否则 0 条链）; ② 链匹配改为**链方向无关**（非链特异性组装本方向任意, 链标记匹配丢弃
> 约半数真实命中, dev1 实测 Jprec 10.0→21.0%）; ③ Jtx 语义锁定为"多外显子链精确复现"
> （单外显子回退只计 Jgene）。T5 考古发现 dev 参照本就是 chr5 引导+全基因组重比对产物,
> real 镜像同源流程。dev sanity Jgene=30.4%（越 30% 线）; real 四配置区分力/稳定性双 PASS。

## 完成定义

- [x] junction.rs 三件套纯函数 + 5 个手算单测全绿
- [x] bench opt-in 旗标工作, 不给旗标行为不变
- [x] dev sanity: Jgene ≥30%
- [x] real 参照构建完成且入档
- [x] 四配置 + dev 的 junction 对照表入档, 区分力/稳定性两项验收过
- [x] 全程无 git 提交（冻结）
