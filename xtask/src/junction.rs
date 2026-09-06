//! 评估 v2（junction 级）: StringTie GTF 参考模型 + 组装本内含子链 + 判定规则。
//! 全部纯函数, 手算用例锁定。坐标 1-based inclusive（PAF tstart 0-based 在
//! introns_from_cigar 内转换）。
//! 匹配为链方向无关（strand-agnostic）: de novo 非链特异性组装本方向任意,
//! 链方向不参与匹配（证据: dev1 sanity 中链标记匹配丢弃约半数真实命中,
//! Jprec 10.0%→21.0%, 基因链命中 519→1062）。

use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct RefTx {
    pub tx: String,
    pub chrom: String,
    /// 保留链字段仅为解析保真与测试断言; 匹配不使用（strand-agnostic）。
    #[allow(dead_code)]
    pub strand: char,
    pub exons: Vec<(u64, u64)>,
    pub chain: Vec<(u64, u64)>,
}

impl RefTx {
    pub fn span(&self) -> (u64, u64) {
        (
            self.exons.first().map(|e| e.0).unwrap_or(0),
            self.exons.last().map(|e| e.1).unwrap_or(0),
        )
    }
}

/// 链键: (chrom, 完整内含子链)。链方向不参与匹配（见模块注释）。
pub type ChainKey = (String, Vec<(u64, u64)>);

pub struct RefModel {
    pub genes: HashMap<String, Vec<RefTx>>,
}

impl RefModel {
    pub fn intron_set(&self) -> HashSet<(String, (u64, u64))> {
        let mut set = HashSet::new();
        for txs in self.genes.values() {
            for tx in txs {
                for iv in &tx.chain {
                    set.insert((tx.chrom.clone(), *iv));
                }
            }
        }
        set
    }

    /// 多外显子转录本的内含子链集合（chrom+完整链, 链方向不参与）。
    pub fn chain_set(&self) -> HashSet<ChainKey> {
        let mut set = HashSet::new();
        for txs in self.genes.values() {
            for tx in txs {
                if !tx.chain.is_empty() {
                    set.insert((tx.chrom.clone(), tx.chain.clone()));
                }
            }
        }
        set
    }
}

/// GTF 属性列取值: `key "value";` 格式。
fn attr(field: &str, key: &str) -> Option<String> {
    let pat = format!("{key} \"");
    field.split(';').find_map(|seg| {
        let seg = seg.trim();
        seg.strip_prefix(&pat)
            .map(|rest| rest.trim_end_matches('"').to_string())
    })
}

/// 解析 StringTie GTF: 只取 exon 行, 按 transcript_id 聚合, 排序后由相邻外显子
/// 间隔得内含子链, 按 gene_id 分组。
pub fn parse_reference_gtf(text: &str) -> RefModel {
    // transcript_id -> (gene_id, RefTx)
    let mut txs: HashMap<String, (String, RefTx)> = HashMap::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 || f[2] != "exon" {
            continue;
        }
        let (Some(gene), Some(tx)) = (attr(f[8], "gene_id"), attr(f[8], "transcript_id")) else {
            continue;
        };
        let (Ok(start), Ok(end)) = (f[3].parse::<u64>(), f[4].parse::<u64>()) else {
            continue;
        };
        let entry = txs.entry(tx).or_insert_with(|| {
            (
                gene.clone(),
                RefTx {
                    tx: String::new(),
                    chrom: f[0].to_string(),
                    strand: f[6].chars().next().unwrap_or('.'),
                    exons: Vec::new(),
                    chain: Vec::new(),
                },
            )
        });
        entry.1.exons.push((start, end));
    }

    let mut genes: HashMap<String, Vec<RefTx>> = HashMap::new();
    for (tx_id, (gene, mut r)) in txs {
        r.tx = tx_id;
        r.exons.sort_unstable();
        r.chain = r
            .exons
            .windows(2)
            .map(|w| (w[0].1 + 1, w[1].0 - 1))
            .collect();
        genes.entry(gene).or_default().push(r);
    }
    RefModel { genes }
}

pub struct AsmTx {
    pub qid: String,
    pub chrom: String,
    /// 保留链字段仅为解析保真与测试断言; 匹配不使用（strand-agnostic）。
    #[allow(dead_code)]
    pub strand: char,
    pub chain: Vec<(u64, u64)>,
    pub ident: f64,              // mch/aln
    pub blocks: Vec<(u64, u64)>, // 比对块（基因组 1-based inclusive）
    pub qlen_hit: u64,           // query 对齐长度（主次判定）
}

/// CIGAR 切成 (op, len) 序列; 遇畸形前缀即止。
fn cigar_ops(cg: &str) -> Vec<(char, u64)> {
    let b = cg.as_bytes();
    let mut ops = Vec::new();
    let mut i = 0;
    while i < b.len() {
        let num_start = i;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
        if i == num_start || i >= b.len() {
            break;
        }
        let n: u64 = cg[num_start..i].parse().unwrap_or(0);
        ops.push((b[i] as char, n));
        i += 1;
    }
    ops
}

/// cg:Z:CIGAR + tstart（0-based, PAF f[7]）→ 内含子链（1-based inclusive）。
/// N 记内含子; M/=/X/D 推进参考游标; I/S/H 不推进。
pub fn introns_from_cigar(cg: &str, tstart: u64) -> Vec<(u64, u64)> {
    let mut out = Vec::new();
    let mut pos = tstart; // 0-based 参考游标
    for (op, n) in cigar_ops(cg) {
        match op {
            'N' => {
                out.push((pos + 1, pos + n));
                pos += n;
            }
            'M' | '=' | 'X' | 'D' => pos += n,
            _ => {} // I/S/H 及未知 op 不推进参考
        }
    }
    out
}

/// cg → 比对块（M/=/X 块, 1-based inclusive）。
fn blocks_from_cigar(cg: &str, tstart: u64) -> Vec<(u64, u64)> {
    let mut out = Vec::new();
    let mut pos = tstart; // 0-based
    for (op, n) in cigar_ops(cg) {
        match op {
            'M' | '=' | 'X' => {
                out.push((pos + 1, pos + n));
                pos += n;
            }
            'N' | 'D' => pos += n,
            _ => {}
        }
    }
    out
}

/// minimap2 PAF → 组装本列表。每查询取主比对（tp:A:P 行, 无则 qlen_hit 最大者）;
/// 无 cg:Z: 的行跳过。
pub fn parse_splice_paf(text: &str) -> Vec<AsmTx> {
    // (AsmTx, 是否 tp:A:P)
    let mut hits: Vec<Option<(AsmTx, bool)>> = Vec::new();
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        let tags = &f[12..];
        let Some(cg) = tags.iter().find_map(|t| t.strip_prefix("cg:Z:")) else {
            continue; // 无 CIGAR 的行跳过
        };
        let (Ok(qs), Ok(qe), Ok(tstart), Ok(mch), Ok(aln)) = (
            f[2].parse::<u64>(),
            f[3].parse::<u64>(),
            f[7].parse::<u64>(),
            f[9].parse::<u64>(),
            f[10].parse::<u64>(),
        ) else {
            continue;
        };
        hits.push(Some((
            AsmTx {
                qid: f[0].to_string(),
                chrom: f[5].to_string(),
                strand: f[4].chars().next().unwrap_or('.'),
                chain: introns_from_cigar(cg, tstart),
                ident: if aln > 0 {
                    mch as f64 / aln as f64
                } else {
                    0.0
                },
                blocks: blocks_from_cigar(cg, tstart),
                qlen_hit: qe.saturating_sub(qs),
            },
            tags.contains(&"tp:A:P"),
        )));
    }

    // 按查询分组（保持首次出现顺序）, 每组选主比对
    let mut first: HashMap<String, usize> = HashMap::new();
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for (i, h) in hits.iter().enumerate() {
        let g = *first
            .entry(h.as_ref().unwrap().0.qid.clone())
            .or_insert_with(|| {
                groups.push(Vec::new());
                groups.len() - 1
            });
        groups[g].push(i);
    }

    let mut out = Vec::with_capacity(groups.len());
    for cand in groups {
        let primary: Vec<usize> = cand
            .iter()
            .copied()
            .filter(|&i| hits[i].as_ref().unwrap().1)
            .collect();
        let pool: &[usize] = if primary.is_empty() { &cand } else { &primary };
        // qlen_hit 最大者; 平局取先出现
        let mut pick = pool[0];
        for &i in &pool[1..] {
            if hits[i].as_ref().unwrap().0.qlen_hit > hits[pick].as_ref().unwrap().0.qlen_hit {
                pick = i;
            }
        }
        out.push(hits[pick].take().unwrap().0);
    }
    out
}

pub struct JunctionMetrics {
    pub j_gene: f64,
    pub j_tx: f64,
    pub j_prec: f64,
    pub novel_j: f64,
    pub single_exon: f64,
    pub n_gene: usize,
    pub n_tx: usize,
    pub n_asm: usize,
}

/// blocks 并集与 span 交集的覆盖长度（区间裁剪 → 排序 → 合并求和）。
fn union_covered(blocks: &[(u64, u64)], span: (u64, u64)) -> u64 {
    let mut segs: Vec<(u64, u64)> = blocks
        .iter()
        .map(|&(s, e)| (s.max(span.0), e.min(span.1)))
        .filter(|&(s, e)| e >= s)
        .collect();
    segs.sort_unstable();
    let mut total = 0u64;
    let mut cur: Option<(u64, u64)> = None;
    for (s, e) in segs {
        match cur {
            Some((cs, ce)) if s <= ce.saturating_add(1) => {
                if e > ce {
                    cur = Some((cs, e));
                }
            }
            _ => {
                if let Some((cs, ce)) = cur {
                    total += ce - cs + 1;
                }
                cur = Some((s, e));
            }
        }
    }
    if let Some((cs, ce)) = cur {
        total += ce - cs + 1;
    }
    total
}

/// - 正确: 组装本链 == 某参考链（chrom+链完全一致; 链方向不参与, 见模块注释）
/// - 基因找回: 任一多外显子参考链被复现, 或任一单外显子参考转录本被
///   同 chrom、ident≥0.95 的组装本 blocks 并集覆盖 ≥70% 基因组跨度
/// - Jtx 分子: 精确复现的多外显子参考（单外显子回退只计入 Jgene; 手算用例锁定）
pub fn eval_junction(asm: &[AsmTx], model: &RefModel) -> JunctionMetrics {
    let introns = model.intron_set();
    let chains = model.chain_set();

    // 组装本侧: 精确链 / 新奇内含子 / 单外显子
    let (mut n_multi, mut n_correct, mut n_novel, mut n_single) = (0usize, 0usize, 0usize, 0usize);
    for a in asm {
        if a.chain.is_empty() {
            n_single += 1;
            continue;
        }
        n_multi += 1;
        if chains.contains(&(a.chrom.clone(), a.chain.clone())) {
            n_correct += 1;
        }
        if a.chain
            .iter()
            .any(|iv| !introns.contains(&(a.chrom.clone(), *iv)))
        {
            n_novel += 1;
        }
    }

    // 参考侧: 转录本精确复现 + 基因找回（链复现或单外显子回退）, 均不比对链方向
    let n_gene = model.genes.len();
    let (mut n_tx, mut tx_exact, mut genes_hit) = (0usize, 0usize, 0usize);
    for txs in model.genes.values() {
        n_tx += txs.len();
        let mut hit = false;
        for r in txs {
            if r.chain.is_empty() {
                continue;
            }
            if asm.iter().any(|a| a.chrom == r.chrom && a.chain == r.chain) {
                hit = true;
                tx_exact += 1;
            }
        }
        if !hit {
            for r in txs.iter().filter(|r| r.chain.is_empty()) {
                let span = r.span();
                let span_len = span.1.saturating_sub(span.0) + 1;
                if span_len == 0 {
                    continue;
                }
                let blocks: Vec<(u64, u64)> = asm
                    .iter()
                    .filter(|a| a.chrom == r.chrom && a.ident >= 0.95)
                    .flat_map(|a| a.blocks.iter().copied())
                    .collect();
                if union_covered(&blocks, span) as f64 >= 0.7 * span_len as f64 {
                    hit = true;
                    break;
                }
            }
        }
        if hit {
            genes_hit += 1;
        }
    }

    let rate = |num: usize, den: usize| {
        if den == 0 {
            0.0
        } else {
            num as f64 / den as f64
        }
    };
    JunctionMetrics {
        j_gene: rate(genes_hit, n_gene),
        j_tx: rate(tx_exact, n_tx),
        j_prec: rate(n_correct, n_multi),
        novel_j: rate(n_novel, n_multi),
        single_exon: rate(n_single, asm.len()),
        n_gene,
        n_tx,
        n_asm: asm.len(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        assert_eq!(g1.len(), 2);
        let tx11 = g1.iter().find(|t| t.tx == "STRG.1.1").unwrap();
        assert_eq!(tx11.chain, vec![(201u64, 300u64), (401, 500)]);
        assert_eq!(tx11.chrom, "chr1");
        assert_eq!(tx11.strand, '+');
        let tx12 = g1.iter().find(|t| t.tx == "STRG.1.2").unwrap();
        assert_eq!(tx12.chain, vec![(201, 500)]);
        let g2 = &model.genes["STRG.2"];
        assert!(g2[0].chain.is_empty());
        assert_eq!(g2[0].span(), (1000, 1500));
    }

    #[test]
    fn intron_set_dedupes_across_transcripts() {
        let model = parse_reference_gtf(GTF);
        let set = model.intron_set();
        assert!(set.contains(&(String::from("chr1"), (201, 300))));
        assert!(set.contains(&(String::from("chr1"), (401, 500))));
        assert!(set.contains(&(String::from("chr1"), (201, 500))));
        assert_eq!(set.len(), 3);
    }

    #[test]
    fn introns_from_cigar() {
        let chain = super::introns_from_cigar("100M200N50M30I20M", 1000);
        assert_eq!(chain, vec![(1101, 1300)]);
        assert!(super::introns_from_cigar("300M", 1000).is_empty());
        assert!(super::introns_from_cigar("100M10D100M", 0).is_empty());
    }

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
        assert!((asm[0].ident - 540.0 / 560.0).abs() < 1e-9);
        assert!(asm[1].chain.is_empty());
    }

    /// 链方向无关匹配: 参考链 chr1 '+' 与组装本 chr1 '-' 同链坐标即算正确;
    /// 单外显子参考回退同样不要求同链; intron_set 不再区分链方向。
    #[test]
    fn strand_agnostic_matching() {
        // 参考: S1.1 = chr1 '+' 双外显子(内含子 201-300); S2.1 = chr4 '+' 单外显子
        let gtf = "chr1\tStringTie\texon\t100\t200\t1000\t+\t.\tgene_id \"S1\"; transcript_id \"S1.1\"; exon_number \"1\";\n\
                   chr1\tStringTie\texon\t301\t400\t1000\t+\t.\tgene_id \"S1\"; transcript_id \"S1.1\"; exon_number \"2\";\n\
                   chr4\tStringTie\texon\t1000\t1500\t1000\t+\t.\tgene_id \"S2\"; transcript_id \"S2.1\"; exon_number \"1\";\n";
        let model = parse_reference_gtf(gtf);
        // 组装本方向与参考相反
        let asm = vec![
            AsmTx {
                qid: "B1".into(),
                chrom: "chr1".into(),
                strand: '-',
                chain: vec![(201, 300)],
                ident: 0.99,
                blocks: vec![(100, 400)],
                qlen_hit: 300,
            },
            AsmTx {
                qid: "B2".into(),
                chrom: "chr4".into(),
                strand: '-',
                chain: vec![],
                ident: 0.99,
                blocks: vec![(1000, 1500)],
                qlen_hit: 501,
            },
        ];
        let m = eval_junction(&asm, &model);
        assert!((m.j_prec - 1.0).abs() < 1e-9, "反向链多外显子链应算正确");
        assert!((m.j_gene - 1.0).abs() < 1e-9, "反向链单外显子参考应被找回");
        assert!((m.j_tx - 0.5).abs() < 1e-9);
        assert!((m.novel_j - 0.0).abs() < 1e-9, "反向链已知内含子不应算新奇");
        assert!((m.single_exon - 0.5).abs() < 1e-9);
        assert_eq!((m.n_gene, m.n_tx, m.n_asm), (2, 2, 2));

        // intron_set: 同 chrom 同坐标、方向相反的参考内含子合并为 1 条
        let model2 = RefModel {
            genes: HashMap::from([(
                String::from("G"),
                vec![
                    RefTx {
                        tx: String::from("G.1"),
                        chrom: String::from("chr1"),
                        strand: '+',
                        exons: vec![(100, 200), (301, 400)],
                        chain: vec![(201, 300)],
                    },
                    RefTx {
                        tx: String::from("G.2"),
                        chrom: String::from("chr1"),
                        strand: '-',
                        exons: vec![(100, 200), (301, 400)],
                        chain: vec![(201, 300)],
                    },
                ],
            )]),
        };
        let set = model2.intron_set();
        assert_eq!(set.len(), 1, "反向链等坐标内含子应去重为 1");
        assert!(set.contains(&(String::from("chr1"), (201, 300))));
    }

    #[test]
    fn eval_junction_hand_case() {
        let model = parse_reference_gtf(GTF);
        let asm = vec![
            AsmTx {
                qid: "A1".into(),
                chrom: "chr1".into(),
                strand: '+',
                chain: vec![(201, 300), (401, 500)],
                ident: 0.99,
                blocks: vec![(100, 600)],
                qlen_hit: 500,
            },
            AsmTx {
                qid: "A2".into(),
                chrom: "chr1".into(),
                strand: '+',
                chain: vec![(201, 999)],
                ident: 0.99,
                blocks: vec![(100, 1000)],
                qlen_hit: 900,
            },
            AsmTx {
                qid: "A3".into(),
                chrom: "chr4".into(),
                strand: '-',
                chain: vec![],
                ident: 0.99,
                blocks: vec![(1000, 1500)],
                qlen_hit: 501,
            },
            AsmTx {
                qid: "A4".into(),
                chrom: "chr1".into(),
                strand: '+',
                chain: vec![],
                ident: 0.99,
                blocks: vec![(5000, 5200)],
                qlen_hit: 201,
            },
        ];
        let m = eval_junction(&asm, &model);
        assert!((m.j_prec - 0.5).abs() < 1e-9);
        assert!((m.j_gene - 1.0).abs() < 1e-9);
        assert!((m.j_tx - 1.0 / 3.0).abs() < 1e-9);
        assert!((m.novel_j - 0.5).abs() < 1e-9);
        assert!((m.single_exon - 0.5).abs() < 1e-9);
        assert_eq!((m.n_gene, m.n_tx, m.n_asm), (2, 3, 4));
    }
}
