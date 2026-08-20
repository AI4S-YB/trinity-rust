//! xtask bench —— F1 基准: 组装(缓存) → minimap2 → R/P/F1。
//! 口径锁定见 M0 计划; 纯函数 compute_metrics 便于单测。

use std::collections::HashMap;

pub const MINIMAP2: &str = "/public/home/senior007/miniconda3/envs/cactus/bin/minimap2";

#[derive(Debug, Default, PartialEq)]
pub struct Metrics {
    pub refs_total: usize,
    pub r_full: f64,     // ≥95%cov & ≥95%id
    pub r70: f64,        // ≥70%cov & ≥90%id
    pub precision: f64,  // 碱基级
    pub f1: f64,
    pub n50: usize,
    pub denovo_total: usize,
}

/// PAF 行: q qlen qs qe strand t tlen ts te mch aln ...
/// ref_lens: 参照名→长度; denovo_lens: de novo 名→长度
pub fn compute_metrics(paf: &str, ref_lens: &HashMap<String, usize>, denovo_lens: &HashMap<String, usize>) -> Metrics {
    let mut best: HashMap<&str, (f64, f64)> = HashMap::new(); // ref -> (best cov, best ident)（按 cov*ident 取）
    let mut qhit: HashMap<&str, usize> = HashMap::new();
    for line in paf.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 11 { continue; }
        let (q, t) = (f[0], f[5]);
        let qs: usize = f[2].parse().unwrap_or(0);
        let qe: usize = f[3].parse().unwrap_or(0);
        let ts: usize = f[7].parse().unwrap_or(0);
        let te: usize = f[8].parse().unwrap_or(0);
        let mch: usize = f[9].parse().unwrap_or(0);
        let aln: usize = f[10].parse().unwrap_or(0);
        let ident = if aln > 0 { mch as f64 / aln as f64 } else { 0.0 };
        if let Some(&tl) = ref_lens.get(t) {
            let cov = te.saturating_sub(ts) as f64 / tl as f64;
            let score = cov * ident;
            let e = best.entry(t).or_insert((-1.0, 0.0));
            if score > e.0 * e.1 { *e = (cov, ident); }
        }
        if ident >= 0.95 {
            let span = qe.saturating_sub(qs);
            let e = qhit.entry(q).or_insert(0);
            if span > *e { *e = span; }
        }
    }
    let refs_total = ref_lens.len();
    let mut r_full = 0usize;
    let mut r70 = 0usize;
    for (_, &(cov, ident)) in best.iter() {
        if cov >= 0.95 && ident >= 0.95 { r_full += 1; }
        if cov >= 0.70 && ident >= 0.90 { r70 += 1; }
    }
    let total_bases: usize = denovo_lens.values().sum();
    let hit_bases: usize = qhit.values().sum();
    let precision = if total_bases > 0 { hit_bases as f64 / total_bases as f64 } else { 0.0 };
    let r = if refs_total > 0 { r_full as f64 / refs_total as f64 } else { 0.0 };
    let f1 = if r + precision > 0.0 { 2.0 * r * precision / (r + precision) } else { 0.0 };
    // N50
    let mut lens: Vec<usize> = denovo_lens.values().copied().collect();
    lens.sort_unstable_by(|a, b| b.cmp(a));
    let half = total_bases as f64 / 2.0;
    let (mut acc, mut n50) = (0usize, 0usize);
    for &l in &lens { acc += l; if acc as f64 >= half { n50 = l; break; } }
    Metrics { refs_total, r_full: r, r70: if refs_total>0 { r70 as f64/refs_total as f64 } else {0.0},
              precision, f1, n50, denovo_total: denovo_lens.len() }
}

fn read_fasta_lens(path: &std::path::Path) -> std::io::Result<HashMap<String, usize>> {
    let mut m = HashMap::new();
    let mut name: Option<String> = None; let mut len = 0usize;
    for line in std::io::BufRead::lines(&mut std::io::BufReader::new(std::fs::File::open(path)?)) {
        let line = line?;
        if line.starts_with('>') {
            if let Some(n) = name { m.insert(n, len); }
            name = Some(line[1..].split_whitespace().next().unwrap_or("").to_string());
            len = 0;
        } else { len += line.trim().len(); }
    }
    if let Some(n) = name { m.insert(n, len); }
    Ok(m)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hm(pairs: &[(&str, usize)]) -> HashMap<String, usize> {
        pairs.iter().map(|(a,b)| (a.to_string(), *b)).collect()
    }

    #[test]
    fn metrics_hand_case() {
        // 参照 r1=1000bp, r2=500bp; de novo d1=900bp, d2=400bp, d3=100bp(垃圾)
        let ref_l = hm(&[("r1",1000),("r2",500)]);
        let dn_l  = hm(&[("d1",900),("d2",400),("d3",100)]);
        // d1 覆盖 r1 960/1000=0.96, id=0.97 → r1 全长
        // d2 覆盖 r2 400/500=0.80, id=0.95 → R70 计入
        let paf = "d1\t900\t0\t960\t+\tr1\t1000\t40\t1000\t941\t970\n\
                   d2\t400\t0\t400\t+\tr2\t500\t100\t500\t390\t400\n\
                   d3\t100\t0\t50\t+\tr1\t1000\t0\t50\t30\t50\n";
        // d3 id=30/50=0.6 <0.95 不计精确度; d1/d2 全量计入
        let m = compute_metrics(paf, &ref_l, &dn_l);
        assert_eq!(m.refs_total, 2);
        assert!((m.r_full - 0.5).abs() < 1e-9);      // r1 全长
        assert!((m.r70 - 1.0).abs() < 1e-9);         // r1+r2
        assert!((m.precision - (960.0+400.0)/1400.0).abs() < 1e-9);
        assert!((m.f1 - 2.0*0.5*(1360.0/1400.0)/(0.5+1360.0/1400.0)).abs() < 1e-9);
    }

    #[test]
    fn n50_case() {
        let dn_l = hm(&[("a",100),("b",900),("c",400)]);
        let m = compute_metrics("", &HashMap::new(), &dn_l);
        assert_eq!(m.n50, 900); // total 1400, half 700, 900 先跨过
        assert!((m.precision - 0.0).abs() < 1e-9);
        assert_eq!(m.denovo_total, 3);
    }
}
