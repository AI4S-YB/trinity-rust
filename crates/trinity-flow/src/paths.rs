//! 路径 → 序列转换 + 相似去冗余 + FASTA 输出。

pub struct PathOutput {
    pub seq: Vec<u8>,
    pub support: u64,
    pub path: Vec<i64>,
}

/// 首节点全 kmer + 后续节点按最长后缀-前缀重叠拼接。
/// 等长 k-mer 图（每节点 24bp 相邻差 1）上退化为「后续各取末 1 碱基」。
pub fn path_to_seq(path: &[i64], kmers: &[(i64, Vec<u8>)], _k: usize) -> Vec<u8> {
    let map: std::collections::HashMap<i64, &[u8]> =
        kmers.iter().map(|(a, b)| (*a, b.as_slice())).collect();
    let mut seq = Vec::new();
    for &id in path {
        let Some(kmer) = map.get(&id).copied() else {
            continue;
        };
        if seq.is_empty() {
            seq.extend_from_slice(kmer);
            continue;
        }
        let max_ov = seq.len().min(kmer.len());
        let mut ov = max_ov;
        while ov > 0 && seq[seq.len() - ov..] != kmer[..ov] {
            ov -= 1;
        }
        seq.extend_from_slice(&kmer[ov..]);
    }
    seq
}

/// 简化 cd-hit 式去冗余, 两级:
/// 1. 包含淘汰: 序列是某条**更长**路径的子串 → 剥离残余（源/汇中途起止的
///    残余流）, 丢弃（长度降序处理, 保长不保 support）;
/// 2. 逐位同一率 ≥ threshold → 保留高 support。
pub fn dedup_similar(paths: Vec<PathOutput>, threshold: f64) -> Vec<PathOutput> {
    // --- 1. 包含淘汰（长度降序）---
    let mut by_len = paths;
    by_len.sort_by_key(|p| std::cmp::Reverse(p.seq.len()));
    let mut survivors: Vec<PathOutput> = Vec::new();
    for p in by_len {
        let contained = survivors.iter().any(|q| {
            p.seq.len() <= q.seq.len()
                && !p.seq.is_empty()
                && q.seq.windows(p.seq.len()).any(|w| w == p.seq.as_slice())
        });
        if !contained {
            survivors.push(p);
        }
    }
    // --- 2. 同一率去冗余（support 降序）---
    let mut sorted = survivors;
    sorted.sort_by_key(|p| std::cmp::Reverse(p.support));
    let mut out: Vec<PathOutput> = Vec::new();
    for p in sorted {
        let dup = out.iter().any(|q| {
            let min_len = p.seq.len().min(q.seq.len());
            if min_len == 0 {
                return false;
            }
            let same = p
                .seq
                .iter()
                .zip(q.seq.iter())
                .filter(|(a, b)| a == b)
                .count();
            (same as f64 / min_len as f64) >= threshold
        });
        if !dup {
            out.push(p);
        }
    }
    out
}

pub fn to_fasta(paths: &[PathOutput], comp: &str, gene: usize, wrap: usize) -> String {
    let mut s = String::new();
    for (i, p) in paths.iter().enumerate() {
        s.push_str(&format!(
            ">{}_g{}_i{} len={} support={} path=[{}]\n",
            comp,
            gene + 1,
            i + 1,
            p.seq.len(),
            p.support,
            p.path
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(",")
        ));
        for chunk in p.seq.chunks(wrap.max(1)) {
            s.push_str(&String::from_utf8_lossy(chunk));
            s.push('\n');
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn path_to_sequence_chain() {
        // K=4: "AAAA" → "AAAC" → "AACC" → 首全串 + 后续各末 1 碱基 = "AAAACC"
        let kmers = vec![
            (1i64, b"AAAA".to_vec()),
            (2, b"AAAC".to_vec()),
            (3, b"AACC".to_vec()),
        ];
        assert_eq!(path_to_seq(&[1, 2, 3], &kmers, 4), b"AAAACC".to_vec());
    }

    #[test]
    fn path_to_sequence_branch() {
        // 菱形: 1→2→4 与 1→3→4（2/3 同为中间节点, 各去前 K-1）
        let kmers = vec![(1i64, b"AAAACCCCGGGG".to_vec()), (2, b"CCCGGGG".to_vec())];
        assert_eq!(path_to_seq(&[1, 2], &kmers, 7), b"AAAACCCCGGGG".to_vec());
        // 注意: 中间节点 kmer 从 graph 解析来时是全 24-mer（不是增量）——
        // path_to_seq 需要用 K-1 前缀跳过。这里 K=7 时节点 1 的 kmer 长度是 12（已被 min100 压缩图的特征）,
        // 实际上 flow 引擎直接消费原始 24-mer 图, 每节点都是 24bp。首节点全长 + 后续各去前 K-1。
    }

    #[test]
    fn dedup_similar() {
        // 两条路径, 序列 95% 相同且支持度不同 → 保留高支持度
        let seqs = vec![
            PathOutput {
                seq: b"AAAACCCCGGGGTTTTAAAACCCC".to_vec(),
                support: 10,
                path: vec![1, 2],
            },
            PathOutput {
                seq: b"AAAACCCCGGGGTTTTAAAACCCG".to_vec(),
                support: 5,
                path: vec![1, 3],
            },
        ];
        let out = super::dedup_similar(seqs, 0.95);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].support, 10);
    }

    #[test]
    fn dedup_drops_contained_peel_residue() {
        // 剥离残余: 短路径是长路径的子串（support 更高也仍丢弃——保长）
        let mk = |seq: &[u8], support: u64| PathOutput {
            seq: seq.to_vec(),
            support,
            path: vec![1, 2],
        };
        let long = b"AAAACCCCGGGGTTTTAAAACCCC";
        let out = super::dedup_similar(
            vec![mk(&long[..10], 5), mk(long, 1)],
            0.95,
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].seq, long);
        assert_eq!(out[0].support, 1);
    }

    #[test]
    fn fasta_format() {
        let po = vec![PathOutput {
            seq: vec![b'A'; 130],
            support: 5,
            path: vec![1, 2, 3],
        }];
        let fa = to_fasta(&po, "c0", 0, 60);
        assert!(fa.starts_with(">c0_g1_i1 len=130"));
        assert!(fa.contains("\n")); // 60 列折行
        let lines: Vec<&str> = fa.lines().collect();
        assert_eq!(lines[1].len(), 60);
        assert_eq!(lines[2].len(), 60);
        assert_eq!(lines[3].len(), 10);
    }
}
