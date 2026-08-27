//! 候选重排数据集构建器（M9, 全长级判别模型）。
//!
//! 输入: bench 缓存（cache 根 / out 目录, 内含 chrysalis/Component_bins 的
//! c*.graph.out + c*.graph.reads）与参照转录本 FASTA。对每个组件以**库方式**
//! 离线重跑 trinity-butterfly 的 classic 引擎（read 支持过滤放宽到
//! min_read_support=1, 管线零改动、不重跑组装）, 收集候选转录本,
//! （特性面精简后 PogFlow 已移除; classic 头无 support 字段 → 该列填 -1）,
//! 提取全长级特征（序列浅特征 + 最长 ORF + 25-mer 覆盖度形状 + read 支持）,
//! 再用 minimap2 锚定参照打 oracle 标签, 每候选一行写 TSV。
//!
//! 标签语义（镜像 bench.rs compute_metrics 同口径）: qcov = query 侧覆盖
//! （ident ≥ 0.95 命中的最长 query 比对段 / qlen, eval 精度方向）;
//! 正例 = qcov ≥ 0.95。ref_cov/ref_ident 为召回方向参照侧 best hit, 仅切片用。
//! （首版曾误用 ref_cov 作标签——把真实片段全判负, 9.5% 正例 vs eval P=84%,
//! 口径错配教训入档。）
//! ref_cov/ref_ident/qcov/label 为 oracle 列——仅训练/切片用, 推理禁用。

use std::collections::HashMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use rayon::prelude::*;
use trinity_butterfly::{run_component, ComponentParams};
use trinity_common::kmer::{get_ds_kmer_val, kmer_to_intval};
use trinity_kmer::counter::CountMap;

const KMER: usize = 25; // 计数表键长（graph KMER_SIZE; graph.out 顶点名长 24）

pub const TSV_HEADER: &str = "comp\tcand\tlen\tsupport\tgc\tent\thom\torf_aa\torf_frac\torf_m\tcov_mean\tcov_median\tcov_min\tcov_ge2_frac\tcov_cv\tref_cov\tref_ident\tqcov\tlabel";

// ---------------------------------------------------------------------------
// 密码子翻译与数据工具（特性面精简后自 protein_extend/export_edges 收编内联）
// ---------------------------------------------------------------------------

pub fn translate_codon(codon: &[u8]) -> u8 {
    let n = |b: u8| match b {
        b'T' => 0,
        b'C' => 1,
        b'A' => 2,
        b'G' => 3,
        _ => 4,
    };
    match codon.len() {
        3 => match (n(codon[0]), n(codon[1]), n(codon[2])) {
            (0, 0, 0) => b'F',
            (0, 0, 1) => b'F',
            (0, 0, 2) => b'L',
            (0, 0, 3) => b'L',
            (0, 1, 0) => b'S',
            (0, 1, 1) => b'S',
            (0, 1, 2) => b'S',
            (0, 1, 3) => b'S',
            (0, 2, 0) => b'Y',
            (0, 2, 1) => b'Y',
            (0, 2, 2) => b'*',
            (0, 2, 3) => b'*',
            (0, 3, 0) => b'C',
            (0, 3, 1) => b'C',
            (0, 3, 2) => b'*',
            (0, 3, 3) => b'W',
            (1, 0, 0) => b'L',
            (1, 0, 1) => b'L',
            (1, 0, 2) => b'L',
            (1, 0, 3) => b'L',
            (1, 1, 0) => b'P',
            (1, 1, 1) => b'P',
            (1, 1, 2) => b'P',
            (1, 1, 3) => b'P',
            (1, 2, 0) => b'H',
            (1, 2, 1) => b'H',
            (1, 2, 2) => b'Q',
            (1, 2, 3) => b'Q',
            (1, 3, 0) => b'R',
            (1, 3, 1) => b'R',
            (1, 3, 2) => b'R',
            (1, 3, 3) => b'R',
            (2, 0, 0) => b'I',
            (2, 0, 1) => b'I',
            (2, 0, 2) => b'I',
            (2, 0, 3) => b'M',
            (2, 1, 0) => b'T',
            (2, 1, 1) => b'T',
            (2, 1, 2) => b'T',
            (2, 1, 3) => b'T',
            (2, 2, 0) => b'N',
            (2, 2, 1) => b'N',
            (2, 2, 2) => b'K',
            (2, 2, 3) => b'K',
            (2, 3, 0) => b'S',
            (2, 3, 1) => b'S',
            (2, 3, 2) => b'R',
            (2, 3, 3) => b'R',
            (3, 0, 0) => b'V',
            (3, 0, 1) => b'V',
            (3, 0, 2) => b'V',
            (3, 0, 3) => b'V',
            (3, 1, 0) => b'A',
            (3, 1, 1) => b'A',
            (3, 1, 2) => b'A',
            (3, 1, 3) => b'A',
            (3, 2, 0) => b'D',
            (3, 2, 1) => b'D',
            (3, 2, 2) => b'E',
            (3, 2, 3) => b'E',
            (3, 3, 0) => b'G',
            (3, 3, 1) => b'G',
            (3, 3, 2) => b'G',
            (3, 3, 3) => b'G',
            _ => b'X',
        },
        _ => b'X',
    }
}

pub fn gc_frac(s: &str) -> f64 {
    let mut gc = 0usize;
    let mut n = 0usize;
    for b in s.bytes() {
        match b {
            b'G' | b'C' => {
                gc += 1;
                n += 1;
            }
            b'A' | b'T' => n += 1,
            _ => {}
        }
    }
    if n == 0 {
        0.0
    } else {
        gc as f64 / n as f64
    }
}

pub fn entropy(s: &str) -> f64 {
    let mut cnt = [0usize; 4];
    let mut n = 0usize;
    for b in s.bytes() {
        let i = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => continue,
        };
        cnt[i] += 1;
        n += 1;
    }
    if n == 0 {
        return 0.0;
    }
    let mut h = 0.0f64;
    for c in cnt {
        if c > 0 {
            let p = c as f64 / n as f64;
            h -= p * p.log2();
        }
    }
    h
}

pub fn max_homopolymer(s: &str) -> usize {
    let mut best = 0;
    let mut run = 0;
    let mut prev = 0u8;
    for b in s.bytes() {
        if b == prev && matches!(b, b'A' | b'C' | b'G' | b'T') {
            run += 1;
        } else {
            run = if matches!(b, b'A' | b'C' | b'G' | b'T') {
                1
            } else {
                0
            };
        }
        best = best.max(run);
        prev = b;
    }
    best
}

/// 定位管线 asm out 目录（含 chrysalis/ 的那层）: 接受 cache 根或 out 本身。
pub fn locate_asm_out(p: &Path) -> Option<PathBuf> {
    [p.to_path_buf(), p.join("out")]
        .into_iter()
        .find(|cand| cand.join("chrysalis").is_dir())
}

/// 取 25-mer DS 计数表: override 路径 > mer_counts bin > norm.fq 重数。
/// 返回 (表, 来源说明)。v1 只支持 K=25 DS（dev/real 全部数据集的形态）。
pub fn counts_for_out_dir(asm_out: &Path, override_path: &str) -> (CountMap, String) {
    let bin = if !override_path.is_empty() {
        PathBuf::from(override_path)
    } else {
        asm_out.join("mer_counts.25.asm.rs.bin")
    };
    if bin.is_file() {
        let m = trinity_kmer::counter::load_counts_bin(&bin, KMER, true)
            .unwrap_or_else(|e| panic!("load {}: {e}", bin.display()));
        return (m, format!("bin {}", bin.display()));
    }
    let norm = asm_out.join("insilico_read_normalization");
    let mut merged = CountMap::default();
    let mut n_files = 0;
    for name in ["left.norm.fq", "right.norm.fq"] {
        let p = norm.join(name);
        if !p.is_file() {
            continue;
        }
        let data = std::fs::read(&p).unwrap_or_else(|e| panic!("read {}: {e}", p.display()));
        // fq → fa（4 行记录取序列行; both.fa 经 seqtk -A 转换, 序列多重集相同）
        let mut fa = Vec::with_capacity(data.len() / 2 + 64);
        for (i, line) in data.split(|&b| b == b'\n').enumerate() {
            if i % 4 == 1 && !line.is_empty() {
                fa.extend_from_slice(b">r\n");
                fa.extend_from_slice(line);
                fa.push(b'\n');
            }
        }
        let m = trinity_kmer::counter::KmerCountTable::count_fasta_data_with_capacity(
            &fa,
            KMER,
            true,
            1 << 20,
        );
        for (k, v) in m {
            *merged.entry(k).or_insert(0) += v;
        }
        n_files += 1;
    }
    assert!(
        n_files > 0,
        "{} 下既无 mer_counts.25.asm.rs.bin 也无 left/right.norm.fq",
        asm_out.display()
    );
    (merged, format!("norm.fq 重数 ×{n_files}"))
}

/// 接受 cache 根（内含 out/）、out 目录、或 chrysalis 目录, 返回 Component_bins。
pub fn locate_component_bins(p: &Path) -> Option<PathBuf> {
    [
        p.join("chrysalis/Component_bins"),
        p.join("out/chrysalis/Component_bins"),
        p.join("Component_bins"),
    ]
    .into_iter()
    .find(|cand| cand.is_dir())
}

pub fn read_fasta(path: &Path) -> Vec<(String, String)> {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let text: String = if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut s = String::new();
        use std::io::Read;
        flate2::read::GzDecoder::new(&bytes[..])
            .read_to_string(&mut s)
            .unwrap_or_else(|e| panic!("gunzip {}: {e}", path.display()));
        s
    } else {
        String::from_utf8(bytes).unwrap_or_else(|e| panic!("utf8 {}: {e}", path.display()))
    };
    let mut out = Vec::new();
    let mut name = String::new();
    let mut seq = String::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if !name.is_empty() {
                out.push((std::mem::take(&mut name), std::mem::take(&mut seq)));
            }
            name = h.split_whitespace().next().unwrap_or("").to_string();
        } else {
            seq.push_str(line.trim());
        }
    }
    if !name.is_empty() {
        out.push((name, seq));
    }
    out
}

// ---------------------------------------------------------------------------
// 候选 FASTA 解析
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
pub struct Candidate {
    pub comp: i32,
    pub name: String,
    pub len: usize,
    pub support: i64,
    pub seq: String,
}

/// 解析 run_component 的 allProbPaths FASTA（classic print_final_paths 产物:
/// header `>c6341_g1_i2 len=1730 path=[...] {...}`, 序列 60 列折行;
/// classic 无 support= 字段 → support 记 -1（未知, 与 0 区分））。
/// header 缺 len= 时回退用拼接后的序列长。
pub fn parse_candidates(comp: i32, fasta: &str) -> Vec<Candidate> {
    let mut out: Vec<Candidate> = Vec::new();
    let mut name = String::new();
    let mut len: Option<usize> = None;
    let mut support = -1i64;
    let mut seq = String::new();
    for line in fasta.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if !name.is_empty() {
                out.push(Candidate {
                    comp,
                    name: std::mem::take(&mut name),
                    len: len.take().unwrap_or(seq.len()),
                    support,
                    seq: std::mem::take(&mut seq),
                });
            }
            let mut toks = h.split_whitespace();
            name = toks.next().unwrap_or("").to_string();
            len = None;
            support = -1;
            for t in toks {
                if let Some(v) = t.strip_prefix("len=") {
                    len = v.parse().ok();
                } else if let Some(v) = t.strip_prefix("support=") {
                    support = v.parse().unwrap_or(-1);
                }
            }
        } else {
            seq.push_str(line.trim());
        }
    }
    if !name.is_empty() {
        out.push(Candidate {
            comp,
            name,
            len: len.unwrap_or(seq.len()),
            support,
            seq,
        });
    }
    out
}

// ---------------------------------------------------------------------------
// 最长 ORF 特征
// ---------------------------------------------------------------------------

/// 最长 ORF 特征（语义镜像 bench.rs translate_longest_orf: 序列先大写化,
/// 3 个正向框架按 '*' 切段, 段内有 M 则从首个 M 起取, 全体取最长段;
/// 等长平局保留先到者, 与 bench.rs 的 `>` 比较一致）。
/// 返回 (orf_aa, orf_frac, orf_m): 最长 ORF aa 数 / orf_aa*3/len / 是否起始于 M(0/1)。
pub fn orf_features(seq: &str) -> (usize, f64, u8) {
    let su = seq.to_ascii_uppercase();
    let b = su.as_bytes();
    let mut best = 0usize;
    let mut best_m = 0u8;
    for frame in 0..3 {
        let mut prot: Vec<u8> = Vec::new();
        let mut i = frame;
        while i + 3 <= b.len() {
            prot.push(translate_codon(&b[i..i + 3]));
            i += 3;
        }
        for seg in prot.split(|&c| c == b'*') {
            let (cand, m) = match seg.iter().position(|&c| c == b'M') {
                Some(p) => (&seg[p..], 1u8),
                None => (seg, 0u8),
            };
            if cand.len() > best {
                best = cand.len();
                best_m = m;
            }
        }
    }
    let frac = if b.is_empty() {
        0.0
    } else {
        best as f64 * 3.0 / b.len() as f64
    };
    (best, frac, best_m)
}

// ---------------------------------------------------------------------------
// 覆盖度形状
// ---------------------------------------------------------------------------

/// 覆盖度形状统计（沿序列的 25-mer 计数分布）。
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct CovStats {
    pub mean: f64,
    pub median: f64,
    pub min: f64,
    /// 计数 ≥2 窗口占比
    pub ge2_frac: f64,
    /// 变异系数（总体标准差/均值; 均值 0 → 0）
    pub cv: f64,
}

/// 沿候选序列滑 25-mer 窗查计数表（含非 ACGT 字符的窗跳过, 表缺 = 0）。
/// 序列短于 25 或无有效窗 → 全 0。中位数: 偶数窗取中间两值均值。
pub fn cov_features(counts: &CountMap, seq: &str) -> CovStats {
    let b = seq.as_bytes();
    if b.len() < KMER {
        return CovStats::default();
    }
    let mut vals: Vec<f64> = Vec::with_capacity(b.len() - KMER + 1);
    for w in b.windows(KMER) {
        if !w.iter().all(|&c| matches!(c, b'A' | b'C' | b'G' | b'T')) {
            continue;
        }
        let Ok(val) = kmer_to_intval(w) else {
            continue;
        };
        let key = get_ds_kmer_val(val, KMER);
        vals.push(counts.get(&key).map(|&c| c as f64).unwrap_or(0.0));
    }
    if vals.is_empty() {
        return CovStats::default();
    }
    let n = vals.len() as f64;
    let mean = vals.iter().sum::<f64>() / n;
    let mut sorted = vals.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = sorted.len();
    let median = if m % 2 == 1 {
        sorted[m / 2]
    } else {
        (sorted[m / 2 - 1] + sorted[m / 2]) / 2.0
    };
    let ge2_frac = vals.iter().filter(|&&v| v >= 2.0).count() as f64 / n;
    let cv = if mean > 0.0 {
        let var = vals.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n;
        var.sqrt() / mean
    } else {
        0.0
    };
    CovStats {
        mean,
        median,
        min: sorted[0],
        ge2_frac,
        cv,
    }
}

// ---------------------------------------------------------------------------
// minimap2 PAF best-hit 与标签
// ---------------------------------------------------------------------------

/// minimap2 PAF → 每个 query（候选）的 (ref_cov, ref_ident, qcov)。
/// ref_cov/ref_ident: 参照侧 best hit（语义镜像 bench.rs compute_metrics:
/// 只计 ref_lens 已知的目标; cov = 参照侧覆盖 (te-ts)/tlen; ident = mch/aln;
/// best 按 cov*ident 取）。
/// qcov: query 侧覆盖——ident ≥ 0.95 的命中里最长 query 比对段 / qlen
/// （镜像 compute_metrics 的 qhit 逻辑; 与 eval 精度（碱基级）同方向,
/// 转录本粒度）。**标签的正确方向**: "这条候选是不是真实序列"= qcov,
/// 不是 ref_cov（ref_cov 是召回方向, 会把真实片段误判负——首版踩坑,
/// 745/7869 = 9.5% 正例 vs eval P=84% 即此口径错配）。
pub fn best_hits(paf: &str, ref_lens: &HashMap<String, usize>) -> HashMap<String, (f64, f64, f64)> {
    let mut best: HashMap<String, (f64, f64)> = HashMap::new();
    let mut qspan: HashMap<String, (usize, usize)> = HashMap::new(); // q -> (qlen, max span)
    for line in paf.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 11 {
            continue;
        }
        let qlen: usize = f[1].parse().unwrap_or(0);
        let qs: usize = f[2].parse().unwrap_or(0);
        let qe: usize = f[3].parse().unwrap_or(0);
        let mch: usize = f[9].parse().unwrap_or(0);
        let aln: usize = f[10].parse().unwrap_or(0);
        let ident = if aln > 0 {
            mch as f64 / aln as f64
        } else {
            0.0
        };
        if ident >= 0.95 {
            let span = qe.saturating_sub(qs);
            let e = qspan.entry(f[0].to_string()).or_insert((qlen, 0));
            if span > e.1 {
                *e = (qlen, span);
            }
        }
        let Some(&tl) = ref_lens.get(f[5]) else {
            continue;
        };
        let ts: usize = f[7].parse().unwrap_or(0);
        let te: usize = f[8].parse().unwrap_or(0);
        let cov = te.saturating_sub(ts) as f64 / tl as f64;
        let score = cov * ident;
        match best.get(f[0]) {
            Some(&(bc, bi)) if bc * bi >= score => {}
            _ => {
                best.insert(f[0].to_string(), (cov, ident));
            }
        }
    }
    let mut out: HashMap<String, (f64, f64, f64)> = HashMap::new();
    for (q, (qlen, span)) in &qspan {
        let qcov = if *qlen > 0 {
            *span as f64 / *qlen as f64
        } else {
            0.0
        };
        let (rc, ri) = best.get(q).copied().unwrap_or((0.0, 0.0));
        out.insert(q.clone(), (rc, ri, qcov));
    }
    // 只有 ref 侧命中（ident < 0.95）的候选也进表（qcov=0）
    for (q, &(rc, ri)) in &best {
        out.entry(q.clone()).or_insert((rc, ri, 0.0));
    }
    out
}

/// 正例阈值: qcov ≥ 0.95（候选 ≥95% 碱基以 ≥95% 同一率比对上参照——
/// eval 精度方向的转录本粒度化; hit 本身已要求 ident ≥ 0.95）。
pub fn is_positive(qcov: f64) -> bool {
    qcov >= 0.95
}

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

/// 用法: cargo xtask export-candidates --out-dir <cache根|out目录> --ref <ref.fa> --tsv <out.tsv>
///       [--threads N] [--min-support 1] [--counts <mer_counts.rs.bin>]
pub fn run(argv: &[String]) -> i32 {
    let mut out_dir = String::new();
    let mut reff = String::new();
    let mut tsv = String::new();
    let mut threads = 8usize;
    let mut min_support = 1i64;
    let mut counts_override = String::new();
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "--out-dir" => {
                out_dir = argv.get(i + 1).cloned().unwrap_or_default();
                i += 1;
            }
            "--ref" => {
                reff = argv.get(i + 1).cloned().unwrap_or_default();
                i += 1;
            }
            "--tsv" => {
                tsv = argv.get(i + 1).cloned().unwrap_or_default();
                i += 1;
            }
            "--threads" => {
                threads = argv.get(i + 1).and_then(|v| v.parse().ok()).unwrap_or(8);
                i += 1;
            }
            "--min-support" => {
                min_support = argv.get(i + 1).and_then(|v| v.parse().ok()).unwrap_or(1);
                i += 1;
            }
            "--counts" => {
                counts_override = argv.get(i + 1).cloned().unwrap_or_default();
                i += 1;
            }
            other => {
                eprintln!("未知参数: {other}");
                return 2;
            }
        }
        i += 1;
    }
    if out_dir.is_empty() || reff.is_empty() || tsv.is_empty() || threads == 0 || min_support < 1 {
        eprintln!("用法: cargo xtask export-candidates --out-dir <cache根|out目录> --ref <ref.fa> --tsv <out.tsv> [--threads N] [--min-support 1] [--counts <mer_counts.rs.bin>]");
        return 2;
    }

    let bins = locate_component_bins(Path::new(&out_dir))
        .unwrap_or_else(|| panic!("在 {out_dir} 下找不到 Component_bins"));
    eprintln!("[export-candidates] Component_bins: {}", bins.display());

    // 收集并按组件号排序（文件名 c<N>.graph.out）
    let mut files: Vec<(i32, PathBuf)> = Vec::new();
    for cb in std::fs::read_dir(&bins).unwrap() {
        let cb = cb.unwrap().path();
        if !cb.is_dir() {
            continue;
        }
        for f in std::fs::read_dir(&cb).unwrap() {
            let f = f.unwrap().path();
            let name = f.file_name().unwrap().to_string_lossy().to_string();
            if let Some(rest) = name.strip_suffix(".graph.out") {
                if let Some(idstr) = rest.strip_prefix('c') {
                    if let Ok(id) = idstr.parse::<i32>() {
                        files.push((id, f));
                    }
                }
            }
        }
    }
    files.sort_by_key(|(id, _)| *id);
    eprintln!("[export-candidates] 组件 {} 个", files.len());

    let asm_out = locate_asm_out(Path::new(&out_dir))
        .unwrap_or_else(|| panic!("在 {out_dir} 下找不到 asm out（含 chrysalis/）"));
    eprintln!("[export-candidates] 构建/装载 25-mer 计数表…");
    let t0 = std::time::Instant::now();
    let (counts, provenance) = counts_for_out_dir(&asm_out, &counts_override);
    eprintln!(
        "[export-candidates] 计数表 {} 键（来源: {}, {:.1?}）",
        counts.len(),
        provenance,
        t0.elapsed()
    );

    let refs = read_fasta(Path::new(&reff));
    let ref_lens: HashMap<String, usize> = refs.iter().map(|(n, s)| (n.clone(), s.len())).collect();
    eprintln!("[export-candidates] 参照转录本 {} 条", refs.len());

    // rayon 并行重跑组件（run_component 内部自带大栈线程隔离;
    // 单组件 Err 记失败数继续, 不中断其余）
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap_or_else(|e| panic!("rayon pool: {e}"));
    let failed = AtomicUsize::new(0);
    let done = AtomicUsize::new(0);
    let n_comps = files.len();
    let t0 = std::time::Instant::now();
    let per_comp: Vec<Vec<Candidate>> = pool.install(|| {
        files
            .par_iter()
            .map(|(id, graph_out_path)| {
                let reads_path = graph_out_path.with_file_name(format!("c{id}.graph.reads"));
                let graph_out = std::fs::read_to_string(graph_out_path)
                    .unwrap_or_else(|e| panic!("read {}: {e}", graph_out_path.display()));
                let graph_reads = match std::fs::read_to_string(&reads_path) {
                    Ok(t) => t,
                    Err(e) => {
                        failed.fetch_add(1, Ordering::Relaxed);
                        eprintln!("[export-candidates] 组件 c{id} reads 读取失败: {e}");
                        return Vec::new();
                    }
                };
                let params = ComponentParams {
                    n: 100_000,
                    min_len: 200,
                    max_pair_distance: 10_000,
                    no_em_reduce: true,
                    min_read_support: min_support,
                    path_reinforcement_distance: Some(25),
                    name: format!("c{id}"),
                    ..Default::default()
                };
                let cands = match run_component(&graph_out, &graph_reads, &params) {
                    Ok(res) => parse_candidates(*id, &res.all_prob_paths_fasta),
                    Err(e) => {
                        failed.fetch_add(1, Ordering::Relaxed);
                        eprintln!("[export-candidates] 组件 c{id} 失败: {e}");
                        Vec::new()
                    }
                };
                let d = done.fetch_add(1, Ordering::Relaxed) + 1;
                if d.is_multiple_of(500) || d == n_comps {
                    eprintln!(
                        "[export-candidates] 进度 {d}/{n_comps}（{:.1?}）",
                        t0.elapsed()
                    );
                }
                cands
            })
            .collect()
    });
    let n_failed = failed.load(Ordering::Relaxed);
    eprintln!(
        "[export-candidates] 组件 {} 个, 失败 {}（{:.2}%）, 引擎重跑 {:.1?}",
        n_comps,
        n_failed,
        100.0 * n_failed as f64 / n_comps.max(1) as f64,
        t0.elapsed()
    );

    let cands: Vec<Candidate> = per_comp.into_iter().flatten().collect();
    eprintln!("[export-candidates] 候选 {} 条", cands.len());

    // minimap2 只跑一次: 全部候选写 fasta 后一次比对。
    // fasta 与 PAF 均落盘保留（操作口径评估——按模型分过滤后重算 R_full/P/F1——
    // 需要序列与逐命中信息, TSV 不含）。
    let paf = if cands.is_empty() {
        String::new()
    } else {
        let tmp_fa = format!("{tsv}.candidates.fa");
        let mut fa = String::with_capacity(1 << 22);
        for c in &cands {
            fa.push('>');
            fa.push_str(&c.name);
            fa.push('\n');
            for chunk in c.seq.as_bytes().chunks(60) {
                fa.push_str(&String::from_utf8_lossy(chunk));
                fa.push('\n');
            }
        }
        std::fs::write(&tmp_fa, &fa).unwrap_or_else(|e| panic!("write {tmp_fa}: {e}"));
        let out = std::process::Command::new(crate::bench::MINIMAP2)
            .args(["-x", "asm5", "-t", &threads.to_string(), "--secondary=no"])
            .arg(&reff)
            .arg(&tmp_fa)
            .output()
            .expect("无法启动 minimap2");
        if !out.status.success() {
            eprintln!(
                "minimap2 失败 ({})\n{}",
                out.status,
                String::from_utf8_lossy(&out.stderr)
            );
            return 1;
        }
        let paf = String::from_utf8_lossy(&out.stdout).to_string();
        std::fs::write(format!("{tsv}.paf"), &paf)
            .unwrap_or_else(|e| panic!("write {tsv}.paf: {e}"));
        paf
    };
    let hits = best_hits(&paf, &ref_lens);
    eprintln!("[export-candidates] minimap2 命中候选 {} 条", hits.len());

    // 统一写 TSV（特征全部推理时可得; ref_cov/ref_ident/label 为 oracle 列）
    let mut out = String::with_capacity(1 << 20);
    out.push_str(TSV_HEADER);
    out.push('\n');
    let (mut n_pos, mut s_unk, mut s1, mut s2, mut s3) = (0usize, 0usize, 0usize, 0usize, 0usize);
    for c in &cands {
        let (orf_aa, orf_frac, orf_m) = orf_features(&c.seq);
        let cov = cov_features(&counts, &c.seq);
        let (rc, ri, qc) = hits.get(&c.name).copied().unwrap_or((0.0, 0.0, 0.0));
        let label = if is_positive(qc) { 1 } else { 0 };
        n_pos += label;
        match c.support {
            -1 => s_unk += 1, // classic 引擎无 support 字段
            1 => s1 += 1,
            2 => s2 += 1,
            _ => s3 += 1,
        }
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{:.4}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}",
            c.comp,
            c.name,
            c.len,
            c.support,
            gc_frac(&c.seq),
            entropy(&c.seq),
            max_homopolymer(&c.seq),
            orf_aa,
            orf_frac,
            orf_m,
            cov.mean,
            cov.median,
            cov.min,
            cov.ge2_frac,
            cov.cv,
            rc,
            ri,
            qc,
            label,
        )
        .unwrap();
    }
    std::fs::write(&tsv, &out).unwrap_or_else(|e| panic!("write {tsv}: {e}"));
    eprintln!("[export-candidates] 写出 {tsv}");
    eprintln!(
        "[export-candidates] 正例 {} / {}（{:.2}%）; support 分布: 未知(classic): {}, =1: {}, =2: {}, ≥3: {}",
        n_pos,
        cands.len(),
        100.0 * n_pos as f64 / cands.len().max(1) as f64,
        s_unk,
        s1,
        s2,
        s3
    );
    0
}

// ---------------------------------------------------------------------------
// 测试（手算用例锁定）
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orf_hand_computed() {
        // 框架0: ATG AAA CCG TAA → M K P *, 段 "MKP" 起始于 M, 长 3;
        // 框架2: GAA ACC GTA → E T V 长 3 无 M, 等长平局保留先到者（M 锚定）。
        let (aa, frac, m) = orf_features("ATGAAACCGTAA");
        assert_eq!(aa, 3);
        assert!((frac - 9.0 / 12.0).abs() < 1e-9);
        assert_eq!(m, 1);
        // 无 M: 框架0 AAA CCG TAA → K P *, 段 "KP" 长 2, 不起始于 M
        let (aa, frac, m) = orf_features("AAACCGTAA");
        assert_eq!(aa, 2);
        assert!((frac - 6.0 / 9.0).abs() < 1e-9);
        assert_eq!(m, 0);
        // 段内 M（M 锚定输给更长非 M 段）: 框架0 GAA ATG AAA AAA TAA → E M K K *,
        // 取 "MKK" 长 3; 框架2 AAT GAA AAA AAA → N E K K 长 4 无 M, 更长 → 胜。
        let (aa, frac, m) = orf_features("GAAATGAAAAAATAA");
        assert_eq!(aa, 4);
        assert!((frac - 12.0 / 15.0).abs() < 1e-9);
        assert_eq!(m, 0);
        // 段内 M（M 锚定等长平局胜出, 先到者保留）: 框架0 GAA ATG AAA AAA AAA(残 AA)
        // → E M K K K 取 "MKKK" 长 4; 框架2 AAT GAA AAA AAA TAA → N E K K * 长 4,
        // 等长不覆盖 → 保留框架0 的 M 锚定段。
        let (aa, frac, m) = orf_features("GAAATGAAAAAAAATAA");
        assert_eq!(aa, 4);
        assert!((frac - 12.0 / 17.0).abs() < 1e-9);
        assert_eq!(m, 1);
        // 无完整密码子
        assert_eq!(orf_features("AC"), (0, 0.0, 0));
    }

    #[test]
    fn cov_stats_hand_computed() {
        // 28nt → 4 个 25-mer 窗; 计数依次 4 / 1 / 缺表(=0) / 3。
        let seq = "TGCAGTCCGATGACAATTCCGGTAAGCA";
        assert_eq!(seq.len(), 28);
        let b = seq.as_bytes();
        let mut counts = CountMap::default();
        let want = [4u64, 1, 0, 3];
        let mut keys = Vec::new();
        for (i, w) in b.windows(25).enumerate() {
            let key = get_ds_kmer_val(kmer_to_intval(w).unwrap(), 25);
            assert!(!keys.contains(&key), "测试序列窗口 ds 键冲突");
            keys.push(key);
            if want[i] > 0 {
                counts.insert(key, want[i] as u32);
            }
        }
        let st = cov_features(&counts, seq);
        // vals = [4,1,0,3]: mean 2.0; sorted [0,1,3,4] median (1+3)/2 = 2.0;
        // min 0; ≥2 的窗 2/4 = 0.5; var = (4+1+4+1)/4 = 2.5, cv = sqrt(2.5)/2
        assert!((st.mean - 2.0).abs() < 1e-9);
        assert!((st.median - 2.0).abs() < 1e-9);
        assert_eq!(st.min, 0.0);
        assert!((st.ge2_frac - 0.5).abs() < 1e-9);
        assert!((st.cv - 2.5f64.sqrt() / 2.0).abs() < 1e-9);
        // 含 N 的窗跳过: 中部插入 N → 含 N 的 25 个窗剔除, 仅剩首尾无 N 窗
        let seq_n = format!("{}N{}", &seq[..14], &seq[14..]);
        let st_n = cov_features(&counts, &seq_n);
        // 29nt, 5 个窗, 全部含 N? 否: 窗 [0..25) 与 [4..29) 含 N? 窗0=[0..25) 含位置14
        // 的 N → 跳过; 同理窗 0..=4 全跨位置 14? 窗 i 覆盖 [i, i+25), i∈0..=4,
        // 位置 14 落在每个窗内 → 全部跳过
        assert_eq!(st_n, CovStats::default());
        // 短于 25 → 全 0
        assert_eq!(cov_features(&counts, "ACGT"), CovStats::default());
    }

    #[test]
    fn parse_candidates_fasta() {
        let fasta = concat!(
            ">c7_g1_i1 len=70 support=5 path=[1,2,3]\n",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
            "GGGGGGGGGG\n",
            ">c7_g1_i2 support=2 path=[1,4]\n",
            "TTTTTTTTTT\n",
        );
        let c = parse_candidates(7, fasta);
        assert_eq!(c.len(), 2);
        assert_eq!(c[0].comp, 7);
        assert_eq!(c[0].name, "c7_g1_i1");
        assert_eq!(c[0].len, 70);
        assert_eq!(c[0].support, 5);
        assert_eq!(c[0].seq.len(), 70);
        assert_eq!(&c[0].seq[..30], &"A".repeat(30));
        // header 缺 len= → 回退序列长
        assert_eq!(c[1].name, "c7_g1_i2");
        assert_eq!(c[1].len, 10);
        assert_eq!(c[1].support, 2);
        assert_eq!(c[1].seq, "TTTTTTTTTT");
    }

    #[test]
    fn paf_best_hit_and_label() {
        // PAF 列: q qlen qs qe strand t tlen ts te mch aln mapq
        let paf = concat!(
            "q1\t1000\t0\t950\t+\tt1\t1000\t0\t950\t940\t1000\t60\n",
            "q1\t1000\t0\t1000\t+\tt2\t2000\t0\t1000\t990\t1000\t60\n",
            "q2\t500\t0\t480\t+\tt3\t500\t0\t480\t470\t480\t60\n",
            "q3\t500\t0\t500\t+\tt1\t1000\t0\t500\t400\t500\t60\n",
        );
        let ref_lens: HashMap<String, usize> =
            [("t1".into(), 1000), ("t2".into(), 2000), ("t3".into(), 500)]
                .into_iter()
                .collect();
        let hits = best_hits(paf, &ref_lens);
        // q1: ref 侧 best = hit1（cov .95 ident .94, score .893 > hit2 .495）;
        // qcov: hit2 ident .99 ≥ .95, span 1000/1000 = 1.0 → 正例
        let (rc, ri, qc) = hits["q1"];
        assert!((rc - 0.95).abs() < 1e-9);
        assert!((ri - 0.94).abs() < 1e-9);
        assert!((qc - 1.0).abs() < 1e-9);
        assert!(is_positive(qc));
        // q2: qcov 480/500 = 0.96 → 正例; ref_cov 0.96
        let (rc, _, qc) = hits["q2"];
        assert!((rc - 0.96).abs() < 1e-9);
        assert!((qc - 0.96).abs() < 1e-9);
        assert!(is_positive(qc));
        // q3: ident 0.8 < 0.95 → qcov 0 → 负例（ref 侧 best 仍记录）
        let (rc, ri, qc) = hits["q3"];
        assert!((rc - 0.5).abs() < 1e-9);
        assert!((ri - 0.8).abs() < 1e-9);
        assert_eq!(qc, 0.0);
        assert!(!is_positive(qc));
        // q0 无命中 → 表缺 → 调用侧缺省 (0,0,0) 负例
        assert!(!hits.contains_key("q0"));
        assert!(!is_positive(0.0));
        // 阈值边界
        assert!(is_positive(0.95));
        assert!(!is_positive(0.9499));
    }
}
