//! Chrysalis 侧 DNAVector / vecDNAVector / DNAStringStreamFast 的最小语义移植
//! （P3-T1 基础层）。
//!
//! 镜像来源：
//! - `Chrysalis/analysis/DNAVector.h:67-118`（plain_table / NucIndex —— A=0,C=1,G=2,T=3,N=4）
//! - `Chrysalis/analysis/DNAVector.cc:856-971`（vecDNAVector::Read 的 FASTA 读入语义）
//! - `Chrysalis/analysis/sequenceUtil.cc:124-177`（revcomp）、`:326-355`（compute_entropy(string)）
//! - `Chrysalis/analysis/GraphFromFasta.cc:70-160`（is_simple_repeat）、
//!   `:197-236`（IsSimple）、`:239-276`（SimpleHalves）
//! - `Chrysalis/analysis/DNAVector.cc:1456-1501`（DNAStringStreamFast 流式读）
//!
//! 已证差异（有意为之，勿"修复"）：
//! - CRLF：本版剥行尾 `\r`（原版 getline 保留，Chrysalis 读自身产物均为 LF）；
//! - 原版 `is_simple_repeat` 的 stringstream 累积 bug（`ss.clear()` 不清内容，
//!   DEBUG 输出的 best_left/right_kmer 越滚越长）不复刻——只影响调试打印；
//! - `bases_to_number` 对 `from+range > seq.len()` 返回 -1（原版越界读是 UB；
//!   调用方循环边界本就保证不越界）。

use std::path::Path;

use trinity_common::error::CommonError;

/// GraphFromFasta.cc:21 `static float MIN_KMER_ENTROPY = 1.3;`
pub const MIN_KMER_ENTROPY: f32 = 1.3;

/// GraphFromFasta.cc:24 `static float MAX_RATIO_INTERNALLY_REPETITIVE = 0.85;`
pub const MAX_RATIO_INTERNALLY_REPETITIVE: f32 = 0.85;

/// name 含 '>' 前缀（原版 `vecDNAVector::Name()` 保留 '>'，`NameClean()` 才去掉）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DnaSeq {
    pub name: String,
    pub seq: Vec<u8>,
}

/// 镜像 `vecDNAVector::Read(file, false, shortName=false, allUpper=true, ...)`。
///
/// - header：以首个 token 首字符为 '>' 判定；shortName=false 时**全部**空白分隔
///   token 用 '_' 连接（`>a1;43 total_counts: 123` → `>a1;43_total_counts:_123`）；
/// - 序列行：原版只取行内**第一个**空白 token（DNAVector.cc:938 `AsString(0)`），
///   多行拼接无分隔，读毕整体 toupper；
/// - 首条 header 之前的裸序列行并入第一条记录（镜像原版 tmpVec/j 状态机）。
pub fn read_fasta(path: impl AsRef<Path>) -> Result<Vec<DnaSeq>, CommonError> {
    read_fasta_impl(path, false)
}

/// 镜像 `Read(file, false, shortName=true, allUpper=true, ...)`（QuantifyGraph.cc:346
/// `seq.Read(aString, false, true, true, 1000)`）——name = header 首 token（含 '>'）。
pub fn read_fasta_short_names(path: impl AsRef<Path>) -> Result<Vec<DnaSeq>, CommonError> {
    read_fasta_impl(path, true)
}

fn read_fasta_impl(path: impl AsRef<Path>, short_name: bool) -> Result<Vec<DnaSeq>, CommonError> {
    let data = std::fs::read(path)?;
    Ok(read_fasta_bytes(&data, short_name))
}

/// 字节版 vecDNAVector::Read（管线库层入参用：内存中的 fasta 文本）。
pub fn read_fasta_bytes(data: &[u8], short_name: bool) -> Vec<DnaSeq> {
    let mut out: Vec<DnaSeq> = Vec::new();
    let mut cur: Vec<u8> = Vec::new(); // tmpVec：未 flush 的序列累计
    let mut last_name = String::new();
    let mut active = false; // pVec != NULL：已见过至少一条 header

    for raw_line in data.split(|&b| b == b'\n') {
        let line = strip_cr(raw_line);
        // FlatFileParser 语义：token = 空白(' ' 与 '\t')分隔的非空段（mutil.cc Tokenize）
        let mut tokens = line
            .split(|&c| c == b' ' || c == b'\t')
            .filter(|t| !t.is_empty());
        let first = match tokens.next() {
            Some(t) => t,
            None => continue, // GetItemCount()==0 → 跳过空行
        };
        if first.first() == Some(&b'>') {
            if active {
                out.push(DnaSeq {
                    name: std::mem::take(&mut last_name),
                    seq: std::mem::take(&mut cur).to_ascii_uppercase(),
                });
            }
            // shortName=false：余下 token 全部 '_' 连接
            let mut name = String::from_utf8_lossy(first).into_owned();
            if !short_name {
                for t in tokens {
                    name.push('_');
                    name.push_str(&String::from_utf8_lossy(t));
                }
            }
            last_name = name;
            active = true;
        } else {
            // 序列行：原版只取第一个 token（quirk，见模块文档）
            cur.extend_from_slice(first);
        }
    }

    if active {
        out.push(DnaSeq {
            name: last_name,
            seq: cur.to_ascii_uppercase(),
        });
    }
    out
}

fn strip_cr(line: &[u8]) -> &[u8] {
    match line.split_last() {
        Some((b'\r', rest)) => rest,
        _ => line,
    }
}

/// DNAStringStreamFast（DNAVector.cc:1456-1501）等价的批量流式读：
/// - 首条 '>' 之前的行丢弃（ReadStream 先定位 header）；
/// - header 行丢弃，多行序列**无分隔**拼接，不做 toupper（调用方做）；
/// - **空序列记录终止整个流**（原版 `Next()` 返回空串 → AddData 的
///   `while(true)` 立即 break，后续记录全部不读）——按原样复刻。
pub fn stream_fasta_records(data: &str) -> Vec<Vec<u8>> {
    let mut out: Vec<Vec<u8>> = Vec::new();
    let mut cur: Vec<u8> = Vec::new();
    let mut in_record = false;
    for raw in data.split('\n') {
        let line = raw.strip_suffix('\r').unwrap_or(raw);
        if line.as_bytes().first() == Some(&b'>') {
            if in_record {
                if cur.is_empty() {
                    return out; // 空记录终止流（原版 quirk）
                }
                out.push(std::mem::take(&mut cur));
            }
            in_record = true;
        } else if in_record {
            cur.extend_from_slice(line.as_bytes());
        }
    }
    if in_record && !cur.is_empty() {
        out.push(cur);
    }
    out
}

/// sequenceUtil.cc:124-177 `revcomp(string)`：A<->T、C<->G（大小写各自互补），
/// 其余（含 N 与非法字符）→ 'N'。
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&c| match c {
            b'g' => b'c',
            b'G' => b'C',
            b'a' => b't',
            b'A' => b'T',
            b't' => b'a',
            b'T' => b'A',
            b'c' => b'g',
            b'C' => b'G',
            _ => b'N',
        })
        .collect()
}

/// DNAVector.h:67-84 plain_table / :115-118 `NucIndex`：
/// **A=0, C=1, G=2, T=3, N=4，其余（含小写）=-1**。
///
/// 与 Inchworm 侧 sequenceUtil 的 `_base_to_int`（G=0,A=1,T=2,C=3，
/// 即 trinity_common::kmer::base_to_int）完全不同——锁定差异的测试见下。
#[inline]
pub fn nuc_index(letter: u8) -> i32 {
    match letter {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        b'N' => 4,
        _ => -1,
    }
}

/// KmerAlignCore.h:46-64 `TranslateBasesToNumberExact::BasesToNumber`：
/// 小端编码 `Σ enc(seq[from+t]) << 2t`（首碱基权重最低位）；
/// 任一碱基 nuc_index 为 -1（非法/小写）或 4（N）→ 返回 -1。
/// `from + range > seq.len()` 返回 -1（原版越界读 UB，见模块文档）。
pub fn bases_to_number(seq: &[u8], from: usize, range: usize) -> i64 {
    if from + range > seq.len() {
        return -1;
    }
    let mut num: i64 = 0;
    let mut shift: i64 = 1;
    for t in 0..range {
        let v = nuc_index(seq[from + t]);
        if v == -1 || v == 4 {
            return -1;
        }
        num += (v as i64) * shift;
        shift <<= 2;
    }
    num
}

/// Chrysalis 侧字符串熵（sequenceUtil.cc:326-355 / GraphFromFasta.cc:166-195，
/// 两份实现逐行同型）：只统计大写 G/A/T/C 四种计数，**分母是全长**（N、小写
/// 等都计入分母但不进分子），`entropy = Σ p·log2(1/p)`。
/// 空串 / 全非 GATC → 0（prob=NaN 不入分支，与原版一致）。
pub fn compute_entropy(kmer: &[u8]) -> f32 {
    let mut entropy = 0f32;
    let len = kmer.len();
    for nuc in [b'G', b'A', b'T', b'C'] {
        let count = kmer.iter().filter(|&&c| c == nuc).count();
        let prob = count as f32 / len as f32;
        if prob > 0. {
            // C++: float val = prob * log(1/prob)/log(2.0f); —— 中间按 double 计
            let val = (prob as f64 * ((1.0 / prob as f64).ln() / std::f64::consts::LN_2)) as f32;
            entropy += val;
        }
    }
    entropy
}

/// GraphFromFasta.cc:197-214 `IsSimple`：`compute_entropy(d) < 1.3`（**严格 <**）。
pub fn is_simple(kmer: &[u8]) -> bool {
    compute_entropy(kmer) < MIN_KMER_ENTROPY
}

/// GraphFromFasta.cc:70-160 `is_simple_repeat`：mid = len/2；i∈[0,mid)、
/// j∈(i,mid]，比较 `kmer[i..i+mid)` 与 `kmer[j..j+mid)`（各 mid 位）的相等率；
/// 早退 `ratio >= 0.85`（**含等号**）、终判 `max > 0.85`（**严格**）。
/// 不复刻 stringstream 累积 bug（仅影响 DEBUG 输出）。
pub fn is_simple_repeat(kmer: &[u8]) -> bool {
    let mid = kmer.len() / 2;
    let mut max_ratio = 0f32;
    for i in 0..mid {
        for j in i + 1..=mid {
            let mut bases_compared = 0;
            let mut bases_common = 0;
            for t in 0..mid {
                bases_compared += 1;
                if kmer[i + t] == kmer[j + t] {
                    bases_common += 1;
                }
            }
            let ratio_same = bases_common as f32 / bases_compared as f32;
            if ratio_same > max_ratio {
                max_ratio = ratio_same;
            }
            if ratio_same >= MAX_RATIO_INTERNALLY_REPETITIVE {
                return true; // !DEBUG 分支：早退
            }
        }
    }
    max_ratio > MAX_RATIO_INTERNALLY_REPETITIVE
}

/// GraphFromFasta.cc:239-276 `SimpleHalves`（DISABLE_REPEAT_CHECK=false 默认语义）：
/// 左半 `[0,len/2)`、右半 `[len/2,len)`；任一半低熵（<1.3）或内部重复即 true。
pub fn simple_halves(s: &[u8]) -> bool {
    simple_halves_with(s, false)
}

/// 带 `-disable_repeat_check` 开关的完整语义（C++ 运算符优先级
/// `(!D && isr(left)) || isr(right)` 在 D=false 下化简为 `isr(left) || isr(right)`）。
pub fn simple_halves_with(s: &[u8], disable_repeat_check: bool) -> bool {
    let len = s.len();
    let mid_pos = len / 2;
    let left = &s[..mid_pos];
    let right = &s[mid_pos..];
    compute_entropy(left) < MIN_KMER_ENTROPY
        || compute_entropy(right) < MIN_KMER_ENTROPY
        || ((!disable_repeat_check && is_simple_repeat(left)) || is_simple_repeat(right))
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---------- nuc_index / bases_to_number：两套编码表差异锁定 ----------

    /// Chrysalis 表 A=0,C=1,G=2,T=3,N=4 —— 注意 'A' 是 0（Inchworm 表里 A=1）。
    #[test]
    fn nuc_index_chrysalis_table() {
        assert_eq!(nuc_index(b'A'), 0);
        assert_eq!(nuc_index(b'C'), 1);
        assert_eq!(nuc_index(b'G'), 2);
        assert_eq!(nuc_index(b'T'), 3);
        assert_eq!(nuc_index(b'N'), 4);
        // 其余（含小写、X、'-'）均为 -1
        for c in [b'a', b'c', b'g', b't', b'n', b'X', b'-', b'*'] {
            assert_eq!(nuc_index(c), -1, "char {}", c as char);
        }
    }

    /// 与 trinity-common（Inchworm 表 G=0,A=1,T=2,C=3）逐字母对拍——除 'T' 外
    /// 数值全部不同，防止两套表混用。
    #[test]
    fn nuc_index_differs_from_inchworm_table() {
        let inchworm = |c: u8| trinity_common::kmer::base_to_int(c).map(i32::from);
        assert_eq!(inchworm(b'A'), Some(1));
        assert_eq!(inchworm(b'G'), Some(0));
        assert_eq!(inchworm(b'T'), Some(2));
        assert_eq!(inchworm(b'C'), Some(3));
        for c in [b'A', b'C', b'G', b'T'] {
            assert_ne!(
                nuc_index(c),
                inchworm(c).unwrap(),
                "两套编码表在此字母上不应相等: {}",
                c as char
            );
        }
    }

    /// 小端权重：首碱基在最低 2bit。"ACGT" = 0 + 1·4 + 2·16 + 3·64 = 228。
    #[test]
    fn bases_to_number_little_endian() {
        assert_eq!(bases_to_number(b"ACGT", 0, 4), 4 + 32 + 192); // A(0)·1+C(1)·4+G(2)·16+T(3)·64
        assert_eq!(bases_to_number(b"AG", 0, 2), 2 * 4); // A(0)·1 + G(2)·4
        assert_eq!(bases_to_number(b"GA", 0, 2), 2); // G(2)·1 + A(0)·4 → 交换后不同
        assert_eq!(bases_to_number(b"AAAA", 0, 4), 0);
        assert_eq!(bases_to_number(b"TTTT", 0, 4), 3 * (1 + 4 + 16 + 64));
    }

    #[test]
    fn bases_to_number_miss_cases() {
        assert_eq!(bases_to_number(b"ACGN", 0, 4), -1); // N → -1
        assert_eq!(bases_to_number(b"acgt", 0, 4), -1); // 小写 → -1
        assert_eq!(bases_to_number(b"ACGT", 1, 4), -1); // 越界窗口 → -1
        assert_eq!(bases_to_number(b"ACGT", 0, 5), -1);
        // from 偏移正常工作：[1..4) = "CGT" = 1 + 2·4 + 3·16 = 57
        assert_eq!(bases_to_number(b"ACGT", 1, 3), 1 + 8 + 48);
    }

    // ---------- revcomp ----------

    #[test]
    fn revcomp_basics() {
        assert_eq!(revcomp(b"AACGT"), b"ACGTT");
        assert_eq!(revcomp(b"A"), b"T");
        assert_eq!(revcomp(b"ACGT"), b"ACGT"); // 自互补
        assert_eq!(revcomp(b""), b"");
    }

    /// 大小写各自互补（原版 switch 显式列出小写分支）；其余字符 → 'N'。
    #[test]
    fn revcomp_case_and_unknown() {
        assert_eq!(revcomp(b"aacg"), b"cgtt");
        assert_eq!(revcomp(b"AN"), b"NT"); // N → N
        assert_eq!(revcomp(b"XYZ"), b"NNN"); // 非法 → N
    }

    // ---------- compute_entropy / is_simple ----------

    #[test]
    fn entropy_golden_values() {
        assert_eq!(compute_entropy(b"AAAAAAAAAAAA"), 0.0); // poly-A
        assert!((compute_entropy(b"ACGTACGTACGT") - 2.0).abs() < 1e-5); // 等比四元
                                                                        // 含 N：分母全长（10）但分子只算 A（4）→ 0.4·log2(2.5) = 0.52877
        assert!((compute_entropy(b"AAAANNNNNN") - 0.528_771_2).abs() < 1e-5);
        // A=8,T=2：0.8·log2(1.25)+0.2·log2(5) = 0.72193
        assert!((compute_entropy(b"AAAATAAAAT") - 0.721_928_1).abs() < 1e-5);
    }

    /// 只统计大写 GATC：小写全不入分子；空串/全非 GATC → 0。
    #[test]
    fn entropy_denominator_quirks() {
        assert_eq!(compute_entropy(b"acgtacgt"), 0.0); // 小写不计分子
        assert_eq!(compute_entropy(b""), 0.0);
        assert_eq!(compute_entropy(b"NNNN"), 0.0);
        // 混大小写：只有 4 个大写 A 计分子，分母 12 → p=1/3 → log2(3)/3
        assert!((compute_entropy(b"AAAAnnnnaaaa") - 0.528_320_3).abs() < 1e-5);
    }

    #[test]
    fn is_simple_boundary_semantics() {
        // 熵 0 / 0.72 / 0.92 → 严格 < 1.3 → simple
        assert!(is_simple(b"AAAAAAAAAAAA"));
        assert!(is_simple(b"AAAATAAAAT")); // 0.7219
        assert!(is_simple(b"AAAACC")); // 0.9183
                                       // 熵 1.459 / 1.585 / 2.0 → not simple
        assert!(!is_simple(b"AAACCG")); // 1.4591
        assert!(!is_simple(b"CAAGGC")); // log2(3) = 1.585
        assert!(!is_simple(b"ACGTACGTACGT")); // 2.0
    }

    // ---------- is_simple_repeat ----------

    #[test]
    fn is_simple_repeat_poly_and_random() {
        assert!(is_simple_repeat(b"AAAAAAAAAA")); // 任意窗口全等 → 1.0
        assert!(!is_simple_repeat(b"ACTG")); // mid=2 各窗口 0 相等
        assert!(!is_simple_repeat(b"AC")); // mid=1：A vs C → 0.0
        assert!(is_simple_repeat(b"AA")); // 1.0 ≥ 0.85
        assert!(!is_simple_repeat(b"")); // 无窗口
    }

    /// 手推全部 15 个 (i,j) 对（mid=5）："AAAATAAAAG"（周期 5 的 poly-A，
    /// 末端 T/G 破坏对齐全等）：逐对比值——(0,5) 与 (1,5) 均 4/5 = 0.8，
    /// 其余 13 对均为 3/5 = 0.6 —— max = 0.8 < 0.85 → false。
    /// （对照：末位也是 T 的 "AAAATAAAAT" 周期 5 完整 → (0,5) 全等 → true。）
    #[test]
    fn is_simple_repeat_max_exactly_0p8_is_not_simple() {
        assert!(!is_simple_repeat(b"AAAATAAAAG"));
        assert!(is_simple_repeat(b"AAAATAAAAT")); // 对齐全等 → ≥ 0.85
    }

    /// 恰等 0.85（17/20）：早退分支用 `>=`，应判 simple。
    /// 构造：X=20×A，Y=X 仅在 t=3,8,18 置 T → (i=0,j=20) 两半窗口 17/20 相等。
    #[test]
    fn is_simple_repeat_exactly_0p85_is_simple() {
        let x = b"A".repeat(20);
        let mut y = x.clone();
        y[3] = b'T';
        y[8] = b'T';
        y[18] = b'T';
        let mut kmer = x;
        kmer.extend_from_slice(&y);
        assert!(is_simple_repeat(&kmer));
    }

    /// 终判分支 `max > 0.85` 在 !DEBUG 路径下实为不可达（任何 ratio ≥ 0.85
    /// 都先被早退拦截）——严格性已由上面的 0.8 案例与恰 0.85 案例共同锁定，
    /// 无需单测。另注：周期性序列（如 "AAAAT"×2 对齐偏移 5 的窗口全等）必为
    /// true，构造 false 案例时须避开小周期对齐。

    // ---------- simple_halves ----------

    #[test]
    fn simple_halves_golden_cases() {
        // 左半 "AAAAT" 熵 0.722 < 1.3 → true
        assert!(simple_halves(b"AAAATAAAAT"));
        // "ACGTTG|CAAGGC"：左熵 1.918、右熵 1.585；两半 isr 均 ≤ 2/3 → false
        assert!(!simple_halves(b"ACGTTGCAAGGC"));
        // 左半 poly-A → true（右半无论如何）
        assert!(simple_halves(b"AAAAAACGTACG"));
        // disable_repeat_check 只关掉 isr 分支，熵分支仍在
        let s = b"AAAATAAAAT";
        assert!(simple_halves_with(s, true)); // 左半熵仍 < 1.3
    }

    /// 左右两半熵均 ≥ 1.3 但左半内部重复 → true；开关只护左半（C++ 优先级
    /// quirk：`(!D && isr(left)) || isr(right)`——isr(right) 永不受开关控制）。
    #[test]
    fn simple_halves_repeat_branch() {
        let left = b"ACGACG"; // 周期 3 → isr=true；熵 log2(3)=1.585 ≥ 1.3
        let right = b"TTGCAA"; // isr=false（逐对手推 max=1/3）；熵 1.918 ≥ 1.3
        assert!(is_simple_repeat(left));
        assert!(!is_simple_repeat(right));
        assert!(compute_entropy(left) >= 1.3);
        assert!(compute_entropy(right) >= 1.3);
        let mut s = left.to_vec();
        s.extend_from_slice(right);
        assert!(simple_halves_with(&s, false)); // 左半 isr 触发
        assert!(!simple_halves_with(&s, true)); // 开关关掉左半 isr，右半本就干净
    }

    /// 反向验证优先级 quirk：右半 isr 不受开关控制——关掉开关仍为 true。
    #[test]
    fn simple_halves_right_repeat_ignores_flag() {
        let left = b"TTGCAA"; // 干净
        let right = b"ACGACG"; // isr=true
        let mut s = left.to_vec();
        s.extend_from_slice(right);
        assert!(simple_halves_with(&s, false));
        assert!(simple_halves_with(&s, true)); // isr(right) 仍生效
    }

    // ---------- read_fasta ----------

    fn write_tmp(name: &str, content: &[u8]) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join("trinity_chrysalis_dnav_t1");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join(name);
        std::fs::write(&p, content).unwrap();
        p
    }

    /// 多 token header 的 '_' 连接、多行序列拼接、大写化。
    #[test]
    fn read_fasta_multi_token_header_and_uppercase() {
        let p = write_tmp(
            "multi.fa",
            b">a1;43 total_counts: 123\nacgt\nTTTT\n>b2\nACGTNN\n",
        );
        let recs = read_fasta(&p).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, ">a1;43_total_counts:_123");
        assert_eq!(recs[0].seq, b"ACGTTTTT");
        assert_eq!(recs[1].name, ">b2");
        assert_eq!(recs[1].seq, b"ACGTNN");
        let _ = std::fs::remove_file(&p);
    }

    /// shortName=true：name = 首 token（含 '>'），余 token 丢弃。
    #[test]
    fn read_fasta_short_names_first_token_only() {
        let p = write_tmp(
            "short.fa",
            b">a1;43 total_counts: 123\nACGT\n>b2 x y\nTTTT\n",
        );
        let recs = read_fasta_short_names(&p).unwrap();
        assert_eq!(recs[0].name, ">a1;43");
        assert_eq!(recs[1].name, ">b2");
        assert_eq!(recs[1].seq, b"TTTT");
        let _ = std::fs::remove_file(&p);
    }

    /// 空行跳过；首条 header 前的裸序列行并入首条记录（原版 tmpVec 状态机）。
    #[test]
    fn read_fasta_leading_and_blank_lines() {
        let p = write_tmp("blank.fa", b"GGGG\n\n>a\nACGT\n\nTT\n>b\nAAA\n");
        let recs = read_fasta(&p).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].seq, b"GGGGACGTTT"); // 领先行并入
        assert_eq!(recs[1].seq, b"AAA");
        let _ = std::fs::remove_file(&p);
    }

    /// 序列行只取第一个空白 token（DNAVector.cc:938 quirk）。
    #[test]
    fn read_fasta_sequence_line_first_token_only() {
        let p = write_tmp("tok.fa", b">a\nACGT TTTT GG\n");
        let recs = read_fasta(&p).unwrap();
        assert_eq!(recs[0].seq, b"ACGT");
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn read_fasta_missing_file_is_error() {
        assert!(read_fasta("/nonexistent/zz.fa").is_err());
    }

    // ---------- stream_fasta_records ----------

    /// header 丢弃、多行无分隔拼接、不做 toupper；首 header 前的行丢弃。
    #[test]
    fn stream_fasta_records_basic() {
        let data = "GGG\n>ignored header\nacgt\nCCCC\n\n>two\nGG\n";
        let recs = stream_fasta_records(data);
        assert_eq!(recs, vec![b"acgtCCCC".to_vec(), b"GG".to_vec()]);
    }

    /// 空序列记录终止整个流（原版 Next()=="" → break）。
    #[test]
    fn stream_fasta_records_empty_record_terminates() {
        let data = ">a\nACGT\n>b\n>c\nTTTT\n";
        let recs = stream_fasta_records(data);
        assert_eq!(recs, vec![b"ACGT".to_vec()]); // >b 的空记录截断后续
    }

    /// EOF 收尾：最后一条非空记录照常产出；末尾空 header 不产记录。
    #[test]
    fn stream_fasta_records_eof_cases() {
        assert_eq!(stream_fasta_records(">a\nAAAA"), vec![b"AAAA".to_vec()]);
        assert_eq!(stream_fasta_records(">a\nAAAA\n>b"), vec![b"AAAA".to_vec()]);
        assert_eq!(stream_fasta_records(""), Vec::<Vec<u8>>::new());
        assert_eq!(stream_fasta_records("AAAA\nCCCC"), Vec::<Vec<u8>>::new()); // 无 header
    }
}
