//! fastaToKmerCoverageStats 移植（Inchworm/src/fastaToKmerCoverageStats.cpp）。
//! 读 jellyfish dump（">count\nKMER"）建计数表 + reads FASTA → 每 read 一行覆盖统计 TSV。
//!
//! 与 counter::for_each_kmer 的关键语义差异（勿合并）:
//! 计数器对含非 gatc 字符的窗口**跳过**（不回调）; 本模块镜像 cpp L300-335——
//! **每个**滑动窗口（共 len-k+1 个）都产生一个覆盖值: 窗口全 gatc → 查表
//! （DS 用 canonical 键，查不到按 1）; 窗口含非 gatc（N 等）→ 计数按 1。
//!
//! 公式逐字镜像（cpp L337-347/L371-387/L389-402）: median = 无符号整数截断的
//! 中位数; mean = long 求和转 f32 除法; stdev = 样本标准差 `sqrtf(Σ(d²)/(n-1))`
//! （n=0 → -0, n=1 → NaN，x86 默认 QNaN 符号位为负 → 原版输出 "-nan"）。

use std::io::{BufReader, Write};

use trinity_common::error::CommonError;
use trinity_common::fasta::FastaReader;
use trinity_common::kmer::{base_to_int, get_ds_kmer_val, kmer_to_intval, KmerId};

use crate::counter::CountMap;

/// 一行统计（cpp L140-148 的 stats_text 内容，tid 列由 write_stats_tsv 固定为 thread:0）。
#[derive(Debug, Clone, PartialEq)]
pub struct CoverageStatsRow {
    pub acc: String,
    pub median: u32,
    pub mean: f32,
    pub stdev: f32,
}

/// cpp L337-347 median_coverage: 排序取中位; 偶数个 (a+b)/2 无符号截断。
/// wrapping_add 镜像 C 无符号回绕（计数实际不会溢出，仅为防御语义一致）。
pub fn median_cov(vals: &mut [u32]) -> u32 {
    let n = vals.len();
    if n == 0 {
        return 0;
    }
    vals.sort_unstable();
    if n % 2 == 1 {
        vals[n / 2]
    } else {
        vals[(n - 1) / 2].wrapping_add(vals[n / 2]) / 2
    }
}

/// cpp L371-387 sum/mean: long 求和后 `(float)sum / size`。
pub fn mean_f32(vals: &[u32]) -> f32 {
    if vals.is_empty() {
        return 0.0;
    }
    let sum: i64 = vals.iter().map(|&v| v as i64).sum();
    sum as f32 / vals.len() as f32
}

/// cpp L389-402 stDev: 样本标准差（除 n-1）; n=0 → -0.0, n=1 → NaN。
pub fn stdev_f32(vals: &[u32]) -> f32 {
    let n = vals.len();
    let avg = mean_f32(vals);
    let mut sum_avg_diffs_sqr = 0.0f32;
    for &v in vals {
        let delta = v as f32 - avg;
        sum_avg_diffs_sqr += delta * delta;
    }
    (sum_avg_diffs_sqr / ((n as i64 - 1) as f32)).sqrt()
}

/// cpp L300-335 compute_kmer_coverage。`seq` 应已过 Fasta_reader 预处理
/// （大写 + 去空白）; 小写 gatc 仍按同码接受，其余字符使所在窗口计数按 1。
/// 序列长度 < k 时返回空向量（并对 stderr 发原版逐字警告——seq 与 "is" 间无空格）。
pub fn kmer_coverage_vector(seq: &[u8], counts: &CountMap, k: usize, ds: bool) -> Vec<u32> {
    if seq.len() < k {
        eprintln!(
            "Sequence: {}is smaller than {} base pairs, skipping",
            String::from_utf8_lossy(seq),
            k
        );
        return Vec::new();
    }
    let mask: KmerId = if k >= 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k as u32)) - 1
    };
    let mut coverage = Vec::with_capacity(seq.len() - k + 1);
    let mut kmer: KmerId = 0;
    let mut valid = 0usize; // 当前位置为止连续 gatc 计数; < k ⇒ 窗口含非 gatc
    for (i, &c) in seq.iter().enumerate() {
        match base_to_int(c) {
            Some(v) => {
                kmer = ((kmer << 2) | v as KmerId) & mask;
                if valid < k {
                    valid += 1;
                }
            }
            None => {
                kmer = 0;
                valid = 0;
            }
        }
        if i + 1 >= k {
            // i+1 = 窗口右端（含）⇒ 窗口 [i+1-k, i]
            let count = if valid == k {
                let key = if ds { get_ds_kmer_val(kmer, k) } else { kmer };
                counts.get(&key).copied().unwrap_or(0)
            } else {
                0 // cpp: 非 gatc 窗口保持 0 → 下面 clamp 到 1
            };
            coverage.push(count.max(1));
        }
    }
    coverage
}

/// cpp L181-228 populate_kmer_counter_from_kmers: 解析 jellyfish dump FASTA。
/// 长度 != k 的条目 stderr 警告并跳过; 首个空序列条目终止（原版 `break`）;
/// 同键计数**累加**（add_kmer 是 `map[kmer] += count`）; DS 键 canonical 化存储。
/// 计数取 header 的 atoi 语义（前导空白/符号/数字，遇非数字截止）。
pub fn load_kmer_dump(data: &[u8], k: usize, ds: bool) -> Result<CountMap, CommonError> {
    let mut counts: CountMap = Default::default();
    let mut reader = FastaReader::new(BufReader::new(data));
    while let Some(rec) = reader.next_record()? {
        if rec.sequence.is_empty() {
            break; // 原版: get_sequence() == "" → break
        }
        if rec.sequence.len() != k {
            eprintln!("ERROR: kmer {} is not of length: {}", rec.sequence, k);
            continue;
        }
        let key = kmer_to_intval(rec.sequence.as_bytes())?; // 非 gatc → Err（原版 throw→exit 1）
        let key = if ds { get_ds_kmer_val(key, k) } else { key };
        *counts.entry(key).or_insert(0) += atoi_prefix(&rec.header);
    }
    Ok(counts)
}

/// C `atoi`（std::atoi）镜像: 跳过前导空白，可选正负号，十进制数字截止。
/// 溢出时 C 为 UB，此处按 i32 饱和（防御，实际 dump 计数远小于此）。
fn atoi_prefix(s: &str) -> u32 {
    let b = s.as_bytes();
    let mut i = 0;
    while i < b.len() && b[i].is_ascii_whitespace() {
        i += 1;
    }
    let neg = i < b.len() && (b[i] == b'+' || b[i] == b'-');
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let mut val: i64 = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        val = (val.saturating_mul(10))
            .saturating_add((b[i] - b'0') as i64)
            .min(i32::MAX as i64);
        i += 1;
    }
    let signed = if neg { -val } else { val };
    signed as i32 as u32 // int → unsigned int 按位回绕（-1 → u32::MAX）
}

/// cpp L122-172 主循环（单线程形态）: 逐 read 计算统计; 空 seq 记录跳过（continue）。
/// 原版对畸形 FASTA throw → exit(1)，此处 Err 对应致命失败。
pub fn coverage_stats_rows(
    reads_fa: &[u8],
    counts: &CountMap,
    k: usize,
    ds: bool,
) -> Result<Vec<CoverageStatsRow>, CommonError> {
    let mut rows = Vec::new();
    let mut reader = FastaReader::new(BufReader::new(reads_fa));
    while let Some(rec) = reader.next_record()? {
        if rec.sequence.is_empty() {
            continue;
        }
        let cov = kmer_coverage_vector(rec.sequence.as_bytes(), counts, k, ds);
        // cpp median_coverage 按值传参——排序发生在副本上; mean/stDev 用**原序**向量
        // （f32 累加对顺序敏感，先排序会改变 stdev 舍入，实测 3 ulp 差异）
        let mut sorted = cov.clone();
        rows.push(CoverageStatsRow {
            acc: rec.accession,
            median: median_cov(&mut sorted),
            mean: mean_f32(&cov),
            stdev: stdev_f32(&cov),
        });
    }
    Ok(rows)
}

/// cpp L116-120/L142-148 输出。tid 列原样 `thread:{tid}`; 原版 tid 是 omp 线程号，
/// 本移植单线程恒 0。mean/stdev 用 format_g6（ostream 默认 6 位有效数字）。
pub fn write_stats_tsv<W: Write>(w: &mut W, rows: &[CoverageStatsRow]) -> std::io::Result<()> {
    writeln!(w, "acc\tmedian_cov\tmean_cov\tstdev\ttid")?;
    for r in rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\tthread:0",
            r.acc,
            r.median,
            crate::fmt::format_g6(r.mean),
            crate::fmt::format_g6(r.stdev)
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn key(s: &str) -> KmerId {
        kmer_to_intval(s.as_bytes()).unwrap()
    }

    #[test]
    fn median_hand_vectors() {
        assert_eq!(median_cov(&mut [2u32, 4]), 3);
        assert_eq!(median_cov(&mut [1, 2]), 1); // 偶数个截断: (1+2)/2=1
        assert_eq!(median_cov(&mut [5]), 5);
        assert_eq!(median_cov(&mut []), 0);
        assert_eq!(median_cov(&mut [1, 2, 3]), 2);
        assert_eq!(median_cov(&mut [10, 4, 7]), 7); // 先排序 → [4,7,10]
                                                    // cpp 无符号截断: (1+2)/2 = 1（非 1.5）
        assert_eq!(median_cov(&mut [2, 5]), 3);
        assert_eq!(median_cov(&mut [1, 4]), 2);
    }

    #[test]
    fn mean_hand_vectors() {
        assert_eq!(mean_f32(&[1, 2]), 1.5);
        assert_eq!(mean_f32(&[]), 0.0);
        assert_eq!(mean_f32(&[3]), 3.0);
        // 大计数走 i64 求和（原版 long）再转 f32: (2^32+1)/2 精确
        assert_eq!(mean_f32(&[u32::MAX, 1]), 2147483648.0f32);
    }

    #[test]
    fn stdev_hand_vectors() {
        // n=1: 0/0 → NaN（x86 默认 QNaN 带符号位 → 原版输出 "-nan"）
        let nan = stdev_f32(&[5]);
        assert!(nan.is_nan());
        assert_eq!(crate::fmt::format_g6(f32::NAN.copysign(-1.0)), "-nan");
        // n=0: 0.0/-1 = -0 → sqrt(-0) = -0
        let neg0 = stdev_f32(&[]);
        assert_eq!(neg0, -0.0);
        assert_eq!(crate::fmt::format_g6(neg0), "-0");
        // [2,4]: avg=3, Σd²=2, 2/(2-1)=2 → sqrt(2)
        assert_eq!(stdev_f32(&[2, 4]), 2.0f32.sqrt());
        assert_eq!(crate::fmt::format_g6(stdev_f32(&[2, 4])), "1.41421");
        // 单调重复值 → stdev = 0（非 NaN）
        assert_eq!(crate::fmt::format_g6(stdev_f32(&[1, 1, 1, 1])), "0");
    }

    #[test]
    fn coverage_vector_semantics() {
        let mut m: CountMap = Default::default();
        m.insert(key("ACGT"), 7);
        // 空 counts + "ACGTACGT"(k=4) → 5 窗口全部按 1
        let empty: CountMap = Default::default();
        assert_eq!(
            kmer_coverage_vector(b"ACGTACGT", &empty, 4, false),
            vec![1, 1, 1, 1, 1]
        );
        // 查表命中 7 / 缺失按 1: ACGT,CGTA,GTAC,TACG,ACGT
        assert_eq!(
            kmer_coverage_vector(b"ACGTACGT", &m, 4, false),
            vec![7, 1, 1, 1, 7]
        );
        // 窗口含 N → 计数 1（不跳过! 与计数器语义不同）: ACGTNACG → 5 窗口
        assert_eq!(
            kmer_coverage_vector(b"ACGTNACG", &m, 4, false),
            vec![7, 1, 1, 1, 1]
        );
        // ACGT 是回文（canonical = 自身）→ DS 查同一键
        assert_eq!(
            kmer_coverage_vector(b"ACGTACGT", &m, 4, true),
            vec![7, 1, 1, 1, 7]
        );
        // DS 差分: 表键 canonical("AAAA")=TTTT（编码序 max）; AAAA 窗口 DS 命中、SS 不中
        let mut ds_m: CountMap = Default::default();
        ds_m.insert(key("TTTT"), 5);
        assert_eq!(
            kmer_coverage_vector(b"TTTTAAAA", &ds_m, 4, true),
            vec![5, 1, 1, 1, 5]
        );
        assert_eq!(
            kmer_coverage_vector(b"TTTTAAAA", &ds_m, 4, false),
            vec![5, 1, 1, 1, 1]
        );
        // 短序列/空序列 → 空向量
        assert_eq!(
            kmer_coverage_vector(b"ACG", &m, 4, false),
            Vec::<u32>::new()
        );
        assert_eq!(kmer_coverage_vector(b"", &m, 4, false), Vec::<u32>::new());
        // 恰 k 长度 → 1 窗口
        assert_eq!(kmer_coverage_vector(b"ACGT", &m, 4, false), vec![7]);
    }

    #[test]
    fn dump_loader_mirrors_populate() {
        let dump = b">3\nACGTACGTACGTACGTACGTACGTA\n>5\nTTTTGGGGCCCCAAAATTTTGGGGC\n";
        let m = load_kmer_dump(dump, 25, false).unwrap();
        assert_eq!(m[&key("ACGTACGTACGTACGTACGTACGTA")], 3);
        assert_eq!(m.len(), 2);
        // 重复键累加（cpp add_kmer: _kmer_counter[kmer_val] += count）
        let dup = b">2\nACGTACGTACGTACGTACGTACGTA\n>3\nACGTACGTACGTACGTACGTACGTA\n";
        let m = load_kmer_dump(dup, 25, false).unwrap();
        assert_eq!(m[&key("ACGTACGTACGTACGTACGTACGTA")], 5);
        // 长度不符 → 跳过
        let bad = b">9\nACGT\n";
        assert!(load_kmer_dump(bad, 25, false).unwrap().is_empty());
        // atoi 语义: 前导数字后缀垃圾忽略; 无数字 → 0（仍入库，查询侧 clamp 到 1）
        let junk = b">12xyz\nACGTACGTACGTACGTACGTACGTA\n>abc\nTTTTGGGGCCCCAAAATTTTGGGGC\n";
        let m = load_kmer_dump(junk, 25, false).unwrap();
        assert_eq!(m[&key("ACGTACGTACGTACGTACGTACGTA")], 12);
        assert_eq!(m[&key("TTTTGGGGCCCCAAAATTTTGGGGC")], 0);
        // DS: dump 的词典序代表键 canonical 化入库（AAAA 条目落 TTTT 键）
        let dsm = load_kmer_dump(b">4\nAAAAAAAAAAAAAAAAAAAAAAAAA\n", 25, true).unwrap();
        assert_eq!(dsm.len(), 1);
        assert_eq!(dsm[&key("TTTTTTTTTTTTTTTTTTTTTTTTT")], 4);
        // 负计数（atoi 负号）按位回绕成大 u32（镜像 C int→unsigned 赋值）
        let neg = load_kmer_dump(b">-1\nAAAAAAAAAAAAAAAAAAAAAAAAA\n", 25, false).unwrap();
        assert_eq!(neg[&key("AAAAAAAAAAAAAAAAAAAAAAAAA")], u32::MAX);
        // 非 gatc kmer → Err（原版 kmer_to_intval throw → exit 1）
        assert!(load_kmer_dump(b">1\nACGTNCGTACGTACGTACGTACGTA\n", 25, false).is_err());
    }

    #[test]
    fn rows_and_tsv_format() {
        let dump = b">3\nACGTACGTACGTACGTACGTACGTA\n";
        let m = load_kmer_dump(dump, 25, false).unwrap();
        // r1 = dump kmer + 前后各延 2 碱基（29bp → 5 窗口 [1,1,3,1,1]，仅中间窗口命中）
        let reads =
            b">r1 desc\nCGACGTACGTACGTACGTACGTACGTACG\n>short\nACGT\n>r3\nTTTTGGGGCCCCAAAATTTTGGGGCTTT\n";
        let rows = coverage_stats_rows(reads, &m, 25, false).unwrap();
        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0].acc, "r1");
        // cov 原序 [1,1,3,1,1] → median 1（排序副本）; mean 7/5=1.4;
        // stdev 按**原序**累加（numpy f32 精确模型）= 0x3f64f92f → "0.894427"
        assert_eq!(rows[0].median, 1);
        assert_eq!(rows[0].mean, 1.4);
        assert_eq!(rows[0].stdev.to_bits(), 0x3f64f92f);
        assert_eq!(rows[1].acc, "short"); // <25 → 空向量
        assert_eq!(rows[1].median, 0);
        assert_eq!(rows[1].mean, 0.0);
        assert_eq!(rows[1].stdev, -0.0);
        // r3: 28bp → 4 窗口全缺失 → [1,1,1,1]
        assert_eq!(rows[2].median, 1);
        assert_eq!(rows[2].mean, 1.0);
        assert_eq!(rows[2].stdev, 0.0);
        let mut buf = Vec::new();
        write_stats_tsv(&mut buf, &rows).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "acc\tmedian_cov\tmean_cov\tstdev\ttid");
        assert_eq!(lines[1], "r1\t1\t1.4\t0.894427\tthread:0");
        assert_eq!(lines[2], "short\t0\t0\t-0\tthread:0");
        assert_eq!(lines[3], "r3\t1\t1\t0\tthread:0");
    }

    #[test]
    fn empty_sequence_records_skip_row() {
        let m: CountMap = Default::default();
        // header 后直接跟下一 header → 空 seq → 原版 continue（无行）
        let reads = b">a\n>b\nTTTTGGGGCCCCAAAATTTTGGGGC\n";
        let rows = coverage_stats_rows(reads, &m, 25, false).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].acc, "b");
    }

    /// 畸形 FASTA（首行非 '>'）→ Err 向上传播（原版 throw → exit(1)）。
    /// 修 T3 审查 Minor 1: 签名从 unwrap/panic 改为 Result。
    #[test]
    fn malformed_reads_fa_is_err() {
        let m: CountMap = Default::default();
        assert!(coverage_stats_rows(b"notafasta\nACGT\n", &m, 25, false).is_err());
    }
}
