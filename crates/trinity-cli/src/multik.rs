//! `--multi-k`（rnaSPAdes 式多 k 迭代）的策略层: auto k 公式、支持域钳制、
//! reads 读长抽样测定。
//!
//! 机制背景（详见 docs/opt-log.md「rnaSPAdes 构图机制调研」）: 实测失败模式是
//! de Bruijn 图"低深度山谷"（转录本中段 kmer 深度 1-3）导致断图; rnaSPAdes 的
//! 解法是多 k 迭代 + trusted contig 注入——每个 k 独立全量构图, 上一 k 的简化
//! contig 并入下一 k 的 kmer 计数流（贡献序列、不贡献覆盖）。本模块只管
//! 「k 列表如何确定」; 迭代骨架与注入点在 orchestrate.rs。

use std::path::Path;

use trinity_common::error::CommonError;
use trinity_common::fasta::FastaReader;
use trinity_common::fastq::FastqReader;
use trinity_common::kmer::MAX_KMER_LENGTH;

use crate::args::{SeqType, TrinityArgs};

/// rnaSPAdes RNA_MIN_K（spades options_storage.py:77）。
pub const RNA_MIN_K: usize = 29;
/// rnaSPAdes RNA_MAX_LOWER_K（options_storage.py:78）。
pub const RNA_MAX_LOWER_K: usize = 55;
/// `--multi-k` 合法下界（KmerId 位宽上限见上界 MAX_KMER_LENGTH=32）。
pub const MULTI_K_MIN: usize = 15;
/// `--multi-k` 合法上界（u64 / 2bit 编码）。
pub const MULTI_K_MAX: usize = MAX_KMER_LENGTH;

/// rnaSPAdes auto k 公式（spades_stage.py `rna_k_values` 直译）:
/// - `upper_k = int(RL/2) - 1`
/// - `lower_k = clamp(int(RL/3), RNA_MIN_K=29, RNA_MAX_LOWER_K=55)`
/// - 两者均向下取奇（偶则 -1）
/// - `upper_k <= lower_k` 时单 k = `[upper_k]`, 否则 `[lower_k, upper_k]`
///
/// 忠实保留公式原值——可能超出本实现支持域 [15, 32]（如 RL=76 → [29, 37]）,
/// 支持域钳制归 [`clamp_k_list`]（对应 SPAdes 公式内的 MAX_K=127 钳制步骤,
/// 我们的位宽上限是 32）。
pub fn rna_k_values(rl: usize) -> Vec<usize> {
    // int(RL/2) - 1: usize 除法即截断（正数 = floor）; saturating 防极小 RL 下溢
    let mut upper_k = (rl / 2).saturating_sub(1);
    if upper_k.is_multiple_of(2) {
        upper_k = upper_k.saturating_sub(1);
    }
    // min(max(int(RL/3), 29), 55); 29/55 均为奇, clamp 后再 -1 不出界
    let mut lower_k = (rl / 3).clamp(RNA_MIN_K, RNA_MAX_LOWER_K);
    if lower_k.is_multiple_of(2) {
        lower_k -= 1;
    }
    if upper_k <= lower_k {
        vec![upper_k]
    } else {
        vec![lower_k, upper_k]
    }
}

/// 公式值 → 支持域 [MULTI_K_MIN, MULTI_K_MAX] 钳制 + 保序去重。
/// 输入递升 → 钳制后仍递升（可能出相邻重复）, `dedup` 收敛; 结果恒非空。
/// 例: RL=76 公式 [29,37] → [29,32]; RL=150 公式 [49,73] → [32]（双双触顶
/// 去重成单 k, 退化为 k=32 的 classic 结构——位宽上限下的诚实结果）。
pub fn clamp_k_list(ks: Vec<usize>) -> Vec<usize> {
    let mut out: Vec<usize> = ks
        .into_iter()
        .map(|k| k.clamp(MULTI_K_MIN, MULTI_K_MAX))
        .collect();
    out.dedup();
    out
}

/// `--multi-k auto` 的读长测定: 从各 reads 文件（原始输入, gz 透明）抽样
/// 前 1000 条记录, 每文件取样本内最大读长, 跨文件取最小——SPAdes
/// `get_reads_length`（file_operations.py:319-331: per-file max → min over
/// files）的语义对齐。所有文件均无记录时报错。
pub fn detect_read_length(args: &TrinityArgs) -> Result<usize, CommonError> {
    /// 每文件抽样上限（任务口径; SPAdes 用 10^4, 1000 对定长测序足够收敛）。
    const SAMPLE: usize = 1000;
    let files: Vec<&Path> = args
        .left
        .iter()
        .chain(&args.right)
        .chain(&args.single)
        .map(|p| p.as_path())
        .collect();
    let mut rl_min: Option<usize> = None;
    for f in files {
        let reader = trinity_common::io_util::open_maybe_gz(f)?;
        let mut max_len = 0usize;
        match args.seq_type {
            SeqType::Fq => {
                let mut r = FastqReader::new(reader);
                for _ in 0..SAMPLE {
                    match r.next_record()? {
                        Some(rec) => max_len = max_len.max(rec.sequence.len()),
                        None => break,
                    }
                }
            }
            SeqType::Fa => {
                let mut r = FastaReader::new(reader);
                for _ in 0..SAMPLE {
                    match r.next_record()? {
                        Some(rec) => max_len = max_len.max(rec.sequence.len()),
                        None => break,
                    }
                }
            }
        }
        if max_len > 0 {
            rl_min = Some(rl_min.map_or(max_len, |m| m.min(max_len)));
        }
    }
    rl_min.ok_or_else(|| {
        CommonError::Parse(
            "--multi-k auto: 无法从 reads 文件测定读长（所有输入均无记录）".to_string(),
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 任务验收边界: RL=76 → [29, 37]; RL=50 → 单 k。
    #[test]
    fn rna_k_values_acceptance_boundaries() {
        // upper = 76/2-1 = 37（奇）; lower = clamp(25, 29, 55) = 29 → [29, 37]
        assert_eq!(rna_k_values(76), vec![29, 37]);
        // upper = 50/2-1 = 24 → 偶取奇 23; lower = clamp(16,29,55) = 29;
        // 23 <= 29 → 单 k
        assert_eq!(rna_k_values(50), vec![23]);
    }

    #[test]
    fn rna_k_values_spades_mirror() {
        // RL=150: upper = 74 → 73; lower = clamp(50,29,55) = 50 → 49 → [49, 73]
        assert_eq!(rna_k_values(150), vec![49, 73]);
        // RL=250: upper = 124 → 123; lower = clamp(83,29,55) = 55 → [55, 123]
        assert_eq!(rna_k_values(250), vec![55, 123]);
        // RL=76 奇偶降档: lower 公式值 25 被 clamp 到 29（奇, 不再降）
        // RL=77: upper = 38-1 = 37; lower = clamp(25,29,55) = 29 → [29, 37]
        assert_eq!(rna_k_values(77), vec![29, 37]);
        // RL=60: upper = 29; lower = 29 → 29 <= 29 → 单 k [29]
        assert_eq!(rna_k_values(60), vec![29]);
        // 极小读长不下溢（防御; clamp_k_list 会再兜到 15）
        assert_eq!(rna_k_values(2), vec![0]);
    }

    #[test]
    fn clamp_k_list_to_supported_domain() {
        // RL=76 的公式值 [29,37] → 37 超位宽上限 → 钳到 32
        assert_eq!(clamp_k_list(vec![29, 37]), vec![29, 32]);
        // RL=150: [49,73] 双双触顶 → 去重成单 k=32
        assert_eq!(clamp_k_list(vec![49, 73]), vec![32]);
        // 下界: 公式给出 <15 的（短读长）→ 钳到 15
        assert_eq!(clamp_k_list(vec![11]), vec![15]);
        // 域内不变
        assert_eq!(clamp_k_list(vec![21, 25]), vec![21, 25]);
    }

    fn write(p: &Path, bytes: &[u8]) {
        std::fs::write(p, bytes).unwrap();
    }

    fn args_for(files: &[&Path], seq_type: SeqType) -> TrinityArgs {
        let argv: Vec<String> = [
            "--seqType",
            match seq_type {
                SeqType::Fq => "fq",
                SeqType::Fa => "fa",
            },
            "--single",
            files[0].to_str().unwrap(),
            "--max_memory",
            "1G",
            "--output",
            "out",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        crate::args::parse_args(&argv).unwrap()
    }

    #[test]
    fn detect_read_length_fq() {
        let dir = std::env::temp_dir().join(format!("multik_rl_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fq = dir.join("a.fq");
        write(&fq, b"@r1\nACGTACGTAC\n+\nIIIIIIIIII\n@r2\nACGT\n+\nIIII\n");
        let a = args_for(&[&fq], SeqType::Fq);
        assert_eq!(detect_read_length(&a).unwrap(), 10); // 抽样内最大读长
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn detect_read_length_fa_and_gz() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write as _;
        let dir = std::env::temp_dir().join(format!("multik_rlgz_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("a.fa.gz");
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(b">r1\nACGTAC\n>r2\nACGTACGTACGTACGT\n")
            .unwrap();
        std::fs::write(&fa, enc.finish().unwrap()).unwrap();
        let a = args_for(&[&fa], SeqType::Fa);
        assert_eq!(detect_read_length(&a).unwrap(), 16); // r2 序列 16bp（gzip 透明读取）
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn detect_read_length_empty_is_error() {
        let dir = std::env::temp_dir().join(format!("multik_rle_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fq = dir.join("empty.fq");
        write(&fq, b"");
        let a = args_for(&[&fq], SeqType::Fq);
        assert!(detect_read_length(&a).is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
