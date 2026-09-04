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
/// `--multi-k` 合法下界（KmerId 位宽上限见上界 MAX_KMER_LENGTH=128）。
pub const MULTI_K_MIN: usize = 15;
/// `--multi-k` 合法上界: 127 = MAX_KMER_LENGTH(128) - 1（KmerId=U256 / 2bit,
/// 2026-09-02 扩容后由 63 放宽, 对齐 SPAdes MAX_K=127）。保守留 1 碱基余量:
/// k=128 物理可编码（编码层 MAX_KMER_LENGTH 允许）, 但 CLI 支持域取 127——
/// 让 `(1 << 2k)` 族表达式恒落在位宽内（2k ≤ 254）, 不依赖满宽边界的自限回绕。
/// rnaSPAdes 公式域内 RL=150 的 [49,73] 与 RL=250 的 [55,123] 均原值落地。
pub const MULTI_K_MAX: usize = MAX_KMER_LENGTH - 1;

/// auto 模式的顶 k 钳制（默认不回退护栏, 2026-09-02）: RL=150 的公式值
/// upper_k=73 在 QM4 实测**劣于** 63——[29,49,73] R=44.4% vs [29,49,63]
/// R=45.3%（−0.9pp, 超波动带）, [49,73] R=43.3% vs [49,63] R=43.9%（−0.6pp）。
/// 机制: kmer 深度 ∝ (RL−k+1), 大 k 对低表达尾巴的断图损失 > 解纠缠收益
/// （rnaSPAdes 公式的迭代简化假设与我们 trusted 注入 + 末轮完整管线的架构
/// 不共享）。长读（RL≥200）公式值 [55,123] 未实测, 默认同样钳到 63;
/// 显式 `--multi-k` 列表不受此限（支持域 [15,127]）。
pub const AUTO_K_CAP: usize = 63;

/// auto 公式的实测钳制: 各公式值 min AUTO_K_CAP 后走支持域钳制。
pub fn clamp_auto_k_list(ks: Vec<usize>) -> Vec<usize> {
    clamp_k_list(ks.into_iter().map(|k| k.min(AUTO_K_CAP)).collect())
}

/// rnaSPAdes auto k 公式（spades_stage.py `rna_k_values` 直译）:
/// - `upper_k = int(RL/2) - 1`
/// - `lower_k = clamp(int(RL/3), RNA_MIN_K=29, RNA_MAX_LOWER_K=55)`
/// - 两者均向下取奇（偶则 -1）
/// - `upper_k <= lower_k` 时单 k = `[upper_k]`, 否则 `[lower_k, upper_k]`
///
/// 忠实保留公式原值——可能超出本实现支持域 [15, 127]（如 RL=300 → [55,149]）,
/// 支持域钳制归 [`clamp_k_list`]（对应 SPAdes 公式内的 MAX_K=127 钳制步骤,
/// 我们的位宽上限经 U256 扩容后为 127）。u64 时代的 [29,37]→[29,32] 钳制已解除。
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
/// 例: RL=76 公式 [29,37] 原值落地; RL=150 公式 [49,73] **U256 扩容后原值
/// 落地**（u128 时代被钳到 [49,63]）; RL=250 [55,123] 同样原值落地。
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

/// `--weld-min-glue` auto 的读长分界（2026-09-04 转正）: RL ≥ 150 → 1（mg1,
/// QM4 组合验证 R+0.6）; 短读 → 2（原版, dev1 76bp 上单 read 证据 P−1.7 系噪声）。
pub const WELD_MG1_READ_LEN_MIN: usize = 150;

/// `--weld-min-glue` auto（哨兵 0）resolve: 按抽样读长给 mg1/原版阈值;
/// 显式值（≠0）原样返回、不读文件。`detected_rl` 传入已测读长避免重复抽样
/// （None 且需要时才现场测）。
pub fn resolve_weld_min_glue(
    args: &TrinityArgs,
    detected_rl: Option<usize>,
) -> Result<u32, CommonError> {
    if args.weld_min_glue != 0 {
        return Ok(args.weld_min_glue);
    }
    let rl = match detected_rl {
        Some(rl) => rl,
        None => detect_read_length(args)?,
    };
    Ok(if rl >= WELD_MG1_READ_LEN_MIN { 1 } else { 2 })
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
        // U256 扩容（[15,127]）: RL=76/150/250 的公式值全部原值落地, 不再被钳
        assert_eq!(clamp_k_list(vec![29, 37]), vec![29, 37]);
        assert_eq!(clamp_k_list(vec![49, 73]), vec![49, 73]);
        assert_eq!(clamp_k_list(vec![55, 123]), vec![55, 123]);
        // 双侧触顶去重成单 k（构造超出 [15,127]: [130,140] → [127]）
        assert_eq!(clamp_k_list(vec![130, 140]), vec![127]);
        // 下界: 公式给出 <15 的（短读长）→ 钳到 15
        assert_eq!(clamp_k_list(vec![11]), vec![15]);
        // 域内不变（含旧上界 63 以上区域与新上界 127）
        assert_eq!(clamp_k_list(vec![21, 25]), vec![21, 25]);
        assert_eq!(clamp_k_list(vec![29, 55, 63]), vec![29, 55, 63]);
        assert_eq!(clamp_k_list(vec![29, 49, 73]), vec![29, 49, 73]);
    }

    /// auto 实测钳制: 公式值超 AUTO_K_CAP 的部分压到 63（RL=150 → [49,63];
    /// RL=250 [55,123] → [55,63]）; ≤63 的公式值（RL=76 [29,37]）原值落地;
    /// 双双压到 63 时去重成单 k。
    #[test]
    fn clamp_auto_k_list_caps_top() {
        assert_eq!(clamp_auto_k_list(vec![49, 73]), vec![49, 63]);
        assert_eq!(clamp_auto_k_list(vec![29, 37]), vec![29, 37]);
        assert_eq!(clamp_auto_k_list(vec![55, 123]), vec![55, 63]);
        assert_eq!(clamp_auto_k_list(vec![73, 123]), vec![63]);
    }

    /// weld-min-glue auto 分界: RL≥150 → 1（mg1）, 短读 → 2（原版）;
    /// 显式值直通（不触发读长抽样）。
    #[test]
    fn resolve_weld_min_glue_boundaries() {
        fn args_with(mg: u32) -> TrinityArgs {
            parse_with_weld_min_glue(mg)
        }
        // 显式值直通（哨兵之外的任何值, 含 1/2）
        assert_eq!(resolve_weld_min_glue(&args_with(1), None).unwrap(), 1);
        assert_eq!(resolve_weld_min_glue(&args_with(2), None).unwrap(), 2);
        assert_eq!(resolve_weld_min_glue(&args_with(5), Some(150)).unwrap(), 5);
        // auto: 分界 150（≥ → mg1, < → 原版）——detected_rl 直通免抽样
        assert_eq!(resolve_weld_min_glue(&args_with(0), Some(150)).unwrap(), 1);
        assert_eq!(resolve_weld_min_glue(&args_with(0), Some(250)).unwrap(), 1);
        assert_eq!(resolve_weld_min_glue(&args_with(0), Some(149)).unwrap(), 2);
        assert_eq!(resolve_weld_min_glue(&args_with(0), Some(76)).unwrap(), 2);
    }

    /// 测试用 TrinityArgs 构造（--weld-min-glue 显式注入, 其余走 CLI 解析）。
    fn parse_with_weld_min_glue(mg: u32) -> TrinityArgs {
        let argv: Vec<String> = [
            "--seqType",
            "fq",
            "--left",
            "l.fq",
            "--right",
            "r.fq",
            "--max_memory",
            "2G",
            "--output",
            "out",
            "--weld-min-glue",
            &mg.to_string(),
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        crate::args::parse_args(&argv).unwrap()
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
