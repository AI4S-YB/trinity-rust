//! DigiNorm 编排 — insilico_read_normalization.pl 的 together 模式镜像:
//! fq→fa（名字规范 + SS 'R' 端 revcomp）→ 计数（dump -L 2 语义）→ 双端 stats → merge → nbkc 选择
//! → 按名单扫描原始文件抽取记录（原格式原样写出）。
//!
//! 原版管线（v2.15.2，行号对照）: prep L318-368/L660-719（fq 走 `seqtk-trinity seq -A -R <1|2> [-r]`，
//! fa 走 `cat`——**名字已规范假定，不过 seqtk**; SS 端类型恰为 "R" 时 fa 走 revcomp_fasta.pl），
//! both.fa = `cat left.fa right.fa` L350; jellyfish count [--canonical] + `dump -L 2` L588-654;
//! fastaToKmerCoverageStats per-side L816-898 + `sort -k1,1`; PE merge L950-985
//! （nbkc_merge_left_right_stats.pl --sorted）; nbkc_normalize.pl（srand(12345)）选名单;
//! 抽取 L505-584（Fastq_reader/Fasta_reader.pm 逐记录，命中原样写出）;
//! 输出命名 L398-410 + left/right/single.norm.{fq,fa} 便捷名 L442-453（原版为 symlink，此处双写）。
//!
//! SS_lib_type 语义（读原版确认，与任务稿的假设不同处以源码为准）:
//! - 空 = DS（jellyfish --canonical + stats --DS）;
//! - 单端 'R' 在 L200-202 被改写为 'F'（"just treat it as F ... waste of time"）→ **单端恒不 revcomp**;
//! - PE 按 `split(//, $SS_lib_type)` 逐字符: 'RF' → left='R' → **left revcomp**;
//!   'FR' → right='R' → **right revcomp**（不是"无 revcomp"!）。
//!
//! **SS 模式 stats 阶段恒 canonical**（两处源码互证）:
//! - fastaToKmerCoverageStats.cpp:75 `bool is_DS = (! args.isArgSet("--SS"))` → DS 是默认;
//! - insilico_read_normalization.pl:838 只在**非** SS 时追加 ` --DS `、从不传 `--SS`
//!   → stats 工具的 is_DS 恒真，与 SS_lib_type 无关。
//!
//! 故 SS 链路是: jellyfish 无 --canonical 计数（literal 键）→ dump -L 过滤（literal 计数）→
//! fastaToKmerCoverageStats 以 DS 装载 dump（KmerCounter.cpp:476-487 add_kmer: 键 canonical 化 +
//! `_kmer_counter[kmer_val] += count` 求和）→ canonical 键查询。SS_lib_type 只影响计数阶段，
//! 触发场景: SS F PE 左右端互补（典型 dUTP 文库）与单端混链输入。
//!
//! 双次文本往返是本模块的关键保真点（T4 审查 Important）:
//! 1. stats → StatsRow: mean/stdev 经 `format_g6 → parse::<f64>()`（镜像 stats TSV 落盘再被
//!    nbkc/merge 脚本数值化解析; 直接 f32 as f64 会在 CV/ratio 边界分歧）。median 是整数文本 → 精确直转。
//! 2. PE merge → StatsRow: MergedRow 的 %.1f 字符串再 parse 回 f64（镜像 pairs.stats 落盘文本）。

use std::collections::HashSet;
use std::io::Read;
use std::path::{Path, PathBuf};

use trinity_common::error::CommonError;
use trinity_common::io_util::open_maybe_gz;
use trinity_common::kmer::get_ds_kmer_val;

use crate::counter::{CountMap, KmerCountTable};
use crate::coverage_stats::{coverage_stats_rows, write_stats_tsv, CoverageStatsRow};
use crate::drand48::Drand48;
use crate::nbkc::{merge_pairs, select, MergedRow, NbkcParams, StatsRow};
use crate::read_names::{core_read_name, fq_records_to_fa, ReadType};

/// together 模式参数（默认值 = 原版默认: KMER_SIZE 25、MIN_KMER_COV_CONST 2 "DO NOT CHANGE"、
/// max_cov 由 CLI 必填、min_cov 0→本 CLI 取 1、max_CV 10000）。
#[derive(Debug, Clone)]
pub struct DigiNormParams {
    pub k: usize,
    /// jellyfish `dump -L`（原版 L45 固定 2）: 计数 < 2 的 kmer 不进查表（查表侧缺失按 1）。
    pub min_kmer_cov_const: u32,
    pub max_cov: f64,
    pub min_cov: f64,
    pub max_cv: f64,
    /// SS 'R' 端 revcomp（PE: left 对应 'RF'; 单端原版恒 false，见模块文档）。
    pub ss_revcomp_left: bool,
    /// PE 'FR' → right revcomp。
    pub ss_revcomp_right: bool,
    /// 非 SS 即 true。只影响**计数阶段**（jellyfish --canonical）; stats 阶段恒 DS
    /// （见模块文档），SS 时由 stats_query_table 补 canonical 合并。
    pub ds: bool,
}

impl Default for DigiNormParams {
    fn default() -> Self {
        DigiNormParams {
            k: 25,
            min_kmer_cov_const: 2,
            max_cov: 200.0,
            min_cov: 1.0,
            max_cv: 10000.0,
            ss_revcomp_left: false,
            ss_revcomp_right: false,
            ds: true,
        }
    }
}

/// 输入: 双端（left/right）或单端。格式（fq/fa）按扩展名或魔数嗅探，各文件
/// 与两侧须一致。多文件（逗号列表展开）：按列表序拼接成一个输入流，全部
/// 文件**一起**归一化（单一输出、单一 outdir——原版 prep_list_of_seqs 把
/// 每侧文件列表串成一个 left.fa/right.fa/single.fa，kmer 计数与选取都在
/// 合并流上做；输出 basename = 第一个文件名，多文件时加 `_ext_all_reads`，
/// insilico_read_normalization.pl L399）。
pub enum ReadsInput {
    Single(Vec<PathBuf>),
    Paired(Vec<PathBuf>, Vec<PathBuf>),
}

/// 一侧文件列表：按序读入拼接 + 逐文件格式判定（须全一致）。
fn read_side(paths: &[PathBuf]) -> Result<(Vec<u8>, SeqFormat), CommonError> {
    let mut data = Vec::new();
    let mut fmt_all: Option<SeqFormat> = None;
    for p in paths {
        let raw = read_input(p)?;
        let f = detect_format(p, &raw)?;
        match fmt_all {
            Some(f0) if f0 != f => {
                return Err(CommonError::FastqFormat(format!(
                    "Error, {} is {f0:?} but {} is {f:?}; seqType must agree",
                    paths[0].display(),
                    p.display()
                )))
            }
            _ => fmt_all = Some(f),
        }
        data.extend_from_slice(&raw);
    }
    Ok((data, fmt_all.unwrap()))
}

/// 输出 basename：第一个文件名；多文件时 `<first>_ext_all_reads`（L399）。
fn side_base_name(paths: &[PathBuf]) -> String {
    let first = paths[0]
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "reads".to_string());
    if paths.len() == 1 {
        first
    } else {
        format!("{first}_ext_all_reads")
    }
}

/// 输出文件（长名，镜像原版 L398-410 的 `{basename}.normalized_K.._maxC.._minC.._maxCV..` 命名）。
/// 便捷名 `left.norm.*`/`right.norm.*`/`single.norm.*` 一并写出（内容双写，见 [`run`]）。
pub struct DigiNormOutputs {
    pub left: Option<PathBuf>,
    pub right: Option<PathBuf>,
    pub single: Option<PathBuf>,
}

/// 输入序列格式。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqFormat {
    Fq,
    Fa,
}

/// 运行 together 模式 DigiNorm。out_dir 递归创建；返回各侧输出路径。
///
/// 单端走 ss_revcomp_left 作为该侧标志（库层通用性; 原版管线因 L200-202 的 'R'→'F'
/// 改写，单端实际恒为 false）。
pub fn run(
    params: &DigiNormParams,
    reads: &ReadsInput,
    out_dir: &Path,
) -> Result<DigiNormOutputs, CommonError> {
    std::fs::create_dir_all(out_dir)?;
    match reads {
        ReadsInput::Paired(left_paths, right_paths) => {
            let (left_raw, lf) = read_side(left_paths)?;
            let (right_raw, rf) = read_side(right_paths)?;
            let fmt = if lf != rf {
                return Err(CommonError::FastqFormat(format!(
                    "Error, left is {lf:?} but right is {rf:?}; seqType must agree"
                )));
            } else {
                lf
            };
            let left_fa = prep_side(&left_raw, fmt, ReadType::R1, params.ss_revcomp_left)?;
            let right_fa = prep_side(&right_raw, fmt, ReadType::R2, params.ss_revcomp_right)?;
            // both.fa = left.fa 字节后接 right.fa 字节（原版 `cat left.fa right.fa` L350）
            let mut both = left_fa.clone();
            both.extend_from_slice(&right_fa);
            let table = stats_query_table(&both, params);
            let left_rows = sorted_stats_rows(&left_fa, &table, params)?;
            let right_rows = sorted_stats_rows(&right_fa, &table, params)?;
            let merged = merge_pairs(&left_rows, &right_rows);
            let final_rows = merged_to_stats_rows(&merged);
            let selected = select_rows(&final_rows, params)?;
            let index = index_keys(&selected);
            let out_left = normalized_out_path(out_dir, &side_base_name(left_paths), fmt, params);
            let out_right = normalized_out_path(out_dir, &side_base_name(right_paths), fmt, params);
            extract_to_file(&left_raw, fmt, &index, &out_left)?;
            extract_to_file(&right_raw, fmt, &index, &out_right)?;
            write_convenience(&out_left, out_dir, "left", fmt)?;
            write_convenience(&out_right, out_dir, "right", fmt)?;
            Ok(DigiNormOutputs {
                left: Some(out_left),
                right: Some(out_right),
                single: None,
            })
        }
        ReadsInput::Single(paths) => {
            let (raw, fmt) = read_side(paths)?;
            let fa = prep_side(&raw, fmt, ReadType::R1, params.ss_revcomp_left)?;
            let table = stats_query_table(&fa, params);
            let rows = sorted_stats_rows(&fa, &table, params)?;
            let selected = select_rows(&rows, params)?;
            let index = index_keys(&selected);
            let out = normalized_out_path(out_dir, &side_base_name(paths), fmt, params);
            extract_to_file(&raw, fmt, &index, &out)?;
            write_convenience(&out, out_dir, "single", fmt)?;
            Ok(DigiNormOutputs {
                left: None,
                right: None,
                single: Some(out),
            })
        }
    }
}

// ---------------------------------------------------------------- 读入与格式

/// 读原始输入（gzip 魔数嗅探透明解压——原版经 FIFO+gunzip 等价）。
fn read_input(path: &Path) -> Result<Vec<u8>, CommonError> {
    let mut data = Vec::new();
    open_maybe_gz(path)?.read_to_end(&mut data)?;
    Ok(data)
}

/// 格式判定: 优先扩展名（.gz 先剥再判），否则按首个非空行首字节嗅探（'@'→fq, '>'→fa）。
fn detect_format(path: &Path, data: &[u8]) -> Result<SeqFormat, CommonError> {
    let mut name = path
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    if let Some(stem) = name.strip_suffix(".gz") {
        name = stem.to_string();
    }
    match name.rsplit('.').next() {
        Some("fa" | "fasta" | "fna") => return Ok(SeqFormat::Fa),
        Some("fq" | "fastq") => return Ok(SeqFormat::Fq),
        _ => {}
    }
    let first = data
        .split(|&b| b == b'\n')
        .find(|l| !l.is_empty())
        .and_then(|l| l.first())
        .copied();
    match first {
        Some(b'@') => Ok(SeqFormat::Fq),
        Some(b'>') => Ok(SeqFormat::Fa),
        other => Err(CommonError::FastqFormat(format!(
            "Error, cannot determine seqType (fq or fa) for {}: first byte {:?} not '@' or '>'",
            path.display(),
            other
        ))),
    }
}

/// PE 两侧格式须一致（原版单一 --seqType）。
#[cfg(test)]
fn detect_pair_format(
    left_path: &Path,
    left: &[u8],
    right_path: &Path,
    right: &[u8],
) -> Result<SeqFormat, CommonError> {
    let lf = detect_format(left_path, left)?;
    let rf = detect_format(right_path, right)?;
    if lf != rf {
        return Err(CommonError::FastqFormat(format!(
            "Error, left is {:?} but right is {:?}; seqType must agree",
            lf, rf
        )));
    }
    Ok(lf)
}

// ---------------------------------------------------------------- prep（L660-719）

/// 一侧预处理: fq → `seqtk-trinity seq -A -R <1|2> [-r]`（fq_records_to_fa）;
/// fa → `cat` 语义（字节原样，**名字已规范假定**——原版 fa 不再过 seqtk），
/// SS 端 'R' 时镜像 revcomp_fasta.pl。
pub fn prep_side(
    data: &[u8],
    fmt: SeqFormat,
    read_type: ReadType,
    revcomp: bool,
) -> Result<Vec<u8>, CommonError> {
    match fmt {
        SeqFormat::Fq => fq_records_to_fa(data, read_type, revcomp),
        SeqFormat::Fa => Ok(if revcomp {
            revcomp_fasta_mirror(data)
        } else {
            data.to_vec()
        }),
    }
}

/// revcomp_fasta.pl 逐字镜像（fa 输入 + SS 端 'R'，原版 L698-700 调用）:
/// 按 `>` 分块（$/=">"，**任意位置**的 '>' 都切）; header = 首个 \n 前;
/// 序列 = 从该 \n 起（**含 \n**）整体字节反转 → tr 映射（\n 无替换 → /d 删除;
/// 换行先参与反转再删，因互补是逐字节的，数学上等价于"拼行后 revcomp"）→ 60 列折行。
/// tr 表: ACGTUMRWSYKVHDBN ↔ TGCAAKYWSRMBDHVN（注意 U→A 不是 T、N→N），
/// 同字母表外字符（E/F/I/X/Z/数字等）原样保留。
const REVCOMP_TR: [(u8, u8); 32] = [
    (b'A', b'T'),
    (b'C', b'G'),
    (b'G', b'C'),
    (b'T', b'A'),
    (b'U', b'A'),
    (b'M', b'K'),
    (b'R', b'Y'),
    (b'W', b'S'),
    (b'S', b'S'),
    (b'Y', b'R'),
    (b'K', b'M'),
    (b'V', b'B'),
    (b'H', b'D'),
    (b'D', b'H'),
    (b'B', b'V'),
    (b'N', b'N'),
    (b'a', b't'),
    (b'c', b'g'),
    (b'g', b'c'),
    (b't', b'a'),
    (b'u', b'a'),
    (b'm', b'k'),
    (b'r', b'y'),
    (b'w', b's'),
    (b's', b's'),
    (b'y', b'r'),
    (b'k', b'm'),
    (b'v', b'b'),
    (b'h', b'd'),
    (b'd', b'h'),
    (b'b', b'v'),
    (b'n', b'n'),
];

fn revcomp_fasta_mirror(data: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(data.len() + data.len() / 60 + 8);
    // while(<>) with $/=">" : 首块为首个 '>' 前的内容（正常为空），其后每块无前导 '>'
    for chunk in data.split(|&b| b == b'>') {
        let Some(nl) = chunk.iter().position(|&b| b == b'\n') else {
            continue; // 无 \n: header 为整块（含空块 → if($header) 跳过空块）
        };
        let header = &chunk[..nl];
        // Perl `if($header)`: 空串与 "0" 均为假
        if header.is_empty() || header == b"0" {
            continue;
        }
        out.extend_from_slice(b">");
        out.extend_from_slice(header);
        out.push(b'\n');
        // reverse(substr($_, lineSepPos)): 含首个 \n 起的所有字节反转
        let mut seq: Vec<u8> = chunk[nl..].iter().rev().copied().collect();
        // tr/.../\n/d: 表内映射，\n 无对应 → 删除; 表外字符保留
        for b in seq.iter_mut() {
            if *b == b'\n' {
                *b = 0; // 标记删除（tr 表与 0 不相交，安全）
            } else if let Some(&(_, to)) = REVCOMP_TR.iter().find(|(from, _)| *from == *b) {
                *b = to;
            }
        }
        seq.retain(|&b| b != 0);
        for line in seq.chunks(60) {
            out.extend_from_slice(line);
            out.push(b'\n');
        }
    }
    out
}

// ---------------------------------------------------------------- 计数（dump -L 语义）

/// stats 查表（对齐原版管道 "jellyfish count → dump -L → fastaToKmerCoverageStats 装载" 三步）:
/// 计数 → -L 过滤 → SS 模式再做 canonical 合并求和。DS 计数的表已 canonical，合并为
/// 恒等映射，直接原样返回（等价于只对 literal 计数表做 [`canonical_merge_sum`]）。
/// 时序推理见 [`canonical_merge_sum`] 的注记。
fn stats_query_table(both_fa: &[u8], params: &DigiNormParams) -> CountMap {
    let filtered = filtered_counts(both_fa, params);
    if params.ds {
        filtered // 已是 canonical 键（jellyfish --canonical），add_kmer 合并为恒等映射
    } else {
        canonical_merge_sum(filtered, params.k)
    }
}

/// `jellyfish count [--canonical] both.fa` + `dump -L min_kmer_cov_const` 语义:
/// 完整计数后过滤 count < min_kmer_cov_const（计数 1 的不进表 → 查表侧缺失按 1，
/// 与"全表 + 查表 clamp"数学等价，此处按原版路径走以保证语义保真）。
/// DS 计数得 canonical 键、SS 得 literal 键——过滤均作用在**该表**的计数上。
fn filtered_counts(both_fa: &[u8], params: &DigiNormParams) -> CountMap {
    KmerCountTable::count_fasta_data(both_fa, params.k, params.ds)
        .into_iter()
        .filter(|&(_, c)| c >= params.min_kmer_cov_const)
        .collect()
}

/// KmerCounter.cpp:476-487 add_kmer 的 DS 装载语义: 键 canonical 化（get_DS_kmer_val）
/// 后 `_kmer_counter[kmer_val] += count` 求和——两条互补链的计数在此合并。
///
/// **先 -L 过滤、后合并**（时序承重）: 原版管道序是 `dump -L` 过滤计数表 → stats 工具
/// add_kmer 求和装载。DS 下 canonical 计数就是合并后的计数（jellyfish 计数时已合并），
/// 两种时序等价; SS 下两条链各 literal 计数 1 的 kmer 先被 dump -L 丢弃、**不参与**合并
/// （canonical 键缺失 → 查表侧 clamp 1）; 若先合并再过滤会错得计数 2。
/// fixtures/p1/diginorm 的 ssB1（左 1 份 + 右 2 份互补，正确 median 2）锁定该序。
fn canonical_merge_sum(table: CountMap, k: usize) -> CountMap {
    let mut merged: CountMap = CountMap::default();
    for (key, count) in table {
        let canon = get_ds_kmer_val(key, k);
        *merged.entry(canon).or_insert(0) += count;
    }
    merged
}

// ---------------------------------------------------------------- stats → 排序 → 文本往返

/// 一侧 stats: 行集 → TSV 文本（与 fastaToKmerCoverageStats 落盘字节一致）→
/// 按 acc 字节序排序（镜像 `head -n1 > sort && tail -n +2 | sort -k1,1`）→
/// 逐行 parse 成 StatsRow（**第一次文本往返**: median 整数文本精确、mean/stdev 经 %g6 文本）。
///
/// 排序用字节序（C locale）; 原版 `sort -k1,1` 用环境 locale，非 ASCII acc 时可能分歧
/// （Trinity 管道 acc 均为 ASCII）。tie-break 比较整行，镜像 GNU sort 的 last-resort 比较。
///
/// 查表恒用 canonical 键（ds=true）: fastaToKmerCoverageStats.cpp:75 `is_DS = !isArgSet("--SS")`
/// 且 insilico_read_normalization.pl:838 只在非 SS 时加 `--DS`、**从不传 `--SS`** →
/// stats 阶段恒 DS。SS 模式的 canonical 化已在 [`stats_query_table`] 的表侧完成。
fn sorted_stats_rows(
    reads_fa: &[u8],
    counts: &CountMap,
    params: &DigiNormParams,
) -> Result<Vec<StatsRow>, CommonError> {
    let rows = coverage_stats_rows(reads_fa, counts, params.k, true)?;
    let mut lines = stats_lines(&rows);
    lines.sort_by(|a, b| {
        acc_field(a)
            .cmp(acc_field(b))
            .then_with(|| a.as_bytes().cmp(b.as_bytes()))
    });
    Ok(lines.iter().map(|l| parse_stats_line(l)).collect())
}

/// 渲染 stats 行文本（无表头; 与 write_stats_tsv 的数据行逐字节一致）。
fn stats_lines(rows: &[CoverageStatsRow]) -> Vec<String> {
    let mut buf = Vec::new();
    write_stats_tsv(&mut buf, rows).expect("Vec<u8> 写入不会失败");
    let text = String::from_utf8(buf).expect("stats 行均为 ASCII/UTF-8");
    text.lines().skip(1).map(str::to_string).collect()
}

fn acc_field(line: &str) -> &str {
    line.split('\t').next().unwrap_or("")
}

/// parse 一行 stats 文本（DelimParser 语义: \t 列; acc/median_cov/mean_cov/stdev）。
/// "nan"/"-nan"/"-0" 等文本 Rust f64 parse 均接受; 数值语义与 perl 数值化一致
/// （NaN 比较恒 false、-0 == 0）——除符号: perl 数值化丢弃 -0 的符号
/// （`"-0"+0 == 0`），Rust parse 保号 → 此处归一为 +0，否则 merge 的 %.1f 会打
/// 出 "-0.0" 而原版是 "0.0"（oracle 实测，pair4 短 read 的 stdev 列）。
fn parse_stats_line(line: &str) -> StatsRow {
    let f: Vec<&str> = line.split('\t').collect();
    StatsRow {
        acc: f.first().copied().unwrap_or("").to_string(),
        median: parse_num(f.get(1).copied().unwrap_or("nan")),
        mean: parse_num(f.get(2).copied().unwrap_or("nan")),
        stdev: parse_num(f.get(3).copied().unwrap_or("nan")),
    }
}

fn parse_num(s: &str) -> f64 {
    let v: f64 = s.parse().unwrap_or(f64::NAN);
    if v == 0.0 {
        0.0
    } else {
        v
    }
}

/// PE merge → StatsRow（**第二次文本往返**: MergedRow 的 %.1f 列被 nbkc 当作
/// pairs.stats 文本 parse; "10.5"/"NaN" 等）。单端跳过 merge。
fn merged_to_stats_rows(merged: &[MergedRow]) -> Vec<StatsRow> {
    merged
        .iter()
        .map(|m| StatsRow {
            acc: m.core.clone(),
            median: parse_num(&m.median),
            mean: parse_num(&m.mean),
            stdev: parse_num(&m.stdev),
        })
        .collect()
}

/// nbkc_normalize.pl: srand(12345) 每 run 一次（PE together 是单次进程调用）。
/// 名单为空输入（stats 无数据行）时镜像原版 die "no reads made it..."。
fn select_rows(rows: &[StatsRow], params: &DigiNormParams) -> Result<Vec<String>, CommonError> {
    if rows.is_empty() {
        return Err(CommonError::FastqFormat(
            "Error, no reads made it to the normalization process...".to_string(),
        ));
    }
    let mut rng = Drand48::new(12345);
    Ok(select(
        rows,
        &NbkcParams {
            max_cov: params.max_cov,
            min_cov: params.min_cov,
            max_cv: params.max_cv,
        },
        &mut rng,
    ))
}

/// build_selected_index（L470-501）: 名单键 = 名单行剥尾部 "/"+单字 \w 字符
/// （`s|/\w$||`——比 core_acc 的仅 /1|/2 更宽: "/3" 也剥; \w = [A-Za-z0-9_]）。
fn index_keys(selected: &[String]) -> HashSet<String> {
    selected
        .iter()
        .map(|s| {
            let b = s.as_bytes();
            if b.len() >= 2 && b[b.len() - 2] == b'/' && is_word_byte(b[b.len() - 1]) {
                s[..s.len() - 2].to_string()
            } else {
                s.clone()
            }
        })
        .collect()
}

fn is_word_byte(b: u8) -> bool {
    b.is_ascii_alphanumeric() || b == b'_'
}

// ---------------------------------------------------------------- 抽取（L505-584）

/// 扫描原始输入抽取命中记录写入 out 文件。返回后校验名单全部命中
/// （原版 die "not all specified records have been retrieved"，文件已写出，镜像）。
fn extract_to_file(
    data: &[u8],
    fmt: SeqFormat,
    index: &HashSet<String>,
    out_path: &Path,
) -> Result<(), CommonError> {
    let mut out = Vec::new();
    let mut found: HashSet<String> = HashSet::new();
    match fmt {
        SeqFormat::Fq => extract_fq_records(data, index, &mut out, &mut found)?,
        SeqFormat::Fa => extract_fa_records(data, index, &mut out, &mut found)?,
    }
    if let Some(parent) = out_path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    std::fs::write(out_path, &out)?;
    let missing: Vec<&String> = index.difference(&found).collect();
    if !missing.is_empty() {
        return Err(CommonError::FastqFormat(format!(
            "Error, not all specified records have been retrieved (missing {}): {}",
            missing.len(),
            missing
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>()
                .join(",")
        )));
    }
    Ok(())
}

/// fq 抽取 — Fastq_reader.pm 的字节级 4 行扫描（**不做 UTF-8 校验、原始字节透传**，
/// P5 字节路径的前置）: 每次取 4 条原始行（空行也是一行、\r 保留、末行可无 \n）;
/// 首 token（'@' 后、kseq isspace 前）经 [`core_read_name`] 查名单; 命中输出原始 4 行。
/// 错误镜像 Perl confess: 少于 4 行 / 首行不以 '@' 开始。
fn extract_fq_records(
    data: &[u8],
    index: &HashSet<String>,
    out: &mut Vec<u8>,
    found: &mut HashSet<String>,
) -> Result<(), CommonError> {
    let mut pos = 0usize;
    while pos < data.len() {
        let mut lines: [Option<&[u8]>; 4] = [None; 4];
        let mut p = pos;
        for slot in lines.iter_mut() {
            if p >= data.len() {
                break;
            }
            let end = data[p..]
                .iter()
                .position(|&b| b == b'\n')
                .map_or(data.len(), |i| p + i + 1);
            *slot = Some(&data[p..end]);
            p = end;
        }
        let n_lines = lines.iter().filter(|l| l.is_some()).count();
        let record: Vec<&[u8]> = lines.iter().flatten().copied().collect();
        if record.is_empty() {
            break;
        }
        if n_lines < 4 {
            return Err(CommonError::FastqFormat(format!(
                "Error, fastQ entry doesn't have 4 lines: {}",
                String::from_utf8_lossy(&record.concat())
            )));
        }
        let name_line = trim_line_end(record[0]);
        if !name_line.starts_with(b"@") {
            return Err(CommonError::FastqFormat(format!(
                "Error, cannot identify first line as read name line: {}",
                String::from_utf8_lossy(record[0])
            )));
        }
        // Fastq_record.pm: read_name = split(/\s+/, $name_line)[0] 去 '@'
        let token = first_token(&name_line[1..]);
        let acc = core_read_name(&String::from_utf8_lossy(token));
        if index.contains(&acc) {
            found.insert(acc);
            for l in record {
                out.extend_from_slice(l);
            }
        }
        pos = p;
    }
    Ok(())
}

/// fa 抽取 — Fasta_reader.pm（$/="\n>"）镜像: 块内先 tr 删控制字符
/// （保留 \t \n; \r 及 >=0x7F 全删）→ header = 首行（含空白，原样保留）→
/// sequence = 其余行 join 后去全部空白 → acc = header 首 token 再剥尾部 /1|/2。
/// 命中输出 `>{header}\n{sequence}\n`（fasta_line_len=-1 → 不折行、序列单行化——
/// 多行 fa 的输出与输入**不是**字节一致，这是原版行为）。无 word 字符的块跳过。
fn extract_fa_records(
    data: &[u8],
    index: &HashSet<String>,
    out: &mut Vec<u8>,
    found: &mut HashSet<String>,
) -> Result<(), CommonError> {
    // 分块: "\n>" 为界（首块 = 首个 "\n>" 前）
    let mut chunks: Vec<&[u8]> = Vec::new();
    let mut start = 0usize;
    for i in 0..data.len().saturating_sub(1) {
        if data[i] == b'\n' && data[i + 1] == b'>' {
            chunks.push(&data[start..i]); // 不含 "\n>"（\n 归上一块——原版分隔符含 \n）
            start = i + 2;
        }
    }
    chunks.push(&data[start..]);
    for chunk in chunks {
        // tr/\t\n\000-\037\177-\377/\t\n/d
        let cleaned: Vec<u8> = chunk
            .iter()
            .copied()
            .filter(|&b| b == b'\t' || b == b'\n' || (0x20..0x7F).contains(&b))
            .collect();
        // next(): 块内无 \w → 跳过（前导空白块 / 尾部空块）
        if !cleaned.iter().any(|&b| is_word_byte(b)) {
            continue;
        }
        let Some(nl) = cleaned.iter().position(|&b| b == b'\n') else {
            return Err(CommonError::FastaFormat(format!(
                "Error, no newline in fasta record: {}",
                String::from_utf8_lossy(&cleaned)
            )));
        };
        let header = &cleaned[..nl];
        // sequence = 其余行 join（\n 已被 split 拿掉）后 s/\s//g（残留空白仅 ' ' 与 \t）
        let sequence: Vec<u8> = cleaned[nl + 1..]
            .iter()
            .copied()
            .filter(|&b| b != b' ' && b != b'\t' && b != b'\n')
            .collect();
        // s/^>|>$//g（分块后首块以 '>' 开始时剥之）
        let header = header.strip_prefix(b">").unwrap_or(header);
        let header = header.strip_suffix(b">").unwrap_or(header);
        // acc = split(/\s+/, $header)[0]，再 s|/[12]\s*$||
        let token = first_token(header);
        let acc = strip_pair_suffix(&String::from_utf8_lossy(token));
        if index.contains(&acc) {
            found.insert(acc);
            out.extend_from_slice(b">");
            out.extend_from_slice(header);
            out.push(b'\n');
            out.extend_from_slice(&sequence);
            out.push(b'\n');
        }
    }
    Ok(())
}

/// L539-541 fa 路径: `s|/[12]\s*$||`（比 core_read_name 少 _forward/_reverse 删除）。
fn strip_pair_suffix(acc: &str) -> String {
    acc.strip_suffix("/1")
        .or_else(|| acc.strip_suffix("/2"))
        .unwrap_or(acc)
        .to_string()
}

/// 行内容去尾部 \n（\r 保留——Perl split(/\n/) 语义）。
fn trim_line_end(line: &[u8]) -> &[u8] {
    line.strip_suffix(b"\n").unwrap_or(line)
}

/// 首 token: 到 kseq/C isspace（' ' \t \v \f \r）为止。
fn first_token(s: &[u8]) -> &[u8] {
    let end = s
        .iter()
        .position(|b| matches!(b, b' ' | b'\t' | 0x0B | 0x0C | b'\r'))
        .unwrap_or(s.len());
    &s[..end]
}

// ---------------------------------------------------------------- 输出命名

/// L398-410: `{basename(原输入)}.normalized_K{K}_maxC{max}_minC{min}_maxCV{cv}.{fq|fa}`。
/// 数值格式化镜像 Perl 整数插值（Rust Display 对整值 f64 恰为无小数点形式）。
fn normalized_out_path(
    out_dir: &Path,
    base: &str,
    fmt: SeqFormat,
    params: &DigiNormParams,
) -> PathBuf {
    let ext = match fmt {
        SeqFormat::Fq => "fq",
        SeqFormat::Fa => "fa",
    };
    out_dir.join(format!(
        "{}.normalized_K{}_maxC{}_minC{}_maxCV{}.{}",
        base, params.k, params.max_cov, params.min_cov, params.max_cv, ext
    ))
}

/// 便捷名双写（原版为 `ln -sf`，此处复制内容，便于下游按固定名取用）。
fn write_convenience(
    long: &Path,
    out_dir: &Path,
    side: &str,
    fmt: SeqFormat,
) -> Result<(), CommonError> {
    let ext = match fmt {
        SeqFormat::Fq => "fq",
        SeqFormat::Fa => "fa",
    };
    let dst = out_dir.join(format!("{side}.norm.{ext}"));
    std::fs::copy(long, dst)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    // ---------------- 格式判定 ----------------

    #[test]
    fn detect_format_by_extension_then_magic() {
        let p = Path::new("/tmp/x.l.fq");
        assert_eq!(detect_format(p, b"@r\nAC\n").unwrap(), SeqFormat::Fq);
        assert_eq!(
            detect_format(Path::new("/tmp/x.fa"), b">r\nAC\n").unwrap(),
            SeqFormat::Fa
        );
        assert_eq!(
            detect_format(Path::new("/tmp/x.fasta"), b"").unwrap(),
            SeqFormat::Fa
        );
        // 无扩展名 → 嗅探（跳过前导空行）
        assert_eq!(
            detect_format(Path::new("/tmp/x"), b"\n\n@r\n").unwrap(),
            SeqFormat::Fq
        );
        assert_eq!(
            detect_format(Path::new("/tmp/x"), b">r\n").unwrap(),
            SeqFormat::Fa
        );
        assert!(detect_format(Path::new("/tmp/x"), b"garbage\n").is_err());
        // .gz 剥后再看扩展名（数据已解压）
        assert_eq!(
            detect_format(Path::new("/tmp/x.fq.gz"), b"@r\n").unwrap(),
            SeqFormat::Fq
        );
    }

    #[test]
    fn pair_format_must_agree() {
        let l = Path::new("/tmp/l.fq");
        let r = Path::new("/tmp/r.fq");
        assert_eq!(
            detect_pair_format(l, b"@a\n", r, b"@b\n").unwrap(),
            SeqFormat::Fq
        );
        assert!(detect_pair_format(l, b"@a\n", Path::new("/tmp/r.fa"), b">b\n").is_err());
    }

    // ---------------- revcomp_fasta 镜像 ----------------

    /// revcomp_fasta.pl 实测（od 逐字节）: ">a desc\nACGTacgtNRY\n>b\nTTTT\n" →
    /// ">a desc\nRYNacgtACGT\n>b\nAAAA\n"（U→A、N→N 的 tr 表; 单行序列输出单行）。
    #[test]
    fn revcomp_fasta_oracle_locked() {
        let out = revcomp_fasta_mirror(b">a desc\nACGTacgtNRY\n>b\nTTTT\n");
        assert_eq!(out, b">a desc\nRYNacgtACGT\n>b\nAAAA\n");
    }

    /// tr 表外字符原样保留; 60 列折行（65bp → 60+5）。
    #[test]
    fn revcomp_fasta_wrap_and_passthrough() {
        let seq = "ACGT".repeat(17); // 68bp → revcomp 后折 60+8
        let inp = format!(">x\n{seq}\n");
        let out = revcomp_fasta_mirror(inp.as_bytes());
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], ">x");
        assert_eq!(lines[1].len(), 60);
        assert_eq!(lines[2].len(), 8);
        assert_eq!(lines[1].parse::<String>().unwrap().chars().count(), 60);
        // 表外字符: 数字与 E/X 原样
        let out = revcomp_fasta_mirror(b">y\n123EX\n");
        assert_eq!(out, b">y\nXE321\n");
    }

    /// \n 参与反转后再删除: 逐字节互补与删除可交换，数学上等价于"先拼行再 revcomp"
    /// （"ACGT\nAAAA" → TTTTACGT; 单行序列是特例）。锁多行输入路径。
    #[test]
    fn revcomp_fasta_multiline() {
        assert_eq!(revcomp_fasta_mirror(b">z\nACGT\nAAAA\n"), b">z\nTTTTACGT\n");
        assert_eq!(revcomp_fasta_mirror(b">z\nACGT\nAAAA\n"), {
            // 与"拼行后 revcomp"的一致性（同上论证）
            let joined = "ACGTAAAA".as_bytes();
            let comp: Vec<u8> = joined
                .iter()
                .rev()
                .map(|&b| match b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    o => o,
                })
                .collect();
            let expect = format!(">z\n{}\n", String::from_utf8(comp).unwrap());
            expect.as_bytes().to_vec()
        });
    }

    // ---------------- dump -L 过滤 ----------------

    #[test]
    fn filtered_counts_drops_below_min_kmer_cov() {
        let mut p = DigiNormParams {
            k: 4,
            ds: false, // SS 模式使键与字面 kmer 一致，便于断言
            ..Default::default()
        };
        // AAAA ×2、CCCC ×1（均为单窗序列，窗口数=1）→ dump -L 2 后仅 AAAA
        let fa = b">a\nAAAA\n>b\nAAAA\n>c\nCCCC\n";
        let t = filtered_counts(fa, &p);
        use trinity_common::kmer::kmer_to_intval;
        assert_eq!(t[&kmer_to_intval(b"AAAA").unwrap()], 2);
        assert!(!t.contains_key(&kmer_to_intval(b"CCCC").unwrap()));
        // min_kmer_cov_const=1 → 全表（"count>=1 全在表 + 查表缺失按 1"的等价路径 sanity）
        p.min_kmer_cov_const = 1;
        assert_eq!(
            filtered_counts(fa, &p)[&kmer_to_intval(b"CCCC").unwrap()],
            1
        );
    }

    // ---------------- SS stats 恒 canonical（先 -L 过滤后合并） ----------------

    /// SS 计数(literal) → -L2 过滤 → canonical 合并求和。ACGA=2 保留、其互补 TCGT=1 被
    /// dump -L 丢弃 → 合并后 canonical 键（ACGA/TCGT 的编码序大者）计数 **2**;
    /// 若（错误地）先合并再过滤会得 3。canonical("ACGA") vs ("TCGT"):
    /// G=0,A=1,T=2,C=3 编码下 TCGT(2,1,3,2)=158 > ACGA(1,1,3,1)=93 → 键为 TCGT。
    #[test]
    fn ss_stats_query_table_merges_after_l2_filter() {
        use trinity_common::kmer::kmer_to_intval;
        let p = DigiNormParams {
            k: 4,
            ds: false, // SS: literal 计数
            ..Default::default()
        };
        let fa = b">a\nACGA\n>b\nACGA\n>c\nTCGT\n";
        let table = stats_query_table(fa, &p);
        assert_eq!(table.len(), 1);
        assert_eq!(table[&kmer_to_intval(b"TCGT").unwrap()], 2);
        // 对照: dump -L 表（未合并）仍是 literal 键
        let dump = filtered_counts(fa, &p);
        assert_eq!(dump[&kmer_to_intval(b"ACGA").unwrap()], 2);
        assert!(!dump.contains_key(&kmer_to_intval(b"TCGT").unwrap()));
        // DS 路径: 表已 canonical，stats_query_table 与 filtered_counts 恒等
        let pds = DigiNormParams {
            k: 4,
            ds: true,
            ..Default::default()
        };
        assert_eq!(stats_query_table(fa, &pds), filtered_counts(fa, &pds));
    }

    /// SS stats 查询恒 canonical: 互补链 read（rc）经 ds=true 查表命中合并计数。
    /// ">r\nACGA" 与 ">s\nTCGT"（rc 关系）各 ×2 → canonical 计数 4 → stats median 4。
    #[test]
    fn ss_stats_rows_query_canonically() {
        use trinity_common::kmer::kmer_to_intval;
        let p = DigiNormParams {
            k: 4,
            ds: false,
            ..Default::default()
        };
        let fa = b">r\nACGA\n>r2\nACGA\n>s\nTCGT\n>s2\nTCGT\n";
        let table = stats_query_table(fa, &p);
        assert_eq!(table[&kmer_to_intval(b"TCGT").unwrap()], 4);
        let rows = sorted_stats_rows(fa, &table, &p).unwrap();
        // 两 read 的全部窗口查同一 canonical 键 → median 4（buggy 版: r 命中 literal 2、
        // s 查 literal 缺失 clamp 1）
        for r in &rows {
            assert_eq!(r.median, 4.0, "acc={}", r.acc);
        }
    }

    // ---------------- stats 文本往返与排序 ----------------

    #[test]
    fn stats_lines_match_write_stats_tsv_and_sort_by_acc_bytes() {
        use crate::coverage_stats::CoverageStatsRow;
        let rows = vec![
            CoverageStatsRow {
                acc: "z/1".into(),
                median: 2,
                mean: 2.0,
                stdev: 0.0,
            },
            CoverageStatsRow {
                acc: "a/1".into(),
                median: 1,
                mean: 1.4,
                stdev: 0.8,
            },
        ];
        let lines = stats_lines(&rows);
        assert_eq!(
            lines,
            vec!["z/1\t2\t2\t0\tthread:0", "a/1\t1\t1.4\t0.8\tthread:0"]
        );
        let mut sorted = lines;
        sorted.sort_by(|a, b| {
            acc_field(a)
                .cmp(acc_field(b))
                .then_with(|| a.as_bytes().cmp(b.as_bytes()))
        });
        assert!(sorted[0].starts_with("a/1"));
        // 字节序: 'Z'(0x5A) < 'a'(0x61)
        let mut l = vec!["a/1\t1", "Z/1\t2"];
        l.sort_by(|a, b| acc_field(a).cmp(acc_field(b)));
        assert_eq!(l, vec!["Z/1\t2", "a/1\t1"]);
    }

    /// 第一次文本往返: mean/stdev 经 %g6 文本再 parse（f32 直接 as f64 会保留
    /// 尾随精度，0.894427f32 as f64 = 0.8944269895553589 ≠ parse("0.894427")）。
    #[test]
    fn parse_stats_line_round_trips_g6_text() {
        let r = parse_stats_line("r1\t1\t0.894427\t-0\tthread:0");
        assert_eq!(r.acc, "r1");
        assert_eq!(r.median, 1.0);
        assert_eq!(r.mean, "0.894427".parse::<f64>().unwrap());
        assert_ne!(r.mean, 0.894427f32 as f64);
        assert_eq!(r.stdev, 0.0); // "-0" 归一为 +0（perl 数值化丢符号，见 parse_num 文档）
        assert!(r.stdev.is_sign_positive());
        // nan 文本（原版短 read 的 stdev 列为 "-nan"）
        let r = parse_stats_line("s\t0\t0\t-nan\tthread:0");
        assert!(r.stdev.is_nan());
        assert!(r.median == 0.0);
    }

    #[test]
    fn merged_rows_round_trip_percent_1f_text() {
        let merged = vec![MergedRow {
            core: "p".into(),
            median: "2.5".into(),
            mean: "2.0".into(),
            stdev: "NaN".into(),
        }];
        let rows = merged_to_stats_rows(&merged);
        assert_eq!(rows[0].acc, "p");
        assert_eq!(rows[0].median, 2.5);
        assert_eq!(rows[0].mean, 2.0);
        assert!(rows[0].stdev.is_nan());
    }

    // ---------------- 名单索引 ----------------

    #[test]
    fn index_keys_strip_trailing_slash_wordchar() {
        let idx = index_keys(&["rd".into(), "x/3".into(), "y/12".into(), "/1".into()]);
        assert!(idx.contains("rd"));
        // s|/\w$||: 任意单字 word 字符都剥（比 core_acc 的 [12] 宽）
        assert!(idx.contains("x"));
        // "y/12" 尾二字符 "/2" 前是 '1' 不是 '/' → 不剥
        assert!(idx.contains("y/12"));
        // "/1" 整体被剥 → 空键
        assert!(idx.contains(""));
    }

    // ---------------- fq 字节级抽取 ----------------

    /// 4 行原样透传（含 \r、'+' 行内容）; 空 seq 行也是 4 行之一; 末行无 \n 原样。
    #[test]
    fn extract_fq_byte_transparent() {
        // CRLF 的 header 行 + 末条记录 qual 无换行（原版 get_fastq_record 原文输出）
        let data = b"@r1/1 desc\r\nACGT\n+r1\nIIII\n@r2/1\nTT\n+\nII";
        let idx = index_keys(&["r1".into()]);
        let mut out = Vec::new();
        let mut found = HashSet::new();
        extract_fq_records(data, &idx, &mut out, &mut found).unwrap();
        assert_eq!(out, b"@r1/1 desc\r\nACGT\n+r1\nIIII\n");
        let idx = index_keys(&["r2".into()]);
        let mut out = Vec::new();
        let mut found = HashSet::new();
        extract_fq_records(data, &idx, &mut out, &mut found).unwrap();
        assert_eq!(out, b"@r2/1\nTT\n+\nII"); // 无尾换行原样透传
                                              // 记录内空行: 作为 4 行之一透传（seq 位空行 → 原版无长度校验照样输出）
        let data = b"@r4/1\n\n+\nII\n";
        let idx = index_keys(&["r4".into()]);
        let mut out = Vec::new();
        let mut found = HashSet::new();
        extract_fq_records(data, &idx, &mut out, &mut found).unwrap();
        assert_eq!(out, b"@r4/1\n\n+\nII\n");
    }

    /// 少于 4 行 / 首行非 '@' → Err（镜像 Perl confess）。
    #[test]
    fn extract_fq_errors_mirror_perl() {
        let idx = index_keys(&[]);
        // 尾部多余空行 → 读到 1 行即 EOF
        let mut out = Vec::new();
        let mut found = HashSet::new();
        let e = extract_fq_records(b"@r\nAC\n+\nII\n\n", &idx, &mut out, &mut found);
        assert!(e.is_err());
        let e = extract_fq_records(b"r\nAC\n+\nII\n", &idx, &mut out, &mut found);
        assert!(e.is_err());
        // 恰好 3 行的残尾
        let e = extract_fq_records(b"@r\nAC\n+\n", &idx, &mut out, &mut found);
        assert!(e.is_err());
    }

    /// _forward/_reverse 与 /1|/2 的 core 名匹配（L534-535 + Fastq_reader.pm）。
    #[test]
    fn extract_fq_core_name_forward_reverse() {
        let idx = index_keys(&["xandmore".into()]);
        let data = b"@x_forwardandmore_reverse/1\nAC\n+\nII\n";
        let mut out = Vec::new();
        let mut found = HashSet::new();
        extract_fq_records(data, &idx, &mut out, &mut found).unwrap();
        assert_eq!(out.len(), data.len());
    }

    // ---------------- fa 抽取 ----------------

    /// header 原样保留（含描述）、序列单行化、多行序列拼接、acc 剥 /1|/2。
    #[test]
    fn extract_fa_record_formatting() {
        let idx = index_keys(&["a".into(), "b".into()]);
        let data = b">a desc line\nACGT\nACGT\n>b/1\nTT\n";
        let mut out = Vec::new();
        let mut found = HashSet::new();
        extract_fa_records(data, &idx, &mut out, &mut found).unwrap();
        assert_eq!(out, b">a desc line\nACGTACGT\n>b/1\nTT\n");
        assert!(found.contains("a"));
        // ">" 在序列行中间不切（分隔符是 "\n>"）
        let data = b">c/2 x\nAC>GT\n";
        let idx2 = index_keys(&["c".into()]);
        let mut out = Vec::new();
        let mut found = HashSet::new();
        extract_fa_records(data, &idx2, &mut out, &mut found).unwrap();
        assert_eq!(out, b">c/2 x\nAC>GT\n");
        // 控制字符（\r 等）删除; 前导空块跳过
        let data = b"\n>d/1\r\nAC\rGT\n";
        let idx3 = index_keys(&["d".into()]);
        let mut out = Vec::new();
        let mut found = HashSet::new();
        extract_fa_records(data, &idx3, &mut out, &mut found).unwrap();
        assert_eq!(out, b">d/1\nACGT\n");
    }

    // ---------------- 输出命名 ----------------

    #[test]
    fn normalized_out_path_matches_original_pattern() {
        let p = DigiNormParams::default();
        let path = normalized_out_path(Path::new("/o"), "pe.l.fq", SeqFormat::Fq, &p);
        assert_eq!(
            path.to_str().unwrap(),
            "/o/pe.l.fq.normalized_K25_maxC200_minC1_maxCV10000.fq"
        );
        let p2 = DigiNormParams {
            max_cov: 50.0,
            ..Default::default()
        };
        let path = normalized_out_path(Path::new("/o"), "s.fa", SeqFormat::Fa, &p2);
        assert!(path
            .to_str()
            .unwrap()
            .ends_with("s.fa.normalized_K25_maxC50_minC1_maxCV10000.fa"));
    }
}
