//! prep_seqs 移植（原版 Trinity L2746-2879 + 主流程 L1604-1670）:
//! fq → fa（seqtk-trinity 语义，复用 `trinity_kmer::read_names::fq_records_to_fa`），
//! fa → 字节拼接（SS 端 'R' 时 revcomp_fasta.pl 镜像——复用
//! `trinity_kmer::diginorm::prep_side`），多文件列表序拼接;
//! PE 再 cat left.fa right.fa → both.fa（先左后右，字节数校验）。
//! 各步 `.ok` 断点存在即跳过。

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

use trinity_common::error::CommonError;
use trinity_kmer::diginorm::{prep_side, SeqFormat};
use trinity_kmer::read_names::ReadType;

use crate::args::{SeqType, TrinityArgs};

/// 读入一个输入文件（支持 .gz——原版逐文件 gunzip -c 语义）。
fn read_input(p: &Path) -> Result<Vec<u8>, CommonError> {
    let mut buf = Vec::new();
    trinity_common::io_util::open_maybe_gz(p)?.read_to_end(&mut buf)?;
    Ok(buf)
}

fn fmt_of(t: SeqType) -> SeqFormat {
    match t {
        SeqType::Fq => SeqFormat::Fq,
        SeqType::Fa => SeqFormat::Fa,
    }
}

/// fa 记录数（trinity_target_fa.read_count 等价: 行首 '>' 计数）。
fn fa_record_count(data: &[u8]) -> u64 {
    let mut n = 0u64;
    let mut at_line_start = true;
    for &b in data {
        if at_line_start && b == b'>' {
            n += 1;
        }
        at_line_start = b == b'\n';
    }
    n
}

/// 单侧 prep_seqs: 列表序逐文件转换拼接 → `<prefix>.fa` + `.ok`。
/// 断点: `<prefix>.fa.ok` 存在即跳过（原版 `return if -e "$file_prefix.fa.ok"`）。
fn prep_side_files(
    files: &[PathBuf],
    seq_type: SeqType,
    prefix: &str,
    revcomp: bool,
    outdir: &Path,
) -> Result<(), CommonError> {
    let fa = outdir.join(format!("{prefix}.fa"));
    let ok = outdir.join(format!("{prefix}.fa.ok"));
    // 断点: .ok 存在**且产物在**即跳过（收尾清理删除产物后, resume 自动重建）。
    if crate::checkpoint::checkpoint_exists(&ok) && fa.exists() {
        eprintln!("---- checkpoint found, skipping: {}", ok.display());
        return Ok(());
    }
    let mut data = Vec::new();
    for f in files {
        let raw = read_input(f)?;
        // eval 块语义: 失败时移除半成品 fa。
        match prep_side(&raw, fmt_of(seq_type), read_type_of(prefix), revcomp) {
            Ok(part) => data.extend_from_slice(&part),
            Err(e) => {
                let _ = fs::remove_file(&fa);
                return Err(e);
            }
        }
    }
    fs::write(&fa, &data)?;
    fs::write(&ok, b"")?;
    Ok(())
}

/// prep_seqs L2749: read_type = (prefix == "right") ? 2 : 1（/1、/2 追加）。
fn read_type_of(prefix: &str) -> ReadType {
    if prefix == "right" {
        ReadType::R2
    } else {
        ReadType::R1
    }
}

/// PE SS 端字符（split(//, ss)）: 只有 'R' 才 revcomp（L1611）。
fn side_revcomp(ss: &Option<String>, idx: usize) -> bool {
    ss.as_ref()
        .and_then(|s| s.as_bytes().get(idx).map(|&b| b == b'R'))
        .unwrap_or(false)
}

/// 主入口。返回 `(trinity_target_fa 路径, read count)`。
/// PE: both.fa（= left.fa ++ right.fa，先左后右）; SE: single.fa。
pub fn prep_reads(args: &TrinityArgs, outdir: &Path) -> Result<(PathBuf, u64), CommonError> {
    fs::create_dir_all(outdir)?;

    if !args.single.is_empty() {
        // SE: SS 整串即端类型（'R' → revcomp）。
        prep_side_files(
            &args.single,
            args.seq_type,
            "single",
            args.ss_lib_type.as_deref() == Some("R"),
            outdir,
        )?;
        let target = outdir.join("single.fa");
        let data = fs::read(&target)?;
        return Ok((target, fa_record_count(&data)));
    }

    // PE
    prep_side_files(
        &args.left,
        args.seq_type,
        "left",
        side_revcomp(&args.ss_lib_type, 0),
        outdir,
    )?;
    prep_side_files(
        &args.right,
        args.seq_type,
        "right",
        side_revcomp(&args.ss_lib_type, 1),
        outdir,
    )?;

    let left_fa = outdir.join("left.fa");
    let right_fa = outdir.join("right.fa");
    let target = outdir.join("both.fa");
    let target_ok = outdir.join("both.fa.ok");

    let left_data = fs::read(&left_fa)?;
    let right_data = fs::read(&right_fa)?;
    let expected = left_data.len() + right_data.len();
    let need_rebuild = !crate::checkpoint::checkpoint_exists(&target_ok)
        || fs::metadata(&target).map(|m| m.len() as usize).unwrap_or(0) != expected;
    if need_rebuild {
        let mut both = Vec::with_capacity(expected);
        both.extend_from_slice(&left_data);
        both.extend_from_slice(&right_data);
        fs::write(&target, &both)?;
        // L1645-1648: 字节数校验（失败 die）。
        let got = fs::metadata(&target)?.len() as usize;
        if got != expected {
            return Err(CommonError::Parse(format!(
                "both.fa is smaller ({got} bytes) than the combined size of left.fa and right.fa ({expected} bytes)"
            )));
        }
        fs::write(&target_ok, b"")?;
    }
    Ok((
        target,
        fa_record_count(&left_data) + fa_record_count(&right_data),
    ))
}
