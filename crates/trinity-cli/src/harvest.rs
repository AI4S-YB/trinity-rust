//! 汇总输出（P5-T3）——原版 `print_butterfly_assemblies.pl`（全序列字符串去重
//! 与 min_seq_length 过滤; 缺失 allProbPaths 记录 → stderr + 退出码 1）以及
//! `get_Trinity_gene_to_trans_map.pl`（gene_trans_map 由 orchestrate 层落）。
//!
//! allProbPaths 从每 Cbin 的 `butterfly.pack` 按索引偏移读（IO 优化,
//! 替代 per-component 小文件; 见 trinity-chrysalis pack 模块）。

use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};

use trinity_chrysalis::pack::{read_at, PackIndex, BUTTERFLY_PACK, KIND_ALL_PROB_PATHS};
use trinity_common::error::CommonError;

/// 逐组件收 allProbPaths → `final_path`。
///
/// * header 原样保留, 序列以 60 列折行输出;
/// * `len >= min_len` 才收录; 全序列字符串去重（重复 → stderr 警告并排除）;
/// * 任一组件的 allProbPaths 记录缺失 → 处理完其余后 Err
///   （原版 stderr "Error, no fasta file reported as: ..." + exit 1）。
///
/// 返回收录的转录本数。
pub fn harvest(
    listing: &[(u64, PathBuf)],
    min_len: usize,
    final_path: &Path,
) -> Result<usize, CommonError> {
    let mut seen: HashSet<String> = HashSet::new();
    let mut out = String::new();
    let mut missing: Vec<String> = Vec::new();
    let mut n = 0usize;
    // listing 按 id 升序 → 同 Cbin 连续; 当前 bin 的索引与句柄各开一次。
    let mut cur_bin: Option<(PathBuf, PackIndex, fs::File)> = None;

    for (id, base) in listing {
        let dir = base.parent().map(|d| d.to_path_buf()).unwrap_or_default();
        if cur_bin.as_ref().map(|(d, _, _)| d != &dir).unwrap_or(true) {
            cur_bin = load_bin(&dir).ok();
            if cur_bin.is_none() {
                // 冷路径（包缺失/损坏）: 逐组件报 missing, 同 bin 的下一组件会重试加载
                let logical = format!("{}.graph.allProbPaths.fasta", base.display());
                eprintln!("Error, no fasta file reported as: {logical}");
                missing.push(logical);
                continue;
            }
        }
        let (_, idx, file) = cur_bin.as_ref().unwrap();
        let logical = format!("{}.graph.allProbPaths.fasta", base.display());
        let text = match idx.find(*id, KIND_ALL_PROB_PATHS) {
            Some((offset, len)) => match read_at(file, offset, len) {
                Ok(bytes) => match String::from_utf8(bytes) {
                    Ok(t) => t,
                    Err(_) => {
                        eprintln!("Error, no fasta file reported as: {logical}");
                        missing.push(logical);
                        continue;
                    }
                },
                Err(_) => {
                    eprintln!("Error, no fasta file reported as: {logical}");
                    missing.push(logical);
                    continue;
                }
            },
            None => {
                eprintln!("Error, no fasta file reported as: {logical}");
                missing.push(logical);
                continue;
            }
        };
        for (header, seq) in records(&text) {
            if seq.len() < min_len {
                continue;
            }
            if !seen.insert(seq.clone()) {
                eprintln!("-duplicate sequence detected, excluding it. (component {id})");
                continue;
            }
            out.push('>');
            out.push_str(&header);
            out.push('\n');
            for chunk in seq.as_bytes().chunks(60) {
                out.push_str(std::str::from_utf8(chunk).unwrap_or(""));
                out.push('\n');
            }
            n += 1;
        }
    }

    if !missing.is_empty() {
        return Err(CommonError::Parse(format!(
            "Error, {} butterfly allProbPaths fasta file(s) missing (first: {})",
            missing.len(),
            missing[0]
        )));
    }
    fs::write(final_path, out)?;
    Ok(n)
}

/// 打开一个 Cbin 的 butterfly 包（索引 + 数据句柄）。
pub(crate) fn load_bin(dir: &Path) -> Result<(PathBuf, PackIndex, fs::File), CommonError> {
    let idx = PackIndex::load(&dir.join(format!("{BUTTERFLY_PACK}.idx")))?;
    let file = fs::File::open(dir.join(format!("{BUTTERFLY_PACK}.pack")))?;
    Ok((dir.to_path_buf(), idx, file))
}

/// 简单 FASTA 解析: (header 行原文含 '>', 去空白序列)。
pub(crate) fn records(text: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut header: Option<String> = None;
    let mut seq = String::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if let Some(h0) = header.take() {
                out.push((h0, std::mem::take(&mut seq)));
            }
            header = Some(h.to_string());
        } else if header.is_some() {
            seq.push_str(line.trim());
        }
    }
    if let Some(h0) = header {
        out.push((h0, seq));
    }
    out
}
