//! 汇总输出（P5-T3）——原版 `print_butterfly_assemblies.pl`（全序列字符串去重
//! 与 min_seq_length 过滤; 缺失 allProbPaths 文件 → stderr + 退出码 1）以及
//! `get_Trinity_gene_to_trans_map.pl`（gene_trans_map 由 orchestrate 层落）。

use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};

use trinity_common::error::CommonError;

/// 逐组件收 allProbPaths.fasta → `final_path`。
///
/// * header 原样保留, 序列以 60 列折行输出;
/// * `len >= min_len` 才收录; 全序列字符串去重（重复 → stderr 警告并排除）;
/// * 任一组件的 `<base>.graph.allProbPaths.fasta` 缺失 → 处理完其余后 Err
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

    for (id, base) in listing {
        let p = PathBuf::from(format!("{}.graph.allProbPaths.fasta", base.display()));
        let text = match fs::read_to_string(&p) {
            Ok(t) => t,
            Err(_) => {
                eprintln!("Error, no fasta file reported as: {}", p.display());
                missing.push(p.display().to_string());
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

/// 简单 FASTA 解析: (header 行原文含 '>', 去空白序列)。
fn records(text: &str) -> Vec<(String, String)> {
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
