//! `Chrysalis/analysis/CreateIwormFastaBundle.cc` 的移植——把
//! BubbleUpClustering 的 COMPONENT 块折叠成捆绑 FASTA（`>s_<no> <cov>...` +
//! `X` 连接的多序列一行）。
//!
//! 过滤（组件号留空洞，**不重编号**——下游 ReadsToComponents 按号对齐）：
//! 单 contig 且长 < min；或组件总长 < min。

use trinity_common::error::CommonError;

/// get_iworm_coverage（:112-137）：`[iworm>a1;43_total_counts:_...]` → `43`。
///
/// 前缀必须是 `[iworm>a`（否则 Err）；取首个 `';'` 到首个 `'_'` 之间
/// （start > end 时 Err，镜像原版 exit(5)）。
pub fn get_iworm_coverage(iworm_info: &str) -> Result<String, CommonError> {
    if !iworm_info.starts_with("[iworm>a") {
        return Err(CommonError::Parse(format!(
            "Error, iworm_info: {iworm_info} doesn't start with iworm>a"
        )));
    }
    let start = iworm_info.find(';');
    let end = iworm_info.find('_');
    match (start, end) {
        (Some(s), Some(e)) if s <= e => Ok(iworm_info[s + 1..e].to_string()),
        _ => Err(CommonError::Parse(format!(
            "Error extracting coverage info from: {iworm_info}"
        ))),
    }
}

/// CreateIwormFastaBundle 主流程：返回捆绑 FASTA 全文（写文件内容的等价物）。
#[allow(clippy::while_let_on_iterator)] // 外/内两层消费同一行迭代器
pub fn create_iworm_fasta_bundle(
    component_out: &str,
    min_len: usize,
) -> Result<String, CommonError> {
    let mut out = String::new();
    let mut lines = component_out.lines().peekable();

    // 镜像原版外层 while (ParseLine) + 内层 while (ParseLine) until END 的控制流：
    // token 化按空白切分；空行与 token0 以 '#' 开头的行跳过。
    while let Some(line) = lines.next() {
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.is_empty() || toks[0].starts_with('#') {
            continue;
        }
        if toks[0] != "COMPONENT" {
            continue;
        }
        let component_no: i64 = toks[1].parse().unwrap_or(0); // AsInt
        let _num_iworm_contigs = toks[2]; // 读了不用（原版 quirk）

        let mut seqs: Vec<String> = Vec::new(); // tmpSeq
        let mut cov_vals = String::new();

        while let Some(line) = lines.next() {
            let toks: Vec<&str> = line.split_whitespace().collect();
            if toks.is_empty() || toks[0].starts_with('#') {
                continue;
            }
            if toks[0] == "END" {
                break;
            }
            if toks[0].starts_with('>') {
                // '>' 行：第 4 列 = [iworm>...]
                let iworm_info = toks[3];
                cov_vals.push(' ');
                cov_vals.push_str(&get_iworm_coverage(iworm_info)?);
                seqs.push(String::new());
                continue;
            }
            // 序列行：只取 token0（AsString(0)），追加到当前 iworm
            if let Some(last) = seqs.last_mut() {
                last.push_str(toks[0]);
            }
        }

        // 单 contig 太短 → 丢
        if seqs.len() == 1 && seqs[0].len() < min_len {
            continue;
        }
        // 总长太短 → 丢（组件号留空洞，不重编号）
        if seqs.iter().map(|s| s.len()).sum::<usize>() < min_len {
            continue;
        }

        out.push_str(&format!(">s_{component_no}{cov_vals}\n"));
        out.push_str(&seqs.join("X"));
        out.push('\n');
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coverage_extraction() {
        assert_eq!(
            get_iworm_coverage("[iworm>a1;43_total_counts:_59920_Seed:_57]").unwrap(),
            "43"
        );
        assert_eq!(get_iworm_coverage("[iworm>a17;123.5_x").unwrap(), "123.5");
    }

    #[test]
    fn coverage_errors() {
        assert!(get_iworm_coverage("[iworm>b1;43_x").is_err()); // 非 'a' 前缀
        assert!(get_iworm_coverage("a1;43_x").is_err());
        assert!(get_iworm_coverage("[iworm>a143_x").is_err()); // 无 ';' → start 不存在
    }

    fn block(no: usize, entries: &[(&str, &str)]) -> String {
        let mut s = format!("COMPONENT {no}\t{}\n", entries.len());
        for (name, seq) in entries {
            s.push_str(&format!(
                ">Component_{no} {} 0 [iworm>{name}]\n",
                entries.len()
            ));
            s.push_str(seq);
            s.push('\n');
        }
        s.push_str("#POOL_INFO\nEND\n");
        s
    }

    #[test]
    fn bundle_joins_with_x_and_covs_in_order() {
        let inp = format!(
            "{}{}",
            block(
                0,
                &[("a1;10_len", "AAAACCCCGGGGTTTT"), ("a2;20_len", "AAAA")]
            ),
            block(3, &[("a3;30_len", "CCCC")]),
        );
        let out = create_iworm_fasta_bundle(&inp, 4).unwrap();
        assert_eq!(out, ">s_0 10 20\nAAAACCCCGGGGTTTTXAAAA\n>s_3 30\nCCCC\n");
    }

    #[test]
    fn multi_line_sequences_concatenated() {
        let inp = "COMPONENT 0\t1\n>Component_0 1 0 [iworm>a1;5_l]\nAAAA\nCCCC\nGGGG\nEND\n";
        let out = create_iworm_fasta_bundle(inp, 6).unwrap();
        assert_eq!(out, ">s_0 5\nAAAACCCCGGGG\n");
    }

    #[test]
    fn filters_leave_component_holes() {
        // 组件 0：单 contig 长 4 < 8 → 丢；组件 1：两 contig 总 6 < 8 → 丢；
        // 组件 2：单 contig 长 8 → 保 → 输出仍叫 s_2
        let inp = format!(
            "{}{}{}",
            block(0, &[("a1;1_l", "AAAA")]),
            block(1, &[("a2;1_l", "AAAA"), ("a3;1_l", "CC")]),
            block(2, &[("a4;1_l", "AAAACCCCGGGGTTTT")]),
        );
        let out = create_iworm_fasta_bundle(&inp, 8).unwrap();
        assert_eq!(out, ">s_2 1\nAAAACCCCGGGGTTTT\n");
    }

    #[test]
    fn skips_blank_and_comment_lines() {
        let inp =
            "\n#comment\nCOMPONENT 0\t1\n#mid\n\n>Component_0 1 0 [iworm>a1;7_l]\nAAAA\nEND\n";
        let out = create_iworm_fasta_bundle(inp, 2).unwrap();
        assert_eq!(out, ">s_0 7\nAAAA\n");
    }
}
