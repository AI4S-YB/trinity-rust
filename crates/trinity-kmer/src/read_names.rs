//! seqtk-trinity 的 read 名规范化（update_read_name_for_Trinity）+ fq→fa 预处理移植。
//!
//! oracle: trinity-plugins/seqtk-trinity（seqtk.c，本机已编译）。
//! 管线调用形态（insilico_read_normalization.pl L688-695）:
//! `seqtk-trinity seq -A -R <1|2> [-r] file >> out.fa`，`-r` 仅 SS_lib_type == "R" 时加。
//!
//! 名字规则（seqtk.c:201-290，按原序）:
//! 1. 名字含 "_forward"/"_reverse"（**任意位置**首次出现，_forward 优先）→ **截断**到该处;
//! 2. 未截断 且 尾二字符 "/1"|"/2" → 校验与期望类型一致（不一致 = 原版 exit(2)），原样返回;
//! 3. 未截断 且 comment 以 "1:"|"2:" 开头（Illumina 新式，comment 至少 2 字符）
//!    → 校验后追加 "/"+comment[0];
//! 4. 兜底 → 追加 "/1" 或 "/2"（截断过的名也走这里）。
//!
//! 注意 `_forward` 在 seqtk 是**截断**（`name_len = found - name`），而抽取阶段的
//! core 名（本模块 [`core_read_name`]，Perl 语义）是**删子串再拼接**——两者刻意不同。

use std::io::BufReader;

use trinity_common::error::CommonError;
use trinity_common::fastq::FastqReader;

/// seqtk.c:52 `read_type` 全局变量（-R 参数）: 1 = left，2 = right。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadType {
    R1 = 1,
    R2 = 2,
}

/// exit(2) 的错误消息（seqtk.c:230/256 逐字）: "Error, found read_type %c but expecting read_type %i"。
fn read_type_mismatch(found: char, expected: ReadType) -> CommonError {
    CommonError::FastqFormat(format!(
        "Error, found read_type {found} but expecting read_type {}",
        expected as u8
    ))
}

/// seqtk.c:201-290 update_read_name_for_Trinity 的逐字镜像。
///
/// 返回 Err 对应原版 exit(2)（管线级致命）。与 C 的差异:
/// - C 对长度 0/1 的名读 `name[-1]`（未定义行为）; 这里的 `len >= 2` 守卫是安全化，
///   短名直接落入 comment/兜底分支——与 C 的实际常见路径（越界字节恰为 '/' 且尾字符
///   为 '1'/'2' 才会分叉）行为一致。
/// - C 就地改 `name_copy[1000]`（名字 ≥ 1000 字节会溢出栈缓冲）; 这里返回新 String。
pub fn update_read_name_for_trinity(
    name: &str,
    comment: Option<&str>,
    read_type: ReadType,
) -> Result<String, CommonError> {
    // seqtk.c:209-219: strstr 找 "_forward"（优先）/"_reverse" 的首次出现（任意位置），
    // 找到即把 name_len 截到该处。found != NULL 时规则 2/3 一律跳过。
    let mut name = name.to_string();
    let truncated = if let Some(pos) = name.find("_forward").or_else(|| name.find("_reverse")) {
        name.truncate(pos);
        true
    } else {
        false
    };

    // seqtk.c:222-235 旧式: 倒数第二字符 '/' 且末字符 '1'|'2' → 校验后原样返回。
    // （truncated 时跳过。短名守卫见函数文档的 UB 说明。）
    if !truncated && name.len() >= 2 {
        let b = name.as_bytes();
        let last = b[name.len() - 1];
        if b[name.len() - 2] == b'/' && (last == b'1' || last == b'2') {
            let found = if last == b'1' {
                ReadType::R1
            } else {
                ReadType::R2
            };
            if found != read_type {
                return Err(read_type_mismatch(last as char, read_type));
            }
            return Ok(name);
        }
    }

    // seqtk.c:243-267 新式（Illumina）: comment 长度 > 1、comment[1] == ':'、
    // comment[0] ∈ {'1','2'} → 校验后追加 "/"+comment[0]。
    if !truncated {
        if let Some(c) = comment {
            let cb = c.as_bytes();
            if cb.len() > 1 && cb[1] == b':' && (cb[0] == b'1' || cb[0] == b'2') {
                let found = if cb[0] == b'1' {
                    ReadType::R1
                } else {
                    ReadType::R2
                };
                if found != read_type {
                    return Err(read_type_mismatch(cb[0] as char, read_type));
                }
                return Ok(format!("{name}/{}", cb[0] as char));
            }
        }
    }

    // seqtk.c:276-284 兜底: 追加期望的 "/1"|"/2"（截断过的名也走这里）。
    // 注意是数字字符 '1'|'2'（b'0' + 判别值），不是判别值本身。
    Ok(format!("{name}/{}", char::from(b'0' + read_type as u8)))
}

/// seqtk.c:163-172 comp_tab 逐字镜像: A<->T、C<->G、U->A/a，大小写保留，
/// IUPAC 互补（B<->V、D<->H、K<->M、R<->Y、S/W/N/X 自身），其余按原表（多为恒等）。
const COMP_TAB: [u8; 128] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, //
    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, //
    32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, //
    48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, //
    b'@', b'T', b'V', b'G', b'H', b'E', b'F', b'C', b'D', b'I', b'J', b'M', b'L', b'K', b'N',
    b'O', //
    b'P', b'Q', b'Y', b'S', b'A', b'A', b'B', b'W', b'X', b'R', b'Z', b'[', b'\\', b']', b'^',
    b'_', //
    b'@', b't', b'v', b'g', b'h', b'e', b'f', b'c', b'd', b'i', b'j', b'm', b'l', b'k', b'n',
    b'o', //
    b'p', b'q', b'y', b's', b'a', b'a', b'b', b'w', b'x', b'r', b'z', b'{', b'|', b'}', b'~',
    127, //
];

/// seqtk.c:538-551（flag&4）反向互补: 大小写保留（oracle 实测 ACGTacgt → acgtACGT）。
/// 字节 >= 0x80 原版以 signed char 负下标读表（UB），此处按原样保留。
fn revcomp_seq(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| if b < 128 { COMP_TAB[b as usize] } else { b })
        .collect()
}

/// kseq.h:111-113/182-183: name = header 到首个空白为止，其余为 comment（无分隔符则
/// None）。空白 = C `isspace`（kseq 逐字节调用）: ' '、'\t'、'\v'、'\f'、'\r'
/// （'\n' 不会出现在行内）。comment 存在与否不影响名字规则——只有规则 3 会读它，
/// 且要求长度 > 1。
fn split_name_comment(header: &str) -> (&str, Option<&str>) {
    match header.find([' ', '\t', '\x0b', '\x0c', '\r']) {
        Some(pos) => (&header[..pos], Some(&header[pos + 1..])),
        None => (header, None),
    }
}

/// fq→fa 转换，镜像 `seqtk-trinity seq -A -R <1|2> [-r]`（DigiNorm 预处理，
/// insilico_read_normalization.pl L688-695 的调用）:
/// 名字规范化 + 可选反向互补（-r），输出 `>name\nseq\n`，序列大小写原样保留。
///
/// **输入须为有效 UTF-8**（内部经 [`trinity_common::fastq::FastqReader`] 按 String 读行，
/// 与原版 seqtk 的字节流不同——非 UTF-8 header 会在读入层报错而非由本函数处理;
/// 抽取阶段 [`crate::diginorm`] 的字节级扫描器无此约束）。
///
/// - 类型不匹配 → Err（原版 exit(2)）;
/// - 解析失败 → Err（FastqReader 的 4 行严格校验）;
/// - 一条记录都读不到 → Err（原版 "no records were correctly parsed" exit(5)）。
///   与原版的流式差异: 原版逐条打印、在出错条目处中止（此前的条目已写出）;
///   本函数整体返回，Err 时调用方拿不到部分输出——DigiNorm 中 exit(2)/exit(5)
///   同样直接终止管线，语义等价。
pub fn fq_records_to_fa(
    input: &[u8],
    read_type: ReadType,
    revcomp: bool,
) -> Result<Vec<u8>, CommonError> {
    let mut reader = FastqReader::new(BufReader::new(input));
    let mut out = Vec::new();
    let mut n_seqs = 0usize;
    while let Some(rec) = reader.next_record()? {
        n_seqs += 1;
        let (name, comment) = split_name_comment(&rec.header);
        let new_name = update_read_name_for_trinity(name, comment, read_type)?;
        out.reserve(new_name.len() + rec.sequence.len() + 3);
        out.push(b'>');
        out.extend_from_slice(new_name.as_bytes());
        out.push(b'\n');
        if revcomp {
            out.extend_from_slice(&revcomp_seq(rec.sequence.as_bytes()));
        } else {
            out.extend_from_slice(rec.sequence.as_bytes());
        }
        out.push(b'\n');
    }
    if n_seqs == 0 {
        // seqtk.c:564-567: "Error, no records were correctly parsed from <file>" → exit(5)。
        return Err(CommonError::FastqFormat(
            "Error, no records were correctly parsed from input".to_string(),
        ));
    }
    Ok(out)
}

/// 抽取阶段的 core 名（make_normalized_reads_file 的 fq 路径）:
/// Fastq_reader.pm L113-116 `^(\S+)/([12])$`（剥**尾部** "/1"|"/2"，要求前面至少
/// 1 个字符）+ insilico_read_normalization.pl L534-535 `s/_forward//; s/_reverse//`
/// （两条独立语句，各删**首个**出现的子串后**拼接**——不是截断!）。
///
/// perl 实测: "x_forwardandmore_reverse" → "xandmore"（与 seqtk 的截断语义不同）。
pub fn core_read_name(header_first_token: &str) -> String {
    let mut name = header_first_token.to_string();
    // Fastq_reader.pm L113-116: ^(\S+)/([12])$ —— 剥尾部 "/1"|"/2"（\S+ 要求至少
    // 1 个字符，故 "/1" 整体保留）; strip_suffix 语义即贪婪回溯到最后一节。
    if name.len() > 2 {
        if let Some(core) = name.strip_suffix("/1").or_else(|| name.strip_suffix("/2")) {
            name = core.to_string();
        }
    }
    // insilico_read_normalization.pl L534-535: s/_forward//; s/_reverse// ——
    // 两条独立语句（无 /g），各删**首个**出现的子串后拼接。
    for pat in ["_forward", "_reverse"] {
        if let Some(pos) = name.find(pat) {
            name.replace_range(pos..pos + pat.len(), "");
        }
    }
    name
}

#[cfg(test)]
mod tests {
    use super::*;

    fn err_msg(e: CommonError) -> String {
        e.to_string()
    }

    // ---------------- update_read_name_for_trinity ----------------

    /// oracle 实测（seq -A -R 1，4 种 header 形态）:
    /// rd1 1:N:0:1 → rd1/1（新式）; rd2/1 → rd2/1（旧式原样）;
    /// rd3_forward → rd3/1（截断+兜底）; rd4 → rd4/1（兜底）。
    #[test]
    fn four_header_forms_oracle_locked() {
        assert_eq!(
            update_read_name_for_trinity("rd1", Some("1:N:0:1"), ReadType::R1).unwrap(),
            "rd1/1"
        );
        assert_eq!(
            update_read_name_for_trinity("rd2/1", None, ReadType::R1).unwrap(),
            "rd2/1"
        );
        assert_eq!(
            update_read_name_for_trinity("rd3_forward", None, ReadType::R1).unwrap(),
            "rd3/1"
        );
        assert_eq!(
            update_read_name_for_trinity("rd4", None, ReadType::R1).unwrap(),
            "rd4/1"
        );
    }

    /// 旧式 /1|/2 尾巴与期望类型不符 → exit(2)（消息逐字）。
    #[test]
    fn old_format_type_mismatch_is_error() {
        let e = update_read_name_for_trinity("rd2/1", None, ReadType::R2).unwrap_err();
        assert_eq!(
            err_msg(e),
            "fastq format error: Error, found read_type 1 but expecting read_type 2"
        );
        let e = update_read_name_for_trinity("rd5/2", None, ReadType::R1).unwrap_err();
        assert_eq!(
            err_msg(e),
            "fastq format error: Error, found read_type 2 but expecting read_type 1"
        );
        // 匹配时即便 comment 声明另一类型也不查（规则 2 先于规则 3 返回）
        assert_eq!(
            update_read_name_for_trinity("rd/1", Some("2:N:0:1"), ReadType::R1).unwrap(),
            "rd/1"
        );
    }

    /// 新式 comment（Illumina）类型不符 → exit(2)。oracle: rd5 2:N:0:1 + -R 1。
    #[test]
    fn new_format_type_mismatch_is_error() {
        let e = update_read_name_for_trinity("rd5", Some("2:N:0:1"), ReadType::R1).unwrap_err();
        assert_eq!(
            err_msg(e),
            "fastq format error: Error, found read_type 2 but expecting read_type 1"
        );
    }

    /// _forward/_reverse 在**任意位置**首次出现即截断（_forward 优先），截断后
    /// 跳过规则 2/3 直接兜底追加期望类型。oracle 实测（-R 1 与 -R 2）。
    #[test]
    fn forward_reverse_truncation_anywhere() {
        assert_eq!(
            update_read_name_for_trinity("x_forwardandmore_reverse", None, ReadType::R1).unwrap(),
            "x/1"
        );
        assert_eq!(
            update_read_name_for_trinity("mid_forward_x", None, ReadType::R1).unwrap(),
            "mid/1"
        );
        assert_eq!(
            update_read_name_for_trinity("mid2_reverse_x", None, ReadType::R1).unwrap(),
            "mid2/1"
        );
        assert_eq!(
            update_read_name_for_trinity("rd3_forward", None, ReadType::R2).unwrap(),
            "rd3/2"
        );
        // 截断后旧式尾巴被覆盖、comment 被忽略
        assert_eq!(
            update_read_name_for_trinity("rd_forward/1", None, ReadType::R1).unwrap(),
            "rd/1"
        );
        assert_eq!(
            update_read_name_for_trinity("rd_forward", Some("2:N:0:1"), ReadType::R1).unwrap(),
            "rd/1"
        );
    }

    /// 尾字符非 1|2（如 "/3"）不满足规则 2 → 兜底追加。oracle: rd6/3 → rd6/3/1。
    #[test]
    fn suffix_not_1_or_2_falls_back() {
        assert_eq!(
            update_read_name_for_trinity("rd6/3", None, ReadType::R1).unwrap(),
            "rd6/3/1"
        );
    }

    /// comment 长度边界: 规则 3 要求 comment 长度 > 1（C: comment_length > 1）。
    /// oracle: rd7 x → rd7/1; rd8 1 → rd8/1; rd9 2: → exit(2)（恰 2 字符即触发规则 3）。
    #[test]
    fn comment_length_edges() {
        assert_eq!(
            update_read_name_for_trinity("rd7", Some("x"), ReadType::R1).unwrap(),
            "rd7/1"
        );
        assert_eq!(
            update_read_name_for_trinity("rd8", Some("1"), ReadType::R1).unwrap(),
            "rd8/1"
        );
        assert!(update_read_name_for_trinity("rd9", Some("2:"), ReadType::R1).is_err());
        assert_eq!(
            update_read_name_for_trinity("rd9", Some("2:"), ReadType::R2).unwrap(),
            "rd9/2"
        );
        // 空 comment（header 以空白结尾）等价于无 comment
        assert_eq!(
            update_read_name_for_trinity("rd", Some(""), ReadType::R1).unwrap(),
            "rd/1"
        );
    }

    /// 短名/空名: C 读 name[-1] 属 UB，这里安全化为跳过规则 2（短名走 comment/兜底）。
    /// 空名兜底结果 "/1" 与 C 的实际写入一致（name[0]='/'; name[1]='1'）。
    #[test]
    fn short_and_empty_names_take_fallback() {
        assert_eq!(
            update_read_name_for_trinity("A", None, ReadType::R1).unwrap(),
            "A/1"
        );
        assert_eq!(
            update_read_name_for_trinity("1", None, ReadType::R1).unwrap(),
            "1/1"
        );
        assert_eq!(
            update_read_name_for_trinity("", None, ReadType::R2).unwrap(),
            "/2"
        );
        // 短名仍可走规则 3（长度守卫只在规则 2）
        assert_eq!(
            update_read_name_for_trinity("x", Some("1:N:0:1"), ReadType::R1).unwrap(),
            "x/1"
        );
    }

    // ---------------- revcomp ----------------

    /// -r 大小写行为 oracle 实测: 大小写保留的互补（ACGTacgt → acgtACGT; aCcG → CgGt）。
    #[test]
    fn revcomp_preserves_case_oracle_locked() {
        assert_eq!(revcomp_seq(b"ACGTacgt"), b"acgtACGT");
        assert_eq!(revcomp_seq(b"aCcG"), b"CgGt");
        assert_eq!(revcomp_seq(b"TTTTgggg"), b"ccccAAAA");
        assert_eq!(revcomp_seq(b"CCCC"), b"GGGG");
        assert_eq!(revcomp_seq(b""), b"");
        assert_eq!(revcomp_seq(b"A"), b"T"); // 奇数长度的中间碱基也要互补
    }

    /// comp_tab 全表行为（含 IUPAC 与 U→A）。oracle 实测:
    /// RYKMSWBDHVNrykmswbdhvnUuXx → xXaAnbdhvwskmryNBDHVWSKMRY。
    #[test]
    fn revcomp_iupac_table_oracle_locked() {
        assert_eq!(
            revcomp_seq(b"RYKMSWBDHVNrykmswbdhvnUuXx"),
            b"xXaAnbdhvwskmryNBDHVWSKMRY"
        );
        assert_eq!(revcomp_seq(b"NnRrYyUu"), b"aArRyYnN");
        // 非 ACGT 字符按 comp_tab 恒等（'E'→'E'、'-'→'-'）
        assert_eq!(revcomp_seq(b"E-Z"), b"Z-E");
    }

    // ---------------- fq_records_to_fa ----------------

    /// 冒烟输入 = oracle 冒烟用的 /tmp/p5_clean.fq（4 条、4 种 header 形态）。
    fn p5_clean() -> Vec<u8> {
        b"@rd1 1:N:0:1\nACGTacgt\n+\nIIIIIIII\n\
          @rd2/1\nTTTTgggg\n+\nIIIIIIII\n\
          @rd3_forward\naCcG\n+\nIIII\n\
          @rd4\nCCCC\n+\nIIII\n"
            .to_vec()
    }

    /// 黄金对拍: 输出逐字节 == fixtures/p1/seqtk_names.fa（`seq -A -R 1` 生成）。
    /// 大小写原样保留（无 -U，无质量掩蔽——qual_thres 0 < qual_shift 33）。
    #[test]
    fn fq_to_fa_golden_fixture_r1() {
        let golden = std::fs::read(format!(
            "{}/../../fixtures/p1/seqtk_names.fa",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        assert_eq!(
            fq_records_to_fa(&p5_clean(), ReadType::R1, false).unwrap(),
            golden
        );
    }

    /// `seq -A -R 1 -r` oracle 实测（大小写保留的反向互补）。
    #[test]
    fn fq_to_fa_revcomp_oracle_locked() {
        let out =
            String::from_utf8(fq_records_to_fa(&p5_clean(), ReadType::R1, true).unwrap()).unwrap();
        assert_eq!(
            out,
            ">rd1/1\nacgtACGT\n>rd2/1\nccccAAAA\n>rd3/1\nCgGt\n>rd4/1\nGGGG\n"
        );
    }

    /// 类型不匹配（rd5 2:N:0:1 配 -R 1）→ Err（原版 exit(2)，rd5 之前无错条目已打印——
    /// 本函数整体 Err，见文档说明的流式差异）。
    #[test]
    fn fq_to_fa_type_mismatch_returns_err() {
        let mut input = p5_clean();
        input.extend_from_slice(b"@rd5 2:N:0:1\nAAAA\n+\nIIII\n");
        let e = fq_records_to_fa(&input, ReadType::R1, false).unwrap_err();
        assert_eq!(
            e.to_string(),
            "fastq format error: Error, found read_type 2 but expecting read_type 1"
        );
    }

    /// 空输入/仅空行 → 一条记录都没有 → Err（原版 exit(5) "no records were correctly parsed"）。
    #[test]
    fn fq_to_fa_no_records_is_error() {
        assert!(fq_records_to_fa(b"", ReadType::R1, false).is_err());
        assert!(fq_records_to_fa(b"\n\n", ReadType::R1, false).is_err());
    }

    /// 解析错误（4 行严格校验）向上传播。
    #[test]
    fn fq_to_fa_malformed_propagates_err() {
        // header 未以 '@' 开始
        assert!(fq_records_to_fa(b"rd1\nACGT\n+\nIIII\n", ReadType::R1, false).is_err());
        // seq/qual 长度不一致
        assert!(fq_records_to_fa(b"@rd1\nACGTA\n+\nIIII\n", ReadType::R1, false).is_err());
    }

    /// kseq 的 name/comment 分隔是任意行内空白（空格或制表符）。
    #[test]
    fn fq_to_fa_tab_separated_comment() {
        let out =
            fq_records_to_fa(b"@tabN\t1:N:0:1\nACGT\n+\nIIII\n", ReadType::R1, false).unwrap();
        assert_eq!(out, b">tabN/1\nACGT\n");
    }

    /// kseq isspace 全集: '\v'(0x0B)、'\f'(0x0C)、'\r' 也切 name/comment
    /// （'\n' 不会出现在行内）。修 T5 审查 Issue 1: 此前只认空格与制表符。
    #[test]
    fn fq_to_fa_all_isspace_separators() {
        // comment "1:N:0:1" 被 \v 分隔 → 规则 3 触发
        assert_eq!(
            fq_records_to_fa(b"@vN\x0b1:N:0:1\nACGT\n+\nIIII\n", ReadType::R1, false).unwrap(),
            b">vN/1\nACGT\n"
        );
        // \f: comment "x" 长度 1 不触发规则 3 → 兜底
        assert_eq!(
            fq_records_to_fa(b"@fN\x0cx\nACGT\n+\nIIII\n", ReadType::R1, false).unwrap(),
            b">fN/1\nACGT\n"
        );
        // \r 作为分隔符（CRLF 残留场景）
        assert_eq!(
            fq_records_to_fa(b"@rN\r1:N:0:1\nACGT\n+\nIIII\n", ReadType::R1, false).unwrap(),
            b">rN/1\nACGT\n"
        );
        // 类型不匹配经 \v 分隔的 comment 照样报错
        assert!(fq_records_to_fa(b"@vN\x0b2:N:0:1\nACGT\n+\nIIII\n", ReadType::R1, false).is_err());
    }

    // ---------------- core_read_name ----------------

    /// Fastq_reader.pm L113-116 + insilico_read_normalization.pl L534-535 的 perl 实测
    /// （逐值对拍，含与 seqtk 截断语义的分叉点）。
    #[test]
    fn core_read_name_perl_semantics() {
        assert_eq!(core_read_name("rd/1"), "rd");
        assert_eq!(core_read_name("rd/2"), "rd");
        assert_eq!(core_read_name("rd"), "rd");
        assert_eq!(core_read_name("rd_forward"), "rd");
        assert_eq!(core_read_name("rd_reverse"), "rd");
        // 先剥尾部 /1 再删 _forward
        assert_eq!(core_read_name("rd_forward/1"), "rd");
        // 删子串后**拼接**（与 seqtk 截断到 "x" 不同!）
        assert_eq!(core_read_name("x_forwardandmore_reverse"), "xandmore");
        assert_eq!(core_read_name("a_forwardb_reversec"), "abc");
        // 尾部 /1 剥离去要求前面至少 1 字符: "/1" 原样
        assert_eq!(core_read_name("/1"), "/1");
        // 只剥 [12]: "/3"、"/12" 不剥
        assert_eq!(core_read_name("rd/3"), "rd/3");
        assert_eq!(core_read_name("rd/12"), "rd/12");
        // 贪婪回溯: 尾巴按最后一节剥（perl ^(\S+)/([12])$）
        assert_eq!(core_read_name("rd/1/2"), "rd/1");
    }
}
