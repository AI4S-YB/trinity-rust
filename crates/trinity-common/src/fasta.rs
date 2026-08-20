//! FASTA 读取 — 镜像 Inchworm/src/Fasta_reader.cpp + Fasta_entry.cpp
//! 行为: header 去 '>' 取整行为 _header；首个空白分隔 token 为 _accession（Fasta_entry.cpp:22-30）；
//! 序列大写化并去空白（Fasta_reader.cpp:121-122）。
//! 与原版的已证差异（均为有意为之，勿"修复"）:
//! - 跳过空行；原版读到即用
//! - 首个 '>' 前的非空白垃圾行: 本版报错, 原版静默跳过（FASTQ 误喂时本版更安全）
//! - CRLF: 本版剥 \r, 原版序列中留 \r 字节（下游 nongatc 报错）; is_ascii_whitespace 亦滤 \f
//! - 文件末尾 header 无换行: 本版产出该（空序列）记录, 原版因 eof 标志缺陷丢弃

use std::io::BufRead;

use crate::error::CommonError;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    /// 不含 '>' 的完整 header 行
    pub header: String,
    /// header 首个 token（按空格/制表符切）
    pub accession: String,
    /// 已大写、已去空白的序列
    pub sequence: String,
}

impl FastaRecord {
    /// Fasta_entry 构造逻辑
    pub fn new(header_line: &str, sequence: String) -> Self {
        let header = header_line
            .strip_prefix('>')
            .unwrap_or(header_line)
            .to_string();
        let accession = header
            .split([' ', '\t'])
            .find(|t| !t.is_empty())
            .unwrap_or("")
            .to_string();
        FastaRecord {
            header,
            accession,
            sequence,
        }
    }
}

pub struct FastaReader<R: BufRead> {
    reader: R,
    pending_header: Option<String>,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(reader: R) -> Self {
        FastaReader {
            reader,
            pending_header: None,
        }
    }

    /// 返回下一条记录；EOF 返回 Ok(None)。
    /// 比原版略宽容：跳过空行（原版读到即用，Trinity 管道产物无空行）。
    pub fn next_record(&mut self) -> Result<Option<FastaRecord>, CommonError> {
        let header_line = match self.pending_header.take() {
            Some(h) => h,
            None => match self.read_line_trimmed(true)? {
                Some(l) if l.starts_with('>') => l,
                Some(l) => return Err(CommonError::FastaFormat(format!("记录未以 '>' 开始: {l}"))),
                None => return Ok(None), // EOF
            },
        };

        let mut sequence = String::new();
        loop {
            match self.read_line_trimmed(false)? {
                Some(l) if l.starts_with('>') => {
                    self.pending_header = Some(l);
                    break;
                }
                Some(l) => {
                    for c in l.chars() {
                        if !c.is_ascii_whitespace() {
                            sequence.push(c.to_ascii_uppercase());
                        }
                    }
                }
                None => break, // EOF 结束本记录
            }
        }
        Ok(Some(FastaRecord::new(&header_line, sequence)))
    }

    /// 读一行并去掉 \r\n；skip_empty=true 时跳过空行。
    fn read_line_trimmed(&mut self, skip_empty: bool) -> Result<Option<String>, CommonError> {
        let mut line = String::new();
        loop {
            line.clear();
            let n = self.reader.read_line(&mut line)?;
            if n == 0 {
                return Ok(None);
            }
            let t = line.trim_end_matches(['\n', '\r']);
            if skip_empty && t.trim().is_empty() {
                continue;
            }
            return Ok(Some(t.to_string()));
        }
    }
}

/// IRKE.cpp:386 add_fasta_seq_line_breaks — 每 interval 个字符插入换行，
/// 末尾不追加（调用方负责最后的 '\n'）。interval=0 视为不折行（防御，原版不会传 0）。
pub fn add_fasta_seq_line_breaks(sequence: &[u8], interval: usize) -> String {
    if interval == 0 || sequence.len() <= interval {
        return String::from_utf8_lossy(sequence).into_owned();
    }
    let mut out = Vec::with_capacity(sequence.len() + sequence.len() / interval + 8);
    for (i, &b) in sequence.iter().enumerate() {
        out.push(b);
        if (i + 1) % interval == 0 && i + 1 != sequence.len() {
            out.push(b'\n');
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// 常规 FASTA 记录输出（header 不带 '>' 传入——与 FastaRecord::new 的可选 '>' 宽容相反，本函数不剥）。
pub fn format_fasta_record(header: &str, sequence: &[u8], interval: usize) -> String {
    format!(
        ">{header}\n{}\n",
        add_fasta_seq_line_breaks(sequence, interval)
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    fn read_all(input: &str) -> Vec<FastaRecord> {
        let mut r = FastaReader::new(BufReader::new(input.as_bytes()));
        let mut out = Vec::new();
        while let Some(rec) = r.next_record().unwrap() {
            out.push(rec);
        }
        out
    }

    #[test]
    fn basic_single_record() {
        let recs = read_all(">acc1 some description\nacgt\nACGT\n");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].accession, "acc1");
        assert_eq!(recs[0].header, "acc1 some description");
        assert_eq!(recs[0].sequence, "ACGTACGT"); // 大写化
    }

    #[test]
    fn multiple_records_and_pending_header() {
        let recs = read_all(">a\nAAAA\n>b\ncccc\ngggg\n");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].sequence, "AAAA");
        assert_eq!(recs[1].sequence, "CCCCGGGG"); // 多行拼接 + 大写化
    }

    #[test]
    fn skips_blank_lines_and_crlf() {
        let recs = read_all("\n>a\nAC\r\nGT\n\n>b\nTT\n");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].sequence, "ACGT"); // \r 被剥、空行跳过
        assert_eq!(recs[1].sequence, "TT");
    }

    #[test]
    fn strips_whitespace_in_sequence() {
        let recs = read_all(">a\nAC GT\tAA\n");
        assert_eq!(recs[0].sequence, "ACGTAA"); // 去 ' ' 与 '\t'（原版 remove_whitespace 语义）
    }

    #[test]
    fn error_when_missing_gt() {
        let mut r = FastaReader::new(BufReader::new("notafasta\nACGT\n".as_bytes()));
        assert!(matches!(r.next_record(), Err(CommonError::FastaFormat(_))));
    }

    #[test]
    fn empty_input() {
        assert_eq!(read_all(""), Vec::new());
    }

    #[test]
    fn line_breaks_at_interval() {
        // IRKE.cpp:386: 每 interval 字符换行，最后一组不换
        assert_eq!(
            add_fasta_seq_line_breaks(b"AAAAABBBBBCCCCC", 5),
            "AAAAA\nBBBBB\nCCCCC"
        );
        // 16 字符: 位置 5/10/15 处换行（15 != 16 所以换）
        assert_eq!(
            add_fasta_seq_line_breaks(b"AAAAABBBBBCCCCCD", 5),
            "AAAAA\nBBBBB\nCCCCC\nD"
        );
        // 短于 interval 不折
        assert_eq!(add_fasta_seq_line_breaks(b"ACGT", 60), "ACGT");
        // interval=0 防御: 不折行
        assert_eq!(add_fasta_seq_line_breaks(b"ACGTACGT", 0), "ACGTACGT");
    }

    #[test]
    fn format_record() {
        let out = format_fasta_record("a1;25 total_counts: 30", b"ACGTACGT", 4);
        assert_eq!(out, ">a1;25 total_counts: 30\nACGT\nACGT\n");
    }

    #[test]
    fn accession_skips_leading_whitespace() {
        // 原版 tokenize 跳过前导分隔符（string_util.cpp tokenize）
        let recs = read_all(">  acc1 desc\nACGT\n");
        assert_eq!(recs[0].accession, "acc1");
        assert_eq!(recs[0].header, "  acc1 desc");
        // 全空白 header → accession 仍为空
        let recs2 = read_all("> \t\nACGT\n");
        assert_eq!(recs2[0].accession, "");
    }

    #[test]
    fn header_without_sequence() {
        let recs = read_all(">a\n>b\nTT\n");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].sequence, ""); // 空 seq 记录（与原版实证一致）
    }

    #[test]
    fn iupac_and_digits_passthrough_uppercased() {
        // 只大写不过滤——下游靠 kmer_to_intval 报 nongatc
        let recs = read_all(">a\nacgtnrysw123\n");
        assert_eq!(recs[0].sequence, "ACGTNRYSW123");
    }

    #[test]
    fn trailing_header_without_newline() {
        let recs = read_all(">a\nACGT\n>b");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[1].accession, "b");
        assert_eq!(recs[1].sequence, "");
    }
}
