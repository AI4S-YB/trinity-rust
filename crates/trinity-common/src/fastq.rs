//! FASTQ 读取 — 严格 4 行记录（@/seq/+/qual）。
//! 原版由 Perl 脚本转换，此为 Rust 版统一入口；校验比原版严格（报错带记录号与 header）。

use std::io::BufRead;

use crate::error::CommonError;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqRecord {
    pub accession: String, // header 去 '@' 后首个 token
    pub header: String,    // 去 '@' 的完整 header
    pub sequence: String,  // 原样保留大小写（不强制大写，与原版 fq->fa 转换一致）
    pub plus: String,      // '+' 行内容（可为空或重复 header）
    pub qual: String,
}

pub struct FastqReader<R: BufRead> {
    reader: R,
    record_no: u64,
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        FastqReader {
            reader,
            record_no: 0,
        }
    }

    /// 读一条记录；容忍文件级结尾多余空行；EOF 返回 Ok(None)。
    /// 出错后流位置未定义——调用方收到 Err 应当中止读取，不要尝试继续（错位流会连续报错后静默重同步）。
    pub fn next_record(&mut self) -> Result<Option<FastqRecord>, CommonError> {
        let header_line = match self.next_nonempty_line(true)? {
            Some(l) => l,
            None => return Ok(None),
        };
        self.record_no += 1;
        let where_ = format!("记录 #{} ({})", self.record_no, header_line);

        if !header_line.starts_with('@') {
            return Err(CommonError::FastqFormat(format!(
                "header 未以 '@' 开始: {where_}"
            )));
        }
        let sequence = self
            .next_nonempty_line(false)?
            .ok_or_else(|| CommonError::FastqFormat(format!("缺失序列行: {where_}")))?;
        let plus_line = self
            .next_nonempty_line(false)?
            .ok_or_else(|| CommonError::FastqFormat(format!("缺失 '+' 行: {where_}")))?;
        if !plus_line.starts_with('+') {
            return Err(CommonError::FastqFormat(format!(
                "第三行未以 '+' 开始: {where_}"
            )));
        }
        let qual = self
            .next_nonempty_line(false)?
            .ok_or_else(|| CommonError::FastqFormat(format!("缺失质量行: {where_}")))?;
        if qual.len() != sequence.len() {
            return Err(CommonError::FastqFormat(format!(
                "序列({})与质量({})长度不一致: {where_}",
                sequence.len(),
                qual.len()
            )));
        }

        let header = header_line[1..].to_string();
        let accession = header
            .split([' ', '\t'])
            .find(|t| !t.is_empty())
            .unwrap_or("")
            .to_string();
        Ok(Some(FastqRecord {
            accession,
            header,
            sequence,
            plus: plus_line[1..].to_string(),
            qual,
        }))
    }

    /// 读下一非空行。`at_record_start` 为 true（header 位）时跳过空行——仅记录边界允许空行，
    /// 文件级结尾多余换行由此容忍；为 false（seq/+/qual 位）时空行是损坏信号，直接报错。
    fn next_nonempty_line(&mut self, at_record_start: bool) -> Result<Option<String>, CommonError> {
        let mut line = String::new();
        loop {
            line.clear();
            let n = self.reader.read_line(&mut line)?;
            if n == 0 {
                return Ok(None);
            }
            let t = line.trim_end_matches(['\n', '\r']);
            if t.is_empty() {
                if at_record_start {
                    continue; // 记录边界空行（含文件级结尾多余换行）
                }
                return Err(CommonError::FastqFormat(format!(
                    "记录内出现空行: 记录 #{}",
                    self.record_no
                )));
            }
            return Ok(Some(t.to_string()));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    fn read_all(input: &str) -> Vec<FastqRecord> {
        let mut r = FastqReader::new(BufReader::new(input.as_bytes()));
        let mut out = Vec::new();
        while let Some(rec) = r.next_record().unwrap() {
            out.push(rec);
        }
        out
    }

    #[test]
    fn basic_two_records() {
        let recs = read_all("@r1 desc\nACGT\n+\nIIII\n@r2\nTTTT\n+r2\n!!!!\n");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].accession, "r1");
        assert_eq!(recs[0].header, "r1 desc");
        assert_eq!(recs[0].sequence, "ACGT");
        assert_eq!(recs[0].plus, "");
        assert_eq!(recs[0].qual, "IIII");
        assert_eq!(recs[1].plus, "r2");
        assert_eq!(recs[1].qual, "!!!!");
    }

    #[test]
    fn tolerates_crlf_and_trailing_newline() {
        let recs = read_all("@r1\r\nACGT\r\n+\r\nIIII\r\n");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].sequence, "ACGT");
    }

    #[test]
    fn missing_at_header_is_error() {
        let mut r = FastqReader::new(BufReader::new("r1\nACGT\n+\nIIII\n".as_bytes()));
        let e = r.next_record().unwrap_err();
        assert!(matches!(e, CommonError::FastqFormat(_)));
    }

    #[test]
    fn seq_qual_length_mismatch_is_error() {
        let mut r = FastqReader::new(BufReader::new("@r1\nACGTA\n+\nIIII\n".as_bytes()));
        let e = r.next_record().unwrap_err();
        assert!(matches!(e, CommonError::FastqFormat(_)));
    }

    #[test]
    fn truncated_record_is_error() {
        let mut r = FastqReader::new(BufReader::new("@r1\nACGT\n".as_bytes()));
        assert!(matches!(r.next_record(), Err(CommonError::FastqFormat(_))));
    }

    #[test]
    fn blank_line_inside_record_is_error() {
        // 空行出现在 seq/+ 或 qual 位时是损坏信号，必须报错而非静默吞掉
        let mut r = FastqReader::new(BufReader::new(
            "@r1\nACGT\n+\n\n@r2x\nAAAA\n+\nIIII\n".as_bytes(),
        ));
        assert!(matches!(r.next_record(), Err(CommonError::FastqFormat(_))));
    }

    #[test]
    fn file_level_trailing_blank_lines_tolerated() {
        let recs = read_all("@r1\nACGT\n+\nIIII\n\n\n");
        assert_eq!(recs.len(), 1);
    }

    #[test]
    fn accession_skips_leading_whitespace() {
        // 与 fasta.rs 的 FastaRecord::new 及原版 tokenize 语义一致
        let recs = read_all("@  r1 desc\nACGT\n+\nIIII\n");
        assert_eq!(recs[0].accession, "r1");
        assert_eq!(recs[0].header, "  r1 desc");
    }

    #[test]
    fn lone_at_header_empty_accession() {
        let recs = read_all("@\nACGT\n+\nIIII\n");
        assert_eq!(recs[0].accession, "");
        assert_eq!(recs[0].sequence, "ACGT");
    }

    #[test]
    fn empty_input_ok() {
        assert_eq!(read_all(""), Vec::new());
    }
}
