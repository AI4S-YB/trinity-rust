//! gzip 魔数嗅探读入 — 纯 Rust 解压（flate2/miniz_oxide）。
//! 仅流式入口（Chain/MultiGzDecoder 非 Seek）；.fq.gz 随机访问不在范围内。

use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::error::CommonError;

/// 打开文件，按魔数 (1f 8b) 自动套 gzip 解压层。
pub fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead + Send>, CommonError> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    // 魔数字节已被消费，须前置于剩余文件流（cursor 在前、file 在后），否则顺序错乱。
    let chain = Cursor::new(magic[..n].to_vec()).chain(file);
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(chain))))
    } else {
        Ok(Box::new(BufReader::new(chain)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plain_file_passthrough() {
        let dir = std::env::temp_dir().join("trinity_common_ioutil_t1");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("x.fa");
        std::fs::write(&p, b">a\nACGT\n").unwrap();
        let mut r = open_maybe_gz(&p).unwrap();
        let mut s = String::new();
        r.read_to_string(&mut s).unwrap();
        assert_eq!(s, ">a\nACGT\n");
    }

    #[test]
    fn gz_file_transparent() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;
        let dir = std::env::temp_dir().join("trinity_common_ioutil_t2");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("x.fq.gz");
        let mut enc = GzEncoder::new(File::create(&p).unwrap(), Compression::default());
        enc.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
        enc.finish().unwrap();
        let mut r = open_maybe_gz(&p).unwrap();
        let mut s = String::new();
        r.read_to_string(&mut s).unwrap();
        assert_eq!(s, "@r1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn multi_member_gz_fully_decoded() {
        // 回归锁：将来把 MultiGzDecoder 简化为 GzDecoder 会静默破坏拼接成员
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;
        let dir = std::env::temp_dir().join("trinity_common_ioutil_t3");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("multi.fq.gz");
        let mut buf = Vec::new();
        for chunk in ["@r1\nACGT\n+\nIIII\n", "@r2\nTTTT\n+\n!!!!\n"] {
            let mut enc = GzEncoder::new(Vec::new(), Compression::default());
            enc.write_all(chunk.as_bytes()).unwrap();
            buf.extend_from_slice(&enc.finish().unwrap());
        }
        std::fs::write(&p, &buf).unwrap();
        let mut r = open_maybe_gz(&p).unwrap();
        let mut s = String::new();
        r.read_to_string(&mut s).unwrap();
        assert_eq!(s, "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\n!!!!\n");
    }
}
