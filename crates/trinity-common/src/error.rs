//! trinity-common 统一错误类型。消息格式贴近原版（"error, ..." 前缀风格）。

use thiserror::Error;

#[derive(Debug, Error)]
pub enum CommonError {
    /// sequenceUtil.cpp:260-261 "error, kmer length exceeds 32"
    #[error("error, kmer length exceeds 32: {len}")]
    KmerTooLong { len: usize },

    /// sequenceUtil.cpp:276-280 "error, kmer contains nongatc: {kmer}"
    #[error("error, kmer contains nongatc: {kmer}")]
    NonGatcChar { kmer: String },

    /// IRKE 贪心链 throw(string) 的移植载体（IRKE.cpp:852-854 inchworm 轮数超限）。
    /// 消息保留原版原文。
    #[error("error, {0}")]
    Inchworm(String),

    #[error("fasta format error: {0}")]
    FastaFormat(String),

    #[error("fastq format error: {0}")]
    FastqFormat(String),

    /// Chrysalis 侧 `exit(n)` 硬错误的移植载体（消息保留原版 cerr 文案风格）。
    #[error("error, {0}")]
    Parse(String),

    #[error(transparent)]
    Io(#[from] std::io::Error),
}
