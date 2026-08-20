//! trinity-cli 主控骨架（P5-T2）: 原版同名参数面（[`args::parse_args`]）、
//! Pipeliner.pm 语义的 checkpoint 框架（[`checkpoint`]）、prep_seqs 移植
//! （[`prep::prep_reads`]——SS 端 revcomp/多文件拼接/both.fa）。
//!
//! 主线参数集之外的原版选项一律按未知参数报错（trinity-common cli 的
//! `do not understand option` 文案）。

pub mod args;
pub mod butterfly_pool;
pub mod checkpoint;
pub mod harvest;
pub mod orchestrate;
pub mod prep;

pub use args::{parse_args, SeqType, TrinityArgs};
