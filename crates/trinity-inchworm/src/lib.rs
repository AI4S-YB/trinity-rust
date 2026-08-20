//! P2 trinity-inchworm — IRKE 贪心组装移植（kmer 目录、剪枝、种子、贪心延伸）。
//! 镜像 Inchworm/src/{IRKE.cpp, KmerCounter.cpp, IRKE_run.cpp}。

pub mod counter_sync;
pub mod debruijn;
pub mod glibc_rand;
pub mod irke;
pub mod kmer_counter;
pub mod visitor;
