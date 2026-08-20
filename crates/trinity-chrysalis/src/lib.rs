//! P3 trinity-chrysalis — Chrysalis 六子阶段移植（GraphFromFasta 焊接聚类 /
//! BubbleUpClustering / CreateIwormFastaBundle / ReadsToTranscripts /
//! QuantifyGraph + 分区器）。镜像 `Chrysalis/analysis/*`。
//!
//! 本文件为 T1 基础层：dna_vector（读入/编码/熵/简单性）、kmer_align
//! （12-mer CSR 倒排索引）、nonred_table（排序数组 + 二分计数表）。
//!
//! **编码表警示**：Chrysalis 侧存在两套核苷酸编码——本 crate 的
//! [`dna_vector::nuc_index`]（DNAVector.h plain_table：A=0,C=1,G=2,T=3,N=4）
//! 与 trinity-common::kmer 的 Inchworm 表（G=0,A=1,T=2,C=3）互不相容，
//! 二者分别服务 KmerAlignCore 的 12-mer 桶索引与 k-mer 整数化两条数据流。

pub mod bubble_up;
pub mod bundle;
pub mod dna_vector;
pub mod graph_from_fasta;
pub mod kmer_align;
pub mod nonred_table;
pub mod partition;
pub mod pipeline;
pub mod quantify;
pub mod reads_to_transcripts;
