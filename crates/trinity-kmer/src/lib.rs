//! trinity-kmer — jellyfish 计数替代 + DigiNorm（jellyfish/fastaToKmerCoverageStats/nbkc/seqtk 的镜像移植）。

pub mod counter;
pub mod coverage_stats;
pub mod diginorm;
pub mod drand48;
pub mod dump;
pub mod fmt;
pub mod nbkc;
pub mod read_names;
