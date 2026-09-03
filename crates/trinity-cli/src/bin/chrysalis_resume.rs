//! 一次性实验工具: 对既有 run 目录从 RTT 起重跑 chrysalis 后半段 + butterfly。
//!
//! 用法: chrysalis_resume <iter_dir> <left.fa> <right.fa>
//!         [--multi-kmer] [--no-path-merging] [--cpu N] [--out <fasta>] [--suffix S]
//!
//! <iter_dir> = 含 chrysalis/inchworm.DS.fa.min100 的轮目录。新 chrysalis
//! 产物写到 <iter_dir>/chrysalis_S/（S 默认 rttmk; 不碰原 chrysalis/）。
//! reads 输入先左后右拼接（= both.fa 语义）。

use std::path::{Path, PathBuf};

use trinity_cli::args::PathEngine;
use trinity_cli::butterfly_pool::{run_butterfly_pool, ButterflyPoolParams};
use trinity_cli::harvest::harvest;
use trinity_chrysalis::pipeline::{run_chrysalis_pipeline, ChrysalisParams};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 3 {
        eprintln!(
            "usage: chrysalis_resume <iter_dir> <left.fa> <right.fa> \
             [--multi-kmer] [--no-path-merging] [--cpu N] [--out <fasta>] [--suffix S]"
        );
        std::process::exit(2);
    }
    let iter_dir = PathBuf::from(args.remove(0));
    let left = PathBuf::from(args.remove(0));
    let right = PathBuf::from(args.remove(0));
    let (mut multi_kmer, mut no_path_merging, mut cpu) = (false, false, 32usize);
    let mut out = iter_dir.join("chrysalis_resume.Trinity.fasta");
    let mut suffix = "rttmk".to_string();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--multi-kmer" => multi_kmer = true,
            "--no-path-merging" => no_path_merging = true,
            "--cpu" => {
                i += 1;
                cpu = args[i].parse()?;
            }
            "--out" => {
                i += 1;
                out = PathBuf::from(&args[i]);
            }
            "--suffix" => {
                i += 1;
                suffix = args[i].clone();
            }
            f => return Err(format!("unknown flag {f}").into()),
        }
        i += 1;
    }

    // reads = left ++ right（both.fa 语义, 先左后右）
    let mut reads_data = std::fs::read(&left)?;
    reads_data.extend_from_slice(&std::fs::read(&right)?);
    eprintln!(
        "reads: {} + {} -> {} bytes",
        left.display(),
        right.display(),
        reads_data.len()
    );
    let iworm_data = std::fs::read(iter_dir.join("chrysalis/inchworm.DS.fa.min100"))?;

    let chry_dir = iter_dir.join(format!("chrysalis_{suffix}"));
    let listing = run_chrysalis_pipeline(
        &iworm_data,
        &reads_data,
        &chry_dir,
        &ChrysalisParams {
            strand: false,
            min_contig_length: 200,
            pct_required: 50,
            graphs_per_partition: 1000,
            max_reads: 200_000,
            threads: cpu,
            rtt_multi_kmer: multi_kmer,
            ..Default::default()
        },
    )?;
    eprintln!("chrysalis done: {} components", listing.len());

    run_butterfly_pool(
        &listing,
        &chry_dir,
        &ButterflyPoolParams {
            cpu,
            min_contig_length: 200,
            group_pairs_distance: 500,
            stack_size_mb: 256,
            path_engine: PathEngine::Classic,
            gtasm: None,
            gtasm_visited_walk: false,
            gtasm_bulge_dual: false,
            gtasm_pair_support: false,
            no_path_merging,
            min_read_support: 2,
        },
    )?;
    let n = harvest(&listing, 200, Path::new(&out))?;
    eprintln!("chrysalis_resume: harvested {n} transcripts -> {}", out.display());
    Ok(())
}
