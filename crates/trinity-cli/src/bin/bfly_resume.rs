//! 一次性实验工具: 对既有 run 目录（如 iter_k63）只重跑 butterfly 池 + harvest。
//!
//! 用法: bfly-resume <iter_dir> [--no-path-merging] [--cpu N] [--out <fasta>]
//!   <iter_dir> = 含 chrysalis/Component_bins/ 的轮目录（graphs.pack 必在,
//!   butterfly.pack 会被 truncate 重建——注意先备份原产物）。
//! 默认输出 <iter_dir>/bfly_resume.Trinity.fasta。

use std::path::{Path, PathBuf};

use trinity_cli::args::PathEngine;
use trinity_cli::butterfly_pool::{run_butterfly_pool, ButterflyPoolParams};
use trinity_cli::harvest::harvest;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        eprintln!("usage: bfly-resume <iter_dir> [--no-path-merging] [--cpu N] [--out <fasta>]");
        std::process::exit(2);
    }
    let iter_dir = PathBuf::from(args.remove(0));
    let mut no_path_merging = false;
    let mut cpu = 32usize;
    let mut min_read_support = 2i64;
    let mut min_contig_length = 200usize;
    let mut out = iter_dir.join("bfly_resume.Trinity.fasta");
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--no-path-merging" => no_path_merging = true,
            "--min-read-support" | "-R" => {
                i += 1;
                min_read_support = args[i].parse()?;
            }
            "--min-contig-length" | "-L" => {
                i += 1;
                min_contig_length = args[i].parse()?;
            }
            "--cpu" => {
                i += 1;
                cpu = args[i].parse()?;
            }
            "--out" => {
                i += 1;
                out = PathBuf::from(&args[i]);
            }
            f => return Err(format!("unknown flag {f}").into()),
        }
        i += 1;
    }

    // listing 从 component_base_listing.txt（id 升序）——兼容 <iter>/chrysalis/
    // 与直接传 chrysalis 目录（如 chrysalis_rttmk）两种形态
    let listing_path = {
        let a = iter_dir.join("chrysalis/Component_bins/component_base_listing.txt");
        if a.exists() {
            a
        } else {
            iter_dir.join("Component_bins/component_base_listing.txt")
        }
    };
    let mut listing: Vec<(u64, PathBuf)> = Vec::new();
    for line in std::fs::read_to_string(&listing_path)?.lines() {
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() >= 2 {
            listing.push((toks[0].parse()?, PathBuf::from(toks[1])));
        }
    }
    eprintln!(
        "bfly-resume: {} components from {}",
        listing.len(),
        listing_path.display()
    );

    run_butterfly_pool(
        &listing,
        &iter_dir,
        &ButterflyPoolParams {
            cpu,
            min_contig_length,
            group_pairs_distance: 500,
            stack_size_mb: 256,
            path_engine: PathEngine::Classic,
            gtasm: None,
            gtasm_visited_walk: false,
            gtasm_bulge_dual: false,
            gtasm_pair_support: false,
            no_path_merging,
            min_read_support,
        },
    )?;
    let n = harvest(&listing, min_contig_length, Path::new(&out))?;
    eprintln!("bfly-resume: harvested {n} transcripts -> {}", out.display());
    Ok(())
}
