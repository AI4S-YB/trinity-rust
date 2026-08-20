//! trinity-chrysalis CLI —— Chrysalis 六子程序形态的统一入口 + chrysalis-all
//! 全链路编排（镜像原版 `Chrysalis/bin/*` 各二进制 + Trinity 管线 shell/perl 串联）。
//!
//! 参数名对齐原版二进制（GraphFromFasta/BubbleUpClustering/…），错误处理沿用
//! 本仓库 bin 约定：参数错 exit 2 / 运行错 exit 1（`trinity_common::cli::CliError`）。
//!
//! `sort-welds` / `sort-rtc` 是管线便利命令（原版用 shell `sort` /
//! `sort -k1,1n -k3,3nr -k2,2`），stdin → stdout。

use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use trinity_chrysalis::bubble_up::{bubble_up_clustering, BubbleParams};
use trinity_chrysalis::bundle::create_iworm_fasta_bundle;
use trinity_chrysalis::dna_vector::{read_fasta, read_fasta_short_names};
use trinity_chrysalis::graph_from_fasta::{graph_from_fasta_files, sort_weld_graph, GffParams};
use trinity_chrysalis::partition::{partition, PartParams};
use trinity_chrysalis::quantify::{quantify_graph, reads_ext_filename, QgParams};
use trinity_chrysalis::reads_to_transcripts::{
    reads_to_transcripts, sort_reads_to_components, RttParams,
};
use trinity_common::cli::{Cli, CliError};
use trinity_common::fasta::FastaReader;
use trinity_inchworm::debruijn::graph_per_record;

const USAGE: &str = "用法: trinity-chrysalis <子命令> [参数]
子命令:
  graph-from-fasta    -i <iworm.fa> -r <reads.fa> [-strand] [-k 24] [-kk 48] [-t N]
                      [-glue_factor 0.05] [-min_glue 2] [-max_glue_required -1]
                      [-min_iso_ratio 0.05] [-no_welds] [-no_glue_required]
                      [-disable_repeat_check] [-report_welds]
  sort-welds          (stdin → stdout)
  bubble-up           -i <iworm.fa> -weld_graph <sorted> [-min_contig_length 200]
                      [-max_cluster_size 25] [-debug_weld_all]
  create-bundle       -i <component.out> -o <out.fa> [-min 200]
  reads-to-transcripts -i <reads.fa> -f <bundle.fa> -o <out> [-t N] [-strand] [-p 50]
                      [-min_kmer_entropy 1.5] [-max_mem_reads N]
  sort-rtc            (stdin → stdout)
  fasta-to-debruijn   --fasta <fa> -K 24 --graph_per_record [--SS] [--threads N]
  partition           --deBruijns <f> --componentReads <sorted> [-N 1000] [-L 200]
                      --outdir <dir>
  quantify-graph      -g <graph.tmp> -i <reads.tmp> -o <graph.out> [-k 24] [-strand]
                      [-max_reads 200000] [-no_cleanup]
  chrysalis-all       -i <iworm.fa> -r <reads.fa> -o <outdir> [--SS] [-L 200] [-p 50]
                      [-N 1000] [--max_reads 200000] [-t N]";

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    std::process::exit(run(&args));
}

fn run(args: &[String]) -> i32 {
    let Some((sub, rest)) = args.split_first() else {
        eprintln!("{USAGE}");
        return 2;
    };
    let mut cli = Cli::new(rest);
    let r = match sub.as_str() {
        "graph-from-fasta" => cmd_graph_from_fasta(&mut cli),
        "sort-welds" => cmd_sort_welds(&mut cli),
        "bubble-up" => cmd_bubble_up(&mut cli),
        "create-bundle" => cmd_create_bundle(&mut cli),
        "reads-to-transcripts" => cmd_reads_to_transcripts(&mut cli),
        "sort-rtc" => cmd_sort_rtc(&mut cli),
        "fasta-to-debruijn" => cmd_fasta_to_debruijn(&mut cli),
        "partition" => cmd_partition(&mut cli),
        "quantify-graph" => cmd_quantify_graph(&mut cli),
        "chrysalis-all" => cmd_chrysalis_all(&mut cli),
        _ => {
            eprintln!("Error, unknown subcommand: {sub}\n{USAGE}");
            return 2;
        }
    };
    match r {
        Ok(()) => 0,
        Err(e) => {
            if e.usage {
                eprintln!("{USAGE}\nError, {}", e.msg);
                2
            } else {
                eprintln!("Error, {}", e.msg);
                1
            }
        }
    }
}

fn read_stdin() -> Result<String, CliError> {
    let mut s = String::new();
    std::io::stdin()
        .read_to_string(&mut s)
        .map_err(|e| CliError::run(format!("cannot read stdin: {e}")))?;
    Ok(s)
}

fn write_stdout(text: &str) -> Result<(), CliError> {
    std::io::stdout()
        .write_all(text.as_bytes())
        .map_err(CliError::run)
}

fn read_path(p: &str) -> Result<String, CliError> {
    std::fs::read_to_string(p).map_err(|e| CliError::run(format!("cannot read {p}: {e}")))
}

fn write_file(p: &Path, text: &str) -> Result<(), CliError> {
    std::fs::write(p, text).map_err(|e| CliError::run(format!("cannot write {}: {e}", p.display())))
}

// ---------------------------------------------------------------------------
// 子命令
// ---------------------------------------------------------------------------

fn cmd_graph_from_fasta(cli: &mut Cli) -> Result<(), CliError> {
    let iworm = PathBuf::from(cli.req_flag("-i")?);
    let reads = PathBuf::from(cli.req_flag("-r")?);
    let p = GffParams {
        k: cli.uint_flag("-k", 24usize)?,
        kk: cli.uint_flag("-kk", 48usize)?,
        strand: cli.bool_flag("-strand"),
        glue_factor: cli.float_flag("-glue_factor", 0.05)?,
        min_glue_required: cli.uint_flag("-min_glue", 2u32)?,
        max_glue_required: cli.int_flag("-max_glue_required", -1i64)?,
        min_iso_ratio: cli.float_flag("-min_iso_ratio", 0.05)?,
        no_welds: cli.bool_flag("-no_welds"),
        no_glue_required: cli.bool_flag("-no_glue_required"),
        disable_repeat_check: cli.bool_flag("-disable_repeat_check"),
        report_welds: cli.bool_flag("-report_welds"),
        debug: false,
        threads: cli.uint_flag("-t", 1usize)?.max(1),
    };
    cli.finish()?;
    write_stdout(&graph_from_fasta_files(&iworm, &reads, &p)?)
}

fn cmd_sort_welds(cli: &mut Cli) -> Result<(), CliError> {
    cli.finish()?;
    write_stdout(&sort_weld_graph(&read_stdin()?))
}

fn cmd_bubble_up(cli: &mut Cli) -> Result<(), CliError> {
    let iworm = PathBuf::from(cli.req_flag("-i")?);
    let weld_graph = cli.req_flag("-weld_graph")?;
    let p = BubbleParams {
        min_contig_length: cli.uint_flag("-min_contig_length", 24usize)?,
        max_cluster_size: cli.uint_flag("-max_cluster_size", 25usize)?,
        debug_weld_all: cli.bool_flag("-debug_weld_all"),
    };
    cli.finish()?;
    let seqs = read_fasta(&iworm)?;
    let graph = read_path(&weld_graph)?;
    write_stdout(&bubble_up_clustering(&seqs, &graph, &p)?)
}

fn cmd_create_bundle(cli: &mut Cli) -> Result<(), CliError> {
    let input = cli.req_flag("-i")?;
    let out = PathBuf::from(cli.req_flag("-o")?);
    let min = cli.uint_flag("-min", 200usize)?;
    cli.finish()?;
    let text = create_iworm_fasta_bundle(&read_path(&input)?, min)?;
    write_file(&out, &text)
}

fn cmd_reads_to_transcripts(cli: &mut Cli) -> Result<(), CliError> {
    let reads = PathBuf::from(cli.req_flag("-i")?);
    let bundle = PathBuf::from(cli.req_flag("-f")?);
    let out = PathBuf::from(cli.req_flag("-o")?);
    let p = RttParams {
        strand: cli.bool_flag("-strand"),
        pct_required: cli.uint_flag("-p", 0u32)?,
        min_kmer_entropy: cli.float_flag("-min_kmer_entropy", 1.5f32)?,
        max_mem_reads: cli.uint_flag("-max_mem_reads", usize::MAX)?,
        threads: cli.uint_flag("-t", 1usize)?.max(1),
    };
    cli.finish()?;
    let r = read_fasta_short_names(&reads)?;
    let b = read_fasta(&bundle)?;
    let o = reads_to_transcripts(&r, &b, &p)?;
    write_file(&out, &o.text)?;
    let rcts = PathBuf::from(format!("{}.rcts.out", out.display()));
    write_file(&rcts, &format!("{}\n", o.mapped_count))
}

fn cmd_sort_rtc(cli: &mut Cli) -> Result<(), CliError> {
    cli.finish()?;
    write_stdout(&sort_reads_to_components(&read_stdin()?))
}

fn cmd_fasta_to_debruijn(cli: &mut Cli) -> Result<(), CliError> {
    let fasta = PathBuf::from(cli.req_flag("--fasta")?);
    let kmer_length = cli.uint_flag("-K", 24usize)?;
    let per_record = cli.bool_flag("--graph_per_record");
    let ss = cli.bool_flag("--SS");
    let _ = cli.uint_flag("--threads", 1u32)?;
    cli.finish()?;
    if !per_record {
        // 原版默认 --graph_per_record；其他模式（单全局图）非本管线所需。
        return Err(CliError::usage("--graph_per_record required"));
    }
    let f = std::fs::File::open(&fasta)
        .map_err(|e| CliError::run(format!("cannot read {}: {e}", fasta.display())))?;
    let mut reader = FastaReader::new(std::io::BufReader::new(f));
    let mut recs = Vec::new();
    while let Some(r) = reader.next_record().map_err(CliError::run)? {
        recs.push(r);
    }
    write_stdout(&graph_per_record(&recs, kmer_length, ss)?)
}

fn cmd_partition(cli: &mut Cli) -> Result<(), CliError> {
    let debruijns = cli.req_flag("--deBruijns")?;
    let comp_reads = cli.req_flag("--componentReads")?;
    let outdir = PathBuf::from(cli.req_flag("--outdir")?);
    let p = PartParams {
        graphs_per_partition: cli.uint_flag("-N", 1000usize)?,
        min_contig_length: cli.uint_flag("-L", 200usize)?,
    };
    cli.finish()?;
    let listing = partition(
        &read_path(&debruijns)?,
        &read_path(&comp_reads)?,
        &outdir,
        &p,
    )?;
    for (id, base) in &listing {
        println!("{id}\t{}", base.display());
    }
    Ok(())
}

fn cmd_quantify_graph(cli: &mut Cli) -> Result<(), CliError> {
    let graph_tmp = cli.req_flag("-g")?;
    let reads_tmp = cli.req_flag("-i")?;
    let out_flag = cli.req_flag("-o")?;
    let p = QgParams {
        k: cli.uint_flag("-k", 24usize)?,
        strand: cli.bool_flag("-strand"),
        max_reads: cli.int_flag("-max_reads", 200_000i64)?,
    };
    let no_cleanup = cli.bool_flag("-no_cleanup");
    cli.finish()?;

    // 库层从输入路径推导输出名；CLI -o 优先——跑完后改名对齐（含 .reads 派生）。
    let o = quantify_graph(&graph_tmp, &reads_tmp, &p)?;
    let want_graph = PathBuf::from(&out_flag);
    let derived_graph = PathBuf::from(&o.graph_out);
    let derived_reads = PathBuf::from(&o.reads_out);
    if want_graph != derived_graph {
        let want_reads = PathBuf::from(reads_ext_filename(&out_flag));
        std::fs::rename(&derived_graph, &want_graph)
            .map_err(|e| CliError::run(format!("cannot rename output: {e}")))?;
        std::fs::rename(&derived_reads, &want_reads)
            .map_err(|e| CliError::run(format!("cannot rename output: {e}")))?;
    }
    // 原版默认成功后删输入（QuantifyGraph.cc:489-493 unlink）
    if !no_cleanup {
        let _ = std::fs::remove_file(&graph_tmp);
        let _ = std::fs::remove_file(&reads_tmp);
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// chrysalis-all —— 全链路编排（输出布局镜像原版 Chrysalis 目录）
// ---------------------------------------------------------------------------

fn cmd_chrysalis_all(cli: &mut Cli) -> Result<(), CliError> {
    let iworm = PathBuf::from(cli.req_flag("-i")?);
    let reads = PathBuf::from(cli.req_flag("-r")?);
    let outdir = PathBuf::from(cli.req_flag("-o")?);
    let ss = cli.bool_flag("--SS");
    let min_len = cli.uint_flag("-L", 200usize)?;
    let pct = cli.uint_flag("-p", 50u32)?;
    let graphs_per = cli.uint_flag("-N", 1000usize)?;
    let max_reads = cli.int_flag("--max_reads", 200_000i64)?;
    let threads = cli.uint_flag("-t", 1u32)?;
    cli.finish()?;

    // 全链编排上移到库层 run_chrysalis_pipeline（CLI/测试同入口）
    let iworm_fa = std::fs::read(&iworm)
        .map_err(|e| CliError::run(format!("cannot read {}: {e}", iworm.display())))?;
    let reads_fa = std::fs::read(&reads)
        .map_err(|e| CliError::run(format!("cannot read {}: {e}", reads.display())))?;
    trinity_chrysalis::pipeline::run_chrysalis_pipeline(
        &iworm_fa,
        &reads_fa,
        &outdir,
        &trinity_chrysalis::pipeline::ChrysalisParams {
            strand: ss,
            min_contig_length: min_len,
            pct_required: pct,
            graphs_per_partition: graphs_per,
            max_reads: max_reads.max(0) as u64,
            threads: threads as usize,
        },
    )?;
    Ok(())
}
