//! butterfly 组件池（P5-T3）——原版 `PARAFLY -c butterfly_commands -shuffle
//! -CPU $bflyCPU -failed_cmds failed_butterfly_commands.$$.txt`（L1931-1935）
//! 的 rayon 移植。
//!
//! 原版 L1236-1240: `bflyCalculateCPU` 默认 0 → `bflyCPU` 未定义 →
//! `bflyCPU = CPU`——即组件池线程数 = `--CPU`（本实现 `params.cpu`）。
//!
//! 每组件: 读 `<base>.graph.out` / `<base>.graph.reads` →
//! [`trinity_butterfly::run_component`]（其内部已在 `stack_size` 配置的
//! 专用线程上跑, 池线程只做 IO + 驱动）→ 写 `<base>.graph.allProbPaths.fasta`。
//!
//! 隔离（PARAFLY -failed_cmds 语义）: 单组件 Err 不中断其余; 结束时失败
//! 清单非空 → 写 `failed_butterfly_commands.txt` 并整体 Err。

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use rayon::prelude::*;
use trinity_butterfly::{run_component, ComponentParams};
use trinity_common::error::CommonError;
use trinity_flow::{run_flow_component, FlowParams};

use crate::args::PathEngine;

/// 组件池参数（从 TrinityArgs 收敛）。
#[derive(Debug, Clone)]
pub struct ButterflyPoolParams {
    /// 线程池大小（原版 bflyCPU = CPU）。
    pub cpu: usize,
    /// `-L`: 最短报告序列长度（--min_contig_length）。
    pub min_contig_length: usize,
    /// `-F`: 配对端最大跨度（--group_pairs_distance）。
    pub group_pairs_distance: usize,
    /// 专用组件线程栈（--bfly_stack_mb）。
    pub stack_size_mb: usize,
    /// 输出层引擎（--path-engine; flow 全程迭代/显式栈——DP 剥离为循环,
    /// Tarjan SCC 用显式栈——不依赖大栈线程, 但仍受 stack_size 保守配置）。
    pub path_engine: PathEngine,
    /// kmer 长（flow 引擎 path_to_seq 用; classic 不消费）。
    pub kmer_size: usize,
}

/// 原版 butterfly 命令实参（L2363-2388）:
/// `-N 100000 -L <min_contig_length> -F <group_pairs_distance> -C <base>.graph
/// --path_reinforcement_distance=25 --NO_EM_REDUCE`。
fn component_params(base: &Path, p: &ButterflyPoolParams) -> ComponentParams {
    ComponentParams {
        n: 100_000,
        min_len: p.min_contig_length,
        max_pair_distance: p.group_pairs_distance,
        no_em_reduce: true,
        // 原版 L97-99/1271-1277: path_reinforcement_distance 默认恒 25
        // （PE/SE 同值, 与 -F 无关——非 25% * F）。
        path_reinforcement_distance: Some(25),
        name: base
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "component".to_string()),
        stack_size: p.stack_size_mb.saturating_mul(1024 * 1024),
        ..Default::default()
    }
}

/// 逐组件跑 butterfly。`outdir` 用于落失败清单
/// （原版 `failed_butterfly_commands.$$.txt`; 本实现固定名无 PID）。
pub fn run_butterfly_pool(
    listing: &[(u64, PathBuf)],
    outdir: &Path,
    params: &ButterflyPoolParams,
) -> Result<(), CommonError> {
    let failed: Mutex<Vec<(u64, String)>> = Mutex::new(Vec::new());
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(params.cpu.max(1))
        .build()
        .map_err(|e| CommonError::Parse(format!("cannot build butterfly thread pool: {e}")))?;

    pool.install(|| {
        listing.par_iter().for_each(|(id, base)| {
            if let Err(e) = run_one(base, params) {
                eprintln!("-butterfly component {id} failed: {e}");
                failed
                    .lock()
                    .unwrap()
                    .push((*id, format!("{}: {e}", base.display())));
            }
        });
    });

    let mut failed = failed.into_inner().unwrap();
    if failed.is_empty() {
        return Ok(());
    }
    failed.sort();
    let path = outdir.join("failed_butterfly_commands.txt");
    let text: String = failed.iter().map(|(_, s)| format!("{s}\n")).collect();
    let _ = fs::write(&path, &text);
    Err(CommonError::Parse(format!(
        "Error, butterfly failed for {} component(s): {} (see {})",
        failed.len(),
        failed.len().min(5),
        path.display()
    )))
}

fn run_one(base: &Path, params: &ButterflyPoolParams) -> Result<(), CommonError> {
    let graph_out = PathBuf::from(format!("{}.graph.out", base.display()));
    let graph_reads = PathBuf::from(format!("{}.graph.reads", base.display()));
    let out = PathBuf::from(format!("{}.graph.allProbPaths.fasta", base.display()));
    let g = fs::read_to_string(&graph_out)
        .map_err(|e| CommonError::Parse(format!("cannot read {}: {e}", graph_out.display())))?;
    let r = fs::read_to_string(&graph_reads)
        .map_err(|e| CommonError::Parse(format!("cannot read {}: {e}", graph_reads.display())))?;
    match params.path_engine {
        PathEngine::Classic => {
            let result = run_component(&g, &r, &component_params(base, params))?;
            fs::write(&out, &result.all_prob_paths_fasta)
                .map_err(|e| CommonError::Parse(format!("cannot write {}: {e}", out.display())))?;
        }
        PathEngine::PogFlow => {
            // M2.1: 走 run_component（POG 递归需要 butterfly 的大栈线程），
            // 引擎参数切 PogFlow——create_dag_from_overlap_layout 后流分解。
            let mut cp = component_params(base, params);
            cp.path_engine = trinity_butterfly::PathEngine::PogFlow;
            let result = run_component(&g, &r, &cp)?;
            fs::write(&out, &result.all_prob_paths_fasta)
                .map_err(|e| CommonError::Parse(format!("cannot write {}: {e}", out.display())))?;
        }
        PathEngine::Flow => {
            // 流分解引擎: 纯 Rust 无递归（不需大栈线程）; 输出同为
            // allProbPaths 格式（header `{comp}_g1_i{n}`, 由 to_fasta 产生）。
            let fasta = run_flow_component(
                &g,
                &r,
                &FlowParams {
                    min_len: params.min_contig_length,
                    min_support: 1,
                    dedup_threshold: 0.95,
                    k: params.kmer_size,
                },
            )?;
            fs::write(&out, &fasta)
                .map_err(|e| CommonError::Parse(format!("cannot write {}: {e}", out.display())))?;
        }
    }
    Ok(())
}
