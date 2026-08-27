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
    /// `--flow-cov-fusion`（M8）: 无 read 证据边以 junction k-mer 覆盖度兜底。
    pub flow_cov_fusion: bool,
    /// `--flow-cov-cap`（默认 10）。
    pub flow_cov_cap: u64,
    /// 计数表路径（mer_counts.<k>.<SR>.rs.bin, 收尾清理会删——
    /// 缺失时回退 count_source_fa 重数）。
    pub counts_bin: PathBuf,
    /// 计数源 FASTA（both.fa / single.fa, 重数回退用）。
    pub count_source_fa: PathBuf,
    /// 计数表 DS 模式（load_counts_bin 校验用）。
    pub counts_ds: bool,
}

/// M8 覆盖度表装载: bin → 源 FASTA 重数 → 双缺失优雅降级（None = 现状行为）。
/// 一次性装载, `Arc` 跨组件线程只读共享。
fn load_flow_cov(
    params: &ButterflyPoolParams,
) -> Option<(std::sync::Arc<trinity_kmer::counter::CountMap>, u64)> {
    use trinity_kmer::counter::{CountMap, KmerCountTable};
    let cap = params.flow_cov_cap;
    if params.counts_bin.is_file() {
        match trinity_kmer::counter::load_counts_bin(
            &params.counts_bin,
            params.kmer_size,
            params.counts_ds,
        ) {
            Ok(m) => {
                eprintln!(
                    "[flow-cov] 装载计数表 {} 键（{}）",
                    m.len(),
                    params.counts_bin.display()
                );
                return Some((std::sync::Arc::new(m), cap));
            }
            Err(e) => eprintln!("[flow-cov] 计数表装载失败（{e}）, 尝试源 FASTA 重数"),
        }
    }
    if params.count_source_fa.is_file() {
        let t0 = std::time::Instant::now();
        let data = fs::read(&params.count_source_fa)
            .unwrap_or_else(|e| panic!("read {}: {e}", params.count_source_fa.display()));
        let m: CountMap = KmerCountTable::count_fasta_data_with_capacity(
            &data,
            params.kmer_size,
            params.counts_ds,
            1 << 20,
        );
        eprintln!(
            "[flow-cov] 源 FASTA 重数 {} 键（{:.1?}）",
            m.len(),
            t0.elapsed()
        );
        return Some((std::sync::Arc::new(m), cap));
    }
    eprintln!("[flow-cov] 计数表与计数源均缺失, 融合禁用（行为同 floor 1）");
    None
}

/// 原版 butterfly 命令实参（L2363-2388）:
/// `-N 100000 -L <min_contig_length> -F <group_pairs_distance> -C <base>.graph
/// --path_reinforcement_distance=25 --NO_EM_REDUCE`。
fn component_params(
    base: &Path,
    p: &ButterflyPoolParams,
    flow_cov: &Option<(std::sync::Arc<trinity_kmer::counter::CountMap>, u64)>,
) -> ComponentParams {
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
        flow_cov: flow_cov.clone(),
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
    // M8: fusion 开启才装载（关闭零成本零行为差）
    let flow_cov = if params.flow_cov_fusion {
        load_flow_cov(params)
    } else {
        None
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(params.cpu.max(1))
        .build()
        .map_err(|e| CommonError::Parse(format!("cannot build butterfly thread pool: {e}")))?;

    pool.install(|| {
        listing.par_iter().for_each(|(id, base)| {
            if let Err(e) = run_one(base, params, &flow_cov) {
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

fn run_one(
    base: &Path,
    params: &ButterflyPoolParams,
    flow_cov: &Option<(std::sync::Arc<trinity_kmer::counter::CountMap>, u64)>,
) -> Result<(), CommonError> {
    let graph_out = PathBuf::from(format!("{}.graph.out", base.display()));
    let graph_reads = PathBuf::from(format!("{}.graph.reads", base.display()));
    let out = PathBuf::from(format!("{}.graph.allProbPaths.fasta", base.display()));
    let g = fs::read_to_string(&graph_out)
        .map_err(|e| CommonError::Parse(format!("cannot read {}: {e}", graph_out.display())))?;
    let r = fs::read_to_string(&graph_reads)
        .map_err(|e| CommonError::Parse(format!("cannot read {}: {e}", graph_reads.display())))?;
    match params.path_engine {
        PathEngine::Classic => {
            let result = run_component(&g, &r, &component_params(base, params, flow_cov))?;
            fs::write(&out, &result.all_prob_paths_fasta)
                .map_err(|e| CommonError::Parse(format!("cannot write {}: {e}", out.display())))?;
        }
        PathEngine::PogFlow => {
            // M2.1: 走 run_component（POG 递归需要 butterfly 的大栈线程），
            // 引擎参数切 PogFlow——create_dag_from_overlap_layout 后流分解。
            let mut cp = component_params(base, params, flow_cov);
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

#[cfg(test)]
mod tests {
    use super::*;

    fn tmpdir(tag: &str) -> PathBuf {
        let d =
            std::env::temp_dir().join(format!("trinity-rust-flowcov-{tag}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    fn params(dir: &Path, fusion: bool) -> ButterflyPoolParams {
        ButterflyPoolParams {
            cpu: 1,
            min_contig_length: 100,
            group_pairs_distance: 500,
            stack_size_mb: 16,
            path_engine: PathEngine::PogFlow,
            kmer_size: 25,
            flow_cov_fusion: fusion,
            flow_cov_cap: 10,
            counts_bin: dir.join("mer_counts.25.asm.rs.bin"),
            count_source_fa: dir.join("both.fa"),
            counts_ds: true,
        }
    }

    #[test]
    fn flow_cov_missing_everything_degrades_to_none() {
        let dir = tmpdir("none");
        assert!(load_flow_cov(&params(&dir, true)).is_none());
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn flow_cov_loads_bin_when_present() {
        let dir = tmpdir("bin");
        let mut m = trinity_kmer::counter::CountMap::default();
        m.insert(42u64, 7u32);
        trinity_kmer::counter::save_counts_bin(&dir.join("mer_counts.25.asm.rs.bin"), &m, 25, true)
            .unwrap();
        let (table, cap) = load_flow_cov(&params(&dir, true)).unwrap();
        assert_eq!(cap, 10);
        assert_eq!(table.get(&42), Some(&7));
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn flow_cov_recounts_from_source_fa() {
        let dir = tmpdir("fa");
        // 无 bin, 只有源 FASTA: 26 个 A → 25-mer "A"*25 计数 2
        fs::write(
            dir.join("both.fa"),
            ">r1\nAAAAAAAAAAAAAAAAAAAAAAAAAA\n>r2\nAAAAAAAAAAAAAAAAAAAAAAAAAA\n",
        )
        .unwrap();
        let (table, _) = load_flow_cov(&params(&dir, true)).unwrap();
        use trinity_common::kmer::{get_ds_kmer_val, kmer_to_intval};
        let key = get_ds_kmer_val(kmer_to_intval(b"AAAAAAAAAAAAAAAAAAAAAAAAA").unwrap(), 25);
        assert_eq!(table.get(&key), Some(&4)); // 每条 26nt 序列 2 窗 × 2 序列
        let _ = fs::remove_dir_all(&dir);
    }
}
