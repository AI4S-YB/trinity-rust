//! butterfly 组件池（P5-T3）——原版 `PARAFLY -c butterfly_commands -shuffle
//! -CPU $bflyCPU -failed_cmds failed_butterfly_commands.$$.txt`（L1931-1935）
//! 的 rayon 移植。
//!
//! 原版 L1236-1240: `bflyCalculateCPU` 默认 0 → `bflyCPU` 未定义 →
//! `bflyCPU = CPU`——即组件池线程数 = `--CPU`（本实现 `params.cpu`）。
//!
//! 每组件: 从 Cbin 的 `graphs.pack` 按偏移读 graph.out / graph.reads 记录 →
//! [`trinity_butterfly::run_component`]（其内部已在 `stack_size` 配置的
//! 专用线程上跑, 池线程只做 IO + 驱动）→ allProbPaths 追加进 Cbin 的
//! `butterfly.pack`（IO 优化: 替代 per-component 小文件, 见 chrysalis pack 模块）。
//!
//! 隔离（PARAFLY -failed_cmds 语义）: 单组件 Err 不中断其余; 结束时失败
//! 清单非空 → 写 `failed_butterfly_commands.txt` 并整体 Err。

use std::cell::RefCell;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use rayon::prelude::*;
use trinity_butterfly::{
    run_component, run_component_gtasm, run_component_rnaspades, ComponentParams, EdgeSupport,
    GtModel,
};
use trinity_chrysalis::pack::{
    read_at, PackIndex, PackWriter, BUTTERFLY_PACK, GRAPHS_PACK, KIND_ALL_PROB_PATHS,
    KIND_GRAPH_OUT, KIND_GRAPH_READS,
};
use trinity_common::error::CommonError;

use crate::args::PathEngine;

/// gtasm 引擎的池级共享上下文（模型 + 原始 reads 边支持表 + read pair 分支
/// 支持表; 逐 comp 只读）。
#[derive(Debug)]
pub struct GtasmShared {
    pub model: GtModel,
    pub support: EdgeSupport,
    /// read pair 分支支持（--gtasm-pair-support 关时为空表, 解码零开销）。
    pub pair_support: trinity_butterfly::PairSupport,
}

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
    /// 路径引擎（classic 默认 = 现状逐字节; gtasm = GT 打分 + GTasm 解码;
    /// rnaspades = rnaSPAdes MultiExtender 核心规则复刻）。
    pub path_engine: PathEngine,
    /// gtasm 共享上下文（path_engine == Gtasm 时必有）。
    pub gtasm: Option<Arc<GtasmShared>>,
    /// gtasm 解码的 visited 穿越（--no-gtasm-visited-walk 关闭; 默认开,
    /// 机制见 trinity-butterfly gtdecode 模块文档）。仅 gtasm 引擎消费。
    pub gtasm_visited_walk: bool,
    /// gtasm 解码的简单 bulge 双臂产出（--gtasm-bulge-dual 开; 默认关,
    /// 机制见 gtdecode 模块文档）。仅 gtasm 引擎消费。
    pub gtasm_bulge_dual: bool,
    /// gtasm 解码的 read pair 分支证据乘性消费（--gtasm-pair-support 开;
    /// 默认关）。仅 gtasm 引擎消费。
    pub gtasm_pair_support: bool,
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

// graphs.pack 的 per-thread 打开句柄缓存（同一线程连续命中同 Cbin 时省
// open/close; NFS 上每次 open 都是一次元数据往返）。
thread_local! {
    static PACK_FILES: RefCell<HashMap<PathBuf, fs::File>> = RefCell::new(HashMap::new());
}

fn with_pack_file<R>(
    path: &Path,
    f: impl FnOnce(&fs::File) -> Result<R, CommonError>,
) -> Result<R, CommonError> {
    PACK_FILES.with(|m| {
        let mut m = m.borrow_mut();
        if !m.contains_key(path) {
            let file = fs::File::open(path)
                .map_err(|e| CommonError::Parse(format!("cannot read {}: {e}", path.display())))?;
            m.insert(path.to_path_buf(), file);
        }
        f(&m[path])
    })
}

/// 逐组件跑 butterfly。`outdir` 用于落失败清单
/// （原版 `failed_butterfly_commands.$$.txt`; 本实现固定名无 PID）。
pub fn run_butterfly_pool(
    listing: &[(u64, PathBuf)],
    outdir: &Path,
    params: &ButterflyPoolParams,
) -> Result<(), CommonError> {
    // 索引预载（每 Cbin 一次）+ butterfly 包 truncate 重建
    //（本阶段由 .butterfly.ok 整段 checkpoint 保护, 重跑 = 重建, 幂等）。
    let mut indexes: HashMap<PathBuf, PackIndex> = HashMap::new();
    let mut writers: HashMap<PathBuf, Mutex<PackWriter>> = HashMap::new();
    for (_, base) in listing {
        let dir = base.parent().map(|d| d.to_path_buf()).unwrap_or_default();
        if !indexes.contains_key(&dir) {
            indexes.insert(
                dir.clone(),
                PackIndex::load(&dir.join(format!("{GRAPHS_PACK}.idx")))?,
            );
            writers.insert(
                dir.clone(),
                Mutex::new(PackWriter::create(&dir, BUTTERFLY_PACK)?),
            );
        }
    }
    let indexes = indexes;
    let writers = writers;

    let failed: Mutex<Vec<(u64, String)>> = Mutex::new(Vec::new());
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(params.cpu.max(1))
        .build()
        .map_err(|e| CommonError::Parse(format!("cannot build butterfly thread pool: {e}")))?;

    pool.install(|| {
        listing.par_iter().for_each(|(id, base)| {
            let dir = base.parent().map(|d| d.to_path_buf()).unwrap_or_default();
            let idx = &indexes[&dir];
            let w = &writers[&dir];
            if let Err(e) = run_one(*id, base, idx, w, params) {
                eprintln!("-butterfly component {id} failed: {e}");
                failed
                    .lock()
                    .unwrap()
                    .push((*id, format!("{}: {e}", base.display())));
            }
        });
    });
    for w in writers.values() {
        w.lock().unwrap().finish()?;
    }

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
    id: u64,
    base: &Path,
    idx: &PackIndex,
    writer: &Mutex<PackWriter>,
    params: &ButterflyPoolParams,
) -> Result<(), CommonError> {
    let dir = base.parent().map(|d| d.to_path_buf()).unwrap_or_default();
    let pack = dir.join(format!("{GRAPHS_PACK}.pack"));
    let read_rec = |kind: u64| -> Result<String, CommonError> {
        let (offset, len) = idx.find(id, kind).ok_or_else(|| {
            CommonError::Parse(format!(
                "graphs pack {} has no record (id={id}, kind={kind})",
                pack.display()
            ))
        })?;
        let bytes = with_pack_file(&pack, |f| read_at(f, offset, len))?;
        String::from_utf8(bytes)
            .map_err(|e| CommonError::Parse(format!("non-UTF8 record in {}: {e}", pack.display())))
    };
    let fasta = match params.path_engine {
        PathEngine::Classic => {
            let g = read_rec(KIND_GRAPH_OUT)?;
            let r = read_rec(KIND_GRAPH_READS)?;
            run_component(&g, &r, &component_params(base, params))?.all_prob_paths_fasta
        }
        PathEngine::Gtasm => {
            // gtasm 只需 graph.out（支持来自原始 reads 支持表, 不读 graph.reads）
            let g = read_rec(KIND_GRAPH_OUT)?;
            let shared = params.gtasm.as_ref().ok_or_else(|| {
                CommonError::Parse("gtasm 模式缺少共享上下文（模型/支持表）".to_string())
            })?;
            let name = base
                .file_name()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| "component".to_string());
            run_component_gtasm(
                &g,
                id,
                &name,
                &shared.model,
                &shared.support,
                &shared.pair_support,
                &trinity_butterfly::GtasmParams {
                    min_len: params.min_contig_length,
                    seed: trinity_butterfly::GTASM_DEFAULT_SEED,
                    allow_visited_walk: params.gtasm_visited_walk,
                    bulge_dual: params.gtasm_bulge_dual,
                    pair_boost: params.gtasm_pair_support,
                },
            )
            .map_err(|e| CommonError::Parse(format!("gtasm 组件 {id} 失败: {e}")))?
        }
        PathEngine::RnaSpades => {
            // rnaspades 只需 graph.out（简化点: 不做 PE 权重, 不读 graph.reads;
            // 长度过滤走自身 path_cleaning 110bp/cov≥2, harvest 再按
            // --min_contig_length 收口）
            let g = read_rec(KIND_GRAPH_OUT)?;
            let name = base
                .file_name()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| "component".to_string());
            run_component_rnaspades(&g, &name)
                .map_err(|e| CommonError::Parse(format!("rnaspades 组件 {id} 失败: {e}")))?
        }
    };
    writer
        .lock()
        .unwrap()
        .append(id, KIND_ALL_PROB_PATHS, fasta.as_bytes())?;
    Ok(())
}
