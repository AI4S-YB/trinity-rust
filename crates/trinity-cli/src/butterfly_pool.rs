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
use std::sync::Mutex;

use rayon::prelude::*;
use trinity_butterfly::{run_component, ComponentParams};
use trinity_chrysalis::pack::{
    read_at, PackIndex, PackWriter, BUTTERFLY_PACK, GRAPHS_PACK, KIND_ALL_PROB_PATHS,
    KIND_GRAPH_OUT, KIND_GRAPH_READS,
};
use trinity_common::error::CommonError;

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
    let g = read_rec(KIND_GRAPH_OUT)?;
    let r = read_rec(KIND_GRAPH_READS)?;
    let result = run_component(&g, &r, &component_params(base, params))?;
    writer.lock().unwrap().append(
        id,
        KIND_ALL_PROB_PATHS,
        result.all_prob_paths_fasta.as_bytes(),
    )?;
    Ok(())
}
