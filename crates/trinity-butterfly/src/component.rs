//! 单组件全链驱动（库层）——原版 `java -jar Butterfly.jar -C <prefix>` 的
//! main 编排（TransAssembly_allProbPaths.java L726-912）上移。
//!
//! 输入是 **文本**（graph.out / graph.reads 内容），输出是 allProbPaths FASTA
//! 文本——CLI、集成测试（c0/c1/c2）与 xtask 对拍共用同一入口，防漂移。
//!
//! 递归保护：POG 的 DFS_add_path 与 threading 递归在深组件上可能爆默认
//! 栈（Java 默认虚拟机栈可增长到数百 MB），因此在 `stack_size` 配置的
//! 专用线程上运行（默认 256MB，Trinity 主线 JVM 经验值）。

use std::io::Write;

use rustc_hash::FxHashMap;

use crate::postprocess::PostProcessParams;
use crate::PathSearchParams;
use trinity_common::error::CommonError;

use crate::{
    build_new_graph_use_kmers, calc_sub_components_stats, compact_linear_paths,
    create_dag_from_overlap_layout, fix_extremely_high_single_edges,
    get_original_ver_ids_mapping_hash, get_read_starts_ordered, get_suff_stats_w_pairs,
    parse_graph_reads, pre_process_graph_file, print_final_paths, remove_light_edges,
    remove_single_nt_bubbles, run_butterfly_all_prob_paths_with, run_dfs2,
};

/// 组件驱动参数（镜像 Butterfly.jar getopt 面的库层子集）。
#[derive(Debug, Clone)]
pub struct ComponentParams {
    /// `-N`：total reads / fragment pairs（Java printUsage 只要求非零）。
    pub n: u64,
    /// `-L`：最短报告序列长度。
    pub min_len: usize,
    /// `-F`：配对端最大跨度（max_pair_distance）。
    pub max_pair_distance: usize,
    /// `--NO_EM_REDUCE`（Trinity 主线传 true；getopt 默认 false）。
    pub no_em_reduce: bool,
    /// `--no_path_merging`。
    pub no_path_merging: bool,
    /// `-R`：最小 read 支持（默认 2）。
    pub min_read_support: i64,
    /// 显式 reinforcement 距离（覆盖 `-O` 百分比推算；None → 25% * max_pair_distance）。
    pub path_reinforcement_distance: Option<usize>,
    /// FASTA header 里的组件名（原版 = -C 值最后一段；printFinalPaths 去后缀）。
    pub name: String,
    /// 专用线程栈大小（递归逃生口；默认 256MB）。
    pub stack_size: usize,
}

impl Default for ComponentParams {
    fn default() -> Self {
        ComponentParams {
            n: 100_000,
            min_len: 200,
            max_pair_distance: 10_000,
            no_em_reduce: false,
            no_path_merging: false,
            min_read_support: 2,
            path_reinforcement_distance: None,
            name: "component".to_string(),
            stack_size: 256 * 1024 * 1024,
        }
    }
}

/// 单组件运行结果。
pub struct ComponentResult {
    /// `<prefix>.allProbPaths.fasta` 全文（调用方决定写哪/打 stderr）。
    pub all_prob_paths_fasta: String,
    /// 最终路径（转录本）数。
    pub num_paths: usize,
}

/// 单组件全链：graph.out 文本 + graph.reads 文本 → allProbPaths FASTA 文本。
///
/// 链路 = graph 解析 → 剪枝链（B/C/H/D + dfs2 + bubble + stats）→ 穿线
/// → SuffStats → POG → 路径搜索 → T9 后处理 → printFinalPaths。
/// 在 `stack_size` 配置的专用线程上运行（递归逃生口）；失败返回 Err
/// （组件级隔离由调用方线程池做）。
pub fn run_component(
    graph_out: &str,
    graph_reads: &str,
    params: &ComponentParams,
) -> Result<ComponentResult, CommonError> {
    if params.n == 0 || params.max_pair_distance == 0 || params.min_read_support < 1 {
        return Err(CommonError::Parse(
            "缺少必选参数（-N/-F/-R 无效值）".to_string(),
        ));
    }
    // Java L700+：显式 --path_reinforcement_distance 优先，否则 O% * F（O 默认 25）
    let prd = params
        .path_reinforcement_distance
        .unwrap_or(params.max_pair_distance / 4);

    let graph_out = graph_out.to_string();
    let graph_reads = graph_reads.to_string();
    let params = params.clone();

    // 递归逃生口：POG DFS_add_path / threading 递归在大组件上会爆默认 8MB 栈。
    let handle = std::thread::Builder::new()
        .stack_size(params.stack_size)
        .spawn(move || run_component_inner(&graph_out, &graph_reads, &params, prd))
        .map_err(|e| CommonError::Parse(format!("cannot spawn component thread: {e}")))?;
    match handle.join() {
        Ok(r) => r,
        Err(_) => Err(CommonError::Parse(
            "butterfly 组件线程 panic（深递归/数据异常）".to_string(),
        )),
    }
}

fn run_component_inner(
    graph_text: &str,
    reads_text: &str,
    params: &ComponentParams,
    path_reinforcement_distance: usize,
) -> Result<ComponentResult, CommonError> {
    // 原始 kmer 图（边权重重标定基准）
    let br0 = build_new_graph_use_kmers(graph_text)
        .map_err(|e| CommonError::Parse(format!("graph 解析失败: {e}")))?;
    let mut orig_kmer_to_node: FxHashMap<String, i32> = FxHashMap::default();
    for &vid in br0.graph.vertex_ids() {
        orig_kmer_to_node.insert(br0.graph.get_vertex(vid).unwrap().name.clone(), vid);
    }

    // 主链 B/C/H/D + 穿线 + SuffStats + POG
    let br = build_new_graph_use_kmers(graph_text)
        .map_err(|e| CommonError::Parse(format!("graph 解析失败: {e}")))?;
    let (mut graph, mut ctx) = (br.graph, br.ctx);
    let (in_flow, out_flow, _) = pre_process_graph_file(graph_text);
    fix_extremely_high_single_edges(&mut graph, &out_flow, &in_flow);
    remove_light_edges(&mut graph, &ctx);
    compact_linear_paths(&mut graph, &ctx);
    run_dfs2(&mut graph);
    remove_single_nt_bubbles(&mut graph, &mut ctx);
    calc_sub_components_stats(&mut graph, params.min_len);

    let orig_id_map = get_original_ver_ids_mapping_hash(&mut graph, ctx.last_real_id);
    let raw_reads = parse_graph_reads(reads_text, ctx.kmer_size, false)
        .map_err(|e| CommonError::Parse(format!("reads 解析失败: {e}")))?;
    let ordered =
        get_read_starts_ordered(&graph, &ctx, &orig_id_map, &orig_kmer_to_node, &raw_reads);
    let suff = get_suff_stats_w_pairs(&graph, &ordered);
    let result = create_dag_from_overlap_layout(&graph, &mut ctx, &suff);

    let sp = PathSearchParams {
        kmer_size: ctx.kmer_size,
        min_output_seq: params.min_len,
        min_read_support_thr: params.min_read_support,
        path_reinforcement_distance,
        ..PathSearchParams::default()
    };
    let post = PostProcessParams {
        no_em_reduce: params.no_em_reduce,
        no_path_merging: params.no_path_merging,
        ..PostProcessParams::default()
    };
    let bfly = run_butterfly_all_prob_paths_with(
        &result.seqvertex_graph,
        &br0.graph,
        &graph,
        &result.combined_read_hash,
        &sp,
        &post,
    );

    let fasta = print_final_paths(
        &graph,
        &bfly.final_paths,
        &bfly.gene_ids,
        &params.name,
        ctx.kmer_size,
    );
    let _ = std::io::stderr().flush();
    Ok(ComponentResult {
        all_prob_paths_fasta: fasta,
        num_paths: bfly.final_paths.len(),
    })
}
