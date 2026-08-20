//! Butterfly 全局状态（TransAssembly_allProbPaths.java L49-170 的实例级归置，
//! 替代 Java 静态字段与 SeqVertex 静态 tracker）。

use rustc_hash::FxHashMap;

/// 运行期上下文：LAST_ID/LAST_REAL_ID、阈值参数、节点 tracker。
#[derive(Debug, Clone)]
pub struct BflyContext {
    pub last_id: i32,
    pub last_real_id: i32,
    /// EDGE_THR，默认 0.02。
    pub edge_thr: f64,
    /// FLOW_THR，默认 0.02。
    pub flow_thr: f64,
    /// MAX_READ_SEQ_DIVERGENCE，默认 0.05。
    pub max_read_seq_divergence: f64,
    /// MAX_READ_LOCAL_SEQ_DIVERGENCE，默认 0.1。
    pub max_read_local_seq_divergence: f64,
    /// READ_END_PATH_TRIM_LENGTH，默认 0。
    pub read_end_path_trim_length: usize,
    /// USE_DP_READ_TO_VERTEX_ALIGN，默认 true。
    pub use_dp_read_to_vertex_align: bool,
    /// KMER_SIZE（由 graph 文件首 kmer 推断）。
    pub kmer_size: usize,
    /// nodeTracker：id → 是否仍在图中（Java: Map<Integer,SeqVertex>，含已从图移除者）。
    pub node_tracker: FxHashMap<i32, i32>,
    /// origIDnodeTracker：orig_butterfly_id → 图中该 orig id 的所有节点。
    pub orig_id_node_tracker: FxHashMap<i32, Vec<i32>>,
}

impl Default for BflyContext {
    fn default() -> Self {
        Self::new()
    }
}

impl BflyContext {
    pub fn new() -> Self {
        Self {
            last_id: -1,
            last_real_id: -1,
            edge_thr: 0.02,
            flow_thr: 0.02,
            max_read_seq_divergence: 0.05,
            max_read_local_seq_divergence: 0.1,
            read_end_path_trim_length: 0,
            use_dp_read_to_vertex_align: true,
            kmer_size: 0,
            node_tracker: FxHashMap::default(),
            orig_id_node_tracker: FxHashMap::default(),
        }
    }

    /// `getNextID()`：++LAST_ID。
    pub fn get_next_id(&mut self) -> i32 {
        self.last_id += 1;
        self.last_id
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_match_java_globals() {
        let ctx = BflyContext::new();
        assert_eq!(ctx.last_id, -1);
        assert_eq!(ctx.last_real_id, -1);
        assert_eq!(ctx.edge_thr, 0.02);
        assert_eq!(ctx.flow_thr, 0.02);
        assert_eq!(ctx.max_read_seq_divergence, 0.05);
        assert_eq!(ctx.max_read_local_seq_divergence, 0.1);
        assert_eq!(ctx.read_end_path_trim_length, 0);
        assert!(ctx.use_dp_read_to_vertex_align);
        assert_eq!(ctx.kmer_size, 0);
    }

    #[test]
    fn next_id_increments() {
        let mut ctx = BflyContext::new();
        ctx.last_id = 10;
        assert_eq!(ctx.get_next_id(), 11);
        assert_eq!(ctx.get_next_id(), 12);
        assert_eq!(ctx.last_id, 12);
    }
}
