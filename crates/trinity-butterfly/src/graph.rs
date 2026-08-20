//! 镜像 JUNG `DirectedSparseGraph<SeqVertex, SimpleEdge>` + `SeqVertex.java` / `SimpleEdge.java`。
//!
//! 语义要点：
//! - 节点以 id 为 key（Java SeqVertex.hashCode == id），禁止平行边（JUNG
//!   DirectedSparseGraph.addEdge 对重复 (u,v) 抛异常；这里重复 add 视为替换）。
//! - 邻接表保留插入序（方便确定性遍历/测试；JUNG 本身是 HashSet 无序）。
//! - 特殊节点 id：ROOT = -1（名字 "S"），T_VERTEX = -2（名字 "E"）。

use rustc_hash::FxHashMap;

pub const VERTEX_ROOT_ID: i32 = -1;
pub const T_VERTEX_ID: i32 = -2;

// ---------------------------------------------------------------------------
// SimpleEdge (SimpleEdge.java)
// ---------------------------------------------------------------------------

/// 有向边；镜像 `SimpleEdge.java`（权重 + 圈标记 + loop 计数 + repeat unroll 权重）。
#[derive(Debug, Clone, PartialEq)]
pub struct SimpleEdge {
    pub weight: f64,
    pub is_in_circle: bool,
    pub num_loops_involved: i32,
    pub from_vertex_id: i32,
    pub to_vertex_id: i32,
    pub repeat_unroll_weight: f64,
}

impl SimpleEdge {
    pub fn new(weight: f64, from_vertex_id: i32, to_vertex_id: i32) -> Self {
        Self {
            weight,
            is_in_circle: false,
            num_loops_involved: 0,
            from_vertex_id,
            to_vertex_id,
            repeat_unroll_weight: 0.0,
        }
    }
}

// ---------------------------------------------------------------------------
// SeqVertex (SeqVertex.java)
// ---------------------------------------------------------------------------

/// 图节点；镜像 `SeqVertex.java`（静态 tracker 移入 `BflyContext`）。
#[derive(Debug, Clone, PartialEq)]
pub struct SeqVertex {
    pub id: i32,
    /// kmer 或压缩合并后的全长序列。
    pub name: String,
    /// 被压缩进本节点的所有段权重（root 初始为每碱基一份）。
    pub weights: Vec<f64>,
    /// 合并历史：每次 concatVertex 追加一个 Vec<i32>（嵌套保留 v 的历史）。
    pub prev_vertices_id: Vec<Vec<i32>>,
    pub depth: i32,
    pub node_depth: i32,
    pub dfs_discovery_time: i32,
    pub dfs_finish_time: i32,
    pub is_in_circle: bool,
    pub to_be_deleted: bool,
    pub orig_butterfly_id: i32,
    // degenerate 三元组（SNP 简并码，占位；COLLAPSE_SNPs 后续任务使用）
    pub degenerative_freq: Vec<Vec<i32>>,
    pub degenerative_letters: Vec<Vec<String>>,
    pub degenerative_locations: Vec<i32>,
}

impl SeqVertex {
    /// `new SeqVertex(id, name)`：空 weights、dfs 时间 -1、node_depth -1。
    pub fn new(id: i32, name: impl Into<String>) -> Self {
        Self {
            id,
            name: name.into(),
            weights: Vec::new(),
            prev_vertices_id: Vec::new(),
            depth: 0,
            node_depth: -1,
            dfs_discovery_time: -1,
            dfs_finish_time: -1,
            is_in_circle: false,
            to_be_deleted: false,
            orig_butterfly_id: id,
            degenerative_freq: Vec::new(),
            degenerative_letters: Vec::new(),
            degenerative_locations: Vec::new(),
        }
    }

    /// `new SeqVertex(id, name, wei)`：每个碱基一份权重 wei。
    pub fn with_per_base_weight(id: i32, name: impl Into<String>, wei: f64) -> Self {
        let mut v = Self::new(id, name);
        let len = v.name.len();
        v.weights = vec![wei; len];
        v
    }

    /// `getWeightAvg()`：Java `(int) Math.round(sum/size)` == `floor(x + 0.5)`；空 → -1。
    pub fn get_weight_avg(&self) -> i64 {
        if self.weights.is_empty() {
            return -1;
        }
        let sum: f64 = self.weights.iter().sum();
        let avg = sum / self.weights.len() as f64;
        (avg + 0.5).floor() as i64
    }

    /// `getWeightSum()`。
    pub fn get_weight_sum(&self) -> f64 {
        self.weights.iter().sum()
    }

    /// `getNameKmerAdj()`（需要调用方给出"是否有前驱"，Java 通过静态 _graph 查询）：
    /// - id < 0（ROOT/T）→ ""
    /// - 无前驱（root）→ 全名
    /// - 有前驱 → 去掉与前驱共享的前 K-1 个碱基
    ///
    /// 名字短于一个 kmer 时 panic（Java 抛 RuntimeException）。
    pub fn get_name_kmer_adj(&self, kmer_length: usize, has_predecessors: bool) -> &str {
        if self.id < 0 {
            return "";
        }
        if !has_predecessors {
            return &self.name;
        }
        assert!(
            self.name.len() >= kmer_length,
            "ERROR, Node: {} has length shorter than a kmer",
            self.name
        );
        &self.name[kmer_length - 1..]
    }

    /// `concatVertex(vertex, w, lastRealID)`：把 vertex 压缩进本节点。
    pub fn concat_vertex(
        &mut self,
        vertex: SeqVertex,
        w: f64,
        last_real_id: i32,
        kmer_length: usize,
    ) {
        // name += v.name[K-1..]
        let suffix = vertex.get_name_kmer_adj(kmer_length, true).to_owned();
        self.name.push_str(&suffix);
        // weights: 当前边权重 + v 的所有段权重
        self.weights.push(w);
        self.weights.extend_from_slice(&vertex.weights);
        // prev id 历史（嵌套保留）
        if vertex.id <= last_real_id {
            self.prev_vertices_id.push(vec![vertex.id]);
        }
        self.prev_vertices_id.extend(vertex.prev_vertices_id);
        // degenerate 信息（占位镜像）
        if !vertex.degenerative_freq.is_empty() {
            let prev_len = self.name.len() - suffix.len();
            for &loc in &vertex.degenerative_locations {
                self.degenerative_locations.push(prev_len as i32 + loc);
            }
            self.degenerative_freq.extend(vertex.degenerative_freq);
            self.degenerative_letters
                .extend(vertex.degenerative_letters);
        }
    }

    /// `copyTheRest(v)`（SeqVertex.java L669）：拷贝 prev 历史 / dfsFinish / depth / weights。
    /// Java 拷贝的是引用（共享可变列表）；这里值拷贝，语义等价（原节点随后即被删除）。
    pub fn copy_the_rest(&mut self, v: &SeqVertex) {
        self.prev_vertices_id = v.prev_vertices_id.clone();
        self.dfs_finish_time = v.dfs_finish_time;
        self.depth = v.depth;
        self.weights = v.weights.clone();
    }

    /// `addToPrevIDs(vToKeep, vToRemove, lastRealID)`（SeqVertex.java L738）。
    ///
    /// **逐字镜像 Java quirk**：只在 `_prevVerticesID` 为空时执行；否则 no-op
    /// （Java else 分支是 `assert(true)`）。因此 `copyTheRest` 先行且 vToKeep 已有
    /// prev 历史时（压缩后的常态），本方法什么都不加。
    pub fn add_to_prev_ids(&mut self, v_keep: &SeqVertex, v_remove: &SeqVertex, last_real_id: i32) {
        if !self.prev_vertices_id.is_empty() {
            return;
        }
        let mut this_v: Vec<i32> = Vec::new();
        if v_keep.id >= last_real_id {
            this_v.push(v_keep.id);
        }
        if v_remove.id >= last_real_id {
            this_v.push(v_remove.id);
        }
        if !v_keep.prev_vertices_id.is_empty() {
            this_v.extend_from_slice(&v_keep.prev_vertices_id[0]);
        }
        if !v_remove.prev_vertices_id.is_empty() {
            this_v.extend_from_slice(&v_remove.prev_vertices_id[0]);
        }
        self.prev_vertices_id.push(this_v);
    }

    /// `clearDoubleEntriesToPrevIDs()`（SeqVertex.java L714）：删除与后一段
    /// 完全相等的相邻 prevID 段（保留后者，倒序删除保持索引有效）。
    pub fn clear_double_entries_to_prev_ids(&mut self) {
        let mut remove_indices: Vec<usize> = Vec::new();
        let n = self.prev_vertices_id.len();
        for i in 0..n.saturating_sub(1) {
            if self.prev_vertices_id[i] == self.prev_vertices_id[i + 1] {
                remove_indices.push(i);
            }
        }
        for &i in remove_indices.iter().rev() {
            self.prev_vertices_id.remove(i);
        }
    }

    /// `getLastKmer()`：不足 K 时返回全名。
    pub fn get_last_kmer(&self, kmer_length: usize) -> &str {
        let n = self.name.len();
        if n < kmer_length {
            return &self.name;
        }
        &self.name[n - kmer_length..]
    }

    /// `getFirstKmer()`。
    pub fn get_first_kmer(&self, kmer_length: usize) -> &str {
        if self.name.len() < kmer_length {
            return &self.name;
        }
        &self.name[..kmer_length]
    }

    /// `getShortSeq()`：kmerAdj 截断（>30 → 前 10 + "..." + 后 10）+ ":W{avg}"。
    pub fn get_short_seq(&self, kmer_length: usize, has_predecessors: bool) -> String {
        let adj = self.get_name_kmer_adj(kmer_length, has_predecessors);
        let mut res = if adj.len() > 30 {
            format!("{}...{}", &adj[..10], &adj[adj.len() - 10..])
        } else {
            adj.to_string()
        };
        res.push_str(&format!(":W{}", self.get_weight_avg()));
        res
    }

    /// `getShortSeqWID()`：shortSeq + "(V{id}[_{orig}]_D{node_depth})"。
    pub fn get_short_seq_wid(&self, kmer_length: usize, has_predecessors: bool) -> String {
        let mut label = format!(
            "{}(V{}",
            self.get_short_seq(kmer_length, has_predecessors),
            self.id
        );
        if self.orig_butterfly_id != self.id {
            label.push_str(&format!("_{}", self.orig_butterfly_id));
        }
        label.push_str(&format!("_D{})", self.node_depth));
        label
    }

    /// `getLongtSeqWID()`：kmerAdj + "(V{id})"。
    pub fn get_long_seq_wid(&self, kmer_length: usize, has_predecessors: bool) -> String {
        format!(
            "{}(V{})",
            self.get_name_kmer_adj(kmer_length, has_predecessors),
            self.id
        )
    }
}

// ---------------------------------------------------------------------------
// DiGraph (JUNG DirectedSparseGraph 语义)
// ---------------------------------------------------------------------------

/// 有向图：节点 id 即 key（禁平行边），邻接表保留插入序。
#[derive(Debug, Clone, Default)]
pub struct DiGraph {
    vertices: FxHashMap<i32, SeqVertex>,
    vertex_order: Vec<i32>,
    out_edges: FxHashMap<i32, Vec<i32>>,
    in_edges: FxHashMap<i32, Vec<i32>>,
    edges: FxHashMap<(i32, i32), SimpleEdge>,
}

impl DiGraph {
    pub fn new() -> Self {
        Self::default()
    }

    /// JUNG addVertex：已存在则 no-op（返回 false）。
    pub fn add_vertex(&mut self, v: SeqVertex) -> bool {
        if self.vertices.contains_key(&v.id) {
            return false;
        }
        self.vertex_order.push(v.id);
        self.vertices.insert(v.id, v);
        true
    }

    /// 不存在则创建（new(id,name)）并返回引用。
    pub fn get_or_insert_vertex(&mut self, id: i32, name: impl Into<String>) -> &mut SeqVertex {
        if !self.vertices.contains_key(&id) {
            self.add_vertex(SeqVertex::new(id, name));
        }
        self.vertices.get_mut(&id).expect("just inserted")
    }

    pub fn get_vertex(&self, id: i32) -> Option<&SeqVertex> {
        self.vertices.get(&id)
    }

    pub fn get_vertex_mut(&mut self, id: i32) -> Option<&mut SeqVertex> {
        self.vertices.get_mut(&id)
    }

    pub fn contains_vertex(&self, id: i32) -> bool {
        self.vertices.contains_key(&id)
    }

    /// JUNG removeVertex：同时清理邻接与关联边。
    pub fn remove_vertex(&mut self, id: i32) -> bool {
        if !self.vertices.contains_key(&id) {
            return false;
        }
        // 清出边
        if let Some(succs) = self.out_edges.remove(&id) {
            for s in succs {
                self.edges.remove(&(id, s));
                if let Some(preds) = self.in_edges.get_mut(&s) {
                    preds.retain(|&p| p != id);
                }
            }
        }
        // 清入边
        if let Some(preds) = self.in_edges.remove(&id) {
            for p in preds {
                self.edges.remove(&(p, id));
                if let Some(succs) = self.out_edges.get_mut(&p) {
                    succs.retain(|&s| s != id);
                }
            }
        }
        self.vertices.remove(&id);
        self.vertex_order.retain(|&v| v != id);
        true
    }

    /// JUNG addEdge：需两端节点存在；已有 (u,v) 边则替换（JUNG 会抛异常，这里宽松处理）。
    pub fn add_edge(&mut self, u: i32, v: i32, e: SimpleEdge) -> bool {
        if !self.contains_vertex(u) || !self.contains_vertex(v) {
            return false;
        }
        debug_assert!(
            !self.edges.contains_key(&(u, v)),
            "parallel edge {u}->{v}: JUNG DirectedSparseGraph would throw"
        );
        if !self.edges.contains_key(&(u, v)) {
            self.out_edges.entry(u).or_default().push(v);
            self.in_edges.entry(v).or_default().push(u);
        }
        self.edges.insert((u, v), e);
        true
    }

    pub fn remove_edge(&mut self, u: i32, v: i32) -> bool {
        if self.edges.remove(&(u, v)).is_none() {
            return false;
        }
        if let Some(succs) = self.out_edges.get_mut(&u) {
            succs.retain(|&s| s != v);
        }
        if let Some(preds) = self.in_edges.get_mut(&v) {
            preds.retain(|&p| p != u);
        }
        true
    }

    pub fn find_edge(&self, u: i32, v: i32) -> Option<&SimpleEdge> {
        self.edges.get(&(u, v))
    }

    pub fn find_edge_mut(&mut self, u: i32, v: i32) -> Option<&mut SimpleEdge> {
        self.edges.get_mut(&(u, v))
    }

    /// 后继（插入序）。
    pub fn get_successors(&self, id: i32) -> &[i32] {
        self.out_edges.get(&id).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// 前驱（插入序）。
    pub fn get_predecessors(&self, id: i32) -> &[i32] {
        self.in_edges.get(&id).map(|v| v.as_slice()).unwrap_or(&[])
    }

    pub fn out_degree(&self, id: i32) -> usize {
        self.out_edges.get(&id).map_or(0, |v| v.len())
    }

    pub fn in_degree(&self, id: i32) -> usize {
        self.in_edges.get(&id).map_or(0, |v| v.len())
    }

    /// 所有节点 id（插入序，已删节点除外）。
    pub fn vertex_ids(&self) -> &[i32] {
        &self.vertex_order
    }

    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// 所有边 ((u,v), &SimpleEdge)（按 (u,v) 排序保证确定性）。
    pub fn edges_sorted(&self) -> Vec<((i32, i32), &SimpleEdge)> {
        let mut es: Vec<_> = self.edges.iter().map(|(k, v)| (*k, v)).collect();
        es.sort_by_key(|(k, _)| *k);
        es
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weight_avg_java_round_semantics() {
        // 空 → -1
        assert_eq!(SeqVertex::new(1, "ACGT").get_weight_avg(), -1);
        // floor(x+0.5)：Java Math.round(0.5)==1, Math.round(2.5)==3, Math.round(2.4)==2
        let mut v = SeqVertex::new(1, "AC");
        v.weights = vec![0.5];
        assert_eq!(v.get_weight_avg(), 1);
        v.weights = vec![2.5];
        assert_eq!(v.get_weight_avg(), 3);
        v.weights = vec![2.4];
        assert_eq!(v.get_weight_avg(), 2);
        v.weights = vec![1.0, 2.0];
        assert_eq!(v.get_weight_avg(), 2); // 1.5 → 2
    }

    #[test]
    fn name_kmer_adj_root_vs_nonroot_vs_sink() {
        let k = 4;
        let root = SeqVertex::new(5, "ACGT");
        assert_eq!(root.get_name_kmer_adj(k, false), "ACGT");
        let inner = SeqVertex::new(6, "ACGTA");
        assert_eq!(inner.get_name_kmer_adj(k, true), "TA");
        // sink/ROOT/T
        let sink = SeqVertex::new(VERTEX_ROOT_ID, "S");
        assert_eq!(sink.get_name_kmer_adj(k, false), "");
        let t = SeqVertex::new(T_VERTEX_ID, "E");
        assert_eq!(t.get_name_kmer_adj(k, false), "");
    }

    #[test]
    fn concat_vertex_merge_order() {
        // u(name=ACGT) concat v(name=TTAG, weights=[3,3,3,3], prev=[[7]]) with edge weight 5
        let k = 4;
        let mut u = SeqVertex::new(1, "ACGT");
        u.weights = vec![1.0; 4];
        let mut v = SeqVertex::new(2, "TTAG");
        v.weights = vec![3.0; 4];
        v.prev_vertices_id = vec![vec![7]];
        u.concat_vertex(v, 5.0, 100, k);
        assert_eq!(u.name, "ACGTG"); // + v.name[K-1..] = "G"
                                     // weights: [1,1,1,1] ++ [5] ++ [3,3,3,3]
        let mut u2 = SeqVertex::new(1, "ACGT");
        u2.weights = vec![1.0; 4];
        let mut v2 = SeqVertex::new(2, "TTAG");
        v2.weights = vec![3.0; 4];
        v2.prev_vertices_id = vec![vec![7]];
        u2.concat_vertex(v2, 5.0, 100, k);
        assert_eq!(
            u2.weights,
            vec![1.0, 1.0, 1.0, 1.0, 5.0, 3.0, 3.0, 3.0, 3.0]
        );
        // v.id=2 <= last_real_id=100 → [[2]] ++ v.prev [[7]]
        assert_eq!(u2.prev_vertices_id, vec![vec![2], vec![7]]);
        // v.id > last_real_id → 只保留 v 的历史
        let mut u3 = SeqVertex::new(1, "ACGT");
        let mut v3 = SeqVertex::new(200, "TTAG");
        v3.prev_vertices_id = vec![vec![7]];
        u3.concat_vertex(v3, 5.0, 100, k);
        assert_eq!(u3.prev_vertices_id, vec![vec![7]]);
    }

    #[test]
    fn digraph_insertion_order_and_removal() {
        let mut g = DiGraph::new();
        for id in [1, 2, 3, 4] {
            g.add_vertex(SeqVertex::new(id, "ACGT"));
        }
        assert_eq!(g.vertex_ids(), &[1, 2, 3, 4]);
        g.add_edge(1, 2, SimpleEdge::new(1.0, 1, 2));
        g.add_edge(1, 3, SimpleEdge::new(2.0, 1, 3));
        g.add_edge(2, 3, SimpleEdge::new(3.0, 2, 3));
        assert_eq!(g.get_successors(1), &[2, 3]);
        assert_eq!(g.get_predecessors(3), &[1, 2]);
        assert_eq!(g.out_degree(1), 2);
        assert_eq!(g.in_degree(3), 2);
        // 重复 add 顶点 no-op
        assert!(!g.add_vertex(SeqVertex::new(1, "TTTT")));
        assert_eq!(g.get_vertex(1).unwrap().name, "ACGT");
        // 删点清邻接/边
        assert!(g.remove_vertex(3));
        assert!(!g.contains_vertex(3));
        assert_eq!(g.get_successors(1), &[2]);
        assert!(g.find_edge(1, 3).is_none());
        assert!(g.find_edge(2, 3).is_none());
        assert_eq!(g.edge_count(), 1);
        assert_eq!(g.vertex_ids(), &[1, 2, 4]);
        // 删边
        assert!(g.remove_edge(1, 2));
        assert_eq!(g.get_successors(1), &[] as &[i32]);
        assert_eq!(g.get_predecessors(2), &[] as &[i32]);
        assert!(!g.remove_edge(1, 2));
    }

    #[test]
    fn parallel_edge_replaces() {
        let mut g = DiGraph::new();
        g.add_vertex(SeqVertex::new(1, "ACGT"));
        g.add_vertex(SeqVertex::new(2, "CGTA"));
        assert!(g.add_edge(1, 2, SimpleEdge::new(1.0, 1, 2)));
        // 平行边在 JUNG 会抛异常；这里 debug_assert 提示 + 替换语义
        let dup = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            g.add_edge(1, 2, SimpleEdge::new(9.0, 1, 2));
        }));
        if cfg!(debug_assertions) {
            assert!(dup.is_err(), "duplicate add_edge should debug_assert");
        } else {
            assert_eq!(g.edge_count(), 1);
            assert_eq!(g.find_edge(1, 2).unwrap().weight, 9.0);
            assert_eq!(g.out_degree(1), 1);
        }
        // 先删再加：合法路径，权重生效
        assert!(g.remove_edge(1, 2));
        assert!(g.add_edge(1, 2, SimpleEdge::new(9.0, 1, 2)));
        assert_eq!(g.edge_count(), 1);
        assert_eq!(g.find_edge(1, 2).unwrap().weight, 9.0);
        assert_eq!(g.out_degree(1), 1);
    }

    #[test]
    fn first_last_kmer() {
        let v = SeqVertex::new(1, "ACGTT");
        assert_eq!(v.get_first_kmer(3), "ACG");
        assert_eq!(v.get_last_kmer(3), "GTT");
        assert_eq!(v.get_first_kmer(9), "ACGTT");
        assert_eq!(v.get_last_kmer(9), "ACGTT");
    }
}
