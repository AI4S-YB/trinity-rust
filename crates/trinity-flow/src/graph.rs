//! graph.out 解析 → 流图 + 弱连通分量 + 拓扑排序。

use std::collections::HashMap;
use trinity_common::error::CommonError;

#[derive(Debug, Clone)]
pub struct FlowNode {
    pub kmer: Vec<u8>,
}

#[derive(Debug, Default)]
pub struct FlowGraph {
    pub nodes: HashMap<i64, FlowNode>,
    /// (prev, cur) → weight; prev=-1 表示虚拟源边（不参与分量但参与拓扑起点判定）
    pub edges: HashMap<(i64, i64), u64>,
    pub component_id: u64,
}

impl FlowGraph {
    pub fn parse(data: &[u8]) -> Result<Self, CommonError> {
        let text = std::str::from_utf8(data)
            .map_err(|e| CommonError::Parse(format!("graph.out UTF-8: {e}")))?;
        let mut g = FlowGraph::default();
        for (i, line) in text.lines().enumerate() {
            if line.starts_with("Component") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                g.component_id = parts.get(1).and_then(|s| s.parse().ok()).unwrap_or(0);
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 4 {
                if i > 0 && !line.is_empty() {
                    return Err(CommonError::Parse(format!("graph.out 行格式: {line}")));
                }
                continue;
            }
            let node: i64 = f[0]
                .parse()
                .map_err(|_| CommonError::Parse(format!("node id: {}", f[0])))?;
            let prev: i64 = f[1]
                .parse()
                .map_err(|_| CommonError::Parse(format!("prev id: {}", f[1])))?;
            let w: u64 = f[2].parse().unwrap_or(0);
            g.nodes.entry(node).or_insert(FlowNode {
                kmer: f[3].as_bytes().to_vec(),
            });
            if prev >= 0 {
                g.nodes.entry(prev).or_insert(FlowNode { kmer: Vec::new() });
                g.edges.insert((prev, node), w);
            }
        }
        // 填充 prev 节点的 kmer（其在以 node 身份出现时有值）
        let kmers: HashMap<i64, Vec<u8>> = g
            .nodes
            .iter()
            .filter(|(_, n)| !n.kmer.is_empty())
            .map(|(&id, n)| (id, n.kmer.clone()))
            .collect();
        for (id, n) in g.nodes.iter_mut() {
            if n.kmer.is_empty() {
                n.kmer = kmers.get(id).cloned().unwrap_or_default();
            }
        }
        Ok(g)
    }

    pub fn edge_weight(&self, prev: i64, cur: i64) -> Option<u64> {
        self.edges.get(&(prev, cur)).copied()
    }

    /// 邻接表（前向）
    pub fn successors(&self, id: i64) -> Vec<i64> {
        let mut v: Vec<i64> = self
            .edges
            .keys()
            .filter(|&&(p, _)| p == id)
            .map(|&(_, c)| c)
            .collect();
        v.sort_unstable();
        v
    }

    pub fn predecessors(&self, id: i64) -> Vec<i64> {
        let mut v: Vec<i64> = self
            .edges
            .keys()
            .filter(|&(_, c)| c == id)
            .map(|&(p, _)| p)
            .collect();
        v.sort_unstable();
        v
    }

    /// 弱连通分量（并查集, 无向）
    pub fn weak_components(&self) -> Vec<Vec<i64>> {
        let mut parent: HashMap<i64, i64> = self.nodes.keys().map(|&id| (id, id)).collect();
        fn find(p: &mut HashMap<i64, i64>, x: i64) -> i64 {
            let mut r = x;
            while p[&r] != r {
                r = p[&r];
            }
            let mut c = x;
            while p[&c] != r {
                let nxt = p[&c];
                p.insert(c, r);
                c = nxt;
            }
            r
        }
        for &(u, v) in self.edges.keys() {
            let (ru, rv) = (find(&mut parent, u), find(&mut parent, v));
            if ru != rv {
                parent.insert(ru, rv);
            }
        }
        let mut groups: HashMap<i64, Vec<i64>> = HashMap::new();
        for &id in self.nodes.keys() {
            groups.entry(find(&mut parent, id)).or_default().push(id);
        }
        let mut out: Vec<Vec<i64>> = groups.into_values().collect();
        for v in out.iter_mut() {
            v.sort_unstable();
        }
        out.sort_by_key(|v| v[0]);
        out
    }

    /// 分量内 Kahn 拓扑排序（有环返回 None）
    pub fn toposort_component(&self, comp: &[i64]) -> Option<Vec<i64>> {
        let set: std::collections::HashSet<i64> = comp.iter().copied().collect();
        let mut indeg: HashMap<i64, usize> = comp.iter().map(|&id| (id, 0usize)).collect();
        for &(u, v) in self.edges.keys() {
            if set.contains(&u) && set.contains(&v) {
                *indeg.entry(v).or_default() += 1;
            }
        }
        let mut queue: std::collections::BinaryHeap<std::cmp::Reverse<i64>> = indeg
            .iter()
            .filter(|(_, &d)| d == 0)
            .map(|(&id, _)| std::cmp::Reverse(id))
            .collect();
        let mut order = Vec::with_capacity(comp.len());
        while let Some(std::cmp::Reverse(id)) = queue.pop() {
            order.push(id);
            for succ in self.successors(id) {
                if let Some(d) = indeg.get_mut(&succ) {
                    *d -= 1;
                    if *d == 0 {
                        queue.push(std::cmp::Reverse(succ));
                    }
                }
            }
        }
        if order.len() == comp.len() {
            Some(order)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINI: &str = "Component 0\n\
        1\t-1\t5\tAAAACCCCGGGGTTTTAAAACCCC\t1\n\
        2\t1\t3\tAAACCCCGGGGTTTTAAAACCCCA\t1\n\
        3\t1\t2\tAAACCCCGGGGTTTTAAAACCCCG\t1\n\
        4\t2\t3\tAACCCCGGGGTTTTAAAACCCCAA\t1\n\
        5\t3\t2\tACCCCGGGGTTTTAAAACCCCAAA\t1\n\
        10\t-1\t1\tTTTTAAAACCCCAAAGGGGCCCC\t1\n\
        11\t10\t1\tTTTAAAACCCCAAAGGGGCCCAA\t1\n";

    #[test]
    fn parse_and_components() {
        let g = FlowGraph::parse(MINI.as_bytes()).unwrap();
        assert_eq!(g.nodes.len(), 7);
        assert_eq!(g.edge_weight(1, 2), Some(3));
        assert_eq!(g.edge_weight(1, 3), Some(2));
        // 两个弱连通分量: {1..5} 和 {10,11}
        let comps = g.weak_components();
        assert_eq!(comps.len(), 2);
        assert!(comps.iter().any(|c| c.len() == 5));
        assert!(comps.iter().any(|c| c.len() == 2));
        // 节点 1 与 4 的 kmer 长度 = 24
        assert_eq!(g.nodes[&1].kmer.len(), 24);
    }

    #[test]
    fn topological_order_dag() {
        let g = FlowGraph::parse(MINI.as_bytes()).unwrap();
        let order = g.toposort_component(&g.weak_components()[0]).unwrap();
        assert_eq!(order.len(), 5);
        // 1 在 2/3 之前; 2 在 4 之前; 3 在 5 之前
        let pos = |id: i64| order.iter().position(|&x| x == id).unwrap();
        assert!(pos(1) < pos(2) && pos(1) < pos(3));
        assert!(pos(2) < pos(4) && pos(3) < pos(5));
    }
}
