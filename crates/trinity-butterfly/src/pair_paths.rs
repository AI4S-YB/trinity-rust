//! P4 Butterfly：PairPath 兼容性族 + getSuffStats_wPairs 配对组装。
//!
//! 镜像 trinityrnaseq v2.15.2 `Butterfly/Butterfly/src/src/PairPath.java` 与
//! `TransAssembly_allProbPaths.java` 的 `getSuffStats_wPairs`（L9294-9410）。
//!
//! Java `PairPath.equals/hashCode` = 两条路径列表的整体相等 —— Rust 直接派生
//! `PartialEq/Eq/Hash`（HashMap 键语义等价于 Java 的 `path1|path2` 字符串键）。
//!
//! `node_is_contained_or_possibly_in_gap` 的 Dijkstra 部分（isAncestral 可达性）
//! 以闭包注入：本任务只实现主干（containsID 短路），gap 检查由 T8 的
//! DijkstraDistanceWoVer 供参。

use rustc_hash::FxHashMap;

use crate::threading::ReadPath;

/// Java `PairPath`：双端 read 的两条路径（单端只占 path1）。
#[derive(Debug, Clone, Default, PartialEq, Eq, Hash)]
pub struct PairPath {
    pub path1: Vec<i32>,
    pub path2: Vec<i32>,
    is_circular: bool,
}

impl PairPath {
    /// Java `PairPath(path1)`：path2 置空。
    pub fn new(path1: Vec<i32>) -> Self {
        Self {
            path1,
            path2: Vec::new(),
            is_circular: false,
        }
    }

    /// Java `PairPath(path1, path2)`。
    pub fn with_pair(path1: Vec<i32>, path2: Vec<i32>) -> Self {
        Self {
            path1,
            path2,
            is_circular: false,
        }
    }

    pub fn is_circular(&self) -> bool {
        self.is_circular
    }

    pub fn set_circular(&mut self) {
        self.is_circular = true;
    }

    pub fn is_empty(&self) -> bool {
        self.path1.is_empty() && self.path2.is_empty()
    }

    pub fn has_second_path(&self) -> bool {
        !self.path2.is_empty()
    }

    /// Java `getFirstID`：path1 首节点；空则 -10。
    pub fn get_first_id(&self) -> i32 {
        if self.is_empty() {
            -10
        } else {
            self.path1[0]
        }
    }

    /// Java `getLastID`：有 path2 取 path2 末节点，否则 path1 末节点（空则 -10）。
    pub fn get_last_id(&self) -> i32 {
        if !self.has_second_path() {
            self.get_last_id_path1()
        } else {
            *self.path2.last().unwrap()
        }
    }

    /// Java `getLastID_path1`（空则 -10）。
    pub fn get_last_id_path1(&self) -> i32 {
        if self.is_empty() {
            -10
        } else {
            *self.path1.last().unwrap()
        }
    }

    /// Java `getFirstID_path2`（无 path2 则 None）。
    pub fn get_first_id_path2(&self) -> Option<i32> {
        if self.has_second_path() {
            Some(self.path2[0])
        } else {
            None
        }
    }

    /// Java `containsID`：**未 trim** 的两条列表直接 contains。
    pub fn contains_id(&self, id: i32) -> bool {
        self.path1.contains(&id) || self.path2.contains(&id)
    }

    /// Java `getMaxPathLength`。
    pub fn get_max_path_length(&self) -> usize {
        self.path1.len().max(self.path2.len())
    }

    /// Java `trimSinkNodes()`（实例版）：sink 节点（id < 0）只 trim 首尾。
    pub fn trim_sink_nodes(&self) -> PairPath {
        PairPath::with_pair(trim_sink_nodes(&self.path1), trim_sink_nodes(&self.path2))
    }

    /// Java 实例版 `haveAnyNodeInCommon(List<Integer> path)`（sink 不计）。
    pub fn have_any_node_in_common_path(&self, path: &[i32]) -> bool {
        if path.is_empty() || self.is_empty() {
            return false;
        }
        let path_t = trim_sink_nodes(path);
        for &id in &trim_sink_nodes(&self.path1) {
            if path_t.contains(&id) {
                return true;
            }
        }
        if self.has_second_path() {
            for &id in &trim_sink_nodes(&self.path2) {
                if path_t.contains(&id) {
                    return true;
                }
            }
        }
        false
    }

    /// Java 实例版 `haveAnyNodeInCommon(PairPath other)`（sink 不计）。
    pub fn have_any_node_in_common_pair(&self, other: &PairPath) -> bool {
        let mut path: Vec<i32> = trim_sink_nodes(&other.path1);
        if other.has_second_path() {
            path.extend(trim_sink_nodes(&other.path2));
        }
        for &id in &trim_sink_nodes(&self.path1) {
            if path.contains(&id) {
                return true;
            }
        }
        if self.has_second_path() {
            for &id in &trim_sink_nodes(&self.path2) {
                if path.contains(&id) {
                    return true;
                }
            }
        }
        false
    }

    /// Java `isCompatible(List<Integer>)`：包装为 path2 空的 PairPath 再比较。
    pub fn is_compatible_path(&self, path: &[i32]) -> bool {
        self.is_compatible(&PairPath::new(path.to_vec()))
    }

    /// Java `isCompatible(PairPath)`（L361）：任一对路径有公共节点即要求逐对
    /// compatible；完全无公共节点则返回 false（have_overlap 语义）。
    pub fn is_compatible(&self, other: &PairPath) -> bool {
        let mut have_overlap = false;
        for my_path in [&self.path1, &self.path2] {
            for other_path in [&other.path1, &other.path2] {
                if have_any_node_in_common_paths(my_path, other_path) {
                    have_overlap = true;
                    if !individual_paths_are_compatible(my_path, other_path) {
                        return false;
                    }
                }
            }
        }
        have_overlap
    }

    /// Java `isCompatibleAndContainedByPairPath(other_pp)`：本 pairpath 的
    /// path1/path2 都是 other 的 subPath。
    pub fn is_compatible_and_contained_by_pair_path(&self, other_pp: &PairPath) -> bool {
        if !other_pp.contains_sub_path(&self.path1) {
            return false;
        }
        if self.has_second_path() && !other_pp.contains_sub_path(&self.path2) {
            return false;
        }
        true
    }

    /// Java `isCompatibleAndContainedBySinglePath`（L472）：
    /// read（path1/path2）完整落在 path 内，且从 read 首节点开始逐节点相等。
    pub fn is_compatible_and_contained_by_single_path(&self, path: &[i32]) -> bool {
        let path = trim_sink_nodes(path);
        let first_path = trim_sink_nodes(&self.path1);

        if !self.have_any_node_in_common_path(&path) {
            return false;
        }

        let Some(first_common_node) = get_first_common_id(&path, &first_path) else {
            return false; // first path 不在 path 内
        };
        let mut i = path.iter().position(|&x| x == first_common_node).unwrap();
        let mut j = first_path
            .iter()
            .position(|&x| x == first_common_node)
            .unwrap();
        if j != 0 {
            return false; // 必须从 read 首节点开始
        }
        while i < path.len() && j < first_path.len() {
            if path[i] != first_path[j] {
                return false;
            }
            i += 1;
            j += 1;
        }
        if j != first_path.len() {
            return false; // path 未完整包含 read
        }

        // 第二条路径（若存在）
        if self.has_second_path() {
            let second_path = trim_sink_nodes(&self.path2);
            let Some(common_node) = get_first_common_id(&path, &second_path) else {
                return false; // 第二条路径不在 path 内
            };
            let mut i = path.iter().position(|&x| x == common_node).unwrap();
            let mut j = second_path.iter().position(|&x| x == common_node).unwrap();
            if j != 0 {
                return false;
            }
            while i < path.len() && j < second_path.len() {
                if path[i] != second_path[j] {
                    return false;
                }
                i += 1;
                j += 1;
            }
            if j != second_path.len() {
                return false;
            }
        }
        true
    }

    /// Java `containsSubPath`（L579）：path（sub）被本 pairpath 的 path1（或
    /// path1→path2 接力，允许 discontiguous）覆盖。
    pub fn contains_sub_path(&self, path: &[i32]) -> bool {
        let path = trim_sink_nodes(path);
        let first_path = trim_sink_nodes(&self.path1);

        if !self.have_any_node_in_common_path(&path) {
            return false;
        }

        let first_common_node = get_first_common_id(&path, &first_path);

        let mut last_index_of_path_covered: isize = -1;

        if let Some(first_common_node) = first_common_node {
            let mut i = path.iter().position(|&x| x == first_common_node).unwrap() as isize;
            let mut j = first_path
                .iter()
                .position(|&x| x == first_common_node)
                .unwrap() as isize;

            if i != 0 {
                return false; // 必须从 path 首节点开始
            }
            while i < path.len() as isize && j < first_path.len() as isize {
                if path[i as usize] != first_path[j as usize] {
                    return false;
                }
                i += 1;
                j += 1;
            }
            if !(i == path.len() as isize || j == first_path.len() as isize) {
                return false; // 未在 read 范围内走完 subPath
            } else if i == path.len() as isize {
                return true; // path1 已完整覆盖
            } else {
                last_index_of_path_covered = i - 1;
            }
        }

        // 第二条路径（若存在）
        if self.has_second_path() {
            let second_path = trim_sink_nodes(&self.path2);
            if let Some(common_node) = get_first_common_id(&path, &second_path) {
                let mut i = path.iter().position(|&x| x == common_node).unwrap() as isize;
                let mut j = second_path.iter().position(|&x| x == common_node).unwrap() as isize;

                if i > last_index_of_path_covered + 1 {
                    return false; // path 中间有未覆盖段
                }
                if i > 0 && j != 0 {
                    return false; // 已在 path 上行走时必须从第二 read 首节点起
                }
                while i < path.len() as isize && j < second_path.len() as isize {
                    if path[i as usize] != second_path[j as usize] {
                        return false;
                    }
                    i += 1;
                    j += 1;
                }
                if i == path.len() as isize {
                    return true; // 走完 path 其余部分
                }
            }
        }
        false
    }

    /// Java `node_is_contained_or_possibly_in_gap`（L854）。
    ///
    /// 主干：containsID 命中即 true。gap 检查（last(path1) → node → first(path2)
    /// 双向可达）通过 `is_ancestral(a, b)` 闭包注入（Java 传 DijkstraDistanceWoVer，
    /// `isAncestral > 0` = a→b 正向可达）。本任务主线默认不依赖 gap 分支
    /// （subPathHasEnoughReadSupport 的调用方在 T8 接 Dijkstra）。
    pub fn node_is_contained_or_possibly_in_gap(
        &self,
        node_id: i32,
        is_ancestral: &dyn Fn(i32, i32) -> i32,
    ) -> bool {
        if self.contains_id(node_id) {
            return true;
        }

        if self.has_second_path() {
            let last_v = self.get_last_id_path1();
            let first_v = self.get_first_id_path2().unwrap();
            // lastV --> node --> firstV
            if is_ancestral(last_v, node_id) > 0 && is_ancestral(node_id, first_v) > 0 {
                return true;
            }
        }
        false
    }

    /// Java `toString` 键语义（HashMap 键显示）：`PairPath [_paths=[[p1],[p2]]]`。
    pub fn to_key_string(&self) -> String {
        format!(
            "PairPath [_paths=[{},{}]]",
            i32_list_java_repr(&self.path1),
            i32_list_java_repr(&self.path2)
        )
    }
}

fn i32_list_java_repr(v: &[i32]) -> String {
    format!(
        "[{}]",
        v.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(", ")
    )
}

/// Java 静态 `getFirstCommonID`（L226）：取**较短**路径（等长取 pathA）中第一个
/// 出现在另一条路径中的节点。
pub fn get_first_common_id(path_a: &[i32], path_b: &[i32]) -> Option<i32> {
    let (longer, shorter) = if path_a.len() > path_b.len() {
        (path_a, path_b)
    } else {
        (path_b, path_a)
    };
    shorter.iter().find(|&&i| longer.contains(&i)).copied()
}

/// Java 静态 `haveAnyNodeInCommon(pathA, pathB)`（L250，先 trimSinkNodes）。
pub fn have_any_node_in_common_paths(path_a: &[i32], path_b: &[i32]) -> bool {
    let a = trim_sink_nodes(path_a);
    let b = trim_sink_nodes(path_b);
    if a.is_empty() || b.is_empty() {
        return false;
    }
    a.iter().any(|&id| b.contains(&id))
}

/// Java 静态 `trimSinkNodes(List<Integer>)`（L678）：首尾 sink 节点（id<0）截除。
pub fn trim_sink_nodes(path: &[i32]) -> Vec<i32> {
    let mut p = path.to_vec();
    if p.first().is_some_and(|&x| x < 0) {
        p.remove(0);
    }
    if p.last().is_some_and(|&x| x < 0) {
        p.pop();
    }
    p
}

/// Java 静态 `individual_paths_are_compatible`（L401）：
/// 找第一个公共节点，要求至少一方从重叠区首节点起（i==0 || j==0），
/// 随后逐节点相等走到某一方尽头。
pub fn individual_paths_are_compatible(path_a_in: &[i32], path_b_in: &[i32]) -> bool {
    let path_a = trim_sink_nodes(path_a_in);
    let path_b = trim_sink_nodes(path_b_in);

    if !have_any_node_in_common_paths(&path_a, &path_b) {
        return false;
    }

    let first_common_node = get_first_common_id(&path_a, &path_b).unwrap();

    let mut i = path_a.iter().position(|&x| x == first_common_node).unwrap();
    let mut j = path_b.iter().position(|&x| x == first_common_node).unwrap();

    if !(i == 0 || j == 0) {
        return false; // 一方须从重叠区首节点开始
    }

    while i < path_a.len() && j < path_b.len() {
        if path_a[i] != path_b[j] {
            return false;
        }
        i += 1;
        j += 1;
    }
    true
}

/// `getSuffStats_wPairs` 输出（Java 静态全局 LONG_READ_* 一并收进结构体）。
#[derive(Debug, Default)]
pub struct SuffStats {
    /// start vertex（= PairPath path1 首节点）→ PairPath → 计数。
    pub combined_read_hash: FxHashMap<i32, FxHashMap<PairPath, i64>>,
    pub num_singletons: i64,
    pub num_pairs: i64,
    pub num_pairs_discarded: i64,
    pub num_reads_used: i64,
    /// Java `LONG_READ_NAME_TO_PPath`（"LR$|" 前缀长读）。
    pub long_read_name_to_ppath: FxHashMap<String, PairPath>,
    /// Java `LONG_READ_PATH_MAP`。
    pub long_read_path_map: FxHashMap<PairPath, Vec<String>>,
}

impl SuffStats {
    /// combined_read_hash 的总计数（= num_reads_used，自洽校验用）。
    pub fn total_count(&self) -> i64 {
        self.combined_read_hash
            .values()
            .flat_map(|m| m.values())
            .sum()
    }

    /// 不同 start vertex 数。
    pub fn num_start_vertices(&self) -> usize {
        self.combined_read_hash.len()
    }
}

/// Java `getSuffStats_wPairs`（L9294）：按 read 名分组（同 T5 的 occurrence 序），
/// 名下 1 条 → 单端 PairPath(path1)；≥2 条 → 取**前两条**（Java curList.get(0/1)，
/// 更多 occurrence 忽略）。计数按 (firstV, PairPath) 累积。
/// "LR$|" 前缀名单读记入长读映射。
///
/// 注：v2.15.2 源码里 combinePaths 调用被注释（"move this to after repeat
/// unrolling"），但**发布的 Butterfly.jar 实际仍执行合并**——c1 对拍实证：
/// jar Init 的 `[3768]=5030` = 我们未合并的 `[3768]=1357 + [3768|3768]=3673`、
/// `[3280,3298]=21` = `6 + [3280,3298|3768]=15`、`[3768|3298]=254` =
/// `[3298|3768]=254`（path2→path1 交换）。不合并会使后续 POG 重映射把配对
/// 全部落进"多映射拆分 support=1"分支，读支持塌缩、路径搜索截断（c1 症状）。
/// 合并失败（空 PairPath）→ num_pairs_discarded，丢弃。
/// EXPERIMENT B helper：p1 在前、p2 在后，尝试唯一后继 imputation 合并；
/// 失败则保留配对。
fn impute_merge(graph: &crate::graph::DiGraph, p1: Vec<i32>, p2: Vec<i32>) -> PairPath {
    use crate::paths::is_ancestral;
    let (mut p1, p2) = (p1, p2);
    let l1 = p1[p1.len() - 1];
    let f2 = p2[0];
    if is_ancestral(graph, l1, f2) > 0 {
        let mut impute = true;
        let mut v = l1;
        loop {
            let mut next: Option<i32> = None;
            let mut count = 0;
            for &succ in graph.get_successors(v) {
                if is_ancestral(graph, succ, f2) > 0 {
                    count += 1;
                    next = Some(succ);
                }
            }
            match next {
                Some(n) if count == 1 => {
                    if n == f2 {
                        break;
                    }
                    p1.push(n);
                    v = n;
                }
                _ => {
                    impute = false;
                    break;
                }
            }
        }
        if impute {
            p1.extend(p2);
            return PairPath::new(p1);
        }
    }
    PairPath::with_pair(p1, p2)
}

/// Java `combinePaths`（L9417）的完整镜像（jar 实际在 getSuffStats_wPairs 中
/// 执行——源码树该调用被注释但发布版未注释，见上方对拍记录）。
/// 返回 `None` = 空 PairPath（两端无一致方向且互不包含）→ 调用方丢弃。
///
/// else-if 链的顺序与条件逐条对齐：
/// 1. path1 ⊇ path2 → 单端 path1；2. path2 ⊇ path1 → 单端 path2；
/// 3. last1→first2 可达且不等 → (p1,p2)；4. last2→first1 可达且不等 → (p2,p1)；
/// 5. first2 与 first1 无向、last2 与 last1 无向 → 空（不合并，overlap 分支
///    也被跳过）；6/7. overlap：first1→first2 可达且 p1 含 first2 →
///    p1[..i]+p2；反向同理。最后统一 imputation（见 impute_merge）。
pub fn combine_paths(
    graph: &crate::graph::DiGraph,
    path1: &[i32],
    path2: &[i32],
) -> Option<PairPath> {
    use crate::paths::is_ancestral;
    let contains_all = |a: &[i32], b: &[i32]| b.iter().all(|x| a.contains(x));

    let (p1, p2): (Vec<i32>, Vec<i32>);
    let mut has_pair = false;

    let l1 = path1[path1.len() - 1];
    let f1 = path1[0];
    let l2 = path2[path2.len() - 1];
    let f2 = path2[0];

    if contains_all(path1, path2) {
        p1 = path1.to_vec();
        p2 = Vec::new();
    } else if contains_all(path2, path1) {
        // Java setPath2 → 随后 movePath2To1
        p1 = path2.to_vec();
        p2 = Vec::new();
    } else if is_ancestral(graph, l1, f2) > 0 && l1 != f2 {
        p1 = path1.to_vec();
        p2 = path2.to_vec();
        has_pair = true;
    } else if is_ancestral(graph, l2, f1) > 0 && l2 != f1 {
        p1 = path2.to_vec();
        p2 = path1.to_vec();
        has_pair = true;
    } else if is_ancestral(graph, f2, f1) == 0 && is_ancestral(graph, l2, l1) == 0 {
        // 两端无一致方向：combinePaths 返回空 → 丢弃
        return None;
    } else if is_ancestral(graph, f1, f2) > 0 && path1.contains(&f2) {
        let i = path1.iter().position(|&v| v == f2).unwrap();
        p1 = [&path1[..i], path2].concat();
        p2 = Vec::new();
    } else if is_ancestral(graph, f2, f1) > 0 && path2.contains(&f1) {
        let i = path2.iter().position(|&v| v == f1).unwrap();
        p1 = [&path2[..i], path1].concat();
        p2 = Vec::new();
    } else {
        // 无分支命中：Java 返回空 PairPath → 丢弃
        return None;
    }

    if has_pair {
        Some(impute_merge(graph, p1, p2))
    } else {
        Some(PairPath::new(p1))
    }
}

pub fn get_suff_stats_w_pairs(
    graph: &crate::graph::DiGraph,
    read_paths_ordered: &[(String, ReadPath)],
) -> SuffStats {
    // 名 → 该名下所有成功穿线路径（occurrence 序，镜像 readNameHash）
    let mut read_name_hash: FxHashMap<&str, Vec<&ReadPath>> = FxHashMap::default();
    for (name, p) in read_paths_ordered {
        read_name_hash.entry(name.as_str()).or_default().push(p);
    }

    let mut stats = SuffStats::default();

    for (name, cur_list) in read_name_hash {
        let path: PairPath = if cur_list.len() == 1 {
            // 单端 read
            let path = PairPath::new(cur_list[0].path.clone());
            stats.num_singletons += 1;

            // 长读特判（LR$| 前缀）
            if let Some(lr_name) = name.strip_prefix("LR$|") {
                let _ = lr_name;
                stats
                    .long_read_name_to_ppath
                    .insert(name.to_string(), path.clone());
                stats
                    .long_read_path_map
                    .entry(path.clone())
                    .or_default()
                    .push(name.to_string());
            }
            path
        } else {
            // 配对 read：取前两条（Java curList.get(0)/get(1)），
            // 再 combinePaths 合并（见函数文档——jar 实际行为）
            let path1 = cur_list[0].path.clone();
            let path2 = cur_list[1].path.clone();
            match combine_paths(graph, &path1, &path2) {
                Some(p) => {
                    stats.num_pairs += 1;
                    p
                }
                None => {
                    // jar（合并生效版）：空 PairPath → numPairsDiscarded++，跳过
                    stats.num_pairs_discarded += 1;
                    continue;
                }
            }
        };
        if std::env::var_os("TR_DUMP_MERGES").is_some() && cur_list.len() > 1 {
            eprintln!(
                "MERGE p1: {:?}, p2: {:?} => {}",
                cur_list[0].path,
                cur_list[1].path,
                path.to_key_string()
            );
        }

        let first_v = path.get_first_id();
        *stats
            .combined_read_hash
            .entry(first_v)
            .or_default()
            .entry(path)
            .or_insert(0) += 1;
        stats.num_reads_used += 1;
    }

    stats
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pp(p1: &[i32], p2: &[i32]) -> PairPath {
        PairPath::with_pair(p1.to_vec(), p2.to_vec())
    }

    // ---------------- trim / common ----------------

    #[test]
    fn trim_sink_nodes_strips_leading_trailing_negatives() {
        assert_eq!(trim_sink_nodes(&[-1, 3, 5, -2]), vec![3, 5]);
        assert_eq!(trim_sink_nodes(&[-1, -2]), Vec::<i32>::new());
        assert_eq!(trim_sink_nodes(&[3, -1, 5]), vec![3, -1, 5]); // 只 trim 首尾
        assert_eq!(trim_sink_nodes(&[]), Vec::<i32>::new());
    }

    #[test]
    fn get_first_common_id_picks_from_shorter_path() {
        assert_eq!(get_first_common_id(&[1, 2, 3], &[4, 2, 5, 6]), Some(2));
        // 较短者（a）中第一个出现在 b 的：7 不在 b，2 在 → 2
        assert_eq!(get_first_common_id(&[7, 2], &[1, 2, 3]), Some(2));
        // 等长：shorter = a，取 a 中第一个出现在 b 的 → 9
        assert_eq!(get_first_common_id(&[9, 2], &[2, 9]), Some(9));
        assert_eq!(get_first_common_id(&[1, 2], &[3, 4]), None);
    }

    #[test]
    fn have_any_node_in_common_ignores_sinks() {
        assert!(have_any_node_in_common_paths(&[-1, 3, 5], &[5, 6]));
        assert!(have_any_node_in_common_paths(&[-1, 3], &[3, -2])); // trim 后 [3] vs [3]
        assert!(!have_any_node_in_common_paths(&[-1, -2], &[-1, 3])); // trim 后 a 空
        assert!(!have_any_node_in_common_paths(&[1], &[]));
    }

    // -------- individual_paths_are_compatible --------

    #[test]
    fn ipac_read_inside_path() {
        // read: --- ; path: ---------
        assert!(individual_paths_are_compatible(&[2, 3], &[1, 2, 3, 4]));
        // 镜像：read 更长
        assert!(individual_paths_are_compatible(&[1, 2, 3, 4], &[2, 3]));
    }

    #[test]
    fn ipac_partial_overlap_both_sides() {
        // read: ------- ; path:    --------
        assert!(individual_paths_are_compatible(&[1, 2, 3], &[3, 4, 5]));
        // read:    -------- ; path: -------
        assert!(individual_paths_are_compatible(&[3, 4, 5], &[1, 2, 3]));
    }

    #[test]
    fn ipac_rejects_mismatch_in_overlap() {
        // 重叠区起点条件满足（j==0）但节点不等
        assert!(!individual_paths_are_compatible(&[1, 2, 3], &[2, 9, 3]));
    }

    #[test]
    fn ipac_requires_one_starting_at_overlap() {
        // 公共节点 3，但两方都不从 3 开始（i=2, j=1）→ false
        assert!(!individual_paths_are_compatible(&[1, 2, 3, 4], &[9, 3, 4]));
    }

    #[test]
    fn ipac_no_common_node() {
        assert!(!individual_paths_are_compatible(&[1, 2], &[3, 4]));
    }

    #[test]
    fn ipac_sink_nodes_do_not_count_as_overlap() {
        // 公共节点只有 sink（-1）——trim 后无公共节点
        assert!(!individual_paths_are_compatible(&[-1, 2], &[-1, 3]));
    }

    // ---------------- is_compatible ----------------

    #[test]
    fn compatible_four_position_relationships() {
        // read 相对 path：内 / 外（read 更长）/ 部分重叠左 / 部分重叠右
        let path = pp(&[1, 2, 3, 4, 5], &[]);
        assert!(path.is_compatible_path(&[2, 3, 4])); // 内
        assert!(path.is_compatible_path(&[0, 1, 2, 3, 4, 5, 6])); // 外
        assert!(path.is_compatible_path(&[0, 1, 2])); // 左重叠
        assert!(path.is_compatible_path(&[4, 5, 6, 7])); // 右重叠
        assert!(!path.is_compatible_path(&[1, 9, 3])); // 重叠区不匹配
        assert!(!path.is_compatible_path(&[8, 9])); // 无公共节点
    }

    #[test]
    fn compatible_cross_path_via_path2() {
        // 单条路径与 path2 无公共但与 path1 有 → compatible；
        // 与 path2 有公共且不匹配 → false
        let pair = pp(&[1, 2, 3], &[7, 8]);
        assert!(pair.is_compatible_path(&[2, 3, 4]));
        // 与 path2 部分重叠（[7,8] vs [8,9]：j==0 从重叠首节点起）→ compatible
        assert!(pair.is_compatible_path(&[8, 9]));
        // 与 path2 重叠区错位（[7,8] vs [7,9]：j==0 但第二节点不等）→ false
        assert!(!pair.is_compatible_path(&[7, 9]));
    }

    #[test]
    fn compatible_pair_vs_pair() {
        let a = pp(&[1, 2, 3], &[7, 8, 9]);
        let b = pp(&[2, 3, 4], &[9, 10]);
        assert!(a.is_compatible(&b));

        let c = pp(&[2, 9, 3], &[]); // path1 重叠区错位
        assert!(!a.is_compatible(&c));

        let d = pp(&[50, 51], &[60, 61]); // 完全无公共
        assert!(!a.is_compatible(&d));
    }

    // -------- contained / contains_sub_path --------

    #[test]
    fn contained_by_single_path_happy_and_unhappy() {
        let read = pp(&[2, 3], &[4, 5]);
        assert!(read.is_compatible_and_contained_by_single_path(&[1, 2, 3, 4, 5, 6]));
        // path1 未完整包含
        assert!(!read.is_compatible_and_contained_by_single_path(&[1, 2, 3, 4]));
        // path2 在 path 中错位（[4,5] vs path 尾部 [5,4]）
        assert!(!read.is_compatible_and_contained_by_single_path(&[1, 2, 3, 5, 4]));
        // read 必须从自身首节点对齐（j != 0）：公共节点 3 在 read 中下标 1
        assert!(!read.is_compatible_and_contained_by_single_path(&[1, 3, 4]));
        // 首节点对齐（2 在 read 下标 0）即使 path 前有其他节点也可
        assert!(read.is_compatible_and_contained_by_single_path(&[3, 2, 3, 4, 5]));
        // 单端 read
        assert!(pp(&[2, 3], &[]).is_compatible_and_contained_by_single_path(&[1, 2, 3]));
        // 无公共节点
        assert!(!read.is_compatible_and_contained_by_single_path(&[8, 9]));
    }

    #[test]
    fn contained_by_single_path_discontiguous_ok() {
        // read 的两条 path 分别落在 path 的两段即可（不必连续）
        let read = pp(&[1, 2], &[5, 6]);
        assert!(read.is_compatible_and_contained_by_single_path(&[1, 2, 3, 4, 5, 6]));
    }

    #[test]
    fn contains_sub_path_single_and_pair_coverage() {
        // path 完整落在 path1
        assert!(pp(&[1, 2, 3, 4], &[]).contains_sub_path(&[2, 3, 4]));
        // path 完整落在 path2
        assert!(pp(&[1, 2], &[7, 8, 9]).contains_sub_path(&[8, 9]));
        // path1→path2 接力（discontiguous：path1=[1,2] 先耗尽，path2 接上 5,6）
        assert!(pp(&[1, 2], &[5, 6]).contains_sub_path(&[1, 2, 5, 6]));
        // path1 未耗尽就与 path 错位（[1,2,3] vs [1,2,5,6] 重叠区 3≠5）→ false
        assert!(!pp(&[1, 2, 3], &[5, 6]).contains_sub_path(&[1, 2, 5, 6]));
        // path1 覆盖后 path 中间断档（last_covered=1，path2 覆盖起点 4 > 2）
        assert!(!pp(&[1, 2], &[4, 6]).contains_sub_path(&[1, 2, 3, 4]));
        // 必须从 path 首节点开始（公共节点 3 在 path [5,3] 中下标 1）
        assert!(!pp(&[1, 2, 3], &[]).contains_sub_path(&[5, 3]));
        // 重叠区不匹配
        assert!(!pp(&[1, 2, 3], &[]).contains_sub_path(&[2, 9, 3]));
        // 无公共节点
        assert!(!pp(&[1, 2], &[]).contains_sub_path(&[8, 9]));
    }

    #[test]
    fn compatible_and_contained_by_pair_path() {
        let outer = pp(&[1, 2, 3, 4], &[8, 9, 10]);
        let inner = pp(&[2, 3], &[9, 10]);
        assert!(inner.is_compatible_and_contained_by_pair_path(&outer));
        assert!(!outer.is_compatible_and_contained_by_pair_path(&inner));
        let not_inner = pp(&[2, 3], &[11, 12]);
        assert!(!not_inner.is_compatible_and_contained_by_pair_path(&outer));
    }

    // -------- node_is_contained_or_possibly_in_gap --------

    #[test]
    fn node_in_gap_trunk_contains_id() {
        let pair = pp(&[1, 2, 3], &[6, 7]);
        // 主干：containsID 直接命中
        assert!(pair.node_is_contained_or_possibly_in_gap(2, &|_, _| 0));
        assert!(pair.node_is_contained_or_possibly_in_gap(7, &|_, _| 0));
        assert!(!pair.node_is_contained_or_possibly_in_gap(5, &|_, _| 0));
    }

    #[test]
    fn node_in_gap_dijkstra_closure_branch() {
        let pair = pp(&[1, 2, 3], &[6, 7]);
        // gap 检查：last(path1)=3 → 5 → first(path2)=6 双向可达
        let reach = |a: i32, b: i32| if a < b { 1 } else { -1 };
        assert!(pair.node_is_contained_or_possibly_in_gap(5, &reach));
        let none = |_: i32, _: i32| 0;
        assert!(!pair.node_is_contained_or_possibly_in_gap(5, &none));
        // 单端 read（无 path2）不走 gap 分支
        assert!(!pp(&[1, 2], &[]).node_is_contained_or_possibly_in_gap(9, &reach));
    }

    // ---------------- misc accessors ----------------

    #[test]
    fn accessors_and_key_string() {
        let pair = pp(&[1, 2], &[7, 8]);
        assert_eq!(pair.get_first_id(), 1);
        assert_eq!(pair.get_last_id(), 8); // 有 path2 取 path2 末节点
        assert_eq!(pair.get_last_id_path1(), 2);
        assert_eq!(pair.get_first_id_path2(), Some(7));
        assert_eq!(pp(&[1, 2], &[]).get_last_id(), 2);
        assert_eq!(PairPath::default().get_first_id(), -10);
        assert_eq!(pair.get_max_path_length(), 2);
        assert_eq!(pair.to_key_string(), "PairPath [_paths=[[1, 2],[7, 8]]]");
        // HashMap 键语义：同 path1/path2 的 PairPath 相等
        assert_eq!(pp(&[1, 2], &[7, 8]), pp(&[1, 2], &[7, 8]));
        assert_ne!(pp(&[1, 2], &[7, 8]), pp(&[7, 8], &[1, 2]));
    }

    // ---------------- getSuffStats_wPairs ----------------

    fn rp(path: &[i32]) -> ReadPath {
        ReadPath {
            mismatch_count: 0,
            path: path.to_vec(),
            positions: Vec::new(),
        }
    }

    /// combinePaths 需要的可达性图（配对两端在图上前后相接 → 保留配对）。
    fn reach_graph(edges: &[(i32, i32)]) -> crate::graph::DiGraph {
        use crate::graph::{DiGraph, SeqVertex, SimpleEdge};
        let mut g = DiGraph::new();
        let mut vids: Vec<i32> = edges.iter().flat_map(|&(a, b)| [a, b]).collect();
        vids.sort();
        vids.dedup();
        for id in vids {
            g.add_vertex(SeqVertex::new(id, "A".to_string()));
        }
        for &(a, b) in edges {
            g.add_edge(a, b, SimpleEdge::new(1.0, a, b));
        }
        g
    }

    #[test]
    fn suff_stats_single_and_paired() {
        // 11 的两个后继（12/13）都可达 20 → imputation 失败（count_connectable=2）
        // → 保留配对（注意 isAncestral(f2,f2)>0，若 11 直连 20 会反被 impute 合并）
        let g = reach_graph(&[
            (10, 11),
            (11, 12),
            (11, 13),
            (12, 20),
            (13, 20),
            (20, 21),
            (30, 31),
        ]);
        let reads = vec![
            ("r1".to_string(), rp(&[10, 11])),
            ("r1".to_string(), rp(&[20, 21])), // r1 配对（11→20 可达）
            ("r2".to_string(), rp(&[10, 11])), // r2 单端，路径同 r1.path1
            ("r3".to_string(), rp(&[30, 31])),
        ];
        let s = get_suff_stats_w_pairs(&g, &reads);
        assert_eq!(s.num_singletons, 2);
        assert_eq!(s.num_pairs, 1);
        assert_eq!(s.num_reads_used, 3);
        assert_eq!(s.num_pairs_discarded, 0);
        assert_eq!(s.total_count(), 3);

        // start vertex 分桶
        let v10 = &s.combined_read_hash[&10];
        for k in v10.keys() {
            eprintln!("KEY: {}", k.to_key_string());
        }
        assert_eq!(v10.len(), 2); // PairPath([10,11],[]) 和 PairPath([10,11],[20,21])
        assert_eq!(v10[&pp(&[10, 11], &[])], 1);
        assert_eq!(v10[&pp(&[10, 11], &[20, 21])], 1);
        assert_eq!(s.combined_read_hash[&30][&pp(&[30, 31], &[])], 1);
    }

    #[test]
    fn suff_stats_containment_pair_merges_to_single() {
        // jar 实际行为：path1 ⊇ path2 → 合并为单端（c1 的 [3768|3768]→[3768]）
        let g = reach_graph(&[(1, 2)]);
        let reads = vec![
            ("a".to_string(), rp(&[1, 2])),
            ("a".to_string(), rp(&[2])),
            ("b".to_string(), rp(&[1, 2])),
        ];
        let s = get_suff_stats_w_pairs(&g, &reads);
        assert_eq!(s.num_pairs, 1);
        assert_eq!(s.num_pairs_discarded, 0);
        // 合并后 [1,2] 单端累积计数 2
        assert_eq!(s.combined_read_hash[&1][&pp(&[1, 2], &[])], 2);
    }

    #[test]
    fn suff_stats_discards_inconsistent_pair() {
        // 两端无一致方向且互不包含 → combinePaths 返回空 → 丢弃（jar 行为）
        let g = reach_graph(&[(1, 2), (5, 6)]);
        let reads = vec![
            ("a".to_string(), rp(&[1, 2])),
            ("a".to_string(), rp(&[5, 6])),
        ];
        let s = get_suff_stats_w_pairs(&g, &reads);
        assert_eq!(s.num_pairs, 0);
        assert_eq!(s.num_reads_used, 0);
        assert_eq!(s.num_pairs_discarded, 1);
        assert!(s.combined_read_hash.is_empty());
    }

    #[test]
    fn suff_stats_same_pair_accumulates_count() {
        let g = reach_graph(&[(1, 2), (2, 3), (2, 4), (3, 5), (4, 5), (5, 6)]);
        let reads = vec![
            ("a".to_string(), rp(&[1, 2])),
            ("a".to_string(), rp(&[5, 6])),
            ("b".to_string(), rp(&[1, 2])),
            ("b".to_string(), rp(&[5, 6])),
        ];
        let s = get_suff_stats_w_pairs(&g, &reads);
        assert_eq!(s.num_pairs, 2);
        assert_eq!(s.combined_read_hash[&1][&pp(&[1, 2], &[5, 6])], 2);
    }

    #[test]
    fn suff_stats_multi_occurrence_uses_first_two() {
        // 同名三次出现：Java 只取 curList.get(0)/get(1)，第三次忽略
        let g = reach_graph(&[(1, 2), (2, 7), (2, 8), (7, 3), (8, 3), (3, 4)]);
        let reads = vec![
            ("m".to_string(), rp(&[1, 2])),
            ("m".to_string(), rp(&[3, 4])),
            ("m".to_string(), rp(&[9, 9])),
        ];
        let s = get_suff_stats_w_pairs(&g, &reads);
        assert_eq!(s.num_pairs, 1);
        assert!(s.combined_read_hash[&1].contains_key(&pp(&[1, 2], &[3, 4])));
    }

    #[test]
    fn suff_stats_long_read_prefix() {
        let g = reach_graph(&[]);
        let reads = vec![
            ("LR$|long1".to_string(), rp(&[1, 2, 3])),
            ("plain".to_string(), rp(&[4, 5])),
        ];
        let s = get_suff_stats_w_pairs(&g, &reads);
        assert_eq!(s.long_read_name_to_ppath.len(), 1);
        assert_eq!(s.long_read_name_to_ppath["LR$|long1"], pp(&[1, 2, 3], &[]));
        assert_eq!(
            s.long_read_path_map[&pp(&[1, 2, 3], &[])],
            vec!["LR$|long1"]
        );
        // "LR$" 但无 "|" 分隔不触发
        let s2 = get_suff_stats_w_pairs(&g, &[("LR$x".to_string(), rp(&[1, 2]))]);
        assert!(s2.long_read_name_to_ppath.is_empty());
    }
}
