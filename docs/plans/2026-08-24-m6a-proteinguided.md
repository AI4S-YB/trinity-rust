# M6a 蛋白引导组装（探路阶段）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **状态（2026-08-24 收档）:** 8 任务全部完成; probe 四轮收敛（F1 15.1~15.3 ≈ 基线,
> 详见 opt-log「M6a 三种验收域的收敛」）——**中性偏负结果**, `--protein-guided`
> 留默认 off 的 opt-in flag。后续迭代（PE 跨接门 M6a''）与深度诊断均在计划执行期
> 追加, 未另立计划。M6b 转向 CDS 完整性度量（另行立项）。

> **提交冻结（2026-08-24 用户指示）:** 执行本计划期间**跳过所有 Commit 步骤**（本地也不提交），直到用户解除冻结。其余步骤照常。

**Goal:** `--protein-guided` 探路版：PogFlow 初步转录本经 DIAMOND blastx 选蛋白模板 → 缺口轮廓 → 原始 k-mer 图上的读框感知引导延长 → 重比对验收 → dev1 现库 probe 出信号强弱数据（Gate 数字据此另定，不在本计划内）。

**Architecture:** butterfly 池与 harvest 之间插入一个 stage：批量 blastx（一次调用）选每转录本最佳模板，按 HSP 读框把延长中的密码子映射到蛋白坐标，用「密码子匹配 × 边权重」的 3 步前瞻贪心在原始 k-mer 图（`build_new_graph_use_kmers` 重建，顶点名 = k-mer，转录本末端可直接锚定）上延长，验收不过则回退原转录本。域语法表（Pfam）不在 M6a 范围。

**Tech Stack:** 既有组件 + DIAMOND 2.2.4（外部进程，`/public/home/senior007/miniconda3/envs/diamond2go/bin/diamond`，env `TRINITY_DIAMOND` 可覆盖）。

**对 spec 的两处落地偏差（均为简化，机制不变）:**
1. 翻译不写独立 ORF finder——blastx 内部完成翻译；只写按 HSP 读框逐密码子的最小翻译器（用于走图打分）
2. `--domain-table` 旗标 M6a 不加（无表即默认行为），留 M6b
3. 模板选择用 bitscore（DIAMOND 标准最佳命中准则，隐含 coverage×identity），缺口用 N/C 两端缺失跨度（单 HSP 天然只缺两端；M/C 三分轮廓留 M6b 域语法）

---

## 文件结构

```
crates/trinity-butterfly/src/protein_extend.rs   # 新建: 密码子翻译 + 读框映射 + 引导走图
crates/trinity-butterfly/src/lib.rs              # +1 行 mod
crates/trinity-cli/src/protein_guided.rs         # 新建: DIAMOND 接线 + 模板/缺口 + stage 编排
crates/trinity-cli/src/lib.rs                    # +1 行 mod
crates/trinity-cli/src/args.rs                   # --protein-guided / --protein-db
crates/trinity-cli/src/orchestrate.rs            # butterfly 与 harvest 之间插 stage（1 个 if）
crates/trinity-cli/tests/protein_guided_e2e.rs   # 新建: e2e
crates/trinity-cli/tests/cli.rs                  # args 校验测试
docs/opt-log.md                                  # probe 入档
```

## 核心数学（全计划共用，先钉死）

HSP 坐标全部 1-based inclusive；`plus = qframe > 0`。绝对查询坐标允许 ≤0 / >qlen（延长侧越界），用 i64：

- **plus 框**: 氨基酸 j（0-based，自 sstart）↔ 查询区间 `[qs + 3j, qs + 3j + 2]`，顺读直译
- **minus 框**: 氨基酸 j ↔ 查询区间 `[qs − 3j − 2, qs − 3j]`，反读（revcomp 后直译）
- 延长 base 的绝对坐标: 右延长第 m 个 base = `qlen + m`；左延长第 m 个 = `1 − m`（即第 i 个到达序 base ↔ 位置 `1 − (i+1) + ... `）——统一按索引: **到达序 ext[i] ↔ 绝对位置 `qlen + 1 + i`（右）/ `−i`（左，ext[0] ↔ 0）**
- 密码子窗口对齐: plus 框窗口 `[lo, lo+2]` 要求 `(lo − qs) % 3 == 0`（rem_euclid）；aa = `sstart + (lo − qs)/3`。minus 框窗口 `[hi−2, hi]` 要求 `(qs − hi) % 3 == 0`；aa = `sstart + (qs − hi)/3`
- 图上一步 = 一个新 base（顶点名是 k-mer，相邻顶点重叠 k−1，右行取后继顶点名的末字符，左行取前驱顶点名的首字符）

---

### Task 1: 最小翻译器（trinity-butterfly/protein_extend.rs）

**Files:** Create `crates/trinity-butterfly/src/protein_extend.rs`; Modify `crates/trinity-butterfly/src/lib.rs`

- [x] **Step 1: 写失败测试**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codon_table_basics() {
        assert_eq!(translate_codon(b"ATG"), b'M');
        assert_eq!(translate_codon(b"TAA"), b'*');
        assert_eq!(translate_codon(b"TAG"), b'*');
        assert_eq!(translate_codon(b"TGA"), b'*');
        assert_eq!(translate_codon(b"GGC"), b'G');
        assert_eq!(translate_codon(b"TRY"), b'X'); // 歧义
    }

    #[test]
    fn revcomp_works() {
        assert_eq!(revcomp(b"ATGC"), b"GCAT");
    }
}
```

- [x] **Step 2: 跑测试确认失败**

Run: `cargo test -p trinity-butterfly protein_extend`
Expected: FAIL（函数未定义，编译错误）

- [x] **Step 3: 最小实现**

```rust
//! M6a 蛋白引导延长：HSP 读框感知的密码子翻译 + 原始 k-mer 图引导走图。

/// 标准密码子表（歧义 → X，stop → *）。
pub fn translate_codon(codon: &[u8]) -> u8 {
    let n = |b: u8| match b {
        b'T' => 0,
        b'C' => 1,
        b'A' => 2,
        b'G' => 3,
        _ => return b'X',
    };
    match codon.len() {
        3 => match (n(codon[0]), n(codon[1]), n(codon[2])) {
            (0, 0, 0) => b'F', (0, 0, 1) => b'L', (0, 0, 2) => b'Y', (0, 0, 3) => b'*',
            (0, 1, 0) => b'C', (0, 1, 3) => b'W',
            (0, 2, 0) => b'L', (0, 2, 1) => b'L', (0, 2, 2) => b'Q', (0, 2, 3) => b'R',
            (0, 3, 0) => b'L', (0, 3, 1) => b'L', (0, 3, 2) => b'R', (0, 3, 3) => b'R',
            (1, 0, 0) => b'I', (1, 0, 1) => b'I', (1, 0, 2) => b'I', (1, 0, 3) => b'M',
            (1, 1, 0) => b'T', (1, 1, 1) => b'T', (1, 1, 2) => b'T', (1, 1, 3) => b'T',
            (1, 2, 0) => b'N', (1, 2, 1) => b'K',
            (1, 3, 0) => b'S', (1, 3, 1) => b'S', (1, 3, 2) => b'R', (1, 3, 3) => b'S',
            (2, 0, 0) => b'I', (2, 0, 1) => b'M',
            (2, 1, 0) => b'T', (2, 1, 1) => b'T', (2, 1, 2) => b'T', (2, 1, 3) => b'T',
            (2, 2, 0) => b'N', (2, 2, 1) => b'K',
            (2, 3, 0) => b'S', (2, 3, 1) => b'R',
            (3, 0, 0) => b'V', (3, 0, 1) => b'V', (3, 0, 2) => b'V', (3, 0, 3) => b'V',
            (3, 1, 0) => b'A', (3, 1, 1) => b'A', (3, 1, 2) => b'A', (3, 1, 3) => b'A',
            (3, 2, 0) => b'D', (3, 2, 1) => b'E',
            (3, 3, 0) => b'G', (3, 3, 1) => b'G', (3, 3, 2) => b'G', (3, 3, 3) => b'G',
            _ => b'X',
        },
        _ => b'X',
    }
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}
```

`crates/trinity-butterfly/src/lib.rs` 加一行（与现有 mod 声明并列）：

```rust
pub mod protein_extend;
```

- [x] **Step 4: 跑测试确认通过**

Run: `cargo test -p trinity-butterfly protein_extend`
Expected: PASS (2 tests)

> **执行偏差（2026-08-24）:** 计划草稿里的密码子表本身有错（ATG 落 `_=>X`、TAA→Q 等）——
> 实现者按"测试即规格"换成标准遗传密码全表; 控制者复核时又抓到 TTA/TTG 误为 Y 并修正,
> 补了全 64 密码子对照 NCBI 表 1 的完备测试（codon_table_complete_standard_code）。
> 最终 3 测试全绿。

- [x] **Step 5: Commit（冻结期间跳过）**

```bash
git add crates/trinity-butterfly/src/protein_extend.rs crates/trinity-butterfly/src/lib.rs
git commit -m "feat(butterfly): M6a 密码子翻译器 + revcomp"
```

---

### Task 2: 读框坐标映射（同文件）

**Files:** Modify `crates/trinity-butterfly/src/protein_extend.rs`

- [x] **Step 1: 写失败测试**

在 `mod tests` 内追加：

```rust
    /// plus 框 HSP: qs=2, sstart=5。
    /// aa j=0 ↔ [2,4], j=1 ↔ [5,7]; 窗口 [2,4] → aa 5; [8,10] → aa 7。
    #[test]
    fn plus_frame_windows() {
        let hsp = AnchorHsp { qstart: 2, sstart: 5, plus: true };
        assert_eq!(hsp.aa_of_window(2, 4), Some(5));
        assert_eq!(hsp.aa_of_window(5, 7), Some(6));
        assert_eq!(hsp.aa_of_window(8, 10), Some(7));
        assert_eq!(hsp.aa_of_window(3, 5), None, "未对齐窗口");
        assert_eq!(hsp.aa_of_window(-1, 1), Some(4), "左侧越界 j=-1");
    }

    /// minus 框 HSP: qs=10, sstart=3。
    /// aa j=0 ↔ [8,10], j=1 ↔ [5,7]; 窗口 [11,13] → j=-1 → aa 2。
    #[test]
    fn minus_frame_windows() {
        let hsp = AnchorHsp { qstart: 10, sstart: 3, plus: false };
        assert_eq!(hsp.aa_of_window(8, 10), Some(3));
        assert_eq!(hsp.aa_of_window(5, 7), Some(4));
        assert_eq!(hsp.aa_of_window(11, 13), Some(2), "右侧越界 j=-1");
        assert_eq!(hsp.aa_of_window(7, 9), None, "未对齐窗口");
    }

    /// 翻译一个窗口: plus 顺读, minus 反读。
    #[test]
    fn window_translation_respects_strand() {
        let plus = AnchorHsp { qstart: 1, sstart: 1, plus: true };
        assert_eq!(plus.translate_window(b"ATG"), b'M');
        let minus = AnchorHsp { qstart: 10, sstart: 1, plus: false };
        // 反读: CAT → revcomp ATG → M
        assert_eq!(minus.translate_window(b"CAT"), b'M');
    }
```

- [x] **Step 2: 跑测试确认失败**

Run: `cargo test -p trinity-butterfly protein_extend`
Expected: FAIL（AnchorHsp 未定义）

- [x] **Step 3: 实现**

```rust
/// HSP 锚点（只留走图需要的三个量）。坐标 1-based inclusive。
#[derive(Debug, Clone, Copy)]
pub struct AnchorHsp {
    /// HSP 对齐起点的查询坐标（plus 框 = 低侧, minus 框 = 高侧）。
    pub qstart: i64,
    /// 与 qstart 对应的蛋白坐标（1-based）。
    pub sstart: i64,
    /// 读框正负（blastx qframe 符号）。
    pub plus: bool,
}

impl AnchorHsp {
    /// 密码子窗口 [lo, hi]（lo ≤ hi, 3 长）→ 蛋白 1-based 坐标; 未对齐返回 None。
    pub fn aa_of_window(&self, lo: i64, hi: i64) -> Option<usize> {
        let aa = if self.plus {
            if (lo - self.qstart).rem_euclid(3) != 0 {
                return None;
            }
            self.sstart + (lo - self.qstart) / 3
        } else {
            if (self.qstart - hi).rem_euclid(3) != 0 {
                return None;
            }
            self.sstart + (self.qstart - hi) / 3
        };
        aa.try_into().ok()
    }

    /// 窗口三碱基 → 氨基酸（plus 顺读, minus 反读）。
    pub fn translate_window(&self, w: &[u8]) -> u8 {
        if w.len() != 3 {
            return b'X';
        }
        if self.plus {
            translate_codon(w)
        } else {
            let rc = revcomp(w);
            translate_codon(&rc)
        }
    }

    /// 覆盖绝对坐标 frontier 的密码子窗口（plus: [lo,lo+2] 对齐 qs;
    /// minus: [hi-2,hi] 对齐 qs）。
    pub fn window_of_frontier(&self, frontier: i64) -> (i64, i64) {
        if self.plus {
            let lo = frontier - (frontier - self.qstart).rem_euclid(3);
            (lo, lo + 2)
        } else {
            let hi = frontier + (self.qstart - frontier).rem_euclid(3);
            (hi - 2, hi)
        }
    }
}
```

- [x] **Step 4: 跑测试确认通过**

Run: `cargo test -p trinity-butterfly protein_extend`
Expected: PASS (5 tests)

> **执行偏差（2026-08-24）:** 计划测试里 `window_of_frontier(5)` 期望 `(3,5)` 是错的
> （3..5 不在 minus 框对齐位上, `aa_of_window(3,5)`=None 自相矛盾）——实现者按实现公式
> 的正确语义改为 `(5,7)` 并注明。7 测试全绿。

- [x] **Step 5: Commit（冻结期间跳过）**

```bash
git add crates/trinity-butterfly/src/protein_extend.rs
git commit -m "feat(butterfly): M6a HSP 读框窗口映射"
```

---

### Task 3: 引导走图 extend_transcript（同文件）

**Files:** Modify `crates/trinity-butterfly/src/protein_extend.rs`（引用同 crate 的 `crate::graph::{DiGraph, SeqVertex, SimpleEdge}`，无需改 pub）

- [x] **Step 1: 写失败测试**

`mod tests` 追加：

```rust
    use crate::graph::{DiGraph, SeqVertex, SimpleEdge};

    /// k=6 分叉图。末端锚 "AAAAAC"(v2):
    ///   主干 v2→v3("AAAACT")→v4("AAACTG")→v5("AACTGG")  重边 10, 碱基 T,G,G
    ///   分支 v2→v9("AAAAAG")→v8("AAAAGT")→v10("AAAGTC") 轻边 2,  碱基 G,T,C
    /// 窗口 [7,9]: 主干 "TGG"→W, 分支 "GTC"→V。
    fn branch_graph() -> DiGraph {
        let mut g = DiGraph::new();
        for (id, name) in [
            (1, "AAAAAA"),
            (2, "AAAAAC"),
            (3, "AAAACT"),
            (4, "AAACTG"),
            (5, "AACTGG"),
            (9, "AAAAAG"),
            (8, "AAAAGT"),
            (10, "AAAGTC"),
        ] {
            g.add_vertex(SeqVertex::new(id, name));
        }
        g.add_edge(1, 2, SimpleEdge::new(10.0, 1, 2));
        g.add_edge(2, 3, SimpleEdge::new(10.0, 2, 3)); // 重
        g.add_edge(3, 4, SimpleEdge::new(10.0, 3, 4));
        g.add_edge(4, 5, SimpleEdge::new(10.0, 4, 5));
        g.add_edge(2, 9, SimpleEdge::new(2.0, 2, 9)); // 轻
        g.add_edge(9, 8, SimpleEdge::new(2.0, 9, 8));
        g.add_edge(8, 10, SimpleEdge::new(2.0, 8, 10));
        g
    }

    #[test]
    fn weight_only_walk_takes_heavy_branch() {
        let g = branch_graph();
        // 无引导（空模板）: 走重边 → 追加 "TGG"（走到 v5 尽头）
        let ext = extend_transcript(&g, 6, b"AAAAAC", Direction::Right,
            &AnchorHsp { qstart: 1, sstart: 1, plus: true }, b"", 10);
        assert_eq!(ext, b"TGG");
    }

    #[test]
    fn guided_walk_overrides_weight_on_codon_match() {
        let g = branch_graph();
        // qs=1 → 窗口 [7,9] = aa 3。模板 aa3='V': 分支 GTC→V 匹配(×3.0×2=6),
        // 主干 TGG→W 失配(×0.1×10=1) → 轻分支胜出。
        let template = b"MMVMM";
        let ext = extend_transcript(&g, 6, b"AAAAAC", Direction::Right,
            &AnchorHsp { qstart: 1, sstart: 1, plus: true }, template, 10);
        assert_eq!(ext, b"GTC", "guided ext = {:?}", String::from_utf8_lossy(&ext));
    }

    #[test]
    fn budget_caps_extension() {
        let g = branch_graph();
        let ext = extend_transcript(&g, 6, b"AAAAAC", Direction::Right,
            &AnchorHsp { qstart: 1, sstart: 1, plus: true }, b"", 2);
        assert_eq!(ext.len(), 2, "预算在密码子中途截断");
    }

    #[test]
    fn left_walk_uses_predecessors() {
        let g = branch_graph();
        // 从 "AAAAAC" 左延: 前驱 v1 "AAAAAA", 追加首字符 'A', v1 无前驱停。
        let ext = extend_transcript(&g, 6, b"AAAAAC", Direction::Left,
            &AnchorHsp { qstart: 6, sstart: 1, plus: true }, b"", 5);
        assert_eq!(ext, b"A");
    }

    #[test]
    fn anchor_kmer_missing_returns_empty() {
        let g = branch_graph();
        let ext = extend_transcript(&g, 6, b"TTTTTT", Direction::Right,
            &AnchorHsp { qstart: 1, sstart: 1, plus: true }, b"", 5);
        assert!(ext.is_empty());
    }
```

- [x] **Step 2: 跑测试确认失败**

Run: `cargo test -p trinity-butterfly protein_extend`
Expected: FAIL（extend_transcript/Direction 未定义）

- [x] **Step 3: 实现**

```rust
use crate::graph::DiGraph;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    Right,
    Left,
}

/// 延长一条转录本，返回**延长区字节**。
/// Right = 追加序（拼法 `t + ext`）; Left = 前插序（拼法 `ext + t`）。
///
/// 密码子窗口驱动的 3 步前瞻贪心：每次填满一个对齐窗口
/// （plus: [lo,lo+2]; minus: [hi-2,hi]），候选路径按
/// `min 边权 × (aa 匹配 3.0 / 失配 0.1 / 模板无信息 1.0)` 打分。
/// 终止：无路 / 预算尽 / 顶点已访问（环守卫）/ 窗口译出 '*' 且 aa 已到模板末位。
pub fn extend_transcript(
    graph: &DiGraph,
    kmer_size: usize,
    transcript: &[u8],
    dir: Direction,
    anchor: &AnchorHsp,
    template_aa: &[u8],
    budget_nt: usize,
) -> Vec<u8> {
    if transcript.len() < kmer_size || budget_nt == 0 {
        return Vec::new();
    }
    let qlen = transcript.len() as i64;
    let anchor_kmer: &[u8] = match dir {
        Direction::Right => &transcript[transcript.len() - kmer_size..],
        Direction::Left => &transcript[..kmer_size],
    };
    let mut by_name: HashMap<&str, i32> = HashMap::new();
    for &vid in graph.vertex_ids() {
        if let Some(v) = graph.get_vertex(vid) {
            by_name.insert(v.name.as_str(), vid);
        }
    }
    let Some(&cur) = by_name.get(std::str::from_utf8(anchor_kmer).unwrap_or("")) else {
        return Vec::new();
    };
    let mut cur = cur;

    // 绝对坐标 → 碱基。buf = 到达序延长区: 右 buf[i]↔qlen+1+i, 左 buf[i]↔−i。
    let base_at = |buf: &[u8], pos: i64| -> Option<u8> {
        let idx = match dir {
            Direction::Right => pos - qlen - 1,
            Direction::Left => -pos,
        };
        if idx >= 0 && (idx as usize) < buf.len() {
            Some(buf[idx as usize])
        } else if (1..=qlen).contains(&pos) {
            transcript[(pos - 1) as usize].to_ascii_upper()
        } else {
            None
        }
    };

    let mut ext: Vec<u8> = Vec::new(); // 到达序
    let mut visited: HashSet<i32> = HashSet::new();
    visited.insert(cur);
    let mut added = 0usize;
    while added < budget_nt {
        let frontier = match dir {
            // ext[i] ↔ qlen+1+i → 下一个未填位 = qlen + len + 1
            Direction::Right => qlen + ext.len() as i64 + 1,
            // ext[i] ↔ −i → 已填 0..=−(len−1), 下一个未填位 = −len
            Direction::Left => -ext.len() as i64,
        };
        let (lo, hi) = anchor.window_of_frontier(frontier);
        let need: Vec<i64> = (lo..=hi).filter(|&p| base_at(&ext, p).is_none()).collect();
        if need.is_empty() {
            break; // 防御: 理论不可达
        }
        let Some((bases, next_v)) = dfs_codon_paths(
            graph, cur, &need, dir, &ext, &base_at, anchor,
            template_aa, &visited, (lo, hi),
        ) else {
            break;
        };
        for b in bases {
            ext.push(b);
            added += 1;
            if added >= budget_nt {
                break;
            }
        }
        cur = next_v;
        visited.insert(cur);
        // 蛋白完整: 末位窗口译出 stop
        if let Some(aa) = anchor.aa_of_window(lo, hi) {
            if aa == template_aa.len() {
                let w: Vec<u8> = (lo..=hi).filter_map(|p| base_at(&ext, p)).collect();
                if w.len() == 3 && anchor.translate_window(&w) == b'*' {
                    break;
                }
            }
        }
    }
    match dir {
        Direction::Right => ext,
        // 到达序第一个最贴近原转录本, 前插后应为最靠右 → 逆序
        Direction::Left => ext.into_iter().rev().collect(),
    }
}

/// 从 cur 出发枚举覆盖 need（1..=3 个位置）的路径，返回得分最高者
/// `(到达序碱基, 结束顶点)`；无可行路径返回 None。
#[allow(clippy::too_many_arguments)]
fn dfs_codon_paths(
    graph: &DiGraph,
    cur: i32,
    need: &[i64],
    dir: Direction,
    ext: &[u8],
    base_at: &dyn Fn(&[u8], i64) -> Option<u8>,
    anchor: &AnchorHsp,
    template_aa: &[u8],
    visited: &HashSet<i32>,
    window: (i64, i64),
) -> Option<(Vec<u8>, i32)> {
    let mut best: Option<(f64, Vec<u8>, i32)> = None;
    // 迭代式小 DFS（深度 ≤ 3; 组件出度有限无爆炸）
    let mut stack: Vec<(i32, Vec<u8>, f64)> = vec![(cur, Vec::new(), f64::INFINITY)];
    while let Some((v, bases, min_w)) = stack.pop() {
        if bases.len() == need.len() {
            let mut all = ext.to_vec();
            all.extend_from_slice(&bases);
            let w: Vec<u8> = (window.0..=window.1)
                .filter_map(|p| base_at(&all, p))
                .collect();
            let score = if w.len() == 3 {
                match anchor.aa_of_window(window.0, window.1) {
                    Some(aa) if (1..=template_aa.len()).contains(&aa) => {
                        if anchor.translate_window(&w) == template_aa[aa - 1] {
                            min_w * 3.0
                        } else {
                            min_w * 0.1
                        }
                    }
                    _ => min_w,
                }
            } else {
                min_w
            };
            let better = match &best {
                Some((s, _, _)) => score > *s,
                None => true,
            };
            if better {
                best = Some((score, bases, v));
            }
            continue;
        }
        let succs: Vec<i32> = match dir {
            Direction::Right => graph.get_successors(v).to_vec(),
            Direction::Left => graph.get_predecessors(v).to_vec(),
        };
        // 逆序压栈 → 弹出按顶点 id 升序（平局确定性）
        for nv in succs.into_iter().rev() {
            if visited.contains(&nv) {
                continue;
            }
            let Some(vv) = graph.get_vertex(nv) else { continue };
            let name = vv.name.as_bytes();
            if name.is_empty() {
                continue;
            }
            let b = match dir {
                Direction::Right => name[name.len() - 1].to_ascii_upper(),
                Direction::Left => name[0].to_ascii_upper(),
            };
            let w = graph.find_edge(v, nv).map(|e| e.weight).unwrap_or(0.0).max(0.0);
            let mut nb = bases.clone();
            nb.push(b);
            stack.push((nv, nb, min_w.min(w)));
        }
    }
    best.map(|(_, bases, v)| (bases, v))
}
```

- [x] **Step 4: 跑测试确认通过**

Run: `cargo test -p trinity-butterfly protein_extend`
Expected: PASS (10 tests)

调试提示: `guided_walk_overrides_weight_on_codon_match` 若失败——核对窗口对齐（6nt 转录本 + qs=1 → 窗口 [7,9]，`(7−1)%3==0` ✓）、主干/分支碱基提取（右行取名字末字符）、评分（匹配 3.0/失配 0.1 乘 min 边权）。

> **执行记录（2026-08-24）:** 质量审查后修 3 项——锚 k-mer 大写化、
> `dfs_codon_paths` 加 kmer_size 参数并跳过名字长度≠k 的后继、补 3 个测试
> （左延多碱基次序 / minus 框引导 / stop 终止; 左延测试的图设置经修复者
> 改为独立 4 前驱链 + budget 3, 断言 b"CGT" 受保护未动）。复审通过, 15 测试全绿。
> 顺带修 `-ext.len() as i64` 的 Rust 解析问题（一元负号优先于 as）。

- [x] **Step 5: clippy**

Run: `cargo clippy -p trinity-butterfly --tests`
Expected: 无 warning

- [x] **Step 6: Commit（冻结期间跳过）**

```bash
git add crates/trinity-butterfly/src/protein_extend.rs
git commit -m "feat(butterfly): M6a 引导走图 extend_transcript"
```

---

### Task 4: DIAMOND 接线 + 解析（trinity-cli/protein_guided.rs）

**Files:** Create `crates/trinity-cli/src/protein_guided.rs`; Modify `crates/trinity-cli/src/lib.rs`（`pub mod protein_guided;`）

- [x] **Step 1: 写失败测试**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_tab_roundtrip() {
        let text = "t1\tP1\t90.1\t100\t300\t200\t1\t300\t1\t100\t1e-40\t185.2\t1\n\
                    t1\tP2\t50.0\t50\t300\t80\t10\t160\t5\t55\t1e-9\t45.0\t-1\n";
        let hits = parse_diamond_tab(text);
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].qid, "t1");
        assert_eq!(hits[0].sid, "P1");
        assert_eq!(hits[0].qlen, 300);
        assert_eq!(hits[0].slen, 200);
        assert_eq!(hits[0].qframe, 1);
        assert_eq!(hits[1].qframe, -1);
        assert!((hits[0].bitscore - 185.2).abs() < 1e-9);
    }

    #[test]
    fn best_hit_per_query_max_bitscore() {
        let text = "a\tP1\t90\t100\t300\t200\t1\t300\t1\t100\t1e-40\t185.2\t1\n\
                    a\tP2\t99\t300\t300\t310\t1\t300\t1\t300\t1e-99\t520.0\t1\n\
                    b\tP3\t88\t90\t100\t120\t5\t95\t8\t98\t1e-20\t95.0\t1\n";
        let hits = parse_diamond_tab(text);
        let best = best_hit_per_query(&hits);
        assert_eq!(best["a"].sid, "P2");
        assert_eq!(best["b"].sid, "P3");
    }
}
```

- [x] **Step 2: 跑测试确认失败**

Run: `cargo test -p trinity-cli protein_guided`
Expected: FAIL（模块未建）

- [x] **Step 3: 实现**

```rust
//! M6a 蛋白引导组装（探路版）: DIAMOND blastx 模板选择 + 缺口轮廓 +
//! 引导延长 stage（插在 butterfly 池与 harvest 之间）。

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use trinity_common::error::CommonError;

/// DIAMOND 可执行文件（env `TRINITY_DIAMOND` 可覆盖）。
fn diamond_bin() -> &'static str {
    static DEFAULT: &str = "/public/home/senior007/miniconda3/envs/diamond2go/bin/diamond";
    match std::env::var("TRINITY_DIAMOND") {
        Ok(v) if !v.is_empty() => Box::leak(v.into_boxed_str()),
        _ => DEFAULT,
    }
}

pub struct DiamondHit {
    pub qid: String,
    pub sid: String,
    pub pident: f64,
    pub qlen: i64,
    pub slen: i64,
    pub qstart: i64,
    pub qend: i64,
    pub sstart: i64,
    pub send: i64,
    pub evalue: f64,
    pub bitscore: f64,
    pub qframe: i32,
}

/// outfmt 6: qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qframe
pub fn parse_diamond_tab(text: &str) -> Vec<DiamondHit> {
    let mut out = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 13 {
            continue;
        }
        let (Ok(qlen), Ok(slen)) = (f[4].parse::<i64>(), f[5].parse::<i64>()) else {
            continue;
        };
        let (Ok(qs), Ok(qe), Ok(ss), Ok(se)) = (
            f[6].parse::<i64>(),
            f[7].parse::<i64>(),
            f[8].parse::<i64>(),
            f[9].parse::<i64>(),
        ) else {
            continue;
        };
        out.push(DiamondHit {
            qid: f[0].to_string(),
            sid: f[1].to_string(),
            pident: f[2].parse().unwrap_or(0.0),
            qlen,
            slen,
            qstart: qs,
            qend: qe,
            sstart: ss,
            send: se,
            evalue: f[10].parse().unwrap_or(f64::MAX),
            bitscore: f[11].parse().unwrap_or(0.0),
            qframe: f[12].parse().unwrap_or(0),
        });
    }
    out
}

/// 每查询最佳命中（bitscore 最高; 平局取先出现者）。
pub fn best_hit_per_query(hits: &[DiamondHit]) -> HashMap<String, &DiamondHit> {
    let mut best: HashMap<String, &DiamondHit> = HashMap::new();
    for h in hits {
        match best.get(&h.qid) {
            Some(b) if b.bitscore >= h.bitscore => {}
            _ => {
                best.insert(h.qid.clone(), h);
            }
        }
    }
    best
}

/// makedb（存在即跳过）→ 返回 dmnd 路径。
pub fn diamond_makedb_if_needed(db_fa: &Path, workdir: &Path) -> Result<PathBuf, CommonError> {
    fs::create_dir_all(workdir)?;
    let dmnd = workdir.join("protein.dmnd");
    if dmnd.exists() {
        return Ok(dmnd);
    }
    let out = std::process::Command::new(diamond_bin())
        .arg("makedb")
        .arg("--in")
        .arg(db_fa)
        .arg("-d")
        .arg(&dmnd)
        .stderr(std::process::Stdio::null())
        .output()
        .map_err(|e| {
            CommonError::Parse(format!("cannot spawn diamond ({}): {e}", diamond_bin()))
        })?;
    if !out.status.success() {
        return Err(CommonError::Parse("diamond makedb failed".to_string()));
    }
    Ok(dmnd)
}

/// blastx（DNA 查询 vs 蛋白库）→ tabular 全文（临时输出文件读回后即留盘）。
pub fn diamond_blastx(query_fa: &Path, dmnd: &Path, threads: usize, tag: &str) -> Result<String, CommonError> {
    let out_path = std::env::temp_dir().join(format!("diamond_{tag}_{}", std::process::id()));
    let out = std::process::Command::new(diamond_bin())
        .arg("blastx")
        .arg("-q")
        .arg(query_fa)
        .arg("-d")
        .arg(dmnd)
        .arg("--threads")
        .arg(threads.to_string())
        .arg("-e")
        .arg("1e-5")
        .arg("--max-target-seqs")
        .arg("5")
        .arg("--outfmt")
        .arg("6")
        .arg("-o")
        .arg(&out_path)
        .stderr(std::process::Stdio::null())
        .output()
        .map_err(|e| {
            CommonError::Parse(format!("cannot spawn diamond ({}): {e}", diamond_bin()))
        })?;
    if !out.status.success() {
        return Err(CommonError::Parse("diamond blastx failed".to_string()));
    }
    fs::read_to_string(&out_path).map_err(|e| {
        CommonError::Parse(format!(
            "cannot read diamond output {}: {e}",
            out_path.display()
        ))
    })
}
```

`crates/trinity-cli/src/lib.rs` 加：`pub mod protein_guided;`

- [x] **Step 4: 跑测试确认通过**

Run: `cargo test -p trinity-cli protein_guided`
Expected: PASS (2 tests)

- [x] **Step 5: Commit（冻结期间跳过）**

```bash
git add crates/trinity-cli/src/protein_guided.rs crates/trinity-cli/src/lib.rs
git commit -m "feat(cli): M6a DIAMOND 接线 + tabular 解析"
```

---

### Task 5: 蛋白库索引 + 缺口轮廓（同文件）

**Files:** Modify `crates/trinity-cli/src/protein_guided.rs`

- [x] **Step 1: 写失败测试**

`mod tests` 追加：

```rust
    #[test]
    fn protein_db_index() {
        let dir = std::env::temp_dir().join(format!("pg_idx_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("p.fa");
        fs::write(&fa, ">P1\nMKVL\n>P2\nAAA\n").unwrap();
        let idx = index_protein_db(&fa).unwrap();
        assert_eq!(idx["P1"], b"MKVL");
        assert_eq!(idx.len(), 2);
    }

    /// 转录本 300nt 命中 P1(200aa) 的 aa 60..160（plus 框）:
    /// N 缺 59 > C 缺 39 → 补 N → plus 框向左; budget = 59 × 3 × 1.5。
    #[test]
    fn gap_profile_missing_nterm() {
        let h = DiamondHit {
            qid: "t".into(), sid: "P1".into(), pident: 90.0,
            qlen: 300, slen: 200, qstart: 180, qend: 300,
            sstart: 60, send: 160, evalue: 1e-40, bitscore: 200.0, qframe: 1,
        };
        let gp = gap_profile(&h);
        assert!(gp.plus);
        assert!(!gp.dir_right, "缺 N 端更大 → plus 框向左");
        assert!(gp.missing_nterm());
        assert!(!gp.missing_cterm());
        assert_eq!(gp.budget_nt(), 59 * 3 * 3 / 2);
    }

    /// 只缺 C 端（覆盖 aa 1..160）→ 向右。
    #[test]
    fn gap_profile_missing_cterm_only() {
        let h = DiamondHit {
            qid: "t".into(), sid: "P1".into(), pident: 90.0,
            qlen: 483, slen: 200, qstart: 4, qend: 483,
            sstart: 1, send: 160, evalue: 1e-40, bitscore: 200.0, qframe: 1,
        };
        let gp = gap_profile(&h);
        assert!(gp.dir_right);
        assert_eq!(gp.budget_nt(), 39 * 3 * 3 / 2);
    }

    /// minus 框方向翻转: 只缺 C 端 → 向左。
    #[test]
    fn gap_profile_minus_directions() {
        let h = DiamondHit {
            qid: "t".into(), sid: "P1".into(), pident: 90.0,
            qlen: 300, slen: 200, qstart: 480, qend: 1,
            sstart: 160, send: 1, evalue: 1e-40, bitscore: 200.0, qframe: -1,
        };
        let gp = gap_profile(&h);
        assert!(!gp.plus);
        assert!(!gp.dir_right, "minus + 只缺 C → 向左");
    }
```

- [x] **Step 2: 跑测试确认失败**

Run: `cargo test -p trinity-cli protein_guided`
Expected: FAIL

- [x] **Step 3: 实现**

```rust
/// 蛋白库 fasta → sseqid → 氨基酸序列（header 取首空白 token）。
pub fn index_protein_db(path: &Path) -> Result<HashMap<String, Vec<u8>>, CommonError> {
    let text = fs::read_to_string(path).map_err(|e| {
        CommonError::Parse(format!("cannot read protein db {}: {e}", path.display()))
    })?;
    let mut out = HashMap::new();
    let mut cur: Option<(String, Vec<u8>)> = None;
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if let Some((k, v)) = cur.take() {
                out.insert(k, v);
            }
            cur = Some((
                h.split_whitespace().next().unwrap_or_default().to_string(),
                Vec::new(),
            ));
        } else if let Some((_, seq)) = cur.as_mut() {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if let Some((k, v)) = cur.take() {
        out.insert(k, v);
    }
    Ok(out)
}

/// 缺口轮廓: 延长方向 + 预算（单 HSP 版）。
#[derive(Debug, Clone, Copy)]
pub struct GapProfile {
    pub plus: bool,
    /// true = 向转录本 3' 延长。
    pub dir_right: bool,
    n_missing: usize,
    c_missing: usize,
    budget_nt: usize,
}

impl GapProfile {
    pub fn missing_nterm(&self) -> bool {
        self.n_missing > 0
    }
    pub fn missing_cterm(&self) -> bool {
        self.c_missing > 0
    }
    pub fn budget_nt(&self) -> usize {
        self.budget_nt
    }
}

/// 单 HSP: 覆盖区间 [min(sstart,send), max(sstart,send)]。
/// 补缺失跨度大的一侧（平局取 C 侧）; plus: C→Right / N→Left; minus 反转。
/// budget = 该侧缺失 aa × 3 × 1.5, 上限 5000。
pub fn gap_profile(h: &DiamondHit) -> GapProfile {
    let (lo, hi) = (h.sstart.min(h.send), h.sstart.max(h.send));
    let n_missing = (lo - 1).max(0) as usize;
    let c_missing = (h.slen - hi).max(0) as usize;
    let plus = h.qframe > 0;
    let fill_c = c_missing >= n_missing;
    let (dir_right, missing) = if fill_c {
        (plus, c_missing)
    } else {
        (!plus, n_missing)
    };
    GapProfile {
        plus,
        dir_right,
        n_missing,
        c_missing,
        budget_nt: (missing * 3 * 3 / 2).min(5000),
    }
}
```

- [x] **Step 4: 跑测试确认通过**

Run: `cargo test -p trinity-cli protein_guided`
Expected: PASS (6 tests)

- [x] **Step 5: Commit（冻结期间跳过）**

```bash
git add crates/trinity-cli/src/protein_guided.rs
git commit -m "feat(cli): M6a 蛋白库索引 + 缺口轮廓"
```

---

### Task 6: stage 编排 + args + orchestrate 接线

**Files:** Modify `crates/trinity-cli/src/protein_guided.rs`, `args.rs`, `orchestrate.rs`, `tests/cli.rs`

- [x] **Step 1: 写失败测试（args 校验, tests/cli.rs 追加）**

```rust
fn pg_base_argv(r1: &str, r2: &str, out: &str) -> Vec<String> {
    ["--seqType", "fq", "--left", r1, "--right", r2, "--max_memory", "1G",
     "--output", out, "--CPU", "2"]
        .iter().map(|s| s.to_string()).collect()
}

#[test]
fn protein_guided_requires_protein_db() {
    let mut argv = pg_base_argv("/tmp/a.fq", "/tmp/b.fq", "/tmp/out");
    argv.push("--protein-guided".into());
    let err = parse_args(&argv).unwrap_err();
    assert!(err.msg.contains("--protein-guided requires --protein-db"));
}

#[test]
fn protein_guided_db_must_exist() {
    let mut argv = pg_base_argv("/tmp/a.fq", "/tmp/b.fq", "/tmp/out");
    argv.push("--protein-guided".into());
    argv.push("--protein-db".into());
    argv.push("/nonexistent/P.fa".into());
    let err = parse_args(&argv).unwrap_err();
    assert!(err.msg.contains("not found"));
}
```

（`tests/cli.rs` 已 `use trinity_cli::args::parse_args` 的话直接用; read 路径文件不存在没关系——校验只查 `--protein-db` 自身存在性。）

- [x] **Step 2: 跑测试确认失败**

Run: `cargo test -p trinity-cli --test cli protein_guided`
Expected: FAIL (0 passed)

- [x] **Step 3: args.rs 三处修改**

结构体（`pub model: Option<PathBuf>,` 之后）：

```rust
    /// `--protein-guided`（M6a）: DIAMOND 蛋白模板引导延长。
    pub protein_guided: bool,
    /// `--protein-db <fa>`: 蛋白库（--protein-guided 必需）。
    pub protein_db: Option<PathBuf>,
```

解析（`let model = ...` 之后）：

```rust
    let protein_guided = cli.bool_flag("--protein-guided");
    let protein_db = cli.value_flag("--protein-db")?.map(PathBuf::from);
```

校验（`--model` 校验块之后）：

```rust
    if protein_guided {
        let Some(db) = &protein_db else {
            return Err(CliError::usage(
                "Error, --protein-guided requires --protein-db <fa>",
            ));
        };
        if !db.exists() {
            return Err(CliError::usage(format!(
                "Error, --protein-db file not found: {}",
                db.display()
            )));
        }
    }
```

构造器 `Ok(TrinityArgs { ... })` 末尾加 `protein_guided, protein_db,`。

- [x] **Step 4: 跑 args 测试确认通过**

Run: `cargo test -p trinity-cli --test cli protein_guided`
Expected: PASS (2 tests)

- [x] **Step 5: 实现 stage 主体（protein_guided.rs 追加）**

```rust
use trinity_butterfly::graph_io::build_new_graph_use_kmers;
use trinity_butterfly::protein_extend::{extend_transcript, AnchorHsp, Direction};

/// 延长后验收门槛: 新 bitscore 须超旧值 × (1 + 该系数)。
const ACCEPT_MARGIN: f64 = 0.05;

/// stage 入口: 插在 butterfly 池之后、harvest 之前。
/// 逐组件改写 `<base>.graph.allProbPaths.fasta`（原件备份 `.orig`, 一次性）。
pub fn run_protein_guided(
    listing: &[(u64, PathBuf)],
    args: &crate::args::TrinityArgs,
    outdir: &Path,
) -> Result<(), CommonError> {
    use crate::checkpoint::checkpoint_exists;
    let ok = outdir.join(".protein_guided.ok");
    if checkpoint_exists(&ok) {
        eprintln!(
            "---- Trinity (rust) checkpoint found, skipping: {}",
            ok.display()
        );
        return Ok(());
    }
    let Some(db_fa) = args.protein_db.clone() else {
        return Err(CommonError::Parse(
            "--protein-guided requires --protein-db".to_string(),
        ));
    };
    let work = outdir.join("protein_guided");
    let dmnd = diamond_makedb_if_needed(&db_fa, &work)?;
    let threads = args.cpu.max(1);

    // ---- 1. 收集全部初步转录本 → queries.fa ----
    let mut queries: Vec<(String, PathBuf, String, String)> = Vec::new(); // (qid, base, header, seq)
    let mut qfa = String::new();
    for (i, (_id, base)) in listing.iter().enumerate() {
        let fa = PathBuf::from(format!("{}.graph.allProbPaths.fasta", base.display()));
        let Ok(text) = fs::read_to_string(&fa) else { continue };
        for (header, seq) in parse_fasta_pairs(&text) {
            let acc = header.split_whitespace().next().unwrap_or("t");
            let qid = format!("c{i}_{acc}");
            qfa.push_str(&format!(">{qid}\n{seq}\n"));
            queries.push((qid, base.clone(), header, seq));
        }
    }
    let qfa_path = work.join("queries.fa");
    fs::write(&qfa_path, &qfa)?;

    // ---- 2. blastx + 模板 ----
    let hits = parse_diamond_tab(&diamond_blastx(&qfa_path, &dmnd, threads, "q")?);
    let best = best_hit_per_query(&hits);
    let proteins = index_protein_db(&db_fa)?;

    // ---- 3. 逐转录本延长（有缺口才动）----
    struct ExtOut {
        qid: String,
        base: PathBuf,
        header: String,
        new_seq: String,
        old_bitscore: f64,
    }
    let mut extended: Vec<ExtOut> = Vec::new();
    let mut with_t = 0usize;
    for (qid, base, header, seq) in &queries {
        let Some(h) = best.get(qid) else { continue };
        with_t += 1;
        let gp = gap_profile(h);
        if gp.budget_nt() == 0 {
            continue;
        }
        let gout = PathBuf::from(format!("{}.graph.out", base.display()));
        let Ok(gtext) = fs::read_to_string(&gout) else {
            eprintln!("-protein-guided: cannot reload {}", gout.display());
            continue;
        };
        let Ok(br) = build_new_graph_use_kmers(&gtext) else {
            eprintln!("-protein-guided: cannot parse {}", gout.display());
            continue;
        };
        let anchor = AnchorHsp {
            qstart: if h.qframe > 0 {
                h.qstart.min(h.qend)
            } else {
                h.qstart.max(h.qend)
            },
            sstart: h.sstart,
            plus: h.qframe > 0,
        };
        let dir = if gp.dir_right {
            Direction::Right
        } else {
            Direction::Left
        };
        let template = proteins.get(&h.sid).cloned().unwrap_or_default();
        let ext = extend_transcript(
            &br.graph,
            args.kmer_size,
            seq.as_bytes(),
            dir,
            &anchor,
            &template,
            gp.budget_nt(),
        );
        if ext.is_empty() {
            continue;
        }
        let new_seq = match dir {
            Direction::Right => format!("{seq}{}", String::from_utf8_lossy(&ext)),
            Direction::Left => format!("{}{seq}", String::from_utf8_lossy(&ext)),
        };
        extended.push(ExtOut {
            qid: qid.clone(),
            base: base.clone(),
            header: header.clone(),
            new_seq,
            old_bitscore: h.bitscore,
        });
    }

    // ---- 4. 验收: 延长产物重 blastx ----
    let accepted: Vec<&ExtOut> = if extended.is_empty() {
        Vec::new()
    } else {
        let mut xfa = String::new();
        for (i, e) in extended.iter().enumerate() {
            xfa.push_str(&format!(">x{i}\n{}\n", e.new_seq));
        }
        let xfa_path = work.join("extended.fa");
        fs::write(&xfa_path, &xfa)?;
        let xb = best_hit_per_query(&parse_diamond_tab(&diamond_blastx(
            &xfa_path, &dmnd, threads, "x",
        )?));
        extended
            .iter()
            .enumerate()
            .filter(|(i, e)| {
                let new = xb.get(&format!("x{i}")).map(|h| h.bitscore).unwrap_or(0.0);
                new > e.old_bitscore * (1.0 + ACCEPT_MARGIN)
            })
            .map(|(_, e)| e)
            .collect()
    };

    // ---- 5. 改写组件 fasta（原件备份一次）----
    let mut by_comp: HashMap<PathBuf, Vec<&ExtOut>> = HashMap::new();
    for e in &accepted {
        by_comp.entry(e.base.clone()).or_default().push(e);
    }
    let mut total_nt = 0usize;
    for (base, exts) in &by_comp {
        let fa = PathBuf::from(format!("{}.graph.allProbPaths.fasta", base.display()));
        let backup = PathBuf::from(format!("{}.graph.allProbPaths.orig", base.display()));
        let text = fs::read_to_string(&fa)?;
        if !backup.exists() {
            fs::write(&backup, &text)?;
        }
        let mut out = String::new();
        for (header, seq) in parse_fasta_pairs(&text) {
            match exts.iter().find(|e| e.header == header) {
                Some(e) => {
                    total_nt += e.new_seq.len().saturating_sub(seq.len());
                    out.push_str(&format!(">{header}\n{}\n", e.new_seq));
                }
                None => out.push_str(&format!(">{header}\n{seq}\n")),
            }
        }
        fs::write(&fa, &out)?;
    }
    eprintln!(
        "-- protein-guided: {with_t} transcripts with template, {} extended ({} accepted, {} reverted), +{total_nt} nt",
        extended.len(),
        accepted.len(),
        extended.len() - accepted.len()
    );
    fs::write(&ok, b"")?;
    Ok(())
}

fn parse_fasta_pairs(text: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut header: Option<String> = None;
    let mut seq = String::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if let Some(h0) = header.take() {
                out.push((h0, std::mem::take(&mut seq)));
            }
            header = Some(h.to_string());
        } else if header.is_some() {
            seq.push_str(line.trim());
        }
    }
    if let Some(h0) = header {
        out.push((h0, seq));
    }
    out
}
```

- [x] **Step 6: orchestrate.rs 接线**

在 butterfly checkpoint 块结束（`run_with_checkpoint(... run_butterfly_pool ...)  ?;` 之后）与 `// ---- 7. 汇总` 之间插入：

```rust
    // ---- 6b. 蛋白引导延长（M6a; checkpoint .protein_guided.ok）----
    if args.protein_guided {
        crate::protein_guided::run_protein_guided(&listing, args, &outdir)?;
    }
```

- [x] **Step 7: 编译 + 全量测试 + lint**

Run: `cargo test -p trinity-cli && cargo test -p trinity-butterfly`
Expected: 全 PASS（旗标 off 不进新分支, 现有测试零回归）

Run: `cargo clippy -p trinity-cli -p trinity-butterfly --tests && cargo fmt -p trinity-cli -p trinity-butterfly`
Expected: 干净

> **执行记录（2026-08-24）:** 双审通过后按质量审查修 7 项——组件 fasta 原子改写
> （.pg_tmp + rename）、改写/验收逻辑抽纯函数并补 4 单测、按组件分组使
> `build_new_graph_use_kmers` 每组件只解析一次、延长字节 from_utf8 显式检查、
> DIAMOND stderr 并入错误信息、临时输出清理、空 queries 早退。
> 复审通过, protein_guided 10 单测全绿。
> 计划测试草稿的 c_missing 差一（inclusive 坐标 200−160=40 非 39）在 Task 5 已修。

- [x] **Step 8: Commit（冻结期间跳过）**

```bash
git add crates/trinity-cli
git commit -m "feat(cli): M6a --protein-guided stage 接线（butterfly→harvest 之间）"
```

---

### Task 7: e2e（protein_guided_e2e.rs）

**Files:** Create `crates/trinity-cli/tests/protein_guided_e2e.rs`

- [x] **Step 1: 写测试**

```rust
//! M6a e2e: `--protein-guided`——合成 CDS 基因中段低覆盖 → vanilla 截断,
//! 蛋白库含真蛋白 → 引导延长恢复全长; 旗标关零变化。

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use trinity_cli::args::parse_args;
use trinity_cli::orchestrate::run_trinity;

fn tmpdir(tag: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!("trinity-cli-pg-{}-{}", tag, std::process::id()));
    let _ = fs::remove_dir_all(&p);
    fs::create_dir_all(&p).unwrap();
    p
}

fn xorshift_nt(x: &mut u64) -> usize {
    *x ^= *x << 13;
    *x ^= *x >> 7;
    *x ^= *x << 17;
    (*x % 20) as usize
}

fn revcomp(s: &str) -> String {
    s.bytes().rev().map(|b| match b {
        b'A' => b'T', b'C' => b'G', b'G' => b'C', _ => b'A',
    }).collect()
}

/// 密码子表（每 aa 一码, 避开 stop）。
const CODONS: [&str; 20] = [
    "GCT", "TGT", "GAT", "GAA", "TTC", "GGC", "CAC", "ATC", "AAA", "CTG",
    "ATG", "AAC", "CCC", "CAA", "CGT", "AGC", "ACC", "GTT", "TGG", "TAC",
];
const AAS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

/// 合成 CDS: ATG + n 个非 stop 密码子 + TAA（3n+3 nt, 无内部 stop）。
fn synth_cds(seed: u64, n_codons: usize) -> (String, String) {
    let mut x = seed | 1;
    let mut dna = String::from("ATG");
    let mut aa = String::from("M");
    for _ in 0..n_codons {
        let i = xorshift_nt(&mut x);
        dna.push_str(CODONS[i]);
        aa.push(AAS[i] as char);
    }
    dna.push_str("TAA");
    aa.push('*');
    (dna, aa)
}

/// 平铺 PE reads（150nt, 插入 250, 步进 8）; middle 窗口 [lo,hi) 只放 sparse_pairs 对。
fn write_reads(
    dir: &Path, cds: &str, mid_lo: usize, mid_hi: usize, sparse_pairs: usize,
) -> (PathBuf, PathBuf) {
    let b = cds.as_bytes();
    let mut r1 = String::new();
    let mut r2 = String::new();
    let mut n = 0usize;
    let mut s = 0usize;
    while s + 250 <= b.len() {
        if s < mid_lo || s + 250 > mid_hi {
            let f1 = std::str::from_utf8(&b[s..s + 150]).unwrap();
            let f2 = std::str::from_utf8(&b[s + 100..s + 250]).unwrap();
            r1.push_str(&format!("@r{n}\n{f1}\n+\n{}\n", "I".repeat(150)));
            r2.push_str(&format!("@r{n}\n{}\n+\n{}\n", revcomp(f2), "I".repeat(150)));
            n += 1;
            s += 8;
        } else {
            s += 4; // 跳过中段
        }
    }
    for i in 0..sparse_pairs {
        let step = (mid_hi - mid_lo - 250).max(1) / sparse_pairs.max(1);
        let s = mid_lo + i * step;
        let f1 = std::str::from_utf8(&b[s..s + 150]).unwrap();
        let f2 = std::str::from_utf8(&b[s + 100..s + 250]).unwrap();
        r1.push_str(&format!("@m{i}\n{f1}\n+\n{}\n", "I".repeat(150)));
        r2.push_str(&format!("@m{i}\n{}\n+\n{}\n", revcomp(f2), "I".repeat(150)));
    }
    let p1 = dir.join("r1.fq.gz");
    let p2 = dir.join("r2.fq.gz");
    for (p, body) in [(&p1, &r1), (&p2, &r2)] {
        let mut w = flate2::write::GzEncoder::new(
            fs::File::create(p).unwrap(), flate2::Compression::fast());
        w.write_all(body.as_bytes()).unwrap();
        w.finish().unwrap();
    }
    (p1, p2)
}

fn base_args(r1: &Path, r2: &Path, outdir: &Path) -> Vec<String> {
    ["--seqType", "fq", "--left", r1.to_str().unwrap(), "--right",
     r2.to_str().unwrap(), "--max_memory", "1G", "--output",
     outdir.to_str().unwrap(), "--no_normalize_reads", "--CPU", "2",
     "--min_contig_length", "100", "--path-engine", "pogflow"]
        .iter().map(|s| s.to_string()).collect()
}

fn diamond_available() -> bool {
    let bin = std::env::var("TRINITY_DIAMOND").ok()
        .filter(|v| !v.is_empty())
        .unwrap_or_else(|| {
            "/public/home/senior007/miniconda3/envs/diamond2go/bin/diamond".to_string()
        });
    std::process::Command::new(bin).arg("--version")
        .stdout(std::process::Stdio::null()).stderr(std::process::Stdio::null())
        .status().map(|s| s.success()).unwrap_or(false)
}

fn max_len(fa: &Path) -> usize {
    fs::read_to_string(fa).unwrap()
        .lines().skip(1).fold(0usize, |a, l| a.max(l.trim().len()))
}

#[test]
fn protein_guided_e2e_extends_low_coverage_gap() {
    if !diamond_available() {
        eprintln!("skip: diamond not available");
        return;
    }
    let dir = tmpdir("e2e");
    // 900nt CDS（299 密码子+起止）; 中段 [300,600) 只给 1 对 reads
    let (cds, aa) = synth_cds(0xA5A5_1234_DEAD_BEEF, 299);
    assert_eq!(cds.len(), 900);
    let (r1, r2) = write_reads(&dir, &cds, 300, 600, 1);

    let pdb = dir.join("prot.fa");
    fs::write(&pdb, format!(">REAL\n{aa}\n")).unwrap();

    // ---- off: vanilla ----
    let out_off = dir.join("trinity_off");
    let f_off = run_trinity(&parse_args(&base_args(&r1, &r2, &out_off)).unwrap()).unwrap();
    let off_max = max_len(&f_off);

    // ---- on ----
    let out_on = dir.join("trinity_on");
    let mut argv = base_args(&r1, &r2, &out_on);
    argv.push("--protein-guided".into());
    argv.push("--protein-db".into());
    argv.push(pdb.display().to_string());
    let f_on = run_trinity(&parse_args(&argv).unwrap()).unwrap();
    let on_max = max_len(&f_on);

    assert!(out_on.join(".protein_guided.ok").exists());
    assert!(!out_off.join(".protein_guided.ok").exists());
    assert!(!out_off.join("protein_guided").exists());
    assert!(on_max >= off_max, "延长不应变短: {on_max} < {off_max}");
    // 主断言: 引导后接近全长（允许 10% 缺口）
    assert!(on_max >= 810, "on_max={on_max} off_max={off_max}");
}
```

- [x] **Step 2: 跑测试**

Run: `cargo test -p trinity-cli --test protein_guided_e2e`
Expected: PASS。

> **执行记录（2026-08-24）:** e2e 揭出真 bug——`--outfmt 6` 不带字段列表时 DIAMOND
> 输出默认 12 列（无 qlen/slen/qframe），解析器要求 13 列 → 全部丢弃、0 模板命中。
> 已修（字段逐个传参, 单测锁定 12 列被拒）。
> **结构性发现**: 单基因合成数据上管线在 k-mer 层全耦合——稀疏 read 与侧翼共享
> k-mer 才进组件图, 而共享即链通, butterfly 必产出全长; 不共享则图也无路可走。
> 即"图有路 ⟺ vanilla 全长", e2e 无法构造"图有路但枚举丢弃"场景（该场景需真实
> 数据的分支/重复复杂度）。最终 e2e 定位为**接线级**: DIAMOND 真实调用 + 产物
> （queries.fa/protein.dmnd）+ 旗标语义 + on≥off≥810 + off 零残留; 机制级正确性
> 由 Task 3 单测（分叉图引导翻转/minus 框/stop 终止）覆盖。机制真伪由 Task 8
> probe 的 extended/accepted 统计实证。

- [x] **Step 3: 全量回归**

Run: `cargo test -p trinity-cli`
Expected: 全 PASS

- [x] **Step 4: Commit（冻结期间跳过）**

```bash
git add crates/trinity-cli/tests/protein_guided_e2e.rs
git commit -m "test(cli): M6a 蛋白引导 e2e——低覆盖缺口延长"
```

---

### Task 8: dev1 现库 probe + 入档

**Files:** Modify `docs/opt-log.md`

- [x] **Step 1: 构建 release**

Run: `cargo build --release -p trinity-cli && cargo build --release -p xtask`
Expected: 成功

- [x] **Step 2: dev1 probe（现库, 含拟南芥——循环上限估计）**

```bash
A=/storage/home/senior007/test/athal/bench
cargo run -q --release -p xtask -- bench --tag m6-probe \
  --reads $A/reads/dev1_1.fq.gz,$A/reads/dev1_2.fq.gz --ref $A/refs/dev1.fa \
  -- --path-engine pogflow --protein-guided \
     --protein-db /storage/home/senior007/database/protein.fa \
     --max_memory 16G --CPU 8
```

Expected: R/P/F1 + N50/n; stderr 有 `-- protein-guided: ...` 统计行（extended / accepted / +nt）。

- [x] **Step 3: 对照基线解读**

对照 M5 最强（R=8.6 P=83.6 F1=15.52 N50=1735 n=8165）与 M3 PogFlow（R=8.3 P=84.5 F1=15.11）:
- R 升、P 掉幅 <2pp → 信号正, 报用户定 Gate
- R 平/降 → 用延长统计区分「没延长」（extended≈0: 模板命中或锚定问题）vs「延长了没用/有害」（accepted 高但 R 不动）, 报用户
- R 微升 P 掉 2-4pp → 列 ACCEPT_MARGIN / budget 旋钮敏感性建议, 报用户

- [x] **Step 4: opt-log 入档（新表 + 分析段; 不下 Gate 结论——那是用户的决定）**

在 opt-log.md 追加:

```markdown
## M6 蛋白引导组装——M6a 探路（2026-08-XX, --protein-guided）

| 日期 | tag | 指标 | 其他 | 参数 |
|---|---|---|---|---|
| 2026-08-XX | m6-probe | R=?% P=?% F1=?% | N50=? n=? | --protein-guided --protein-db 现库(含拟南芥, 循环上限估计) |

### 延长统计
- X/Y 转录本有模板, Z 尝试延长, W 验收通过, +N nt
- （对照分析按 Step 3 实测填写）
```

- [x] **Step 5: Commit（冻结期间跳过）**

```bash
git add docs/opt-log.md
git commit -m "docs: M6a 蛋白引导 dev1 现库 probe 入档"
```

---

## 完成定义

- [x] workspace 测试全绿 + clippy/fmt 干净
- [x] e2e 主断言（on_max ≥ 810）通过且参数定稿
- [x] dev1 probe 数字 + 延长统计入档 opt-log
- [x] 未做任何 git 提交（冻结）——等待用户定 Gate 与解冻
