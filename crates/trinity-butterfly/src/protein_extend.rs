//! M6a 蛋白引导延长：HSP 读框感知的密码子翻译 + 原始 k-mer 图引导走图。

/// 标准密码子表（歧义 → X，stop → *）。
pub fn translate_codon(codon: &[u8]) -> u8 {
    let n = |b: u8| match b {
        b'T' => 0,
        b'C' => 1,
        b'A' => 2,
        b'G' => 3,
        _ => 4,
    };
    match codon.len() {
        3 => match (n(codon[0]), n(codon[1]), n(codon[2])) {
            (0, 0, 0) => b'F',
            (0, 0, 1) => b'F',
            (0, 0, 2) => b'L',
            (0, 0, 3) => b'L',
            (0, 1, 0) => b'S',
            (0, 1, 1) => b'S',
            (0, 1, 2) => b'S',
            (0, 1, 3) => b'S',
            (0, 2, 0) => b'Y',
            (0, 2, 1) => b'Y',
            (0, 2, 2) => b'*',
            (0, 2, 3) => b'*',
            (0, 3, 0) => b'C',
            (0, 3, 1) => b'C',
            (0, 3, 2) => b'*',
            (0, 3, 3) => b'W',
            (1, 0, 0) => b'L',
            (1, 0, 1) => b'L',
            (1, 0, 2) => b'L',
            (1, 0, 3) => b'L',
            (1, 1, 0) => b'P',
            (1, 1, 1) => b'P',
            (1, 1, 2) => b'P',
            (1, 1, 3) => b'P',
            (1, 2, 0) => b'H',
            (1, 2, 1) => b'H',
            (1, 2, 2) => b'Q',
            (1, 2, 3) => b'Q',
            (1, 3, 0) => b'R',
            (1, 3, 1) => b'R',
            (1, 3, 2) => b'R',
            (1, 3, 3) => b'R',
            (2, 0, 0) => b'I',
            (2, 0, 1) => b'I',
            (2, 0, 2) => b'I',
            (2, 0, 3) => b'M',
            (2, 1, 0) => b'T',
            (2, 1, 1) => b'T',
            (2, 1, 2) => b'T',
            (2, 1, 3) => b'T',
            (2, 2, 0) => b'N',
            (2, 2, 1) => b'N',
            (2, 2, 2) => b'K',
            (2, 2, 3) => b'K',
            (2, 3, 0) => b'S',
            (2, 3, 1) => b'S',
            (2, 3, 2) => b'R',
            (2, 3, 3) => b'R',
            (3, 0, 0) => b'V',
            (3, 0, 1) => b'V',
            (3, 0, 2) => b'V',
            (3, 0, 3) => b'V',
            (3, 1, 0) => b'A',
            (3, 1, 1) => b'A',
            (3, 1, 2) => b'A',
            (3, 1, 3) => b'A',
            (3, 2, 0) => b'D',
            (3, 2, 1) => b'D',
            (3, 2, 2) => b'E',
            (3, 2, 3) => b'E',
            (3, 3, 0) => b'G',
            (3, 3, 1) => b'G',
            (3, 3, 2) => b'G',
            (3, 3, 3) => b'G',
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
///
/// `kmer_size` 此处 = **图顶点名长度**（chrysalis graph.out 为 K−1-mer,
/// 取 `BflyContext::kmer_size`, 不要传 `--KMER_SIZE`——名字不等长时锚定与
/// 后继校验都会失败）。
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
    let anchor_kmer: Vec<u8> = match dir {
        Direction::Right => transcript[transcript.len() - kmer_size..].to_ascii_uppercase(),
        Direction::Left => transcript[..kmer_size].to_ascii_uppercase(),
    };
    let mut by_name: HashMap<&str, i32> = HashMap::new();
    for &vid in graph.vertex_ids() {
        if let Some(v) = graph.get_vertex(vid) {
            by_name.insert(v.name.as_str(), vid);
        }
    }
    let Some(cur) = by_name.get(std::str::from_utf8(&anchor_kmer).unwrap_or("")) else {
        return Vec::new();
    };
    let mut cur = *cur;

    // 绝对坐标 → 碱基。buf = 到达序延长区: 右 buf[i]↔qlen+1+i, 左 buf[i]↔−i。
    let base_at = |buf: &[u8], pos: i64| -> Option<u8> {
        let idx = match dir {
            Direction::Right => pos - qlen - 1,
            Direction::Left => -pos,
        };
        if idx >= 0 && (idx as usize) < buf.len() {
            Some(buf[idx as usize])
        } else if (1..=qlen).contains(&pos) {
            Some(transcript[(pos - 1) as usize].to_ascii_uppercase())
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
            Direction::Left => -(ext.len() as i64),
        };
        let (lo, hi) = anchor.window_of_frontier(frontier);
        let need: Vec<i64> = (lo..=hi).filter(|&p| base_at(&ext, p).is_none()).collect();
        if need.is_empty() {
            break; // 防御: 理论不可达
        }
        let Some((bases, next_v)) = dfs_codon_paths(
            graph,
            kmer_size,
            cur,
            &need,
            dir,
            &ext,
            &base_at,
            anchor,
            template_aa,
            &visited,
            (lo, hi),
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
        // 环守卫按密码子端点记; 密码子中点可复访（de Bruijn 走图语义合法）
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
    kmer_size: usize,
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
        // 逆序压栈 → 同层弹出按顶点 id 升序; LIFO 深探使全局平局取首个完整路径（确定性, 但非全局 id 升序）
        for nv in succs.into_iter().rev() {
            if visited.contains(&nv) {
                continue;
            }
            let Some(vv) = graph.get_vertex(nv) else {
                continue;
            };
            if vv.name.len() != kmer_size {
                continue; // 非 k-mer 名（压缩合并图不适用本走图的前提）
            }
            let name = vv.name.as_bytes();
            if name.is_empty() {
                continue;
            }
            let b = match dir {
                Direction::Right => name[name.len() - 1].to_ascii_uppercase(),
                Direction::Left => name[0].to_ascii_uppercase(),
            };
            let w = graph
                .find_edge(v, nv)
                .map(|e| e.weight)
                .unwrap_or(0.0) // 后继表与边表应一致; 防御性 0.0 不 panic
                .max(0.0);
            let mut nb = bases.clone();
            nb.push(b);
            stack.push((nv, nb, min_w.min(w)));
        }
    }
    best.map(|(_, bases, v)| (bases, v))
}

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
        assert_eq!(translate_codon(b"TTA"), b'L');
        assert_eq!(translate_codon(b"TTG"), b'L');
        assert_eq!(translate_codon(b"TAT"), b'Y');
        assert_eq!(translate_codon(b"TRY"), b'X'); // 歧义
    }

    #[test]
    fn revcomp_works() {
        assert_eq!(revcomp(b"ATGC"), b"GCAT");
    }

    /// 全 64 密码子对照标准遗传密码（NCBI 表 1）。
    #[test]
    fn codon_table_complete_standard_code() {
        const STD: [(&str, u8); 64] = [
            ("TTT", b'F'),
            ("TTC", b'F'),
            ("TTA", b'L'),
            ("TTG", b'L'),
            ("CTT", b'L'),
            ("CTC", b'L'),
            ("CTA", b'L'),
            ("CTG", b'L'),
            ("ATT", b'I'),
            ("ATC", b'I'),
            ("ATA", b'I'),
            ("ATG", b'M'),
            ("GTT", b'V'),
            ("GTC", b'V'),
            ("GTA", b'V'),
            ("GTG", b'V'),
            ("TCT", b'S'),
            ("TCC", b'S'),
            ("TCA", b'S'),
            ("TCG", b'S'),
            ("CCT", b'P'),
            ("CCC", b'P'),
            ("CCA", b'P'),
            ("CCG", b'P'),
            ("ACT", b'T'),
            ("ACC", b'T'),
            ("ACA", b'T'),
            ("ACG", b'T'),
            ("GCT", b'A'),
            ("GCC", b'A'),
            ("GCA", b'A'),
            ("GCG", b'A'),
            ("TAT", b'Y'),
            ("TAC", b'Y'),
            ("TAA", b'*'),
            ("TAG", b'*'),
            ("CAT", b'H'),
            ("CAC", b'H'),
            ("CAA", b'Q'),
            ("CAG", b'Q'),
            ("AAT", b'N'),
            ("AAC", b'N'),
            ("AAA", b'K'),
            ("AAG", b'K'),
            ("GAT", b'D'),
            ("GAC", b'D'),
            ("GAA", b'E'),
            ("GAG", b'E'),
            ("TGT", b'C'),
            ("TGC", b'C'),
            ("TGA", b'*'),
            ("TGG", b'W'),
            ("CGT", b'R'),
            ("CGC", b'R'),
            ("CGA", b'R'),
            ("CGG", b'R'),
            ("AGT", b'S'),
            ("AGC", b'S'),
            ("AGA", b'R'),
            ("AGG", b'R'),
            ("GGT", b'G'),
            ("GGC", b'G'),
            ("GGA", b'G'),
            ("GGG", b'G'),
        ];
        for (codon, aa) in STD {
            assert_eq!(
                translate_codon(codon.as_bytes()),
                aa,
                "codon {codon} should be {}",
                aa as char
            );
        }
    }

    /// plus 框 HSP: qs=2, sstart=5。
    /// aa j=0 ↔ [2,4], j=1 ↔ [5,7]; 窗口 [2,4] → aa 5; [8,10] → aa 7。
    #[test]
    fn plus_frame_windows() {
        let hsp = AnchorHsp {
            qstart: 2,
            sstart: 5,
            plus: true,
        };
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
        let hsp = AnchorHsp {
            qstart: 10,
            sstart: 3,
            plus: false,
        };
        assert_eq!(hsp.aa_of_window(8, 10), Some(3));
        assert_eq!(hsp.aa_of_window(5, 7), Some(4));
        assert_eq!(hsp.aa_of_window(11, 13), Some(2), "右侧越界 j=-1");
        assert_eq!(hsp.aa_of_window(7, 9), None, "未对齐窗口");
    }

    /// 翻译一个窗口: plus 顺读, minus 反读。
    #[test]
    fn window_translation_respects_strand() {
        let plus = AnchorHsp {
            qstart: 1,
            sstart: 1,
            plus: true,
        };
        assert_eq!(plus.translate_window(b"ATG"), b'M');
        let minus = AnchorHsp {
            qstart: 10,
            sstart: 1,
            plus: false,
        };
        // 反读: CAT → revcomp ATG → M
        assert_eq!(minus.translate_window(b"CAT"), b'M');
    }

    #[test]
    fn window_of_frontier_aligned() {
        let plus = AnchorHsp {
            qstart: 2,
            sstart: 5,
            plus: true,
        };
        // frontier 7: 7-2=5, 5%3=2 → lo=7-2=5 → 窗口 [5,7]（覆盖 7）
        assert_eq!(plus.window_of_frontier(7), (5, 7));
        // frontier 8: 8-2=6 → lo=8 → [8,10]
        assert_eq!(plus.window_of_frontier(8), (8, 10));
        let minus = AnchorHsp {
            qstart: 10,
            sstart: 3,
            plus: false,
        };
        // frontier 6: 10-6=4, 4%3=1 → hi=6+1=7 → [5,7]
        assert_eq!(minus.window_of_frontier(6), (5, 7));
        // frontier 5: 10-5=5, 5%3=2 → hi=5+2=7 → [5,7]（[3,5] 不在框上）
        assert_eq!(minus.window_of_frontier(5), (5, 7));
    }

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
        let ext = extend_transcript(
            &g,
            6,
            b"AAAAAC",
            Direction::Right,
            &AnchorHsp {
                qstart: 1,
                sstart: 1,
                plus: true,
            },
            b"",
            10,
        );
        assert_eq!(ext, b"TGG");
    }

    #[test]
    fn guided_walk_overrides_weight_on_codon_match() {
        let g = branch_graph();
        // qs=1 → 窗口 [7,9] = aa 3。模板 aa3='V': 分支 GTC→V 匹配(2×3.0=6),
        // 主干 TGG→W 失配(10×0.1=1) → 轻分支胜出。
        let template = b"MMVMM";
        let ext = extend_transcript(
            &g,
            6,
            b"AAAAAC",
            Direction::Right,
            &AnchorHsp {
                qstart: 1,
                sstart: 1,
                plus: true,
            },
            template,
            10,
        );
        assert_eq!(
            ext,
            b"GTC",
            "guided ext = {:?}",
            String::from_utf8_lossy(&ext)
        );
    }

    #[test]
    fn budget_caps_extension() {
        let g = branch_graph();
        let ext = extend_transcript(
            &g,
            6,
            b"AAAAAC",
            Direction::Right,
            &AnchorHsp {
                qstart: 1,
                sstart: 1,
                plus: true,
            },
            b"",
            2,
        );
        assert_eq!(ext.len(), 2, "预算在密码子中途截断");
    }

    #[test]
    fn left_walk_uses_predecessors() {
        let g = branch_graph();
        // 从 "AAAAAC" 左延: 前驱 v1 "AAAAAA", 追加首字符 'A', v1 无前驱停。
        let ext = extend_transcript(
            &g,
            6,
            b"AAAAAC",
            Direction::Left,
            &AnchorHsp {
                qstart: 6,
                sstart: 1,
                plus: true,
            },
            b"",
            5,
        );
        assert_eq!(ext, b"A");
    }

    /// 左延多碱基次序: 前驱链 p1="TAAAAA"→p2="GTAAAA"→p3="CGTAAA"→p4="ACGTAA"
    /// (各自与前一个重叠 5)。窗口序 [0,2](1 碱基 T) + [-3,-1](3 碱基 G,C,A),
    /// 预算 3 截到到达序 [T,G,C] → 返回前插序 "CGT"。
    #[test]
    fn left_walk_multi_base_order() {
        // 注: 不复用 branch_graph —— 其 v1("AAAAAA")→v2 重边(10.0)会在首个
        // 无模板信息窗口按纯权重夺走前驱, 干扰本测试的前驱链。
        let mut g = DiGraph::new();
        for (id, name) in [
            (2, "AAAAAC"),
            (21, "TAAAAA"),
            (22, "GTAAAA"),
            (23, "CGTAAA"),
            (24, "ACGTAA"),
        ] {
            g.add_vertex(SeqVertex::new(id, name));
        }
        g.add_edge(21, 2, SimpleEdge::new(5.0, 21, 2));
        g.add_edge(22, 21, SimpleEdge::new(5.0, 22, 21));
        g.add_edge(23, 22, SimpleEdge::new(5.0, 23, 22));
        g.add_edge(24, 23, SimpleEdge::new(5.0, 24, 23));
        let ext = extend_transcript(
            &g,
            6,
            b"AAAAAC",
            Direction::Left,
            &AnchorHsp {
                qstart: 6,
                sstart: 1,
                plus: true,
            },
            b"",
            3,
        );
        assert_eq!(ext, b"CGT", "前插序应为到达序的逆序");
    }

    /// minus 框右延引导: qs=6(高侧), sstart=2 → 窗口 [7,9] = aa 1,
    /// 密码子反读。分支 "GTC" revcomp="GAC"→D 匹配模板 aa1 → 轻分支胜。
    #[test]
    fn minus_frame_guided_walk() {
        let g = branch_graph();
        let ext = extend_transcript(
            &g,
            6,
            b"AAAAAC",
            Direction::Right,
            &AnchorHsp {
                qstart: 6,
                sstart: 2,
                plus: false,
            },
            b"DMM",
            10,
        );
        assert_eq!(
            ext,
            b"GTC",
            "minus 框引导 ext = {:?}",
            String::from_utf8_lossy(&ext)
        );
    }

    /// stop 终止: 模板长 3, aa3 窗口 [7,9] 译出 TAA=*
    /// → 停止, 不再走 c4(高权后续)。
    #[test]
    fn stop_codon_terminates() {
        let mut g = DiGraph::new();
        for (id, name) in [
            (2, "AAAAAC"),
            (31, "AAAACT"), // +T
            (32, "AAACTA"), // +A
            (33, "AACTAA"), // +A → "TAA" 完成密码子
            (34, "ACTAAG"), // +G (若无 stop 规则会继续)
        ] {
            g.add_vertex(SeqVertex::new(id, name));
        }
        for (u, v) in [(2, 31), (31, 32), (32, 33), (33, 34)] {
            g.add_edge(u, v, SimpleEdge::new(5.0, u, v));
        }
        let ext = extend_transcript(
            &g,
            6,
            b"AAAAAC",
            Direction::Right,
            &AnchorHsp {
                qstart: 1,
                sstart: 1,
                plus: true,
            },
            b"MMA",
            10,
        );
        assert_eq!(
            ext,
            b"TAA",
            "模板末位密码子为 stop → 截停: {:?}",
            String::from_utf8_lossy(&ext)
        );
    }

    /// 模板末位无 stop → 不因越过模板末端而截停（M6a'' 撤回了 span-stop
    /// 规则, 真实 UTR 会延伸到模板蛋白之外）: 模板长 3, aa3 窗口译 "TCA"=S
    /// 后继续走图, 直到无路（图末）或预算尽。
    #[test]
    fn walk_continues_past_template_end_without_stop() {
        let mut g = DiGraph::new();
        for (id, name) in [
            (2, "AAAAAC"),
            (41, "AAAACT"), // +T
            (42, "AAACTC"), // +C
            (43, "AACTCA"), // +A → "TCA" 完成模板末位密码子
            (44, "ACTCAG"), // +G (模板外, 仍可行)
            (45, "CTCAGC"), // +C
            (46, "TCAGCT"), // +T (完成窗口 [10,12], 图末)
        ] {
            g.add_vertex(SeqVertex::new(id, name));
        }
        for (u, v) in [(2, 41), (41, 42), (42, 43), (43, 44), (44, 45), (45, 46)] {
            g.add_edge(u, v, SimpleEdge::new(5.0, u, v));
        }
        let ext = extend_transcript(
            &g,
            6,
            b"AAAAAC",
            Direction::Right,
            &AnchorHsp {
                qstart: 1,
                sstart: 1,
                plus: true,
            },
            b"MMS",
            10,
        );
        assert_eq!(
            ext,
            b"TCAGCT",
            "无 stop → 走到图末/预算尽: {:?}",
            String::from_utf8_lossy(&ext)
        );
    }

    #[test]
    fn anchor_kmer_missing_returns_empty() {
        let g = branch_graph();
        let ext = extend_transcript(
            &g,
            6,
            b"TTTTTT",
            Direction::Right,
            &AnchorHsp {
                qstart: 1,
                sstart: 1,
                plus: true,
            },
            b"",
            5,
        );
        assert!(ext.is_empty());
    }
}
