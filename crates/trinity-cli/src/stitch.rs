//! 输出层 split-read 缝合 pass（2026-09-02, 断点解剖驱动）。
//!
//! 机制背景: QM4 蛋白层诊断（docs/opt-log.md「碎片缝合路线」）——缺失蛋白的
//! 74% 是"碎片紧邻"型: Butterfly/图简化在分支处把同一转录本拆成首尾相接的
//! 独立 contig（82% 无 isoform 分歧证据, 断点处读深 95~2055×）。碎片按序
//! 缝合即全长（201/201 ≥95% 覆盖）, 序列 identity 无虞（C 类=0）。
//!
//! 判据（不依赖任何参照, 无评估循环）: contig A 的**末端锚**与 contig B 的
//! **首端锚**（各 30bp）在**同一条 read** 内共现（正反链皆收）→ 桥证据。
//! 防误连三件套:
//! 1. 锚唯一性——锚 30bp 在全部 contigs 中出现次数 ==1, 否则弃;
//! 2. 桥 read 数下限（默认 ≥2）;
//! 3. 连接约束——每端至多一桥（多候选取桥数最大者, 平局按名字确定性）,
//!    并查集禁环 → 缝合关系必为链（DAG 退化为一组路径）。
//!
//! 接缝: A 尾与 B 头取 0..=60bp 内**最大精确重叠**（suffix==prefix）去重拼接,
//! 无重叠则直接相连（紧邻型主形态）。缝合链以链首 contig 命名输出。

use std::collections::HashMap;
use std::path::Path;

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use trinity_common::error::CommonError;
use trinity_common::io_util::open_maybe_gz;

use crate::args::SeqType;

#[derive(Debug, Clone)]
pub struct StitchParams {
    /// 断点两侧锚长（bp）; 30 = split-read 共现的最小可靠窗
    pub anchor_len: usize,
    /// 桥成立的最低 read 数
    pub min_bridge_reads: u32,
    /// 同 read 内尾锚结束与头锚开始之间允许的最大间距（bp）
    pub max_gap: usize,
    /// 接缝去重检查的最大重叠（bp）
    pub max_seam_overlap: usize,
}

impl Default for StitchParams {
    fn default() -> Self {
        StitchParams { anchor_len: 30, min_bridge_reads: 2, max_gap: 40, max_seam_overlap: 60 }
    }
}

#[derive(Debug, Default)]
pub struct StitchStats {
    pub contigs_in: usize,
    pub contigs_out: usize,
    /// 锚构建后存活的（尾/头）锚数
    pub anchors_live: usize,
    /// read 扫描产生的原始桥对数（去重前）
    pub bridges_total: usize,
    /// 达 read 数阈值的桥对数
    pub bridges_kept: usize,
    /// 端唯一化 + 禁环后被采纳的桥数
    pub bridges_applied: usize,
    /// 缝合出的链数（长度 ≥2）
    pub chains: usize,
}

/// 2-bit 打包的锚键（anchor_len ≤ 32 时用 u64; N 等非 gatc 断键返回 None）。
fn pack_key(seq: &[u8]) -> Option<u64> {
    let mut v: u64 = 0;
    for &c in seq {
        let x = match c {
            b'G' | b'g' => 0,
            b'A' | b'a' => 1,
            b'T' | b't' => 2,
            b'C' | b'c' => 3,
            _ => return None,
        };
        v = (v << 2) | x as u64;
    }
    Some(v)
}
/// 对 `fasta`（harvest 产物）就地执行缝合, 就地重写文件。
pub fn stitch_fasta(
    fasta: &Path,
    read_files: &[std::path::PathBuf],
    seq_type: SeqType,
    params: &StitchParams,
) -> Result<StitchStats, CommonError> {
    let raw = std::fs::read(fasta)?;
    let contigs = parse_fasta(&raw)?;
    let mut stats = StitchStats { contigs_in: contigs.len(), ..Default::default() };
    if contigs.len() < 2 {
        stats.contigs_out = contigs.len();
        std::fs::write(fasta, &raw)?;
        return Ok(stats);
    }

    // ---- 1. 锚构建: 每条 contig 的尾/头锚（长度不足 anchor*2 的 contig 不参与）----
    // v1 端位唯一: 同串只允许作为**一个端**的锚出现（端位计数==1）; 串在别的
    // contig **内部**出现不再剔除（v0 全局唯一砍掉 60% 锚——isoform 家族共享
    // 末端使端锚串普遍在兄弟 contig 内部出现, 而紧邻断点恰密集于共享区）。
    // 内部命中混入的风险由桥竞争仲裁（端唯一化取 read 数最大）+ 接缝几何兜底。
    let k = params.anchor_len;
    // 键 -> (contig 下标, 是否尾锚); 端位出现计数
    let mut anchors: FxHashMap<u64, (u32, bool)> = FxHashMap::default();
    let mut end_seen: FxHashMap<u64, u32> = FxHashMap::default();
    for (i, (_, seq)) in contigs.iter().enumerate() {
        if seq.len() < k * 2 {
            continue;
        }
        for (key, is_tail) in [
            (pack_key(&seq[seq.len() - k..]), true),
            (pack_key(&seq[..k]), false),
        ] {
            let Some(key) = key else { continue };
            *end_seen.entry(key).or_insert(0) += 1;
            anchors.insert(key, (i as u32, is_tail));
        }
    }
    anchors.retain(|key, _| end_seen.get(key).copied().unwrap_or(0) == 1);
    stats.anchors_live = anchors.len();
    if anchors.is_empty() {
        stats.contigs_out = contigs.len();
        std::fs::write(fasta, &raw)?;
        return Ok(stats);
    }

    // ---- 3. 扫 reads: 同 read 内 (尾锚, 头锚) 共现 → 桥计数 ----
    // bridges[(a, b)] = a 尾接 b 头的 read 数。
    let mut bridges: FxHashMap<(u32, u32), u32> = FxHashMap::default();
    for f in read_files {
        scan_reads(f, seq_type, k, &anchors, &mut bridges, params.max_gap)?;
    }
    stats.bridges_total = bridges.len();
    bridges.retain(|_, n| *n >= params.min_bridge_reads);
    stats.bridges_kept = bridges.len();

    // ---- 4. 端唯一化: 每尾至多一出桥 / 每头至多一入桥 ----
    // 出桥: (a, b, n); 平局按 b 序列字典序（确定性, 与 rayon 归约序无关）。
    let mut best_out: HashMap<u32, (u32, u32)> = HashMap::new(); // a -> (b, n)
    let mut best_in: HashMap<u32, (u32, u32)> = HashMap::new(); // b -> (a, n)
    let mut cand: Vec<((u32, u32), u32)> = bridges.into_iter().collect();
    cand.sort_by(|x, y| y.1.cmp(&x.1).then_with(|| x.0.cmp(&y.0)));
    let mut kept: Vec<((u32, u32), u32)> = Vec::new();
    for ((a, b), n) in cand {
        if best_out.contains_key(&a) || best_in.contains_key(&b) {
            continue;
        }
        best_out.insert(a, (b, n));
        best_in.insert(b, (a, n));
        kept.push(((a, b), n));
    }

    // ---- 5. 并查集禁环: 按 n 降序（kept 已序）贪心加边 ----
    let mut uf = UnionFind::new(contigs.len());
    let mut applied: Vec<((u32, u32), u32)> = Vec::new();
    for &e @ ((a, b), _) in &kept {
        if uf.union(a, b) {
            applied.push(e);
        }
    }
    stats.bridges_applied = applied.len();

    // ---- 6. 链合并 + 重写 ----
    // 出边表（端唯一化保证每 a 至多一）+ 前驱表
    let mut next: HashMap<u32, u32> = HashMap::new();
    for &((a, b), _) in &applied {
        next.insert(a, b);
    }
    let mut prev = vec![u32::MAX; contigs.len()];
    for (&a, &b) in next.iter() {
        prev[b as usize] = a;
    }

    let mut removed = vec![false; contigs.len()];
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    for i in 0..contigs.len() {
        if prev[i] != u32::MAX || removed[i] {
            continue; // 非链首（或已被并入前序链——后者不可达, 防御）
        }
        let is_chain = next.contains_key(&(i as u32));
        let mut seq = contigs[i].1.clone();
        let mut cur = i as u32;
        while let Some(&nxt) = next.get(&cur) {
            let n = nxt as usize;
            let seam = max_exact_overlap(&seq, &contigs[n].1, params.max_seam_overlap);
            // 重叠 = B 头部 r bp（suffix(seq)==prefix(B)）→ 拼接 B[r..]
            seq.extend_from_slice(&contigs[n].1[seam..]);
            removed[n] = true;
            cur = nxt;
        }
        if is_chain {
            stats.chains += 1;
        }
        out.push((contigs[i].0.clone(), seq));
    }
    stats.contigs_out = out.len();
    let mut buf = String::new();
    for (name, seq) in &out {
        buf.push('>');
        buf.push_str(name);
        buf.push('\n');
        for chunk in seq.chunks(60) {
            buf.push_str(std::str::from_utf8(chunk).map_err(|_| CommonError::Parse(
                "stitch: 非 UTF-8 序列".to_string(),
            ))?);
            buf.push('\n');
        }
    }
    std::fs::write(fasta, buf)?;
    Ok(stats)
}

/// fasta 字节解析为 (name, seq)（大写化; 空序列跳过）。
fn parse_fasta(data: &[u8]) -> Result<Vec<(String, Vec<u8>)>, CommonError> {
    let mut out = Vec::new();
    let mut name: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();
    for line in data.split(|&b| b == b'\n') {
        let line = if line.last() == Some(&b'\r') { &line[..line.len() - 1] } else { line };
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if let Some(n) = name.take() {
                if !seq.is_empty() {
                    out.push((n, std::mem::take(&mut seq)));
                }
            }
            name = Some(String::from_utf8_lossy(&line[1..])
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string());
        } else {
            seq.extend(line.iter().map(|&b| b.to_ascii_uppercase()));
        }
    }
    if let Some(n) = name {
        if !seq.is_empty() {
            out.push((n, seq));
        }
    }
    Ok(out)
}

/// 扫单个 reads 文件（v1 并行版）: 整读入内存 → 切序列行区间 → rayon 分块
/// 并行扫描 → 局部桥表归并。命中 (尾锚@pos, 头锚@pos') 且接缝几何成立 → 计数。
fn scan_reads(
    path: &Path,
    seq_type: SeqType,
    k: usize,
    anchors: &FxHashMap<u64, (u32, bool)>,
    bridges: &mut FxHashMap<(u32, u32), u32>,
    max_gap: usize,
) -> Result<(), CommonError> {
    let mut reader = open_maybe_gz(path)?;
    let mut data = Vec::new();
    std::io::Read::read_to_end(&mut reader, &mut data)?;

    // 序列字节区间收集: fq 按 4 行组（第 2 行 = 序列）; fa 按记录（序列行拼接,
    // 多行时退化为逐行区间——跨行 k-mer 会断窗, 但 fa 输入在此路径为归一化
    // 产物（单行）, 实际不触发）。
    let mut seqs: Vec<(usize, usize)> = Vec::new();
    let mut pos = 0usize;
    let mut line_no = 0usize;
    while pos < data.len() {
        let end = memchr(b'\n', &data[pos..]).map(|i| pos + i).unwrap_or(data.len());
        let content_end = if end > pos && data[end - 1] == b'\r' { end - 1 } else { end };
        let is_hdr = content_end > pos
            && ((seq_type == SeqType::Fq && data[pos] == b'@')
                || (seq_type == SeqType::Fa && data[pos] == b'>'));
        let is_blank = content_end == pos;
        match seq_type {
            SeqType::Fq => {
                let m = line_no % 4;
                if m == 0 && !is_hdr && !is_blank {
                    return Err(CommonError::Parse(format!(
                        "stitch: 非 fq 记录边界（{} 第 {} 行）",
                        path.display(),
                        line_no + 1
                    )));
                }
                if m == 1 && !is_blank {
                    seqs.push((pos, content_end));
                }
                if m == 0 && !is_hdr {
                    line_no += 1; // 空行宽容: 只推进不校验
                    pos = end + 1;
                    continue;
                }
            }
            SeqType::Fa => {
                if !is_blank && !is_hdr {
                    seqs.push((pos, content_end));
                }
            }
        }
        line_no += 1;
        pos = end + 1;
    }

    // 并行扫描: 每 4096 条序列一块, 局部桥表 + 归并（小并大）。
    let partials: Vec<FxHashMap<(u32, u32), u32>> = seqs
        .par_chunks(4096)
        .map(|chunk| {
            let mut local: FxHashMap<(u32, u32), u32> = FxHashMap::default();
            let mut rc_buf: Vec<u8> = Vec::new();
            for &(s, e) in chunk {
                scan_one_read(&data[s..e], k, anchors, &mut local, max_gap, &mut rc_buf);
            }
            local
        })
        .collect();
    for p in partials {
        for (key, c) in p {
            *bridges.entry(key).or_insert(0) += c;
        }
    }
    Ok(())
}

/// 单条 read 扫描: 正反链各滚一遍（rc 复用缓冲免逐条分配）。
fn scan_one_read(
    seq: &[u8],
    k: usize,
    anchors: &FxHashMap<u64, (u32, bool)>,
    bridges: &mut FxHashMap<(u32, u32), u32>,
    max_gap: usize,
    rc_buf: &mut Vec<u8>,
) {
    if seq.len() < 2 * k {
        return;
    }
    let mut hits: Vec<(usize, u32, bool)> = Vec::new();
    for flip in [false, true] {
        hits.clear();
        let strand: &[u8] = if flip {
            rc_buf.clear();
            rc_buf.extend(seq.iter().rev().map(|&c| match c {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                _ => b'N',
            }));
            rc_buf
        } else {
            seq
        };
        let mask = (1u64 << (2 * k as u32)) - 1;
        let mut valid = 0usize;
        let mut v: u64 = 0;
        for i in 0..strand.len() {
            let x = match strand[i] {
                b'G' | b'g' => 0,
                b'A' | b'a' => 1,
                b'T' | b't' => 2,
                b'C' | b'c' => 3,
                _ => {
                    valid = 0;
                    v = 0;
                    continue;
                }
            };
            v = ((v << 2) | x as u64) & mask;
            valid += 1;
            if valid >= k {
                if let Some(&(c, is_tail)) = anchors.get(&v) {
                    hits.push((i + 1, c, is_tail)); // i+1 = 窗口右端(含)
                }
            }
        }
        // 桥判定: 尾锚区间 [e-k, e), 头锚区间 [s2, s2+k), s2 = e2-k。
        // 顺序: 头锚起点不早于尾锚起点（允许锚重叠）; 间距: e - s2 <= k + max_gap。
        for &(e, ca, ta) in &hits {
            if !ta {
                continue;
            }
            for &(e2, cb, hb) in &hits {
                if hb {
                    continue;
                }
                let s2 = e2.saturating_sub(k);
                if ca != cb && s2 >= e.saturating_sub(k) && e <= s2 + k + max_gap {
                    *bridges.entry((ca, cb)).or_insert(0) += 1;
                }
            }
        }
    }
}

fn memchr(needle: u8, hay: &[u8]) -> Option<usize> {
    hay.iter().position(|&b| b == needle)
}

/// suffix(a) 与 prefix(b) 在 0..=max 内的最大精确重叠长度。
fn max_exact_overlap(a: &[u8], b: &[u8], max: usize) -> usize {
    let hi = max.min(a.len()).min(b.len());
    for r in (1..=hi).rev() {
        if a[a.len() - r..] == b[..r] {
            return r;
        }
    }
    0
}

struct UnionFind {
    parent: Vec<u32>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind { parent: (0..n as u32).collect() }
    }
    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            self.parent[x as usize] = self.parent[self.parent[x as usize] as usize];
            x = self.parent[x as usize];
        }
        x
    }
    /// 合并成功 true; 已同连通分量（会成环）false。
    fn union(&mut self, a: u32, b: u32) -> bool {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return false;
        }
        self.parent[rb as usize] = ra;
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmpdir(tag: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("stitch_{tag}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    /// xorshift31 伪随机 ACGT 序列（避免重复串锚被唯一性剔除）。
    fn rand_seq(seed: u32, n: usize) -> Vec<u8> {
        let mut x = seed | 1;
        let mut s = Vec::with_capacity(n);
        for _ in 0..n {
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;
            s.push(b"ACGT"[(x % 4) as usize]);
        }
        s
    }

    /// 3 contigs: A→B 紧邻(0 重叠), B→C 20bp 重叠; 横跨 read 各 2 条(正反链)。
    #[test]
    fn stitch_chain_with_overlap_seam() {
        let d = tmpdir("chain");
        // A 尾 20bp = B 头 20bp（接缝重叠）; 其余随机不相交
        let head_a = rand_seq(11, 100);
        let seam: Vec<u8> = rand_seq(23, 20); // 20bp 重叠
        let mid_b = rand_seq(35, 100);
        // 挑一个与 B 尾无偶然重叠的 C（排除启发式接缝剪头的干扰, 聚焦 0 重叠形态）
        let mut seed_c = 47;
        let contig_b = [seam.clone(), mid_b.clone()].concat();
        let tail_c = loop {
            let cand = rand_seq(seed_c, 100);
            if max_exact_overlap(&contig_b, &cand, 60) == 0 {
                break cand;
            }
            seed_c += 1;
        };
        let contig_a = [head_a.clone(), seam.clone()].concat();
        let contig_c = [tail_c.clone()].concat();
        // 真实转录本: A + B(去 20 重叠) + C(紧邻 0 重叠)
        let fa = d.join("t.fasta");
        let mut txt = String::new();
        for (n, s) in [("A", &contig_a), ("B", &contig_b), ("C", &contig_c)] {
            txt.push_str(&format!(">{n}\n{}\n", String::from_utf8_lossy(s).into_owned()));
        }
        std::fs::write(&fa, txt).unwrap();
        // reads: 横跨 A/B 断点（A 尾 30 + B 头 30 紧接）×2; 横跨 B/C ×2
        let ab = [contig_a[contig_a.len() - 40..].to_vec(), contig_b[..40].to_vec()].concat();
        let bc = [contig_b[contig_b.len() - 40..].to_vec(), contig_c[..40].to_vec()].concat();
        let mut fq = String::new();
        let mut i = 0;
        for s in [&ab, &ab, &bc, &bc] {
            i += 1;
            fq.push_str(&format!("@r{i}\n{}\n+\n{}\n", String::from_utf8_lossy(s).into_owned(), "I".repeat(s.len())));
        }
        let rf = d.join("r.fq");
        std::fs::write(&rf, fq).unwrap();
        let p = StitchParams::default();
        let st = stitch_fasta(&fa, &[rf], SeqType::Fq, &p).unwrap();
        assert_eq!(st.contigs_in, 3);
        assert_eq!(st.contigs_out, 1, "应缝合成单链: {st:?}");
        assert_eq!(st.bridges_applied, 2);
        // 校验缝合序列 = A + B(去 20 重叠) + C
        let expect = [contig_a.clone(), contig_b[20..].to_vec(), contig_c.clone()].concat();
        let data = std::fs::read(&fa).unwrap();
        let out = parse_fasta(&data).unwrap();
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].1, expect, "接缝去重与拼接序列");
        let _ = std::fs::remove_dir_all(&d);
    }

    /// 环防护: A→B 与 B→A 双桥只留一（并查集拒绝后加的）。
    #[test]
    fn cycle_rejected() {
        let d = tmpdir("cycle");
        let s1: Vec<u8> = rand_seq(101, 100);
        let s2: Vec<u8> = rand_seq(203, 100);
        let fa = d.join("t.fasta");
        std::fs::write(
            &fa,
            format!(">A\n{}\n>B\n{}\n", String::from_utf8_lossy(&s1).into_owned(), String::from_utf8_lossy(&s2).into_owned()),
        )
        .unwrap();
        // read1 含 A 尾30 + B 头30; read2 含 B 尾30 + A 头30 → 双向桥
        let r1 = [s1[70..].to_vec(), s2[..30].to_vec()].concat();
        let r2 = [s2[70..].to_vec(), s1[..30].to_vec()].concat();
        let rf = d.join("r.fq");
        let mut fq = String::new();
        for (i, s) in [&r1, &r1, &r2, &r2].iter().enumerate() {
            fq.push_str(&format!("@r{}\n{}\n+\n{}\n", i + 1, String::from_utf8_lossy(s).into_owned(), "I".repeat(s.len())));
        }
        std::fs::write(&rf, fq).unwrap();
        let st = stitch_fasta(&fa, &[rf], SeqType::Fq, &StitchParams::default()).unwrap();
        assert_eq!(st.bridges_kept, 2);
        assert_eq!(st.bridges_applied, 1, "环必须被拒: {st:?}");
        assert_eq!(st.contigs_out, 1);
        let _ = std::fs::remove_dir_all(&d);
    }

    /// 歧义锚（锚串在其他 contig 内部出现）→ 桥被唯一性剔除。
    #[test]
    fn ambiguous_anchor_dropped() {
        let d = tmpdir("amb");
        let amb = rand_seq(57, 30); // 将同时出现在 A 尾与 B 头（歧义锚）
        let mut a: Vec<u8> = rand_seq(61, 100);
        a.extend_from_slice(&amb);
        let mut bseq = amb.clone();
        bseq.extend_from_slice(&rand_seq(67, 100));
        let fa = d.join("t.fasta");
        std::fs::write(
            &fa,
            format!(">A\n{}\n>B\n{}\n", String::from_utf8_lossy(&a).into_owned(), String::from_utf8_lossy(&bseq).into_owned()),
        )
        .unwrap();
        let r = [a[a.len() - 40..].to_vec(), bseq[..40].to_vec()].concat();
        let rf = d.join("r.fq");
        let fq = format!("@r1\n{}\n+\n{}\n@r2\n{}\n+\n{}\n",
            String::from_utf8_lossy(&r).into_owned(), "I".repeat(r.len()),
            String::from_utf8_lossy(&r).into_owned(), "I".repeat(r.len()));
        std::fs::write(&rf, fq).unwrap();
        let st = stitch_fasta(&fa, &[rf], SeqType::Fq, &StitchParams::default()).unwrap();
        assert_eq!(st.bridges_kept, 0, "歧义锚不得产生桥: {st:?}");
        assert_eq!(st.contigs_out, 2);
        let _ = std::fs::remove_dir_all(&d);
    }

    /// v1 端位唯一: 锚串在别的 contig **内部**出现不再剔除（可缝合）;
    /// 端位出现 2 次（另一端同串）仍是歧义（ambiguous 测试锁定）。
    #[test]
    fn internal_anchor_hit_tolerated() {
        let d = tmpdir("internal");
        let s1: Vec<u8> = rand_seq(501, 100);
        let s2: Vec<u8> = rand_seq(601, 100);
        // C = 别的 contig, 其内部恰好包含 s1 的尾锚 30bp（非端位）
        let a_tail30 = s1[s1.len() - 30..].to_vec();
        let mut c = rand_seq(701, 40);
        c.extend_from_slice(&a_tail30);
        c.extend_from_slice(&rand_seq(801, 40));
        let fa = d.join("t.fasta");
        std::fs::write(
            &fa,
            format!(
                ">A\n{}\n>B\n{}\n>C\n{}\n",
                String::from_utf8_lossy(&s1),
                String::from_utf8_lossy(&s2),
                String::from_utf8_lossy(&c)
            ),
        )
        .unwrap();
        let r = [s1[s1.len() - 40..].to_vec(), s2[..40].to_vec()].concat();
        let rf = d.join("r.fq");
        let mut fq = String::new();
        for i in 0..2 {
            fq.push_str(&format!(
                "@r{}\n{}\n+\n{}\n",
                i + 1,
                String::from_utf8_lossy(&r),
                "I".repeat(r.len())
            ));
        }
        std::fs::write(&rf, fq).unwrap();
        let st = stitch_fasta(&fa, &[rf], SeqType::Fq, &StitchParams::default()).unwrap();
        assert_eq!(st.bridges_applied, 1, "内部命中不阻断缝合: {st:?}");
        assert_eq!(st.contigs_out, 2);
        let _ = std::fs::remove_dir_all(&d);
    }

    /// 桥 read 数不足阈值 → 不缝。
    #[test]
    fn below_min_reads_not_stitched() {
        let d = tmpdir("min");
        let s1: Vec<u8> = rand_seq(301, 100);
        let s2: Vec<u8> = rand_seq(409, 100);
        let fa = d.join("t.fasta");
        std::fs::write(
            &fa,
            format!(">A\n{}\n>B\n{}\n", String::from_utf8_lossy(&s1).into_owned(), String::from_utf8_lossy(&s2).into_owned()),
        )
        .unwrap();
        let r = [s1[70..].to_vec(), s2[..30].to_vec()].concat(); // 仅 1 条
        let rf = d.join("r.fq");
        std::fs::write(&rf, format!("@r1\n{}\n+\n{}\n", String::from_utf8_lossy(&r).into_owned(), "I".repeat(r.len()))).unwrap();
        let st = stitch_fasta(&fa, &[rf], SeqType::Fq, &StitchParams::default()).unwrap();
        assert_eq!(st.bridges_kept, 0);
        assert_eq!(st.contigs_out, 2);
        let _ = std::fs::remove_dir_all(&d);
    }
}
