//! `Chrysalis/analysis/ReadsToTranscripts.cc` 移植——read → Inchworm bundle
//! 的 k-mer 投票分配（谱对齐）。
//!
//! 镜像语义（以原版 413 行逐行为准）：
//! - k=25 固定；bundle 索引 `NonRedKmerTable::set_up_templates(no_ns=true)`
//!   （跨 'X' / 含非大写 ACGT 的窗口不入表），`set_all_counts(-1)` 后逐
//!   bundle 逐位置 `set_count(seq, j, i)` —— **count 字段存 bundle 行下标**，
//!   共享 25-mer 按"后写者赢"（单线程 `#pragma omp parallel for` +
//!   本实现顺序写，i 升序 → bundle i 越大越赢；原版多线程时竞争不确定，
//!   `-t 1` 下与顺序 i 升序一致）；
//! - 每 read toupper 后正向枚举：窗口熵 < min_kmer_entropy 跳过，
//!   `get_count_real` miss(-1) 跳过，命中 push；!strand 时整条 revcomp
//!   再枚举一遍 push 进**同一** comp（DS 分子含两个方向的命中）；
//! - **最长 run 怪癖**（cc:254-268，逐位复刻勿"修复"）：comp 升序排序后
//!   `if (comp[j]!=comp[j-1] || j+1==len)` 才结算 run——非末组 run=m-1、
//!   末组 run=m-2、大小 1 的组永不成 best（见 [`best_component`]）；
//! - `num_kmer_pos = read.len()-25+1`（**仅正向**位置数，DS 分子含 revcomp
//!   命中而分母不变）；`pct = (int)(max/pos*100 + 0.5)`（f32 加 0.5 截断
//!   —— 四舍五入）；
//! - 接受 `best != -1 && pct >= pct_required`；输出行
//!   `<component_id>\t<read_name>\t<pct>%\t<seq>\n`，component_id 从
//!   bundle 名 `>s_N_...` 提取 N（atoi，短止于首个非数字），read_name
//!   保留 '>'、内部空格→'_'（`formatReadNameString`）；
//! - 行序 = `multimap<int,int>(best, read_idx)`：best 升序、同 best 内
//!   read 下标升序（multimap 等键保插入序）。

use rayon::prelude::*;

use crate::dna_vector::{compute_entropy, revcomp, DnaSeq};
use crate::nonred_table::NonRedKmerTable;
use trinity_common::error::CommonError;

/// k=25 固定（cc:130）。
pub const K: usize = 25;

/// 默认 `min_kmer_entropy = 1.5`（cc:76）。
pub const DEFAULT_MIN_KMER_ENTROPY: f32 = 1.5;

#[derive(Debug, Clone)]
pub struct RttParams {
    /// `-strand`：strand specific（true 时不做 revcomp 二次枚举）。
    pub strand: bool,
    /// `-p`：要求的最小映射百分比（Trinity 传 50）。
    pub pct_required: u32,
    /// `-min_kmer_entropy`（默认 1.5）。
    pub min_kmer_entropy: f32,
    /// `-max_mem_reads`（原版仅控制流式分块，不影响输出——本版单块）。
    pub max_mem_reads: usize,
    /// `-t`：read 查询循环的线程数。索引构建（set_count"后写者赢"）
    /// 恒串行保序；查询表只读、每 read 独立 → 显式 rayon 池并行。
    /// 输出按 (best, read_idx) 排序收集 → 与线程数无关（确定性）。
    pub threads: usize,
}

impl Default for RttParams {
    fn default() -> Self {
        RttParams {
            strand: false,
            pct_required: 0,
            min_kmer_entropy: DEFAULT_MIN_KMER_ENTROPY,
            max_mem_reads: usize::MAX,
            threads: 1,
        }
    }
}

/// `reads_to_transcripts` 返回：映射输出全文 + 成功映射 read 总数
/// （原版写到 `<out>.rcts.out` 的 readCount）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RttOutput {
    pub text: String,
    pub mapped_count: u64,
}

/// `atoi` 语义：跳过前导空白，可选符号，取最长数字前缀；无数字 → 0。
fn atoi(s: &str) -> i32 {
    let b = s.as_bytes();
    let mut i = 0;
    while i < b.len() && b[i].is_ascii_whitespace() {
        i += 1;
    }
    let neg = match b.get(i) {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let mut v: i64 = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        v = (v * 10 + (b[i] - b'0') as i64).min(i32::MAX as i64 + 1);
        i += 1;
    }
    if neg {
        v = -v;
    }
    v.clamp(i32::MIN as i64, i32::MAX as i64) as i32
}

/// bundle 名（`read_fasta` 的 '_' 连接全名，如 `>s_0_96_13_17_11`）→
/// 组件号：`substr(3)`（去 `>s_`）后 atoi（短止于首个非数字）。
fn component_number(bundle_acc: &str) -> i32 {
    atoi(&bundle_acc[3.min(bundle_acc.len())..])
}

/// `DNAStringStreamFast::formatReadNameString`（DNAVector.cc:1504-1514）：
/// 去首部空格、内部空格→'_'、去尾部空格（'>' 保留）。
pub fn format_read_name_string(name: &str) -> String {
    let trimmed = name.trim_matches(' ');
    trimmed.replace(' ', "_")
}

/// cc:254-268 最长 run 怪癖——**逐位复刻**，返回 (best, max)。
///
/// 组大小 m 的 run 计数：非末组 m-1、末组 m-2（末组的最后一个元素在
/// `j+1==len` 分支结算时 run 少加一次）；大小 1 的组 run=0 永不 `> max`
/// （max 初值 0）→ 不成 best。空 comp / 全大小 1 → best=-1。
pub fn best_component(comp: &[i32]) -> (i32, i32) {
    let mut best: i32 = -1;
    let mut max: i32 = 0;
    let mut run = 0;
    for j in 1..comp.len() {
        if comp[j] != comp[j - 1] || j + 1 == comp.len() {
            if run > max {
                max = run;
                best = comp[j - 1];
            }
            run = 0;
        } else {
            run += 1;
        }
    }
    (best, max)
}

/// cc:271 `pct_read_mapped = (int)((float)max/num_kmer_pos*100 + 0.5)`
/// （f32 计算 + C 强转截断 = 四舍五入）。
fn pct_read_mapped(max: i32, num_kmer_pos: i32) -> i32 {
    ((max as f32 / num_kmer_pos as f32) * 100.0 + 0.5) as i32
}

/// 主流程：`reads_to_transcripts(reads, bundles, params) -> RttOutput`。
///
/// - `bundles`：`read_fasta` 读入（'_' 连接全名 + toupper）——名字用于
///   提取组件号，序列用于建 25-mer 索引（count 存 bundle 行下标）；
/// - `reads`：`read_fasta_short_names` 读入（name = 首 token 含 '>'；
///   本函数内部再对 name 做 formatReadNameString 变换）。
pub fn reads_to_transcripts(
    reads: &[DnaSeq],
    bundles: &[DnaSeq],
    p: &RttParams,
) -> Result<RttOutput, CommonError> {
    let templates: Vec<Vec<u8>> = bundles.iter().map(|b| b.seq.clone()).collect();
    let mut kt = NonRedKmerTable::set_up_templates(&templates, K, true);
    kt.set_all_counts(-1);

    let component_number_mapping: Vec<i32> =
        bundles.iter().map(|b| component_number(&b.name)).collect();

    // 单线程顺序写 → 共享 25-mer 后写者（更大 bundle 下标）赢
    for (i, b) in bundles.iter().enumerate() {
        if b.seq.len() >= K {
            for j in 0..=b.seq.len() - K {
                kt.set_count(&b.seq, j, i as i32);
            }
        }
    }

    // read 查询循环：表只读、每 read 独立 → 显式线程池并行（`-t`）。
    // 每 read 返回 Option<(read_idx, best, pct)>，随后按 (best, read_idx)
    // 排序输出 —— 与原版 multimap(best, read_idx) 行序等价，且与线程数无关。
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(p.threads.max(1))
        .build()
        .map_err(|e| CommonError::Parse(format!("cannot build RTT thread pool: {e}")))?;
    let assignments: Vec<Option<(usize, i32, i32)>> = pool.install(|| {
        reads
            .par_iter()
            .enumerate()
            .map(|(i, r)| query_read(i, r, &kt, p))
            .collect()
    });

    let mut mapped: Vec<(usize, i32, i32)> = assignments.into_iter().flatten().collect();
    let mapped_count = mapped.len() as u64;
    // multimap<int,int>(best, read_idx)：best 升序、同 best 内 read 下标升序。
    mapped.sort_unstable_by_key(|&(i, best, _)| (best, i));

    let mut text = String::new();
    for &(i, best, pct) in &mapped {
        text.push_str(&component_number_mapping[best as usize].to_string());
        text.push('\t');
        text.push_str(&format_read_name_string(&reads[i].name));
        text.push('\t');
        text.push_str(&format!("{}%\t", pct));
        text.push_str(&String::from_utf8_lossy(&reads[i].seq));
        text.push('\n');
    }

    Ok(RttOutput { text, mapped_count })
}

/// 单 read 查询（cc:225-268 逐位复刻，与原串行版逐行等价）：
/// 正向枚举 + (!strand) revcomp 二次枚举 → 最长 run 怪癖 → pct 四舍五入。
fn query_read(
    i: usize,
    r: &DnaSeq,
    kt: &NonRedKmerTable,
    p: &RttParams,
) -> Option<(usize, i32, i32)> {
    let d: Vec<u8> = r.seq.iter().map(|&c| c.to_ascii_uppercase()).collect();
    let mut comp: Vec<i32> = Vec::with_capacity(4000);
    // 原版 int num_kmer_pos = d.size()-k+1（len<k 时下溢回绕；此时无窗口、
    // best 必为 -1，pct 值无意义 → 安全等价）
    let num_kmer_pos = d.len() as i64 - K as i64 + 1;

    if d.len() >= K {
        for j in 0..=d.len() - K {
            let w = &d[j..j + K];
            if compute_entropy(w) < p.min_kmer_entropy {
                continue;
            }
            let c = kt.get_count_real(&d, j);
            if c >= 0 {
                comp.push(c);
            }
        }
    }
    if !p.strand && d.len() >= K {
        let dd = revcomp(&d);
        for j in 0..=dd.len() - K {
            let w = &dd[j..j + K];
            if compute_entropy(w) < p.min_kmer_entropy {
                continue;
            }
            let c = kt.get_count_real(&dd, j);
            if c >= 0 {
                comp.push(c);
            }
        }
    }

    comp.sort_unstable();
    let (best, max) = best_component(&comp);
    let pct = pct_read_mapped(
        max,
        num_kmer_pos.clamp(i32::MIN as i64, i32::MAX as i64) as i32,
    );
    (best != -1 && pct >= p.pct_required as i32).then_some((i, best, pct))
}

/// 镜像下游 `sort -k1,1n -k3,3nr -k2,2`（GNU sort，无 -s）：
/// 第 1 列数值升序、第 3 列（`45%` 取前导数字）数值降序、第 2 列字典序
/// 升序；全键相等时 last-resort 整行字节序升序。
///
/// 数值键语义（GNU numeric）：跳过前导空白、可选符号、最长数字前缀；
/// 无数字视为 0（数值相等不回退字符串比较，直接进下一键）。
pub fn sort_reads_to_components(text: &str) -> String {
    let mut lines: Vec<&str> = text.lines().collect();
    lines.sort_by(|a, b| {
        let fa: Vec<&str> = a.split('\t').collect();
        let fb: Vec<&str> = b.split('\t').collect();
        let ka = num_key(fa.first().copied().unwrap_or(""));
        let kb = num_key(fb.first().copied().unwrap_or(""));
        if ka != kb {
            return ka.cmp(&kb);
        }
        let na = num_key(fa.get(2).copied().unwrap_or(""));
        let nb = num_key(fb.get(2).copied().unwrap_or(""));
        if na != nb {
            return nb.cmp(&na);
        }
        let sa = fa.get(1).copied().unwrap_or("");
        let sb = fb.get(1).copied().unwrap_or("");
        if sa != sb {
            return sa.cmp(sb);
        }
        a.cmp(b)
    });
    let mut out = lines.join("\n");
    if !out.is_empty() {
        out.push('\n');
    }
    out
}

/// GNU sort 数值键：`[blank]*[-+]?digit*` → (值, )；与 i64 对比较即可
/// （无数字 → 0）。返回 Option 区分正负无效？GNU 对 "+"/"-" 无数字按 0 计。
fn num_key(field: &str) -> i64 {
    let b = field.as_bytes();
    let mut i = 0;
    while i < b.len() && (b[i] == b' ' || b[i] == b'\t') {
        i += 1;
    }
    let neg = match b.get(i) {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let mut v: i64 = 0;
    let mut any = false;
    while i < b.len() && b[i].is_ascii_digit() {
        any = true;
        v = v.saturating_mul(10).saturating_add((b[i] - b'0') as i64);
        i += 1;
    }
    let _ = any;
    if neg {
        -v
    } else {
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---------- 最长 run 怪癖（逐场景手推） ----------

    /// [5,5,5]：末组 run = m-2 = 1 → best=5, max=1。
    #[test]
    fn run_quirk_last_group() {
        assert_eq!(best_component(&[5, 5, 5]), (5, 1));
    }

    /// [3,3,3,7,7]：组 3 非末 run=2，组 7 末 run=0 → best=3, max=2。
    #[test]
    fn run_quirk_nonlast_beats_last() {
        assert_eq!(best_component(&[3, 3, 3, 7, 7]), (3, 2));
    }

    /// 大小 1 的组永不成 best；[1] 与 [1,2] 均 best=-1。
    #[test]
    fn run_quirk_singletons_never_best() {
        assert_eq!(best_component(&[1]), (-1, 0));
        assert_eq!(best_component(&[1, 2]), (-1, 0));
        assert_eq!(best_component(&[]), (-1, 0));
    }

    /// [1,1,1,2,2,2]：组 1（非末）run=2、组 2（末）run=1 → best=1, max=2。
    #[test]
    fn run_quirk_nonlast_group() {
        assert_eq!(best_component(&[1, 1, 1, 2, 2, 2]), (1, 2));
    }

    /// 平局取先到（run 相同不覆盖，`>` 严格）：[1,1,1,2,2,2,2] 组 1 run=2、
    /// 组 2（末，m=4）run=2 → 平局，先结算的组 1 赢。
    #[test]
    fn run_quirk_tie_keeps_first() {
        assert_eq!(best_component(&[1, 1, 1, 2, 2, 2, 2]), (1, 2));
    }

    /// 单组也是"末组"：run = m-2。[0]×11 → max=9（非 10）。
    #[test]
    fn run_quirk_single_group_is_last_group() {
        let comp = vec![0; 11];
        assert_eq!(best_component(&comp), (0, 9));
    }

    /// 末组恰比非末组多 1 个元素时（真实计数 3 vs 2）run 相同（2 vs 2）→
    /// 先结算的非末组赢，即使末组真实命中更多——怪癖核心。
    #[test]
    fn run_quirk_off_by_one_real_counts() {
        // 真实计数：2 出现 3 次、5 出现 4 次；run：组2→2，组5（末）→2
        assert_eq!(best_component(&[2, 2, 2, 5, 5, 5, 5]), (2, 2));
    }

    // ---------- pct 四舍五入 ----------

    /// max=10, pos=20 → 50；max=1, pos=200 → 0.5+0.5=1.0 → 1；
    /// max=1, pos=300 → 1/3%+0.5=0.83 → 0（截断）。
    #[test]
    fn pct_rounding() {
        assert_eq!(pct_read_mapped(10, 20), 50);
        assert_eq!(pct_read_mapped(1, 200), 1);
        assert_eq!(pct_read_mapped(1, 300), 0);
        assert_eq!(pct_read_mapped(0, 100), 0); // (0+0.5) 截断 → 0
    }

    // ---------- atoi / component_number / formatReadNameString ----------

    #[test]
    fn atoi_semantics() {
        assert_eq!(atoi("0_96_13_17_11"), 0);
        assert_eq!(atoi("12_34"), 12);
        assert_eq!(atoi(">s_1"), 0); // '>' 非数字 → 0
        assert_eq!(atoi("  -7abc"), -7);
        assert_eq!(atoi(""), 0);
    }

    #[test]
    fn component_number_from_bundle_name() {
        assert_eq!(component_number(">s_0_96_13_17_11"), 0);
        assert_eq!(component_number(">s_1_142_9"), 1);
        assert_eq!(component_number(">s_23_x"), 23);
    }

    /// formatReadNameString：去首/尾空格、内部空格→'_'、保留 '>'。
    #[test]
    fn read_name_transform() {
        assert_eq!(format_read_name_string(">r1 foo bar"), ">r1_foo_bar");
        assert_eq!(format_read_name_string("  >r2  x  "), ">r2__x");
        assert_eq!(format_read_name_string(">plain"), ">plain");
    }

    // ---------- 端到端（手造小型 25-mer 案例） ----------

    /// 生成确定性高熵 ACGT 序列（LCG），确保每个 25-mer 窗口熵 ≥ 1.5。
    fn high_entropy(seed: u64, len: usize) -> Vec<u8> {
        let mut s = seed;
        let mut v = Vec::with_capacity(len);
        for _ in 0..len {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            v.push(b"ACGT"[(s >> 33) as usize % 4]);
        }
        // 过滤：理论上不应出现低熵窗口，保守断言
        for j in 0..=v.len().saturating_sub(K) {
            assert!(
                compute_entropy(&v[j..j + K]) >= 1.5,
                "seed {seed} 窗口 {j} 熵 < 1.5"
            );
        }
        v
    }

    fn seq(name: &str, seq: Vec<u8>) -> DnaSeq {
        DnaSeq {
            name: name.to_string(),
            seq,
        }
    }

    /// 基本正向映射：read = bundle0 中段 → pct=100 → 行含组件号 7。
    #[test]
    fn assigns_read_to_bundle() {
        let b0 = high_entropy(1, 60);
        let read = b0[10..45].to_vec(); // 35nt → num_kmer_pos=11
        let bundles = vec![seq(">s_7_1", b0.clone())];
        let reads = vec![seq(">r1", read)];
        let out = reads_to_transcripts(
            &reads,
            &bundles,
            &RttParams {
                pct_required: 50,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(out.mapped_count, 1);
        let line = out.text.lines().next().unwrap();
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f[0], "7");
        assert_eq!(f[1], ">r1");
        // 11 命中单组 = 末组 → run=9；9/11 → 81.8+0.5 → 82%（怪癖）
        assert_eq!(f[2], "82%");
        assert_eq!(f[3].as_bytes(), &b0[10..45]);
    }

    /// 跨 'X' 的 k-mer 不入表：read 与 bundle 仅在含 X 的段重叠 → 不映射。
    #[test]
    fn kmers_spanning_x_are_excluded() {
        let left = high_entropy(2, 30);
        let right = high_entropy(3, 30);
        let mut joined = left.clone();
        joined.push(b'X');
        joined.extend_from_slice(&right);
        // read 横跨 X（两侧各 20nt）：任何 25-mer 都含 X → 表 miss
        let read = joined[15..45].to_vec();
        let bundles = vec![seq(">s_0_1", joined)];
        let reads = vec![seq(">r1", read)];
        let out = reads_to_transcripts(
            &reads,
            &bundles,
            &RttParams {
                pct_required: 50,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(out.mapped_count, 0);
        assert!(out.text.is_empty());
    }

    /// DS（!strand）：revcomp 命中合并进同一 comp；strand=true 时同 read 不映射。
    #[test]
    fn ds_revcomp_hits_merge() {
        let b0 = high_entropy(4, 60);
        let read = revcomp(&b0[10..45]); // 反向互补 → 仅 revcomp 枚举命中
        let bundles = vec![seq(">s_3_1", b0.clone())];
        let reads = vec![seq(">r1", read.clone())];
        let params = RttParams {
            pct_required: 50,
            ..Default::default()
        };
        let ds = reads_to_transcripts(&reads, &bundles, &params).unwrap();
        assert_eq!(ds.mapped_count, 1, "DS 应经 revcomp 命中映射");
        assert!(ds.text.starts_with("3\t>r1\t"));

        let ss = reads_to_transcripts(
            &reads,
            &bundles,
            &RttParams {
                strand: true,
                ..params.clone()
            },
        )
        .unwrap();
        assert_eq!(ss.mapped_count, 0, "SS 下无正向命中 → 不映射");
    }

    /// DS 分子含 revcomp 命中、分母仅正向位置数：read = rc(b0[10..45])
    /// （35nt）—— 正向 0 命中、revcomp 枚举 11 命中 → 单末组 run=9，
    /// 9/11 → 82%（若分母含双向 22 位置则 41% < 50 被拒——区分点）。
    #[test]
    fn ds_denominator_forward_only() {
        let b0 = high_entropy(5, 70);
        let read = revcomp(&b0[10..45]); // 35nt, revcomp 枚举 11 命中
        let bundles = vec![seq(">s_0_1", b0)];
        let out = reads_to_transcripts(
            &[seq(">r1", read)],
            &bundles,
            &RttParams {
                pct_required: 50,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(out.mapped_count, 1);
        assert!(out.text.contains("\t82%\t"), "9/11 → 82%（分母仅正向）");
    }

    /// 共享 25-mer 后写者赢：两 bundle 共享 read 的唯一命中窗口 → 归 bundle1。
    #[test]
    fn shared_kmer_last_writer_wins() {
        // shared 25-mer 唯一窗口 + 各自独有背景
        let shared = high_entropy(6, K);
        let mut b0 = high_entropy(7, K); // 独有段（不同 seed）
        b0.extend_from_slice(&shared);
        let mut b1 = high_entropy(8, K);
        b1.extend_from_slice(&shared);
        let read = shared.clone(); // 25nt：唯一窗口 = shared
        let bundles = vec![seq(">s_0_1", b0), seq(">s_9_1", b1)];
        let out = reads_to_transcripts(
            &[seq(">r1", read)],
            &bundles,
            &RttParams {
                pct_required: 50,
                ..Default::default()
            },
        )
        .unwrap();
        // comp = [1]（bundle0 写 0 后被 bundle1 覆写）；大小 1 组不成 best！
        // → best=-1 不映射——恰好复刻原版怪癖
        assert_eq!(out.mapped_count, 0);
    }

    /// 3 连命中（comp=[1,1,1] 末组 run=1）→ 映射 bundle1
    /// （组件 5），证明共享窗口被 bundle1 覆写。
    #[test]
    fn shared_kmer_overwritten_assigns_to_second_bundle() {
        let shared = high_entropy(6, K);
        let unique_a = high_entropy(9, K); // 仅 bundle1
        let unique_b = high_entropy(10, K); // 仅 bundle1
        let mut b0 = high_entropy(7, K);
        b0.extend_from_slice(&shared);
        let mut b1 = unique_a.clone();
        b1.extend_from_slice(&shared);
        b1.extend_from_slice(&unique_b);
        // read = uniqueA + shared + uniqueB 是 bundle1 的精确子串 → 51 命中
        let mut read = unique_a;
        read.extend_from_slice(&shared);
        read.extend_from_slice(&unique_b);
        let bundles = vec![seq(">s_0_1", b0), seq(">s_5_1", b1)];
        let out = reads_to_transcripts(
            &[seq(">r1", read)],
            &bundles,
            &RttParams {
                pct_required: 50,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(out.mapped_count, 1);
        // 单末组 run=49；49/51 → 96.6+0.5 → 96%
        assert!(
            out.text.starts_with("5\t>r1\t96%\t"),
            "共享 k-mer 归 bundle1"
        );
    }

    /// 行序：multimap(best, read_idx) —— best 升序、同 best 内 read 序升序。
    #[test]
    fn output_ordering_by_best_then_read_index() {
        let b0 = high_entropy(13, 40);
        let b1 = high_entropy(14, 40);
        let bundles = vec![seq(">s_0_1", b0.clone()), seq(">s_2_1", b1.clone())];
        // read1/2 → bundle0，read0 → bundle1：输出 b0 行在前、行内 r1 先于 r2
        let reads = vec![
            seq(">r0", b1[5..35].to_vec()),
            seq(">r1", b0[5..35].to_vec()),
            seq(">r2", b0[8..38].to_vec()),
        ];
        let out = reads_to_transcripts(
            &reads,
            &bundles,
            &RttParams {
                pct_required: 50,
                ..Default::default()
            },
        )
        .unwrap();
        let names: Vec<&str> = out
            .text
            .lines()
            .map(|l| l.split('\t').nth(1).unwrap())
            .collect();
        assert_eq!(names, vec![">r1", ">r2", ">r0"]);
    }

    // ---------- 并行查询（-t）: 确定性 + 宽松加速 ----------

    /// 大输入下 -t1 与 -t4 输出逐字节一致（确定性）, 且 -t4 不劣于 -t1×1.5
    /// （宽松时间断言; 读查询为主导工作量时应有净加速）。
    #[test]
    fn parallel_threads_match_serial_and_speed_up() {
        // 种子逐个尝试直到所有窗口熵达标（high_entropy 内部有断言）。
        fn gen(seed: u64, len: usize) -> Vec<u8> {
            let mut s = seed;
            loop {
                match std::panic::catch_unwind(|| high_entropy(s, len)) {
                    Ok(v) => return v,
                    Err(_) => s += 1,
                }
            }
        }
        let mut bundles = Vec::new();
        for b in 0..4u64 {
            bundles.push(seq(&format!(">s_{b}_1"), gen(100 + b, 2000)));
        }
        let mut reads = Vec::new();
        for i in 0..2000usize {
            let b = &bundles[i % 4].seq;
            let off = (i * 37) % (b.len() - 300);
            let mut r = b[off..off + 300].to_vec();
            if i % 3 == 0 {
                r = revcomp(&r); // 混入反向互补读（DS 路径）
            }
            reads.push(seq(&format!(">r{i}"), r));
        }
        let base = RttParams {
            pct_required: 50,
            ..Default::default()
        };
        let t0 = std::time::Instant::now();
        let s1 = reads_to_transcripts(
            &reads,
            &bundles,
            &RttParams {
                threads: 1,
                ..base.clone()
            },
        )
        .unwrap();
        let d1 = t0.elapsed();
        let t0 = std::time::Instant::now();
        let s4 = reads_to_transcripts(&reads, &bundles, &RttParams { threads: 4, ..base }).unwrap();
        let d4 = t0.elapsed();
        assert_eq!(s1, s4, "-t1 与 -t4 输出必须逐字节一致");
        assert!(s1.mapped_count > 0);
        eprintln!(
            "rtt bench: -t1 {d1:?} vs -t4 {d4:?} ({} reads)",
            reads.len()
        );
        assert!(
            d4 <= d1 * 3 / 2,
            "parallel should not be much slower: -t1 {d1:?} vs -t4 {d4:?}"
        );
    }

    // ---------- sort_reads_to_components ----------

    #[test]
    fn sort_numeric_then_pct_desc_then_name() {
        let text = "10\t>b\t50%\tAAA\n2\t>a\t50%\tAAA\n2\t>c\t80%\tAAA\n2\t>a\t99%\tAAA\n";
        assert_eq!(
            sort_reads_to_components(text),
            "2\t>a\t99%\tAAA\n2\t>c\t80%\tAAA\n2\t>a\t50%\tAAA\n10\t>b\t50%\tAAA\n"
        );
    }

    /// 全键相等 → last-resort 整行字节序（GNU 无 -s）。
    #[test]
    fn sort_last_resort_whole_line() {
        let text = "1\t>a\t50%\tTTT\n1\t>a\t50%\tAAA\n";
        assert_eq!(
            sort_reads_to_components(text),
            "1\t>a\t50%\tAAA\n1\t>a\t50%\tTTT\n"
        );
    }

    /// 数值键：`45%` 取前导 45；无数字字段按 0；前导空白/符号容忍。
    #[test]
    fn sort_numeric_key_quirks() {
        assert_eq!(num_key("45%"), 45);
        assert_eq!(num_key("100%"), 100);
        assert_eq!(num_key(""), 0);
        assert_eq!(num_key("-3"), -3);
        assert_eq!(num_key("  7x"), 7);
        // 两行第 1 键数值均为 7（"07" 与 "7"）、第 2 键 ">a" < ">b" → a 行在前；
        // 数值相等的不同表示不回退字符串比较
        let text = "07\t>b\t50%\tA\n7\t>a\t50%\tA\n";
        assert_eq!(
            sort_reads_to_components(text),
            "7\t>a\t50%\tA\n07\t>b\t50%\tA\n"
        );
    }

    #[test]
    fn sort_empty_text() {
        assert_eq!(sort_reads_to_components(""), "");
    }
}
