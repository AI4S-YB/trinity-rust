//! GraphFromFasta 移植（P3-T2）——焊接聚类：把共享 24-mer 且 read 支持
//! 跨越接缝的 inchworm contig 连成 weld 图。
//! 镜像 `Chrysalis/analysis/GraphFromFasta.cc`（Phase 1 :1304-1403、
//! weldmer 计数 :1405-1427、Phase 2 :1435-1754、report :672-747）。
//!
//! **输出序契约**：原版 report 用非稳定 `std::sort`（仅按 pool size 升序）+
//! OMP 并行插入 map（`map<int,Pool>` 本身按 id 有序）——本版取确定性序：
//! 按 (pool_size 升序, pool_id 升序, 成员插入序)。行序不是下游契约
//! （BubbleUpClustering 消费前会过 `sort -k9,9gr`，见 [`sort_weld_graph`]）。
//!
//! **Phase 2 单线程决策**：原版 Phase 2 OMP 并行下 `toasted` 无锁读 + 中途
//! toast 会影响同轮其它线程的判定（数据竞争 + 结果与调度相关）。本版按
//! 串行语义实现（i 升序 → j 升序 → FW 后 RC → 命中序），结果 = 原版单线程
//! 运行的确定值；rayon 并行化留作 T8 基准后的优化项。
//!
//! **add_scaffolds_to_clusters（:903-995）不移植**——P3 阶段无 PE scaffolding
//! 输入，`scaff_pairs` 恒 0（占位注释，语义与原版无 scaffolding 文件时一致）。

use std::collections::{BTreeMap, BTreeSet, HashMap};

use rayon::prelude::*;

use crate::dna_vector::{is_simple, read_fasta, revcomp, simple_halves_with, DnaSeq};
use crate::kmer_align::KmerAlignCore;
use crate::nonred_table::NonRedKmerTable;
use trinity_common::error::CommonError;

/// GraphFromFasta.cc:29 `static int TOO_SIMILAR = 97;`
pub const TOO_SIMILAR: f64 = 97.;

/// 参数束（Trinity:2180 的实际调用形态：-min_contig_length 200 -min_glue 2
/// -glue_factor 0.05 -min_iso_ratio 0.05 -k 24 -kk 48，非 SS 无 -strand）。
#[derive(Debug, Clone)]
pub struct GffParams {
    /// 24-mer 池化窗口（原版 :1252-1258 强制 k=24——12-mer 索引的 2×12 机制）
    pub k: usize,
    /// weldmer 长（k + 两侧各 (kk-k)/2 flank）
    pub kk: usize,
    /// strand-specific 模式（SS 时为 true：不做 RC 匹配）
    pub strand: bool,
    pub glue_factor: f64,
    pub min_glue_required: u32,
    /// -1 关闭上限（默认）
    pub max_glue_required: i64,
    pub min_iso_ratio: f64,
    pub no_welds: bool,
    pub no_glue_required: bool,
    pub disable_repeat_check: bool,
    pub report_welds: bool,
    pub debug: bool,
    /// `-t`：Phase1 候选收集的显式 rayon 池大小（原版 omp_set_num_threads）。
    /// 收集按 contig 索引序合并 → 与线程数无关（确定性）。
    pub threads: usize,
}

impl Default for GffParams {
    fn default() -> Self {
        GffParams {
            k: 24,
            kk: 48,
            strand: false,
            glue_factor: 0.05,
            min_glue_required: 2,
            max_glue_required: -1,
            min_iso_ratio: 0.05,
            no_welds: false,
            no_glue_required: false,
            disable_repeat_check: false,
            report_welds: false,
            debug: false,
            threads: 1,
        }
    }
}

// ---------------------------------------------------------------------------
// Coverage（:396-408）与 IsGoodCoverage（:410-425）
// ---------------------------------------------------------------------------

/// strtod 最长数字前缀（十进制 + 可选指数；不接受 inf/nan/hex——iworm 名中
/// 不可能出现，且 C++ 语义下 `nan < 1.0` 为假会原样放行 NaN，无业务意义）。
fn strtod_prefix(s: &str) -> Option<f64> {
    let b = s.as_bytes();
    let mut i = 0;
    while i < b.len() && (b[i] as char).is_whitespace() {
        i += 1;
    }
    let start_num = i;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let mut saw_digit = false;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
        saw_digit = true;
    }
    if i < b.len() && b[i] == b'.' {
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
            saw_digit = true;
        }
    }
    if !saw_digit {
        return None;
    }
    let mut end = i;
    if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        let mut j = i + 1;
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
            j += 1;
        }
        let mut exp_digits = false;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
            exp_digits = true;
        }
        if exp_digits {
            end = j;
        }
    }
    let _ = start_num;
    s[..end].parse::<f64>().ok()
}

/// GraphFromFasta.cc:396-408 `Coverage(name)`：name 形如
/// `>a1;43_total_counts:_123_...`——首个 `';'` 之后 strtod；未命中或 `< 1.0`
/// → 1.0。
pub fn coverage(name: &str) -> f64 {
    let Some(semi) = name.find(';') else {
        return 1.0;
    };
    match strtod_prefix(&name[semi + 1..]) {
        Some(v) if v >= 1.0 => v,
        _ => 1.0,
    }
}

/// GraphFromFasta.cc:410-425 `IsGoodCoverage`：min/max > min_iso_ratio（**严格 >**）。
pub fn is_good_coverage(a: f64, b: f64, min_iso_ratio: f64) -> bool {
    let (lo, hi) = if a > b { (b, a) } else { (a, b) };
    lo / hi > min_iso_ratio
}

// ---------------------------------------------------------------------------
// IsShadow（:283-324）/ encapsulates（:328-347）/ align_get_per_id（:351-386）
// ---------------------------------------------------------------------------

/// 阴影 contig 判定（近重复、源于含测序错误的 k-mer）。
///
/// 逐字符比对 `a[startA..]` 与 `b[startB..]` 的错配；连续 ≥3 错配（当前
/// 错配的前两位也错）→ break；`dist==k+1(25)` 的错配计 n、其它计 nn（首个
/// 错配只记 last）；`expect = (int)(0.9 * (len/(k+1) - 1))`——**len/(k+1)
/// 是整数除法**；判 `n >= expect && n > 4 && nn < n/5`（**n/5 整数除法**）。
pub fn is_shadow(a: &[u8], b: &[u8], start_a: usize, start_b: usize, k: usize) -> bool {
    let mut n = 0i64;
    let mut nn = 0i64;
    let mut last: i64 = -1;
    let mut len = 0i64;
    let mut i = start_a;
    while i < a.len() {
        let x = i - start_a + start_b;
        if x >= b.len() {
            break;
        }
        len += 1;
        if a[i] != b[x] {
            if last >= 0 {
                let dist = i as i64 - last;
                // 原版 x>3 && i>3（int 比较，恒真于 startA≥0 的正常路径）
                if x > 3 && i > 3 && a[i - 1] != b[x - 1] && a[i - 2] != b[x - 2] {
                    break;
                }
                if dist == (k as i64 + 1) {
                    n += 1;
                } else {
                    nn += 1;
                }
            }
            last = i as i64;
        }
        i += 1;
    }
    // (int)(0.9 * (double)(len/(k+1) - 1))：整数除法先于 -1 与乘法
    let expect = (0.9f64 * (len / (k as i64 + 1) - 1) as f64) as i64;
    n >= expect && n > 4 && nn < n / 5
}

/// GraphFromFasta.cc:328-347 `encapsulates(largerA, smallerB, startA, startB)`：
/// 以 k-mer 匹配锚定的线性包含判定 `startA > startB &&
/// (startA - startB) + smaller.len() < larger.len()`。
pub fn encapsulates(larger: &[u8], smaller: &[u8], start_a: i64, start_b: i64) -> bool {
    start_a > start_b && ((start_a - start_b) + smaller.len() as i64) < larger.len() as i64
}

/// GraphFromFasta.cc:351-386 `align_get_per_id`：锚定对齐（无 gap 对角线）
/// 的百分比一致率。原版用 **float** 计算（`(len-mismatch)/float(len) * 100`）
/// ——恰 97%（3/100）在 f32 下 == 97.0f 严格 > 判假；本版逐位复刻 f32 语义。
/// len==0 时 C++ 得 NaN（`> 97` 为假）——Rust 0/0 同为 NaN，比较一致为 false。
pub fn align_get_per_id(a: &[u8], b: &[u8], start_a: usize, start_b: usize, _k: usize) -> f64 {
    let (mut start_a, mut start_b) = (start_a as i64, start_b as i64);
    if start_a < start_b {
        start_b -= start_a;
        start_a = 0;
    } else {
        start_a -= start_b;
        start_b = 0;
    }
    let mut len = 0i64;
    let mut mismatch = 0i64;
    let mut i = start_a;
    while i < a.len() as i64 {
        let x = (i - start_a + start_b) as usize;
        if x >= b.len() {
            break;
        }
        len += 1;
        if a[i as usize] != b[x] {
            mismatch += 1;
        }
        i += 1;
    }
    ((len - mismatch) as f32 / len as f32 * 100.0f32) as f64
}

// ---------------------------------------------------------------------------
// Welder（:466-595）
// ---------------------------------------------------------------------------

/// GraphFromFasta.cc:481-518 `WeldableKmer`：flank=(kk-k)/2；
/// `startA=one-flank; stopA=one+k; startB=two+k; stopB=startB+flank`；
/// `startA<0 || stopB >= b.len()`（**>= 非 >**）→ 空（越界拒绝）；
/// 否则 `a[startA..stopA] ++ b[startB..stopB]`（kk=48）。
pub fn weldable_kmer(a: &[u8], one: i64, b: &[u8], two: i64, k: usize, kk: usize) -> Vec<u8> {
    let flank = ((kk - k) / 2) as i64;
    let start_a = one - flank;
    let stop_a = one + k as i64;
    let start_b = two + k as i64;
    let stop_b = start_b + flank;
    if start_a < 0 || stop_b >= b.len() as i64 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(kk);
    out.extend_from_slice(&a[start_a as usize..stop_a as usize]);
    out.extend_from_slice(&b[start_b as usize..stop_b as usize]);
    out
}

/// GraphFromFasta.cc:521-577 `Weldable`：构造 weldmer，越界 → false；
/// low-complexity（IsSimple/SimpleHalves，后者受 disable_repeat_check 控制）
/// → false；`thresh==0` → true（只要 k-mer 匹配）；否则表计数 ≥ thresh。
/// 返回 (通过与否, weldmer 字符串, 计数)——后两者镜像原版出参（即便拒绝
/// 也带着最后一次构造的值返回，供 debug 输出）。
#[allow(clippy::too_many_arguments)] // 镜像原版单函数 8 参形态
pub fn weldable(
    table: &NonRedKmerTable,
    a: &[u8],
    one: i64,
    b: &[u8],
    two: i64,
    thresh: i32,
    k: usize,
    kk: usize,
    disable_repeat_check: bool,
) -> (bool, String, u32) {
    let d = weldable_kmer(a, one, b, two, k, kk);
    if d.is_empty() {
        return (false, String::new(), 0);
    }
    let welding_kmer = String::from_utf8(d.clone()).unwrap();
    if is_simple(&d) || simple_halves_with(&d, disable_repeat_check) {
        return (false, welding_kmer, 0);
    }
    let count = table.get_count(&d, 0);
    if thresh == 0 {
        return (true, welding_kmer, count as u32);
    }
    (count >= thresh, welding_kmer, count as u32)
}

// ---------------------------------------------------------------------------
// 主流程
// ---------------------------------------------------------------------------

/// `make_iworm_pair_token`（:1017-1033）：`"min^max"`。
fn pair_token(a: usize, b: usize) -> String {
    format!("{}^{}", a.min(b), a.max(b))
}

/// 文件入口便捷封装：读 iworm fasta（vecDNAVector::Read 语义）。
pub fn graph_from_fasta_files(
    iworm_fasta: &std::path::Path,
    reads_fasta: &std::path::Path,
    p: &GffParams,
) -> Result<String, CommonError> {
    let iworm = read_fasta(iworm_fasta)?;
    let text = std::fs::read_to_string(reads_fasta)?;
    // DNAStringStreamFast 语义：首 header 前丢弃、空记录截断（dna_vector 模块）
    let reads = crate::dna_vector::stream_fasta_records(&text);
    graph_from_fasta(&iworm, &reads, p)
}

/// 焊接聚类主入口。`iworm_seqs`：vecDNAVector::Read 语义读入的 contig
/// （name 含 '>' 前缀）；`reads`：DNAStringStreamFast 流式语义的 read 集。
/// 返回 weld 图文本（stdout 等价，行序见模块文档的确定性契约）。
pub fn graph_from_fasta(
    iworm_seqs: &[DnaSeq],
    reads: &[Vec<u8>],
    p: &GffParams,
) -> Result<String, CommonError> {
    let k = 24usize; // :1252-1258 强制 k=24
    let kk = p.kk;
    let mut min_glue_required = p.min_glue_required;
    let mut glue_factor = p.glue_factor;
    // :1217-1221 min_glue < 1 → 关闭 read 需求
    if min_glue_required < 1 {
        min_glue_required = 0;
        glue_factor = 0.0;
    }

    // :1282-1287 KmerAlignCore 解码
    let core = KmerAlignCore::build(iworm_seqs);
    let n = iworm_seqs.len();

    // ---- Phase 1（:1304-1403）：收集候选 weldmers ----
    // 原版 omp critical push_back——本版 per-contig 局部收集 + 按索引序合并
    // （表构造前会排序去重，顺序只影响日志；仍取确定性）。
    let crossover: Vec<Vec<u8>> = if p.no_welds {
        Vec::new()
    } else {
        // 显式 `-t` 池（原版 -t1 时也不该吃满全局核）。
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(p.threads.max(1))
            .build()
            .map_err(|e| CommonError::Parse(format!("cannot build GFF thread pool: {e}")))?;
        let per_contig: Vec<Vec<Vec<u8>>> = pool.install(|| {
            (0..n)
                .into_par_iter()
                .map(|i| {
                    let d = &iworm_seqs[i].seq;
                    let mut local: Vec<Vec<u8>> = Vec::new();
                    if d.len() >= k {
                        for j in 0..=d.len() - k {
                            let sub = &d[j..j + k];
                            // :1325-1326 Phase 1 无条件跳过 low-complexity
                            if is_simple(sub) {
                                continue;
                            }
                            let matches_fw = core.get_matches(sub);
                            let matches_rc = if p.strand {
                                Vec::new()
                            } else {
                                core.get_matches(&revcomp(sub))
                            };
                            for (rc, matches) in [(false, &matches_fw), (true, &matches_rc)] {
                                for m in matches.iter() {
                                    let c = m.contig as usize;
                                    if c == i {
                                        continue;
                                    }
                                    // RC 命中：坐标换算到 revcomp(dna[c]) 系
                                    // （:1370-1374；FW 直接用 hit pos）
                                    let (dd, start): (Vec<u8>, i64) = if rc {
                                        let dd = revcomp(&iworm_seqs[c].seq);
                                        let start = dd.len() as i64 - m.pos as i64 - k as i64;
                                        (dd, start)
                                    } else {
                                        (iworm_seqs[c].seq.clone(), m.pos as i64)
                                    };
                                    // 双向各试一次（:1358-1365）
                                    for (a, one, b, two) in [
                                        (d.as_slice(), j as i64, dd.as_slice(), start),
                                        (dd.as_slice(), start, d.as_slice(), j as i64),
                                    ] {
                                        let add = weldable_kmer(a, one, b, two, k, kk);
                                        if !add.is_empty()
                                            && !is_simple(&add)
                                            && !simple_halves_with(&add, p.disable_repeat_check)
                                        {
                                            local.push(add);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    local
                })
                .collect()
        });
        per_contig.into_iter().flatten().collect()
    };

    // ---- weldmer 计数（:1412-1427）----
    let table = NonRedKmerTable::set_up_templates(&crossover, kk, false);
    let table = table.with_read_counts(reads);

    // ---- Phase 2（:1435-1754）：单线程确定性（见模块文档决策）----
    let mut pools: BTreeMap<usize, Vec<usize>> = BTreeMap::new(); // map<int,Pool>
    let mut toasted: BTreeSet<usize> = BTreeSet::new();
    let mut weld_support: HashMap<String, u32> = HashMap::new();
    let mut iworm_lengths = vec![0i64; n];

    for i in 0..n {
        let d = &iworm_seqs[i].seq;
        // :1457-1458 所有 contig 都记录长度（先于长度判断）
        iworm_lengths[i] = d.len() as i64;
        if d.len() < k {
            continue; // :1460-1467
        }
        let cov_i = coverage(&iworm_seqs[i].name);

        for j in 0..=d.len() - k {
            let sub = &d[j..j + k];
            let matches_fw = core.get_matches(sub);
            let matches_rc = if p.strand {
                Vec::new()
            } else {
                core.get_matches(&revcomp(sub))
            };
            // :1485-1491 豁免：simple 且命中数 > 1（**无条件判断在 Phase 2
            // 变成了合取**——与 Phase 1 不同）
            if is_simple(sub) && matches_fw.len() + matches_rc.len() > 1 {
                continue;
            }

            for (rc, matches) in [(false, &matches_fw), (true, &matches_rc)] {
                for m in matches.iter() {
                    let c = m.contig as usize;
                    if c == i {
                        continue; // :1501-1503 / :1635-1637
                    }
                    if toasted.contains(&c) {
                        continue; // :1505-1508 / :1640-1643
                    }
                    let cov_c = coverage(&iworm_seqs[c].name);
                    // :1518-1525 minCov 截断 + clamp
                    let higher = if cov_i > cov_c { cov_i } else { cov_c };
                    let mut min_cov = (higher * glue_factor) as i32; // (int) 截断
                    if (min_cov as i64) < min_glue_required as i64 {
                        min_cov = min_glue_required as i32;
                    }
                    if p.max_glue_required > 0 && min_cov as i64 > p.max_glue_required {
                        min_cov = p.max_glue_required as i32;
                    }

                    let (dd, start): (Vec<u8>, i64) = if rc {
                        let dd = revcomp(&iworm_seqs[c].seq);
                        let start = dd.len() as i64 - m.pos as i64 - k as i64;
                        (dd, start)
                    } else {
                        (iworm_seqs[c].seq.clone(), m.pos as i64)
                    };

                    // :1538-1607（FW）/ :1676-1742（RC）同型体
                    let mut welding_kmer = String::new();
                    let mut welding_count = 0u32;
                    if p.no_glue_required {
                        add_reciprocal_iworm_link(
                            &mut pools,
                            &mut weld_support,
                            i,
                            c,
                            welding_count,
                        );
                    } else {
                        if !p.no_welds {
                            // `!bNoWeld && !(W1 || W2)` —— 短路语义下 welding
                            // 出参取最后一次实际调用（:1552-1554）
                            let (ok1, k1, c1) = weldable(
                                &table,
                                d,
                                j as i64,
                                &dd,
                                start,
                                min_cov,
                                k,
                                kk,
                                p.disable_repeat_check,
                            );
                            welding_kmer = k1;
                            welding_count = c1;
                            let ok = ok1 || {
                                let (ok2, k2, c2) = weldable(
                                    &table,
                                    &dd,
                                    start,
                                    d,
                                    j as i64,
                                    min_cov,
                                    k,
                                    kk,
                                    p.disable_repeat_check,
                                );
                                welding_kmer = k2;
                                welding_count = c2;
                                ok2
                            };
                            if !ok {
                                continue; // :1561
                            }
                        }
                        if is_shadow(d, &dd, j, start as usize, k) && cov_i > 2.0 * cov_c {
                            // :1570-1577 toast 阴影（较小覆盖方）
                            toasted.insert(c);
                            continue;
                        } else if encapsulates(d, &dd, j as i64, start)
                            && d.len() / 10 > dd.len() // 整数除法（:1581）
                            && align_get_per_id(d, &dd, j, start as usize, k) > TOO_SIMILAR
                        {
                            // :1579-1592 toast 被包含者
                            toasted.insert(c);
                            continue;
                        } else if min_glue_required > 0
                            && !is_good_coverage(cov_i, cov_c, p.min_iso_ratio)
                        {
                            continue; // :1594-1602
                        }
                        add_reciprocal_iworm_link(
                            &mut pools,
                            &mut weld_support,
                            i,
                            c,
                            welding_count,
                        );
                    }
                    let _ = (&welding_kmer, p.report_welds, p.debug);
                }
            }
        }
    }

    // ---- report（:672-747）----
    Ok(report_iworm_graph(
        &pools,
        &toasted,
        &iworm_lengths,
        &weld_support,
    ))
}

/// `add_iworm_link`（:998-1014）+ `add_reciprocal_iworm_link`（:1036-1058）：
/// 邻接表双向去重追加 + `pair_token` 历史最大计数。
fn add_reciprocal_iworm_link(
    pools: &mut BTreeMap<usize, Vec<usize>>,
    weld_support: &mut HashMap<String, u32>,
    a: usize,
    b: usize,
    count: u32,
) {
    for (x, y) in [(a, b), (b, a)] {
        pools.entry(x).or_default().push(y);
        // add_iworm_link 的 contains 去重（:1010-1012）——push 前查
        let v = pools.get_mut(&x).unwrap();
        if v[..v.len() - 1].contains(&y) {
            v.pop();
        }
    }
    let token = pair_token(a, b);
    let e = weld_support.entry(token).or_insert(0);
    if *e < count {
        *e = count;
    }
}

/// report_iworm_graph（:672-747）。scaff_pairs 恒 0——P3 无 add_scaffolds
/// 输入（见模块文档占位说明）。
fn report_iworm_graph(
    pools: &BTreeMap<usize, Vec<usize>>,
    toasted: &BTreeSet<usize>,
    iworm_lengths: &[i64],
    weld_support: &HashMap<String, u32>,
) -> String {
    // map<int,Pool> 迭代序 = pool_id 升序；sort_pool_sizes_ascendingly 的
    // 非稳定 sort → 确定性化： (size 升序, pool_id 升序)
    let mut pool_vec: Vec<(&usize, &Vec<usize>)> = pools.iter().collect();
    pool_vec.sort_by_key(|(id, members)| (members.len(), **id));

    let mut out = String::new();
    for (pool_id, members) in pool_vec {
        if toasted.contains(pool_id) {
            continue; // :697-700 跳过 toasted 池主
        }
        // :704-710 成员里 toasted 的剔除（保持剩余序）
        for &member in members {
            if toasted.contains(&member) {
                continue;
            }
            let weld = weld_support
                .get(&pair_token(*pool_id, member))
                .copied()
                .unwrap_or(0);
            let scaff = 0u32; // 无 scaffolding（占位）
            let total = weld + scaff;
            let min_len = iworm_lengths[*pool_id].min(iworm_lengths[member]);
            out.push_str(&format!(
                "{a} -> {b} weldmers: {w} scaff_pairs: {s} total: {t} min_len: {m}\n",
                a = pool_id,
                b = member,
                w = weld,
                s = scaff,
                t = total,
                m = min_len
            ));
        }
    }
    out
}

/// 镜像 `sort -k9,9gr`（Trinity:2192 对 weld 图的下游排序）：按第 9 字段
/// （total 的数值）降序数值排序；**tie 用 GNU sort 缺省的 last-resort
/// 整行字节序升序**（无 `-s` 时 GNU sort 对 key 相等的行回退到整行比较，
/// C locale 即字节序）——BubbleUpClustering 的簇序（进而 component 编号）
/// 逐序依赖这一 tie-break。缺第 9 字段的行按 0 处理（GNU 数值 sort 对
/// 缺失 key 的行为同型）。
pub fn sort_weld_graph(text: &str) -> String {
    let mut lines: Vec<&str> = text.split('\n').collect();
    let trailing_nl = lines.last().is_some_and(|l| l.is_empty());
    if trailing_nl {
        lines.pop();
    }
    lines.sort_by(|x, y| {
        let key = |s: &str| -> f64 {
            s.split_whitespace()
                .nth(8)
                .and_then(|f| f.parse::<f64>().ok())
                .unwrap_or(0.0)
        };
        match key(y).partial_cmp(&key(x)) {
            Some(std::cmp::Ordering::Equal) | None => x.cmp(y), // GNU last-resort 整行字节序
            Some(o) => o,
        }
    });
    let mut out = lines.join("\n");
    if !lines.is_empty() {
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---------- Coverage ----------

    #[test]
    fn coverage_name_forms() {
        // '>a1;43_total_counts:_123' → 43
        assert_eq!(coverage(">a1;43_total_counts:_123"), 43.0);
        // 无 ';' → 1
        assert_eq!(coverage(">a1_total_counts:_123"), 1.0);
        // ';' 后非数字 → 1
        assert_eq!(coverage(">a1;abc"), 1.0);
        // < 1.0 → 1
        assert_eq!(coverage(">a1;0.5"), 1.0);
        // 恰 1.0 放行（< 严格）
        assert_eq!(coverage(">a1;1.0"), 1.0);
        // 小数覆盖
        assert_eq!(coverage(">a1;123.75 xxx"), 123.75);
        // 首个 ';' 即止（后续 ';' 不影响）
        assert_eq!(coverage(">a;7;99"), 7.0);
    }

    // ---------- IsGoodCoverage ----------

    #[test]
    fn is_good_coverage_strict_ratio() {
        // 1:20 = 0.05 恰等 → 拒绝（严格 >）
        assert!(!is_good_coverage(1.0, 20.0, 0.05));
        // 1.01:20 略高 → 通过
        assert!(is_good_coverage(1.01, 20.0, 0.05));
        // 参数序无关（内部 swap）
        assert!(!is_good_coverage(20.0, 1.0, 0.05));
    }

    // ---------- IsShadow ----------

    /// 构造：a 长片段，b 与 a 在对角线上每 25bp 一个 SNP（dist==k+1 → 计 n）。
    /// len = 25×M；n = M-1（首个错配只记 last）；expect = (int)(0.9×(M-1))。
    /// M=30：n=29 > 4、expect=26 → n>=expect、nn=0 < 29/5 → shadow。
    #[test]
    fn is_shadow_periodic_snps() {
        let m = 30usize;
        let a = vec![b'A'; m * 25];
        let mut b = a.clone();
        for t in 0..m {
            b[t * 25] = b'C'; // 错配间距恰 25
        }
        assert!(is_shadow(&a, &b, 0, 0, 24));
    }

    /// 错配间距 26（dist != 25）→ 全计入 nn → 非阴影（nn < n/5 失败，n=0）。
    #[test]
    fn is_shadow_wrong_period_counts_nn() {
        let m = 30usize;
        let a = vec![b'A'; m * 26];
        let mut b = a.clone();
        for t in 0..m {
            b[t * 26] = b'C';
        }
        assert!(!is_shadow(&a, &b, 0, 0, 24));
    }

    /// 连续 ≥3 错配中止：b 尾部三连错后即便周期 SNP 达标也 break（len 截短，
    /// n 不再增长——此场景 n 本就不足）。
    #[test]
    fn is_shadow_triple_mismatch_breaks() {
        let m = 30usize;
        let a = vec![b'A'; m * 25 + 50];
        let mut b = a.clone();
        for t in 0..m {
            b[t * 25] = b'C';
        }
        // 在尾部插入三连错配（位置 m*25 起）
        b[m * 25] = b'G';
        b[m * 25 + 1] = b'G';
        b[m * 25 + 2] = b'G';
        // break 发生在第三连错（前两位也错），n 仍是 29 → 仍判 shadow，
        // 但 len 停在 m*25+3；验证不 panic 且判定与手推一致
        assert!(is_shadow(&a, &b, 0, 0, 24));
    }

    /// 零错配 → n=0 不满足 n>4 → 非阴影。
    #[test]
    fn is_shadow_identical_is_not_shadow() {
        let a = vec![b'A'; 1000];
        assert!(!is_shadow(&a, &a, 0, 0, 24));
    }

    /// 短对齐：len < 25 → len/25 = 0 → expect = (int)(0.9×(-1)) = 0
    /// （向零截断）→ n>=0 恒真但 n>4 需 5 个周期错配 → false。
    #[test]
    fn is_shadow_short_alignment_negative_expect() {
        let a = b"AAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let mut b = a.to_vec();
        b[0] = b'C';
        assert!(!is_shadow(a, &b, 0, 0, 24));
    }

    // ---------- align_get_per_id ----------

    /// 100bp 对齐 3 错配 = 97.0（f32 下恰等，见函数文档）→ 严格 > 97 为假；
    /// 2 错配 = 98 > 97 为真。
    #[test]
    fn per_id_strict_boundary() {
        let a = vec![b'A'; 100];
        let mut b3 = a.clone();
        for p in [0usize, 10, 20] {
            b3[p] = b'C';
        }
        let mut b2 = a.clone();
        for p in [0usize, 10] {
            b2[p] = b'C';
        }
        assert_eq!(align_get_per_id(&a, &b3, 0, 0, 24), 97.0);
        assert!(align_get_per_id(&a, &b3, 0, 0, 24) <= TOO_SIMILAR); // 恰等拒绝
        assert!(align_get_per_id(&a, &b2, 0, 0, 24) > TOO_SIMILAR); // 98 通过
    }

    /// 锚定偏移：b 右移 5（startB=5）→ 对齐去掉 a 尾 5bp，len=95。
    #[test]
    fn per_id_offset_anchors() {
        let a = vec![b'A'; 100];
        let mut b = vec![b'A'; 100];
        b[5] = b'G';
        b[6] = b'G'; // 对齐后 pos1,2 错配
                     // f32 语义（见函数文档）
        assert_eq!(
            align_get_per_id(&a, &b, 0, 5, 24),
            (93f32 / 95f32 * 100f32) as f64
        );
    }

    /// 零长度对齐（b 空）→ NaN > 97 为假（与 C++ 0/0=NaN 一致）。
    #[test]
    fn per_id_zero_len_is_nan_false() {
        let per = align_get_per_id(b"AAAA", b"", 0, 0, 24);
        assert!(per.is_nan());
        assert!(per <= TOO_SIMILAR || per.is_nan());
    }

    // ---------- encapsulates ----------

    #[test]
    fn encapsulates_geometry() {
        // larger=100, smaller=50, startA=40 > startB=10, 30+50=80 < 100 → true
        assert!(encapsulates(&[b'A'; 100], &[b'C'; 50], 40, 10));
        // (startA-startB)+smaller == larger 恰等 → 严格 < 失败 → false
        assert!(!encapsulates(&[b'A'; 100], &[b'C'; 50], 60, 10));
        // startA <= startB → false
        assert!(!encapsulates(&[b'A'; 100], &[b'C'; 50], 10, 40));
    }

    // ---------- WeldableKmer ----------

    /// 48 拼接：a 的 one-12 起 36bp + b 的 two+24 起 12bp。
    #[test]
    fn weldable_kmer_composition() {
        let mut a = vec![b'A'; 24];
        a.extend_from_slice(b"ACGTACGTACGTACGTACGTACGT"); // j=24 的 24-mer
        let mut b = b"TTTTTTTTTTTTTTTTTTTTTTTT".to_vec();
        b.extend_from_slice(&[b'G'; 13]); // len 37：stopB=36 = len-1 恰好通过
        let w = weldable_kmer(&a, 24, &b, 0, 24, 48);
        assert_eq!(w.len(), 48);
        // a[12..48] = 12×A + 24-mer；b[24..36] = G×12
        assert_eq!(&w[..12], &[b'A'; 12]);
        assert_eq!(&w[12..36], b"ACGTACGTACGTACGTACGTACGT");
        assert_eq!(&w[36..], &[b'G'; 12]);
    }

    /// 越界：stopB == b.len() 恰等 → 拒绝（>= 语义）；stopB == b.len()-1 通过。
    #[test]
    fn weldable_kmer_boundaries() {
        let a = vec![b'A'; 48];
        // b 使 stopB = two+24+12 = b.len() 恰等 → 空
        let b_short = vec![b'C'; 36];
        assert!(weldable_kmer(&a, 12, &b_short, 0, 24, 48).is_empty());
        // b.len() = 37 → stopB=36 = len-1 → 通过
        let mut b_ok = vec![b'C'; 37];
        b_ok[35] = b'G'; // 落在 b[24..36) 窗口末位 → w[47]
        let w = weldable_kmer(&a, 12, &b_ok, 0, 24, 48);
        assert_eq!(w.len(), 48);
        assert_eq!(w[47], b'G');
        // startA < 0（one < 12）→ 空
        assert!(weldable_kmer(&a, 11, &b_ok, 0, 24, 48).is_empty());
        // one=12 恰好 startA=0 → 通过
        assert_eq!(weldable_kmer(&a, 12, &b_ok, 0, 24, 48).len(), 48);
    }

    // ---------- minCov 截断与 clamp（借 graph 内部逻辑的独立复算） ----------

    #[test]
    fn min_cov_truncation_and_clamp() {
        let higher = 43.0f64;
        let glue = 0.05f64;
        // (int)(43*0.05) = (int)2.15 = 2（截断非四舍五入）
        assert_eq!((higher * glue) as i32, 2);
        // clamp 下限：higher=10 → 0 → min_glue=2
        assert_eq!((10.0f64 * 0.05) as i32, 0);
        // 上限：higher=1000 → 50 → max_glue=10
        let mut mc = (1000.0f64 * 0.05) as i32;
        if mc > 10 {
            mc = 10;
        }
        assert_eq!(mc, 10);
    }

    // ---------- 端到端小场景 ----------

    fn mkseq(name: &str, seq: &[u8]) -> DnaSeq {
        DnaSeq {
            name: name.to_string(),
            seq: seq.to_vec(),
        }
    }

    /// 3 contig：c0/c1 共享 24-mer 且两侧 flank 齐备的 weldmer 在 reads 出现
    /// 2 次（≥ minCov=2）；c2 独立。断言边集合 = {0-1}。
    #[test]
    fn end_to_end_small_weld() {
        // 随机 48bp 母序列（熵 2.0、两半非重复）
        let head: &[u8] = b"ACGTTGCAACGTTGCAGGCATTAC"; // 24
        let mid: &[u8] = b"TTGCAAGGCATTACACGTTGCAAG"; // 24：c0[24..48] = c1[0..24]
        let tail: &[u8] = b"CAAGGCATTACGTTGCAACGTTGC"; // 24
        let c0 = [head, mid].concat(); // 48
        let c1 = [mid, tail].concat(); // 48
        let mut c2 = vec![b'G'; 48];
        c2[0] = b'A';
        c2[24] = b'T';
        let iworm = vec![
            mkseq(">a0;20", &c0),
            mkseq(">a1;20", &c1),
            mkseq(">a2;20", &c2),
        ];
        // 共享 24-mer = mid（c0 j=24 / c1 two=0）。
        // minCov = (int)(20*0.05) = 1 → clamp 到 min_glue=2。
        // weldmer = c0[12..48] ++ c1[24..36]（startA=12, stopB=36 < 48 不越界）
        let weldmer = [&c0[12..48], &c1[24..36]].concat();
        assert_eq!(weldmer.len(), 48);
        let reads = vec![weldmer.clone(), weldmer.clone(), c2.clone()];
        let out = graph_from_fasta(&iworm, &reads, &GffParams::default()).unwrap();
        // 期望边：0-1 双向（邻接表互指）
        let lines: Vec<&str> = out.lines().collect();
        let edges: std::collections::BTreeSet<(usize, usize)> = lines
            .iter()
            .filter_map(|l| {
                let f: Vec<&str> = l.split_whitespace().collect();
                Some((f[0].parse().ok()?, f[2].parse().ok()?))
            })
            .collect();
        assert_eq!(edges, [(0usize, 1usize), (1, 0)].into_iter().collect());
        // c2 无边
        assert!(!edges.iter().any(|e| e.0 == 2 || e.1 == 2));
        // weld 计数 = 2（表计数 ≥ minCov=2 的 weldmer）
        assert!(lines
            .iter()
            .all(|l| l.contains("weldmers: 2") || l.contains("weldmers: 0")));
        assert!(lines.iter().any(|l| l.contains("weldmers: 2")));
        // min_len = 48
        assert!(lines.iter().all(|l| l.ends_with("min_len: 48")));
    }

    /// weldmer 只出现 1 次（< minCov=2）→ 无边。
    #[test]
    fn end_to_end_below_min_glue_no_edge() {
        let head: &[u8] = b"ACGTTGCAACGTTGCAGGCATTAC";
        let mid: &[u8] = b"TTGCAAGGCATTACACGTTGCAAG";
        let tail: &[u8] = b"CAAGGCATTACGTTGCAACGTTGC";
        let c0 = [head, mid].concat();
        let c1 = [mid, tail].concat();
        let iworm = vec![mkseq(">a0;20", &c0), mkseq(">a1;20", &c1)];
        let weldmer = [&c0[12..48], &c1[24..36]].concat();
        let reads = vec![weldmer]; // 只 1 次
        let out = graph_from_fasta(&iworm, &reads, &GffParams::default()).unwrap();
        assert_eq!(out, "");
    }

    /// no_glue_required：只要共享 24-mer 即连边（无 read）。
    #[test]
    fn end_to_end_no_glue_required() {
        let head: &[u8] = b"ACGTTGCAACGTTGCAGGCATTAC";
        let mid: &[u8] = b"TTGCAAGGCATTACACGTTGCAAG";
        let tail: &[u8] = b"CAAGGCATTACGTTGCAACGTTGC";
        let c0 = [head, mid].concat();
        let c1 = [mid, tail].concat();
        let iworm = vec![mkseq(">a0;100", &c0), mkseq(">a1;100", &c1)];
        let out = graph_from_fasta(
            &iworm,
            &[],
            &GffParams {
                no_glue_required: true,
                ..GffParams::default()
            },
        )
        .unwrap();
        assert_eq!(out.lines().count(), 2); // 0->1 与 1->0，weldmers: 0
        assert!(out.lines().all(|l| l.contains("weldmers: 0")));
    }

    /// 覆盖比悬殊（1 vs 100 → 0.01 < 0.05）→ IsGoodCoverage 拒绝。
    #[test]
    fn end_to_end_coverage_ratio_rejects() {
        let head: &[u8] = b"ACGTTGCAACGTTGCAGGCATTAC";
        let mid: &[u8] = b"TTGCAAGGCATTACACGTTGCAAG";
        let tail: &[u8] = b"CAAGGCATTACGTTGCAACGTTGC";
        let c0 = [head, mid].concat();
        let c1 = [mid, tail].concat();
        // minCov = (int)(20*0.05)=1 → clamp 到 2 → weldmer×2 通过 glue；
        // 但 IsGoodCoverage(20, 1) = 0.05 > 0.05 为假 → 拒绝
        let iworm = vec![mkseq(">a0;20", &c0), mkseq(">a1;1", &c1)];
        let weldmer = [&c0[12..48], &c1[24..36]].concat();
        let reads = vec![weldmer.clone(), weldmer.clone()];
        let out = graph_from_fasta(&iworm, &reads, &GffParams::default()).unwrap();
        assert_eq!(out, "");
    }

    // ---------- sort_weld_graph ----------

    #[test]
    fn sort_weld_graph_numeric_descending() {
        let text = "0 -> 1 weldmers: 5 scaff_pairs: 0 total: 10 min_len: 300\n\
                    2 -> 3 weldmers: 2 scaff_pairs: 0 total: 9 min_len: 500\n\
                    4 -> 5 weldmers: 1 scaff_pairs: 0 total: 100 min_len: 50\n";
        let sorted = sort_weld_graph(text);
        let totals: Vec<f64> = sorted
            .lines()
            .map(|l| l.split_whitespace().nth(8).unwrap().parse().unwrap())
            .collect();
        // 数值降序：100 > 10 > 9（字典序会是 10 > 100 > 9 的陷阱）
        assert_eq!(totals, vec![100.0, 10.0, 9.0]);
    }

    /// tie → GNU last-resort 整行字节序升序（无 `-s`）：`0 ->` < `2 ->`。
    #[test]
    fn sort_weld_graph_gnu_last_resort_ties() {
        let text = "2 -> 3 weldmers: 2 scaff_pairs: 0 total: 7 min_len: 500\n\
                    0 -> 1 weldmers: 5 scaff_pairs: 0 total: 7 min_len: 300\n";
        let sorted = sort_weld_graph(text);
        assert_eq!(
            sorted.lines().next().unwrap(),
            "0 -> 1 weldmers: 5 scaff_pairs: 0 total: 7 min_len: 300"
        );
    }

    #[test]
    fn sort_weld_graph_empty() {
        assert_eq!(sort_weld_graph(""), "");
    }
}
