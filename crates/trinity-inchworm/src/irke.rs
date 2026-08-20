//! IRKE 装载 / 剪枝 / 种子列表 / 贪心延伸核心 — 直译 IRKE.cpp 的 populate_*、
//! KmerCounter.cpp 的 prune_some_kmers / get_kmers_sort_descending_counts，
//! 以及 IRKE.cpp 的贪心组装核心（is_good_seed_kmer / build_inchworm_contig_from_seed /
//! inchworm / inchworm_step / extract_best_seed / reconstruct_path_sequence）。
//!
//! 语义要点（对照原版行号）:
//! - populate_from_kmers（IRKE.cpp:81-154）: jellyfish dump 装载。逐记录，
//!   **空序列记录即终止解析**（IRKE.cpp:107-108 `get_sequence() == ""` break——
//!   原版同时用它实现 EOF 检测，中文件空记录也会提前终止且不计入 parsed）；
//!   `len != K` 严格不等跳过但计入 parsed（record_counter 在校验前递增，
//!   IRKE.cpp:110 → 124）；header 按 C atoi 语义取计数（IRKE.cpp:129）；
//!   非 gatc → Err（原版 kmer_to_intval throw → main catch 退出，IRKE.cpp:128）。
//! - populate_from_reads（IRKE.cpp:157-280，非 reassemble 路径）: hasNext 循环
//!   **无空记录终止**（IRKE.cpp:183-184，与 kmers 模式不同！空记录照常计数，
//!   len < K 跳过，IRKE.cpp:207-209 严格 <——与 kmers 模式的 != 不同）；
//!   每窗口 add_kmer(kmer, 1)（add_sequence 默认 cov=1，KmerCounter.hpp:74）。
//! - add_sequence（KC:34-44 + KC:493-505）: 逐窗 substr(i, K)，含非 gatc 的窗口
//!   被 add_kmer(string) 的 contains_non_gatc 前置 `return(false)` **静默跳过**
//!   （跳过点不在 add_sequence 本身！），窗口继续前滑 1 碱基——滑动语义，
//!   与 P1 counter::for_each_kmer 枚举同一窗口集合。
//! - prune_some_kmers（KC:135-280）: 三段剪枝，见函数文档。
//! - 种子列表（KC:701-752 + IRKE.cpp:1332-1346）: count>0 收集；sort=true 用
//!   _DEBUG 比较器（count 降序 + 平局 kmer 值降序）——原版发布版比较器只有
//!   count 降序 + std::sort 不稳定（平局为 introsort 实现定义序），此处显式化
//!   为确定性全序（计划的选择）。

use std::fs;
use std::io::{BufReader, Write};
use std::path::Path;
use std::time::Instant;

use rayon::prelude::*;
use rustc_hash::FxHashSet;
use trinity_common::error::CommonError;
use trinity_common::fasta::{add_fasta_seq_line_breaks, FastaReader};
use trinity_common::kmer::{compute_entropy, kmer_to_intval, revcomp_val, KmerId};
use trinity_common::seq_hash::generate_hash;

use crate::counter_sync::SyncKmerCounter;
use crate::glibc_rand::GlibcRand;
use crate::kmer_counter::{KmerCatalog, KmerCounter};
use crate::visitor::KmerVisitor;

/// C `atoi`（glibc 实现 = `(int)strtol(s, NULL, 10)`）语义，用于 dump header
/// 计数（IRKE.cpp:129 `unsigned int count = atoi(header.c_str())`）:
/// 跳过前导 C 空白（空格 \t \n \v \f \r）、可选 +/-、十进制数字串在首个非数字
/// 处截断、无数字 = 0；数值按 strtol 的 long 饱和（正 LONG_MAX / 负 LONG_MIN），
/// 再 `(int)` 截断回绕、按位解释为 u32。黄金值由 C harness（glibc atoi）锁定。
pub fn c_atoi(bytes: &[u8]) -> u32 {
    let mut i = 0;
    // C isspace: ' ' \t \n \v \f \r
    while i < bytes.len() && matches!(bytes[i], b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r') {
        i += 1;
    }
    let negative = matches!(bytes.get(i), Some(b'-'));
    if matches!(bytes.get(i), Some(b'+') | Some(b'-')) {
        i += 1;
    }
    let mut magnitude: u128 = 0;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        // 逐步封顶防 u128 溢出（超长数字串）; 语义上只需判"是否超过 long 上限"
        magnitude = (magnitude * 10 + (bytes[i] - b'0') as u128).min(u64::MAX as u128);
        i += 1;
    }
    // strtol 饱和: 正数上限 LONG_MAX = 2^63-1，负数上限 LONG_MIN 幅值 = 2^63
    let limit: u128 = if negative { 1 << 63 } else { (1 << 63) - 1 };
    let saturated = magnitude.min(limit);
    let signed = if negative {
        -(saturated as i128) as i64
    } else {
        saturated as i64
    };
    signed as i32 as u32 // (int) 截断回绕 + 赋值 unsigned int（位型保留）
}

/// populate_Kmers_from_kmers（IRKE.cpp:81-154）—— `--kmers` jellyfish dump 装载。
/// 返回 (counter, parsed)；parsed 含被长度校验跳过的记录（镜像 record_counter
/// 在校验前递增），不含触发终止的空序列记录。非 gatc 碱基 → Err（镜像原版 throw）。
/// 读入后、剪枝前由 CLI 调 write_kmer_count_report（镜像 IRKE.cpp:144-147）。
pub fn populate_from_kmers(
    data: &[u8],
    k: usize,
    ds: bool,
) -> Result<(KmerCounter, usize), CommonError> {
    let mut kcounter = KmerCounter::new(k, ds);
    let mut parsed = 0usize;
    let mut reader = FastaReader::new(BufReader::new(data));
    while let Some(record) = reader.next_record()? {
        // IRKE.cpp:107-108: 空序列记录终止解析（原版 getNext 的 EOF 哨兵同款判定，
        // 中文件空记录同样提前终止，且该记录不计入 parsed——计数在其之后）
        if record.sequence.is_empty() {
            break;
        }
        // IRKE.cpp:110: record_counter 在长度校验前递增（跳过的记录也计数）
        parsed += 1;
        // IRKE.cpp:124-126: 长度严格 != 才跳过（恰好 == K 的记录继续）
        if record.sequence.len() != k {
            continue;
        }
        // IRKE.cpp:128-130: get_kmer_intval（非 gatc → Err，镜像原版 throw）+ atoi 计数；
        // DS 折叠在 add_kmer（KmerCounter.cpp:479）——与原版调用序一致
        let kmer = kcounter.get_kmer_intval(record.sequence.as_bytes())?;
        let count = c_atoi(record.header.as_bytes());
        kcounter.add_kmer(kmer, count);
    }
    Ok((kcounter, parsed))
}

/// populate_Kmers_from_fasta（IRKE.cpp:157-280，非 reassemble 路径）—— `--reads` 模式。
/// 返回 (counter, parsed)；parsed 含空序列与短记录（hasNext 循环无终止哨兵）。
/// reassembleIworm 分支（IRKE.cpp:211-249）与 PRUNE_SINGLETON_READ_INTERVAL
/// （IRKE.cpp:255-264，默认 0 且 CLI 无入口）不在主线，未移植。
pub fn populate_from_reads(
    data: &[u8],
    k: usize,
    ds: bool,
) -> Result<(KmerCounter, usize), CommonError> {
    let mut kcounter = KmerCounter::new(k, ds);
    let mut parsed = 0usize;
    let mut reader = FastaReader::new(BufReader::new(data));
    while let Some(record) = reader.next_record()? {
        // IRKE.cpp:183-189: hasNext/getNext 循环——无空记录终止，空记录照常计数
        parsed += 1;
        // IRKE.cpp:207-209: 长度严格 < 才跳过（与 kmers 模式的 != 不同）
        if record.sequence.len() < k {
            continue;
        }
        add_sequence(&mut kcounter, record.sequence.as_bytes(), 1);
    }
    Ok((kcounter, parsed))
}

/// IRKE.cpp:144-147: 写单行 kcounter.size()（键总数，含 0 值键；读入后剪枝前）。
/// 原版固定写 CWD 的 "inchworm.kmer_count"，此处路径参数化（CLI 层决定落点）。
pub fn write_kmer_count_report(path: &Path, kmer_count: usize) -> Result<(), CommonError> {
    fs::write(path, format!("{kmer_count}\n"))?;
    Ok(())
}

/// prune_some_kmers（KmerCounter.cpp:135-280）—— 三段剪枝，返回 count_pruned
/// （min_count 段 + entropy 段 + 错误 kmer 段之和；错误 kmer 可重复入列，每次都计）。
/// 原版返回 bool(count_pruned > 0) 并打印 "Pruned N kmers from catalog."——
/// 这里返回计数，monitor 日志由 CLI 层输出。
///
/// 遍历域 = 全键含 0 值（C++ map 迭代域），序 = FxHashMap 迭代序（与原版哈希序
/// 不同）。min_count/entropy 段**即时置 0** 并 continue——会影响本趟后续 kmer 的
/// 候选计数查询，该行为原样镜像；遍历序差异只在剪枝集合于遍历中途相交的极端
/// 情况下使最终态偏离原版（两版哈希序本就不同，by-design）。错误 kmer 段的
/// deletion_list 统一在遍历结束后置 0（幂等），镜像 KC:267-273。
pub fn prune_some_kmers(
    counter: &mut KmerCounter,
    min_count: u32,
    min_entropy: f32,
    prune_error_kmers: bool,
    min_ratio_non_error: f32,
) -> usize {
    let k = counter.get_kmer_length();
    // 遍历域 = 全键含 0 值（C++ map 迭代域，KC:143）。剪枝只置值不增删键，
    // 先快照键序等价于原版的活迭代；序为 FxHashMap 迭代序（≠原版哈希序）。
    let keys: Vec<KmerId> = counter.iter_all().map(|(kmer, _)| kmer).collect();
    let mut deletion_list: Vec<KmerId> = Vec::new();
    let mut count_pruned: usize = 0;

    for kmer in keys {
        let count = counter.get_kmer_count(kmer);
        if count == 0 {
            continue; // KC:148-149
        }
        if count < min_count {
            counter.clear_kmer(kmer); // KC:156 即时置 0（影响本趟后续候选查询）
            count_pruned += 1;
            continue;
        }
        if compute_entropy(kmer, k) < min_entropy {
            counter.clear_kmer(kmer); // KC:164 即时置 0
            count_pruned += 1;
            continue;
        }
        if !prune_error_kmers {
            continue;
        }
        // KC:172-258: 前向与反向各一遍，结构对称；dominant = 排序后 candidates[0]
        for candidates in [
            counter.get_forward_kmer_candidates(kmer),
            counter.get_reverse_kmer_candidates(kmer),
        ] {
            if candidates.len() <= 1 {
                continue; // KC:175/221
            }
            let dominant_count = candidates[0].1 as i32; // u32→int（g++ 位型保留）
            for &(candidate_key, candidate_u32) in &candidates[1..] {
                if candidate_u32 == 0 {
                    continue; // KC:183/227（count>0 过滤后恒真，保留镜像结构）
                }
                let candidate_count = candidate_u32 as i32;
                let ratio_dominant_count = candidate_count as f32 / dominant_count as f32;
                let ratio_curr_count = candidate_count as f32 / count as f32;
                if dominant_count > 0
                    && ratio_dominant_count < min_ratio_non_error
                    && ratio_curr_count < min_ratio_non_error
                {
                    // KC:206/250 入列（不即时置 0——207 行的即时置 0 在原版被注释），
                    // 可重复入列（幂等置 0）；DS 下 candidate_key 是原始位运算值，
                    // clear_kmer 内部 canonical 与原版 find_kmer 折叠到同一键
                    deletion_list.push(candidate_key);
                    count_pruned += 1;
                }
            }
        }
    }

    if count_pruned > 0 {
        for kmer in deletion_list {
            counter.clear_kmer(kmer); // KC:267-273: 统一置 0（不缩表）
        }
    }
    count_pruned
}

/// 种子列表（KC:701-752 get_kmers_sort_descending_counts + IRKE.cpp:1332-1338）。
/// 收集 count > 0 的 (kmer, count)。sort=false（PARALLEL，__DEVEL_no_kmer_sort，
/// KC:732-737）保持容器迭代序返回；sort=true（单线程）按 count 降序 + 平局
/// kmer 值降序（KC:823-831 _DEBUG 比较器的显式化——发布版平局为不稳定排序的
/// 实现定义序，此处选定确定性全序；平局键唯一 → 结果唯一）。
pub fn sorted_seed_list(counter: &KmerCounter, sort: bool) -> Vec<(KmerId, u32)> {
    let mut kmer_list = collect_nonzero_seeds(counter);
    if sort {
        sort_seeds_desc(&mut kmer_list);
    }
    kmer_list
}

/// KC:719-728 收集阶段（get_kmers_sort_descending_counts 前半）: count>0 快照，
/// 迭代序。compute_sequence_assemblies 需要在收集与排序之间插入 monitor 行，
/// 故拆出。
pub fn collect_nonzero_seeds(counter: &KmerCounter) -> Vec<(KmerId, u32)> {
    counter.iter_nonzero().collect()
}

/// KC:739-743 排序阶段 + KC:823-831（_DEBUG 比较器）: count 降序，平局 kmer 值降序。
fn sort_seeds_desc(kmer_list: &mut [(KmerId, u32)]) {
    kmer_list.sort_by(|&(kmer_a, count_a), &(kmer_b, count_b)| {
        count_b.cmp(&count_a).then(kmer_b.cmp(&kmer_a))
    });
}

/// KC:34-44 add_sequence + KC:493-505 add_kmer(string): 逐窗 substr(i, K)，
/// 含非 gatc 的窗口 contains_non_gatc → return(false) 静默跳过，窗口继续前滑。
/// 调用方已保证 len >= K（原版此处 unsigned 下溢为 UB，防御性返回）。
fn add_sequence(kcounter: &mut KmerCounter, seq: &[u8], cov: u32) {
    let k = kcounter.get_kmer_length();
    if seq.len() < k {
        return;
    }
    for i in 0..=seq.len() - k {
        // 唯一 Err 原因是非 gatc（窗口长度 == K ≤ 32）→ 镜像 add_kmer(string) 跳过
        if let Ok(kmer) = kmer_to_intval(&seq[i..i + k]) {
            kcounter.add_kmer(kmer, cov);
        }
    }
}

// ===========================================================================
// 贪心延伸核心（IRKE.cpp:719-818, 821-918, 933-1150, 1159-1213, 1348-1371）
// ===========================================================================

/// IRKE_run.cpp:79 `MAX_RECURSION = 1`（inchworm step size，CLI -R 可调）
pub const MAX_RECURSION: u32 = 1;

/// IRKE.cpp:29 `const unsigned int MAX_RECURSION_HARD_STOP = 50;`
pub const MAX_RECURSION_HARD_STOP: u32 = 50;

/// IRKE.hpp:102-105 成员参数的显式化。默认值取 IRKE_run.cpp:79/101/102。
#[derive(Debug, Clone, Copy)]
pub struct IrkeParams {
    /// MIN_CONNECTIVITY_RATIO（IRKE_run.cpp 默认 0 → exceeds_min_connectivity 恒过）
    pub min_connectivity: f32,
    /// MIN_SEED_ENTROPY = 1.5（IRKE_run.cpp:101）
    pub min_seed_entropy: f32,
    /// MIN_SEED_COVERAGE = 2（IRKE_run.cpp:102）
    pub min_seed_coverage: u32,
}

impl Default for IrkeParams {
    fn default() -> Self {
        IrkeParams {
            min_connectivity: 0.0,
            min_seed_entropy: 1.5,
            min_seed_coverage: 2,
        }
    }
}

/// IRKE.hpp:12 `typedef pair<vector<kmer_int_type_t>,int> Path_n_count_pair`
/// —— count 是 **int(i32)**（C++ int 的回绕语义用 wrapping_add 镜像）。
#[derive(Debug, Clone, Default)]
pub struct PathNCount {
    /// 深→浅压栈：path[0] = 最远端，path[len-1] = 紧邻上一步的 kmer
    pub path: Vec<KmerId>,
    pub count: i32,
}

/// inchworm 的延伸方向（原版 char 'F'/'R'，IRKE.cpp:822 等）
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    Forward, // 'F'
    Reverse, // 'R'
}

/// IRKE.cpp:719-760 is_good_seed_kmer。
/// 4 段判定全部按原版顺序、严格比较：count==0 → 回文（kmer == revcomp_val(kmer, K)
/// 精确等值；DS 下调用方传入的已是 canonical 值，回文自映射不变）→
/// count < MIN_SEED_COVERAGE（严格 <）→ entropy < MIN_SEED_ENTROPY（严格 <）。
/// 原版第 4 参 float min_connectivity 未使用（定义处形参匿名），不移植。
///
/// 泛型 `&C: KmerCatalog`——单线程 KmerCounter 与 PARALLEL 的 SyncKmerCounter
/// （dashmap）共用本直译（原版只有一份 KmerCounter& 代码）。
pub fn is_good_seed_kmer<C: KmerCatalog>(
    counter: &C,
    kmer: KmerId,
    kmer_count: u32,
    p: &IrkeParams,
) -> bool {
    let kmer_length = counter.get_kmer_length();

    if kmer_count == 0 {
        return false; // IRKE.cpp:725-727
    }
    if kmer == revcomp_val(kmer, kmer_length) {
        // 回文 kmer 不作种子（IRKE.cpp:730-739）
        return false;
    }
    if kmer_count < p.min_seed_coverage {
        return false; // IRKE.cpp:742-748
    }
    if compute_entropy(kmer, kmer_length) < p.min_seed_entropy {
        return false; // IRKE.cpp:750-758
    }
    true
}

/// IRKE.cpp:764-818 build_inchworm_contig_from_seed。
/// 返回 (joined_path, total_counts)——原版 total_counts 是 `unsigned int&` 出参，
/// int 求和后转 unsigned 的位型与 i32 一致，此处直接给 i32。
///
/// 流程：前向 inchworm('F') → visitor 清空后回放前向路径、再 erase 种子 →
/// 反向 inchworm('R') → total = fwd + rev + 种子 count（种子不在两条 path 中）→
/// joined = reverse_path 逆序 + seed + forward_path（IRKE.cpp:1150）。
///
/// 泛型 `&C: KmerCatalog`——单线程与 PARALLEL（SyncKmerCounter）共用；
/// 原版末参 `bool PARALLEL_IWORM` 定义处即匿名未用（IRKE.cpp:766），不移植。
pub fn build_inchworm_contig_from_seed<C: KmerCatalog>(
    seed: KmerId,
    counter: &C,
    p: &IrkeParams,
    rng: &mut GlibcRand,
) -> Result<(Vec<KmerId>, i32), CommonError> {
    // IRKE.cpp:769: 种子实时计数
    let kmer_count = counter.get_kmer_count(seed);

    let mut visitor = KmerVisitor::new(counter.get_kmer_length(), counter.is_double_stranded());

    /* Extend to the right */
    let forward = inchworm(counter, Direction::Forward, seed, &mut visitor, p, rng)?;

    // IRKE.cpp:775-792: 清空 visitor 后把前向路径全部标记，使反向延伸不会复用
    visitor.clear();
    for &kmer in &forward.path {
        visitor.add(kmer);
    }

    /* Extend to the left */
    visitor.erase(seed); // IRKE.cpp:795: reset the seed
    let reverse = inchworm(counter, Direction::Reverse, seed, &mut visitor, p, rng)?;

    // IRKE.cpp:801: int + int + unsigned → unsigned（回绕加，位型保留）
    let total_counts = forward
        .count
        .wrapping_add(reverse.count)
        .wrapping_add(kmer_count as i32);

    let joined = join_forward_n_reverse_paths(&reverse.path, seed, &forward.path);
    Ok((joined, total_counts))
}

/// IRKE.cpp:821-918 inchworm —— 从 kmer 出发单方向贪心延伸的外层循环。
/// 每轮：eliminator 清空（PACMAN=false 时恒空）→ 轮数超限检查（原版 throw）→
/// visitor.erase(当前 kmer) → inchworm_step(depth=0, max_recurse=MAX_RECURSION) →
/// `best.count > 0`（严格 >）才延伸：step 返回的 path 是深→浅压栈，**倒序**消费
/// （CRAWL=0 → 全量），逐个 visitor.add；当前 kmer ← entire.path 末尾；
/// entire.count += best.count。
fn inchworm<C: KmerCatalog>(
    counter: &C,
    direction: Direction,
    mut kmer: KmerId,
    visitor: &mut KmerVisitor,
    p: &IrkeParams,
    rng: &mut GlibcRand,
) -> Result<PathNCount, CommonError> {
    let mut entire = PathNCount {
        path: Vec::new(),
        count: 0, // init cumulative path coverage（IRKE.cpp:826）
    };

    let mut inchworm_round: u32 = 0;

    // IRKE.cpp:829: 进入时快照（size 含 count=0 键——惰性删除不缩表）
    let num_total_kmers = counter.size();

    // IRKE.cpp:831: eliminator 在循环外构造、每轮 clear（PACMAN=false 恒空，
    // 保留语义位——Task 6 PARALLEL 版同构）
    let mut eliminator = KmerVisitor::new(counter.get_kmer_length(), counter.is_double_stranded());

    loop {
        inchworm_round += 1;
        eliminator.clear();

        // IRKE.cpp:851-853: 原版 throw(string) → 此处 Err
        if inchworm_round as usize > num_total_kmers {
            return Err(CommonError::Inchworm(
                "Error, inchworm rounds have exceeded the number of possible seed kmers"
                    .to_string(),
            ));
        }

        visitor.erase(kmer); // IRKE.cpp:864: 种子/当前 kmer 必须未访问

        let kmer_pair = (kmer, counter.get_kmer_count(kmer)); // Kmer_Occurence_Pair
        let best = inchworm_step(
            counter,
            direction,
            kmer_pair,
            visitor,
            &mut eliminator,
            inchworm_round,
            0,
            p,
            MAX_RECURSION,
            rng,
        );

        // IRKE.cpp:878: (__DEVEL_zero_kmer_on_use=false) ⇒ 条件就是 best.second > 0
        if best.count > 0 {
            // IRKE.cpp:884-897: 从 last_index(=0, CRAWL=0) 到 first_index 倒序消费
            for &kmer_extend in best.path.iter().rev() {
                entire.path.push(kmer_extend);
                visitor.add(kmer_extend);
            }
            // count>0 ⇒ path 非空（每个 depth>0 的贡献节点各 push 一次，
            // 候选 count!=0；u32>i32::MAX 回绕成负数会被上面的 >0 拦下）
            kmer = *entire.path.last().unwrap(); // IRKE.cpp:898
            entire.count = entire.count.wrapping_add(best.count); // IRKE.cpp:901
        } else {
            // no extension possible（IRKE.cpp:905-907）
            break;
        }
    }

    Ok(entire)
}

/// IRKE.cpp:933-1125 inchworm_step —— 递归贪心体。
///
/// - 已访问/已排除 → 返回空（count==0 的拦截在候选过滤，不在此——原版 953 行
///   的 `!kmer.second` 被注释掉）
/// - `depth < max_recurse` 时：取候选（count>0、!visitor、connectivity 恒过），
///   逐个递归（传 depth+1 与**本层 recurse_cap**），返回后 visitor.erase 回溯
///   （深探完全展开后各层把自己的候选 erase 掉——全量回溯）
/// - paths（size≥1 者）分析：
///   - size>1 → 排序（compare: 发布版仅 count 降序；std::sort 不稳定 → 此处
///     确定性化为稳定排序，平局保持收集序即候选 G,A,T,C 序）
///   - 同分且 `path[0]`（**最远端**，深→浅压栈）不同 → 真 tie：
///     `recurse_cap >= MAX_RECURSION_HARD_STOP(50)` → **rand()%2 二选一**（唯一
///     rand 调用点）；否则 `paths[0].path.len() > best_path_length` → recurse_cap++
///     继续循环（唯一保持 tie 的分支）；否则取 paths[0]
///   - 同分同最远端 / 无 tie → 取 paths[0]；单路径 → 它；零路径 → 空
/// - 尾部：`depth > 0` 时 push 自身 kmer + count += 进入时的 count 快照
///   （**种子 depth=0 不入 path 不计 count**）
// 原版签名（IRKE.cpp:933-937）即 9 参，加 p/rng 后超过 clippy 默认阈值——
// 直译优先，保持参数一一对应
#[allow(clippy::too_many_arguments)]
fn inchworm_step<C: KmerCatalog>(
    counter: &C,
    direction: Direction,
    kmer: (KmerId, u32),
    visitor: &mut KmerVisitor,
    eliminator: &mut KmerVisitor,
    _inchworm_round: u32, // 原版仅 MONITOR 日志用（IRKE.cpp:941-945）
    depth: usize,
    p: &IrkeParams,
    max_recurse: u32,
    rng: &mut GlibcRand,
) -> PathNCount {
    let mut best = PathNCount {
        path: Vec::new(),
        count: 0, // IRKE.cpp:951
    };

    // IRKE.cpp:953-968: visited / eliminated → 空路径
    if visitor.exists(kmer.0) || eliminator.exists(kmer.0) {
        return best;
    }

    visitor.add(kmer.0); // IRKE.cpp:970

    // IRKE.cpp:972-975: PACMAN && depth>0 → eliminator.add —— 本移植 PACMAN=false，
    // eliminator 恒空（inchworm 每轮 clear 后从未 add），保留调用位不保留行为

    if (depth as u32) < max_recurse {
        // IRKE.cpp:977-986: 候选在 while(tie) 之外取一次
        let kmer_candidates = match direction {
            Direction::Forward => counter.get_forward_kmer_candidates(kmer.0),
            Direction::Reverse => counter.get_reverse_kmer_candidates(kmer.0),
        };

        // IRKE.cpp:994-996
        let mut tie = true;
        let mut recurse_cap = max_recurse;
        let mut best_path_length = 0usize;

        while tie {
            // 继续加深递归以打破 tie: recurse_cap 每次 +1 直到 tie 消失或硬停
            let mut paths: Vec<PathNCount> = Vec::new();

            for &kmer_candidate in &kmer_candidates {
                // IRKE.cpp:1007-1010: count!=0（候选已过滤，恒真——保留镜像结构）、
                // 未访问、connectivity 过滤（默认恒过）
                if kmer_candidate.1 != 0
                    && !visitor.exists(kmer_candidate.0)
                    && exceeds_min_connectivity(counter, kmer, kmer_candidate, p.min_connectivity)
                {
                    // IRKE.cpp:1016: 递归传 depth+1 与本层 recurse_cap
                    let sub = inchworm_step(
                        counter,
                        direction,
                        kmer_candidate,
                        visitor,
                        eliminator,
                        _inchworm_round,
                        depth + 1,
                        p,
                        recurse_cap,
                        rng,
                    );
                    if !sub.path.is_empty() {
                        // 只保留实际包含节点的路径（IRKE.cpp:1018-1020 size() >= 1）
                        paths.push(sub);
                    }
                    visitor.erase(kmer_candidate.0); // IRKE.cpp:1022: un-visiting 回溯
                }
            }

            if paths.len() > 1 {
                // IRKE.cpp:1040 sort + compare（IRKE.cpp:921-929 发布版仅 second 降序）。
                // std::sort 不稳定 → 平局为 introsort 实现定义序；确定性化为稳定排序
                paths.sort_by_key(|pn| std::cmp::Reverse(pn.count));

                if paths[0].count == paths[1].count
                    // IRKE.cpp:1049: 检查的是 path[0] —— 路径**最远端**
                    //（深→浅压栈；同端点的同分路径不值得打破）
                    && paths[0].path[0] != paths[1].path[0]
                {
                    // 真 tie: 同分、不同端点
                    if recurse_cap >= MAX_RECURSION_HARD_STOP {
                        // IRKE.cpp:1057-1067: __DEVEL_no_tie_breaking=false，唯一
                        // 存活分支是硬停 → rand()%2（glibc random 复刻，唯一调用点）
                        tie = false;
                        let rand_index = (rng.next() % 2) as usize;
                        best = paths[rand_index].clone(); // 原版 paths[rand_index] 拷贝
                    } else if paths[0].path.len() > best_path_length {
                        // IRKE.cpp:1069-1074: 还在加深（还能延伸）→ 唯一保持
                        // tie=true 的分支
                        recurse_cap += 1;
                        best_path_length = paths[0].path.len();
                    } else {
                        // IRKE.cpp:1075-1079: 深不下去 → 取第一条
                        tie = false;
                        best = paths[0].clone();
                    }
                } else if paths[0].count == paths[1].count && paths[0].path[0] == paths[1].path[0] {
                    // IRKE.cpp:1081-1093: 同分同端点——两条路径汇于同一 kmer，任取
                    tie = false;
                    best = paths[0].clone();
                } else {
                    // IRKE.cpp:1095-1098: 无 tie
                    tie = false;
                    best = paths[0].clone();
                }
            } else if paths.len() == 1 {
                // IRKE.cpp:1100-1103
                tie = false;
                best = paths.pop().unwrap();
            } else {
                // IRKE.cpp:1104-1106: 无延伸可能
                tie = false;
            }
        }
    }

    // IRKE.cpp:1112-1117: 追加当前 kmer——只要不是原始种子（depth=0 不入 path、
    // 不计 count；count 用进入本层时的快照，u32→int 位型保留）
    if depth > 0 {
        best.path.push(kmer.0);
        best.count = best.count.wrapping_add(kmer.1 as i32);
    }

    best
}

/// IRKE.cpp:1128-1150 _join_forward_n_reverse_paths。
/// reverse_path 逆序 + 种子 + forward_path 原序。
fn join_forward_n_reverse_paths(
    reverse_path: &[KmerId],
    seed_kmer_val: KmerId,
    forward_path: &[KmerId],
) -> Vec<KmerId> {
    let mut joined = Vec::with_capacity(reverse_path.len() + 1 + forward_path.len());
    for &kmer in reverse_path.iter().rev() {
        joined.push(kmer);
    }
    joined.push(seed_kmer_val);
    joined.extend_from_slice(forward_path);
    joined
}

/// IRKE.cpp:1159-1176 reconstruct_path_sequence。
/// 首 kmer 全串 + 后续各末碱基；同时快照每 kmer 的计数（原版 cov_counter 出参
/// push_back，此处改为返回值——调用方语义等价）。
///
/// 泛型 `&C: KmerCatalog`——PARALLEL 主循环在 SyncKmerCounter 上实时快照
/// （IRKE.cpp:566 同样在（并发读写的）kcounter 上取计数）。
pub fn reconstruct_path_sequence<C: KmerCatalog>(
    counter: &C,
    path: &[KmerId],
) -> (String, Vec<u32>) {
    if path.is_empty() {
        return (String::new(), Vec::new());
    }
    let kmer_length = counter.get_kmer_length();
    let mut cov_counter: Vec<u32> = Vec::with_capacity(path.len());
    // get_kmer_string 输出恒为 gatc ASCII（INT_TO_BASE），from_utf8 不可能失败
    let first = counter.get_kmer_string(path[0]);
    let mut seq = String::with_capacity(kmer_length + path.len() - 1);
    seq.push_str(std::str::from_utf8(&first).unwrap());
    cov_counter.push(counter.get_kmer_count(path[0]));

    for &kmer in &path[1..] {
        let kmer_str = counter.get_kmer_string(kmer);
        seq.push(kmer_str[kmer_length - 1] as char); // substr(len-1, 1)
        cov_counter.push(counter.get_kmer_count(kmer));
    }
    (seq, cov_counter)
}

/// IRKE.cpp:1179-1213 exceeds_min_connectivity。
/// 第一行短路 `min_connectivity < 1e5 → true`（consider test off）——默认 0 恒过，
/// 后续死代码仅为逐行镜像保留。counter 形参原版即未使用（匿名 `KmerCounter&`）。
fn exceeds_min_connectivity<C: KmerCatalog>(
    _counter: &C,
    kmer_a: (KmerId, u32),
    kmer_b: (KmerId, u32),
    min_connectivity: f32,
) -> bool {
    if min_connectivity < 1e5 {
        return true; // consider test off
    }

    let kmer_a_count = kmer_a.1;
    if kmer_a_count == 0 {
        return false;
    }
    let kmer_b_count = kmer_b.1;
    if kmer_b_count == 0 {
        return false;
    }

    let (min_val, max_val) = if kmer_a_count < kmer_b_count {
        (kmer_a_count, kmer_b_count)
    } else {
        (kmer_b_count, kmer_a_count)
    };

    let connectivity_ratio = min_val as f32 / max_val as f32;
    connectivity_ratio >= min_connectivity
}

// ===========================================================================
// 组装主循环（IRKE.cpp:426-716 compute_sequence_assemblies，单线程路径）
// ===========================================================================

/// IRKE_COMMON::MONITOR（IRKE_common.hpp 的进程级全局）的移植载体——显式传递。
/// level 语义同原版: 0 = 关，>=1 常规，>=2 详细（SEED kmer 行等）。
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Monitor {
    pub level: u32,
}

impl Monitor {
    pub fn new(level: u32) -> Self {
        Monitor { level }
    }

    /// 原版 `if (IRKE_COMMON::MONITOR)`（IRKE.cpp:590）——非零即真。
    pub fn enabled(&self) -> bool {
        self.level >= 1
    }
}

/// IRKE.cpp:417-424 compute_sequence_assemblies 的三个标量参数。
/// min_connectivity 是**组装级**参数——原版 IRKE 对象不持有它，由本函数形参
/// 逐层下传（IRKE.cpp:541/552），故此处覆盖 IrkeParams 的同名字段。
#[derive(Debug, Clone, Copy)]
pub struct AssemblyParams {
    pub min_connectivity: f32,
    /// -L，默认 25。注意 IRKE_run.cpp:90 `MIN_ASSEMBLY_LENGTH = kmer_length` 在
    /// -K（IRKE_run.cpp:190）解析**之前**定值——默认恒 25，不随 -K 变（CLI 层
    /// 负责该默认，不在此处绑定 kmer_length）。
    pub min_assembly_length: usize,
    /// --min_assembly_coverage，默认 2（IRKE_run.cpp:91）
    pub min_assembly_coverage: u32,
}

impl Default for AssemblyParams {
    fn default() -> Self {
        AssemblyParams {
            min_connectivity: 0.0,
            min_assembly_length: 25,
            min_assembly_coverage: 2,
        }
    }
}

/// IRKE.cpp:426-716 compute_sequence_assemblies —— 单线程路径（无 PARALLEL_IWORM/
/// TWO_PHASE——那是 Task 6 的分支，IRKE.cpp:549 只在 PARALLEL 下启用）。
///
/// - 种子列表在函数内构建（镜像 IRKE_run.cpp:521-523 的 populate_sorted_kmers_
///   list + KC:701-752 的 stderr 文案）
/// - 逐种子: 实时查 count（排序快照已过期，IRKE.cpp:546 注释）→ is_good_seed_
///   kmer → build_inchworm_contig_from_seed → 记录条件
///   `seq.len() >= MIN_ASSEMBLY_LENGTH && avg_cov >= MIN_ASSEMBLY_COVERAGE`
///   （IRKE.cpp:569-573）→ **无论是否记录**，joined_path 全部 clear_kmer
///   （IRKE.cpp:564-574 的 `__DEVEL_zero_kmer_on_use=false` 分支）
/// - 去重输出: key = generateHash 低 32 位（原版 `unsigned int contig_hash =
///   generateHash(seq)` 截断后回宽 u64 作 map key，等价于低 32 位判等），
///   首见才写，header 逐字段镜像 IRKE.cpp:662-670，60 列折行
/// - rng = GlibcRand::new(1): 原版从不 srandom → srand(1) 序列，**整个组装共用
///   一个实例**（进程级全局 state 的移植）
///
/// 记录类型 (total_counts: i32, avg_cov, seed 快照 count, sequence)——total_counts
/// 位型同原版 unsigned int（i32 回绕一致），输出时 `as u32`。原版经
/// tmp.iworm.fa 临时文件中转（IRKE.cpp:505-517），单线程下内存顺序 = tmp 文件序，
/// 故不再落盘（--keep_tmp_files 由 CLI 收下为 no-op）。
///
/// 写入 w（stdout 由调用方传入），返回产出 contig 数（INCHWORM_ASSEMBLY_COUNTER）。
pub fn compute_sequence_assemblies<W: Write>(
    counter: &mut KmerCounter,
    params: &IrkeParams,
    aparams: &AssemblyParams,
    monitor: &Monitor,
    sort_seeds: bool,
    w: &mut W,
) -> Result<usize, CommonError> {
    // ---- 种子列表（IRKE_run.cpp:521-523 + KC:701-752 的文案） ----
    eprintln!("-populating the kmer seed candidate list."); // IRKE_run.cpp:521
    eprintln!("Kcounter hash size: {}", counter.size()); // KC:716
    let mut kmers = collect_nonzero_seeds(counter); // KC:719-728
    eprintln!(
        "Processed {} non-zero abundance kmers in kcounter.",
        kmers.len()
    ); // KC:729
    if sort_seeds {
        eprint!("Sorting {} kmers ...", counter.size()); // KC:739（无换行）
        let sort_start = Instant::now();
        sort_seeds_desc(&mut kmers); // KC:741 + KC:823-831
        eprintln!(
            "Done sorting {} kmers, taking {} seconds.",
            counter.size(),
            sort_start.elapsed().as_secs()
        ); // KC:746-752
    } else {
        // KC:732-737（__DEVEL_no_kmer_sort，PARALLEL 分支——Task 6 走到）
        eprintln!("-Not sorting list of kmers, given parallel mode in effect.");
    }
    compute_sequence_assemblies_from_seeds(counter, params, aparams, monitor, &kmers, w)
}

/// compute_sequence_assemblies 的主循环体（种子列表已由调用方构建）。
/// 拆出公开: (a) T6 PARALLEL 用未排序列表调同一循环; (b) 对拍取证——用原版
/// --monitor 2 抓到的种子序重放，隔离「平局排序序」与「贪心核心」两类分歧
/// （seeds 的 count 字段不参与主循环——实时重查，IRKE.cpp:546）。
pub fn compute_sequence_assemblies_from_seeds<W: Write>(
    counter: &mut KmerCounter,
    params: &IrkeParams,
    aparams: &AssemblyParams,
    monitor: &Monitor,
    kmers: &[(KmerId, u32)],
    w: &mut W,
) -> Result<usize, CommonError> {
    // ---- IRKE.cpp:441-449: 组装中哈希只减不增的护栏基线 ----
    let init_size = counter.size();
    eprintln!(
        "Total kcounter hash size: {init_size} vs. sorted list size: {}",
        kmers.len()
    ); // IRKE.cpp:449

    // IRKE.cpp:494-501: 单线程 omp_set_num_threads(1) → omp_get_max_threads()==1
    eprintln!("num threads set to: 1");

    let kmer_length = counter.get_kmer_length();
    // 组装级 min_connectivity 下传（覆盖 IrkeParams 同名字段——镜像原版数据流）
    let effective = IrkeParams {
        min_connectivity: aparams.min_connectivity,
        ..*params
    };
    // 原版 rand() 是进程级 glibc 状态（从不 srandom → srand(1)）——整个组装共用
    let mut rng = GlibcRand::new(1);

    // 原版写 tmp.iworm.fa.pid_X.thread_0 再回读去重（IRKE.cpp:505-517 + 624-686）;
    // 单线程下生成序 = tmp 文件序，内存等价（无 tmp 文件）
    let mut records: Vec<(i32, u32, u32, String)> = Vec::new();

    for &(kmer, _cached_count) in kmers {
        // IRKE.cpp:530-539: size 增长即原版 throw 的护栏——本实现只置 0 不插键，
        // 天然满足；保留 debug_assert 防将来回归
        debug_assert!(
            counter.size() <= init_size,
            "Error, Kcounter size has grown from {init_size} to {}",
            counter.size()
        );

        // IRKE.cpp:546: 实时计数（NOT the sorted snapshot——大概率已被清零）
        let kmer_count = counter.get_kmer_count(kmer);
        if !is_good_seed_kmer(counter, kmer, kmer_count, &effective) {
            continue; // IRKE.cpp:548-551
        }

        if monitor.level >= 2 {
            // IRKE.cpp:556-558（"Seed for thread" 行是 omp 风味，省略）
            eprintln!(
                "SEED kmer: {}, count: {kmer_count}",
                String::from_utf8_lossy(&counter.get_kmer_string(kmer))
            );
        }

        // IRKE.cpp:553-554（无 PARALLEL 分支——TWO_PHASE 只在 PARALLEL 下）
        let (joined_path, total_counts) =
            build_inchworm_contig_from_seed(kmer, counter, &effective, &mut rng)?;

        // IRKE.cpp:560-561
        let (sequence, _assembly_base_coverage) = reconstruct_path_sequence(counter, &joined_path);
        // IRKE.cpp:563: (float)total/(len-K+1) + 0.5 → unsigned。
        // 除法两操作数均转 f32（C++ 整型遇 float 提升），+0.5 在 double 完成
        // （0.5 是 double 字面量）；分母 len-K+1 = kmer 数 = joined_path.len()。
        // total_counts 先经 u32 位型再转 f32——镜像原版 unsigned int 形参
        let avg_cov =
            (((total_counts as u32) as f32 / joined_path.len() as f32) as f64 + 0.5) as u32;

        // IRKE.cpp:569-573
        if sequence.len() >= aparams.min_assembly_length && avg_cov >= aparams.min_assembly_coverage
        {
            records.push((total_counts, avg_cov, kmer_count, sequence));
        }

        // IRKE.cpp:564-574: 无论是否记录都清零路径 kmer（else 分支 = 发布版行为）
        for path_kmer in &joined_path {
            counter.clear_kmer(*path_kmer);
        }
    }

    if monitor.enabled() {
        eprintln!(); // IRKE.cpp:590-592
    }

    // ---- IRKE.cpp:624-686: 回读 tmp 文件、去重、输出 ----
    let mut seen_contig_already = FxHashSet::default();
    let mut inchworm_assembly_counter = 0usize;
    for (total_counts, avg_cov, kmer_count, sequence) in records {
        // unsigned int 接收 generateHash 的低 32 位（原版隐式截断）
        let contig_hash = generate_hash(sequence.as_bytes()) as u32;
        if seen_contig_already.insert(contig_hash) {
            inchworm_assembly_counter += 1; // INCHWORM_ASSEMBLY_COUNTER++
                                            // IRKE.cpp:662-670: ">a{i};{avg} total_counts: {tc} Seed: {seed} K: {K} length: {len}"
            writeln!(
                w,
                ">a{inchworm_assembly_counter};{avg_cov} total_counts: {} Seed: {kmer_count} K: {kmer_length} length: {}",
                total_counts as u32, // tmp 文件 unsigned int 位型（i32 回绕一致）
                sequence.len()
            )?;
            writeln!(w, "{}", add_fasta_seq_line_breaks(sequence.as_bytes(), 60))?;
            // IRKE.cpp:671
        }
    }

    Ok(inchworm_assembly_counter)
}

// ===========================================================================
// PARALLEL_IWORM 组装（IRKE.cpp:460-686 的并行分支，Task 6）
// ===========================================================================

/// IRKE.cpp:504 `#pragma omp parallel for ... schedule(dynamic, 1000)` 的
/// chunk 大小——每 1000 个种子为一个工作单元，单元内顺序处理、单元间并行。
pub const PARALLEL_CHUNK_SIZE: usize = 1000;

/// 主循环产出记录（total_counts, avg_cov, 主循环快照的种子 count, sequence）。
/// 原版写线程临时文件 `tmp.iworm.fa.pid_X.thread_N`（IRKE.cpp:474-489，4 字段
/// 逐行），回读去重输出（IRKE.cpp:624-686）；此处内存中转，字段一致。
type IwormRecord = (i32, u32, u32, String);

/// IRKE.cpp:426-686 compute_sequence_assemblies 的 PARALLEL_IWORM 分支。
///
/// 与单线程路径（[`compute_sequence_assemblies`]）的差异全部来自
/// IRKE.cpp:460-620 的并行分支，逐点镜像：
/// - **种子列表不排序**（KC:732-737 `__DEVEL_no_kmer_sort`，收集序即容器迭代序
///   ——本移植为 FxHashMap 迭代序；原版为 hash_map 迭代序，同为任意序）
/// - **chunk 并行**: `par_chunks(1000)` ≙ `schedule(dynamic, 1000)`——chunk 为
///   rayon 工作单元，动态分发到池线程，chunk 内顺序；结果按**全局种子序**
///   收集（原版按 thread 0..N 临时文件序输出去重——两序都源自 nondeterministic
///   的并行过程，去重 key（generateHash 低 32 位）保证内容一致，序为已声明的
///   实现差异）
/// - **目录弱一致**: 组装期读/清零走 [`SyncKmerCounter`]（dashmap 单键原子，
///   跨键无一致性——镜像原版对同一 hash_map 的无锁并发读写语义）
/// - **TWO_PHASE**（默认开，IRKE.cpp:59；`--SINGLE_PHASE` 关闭）: draft path →
///   extract_best_seed（实时查）→ new_seed 已被其它线程清零（zapped）则放弃
///   本种子（不记录也**不清零**，IRKE.cpp:552-557 的 continue）→ 否则从
///   new_seed 重建
/// - **rand**: 原版 rand() 是进程级 glibc 状态，多线程调用交错为 nondeterministic。
///   本移植**每 chunk 独立 `GlibcRand::new(1)`**（= srand(1)，与单线程版同种
///   子）——chunk 边界与 rayon 工作单元对齐 → 同输入同 chunk 划分 → 每 chunk
///   的 tie 打破序列确定；在"rand 仅用于 50 层硬停后的二选一平局打破"这一用途
///   上与原版语义等价（原版交错的随机源也只是任意打破平局）
/// - **线程数**: `omp_set_num_threads(NUM_THREADS)`（IRKE.cpp:464-465）→ 自建
///   rayon 线程池；`num_threads=None` 用 rayon 默认（RAYON_NUM_THREADS 或核数，
///   对齐原版未给 --num_threads 时取 omp_get_max_threads()）
///
/// 不镜像的部分: 原版 per-thread 临时文件的 "Done opening file." stderr 行与
/// 文件落盘/删除（IRKE.cpp:479-489, 663-685）——本实现内存中转，不落盘
/// （`--keep_tmp_files` 对 PARALLEL 同为 no-op）。
///
/// `counter` 按值接收（装载+剪枝已在单线程完成，整表转入 SyncKmerCounter）。
pub fn compute_sequence_assemblies_parallel<W: Write>(
    counter: KmerCounter,
    params: &IrkeParams,
    aparams: &AssemblyParams,
    monitor: &Monitor,
    two_phase: bool,
    num_threads: Option<usize>,
    w: &mut W,
) -> Result<usize, CommonError> {
    // ---- 种子列表（不排序: KC:732-737，收集序 = FxHashMap 迭代序） ----
    eprintln!("-populating the kmer seed candidate list."); // IRKE_run.cpp:521
    eprintln!("Kcounter hash size: {}", counter.size()); // KC:716
    let kmers = collect_nonzero_seeds(&counter); // KC:719-728
    eprintln!(
        "Processed {} non-zero abundance kmers in kcounter.",
        kmers.len()
    ); // KC:729
    eprintln!("-Not sorting list of kmers, given parallel mode in effect."); // KC:732-737

    // ---- IRKE.cpp:441-449: 组装中哈希只减不增的护栏基线 ----
    let init_size = counter.size();
    eprintln!(
        "Total kcounter hash size: {init_size} vs. sorted list size: {}",
        kmers.len()
    ); // IRKE.cpp:449

    // IRKE.cpp:464-465: omp_set_num_threads(NUM_THREADS) → 自建 rayon 池
    let mut builder = rayon::ThreadPoolBuilder::new();
    if let Some(n) = num_threads {
        builder = builder.num_threads(n);
    }
    let pool = builder
        .build()
        .map_err(|e| CommonError::Inchworm(format!("rayon thread pool build failed: {e}")))?;
    eprintln!("num threads set to: {}", pool.current_num_threads()); // IRKE.cpp:471-472（omp_get_max_threads）

    // 目录转入并发结构（此后组装期只读 + 原子置 0，不增不删键）
    let sync = SyncKmerCounter::from_counter(counter);
    let kmer_length = sync.get_kmer_length();
    // 组装级 min_connectivity 下传（覆盖 IrkeParams 同名字段——镜像原版数据流）
    let effective = IrkeParams {
        min_connectivity: aparams.min_connectivity,
        ..*params
    };

    // ---- IRKE.cpp:504-620: chunk 并行主循环 ----
    // par_chunks 为带索引的并行迭代器 → collect 保序（全局种子序）；任一 chunk
    // 的 Err（inchworm 轮数超限 throw 的移植）经 Result 收集提前冒出。
    let per_chunk: Result<Vec<Vec<Option<IwormRecord>>>, CommonError> = pool.install(|| {
        kmers
            .par_chunks(PARALLEL_CHUNK_SIZE)
            .map(|chunk| {
                // 每 chunk 独立 srand(1)（见函数文档 rand 决策）
                let mut rng = GlibcRand::new(1);
                let mut buf: Vec<Option<IwormRecord>> = Vec::with_capacity(chunk.len());
                for &(kmer, _cached_count) in chunk {
                    buf.push(parallel_assemble_seed(
                        kmer, &sync, &effective, aparams, two_phase, monitor, init_size, &mut rng,
                    )?);
                }
                Ok(buf)
            })
            .collect()
    });
    let records: Vec<IwormRecord> = per_chunk?
        .into_iter()
        .flat_map(|chunk| chunk.into_iter().flatten())
        .collect();

    if monitor.enabled() {
        eprintln!(); // IRKE.cpp:622-624
    }

    // ---- IRKE.cpp:624-686: 去重输出（全局种子序;key = generateHash 低 32 位） ----
    let mut seen_contig_already = FxHashSet::default();
    let mut inchworm_assembly_counter = 0usize;
    for (total_counts, avg_cov, kmer_count, sequence) in records {
        // unsigned int 接收 generateHash 的低 32 位（原版隐式截断）
        let contig_hash = generate_hash(sequence.as_bytes()) as u32;
        if seen_contig_already.insert(contig_hash) {
            inchworm_assembly_counter += 1;
            // IRKE.cpp:662-670: ">a{i};{avg} total_counts: {tc} Seed: {seed} K: {K} length: {len}"
            writeln!(
                w,
                ">a{inchworm_assembly_counter};{avg_cov} total_counts: {} Seed: {kmer_count} K: {kmer_length} length: {}",
                total_counts as u32, // tmp 文件 unsigned int 位型（i32 回绕一致）
                sequence.len()
            )?;
            writeln!(w, "{}", add_fasta_seq_line_breaks(sequence.as_bytes(), 60))?;
        }
    }

    Ok(inchworm_assembly_counter)
}

/// PARALLEL 主循环体——单种子处理（IRKE.cpp:504-620 循环内的全部逻辑，含
/// TWO_PHASE 二段确认）。返回 None = 不产出记录（种子不合格 / zapped）。
// 循环体依赖的 8 个量各自独立（目录/参数/开关/护栏/随机源），打包成上下文
// 结构反而遮蔽与原版循环变量的对应——沿用 inchworm_step 的直译豁免
#[allow(clippy::too_many_arguments)]
fn parallel_assemble_seed(
    kmer: KmerId,
    sync: &SyncKmerCounter,
    p: &IrkeParams,
    aparams: &AssemblyParams,
    two_phase: bool,
    monitor: &Monitor,
    init_size: usize,
    rng: &mut GlibcRand,
) -> Result<Option<IwormRecord>, CommonError> {
    // IRKE.cpp:513-523: size 增长即原版 throw 的护栏——本实现只置 0 不插键，
    // 天然满足；保留 debug_assert 防将来回归
    debug_assert!(
        sync.size() <= init_size,
        "Error, Kcounter size has grown from {init_size} to {}",
        sync.size()
    );

    // IRKE.cpp:530-546: 实时计数（NOT the collected snapshot——大概率已被清零）
    let kmer_count = sync.get_kmer_count(kmer);
    if !is_good_seed_kmer(sync, kmer, kmer_count, p) {
        return Ok(None); // IRKE.cpp:548-551
    }

    if monitor.level >= 2 {
        eprintln!(
            "SEED kmer: {}, count: {kmer_count}",
            String::from_utf8_lossy(&sync.get_kmer_string(kmer))
        ); // IRKE.cpp:556-558
        eprintln!(
            "Seed for thread: {} is {} with count: {kmer_count}",
            rayon::current_thread_index().unwrap_or(0),
            String::from_utf8_lossy(&sync.get_kmer_string(kmer))
        ); // IRKE.cpp:559-561（omp critical 的镜像——rayon 线程索引）
    }

    // IRKE.cpp:562-564: draft contig
    let (mut joined_path, mut total_counts) = build_inchworm_contig_from_seed(kmer, sync, p, rng)?;

    // IRKE.cpp:566-577: TWO_PHASE —— draft path 中选实时计数最高的合格种子重建
    if two_phase {
        let new_seed = extract_best_seed(&joined_path, sync, p);
        if sync.get_kmer_count(new_seed) == 0 {
            // must have been zapped by another thread——放弃本种子:
            // 不记录也**不清零** draft path（原版 continue 跳过其后的全部逻辑）
            return Ok(None);
        }
        let rebuilt = build_inchworm_contig_from_seed(new_seed, sync, p, rng)?;
        joined_path = rebuilt.0;
        total_counts = rebuilt.1;
    }

    // IRKE.cpp:579-581: 序列重建
    let (sequence, _assembly_base_coverage) = reconstruct_path_sequence(sync, &joined_path);
    // IRKE.cpp:583: (float)total/(len-K+1) + 0.5 → unsigned（同单线程版注释）
    let avg_cov = (((total_counts as u32) as f32 / joined_path.len() as f32) as f64 + 0.5) as u32;

    // IRKE.cpp:585-595: 记录条件（原版先写线程临时文件，此处收集后统一去重）
    let record = if sequence.len() >= aparams.min_assembly_length
        && avg_cov >= aparams.min_assembly_coverage
    {
        Some((total_counts, avg_cov, kmer_count, sequence))
    } else {
        None
    };

    // IRKE.cpp:597-618: 无论是否记录都清零路径 kmer（dashmap 原子置 0——弱一致，
    // 并发读者可能仍见非零，镜像原版无锁竞态）
    for path_kmer in &joined_path {
        sync.clear_kmer(*path_kmer);
    }

    Ok(record)
}

/// IRKE.cpp:1348-1371 extract_best_seed。
/// `count > best_kmer_count && is_good_seed_kmer(...)`（&& 短路，左边先行）——
/// 平局取路径中更早出现者；找不到返回 0。
///
/// PARALLEL TWO_PHASE 在 draft path 上实时查（可能正被其它线程清零——弱一致，
/// 调用方随后自查 get_kmer_count(new_seed)==0 的 zapped 分支，IRKE.cpp:552-556）。
pub fn extract_best_seed<C: KmerCatalog>(path: &[KmerId], counter: &C, p: &IrkeParams) -> KmerId {
    let mut best_kmer_count: u32 = 0;
    let mut best_seed: KmerId = 0;

    for &kmer in path {
        let count = counter.get_kmer_count(kmer);
        if count > best_kmer_count && is_good_seed_kmer(counter, kmer, count, p) {
            best_kmer_count = count;
            best_seed = kmer;
        }
    }

    best_seed
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    /// K=2 编码表（首碱基 = v>>2，次碱基 = v&3；G=0,A=1,T=2,C=3）:
    /// GA=1 AA=5 AC=7 TG=8 TA=9 TT=10 CG=12 CA=13 CT=14 CC=15
    fn enc(s: &[u8]) -> KmerId {
        kmer_to_intval(s).unwrap()
    }

    // ---- c_atoi（黄金值来自 glibc atoi 的 C harness 实测） ----------------

    #[test]
    fn c_atoi_matches_glibc_harness() {
        // (输入, (int)atoi 赋值 unsigned int 后的位型)
        let cases: &[(&[u8], u32)] = &[
            (b"5x", 5),                              // 首非数字截断
            (b" 7", 7),                              // 跳前导空白
            (b"x", 0),                               // 非数字开头 = 0
            (b"-3", 4294967293),                     // -3 → u32 回绕
            (b"", 0),                                // 空串
            (b"4294967296", 0),                      // long 容纳, (int) 截断回绕
            (b"9999999999999999999999", 4294967295), // long 饱和 LONG_MAX → (int) = -1
            (b"-9999999999999999999999", 0),         // long 饱和 LONG_MIN → 低 32 位 0
            (b"+12", 12),                            // 显式正号
            (b"  -42abc", 4294967254),               // 空白 + 负号 + 截断
            (b"3.9", 3),                             // '.' 截断
            (b"2147483648", 2147483648),             // i32 回绕为负, 位型保留
            (b"-2147483649", 2147483647),            // i64 值 (int) 截断
            (b"0", 0),
            (b"007", 7),     // 前导零
            (b"\t9", 9),     // \t 是 C 空白
            (b"  +0x10", 0), // '0' 后 'x' 截断（无 16 进制）
        ];
        for &(input, expected) in cases {
            assert_eq!(c_atoi(input), expected, "c_atoi({:?})", input);
        }
    }

    // ---- populate_from_kmers ----------------------------------------------

    const K25: &[u8] = b"ACGTACGTACGTACGTACGTACGTA"; // 6*ACGT + A = 25

    #[test]
    fn populate_from_kmers_accumulates_lower_equivalent_and_counts_skipped() {
        // 两条同 kmer 记录（大写 + 小写）累加 5；短(4)/长(26)记录跳过但计入 parsed
        let data = format!(
            ">3\n{}\n>2\n{}\n>9\nACGT\n>4\n{}C\n",
            String::from_utf8_lossy(K25),
            String::from_utf8_lossy(K25).to_lowercase(),
            String::from_utf8_lossy(K25)
        );
        let (counter, parsed) = populate_from_kmers(data.as_bytes(), 25, false).unwrap();
        assert_eq!(parsed, 4);
        assert_eq!(counter.size(), 1);
        assert_eq!(counter.get_kmer_count(enc(K25)), 5);
    }

    #[test]
    fn populate_from_kmers_header_atoi_variants() {
        // "5x"→5、" 7"→7、"x"→0（add 0 仍建键）、"-3"→u32 回绕
        let data = ">5x\nAC\n> 7\nAC\n>x\nAA\n>-3\nAT\n";
        let (counter, parsed) = populate_from_kmers(data.as_bytes(), 2, false).unwrap();
        assert_eq!(parsed, 4);
        assert_eq!(counter.size(), 3);
        assert_eq!(counter.get_kmer_count(enc(b"AC")), 12); // 5 + 7
        assert_eq!(counter.get_kmer_count(enc(b"AA")), 0); // atoi("x")=0，键仍存在
        assert_eq!(counter.get_kmer_count(enc(b"AT")), 4294967293);
    }

    #[test]
    fn populate_from_kmers_ds_mode_folds_canonical() {
        // revcomp("ACGG") = "CCGT" → DS 下同键累加
        let data = ">3\nACGG\n>2\nCCGT\n";
        let (counter, parsed) = populate_from_kmers(data.as_bytes(), 4, true).unwrap();
        assert_eq!(parsed, 2);
        assert_eq!(counter.size(), 1);
        assert_eq!(counter.get_kmer_count(enc(b"ACGG")), 5);
        assert_eq!(counter.get_kmer_count(enc(b"CCGT")), 5);
    }

    #[test]
    fn populate_from_kmers_terminates_on_empty_sequence_record() {
        // IRKE.cpp:107-108: 空序列记录 break——后续 ">7 TT" 不再解析，且该记录不计 parsed
        let data = ">5\nAC\n>x\n>7\nTT\n";
        let (counter, parsed) = populate_from_kmers(data.as_bytes(), 2, false).unwrap();
        assert_eq!(parsed, 1);
        assert_eq!(counter.size(), 1);
        assert_eq!(counter.get_kmer_count(enc(b"AC")), 5);
        assert_eq!(counter.get_kmer_count(enc(b"TT")), 0);
    }

    #[test]
    fn populate_from_kmers_rejects_non_gatc() {
        // 原版 kmer_to_intval throw（IRKE.cpp:128）→ 这里 Err
        let data = ">1\nACGN\n";
        let err = populate_from_kmers(data.as_bytes(), 4, false).unwrap_err();
        assert!(matches!(err, CommonError::NonGatcChar { .. }));
    }

    // ---- populate_from_reads ----------------------------------------------

    #[test]
    fn populate_from_reads_sliding_windows_cov1() {
        // "ACGTA" K=4 → 窗口 ACGT、CGTA 各 +1
        let (counter, parsed) = populate_from_reads(b">r1\nACGTA\n", 4, false).unwrap();
        assert_eq!(parsed, 1);
        assert_eq!(counter.size(), 2);
        assert_eq!(counter.get_kmer_count(enc(b"ACGT")), 1);
        assert_eq!(counter.get_kmer_count(enc(b"CGTA")), 1);
    }

    #[test]
    fn populate_from_reads_skips_windows_with_n_but_keeps_sliding() {
        // "ACGTNCGTA" K=4: 含 N 的窗口跳过，N 后窗口恢复即计 CGTA（断整条记录会漏）
        let (counter, _) = populate_from_reads(b">r1\nACGTNCGTA\n", 4, false).unwrap();
        assert_eq!(counter.size(), 2);
        assert_eq!(counter.get_kmer_count(enc(b"ACGT")), 1);
        assert_eq!(counter.get_kmer_count(enc(b"CGTA")), 1);
    }

    #[test]
    fn populate_from_reads_counts_short_and_empty_records_no_terminate() {
        // reads 模式 hasNext 循环（IRKE.cpp:183）无空记录终止: 空记录/短记录照常计数
        let data = ">r1\nAC\n>x\n>7\nTT\n"; // 记录2 = "x"/""（空序列）
        let (counter, parsed) = populate_from_reads(data.as_bytes(), 2, false).unwrap();
        assert_eq!(parsed, 3);
        assert_eq!(counter.get_kmer_count(enc(b"AC")), 1);
        assert_eq!(counter.get_kmer_count(enc(b"TT")), 1);
        assert_eq!(counter.size(), 2);
        // 长度 < K: 跳过（严格 <，与 kmers 模式的 != 不同——恰好 == K 的记录不跳）
        let (counter2, parsed2) = populate_from_reads(b">r0\nACG\n>r1\nAC\n", 4, false).unwrap();
        assert_eq!(parsed2, 2);
        assert_eq!(counter2.size(), 0);
    }

    #[test]
    fn populate_from_reads_accumulates_and_ds_folds() {
        // 两条同读 → 每窗口 2; DS: "ACGG" 与 revcomp "CCGT" 折同键
        let (c1, _) = populate_from_reads(b">r1\nACGTA\n>r2\nACGTA\n", 4, false).unwrap();
        assert_eq!(c1.get_kmer_count(enc(b"ACGT")), 2);
        assert_eq!(c1.get_kmer_count(enc(b"CGTA")), 2);
        let (c2, _) = populate_from_reads(b">r1\nACGG\n>r2\nCCGT\n", 4, true).unwrap();
        assert_eq!(c2.size(), 1);
        assert_eq!(c2.get_kmer_count(enc(b"ACGG")), 2);
    }

    // ---- write_kmer_count_report ------------------------------------------

    #[test]
    fn write_kmer_count_report_single_line() {
        let path = std::env::temp_dir().join("inchworm_irke_report_t1");
        write_kmer_count_report(&path, 2).unwrap();
        assert_eq!(fs::read_to_string(&path).unwrap(), "2\n");
        let _ = fs::remove_file(&path);
    }

    // ---- prune_some_kmers ---------------------------------------------------

    #[test]
    fn prune_min_count_strict_less() {
        // count == min_count 不剪（严格 <）；count < min 即时置 0
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(7, 5);
        kc.add_kmer(12, 1);
        kc.add_kmer(13, 2); // 恰好等于 min_count → 保留
        let pruned = prune_some_kmers(&mut kc, 2, 0.0, false, 0.0);
        assert_eq!(pruned, 1);
        assert_eq!(kc.get_kmer_count(12), 0);
        assert_eq!(kc.get_kmer_count(7), 5);
        assert_eq!(kc.get_kmer_count(13), 2);
        assert_eq!(kc.size(), 3); // 惰性删除: 键仍占位
    }

    #[test]
    fn prune_low_entropy() {
        // K=4: AAAA 熵 0 < 0.5 → 剪; ACGT 熵 2.0 → 留
        let mut kc = KmerCounter::new(4, false);
        kc.add_kmer(enc(b"AAAA"), 3);
        kc.add_kmer(enc(b"ACGT"), 2);
        let pruned = prune_some_kmers(&mut kc, 1, 0.5, false, 0.0);
        assert_eq!(pruned, 1);
        assert_eq!(kc.get_kmer_count(enc(b"AAAA")), 0);
        assert_eq!(kc.get_kmer_count(enc(b"ACGT")), 2);
    }

    /// 手推场景（K=2, SS, r=0.03, min_count=1, min_entropy=0）。
    /// 前向候选前缀 = (s & 3) << 2（种子末碱基移到首位）: 末碱基为 C 的种子
    /// （AC/TC）前向候选集同为 {CG, CA, CT, CC}。
    /// 键 AC(7):50, TC(11):50, CG(12):100, CA(13):1:
    /// - AC 前向 [CG:100, CA:1]: CA 比值 1/100=0.01 < 0.03 且 1/50=0.02 < 0.03 → 入列
    /// - TC 前向同集 [CG:100, CA:1]: 同上 → **再次**入列（重复入列，幂等置 0）
    /// - AC 反向 [CA:1]、TC 反向空、CG 前向空、CA 前向 [AC:50] 均单元素 → 不处理
    /// - CG 反向 [AC:50, TC:50]: dominant=50, 50/50=1.0 ≮ 0.03 → 不入列（dominant 自身段）
    /// - CA 反向 [AC:50, TC:50]: 同上不入列; CA 前向比值 1/1=1.0 ≮ r → 自保护
    ///
    /// → count_pruned = 2（CA 被两个不同种子各计一次），最终仅 CA 被清。
    #[test]
    fn prune_error_kmers_double_ratio_and_duplicate_push() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AC"), 50);
        kc.add_kmer(enc(b"TC"), 50);
        kc.add_kmer(enc(b"CG"), 100);
        kc.add_kmer(enc(b"CA"), 1);
        let pruned = prune_some_kmers(&mut kc, 1, 0.0, true, 0.03);
        assert_eq!(pruned, 2);
        assert_eq!(kc.get_kmer_count(enc(b"CA")), 0);
        assert_eq!(kc.get_kmer_count(enc(b"AC")), 50);
        assert_eq!(kc.get_kmer_count(enc(b"TC")), 50);
        assert_eq!(kc.get_kmer_count(enc(b"CG")), 100);
        assert_eq!(kc.size(), 4); // 惰性删除: CA 键仍占位
    }

    /// 恰好等于 r 时不剪（严格 <）。r=0.5: AC 前向 [CG:2, CA:1] → 1/2 = 0.5 == r
    ///（f32 精确可表示）→ 不入列; 其余扫描均单元素。
    #[test]
    fn prune_error_kmers_exact_ratio_not_pruned() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AC"), 2);
        kc.add_kmer(enc(b"CG"), 2);
        kc.add_kmer(enc(b"CA"), 1);
        let pruned = prune_some_kmers(&mut kc, 1, 0.0, true, 0.5);
        assert_eq!(pruned, 0);
        assert_eq!(kc.get_kmer_count(enc(b"CA")), 1); // 未被剪
                                                      // 对照: dominant=3 时 1/3 ≈ 0.333 < 0.5 且 1/9 < 0.5 → CA 入列
                                                      //（AC 前向是唯一多元素扫描 → pruned=1）
        let mut kc2 = KmerCounter::new(2, false);
        kc2.add_kmer(enc(b"AC"), 9);
        kc2.add_kmer(enc(b"CG"), 3);
        kc2.add_kmer(enc(b"CA"), 1);
        let pruned2 = prune_some_kmers(&mut kc2, 1, 0.0, true, 0.5);
        assert_eq!(pruned2, 1);
        assert_eq!(kc2.get_kmer_count(enc(b"CA")), 0);
        assert_eq!(kc2.get_kmer_count(enc(b"AC")), 9);
        assert_eq!(kc2.get_kmer_count(enc(b"CG")), 3);
    }

    #[test]
    fn prune_error_kmers_disabled() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AC"), 9);
        kc.add_kmer(enc(b"CG"), 3);
        kc.add_kmer(enc(b"CA"), 1);
        let pruned = prune_some_kmers(&mut kc, 1, 0.0, false, 0.5);
        assert_eq!(pruned, 0);
        assert_eq!(kc.get_kmer_count(enc(b"CA")), 1);
    }

    #[test]
    fn prune_second_pass_zero_keys_not_reprocessed() {
        // 已置 0 的键在下一趟跳过（count==0 continue）→ 返回 0
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(7, 5);
        kc.add_kmer(12, 1);
        assert_eq!(prune_some_kmers(&mut kc, 2, 0.0, false, 0.0), 1);
        assert_eq!(prune_some_kmers(&mut kc, 2, 0.0, false, 0.0), 0);
        assert_eq!(kc.get_kmer_count(12), 0);
    }

    // ---- sorted_seed_list ---------------------------------------------------

    #[test]
    fn seed_list_sorted_count_desc_kmer_value_desc_tiebreak() {
        // 7:5 与 12:5 平局 → kmer 值降序（12 在前）; 13:9 最前
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(7, 5);
        kc.add_kmer(12, 5);
        kc.add_kmer(13, 9);
        assert_eq!(sorted_seed_list(&kc, true), vec![(13, 9), (12, 5), (7, 5)]);
    }

    #[test]
    fn seed_list_excludes_zero_counts() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(7, 5);
        kc.add_kmer(12, 3);
        kc.clear_kmer(12); // 惰性删除: 键在表中但 count=0
        assert_eq!(kc.size(), 2);
        assert_eq!(sorted_seed_list(&kc, true), vec![(7, 5)]);
        assert_eq!(sorted_seed_list(&kc, false), vec![(7, 5)]);
    }

    #[test]
    fn seed_list_unsorted_preserves_map_order() {
        // sort=false（PARALLEL）: 保持容器迭代序 = iter_nonzero 快照
        let mut kc = KmerCounter::new(4, false);
        for i in 0..8 {
            kc.add_kmer(100 + i as KmerId, (i % 3) as u32 + 1);
        }
        kc.clear_kmer(103);
        let snapshot: Vec<(KmerId, u32)> = kc.iter_nonzero().collect();
        let unsorted = sorted_seed_list(&kc, false);
        assert_eq!(unsorted, snapshot);
        assert_eq!(unsorted.len(), 7);
        // 同一 map 状态两次调用序稳定
        assert_eq!(sorted_seed_list(&kc, false), unsorted);
    }

    // ======================================================================
    // 贪心延伸核心（IRKE.cpp:719-1371）
    // ======================================================================

    /// 场景 1（线性链）: K=2, AC(5)→CG(5)→GT(5)，种子 AC。
    /// 前向逐轮各延伸 1 窗（MAX_RECURSION=1）: path=[CG,GT]、count=10；
    /// 反向候选 GA/AA/TA/CA 全无 → path=[]。
    /// total = 10 + 0 + 5 = 15；joined = [AC,CG,GT]；序列 "ACGT"（首全串+末碱基）。
    #[test]
    fn linear_chain_k2_full_trace() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AC"), 5);
        kc.add_kmer(enc(b"CG"), 5);
        kc.add_kmer(enc(b"GT"), 5);
        let mut rng = GlibcRand::new(1);
        let (joined, total) =
            build_inchworm_contig_from_seed(enc(b"AC"), &kc, &IrkeParams::default(), &mut rng)
                .unwrap();
        assert_eq!(joined, vec![enc(b"AC"), enc(b"CG"), enc(b"GT")]);
        assert_eq!(total, 15);
        let (seq, covs) = reconstruct_path_sequence(&kc, &joined);
        assert_eq!(seq, "ACGT");
        assert_eq!(covs, vec![5, 5, 5]);
    }

    /// 场景 2（贪心选主）: AC(10) 前向分叉 CG(8) / CA(2) → 选 CG。
    /// 注意反向: AC 的反向候选恰有 CA(2) → 反向 path=[CA]、count=2
    /// （CA 的反向候选 GC/AC/TC/CC 中只有 AC，但 AC 已在前向 visitor 中被跳过）。
    /// total = 8 + 2 + 10 = 20；joined = [CA,AC,CG]；序列 "CACG"。
    #[test]
    fn greedy_selects_dominant_branch() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AC"), 10);
        kc.add_kmer(enc(b"CG"), 8);
        kc.add_kmer(enc(b"CA"), 2);
        let mut rng = GlibcRand::new(1);
        let (joined, total) =
            build_inchworm_contig_from_seed(enc(b"AC"), &kc, &IrkeParams::default(), &mut rng)
                .unwrap();
        assert_eq!(joined, vec![enc(b"CA"), enc(b"AC"), enc(b"CG")]);
        assert_eq!(total, 20);
        let (seq, _) = reconstruct_path_sequence(&kc, &joined);
        assert_eq!(seq, "CACG");
    }

    /// 场景 3（死端平局，rand 不参与）: AC(5) 前向 CG(5)/CA(5) 同分死端。
    /// cap=1: 平局（端点 CG≠CA）、len 1 > 0 → cap=2; cap=2: 两子路径仍各 len 1
    /// （CG/CA 无后续）、len 1 > best_path_length 1 不成立 → 取 paths[0]=CG
    /// （候选 count 平局保持 G,A,T,C 收集序，CG 的 i=0 在 CA 的 i=1 前）。
    /// 无 rand 调用 → 换 rng 种子结果必须相同（两种子首 rand()%2 已证不同）。
    /// 反向仍吃到 CA: total = 5+5+5 = 15; joined = [CA,AC,CG]。
    #[test]
    fn dead_end_tie_resolves_to_first_candidate_without_rand() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AC"), 5);
        kc.add_kmer(enc(b"CG"), 5);
        kc.add_kmer(enc(b"CA"), 5);
        // 前置: 找到首 rand()%2 奇偶不同的种子对，否则"与种子无关"断言无鉴别力
        let mut seed_b = 0u32;
        for s in 2..=8u32 {
            if GlibcRand::new(1).next() % 2 != GlibcRand::new(s).next() % 2 {
                seed_b = s;
                break;
            }
        }
        assert_ne!(seed_b, 0, "未找到首 rand()%2 不同的种子对");
        let p = IrkeParams::default();
        let (j1, t1) =
            build_inchworm_contig_from_seed(enc(b"AC"), &kc, &p, &mut GlibcRand::new(1)).unwrap();
        let (j2, t2) =
            build_inchworm_contig_from_seed(enc(b"AC"), &kc, &p, &mut GlibcRand::new(seed_b))
                .unwrap();
        assert_eq!(j1, j2); // rand 未被调用（或至少不影响结果）
        assert_eq!(t1, t2);
        assert_eq!(j1, vec![enc(b"CA"), enc(b"AC"), enc(b"CG")]);
        assert_eq!(t1, 15);
    }

    /// 场景 4（recurse_cap 递增打破平局）: AC(10) 分叉 CG(5)→GT(9) 与 CA(5)→AG(1)。
    /// cap=1: [CG]:5 vs [CA]:5 平局 → cap=2; cap=2: [GT,CG]:14 vs [AG,CA]:6
    /// → 无平局 → 选深者（无 rand）。
    /// 外层倒序消费深→浅压栈的 [GT,CG] → entire=[CG,GT]。
    /// 反向: CA(5)（AC 已访问跳过）→ total = 14 + 5 + 10 = 29;
    /// joined = [CA,AC,CG,GT]，序列 "CACGT"，counts [5,10,5,9]。
    #[test]
    fn tie_resolved_by_deeper_recursion_no_rand() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AC"), 10);
        kc.add_kmer(enc(b"CG"), 5);
        kc.add_kmer(enc(b"GT"), 9);
        kc.add_kmer(enc(b"CA"), 5);
        kc.add_kmer(enc(b"AG"), 1);
        let mut rng = GlibcRand::new(1);
        let (joined, total) =
            build_inchworm_contig_from_seed(enc(b"AC"), &kc, &IrkeParams::default(), &mut rng)
                .unwrap();
        assert_eq!(joined, vec![enc(b"CA"), enc(b"AC"), enc(b"CG"), enc(b"GT")]);
        assert_eq!(total, 29);
        let (seq, covs) = reconstruct_path_sequence(&kc, &joined);
        assert_eq!(seq, "CACGT");
        assert_eq!(covs, vec![5, 10, 5, 9]);
    }

    /// 场景 5（50 硬停 + rand 打破）——K=7 构造性场景。
    /// 种子 S = s0 + P6（P6 为 6 位数字），两条 50 窗链自 P6 分叉（第 7 位 a7/b7
    /// 不同），链上每窗 count=3、种子 count=10。构造不变量（xorshift 确定性重试
    /// 直到落成，落成后用真实 KmerCounter 复核）:
    /// - 每个 6 位上下文至多被 1 个在表窗口占据（分叉 P6 恰被两链头双占）
    ///   → 链内每窗前向候选恰为下一窗、链尾无候选（无分叉，C3）
    /// - S 的前向候选恰为两链头（C4）；无窗口的末 6 位等于 S 的前 6 位
    ///   → S 无反向候选（C5）；所有窗口值互异且 ≠ S（C6）。
    ///
    /// 追踪: 每 cap c ∈ [1,50] 两条子路径各长 c、count 3c → 每轮真平局且长度
    /// 递增 → recurse_cap 一路 +1 到 50 → `recurse_cap >= MAX_RECURSION_HARD_STOP`
    /// → rand()%2。srand(1) 首值 1804289383（奇）→ idx=1 → 排序后第二条（候选
    /// 平局保持 G,A,T,C 收集序 → 末位编码小者是 paths[0]）→ **末位 digit 较大**
    /// 的链。整链 50 窗全消费（count 150），链尾无候选 → 停。total = 160。
    #[test]
    fn hard_stop_at_50_broken_by_rand() {
        // 测试构造专用的确定性 PRNG（与 glibc_rand 无关）
        struct XorShift(u64);
        impl XorShift {
            fn next(&mut self) -> u64 {
                let mut x = self.0;
                x ^= x >> 12;
                x ^= x << 25;
                x ^= x >> 27;
                self.0 = x;
                x.wrapping_mul(0x2545_F491_4F6C_DD1D)
            }
            fn digit(&mut self) -> u8 {
                (self.next() % 4) as u8
            }
        }
        // digit 0..3 ↔ G,A,T,C（与候选生成 i=0..3 的 G,A,T,C 序一致）
        const BASES: [u8; 4] = [b'G', b'A', b'T', b'C'];
        let val = |w: &[u8]| {
            kmer_to_intval(&w.iter().map(|&d| BASES[d as usize]).collect::<Vec<_>>()).unwrap()
        };
        let pack6 =
            |w: &[u8]| -> u32 { (0..6).fold(0u32, |acc, m| acc | (w[m] as u32) << (2 * m)) };

        let k = 7usize;
        let chain_len = 50usize;
        let mut rng = XorShift(0x9E37_79B9_7F4A_7C15);
        let mut seed_digits: Vec<u8> = Vec::new();
        let mut chains: Vec<Vec<Vec<u8>>> = Vec::new();

        'attempt: for _ in 0..1000 {
            let p6: Vec<u8> = (0..6).map(|_| rng.digit()).collect();
            let a7 = rng.digit();
            let mut b7 = rng.digit();
            while b7 == a7 {
                b7 = rng.digit();
            }
            let s0 = rng.digit();
            let seed = [vec![s0], p6.clone()].concat();
            let seed_first6 = seed[..6].to_vec();

            let mut used: std::collections::HashSet<u64> = std::collections::HashSet::new();
            let mut claimed: std::collections::HashSet<u32> = std::collections::HashSet::new();
            let mut blocked: std::collections::HashSet<u32> = std::collections::HashSet::new();
            used.insert(val(&seed));
            claimed.insert(pack6(&p6)); // 分叉上下文: 恰被两链头双占

            let mut built: Vec<Vec<Vec<u8>>> = Vec::new();
            for &head_d in &[a7, b7] {
                let head = [p6.clone(), vec![head_d]].concat();
                // C6（值不撞）+ C5（链头末 6 位 ≠ S 前 6 位）
                if used.contains(&val(&head)) || head[1..] == seed_first6[..] {
                    continue 'attempt;
                }
                used.insert(val(&head));
                let mut chain: Vec<Vec<u8>> = vec![head];
                while chain.len() < chain_len {
                    let ctx = chain.last().unwrap()[1..].to_vec(); // 后继的 6 位上下文
                    let ctx_key = pack6(&ctx);
                    // 上下文已被其它窗口占据或已被先完成的链尾封锁 → 无延伸
                    if claimed.contains(&ctx_key) || blocked.contains(&ctx_key) {
                        continue 'attempt;
                    }
                    let mut cands = [0u8, 1, 2, 3];
                    for i in (1..4).rev() {
                        let j = (rng.next() % (i as u64 + 1)) as usize;
                        cands.swap(i, j);
                    }
                    let mut extended = false;
                    for d in cands {
                        let next = [ctx.clone(), vec![d]].concat();
                        if used.contains(&val(&next)) || next[1..] == seed_first6[..] {
                            continue;
                        }
                        used.insert(val(&next));
                        claimed.insert(ctx_key);
                        chain.push(next);
                        extended = true;
                        break;
                    }
                    if !extended {
                        continue 'attempt;
                    }
                }
                // 链尾上下文不得已被任何窗口占据（含周期性回撞/先建链的内部上下文）
                if claimed.contains(&pack6(&chain.last().unwrap()[1..])) {
                    continue 'attempt;
                }
                blocked.insert(pack6(&chain.last().unwrap()[1..])); // 封锁: 无后继
                built.push(chain);
            }
            if built.len() == 2 {
                seed_digits = seed;
                chains = built;
                break;
            }
        }
        assert_eq!(chains.len(), 2, "构造重试 1000 次未成功");

        // ---- 用真实机器复核构造不变量（防御性，构造成立时必然通过） ----
        let sval = val(&seed_digits);
        let mut kc = KmerCounter::new(k, false);
        kc.add_kmer(sval, 10);
        for chain in &chains {
            for w in chain {
                kc.add_kmer(val(w), 3);
            }
        }
        assert_eq!(kc.get_kmer_count(sval), 10); // C6
        let mut expect_fwd = vec![(val(&chains[0][0]), 3), (val(&chains[1][0]), 3)];
        if chains[0][0][6] > chains[1][0][6] {
            expect_fwd.swap(0, 1); // count 平局 → G,A,T,C 序
        }
        assert_eq!(kc.get_forward_kmer_candidates(sval), expect_fwd); // C4
        assert!(kc.get_reverse_kmer_candidates(sval).is_empty()); // C5
        for chain in &chains {
            for m in 0..chain_len {
                let expect: Vec<(KmerId, u32)> = if m + 1 < chain_len {
                    vec![(val(&chain[m + 1]), 3)]
                } else {
                    vec![]
                };
                assert_eq!(kc.get_forward_kmer_candidates(val(&chain[m])), expect);
                // C3
            }
        }

        // ---- 运行与断言 ----
        // srand(1) 首值 1804289383（奇）→ rand()%2 = 1 → paths[1] = 末位 digit 大者
        let (hi, lo) = if chains[0][0][6] > chains[1][0][6] {
            (0, 1)
        } else {
            (1, 0)
        };
        let mut rng = GlibcRand::new(1);
        let (joined, total) =
            build_inchworm_contig_from_seed(sval, &kc, &IrkeParams::default(), &mut rng).unwrap();
        let mut expect = vec![sval];
        expect.extend(chains[hi].iter().map(|w| val(w)));
        assert_eq!(joined, expect);
        assert_eq!(total, 150 + 10); // 链 50*3 + 种子
                                     // 鉴别力: 选中的恰不是"未走 rand 时的默认 paths[0]"（末位 digit 小者）
        assert_ne!(joined[1], val(&chains[lo][0]));
        let (seq, covs) = reconstruct_path_sequence(&kc, &joined);
        assert_eq!(seq.len(), k + chain_len); // 首窗全串 + 每窗末碱基
        assert_eq!(covs.len(), 1 + chain_len);
        assert_eq!(covs[0], 10);
        assert!(covs[1..].iter().all(|&c| c == 3));
    }

    /// 轮数守卫（IRKE.cpp:851-853 throw 的移植）: num_total_kmers = counter.size()
    /// 进入时快照，round 递增后 > 快照即 Err。空计数器在 round 1 即触发
    /// （原版日常图不可达——visitor 每轮净增 ≥1 永久节点；CRAWL 模式才会实际触发）。
    #[test]
    fn inchworm_round_exceeding_total_kmers_is_error() {
        let kc = KmerCounter::new(2, false);
        let mut rng = GlibcRand::new(1);
        let err =
            build_inchworm_contig_from_seed(enc(b"AC"), &kc, &IrkeParams::default(), &mut rng)
                .unwrap_err();
        match err {
            CommonError::Inchworm(msg) => assert_eq!(
                msg,
                "Error, inchworm rounds have exceeded the number of possible seed kmers"
            ),
            other => panic!("unexpected error: {other:?}"),
        }
    }

    /// 场景 8（深→浅压栈）: 直接调 inchworm_step，max_recurse=2 深入两层。
    /// 递归返回时每层把**自己**追加在子路径之后 → path = [GT, CG]（最远端在前），
    /// count = 5+5 = 10；种子 depth=0 不入 path 不计 count。
    #[test]
    fn inchworm_step_stacks_deep_first_and_seed_depth0_free() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AC"), 5);
        kc.add_kmer(enc(b"CG"), 5);
        kc.add_kmer(enc(b"GT"), 5);
        let mut visitor = KmerVisitor::new(2, false);
        let mut eliminator = KmerVisitor::new(2, false);
        let mut rng = GlibcRand::new(1);
        let best = inchworm_step(
            &kc,
            Direction::Forward,
            (enc(b"AC"), 5),
            &mut visitor,
            &mut eliminator,
            1,
            0,
            &IrkeParams::default(),
            2,
            &mut rng,
        );
        assert_eq!(best.path, vec![enc(b"GT"), enc(b"CG")]);
        assert_eq!(best.count, 10);
    }

    /// 场景 6/7: is_good_seed_kmer 四段闸门。
    /// 注意熵是 f32 浮点路径——精确边界（熵恰 == min）不做断言（libm 实现差异），
    /// 只测有安全裕度的两侧；整型闸门（coverage）做精确边界（严格 <）。
    #[test]
    fn is_good_seed_kmer_gates() {
        let kc2 = KmerCounter::new(2, false);
        let kc4 = KmerCounter::new(4, false);
        let loose = IrkeParams {
            min_seed_entropy: 0.5,
            ..Default::default() // min_seed_coverage = 2
        };
        // 回文: "AT" K=2 → 6, revcomp(6) = 6（精确等值）
        assert_eq!(enc(b"AT"), 6);
        assert_eq!(revcomp_val(6, 2), 6);
        assert!(!is_good_seed_kmer(&kc2, enc(b"AT"), 5, &loose));
        // K=4 回文 "ACGT"（熵 2.0 也拦）
        assert_eq!(revcomp_val(enc(b"ACGT"), 4), enc(b"ACGT"));
        assert!(!is_good_seed_kmer(&kc4, enc(b"ACGT"), 5, &loose));
        // count == 0
        assert!(!is_good_seed_kmer(&kc2, enc(b"AC"), 0, &loose));
        // 覆盖不足（严格 <）: 1 < 2
        assert!(!is_good_seed_kmer(&kc2, enc(b"AC"), 1, &loose));
        // 恰好 == min_seed_coverage 通过; "AC" 熵 1.0 ≥ 0.5
        assert!(is_good_seed_kmer(&kc2, enc(b"AC"), 2, &loose));
        // 默认 min_seed_entropy=1.5: "AAAA" 熵 0 → 拒; K=2 最高熵 1.0 → 全拒
        let def = IrkeParams::default();
        assert!(!is_good_seed_kmer(&kc4, enc(b"AAAA"), 5, &def));
        assert!(!is_good_seed_kmer(&kc2, enc(b"AC"), 5, &def));
        // 非回文 + 熵 2.0（"ATCG" 的 revcomp = "CGAT" ≠ 自身）→ 过
        assert_ne!(revcomp_val(enc(b"ATCG"), 4), enc(b"ATCG"));
        assert!(is_good_seed_kmer(&kc4, enc(b"ATCG"), 5, &def));
    }

    /// 场景 9: extract_best_seed——`count > best && is_good_seed`（&& 短路、严格 >）。
    #[test]
    fn extract_best_seed_rules() {
        let p = IrkeParams {
            min_seed_entropy: 1.0,
            ..Default::default()
        };
        let mut kc = KmerCounter::new(4, false);
        kc.add_kmer(enc(b"ATAT"), 9); // 回文 → 不合格（count 最高也没用）
        kc.add_kmer(enc(b"ACGA"), 5); // 熵 ≈1.5 ≥ 1.0 → 合格
        kc.add_kmer(enc(b"AAAC"), 7); // 熵 ≈0.81 < 1.0 → 不合格
        assert_eq!(
            extract_best_seed(&[enc(b"ATAT"), enc(b"ACGA"), enc(b"AAAC")], &kc, &p),
            enc(b"ACGA")
        );
        // 全不合格 → 0
        assert_eq!(extract_best_seed(&[enc(b"ATAT"), enc(b"AAAC")], &kc, &p), 0);
        assert_eq!(extract_best_seed(&[], &kc, &p), 0);
        // 平局（严格 > 才更新）→ 路径中更早者
        let mut kc2 = KmerCounter::new(4, false);
        kc2.add_kmer(enc(b"ACGA"), 5);
        kc2.add_kmer(enc(b"AGTA"), 5); // 熵 ≈1.5 ≥ 1.0
        assert_eq!(
            extract_best_seed(&[enc(b"ACGA"), enc(b"AGTA")], &kc2, &p),
            enc(b"ACGA")
        );
        // count==0 的键（0 > 0 不成立）
        let mut kc3 = KmerCounter::new(4, false);
        kc3.add_kmer(enc(b"ACGA"), 0);
        assert_eq!(extract_best_seed(&[enc(b"ACGA")], &kc3, &p), 0);
    }

    /// 场景 10: reconstruct_path_sequence——首 kmer 全串 + 后续各末碱基 + counts 快照。
    #[test]
    fn reconstruct_path_sequence_rules() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AC"), 3);
        kc.add_kmer(enc(b"CG"), 7);
        assert_eq!(
            reconstruct_path_sequence(&kc, &[enc(b"AC"), enc(b"CG")]),
            ("ACG".to_string(), vec![3, 7])
        );
        assert_eq!(reconstruct_path_sequence(&kc, &[]), (String::new(), vec![]));
    }

    /// DS 模式端到端: visitor/candidate 计数全走 canonical。
    /// AC(5)（canonical 7）前向候选 12..15 中只有 CT(14)（canonical max(14,
    /// revcomp="AG"=4) = 14）有计数 → path=[14]; 反向 GA/AA/TA/CA 的 canonical
    /// 均无计数 → 空。joined 用**原始位运算值** [7,14]，GT 与 AC 同键计数 5。
    #[test]
    fn ds_mode_build_walks_canonical_keys() {
        let mut kc = KmerCounter::new(2, true);
        kc.add_kmer(enc(b"AC"), 5);
        kc.add_kmer(enc(b"CT"), 5);
        let mut rng = GlibcRand::new(1);
        let (joined, total) =
            build_inchworm_contig_from_seed(enc(b"AC"), &kc, &IrkeParams::default(), &mut rng)
                .unwrap();
        assert_eq!(joined, vec![enc(b"AC"), enc(b"CT")]);
        assert_eq!(total, 10);
        assert_eq!(kc.get_kmer_count(enc(b"GT")), 5); // 互补链同键
        assert_eq!(kc.get_kmer_count(enc(b"AG")), 5); // CT 的互补链同键
    }

    // ======================================================================
    // 组装主循环（IRKE.cpp:426-716 compute_sequence_assemblies）
    // ======================================================================

    /// K=2 种子参数: 熵上限 1.0 < 默认 1.5 → 放宽到 0.5（K=2 测试的标准宽松参数，
    /// 同 is_good_seed_kmer_gates 的 loose）。
    fn loose_k2_params() -> IrkeParams {
        IrkeParams {
            min_seed_entropy: 0.5,
            ..Default::default()
        }
    }

    /// 主循环端到端（K=2 线性链 AC(5)→CG(5)→GT(5)，种子 AC）:
    /// joined=[AC,CG,GT]、total=15、seq="ACGT"（len 4）→ avg_cov = 15/3+0.5 = 5.5
    /// → 5（double 加 0.5 后截断）。过关（len 4 >= 3 且 5 >= 2）→ 产出 1 条；
    /// 路径全清零 → 其余种子 count=0 被 is_good_seed_kmer 拒 → 无第二条。
    #[test]
    fn compute_assemblies_linear_chain_record_and_clear() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AC"), 5);
        kc.add_kmer(enc(b"CG"), 5);
        kc.add_kmer(enc(b"GT"), 5);
        let aparams = AssemblyParams {
            min_assembly_length: 3,
            ..Default::default()
        };
        let mut out = Vec::new();
        let n = compute_sequence_assemblies(
            &mut kc,
            &loose_k2_params(),
            &aparams,
            &Monitor::default(),
            true,
            &mut out,
        )
        .unwrap();
        assert_eq!(n, 1);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            ">a1;5 total_counts: 15 Seed: 5 K: 2 length: 4\nACGT\n"
        );
        // 路径全部清零（惰性删除: 键在表中、count=0）
        assert_eq!(kc.get_kmer_count(enc(b"AC")), 0);
        assert_eq!(kc.get_kmer_count(enc(b"CG")), 0);
        assert_eq!(kc.get_kmer_count(enc(b"GT")), 0);
        assert_eq!(kc.size(), 3);
    }

    /// 不过关（seq len 4 < min_assembly_length 5）→ 不记录，但**照样清零**：
    /// 后续种子（CG/GT）count=0 → 无输出（IRKE.cpp:564-574 else 分支的清零
    /// 不依赖记录条件）。
    #[test]
    fn compute_assemblies_short_contig_not_recorded_but_cleared() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AC"), 5);
        kc.add_kmer(enc(b"CG"), 5);
        kc.add_kmer(enc(b"GT"), 5);
        let aparams = AssemblyParams {
            min_assembly_length: 5,
            ..Default::default()
        };
        let mut out = Vec::new();
        let n = compute_sequence_assemblies(
            &mut kc,
            &loose_k2_params(),
            &aparams,
            &Monitor::default(),
            true,
            &mut out,
        )
        .unwrap();
        assert_eq!(n, 0);
        assert!(out.is_empty());
        assert_eq!(kc.get_kmer_count(enc(b"AC")), 0); // 仍被清零
        assert_eq!(kc.get_kmer_count(enc(b"GT")), 0);
    }

    /// avg_cov 边界（IRKE.cpp:563 四舍五入）: AC(6)+CG(5) → total 11、kmer 数 2
    /// → 11/2 = 5.5、+0.5 = 6.0 → 6（半数进位）。avg_cov 门槛 = 6 恰好通过
    /// （严格 >=，IRKE.cpp:573）。
    #[test]
    fn compute_assemblies_avg_cov_rounds_half_up_and_gate_is_ge() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AC"), 6);
        kc.add_kmer(enc(b"CG"), 5);
        let aparams = AssemblyParams {
            min_assembly_length: 3,
            min_assembly_coverage: 6,
            ..Default::default()
        };
        let mut out = Vec::new();
        let n = compute_sequence_assemblies(
            &mut kc,
            &loose_k2_params(),
            &aparams,
            &Monitor::default(),
            true,
            &mut out,
        )
        .unwrap();
        assert_eq!(n, 1);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            ">a1;6 total_counts: 11 Seed: 6 K: 2 length: 3\nACG\n"
        );
        // min_assembly_coverage 提到 7 → avg_cov 6 < 7 → 不产出
        let aparams7 = AssemblyParams {
            min_assembly_length: 3,
            min_assembly_coverage: 7,
            ..Default::default()
        };
        let mut out7 = Vec::new();
        let n7 = compute_sequence_assemblies(
            &mut kc,
            &loose_k2_params(),
            &aparams7,
            &Monitor::default(),
            true,
            &mut out7,
        )
        .unwrap();
        assert_eq!(n7, 0);
        assert!(out7.is_empty());
    }

    /// 多 contig 顺序 + Seed 字段是**主循环快照**（IRKE.cpp:546 的 kmer_count，
    /// 非 build 内部重查——本场景两者一致，锁字段位置与序）。
    /// 两条**互不链接**的链（K=2 平面小、回文/邻接到处相碰——按下表手挑）:
    /// - 链 1 "ACT": AC(7):9 - CT(14):9。CT 前向(TG/TA/TT/TC)与反向(GG/TC/CC)全空，
    ///   AC 反向(GA/AA/TA/CA)全空 → 封闭。
    /// - 链 2 "AGT": AG(1):4 - GT(2):4。AG 前向含 GT(链内)，GT 反向含 AG(链内)，
    ///   其余邻接全空。
    /// - 跨链: AC/CT 的全部邻接 ∩ {AG, GT} = ∅（AG/GT 不在 AC/CT 的候选位上）。
    ///
    /// 种子序: count 9 平局 kmer 值降序 → CT(14) 先 → 反向吃到 AC → "ACT" 为 a1
    /// （total 18、avg 9.5→9）;count 4 平局 → GT(2) 先 → "AGT" 为 a2（8、4）。
    #[test]
    fn compute_assemblies_two_contigs_in_seed_count_order() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AC"), 9);
        kc.add_kmer(enc(b"CT"), 9);
        kc.add_kmer(enc(b"AG"), 4);
        kc.add_kmer(enc(b"GT"), 4);
        let aparams = AssemblyParams {
            min_assembly_length: 3,
            ..Default::default()
        };
        let mut out = Vec::new();
        let n = compute_sequence_assemblies(
            &mut kc,
            &loose_k2_params(),
            &aparams,
            &Monitor::default(),
            true,
            &mut out,
        )
        .unwrap();
        assert_eq!(n, 2);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            ">a1;9 total_counts: 18 Seed: 9 K: 2 length: 3\nACT\n\
             >a2;4 total_counts: 8 Seed: 4 K: 2 length: 3\nAGT\n"
        );
    }

    /// 60 列折行（IRKE.cpp:671 add_fasta_seq_line_breaks(sequence, 60)）:
    /// 确定性 LCG 生成 100 碱基（**不能用周期序列**——GATC 周期 4 下 K=25 滑窗
    /// 每 4 位坍缩同一 kmer，链只剩 4 节点），3 条相同读 → 76 窗各 count 3
    /// → total 228、avg 3.5 → 3、100bp 折为 60+40 两行。
    #[test]
    fn compute_assemblies_wraps_sequence_at_60_columns() {
        let mut state: u64 = 0x853c_49e6_748f_ea9b;
        let mut next_base = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"GATC"[(state >> 33) as usize % 4]
        };
        let seq: Vec<u8> = (0..100).map(|_| next_base()).collect();
        // 结构不变量（确定性——失败即换 LCG 种子）: 76 个 25-窗两两不同（无重复
        // kmer）、77 个 24-mer 两两不同（每窗恰 1 个前向/反向候选 → 唯一线性链）、
        // 无回文窗（种子闸门）
        let uniq25: HashSet<Vec<u8>> = (0..=75).map(|i| seq[i..i + 25].to_vec()).collect();
        assert_eq!(uniq25.len(), 76);
        let uniq24: HashSet<Vec<u8>> = (0..=76).map(|i| seq[i..i + 24].to_vec()).collect();
        assert_eq!(uniq24.len(), 77);
        let no_palindrome = (0..=75).all(|i| {
            let w = &seq[i..i + 25];
            *w != w
                .iter()
                .rev()
                .map(|&b| match b {
                    b'G' => b'C',
                    b'A' => b'T',
                    b'T' => b'A',
                    _ => b'G',
                })
                .collect::<Vec<_>>()
        });
        assert!(no_palindrome);

        let text = String::from_utf8(seq.clone()).unwrap();
        let data = format!(">r1\n{text}\n>r2\n{text}\n>r3\n{text}\n");
        let (mut kc, _) = populate_from_reads(data.as_bytes(), 25, false).unwrap();
        assert_eq!(kc.iter_nonzero().count(), 76);
        assert!(kc.iter_nonzero().all(|(_, c)| c == 3));

        let mut out = Vec::new();
        let n = compute_sequence_assemblies(
            &mut kc,
            &IrkeParams::default(), // 熵 ≈2.0 过默认 1.5
            &AssemblyParams::default(),
            &Monitor::default(),
            true,
            &mut out,
        )
        .unwrap();
        assert_eq!(n, 1);
        let out_text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = out_text.lines().collect();
        assert_eq!(lines.len(), 3); // header + 60 列 + 40 列
        assert_eq!(lines[1].len(), 60);
        assert_eq!(lines[2].len(), 40);
        assert_eq!(lines[1..].concat(), text); // 折行拼回原序列
        assert_eq!(
            lines[0],
            ">a1;3 total_counts: 228 Seed: 3 K: 25 length: 100"
        );
    }

    /// sort_seeds=false（PARALLEL 预留）: 不排序保持容器迭代序——与 sorted_seed_list
    /// (counter,false) 一致（此处仅锁定主循环对 sort 开关的透传）。
    #[test]
    fn compute_assemblies_unsorted_passes_through() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AC"), 5);
        kc.add_kmer(enc(b"CG"), 5);
        let aparams = AssemblyParams {
            min_assembly_length: 3,
            ..Default::default()
        };
        let mut out = Vec::new();
        let n = compute_sequence_assemblies(
            &mut kc,
            &loose_k2_params(),
            &aparams,
            &Monitor::default(),
            false,
            &mut out,
        )
        .unwrap();
        assert_eq!(n, 1); // 无论种子序，链上唯一种子（AC/CG 二者其一先到即整链产出）
        assert!(String::from_utf8(out).unwrap().ends_with("ACG\n"));
    }

    /// 组装级 min_connectivity 覆盖 IrkeParams 同名字段（IRKE.cpp:426 形参下传）。
    /// exceeds_min_connectivity 首行 `mc < 1e5 → true` 短路——只有 mc >= 1e5
    /// 才真正启用比值测试（IRKE.cpp:1185-1188）。
    #[test]
    fn compute_assemblies_min_connectivity_from_aparams() {
        // AC(5) → CA(1): ratio min(5,1)/max(5,1) = 0.2
        let build = |mc: f32| -> usize {
            let mut kc = KmerCounter::new(2, false);
            kc.add_kmer(enc(b"AC"), 5);
            kc.add_kmer(enc(b"CA"), 1);
            // params.min_connectivity 恒 0（默认）——只由 aparams 覆盖
            let aparams = AssemblyParams {
                min_connectivity: mc,
                min_assembly_length: 3,
                ..Default::default()
            };
            let mut out = Vec::new();
            compute_sequence_assemblies(
                &mut kc,
                &loose_k2_params(),
                &aparams,
                &Monitor::default(),
                true,
                &mut out,
            )
            .unwrap()
        };
        // 关（0 < 1e5 短路恒过）: CA 是 AC 的**前向**候选（前缀 C?）→ joined=[AC,CA]、
        // seq "ACA" len 3 → 产出 1
        assert_eq!(build(0.0), 1);
        // 开（1e5）: AC↔CA 比值 0.2 < 1e5 双向被拒 → seq "AC" len 2 < 3 → 0
        assert_eq!(build(1e5), 0);
    }

    // ======================================================================
    // PARALLEL_IWORM（IRKE.cpp:504-686 并行分支，Task 6）
    // ======================================================================

    /// 单 chunk（3 种子 < 1000）→ 无并发竞态窗口: PARALLEL + TWO_PHASE 的行为
    /// 完全确定。线性链 AC(5)→CG(5)→GT(5)（同 compute_assemblies_linear_chain_
    /// record_and_clear 场景）: draft joined=[AC,CG,GT] → extract_best_seed 平局
    /// 取路径首现（AC，严格 > 才更新）→ new_seed 即原种子、未被 zapped → 重建
    /// 同路径 → 记录（total 15、avg 5）→ 路径全清零 → 其余种子 count=0 被闸门拒。
    /// 输出与单线程版**逐字节一致**（线性图上 TWO_PHASE 是不动点）。
    #[test]
    fn parallel_single_chunk_two_phase_linear_chain_deterministic() {
        let build = || {
            let mut kc = KmerCounter::new(2, false);
            kc.add_kmer(enc(b"AC"), 5);
            kc.add_kmer(enc(b"CG"), 5);
            kc.add_kmer(enc(b"GT"), 5);
            let aparams = AssemblyParams {
                min_assembly_length: 3,
                ..Default::default()
            };
            let mut out = Vec::new();
            let n = compute_sequence_assemblies_parallel(
                kc,
                &loose_k2_params(),
                &aparams,
                &Monitor::default(),
                true,    // TWO_PHASE（默认）
                Some(4), // 线程数无关——单 chunk 只有一个工作单元
                &mut out,
            )
            .unwrap();
            (n, out)
        };
        let (n, out) = build();
        assert_eq!(n, 1);
        let expected = ">a1;5 total_counts: 15 Seed: 5 K: 2 length: 4\nACGT\n";
        assert_eq!(String::from_utf8(out.clone()).unwrap(), expected);
        // 确定性: 重复运行逐字节一致（无竞态窗口 + 每 chunk srand(1)）
        let (n2, out2) = build();
        assert_eq!(n2, 1);
        assert_eq!(out2, out.as_slice());
    }

    /// 线程数无关: 单 chunk 时 1/2/8 线程与 None（rayon 默认）输出一致。
    #[test]
    fn parallel_thread_count_irrelevant_single_chunk() {
        let run = |n: Option<usize>| {
            let mut kc = KmerCounter::new(2, false);
            kc.add_kmer(enc(b"AC"), 5);
            kc.add_kmer(enc(b"CG"), 5);
            kc.add_kmer(enc(b"GT"), 5);
            let mut out = Vec::new();
            compute_sequence_assemblies_parallel(
                kc,
                &loose_k2_params(),
                &AssemblyParams {
                    min_assembly_length: 3,
                    ..Default::default()
                },
                &Monitor::default(),
                true,
                n,
                &mut out,
            )
            .unwrap();
            out
        };
        let base = run(Some(1));
        for n in [Some(2), Some(8), None] {
            assert_eq!(run(n), base, "单 chunk 下线程数不得影响结果");
        }
    }

    /// TWO_PHASE 语义直测（直接调主循环体，绕开 FxHashMap 种子序）:
    /// 链 AG(4)→GT(4)→TC(9)。draft 自 AG: 前向 [GT,TC]（贪心选 TC 9 > TA 无）、
    /// 反向无 → joined=[AG,GT,TC]、total=17。extract_best_seed 取 TC(9) → 重建
    /// 自 TC（反向 [GT,AG]）→ joined 仍 [AG,GT,TC]、total 仍 17（线性图不动点）。
    /// 断言三件事:
    /// - **Seed 字段 = 主循环快照的原种子 count（4）**——TWO_PHASE 重建自 TC(9)
    ///   但不更新该字段（IRKE.cpp:546 捕获点在 draft 之前）
    /// - total_counts/avg_cov/sequence 来自重建（17、17/3+0.5→6、"AGTC"）
    /// - 记录后 joined_path 全部清零（dashmap 原子置 0）
    #[test]
    fn parallel_two_phase_records_original_seed_count() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AG"), 4);
        kc.add_kmer(enc(b"GT"), 4);
        kc.add_kmer(enc(b"TC"), 9);
        let sync = crate::counter_sync::SyncKmerCounter::from_counter(kc);
        let init_size = sync.size();
        let mut rng = GlibcRand::new(1);
        let rec = parallel_assemble_seed(
            enc(b"AG"),
            &sync,
            &loose_k2_params(),
            &AssemblyParams {
                min_assembly_length: 3,
                ..Default::default()
            },
            true, // TWO_PHASE
            &Monitor::default(),
            init_size,
            &mut rng,
        )
        .unwrap()
        .expect("须产出记录");
        assert_eq!(rec, (17, 6, 4, "AGTC".to_string()));
        // 路径全部清零（惰性删除: 键在、count=0）
        assert_eq!(sync.get_kmer_count(enc(b"AG")), 0);
        assert_eq!(sync.get_kmer_count(enc(b"GT")), 0);
        assert_eq!(sync.get_kmer_count(enc(b"TC")), 0);
        assert_eq!(sync.size(), init_size);
    }

    /// 主循环体对不合格种子（count=0，已被清零）返回 None 且不清零任何键。
    #[test]
    fn parallel_assemble_seed_rejects_zapped_seed() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AC"), 5);
        kc.add_kmer(enc(b"CG"), 5);
        kc.clear_kmer(enc(b"AC")); // 主种子已被"其它线程"清零
        let sync = crate::counter_sync::SyncKmerCounter::from_counter(kc);
        let mut rng = GlibcRand::new(1);
        let rec = parallel_assemble_seed(
            enc(b"AC"),
            &sync,
            &loose_k2_params(),
            &AssemblyParams::default(),
            true,
            &Monitor::default(),
            sync.size(),
            &mut rng,
        )
        .unwrap();
        assert!(rec.is_none());
        assert_eq!(sync.get_kmer_count(enc(b"CG")), 5); // 未被触碰
    }

    /// --SINGLE_PHASE（two_phase=false）: 跳过重建直接记录 draft。
    /// 线性链 AG(4)→GT(4)→TC(9) 无论哪个种子先到（FxHashMap 序任意），draft
    /// 都覆盖全链 → 恰 1 条记录、total 恒 17（链上计数和）、序列 "AGTC";
    /// Seed 字段随首至种子（AG/GT=4 或 TC=9）——不做逐字节断言，锁结构性不变量。
    #[test]
    fn parallel_single_phase_runs() {
        let mut kc = KmerCounter::new(2, false);
        kc.add_kmer(enc(b"AG"), 4);
        kc.add_kmer(enc(b"GT"), 4);
        kc.add_kmer(enc(b"TC"), 9);
        let mut out = Vec::new();
        let n = compute_sequence_assemblies_parallel(
            kc,
            &loose_k2_params(),
            &AssemblyParams {
                min_assembly_length: 3,
                ..Default::default()
            },
            &Monitor::default(),
            false, // --SINGLE_PHASE
            Some(2),
            &mut out,
        )
        .unwrap();
        assert_eq!(n, 1);
        let text = String::from_utf8(out).unwrap();
        assert!(text.starts_with(">a1;"), "{text}");
        assert!(text.contains("total_counts: 17"), "{text}");
        assert!(text.contains("K: 2 length: 4"), "{text}");
        assert!(text.contains("\nAGTC\n"), "{text}");
    }

    /// 多 chunk（种子数 > PARALLEL_CHUNK_SIZE=1000）: rayon 收集保全局种子序、
    /// 去重输出正常。**唯一合格种子** AC(5)（CG/GT 与 1197 个干扰键 count=1 <
    /// min_seed_coverage=2，永不合格 → 永不清零 → 无跨 chunk 竞态窗口）:
    /// AC 的 draft 仍延伸过 CG/GT（延伸只要求 count != 0，无覆盖闸门），
    /// total = 5+1+1 = 7、avg = 7/3+0.5 → 2 → 恰好过 min_assembly_coverage=2。
    /// TWO_PHASE: extract_best 只有 AC 合格 → 重建自 AC → 不动点。
    /// → 输出逐字节确定，两次运行一致。
    #[test]
    fn parallel_multi_chunk_deterministic_when_no_contention() {
        let build = || {
            let mut kc = KmerCounter::new(2, false);
            kc.add_kmer(enc(b"AC"), 5);
            kc.add_kmer(enc(b"CG"), 1);
            kc.add_kmer(enc(b"GT"), 1);
            // count=1 干扰键撑过 chunk 边界（避开 K=2 低值区防与链撞键）
            let mut i: KmerId = 0;
            while kc.size() < 1200 {
                if !matches!(i, 0..=15) {
                    kc.add_kmer(1u64 << 40 | i, 1);
                }
                i += 1;
            }
            assert!(kc.size() >= 1200);
            let mut out = Vec::new();
            let n = compute_sequence_assemblies_parallel(
                kc,
                &loose_k2_params(),
                &AssemblyParams {
                    min_assembly_length: 3,
                    ..Default::default()
                },
                &Monitor::default(),
                true,
                Some(4),
                &mut out,
            )
            .unwrap();
            (n, out)
        };
        let (n, out) = build();
        assert_eq!(n, 1, "仅唯一合格种子产出");
        assert_eq!(
            String::from_utf8(out.clone()).unwrap(),
            ">a1;2 total_counts: 7 Seed: 5 K: 2 length: 4\nACGT\n"
        );
        // 无竞态窗口 → 完全确定
        let (n2, out2) = build();
        assert_eq!(n2, n);
        assert_eq!(out2, out);
    }
}
