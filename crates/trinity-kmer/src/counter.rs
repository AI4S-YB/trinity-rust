//! 并行 k-mer 计数 — jellyfish count 的 Rust 替代（多重集等价，见 xcheck-kmer）
//! 语义: 小写 gatc 视同大写; 窗口遇非 gatc 跳过该 kmer 并继续滑 1 碱基;
//! DS 模式表键 = max(kmer, revcomp(kmer))（sequenceUtil.cpp:376）。
//! 记录切分对齐 jellyfish: 序列段内滤除 `\n`/`\r` 跨行拼接（T2 实测确认）;
//! 行内空白等其他非 gatc 字符仍按非法字符断窗（Trinity 管道产物无行内空白）。

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use trinity_common::kmer::{KmerId, MAX_KMER_LENGTH, RollingKey};

pub type CountMap = FxHashMap<KmerId, u32>;

/// 对一条序列（原始字节）枚举有效 kmer，回调规范化后的键。
/// 镜像 jellyfish: 非 gatc 字符使包含它的窗口失效，窗口继续前滑（不是重启）。
///
/// 性能分层（2026-09-02）: 泛型核心 [`for_each_kmer_of`] 按键类型实例化——
/// u128 紧凑键（k≤64, diginorm 热路径, 24B/槽 vs U256 40B/槽）与 KmerId 全宽。
/// 本函数 = KmerId 实例化（管线对外签名不变）。
pub fn for_each_kmer(seq: &[u8], k: usize, ds: bool, f: impl FnMut(KmerId)) {
    for_each_kmer_of::<KmerId>(seq, k, ds, f)
}

/// 滚动枚举泛型核心; 语义与 jellyfish/原 for_each_kmer 全域一致
/// （DS canonical = max(kmer, revcomp), 非滚动兜底走 K::revcomp 等价式）。
pub fn for_each_kmer_of<K: RollingKey>(seq: &[u8], k: usize, ds: bool, mut f: impl FnMut(K)) {
    // 满宽兜底自限正确（超窗高位随左移自然掉出）; u128 实例的 low_mask 自带 k≤64 断言。
    let mask: K = K::low_mask(k);
    let mut kmer: K = K::zero();
    let mut valid: usize = 0;
    // M7-norm: 滚动 revcomp（revcomp(km<<2|v) == (revcomp(km)>>2) | (v^3)<<(2k-2)），
    // 免去每 kmer 的 k 迭代 revcomp_val 循环; 与 canonical max 语义数值恒等。
    // rev_shift = 2k-2 ≤ 254（k ≤ MAX_KMER_LENGTH=128）恒安全。
    if ds && k <= MAX_KMER_LENGTH {
        let rev_shift: u32 = (2 * k - 2) as u32;
        let mut rev: K = K::zero();
        for &c in seq {
            let v = match c {
                b'G' | b'g' => 0u8,
                b'A' | b'a' => 1,
                b'T' | b't' => 2,
                b'C' | b'c' => 3,
                _ => {
                    valid = 0;
                    kmer = K::zero();
                    rev = K::zero();
                    continue;
                }
            };
            kmer = kmer.roll_forward(v, mask);
            rev = rev.roll_rev(v, rev_shift, mask);
            if valid < k {
                valid += 1;
            }
            if valid == k {
                f(if rev > kmer { rev } else { kmer });
            }
        }
        return;
    }
    for &c in seq {
        let v = match c {
            b'G' | b'g' => 0u8,
            b'A' | b'a' => 1,
            b'T' | b't' => 2,
            b'C' | b'c' => 3,
            _ => {
                valid = 0;
                kmer = K::zero();
                continue;
            }
        };
        kmer = kmer.roll_forward(v, mask);
        if valid < k {
            valid += 1;
        }
        if valid == k {
            let key = if ds {
                // get_ds_kmer_val 的 K 泛型等价式（非滚动兜底路径）
                let rc = kmer.revcomp(k);
                if rc > kmer { rc } else { kmer }
            } else {
                kmer
            };
            f(key);
        }
    }
}

/// 按记录切分 FASTA 字节流: 每条记录 = `>` header 行之后的序列字节段
/// （到下一个行首 `>` 或文件尾），段内滤除 `\n`/`\r`（jellyfish 跨行拼接，
/// 实测对 ">r1\nACGT\nACGTACGT\n" 按连续 12 字符滑窗）。
/// 其余非 gatc 字符不过滤，由 for_each_kmer 按非法字符断窗。
/// 返回 owned 记录（滤行需脱离原缓冲区）。
pub fn split_fasta_sequences(data: &[u8]) -> Vec<Vec<u8>> {
    let mut records: Vec<Vec<u8>> = Vec::new();
    let mut pos = 0;
    while pos < data.len() && data[pos] == b'>' {
        // header 行结束于 \n; 未闭合（无 \n）则该记录无序列
        let Some(hdr_nl) = data[pos..].iter().position(|&b| b == b'\n') else {
            break;
        };
        let seq_start = pos + hdr_nl + 1;
        // 序列段终点 = 下一个行首 '>'（即 (\n, '>') 字节对）或文件尾
        let end = data
            .windows(2)
            .skip(seq_start - 1)
            .position(|w| w[0] == b'\n' && w[1] == b'>')
            .map_or(data.len(), |p| seq_start + p);
        // 段内滤除 \n 与 \r（对齐 jellyfish 跨行拼接），其余字节保留
        let rec: Vec<u8> = data[seq_start..end]
            .iter()
            .copied()
            .filter(|&b| b != b'\n' && b != b'\r')
            .collect();
        records.push(rec);
        pos = end;
    }
    records
}

/// 借用视图版切分（性能）: 返回各记录的 (start, end) 区间 + "是否存在多行
/// 序列段"标记。单行序列段（段内 \n 仅在末尾）下, 对区间切片直接喂
/// for_each_kmer 与滤除 \n 后滑窗完全等价（段内无 \n 可断窗）;
/// 多行段存在时调用方须回退 [`split_fasta_sequences`]（跨行拼接语义）。
pub fn split_fasta_ranges(data: &[u8]) -> (Vec<(usize, usize)>, bool) {
    let mut ranges = Vec::new();
    let mut multiline = false;
    let mut pos = 0;
    while pos < data.len() && data[pos] == b'>' {
        let Some(hdr_nl) = data[pos..].iter().position(|&b| b == b'\n') else {
            break;
        };
        let seq_start = pos + hdr_nl + 1;
        let end = data
            .windows(2)
            .skip(seq_start - 1)
            .position(|w| w[0] == b'\n' && w[1] == b'>')
            .map_or(data.len(), |p| seq_start + p);
        // 段内换行仅允许出现在末尾（单行序列）; 段内 \r 一并滤除语义
        // 由调用方判断——此处只报多行
        let body = &data[seq_start..end];
        let trimmed = body.strip_suffix(b"\n").unwrap_or(body);
        if trimmed.contains(&b'\n') || trimmed.contains(&b'\r') {
            multiline = true;
        }
        ranges.push((seq_start, end));
        pos = end;
    }
    (ranges, multiline)
}

pub struct KmerCountTable;

impl KmerCountTable {
    /// 并行统计 FASTA 字节流的 k-mer 计数（KmerId 口径; 泛型核心见
    /// [`count_fasta_data_of`]——diginorm 热路径按 k 分层用 u64/u128
    /// 紧凑键实例化）。
    pub fn count_fasta_data(data: &[u8], k: usize, ds: bool) -> CountMap {
        Self::count_fasta_data_with_capacity(data, k, ds, 0)
    }

    /// 带每线程预 reserve 的计数（`estimate_hash_size` 接线）:
    /// `per_thread_capacity` 为每线程局部 HashMap 的初始容量（0 = 不预分配,
    /// 即旧行为）。调用方需保证 capacity 总量（× 线程数 × 每项开销）在
    /// max_memory 预算内——orchestrate 侧按 estimate/threads 截断。
    ///
    /// 性能（2026-09-01）: 单行序列走借用 range 视图（免逐记录拷贝,
    /// 切片无内部 \n 时 for_each_kmer 断窗行为与滤除版逐字节等价）;
    /// 检测到多行序列段自动回退 owned 拷贝路径（跨行拼接语义）。
    pub fn count_fasta_data_with_capacity(
        data: &[u8],
        k: usize,
        ds: bool,
        per_thread_capacity: usize,
    ) -> CountMap {
        count_fasta_data_of::<KmerId>(data, k, ds, per_thread_capacity)
    }
}

/// 分片位宽（2^9 = 512 片）: 每片表 ≈ D/512 条（D = distinct 键数）,
/// 大表的随机插入变为片内局部插入; 片表 Mutex 下按缓冲批量冲刷摊薄锁开销。
/// 片数 256→512（2026-09-06 20GB 战役）: hashbrown 桶数 pow2 粒度下,
/// 每片条目越小、容量浪费比越低（QM4 侧表 D/256 ≈ 52 万条 → 2^20 桶
/// 浪费 2.2×; D/512 ≈ 26 万条 → 2^19 桶 1.77×）——片表常驻省 ~20%。
const COUNT_SHARD_BITS: usize = 9;
/// 单任务分片缓冲的冲刷阈值（**字节**）: flush_keys = 该字节数 / 键宽,
/// 下限 [`COUNT_FLUSH_MIN_KEYS`]。缓冲内存上界 = 任务数 × 阈值字节:
/// 经 with_min_len 粒度控制后 QM4 级（~21.5M 记录 → ≤2.6K 任务）≈ 1.3GB
/// （2026-09-06 20GB 战役 M11 自 1MB 再减半; 锁次数经 wide 对照验证无回退）。
/// 2026-09-05 旧值 2M 键下细粒度任务（单任务 <2M 键）**永不触发冲刷**,
/// 键滞留在任务缓冲里等 collect 终点统一冲刷——整个 k-mer 键流被隐式
/// 缓存在内存（QM4 63 亿键 × 8B + Vec 倍增 ≈ 60GB 瞬时, 即实测阶段
/// 峰值 86GB 的主体; 当时的"峰值 ≈ D + 缓冲"推断未经实测）。
const COUNT_FLUSH_BYTES: usize = 1 << 19;
/// 冲刷阈值键数下限（防宽键下字节数除出过小值——锁次数过多）。
const COUNT_FLUSH_MIN_KEYS: usize = 1 << 14;
/// 每任务最少处理记录数（with_min_len）: 任务数上界 = 记录数/该值,
/// 与冲刷字节阈值共同构成缓冲内存界（2026-09-05 内存战役）。
const COUNT_MIN_RANGES: usize = 8192;

/// 当前键宽下的冲刷键数阈值。
fn count_flush_keys_for<K>() -> usize {
    (COUNT_FLUSH_BYTES / std::mem::size_of::<K>()).max(COUNT_FLUSH_MIN_KEYS)
}

fn cnt_profile_enabled() -> bool {
    static P: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *P.get_or_init(|| std::env::var_os("CNT_PROFILE").is_some())
}

/// CNT_PROFILE 行尾的进程峰值 RSS（MB; /proc 不可读 → x）。
fn vm_hwm_mb() -> String {
    trinity_common::timer::vm_field_kb("VmHWM:")
        .map(|kb| format!("{}", kb / 1024))
        .unwrap_or_else(|| "x".into())
}

/// 每线程分片缓冲: 键按片追加, 攒满上限（调用方传入）整批冲刷进片表。
/// drain 清空保留容量, 后续复用免再分配。
struct ShardBufs<K> {
    bufs: Vec<Vec<K>>,
    buffered: usize,
}

impl<K: RollingKey> ShardBufs<K> {
    fn new(shards: usize) -> Self {
        ShardBufs {
            bufs: vec![Vec::new(); shards],
            buffered: 0,
        }
    }

    #[inline]
    fn push(&mut self, shard: usize, key: K) {
        self.bufs[shard].push(key);
        self.buffered += 1;
    }

    fn flush_into(&mut self, tables: &[std::sync::Mutex<FxHashMap<K, u32>>]) {
        for (i, buf) in self.bufs.iter_mut().enumerate() {
            if buf.is_empty() {
                continue;
            }
            let mut t = tables[i].lock().unwrap();
            for key in buf.drain(..) {
                *t.entry(key).or_insert(0) += 1;
            }
        }
        self.buffered = 0;
    }
}

/// 键 → 片号: FxHash 高位（canonical 键的顶部位分布偏斜, 散列混合消除;
/// 低 k 小时键高位恒 0, 散列同样均匀）。
#[inline]
fn shard_of<K: std::hash::Hash>(key: &K, shard_bits: usize) -> usize {
    use std::hash::Hasher;
    let mut h = rustc_hash::FxHasher::default();
    key.hash(&mut h);
    (h.finish() >> (64 - shard_bits)) as usize
}

/// 并行计数泛型核心（分片共享表, 压时间战役 2026-09-05 刀2; 2026-09-05
/// 内存战役缓冲修复 + 轮计数并入）:
/// rayon 按记录分片（with_min_len 钳任务数）→ 每任务分片缓冲攒键
/// （字节定界冲刷, 内存有界）→ 批量冲刷进 256 张 Mutex 片表
/// （片表 ≈ D/256 条, 局部性好）→ 终点按精确容量把片表一次折叠回单表。
///
/// **适用面沿革**: 刀2 交错 A/B（缓冲病未修时分片 k49 240.5→303.6s）
/// 曾裁定轮计数维持局部表+reduce; 2026-09-05 内存战役定位该 A/B 的分片侧
/// 带着未修的键流隐式全量缓存（见 [`COUNT_FLUSH_BYTES`]）, 且局部表形态
/// 峰值 = 线程数 × 0.4~0.7×D（QM4 数十 GB）不可接受 → 缓冲修复后分层
/// 入口（diginorm 与轮计数）统一走本路径, 旧 A/B 作废重裁。
///
/// 相对局部表+reduce 的机制收益: reads 随机分布下各线程局部表收敛到近
/// 全量 distinct 集（合并 ≈ 线程数×D 次随机插入, 局部表叠加把 RSS 顶到
/// 2~3×D×线程数）; 分片后合并工作量恰为 D 次精确容量插入, 峰值内存
/// ≈ 2D + 有界缓冲。键值多重集与任何计数路径等价（片间键不交,
/// 片内 += 可交换）。
pub fn count_fasta_data_of_sharded<K>(
    data: &[u8],
    k: usize,
    ds: bool,
    per_thread_capacity: usize,
) -> FxHashMap<K, u32>
where
    K: RollingKey + Send,
{
    count_fasta_data_sharded::<K>(
        data,
        k,
        ds,
        per_thread_capacity,
        COUNT_SHARD_BITS,
        count_flush_keys_for::<K>(),
        1,
    )
}

/// [`count_fasta_data_of_sharded`] 的折叠过滤变体（刀M12, 2026-09-06
/// 20GB 战役）: 终点折叠只落 count ≥ `fold_min` 的键——与"全量折叠后
/// retain_min"产出**恒等**（同一过滤谓词作用于同一多重集）, 但输出表
/// 只为保留键付出容量（QM4 diginorm m=2 过滤掉单例错误 k-mer ~30-40%
/// distinct, 侧折叠输出 ≈ -40%）, 且折叠插入次数同步减少。
pub fn count_fasta_data_of_sharded_min<K>(
    data: &[u8],
    k: usize,
    ds: bool,
    per_thread_capacity: usize,
    fold_min: u32,
) -> FxHashMap<K, u32>
where
    K: RollingKey + Send,
{
    count_fasta_data_sharded::<K>(
        data,
        k,
        ds,
        per_thread_capacity,
        COUNT_SHARD_BITS,
        count_flush_keys_for::<K>(),
        fold_min,
    )
}

/// 并行计数泛型核心（局部表+reduce, jellyfish 同型策略）: rayon 按记录
/// 分片 → 每任务局部 FxHashMap<K,u32> → reduce 合并。K = u64（k≤32）/
/// u128（k≤64 紧凑键）或 KmerId。**峰值内存 = 在飞任务数 × 局部表大小
/// （随机 reads 下每表收敛到 ~0.4-0.7×D）——QM4 k49 实测数十 GB**, 故
/// 管线热路径（diginorm 与轮计数分层入口）2026-09-05 起统一走
/// [`count_fasta_data_of_sharded`]; 本路径保留给 k>64 全宽层与公共 API
/// （等价性对拍基准）。
pub fn count_fasta_data_of<K>(
    data: &[u8],
    k: usize,
    ds: bool,
    per_thread_capacity: usize,
) -> FxHashMap<K, u32>
where
    K: RollingKey + Send,
{
    let cap = per_thread_capacity.min(1 << 22); // 单线程上限 ~4M 槽，防 rehash 峰值即可
    let t0 = std::time::Instant::now();
    let (ranges, multiline) = split_fasta_ranges(data);
    // fold（非 map）: 每个 rayon 任务一张局部表跨 range 复用——map+collect
    // 会给每个 range 各建一张 cap 槽初始表, orchestrate 传入非零 per_thread
    // （~4M 槽 × 数百 range）时内存爆至数百 GB（2026-09-05 管线 OOM 实证）。
    let merged = if multiline {
        split_fasta_sequences(data)
            .par_iter()
            .fold(
                || FxHashMap::<K, u32>::with_capacity_and_hasher(cap, Default::default()),
                |mut local, rec| {
                    for_each_kmer_of::<K>(rec, k, ds, |key| {
                        *local.entry(key).or_insert(0) += 1;
                    });
                    local
                },
            )
            .reduce(FxHashMap::default, merge_counts)
    } else {
        ranges
            .par_iter()
            .fold(
                || FxHashMap::<K, u32>::with_capacity_and_hasher(cap, Default::default()),
                |mut local, &(s, e)| {
                    for_each_kmer_of::<K>(&data[s..e], k, ds, |key| {
                        *local.entry(key).or_insert(0) += 1;
                    });
                    local
                },
            )
            .reduce(FxHashMap::default, merge_counts)
    };
    if cnt_profile_enabled() {
        eprintln!(
            "[cnt] k={k} ds={ds}: locals+reduce 合计 {:.1}s, VmHWM {}MB",
            t0.elapsed().as_secs_f64(),
            vm_hwm_mb()
        );
    }
    merged
}

/// 分片计数实现（flush_keys 参数化供测试强制多轮冲刷; fold_min = 终点
/// 折叠的计数下限——1 即无过滤, 与全量折叠+retain_min(min) 恒等）。
fn count_fasta_data_sharded<K>(
    data: &[u8],
    k: usize,
    ds: bool,
    per_thread_capacity: usize,
    shard_bits: usize,
    flush_keys: usize,
    fold_min: u32,
) -> FxHashMap<K, u32>
where
    K: RollingKey + Send,
{
    use std::sync::Mutex;
    let t0 = std::time::Instant::now();
    let shards = 1usize << shard_bits;
    // per_thread_capacity 原语义（每线程局部表初始容量）均摊到片表作初始提示
    let per_shard_cap = (per_thread_capacity.min(1 << 22) / shards).max(64);
    let tables: Vec<Mutex<FxHashMap<K, u32>>> = (0..shards)
        .map(|_| {
            Mutex::new(FxHashMap::with_capacity_and_hasher(
                per_shard_cap,
                Default::default(),
            ))
        })
        .collect();

    let (ranges, multiline) = split_fasta_ranges(data);
    let run_range = |buf: &mut ShardBufs<K>, slice: &[u8]| {
        for_each_kmer_of::<K>(slice, k, ds, |key| {
            buf.push(shard_of(&key, shard_bits), key);
        });
        if buf.buffered >= flush_keys {
            buf.flush_into(&tables);
        }
    };
    // fold 携带每 rayon 任务的缓冲跨记录复用, collect 收回残余（collect 后
    // 逐个冲刷, 不丢尾部）; multiline 记录粒度同型。with_min_len 钳任务数
    // 上界 → 缓冲内存界 = 任务数 × COUNT_FLUSH_KEYS（2026-09-05 内存战役）
    let bufs: Vec<ShardBufs<K>> = if multiline {
        split_fasta_sequences(data)
            .par_iter()
            .with_min_len(COUNT_MIN_RANGES)
            .fold(
                || ShardBufs::new(shards),
                |mut buf, rec| {
                    run_range(&mut buf, rec);
                    buf
                },
            )
            .collect()
    } else {
        ranges
            .par_iter()
            .with_min_len(COUNT_MIN_RANGES)
            .fold(
                || ShardBufs::new(shards),
                |mut buf, &(s, e)| {
                    run_range(&mut buf, &data[s..e]);
                    buf
                },
            )
            .collect()
    };
    for mut buf in bufs {
        buf.flush_into(&tables);
    }
    let t_count = t0.elapsed();

    let t1 = std::time::Instant::now();
    // 终点折叠: 片间键不交（分片由键哈希决定）→ 每键恰好一次插入;
    // 容量按全量预留（pow2 对齐下与过滤后实际差异通常落同一桶数档）。
    // fold_min > 1 时只落达限键——与全量折叠后 retain_min(fold_min) 恒等
    // （刀M12: 过滤谓词同、多重集同; 输出表免为被滤键付容量）。
    let total: usize = tables.iter().map(|t| t.lock().unwrap().len()).sum();
    let mut out = FxHashMap::with_capacity_and_hasher(total, Default::default());
    let mut kept = 0usize;
    for t in tables {
        for (key, cnt) in t.into_inner().unwrap() {
            if cnt >= fold_min {
                out.insert(key, cnt);
                kept += 1;
            }
        }
    }
    if cnt_profile_enabled() {
        if fold_min > 1 {
            eprintln!(
                "[cnt] k={k} ds={ds}: count(shard) {:.1}s, merge(final-fold {total} 键, ≥{fold_min} 保留 {kept}) {:.1}s, VmHWM {}MB",
                t_count.as_secs_f64(),
                t1.elapsed().as_secs_f64(),
                vm_hwm_mb()
            );
        } else {
            eprintln!(
                "[cnt] k={k} ds={ds}: count(shard) {:.1}s, merge(final-fold {total} 键) {:.1}s, VmHWM {}MB",
                t_count.as_secs_f64(),
                t1.elapsed().as_secs_f64(),
                vm_hwm_mb()
            );
        }
    }
    out
}

/// 窄键表 → CountMap（KmerId）一次性 widen（压时间战役 2026-09-04 刀1）:
/// 精确预留容量（无生长 rehash）, 键值多重集原样平移——窄键域内
/// `K::to_wide()` 与 KmerId 低 limb 数值恒等（kmer.rs 数值一致性单测）。
fn widen_to_count_map<K: RollingKey>(m: FxHashMap<K, u32>) -> CountMap {
    let mut out = FxHashMap::with_capacity_and_hasher(m.len(), Default::default());
    for (key, cnt) in m {
        out.insert(key.to_wide(), cnt);
    }
    out
}

/// inchworm 轮计数的分层入口: k≤32 用 u64 键、k≤64 用 u128 键在窄键空间
/// 完成分片共享表计数（终点一次 widen 回 CountMap）——管线既有消费方
/// （trusted 注入 / counts_bin 落盘 / KmerCounter::from_count_map）接口与
/// 语义不变。等价性: 各窄键层与 KmerId 层走同一泛型代码, 键数值域内恒等
/// → 计数多重集相等（单测 count_layered_matches_wide 全 k 对拍）。
/// 迭代序差异不构成回归: 并行合并序本就非确定（inchworm 并行平局序
/// 非确定, backlog 已立案）。
///
/// 路径选择（2026-09-05 内存战役重裁）: 刀2 曾按"32 线程单调用交错 A/B
/// k49 240.5→303.6s 回退"裁定轮计数维持局部表+reduce——但该 A/B 的分片
/// 侧带着未修的缓冲病（细粒度任务永不冲刷, 键流隐式全量缓存, 分配器
/// 狂奔——见 [`COUNT_FLUSH_BYTES`] 注记）, 回退归因被污染。局部表形态
/// 的峰值 = 线程数 × 0.4~0.7×D（QM4 k49 实测数十 GB, 内存战役不可接受）;
/// 缓冲修复后的分片形态峰值 ≈ D×2 + 有界缓冲。重 A/B 见 opt-log。
pub fn count_fasta_data_layered(data: &[u8], k: usize, ds: bool, per_thread_capacity: usize) -> CountMap {
    if k <= 32 {
        widen_to_count_map(count_fasta_data_of_sharded::<u64>(
            data,
            k,
            ds,
            per_thread_capacity,
        ))
    } else if k <= 64 {
        widen_to_count_map(count_fasta_data_of_sharded::<u128>(
            data,
            k,
            ds,
            per_thread_capacity,
        ))
    } else {
        count_fasta_data_of::<KmerId>(data, k, ds, per_thread_capacity)
    }
}

/// 33≤k≤64 的轮计数窄键直通（刀M8, 2026-09-06 内存战役）: 返回 u128 计数表
/// 本体, 不做终点 widen——管线窄链（save_counts_bin128 流式落盘 + trusted
/// 注入 + prune + SyncKmerCounter128 PARALLEL 目录）全程 u128, KmerId 宽表
/// （QM4 ~2 亿键 × 40B/槽, hashbrown pow2 容量放大后 ~10.7GB）不再被物化。
/// 键值多重集与 count_fasta_data_layered 的 widen 前身恒等（同一 sharded
/// 计数, 仅是否 widen 的差别）; 域外 k 返回 None（调用方走宽链）。
pub fn count_fasta_data_layered128(
    data: &[u8],
    k: usize,
    ds: bool,
    per_thread_capacity: usize,
) -> Option<FxHashMap<u128, u32>> {
    if (33..=64).contains(&k) {
        Some(count_fasta_data_of_sharded::<u128>(
            data,
            k,
            ds,
            per_thread_capacity,
        ))
    } else {
        None
    }
}

/// save_counts_bin 的 u128 窄表流式变体（刀M8）: 逐条 to_wide 后同格式写出
/// ——与宽链落盘字节恒等（k≤64 域内 widen 是纯零扩展平移, 表迭代序同一）。
pub fn save_counts_bin128(
    path: &std::path::Path,
    counts: &FxHashMap<u128, u32>,
    k: usize,
    ds: bool,
) -> Result<(), trinity_common::error::CommonError> {
    use std::io::Write;
    let fh = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::with_capacity(1 << 23, fh);
    w.write_all(b"TRM3")?;
    w.write_all(&(k as u32).to_le_bytes())?;
    w.write_all(&[ds as u8])?;
    let mut kb = [0u8; 32];
    let mut cb = [0u8; 4];
    for (&key, &c) in counts.iter() {
        kb.copy_from_slice(&KmerId { hi: 0, lo: key }.to_le_bytes());
        w.write_all(&kb)?;
        cb.copy_from_slice(&c.to_le_bytes());
        w.write_all(&cb)?;
    }
    w.flush()?;
    Ok(())
}

/// inject_trusted_kmers 的 u128 窄表变体（刀M8 窄链注入点）: 语义同宽版
/// ——表中已存在的键不动计数, 缺席键 count=1 presence-only; k≤64 域内
/// u128 键与 KmerId 键数值恒等。
pub fn inject_trusted_kmers128(
    counts: &mut FxHashMap<u128, u32>,
    trusted_fasta: &[u8],
    k: usize,
    ds: bool,
) -> usize {
    let mut inserted = 0usize;
    for rec in split_fasta_sequences(trusted_fasta) {
        for_each_kmer_of::<u128>(&rec, k, ds, |key| {
            if let std::collections::hash_map::Entry::Vacant(e) = counts.entry(key) {
                e.insert(1);
                inserted += 1;
            }
        });
    }
    inserted
}

/// reduce 合并: b 并入 a（计数累加）。小表并入大表（rehash 次数最少化）。
pub fn merge_counts<K: RollingKey>(
    mut a: FxHashMap<K, u32>,
    mut b: FxHashMap<K, u32>,
) -> FxHashMap<K, u32> {
    if a.len() < b.len() {
        std::mem::swap(&mut a, &mut b);
    }
    for (key, cnt) in b {
        *a.entry(key).or_insert(0) += cnt;
    }
    a
}

/// trusted contig 注入（rnaSPAdes construction.cpp:233 的 Trinity 移植——
/// 上一轮组装的简化 contig 并入本轮 k-mer 计数流: **贡献序列、不贡献覆盖**）。
/// 表中已存在的键**不动其计数**; 缺席的键以 count=1 插入（presence-only:
/// 能被 min_count=1 的剪枝保留, 且在 inchworm 贪心延伸中优先级最低——
/// 恰是"只补序列不扰覆盖"）。注入 k-mer 不获任何豁免: 后续熵过滤等剪枝
/// 对它们与 reads k-mer 同等处理（忠实接受低复杂度注入 kmer 被剪掉）。
/// 调用方须保证 `k`/`ds` 与计数表的构建参数一致（注入发生在
/// prune_some_kmers 之前的最终表上）。返回新插入的键数（日志用）。
pub fn inject_trusted_kmers(counts: &mut CountMap, trusted_fasta: &[u8], k: usize, ds: bool) -> usize {
    let mut inserted = 0usize;
    for rec in split_fasta_sequences(trusted_fasta) {
        for_each_kmer(&rec, k, ds, |key| {
            // entry().or_insert(1): 已存在（reads 已覆盖或上轮已注入）→ 不动;
            // 缺席 → count=1 占位（presence-only）。
            if let std::collections::hash_map::Entry::Vacant(e) = counts.entry(key) {
                e.insert(1);
                inserted += 1;
            }
        });
    }
    inserted
}

/// mer_counts 二进制记录布局: `U256 key, u32 count`（小端, 36B/条）, 头部魔数
/// `TRM3` + u32 k + u8 ds（读回时校验）。M7 消 dump 直桥的落盘格式,
/// 管线（orchestrate）与离线工具（xtask export-edges）共用。
///
/// 格式沿革: u64 时代 `TRMC` + 8B key（12B/条）; u128 扩容（2026-08-30）升
/// `TRM2` + 16B key（20B/条）; U256 扩容（2026-09-02）升 `TRM3` + 32B key。
/// 键宽变化后旧文件静默错位解析的风险不可接受 → 每次魔数升版硬拒绝
/// （"bad header" 报错）, 不做双格式读（counts_bin 是易失缓存, 重算即可）。
///
/// 流式写（8MB BufWriter）: QM4 级 ~2 亿键 × 36B ≈ 7GB——旧整表 Vec 物化
/// 与计数表叠加成 15GB 瞬时（2026-09-05 内存战役）; 字节序与逐条顺序不变。
pub fn save_counts_bin(
    path: &std::path::Path,
    counts: &CountMap,
    k: usize,
    ds: bool,
) -> Result<(), trinity_common::error::CommonError> {
    use std::io::Write;
    let fh = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::with_capacity(1 << 23, fh);
    w.write_all(b"TRM3")?;
    w.write_all(&(k as u32).to_le_bytes())?;
    w.write_all(&[ds as u8])?;
    let mut kb = [0u8; 32];
    let mut cb = [0u8; 4];
    for (&key, &c) in counts.iter() {
        kb.copy_from_slice(&key.to_le_bytes());
        w.write_all(&kb)?;
        cb.copy_from_slice(&c.to_le_bytes());
        w.write_all(&cb)?;
    }
    w.flush()?;
    Ok(())
}

/// save_counts_bin 的读回侧; k/ds 不匹配或头损坏报错（`TRMC`/`TRM2` 等
/// 旧格式同走 bad header 硬拒绝——见 save_counts_bin 的格式沿革注记）。
pub fn load_counts_bin(
    path: &std::path::Path,
    k: usize,
    ds: bool,
) -> Result<CountMap, trinity_common::error::CommonError> {
    use trinity_common::error::CommonError;
    let data = std::fs::read(path)?;
    if data.len() < 9 || &data[..4] != b"TRM3" {
        return Err(CommonError::Parse(format!(
            "corrupt k-mer table {}: bad header",
            path.display()
        )));
    }
    let k_in = u32::from_le_bytes([data[4], data[5], data[6], data[7]]) as usize;
    let ds_in = data[8] != 0;
    if k_in != k || ds_in != ds {
        return Err(CommonError::Parse(format!(
            "k-mer table {} was built with different k/DS settings",
            path.display()
        )));
    }
    let mut counts = CountMap::default();
    let mut i = 9;
    while i + 36 <= data.len() {
        let key = KmerId::from_le_bytes(data[i..i + 32].try_into().unwrap());
        let c = u32::from_le_bytes(data[i + 32..i + 36].try_into().unwrap());
        counts.insert(key, c);
        i += 36;
    }
    Ok(counts)
}

/// Trinity:2598-2604 hash 大小护栏（信息性容量提示用，语义与原版一致）
pub fn estimate_hash_size(max_memory: u64, read_file_size: u64) -> u64 {
    let mut h = if max_memory > read_file_size {
        (max_memory - read_file_size) / 7
    } else {
        0
    };
    if h < 100_000_000 || read_file_size < 5_000_000_000 {
        h = 100_000_000;
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;
    use trinity_common::kmer::{decode_kmer_from_intval, get_ds_kmer_val, kmer_to_intval};

    fn counts_of(seq: &str, k: usize, ds: bool) -> Vec<(String, u32)> {
        let mut m: CountMap = FxHashMap::default();
        for_each_kmer(seq.as_bytes(), k, ds, |key| *m.entry(key).or_insert(0) += 1);
        let mut v: Vec<_> = m
            .into_iter()
            .map(|(key, c)| {
                (
                    String::from_utf8(decode_kmer_from_intval(key, k)).unwrap(),
                    c,
                )
            })
            .collect();
        v.sort();
        v
    }

    /// u128 紧凑键计数与 KmerId 全宽计数多重集等价（k≤64 域, hi 恒 0;
    /// 性能分层 2026-09-02 的数值正确性锁）。
    #[test]
    fn u128_compact_counting_matches_kmerid() {
        let seq = b">x\nACGTACGTNNAGGCACGTACGTACGTagct\n>y\nAAAAACCCCCGGGGGTTTTTACGTACGTA";
        for k in [15usize, 25, 31, 63, 64] {
            for ds in [true, false] {
                let wide = KmerCountTable::count_fasta_data(seq, k, ds);
                let narrow = count_fasta_data_of::<u128>(seq, k, ds, 0);
                assert_eq!(wide.len(), narrow.len(), "k={k} ds={ds}");
                for (key, c) in &narrow {
                    assert_eq!(wide.get(&key.to_wide()), Some(c), "k={k} ds={ds}");
                }
            }
        }
    }

    /// 分片计数（刀2, 2026-09-05）: 强制多轮冲刷（flush_keys=64）与常规路径
    /// （大上限, 单轮冲刷）在 u64/u128/KmerId 三层键上多重集等价——冲刷轮数
    /// 只改插入时序, 片内 += 可交换、片间键不交。
    #[test]
    fn count_sharded_forced_flush_matches() {
        let fa = b">x\nACGTACGTNNAGGCACGTACGTACGTagctACGTACGTACGT\n>y\nAAAAACCCCCGGGGGTTTTTACGTACGTACGTACGT";
        let multiline: Vec<u8> = {
            let mut v = b">m\n".to_vec();
            for _ in 0..12 {
                v.extend_from_slice(b"ACGTACGT");
                v.push(b'\n');
            }
            v.extend_from_slice(b"ACGTACGTNNAGG");
            v
        };
        for (tag, seq) in [("single", &fa[..]), ("multiline", &multiline[..])] {
            for k in [7usize, 25, 33, 63] {
                for ds in [true, false] {
                    let normal = count_fasta_data_of::<KmerId>(seq, k, ds, 0);
                    let tiny_flush = count_fasta_data_sharded::<KmerId>(
                        seq, k, ds, 0, COUNT_SHARD_BITS, 64, 1,
                    );
                    assert_eq!(tiny_flush.len(), normal.len(), "{tag} k={k} ds={ds}");
                    for (key, c) in &normal {
                        assert_eq!(tiny_flush.get(key), Some(c), "{tag} k={k} ds={ds}");
                    }
                    // 刀M12: 折叠过滤与"全量折叠后 retain_min"恒等
                    let fold_min2 = count_fasta_data_sharded::<KmerId>(
                        seq, k, ds, 0, COUNT_SHARD_BITS, 64, 2,
                    );
                    let mut retained = tiny_flush.clone();
                    retained.retain(|_, c| *c >= 2);
                    assert_eq!(fold_min2.len(), retained.len(), "{tag} k={k} ds={ds}");
                    for (key, c) in &retained {
                        assert_eq!(fold_min2.get(key), Some(c), "{tag} k={k} ds={ds}");
                    }
                    if k <= 32 {
                        let narrow = count_fasta_data_sharded::<u64>(
                            seq, k, ds, 0, COUNT_SHARD_BITS, 64, 1,
                        );
                        for (key, c) in &normal {
                            debug_assert!(key.hi == 0);
                            assert_eq!(
                                narrow.get(&(key.lo as u64)),
                                Some(c),
                                "{tag} k={k} ds={ds}"
                            );
                        }
                    }
                }
            }
        }
    }

    /// 分层计数入口（count_fasta_data_layered, 压时间战役 2026-09-04 刀1）与
    /// KmerId 直计在全部键层（u64 k≤32 / u128 k≤64 / KmerId 直通）多重集等价。
    /// 多行记录 fixture 额外覆盖 split_fasta_ranges 的 multiline 回退分支。
    #[test]
    fn count_layered_matches_wide() {
        let single = b">x\nACGTACGTNNAGGCACGTACGTACGTagct\n>y\nAAAAACCCCCGGGGGTTTTTACGTACGTA";
        let multiline: Vec<u8> = {
            let mut v = b">m\n".to_vec();
            for _ in 0..10 {
                v.extend_from_slice(b"ACGTACGT"); // 80 base, 分两行
                v.push(b'\n');
            }
            v.extend_from_slice(b"ACGTACGTNNAGG");
            v
        };
        for (tag, seq) in [("single", &single[..]), ("multiline", &multiline[..])] {
            for k in [7usize, 25, 31, 32, 33, 49, 63, 64, 65] {
                for ds in [true, false] {
                    let wide = KmerCountTable::count_fasta_data(seq, k, ds);
                    let layered = count_fasta_data_layered(seq, k, ds, 0);
                    assert_eq!(layered.len(), wide.len(), "{tag} k={k} ds={ds}");
                    for (key, c) in &wide {
                        assert_eq!(layered.get(key), Some(c), "{tag} k={k} ds={ds}");
                    }
                }
            }
        }
    }

    #[test]
    fn basic_ss_counting() {
        assert_eq!(
            counts_of("ACGT", 2, false),
            vec![("AC".into(), 1), ("CG".into(), 1), ("GT".into(), 1)]
        );
    }

    #[test]
    fn ds_collapses_revcomp() {
        assert_eq!(counts_of("AA TT", 2, true), vec![("TT".into(), 2)]);
        assert_eq!(counts_of("GT", 2, true), vec![("AC".into(), 1)]);
    }

    #[test]
    fn lowercase_and_n_window_slide() {
        assert_eq!(counts_of("acgt", 2, false)[0].1, 1);
        assert_eq!(counts_of("ACNT", 2, false), vec![("AC".into(), 1)]);
        // 手推 ACGNNT: 无 N 窗口 = AC(0,1), CG(1,2) → 计两个（计划稿漏了 CG，已修正;
        // 只有包含非法字符的窗口失效，G 在 N 前的 CG 窗口有效）
        assert_eq!(
            counts_of("ACGNNT", 2, false),
            vec![("AC".into(), 1), ("CG".into(), 1)]
        );
        // 滑动 vs 重启判别: N 后窗口恢复即计 CA（重启实现会漏）
        // 手推 AANCAA: 无 N 的 2-mer 窗口 = AA(0,1), CA(3,4), AA(4,5) → AA×2, CA×1
        assert_eq!(
            counts_of("AANCAA", 2, false),
            vec![("AA".into(), 2), ("CA".into(), 1)]
        );
    }

    #[test]
    fn short_seq_no_kmer() {
        assert_eq!(counts_of("AC", 2, false).len(), 1);
        assert!(counts_of("A", 2, false).is_empty());
    }

    #[test]
    fn estimate_hash_size_mirrors_trinity() {
        // Trinity:2601 护栏是 `||`: read_file_size < 5e9 即触发 100e6 下限
        // 手推: 1e9 < 5e9 → 护栏触发（计划稿的 9e9/7 期望未过护栏，已修正）
        assert_eq!(
            estimate_hash_size(10_000_000_000, 1_000_000_000),
            100_000_000
        );
        assert_eq!(
            estimate_hash_size(10_000_000_000, 4_000_000_000),
            100_000_000
        );
        // 大输入文件（≥5e9）且估计值 ≥100e6 时 /7 估计保留
        assert_eq!(
            estimate_hash_size(10_000_000_000, 5_000_000_000),
            5_000_000_000 / 7
        );
        // 下溢防护分支: max_memory ≤ read_file_size → 0 → 护栏
        assert_eq!(estimate_hash_size(1_000, 2_000), 100_000_000);
    }

    #[test]
    fn count_fasta_data_end_to_end() {
        let fa = b">r1\nACGTACGT\n>r2\nACGTACGT\n>r3\nTTTT\n";
        let m = KmerCountTable::count_fasta_data(fa, 4, false);
        let key = |s: &str| kmer_to_intval(s.as_bytes()).unwrap();
        assert_eq!(m[&key("ACGT")], 4);
        assert_eq!(m[&key("CGTA")], 2);
        assert_eq!(m[&key("TTTT")], 1);
        let mds = KmerCountTable::count_fasta_data(fa, 4, true);
        assert_eq!(mds[&key("TTTT")], 1);
    }

    #[test]
    fn count_fasta_multiline_and_crlf() {
        // jellyfish 跨行拼接: 序列段内 \n/\r 滤除后滑窗（实测 jellyfish 2.3.x 对
        // ">r1\nACGT\nACGTACGT\n" 输出 ACGT=3 CGTA=2 GTAC=2 TACG=2，即 12 字符连续滑窗）。
        // 其余非 gatc 字符（N、行内空格等）仍断窗（T1 滑动语义测试锁定）。
        let fa = b">r1 desc\nACGT\nACGT\r\nTT\n";
        let m = KmerCountTable::count_fasta_data(fa, 4, false);
        let key = |s: &str| kmer_to_intval(s.as_bytes()).unwrap();
        // 滤 \n\r 后序列 = ACGTACGTTT（10 字符）→ 7 个窗口:
        // ACGT CGTA GTAC TACG ACGT CGTT GTTT
        assert_eq!(m[&key("ACGT")], 2);
        assert_eq!(m[&key("CGTA")], 1);
        assert_eq!(m[&key("GTAC")], 1);
        assert_eq!(m[&key("TACG")], 1);
        assert_eq!(m[&key("CGTT")], 1);
        assert_eq!(m[&key("GTTT")], 1);
        assert_eq!(m.len(), 6);
    }

    /// M7-norm oracle: 滚动 revcomp 路径产出的 DS 键与逐 kmer get_ds_kmer_val 恒等
    /// （含 N 断窗、小写、滑窗恢复、k=1/k=2 边界; u128 扩容后加 k=33/37——
    /// 滚动路径在 u64 时代对 k>32 回落, 现对 k≤64 全程生效）。
    #[test]
    fn rolling_revcomp_matches_get_ds_kmer_val() {
        let seqs = [
            "ACGTACGTACGTACGTACGTACGTACG",
            "acgtNNNNacgtacgtacgtacgtacg",
            "GGGGG",
            "TTTTGGGGCCCCAAAATTTTGGGGCACGTN",
            // 70bp 长序列覆盖 k=33/37 的多窗滑窗（含首尾窗口的高位碱基）
            "ACGTACGTTTGCAGTCGATCGTACGATCGATGCATCGATCGTACGATGCTAGCTAGCATCGA",
        ];
        for &seq in &seqs {
            for k in [1usize, 2, 4, 25, 33, 37] {
                if k > seq.len() {
                    continue;
                }
                let mut rolling: Vec<KmerId> = Vec::new();
                for_each_kmer(seq.as_bytes(), k, true, |key| rolling.push(key));
                // naive: 逐窗口 kmer_to_intval + get_ds_kmer_val
                let b = seq.as_bytes();
                let mut naive: Vec<KmerId> = Vec::new();
                let mut i = 0;
                while i + k <= b.len() {
                    if let Ok(v) = trinity_common::kmer::kmer_to_intval(&b[i..i + k]) {
                        naive.push(get_ds_kmer_val(v, k));
                    }
                    i += 1;
                }
                assert_eq!(rolling, naive, "seq={seq} k={k}");
            }
        }
    }

    /// counts_bin TRM2（u128 键）落盘/读回往返, 含 k=37 的高位键;
    /// u64 时代 TRMC 旧格式被硬拒绝（防静默错位解析）。
    #[test]
    fn counts_bin_trm2_roundtrip_and_legacy_rejected() {
        let dir = std::env::temp_dir().join(format!("trmc2_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bin = dir.join("m.bin");
        let mut m: CountMap = FxHashMap::default();
        let hi = kmer_to_intval(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTA").unwrap(); // 37-mer
        assert!(hi > KmerId::from(u64::MAX as u128)); // 确为 u64 放不下的键
        m.insert(hi, 7);
        m.insert(kmer_to_intval(b"AAAA").unwrap(), 3);
        save_counts_bin(&bin, &m, 37, true).unwrap();
        let back = load_counts_bin(&bin, 37, true).unwrap();
        assert_eq!(back.len(), 2);
        assert_eq!(back[&hi], 7);
        // k/ds 不匹配 → 报错
        assert!(load_counts_bin(&bin, 25, true).is_err());
        assert!(load_counts_bin(&bin, 37, false).is_err());
        // 旧 TRMC 魔数（u64 版）→ bad header 硬拒绝
        let mut legacy = b"TRMC".to_vec();
        legacy.extend_from_slice(&37u32.to_le_bytes());
        legacy.push(1);
        std::fs::write(&bin, &legacy).unwrap();
        assert!(load_counts_bin(&bin, 37, true).is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn count_fasta_joins_lines_like_jellyfish() {
        // 与 jellyfish 对拍原例（T2 冒烟固化）: 跨行窗口 CGTA/GTAC/TACG 只在拼接语义下出现
        let fa = b">r1\nACGT\nACGTACGT\n";
        let m = KmerCountTable::count_fasta_data(fa, 4, false);
        let key = |s: &str| kmer_to_intval(s.as_bytes()).unwrap();
        assert_eq!(m[&key("ACGT")], 3);
        assert_eq!(m[&key("CGTA")], 2);
        assert_eq!(m[&key("GTAC")], 2);
        assert_eq!(m[&key("TACG")], 2);
        assert_eq!(m.len(), 4);
    }
}

#[cfg(test)]
mod capacity_tests {
    use super::*;

    /// 预 reserve 版本与默认版本多重集等价。
    #[test]
    fn with_capacity_matches_default() {
        let fa = b">r1\nACGTACGTACGTACGTACGTACGTACGTACGT\n>r2\nTTTTACGTACGTACGTACGTACGTACGTACG\n";
        let a = KmerCountTable::count_fasta_data(fa, 25, false);
        let b = KmerCountTable::count_fasta_data_with_capacity(fa, 25, false, 1 << 20);
        assert_eq!(a.len(), b.len());
        assert!(a.iter().all(|(k, v)| b.get(k) == Some(v)));
    }
}

#[cfg(test)]
mod trusted_inject_tests {
    use super::*;
    use trinity_common::kmer::kmer_to_intval;

    fn key(s: &str) -> KmerId {
        kmer_to_intval(s.as_bytes()).unwrap()
    }

    /// 核心语义: 已存在的键不动计数, 缺席的键以 count=1 插入。
    #[test]
    fn existing_untouched_absent_inserted_one() {
        // reads 表: ACG=7; trusted 携带 ACG/CGT/GTA/TAC（k=3, SS）
        let mut m: CountMap = FxHashMap::default();
        m.insert(key("ACG"), 7);
        let n = inject_trusted_kmers(&mut m, b">c1\nACGTAC\n", 3, false);
        assert_eq!(n, 3); // CGT/GTA/TAC 新插入
        assert_eq!(m[&key("ACG")], 7); // 已存在 → 不动（不累加、不覆盖）
        assert_eq!(m[&key("CGT")], 1);
        assert_eq!(m[&key("GTA")], 1);
        assert_eq!(m[&key("TAC")], 1);
        assert_eq!(m.len(), 4);
    }

    /// 幂等: 重复注入（resume 两条构建路径汇合后）不改变表。
    #[test]
    fn injection_is_idempotent() {
        let mut m: CountMap = FxHashMap::default();
        let fa = b">c1\nACGTAC\n";
        assert_eq!(inject_trusted_kmers(&mut m, fa, 3, false), 4);
        // 第二次注入: 全部已存在（含上轮注入的 count=1）→ 一个都不动
        assert_eq!(inject_trusted_kmers(&mut m, fa, 3, false), 0);
        assert_eq!(m.len(), 4);
        assert!(m.values().all(|&c| c == 1));
    }

    /// DS canonical 折叠与计数路径同键: reads 里的 CTT 与 trusted 里的 AAG
    /// 折叠为同一 DS 键 → 已存在不动。
    #[test]
    fn ds_canonical_folding_matches_reads_table() {
        let mut m: CountMap = FxHashMap::default();
        for_each_kmer(b"CTT", 3, true, |kk| *m.entry(kk).or_insert(0) += 5);
        assert_eq!(m.len(), 1);
        let n = inject_trusted_kmers(&mut m, b">c1\nAAG\n", 3, true);
        assert_eq!(n, 0); // AAG 的 DS canonical == CTT 的 canonical → 已存在
        assert_eq!(m.values().copied().collect::<Vec<_>>(), vec![5]);
    }

    /// 注入按**本轮** k 枚举（trusted 由上一轮较小 k 产出）;
    /// 记录短于 k → 无 k-mer 可注入。
    #[test]
    fn inject_uses_current_k_not_contig_build_k() {
        let mut m: CountMap = FxHashMap::default();
        assert_eq!(inject_trusted_kmers(&mut m, b">c1\nACGT\n", 5, false), 0);
        assert!(m.is_empty());
        assert_eq!(inject_trusted_kmers(&mut m, b">c1\nACGT\n", 3, false), 2);
    }

    /// 多记录 + N 断窗语义与计数路径一致（同一 for_each_kmer）。
    #[test]
    fn multi_record_and_n_breaks() {
        let mut m: CountMap = FxHashMap::default();
        let fa = b">c1\nACG\n>c2 desc\nTTNAC\n";
        let n = inject_trusted_kmers(&mut m, fa, 2, false);
        // c1: AC CG; c2 "TTNAC": TT 有效, TN/NA 含 N 失效, N 后窗口恢复 → AC 有效
        assert_eq!(n, 3);
        assert_eq!(m[&key("AC")], 1);
        assert_eq!(m[&key("CG")], 1);
        assert_eq!(m[&key("TT")], 1);
    }
}
