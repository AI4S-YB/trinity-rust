//! kmer_count 分层窄键计数微基准（压时间战役 2026-09-04 刀1）: 对
//! both.fa（管线真实输入）以 orchestrate 同参（estimate_hash_size 口径）
//! A/B 计时 KmerId 直计 vs 分层窄键（k≤32 u64 / k≤64 u128 + 终点 widen）,
//! 并做全表逐键对拍。
//!
//! 用法: count_bench <both.fa> <k> [threads=32]

use std::time::Instant;

use trinity_kmer::counter::{count_fasta_data_layered, estimate_hash_size, KmerCountTable};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let fa = args.get(1).expect("用法: count_bench <both.fa> <k> [threads]");
    let k: usize = args.get(2).expect("缺 k").parse().unwrap();
    let threads: usize = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(32);
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global();
    let data = std::fs::read(fa).expect("读 fa 失败");
    eprintln!("输入 {:.2}GB, k={k}, {threads} 线程", data.len() as f64 / 1e9);
    // orchestrate 同参: estimate / 36B 条目 / 线程数（函数内部再截 4M 槽上限）
    let per_thread =
        (estimate_hash_size(64 << 30, data.len() as u64) / 36 / threads as u64) as usize;

    let t = Instant::now();
    let wide = KmerCountTable::count_fasta_data_with_capacity(&data, k, true, per_thread);
    let t_wide = t.elapsed().as_secs_f64();

    let t = Instant::now();
    let layered = count_fasta_data_layered(&data, k, true, per_thread);
    let t_layered = t.elapsed().as_secs_f64();

    eprintln!("KmerId 直计: {t_wide:.1}s ({})", wide.len());
    eprintln!("分层窄键+widow: {t_layered:.1}s ({})", layered.len());
    eprintln!(
        "加速 {:.2}×, 净省 {:.1}s",
        t_wide / t_layered,
        t_wide - t_layered
    );

    // 全表对拍: 键+计数多重集相等
    assert_eq!(wide.len(), layered.len(), "表长不一致");
    let t = Instant::now();
    for (key, c) in &wide {
        assert_eq!(layered.get(key), Some(c), "键 {key:?} 计数不一致");
    }
    eprintln!("全表逐键对拍通过（{} 键, {:.1}s）", wide.len(), t.elapsed().as_secs_f64());
    println!("k={k}\twide={t_wide:.1}s\tlayered={t_layered:.1}s\tspeedup={:.2}", t_wide / t_layered);
}
