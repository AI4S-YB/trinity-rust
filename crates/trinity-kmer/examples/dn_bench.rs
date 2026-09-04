//! diginorm 归一化性能微基准（压时间战役 2026-09-04）: 直接对 QM4 原始
//! PE fq.gz 跑 diginorm::run（管线同参: DS + fa_direct + 默认覆盖参数）,
//! 配合内置 DN_PROFILE=1 分段计时定位热点。
//!
//! 用法: DN_PROFILE=1 dn_bench <left.fq.gz> <right.fq.gz> <out_dir>

use std::time::Instant;

use trinity_kmer::diginorm::{run, DigiNormParams, ReadsInput};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let (l, r, out) = (
        args.get(1).expect("用法: dn_bench <left.fq.gz> <right.fq.gz> <out_dir>"),
        args.get(2).expect("缺 right"),
        args.get(3).expect("缺 out_dir"),
    );
    // 对齐 bench 口径: 32 线程（不设则 rayon 默认吃满 128 核, 计时不可比）
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(32)
        .build_global();
    let params = DigiNormParams {
        ds: true,
        fa_direct: true,
        ..Default::default()
    };
    let reads = ReadsInput::Paired(vec![l.clone().into()], vec![r.clone().into()]);
    let t = Instant::now();
    let outs = run(&params, &reads, std::path::Path::new(out)).expect("diginorm 失败");
    println!(
        "diginorm 总耗时 {:.1}s → left {} / right {}",
        t.elapsed().as_secs_f64(),
        outs.left.unwrap_or_default().display(),
        outs.right.unwrap_or_default().display(),
    );
}
