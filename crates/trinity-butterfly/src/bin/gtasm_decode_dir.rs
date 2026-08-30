//! gtasm 目录级解码工具——开发期对拍用（不进产线）。
//!
//! 用法: `gtasm_decode_dir <Component_bins 目录> <weights.gtw> <support.gts1> <out.fa> [min_len] [visited_walk=1]`
//! 遍历 `<dir>/Cbin*/c*.graph.out`（cid 取文件名 c<N> 段, 同 python load_comps）,
//! 逐 comp 跑 gtasm 全链, 汇总写 FASTA。`visited_walk=0` 回退原型语义
//! （visited 硬排除）, 用于与 python 原型逐字节对拍。

use trinity_butterfly::{run_component_gtasm, EdgeSupport, GtModel};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 4 || args.len() > 6 {
        eprintln!(
            "用法: gtasm_decode_dir <Component_bins 目录> <weights.gtw> <support.gts1> <out.fa> [min_len] [visited_walk=1]"
        );
        std::process::exit(2);
    }
    let min_len: usize = args.get(4).map(|s| s.parse().unwrap()).unwrap_or(200);
    let allow_visited_walk = args.get(5).map(|s| s != "0").unwrap_or(true);
    let model = GtModel::load(std::path::Path::new(&args[1]))
        .unwrap_or_else(|e| panic!("load weights: {e}"));
    let support = EdgeSupport::load(std::path::Path::new(&args[2]))
        .unwrap_or_else(|e| panic!("load support: {e}"));

    let mut files: Vec<_> = std::fs::read_dir(&args[0])
        .unwrap()
        .flat_map(|d| {
            let d = d.unwrap().path();
            if d.is_dir() {
                std::fs::read_dir(d)
                    .unwrap()
                    .map(|f| f.unwrap().path())
                    .collect::<Vec<_>>()
            } else {
                Vec::new()
            }
        })
        .filter(|p| p.file_name().unwrap().to_string_lossy().ends_with(".graph.out"))
        .collect();
    files.sort();

    let mut out = String::new();
    let mut n_contigs = 0usize;
    for f in &files {
        let fname = f.file_name().unwrap().to_string_lossy();
        // python: cid = int(f.rsplit("/c",1)[1].split(".")[0])
        let cid: u64 = fname
            .strip_prefix('c')
            .and_then(|s| s.split('.').next())
            .and_then(|s| s.parse().ok())
            .unwrap_or_else(|| panic!("bad graph file name: {fname}"));
        let text = std::fs::read_to_string(f).unwrap();
        let name = fname.trim_end_matches(".graph.out").to_string();
        let fasta = run_component_gtasm(&text, cid, &name, min_len, &model, &support, allow_visited_walk)
            .unwrap_or_else(|e| panic!("comp {cid}: {e}"));
        n_contigs += fasta.matches('>').count();
        out.push_str(&fasta);
    }
    std::fs::write(&args[3], &out).unwrap();
    eprintln!("解码完成: {} comps → {} contigs → {}", files.len(), n_contigs, args[3]);
}
