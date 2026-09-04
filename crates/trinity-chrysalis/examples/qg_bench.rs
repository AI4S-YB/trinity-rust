//! QuantifyGraph 性能微基准（压时间战役 2026-09-04）: 真实 QM4 reads 抽样
//! 构造单组件大图, 计时 quantify_graph_mem 整体耗时; 对拍模式（--check 目录）
//! 落盘产物 sha256, 供实现改动前后逐字节一致性验证。
//!
//! 用法:
//!   qg_bench <reads.fq.gz> [read 数, 默认 200000] [边数, 默认 20000]
//!   qg_bench --check <dir> <reads.fq.gz> [...]  // 产物写 <dir>/g.out <dir>/r.out

use std::time::Instant;

use trinity_chrysalis::quantify::{quantify_graph_mem, QgParams};

/// 简易确定性 PRNG（xorshift64*）, 保证两次运行构造相同输入。
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let (check_dir, rest): (Option<String>, Vec<String>) = if args.get(1).map(|s| s.as_str()) == Some("--check") {
        (Some(args[2].clone()), args[3..].to_vec())
    } else {
        (None, args[1..].to_vec())
    };
    let reads_path = rest.first().expect("用法: qg_bench [--check DIR] <reads.fq.gz> [n_reads] [n_edges]");
    let n_reads: usize = rest.get(1).map(|s| s.parse().unwrap()).unwrap_or(200_000);
    let n_edges: usize = rest.get(2).map(|s| s.parse().unwrap()).unwrap_or(20_000);

    // 1. 读 reads（gz 透明）, 截前 n_reads 条
    let t0 = Instant::now();
    let reader = trinity_common::io_util::open_maybe_gz(std::path::Path::new(reads_path)).unwrap();
    let mut fq = trinity_common::fastq::FastqReader::new(reader);
    let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(n_reads);
    let mut fasta = String::new();
    while seqs.len() < n_reads {
        match fq.next_record().unwrap() {
            Some(rec) => {
                fasta.push_str(&format!(">read{}\t100%\n", seqs.len()));
                fasta.push_str(&rec.sequence);
                fasta.push('\n');
                seqs.push(rec.sequence.clone().into_bytes());
            }
            None => break,
        }
    }
    println!(
        "reads: {} 条（{} bp 均长）, fasta {} MB, 读入 {:.2}s",
        seqs.len(),
        seqs.iter().map(|s| s.len()).sum::<usize>() / seqs.len().max(1),
        fasta.len() / 1024 / 1024,
        t0.elapsed().as_secs_f64()
    );

    // 2. 构造组件图: 每条边从随机 read 随机位置取 25-mer（保证高命中）,
    //    拆 first[prev]（1 nt）+ kmer 列（24 nt）; 偶尔注入 N（次群路径）。
    let mut rng = Rng(0x9E3779B97F4A7C15);
    let mut graph = String::from("Component 0\n");
    graph.push_str("0\t-1\t1\tAAAAAAAAAAAAAAAAAAAAAAAA\n");
    let mut node = 1u64;
    for _ in 0..n_edges {
        let ri = rng.below(seqs.len());
        let s = &seqs[ri];
        if s.len() < 40 {
            continue;
        }
        let pos = rng.below(s.len() - 39);
        let sub: Vec<u8> = s[pos..pos + 25].to_vec();
        let kmer: Vec<u8> = sub[1..].to_vec();
        let mut kmer = kmer;
        if rng.below(200) == 0 {
            let i = rng.below(kmer.len());
            kmer[i] = b'N'; // 0.5% 注 N
        }
        graph.push_str(&format!(
            "{}\t{}\t7\t{}\n",
            node,
            node - 1,
            String::from_utf8_lossy(&kmer)
        ));
        node += 1;
    }
    println!(
        "graph: {} 行（{} MB）",
        node,
        graph.len() / 1024 / 1024
    );

    // 3. 跑 quantify_graph_mem
    let t1 = Instant::now();
    let out = quantify_graph_mem(
        &graph,
        fasta.as_bytes(),
        &QgParams { k: 24, strand: false, max_reads: -1 },
    )
    .unwrap();
    let dt = t1.elapsed().as_secs_f64();
    println!(
        "quantify: {:.2}s → graph.out {} MB, reads.out {} MB",
        dt,
        out.0.len() / 1024 / 1024,
        out.1.len() / 1024 / 1024
    );

    // 4. 指纹（对拍）
    let gh = fnv(out.0.as_bytes());
    let rh = fnv(out.1.as_bytes());
    println!("FNV1a64 graph.out = {gh:016x}, reads.out = {rh:016x}");
    if let Some(dir) = check_dir {
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(std::path::Path::new(&dir).join("g.out"), &out.0).unwrap();
        std::fs::write(std::path::Path::new(&dir).join("r.out"), &out.1).unwrap();
        println!("产物已落盘 {dir}/g.out {dir}/r.out");
    }
}

fn fnv(b: &[u8]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for &c in b {
        h ^= c as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}
