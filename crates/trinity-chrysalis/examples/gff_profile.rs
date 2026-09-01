//! gff 阶段性能 profile harness（GFF_PROFILE=1 打印分段计时）。
//! 用法: cargo run --release -p trinity-chrysalis --example gff_profile -- <iworm.fa> <reads1.fa> [reads2.fa...]
use std::time::Instant;

use trinity_chrysalis::dna_vector::{read_fasta, stream_fasta_records};
use trinity_chrysalis::graph_from_fasta::{graph_from_fasta, GffParams};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let iworm = read_fasta(&args[1]).unwrap();
    eprintln!("iworm contigs: {}", iworm.len());
    let mut reads = Vec::new();
    for p in &args[2..] {
        let t = Instant::now();
        let text = std::fs::read_to_string(p).unwrap();
        reads.extend(stream_fasta_records(&text));
        eprintln!("reads +{} ({:?})", p, t.elapsed());
    }
    eprintln!("reads total: {}", reads.len());
    let t = Instant::now();
    let p = GffParams {
        threads: 32,
        ..GffParams::default()
    };
    let out = graph_from_fasta(&iworm, &reads, &p).unwrap();
    eprintln!("gff total: {:?}, out {} bytes", t.elapsed(), out.len());
}
