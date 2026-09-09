//! regroup 输出层探针: 组装产物 fasta 直接跑严格家族拆分（不经主线管线,
//! 机制见 trinity-cli `regroup` 模块文档）。
//! 用法: regroup <in.fa> [cov] [--min-len N] [--min-supp M] [reads1.fa reads2.fq.gz ...]
//! 产物落输入旁同名侧车: 剥最后扩展名接 `.principal.fa` / `.family_map`
//! （mem9.Trinity.fasta → mem9.Trinity.principal.fa + mem9.Trinity.family_map）。
//! 给了 reads 文件即走支持度通道（代表 = argmax(支持度, 长度, id);
//! fq/fa、明文与 gzip 均可, 第二参数能解析为浮点则视为 cov, 否则视为 reads）。
//! --min-len/--min-supp = principal 输出截断门（探针默认 0 = 关, 传 0 值关回
//! 旧行为; 主线管线转正默认 500/10; min-supp 需要 reads 通道, 语义见 regroup
//! 模块文档）。
use std::path::PathBuf;
use std::time::Instant;

use trinity_cli::regroup::{regroup_fasta, sidecar_path, stats_line, RegroupParams};

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    if argv.len() < 2 {
        eprintln!(
            "用法: regroup <in.fa> [cov] [--min-len N] [--min-supp M] [reads1.fa reads2.fq.gz ...]"
        );
        std::process::exit(2);
    }
    let in_fa = PathBuf::from(&argv[1]);
    // 消费任意位置的 --min-len <N> / --min-supp <M>（成对取值, 移出 argv）
    let mut rest: Vec<String> = Vec::new();
    let (mut min_len, mut min_supp) = (0usize, 0u64);
    let mut i = 2;
    while i < argv.len() {
        match argv[i].as_str() {
            "--min-len" | "--min-supp" => {
                let v = argv.get(i + 1).unwrap_or_else(|| {
                    eprintln!("{} 需要一个数值参数", argv[i]);
                    std::process::exit(2);
                });
                if argv[i] == "--min-len" {
                    min_len = v.parse().unwrap_or_else(|_| {
                        eprintln!("--min-len 需要非负整数, 得到 {v}");
                        std::process::exit(2);
                    });
                } else {
                    min_supp = v.parse().unwrap_or_else(|_| {
                        eprintln!("--min-supp 需要非负整数, 得到 {v}");
                        std::process::exit(2);
                    });
                }
                i += 2;
            }
            _ => {
                rest.push(argv[i].clone());
                i += 1;
            }
        }
    }
    if min_supp > 0 && rest.is_empty() {
        eprintln!("WARNING: --min-supp 无 reads 文件时票数全 0, principal 将被剔空");
    }
    let (cov, reads): (f32, Vec<String>) = match rest.first() {
        Some(v) => match v.parse::<f32>() {
            Ok(c) => (c, rest[1..].to_vec()),
            Err(_) => (0.7, rest.to_vec()),
        },
        None => (0.7, Vec::new()),
    };
    let reads: Vec<PathBuf> = reads.into_iter().map(PathBuf::from).collect();
    let params = RegroupParams {
        cov,
        min_principal_len: min_len,
        min_principal_supp: min_supp,
        ..RegroupParams::default()
    };
    let principal = sidecar_path(&in_fa, "principal.fa");
    let map = sidecar_path(&in_fa, "family_map");
    let t = Instant::now();
    let st = regroup_fasta(&in_fa, &principal, &map, &reads, &params).expect("regroup 失败");
    eprintln!("{}", stats_line(&st));
    eprintln!(
        "principal {} / family_map {}（cov {cov}, min-len {min_len}, min-supp {min_supp}, reads 文件 {}, 耗时 {:?}）",
        principal.display(),
        map.display(),
        reads.len(),
        t.elapsed()
    );
}
