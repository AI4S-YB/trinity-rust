//! P5-T3 集成测试: 编排主线小端到端（合成平铺 read → 全程 → Trinity.fasta
//! 与 gene_trans_map、.ok 齐全）、断点续跑（删产物留 .ok → 阶段跳过）、
//! butterfly 组件池失败隔离与参数别名。

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use trinity_cli::args::parse_args;
use trinity_cli::butterfly_pool::{run_butterfly_pool, ButterflyPoolParams};
use trinity_cli::orchestrate::run_trinity;

fn tmpdir(tag: &str) -> PathBuf {
    static N: AtomicU64 = AtomicU64::new(0);
    let p = std::env::temp_dir().join(format!(
        "trinity-cli-t3-{}-{}-{}",
        tag,
        std::process::id(),
        N.fetch_add(1, Ordering::SeqCst)
    ));
    fs::create_dir_all(&p).unwrap();
    p
}

/// 确定性合成输入: 600nt 参考（xorshift 位串）→ 25 条 100bp 平铺 read
/// （步长 20）——inchworm 产出 ~580nt contig → 单组件 → 单转录本。
fn write_synthetic_reads(dir: &Path) -> PathBuf {
    let mut x: u64 = 0x243F6A8885A308D3;
    let mut ref_seq = String::with_capacity(600);
    for _ in 0..600 {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        ref_seq.push(match x % 4 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        });
    }
    let fa = dir.join("synth.fa");
    let mut text = String::new();
    for (i, off) in (0..500usize).step_by(20).enumerate() {
        text.push_str(&format!(
            ">read{i} component=c0\n{}\n",
            &ref_seq[off..off + 100]
        ));
    }
    fs::write(&fa, text).unwrap();
    fa
}

fn synth_args(fa: &Path, outdir: &Path) -> trinity_cli::TrinityArgs {
    let argv: Vec<String> = [
        "--seqType",
        "fa",
        "--single",
        fa.to_str().unwrap(),
        "--max_memory",
        "1G",
        "--output",
        outdir.to_str().unwrap(),
        "--no_normalize_reads",
        "--CPU",
        "2",
        "--min_contig_length",
        "100",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    parse_args(&argv).unwrap()
}

fn synth_args_opt(fa: &Path, outdir: &Path, no_cleanup: bool) -> trinity_cli::TrinityArgs {
    let mut a = synth_args(fa, outdir);
    a.no_cleanup = no_cleanup;
    a
}

fn fa_record_count(p: &Path) -> usize {
    fs::read_to_string(p)
        .unwrap()
        .lines()
        .filter(|l| l.starts_with('>'))
        .count()
}

#[test]
fn end_to_end_small() {
    let dir = tmpdir("e2e");
    let fa = write_synthetic_reads(&dir);
    let outdir = dir.join("trinity_out");
    let final_fa = run_trinity(&synth_args_opt(&fa, &outdir, true)).unwrap();

    // 最终 fasta 存在非空且含转录本
    assert!(
        final_fa.ends_with("trinity_out.Trinity.fasta"),
        "{final_fa:?}"
    );
    let n = fa_record_count(&final_fa);
    assert!(n >= 1, "no transcripts in final fasta");

    // gene_trans_map 行数 == fasta 记录数; 每行 gene\ttrans
    let map = final_fa
        .parent()
        .unwrap()
        .join("trinity_out.Trinity.fasta.gene_trans_map");
    let map_text = fs::read_to_string(&map).unwrap();
    let map_lines: Vec<&str> = map_text.lines().filter(|l| !l.is_empty()).collect();
    assert_eq!(map_lines.len(), n);
    for l in &map_lines {
        let mut f = l.split('\t');
        let gene = f.next().unwrap();
        let trans = f.next().unwrap();
        assert!(gene.ends_with("_g1") || gene.contains("_g"));
        assert!(trans.starts_with(gene), "{gene} vs {trans}");
        assert!(trans.strip_prefix(gene).unwrap().starts_with("_i"));
    }

    // .ok / finished 齐全（no_normalize → 无 normalization.ok）
    for ok in [
        ".jellyfish_count.25.mincov1.asm.ok",
        ".jellyfish_dump.25.mincov1.asm.ok",
        ".iworm.25.asm.ok",
        "iworm_renamed.25.asm.ok",
        ".quantify_graph.ok",
        ".butterfly.ok",
        "inchworm.DS.fa.finished",
        "single.fa.ok",
    ] {
        assert!(outdir.join(ok).exists(), "missing checkpoint: {ok}");
    }
    assert!(!outdir.join("insilico_read_normalization").exists());
    // 单端 → single.fa（非 both.fa）
    assert!(outdir.join("single.fa").exists());
}

#[test]
fn resume_skips_completed_stage() {
    let dir = tmpdir("resume");
    let fa = write_synthetic_reads(&dir);
    let outdir = dir.join("trinity_out");
    run_trinity(&synth_args_opt(&fa, &outdir, true)).unwrap();

    // 删中间产物（保留 .ok）→ 重跑: dump 阶段被跳过（文件不重建）,
    // inchworm 之后各阶段亦全部跳过, 最终仍成功产出。
    let kmer_fa = outdir.join("jellyfish.kmers.25.asm.fa");
    let iworm_fa = outdir.join("inchworm.DS.fa");
    fs::remove_file(&kmer_fa).unwrap();
    let iworm_mtime = fs::metadata(&iworm_fa).unwrap().modified().unwrap();

    let final_fa = run_trinity(&synth_args_opt(&fa, &outdir, true)).unwrap();
    assert!(final_fa.exists());
    assert!(!kmer_fa.exists(), "dump stage should be skipped on resume");
    assert_eq!(
        fs::metadata(&iworm_fa).unwrap().modified().unwrap(),
        iworm_mtime,
        "inchworm stage should be skipped on resume"
    );
}

#[test]
fn butterfly_pool_isolates_failures() {
    let dir = tmpdir("bflypool");
    // 两个坏组件（graph.out 图文本非法）: 池应两者都跑完、都记录、聚合报错,
    // 而不是第一个 Err 即中断。
    let comp_dir = dir.join("Cbin0");
    fs::create_dir_all(&comp_dir).unwrap();
    for id in [1u64, 2u64] {
        fs::write(comp_dir.join(format!("c{id}.graph.out")), "not a graph\n").unwrap();
        fs::write(comp_dir.join(format!("c{id}.graph.reads")), ">r\nACGT\n").unwrap();
    }
    let listing = vec![(1u64, comp_dir.join("c1")), (2u64, comp_dir.join("c2"))];
    let err = run_butterfly_pool(
        &listing,
        &dir,
        &ButterflyPoolParams {
            cpu: 2,
            min_contig_length: 100,
            group_pairs_distance: 500,
            stack_size_mb: 16,
        },
    )
    .unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("2 component(s)"), "{msg}");
    // 失败清单含两个 base 名
    let failed = fs::read_to_string(dir.join("failed_butterfly_commands.txt")).unwrap();
    assert!(failed.contains("c1"));
    assert!(failed.contains("c2"));
}

// ---------------------------------------------------------------- 参数别名

#[test]
fn args_aliases_kmer_size_and_min_kmer_cov() {
    let mk = |extra: &[&str]| -> trinity_cli::TrinityArgs {
        let mut argv: Vec<String> = vec![
            "--seqType".into(),
            "fa".into(),
            "--single".into(),
            "s.fa".into(),
            "--max_memory".into(),
            "1G".into(),
            "--output".into(),
            "out".into(),
        ];
        argv.extend(extra.iter().map(|s| s.to_string()));
        parse_args(&argv).unwrap()
    };
    assert_eq!(mk(&[]).kmer_size, 25);
    assert_eq!(mk(&["--KMER_SIZE", "31"]).kmer_size, 31);
    assert_eq!(mk(&["__KMER_SIZE", "31"]).kmer_size, 31);
    assert_eq!(mk(&[]).min_kmer_count, 1);
    assert_eq!(mk(&["--min_kmer_count", "3"]).min_kmer_count, 3);
    assert_eq!(mk(&["--min_kmer_cov", "3"]).min_kmer_count, 3);
    // 主名优先于别名
    assert_eq!(
        mk(&["--min_kmer_count", "2", "--min_kmer_cov", "9"]).min_kmer_count,
        2
    );
    // 非法值报错
    let argv: Vec<String> = vec![
        "--seqType".into(),
        "fa".into(),
        "--single".into(),
        "s.fa".into(),
        "--max_memory".into(),
        "1G".into(),
        "--output".into(),
        "out".into(),
        "__KMER_SIZE".into(),
        "x".into(),
    ];
    assert!(parse_args(&argv).is_err());
}

/// 默认收尾清理: single.fa / jellyfish 产物删除, .ok 断点与最终产物保留。
#[test]
fn cleanup_removes_intermediates_by_default() {
    let dir = tmpdir("cleanup");
    let fa = write_synthetic_reads(&dir);
    let outdir = dir.join("trinity_out");
    let final_fa = run_trinity(&synth_args(&fa, &outdir)).unwrap();
    assert!(final_fa.exists());
    for gone in [
        "single.fa",
        "jellyfish.kmers.25.asm.fa",
        "mer_counts.25.asm.rs.bin",
    ] {
        assert!(!outdir.join(gone).exists(), "{gone} should be cleaned up");
    }
    // .ok 断点保留 → resume 语义可用
    assert!(outdir.join("single.fa.ok").exists());
    assert!(outdir.join(".jellyfish_dump.25.mincov1.asm.ok").exists());
}

/// resume 在清理后仍可用（产物缺失 + .ok 在 → prep 重建, 各阶段跳过）。
#[test]
fn resume_works_after_cleanup() {
    let dir = tmpdir("resume-clean");
    let fa = write_synthetic_reads(&dir);
    let outdir = dir.join("trinity_out");
    run_trinity(&synth_args(&fa, &outdir)).unwrap();
    assert!(!outdir.join("single.fa").exists());
    let final_fa = run_trinity(&synth_args(&fa, &outdir)).unwrap();
    assert!(final_fa.exists());
}
