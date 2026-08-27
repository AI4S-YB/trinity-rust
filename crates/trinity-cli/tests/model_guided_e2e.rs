//! M5 e2e: `--model-guided`（外部模型）合成集——3 转录本（T0/T1 共享
//! 500nt 前缀 → 同基因, T2 独立）→ 2 基因组, 150 条 read 全部入组
//! （无 orphan 组产物）, 产出转录本 + gene_trans_map + 断点;
//! 关闭旗标时 classic 路径无 model-guided 产物。

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use trinity_cli::args::parse_args;
use trinity_cli::orchestrate::run_trinity;

fn tmpdir(tag: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!(
        "trinity-cli-modelguided-{}-{}",
        tag,
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&p);
    fs::create_dir_all(&p).unwrap();
    p
}

fn xorshift_nt(x: &mut u64) -> char {
    *x ^= *x << 13;
    *x ^= *x >> 7;
    *x ^= *x << 17;
    match *x % 4 {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        _ => 'T',
    }
}

fn randseq(seed: u64, n: usize) -> String {
    let mut x = seed | 1;
    (0..n).map(|_| xorshift_nt(&mut x)).collect()
}

fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect()
}

/// 与 model_guided.rs 单测同款模型: T0/T1 共享 500nt 前缀（同基因）,
/// T2 独立 600nt。
fn synth_model() -> (String, String, String) {
    let common = randseq(0x243F6A8885A308D3, 500);
    let t0 = format!("{}{}", common, randseq(11, 100));
    let t1 = format!("{}{}", common, randseq(22, 100));
    let t2 = randseq(0x9E3779B97F4A7C15, 600);
    (t0, t1, t2)
}

/// 每条转录本平铺 n_pairs 对 PE reads（100nt, FR 插入 250）→ fq.gz 对。
fn write_tiled_pe(name: &str, seq: &str, n_pairs: usize, r1: &mut String, r2: &mut String) {
    let b = seq.as_bytes();
    let step = (b.len() - 250) / n_pairs.max(1);
    for i in 0..n_pairs {
        let s = i * step;
        let f1 = std::str::from_utf8(&b[s..s + 100]).unwrap();
        let f2 = std::str::from_utf8(&b[s + 150..s + 250]).unwrap();
        r1.push_str(&format!("@{name}{i}\n{f1}\n+\n{}\n", "I".repeat(100)));
        r2.push_str(&format!(
            "@{name}{i}\n{}\n+\n{}\n",
            revcomp(f2),
            "I".repeat(100)
        ));
    }
}

/// 模型 fasta + 25 对/转录本 PE fq.gz（3 × 25 × 2 = 150 条 read）。
fn write_synthetic(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let (t0, t1, t2) = synth_model();
    let model = dir.join("model.fa");
    fs::write(&model, format!(">T0\n{t0}\n>T1\n{t1}\n>T2\n{t2}\n")).unwrap();

    let mut r1 = String::new();
    let mut r2 = String::new();
    write_tiled_pe("a", &t0, 25, &mut r1, &mut r2);
    write_tiled_pe("b", &t1, 25, &mut r1, &mut r2);
    write_tiled_pe("c", &t2, 25, &mut r1, &mut r2);
    let p1 = dir.join("r1.fq.gz");
    let p2 = dir.join("r2.fq.gz");
    for (p, body) in [(&p1, &r1), (&p2, &r2)] {
        let mut w = flate2::write::GzEncoder::new(
            fs::File::create(p).unwrap(),
            flate2::Compression::fast(),
        );
        w.write_all(body.as_bytes()).unwrap();
        w.finish().unwrap();
    }
    (p1, p2, model)
}

fn base_args(r1: &Path, r2: &Path, outdir: &Path) -> Vec<String> {
    [
        "--seqType",
        "fq",
        "--left",
        r1.to_str().unwrap(),
        "--right",
        r2.to_str().unwrap(),
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
    .collect()
}

fn parse_records(text: &str) -> Vec<(String, usize)> {
    let mut out = Vec::new();
    let mut header: Option<String> = None;
    let mut len = 0usize;
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if let Some(h0) = header.take() {
                out.push((h0, len));
            }
            header = Some(h.split_whitespace().next().unwrap_or_default().to_string());
            len = 0;
        } else {
            len += line.trim().len();
        }
    }
    if let Some(h0) = header {
        out.push((h0, len));
    }
    out
}

/// minimap2 缺失时跳过（模块单测同样硬依赖, 此处仅免于环境噪音）。
fn minimap2_available() -> bool {
    let bin = std::env::var("TRINITY_MINIMAP2")
        .ok()
        .filter(|v| !v.is_empty())
        .unwrap_or_else(|| {
            "/public/home/senior007/miniconda3/envs/cactus/bin/minimap2".to_string()
        });
    std::process::Command::new(bin)
        .arg("--version")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

#[test]
fn model_guided_e2e_external_model_and_flag_off() {
    if !minimap2_available() {
        eprintln!("skip: minimap2 not available");
        return;
    }
    let dir = tmpdir("e2e");
    let (r1, r2, model) = write_synthetic(&dir);

    // ---- on: 外部模型 → 2 基因组, 全部 reads 入组, 产物齐备 ----
    let out_on = dir.join("trinity_on");
    let mut argv = base_args(&r1, &r2, &out_on);
    argv.push("--model-guided".to_string());
    argv.push("--model".to_string());
    argv.push(model.display().to_string());
    let final_on = run_trinity(&parse_args(&argv).unwrap()).unwrap();

    let groups_dir = out_on.join("model_guided_groups");
    assert!(out_on.join(".model_guided.ok").exists(), "断点缺失");
    assert!(
        groups_dir.join("g0.mini.fa").exists() && groups_dir.join("g1.mini.fa").exists(),
        "应有 2 个非空基因组的 mini 产物"
    );
    assert!(
        !groups_dir.join("orphan.mini.fa").exists(),
        "reads 全部入组, 不应有 orphan 组"
    );
    assert!(
        !out_on.join("first_pass").exists(),
        "外部模型模式不应跑自引导 first pass"
    );

    let records = {
        let text = fs::read_to_string(&final_on).unwrap();
        assert!(!text.is_empty() && text.contains('>'), "产出为空");
        parse_records(&text)
    };
    assert!(records.len() >= 2, "两组应各产出 ≥1 条: {records:?}");
    let longest = records.iter().map(|(_, l)| *l).max().unwrap_or(0);
    assert!(longest >= 300, "最长转录本 {longest}bp < 300");

    let map_path = PathBuf::from(format!("{}.gene_trans_map", final_on.display()));
    let map_text = fs::read_to_string(&map_path).unwrap();
    let map_lines: Vec<&str> = map_text.lines().filter(|l| !l.is_empty()).collect();
    assert_eq!(
        map_lines.len(),
        records.len(),
        "gene_trans_map 行数应与转录本数一致"
    );
    for line in &map_lines {
        let mut f = line.split('\t');
        let (Some(gene), Some(trans), None) = (f.next(), f.next(), f.next()) else {
            panic!("gene_trans_map 行格式错误: {line}");
        };
        assert!(
            trans.starts_with(gene),
            "transcript 应以 gene 为前缀: {trans} !~ {gene}"
        );
    }

    // ---- off: classic 路径无 model-guided 产物 ----
    let out_off = dir.join("trinity_off");
    let final_off = run_trinity(&parse_args(&base_args(&r1, &r2, &out_off)).unwrap()).unwrap();
    assert!(fs::read_to_string(&final_off).unwrap().contains('>'));
    assert!(!out_off.join(".model_guided.ok").exists());
    assert!(!out_off.join("model_guided_groups").exists());
    assert!(!out_off.join("first_pass").exists());
}
