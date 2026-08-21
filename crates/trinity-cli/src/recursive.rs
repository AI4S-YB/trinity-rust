//! 递归迷你管线（M1-T2/T3）——原版 `--trinity_complete` 递归架构（butterfly_cmds
//! == "RECURSIVE_TRINITY_COMPLETE", Trinity L1903）的进程内移植:
//! phase-1 止于分区（`skip_quantify`）, 逐分区在**内存中**重跑一条迷你
//! Trinity（inchworm-from-reads → min100 → chrysalis 全链 → butterfly）,
//! 分区产物落 `read_partitions/c<id>.mini.fa`（存在即跳过 = 断点）。

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use rayon::prelude::*;
use trinity_butterfly::{run_component, ComponentParams};
use trinity_chrysalis::pipeline::{run_chrysalis_pipeline, ChrysalisParams};
use trinity_common::error::CommonError;
use trinity_inchworm::irke::{
    compute_sequence_assemblies, populate_from_reads, prune_some_kmers, AssemblyParams, IrkeParams,
    Monitor,
};

use crate::args::TrinityArgs;
use crate::orchestrate::filter_iworm_min_len_or_cov;

/// 单分区迷你 Trinity 参数。
#[derive(Debug, Clone)]
pub struct MiniParams {
    /// 外层 K（原版迷你用外层 K, 默认 25）。
    pub k: usize,
    /// --SS。
    pub strand: bool,
    /// -L。
    pub min_contig_length: usize,
    /// --path_reinforcement_distance（原版主线恒 25）。
    pub path_reinforcement_distance: usize,
    /// -F（--group_pairs_distance）。
    pub max_pair_distance: usize,
    /// --minKmerCount（phase-2 仅 >1 时生效; ≤1 时等价 inchworm 默认 1）。
    pub min_kmer_count: u32,
    /// butterfly 专用线程栈（--bfly_stack_mb）。
    pub stack_size_mb: usize,
    /// 组件名前缀（分区 id, 保证跨分区转录本名唯一）。
    pub name: String,
}

/// 单分区运行结果。
pub struct MiniOut {
    /// butterfly allProbPaths fasta 全文（可含多条记录）。
    pub fasta: Vec<u8>,
    pub n_partitions_run: usize,
}

/// 单分区迷你 Trinity（TRINITY_COMPLETE 语义）:
/// `reads_tmp` 字节（">name/1 pct\nseq" 格式, /1 /2 保留——PE 配对靠它传到
/// butterfly）→ inchworm(内置 kmer 目录, --reads 装载) → min100 →
/// chrysalis 全链 → butterfly 每组件 → fasta 片段。
///
/// 剪枝开关（已核实原版 Trinity v2.15.2, TRINITY_COMPLETE 分支）:
/// - `--min_any_entropy 1.0` **只加在 phase-1**（L1050 的 else 分支）——
///   phase-2（TRINITY_COMPLETE）不加 → min_entropy = 0.0（不过滤）;
/// - `--no_prune_error_kmers` 仅在 `!TRINITY_COMPLETE && !super_read_mode`
///   时加（L2705-2708）→ **phase-2 保留错误 kmer 剪枝**（prune = true）;
/// - `--minKmerCount` 在 phase-2 的 `--reads` 装载路径仅当 min_kmer_cov > 1
///   时追加（L2670-2672）; ≤1 时 inchworm 自身默认即 1, 等价。
pub fn run_mini(reads_tmp: &[u8], p: &MiniParams) -> Result<MiniOut, CommonError> {
    // 分区池并行调用 → 临时目录须逐调用唯一（pid + 单调序号）。
    static SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let n = SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let tmp = std::env::temp_dir().join(format!("mini_{}_{n}", std::process::id()));
    let _ = fs::remove_dir_all(&tmp);
    fs::create_dir_all(&tmp)?;
    // catch_unwind: panic（深递归/数据异常）时也能落到下方 remove_dir_all
    // 清理临时目录（栈修复后分区 job 不再 abort, 该路径真正生效）。
    let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(
        || -> Result<MiniOut, CommonError> {
            // ① inchworm from reads（phase-2: --reads 装载, 内置 kmer 目录）
            let (mut counter, _parsed) = populate_from_reads(reads_tmp, p.k, !p.strand)?;
            // phase-2 剪枝实参（见函数文档）: min_count=max(min_kmer_count,1),
            // min_entropy=0.0, prune_error_kmers=true, min_ratio_non_error=0.005
            prune_some_kmers(&mut counter, p.min_kmer_count.max(1), 0.0, true, 0.005);
            // ② 单线程排序种子组装（NO_PARALLEL_IWORM——phase-2 强制, inchworm_cpu=1）
            let mut buf = Vec::new();
            compute_sequence_assemblies(
                &mut counter,
                &IrkeParams::default(),
                &AssemblyParams {
                    min_assembly_length: p.k,
                    ..Default::default()
                },
                &Monitor::new(0),
                true, // sort_seeds
                &mut buf,
            )?;
            // ③ min100（与主线 phase-1 相同: filter_iworm_by_min_length_or_cov.pl 100 10）
            let iworm = filter_iworm_min_len_or_cov(&buf, 100, 10)?;
            if iworm.is_empty() {
                return Ok(MiniOut {
                    fasta: Vec::new(),
                    n_partitions_run: 1,
                });
            }
            // ④ chrysalis 全链（phase-2 pct→50 钳制 = 默认; threads 1）
            let listing = run_chrysalis_pipeline(
                &iworm,
                reads_tmp,
                &tmp.join("chry"),
                &ChrysalisParams {
                    strand: p.strand,
                    min_contig_length: p.min_contig_length,
                    pct_required: 50,
                    threads: 1,
                    ..Default::default()
                },
            )?;
            // ⑤ butterfly 每组件（run_component 自带大栈线程）; 组件名带分区前缀
            //    → header `{name}_c{comp}_g.._i..`, 跨分区唯一且匹配 gene_trans_map 正则
            let mut fasta = Vec::new();
            for (_id, base) in &listing {
                let g = fs::read_to_string(base.with_extension("graph.out"))
                    .map_err(|e| CommonError::Parse(format!("mini graph.out read: {e}")))?;
                let r = fs::read_to_string(base.with_extension("graph.reads"))
                    .map_err(|e| CommonError::Parse(format!("mini graph.reads read: {e}")))?;
                let comp_name = format!(
                    "{}_{}",
                    p.name,
                    base.file_name().unwrap_or_default().to_string_lossy()
                );
                let out = run_component(
                    &g,
                    &r,
                    &ComponentParams {
                        n: 100_000,
                        min_len: p.min_contig_length,
                        max_pair_distance: p.max_pair_distance,
                        no_em_reduce: true,
                        path_reinforcement_distance: Some(p.path_reinforcement_distance),
                        name: comp_name,
                        stack_size: p.stack_size_mb.saturating_mul(1024 * 1024),
                        ..Default::default()
                    },
                )?;
                fasta.extend_from_slice(out.all_prob_paths_fasta.as_bytes());
            }
            Ok(MiniOut {
                fasta,
                n_partitions_run: 1,
            })
        },
    ))
    .unwrap_or_else(|_| {
        Err(CommonError::Parse(
            "mini pipeline panic（深递归/数据异常）".to_string(),
        ))
    });
    let _ = fs::remove_dir_all(&tmp);
    r
}

/// TrinityArgs → MiniParams（phase-2 语义收敛: K/SS/-L/-F/prd25/min_kmer_count）。
pub fn mini_params_from_args(args: &TrinityArgs, part_id: u64) -> MiniParams {
    MiniParams {
        k: args.kmer_size,
        strand: args.ss_lib_type.is_some(),
        min_contig_length: args.min_contig_length,
        path_reinforcement_distance: 25,
        max_pair_distance: args.group_pairs_distance,
        min_kmer_count: args.min_kmer_count,
        stack_size_mb: args.bfly_stack_mb,
        name: format!("c{part_id}"),
    }
}

/// 分区池（M1-T3）: 逐分区 `run_mini` → `outdir/c<id>.mini.fa`。
///
/// 断点: mini.fa **存在即跳过**该分区（断点续跑只需重算缺失分区）。
/// 失败隔离（镜像 butterfly_pool）: 单分区 Err 不中断其余, 结束时失败清单
/// 非空 → 写 `failed_read_partitions.txt` 并整体 Err。
///
/// `chrysalis_dir` 保留在签名中以与 butterfly_pool 的调用形状对齐
/// （listing 的 base 路径已含 Component_bins 前缀, 实际不单独使用）。
///
/// 返回**全部分区**的 mini.fa 路径（comp id 升序, 含本次跳过的既有文件）。
pub fn run_partition_pool(
    listing: &[(u64, PathBuf)],
    _chrysalis_dir: &Path,
    args: &TrinityArgs,
    outdir: &Path,
) -> Result<Vec<PathBuf>, CommonError> {
    fs::create_dir_all(outdir)?;
    let failed: Mutex<Vec<(u64, String)>> = Mutex::new(Vec::new());
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.cpu.max(1))
        .build()
        .map_err(|e| CommonError::Parse(format!("cannot build partition pool: {e}")))?;

    // (id, mini.fa 路径, reads.tmp 路径)——reads.tmp = base.with_extension
    // （照抄 pipeline.rs quantify 段的既有用法）
    let jobs: Vec<(u64, PathBuf, PathBuf)> = listing
        .iter()
        .map(|(id, base)| {
            (
                *id,
                outdir.join(format!("c{id}.mini.fa")),
                base.with_extension("reads.tmp"),
            )
        })
        .collect();

    // 大分区优先（reads.tmp 字节数降序）: 最重/最易失败的任务先暴露问题,
    // 也让长尾分区填充尾部空闲槽位。
    let mut jobs = jobs;
    jobs.sort_by_key(|(_, _, reads_tmp)| {
        std::cmp::Reverse(fs::metadata(reads_tmp).map(|m| m.len()).unwrap_or(0))
    });

    pool.install(|| {
        jobs.par_iter().for_each(|(id, mini_fa, reads_tmp)| {
            if mini_fa.exists() {
                return; // 断点: 已完成分区
            }
            let data = match fs::read(reads_tmp) {
                Ok(d) => d,
                Err(e) => {
                    eprintln!("-read partition {id} failed: cannot read {e}");
                    failed
                        .lock()
                        .unwrap()
                        .push((*id, format!("{}: {e}", reads_tmp.display())));
                    return;
                }
            };
            // 大栈线程（对齐 run_component）: rayon worker 默认 ~2MB 栈,
            // inchworm/chrysalis 的排序与种子扩展递归在真数据分区上会溢出。
            // rayon 池保持 args.cpu 路由; 每分区 job 由专用线程执行 run_mini。
            let params = mini_params_from_args(args, *id);
            let handle = std::thread::Builder::new()
                .stack_size(256 * 1024 * 1024)
                .spawn(move || run_mini(&data, &params));
            let res = match handle {
                Ok(h) => match h.join() {
                    Ok(r) => r,
                    Err(_) => Err(CommonError::Parse(format!(
                        "c{id}: mini job thread panic（深递归/数据异常）"
                    ))),
                },
                Err(e) => Err(CommonError::Parse(format!(
                    "c{id}: cannot spawn mini job thread: {e}"
                ))),
            };
            match res {
                Ok(out) => {
                    if let Err(e) = fs::write(mini_fa, &out.fasta) {
                        eprintln!("-read partition {id} failed: {e}");
                        failed
                            .lock()
                            .unwrap()
                            .push((*id, format!("{}: {e}", mini_fa.display())));
                    }
                }
                Err(e) => {
                    eprintln!("-read partition {id} failed: {e}");
                    failed.lock().unwrap().push((*id, format!("c{id}: {e}")));
                }
            }
        });
    });

    let mut failed = failed.into_inner().unwrap();
    if !failed.is_empty() {
        failed.sort();
        let path = outdir.join("failed_read_partitions.txt");
        let text: String = failed.iter().map(|(_, s)| format!("{s}\n")).collect();
        let _ = fs::write(&path, &text);
        return Err(CommonError::Parse(format!(
            "Error, mini pipeline failed for {} partition(s) (see {})",
            failed.len(),
            path.display()
        )));
    }
    Ok(jobs.into_iter().map(|(_, mini_fa, _)| mini_fa).collect())
}

/// 递归模式汇总: mini.fa 拼接 → min_len 过滤 + 全序列去重
/// （print_butterfly_assemblies.pl 语义, 与 harvest 相同; gene_trans_map 由
/// orchestrate 复用正则落盘）。返回收录转录本数。
pub fn harvest_mini(
    mini_files: &[PathBuf],
    min_len: usize,
    final_path: &Path,
) -> Result<usize, CommonError> {
    use std::collections::HashSet;
    let mut seen: HashSet<String> = HashSet::new();
    let mut out = String::new();
    let mut n = 0usize;
    for p in mini_files {
        let text = fs::read_to_string(p).map_err(|e| {
            CommonError::Parse(format!(
                "Error, no fasta file reported as: {} ({e})",
                p.display()
            ))
        })?;
        for (header, seq) in parse_records(&text) {
            if seq.len() < min_len {
                continue;
            }
            if !seen.insert(seq.clone()) {
                eprintln!(
                    "-duplicate sequence detected, excluding it. ({})",
                    p.display()
                );
                continue;
            }
            out.push('>');
            out.push_str(&header);
            out.push('\n');
            for chunk in seq.as_bytes().chunks(60) {
                out.push_str(std::str::from_utf8(chunk).unwrap_or(""));
                out.push('\n');
            }
            n += 1;
        }
    }
    fs::write(final_path, out)?;
    Ok(n)
}

/// 简单 FASTA 解析（harvest::records 同款, 递归路径内联保底）。
fn parse_records(text: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut header: Option<String> = None;
    let mut seq = String::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if let Some(h0) = header.take() {
                out.push((h0, std::mem::take(&mut seq)));
            }
            header = Some(h.to_string());
        } else if header.is_some() {
            seq.push_str(line.trim());
        }
    }
    if let Some(h0) = header {
        out.push((h0, seq));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

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

    /// 合成分区 reads（reads.tmp 格式）: 600nt 参考（xorshift）× 25 对平铺
    /// PE reads（100bp, 步进 14）; R1 = ref[s..s+100],
    /// R2 = revcomp(ref[s+150..s+250])（模拟 FR 文库插入 250）。
    fn synth_partition_reads() -> Vec<u8> {
        let mut x: u64 = 0x243F6A8885A308D3;
        let reference: String = (0..600).map(|_| xorshift_nt(&mut x)).collect();
        let rb = reference.as_bytes();
        let mut text = String::new();
        for i in 0..25usize {
            let s = i * 14;
            text.push_str(&format!(">r{i}/1 100%\n{}\n", &reference[s..s + 100]));
            text.push_str(&format!(
                ">r{i}/2 100%\n{}\n",
                revcomp(std::str::from_utf8(&rb[s + 150..s + 250]).unwrap())
            ));
        }
        text.into_bytes()
    }

    #[test]
    fn mini_pipeline_recovers_reference() {
        let reads = synth_partition_reads();
        let out = run_mini(
            &reads,
            &MiniParams {
                k: 25,
                strand: false,
                min_contig_length: 100,
                path_reinforcement_distance: 25,
                max_pair_distance: 500,
                min_kmer_count: 1,
                stack_size_mb: 64,
                name: "c0".to_string(),
            },
        )
        .unwrap();
        let s = String::from_utf8(out.fasta.clone()).unwrap();
        assert!(
            s.contains('>') && s.len() > 500,
            "应产出至少一条 ≥500bp 转录本, got {} bytes",
            s.len()
        );
        // 转录本应与参考 600nt 高度重叠: 简化断言——最长记录（折行合并）≥ 400
        let longest = parse_records(&s)
            .iter()
            .map(|(_, seq)| seq.len())
            .max()
            .unwrap_or(0);
        assert!(longest >= 400, "最长转录本 {longest}bp < 400");
        // 组件名前缀进入 header（gene_trans_map 正则 `c\d+_g\d+_i\d+` 可解析）
        assert!(
            s.contains(">c0_c"),
            "header 应带分区前缀 c0_c*: {}",
            &s[..s.len().min(120)]
        );
    }

    /// 真实大分区单跑计时（dev1sub c1428, ~9.4MB reads.tmp）:
    /// `cargo test -p trinity-cli --release -- --ignored real_partition_timing --nocapture`
    /// 分阶段计时复刻 run_mini 主体（populate → assembly → chrysalis → butterfly）。
    #[test]
    #[ignore]
    fn real_partition_timing() {
        let path = std::path::PathBuf::from(
            "/tmp/rec_smoke/trinity_out/chrysalis/Component_bins/Cbin1/c1428.reads.tmp",
        );
        let data = match fs::read(&path) {
            Ok(d) => d,
            Err(e) => {
                eprintln!("skip: cannot read {}: {e}", path.display());
                return;
            }
        };
        eprintln!("reads.tmp: {} bytes", data.len());
        // 全程在 256MB 栈线程内（与 run_partition_pool job 一致）——测试线程
        // 默认 2MB 栈会复现被修复的栈溢出。
        let data0 = data.clone();
        let staged = std::thread::Builder::new()
            .stack_size(256 * 1024 * 1024)
            .spawn(move || {
                let t0 = std::time::Instant::now();
                let (mut counter, _) =
                    populate_from_reads(&data0, 25, true).expect("populate_from_reads");
                prune_some_kmers(&mut counter, 1, 0.0, true, 0.005);
                let t1 = std::time::Instant::now();
                let mut buf = Vec::new();
                compute_sequence_assemblies(
                    &mut counter,
                    &IrkeParams::default(),
                    &AssemblyParams {
                        min_assembly_length: 25,
                        ..Default::default()
                    },
                    &Monitor::new(0),
                    true,
                    &mut buf,
                )
                .expect("inchworm");
                let t2 = std::time::Instant::now();
                let iworm = filter_iworm_min_len_or_cov(&buf, 100, 10).expect("min100");
                eprintln!(
                    "inchworm: raw {} contigs ({} bytes) → min100 {} bytes",
                    buf.iter().filter(|&&b| b == b'>').count(),
                    buf.len(),
                    iworm.len()
                );
                (t0, t1, t2)
            })
            .expect("spawn staged")
            .join()
            .expect("join staged");
        let (t0, t1, t2) = staged;
        // 与 run_partition_pool 的 job 线程一致（256MB 栈）——测试线程默认 2MB
        // 会复现被修复的栈溢出。
        let data2 = data.clone();
        let out = std::thread::Builder::new()
            .stack_size(256 * 1024 * 1024)
            .spawn(move || {
                run_mini(
                    &data2,
                    &MiniParams {
                        k: 25,
                        strand: false,
                        min_contig_length: 100,
                        path_reinforcement_distance: 25,
                        max_pair_distance: 500,
                        min_kmer_count: 1,
                        stack_size_mb: 256,
                        name: "c1428".to_string(),
                    },
                )
            })
            .expect("spawn")
            .join()
            .expect("join")
            .expect("run_mini");
        let t3 = std::time::Instant::now();
        eprintln!(
            "TIMING populate+prune={:?} inchworm_asm={:?} run_mini_total={:?} (chrysalis+butterfly={:?}) fasta={} bytes n_records={}",
            t1 - t0,
            t2 - t1,
            t3 - t0,
            (t3 - t0) - (t2 - t0),
            out.fasta.len(),
            out.fasta.iter().filter(|&&b| b == b'>').count()
        );
    }

    #[test]
    fn mini_empty_reads_produce_empty_fasta() {
        let out = run_mini(
            b">r/1 100%\nACGT\n",
            &MiniParams {
                k: 25,
                strand: false,
                min_contig_length: 100,
                path_reinforcement_distance: 25,
                max_pair_distance: 500,
                min_kmer_count: 1,
                stack_size_mb: 64,
                name: "c0".to_string(),
            },
        )
        .unwrap();
        assert!(out.fasta.is_empty());
        assert_eq!(out.n_partitions_run, 1);
    }

    #[test]
    fn harvest_mini_filters_and_dedups() {
        let dir = std::env::temp_dir().join(format!("mini_harvest_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let seq = "ACGT".repeat(60); // 240bp
        fs::write(
            dir.join("c0.mini.fa"),
            format!(
                ">c0_c0_g1_i1 len=240 path=[1] [1]\n{seq}\n>c0_c0_g1_i2 len=240 path=[2] [2]\n{seq}\n"
            ),
        )
        .unwrap();
        let seq2 = format!("{}A", "ACGT".repeat(60)); // 241bp, 与 seq 不同（跨分区不去重）
        fs::write(
            dir.join("c1.mini.fa"),
            format!(
                ">c1_c0_g1_i1 len=241 path=[3] [3]\n{seq2}\n>c1_c0_g1_i2 len=5 path=[4] [4]\nACGTA\n"
            ),
        )
        .unwrap();
        let final_path = dir.join("Trinity.fasta");
        let n = harvest_mini(
            &[dir.join("c0.mini.fa"), dir.join("c1.mini.fa")],
            100,
            &final_path,
        )
        .unwrap();
        // c0_i1 收录; c0_i2 重复序列排除; c1_i1 收录（序列不同）; c1_i2 过短排除
        assert_eq!(n, 2);
        let text = fs::read_to_string(&final_path).unwrap();
        assert!(text.contains(">c0_c0_g1_i1"));
        assert!(text.contains(">c1_c0_g1_i1"));
        assert!(!text.contains("_i2"));
    }
}
