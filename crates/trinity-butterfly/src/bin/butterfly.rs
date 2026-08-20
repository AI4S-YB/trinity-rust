//! butterfly CLI —— TransAssembly_allProbPaths 驱动的薄壳（镜像原版 getopt 面）。
//!
//! 原版：`java -jar Butterfly.jar -N <int> -L <int> -F <int> -C <prefix>
//! [-V <int>] [-R <int>] [-O <int>] [--NO_EM_REDUCE]
//! [--path_reinforcement_distance=<int>] [--no_path_merging]`。
//!
//! 全链编排上移到库层 [`trinity_butterfly::run_component`]（CLI / 集成测试 /
//! xtask 对拍同入口）；本 bin 只做：解析参数 → 读文件 → run_component → 写文件。
//!
//! 文件拼装（Java main L726-912）：
//! * `-C <prefix>` → 读 `<prefix>.out`（图）与 `<prefix>.reads`（read 路径）
//! * 输出 `<prefix>.allProbPaths.fasta`（`--stderr` 改打 stderr；
//!   `--log_stderr` 另写 `<prefix>.err`，本实现只写空占位——Java 的 err
//!   是进度日志，无对拍价值）
//! * 组件名 = `-C` 值的最后一个 `/` 段（L983 `file.split("/")`，
//!   printFinalPaths 再去掉 `.graph` 后缀）
//!
//! `--NO_EM_REDUCE` 的 CLI 默认对齐 Java getopt 默认 = **false**（镜像原则：
//! Trinity 主线脚本显式传该旗标；xcheck 对拍时显式区分两种形态）。

use std::io::Write;

use trinity_butterfly::{run_component, ComponentParams};
use trinity_common::cli::{Cli, CliError};

const USAGE: &str = "用法: butterfly -N <int> -L <int> -F <int> -C <prefix> [选项]
必选:
  -N  <int>     total number of reads or fragment pairs（仅校验非零）
  -L  <int>     min length for an assembled sequence to be reported
  -F  <int>     maximum fragment length (extreme dist between paired ends)
  -C  <string>  prefix for component/reads file（读 <prefix>.out 与 <prefix>.reads，
                写 <prefix>.allProbPaths.fasta）
可选:
  -R  <int>     minimum read support threshold（默认 2）
  -O  <int>     path reinforcement percent of -F（默认 25）
  -V  <int>     verbosity（默认 10；本实现接受但不输出逐级日志）
  --path_reinforcement_distance=<int>  直接指定 reinforcement 距离（覆盖 -O 推算）
  --NO_EM_REDUCE      跳过 EM 削减（Trinity 主线形态；默认 false = 跑 EM，镜像 Java）
  --no_path_merging   禁用 cd-hit 式路径合并
  --stderr            输出改打 stderr
  --log_stderr        写 <prefix>.err（占位）";

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    std::process::exit(run(&args));
}

fn run(args: &[String]) -> i32 {
    let mut cli = Cli::new(args);
    match run_inner(&mut cli) {
        Ok(()) => 0,
        Err(e) => {
            if e.usage {
                eprintln!("{USAGE}\nError, {}", e.msg);
                2
            } else {
                eprintln!("Error, {}", e.msg);
                1
            }
        }
    }
}

fn run_inner(cli: &mut Cli) -> Result<(), CliError> {
    // ---- 必选参数（Java printUsage 条件：N==0 || F==0 || C 空 || R<1）----
    let total_reads: i64 = cli.typed_flag("-N", 0)?;
    let min_output_seq = cli.uint_flag("-L", 0usize)?;
    let max_pair_distance = cli.uint_flag("-F", 0usize)?;
    let prefix = cli.req_flag("-C")?;
    let min_read_support = cli.int_flag("-R", 2i64)?;
    let verbose = cli.int_flag("-V", 10i32)?;
    let reinforcement_pct = cli.int_flag("-O", 25i64)?;
    let prd_flag = cli.value_flag("--path_reinforcement_distance")?;
    let no_em_reduce = cli.bool_flag("--NO_EM_REDUCE");
    let no_path_merging = cli.bool_flag("--no_path_merging");
    let use_stderr = cli.bool_flag("--stderr");
    let log_stderr = cli.bool_flag("--log_stderr");
    cli.finish()?;

    if total_reads == 0 || max_pair_distance == 0 || prefix.is_empty() || min_read_support < 1 {
        return Err(CliError::usage("缺少必选参数（-N/-L/-F/-C）"));
    }

    // Java L700+：显式 --path_reinforcement_distance 优先，否则 O% * F
    let path_reinforcement_distance = match prd_flag {
        Some(v) => v
            .parse::<usize>()
            .map_err(|e| CliError::usage(format!("--path_reinforcement_distance={v}: {e}")))?,
        None => (reinforcement_pct.max(0) as usize * max_pair_distance) / 100,
    };

    // ---- 读输入（<prefix>.out / <prefix>.reads）----
    let graph_text = std::fs::read_to_string(format!("{prefix}.out"))
        .map_err(|e| CliError::run(format!("cannot read {prefix}.out: {e}")))?;
    let reads_text = std::fs::read_to_string(format!("{prefix}.reads"))
        .map_err(|e| CliError::run(format!("cannot read {prefix}.reads: {e}")))?;
    if log_stderr {
        std::fs::write(format!("{prefix}.err"), "").ok();
    }

    // 组件名 = -C 值最后一段（Java file.split("/")；printFinalPaths 去掉 .graph）
    let comp_name = prefix.rsplit('/').next().unwrap_or(&prefix).to_string();

    let result = run_component(
        &graph_text,
        &reads_text,
        &ComponentParams {
            n: total_reads.max(0) as u64,
            min_len: min_output_seq,
            max_pair_distance,
            no_em_reduce,
            no_path_merging,
            min_read_support,
            path_reinforcement_distance: Some(path_reinforcement_distance),
            name: comp_name.clone(),
            ..ComponentParams::default()
        },
    )?;

    if verbose >= 10 {
        eprintln!(
            "Reported {} final paths (component {comp_name})",
            result.num_paths
        );
    }
    if use_stderr {
        std::io::stderr()
            .write_all(result.all_prob_paths_fasta.as_bytes())
            .map_err(CliError::run)?;
    } else {
        std::fs::write(
            format!("{prefix}.allProbPaths.fasta"),
            &result.all_prob_paths_fasta,
        )
        .map_err(|e| CliError::run(format!("cannot write {prefix}.allProbPaths.fasta: {e}")))?;
    }
    Ok(())
}
