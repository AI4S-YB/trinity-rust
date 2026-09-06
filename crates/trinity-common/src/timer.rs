//! 阶段计时探针（M7）: `timed` 包裹任意阶段闭包, 耗时追加到
//! `<outdir>/.stage_timings`（`阶段名\t秒`, OpenOptions append）。
//! 同边界另落 `<outdir>/.stage_memory`（`阶段名\tVmHWM\tVmRSS`, kB）——
//! VmHWM 进程级单调不降, 阶段归属看跳变点; PURGE_DELAY=0 下段间释放
//! 体现为 VmRSS 回落。

use std::fs::OpenOptions;
use std::io::Write;
use std::path::Path;

/// /proc/self/status 的指定字段值（kB）。非 Linux 或字段缺失 → None。
pub fn vm_field_kb(field: &str) -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    vm_kb_of(&status, field)
}

fn vm_kb_of(status: &str, field: &str) -> Option<u64> {
    let line = status.lines().find(|l| l.starts_with(field))?;
    line.split_whitespace().nth(1)?.parse().ok()
}

/// 运行 `f` 并把 `(name, 耗时秒)` 追加到 `<outdir>/.stage_timings`;
/// 同边界追加 `<outdir>/.stage_memory`（kB; /proc 不可读时该文件不写）。
pub fn timed<T>(outdir: &Path, name: &str, f: impl FnOnce() -> T) -> T {
    let start = std::time::Instant::now();
    let out = f();
    let secs = start.elapsed().as_secs_f64();
    if let Some(parent) = outdir.parent() {
        let _ = std::fs::create_dir_all(parent);
    }
    if let Ok(mut fh) = OpenOptions::new()
        .create(true)
        .append(true)
        .open(outdir.join(".stage_timings"))
    {
        let _ = writeln!(fh, "{name}\t{secs:.3}");
    }
    if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
        if let (Some(hwm), Some(rss)) = (
            vm_kb_of(&status, "VmHWM:"),
            vm_kb_of(&status, "VmRSS:"),
        ) {
            if let Ok(mut fh) = OpenOptions::new()
                .create(true)
                .append(true)
                .open(outdir.join(".stage_memory"))
            {
                let _ = writeln!(fh, "{name}\t{hwm}\t{rss}");
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn timer_writes_stage_line() {
        let dir = std::env::temp_dir().join(format!("stage_timer_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let out = timed(&dir, "demo", || 40 + 2);
        assert_eq!(out, 42);
        let text = std::fs::read_to_string(dir.join(".stage_timings")).unwrap();
        let line = text.lines().last().unwrap();
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f.len(), 2);
        assert_eq!(f[0], "demo");
        assert!(f[1].parse::<f64>().unwrap() >= 0.0);
    }

    /// Linux 下 .stage_memory 边界行: 3 列, VmHWM ≥ VmRSS > 0。
    #[test]
    fn timer_writes_stage_memory_line() {
        let dir = std::env::temp_dir().join(format!("stage_timer_mem_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        timed(&dir, "demo", || ());
        if let Ok(text) = std::fs::read_to_string(dir.join(".stage_memory")) {
            let f: Vec<&str> = text.lines().last().unwrap().split('\t').collect();
            assert_eq!(f.len(), 3);
            assert_eq!(f[0], "demo");
            let (hwm, rss): (u64, u64) = (f[1].parse().unwrap(), f[2].parse().unwrap());
            assert!(rss > 0 && hwm >= rss);
        }
    }
}
