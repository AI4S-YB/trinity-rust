//! 阶段计时探针（M7）: `timed` 包裹任意阶段闭包, 耗时追加到
//! `<outdir>/.stage_timings`（`阶段名\t秒`, OpenOptions append）。

use std::fs::OpenOptions;
use std::io::Write;
use std::path::Path;

/// 运行 `f` 并把 `(name, 耗时秒)` 追加到 `<outdir>/.stage_timings`。
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
}
