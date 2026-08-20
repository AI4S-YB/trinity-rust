//! Pipeliner.pm 语义的 checkpoint 框架。
//!
//! 每个步骤: 执行前查 `checkpoint` 存在 → 存在即跳过（返回 `Ok(false)`）;
//! 执行闭包成功 → touch checkpoint（空文件）返回 `Ok(true)`。stderr 进度行
//! 采用原版 `-- CMD` 回显风格的自定变体（`-- ...` / `-p ...`）。

use std::fs;
use std::path::Path;

use trinity_common::error::CommonError;

pub fn checkpoint_exists(p: &Path) -> bool {
    p.exists()
}

/// Pipeliner.pm 等价: exists → skip（Ok(false)）; f 成功 → touch → Ok(true)。
/// f 失败 → Err 原样上抛（不写 checkpoint——原版失败时不 touch）。
pub fn run_with_checkpoint<F>(ckpt: &Path, msg: &str, f: F) -> Result<bool, CommonError>
where
    F: FnOnce() -> Result<(), CommonError>,
{
    if checkpoint_exists(ckpt) {
        eprintln!(
            "---- Trinity (rust) checkpoint found, skipping: {}",
            ckpt.display()
        );
        return Ok(false);
    }
    eprintln!("-- Trinity (rust) | {msg}");
    f()?;
    fs::write(ckpt, b"")?;
    Ok(true)
}
