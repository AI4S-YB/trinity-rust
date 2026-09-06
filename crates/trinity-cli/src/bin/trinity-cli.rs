//! trinity-cli 主程序入口: 参数解析（原版同名参数面）→ 编排主线。
//! 非零退出码 = CommonError/CliError 文案直写 stderr（原版 confess 风格）。

// 全局 mimalloc（分配器 A/B 实测见 docs/opt-log 或提交记录; 小集上 minor
// page fault 与 sys time 显著低于 glibc malloc）。
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// vendored mimalloc v2（libmimalloc-sys 0.1.49）的 `mi_option_purge_delay`
/// 枚举位次（mimalloc.h options 枚举 0-based: eager_commit_delay=14,
/// **purge_delay=15**, use_numa_nodes=16）。升级 mimalloc/libmimalloc-sys
/// 时须对照头文件复核此值。
const MI_OPTION_PURGE_DELAY: std::ffi::c_int = 15;

fn main() {
    // mimalloc 滞留页治理（2026-09-05 刀2 定位）: 分片归一化（diginorm）
    // 的大表释放后, mimalloc 默认 PURGE_DELAY=10ms 滞留已释放页; 与后续
    // multi-k 轮计数（每线程局部表翻倍增长）叠加, 全量 QM4 下 RSS 在 45s
    // 内 35→455GB 被 OOM-kill（3/3 复现; 同二进制同输入的 glibc 构建不
    // 复现）。置 0 即时归还: 全程 34:57 通过, 峰值 RSS 85GB。程序内设置
    // 会覆盖 MIMALLOC_PURGE_DELAY 环境变量（选项置位即初始化）。
    unsafe {
        libmimalloc_sys::mi_option_set(MI_OPTION_PURGE_DELAY, 0);
    }
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let args = match trinity_cli::parse_args(&argv) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("{}", e.msg);
            std::process::exit(1);
        }
    };
    match trinity_cli::orchestrate::run_trinity(&args) {
        Ok(final_fa) => println!("Trinity assembly written to {}", final_fa.display()),
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    }
}
