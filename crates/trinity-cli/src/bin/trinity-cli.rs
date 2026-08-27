//! trinity-cli 主程序入口: 参数解析（原版同名参数面）→ 编排主线。
//! 非零退出码 = CommonError/CliError 文案直写 stderr（原版 confess 风格）。

// 全局 mimalloc（分配器 A/B 实测见 docs/opt-log 或提交记录; 小集上 minor
// page fault 与 sys time 显著低于 glibc malloc）。
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() {
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
