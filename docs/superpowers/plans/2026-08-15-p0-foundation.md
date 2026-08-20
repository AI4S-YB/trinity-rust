# P0 基础工程 — Trinity Rust 复刻实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 搭建 Cargo workspace 脚手架（6 crates + xtask）、实现 `trinity-common`（2-bit k-mer 编码/revcomp/熵 + FASTA/FASTQ IO）、编译原版 Trinity 供交叉验证，并通过黄金向量测试证明 Rust 与原版 C++ 行为位级一致。

**Architecture:** 镜像式 Cargo workspace。本计划只覆盖 spec（`docs/superpowers/specs/2026-08-15-trinity-rust-design.md`）的 P0 阶段；P1（trinity-kmer）在 P0 验证门通过后另写计划。trinity-common 的每个函数直译原版 `Inchworm/src/sequenceUtil.cpp`，注释标注原版行号。

**Tech Stack:** Rust (edition 2021, stable)、thiserror、flate2（rust backend = miniz_oxide，纯 Rust）、原版 Trinity v2.15.2（g++ 编译，仅供交叉验证）。

---

## 环境约定

- 新项目根目录：`/storage/home/senior007/test/trinity_rust/trinity-rust`（已是 git 仓库）
- 原版 Trinity 源码（FULL tarball 解压，含子模块源码）：`/storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2`
- 原版路径可被环境变量 `TRINITY_SRC` 覆盖（默认取上述绝对路径）
- 下文所有 `cargo` 命令默认在新项目根目录执行；所有引用原版的命令用 `$TRINITY_SRC` 或上述绝对路径
- 黄金向量等价性核心：原版编码 **G=0, A=1, T=2, C=3**（互补=按位取反），`_base_to_int` 表同时接受小写 gatc（sequenceUtil.cpp:10-24）

## 文件结构总览（本计划创建）

```
trinity-rust/
├── Cargo.toml                     # workspace 定义
├── .cargo/config.toml             # cargo xtask 别名
├── crates/
│   ├── trinity-common/            # ★ 本计划实现
│   │   ├── Cargo.toml
│   │   ├── src/lib.rs
│   │   ├── src/error.rs           # CommonError
│   │   ├── src/kmer.rs            # 编码/解码/revcomp/熵
│   │   ├── src/fasta.rs           # FastaRecord/FastaReader
│   │   ├── src/fastq.rs           # FastqRecord/FastqReader
│   │   └── src/io_util.rs         # gzip 魔数嗅探
│   │   └── tests/kmer_golden.rs   # 原版 C++ 黄金向量测试
│   ├── trinity-kmer/              # 空 crate（P1 填充）
│   ├── trinity-inchworm/          # 空 crate（P2 填充）
│   ├── trinity-chrysalis/         # 空 crate（P3 填充）
│   ├── trinity-butterfly/         # 空 crate（P4 填充）
│   └── trinity-cli/               # 空 crate（P5 填充）
├── xtask/
│   ├── Cargo.toml
│   ├── src/main.rs                # gen-fixtures 命令
│   └── fixtures-src/dump_kmer_golden.cpp  # 链接原版源码的向量生成器
├── fixtures/
│   ├── kmer_golden_input.txt      # 黄金向量输入（checked in）
│   └── kmer_golden.tsv            # 生成的期望输出（checked in）
└── docs/setup.md                  # 原版构建记录
```

---

### Task 1: 编译原版 Trinity（交叉验证基础设施）

**Files:**
- Create: `docs/setup.md`

- [x] **Step 1: 检查工具链**

```bash
which cmake g++ java || true
cmake --version | head -1
g++ --version | head -1
java -version 2>&1 | head -1
```
Expected: 四者都存在。缺任何一样先安装（HPC 上用 module load / conda）。java 只在运行 Butterfly 时需要，编译不需要。

- [x] **Step 2: 构建 Inchworm 与 Chrysalis（不需要 plugins）**

```bash
cd /storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2
make inchworm_target chrysalis_target 2>&1 | tail -20
```
Expected: 无 error 退出。Inchworm 走 cmake 包装的 Makefile，Chrysalis 直接 make。若报错，读错误信息——常见问题是 cmake 版本过旧或缺 OpenMP 头（g++ 自带）。

- [x] **Step 3: 验证关键二进制存在**

```bash
cd /storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2
ls -la Inchworm/bin/inchworm Inchworm/bin/FastaToDeBruijn \
  Chrysalis/bin/GraphFromFasta Chrysalis/bin/BubbleUpClustering \
  Chrysalis/bin/CreateIwormFastaBundle Chrysalis/bin/ReadsToTranscripts \
  Chrysalis/bin/QuantifyGraph Butterfly/Butterfly.jar
```
Expected: 8 个文件全部存在。`Butterfly.jar` 是源码树预编译的，无需构建。

- [x] **Step 4: smoke 测试 inchworm**

```bash
cd /storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2
printf ">r1\nACGTACGTACGTACGTACGTACGTACGT\n>r1\nACGTACGTACGTACGTACGTACGTACGT\n" > /tmp/tiny.fa
Inchworm/bin/inchworm --reads /tmp/tiny.fa --run_inchworm -K 25 --monitor 1
```
Expected: stdout 输出一条 FASTA contig。注：read 为 28bp；必须写两份才能过默认 `min_seed_coverage=2`（单份会被 coverage 过滤拒绝，exit 0 但无输出）。实测输出：
```
>a1;2 total_counts: 8 Seed: 2 K: 25 length: 28
TACGTACGTACGTACGTACGTACGTACG
```
（周期 4 序列的 de Bruijn 图成环，contig 是 read 的循环移位，起点由被选种子决定——这是正常行为。）header 格式 `>a<N>;<avg_cov> total_counts: <tc> Seed: <n> K: <K> length: <len>` 供 P2 对齐。

- [x] **Step 5: 写 docs/setup.md 并提交**

`docs/setup.md` 内容：

```markdown
# 环境搭建

## 原版 Trinity（仅供交叉验证）

- 源码: /storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2 (v2.15.2 FULL tarball)
- 构建: `make inchworm_target chrysalis_target`（Butterfly.jar 预编译）
- 关键二进制:
  - Inchworm/bin/inchworm, Inchworm/bin/FastaToDeBruijn
  - Chrysalis/bin/{GraphFromFasta,BubbleUpClustering,CreateIwormFastaBundle,ReadsToTranscripts,QuantifyGraph}
  - Butterfly/Butterfly.jar
- 环境变量 `TRINITY_SRC` 可覆盖源码路径（默认上述绝对路径），xtask 使用
```

```bash
cd /storage/home/senior007/test/trinity_rust/trinity-rust
git add docs/setup.md && git commit -m "docs: 原版 Trinity 构建记录（交叉验证基础设施）"
```

---

### Task 2: Cargo workspace 脚手架

**Files:**
- Create: `Cargo.toml`、`.cargo/config.toml`
- Create: `crates/{trinity-common,trinity-kmer,trinity-inchworm,trinity-chrysalis,trinity-butterfly,trinity-cli}/{Cargo.toml,src/lib.rs}`
- Create: `xtask/Cargo.toml`、`xtask/src/main.rs`

- [x] **Step 1: 创建 workspace 根 Cargo.toml**

`Cargo.toml`：
```toml
[workspace]
resolver = "2"
members = [
    "crates/trinity-common",
    "crates/trinity-kmer",
    "crates/trinity-inchworm",
    "crates/trinity-chrysalis",
    "crates/trinity-butterfly",
    "crates/trinity-cli",
    "xtask",
]

[workspace.package]
version = "0.1.0"
edition = "2021"

[workspace.dependencies]
thiserror = "2"
flate2 = "1"
```

- [x] **Step 2: 创建 6 个 crate**

对每个 crate（以 trinity-common 为例，其余 5 个把 name 换成对应名字、dependencies 段为空）：

`crates/trinity-common/Cargo.toml`：
```toml
[package]
name = "trinity-common"
version.workspace = true
edition.workspace = true

[dependencies]
thiserror = { workspace = true }
flate2 = { workspace = true }
```

其余 5 个（trinity-kmer / trinity-inchworm / trinity-chrysalis / trinity-butterfly / trinity-cli）的 Cargo.toml：
```toml
[package]
name = "trinity-kmer"        # 换成对应名字
version.workspace = true
edition.workspace = true

[dependencies]
```

每个 crate 的 `src/lib.rs`（trinity-common 除外，它在 Task 3 重写）：
```rust
//! P0 占位 — 将在后续阶段填充（见 docs/superpowers/specs/2026-08-15-trinity-rust-design.md §5）。
```

trinity-common 的 `src/lib.rs` 初始为：
```rust
//! 共享基础库: 2bit k-mer 编码、FASTA/FASTQ IO — 直译原版 Inchworm/src/sequenceUtil.cpp 等。

#[cfg(test)]
mod tests {
    #[test]
    fn scaffold_smoke() {
        assert!(true);
    }
}
```

- [x] **Step 3: 创建 xtask**

`xtask/Cargo.toml`：
```toml
[package]
name = "xtask"
version.workspace = true
edition.workspace = true

[dependencies]
```

`xtask/src/main.rs`（本任务先放骨架，Task 10 填充 gen-fixtures）：
```rust
use std::env;

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        None => {
            eprintln!("用法: cargo xtask <gen-fixtures>");
            std::process::exit(2);
        }
        Some(other) => {
            eprintln!("未知任务: {other}\n用法: cargo xtask <gen-fixtures>");
            std::process::exit(2);
        }
    }
}
```

- [x] **Step 4: cargo xtask 别名**

`.cargo/config.toml`：
```toml
[alias]
xtask = "run -p xtask --"
```

- [x] **Step 5: 验证构建与测试全绿**

```bash
cargo build --workspace 2>&1 | tail -3
cargo test --workspace 2>&1 | tail -5
```
Expected: `Finished`；测试 `0 failed`（只有 trinity-common 的 scaffold_smoke）。

- [x] **Step 6: Commit**

```bash
git add Cargo.toml .cargo crates xtask Cargo.lock
git commit -m "chore: Cargo workspace 脚手架（6 crates + xtask）"
```

---

### Task 3: kmer 编码与解码

**Files:**
- Modify: `crates/trinity-common/src/lib.rs`（改为模块声明）
- Create: `crates/trinity-common/src/error.rs`
- Create: `crates/trinity-common/src/kmer.rs`

- [x] **Step 1: 创建 error.rs**

`crates/trinity-common/src/error.rs`：
```rust
//! trinity-common 统一错误类型。消息格式贴近原版（"error, ..." 前缀风格）。

use thiserror::Error;

#[derive(Debug, Error)]
pub enum CommonError {
    /// sequenceUtil.cpp:260 "error, kmer length exceeds 32"
    #[error("error, kmer length exceeds 32: {len}")]
    KmerTooLong { len: usize },

    /// sequenceUtil.cpp:282 "error, kmer contains nongatc: {kmer}"
    #[error("error, kmer contains nongatc: {kmer}")]
    NonGatcChar { kmer: String },

    #[error("fasta format error: {0}")]
    FastaFormat(String),

    #[error("fastq format error: {0}")]
    FastqFormat(String),

    #[error(transparent)]
    Io(#[from] std::io::Error),
}
```

`crates/trinity-common/src/lib.rs` 全量替换为：
```rust
//! 共享基础库: 2bit k-mer 编码、FASTA/FASTQ IO — 直译原版 Inchworm/src/sequenceUtil.cpp 等。

pub mod error;
pub mod kmer;
```

- [x] **Step 2: 写失败测试**

`crates/trinity-common/src/kmer.rs`（实现暂写 `unimplemented!()`，测试先行）：

```rust
//! 2-bit k-mer 编码与操作 — 直译 Inchworm/src/sequenceUtil.cpp
//! 编码: G=0, A=1, T=2, C=3（互补 = 按位取反）; 小写 gatc 同样接受（_base_to_int 表）

use crate::error::CommonError;

/// 原版 kmer_int_type_t = unsigned long long（sequenceUtil.hpp:20）
pub type KmerId = u64;

/// sequenceUtil.cpp:258 — kmer 长度上限（64bit / 2bit）
pub const MAX_KMER_LENGTH: usize = 32;

/// _int_to_base 表（sequenceUtil.cpp:10）
pub const INT_TO_BASE: [u8; 4] = [b'G', b'A', b'T', b'C'];

pub fn base_to_int(c: u8) -> Option<u8> {
    match c {
        b'G' | b'g' => Some(0),
        b'A' | b'a' => Some(1),
        b'T' | b't' => Some(2),
        b'C' | b'c' => Some(3),
        _ => None,
    }
}

pub fn kmer_to_intval(kmer: &[u8]) -> Result<KmerId, CommonError> {
    unimplemented!("Task 3 Step 3")
}

pub fn decode_kmer_from_intval(intval: KmerId, kmer_length: usize) -> Vec<u8> {
    unimplemented!("Task 3 Step 3")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_base_encoding() {
        assert_eq!(base_to_int(b'G'), Some(0));
        assert_eq!(base_to_int(b'A'), Some(1));
        assert_eq!(base_to_int(b'T'), Some(2));
        assert_eq!(base_to_int(b'C'), Some(3));
        // 小写同表（sequenceUtil.cpp:12-24）
        assert_eq!(base_to_int(b'g'), Some(0));
        // 非 gatc（原版 255）
        assert_eq!(base_to_int(b'N'), None);
        assert_eq!(base_to_int(b'*'), None);
    }

    #[test]
    fn kmer_to_intval_hand_vectors() {
        assert_eq!(kmer_to_intval(b"G").unwrap(), 0);
        assert_eq!(kmer_to_intval(b"A").unwrap(), 1);
        assert_eq!(kmer_to_intval(b"T").unwrap(), 2);
        assert_eq!(kmer_to_intval(b"C").unwrap(), 3);
        // GA = (0<<2)|1 = 1（不同长度可同值，长度由调用方另行跟踪——原版同此性质）
        assert_eq!(kmer_to_intval(b"GA").unwrap(), 1);
        // ACGT = ((1<<2|3)<<2|0)<<2|2 = 114
        assert_eq!(kmer_to_intval(b"ACGT").unwrap(), 114);
        // 小写接受
        assert_eq!(kmer_to_intval(b"acgt").unwrap(), 114);
        // AAAA = 85
        assert_eq!(kmer_to_intval(b"AAAA").unwrap(), 85);
    }

    #[test]
    fn kmer_to_intval_errors() {
        assert!(matches!(
            kmer_to_intval(b"ACGN"),
            Err(CommonError::NonGatcChar { .. })
        ));
        let long = vec![b'A'; 33];
        assert!(matches!(
            kmer_to_intval(&long),
            Err(CommonError::KmerTooLong { len: 33 })
        ));
    }

    #[test]
    fn decode_hand_vectors() {
        assert_eq!(decode_kmer_from_intval(114, 4), b"ACGT".to_vec());
        assert_eq!(decode_kmer_from_intval(85, 4), b"AAAA".to_vec());
        assert_eq!(decode_kmer_from_intval(0, 1), b"G".to_vec());
        assert_eq!(decode_kmer_from_intval(1, 1), b"A".to_vec());
        assert_eq!(decode_kmer_from_intval(2, 1), b"T".to_vec());
        assert_eq!(decode_kmer_from_intval(3, 1), b"C".to_vec());
    }
}
```

- [x] **Step 3: 运行测试确认失败**

```bash
cargo test -p trinity-common 2>&1 | tail -8
```
Expected: `single_base_encoding` PASS；`kmer_to_intval_*` 与 `decode_hand_vectors` FAIL（panic: unimplemented）。

- [x] **Step 4: 写实现（替换两个 unimplemented!）**

```rust
/// sequenceUtil.cpp:258 kmer_to_intval — 逐字符 kmer_val<<2 | val；非 gatc 抛错（原版 cerr + throw）
pub fn kmer_to_intval(kmer: &[u8]) -> Result<KmerId, CommonError> {
    if kmer.len() > MAX_KMER_LENGTH {
        return Err(CommonError::KmerTooLong { len: kmer.len() });
    }
    let mut kmer_val: KmerId = 0;
    for &c in kmer {
        let val = base_to_int(c).ok_or_else(|| CommonError::NonGatcChar {
            kmer: String::from_utf8_lossy(kmer).into_owned(),
        })?;
        kmer_val <<= 2;
        kmer_val |= val as KmerId;
    }
    Ok(kmer_val)
}

/// sequenceUtil.cpp:298 decode_kmer_from_intval — 从低位端逐 2-bit 解出，写在逆序位置
pub fn decode_kmer_from_intval(intval: KmerId, kmer_length: usize) -> Vec<u8> {
    let mut kmer = vec![0u8; kmer_length];
    let mut v = intval;
    for i in 1..=kmer_length {
        let base_num = (v & 3) as usize;
        kmer[kmer_length - i] = INT_TO_BASE[base_num];
        v >>= 2;
    }
    kmer
}
```

- [x] **Step 5: 运行测试确认通过**

```bash
cargo test -p trinity-common 2>&1 | tail -4
```
Expected: 全部 PASS。

- [x] **Step 6: Commit**

```bash
git add crates/trinity-common
git commit -m "feat(common): CommonError 与 kmer 2bit 编码/解码（直译 sequenceUtil.cpp:258,298）"
```

---

### Task 4: revcomp 与 canonical

**Files:**
- Modify: `crates/trinity-common/src/kmer.rs`（追加函数与测试）

- [x] **Step 1: 写失败测试（追加到 tests 模块内）**

```rust
    #[test]
    fn revcomp_hand_vectors() {
        // revcomp("AA"=5) = "TT"=10
        assert_eq!(revcomp_val(5, 2), 10);
        // revcomp("AC"=7) = "GT"=2
        assert_eq!(revcomp_val(7, 2), 2);
        // ACGT 是回文: revcomp(ACGT)=ACGT
        assert_eq!(revcomp_val(114, 4), 114);
        // ACGTACGT 也是回文（val=29298）
        assert_eq!(revcomp_val(29298, 8), 29298);
        // 单碱基: A->T, G->C
        assert_eq!(revcomp_val(1, 1), 2);
        assert_eq!(revcomp_val(0, 1), 3);
    }

    #[test]
    fn canonical_hand_vectors() {
        // DS 规则 = max(kmer, revcomp)（sequenceUtil.cpp:376-383）
        assert_eq!(get_ds_kmer_val(5, 2), 10); // AA -> TT
        assert_eq!(get_ds_kmer_val(10, 2), 10); // TT -> TT
        assert_eq!(get_ds_kmer_val(7, 2), 7); // AC(7) > GT(2) -> 7
        assert_eq!(get_ds_kmer_val(114, 4), 114); // 回文不变
    }

    #[test]
    fn revcomp_roundtrip() {
        // 任意 kmer 双取 revcomp 复原
        let k = kmer_to_intval(b"AAAATAAAATAAAATAAAATAAAAT").unwrap();
        assert_eq!(revcomp_val(revcomp_val(k, 25), 25), k);
    }
```

- [x] **Step 2: 运行确认失败**

```bash
cargo test -p trinity-common 2>&1 | tail -6
```
Expected: 编译失败（函数未定义）。**本任务测试先行但实现立即补上，编译失败即视为红灯。**

- [x] **Step 3: 写实现（追加到 kmer.rs，测试模块之前）**

```rust
/// sequenceUtil.cpp:181 revcomp_val — ~kmer 完成互补，循环移位完成 2-bit 组反转。
/// 注意 ~ 会翻转全部 64 位，但循环只提取低 kmer_length 组，高位自然丢弃。
pub fn revcomp_val(mut kmer: KmerId, kmer_length: usize) -> KmerId {
    let mut rev_kmer: KmerId = 0;
    kmer = !kmer;
    for _ in 0..kmer_length {
        let base = kmer & 3;
        rev_kmer <<= 2;
        rev_kmer += base;
        kmer >>= 2;
    }
    rev_kmer
}

/// sequenceUtil.cpp:376 get_DS_kmer_val — canonical 形式 = max(kmer, revcomp(kmer))。
/// DS 模式下所有哈希键/visitor 键都必须先过这一步。
pub fn get_ds_kmer_val(kmer_val: KmerId, kmer_length: usize) -> KmerId {
    let rev_kmer = revcomp_val(kmer_val, kmer_length);
    if rev_kmer > kmer_val {
        rev_kmer
    } else {
        kmer_val
    }
}
```

- [x] **Step 4: 运行确认通过**

```bash
cargo test -p trinity-common 2>&1 | tail -4
```
Expected: 全部 PASS。

- [x] **Step 5: Commit**

```bash
git add crates/trinity-common
git commit -m "feat(common): revcomp/canonical（直译 sequenceUtil.cpp:181,376）"
```

---

### Task 5: 熵计算

**Files:**
- Modify: `crates/trinity-common/src/kmer.rs`（追加函数与测试）

- [x] **Step 1: 写失败测试（追加到 tests 模块内）**

```rust
    #[test]
    fn entropy_hand_vectors() {
        let acgt = kmer_to_intval(b"ACGT").unwrap();
        // 均匀分布 4 碱基: H = 2.0
        assert!((compute_entropy(acgt, 4) - 2.0).abs() < 1e-5);
        // 单一碱基: H = 0
        let aaaa = kmer_to_intval(b"AAAA").unwrap();
        assert!(compute_entropy(aaaa, 4).abs() < 1e-6);
        // AAAT: p(A)=0.75, p(T)=0.25 → 0.75*log2(4/3)+0.25*2 ≈ 0.811278
        let aaat = kmer_to_intval(b"AAAT").unwrap();
        assert!((compute_entropy(aaat, 4) - 0.8112781).abs() < 1e-5);
    }
```

- [x] **Step 2: 运行确认失败（编译错误：函数未定义）**

```bash
cargo test -p trinity-common 2>&1 | tail -4
```

- [x] **Step 3: 写实现**

```rust
/// sequenceUtil.cpp:316 compute_entropy — log2 香农熵。
/// 原版用 float（f32）逐项累加: prob * log(1/prob)/log(2.0f)。
/// 运算顺序与类型必须保持 f32，路径等价判定依赖浮点精确性（见 spec §6）。
pub fn compute_entropy(mut kmer: KmerId, kmer_length: usize) -> f32 {
    let mut counts = [0u32; 4];
    for _ in 0..kmer_length {
        let c = (kmer & 3) as usize;
        kmer >>= 2;
        counts[c] += 1;
    }
    let mut entropy = 0.0f32;
    for &cnt in &counts {
        let prob = cnt as f32 / kmer_length as f32;
        if prob > 0.0 {
            entropy += prob * (1.0f32 / prob).ln() / 2.0f32.ln();
        }
    }
    entropy
}
```

- [x] **Step 4: 运行确认通过**

```bash
cargo test -p trinity-common 2>&1 | tail -4
```

- [x] **Step 5: Commit**

```bash
git add crates/trinity-common
git commit -m "feat(common): kmer 香农熵，f32 运算序镜像原版（sequenceUtil.cpp:316）"
```

---

### Task 6: FASTA 读取

**Files:**
- Modify: `crates/trinity-common/src/lib.rs`（追加模块声明）
- Create: `crates/trinity-common/src/fasta.rs`

- [x] **Step 1: 更新 lib.rs 并写失败测试**

`crates/trinity-common/src/lib.rs` 全量替换为：
```rust
//! 共享基础库: 2bit k-mer 编码、FASTA/FASTQ IO — 直译原版 Inchworm/src/sequenceUtil.cpp 等。

pub mod error;
pub mod kmer;
pub mod fasta;
```
（error.rs 已在 Task 3 创建；fasta/fastq/io_util 模块在后续任务加入，此处先加 fasta。）

`crates/trinity-common/src/fasta.rs`（先测试后实现，实现暂 `unimplemented!()`）：

```rust
//! FASTA 读取 — 镜像 Inchworm/src/Fasta_reader.cpp + Fasta_entry.cpp
//! 行为: header 去 '>' 取整行为 _header；首个空白分隔 token 为 _accession（Fasta_entry.cpp:22-30）；
//! 序列大写化并去空白（Fasta_reader.cpp:121-122）。

use std::io::BufRead;

use crate::error::CommonError;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    /// 不含 '>' 的完整 header 行
    pub header: String,
    /// header 首个 token（按空格/制表符切）
    pub accession: String,
    /// 已大写、已去空白的序列
    pub sequence: String,
}

impl FastaRecord {
    /// Fasta_entry 构造逻辑
    pub fn new(header_line: &str, sequence: String) -> Self {
        let header = header_line.strip_prefix('>').unwrap_or(header_line).to_string();
        let accession = header
            .split([' ', '\t'])
            .next()
            .unwrap_or("")
            .to_string();
        FastaRecord { header, accession, sequence }
    }
}

pub struct FastaReader<R: BufRead> {
    reader: R,
    pending_header: Option<String>,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(reader: R) -> Self {
        FastaReader { reader, pending_header: None }
    }

    pub fn next_record(&mut self) -> Result<Option<FastaRecord>, CommonError> {
        unimplemented!("Task 6 Step 3")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    fn read_all(input: &str) -> Vec<FastaRecord> {
        let mut r = FastaReader::new(BufReader::new(input.as_bytes()));
        let mut out = Vec::new();
        while let Some(rec) = r.next_record().unwrap() {
            out.push(rec);
        }
        out
    }

    #[test]
    fn basic_single_record() {
        let recs = read_all(">acc1 some description\nacgt\nACGT\n");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].accession, "acc1");
        assert_eq!(recs[0].header, "acc1 some description");
        assert_eq!(recs[0].sequence, "ACGTACGT"); // 大写化
    }

    #[test]
    fn multiple_records_and_pending_header() {
        let recs = read_all(">a\nAAAA\n>b\ncccc\ngggg\n");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].sequence, "AAAA");
        assert_eq!(recs[1].sequence, "CCCCGGGG"); // 多行拼接 + 大写化
    }

    #[test]
    fn skips_blank_lines_and_crlf() {
        let recs = read_all("\n>a\nAC\r\nGT\n\n>b\nTT\n");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].sequence, "ACGT"); // \r 被剥、空行跳过
        assert_eq!(recs[1].sequence, "TT");
    }

    #[test]
    fn strips_whitespace_in_sequence() {
        let recs = read_all(">a\nAC GT\tAA\n");
        assert_eq!(recs[0].sequence, "ACGTAA"); // 去 ' ' 与 '\t'（原版 remove_whitespace 语义）
    }

    #[test]
    fn error_when_missing_gt() {
        let mut r = FastaReader::new(BufReader::new("notafasta\nACGT\n".as_bytes()));
        assert!(matches!(r.next_record(), Err(CommonError::FastaFormat(_))));
    }

    #[test]
    fn empty_input() {
        assert_eq!(read_all(""), Vec::new());
    }
}
```

- [x] **Step 2: 运行确认失败**

```bash
cargo test -p trinity-common 2>&1 | tail -6
```
Expected: FASTA 测试 panic（unimplemented），kmer 测试仍绿。

- [x] **Step 3: 写实现（替换 unimplemented!；`new` 已在 Step 1 定义，不动）**

```rust
impl<R: BufRead> FastaReader<R> {
    /// 返回下一条记录；EOF 返回 Ok(None)。
    /// 比原版略宽容：跳过空行（原版读到即用，Trinity 管道产物无空行）。
    pub fn next_record(&mut self) -> Result<Option<FastaRecord>, CommonError> {
        let header_line = match self.pending_header.take() {
            Some(h) => h,
            None => match self.read_line_trimmed(true)? {
                Some(l) if l.starts_with('>') => l,
                Some(l) => {
                    return Err(CommonError::FastaFormat(format!(
                        "记录未以 '>' 开始: {l}"
                    )))
                }
                None => return Ok(None), // EOF
            },
        };

        let mut sequence = String::new();
        loop {
            match self.read_line_trimmed(false)? {
                Some(l) if l.starts_with('>') => {
                    self.pending_header = Some(l);
                    break;
                }
                Some(l) => {
                    for c in l.chars() {
                        if !c.is_ascii_whitespace() {
                            sequence.push(c.to_ascii_uppercase());
                        }
                    }
                }
                None => break, // EOF 结束本记录
            }
        }
        Ok(Some(FastaRecord::new(&header_line, sequence)))
    }

    /// 读一行并去掉 \r\n；skip_empty=true 时跳过空行。
    fn read_line_trimmed(&mut self, skip_empty: bool) -> Result<Option<String>, CommonError> {
        let mut line = String::new();
        loop {
            line.clear();
            let n = self.reader.read_line(&mut line)?;
            if n == 0 {
                return Ok(None);
            }
            let t = line.trim_end_matches(['\n', '\r']);
            if skip_empty && t.trim().is_empty() {
                continue;
            }
            return Ok(Some(t.to_string()));
        }
    }
}
```

- [x] **Step 4: 运行确认通过**

```bash
cargo test -p trinity-common 2>&1 | tail -4
```

- [x] **Step 5: Commit**

```bash
git add crates/trinity-common
git commit -m "feat(common): FASTA 流式读取（镜像 Fasta_reader.cpp 行为：大写化/去空白/首 token accession）"
```

---

### Task 7: FASTA 写出（折行）

**Files:**
- Modify: `crates/trinity-common/src/fasta.rs`（追加函数与测试）

- [x] **Step 1: 写失败测试（追加到 tests 模块内）**

```rust
    #[test]
    fn line_breaks_at_interval() {
        // IRKE.cpp:386: 每 interval 字符换行，最后一组不换
        assert_eq!(add_fasta_seq_line_breaks(b"AAAAABBBBBCCCCC", 5), "AAAAA\nBBBBB\nCCCCC");
        // 16 字符: 位置 5/10/15 处换行（15 != 16 所以换）
        assert_eq!(
            add_fasta_seq_line_breaks(b"AAAAABBBBBCCCCCD", 5),
            "AAAAA\nBBBBB\nCCCCC\nD"
        );
        // 短于 interval 不折
        assert_eq!(add_fasta_seq_line_breaks(b"ACGT", 60), "ACGT");
        // interval=0 防御: 不折行
        assert_eq!(add_fasta_seq_line_breaks(b"ACGTACGT", 0), "ACGTACGT");
    }

    #[test]
    fn format_record() {
        let out = format_fasta_record("a1;25 total_counts: 30", b"ACGTACGT", 4);
        assert_eq!(out, ">a1;25 total_counts: 30\nACGT\nACGT\n");
    }
```

- [x] **Step 2: 运行确认失败（编译错误）**

```bash
cargo test -p trinity-common 2>&1 | tail -4
```

- [x] **Step 3: 写实现（追加到 fasta.rs）**

```rust
/// IRKE.cpp:386 add_fasta_seq_line_breaks — 每 interval 个字符插入换行，
/// 末尾不追加（调用方负责最后的 '\n'）。interval=0 视为不折行（防御，原版不会传 0）。
pub fn add_fasta_seq_line_breaks(sequence: &[u8], interval: usize) -> String {
    if interval == 0 || sequence.len() <= interval {
        return String::from_utf8_lossy(sequence).into_owned();
    }
    let mut out =
        String::with_capacity(sequence.len() + sequence.len() / interval + 8);
    for (i, &b) in sequence.iter().enumerate() {
        out.push(b as char);
        if (i + 1) % interval == 0 && i + 1 != sequence.len() {
            out.push('\n');
        }
    }
    out
}

/// 常规 FASTA 记录输出（header 不带 '>' 传入）。
pub fn format_fasta_record(header: &str, sequence: &[u8], interval: usize) -> String {
    format!(">{header}\n{}\n", add_fasta_seq_line_breaks(sequence, interval))
}
```

- [x] **Step 4: 运行确认通过**

```bash
cargo test -p trinity-common 2>&1 | tail -4
```

- [x] **Step 5: Commit**

```bash
git add crates/trinity-common
git commit -m "feat(common): FASTA 输出折行（镜像 IRKE.cpp:386）"
```

---

### Task 8: FASTQ 读取 + gzip 嗅探

**Files:**
- Modify: `crates/trinity-common/src/lib.rs`（追加 `pub mod fastq; pub mod io_util;`）
- Create: `crates/trinity-common/src/fastq.rs`
- Create: `crates/trinity-common/src/io_util.rs`

- [x] **Step 1: 写失败测试**

`crates/trinity-common/src/fastq.rs`：
```rust
//! FASTQ 读取 — 严格 4 行记录（@/seq/+/qual）。
//! 原版由 Perl 脚本转换，此为 Rust 版统一入口；校验比原版严格（报错带记录号与 header）。

use std::io::BufRead;

use crate::error::CommonError;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqRecord {
    pub accession: String, // header 去 '@' 后首个 token
    pub header: String,    // 去 '@' 的完整 header
    pub sequence: String,  // 原样保留大小写（不强制大写，与原版 fq->fa 转换一致）
    pub plus: String,      // '+' 行内容（可为空或重复 header）
    pub qual: String,
}

pub struct FastqReader<R: BufRead> {
    reader: R,
    record_no: u64,
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        FastqReader { reader, record_no: 0 }
    }

    pub fn next_record(&mut self) -> Result<Option<FastqRecord>, CommonError> {
        unimplemented!("Task 8 Step 3")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    fn read_all(input: &str) -> Vec<FastqRecord> {
        let mut r = FastqReader::new(BufReader::new(input.as_bytes()));
        let mut out = Vec::new();
        while let Some(rec) = r.next_record().unwrap() {
            out.push(rec);
        }
        out
    }

    #[test]
    fn basic_two_records() {
        let recs = read_all("@r1 desc\nACGT\n+\nIIII\n@r2\nTTTT\n+r2\n!!!!\n");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].accession, "r1");
        assert_eq!(recs[0].header, "r1 desc");
        assert_eq!(recs[0].sequence, "ACGT");
        assert_eq!(recs[0].plus, "");
        assert_eq!(recs[0].qual, "IIII");
        assert_eq!(recs[1].plus, "r2");
        assert_eq!(recs[1].qual, "!!!!");
    }

    #[test]
    fn tolerates_crlf_and_trailing_newline() {
        let recs = read_all("@r1\r\nACGT\r\n+\r\nIIII\r\n");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].sequence, "ACGT");
    }

    #[test]
    fn missing_at_header_is_error() {
        let mut r = FastqReader::new(BufReader::new("r1\nACGT\n+\nIIII\n".as_bytes()));
        let e = r.next_record().unwrap_err();
        assert!(matches!(e, CommonError::FastqFormat(_)));
    }

    #[test]
    fn seq_qual_length_mismatch_is_error() {
        let mut r = FastqReader::new(BufReader::new("@r1\nACGTA\n+\nIIII\n".as_bytes()));
        let e = r.next_record().unwrap_err();
        assert!(matches!(e, CommonError::FastqFormat(_)));
    }

    #[test]
    fn truncated_record_is_error() {
        let mut r = FastqReader::new(BufReader::new("@r1\nACGT\n".as_bytes()));
        assert!(matches!(r.next_record(), Err(CommonError::FastqFormat(_))));
    }
}
```

`crates/trinity-common/src/io_util.rs`：
```rust
//! gzip 魔数嗅探读入 — 纯 Rust 解压（flate2/miniz_oxide）。

use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::error::CommonError;

/// 打开文件，按魔数 (1f 8b) 自动套 gzip 解压层。
pub fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>, CommonError> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    let chain = file.chain(Cursor::new(magic[..n].to_vec()));
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(chain))))
    } else {
        Ok(Box::new(BufReader::new(chain)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plain_file_passthrough() {
        let dir = std::env::temp_dir().join("trinity_common_ioutil_t1");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("x.fa");
        std::fs::write(&p, b">a\nACGT\n").unwrap();
        let mut r = open_maybe_gz(&p).unwrap();
        let mut s = String::new();
        r.read_to_string(&mut s).unwrap();
        assert_eq!(s, ">a\nACGT\n");
    }

    #[test]
    fn gz_file_transparent() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;
        let dir = std::env::temp_dir().join("trinity_common_ioutil_t2");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("x.fq.gz");
        let mut enc = GzEncoder::new(File::create(&p).unwrap(), Compression::default());
        enc.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
        enc.finish().unwrap();
        let mut r = open_maybe_gz(&p).unwrap();
        let mut s = String::new();
        r.read_to_string(&mut s).unwrap();
        assert_eq!(s, "@r1\nACGT\n+\nIIII\n");
    }
}
```

- [x] **Step 2: 运行确认失败**

```bash
cargo test -p trinity-common 2>&1 | tail -6
```
Expected: io_util 测试通过（无需实现——io_util 一次写完）；fastq 测试 panic（unimplemented）。

- [x] **Step 3: 写 fastq 实现（替换 unimplemented!）**

```rust
impl<R: BufRead> FastqReader<R> {
    /// 读一条记录；容忍文件级结尾多余空行；EOF 返回 Ok(None)。
    pub fn next_record(&mut self) -> Result<Option<FastqRecord>, CommonError> {
        let header_line = match self.next_nonempty_line()? {
            Some(l) => l,
            None => return Ok(None),
        };
        self.record_no += 1;
        let where_ = format!("记录 #{} ({})", self.record_no, header_line);

        if !header_line.starts_with('@') {
            return Err(CommonError::FastqFormat(format!(
                "header 未以 '@' 开始: {where_}"
            )));
        }
        let sequence = self
            .next_nonempty_line()?
            .ok_or_else(|| CommonError::FastqFormat(format!("缺失序列行: {where_}")))?;
        let plus_line = self
            .next_nonempty_line()?
            .ok_or_else(|| CommonError::FastqFormat(format!("缺失 '+' 行: {where_}")))?;
        if !plus_line.starts_with('+') {
            return Err(CommonError::FastqFormat(format!(
                "第三行未以 '+' 开始: {where_}"
            )));
        }
        let qual = self
            .next_nonempty_line()?
            .ok_or_else(|| CommonError::FastqFormat(format!("缺失质量行: {where_}")))?;
        if qual.len() != sequence.len() {
            return Err(CommonError::FastqFormat(format!(
                "序列({})与质量({})长度不一致: {where_}",
                sequence.len(),
                qual.len()
            )));
        }

        let header = header_line[1..].to_string();
        let accession = header.split([' ', '\t']).next().unwrap_or("").to_string();
        Ok(Some(FastqRecord {
            accession,
            header,
            sequence,
            plus: plus_line[1..].to_string(),
            qual,
        }))
    }

    fn next_nonempty_line(&mut self) -> Result<Option<String>, CommonError> {
        let mut line = String::new();
        loop {
            line.clear();
            let n = self.reader.read_line(&mut line)?;
            if n == 0 {
                return Ok(None);
            }
            let t = line.trim_end_matches(['\n', '\r']);
            if t.is_empty() {
                continue; // 文件级空行（含结尾多余换行）
            }
            return Ok(Some(t.to_string()));
        }
    }
}
```

- [x] **Step 4: 运行确认通过**

```bash
cargo test -p trinity-common 2>&1 | tail -4
```

- [x] **Step 5: Commit**

```bash
git add crates/trinity-common
git commit -m "feat(common): FASTQ 严格读取与 gzip 魔数嗅探"
```

---

### Task 9: 原版 C++ 黄金向量 harness 与位级等价测试

**Files:**
- Create: `xtask/fixtures-src/dump_kmer_golden.cpp`
- Create: `fixtures/kmer_golden_input.txt`（checked in）
- Create: `fixtures/kmer_golden.tsv`（生成后 checked in）
- Create: `crates/trinity-common/tests/kmer_golden.rs`

- [x] **Step 1: 写黄金向量输入**

`fixtures/kmer_golden_input.txt`（每行一个 kmer，共 26 行）：
```
G
A
T
C
AA
AC
AG
AT
CA
CC
CG
CT
GA
GC
GG
GT
TA
TC
TG
TT
ACGT
ACGTACGT
acgt
AAAATAAAATAAAATAAAATAAAAT
GGGGCGGGGCGGGGCGGGGCGGGGC
TTCCTTCCATCCTTACCCTTTTCAA
```

- [x] **Step 2: 写 C++ harness（链接原版源码，保证向量 100% 来自原版实现）**

`xtask/fixtures-src/dump_kmer_golden.cpp`：
```cpp
// 黄金向量生成器: 直接链接原版 sequenceUtil.cpp，从 stdin 逐行读 kmer，
// 输出 TSV: kmer \t intval \t revcomp \t dsval \t entropy
// 编译: g++ -O2 -I$TRINITY_SRC/Inchworm/src 本文件 $TRINITY_SRC/Inchworm/src/sequenceUtil.cpp \
//       $TRINITY_SRC/Inchworm/src/stacktrace.cpp -o dump_kmer_golden
#include <cstdio>
#include <iostream>
#include <string>
#include "sequenceUtil.hpp"

int main() {
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        kmer_int_type_t v = kmer_to_intval(line);
        unsigned int k = (unsigned int)line.length();
        printf("%s\t%llu\t%llu\t%llu\t%.9g\n",
               line.c_str(),
               (unsigned long long)v,
               (unsigned long long)revcomp_val(v, k),
               (unsigned long long)get_DS_kmer_val(v, k),
               compute_entropy(v, k));
    }
    return 0;
}
```

- [x] **Step 3: 手工编译并生成 fixture（一次性，Task 10 会自动化）**

```bash
cd /storage/home/senior007/test/trinity_rust/trinity-rust
SRC=/storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2
mkdir -p target/fixture-tools
g++ -O2 -I"$SRC/Inchworm/src" xtask/fixtures-src/dump_kmer_golden.cpp \
    "$SRC/Inchworm/src/sequenceUtil.cpp" "$SRC/Inchworm/src/stacktrace.cpp" \
    -o target/fixture-tools/dump_kmer_golden
./target/fixture-tools/dump_kmer_golden < fixtures/kmer_golden_input.txt > fixtures/kmer_golden.tsv
head -5 fixtures/kmer_golden.tsv
wc -l fixtures/kmer_golden.tsv
```
Expected: 26 行 TSV。前几行形如：
```
G	0	3	3	0
A	1	2	2	0
...
ACGT	114	114	114	2
```
若链接报 undefined symbol（如 string_util）：把缺失符号所在 .cpp 追加到编译命令（原版各工具这样拼装）。

- [x] **Step 4: 写位级等价测试**

`crates/trinity-common/tests/kmer_golden.rs`：
```rust
//! 黄金向量: 与原版 C++（sequenceUtil.cpp 直链 harness）位级一致。

use trinity_common::kmer::{compute_entropy, decode_kmer_from_intval, get_ds_kmer_val, kmer_to_intval, revcomp_val};

#[test]
fn kmer_ops_match_original_cpp_bit_for_bit() {
    let tsv = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../fixtures/kmer_golden.tsv"
    ))
    .expect("缺少 fixtures/kmer_golden.tsv — 先跑 cargo xtask gen-fixtures");
    let mut n = 0;
    for line in tsv.lines() {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 5, "行格式异常: {line}");
        let kmer = cols[0].as_bytes();
        let intval: u64 = cols[1].parse().unwrap();
        let revcomp: u64 = cols[2].parse().unwrap();
        let dsval: u64 = cols[3].parse().unwrap();
        // %.9g 十进制位足以唯一确定 f32 → 解析回的 f32 与原版逐位一致
        let entropy: f32 = cols[4].parse().unwrap();
        let k = kmer.len();

        assert_eq!(kmer_to_intval(kmer).unwrap(), intval, "intval 不一致: {line}");
        assert_eq!(revcomp_val(intval, k), revcomp, "revcomp 不一致: {line}");
        assert_eq!(get_ds_kmer_val(intval, k), dsval, "dsval 不一致: {line}");
        assert_eq!(compute_entropy(intval, k), entropy, "entropy 不一致: {line}");
        // decode 与输入大写形式互逆（小写输入编码相同、解码为大写）
        assert_eq!(decode_kmer_from_intval(intval, k), kmer.to_ascii_uppercase());

        n += 1;
    }
    assert!(n >= 26, "黄金向量行数异常: {n}");
}
```
注意 `compute_entropy` 断言用 `assert_eq!`（f32 精确相等）——这是 P0 的核心命题：f32 运算序逐位对齐。如果在此失败，说明 Rust 与 glibc 的 logf 有位差异，**不要**改成容差比较，先调查（两者在 Linux 都调用系统 libm logf，正常应当一致）。

- [x] **Step 5: 运行确认通过**

```bash
cargo test -p trinity-common --test kmer_golden 2>&1 | tail -4
```
Expected: PASS（1 个测试，26 行向量全对）。

- [x] **Step 6: Commit**

```bash
git add xtask/fixtures-src fixtures crates/trinity-common/tests
git commit -m "test(common): 原版 C++ 直链黄金向量, kmer 运算位级等价验证"
```

---

### Task 10: xtask gen-fixtures 命令

**Files:**
- Modify: `xtask/src/main.rs`（全量替换）

- [x] **Step 1: 写实现**

```rust
use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("gen-fixtures") => gen_kmer_golden(),
        Some(other) => {
            eprintln!("未知任务: {other}\n用法: cargo xtask <gen-fixtures>");
            std::process::exit(2);
        }
        None => {
            eprintln!("用法: cargo xtask <gen-fixtures>");
            std::process::exit(2);
        }
    }
}

fn trinity_src() -> PathBuf {
    env::var_os("TRINITY_SRC")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from("/storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2")
        })
}

fn workspace_root() -> PathBuf {
    Path::new(&env::var("CARGO_MANIFEST_DIR").unwrap())
        .parent()
        .unwrap()
        .to_path_buf()
}

fn gen_kmer_golden() {
    let root = workspace_root();
    let src = trinity_src();
    assert!(
        src.join("Inchworm/src/sequenceUtil.cpp").exists(),
        "找不到原版源码 {} — 请设置 TRINITY_SRC",
        src.display()
    );
    let out_dir = root.join("target/fixture-tools");
    std::fs::create_dir_all(&out_dir).unwrap();
    let bin = out_dir.join("dump_kmer_golden");
    let harness = root.join("xtask/fixtures-src/dump_kmer_golden.cpp");

    let status = Command::new("g++")
        .arg("-O2")
        .arg("-I")
        .arg(src.join("Inchworm/src"))
        .arg(&harness)
        .arg(src.join("Inchworm/src/sequenceUtil.cpp"))
        .arg(src.join("Inchworm/src/stacktrace.cpp"))
        .arg("-o")
        .arg(&bin)
        .status()
        .expect("无法启动 g++（是否安装?）");
    assert!(status.success(), "g++ 编译 harness 失败");

    let input = root.join("fixtures/kmer_golden_input.txt");
    let output = root.join("fixtures/kmer_golden.tsv");
    let status = Command::new(&bin)
        .stdin(std::fs::File::open(&input).unwrap())
        .stdout(std::fs::File::create(&output).unwrap())
        .status()
        .unwrap();
    assert!(status.success(), "运行 harness 失败");
    println!("已生成 {}", output.display());
}
```

- [x] **Step 2: 运行验证可重复生成**

```bash
cargo xtask gen-fixtures
md5sum fixtures/kmer_golden.tsv
cargo xtask gen-fixtures
md5sum fixtures/kmer_golden.tsv
```
Expected: 两次 md5 相同（确定性输出），且与 Task 9 生成的文件一致（`git status` 不显示变更）。

- [x] **Step 3: 全 workspace 测试**

```bash
cargo test --workspace 2>&1 | tail -6
```
Expected: 全绿。

- [x] **Step 4: Commit**

```bash
git add xtask fixtures
git commit -m "feat(xtask): gen-fixtures 命令——重建原版 C++ 黄金向量"
```

---

### Task 11: P0 验证门

**Files:**
- 无新文件（检查与收尾）

- [x] **Step 1: 验证门清单逐项执行**

```bash
cd /storage/home/senior007/test/trinity_rust/trinity-rust
cargo build --workspace 2>&1 | tail -2          # Finished
cargo test --workspace 2>&1 | tail -3           # 0 failed
cargo xtask gen-fixtures && git status --short fixtures/   # 无变更（可重复生成）
cargo test -p trinity-common --test kmer_golden 2>&1 | tail -3   # PASS
```
Expected: 四项全部通过。spec §11 P0 验证门 = 「黄金向量单测全绿；原版二进制可用」达成。

- [x] **Step 2: 原版二进制可用性复查（Task 1 已建，此处确认仍在）**

```bash
ls /storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2/Inchworm/bin/inchworm \
   /storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2/Chrysalis/bin/GraphFromFasta \
   /storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2/Butterfly/Butterfly.jar
```
Expected: 三个代表性文件存在（完整清单见 docs/setup.md）。

- [x] **Step 3: 收尾提交（如有未提交改动）**

```bash
git status --short
git add -A && git commit -m "chore: P0 验证门通过（workspace + common + 黄金向量 + 原版构建）" || echo "无待提交改动"
git log --oneline
```

---

## P0 完成标准（对应 spec §11 验证门）

1. `cargo build/test --workspace` 全绿
2. 黄金向量测试证明 kmer 编码/revcomp/canonical/熵与原版 C++ **位级一致**（f32 精确相等）
3. `cargo xtask gen-fixtures` 确定性地重建黄金向量
4. 原版 Trinity 二进制可用（inchworm smoke 已跑通）

通过后进入 P1（trinity-kmer：并行计数 + dump + DigiNorm），届时另写 `docs/superpowers/plans/<日期>-p1-kmer.md`。
