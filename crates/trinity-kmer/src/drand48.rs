//! Perl_drand48 兼容 PRNG —— perl >= 5.20 的 rand() 即经典 drand48
//! （perl 5.38.2 实测，见 fixtures/p1/perl_rand_12345.txt）。
//! 目的: 复刻 nbkc_normalize.pl 的 srand(12345)+rand(1) 序列，使选择名单与原版逐字节一致。
//!
//! 算法: 48-bit LCG，X' = (a·X + c) mod 2^48，a = 0x5DEECE66D，c = 0xB;
//! srand48(seed): X = ((seed << 16) | 0x330E) mod 2^48（seed 取低 32 位，perl 的
//! seedDrand01 收 U32）; drand48() = X / 2^48。
//!
//! 位级验证: perl `unpack("Q<", pack("d<", rand(1)))` 前 1000 个值与本实现
//! 全部位一致（X/2^48 是精确的 f64——48 位整数仅移指数，无舍入）。
//! 注意 perl 的 `print rand(1)` 默认 %.15g 输出，**不能**逐字符比对（Rust Display
//! 是最短往返表示），也不能解析回 f64 比位（%.15g 截断后有损）——黄金文件改存
//! %.17g（IEEE 保证往返），测试解析后按位断言。

/// 48-bit LCG 状态（仅低 48 位有效）。
pub struct Drand48 {
    state: u64,
}

impl Drand48 {
    /// srand48(seed) 镜像。seed 按无符号 32 位截断（perl seedDrand01 收 U32，
    /// srand(12345) 常量场景不受影响）。
    pub fn new(seed: u64) -> Self {
        let seed32 = (seed & 0xFFFF_FFFF) as u32;
        Drand48 {
            state: (((seed32 as u64) << 16) | 0x330E) & 0xFFFF_FFFF_FFFF,
        }
    }

    /// drand48() 镜像: 返回 [0, 1) 的 f64，与 perl rand(1) 位级一致。
    pub fn next_f64(&mut self) -> f64 {
        self.state = (self.state.wrapping_mul(0x5DEECE66D).wrapping_add(0xB)) & 0xFFFF_FFFF_FFFF;
        self.state as f64 / (1u64 << 48) as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        format!("{}/../../fixtures/p1/{name}", env!("CARGO_MANIFEST_DIR"))
    }

    /// 黄金序列位级锁定: fixtures/p1/perl_rand_12345.txt 为 perl printf("%.17g")
    /// 输出（IEEE 保证 17 位有效数字往返），解析后逐值按位比对。
    #[test]
    fn golden_sequence_bit_exact() {
        let text = std::fs::read_to_string(fixture("perl_rand_12345.txt")).unwrap();
        let mut rng = Drand48::new(12345);
        let mut n = 0;
        for line in text.lines() {
            let expect: f64 = line.trim().parse().unwrap();
            let got = rng.next_f64();
            assert_eq!(got.to_bits(), expect.to_bits(), "第 {n} 个值位不一致");
            n += 1;
        }
        assert_eq!(n, 20);
    }

    /// 前 5 个值的原始位模式硬编码（perl unpack("Q<",pack("d<",...)) 实测），
    /// 双保险: 即使 fixture 文件被误重生成也能立刻发现。
    #[test]
    fn first_five_raw_bits_locked() {
        let expect = [
            4597286335540658304u64,
            4606454484595142400,
            4596620261818962304,
            4604703454901591616,
            4604770283139736224,
        ];
        let mut rng = Drand48::new(12345);
        for (i, e) in expect.iter().enumerate() {
            assert_eq!(rng.next_f64().to_bits(), *e, "index {i}");
        }
    }

    /// 与 perl 默认 print（%.15g）输出逐字符对照——仅证明序列同源，
    /// 不作为实现断言（%.15g 有损，见模块文档）。
    #[test]
    fn matches_perl_15g_print() {
        let mut rng = Drand48::new(12345);
        // 值 ∈ [0.1,1) 时 %.15g == 15 位定点; perl 实测输出见下方字面量
        for expect in [
            "0.225328512796299",
            "0.919183068533556",
            "0.206841253248182",
        ] {
            assert_eq!(format!("{:.15}", rng.next_f64()), expect);
        }
    }

    #[test]
    fn state_init_and_range() {
        // srand48(12345): X0 = (12345<<16)|0x330E = 0x3039330E（48 位内）
        let mut rng = Drand48::new(12345);
        for _ in 0..1000 {
            let v = rng.next_f64();
            assert!((0.0..1.0).contains(&v));
        }
        // 相同种子 → 相同序列
        let mut a = Drand48::new(7);
        let mut b = Drand48::new(7);
        for _ in 0..10 {
            assert_eq!(a.next_f64().to_bits(), b.next_f64().to_bits());
        }
        // 种子按 32 位截断（镜像 perl U32 seedDrand01）
        let mut lo = Drand48::new(0x1_0000_0002u64);
        let mut hi = Drand48::new(2);
        assert_eq!(lo.next_f64().to_bits(), hi.next_f64().to_bits());
        // X0 的低位 0x330E 保证首个输出不为 0
        let mut s = Drand48::new(0);
        assert!(s.next_f64() > 0.0);
    }
}
