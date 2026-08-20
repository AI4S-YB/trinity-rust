//! glibc random() TYPE_3（additive feedback，degree 31 / separation 3）复刻。
//!
//! 原版 inchworm 从不调 srandom → rand() 即 srand(1) 序列，可逐位复现
//! （inchworm_step 的真平局分支用 rand()%2 二选一——见计划 Task 4）。
//!
//! 算法（对照 glibc stdlib/random_r.c 的 __srandom_r / __random_r）:
//! - 初始化: seed == 0 时取 1；r[0] = seed；i ∈ [1,31):
//!   r[i] = (16807 * r[i-1]) mod 2147483647（glibc 用 Schrage 防 31 位溢出，
//!   数学上等价于 i64 乘模后转 i32，结果恒非负）
//! - 31 ≡ -3 (mod 34) 的别名展开: r[31]=r[0], r[32]=r[1], r[33]=r[2]，
//!   此后递推 r[i] = r[i-31] + r[i-3]（i32 回绕）在 34 元环上滚动，
//!   与 glibc 的 fptr/rptr 双指针（fptr=state+3 起步、越界回绕）逐值等价
//! - 每次 srandom 后丢弃前 10*31 = 310 个输出
//! - 输出: (r[i] as u32) >> 1（丢弃最低随机位）
//!
//! 黄金: fixtures/p2/glibcrand_seed1.txt（C srand(1)+rand() 100 值）与
//! glibcrand_mod2.txt（rand()%2 50 值）位级一致。

/// glibc TYPE_3 状态（34 元环 + 回绕索引）。
pub struct GlibcRand {
    r: [i32; 34],
    idx: usize,
}

impl GlibcRand {
    /// srand(seed) 镜像（含 seed==0 → 1 与丢弃 310 输出）。
    pub fn new(seed: u32) -> Self {
        let seed = if seed == 0 { 1 } else { seed };
        let mut r = [0i32; 34];
        r[0] = seed as i32;
        for i in 1..31 {
            // Schrage 的数学等价式: 16807 * r[i-1] ≤ 16807 * 2^31 < 2^53，i64 直算无损
            r[i] = ((16807i64 * r[i - 1] as i64) % 2147483647) as i32;
        }
        r[31] = r[0];
        r[32] = r[1];
        r[33] = r[2];
        let mut rng = GlibcRand { r, idx: 0 };
        for _ in 0..310 {
            rng.next();
        }
        rng
    }

    /// rand() 镜像: r[i] = r[i-31] + r[i-3]（i32 回绕），返回 (r[i] as u32) >> 1。
    /// 名字对齐 C 侧 rand() 的"每次调用一个值"语义，非迭代器（计划 Task 2 接口）。
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> u32 {
        // 31 ≡ -3 (mod 34): i-31 与 i-3 的槽位分别是 (idx+3) 与 (idx+31) 对 34 取模
        self.r[self.idx] = self.r[(self.idx + 31) % 34].wrapping_add(self.r[(self.idx + 3) % 34]);
        let out = (self.r[self.idx] as u32) >> 1;
        self.idx = (self.idx + 1) % 34;
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        format!("{}/../../fixtures/p2/{name}", env!("CARGO_MANIFEST_DIR"))
    }

    #[test]
    fn golden_seed1_100_values() {
        let text = std::fs::read_to_string(fixture("glibcrand_seed1.txt")).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 100);
        let mut rng = GlibcRand::new(1);
        for (i, line) in lines.iter().enumerate() {
            assert_eq!(rng.next(), line.parse::<u32>().unwrap(), "第 {i} 个 rand()");
        }
    }

    #[test]
    fn golden_seed1_mod2_50_values() {
        let text = std::fs::read_to_string(fixture("glibcrand_mod2.txt")).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 50);
        let mut rng = GlibcRand::new(1);
        for (i, line) in lines.iter().enumerate() {
            assert_eq!(
                rng.next() % 2,
                line.parse::<u32>().unwrap(),
                "第 {i} 个 rand()%2"
            );
        }
    }

    #[test]
    fn famous_first_five() {
        // glibc srand(1) 的公认首 5 值（文献常引；亦被黄金 fixture 首五行锁定）
        let mut rng = GlibcRand::new(1);
        let got = [rng.next(), rng.next(), rng.next(), rng.next(), rng.next()];
        assert_eq!(
            got,
            [1804289383u32, 846930886, 1681692777, 1714636915, 1957747793]
        );
    }

    #[test]
    fn seed_zero_behaves_as_one() {
        // glibc __srandom_r: seed == 0 时取 1
        assert_eq!(GlibcRand::new(0).next(), GlibcRand::new(1).next());
    }
}
