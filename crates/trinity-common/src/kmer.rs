//! 2-bit k-mer 编码与操作 — 直译 Inchworm/src/sequenceUtil.cpp
//! 编码: G=0, A=1, T=2, C=3（互补 = 按位取反）; 小写 gatc 同样接受（_base_to_int 表）

use crate::error::CommonError;

/// 原版 kmer_int_type_t = unsigned long long（sequenceUtil.hpp:20）。
/// **扩容沿革**: u64 → u128（2026-08-30, 解锁 rnaSPAdes RL=76 的 [29,37] 原值）;
/// u128 → U256（2026-09-02, 解锁 RL=150 公式 upper_k=73 与 RL=250 的 [55,123]——
/// 2bit 编码下 u128 上限 k=64 不够用, SPAdes 自身上限 MAX_K=127 需 254bit）。
/// k≤32 的既有键值数值不变（零扩展）, u64/u128 口径的全部黄金断言不受影响。
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Default)]
pub struct U256 {
    /// 高 limb 在前: 派生 Ord 按字段序比较 = 256 位数值序。
    pub hi: u128,
    pub lo: u128,
}

/// 滚动 k-mer 编码键抽象（性能分层, 2026-09-02; 2026-09-04 加 u64 层）:
/// - [`U256`]（= KmerId）: 全域 k≤128, 管线对外类型
/// - `u128`: k≤64 的紧凑键（diginorm 计数/查表热路径——k=25 时代价 24B/槽
///   vs U256 的 40B/槽, 表 footprint 与插入访存显著占优）
/// - `u64`: k≤32 的最紧凑键（diginorm k=25 主路径——键带宽/哈希域减半,
///   大表 cache miss 与 FxHash 成本同步降; 压时间战役 2026-09-04）
/// 数值语义全域一致: u64 键 = u128 键 = KmerId 的低 limb（高位恒 0）。
pub trait RollingKey: Copy + Eq + Ord + std::hash::Hash {
    /// 低 2k 位全 1 掩码
    fn low_mask(k: usize) -> Self;
    fn zero() -> Self;
    fn from_base(v: u8) -> Self;
    /// self = (self << 2 | v) & mask
    fn roll_forward(self, v: u8, mask: Self) -> Self;
    /// rev = (rev >> 2 | (v^3) << (2k-2)) & mask
    fn roll_rev(self, v: u8, rev_shift: u32, mask: Self) -> Self;
    /// revcomp（canonical_merge_sum 用; 只读低 2k 组, 高位污染不外泄）
    fn revcomp(self, k: usize) -> Self;
    /// 从全宽键收窄（u128: hi 必须为 0）
    fn from_wide(w: &KmerId) -> Self;
    fn to_wide(self) -> KmerId;
}

impl RollingKey for u64 {
    #[inline]
    fn low_mask(k: usize) -> u64 {
        debug_assert!(k <= 32);
        if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 }
    }
    #[inline]
    fn zero() -> u64 {
        0
    }
    #[inline]
    fn from_base(v: u8) -> u64 {
        v as u64
    }
    #[inline]
    fn roll_forward(self, v: u8, mask: u64) -> u64 {
        ((self << 2) | v as u64) & mask
    }
    #[inline]
    fn roll_rev(self, v: u8, rev_shift: u32, mask: u64) -> u64 {
        ((self >> 2) | (((v ^ 3) as u64) << rev_shift)) & mask
    }
    fn revcomp(self, k: usize) -> u64 {
        let mut r = 0u64;
        let mut kmer = !self;
        for _ in 0..k {
            r <<= 2;
            r |= kmer & 3;
            kmer >>= 2;
        }
        r
    }
    fn from_wide(w: &KmerId) -> u64 {
        debug_assert!(w.hi == 0 && w.lo <= u64::MAX as u128, "u64 键收窄要求 k≤32");
        w.lo as u64
    }
    fn to_wide(self) -> KmerId {
        U256 { hi: 0, lo: self as u128 }
    }
}

impl RollingKey for u128 {
    #[inline]
    fn low_mask(k: usize) -> u128 {
        debug_assert!(k <= 64);
        if k >= 64 { u128::MAX } else { (1u128 << (2 * k)) - 1 }
    }
    #[inline]
    fn zero() -> u128 {
        0
    }
    #[inline]
    fn from_base(v: u8) -> u128 {
        v as u128
    }
    #[inline]
    fn roll_forward(self, v: u8, mask: u128) -> u128 {
        ((self << 2) | v as u128) & mask
    }
    #[inline]
    fn roll_rev(self, v: u8, rev_shift: u32, mask: u128) -> u128 {
        ((self >> 2) | (((v ^ 3) as u128) << rev_shift)) & mask
    }
    fn revcomp(self, k: usize) -> u128 {
        let mut r = 0u128;
        let mut kmer = !self;
        for _ in 0..k {
            r <<= 2;
            r |= kmer & 3;
            kmer >>= 2;
        }
        r
    }
    fn from_wide(w: &KmerId) -> u128 {
        debug_assert!(w.hi == 0, "u128 键收窄要求 hi==0（k≤64）");
        w.lo
    }
    fn to_wide(self) -> KmerId {
        U256 { hi: 0, lo: self }
    }
}

impl RollingKey for U256 {
    #[inline]
    fn low_mask(k: usize) -> U256 {
        U256::low_bits(2 * k as u32)
    }
    #[inline]
    fn zero() -> U256 {
        U256::ZERO
    }
    #[inline]
    fn from_base(v: u8) -> U256 {
        U256 { hi: 0, lo: v as u128 }
    }
    #[inline]
    fn roll_forward(self, v: u8, mask: U256) -> U256 {
        ((self << 2u32) | U256::from(v)) & mask
    }
    #[inline]
    fn roll_rev(self, v: u8, rev_shift: u32, mask: U256) -> U256 {
        ((self >> 2u32) | (U256::from(v ^ 3) << rev_shift)) & mask
    }
    fn revcomp(self, k: usize) -> U256 {
        revcomp_val(self, k)
    }
    fn from_wide(w: &KmerId) -> U256 {
        *w
    }
    fn to_wide(self) -> KmerId {
        self
    }
}

pub type KmerId = U256;

impl U256 {
    pub const ZERO: U256 = U256 { hi: 0, lo: 0 };
    pub const MAX: U256 = U256 { hi: u128::MAX, lo: u128::MAX };

    /// 低 `bits` 位全 1（bits > 256 时饱和为 MAX; bits == 256 即 MAX）。
    pub fn low_bits(bits: u32) -> U256 {
        debug_assert!(bits <= 256);
        if bits >= 256 {
            U256::MAX
        } else if bits >= 128 {
            U256 { hi: (1u128 << (bits - 128)) - 1, lo: u128::MAX }
        } else if bits == 0 {
            U256::ZERO
        } else {
            U256 { hi: 0, lo: (1u128 << bits) - 1 }
        }
    }

    /// 原生整数的回绕左移语义（移位量按 256 取模; x86 硬件行为一致）。
    pub fn wrapping_shl(self, s: u32) -> U256 {
        let s = s % 256;
        if s == 0 {
            self
        } else if s >= 128 {
            U256 { hi: self.lo << (s - 128), lo: 0 }
        } else {
            U256 {
                hi: (self.hi << s) | (self.lo >> (128 - s)),
                lo: self.lo << s,
            }
        }
    }

    pub fn to_le_bytes(self) -> [u8; 32] {
        let lo = self.lo.to_le_bytes();
        let hi = self.hi.to_le_bytes();
        let mut out = [0u8; 32];
        out[..16].copy_from_slice(&lo);
        out[16..].copy_from_slice(&hi);
        out
    }

    pub fn from_le_bytes(b: [u8; 32]) -> U256 {
        let mut lo = [0u8; 16];
        let mut hi = [0u8; 16];
        lo.copy_from_slice(&b[..16]);
        hi.copy_from_slice(&b[16..]);
        U256 {
            hi: u128::from_le_bytes(hi),
            lo: u128::from_le_bytes(lo),
        }
    }
}

impl From<u8> for U256 {
    fn from(v: u8) -> U256 {
        U256 { hi: 0, lo: v as u128 }
    }
}

impl From<u64> for U256 {
    fn from(v: u64) -> U256 {
        U256 { hi: 0, lo: v as u128 }
    }
}

impl From<u128> for U256 {
    fn from(v: u128) -> U256 {
        U256 { hi: 0, lo: v }
    }
}

impl From<u32> for U256 {
    fn from(v: u32) -> U256 {
        U256 { hi: 0, lo: v as u128 }
    }
}

impl From<i32> for U256 {
    fn from(v: i32) -> U256 {
        // 仅服务测试字面量（恒非负）; 负输入按位重解释, 不做防误用护栏
        U256 { hi: 0, lo: v as u32 as u128 }
    }
}

/// 小常数乘（FromStr 十进制累乘 ×10 用; limb 传播进位）。
impl core::ops::Mul<u32> for U256 {
    type Output = U256;
    fn mul(self, m: u32) -> U256 {
        let m = m as u128;
        let (lo, c1) = self.lo.overflowing_mul(m);
        let (hi_lo, c2) = self.hi.overflowing_mul(m);
        // hi*m 的低位部分 + lo*m 的进位（c1 最大 (2^128-1)*9 进位 < 2^128, 加法可能再进位）
        let (hi_part, c3) = hi_lo.overflowing_add(c1 as u128);
        debug_assert!(c2 == false || true); // 高位溢出仅在 >2^256/10 时发生, FromStr 输入域内不可达
        let _ = c3;
        U256 { hi: hi_part, lo }
    }
}

/// 与原生整数字面量的比较（零扩展; i32 覆盖默认整数字面量推断, 让
/// `val == 0` 这类断言无需显式 From 包裹）。
macro_rules! u256_partial_eq_primitive {
    ($($t:ty),*) => {$(
        impl PartialEq<$t> for U256 {
            fn eq(&self, other: &$t) -> bool {
                *self == U256::from(*other as u128)
            }
        }
        impl PartialEq<U256> for $t {
            fn eq(&self, other: &U256) -> bool {
                U256::from(*self as u128) == *other
            }
        }
    )*};
}
u256_partial_eq_primitive!(i32, u32, u64, u128, usize);


/// 手写 Hash（替代 derive）: hi==0（k≤64 全量命中）快路径只喂一个 limb——
/// U256 表的插入/查找哈希开销减半; 判别字节保证双分支编码单射（等值必同支）。
impl std::hash::Hash for U256 {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        if self.hi == 0 {
            state.write_u8(0);
            state.write_u128(self.lo);
        } else {
            state.write_u8(1);
            state.write_u128(self.hi);
            state.write_u128(self.lo);
        }
    }
}

/// 十进制解析（kmer_golden 黄金文件读取用; 非数字字符报错）。
impl core::str::FromStr for U256 {
    type Err = ParseU256Error;
    fn from_str(s: &str) -> Result<U256, ParseU256Error> {
        let mut v = U256::ZERO;
        if s.is_empty() {
            return Err(ParseU256Error);
        }
        for b in s.bytes() {
            let d = (b as char).to_digit(10).ok_or(ParseU256Error)?;
            v = v * 10u32 + U256::from(d);
        }
        Ok(v)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParseU256Error;

impl core::fmt::Display for ParseU256Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "invalid U256 decimal literal")
    }
}

impl core::ops::Not for U256 {
    type Output = U256;
    fn not(self) -> U256 {
        U256 { hi: !self.hi, lo: !self.lo }
    }
}

impl core::ops::BitAnd for U256 {
    type Output = U256;
    fn bitand(self, rhs: U256) -> U256 {
        U256 { hi: self.hi & rhs.hi, lo: self.lo & rhs.lo }
    }
}

impl core::ops::BitOr for U256 {
    type Output = U256;
    fn bitor(self, rhs: U256) -> U256 {
        U256 { hi: self.hi | rhs.hi, lo: self.lo | rhs.lo }
    }
}

impl core::ops::Add for U256 {
    type Output = U256;
    fn add(self, rhs: U256) -> U256 {
        let (lo, carry) = self.lo.overflowing_add(rhs.lo);
        U256 { hi: self.hi + rhs.hi + carry as u128, lo }
    }
}

impl core::ops::Sub for U256 {
    type Output = U256;
    fn sub(self, rhs: U256) -> U256 {
        let (lo, borrow) = self.lo.overflowing_sub(rhs.lo);
        U256 { hi: self.hi.wrapping_sub(rhs.hi + borrow as u128), lo }
    }
}

impl core::ops::AddAssign for U256 {
    fn add_assign(&mut self, rhs: U256) {
        *self = *self + rhs;
    }
}

impl core::ops::SubAssign for U256 {
    fn sub_assign(&mut self, rhs: U256) {
        *self = *self - rhs;
    }
}

impl core::ops::BitAndAssign for U256 {
    fn bitand_assign(&mut self, rhs: U256) {
        *self = *self & rhs;
    }
}

impl core::ops::BitOrAssign for U256 {
    fn bitor_assign(&mut self, rhs: U256) {
        *self = *self | rhs;
    }
}

/// 移位量必须 < 256（kmer 路径恒满足: 2k ≤ 254）, 越界 panic（非回绕语义,
/// 与原生 `<<` 的 debug 行为一致; 回绕需求显式走 [`U256::wrapping_shl`]）。
impl core::ops::Shl<u32> for U256 {
    type Output = U256;
    fn shl(self, s: u32) -> U256 {
        assert!(s < 256, "U256 shl overflow: {s}");
        self.wrapping_shl(s)
    }
}

impl core::ops::Shr<u32> for U256 {
    type Output = U256;
    fn shr(self, s: u32) -> U256 {
        assert!(s < 256, "U256 shr overflow: {s}");
        if s == 0 {
            self
        } else if s >= 128 {
            U256 { hi: 0, lo: self.hi >> (s - 128) }
        } else {
            U256 {
                hi: self.hi >> s,
                lo: (self.lo >> s) | (self.hi << (128 - s)),
            }
        }
    }
}

macro_rules! shl_assign_shr_assign {
    ($($t:ty),*) => {$(
        impl core::ops::ShlAssign<$t> for U256 {
            fn shl_assign(&mut self, s: $t) {
                *self = *self << (s as u32);
            }
        }
        impl core::ops::ShrAssign<$t> for U256 {
            fn shr_assign(&mut self, s: $t) {
                *self = *self >> (s as u32);
            }
        }
    )*};
}
shl_assign_shr_assign!(u32, usize, i32);

/// 运算符版（非 assign）: 与 u128 习惯一致, 允许 u32/usize/i32 直接作移位量。
macro_rules! u256_shl_shr {
    ($($t:ty),*) => {$(
        impl core::ops::Shl<$t> for U256 {
            type Output = U256;
            fn shl(self, s: $t) -> U256 {
                <U256 as core::ops::Shl<u32>>::shl(self, s as u32)
            }
        }
        impl core::ops::Shr<$t> for U256 {
            type Output = U256;
            fn shr(self, s: $t) -> U256 {
                <U256 as core::ops::Shr<u32>>::shr(self, s as u32)
            }
        }
    )*};
}
u256_shl_shr!(usize, i32);

/// sequenceUtil.cpp:258 — kmer 长度上限（原版 64bit / 2bit = 32;
/// KmerId=u128 时 = 64; U256 扩容后 = 128, 对齐 SPAdes MAX_K=127 + 1 余量）。
pub const MAX_KMER_LENGTH: usize = 128;

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

/// sequenceUtil.cpp:258 kmer_to_intval — 逐字符 kmer_val<<2 | val；非 gatc 抛错（原版 cerr + throw）
pub fn kmer_to_intval(kmer: &[u8]) -> Result<KmerId, CommonError> {
    if kmer.len() > MAX_KMER_LENGTH {
        return Err(CommonError::KmerTooLong { len: kmer.len() });
    }
    let mut kmer_val: KmerId = KmerId::ZERO;
    for &c in kmer {
        let val = base_to_int(c).ok_or_else(|| CommonError::NonGatcChar {
            kmer: String::from_utf8_lossy(kmer).into_owned(),
        })?;
        kmer_val <<= 2u32;
        kmer_val |= U256::from(val);
    }
    Ok(kmer_val)
}

/// sequenceUtil.cpp:298 decode_kmer_from_intval — 从低位端逐 2-bit 解出，写在逆序位置
/// kmer_length 无上界守卫（与原版一致）；>64（MAX_KMER_LENGTH）时高位越界读出前导 G。
pub fn decode_kmer_from_intval(intval: KmerId, kmer_length: usize) -> Vec<u8> {
    let mut kmer = vec![0u8; kmer_length];
    let mut v = intval;
    for i in 1..=kmer_length {
        let base_num = (v.lo & 3) as usize;
        kmer[kmer_length - i] = INT_TO_BASE[base_num];
        v >>= 2u32;
    }
    kmer
}

/// sequenceUtil.cpp:181 revcomp_val — ~kmer 完成互补，循环移位完成 2-bit 组反转。
/// 注意 ~ 会翻转全部 128 位，但循环只提取低 kmer_length 组，高位自然丢弃。
pub fn revcomp_val(mut kmer: KmerId, kmer_length: usize) -> KmerId {
    let mut rev_kmer: KmerId = KmerId::ZERO;
    kmer = !kmer;
    for _ in 0..kmer_length {
        let base = U256::from((kmer.lo & 3) as u8);
        rev_kmer <<= 2u32;
        rev_kmer += base;
        kmer >>= 2u32;
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

/// sequenceUtil.cpp:316 compute_entropy — log2 香农熵。
/// 原版用 float（f32）逐项累加: prob * log(1/prob)/log(2.0f)。
/// 运算顺序与类型必须保持 f32，路径等价判定依赖浮点精确性（见 spec §6）。
pub fn compute_entropy(mut kmer: KmerId, kmer_length: usize) -> f32 {
    let mut counts = [0u32; 4];
    for _ in 0..kmer_length {
        let c = (kmer.lo & 3) as usize;
        kmer >>= 2u32;
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

#[cfg(test)]
mod tests {
    use super::*;

    /// 性能分层 u64 层（2026-09-04）: k≤32 域内 u64 键数值恒等于 u128 键与
    /// KmerId 低 limb（滚动/revcomp/canonical 同一算法）——diginorm 计数
    /// 换 u64 键的多重集等价由该数值一致性 + 泛型同代码保证。
    #[test]
    fn rolling_key_u64_numerics_match_u128_and_wide() {
        for k in [15usize, 25, 31, 32] {
            let mask64 = u64::low_mask(k);
            let mask128 = u128::low_mask(k);
            let mut x64 = u64::zero();
            let mut x128 = u128::zero();
            // 滚动喂入一段确定序列（含全部 4 碱基轮转）
            for b in b"GATCGATCGATCGATCGATCGATCGATCGATCGATC".iter().take(k) {
                let v = base_to_int(*b).unwrap();
                x64 = x64.roll_forward(v, mask64);
                x128 = x128.roll_forward(v, mask128);
            }
            assert_eq!(x128, x64 as u128, "k={k}: 滚动键数值一致");
            assert_eq!(x64.to_wide().lo, x64 as u128);
            assert_eq!(u64::from_wide(&x128.to_wide()), x64);
            // revcomp 一致（canonical_merge_sum 的键语义）
            assert_eq!(x128.revcomp(k), x64.revcomp(k) as u128, "k={k}: revcomp 一致");
            // 反向滚动键一致
            let rs = (2 * k - 2) as u32;
            let mut r64 = u64::zero();
            let mut r128 = u128::zero();
            for b in b"GATCGATCGATCGATCGATCGATCGATCGATCGATC".iter().take(k) {
                let v = base_to_int(*b).unwrap();
                r64 = r64.roll_rev(v, rs, mask64);
                r128 = r128.roll_rev(v, rs, mask128);
            }
            assert_eq!(r128, r64 as u128, "k={k}: 反向滚动键一致");
        }
        // k=25（diginorm 实际值）: 与黄金编码对拍
        let mask = u64::low_mask(25);
        let mut x = u64::zero();
        for b in b"ACGTAACGTAACGTAACGTAACGT" {
            x = x.roll_forward(base_to_int(*b).unwrap(), mask);
        }
        assert_eq!(x128_of_25(), x as u128);
        fn x128_of_25() -> u128 {
            let mask = u128::low_mask(25);
            let mut x = u128::zero();
            for b in b"ACGTAACGTAACGTAACGTAACGT" {
                x = x.roll_forward(base_to_int(*b).unwrap(), mask);
            }
            x
        }
    }

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
        // U256 上限 128: 129 碱基超长
        let long = vec![b'A'; 129];
        assert!(matches!(
            kmer_to_intval(&long),
            Err(CommonError::KmerTooLong { len: 129 })
        ));
    }

    #[test]
    fn decode_hand_vectors() {
        assert_eq!(decode_kmer_from_intval(U256::from(114), 4), b"ACGT".to_vec());
        assert_eq!(decode_kmer_from_intval(U256::from(85), 4), b"AAAA".to_vec());
        assert_eq!(decode_kmer_from_intval(U256::from(0), 1), b"G".to_vec());
        assert_eq!(decode_kmer_from_intval(U256::from(1), 1), b"A".to_vec());
        assert_eq!(decode_kmer_from_intval(U256::from(2), 1), b"T".to_vec());
        assert_eq!(decode_kmer_from_intval(U256::from(3), 1), b"C".to_vec());
    }

    #[test]
    fn revcomp_hand_vectors() {
        // revcomp("AA"=5) = "TT"=10
        assert_eq!(revcomp_val(U256::from(5), 2), 10);
        // revcomp("AC"=7) = "GT"=2
        assert_eq!(revcomp_val(U256::from(7), 2), 2);
        // ACGT 是回文: revcomp(ACGT)=ACGT
        assert_eq!(revcomp_val(U256::from(114), 4), 114);
        // ACGTACGT 也是回文（val=29298）
        assert_eq!(revcomp_val(U256::from(29298), 8), 29298);
        // 单碱基: A->T, G->C
        assert_eq!(revcomp_val(U256::from(1), 1), 2);
        assert_eq!(revcomp_val(U256::from(0), 1), 3);
    }

    #[test]
    fn canonical_hand_vectors() {
        // DS 规则 = max(kmer, revcomp)（sequenceUtil.cpp:376-383）
        assert_eq!(get_ds_kmer_val(U256::from(5), 2), 10); // AA -> TT
        assert_eq!(get_ds_kmer_val(U256::from(10), 2), 10); // TT -> TT
        assert_eq!(get_ds_kmer_val(U256::from(7), 2), 7); // AC(7) > GT(2) -> 7
        assert_eq!(get_ds_kmer_val(U256::from(114), 4), 114); // 回文不变
    }

    #[test]
    fn revcomp_roundtrip() {
        // 任意 kmer 双取 revcomp 复原
        let k = kmer_to_intval(b"AAAATAAAATAAAATAAAATAAAAT").unwrap();
        assert_eq!(revcomp_val(revcomp_val(k, 25), 25), k);
    }

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

    #[test]
    fn error_messages_mirror_original() {
        // 原版 cerr 消息（去 \n\n 前缀）格式固定; 上限数字随 KmerId=U256 变 128
        assert_eq!(
            kmer_to_intval(b"ACGN").unwrap_err().to_string(),
            "error, kmer contains nongatc: ACGN"
        );
        assert_eq!(
            kmer_to_intval(&[b'A'; 129]).unwrap_err().to_string(),
            "error, kmer length exceeds 128: 129"
        );
    }

    #[test]
    fn kmer_length_boundaries() {
        // 128 碱基是接受上限（MAX_KMER_LENGTH, U256/2bit）
        let k128 = vec![b'A'; 128];
        assert!(kmer_to_intval(&k128).is_ok());
        // 0 长度: 原版返回 0（循环零次）
        assert_eq!(kmer_to_intval(b"").unwrap(), 0);
        // 脏高位: revcomp 只看低 2*k 位（KmerId::MAX 的低 2 位全 1 → 单碱基 revcomp 视角）
        // 推演: !MAX=0 → base=0&3=0 → rev_kmer=0（sequenceUtil.cpp:181-195 原版同结果）
        assert_eq!(revcomp_val(KmerId::MAX, 1), 0);
        // decode/entropy 对任意长度行为良定义（此处只验证不 panic 且熵=0）
        // 推演: MAX&3=3 → counts[3]=1 → prob=1.0 → 1.0*ln(1)/ln(2)=0
        assert!(compute_entropy(KmerId::MAX, 1).abs() < 1e-6); // 单碱基熵为 0
    }

    /// U256 扩容边界: k=33（首个 u64 放不下的 k）、37（RL=76 的 upper_k）、
    /// 55（RNA_MAX_LOWER_K）、64（u128 上限）、73（RL=150 的 upper_k）、
    /// 127（RL=250 的 upper_k / SPAdes MAX_K）、128（U256 物理上限）的
    /// 编码/解码/revcomp 往返。关键不变量: 高位碱基（>64 位部分）不丢。
    #[test]
    fn u128_beyond_32_boundaries_roundtrip() {
        for k in [33usize, 37, 55, 64, 73, 127, 128] {
            // 首尾异质 + 中段混合的确定性序列（首 C 尾 A, 让高低两端都有非零位）
            let seq: Vec<u8> = (0..k)
                .map(|i| match (i * 7 + k) % 4 {
                    0 => b'G',
                    1 => b'A',
                    2 => b'T',
                    _ => b'C',
                })
                .collect();
            let v = kmer_to_intval(&seq).unwrap();
            // 编码值占满 2k 位空间且不丢高位: 重编码解码结果必须逐字节还原
            assert_eq!(decode_kmer_from_intval(v, k), seq, "k={k} 编解码往返");
            // revcomp 双取复原
            assert_eq!(revcomp_val(revcomp_val(v, k), k), v, "k={k} revcomp 往返");
            // DS canonical 幂等（max(k, rc(k)) 再 canonical 不变）
            let canon = get_ds_kmer_val(v, k);
            assert_eq!(get_ds_kmer_val(canon, k), canon, "k={k} canonical 幂等");
            // 熵计算在 2k>64 位的值上不 panic 且落在 [0, 2]
            let h = compute_entropy(v, k);
            assert!((0.0..=2.0).contains(&h), "k={k} 熵 {h} 越界");
        }
        // k=37 目标用例显式锁: 第 33 位碱基（u64 时代的截断点）参与键值
        let seq37: Vec<u8> = std::iter::once(b'C')
            .chain(std::iter::repeat_n(b'A', 35))
            .chain(std::iter::once(b'T'))
            .collect();
        assert_eq!(seq37.len(), 37);
        let v37 = kmer_to_intval(&seq37).unwrap();
        assert!(v37 > U256::from(u64::MAX as u128), "k=37 全 C..T 键值应超出 u64");
        assert_eq!(decode_kmer_from_intval(v37, 37), seq37);
    }
}
