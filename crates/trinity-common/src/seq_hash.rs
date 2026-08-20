//! sequenceUtil.cpp:422 generateHash — sdbm 型 32 位哈希 + base_val 求和的 64 位组合。
//!
//! 原版（逐行核对）:
//! - `hash` 为 unsigned int: `hash = 65599 * hash + nucleotide`（32 位回绕）
//! - `base_val = base_to_int_value(nucleotide) + 1` → G=1, A=2, T=3, C=4，
//!   非 gatc = 0（-1 + 1；原版注释称不应遇到非 gatc，但表确实如此）
//! - `combined_hashcode += base_val`（64 位累加，无回绕之虞——序列长度受限）
//! - fold: `hash = hash ^ (hash >> 16)`（unsigned int 逻辑右移）
//! - 返回 `(combined << 32) | hash`
//!
//! 用途: inchworm 输出去重 key（接收方再截 `as u32`，复现原版 unsigned int 接收）。
//! 注: nucleotide 是 C 的 char——字节 ≥ 0x80 时符号扩展；本函数面向 ASCII 序列
//! （与原版调用场景一致），`n as u32` 与之等价。

/// sequenceUtil.cpp:422-445 generateHash
pub fn generate_hash(s: &[u8]) -> u64 {
    let mut hash: u32 = 0;
    let mut combined: u64 = 0;
    for &n in s {
        hash = hash.wrapping_mul(65599).wrapping_add(n as u32);
        let base_val: u64 = match n {
            b'G' | b'g' => 1,
            b'A' | b'a' => 2,
            b'T' | b't' => 3,
            b'C' | b'c' => 4,
            _ => 0, // 非 gatc: base_to_int_value = -1，+1 后为 0
        };
        combined += base_val;
    }
    hash ^= hash >> 16; // unsigned int 逻辑右移
    (combined << 32) | (hash as u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hand_vectors_single_base() {
        // "G": hash = 65599*0 + 71 = 71; 71>>16=0 → fold 不变; base_val=1
        assert_eq!(generate_hash(b"G"), (1u64 << 32) | 71);
        // "A": hash = 65; base_val=2
        assert_eq!(generate_hash(b"A"), (2u64 << 32) | 65);
        // "T": hash = 84（ASCII 'T'，非小写 't'=116）; base_val=3
        assert_eq!(generate_hash(b"T"), (3u64 << 32) | 84);
        // "C": hash = 67; base_val=4
        assert_eq!(generate_hash(b"C"), (4u64 << 32) | 67);
    }

    #[test]
    fn empty_string_is_zero() {
        // 原版空串: hash=0, combined=0 → 0
        assert_eq!(generate_hash(b""), 0);
    }

    #[test]
    fn non_gatc_contributes_hash_but_not_base_val() {
        // "N": hash = 78, base_val = 0（-1+1）
        assert_eq!(generate_hash(b"N"), 78u64);
        // "NNNNN": hash = 65599*(...)*5 累加 5 个 78; base_val 全 0 → 高 32 位为 0
        let mut h: u32 = 0;
        for _ in 0..5 {
            h = 65599u32.wrapping_mul(h).wrapping_add(78);
        }
        let folded = h ^ (h >> 16);
        assert_eq!(generate_hash(b"NNNNN"), folded as u64);
    }

    #[test]
    fn case_insensitive_base_val_but_case_sensitive_hash() {
        // 小写进 base_val（G=1 同大写），但 sdbm 哈希吃原始字节 → 低 32 位不同
        let upper = generate_hash(b"GATC");
        let lower = generate_hash(b"gatc");
        assert_eq!(upper >> 32, lower >> 32); // base_val 和相同（1+2+3+4=10）
        assert_ne!(upper & 0xFFFF_FFFF, lower & 0xFFFF_FFFF);
    }
}
