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
/// kmer_length 无上界守卫（与原版一致）；>32 时高位越界读出前导 G。
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
        // 原版 cerr 消息（去 \n\n 前缀）格式固定
        assert_eq!(
            kmer_to_intval(b"ACGN").unwrap_err().to_string(),
            "error, kmer contains nongatc: ACGN"
        );
        assert_eq!(
            kmer_to_intval(&[b'A'; 33]).unwrap_err().to_string(),
            "error, kmer length exceeds 32: 33"
        );
    }

    #[test]
    fn kmer_length_boundaries() {
        // 32 碱基是接受上限（MAX_KMER_LENGTH）
        let k32 = vec![b'A'; 32];
        assert!(kmer_to_intval(&k32).is_ok());
        // 0 长度: 原版返回 0（循环零次）
        assert_eq!(kmer_to_intval(b"").unwrap(), 0);
        // 脏高位: revcomp 只看低 2*k 位（u64::MAX 的低 4 位全 1 → 单碱基 revcomp 视角）
        // 推演: !MAX=0 → base=0&3=0 → rev_kmer=0（sequenceUtil.cpp:181-195 原版同结果）
        assert_eq!(revcomp_val(u64::MAX, 1), 0);
        // decode/entropy 对任意长度行为良定义（此处只验证不 panic 且熵=0）
        // 推演: MAX&3=3 → counts[3]=1 → prob=1.0 → 1.0*ln(1)/ln(2)=0
        assert!(compute_entropy(u64::MAX, 1).abs() < 1e-6); // 单碱基熵为 0
    }
}
