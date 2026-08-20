//! jellyfish dump 兼容输出: ">count\nKMER" FASTA（大写 ACGT，无折行）。
//! 遍历序为哈希序（与 jellyfish dump 一样不确定）——等价性由多重集比较保证。
//!
//! DS（canonical）代表选择: 计数键 = max(kmer, revcomp)（Trinity 编码 G=0,A=1,T=2,C=3 序，
//! sequenceUtil.cpp:376），而 jellyfish canonical = 其编码（A=0,C=1,G=2,T=3，首碱基最高位）
//! 下的整数值较小者 ⇔ 词典序较小串。两种代表各覆盖 DS 类一次，喂原版 inchworm 等价
//! （inchworm 读入后自行 get_DS_kmer_val，两个代表映到同一内部键）; 但要与 jellyfish dump
//! 逐字节对拍，DS 模式必须输出词典序较小串（实测 50 条 read 样本: 1007/1887 类代表不同，
//! 换代表后完全一致; SS 模式两方本就逐字节一致）。

use std::io::Write;

use trinity_common::kmer::{decode_kmer_from_intval, revcomp_val};

use crate::counter::CountMap;

pub fn write_dump<W: Write>(
    w: &mut W,
    counts: &CountMap,
    k: usize,
    min_count: u32,
    ds: bool,
) -> std::io::Result<()> {
    for (&key, &c) in counts.iter() {
        if c < min_count {
            continue;
        }
        let mut seq = decode_kmer_from_intval(key, k);
        if ds {
            // jellyfish 代表 = 类内词典序较小串（DS 键的 revcomp 即另一个代表）
            let other = decode_kmer_from_intval(revcomp_val(key, k), k);
            if other < seq {
                seq = other;
            }
        }
        writeln!(w, ">{c}")?;
        w.write_all(&seq)?;
        writeln!(w)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustc_hash::FxHashMap;
    use trinity_common::kmer::kmer_to_intval;

    #[test]
    fn dump_format_and_min_count_filter() {
        let mut m: CountMap = FxHashMap::default();
        m.insert(kmer_to_intval(b"ACGT").unwrap(), 3);
        m.insert(kmer_to_intval(b"TTTT").unwrap(), 1);
        let mut buf = Vec::new();
        write_dump(&mut buf, &m, 4, 1, false).unwrap();
        let lines: Vec<&str> = std::str::from_utf8(&buf).unwrap().lines().collect();
        assert!(lines.contains(&">3"));
        assert!(lines.contains(&"ACGT"));
        assert!(lines.contains(&">1"));
        assert!(lines.contains(&"TTTT"));
        let mut buf2 = Vec::new();
        write_dump(&mut buf2, &m, 4, 2, false).unwrap();
        let s = String::from_utf8(buf2).unwrap();
        assert!(!s.contains("TTTT"));
        assert!(s.contains("ACGT"));
    }

    #[test]
    fn ds_dump_emits_jellyfish_lex_min_representative() {
        // DS 键 = max(kmer, revcomp)（Trinity 编码序）; jellyfish canonical = 词典序较小串。
        // {AAAA, TTTT} 类: Trinity 键是 TTTT（enc T=2 > A=1），jellyfish 输出 AAAA。
        let mut m: CountMap = FxHashMap::default();
        m.insert(kmer_to_intval(b"TTTT").unwrap(), 5);
        let mut buf = Vec::new();
        write_dump(&mut buf, &m, 4, 1, true).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("AAAA"));
        assert!(!s.contains("TTTT"));
        // SS 模式不作代表置换: 键原样输出
        let mut buf2 = Vec::new();
        write_dump(&mut buf2, &m, 4, 1, false).unwrap();
        let s2 = String::from_utf8(buf2).unwrap();
        assert!(s2.contains("TTTT"));
        assert!(!s2.contains("AAAA"));
        // 回文类两个代表相同，不受影响: ACGT revcomp = ACGT
        let mut mp: CountMap = FxHashMap::default();
        mp.insert(kmer_to_intval(b"ACGT").unwrap(), 2);
        let mut buf3 = Vec::new();
        write_dump(&mut buf3, &mp, 4, 1, true).unwrap();
        assert!(String::from_utf8(buf3).unwrap().contains("ACGT"));
    }
}
