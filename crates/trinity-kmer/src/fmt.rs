//! C++ `std::ostream <<` 默认浮点格式（≈ printf %g，6 位有效数字，去尾零）。
//! fastaToKmerCoverageStats 的 stats 输出依赖此格式（含 "nan"/"-nan"/"-0" 特例）。
//!
//! 实现要点: 先用 Rust `{:.*e}`（与 glibc 同为对精确二进制值的 correct rounding）取
//! 6 位有效数字的科学计数形式并拿到**舍入后**的十进制指数，再按 %g 规则选型:
//! 指数 < -4 或 >= 6 → 科学计数（指数 C 形态 `e+05`，至少 2 位）; 否则定点
//! `{:.(5-exp)}`（对原值格式化——舍入进位跨数量级时 5-exp 与 %g 一致）。
//! 最后去尾零与孤立小数点、补负号。

/// C++ ostream 默认输出（defaultfloat + precision 6，等价 printf %g）。
pub fn format_g6(v: f32) -> String {
    if v.is_nan() {
        return if v.is_sign_negative() {
            "-nan".into()
        } else {
            "nan".into()
        };
    }
    if v.is_infinite() {
        return if v < 0.0 { "-inf".into() } else { "inf".into() };
    }
    if v == 0.0 {
        return if v.is_sign_negative() {
            "-0".into()
        } else {
            "0".into()
        };
    }
    let neg = v < 0.0;
    let a = v.abs() as f64; // f32→f64 精确，舍入语义不变
    let p = 6usize;
    // 科学计数取 p 位有效数字，读出舍入后的指数（进位跨数量级由 Rust 格式化处理）
    let sci = format!("{:.*e}", p - 1, a);
    let (mant_str, exp_str) = sci.split_once('e').unwrap();
    let exp: i32 = exp_str.parse().unwrap();
    let mut body = if exp < -4 || exp >= p as i32 {
        // 科学计数: 尾数去尾零; 指数 C 形态（符号 + 至少 2 位，Rust {:e} 是 "6" 不是 "+06"）
        let mut m = mant_str.to_string();
        strip_trailing_zeros(&mut m);
        format!("{m}e{}{:02}", if exp < 0 { '-' } else { '+' }, exp.abs())
    } else {
        // 定点: 对原值 a 格式化（不是对尾数——尾数已丢失数量级），小数位 = p-1-exp
        let decimals = (p as i32 - 1 - exp) as usize;
        let mut m = format!("{:.*}", decimals, a);
        strip_trailing_zeros(&mut m);
        m
    };
    if neg {
        body.insert(0, '-');
    }
    body
}

/// 去掉小数部分尾零与孤立的小数点（"1.00000"→"1"，"0.810000"→"0.81"）。
fn strip_trailing_zeros(s: &mut String) {
    if s.contains('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g6_basics() {
        assert_eq!(format_g6(0.0), "0");
        assert_eq!(format_g6(-0.0), "-0");
        assert_eq!(format_g6(2.0), "2");
        assert_eq!(format_g6(0.81127816f32), "0.811278");
        assert_eq!(format_g6(1234567.0f32), "1.23457e+06");
        assert_eq!(format_g6(0.00001234567f32), "1.23457e-05");
        assert_eq!(format_g6(123.456f32), "123.456");
        assert_eq!(format_g6(1.5f32), "1.5");
    }

    #[test]
    fn g6_specials() {
        // x86 SSE 默认 QNaN 符号位为 1 → 原版输出 "-nan"（edge.stats.orig.tsv 实证）
        assert_eq!(format_g6(f32::NAN), "nan");
        assert_eq!(format_g6(-f32::NAN), "-nan");
        assert_eq!(format_g6(f32::INFINITY), "inf");
        assert_eq!(format_g6(f32::NEG_INFINITY), "-inf");
    }

    #[test]
    fn g6_golden_vs_printf() {
        let tsv = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../fixtures/p1/g6_golden.tsv"
        ))
        .unwrap();
        let mut n = 0;
        for line in tsv.lines() {
            let (input, expect) = line.split_once('\t').unwrap();
            let v: f32 = input.parse().unwrap();
            assert_eq!(format_g6(v), expect, "input={input}");
            n += 1;
        }
        assert!(n >= 40);
    }
}
