//! nbkc_normalize.pl（L67-127）+ nbkc_merge_left_right_stats.pl（L89-148）镜像:
//! DigiNorm 的概率选择与 PE stats 合并。
//!
//! **随机数消耗序是核心不变量**: 原版 rand(1) 只对**通过全部过滤**的行调用
//! （三个 continue 都在 rand 之前），被过滤的行不消耗随机数——整个序列的
//! 对齐由此决定，select 的测试专门锁定这一点。
//!
//! NaN 语义: perl 数值比较对 NaN 恒 false，Rust 相同——NaN median 通过
//! min_cov/mean 门、NaN cv 通过 max_CV 门、`rand <= NaN` 为 false（该行被丢），
//! 全部自然镜像，无需特判。

use crate::drand48::Drand48;

/// nbkc_normalize.pl 参数（--max_cov/--min_cov/--max_CV，均数值比较语义）。
pub struct NbkcParams {
    pub max_cov: f64,
    pub min_cov: f64,
    pub max_cv: f64,
}

/// 一行 stats（acc/median_cov/mean_cov/stdev，已从 TSV 文本解析为数值）。
#[derive(Debug, Clone, PartialEq)]
pub struct StatsRow {
    pub acc: String,
    pub median: f64,
    pub mean: f64,
    pub stdev: f64,
}

/// 原版 `s|/[12]$||`: 仅剥结尾 "/1" 或 "/2"（"/3"、"/12" 等不动）。
pub fn core_acc(acc: &str) -> String {
    acc.strip_suffix("/1")
        .or_else(|| acc.strip_suffix("/2"))
        .unwrap_or(acc)
        .to_string()
}

/// nbkc_normalize.pl L81-115: 过滤 + 概率保留。
/// 返回保留的 core acc（顺序 = 输入顺序，镜像原版打印序）。
///
/// - `median < min_cov` → 丢（below min coverage）
/// - `mean <= 0` → 丢（aberrant）
/// - `stdev/mean > max_cv` → 丢（aberrant; NaN cv 比较为 false → 通过）
/// - `rand(1) <= max_cov/median` → 留（rand 只在此消耗; rand ∈ [0,1)，
///   故 `max_cov/median >= 1` 时恒留）
pub fn select(rows: &[StatsRow], p: &NbkcParams, rng: &mut Drand48) -> Vec<String> {
    let mut out = Vec::new();
    for r in rows {
        if r.median < p.min_cov {
            continue;
        }
        if r.mean <= 0.0 {
            continue;
        }
        let cv = r.stdev / r.mean;
        if cv > p.max_cv {
            continue;
        }
        if rng.next_f64() <= p.max_cov / r.median {
            out.push(core_acc(&r.acc));
        }
    }
    out
}

/// 合并 stats 行（原版打印 12 列，此处保留合并后的 4 列:
/// core acc + 三项 `sprintf("%.1f", (l+r)/2)` 合成指标）。
#[derive(Debug, Clone, PartialEq)]
pub struct MergedRow {
    pub core: String,
    pub median: String,
    pub mean: String,
    pub stdev: String,
}

/// nbkc_merge_left_right_stats.pl L89-148（--sorted 路径）: 双指针按 core acc
/// 字节序合并两个已排序 stats; 未配对侧静默跳过（advance 指针，不输出）。
///
/// 与原版的两处实现差异（对 DigiNorm 正常输入等价，记录备查）:
/// 1. 原版 core 提取是 `^(\S+)/\d$`（任意一位数字），此处复用 core_acc（仅 /1、/2）;
///    stats 侧 acc 只会是 /1、/2 后缀。
/// 2. 原版推进比较用完整 acc（`lt`），此处比较 core——两侧后缀统一时序一致。
/// 3. 平均对**已解析的数值**做（原版对 TSV 文本数值化后求平均，同一组数值）。
pub fn merge_pairs(left: &[StatsRow], right: &[StatsRow]) -> Vec<MergedRow> {
    let mut li = 0;
    let mut ri = 0;
    let mut out = Vec::new();
    while li < left.len() && ri < right.len() {
        let lc = core_acc(&left[li].acc);
        let rc = core_acc(&right[ri].acc);
        match lc.as_bytes().cmp(rc.as_bytes()) {
            std::cmp::Ordering::Less => {
                li += 1;
            }
            std::cmp::Ordering::Greater => {
                ri += 1;
            }
            std::cmp::Ordering::Equal => {
                out.push(MergedRow {
                    core: lc,
                    median: avg_1f(left[li].median, right[ri].median),
                    mean: avg_1f(left[li].mean, right[ri].mean),
                    stdev: avg_1f(left[li].stdev, right[ri].stdev),
                });
                li += 1;
                ri += 1;
            }
        }
    }
    out
}

/// `sprintf("%.1f", (l+r)/2)` 镜像: glibc 与 Rust 均对**精确二进制值**做
/// round-half-even（黄金 fixtures/p1/avg1f_golden.txt 逐值锁定）。
/// NaN → "NaN"（perl 实测）。注意 +inf 时 perl 打 "Inf"、Rust 为 "inf"
/// （有限覆盖统计不可能出现，未做映射）。
fn avg_1f(l: f64, r: f64) -> String {
    format!("{:.1}", (l + r) / 2.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(acc: &str, median: f64, mean: f64, stdev: f64) -> StatsRow {
        StatsRow {
            acc: acc.into(),
            median,
            mean,
            stdev,
        }
    }

    fn params(max_cov: f64, min_cov: f64, max_cv: f64) -> NbkcParams {
        NbkcParams {
            max_cov,
            min_cov,
            max_cv,
        }
    }

    #[test]
    fn core_acc_variants() {
        assert_eq!(core_acc("r1/1"), "r1");
        assert_eq!(core_acc("r1/2"), "r1");
        assert_eq!(core_acc("r1"), "r1");
        assert_eq!(core_acc("r1/3"), "r1/3"); // 原版只剥 [12]
        assert_eq!(core_acc("r1/12"), "r1/12"); // 结尾是 "12" 不是 "/1"
        assert_eq!(core_acc("r1/"), "r1/");
        assert_eq!(core_acc("a/1/2"), "a/1"); // 只剥一层
        assert_eq!(core_acc("/1"), "");
    }

    #[test]
    fn select_three_filters_each_isolate() {
        // ratio >= 1 保证任何通过全部过滤的行必被保留 → 输出为空即对应行被该门拦下
        let p = params(10.0, 1.0, 5.0);
        let mut rng = Drand48::new(12345);
        let rows = [
            row("below/1", 0.0, 5.0, 0.0),  // median < min_cov
            row("nomean/1", 5.0, 0.0, 0.0), // mean <= 0
            row("negmean/1", 5.0, -1.0, 0.0),
            row("highcv/1", 5.0, 1.0, 10.0), // cv = 10 > 5
        ];
        assert!(select(&rows, &p, &mut rng).is_empty());
        // 各门恰好卡边界时不拦: median == min_cov（< 而非 <=）、cv == max_cv（> 而非 >=）
        let mut rng = Drand48::new(12345);
        let rows = [
            row("edge_med/1", 1.0, 2.0, 10.0), // cv = 5 == max_cv → 通过
            row("good/1", 2.0, 2.0, 0.0),
        ];
        assert_eq!(select(&rows, &p, &mut rng), vec!["edge_med", "good"]);
    }

    #[test]
    fn select_nan_semantics() {
        let p = params(10.0, 1.0, 5.0);
        // NaN stdev → cv NaN → `NaN > max_cv` 为 false → 通过 CV 门（镜像 perl）
        let mut rng = Drand48::new(12345);
        let rows = [row("nan_sd/1", 5.0, 2.0, f64::NAN)];
        assert_eq!(select(&rows, &p, &mut rng), vec!["nan_sd"]);
        // NaN mean: 通过 mean 门（NaN <= 0 为 false）与 CV 门（cv = NaN）;
        // 概率门的比值用 **median**（10/5 = 2 >= 1）→ 仍被保留（perl 同）
        let mut rng = Drand48::new(12345);
        let rows = [row("nan_mean/1", 5.0, f64::NAN, 1.0)];
        assert_eq!(select(&rows, &p, &mut rng), vec!["nan_mean"]);
        // NaN median: 通过 min_cov 门; ratio = NaN → 不保留
        let mut rng = Drand48::new(12345);
        let rows = [row("nan_med/1", f64::NAN, 2.0, 0.0)];
        assert!(select(&rows, &p, &mut rng).is_empty());
    }

    #[test]
    fn select_ratio_ge1_always_kept() {
        // max_cov/median >= 1（rand ∈ [0,1)）→ 恒保留，与具体 rand 值无关
        let p = params(200.0, 1.0, 1000.0);
        let mut rng = Drand48::new(12345);
        let rows = [
            row("a/1", 200.0, 2.0, 1.0), // ratio = 1.0
            row("b/2", 5.0, 2.0, 1.0),   // ratio = 40
            row("c", 1.0, 2.0, 1.0),
        ];
        assert_eq!(select(&rows, &p, &mut rng), vec!["a", "b", "c"]);
    }

    /// **随机数消耗序锁定**: rand 只对通过全部过滤的行消耗。
    /// seed 12345 序列: r0=0.2253, r1=0.9192, r2=0.2068, r3=0.7248, r4=0.7322。
    /// median=20 → ratio = 10/20 = 0.5: r0、r2 <= 0.5 留，r1、r3、r4 > 0.5 丢。
    /// 若（错误地）每行都消耗 rand: B 用 r1、D 用 r3、E 用 r4 → 全丢 → 输出空，
    /// 与本断言区分。
    #[test]
    fn select_rand_consumed_only_for_passing_rows() {
        let p = params(10.0, 1.0, 100.0);
        let rows = [
            row("A/1", 0.0, 5.0, 0.0),  // 过滤断（min_cov）→ 不消耗 rand
            row("B/1", 20.0, 5.0, 0.0), // 用 r0=0.2253 <= 0.5 → 留
            row("C/1", 5.0, 0.0, 0.0),  // 过滤断（mean）→ 不消耗 rand
            row("D/1", 20.0, 5.0, 0.0), // 用 r1=0.9192 > 0.5 → 丢
            row("E/1", 20.0, 5.0, 0.0), // 用 r2=0.2068 <= 0.5 → 留
        ];
        let mut rng = Drand48::new(12345);
        assert_eq!(select(&rows, &p, &mut rng), vec!["B", "E"]);
    }

    #[test]
    fn select_output_order_follows_input() {
        let p = params(10.0, 1.0, 100.0);
        let rows = [
            row("z/1", 5.0, 5.0, 0.0),
            row("a/2", 5.0, 5.0, 0.0),
            row("m/1", 5.0, 5.0, 0.0),
        ];
        let mut rng = Drand48::new(12345);
        // 原版按文件行序打印，不做排序
        assert_eq!(select(&rows, &p, &mut rng), vec!["z", "a", "m"]);
    }

    #[test]
    fn merge_pairs_basic_and_unpaired_skip() {
        let left = [
            row("a/1", 10.0, 5.0, 0.5),
            row("b/1", 20.0, 6.0, 0.0), // 仅左侧 → 静默跳过
            row("d/1", 1.0, 1.0, 1.0),
        ];
        let right = [
            row("a/2", 11.0, 7.0, 0.5),
            row("c/2", 9.0, 9.0, 9.0), // 仅右侧 → 静默跳过
            row("d/2", 2.0, 2.0, 2.0),
        ];
        let merged = merge_pairs(&left, &right);
        assert_eq!(merged.len(), 2);
        assert_eq!(merged[0].core, "a");
        assert_eq!(merged[0].median, "10.5");
        assert_eq!(merged[0].mean, "6.0");
        assert_eq!(merged[0].stdev, "0.5");
        assert_eq!(merged[1].core, "d");
        assert_eq!(merged[1].median, "1.5");
        assert_eq!(merged[1].mean, "1.5");
        assert_eq!(merged[1].stdev, "1.5");
    }

    #[test]
    fn merge_pairs_byte_order_and_edges() {
        // 字节序: 'Z'(0x5A) < 'a'(0x61) → Z 先配对（perl lt 按字节比较）
        let left = [row("Z/1", 1.0, 1.0, 0.0), row("a/1", 2.0, 2.0, 0.0)];
        let right = [row("Z/2", 3.0, 1.0, 0.0), row("a/2", 4.0, 2.0, 0.0)];
        let merged = merge_pairs(&left, &right);
        assert_eq!(merged[0].core, "Z");
        assert_eq!(merged[1].core, "a");
        // 无后缀 acc 与带后缀 acc 同 core → 可配对（merge 脚本 \d 正则同此场景）
        let left = [row("x", 1.0, 1.0, 0.0)];
        let right = [row("x/2", 3.0, 1.0, 0.0)];
        assert_eq!(merge_pairs(&left, &right)[0].median, "2.0");
        // 一侧耗尽即止（剩余不输出）
        let left = [row("a/1", 1.0, 1.0, 0.0), row("b/1", 1.0, 1.0, 0.0)];
        let right = [row("a/2", 1.0, 1.0, 0.0)];
        assert_eq!(merge_pairs(&left, &right).len(), 1);
        assert!(merge_pairs(&[], &[]).is_empty());
    }

    /// avg_1f 黄金锁定: fixtures/p1/avg1f_golden.txt 由 perl
    /// `printf "%.1f\n", (l+r)/2` 生成，输入对与原命令一致。
    #[test]
    fn avg1f_golden_vs_perl() {
        let pairs = [
            (1.25, 1.35),
            (2.05, 2.15),
            (0.5, 0.5),
            (100.25, 100.35),
            (1.05, 1.1),
            (3.3, 3.4),
            (12.345, 12.355),
            (0.05, 0.15),
            (1.0, 2.0),
            (999.99, 999.99),
        ];
        let text = std::fs::read_to_string(format!(
            "{}/../../fixtures/p1/avg1f_golden.txt",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), pairs.len());
        for ((l, r), expect) in pairs.iter().zip(&lines) {
            assert_eq!(&avg_1f(*l, *r), expect, "input=({l}, {r})");
        }
    }

    /// perl 实测补充值（含精确二进制半值 tie 与 NaN）:
    /// 0.25→"0.2"、2.25→"2.2"（half-even 偶侧）; -0.25→"-0.2"; -0.15→"-0.1";
    /// (0.1+0.2)/2 = 0.15000000000000002 → "0.2"; NaN → "NaN"。
    #[test]
    fn avg1f_ties_and_nan() {
        assert_eq!(avg_1f(0.25, 0.25), "0.2");
        assert_eq!(avg_1f(2.25, 2.25), "2.2");
        assert_eq!(avg_1f(0.75, 0.75), "0.8");
        assert_eq!(avg_1f(-0.25, -0.25), "-0.2");
        assert_eq!(avg_1f(-0.15, -0.15), "-0.1");
        assert_eq!(avg_1f(0.1, 0.2), "0.2");
        assert_eq!(avg_1f(0.35, 0.35), "0.3");
        assert_eq!(avg_1f(0.45, 0.45), "0.5");
        assert_eq!(avg_1f(1.5, 2.0), "1.8");
        // NaN 传染: perl `printf "%.1f", (nan+2)/2` → "NaN"（Rust 同形）
        assert_eq!(avg_1f(f64::NAN, 2.0), "NaN");
        assert_eq!(avg_1f(2.0, f64::NAN), "NaN");
    }

    /// **原版脚本实测对拍**（trinityrnaseq-v2.15.2 nbkc_normalize.pl，perl 5.38.2，
    /// --max_cov 10 --min_cov 1 --max_CV 100）: stdout 名单逐字为 B E F G。
    /// 序列 r0..r5 = 0.2253, 0.9192, 0.2068, 0.7248, 0.7322, 0.9065:
    /// A/C 是过滤行不消耗 rand; B(r0≤0.5) E(r2≤0.5) 留、D(r1>0.5) H(r5>0.5) 丢;
    /// F/G ratio=10/5=2≥1 恒留（G 的 stdev=nan 通过 CV 门）。
    #[test]
    fn select_matches_original_script_output() {
        let rows = [
            row("A/1", 0.0, 5.0, 0.0),
            row("B/1", 20.0, 5.0, 0.0),
            row("C/1", 5.0, 0.0, 0.0),
            row("D/1", 20.0, 5.0, 0.0),
            row("E/1", 20.0, 5.0, 0.0),
            row("F/1", 5.0, 1.0, 10.0),
            row("G/1", 5.0, 2.0, f64::NAN),
            row("H/2", 20.0, 5.0, 0.0),
        ];
        let mut rng = Drand48::new(12345);
        let selected = select(&rows, &params(10.0, 1.0, 100.0), &mut rng);
        assert_eq!(selected, vec!["B", "E", "F", "G"]);
    }

    /// **原版脚本实测对拍**（nbkc_merge_left_right_stats.pl --sorted）:
    /// a/c/d 三对合并，b（仅左）、e（仅右）静默跳过，合成列逐字一致。
    #[test]
    fn merge_matches_original_script_output() {
        let left = [
            row("a/1", 10.0, 5.0, 0.5),
            row("b/1", 20.0, 6.0, 0.0),
            row("c/2", 9.0, 9.0, 9.0), // 两侧同 acc（core 相等）→ 仍配对，perl 同
            row("d/1", 1.0, 1.0, 1.0),
        ];
        let right = [
            row("a/2", 11.0, 7.0, 0.5),
            row("c/2", 9.0, 9.0, 9.0),
            row("d/2", 2.0, 2.0, 2.0),
            row("e/2", 3.0, 3.0, 3.0),
        ];
        let merged = merge_pairs(&left, &right);
        let rendered: Vec<String> = merged
            .iter()
            .map(|m| format!("{}\t{}\t{}\t{}", m.core, m.median, m.mean, m.stdev))
            .collect();
        assert_eq!(
            rendered,
            vec!["a\t10.5\t6.0\t0.5", "c\t9.0\t9.0\t9.0", "d\t1.5\t1.5\t1.5"]
        );
        // NaN median 传染: 原版该行合成 median 列 = "NaN"
        let merged = merge_pairs(
            &[row("x/1", f64::NAN, 1.0, 1.0)],
            &[row("x/2", 2.0, 2.0, 2.0)],
        );
        assert_eq!(merged[0].median, "NaN");
    }
}
