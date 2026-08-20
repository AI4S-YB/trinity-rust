//! 各 bin 共用的极简 argv 解析器（镜像原版 ArgProcessor/CommandLineParser 语义：
//! 未知参数报错）。
//!
//! 消费式接口: `Cli::new(argv)` 之后按名取值——各 `*_flag` 方法把匹配的 token
//! 从待处理列表移除，最后 [`Cli::finish`] 检查无残留（残留 → "do not understand
//! option"，原版 ArgProcessor 文案）。约定 `finish()` 在任何副作用（建输出文件、
//! 跑计算）之前调用，保持既有「参数非法时不留半成品」的行为。
//!
//! 重复给出同一值 flag: 后者覆盖前者（镜像原逐 token 赋值的末个生效语义），
//! 全部出现一并消费、不视为残留。值 flag 缺值在消费时报 `{name} expects a value`。
//!
//! 错误文案与既有两 bin（trinity-kmer/inchworm）逐字一致，迁移零文案变化;
//! 类型化解析失败统一为 `invalid {name} value: {v}`（`-K` 等带专属文案的由
//! bin 在 [`Cli::value_flag`] 之上自行包装）。

use std::collections::HashSet;
use std::str::FromStr;

use crate::error::CommonError;

/// usage=true → 参数问题（各 bin 约定 exit 2）; false → 运行失败（exit 1）。
#[derive(Debug)]
pub struct CliError {
    pub msg: String,
    pub usage: bool,
}

impl CliError {
    pub fn usage(msg: impl Into<String>) -> Self {
        CliError {
            msg: msg.into(),
            usage: true,
        }
    }
    pub fn run(msg: impl std::fmt::Display) -> Self {
        CliError {
            msg: msg.to_string(),
            usage: false,
        }
    }
}

/// 库层错误（CommonError）一律归入运行类（exit 1）。
impl From<CommonError> for CliError {
    fn from(e: CommonError) -> Self {
        CliError::run(e)
    }
}

/// 极简 argv 解析器（消费式，无依赖）。
#[derive(Debug)]
pub struct Cli {
    /// 尚未消费的 token。
    args: Vec<String>,
    /// 已消费到值的 flag 名（[`Cli::was_given`] 查询——"set to: N" 类回显需区分
    /// 「显式给出恰好等于默认值」与「未给出」）。
    given: HashSet<String>,
}

impl Cli {
    pub fn new(args: &[String]) -> Self {
        Cli {
            args: args.to_vec(),
            given: HashSet::new(),
        }
    }

    /// 消费 `--name value` 的**全部**出现并返回值（重复末个生效）。缺值报错。
    /// 长旗标同时支持 `--name=value` 等号形式（GNU getopt 风格——原版 Trinity
    /// 主传给 butterfly 的 `--path_reinforcement_distance=25` 即此形式）。
    pub fn value_flag(&mut self, name: &str) -> Result<Option<String>, CliError> {
        let eq_prefix = format!("{name}=");
        let mut value = None;
        let mut kept = Vec::with_capacity(self.args.len());
        let mut i = 0;
        while i < self.args.len() {
            if self.args[i] == name {
                let v = self
                    .args
                    .get(i + 1)
                    .ok_or_else(|| CliError::usage(format!("{name} expects a value")))?;
                value = Some(v.clone());
                i += 2;
            } else if name.starts_with("--") && self.args[i].starts_with(&eq_prefix) {
                value = Some(self.args[i][eq_prefix.len()..].to_string());
                i += 1;
            } else {
                kept.push(self.args[i].clone());
                i += 1;
            }
        }
        if value.is_some() {
            self.given.insert(name.to_string());
        }
        self.args = kept;
        Ok(value)
    }

    /// 必填值 flag: 缺 → `{name} required`。
    pub fn req_flag(&mut self, name: &str) -> Result<String, CliError> {
        self.value_flag(name)?
            .ok_or_else(|| CliError::usage(format!("{name} required")))
    }

    /// 类型化数值 flag: 缺省取 default;解析失败（含 T 的越界/格式错）报
    /// `invalid {name} value: {v}`。T 按调用点推断（u32/usize/f32/…）——
    /// 原 parse_uint/parse_float 的边界语义由所选 T 完整保留。
    pub fn typed_flag<T: FromStr>(&mut self, name: &str, default: T) -> Result<T, CliError> {
        match self.value_flag(name)? {
            None => Ok(default),
            Some(v) => v
                .parse::<T>()
                .map_err(|_| CliError::usage(format!("invalid {name} value: {v}"))),
        }
    }

    /// 整数 flag（`typed_flag` 的具名封装;u32/u64/usize…按调用点推断）。
    pub fn uint_flag<T: FromStr>(&mut self, name: &str, default: T) -> Result<T, CliError> {
        self.typed_flag(name, default)
    }

    /// 整数 flag（i32/i64…;同 [`Cli::uint_flag`]）。
    pub fn int_flag<T: FromStr>(&mut self, name: &str, default: T) -> Result<T, CliError> {
        self.typed_flag(name, default)
    }

    /// 浮点 flag（f32/f64 按调用点推断;有限性检查由需要的一侧自行追加）。
    pub fn float_flag<T: FromStr>(&mut self, name: &str, default: T) -> Result<T, CliError> {
        self.typed_flag(name, default)
    }

    /// 字符串 flag（缺省取 default）。
    pub fn str_flag(&mut self, name: &str, default: &str) -> Result<String, CliError> {
        Ok(self
            .value_flag(name)?
            .unwrap_or_else(|| default.to_string()))
    }

    /// 消费布尔 flag 的**全部**出现。
    pub fn bool_flag(&mut self, name: &str) -> bool {
        let before = self.args.len();
        self.args.retain(|t| t != name);
        before != self.args.len()
    }

    /// 该值 flag 是否被显式给出（须在本名消费之后查询）。
    pub fn was_given(&self, name: &str) -> bool {
        self.given.contains(name)
    }

    /// 有未消费参数 → Err（列出首个;原版 ArgProcessor "do not understand option"）。
    pub fn finish(&self) -> Result<(), CliError> {
        match self.args.first() {
            None => Ok(()),
            Some(tok) => Err(CliError::usage(format!("do not understand option: {tok}"))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(items: &[&str]) -> Vec<String> {
        items.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn value_flag_takes_pair_and_frees_finish() {
        let a = argv(&["--reads", "in.fa", "-o", "out.fa"]);
        let mut cli = Cli::new(&a);
        assert_eq!(cli.value_flag("--reads").unwrap(), Some("in.fa".into()));
        assert_eq!(cli.value_flag("-o").unwrap(), Some("out.fa".into()));
        assert_eq!(cli.value_flag("--absent").unwrap(), None);
        assert!(cli.finish().is_ok());
    }

    #[test]
    fn value_flag_missing_value_errors() {
        let a = argv(&["--reads"]);
        let mut cli = Cli::new(&a);
        let err = cli.value_flag("--reads").unwrap_err();
        assert!(err.usage);
        assert_eq!(err.msg, "--reads expects a value");
    }

    #[test]
    fn repeated_value_flag_last_wins_and_all_consumed() {
        let a = argv(&["--max_cov", "5", "--max_cov", "10"]);
        let mut cli = Cli::new(&a);
        assert_eq!(cli.value_flag("--max_cov").unwrap(), Some("10".into()));
        assert!(cli.was_given("--max_cov"));
        assert!(cli.finish().is_ok(), "重复出现须全部消费");
    }

    #[test]
    fn flag_value_may_look_like_another_flag() {
        // 原版语义: 值 flag 无条件吃下一个 token（--SS_lib_type --DS 的值即 "--DS"）
        let a = argv(&["--SS_lib_type", "--DS"]);
        let mut cli = Cli::new(&a);
        assert_eq!(
            cli.value_flag("--SS_lib_type").unwrap(),
            Some("--DS".into())
        );
        assert!(cli.finish().is_ok());
    }

    #[test]
    fn bool_flag_consumes_all_occurrences() {
        let a = argv(&["--DS", "--run_inchworm", "--DS"]);
        let mut cli = Cli::new(&a);
        assert!(cli.bool_flag("--DS"));
        assert!(cli.bool_flag("--run_inchworm"));
        assert!(!cli.bool_flag("--PARALLEL_IWORM"));
        assert!(cli.finish().is_ok());
    }

    #[test]
    fn unknown_option_reported_by_finish() {
        let a = argv(&["--reads", "in.fa", "--bogus"]);
        let mut cli = Cli::new(&a);
        cli.value_flag("--reads").unwrap();
        let err = cli.finish().unwrap_err();
        assert!(err.usage);
        assert_eq!(err.msg, "do not understand option: --bogus");
    }

    #[test]
    fn required_flag_missing() {
        let mut cli = Cli::new(&argv(&[]));
        let err = cli.req_flag("--kmers").unwrap_err();
        assert!(err.usage);
        assert_eq!(err.msg, "--kmers required");
    }

    #[test]
    fn uint_flag_default_and_invalid() {
        let mut cli = Cli::new(&argv(&["--monitor", "3"]));
        assert_eq!(cli.uint_flag::<u32>("--monitor", 0).unwrap(), 3);
        assert!(cli.was_given("--monitor"));
        let mut cli = Cli::new(&argv(&[]));
        assert_eq!(cli.uint_flag::<u32>("--monitor", 0).unwrap(), 0);
        assert!(!cli.was_given("--monitor"));
        let mut cli = Cli::new(&argv(&["--num_threads", "x"]));
        let err = cli.uint_flag::<u32>("--num_threads", 1).unwrap_err();
        assert!(err.usage);
        assert_eq!(err.msg, "invalid --num_threads value: x");
    }

    #[test]
    fn uint_flag_respects_type_bounds() {
        // u32 边界: 4294967295 可解析、4294967296 越界报错（原 parse_uint 语义）
        let mut cli = Cli::new(&argv(&["--max", "4294967295"]));
        assert_eq!(cli.uint_flag::<u32>("--max", 0).unwrap(), u32::MAX);
        let mut cli = Cli::new(&argv(&["--max", "4294967296"]));
        let err = cli.uint_flag::<u32>("--max", 0).unwrap_err();
        assert_eq!(err.msg, "invalid --max value: 4294967296");
        // 负数对无符号类型同样是解析失败
        let mut cli = Cli::new(&argv(&["--max", "-1"]));
        assert!(cli.uint_flag::<u32>("--max", 0).is_err());
    }

    #[test]
    fn float_flag_parse_and_message() {
        let mut cli = Cli::new(&argv(&["--min_cov", "0.5", "--max_CV", "1e4"]));
        assert_eq!(cli.float_flag::<f64>("--min_cov", 1.0).unwrap(), 0.5);
        assert_eq!(cli.float_flag::<f64>("--max_CV", 0.0).unwrap(), 10000.0);
        let mut cli = Cli::new(&argv(&["--min_cov", "abc"]));
        let err = cli.float_flag::<f64>("--min_cov", 1.0).unwrap_err();
        assert_eq!(err.msg, "invalid --min_cov value: abc");
        // f32/f64 均可（inchworm 侧 f32、kmer 侧 f64——由 T 推断各自保留精度语义）
        let mut cli = Cli::new(&argv(&["--entropy", "1.5"]));
        assert_eq!(cli.float_flag::<f32>("--entropy", 0.0).unwrap(), 1.5f32);
        assert!(cli.finish().is_ok());
    }

    #[test]
    fn int_and_str_flags() {
        let mut cli = Cli::new(&argv(&["--depth", "-3"]));
        assert_eq!(cli.int_flag::<i32>("--depth", 0).unwrap(), -3);
        let mut cli = Cli::new(&argv(&[]));
        assert_eq!(cli.str_flag("--SS_lib_type", "").unwrap(), "");
        let mut cli = Cli::new(&argv(&["--SS_lib_type", "RF"]));
        assert_eq!(cli.str_flag("--SS_lib_type", "").unwrap(), "RF");
        assert!(cli.finish().is_ok());
    }

    #[test]
    fn common_error_maps_to_run_class() {
        let e: CliError = CommonError::KmerTooLong { len: 40 }.into();
        assert!(!e.usage);
        assert_eq!(e.msg, "error, kmer length exceeds 32: 40");
    }

    #[test]
    fn was_given_distinguishes_explicit_default_from_absent() {
        let mut cli = Cli::new(&argv(&["--monitor", "0"]));
        assert_eq!(cli.uint_flag::<u32>("--monitor", 0).unwrap(), 0);
        assert!(cli.was_given("--monitor"));
        assert!(!cli.was_given("--num_threads"));
    }
}

#[cfg(test)]
mod eq_form_tests {
    use super::*;

    fn argv(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn long_flag_equals_form() {
        let a = argv(&["--path_reinforcement_distance=25", "-C", "c0.graph"]);
        let mut cli = Cli::new(&a);
        assert_eq!(
            cli.value_flag("--path_reinforcement_distance").unwrap(),
            Some("25".into())
        );
        assert_eq!(cli.req_flag("-C").unwrap(), "c0.graph");
        assert!(cli.finish().is_ok());
    }

    #[test]
    fn short_flag_no_equals_match() {
        // 短旗标（-K）不启用等号形式——防止 "-K=25" 被误认
        let a = argv(&["-K=25"]);
        let mut cli = Cli::new(&a);
        assert_eq!(cli.value_flag("-K").unwrap(), None);
        assert!(cli.finish().is_err());
    }

    #[test]
    fn mixed_forms_last_wins() {
        let a = argv(&["--prd", "10", "--prd=20"]);
        let mut cli = Cli::new(&a);
        assert_eq!(cli.value_flag("--prd").unwrap(), Some("20".into()));
        assert!(cli.finish().is_ok());
    }
}
