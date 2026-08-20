//! `partition_chrysalis_graphs_n_reads.pl`（util/support_scripts，249 行）移植
//! ——P3-T6 的 Chrysalis 收尾分区器。
//!
//! 镜像语义：
//! - **图分区**：按 `Component` 行切块（Component_reader::next_component）；
//!   `num_kmers = 块体行数`；`num_kmers + 24 < min_contig_length` 跳过
//!   （假定 25-mer）；通过者从 1 计数，目录 `Cbin<通过数/1000>`，每组件写
//!   `<Cbin>/c<id>.graph.tmp`（完整块文本含头行），登记 `comp_id → base`；
//!   **块体为空的组件终止整个读取循环**（perl `if (@lines) ... else undef`）。
//! - **reads 分区**：输入行 `comp \t acc \t pct \t read`（假定已按 comp 排序，
//!   未排序时同 comp 第二段会重开文件截断——照抄）；comp 变化时开新文件
//!   （登记表有的才写，没有的静默丢）；每条 `>acc pct\nread\n`（acc 含 '>'）。
//! - **component_base_listing.txt**：comp id 数值升序，graph.tmp 与 reads.tmp
//!   **都存在且非空**才输出 `id \t base`，同时作为本函数返回值。

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use trinity_common::error::CommonError;

/// `-N`（每分区图数，默认 1000）与 `-L`（最短 contig 长度，默认 200；
/// 判据 `num_kmers + 24 < L`，25-mer 假定）。
#[derive(Debug, Clone)]
pub struct PartParams {
    pub graphs_per_partition: usize,
    pub min_contig_length: usize,
}

impl Default for PartParams {
    fn default() -> Self {
        PartParams {
            graphs_per_partition: 1000,
            min_contig_length: 200,
        }
    }
}

/// 分区主流程。`debruijn_text` = bundled_iworm_contigs.fasta.deBruijn 全文，
/// `sorted_reads_to_components` = readsToComponents.out.sort 全文（已按 comp
/// 排序），输出目录 `out_dir`（对应原版 `<db 目录>/Component_bins`）。
/// 返回 listing 内容（comp id 升序的 (id, base)，双非空过滤后）。
pub fn partition(
    debruijn_text: &str,
    sorted_reads_to_components: &str,
    out_dir: &Path,
    p: &PartParams,
) -> Result<Vec<(u64, PathBuf)>, CommonError> {
    std::fs::create_dir_all(out_dir)?;

    let mut registry: HashMap<u64, PathBuf> = HashMap::new();
    let mut passed: usize = 0;

    for comp in components(debruijn_text) {
        if comp.body_lines + 24 < p.min_contig_length {
            continue;
        }
        passed += 1;
        let cbin = out_dir.join(format!("Cbin{}", passed / p.graphs_per_partition));
        std::fs::create_dir_all(&cbin)?;
        let graph_tmp = cbin.join(format!("c{}.graph.tmp", comp.id));
        std::fs::write(&graph_tmp, comp.text)?;
        registry.insert(comp.id, cbin.join(format!("c{}", comp.id)));
    }

    // reads 分区：comp 变化时（数值比较）开新文件；未登记的 comp 静默丢弃
    let mut prev_comp: i64 = -1;
    let mut cur: Option<std::fs::File> = None;
    use std::io::Write;
    for line in sorted_reads_to_components.lines() {
        let mut f = line.split('\t');
        let comp = perl_num(f.next().unwrap_or(""));
        if comp != prev_comp {
            cur = None;
            if let Some(base) = registry.get(&(comp as u64)) {
                cur = Some(std::fs::File::create(base.with_extension("reads.tmp"))?);
            }
        }
        if let Some(fh) = cur.as_mut() {
            let acc = f.next().unwrap_or("");
            let pct = f.next().unwrap_or("");
            let read = f.next().unwrap_or("");
            write!(fh, "{acc} {pct}\n{read}\n")?;
        }
        prev_comp = comp;
    }

    // listing：数值升序 + 双非空
    let mut ids: Vec<u64> = registry.keys().copied().collect();
    ids.sort_unstable();
    let mut listing = Vec::new();
    let mut text = String::new();
    for id in ids {
        let base = &registry[&id];
        let g = base.with_extension("graph.tmp");
        let r = base.with_extension("reads.tmp");
        let ok = files_nonempty(&g) && files_nonempty(&r);
        if ok {
            text.push_str(&format!("{id}\t{}\n", base.display()));
            listing.push((id, base.clone()));
        }
    }
    std::fs::write(out_dir.join("component_base_listing.txt"), text)?;
    Ok(listing)
}

fn files_nonempty(p: &Path) -> bool {
    matches!(std::fs::metadata(p), Ok(m) if m.len() > 0)
}

/// perl 数值上下文（`$comp != $prev_comp`）：最长数字前缀，无数字 → 0。
fn perl_num(s: &str) -> i64 {
    let t = s.trim_start();
    let mut chars = t.chars().peekable();
    let neg = matches!(chars.peek(), Some('-'));
    if matches!(chars.peek(), Some('+') | Some('-')) {
        chars.next();
    }
    let mut v: i64 = 0;
    for c in chars {
        if let Some(d) = c.to_digit(10) {
            v = v.saturating_mul(10).saturating_add(d as i64);
        } else {
            break;
        }
    }
    if neg {
        -v
    } else {
        v
    }
}

struct Component {
    id: u64,
    body_lines: usize,
    text: String,
}

/// Component_reader：首行必须是 `Component`；块体 = 到下一 `Component` 行
/// 之前的所有行（含空行）；**空块体终止读取**（perl `if (@lines) ... else
/// return undef` 怪癖——其后组件全部不读）。id 取 `Component (\d+)`。
fn components(text: &str) -> Vec<Component> {
    // perl 行迭代语义：末尾单个换行不产生额外空行，双换行保留一个空行
    let mut lines: Vec<&str> = text.split('\n').collect();
    if lines.last() == Some(&"") {
        lines.pop();
    }
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < lines.len() {
        let header = lines[i];
        if !header.starts_with("Component") {
            break; // perl confess；上游保证首行即 Component，防御性终止
        }
        let id = regex_component_id(header);
        let mut j = i + 1;
        while j < lines.len() && !lines[j].starts_with("Component") {
            j += 1;
        }
        let body = &lines[i + 1..j];
        if body.is_empty() {
            break; // 空块体 → undef → 整个读取终止
        }
        let mut buf = String::with_capacity(header.len() + 1);
        buf.push_str(header);
        buf.push('\n');
        for b in body {
            buf.push_str(b);
            buf.push('\n');
        }
        out.push(Component {
            id,
            body_lines: body.len(),
            text: buf,
        });
        i = j;
    }
    out
}

/// `Component (\d+)`：`Component` 后首个数字段。
fn regex_component_id(line: &str) -> u64 {
    let rest = line.split_once("Component").map(|(_, r)| r).unwrap_or(line);
    let digits: String = rest
        .chars()
        .skip_while(|c| !c.is_ascii_digit())
        .take_while(|c| c.is_ascii_digit())
        .collect();
    digits.parse().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmpdir(name: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(name);
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn comp(id: u64, n_lines: usize) -> String {
        let mut s = format!("Component {id}\n");
        for i in 0..n_lines {
            s.push_str(&format!("{i}\t-1\t1\tACGTACGTACGTACGTACGTACGT\t1\n"));
        }
        s
    }

    #[test]
    fn length_filter_and_cbin_split() {
        // min_contig_length=200：num_kmers+24 < 200 跳过 → 176 行体通过
        let db = comp(1, 10) + &comp(2, 176) + &comp(3, 175);
        let d = tmpdir("p3_t6_part_len");
        let l = partition(
            &db,
            "",
            &d,
            &PartParams {
                graphs_per_partition: 2,
                min_contig_length: 200,
            },
        )
        .unwrap();
        // 只有 comp 2 通过；reads 全空 → listing 空（reads.tmp 不存在）
        assert!(l.is_empty());
        assert!(d.join("Cbin0/c2.graph.tmp").metadata().unwrap().len() > 0);
        assert!(!d.join("Cbin0/c1.graph.tmp").exists());
        // listing 文件本身存在
        assert!(d.join("component_base_listing.txt").exists());
    }

    #[test]
    fn cbin_indexing_by_passed_count() {
        // 5 个 200 行组件、N=2 → Cbin0×2、Cbin1×2、Cbin2×1
        let db = comp(1, 200) + &comp(2, 200) + &comp(3, 200) + &comp(4, 200) + &comp(5, 200);
        let d = tmpdir("p3_t6_part_cbin");
        let reads = "5\t>r9\t50%\tACGT\n5\t>r8\t50%\tACGT\n1\t>r1\t50%\tTTTT\n";
        let l = partition(
            &db,
            reads,
            &d,
            &PartParams {
                graphs_per_partition: 2,
                min_contig_length: 200,
            },
        )
        .unwrap();
        assert_eq!(l.len(), 2);
        assert_eq!(l[0], (1, d.join("Cbin0/c1")));
        assert_eq!(l[1], (5, d.join("Cbin2/c5")));
        // graph.tmp 含头行的完整块
        let g = std::fs::read_to_string(d.join("Cbin0/c1.graph.tmp")).unwrap();
        assert!(g.starts_with("Component 1\n"));
        assert_eq!(g.lines().count(), 201);
        // reads.tmp：acc 含 '>'，pct 分隔；未登记 comp（此处 comp 5 的 reads
        // 有写入）… comp 1 有 1 条
        let r1 = std::fs::read_to_string(d.join("Cbin0/c1.reads.tmp")).unwrap();
        assert_eq!(r1, ">r1 50%\nTTTT\n");
        let r5 = std::fs::read_to_string(d.join("Cbin2/c5.reads.tmp")).unwrap();
        assert_eq!(r5, ">r9 50%\nACGT\n>r8 50%\nACGT\n");
        // listing 内容
        let lst = std::fs::read_to_string(d.join("component_base_listing.txt")).unwrap();
        assert_eq!(
            lst,
            format!(
                "1\t{}\n5\t{}\n",
                d.join("Cbin0/c1").display(),
                d.join("Cbin2/c5").display()
            )
        );
    }

    #[test]
    fn unregistered_comp_reads_silently_dropped() {
        let db = comp(7, 200);
        let reads = "99\t>rA\t50%\tACGT\n7\t>rB\t50%\tACGT\n";
        let d = tmpdir("p3_t6_part_drop");
        let l = partition(&db, reads, &d, &PartParams::default()).unwrap();
        assert_eq!(l, vec![(7, d.join("Cbin0/c7"))]);
        assert!(!d.join("Cbin0/c99.reads.tmp").exists());
        assert_eq!(
            std::fs::read_to_string(d.join("Cbin0/c7.reads.tmp")).unwrap(),
            ">rB 50%\nACGT\n"
        );
    }

    #[test]
    fn listing_requires_nonempty_reads() {
        let db = comp(7, 200) + &comp(8, 200);
        // comp 8 无 reads → 不进 listing；comp 7 有
        let reads = "7\t>rB\t50%\tACGT\n";
        let d = tmpdir("p3_t6_part_listing");
        let l = partition(&db, reads, &d, &PartParams::default()).unwrap();
        assert_eq!(l, vec![(7, d.join("Cbin0/c7"))]);
        // c8.graph.tmp 已写出（非空）但 reads.tmp 缺 → 被过滤
        assert!(d.join("Cbin0/c8.graph.tmp").exists());
        assert!(!d.join("Cbin0/c8.reads.tmp").exists());
    }

    #[test]
    fn empty_component_terminates_reading() {
        // 块体为空的组件之后的所有组件不再读（perl next_component 怪癖）
        let db = comp(1, 200) + "Component 2\n" + &comp(3, 200);
        let d = tmpdir("p3_t6_part_empty");
        let l = partition(&db, "", &d, &PartParams::default()).unwrap();
        assert!(l.is_empty());
        assert!(d.join("Cbin0/c1.graph.tmp").exists());
        assert!(!d.join("Cbin0/c3.graph.tmp").exists());
    }

    #[test]
    fn perl_num_semantics() {
        assert_eq!(perl_num("23"), 23);
        assert_eq!(perl_num(" 7"), 7);
        assert_eq!(perl_num("abc"), 0);
        assert_eq!(perl_num("-1"), -1);
        assert_eq!(perl_num("12abc"), 12);
    }
}
