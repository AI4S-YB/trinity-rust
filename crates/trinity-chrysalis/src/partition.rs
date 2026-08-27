//! `partition_chrysalis_graphs_n_reads.pl`（util/support_scripts，249 行）移植
//! ——P3-T6 的 Chrysalis 收尾分区器。
//!
//! 镜像语义：
//! - **图分区**：按 `Component` 行切块（Component_reader::next_component）；
//!   `num_kmers = 块体行数`；`num_kmers + 24 < min_contig_length` 跳过
//!   （假定 25-mer）；通过者从 1 计数，目录 `Cbin<通过数/1000>`，每组件一条
//!   graph.tmp 记录（完整块文本含头行），登记 `comp_id → base`；
//!   **块体为空的组件终止整个读取循环**（perl `if (@lines) ... else undef`）。
//! - **reads 分区**：输入行 `comp \t acc \t pct \t read`（假定已按 comp 排序，
//!   未排序时同 comp 第二段**截断重写**——由 pack 索引"同 key 后者覆盖"镜像）；
//!   comp 变化时换组件（登记表有的才写，没有的静默丢）；每条 `>acc pct\nread\n`
//!   （acc 含 '>'）。
//! - **component_base_listing.txt**：comp id 数值升序，graph.tmp 与 reads.tmp
//!   记录**都存在且非空**才输出 `id \t base`，同时作为本函数返回值。
//!
//! 输出布局（IO 优化, M8）：per-component 的 cN.graph.tmp / cN.reads.tmp 文件
//! 改为每 Cbin 一个 `tmp.pack`+`tmp.idx` 打包（见 pack 模块）——28444 组件
//! 场景再省 ~5.7 万次 create/unlink。逻辑布局（id → 记录字节）与旧文件布局
//! 逐字节等价（partition 测试的 oracle 对拍锁定）。

use std::borrow::Cow;
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use trinity_common::error::CommonError;

use crate::pack::{PackWriter, KIND_GRAPH_TMP, KIND_READS_TMP, TMP_PACK};

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

/// 分区主流程。`debruijn_text` = bundled_iworm_contigs.fasta.deBruijn 全文,
/// `sorted_reads_to_components` = readsToComponents.out.sort 全文（已按 comp
/// 排序）, 输出目录 `out_dir`（对应原版 `<db 目录>/Component_bins`）——
/// 产物为每 Cbin 一个 `tmp.pack`+`tmp.idx`（打包布局, 见 pack 模块）。
/// 返回 listing 内容（comp id 升序的 (id, base)，双非空过滤后）。
pub fn partition(
    debruijn_text: &str,
    sorted_reads_to_components: &str,
    out_dir: &Path,
    p: &PartParams,
) -> Result<Vec<(u64, PathBuf)>, CommonError> {
    std::fs::create_dir_all(out_dir)?;

    let mut registry: HashMap<u64, PathBuf> = HashMap::new();
    // 每 Cbin 一个 tmp 包写句柄（延迟创建; 数量 = 组件数/N, 很小）。
    let mut writers: HashMap<PathBuf, PackWriter> = HashMap::new();
    let mut passed: usize = 0;

    for comp in components(debruijn_text) {
        if comp.body_lines + 24 < p.min_contig_length {
            continue;
        }
        passed += 1;
        let cbin = out_dir.join(format!("Cbin{}", passed / p.graphs_per_partition));
        let w = match writers.get_mut(&cbin) {
            Some(w) => w,
            None => {
                std::fs::create_dir_all(&cbin)?;
                writers.insert(cbin.clone(), PackWriter::create(&cbin, TMP_PACK)?);
                writers.get_mut(&cbin).unwrap()
            }
        };
        w.append(comp.id, KIND_GRAPH_TMP, comp.text.as_bytes())?;
        registry.insert(comp.id, cbin.join(format!("c{}", comp.id)));
    }

    // reads 分区: comp 变化时（数值比较）flush 上一组件的缓冲为一条记录;
    // 未登记的 comp 静默丢弃。同 comp 再现（未排序输入）= 截断重写——
    // 追加新记录 + reads_len 覆盖, 由 PackIndex "同 key 后者覆盖" 镜像。
    // （原版此时对已截断文件持续写入; 末次内容 = 有效记录。）
    let mut prev_comp: i64 = -1;
    let mut prev_registered = false;
    let mut buf: Vec<u8> = Vec::new();
    let mut reads_len: HashMap<u64, u64> = HashMap::new(); // comp → 末次记录长度
    for line in sorted_reads_to_components.lines() {
        let mut f = line.split('\t');
        let comp = perl_num(f.next().unwrap_or(""));
        if comp != prev_comp {
            flush_reads(
                &mut writers,
                &registry,
                &mut reads_len,
                prev_comp,
                prev_registered,
                &buf,
            )?;
            prev_registered = registry.contains_key(&(comp as u64));
            buf.clear();
        }
        if prev_registered {
            let acc = f.next().unwrap_or("");
            let pct = f.next().unwrap_or("");
            let read = f.next().unwrap_or("");
            write!(buf, "{acc} {pct}\n{read}\n")?;
        }
        prev_comp = comp;
    }
    flush_reads(
        &mut writers,
        &registry,
        &mut reads_len,
        prev_comp,
        prev_registered,
        &buf,
    )?;
    for w in writers.values_mut() {
        w.finish()?;
    }

    // listing：数值升序 + 双非空（graph 记录恒非空; reads 看末次记录长度）
    let mut ids: Vec<u64> = registry.keys().copied().collect();
    ids.sort_unstable();
    let mut listing = Vec::new();
    let mut text = String::new();
    for id in ids {
        let base = &registry[&id];
        let ok = reads_len.get(&id).copied().unwrap_or(0) > 0;
        if ok {
            text.push_str(&format!("{id}\t{}\n", base.display()));
            listing.push((id, base.clone()));
        }
    }
    std::fs::write(out_dir.join("component_base_listing.txt"), text)?;
    Ok(listing)
}

/// 把上一组件的 reads 缓冲 flush 为一条 tmp 包记录（已登记组件才写;
/// 空缓冲也写——镜像原版"comp 变化即创建（截断）文件"的空文件语义）。
fn flush_reads(
    writers: &mut HashMap<PathBuf, PackWriter>,
    registry: &HashMap<u64, PathBuf>,
    reads_len: &mut HashMap<u64, u64>,
    prev_comp: i64,
    prev_registered: bool,
    buf: &[u8],
) -> Result<(), CommonError> {
    if !prev_registered {
        return Ok(());
    }
    let id = prev_comp as u64;
    let cbin = registry[&id]
        .parent()
        .map(|d| d.to_path_buf())
        .unwrap_or_default();
    let w = writers
        .get_mut(&cbin)
        .expect("已登记组件的 Cbin 写句柄必存在");
    w.append(id, KIND_READS_TMP, buf)?;
    reads_len.insert(id, buf.len() as u64);
    Ok(())
}

/// perl 数值上下文（`$comp != $prev_comp`）：最长数字前缀，无数字 → 0。
/// 快路径：首字节即 ASCII 数字（真实数据全部如此），零开销字节扫描；
/// 其余（空白/正负号/前缀垃圾）走完整 char 语义慢路径。
fn perl_num(s: &str) -> i64 {
    let b = s.as_bytes();
    if let Some(&c0) = b.first() {
        if c0.is_ascii_digit() {
            let mut v: i64 = 0;
            let mut i = 0;
            while i < b.len() {
                let d = b[i];
                if !d.is_ascii_digit() {
                    break;
                }
                v = v.saturating_mul(10).saturating_add((d - b'0') as i64);
                i += 1;
            }
            return v;
        }
    }
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

struct Component<'a> {
    id: u64,
    body_lines: usize,
    text: Cow<'a, str>,
}

/// Component_reader：首行必须是 `Component`；块体 = 到下一 `Component` 行
/// 之前的所有行（含空行）；**空块体终止读取**（perl `if (@lines) ... else
/// return undef` 怪癖——其后组件全部不读）。id 取 `Component (\d+)`。
///
/// 实现为原文借用切片（perl 行语义：末尾单个换行不产生额外空行，双换行
/// 保留一个空行——`"\nComponent"` 搜索天然满足）。仅 EOF 无换行的末组件
/// 需要 Owned 补 '\n'。
fn components(text: &str) -> Vec<Component<'_>> {
    let mut out = Vec::new();
    if !text.starts_with("Component") {
        return out; // perl confess；上游保证首行即 Component，防御性终止
    }
    let bytes = text.as_bytes();
    let mut start = 0usize; // 当前组件头行起始
    loop {
        // 头行 = [start, nl)；nl = 行尾换行或文本末尾
        let nl = match bytes[start..].iter().position(|&b| b == b'\n') {
            Some(off) => start + off,
            None => text.len(),
        };
        let header = &text[start..nl];
        // 下一 `Component` 行：搜 "\nComponent"，命中处 +1 即下一头行起始
        let next = text[nl..].find("\nComponent").map(|off| nl + off + 1);
        let end = next.unwrap_or(text.len());
        let body_start = (nl + 1).min(end);
        // 体行数 = 体区域换行数（EOF 无换行的末行 +1）
        let mut body_lines = text[body_start..end].matches('\n').count();
        if next.is_none() && body_start < end && !text.ends_with('\n') {
            body_lines += 1;
        }
        if body_lines == 0 {
            break; // 空块体 → undef → 整个读取终止
        }
        let comp_text = if next.is_some() || text.ends_with('\n') {
            Cow::Borrowed(&text[start..end])
        } else {
            // EOF 无换行：原语义 join+"\n" 需补行尾
            Cow::Owned(format!("{}\n", &text[start..end]))
        };
        out.push(Component {
            id: regex_component_id(header),
            body_lines,
            text: comp_text,
        });
        start = match next {
            Some(nh) => nh,
            None => break, // 已到文本末尾
        };
    }
    out
}

/// `Component (\d+)`：`Component` 后首个数字段（零分配）。
fn regex_component_id(line: &str) -> u64 {
    let rest = line.split_once("Component").map(|(_, r)| r).unwrap_or(line);
    let b = rest.as_bytes();
    let s = b.iter().position(|c| c.is_ascii_digit()).unwrap_or(b.len());
    let e = b[s..]
        .iter()
        .position(|c| !c.is_ascii_digit())
        .map_or(b.len(), |off| s + off);
    rest[s..e].parse().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pack::{self, PackIndex};

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

    /// 打包布局的"逐组件可定位"读接口（等价旧布局的读文件断言）:
    /// 读 `dir/Cbin*? /tmp.pack` 中 (id, kind) 的记录。
    fn tmp_record(cbin: &Path, id: u64, kind: u64) -> Option<Vec<u8>> {
        let idx = PackIndex::load(&cbin.join(format!("{}.idx", pack::TMP_PACK))).unwrap();
        let (off, len) = idx.find(id, kind)?;
        let f = std::fs::File::open(cbin.join(format!("{}.pack", pack::TMP_PACK))).unwrap();
        Some(pack::read_at(&f, off, len).unwrap())
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
        // 只有 comp 2 通过；reads 全空 → listing 空（reads 记录不存在）
        assert!(l.is_empty());
        let cbin0 = d.join("Cbin0");
        assert!(!tmp_record(&cbin0, 2, pack::KIND_GRAPH_TMP)
            .unwrap()
            .is_empty());
        assert!(tmp_record(&cbin0, 1, pack::KIND_GRAPH_TMP).is_none());
        // listing 文件本身存在; per-component 小文件不再产生
        assert!(d.join("component_base_listing.txt").exists());
        assert!(!d.join("Cbin0/c2.graph.tmp").exists());
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
        let g = String::from_utf8(tmp_record(&d.join("Cbin0"), 1, pack::KIND_GRAPH_TMP).unwrap())
            .unwrap();
        assert!(g.starts_with("Component 1\n"));
        assert_eq!(g.lines().count(), 201);
        // reads.tmp：acc 含 '>'，pct 分隔；comp 1 有 1 条
        let r1 = tmp_record(&d.join("Cbin0"), 1, pack::KIND_READS_TMP).unwrap();
        assert_eq!(r1, b">r1 50%\nTTTT\n");
        let r5 = tmp_record(&d.join("Cbin2"), 5, pack::KIND_READS_TMP).unwrap();
        assert_eq!(r5, b">r9 50%\nACGT\n>r8 50%\nACGT\n");
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
        assert!(tmp_record(&d.join("Cbin0"), 99, pack::KIND_READS_TMP).is_none());
        assert_eq!(
            tmp_record(&d.join("Cbin0"), 7, pack::KIND_READS_TMP).unwrap(),
            b">rB 50%\nACGT\n"
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
        // c8 的 graph 记录已写入（非空）但 reads 记录缺 → 被过滤
        assert!(!tmp_record(&d.join("Cbin0"), 8, pack::KIND_GRAPH_TMP)
            .unwrap()
            .is_empty());
        assert!(tmp_record(&d.join("Cbin0"), 8, pack::KIND_READS_TMP).is_none());
    }

    #[test]
    fn empty_component_terminates_reading() {
        // 块体为空的组件之后的所有组件不再读（perl next_component 怪癖）
        let db = comp(1, 200) + "Component 2\n" + &comp(3, 200);
        let d = tmpdir("p3_t6_part_empty");
        let l = partition(&db, "", &d, &PartParams::default()).unwrap();
        assert!(l.is_empty());
        assert!(tmp_record(&d.join("Cbin0"), 1, pack::KIND_GRAPH_TMP).is_some());
        assert!(tmp_record(&d.join("Cbin0"), 3, pack::KIND_GRAPH_TMP).is_none());
    }

    #[test]
    fn perl_num_semantics() {
        assert_eq!(perl_num("23"), 23);
        assert_eq!(perl_num(" 7"), 7);
        assert_eq!(perl_num("abc"), 0);
        assert_eq!(perl_num("-1"), -1);
        assert_eq!(perl_num("12abc"), 12);
    }

    // =====================================================================
    // 身份 oracle：原实现逐字节快照（勿同步改动）——新实现产物必须与
    // reference 完全一致（全目录树字节 + 返回 listing）。
    // =====================================================================

    struct RefComponent {
        id: u64,
        body_lines: usize,
        text: String,
    }

    fn ref_components(text: &str) -> Vec<RefComponent> {
        let mut lines: Vec<&str> = text.split('\n').collect();
        if lines.last() == Some(&"") {
            lines.pop();
        }
        let mut out = Vec::new();
        let mut i = 0usize;
        while i < lines.len() {
            let header = lines[i];
            if !header.starts_with("Component") {
                break;
            }
            let id = ref_regex_component_id(header);
            let mut j = i + 1;
            while j < lines.len() && !lines[j].starts_with("Component") {
                j += 1;
            }
            let body = &lines[i + 1..j];
            if body.is_empty() {
                break;
            }
            let mut buf = String::with_capacity(header.len() + 1);
            buf.push_str(header);
            buf.push('\n');
            for b in body {
                buf.push_str(b);
                buf.push('\n');
            }
            out.push(RefComponent {
                id,
                body_lines: body.len(),
                text: buf,
            });
            i = j;
        }
        out
    }

    fn ref_regex_component_id(line: &str) -> u64 {
        let rest = line.split_once("Component").map(|(_, r)| r).unwrap_or(line);
        let digits: String = rest
            .chars()
            .skip_while(|c| !c.is_ascii_digit())
            .take_while(|c| c.is_ascii_digit())
            .collect();
        digits.parse().unwrap_or(0)
    }

    fn ref_perl_num(s: &str) -> i64 {
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

    fn ref_partition(
        debruijn_text: &str,
        sorted_reads_to_components: &str,
        out_dir: &Path,
        p: &PartParams,
    ) -> Result<Vec<(u64, PathBuf)>, CommonError> {
        std::fs::create_dir_all(out_dir)?;

        let mut registry: HashMap<u64, PathBuf> = HashMap::new();
        let mut passed: usize = 0;

        for comp in ref_components(debruijn_text) {
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

        let mut prev_comp: i64 = -1;
        let mut cur: Option<std::fs::File> = None;
        use std::io::Write;
        for line in sorted_reads_to_components.lines() {
            let mut f = line.split('\t');
            let comp = ref_perl_num(f.next().unwrap_or(""));
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

        let mut ids: Vec<u64> = registry.keys().copied().collect();
        ids.sort_unstable();
        let mut listing = Vec::new();
        let mut text = String::new();
        for id in ids {
            let base = &registry[&id];
            let g = base.with_extension("graph.tmp");
            let r = base.with_extension("reads.tmp");
            let ok = ref_files_nonempty(&g) && ref_files_nonempty(&r);
            if ok {
                text.push_str(&format!("{id}\t{}\n", base.display()));
                listing.push((id, base.clone()));
            }
        }
        std::fs::write(out_dir.join("component_base_listing.txt"), text)?;
        Ok(listing)
    }

    fn ref_files_nonempty(p: &Path) -> bool {
        matches!(std::fs::metadata(p), Ok(m) if m.len() > 0)
    }

    /// 递归收集 relpath → 内容字节。
    fn tree(root: &Path) -> std::collections::BTreeMap<std::path::PathBuf, Vec<u8>> {
        let mut m = std::collections::BTreeMap::new();
        let mut stack = vec![root.to_path_buf()];
        while let Some(d) = stack.pop() {
            for e in std::fs::read_dir(&d).unwrap() {
                let e = e.unwrap();
                let p = e.path();
                if p.is_dir() {
                    stack.push(p);
                } else {
                    let rel = p.strip_prefix(root).unwrap().to_path_buf();
                    m.insert(rel, std::fs::read(&p).unwrap());
                }
            }
        }
        m
    }

    /// 打包布局目录 → 虚拟"逐组件文件树"（key 与旧文件布局一致:
    /// `CbinN/c<id>.graph.tmp` 等）——oracle 对拍在新旧布局间比较同一逻辑内容,
    /// "逐组件可定位"由 pack 索引上的等价位接口（PackIndex::entries + read_at）
    /// 承担。
    fn tree_packed(root: &Path) -> std::collections::BTreeMap<std::path::PathBuf, Vec<u8>> {
        let mut m = std::collections::BTreeMap::new();
        let lst = root.join("component_base_listing.txt");
        if lst.exists() {
            m.insert(
                std::path::PathBuf::from("component_base_listing.txt"),
                std::fs::read(&lst).unwrap(),
            );
        }
        for e in std::fs::read_dir(root).unwrap() {
            let d = e.unwrap().path();
            if !d.is_dir() {
                continue;
            }
            let cbin = d.file_name().unwrap().to_string_lossy().into_owned();
            let idx = PackIndex::load(&d.join(format!("{}.idx", pack::TMP_PACK))).unwrap();
            let f = std::fs::File::open(d.join(format!("{}.pack", pack::TMP_PACK))).unwrap();
            for (id, kind, off, len) in idx.entries() {
                let ext = match kind {
                    pack::KIND_GRAPH_TMP => "graph.tmp",
                    pack::KIND_READS_TMP => "reads.tmp",
                    other => panic!("unexpected kind {other}"),
                };
                m.insert(
                    std::path::PathBuf::from(format!("{cbin}/c{id}.{ext}")),
                    pack::read_at(&f, off, len).unwrap(),
                );
            }
        }
        m
    }

    /// 全目录树字节对拍（a = reference 文件布局, b = 新打包布局）。
    /// `component_base_listing.txt` 内嵌各自 out_dir 的
    /// 绝对路径——比较前把两个根路径规范化为占位符（其余文件原样字节比）。
    fn assert_trees_identical(a: &Path, b: &Path) {
        let (mut ta, mut tb) = (tree(a), tree_packed(b));
        let keys: Vec<_> = ta.keys().collect();
        assert_eq!(
            keys,
            tb.keys().collect::<Vec<_>>(),
            "文件集合不一致（{a:?} vs {b:?}）"
        );
        for (m, root) in [(&mut ta, a), (&mut tb, b)] {
            if let Some(v) = m.get_mut(std::path::Path::new("component_base_listing.txt")) {
                let s = String::from_utf8_lossy(v).replace(&root.display().to_string(), "<OUT>");
                *v = s.into_bytes();
            }
        }
        for k in ta.keys() {
            assert_eq!(ta[k], tb[k], "内容不一致: {k:?}");
        }
    }

    fn assert_listing_identical(la: &[(u64, PathBuf)], lb: &[(u64, PathBuf)], a: &Path, b: &Path) {
        assert_eq!(la.len(), lb.len());
        for ((ia, ba), (ib, bb)) in la.iter().zip(lb.iter()) {
            assert_eq!(ia, ib);
            assert_eq!(ba.strip_prefix(a).unwrap(), bb.strip_prefix(b).unwrap());
        }
    }

    /// 确定性伪随机（PCG 简版），固定种子 → 测试可复现。
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n.max(1)
        }
    }

    /// 生成边界覆盖的合成输入：变长体（0 行=终止）、体行以 "Component"
    /// 开头（切块边界）、重复 comp id、非单调 id、reads 含未登记/畸形
    /// comp 字段、同 comp 二次出现（截断重开）、CRLF、缺列。
    fn synth_inputs(seed: u64, eof_newline: bool) -> (String, String) {
        let mut rng = Lcg(seed);
        let ids: Vec<u64> = vec![30, 4, 7, 4, 0, 12, 9, 255, 1]; // 含重复 4
        let mut db = String::new();
        for (k, &id) in ids.iter().enumerate() {
            db.push_str(&format!("Component {id}\n"));
            match k {
                3 => {}             // 空块体：终止读取（其后 0/12/9/255/1 全不读）
                4 => db.push('\n'), // 单空行体
                _ => {
                    // 170..229：跨过 176 阈值（min_contig_length=200），种子
                    // 决定通过/被滤的混合 → listing 非空路径也被覆盖
                    let n = 170 + rng.below(60) as usize;
                    for i in 0..n {
                        match rng.below(16) {
                            0 => db.push_str("Componentish\tjunk\n"), // 切块边界
                            1 => db.push('\n'),
                            2 => db.push_str("99999999999999999999\t-1\t1\tAAAA\t1\n"),
                            _ => db.push_str(&format!("{i}\t-1\t1\t{}\t1\n", "ACGT".repeat(8))),
                        }
                    }
                }
            }
        }
        if eof_newline {
            db.push('\n');
        } else {
            // 去掉末尾换行 → 末组件需补 '\n' 路径
            if db.ends_with('\n') {
                db.pop();
            }
        }
        // reads：按 comp 分组（个别重复分组模拟未排序重开）
        let groups: Vec<Vec<&str>> = vec![
            vec!["30", ">a0", "50%", "ACGTACGT"],
            vec!["4", ">a1", "50%", "ACGT"],
            vec!["4x", ">a2", "30%", "TTTT"], // perl_num → 4（同组续写）
            vec!["+7", ">a3", "50%", "GGGG"], // → 7
            vec!["abc", ">a4", "50%", "CCCC"], // → 0：已登记 comp 0
            vec!["999", ">a5", "50%", "AAAA"], // 未登记 → 丢弃
            vec!["-1", ">a6", "50%", "TTTT"], // → -1 未登记
            vec!["30", ">a7", "50%", "重开截断"], // 30 二次分组：重开截断
            vec![" 12", ">a8", "50%", "AAAA"], // → 12 未登记（空块体后）
            vec!["12", ">a9", "", ""],        // 空字段
        ];
        let mut reads = String::new();
        for g in &groups {
            let reps = 1 + rng.below(3);
            for r in 0..reps {
                reads.push_str(&format!("{}\t{}\t{}\t{}\n", g[0], g[1], g[2], g[3]));
                if r == 0 && rng.below(2) == 0 {
                    reads.push_str(&format!("{}\t{}\t{}\t{}\r\n", g[0], g[1], g[2], g[3]));
                }
            }
        }
        if !eof_newline && reads.ends_with('\n') {
            reads.pop(); // reads 侧也走 EOF 无换行路径
        }
        (db, reads)
    }

    #[test]
    fn identity_vs_reference_snapshot() {
        for seed in [1u64, 42, 20260826] {
            for &eof_nl in &[true, false] {
                let (db, reads) = synth_inputs(seed, eof_nl);
                for &n in &[1usize, 2, 1000] {
                    let p = PartParams {
                        graphs_per_partition: n,
                        min_contig_length: 200,
                    };
                    let da = tmpdir(&format!("p3_t6_ident_a_{seed}_{eof_nl}_{n}"));
                    let db_ = tmpdir(&format!("p3_t6_ident_b_{seed}_{eof_nl}_{n}"));
                    let la = ref_partition(&db, &reads, &da, &p).unwrap();
                    let lb = partition(&db, &reads, &db_, &p).unwrap();
                    assert_trees_identical(&da, &db_);
                    assert_listing_identical(&la, &lb, &da, &db_);
                    let _ = std::fs::remove_dir_all(&da);
                    let _ = std::fs::remove_dir_all(&db_);
                }
            }
        }
    }

    /// 退化输入的身份对拍：空文本、仅空行、无 Component 头、头行即 EOF。
    #[test]
    fn identity_degenerate_inputs() {
        for db in [
            "",
            "\n",
            "not a component\nComponent 1\nACGT\n",
            "Component 5",
            "Component 5\n",
            "Component\n\n\n",
        ] {
            for reads in ["", "\n", "1\t>a\t50%\tACGT", "1\n", "1\t>a"] {
                let da = tmpdir("p3_t6_degen_a");
                let db_ = tmpdir("p3_t6_degen_b");
                let la = ref_partition(db, reads, &da, &PartParams::default()).unwrap();
                let lb = partition(db, reads, &db_, &PartParams::default()).unwrap();
                assert_trees_identical(&da, &db_);
                assert_listing_identical(&la, &lb, &da, &db_);
                let _ = std::fs::remove_dir_all(&da);
                let _ = std::fs::remove_dir_all(&db_);
            }
        }
    }

    /// 性能探针（`--ignored` 手跑）：合成 ~1/3 dev1 规模输入，身份对拍 +
    /// 双方耗时。dev1 量级：6.4k 组件图文本 + 7.6M 行 reads（~810MB 输出）。
    #[test]
    #[ignore = "大输入性能探针：cargo test -p trinity-chrysalis --release -- --ignored partition_perf_probe"]
    fn partition_perf_probe() {
        let t0 = std::time::Instant::now();
        let mut rng = Lcg(7);
        // db：2100 组件，体 170..249 行（~91% 过 176 行阈值 → 与 dev1 一样
        // 绝大多数组件登记在册），每行 ~26B
        let mut db = String::new();
        let ids: Vec<u64> = (0..2100).map(|i| (i * 37) % 2100).collect();
        for &id in &ids {
            db.push_str(&format!("Component {id}\n"));
            let n = 170 + rng.below(80) as usize;
            for i in 0..n {
                db.push_str(&format!("{i}\t-1\t1\t{}\t1\n", "ACGTACGTAC".repeat(2)));
            }
        }
        // reads：2.5M 行 × ~250B ≈ 600MB，**按 comp 分组**（镜像
        // readsToComponents.out.sort 的真实契约：dev1 均值 ~1200 行/comp）
        let seq = "A".repeat(220);
        let mut reads = String::with_capacity(2_500_000 * 260);
        let mut rn = 0u64;
        for &c in &ids {
            let reps = 800 + rng.below(800);
            for _ in 0..reps {
                reads.push_str(&format!("{c}\t>read{rn}/1\t50%\t{seq}\n"));
                rn += 1;
            }
        }
        println!(
            "gen {}s (db {}MB, reads {}MB)",
            t0.elapsed().as_secs_f32(),
            db.len() / 1_000_000,
            reads.len() / 1_000_000
        );

        let p = PartParams::default();
        let da = tmpdir("p3_t6_probe_ref");
        let t1 = std::time::Instant::now();
        let la = ref_partition(&db, &reads, &da, &p).unwrap();
        let t_ref = t1.elapsed();
        let db_ = tmpdir("p3_t6_probe_new");
        let t2 = std::time::Instant::now();
        let lb = partition(&db, &reads, &db_, &p).unwrap();
        let t_new = t2.elapsed();
        assert_trees_identical(&da, &db_);
        assert_listing_identical(&la, &lb, &da, &db_);
        println!(
            "reference: {:.2}s | new: {:.2}s | speedup {:.1}x | files {}",
            t_ref.as_secs_f32(),
            t_new.as_secs_f32(),
            t_ref.as_secs_f32() / t_new.as_secs_f32().max(1e-6),
            la.len()
        );
        let _ = std::fs::remove_dir_all(&da);
        let _ = std::fs::remove_dir_all(&db_);
    }
}
