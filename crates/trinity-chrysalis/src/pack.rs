//! per-Cbin 打包存储（IO 优化）——替代 per-component 小文件。
//!
//! 背景: classic 管线在 `Component_bins/Cbin*/` 下每组件写 3 个持久文件
//! （`cN.graph.out` / `cN.graph.reads` / `cN.graph.allProbPaths.fasta`），
//! 2.8 万组件 ≈ 8.5 万个小文件；输出目录在 NFS 上时每个 create/open/close
//! 都是元数据往返，实测贡献了大量 sys time。
//!
//! 设计（从简）:
//! - 每个 Cbin 两个逻辑包: `graphs`（QuantifyGraph 产物 graph.out+graph.reads）
//!   与 `butterfly`（allProbPaths.fasta）; 各为 `<name>.pack`（记录顺序追加的
//!   裸字节）+ `<name>.idx`（定长 32B 记录: id u64 / kind u64 / offset u64 /
//!   len u64, 全小端）。
//! - 索引读入为 `Vec<PackEntry>` 按 (id, kind) 排序后二分——u64 数组级别,
//!   不建 HashMap<String, _>; 写侧全程追加（顺序写）。
//! - 幂等: 各阶段（chrysalis-quantify / butterfly 池）都有整段 checkpoint,
//!   阶段重跑 =  truncate 重建对应的包, 与原来"重写同名文件"语义一致。
//! - partition 的 `.graph.tmp`/`.reads.tmp` 仍是逐组件文件（quantify 成功即删,
//!   生命周期限于 chrysalis 阶段内; --scratch 下本就在本地盘）——有意不打包。

use std::io::{BufWriter, Write};
use std::path::Path;

use trinity_common::error::CommonError;

/// graph.out 记录（QuantifyGraph 的图输出）。
pub const KIND_GRAPH_OUT: u64 = 1;
/// graph.reads 记录（QuantifyGraph 的 reads 输出）。
pub const KIND_GRAPH_READS: u64 = 2;
/// allProbPaths.fasta 记录（Butterfly 产物）。
pub const KIND_ALL_PROB_PATHS: u64 = 3;

/// QuantifyGraph 产物包名（每个 Cbin 一个）。
pub const GRAPHS_PACK: &str = "graphs";
/// Butterfly 产物包名（每个 Cbin 一个）。
pub const BUTTERFLY_PACK: &str = "butterfly";

/// 索引记录字节数: id / kind / offset / len 各 u64 LE。
const IDX_REC_LEN: usize = 32;

/// 追加写包: 数据顺序写 `<name>.pack`, 索引顺序写 `<name>.idx`（定长记录）。
pub struct PackWriter {
    data: BufWriter<std::fs::File>,
    idx: BufWriter<std::fs::File>,
    offset: u64,
}

impl PackWriter {
    /// 在 `dir` 下创建（truncate）`<name>.pack` / `<name>.idx`。
    /// 阶段级幂等点: 重跑即重建, 等价于旧布局的重写同名文件。
    pub fn create(dir: &Path, name: &str) -> Result<Self, CommonError> {
        let data = std::fs::File::create(dir.join(format!("{name}.pack")))?;
        let idx = std::fs::File::create(dir.join(format!("{name}.idx")))?;
        Ok(PackWriter {
            data: BufWriter::new(data),
            idx: BufWriter::new(idx),
            offset: 0,
        })
    }

    /// 追加一条记录（id, kind）→ 当前偏移。
    pub fn append(&mut self, id: u64, kind: u64, bytes: &[u8]) -> Result<(), CommonError> {
        self.data.write_all(bytes)?;
        let mut rec = [0u8; IDX_REC_LEN];
        rec[0..8].copy_from_slice(&id.to_le_bytes());
        rec[8..16].copy_from_slice(&kind.to_le_bytes());
        rec[16..24].copy_from_slice(&self.offset.to_le_bytes());
        rec[24..32].copy_from_slice(&(bytes.len() as u64).to_le_bytes());
        self.idx.write_all(&rec)?;
        self.offset += bytes.len() as u64;
        Ok(())
    }

    /// 冲刷两侧缓冲（drop 前显式调用以传播写错误）。
    pub fn finish(&mut self) -> Result<(), CommonError> {
        self.data.flush()?;
        self.idx.flush()?;
        Ok(())
    }
}

/// 内存索引: 定长记录读入 Vec 后按 (id, kind) 排序, 二分查找。
pub struct PackIndex {
    entries: Vec<PackEntry>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct PackEntry {
    id: u64,
    kind: u64,
    offset: u64,
    len: u64,
}

impl PackIndex {
    /// 读入 `<name>.idx` 全文并排序（32B/条; 2.8 万组件 × 3 条 ≈ 2.7MB）。
    pub fn load(path: &Path) -> Result<Self, CommonError> {
        let raw = std::fs::read(path)?;
        if raw.len() % IDX_REC_LEN != 0 {
            return Err(CommonError::Parse(format!(
                "pack index {} has trailing {} byte(s) (corrupt?)",
                path.display(),
                raw.len() % IDX_REC_LEN
            )));
        }
        let mut entries = Vec::with_capacity(raw.len() / IDX_REC_LEN);
        for rec in raw.chunks_exact(IDX_REC_LEN) {
            entries.push(PackEntry {
                id: u64::from_le_bytes(rec[0..8].try_into().unwrap()),
                kind: u64::from_le_bytes(rec[8..16].try_into().unwrap()),
                offset: u64::from_le_bytes(rec[16..24].try_into().unwrap()),
                len: u64::from_le_bytes(rec[24..32].try_into().unwrap()),
            });
        }
        entries.sort_unstable_by_key(|e| (e.id, e.kind));
        Ok(PackIndex { entries })
    }

    /// 查 (id, kind) → (offset, len); 同 id+kind 多条（不应出现）取首条。
    pub fn find(&self, id: u64, kind: u64) -> Option<(u64, u64)> {
        self.entries
            .binary_search_by_key(&(id, kind), |e| (e.id, e.kind))
            .ok()
            .map(|i| (self.entries[i].offset, self.entries[i].len))
    }
}

/// 从已打开的 pack 文件按偏移读一条记录（read_exact_at, 与并发 append 兼容）。
pub fn read_at(file: &std::fs::File, offset: u64, len: u64) -> Result<Vec<u8>, CommonError> {
    use std::os::unix::fs::FileExt;
    let mut buf = vec![0u8; len as usize];
    file.read_exact_at(&mut buf, offset)?;
    Ok(buf)
}

/// 一次性读取便捷入口（递归模式 mini 管线等小组件数场景）:
/// 加载 `dir/<name>.idx` + 打开 `<name>.pack` 读一条记录。记录缺失 → Err。
pub fn read_entry(dir: &Path, name: &str, id: u64, kind: u64) -> Result<Vec<u8>, CommonError> {
    let idx = PackIndex::load(&dir.join(format!("{name}.idx")))?;
    let (offset, len) = idx.find(id, kind).ok_or_else(|| {
        CommonError::Parse(format!(
            "pack {}/{} has no record (id={id}, kind={kind})",
            dir.display(),
            name
        ))
    })?;
    let file = std::fs::File::open(dir.join(format!("{name}.pack")))?;
    read_at(&file, offset, len)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn tmpdir(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("pack_{name}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn write_load_find_roundtrip() {
        let d = tmpdir("roundtrip");
        let mut w = PackWriter::create(&d, "graphs").unwrap();
        // 乱序 id + 多 kind + 空记录
        w.append(7, KIND_GRAPH_OUT, b"g7").unwrap();
        w.append(3, KIND_GRAPH_OUT, b"g33").unwrap();
        w.append(7, KIND_GRAPH_READS, b"r7!!").unwrap();
        w.append(3, KIND_GRAPH_READS, b"").unwrap();
        w.finish().unwrap();
        drop(w);

        let idx = PackIndex::load(&d.join("graphs.idx")).unwrap();
        assert_eq!(idx.find(7, KIND_GRAPH_OUT), Some((0, 2)));
        assert_eq!(idx.find(3, KIND_GRAPH_OUT), Some((2, 3)));
        assert_eq!(idx.find(7, KIND_GRAPH_READS), Some((5, 4)));
        assert_eq!(idx.find(3, KIND_GRAPH_READS), Some((9, 0)));
        assert_eq!(idx.find(9, KIND_GRAPH_OUT), None);

        let f = std::fs::File::open(d.join("graphs.pack")).unwrap();
        assert_eq!(read_at(&f, 0, 2).unwrap(), b"g7");
        assert_eq!(read_at(&f, 2, 3).unwrap(), b"g33");
        assert_eq!(read_at(&f, 5, 4).unwrap(), b"r7!!");
        assert_eq!(read_at(&f, 9, 0).unwrap(), b"");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn create_truncates_existing_pack() {
        let d = tmpdir("truncate");
        let mut w = PackWriter::create(&d, "butterfly").unwrap();
        w.append(1, KIND_ALL_PROB_PATHS, b"aaaa").unwrap();
        w.finish().unwrap();
        drop(w);
        // 重建（阶段重跑语义）: 旧记录不可见
        let mut w = PackWriter::create(&d, "butterfly").unwrap();
        w.append(2, KIND_ALL_PROB_PATHS, b"b").unwrap();
        w.finish().unwrap();
        drop(w);
        let idx = PackIndex::load(&d.join("butterfly.idx")).unwrap();
        assert_eq!(idx.find(1, KIND_ALL_PROB_PATHS), None);
        assert_eq!(idx.find(2, KIND_ALL_PROB_PATHS), Some((0, 1)));
        assert_eq!(
            std::fs::metadata(d.join("butterfly.pack")).unwrap().len(),
            1
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn read_entry_helper_and_missing_record() {
        let d = tmpdir("entry");
        let mut w = PackWriter::create(&d, "graphs").unwrap();
        w.append(42, KIND_GRAPH_OUT, b"graph-bytes").unwrap();
        w.finish().unwrap();
        drop(w);
        assert_eq!(
            read_entry(&d, "graphs", 42, KIND_GRAPH_OUT).unwrap(),
            b"graph-bytes"
        );
        assert!(read_entry(&d, "graphs", 42, KIND_GRAPH_READS).is_err());
        let _ = std::fs::remove_dir_all(&d);
    }
}
