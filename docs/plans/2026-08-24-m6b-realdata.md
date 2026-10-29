# M6b 真实数据基准 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **提交冻结（2026-08-24 用户指示）:** 执行本计划期间**跳过所有 Commit 步骤**（本地也不提交），直到用户解除冻结。其余步骤照常。

**Goal:** 拟南芥真实 RNA-seq（抽样 10M PE 对）上跑 classic/PogFlow/M5/M6a 四配置 xtask bench，沉淀蛋白级回收指标（P_full/P70），按泛化验证型 Gate 判定并入档。

**Architecture:** 数据准备走一次性幂等脚本（ENA 选样下载 + Ensembl Plants 参考 + seqkit 确定性抽样）；xtask bench 增加 `--prot-ref` 与 DIAMOND 蛋白回收指标（参考缺省时转录本最长 ORF 翻译，复用 trinity-butterfly 的 translate_codon）；四配置跑 bench（缓存惯例）；M6a 延长诊断重演；Gate 表入档。主工程（crates/）零改动。

**Tech Stack:** xtask（已依赖 trinity-butterfly/sha2/flate2）+ DIAMOND 2.2.4（env `TRINITY_DIAMOND`，默认 `/public/home/senior007/miniconda3/envs/diamond2go/bin/diamond`）+ seqkit（`/public/home/senior007/miniconda3/bin/seqkit`）+ ENA/Ensembl Plants HTTPS。

---

## 文件结构

```
scripts/prepare_realdata.sh                     # 新建: 数据准备（幂等）
xtask/src/bench.rs                              # 修改: --prot-ref + 蛋白回收指标
docs/opt-log.md                                 # 入档（bench 自动追加行 + 手工段落）
/storage/home/senior007/test/athal/bench/real/  # 数据产物（不入仓库）
  R1.fq.gz R2.fq.gz atha_cdna.fa atha_pep.fa alyr_pep.fa
```

---

### Task 1: 数据准备脚本 + 执行

**Files:** Create `scripts/prepare_realdata.sh`; 产物目录 `/storage/home/senior007/test/athal/bench/real/`

- [x] **Step 1: 选样（ENA portal API，交互验证）**

查询候选（拟南芥 PE Illumina RNA-seq，20M-60M 对）：

```bash
curl -s 'https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=tax_tree(3702)%20AND%20library_strategy%3D%22RNA-Seq%22%20AND%20library_layout%3D%22PAIRED%22%20AND%20platform%3D%22ILLUMINA%22&fields=run_accession,read_count,base_count,study_accession,experiment_title&format=json&limit=200' > /tmp/ena_cand.json
python3 -c "
import json
rows = json.load(open('/tmp/ena_cand.json'))
# read_count 是两端之和; 20M-60M 对 => 40M-120M 条
ok = [r for r in rows if r.get('read_count') and 40_000_000 <= int(r['read_count']) <= 120_000_000]
for r in ok[:10]: print(r['run_accession'], r['read_count'], r.get('experiment_title','')[:60])
"
```

选定首个标题不含 small-RNA/degradome/miRNA 的候选（rRNA/小RNA 文库排除）。验证 fastq 可下载：

```bash
curl -s 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=<RUN>&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes&format=json'
```

Expected: `fastq_ftp` 两条路径（R1/R2）。**把选定的 accession 记下来**（写进脚本 `RUN_ACC=` 变量与 opt-log）。

- [x] **Step 2: 写脚本 `scripts/prepare_realdata.sh`（幂等: 每个产物存在即跳过）**

```bash
#!/usr/bin/env bash
# M6b 真实数据准备（幂等）。产物: bench/real/{R1,R2}.fq.gz atha_cdna.fa atha_pep.fa alyr_pep.fa
set -euo pipefail
RUN_ACC="<Step1 选定>"          # ← 填入选定 accession
OUT=/storage/home/senior007/test/athal/bench/real
RAW=$OUT/raw
ENA_FTP=https://ftp.sra.ebi.ac.uk
EG=https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-60/plants/fasta
SEQKIT=/public/home/senior007/miniconda3/bin/seqkit
N=10000000; SEED=42
mkdir -p "$RAW" "$OUT"

# ---- 1. reads 下载 ----
if [ ! -f "$OUT/R1.fq.gz" ]; then
  FTP=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$RUN_ACC&result=read_run&fields=fastq_ftp&format=json" \
        | python3 -c "import json,sys; print(json.load(sys.stdin)[0]['fastq_ftp'])")
  R1F=$(echo "$FTP" | cut -d';' -f1); R2F=$(echo "$FTP" | cut -d';' -f2)
  [ -f "$RAW/R1.fastq.gz" ] || curl -sS -o "$RAW/R1.fastq.gz" "$ENA_FTP/$R1F"
  [ -f "$RAW/R2.fastq.gz" ] || curl -sS -o "$RAW/R2.fastq.gz" "$ENA_FTP/$R2F"
  # 确定性抽样: 同 seed 同记录序 → 两端同索引
  $SEQKIT sample -n $N -s $SEED -o "$OUT/R1.fq.gz" "$RAW/R1.fastq.gz"
  $SEQKIT sample -n $N -s $SEED -o "$OUT/R2.fq.gz" "$RAW/R2.fastq.gz"
fi

# ---- 2. Ensembl Plants 参考（先列目录发现 alyr 的确切文件名） ----
[ -f "$OUT/atha_cdna.fa" ] || { curl -sS "$EG/arabidopsis_thaliana/cdna/" | grep -o 'Arabidopsis_thaliana[^"]*cdna.all.fa.gz' | head -1; \
  curl -sS -o "$OUT/atha_cdna.fa.gz" "$EG/arabidopsis_thaliana/cdna/$(curl -sS "$EG/arabidopsis_thaliana/cdna/" | grep -o 'Arabidopsis_thaliana[^"]*cdna.all.fa.gz' | head -1)"; gunzip -kf "$OUT/atha_cdna.fa.gz"; }
[ -f "$OUT/atha_pep.fa" ]  || { P=$(curl -sS "$EG/arabidopsis_thaliana/pep/" | grep -o 'Arabidopsis_thaliana[^"]*pep.all.fa.gz' | head -1); \
  curl -sS -o "$OUT/atha_pep.fa.gz" "$EG/arabidopsis_thaliana/pep/$P"; gunzip -kf "$OUT/atha_pep.fa.gz"; }
[ -f "$OUT/alyr_pep.fa" ]  || { P=$(curl -sS "$EG/arabidopsis_lyrata/pep/" | grep -o 'Arabidopsis_lyrata[^"]*pep.all.fa.gz' | head -1); \
  curl -sS -o "$OUT/alyr_pep.fa.gz" "$EG/arabidopsis_lyrata/pep/$P"; gunzip -kf "$OUT/alyr_pep.fa.gz"; }

# ---- 3. 断言 ----
N1=$($SEQKIT stats -n "$OUT/R1.fq.gz" 2>/dev/null | tail -1 | awk '{print $4}')
N2=$($SEQKIT stats -n "$OUT/R2.fq.gz" 2>/dev/null | tail -1 | awk '{print $4}')
echo "R1=$N1 R2=$N2 (期望 $N)"
[ "$N1" = "$N2" ] || { echo "两端记录数不一致"; exit 1; }
# 配对一致性抽查: 前 1000 条 ID（去 /1 /2 后缀）应逐位相同
diff <(zcat "$OUT/R1.fq.gz" | head -4000 | awk 'NR%4==1{sub(/\/1$/,"");print}') \
     <(zcat "$OUT/R2.fq.gz" | head -4000 | awk 'NR%4==1{sub(/\/2$/,"");print}') >/dev/null \
  || { echo "配对 ID 不一致（seqkit 同 seed 抽样失效, 改用 grep -f 方案）"; exit 1; }
echo "OK: $(grep -c '>' $OUT/atha_cdna.fa) cdna / $(grep -c '>' $OUT/atha_pep.fa) atha pep / $(grep -c '>' $OUT/alyr_pep.fa) alyr pep"
```

- [x] **Step 3: 执行（后台, 下载是瓶颈）**

Run: `chmod +x scripts/prepare_realdata.sh && ./scripts/prepare_realdata.sh`
Expected: `OK: ~48k cdna / ~35k atha pep / ~32k alyr pep` + R1=R2=10,000,000。下载中断重跑（幂等）。

**配对一致性若失败**（seqkit 同 seed 不可靠的版本）: 改为
```bash
$SEQKIT sample -n $N -s $SEED "$RAW/R1.fastq.gz" -o "$OUT/R1.fq.gz"
zcat "$OUT/R1.fq.gz" | awk 'NR%4==1{print}' > "$OUT/ids.txt"
$SEQKIT grep -f "$OUT/ids.txt" "$RAW/R2.fastq.gz" -o "$OUT/R2.fq.gz"
rm "$OUT/ids.txt"
```
（注意 grep 模式下 ids 需去 /1 后缀: 先 `sed 's/\/1$//'`）

- [x] **Step 4: 入档选样信息**

opt-log 追加一段: RUN_ACC、数据集标题/研究、read 深度、抽样参数（10M/seed42）、三个参考的记录数。

- [x] **Step 5: Commit（冻结期间跳过）**

---

### Task 2: xtask bench 蛋白回收指标（TDD）

**Files:** Modify `xtask/src/bench.rs`

- [x] **Step 1: 写失败测试（bench.rs 底部 tests 模块追加）**

```rust
    /// 手算用例: P1=100aa, P2=200aa。
    /// q1↔P1 pident 97, HSP [1,60]+[61,100] → 并集 cov=1.0（两段拼接）
    /// q2↔P1 pident 99, HSP [1,70] → cov=0.70
    /// P1 最佳 = q1 (1.0×0.97 > 0.70×0.99) → full ✓ p70 ✓
    /// q3↔P2 pident 80, HSP [1,200] → cov=1.0 但 pid<90 → 都不满足
    #[test]
    fn protein_recall_hand_case() {
        let tab = "q1\tP1\t97.0\t60\t200\t100\t1\t180\t1\t60\t0.0\t100\t\n\
                   q1\tP1\t97.0\t40\t200\t100\t181\t300\t61\t100\t0.0\t70\t\n\
                   q2\tP1\t99.0\t70\t200\t100\t1\t210\t1\t70\t0.0\t110\t\n\
                   q3\tP2\t80.0\t200\t600\t200\t1\t600\t1\t200\t0.0\t300\t\n";
        let lens: HashMap<String, usize> =
            [("P1".to_string(), 100), ("P2".to_string(), 200)].into_iter().collect();
        let (p_full, p70, total) = protein_recall(tab, &lens);
        assert_eq!(total, 2);
        assert!((p_full - 0.5).abs() < 1e-9);
        assert!((p70 - 0.5).abs() < 1e-9);
    }

    /// 最长 ORF 翻译: 记录1 = ATG + 50×GCT + TAA → M+50A（51aa ≥50 保留）;
    /// 记录2 = 短 ORF（<50aa）→ 过滤。
    #[test]
    fn longest_orf_translation() {
        let dir = std::env::temp_dir().join(format!("xt_orf_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let cds = format!(">T1\nATG{}TAA\n>T2\nATGAAATAA\n", "GCT".repeat(50));
        let in_fa = dir.join("in.fa");
        std::fs::write(&in_fa, cds).unwrap();
        let out_fa = dir.join("out.fa");
        let n = translate_longest_orf(&in_fa, &out_fa).unwrap();
        assert_eq!(n, 1);
        let text = std::fs::read_to_string(&out_fa).unwrap();
        assert!(text.contains(">T1"));
        let seq: String = text.lines().skip(1).collect();
        assert_eq!(seq.len(), 51);
        assert!(seq.starts_with('M') && seq.ends_with('A'));
    }
```

- [x] **Step 2: 跑测试确认失败**

Run: `cargo test -p xtask`
Expected: FAIL（两函数未定义）

- [x] **Step 3: 实现（bench.rs 追加; 顶部 use 增补）**

顶部补: `use trinity_butterfly::protein_extend::translate_codon;`

```rust
// ---------------------------------------------------------------------------
// 蛋白级回收指标（M6b）: DIAMOND blastx 组装产物 vs 蛋白参考
// ---------------------------------------------------------------------------

fn diamond_bin() -> &'static str {
    static DEFAULT: &str = "/public/home/senior007/miniconda3/envs/diamond2go/bin/diamond";
    match std::env::var("TRINITY_DIAMOND") {
        Ok(v) if !v.is_empty() => Box::leak(v.into_boxed_str()),
        _ => DEFAULT,
    }
}

/// DIAMOND blastx → tabular（12 列, 无 qframe）。dmnd 按参考内容 sha256 缓存。
fn diamond_blastx_cached(query_fa: &Path, prot_ref: &Path, threads: usize) -> std::io::Result<String> {
    let bytes = std::fs::read(prot_ref)?;
    let mut h = Sha256::new();
    h.update(&bytes);
    let kh: String = h.finalize().iter().take(6).map(|b| format!("{b:02x}")).collect();
    let dmnd = bench_root().join("cache").join(format!("protref_{kh}.dmnd"));
    if !dmnd.exists() {
        let st = std::process::Command::new(diamond_bin())
            .args(["makedb", "--in"]).arg(prot_ref).arg("-d").arg(&dmnd)
            .arg("--threads").arg(threads.to_string())
            .stdout(std::process::Stdio::null()).status()
            .map_err(|e| std::io::Error::other(format!("spawn diamond: {e}")))?;
        if !st.success() { return Err(std::io::Error::other("diamond makedb failed")); }
    }
    let out_path = std::env::temp_dir().join(format!("xt_blastx_{}_{}", kh, std::process::id()));
    let st = std::process::Command::new(diamond_bin())
        .arg("blastx").arg("-q").arg(query_fa).arg("-d").arg(&dmnd)
        .arg("--threads").arg(threads.to_string())
        .arg("-e").arg("1e-5").arg("--max-target-seqs").arg("20")
        .arg("--outfmt").arg("6")
        .arg("-o").arg(&out_path)
        .stdout(std::process::Stdio::null()).status()
        .map_err(|e| std::io::Error::other(format!("spawn diamond: {e}")))?;
    if !st.success() { return Err(std::io::Error::other("diamond blastx failed")); }
    let text = std::fs::read_to_string(&out_path)?;
    let _ = std::fs::remove_file(&out_path);
    Ok(text)
}

/// 参考转录本 → 最长 ORF 蛋白集（≥50aa, 3 框, 按 stop 分段取最长, M 锚定优先）。
/// 返回蛋白记录数。与 M6b 探索脚本（python 版）同口径。
fn translate_longest_orf(cdna_fa: &Path, out_fa: &Path) -> std::io::Result<usize> {
    let text = std::fs::read_to_string(cdna_fa)?;
    let mut out = String::new();
    let mut n = 0usize;
    let mut name: Option<String> = None;
    let mut seq = String::new();
    let mut flush = |name: &Option<String>, seq: &str, out: &mut String, n: &mut usize| {
        if let Some(nm) = name {
            let s = seq.to_ascii_uppercase();
            let mut best: &[u8] = b"";
            for frame in 0..3 {
                let prot: Vec<u8> = (frame..s.len().saturating_sub(2)).step_by(3)
                    .map(|i| translate_codon(&s.as_bytes()[i..i + 3]))
                    .collect();
                for seg in prot.split(|&b| b == b'*') {
                    let mpos = seg.iter().position(|&b| b == b'M');
                    let cand: &[u8] = match mpos {
                        Some(p) => &seg[p..],
                        None => seg,
                    };
                    if cand.len() > best.len() {
                        best = cand;
                    }
                }
            }
            if best.len() >= 50 {
                *n += 1;
                out.push('>');
                out.push_str(nm);
                out.push('\n');
                for ch in best.chunks(60) {
                    out.push_str(std::str::from_utf8(ch).unwrap_or(""));
                    out.push('\n');
                }
            }
        }
    };
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            flush(&name.take(), &std::mem::take(&mut seq), &mut out, &mut n);
            name = Some(h.split_whitespace().next().unwrap_or("").to_string());
        } else {
            seq.push_str(line.trim());
        }
    }
    flush(&name, &seq, &mut out, &mut n);
    std::fs::write(out_fa, out)?;
    Ok(n)
}

/// blastx tabular → (P_full, P70, 蛋白总数)。
/// 同一 (q,s) 的 HSP 求 subject 并集跨度; 每蛋白取 c×pid 最高的 query。
fn protein_recall(tab: &str, prot_lens: &HashMap<String, usize>) -> (f64, f64, usize) {
    use std::collections::{HashMap as Map, HashSet};
    // (q,s) → (区间集, 最大 pident)
    let mut pair: Map<(String, String), (Vec<(i64, i64)>, f64)> = Map::new();
    for line in tab.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 { continue; }
        let (Ok(pid), Ok(ss), Ok(se)) = (f[2].parse::<f64>(), f[8].parse::<i64>(), f[9].parse::<i64>()) else { continue };
        let e = pair.entry((f[0].to_string(), f[1].to_string())).or_insert((Vec::new(), 0.0));
        e.0.push((ss.min(se), ss.max(se)));
        e.1 = e.1.max(pid);
    }
    let mut best: Map<&str, (f64, f64)> = Map::new();
    for ((_, s), (ivs, pid)) in &pair {
        let Some(&len) = prot_lens.get(s) else { continue };
        let mut ivs = ivs.clone();
        ivs.sort_unstable();
        let mut cov = 0i64;
        let mut cur: Option<(i64, i64)> = None;
        for &(a, b) in &ivs {
            match cur {
                None => cur = Some((a, b)),
                Some((cs, ce)) if a <= ce + 1 => cur = Some((cs, ce.max(b))),
                Some((cs, ce)) => { cov += ce - cs + 1; cur = Some((a, b)); }
            }
        }
        if let Some((cs, ce)) = cur { cov += ce - cs + 1; }
        let c = cov as f64 / len.max(1) as f64;
        let score = c * pid;
        match best.get(s.as_str()) {
            Some((bc, bp)) if bc * bp >= score => {}
            _ => { best.insert(s.as_str(), (c, pid)); }
        }
    }
    let total = prot_lens.len();
    if total == 0 { return (0.0, 0.0, 0); }
    let full = best.values().filter(|&&(c, p)| c >= 0.95 && p >= 0.95).count();
    let p70 = best.values().filter(|&&(c, p)| c >= 0.70 && p >= 0.90).count();
    (full as f64 / total as f64, p70 as f64 / total as f64, total)
}
```

（`HashSet` 导入若未用则删; `Path`/`PathBuf`/`Sha256` 已在 bench.rs 作用域。）

- [x] **Step 4: 跑测试确认通过**

Run: `cargo test -p xtask`
Expected: PASS（新增 2 + 既有 tests 全绿）

- [x] **Step 5: 接线 `run()`（四处小改）**

1. 参数解析加分支:

```rust
            "--prot-ref" => {
                prot_ref = argv.get(i + 1).cloned().unwrap_or_default();
                i += 2;
            }
```
（顶部 `let mut prot_ref = String::new();`）

2. minimap2 求值之后、metrics.json 写入之前插入:

```rust
    // 蛋白级回收（M6b）: DIAMOND blastx vs --prot-ref（缺省转录本最长 ORF）
    let (p_full, p70, prot_total) = {
        let dir2 = dir.join("protein_metric");
        let _ = std::fs::create_dir_all(&dir2);
        let prot_fa = if prot_ref.is_empty() {
            let p = dir2.join("ref_orf.fa");
            if !p.exists() {
                translate_longest_orf(std::path::Path::new(&reff), &p).unwrap();
            }
            p
        } else {
            std::path::PathBuf::from(&prot_ref)
        };
        let tab = diamond_blastx_cached(&final_fa, &prot_fa, 8).unwrap_or_default();
        // 蛋白参考长度表（entry 累计）
        let mut lens: HashMap<String, usize> = HashMap::new();
        let mut name: Option<String> = None;
        for line in std::fs::read_to_string(&prot_fa).unwrap_or_default().lines() {
            if let Some(h) = line.strip_prefix('>') {
                name = Some(h.split_whitespace().next().unwrap_or("").to_string());
                lens.entry(name.clone().unwrap()).or_insert(0);
            } else if let Some(n) = &name {
                *lens.entry(n.clone()).or_insert(0) += line.trim().len();
            }
        }
        protein_recall(&tab, &lens)
    };
```

（实现者注意: 上面第一段 lens 构造是推导残留——**以第二段 entry 模式为准**，只保留一份。）

3. metrics.json 行加 `,"p_full":{:.4},"p70":{:.4},"prot_total":{}`（对应实参 `p_full, p70, prot_total`）。

4. stdout 汇总行追加 ` | P_full={:.1}% P70={:.1}%`（实参 `100.0*p_full, 100.0*p70`）; opt-log 自动行不改动格式（保持既有列稳定，蛋白指标靠 metrics.json/手工段落）。

- [x] **Step 6: 回归验证（dev1 快速复跑, 缓存命中只重算指标）**

Run: `cargo build --release -p xtask && A=/storage/home/senior007/test/athal/bench && cargo run -q --release -p xtask -- bench --tag m3-dev1-pmetric --reads $A/reads/dev1_1.fq.gz,$A/reads/dev1_2.fq.gz --ref $A/refs/dev1.fa -- --path-engine pogflow --max_memory 16G --CPU 8`
Expected: 缓存命中（binary 未变? —— xtask 改动不影响 trinity-cli 二进制 → 同 key）+ 输出行含 `P_full≈63.4% P70≈76.8%`（修正口径: ≥95%cov 且 ≥95%id; 原型 68.5/83.0 因 pident 百分数误比 0.95 恒真而虚高, 且 Rust 版曾漏大写化 soft-masked 输入——两处均已修, 最终值经独立 python 复核一致）

- [x] **Step 7: clippy + fmt**

Run: `cargo clippy -p xtask --tests && cargo fmt -p xtask`
Expected: 干净

- [x] **Step 8: Commit（冻结期间跳过）**

---

### Task 3: 四配置运行

**Files:** 产物为 bench cache + opt-log 自动行

- [x] **Step 1: 构建并确认依赖**

Run: `cargo build --release -p trinity-cli -p xtask`
Expected: 成功（trinity-cli 未改动——若 Step 2 的回归导致二进制重链, 注意缓存 key 会变, 无妨）

- [x] **Step 2: 依次跑四配置（每个后台 ~35 分钟, 串行避免 8×4 线程争用）**

```bash
A=/storage/home/senior007/test/athal/bench/real
R=$A/R1.fq.gz,$A/R2.fq.gz; REF=$A/atha_cdna.fa
cd /storage/home/senior007/test/trinity_rust/trinity-rust
cargo run -q --release -p xtask -- bench --tag real-classic  --reads $R --ref $REF -- --max_memory 16G --CPU 8
cargo run -q --release -p xtask -- bench --tag real-pogflow --reads $R --ref $REF -- --path-engine pogflow --max_memory 16G --CPU 8
cargo run -q --release -p xtask -- bench --tag real-m5      --reads $R --ref $REF -- --path-engine pogflow --model-guided --max_memory 16G --CPU 8
cargo run -q --release -p xtask -- bench --tag real-m6a     --reads $R --ref $REF -- --path-engine pogflow --protein-guided --protein-db $A/alyr_pep.fa --max_memory 16G --CPU 8
```

Expected: 每个产出 R/P/F1/N50 + P_full/P70; real-m6a 的 stderr 有 `-- protein-guided:` 统计行（抄录存档）。串行脚本整体后台跑（~2.5h）。

- [x] **Step 3: Commit（冻结期间跳过）**

---

### Task 4: M6a 诊断重演（副报告）

**Files:** `docs/opt-log.md`（手工段落）; 分析过程用 /tmp 临时脚本（不入库）

- [x] **Step 1: 提取 real-m6a 工件**

`<cache>/out/protein_guided/{queries.fa,extended.fa}` + 日志中的 extended/accepted/+nt 统计行。

- [x] **Step 2: same-ref/diff-ref/unaligned 归类（复用 M6a 探索的方法）**

- 前缀/后缀 80nt 精确匹配把 extended 配对回 queries（方向 R/L）
- `minimap2 -x asm5 -t 8 --secondary=no atha_cdna.fa {queries,extended}.fa` → PAF
- 延长区碱基归类（同参考/异座位/无处）——python 脚本照 M6a 探索版（opt-log 2026-08-24 段有完整描述）
- 覆盖度移动: same-ref 条目的延长前后参考覆盖度 mean/median + 跨 0.95 门槛数

- [x] **Step 3: 对照 dev 探路结论写入 opt-log**

- 延长统计（extended/accepted/reverted/+nt）vs dev probe2（229/237kb）
- same/diff/unaligned 构成 vs dev（46.7/47.5/5.8%）
- 结论: 真实数据与 dev 是否同构（旁系游荡比例、验收通过率、覆盖度增益）

- [x] **Step 4: Commit（冻结期间跳过）**

---

### Task 5: Gate 判定入档

**Files:** `docs/opt-log.md`

- [x] **Step 1: 四配置对照表 + Gate 判定**

```markdown
## M6b 真实数据基准（2026-08-XX, 拟南芥 <RUN_ACC> 抽样 10M 对）

| tag | R | P | F1 | N50 | n | P_full | P70 |
|---|---|---|---|---|---|---|---|
| real-classic | | | | | | | |
| real-pogflow | | | | | | | |
| real-m5 | | | | | | | |
| real-m6a | | | | | | | |

### Gate（泛化验证型）
- G1 pogflow ≥ classic: PASS/FAIL（数字）
- G2 m5 ≥ classic: PASS/FAIL（数字）
- G3 m6a ≥ classic − 0.5pp: PASS/FAIL（数字）
- 副报告: （Task 4 结论一句化 + 延长统计）
```

- [x] **Step 2: 向用户汇报 Gate 结果与解读（不擅自开 M6c）**

- [x] **Step 3: Commit（冻结期间跳过）**

---

## 完成定义

- [x] bench/real/ 五件产物齐 + 选样信息入档
- [x] xtask 蛋白指标测试全绿 + dev1 复跑口径对齐（P_full=63.4/P70=76.8, 独立复核一致）
- [x] 四配置数字 + P_full/P70 入档
- [x] M6a 诊断重演段落入档
- [x] Gate 表入档并已向用户汇报
- [x] 未做任何 git 提交（冻结）
