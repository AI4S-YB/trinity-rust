#!/usr/bin/env bash
# 评估 v2 real 轨道 StringTie 参照构建（幂等）。
# 产物: bench/real/{guide.gtf, hisat2.bam, hisat2.log, real_stringtie.gtf, real_stringtie.abund}
#
# ---- 考古（Task5 Step1, 2026-08-26）----
# 当年 holdout 轨道参照 athal/stringtie.gtf 头部:
#   stringtie -G chr5.renamed.gtf -p 8 -o stringtie.gtf -A stringtie.abund chr5.bam
# 其中的 chr1(3行)/chr4(5行) 杂染色体记录来自 chr5.bam 自身的散逸比对
# （samtools idxstats chr5.bam: chr1=153, chr2=108, chr3=172, chr4=382 条）——
# 头部命令行诚实, 并不存在「用 all.bam 跑全基因组 StringTie」的一次运行。
# chr5.renamed.gtf 的生成: zcat TAIR10.gff.gz 取 NC_003076.8 子集 → `gffread -E - -o chr5.gtf`
# → NC_003076.8 改名 chr5（athal5.fa.fai 只有 chr1..chr5, 尺寸与 RefSeq 五条核染色体一一对应:
# NC_003070.9=chr1/30427671, NC_003071.7=chr2/19698289, NC_003074.8=chr3/23459830,
# NC_003075.7=chr4/18585056, NC_003076.8=chr5/26975502; 叶绿体 NC_000932.1/线粒体 NC_037304.1 不在 athal5.fa, 剔除）。
# dev2/dev3 轨道同构（bench/reads/devN.gtf 头部 `gffread -E - -o devN.gtf`; devN.st.gtf 头部
# `-G devN.renamed.gtf -p 8`, 无 -A）。
# → real 轨道全基因组引导 GTF = 同一 TAIR10.gff.gz 五条核染色体 → gffread → 同一改名约定。
# hisat2: 当年 hd_pipeline.sh L18/L38 为 `hisat2 -x athal5 -p 8 --dta -1 .. -2 .. | samtools sort -@ 4 -m 2G`;
# 本脚本按计划书给定形态执行（无 --dta）——差异在此注明, 其余参数一致。
set -euo pipefail

ATHAL=/storage/home/senior007/test/athal
R=$ATHAL/bench/real
H2=/public/home/senior007/miniconda3/envs/rnaseq/bin/hisat2
ST=/public/home/senior007/miniconda3/envs/stringtie-test/bin/stringtie
SAM=/public/home/senior007/miniconda3/envs/stringtie-test/bin/samtools
GFFREAD=/public/home/senior007/miniconda3/envs/busco/bin/gffread
mkdir -p "$R"

# ---- 1. 引导 GTF（镜像 chr5.renamed.gtf 的全基因组版）----
if [ ! -f "$R/guide.gtf" ]; then
  echo "[guide] gffread 转换 TAIR10.gff.gz 五条核染色体 ..."
  zcat "$ATHAL/TAIR10.gff.gz" \
    | awk '$1 ~ /^#/ || $1=="NC_003070.9"||$1=="NC_003071.7"||$1=="NC_003074.8"||$1=="NC_003075.7"||$1=="NC_003076.8"' \
    | "$GFFREAD" -E - -o "$R/guide.norename.gtf"
  # 改名对齐 athal5.fa（仅行首 seqid, 不用 sed -i: NFS 原地替换有前科）
  sed -e 's/^NC_003070\.9\t/chr1\t/' -e 's/^NC_003071\.7\t/chr2\t/' \
      -e 's/^NC_003074\.8\t/chr3\t/' -e 's/^NC_003075\.7\t/chr4\t/' \
      -e 's/^NC_003076\.8\t/chr5\t/' "$R/guide.norename.gtf" > "$R/guide.gtf"
  # 断言: guide 染色体名 ⊆ athal5.fa.fai
  BAD=$(grep -v '^#' "$R/guide.gtf" | cut -f1 | sort -u | comm -23 - <(cut -f1 "$ATHAL/athal5.fa.fai" | sort -u) | tr '\n' ',')
  [ -z "$BAD" ] || { echo "[guide] 染色体名不在 athal5.fa: $BAD"; exit 1; }
  echo "[guide] OK: $(grep -vc '^#' "$R/guide.gtf") 行; 特征类型分布: $(grep -v '^#' "$R/guide.gtf" | cut -f3 | sort | uniq -c | sort -rn | tr '\n' ';')"
fi

# ---- 2. HISAT2 比对（10M PE 对, ~10-20 min）----
if [ ! -f "$R/hisat2.bam" ]; then
  echo "[hisat2] 比对 real R1/R2 → athal5 索引 ..."
  if ! "$H2" -x "$ATHAL/athal5" -p 8 -1 "$R/R1.fq.gz" -2 "$R/R2.fq.gz" 2> "$R/hisat2.log" \
      | "$SAM" sort -@ 4 -T "$R/hisat2.sort.tmp" -o "$R/hisat2.tmp.bam" -; then
    rm -f "$R/hisat2.tmp.bam" "$R/hisat2.bam"
    echo "[hisat2] 失败, stderr 见 $R/hisat2.log"; exit 1
  fi
  mv "$R/hisat2.tmp.bam" "$R/hisat2.bam"
fi

# ---- 3. StringTie 引导组装（参照本体, ~5-10 min）----
if [ ! -f "$R/real_stringtie.gtf" ]; then
  echo "[stringtie] 引导组装 ..."
  if ! "$ST" -G "$R/guide.gtf" -p 8 -o "$R/real_stringtie.gtf" -A "$R/real_stringtie.abund" \
      "$R/hisat2.bam" 2> "$R/real_stringtie.log"; then
    rm -f "$R/real_stringtie.gtf"
    echo "[stringtie] 失败, stderr 见 $R/real_stringtie.log"; exit 1
  fi
fi

# ---- 4. sanity ----
echo "== sanity =="
echo "guide: $(grep -v '^#' "$R/guide.gtf" | grep -c $'\texon\t') exon 行, 染色体 $(grep -v '^#' "$R/guide.gtf" | cut -f1 | sort -u | tr '\n' ',')"
echo "hisat2 映射率: $(grep 'overall alignment rate' "$R/hisat2.log")"
echo "首条比对: $("$SAM" view "$R/hisat2.bam" | head -1 | cut -f1-6)"
echo "stringtie: $(grep -c $'\ttranscript\t' "$R/real_stringtie.gtf") 转录本 / $(grep $'\ttranscript\t' "$R/real_stringtie.gtf" | sed 's/ref_gene_id/XX/' | grep -o 'gene_id "[^"]*"' | sort -u | wc -l) 基因"
echo "stringtie stderr: $(cat "$R/real_stringtie.log" 2>/dev/null | head -5)"
echo "染色体分布:"; grep -v '^#' "$R/real_stringtie.gtf" | cut -f1 | sort | uniq -c
