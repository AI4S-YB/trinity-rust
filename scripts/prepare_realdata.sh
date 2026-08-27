#!/usr/bin/env bash
# M6b 真实数据准备（幂等）。产物: bench/real/{R1,R2}.fq.gz atha_cdna.fa atha_pep.fa alyr_pep.fa
#
# 网络现实: 出网带宽被限 (~16KB/s 单连接, ~150KB/s 并行总量), 全量 6.4GB 不可行。
# 策略: 并行分块只下载 fastq 前缀 (覆盖 >N*1.2 read 对), 截断 gzip 可完整解出前缀内
# 所有 record, 再对前缀做确定性抽样。对 benchmark 而言等价于真实数据子集。
set -euo pipefail
RUN_ACC="DRR095137"
OUT=/storage/home/senior007/test/athal/bench/real
RAW=$OUT/raw
ENA_FTP=https://ftp.sra.ebi.ac.uk
EG=https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-60/plants/fasta
SEQKIT=/public/home/senior007/miniconda3/bin/seqkit
N=10000000; SEED=42
THREADS=16; CHUNK=$((8*1024*1024))
mkdir -p "$RAW" "$OUT"

# ---- 0. 从 ENA API 拿路径/大小 ----
read -r R1F R2F B1 B2 RC < <(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$RUN_ACC&result=read_run&fields=fastq_ftp,fastq_bytes,read_count&format=json" \
  | python3 -c "
import json,sys
d=json.load(sys.stdin)[0]
f=d['fastq_ftp'].split(';'); b=d['fastq_bytes'].split(';')
strip=lambda p: p[len('ftp.sra.ebi.ac.uk/'):] if p.startswith('ftp.sra.ebi.ac.uk/') else p
print(strip(f[0]), strip(f[1]), b[0], b[1], d['read_count'])")

# ---- 1. 并行分块下载前缀 (只为覆盖 N 对 read, 不下全量) ----
# prefix_bytes = N * 1.2 * (bytes/read) ; read_count 为两条 fastq 的 read 总数
dl_prefix() { # $1=remote_path $2=total_bytes $3=dest.gz
  local path=$1 total=$2 dest=$3 base avg prefix nch i
  base=$(basename "$3")
  [ -f "$dest.ok" ] && { echo "$base: 已存在, 跳过"; return; }
  # 实测: read_count 为 spot 数 (R1 内 read 数 = RC), 每条 read 压缩 ~56B; 留 20% 余量
  prefix=$(( N * 12 / 10 * total / RC ))
  [ $prefix -gt "$total" ] && prefix=$total
  nch=$(( (prefix + CHUNK - 1) / CHUNK ))
  echo "$base: 需要前缀 $(($prefix/1024/1024))MB / 全量 $(($total/1024/1024))MB, $nch 块"
  mkdir -p "$RAW/chunks_$base"
  python3 - "$ENA_FTP/$path" "$RAW/chunks_$base" "$prefix" "$CHUNK" "$THREADS" <<'PYEOF'
import sys, os, time, urllib.request, concurrent.futures as cf
url, cdir, prefix, chunk, nthr = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
nch = (prefix + chunk - 1)//chunk
def get(i):
    ok=f"{cdir}/{i:05d}.ok"; dst=f"{cdir}/{i:05d}.part"
    s=i*chunk; e=min(s+chunk, prefix)-1; want=e-s+1
    if os.path.exists(ok) and os.path.exists(dst) and os.path.getsize(dst)==want: return 0
    for attempt in range(30):
        try:
            req=urllib.request.Request(url, headers={"Range": f"bytes={s}-{e}"})
            d=urllib.request.urlopen(req, timeout=300).read()
            if len(d)==want:
                with open(dst,"wb") as f: f.write(d)
                open(ok,"w").close()
                return len(d)
        except Exception:
            time.sleep(min(2**attempt, 60))
    raise RuntimeError(f"chunk {i} failed after 30 attempts")
done=0
with cf.ThreadPoolExecutor(nthr) as ex:
    for r in ex.map(get, range(nch)):
        done+=r
        print(f"  {done/1e6:.0f}MB / {prefix/1e6:.0f}MB", flush=True)
print("chunks done")
PYEOF
  # 按块号顺序拼接 (文件名 5 位零填充, 字典序即块序)
  cat $(ls "$RAW/chunks_$base"/*.part | sort) > "$dest"
  touch "$dest.ok"
  echo "$base: 前缀下载完成 $(($(stat -c%s "$dest")/1024/1024))MB"
}

if [ ! -f "$OUT/R1.fq.gz" ]; then
  dl_prefix "$R1F" "$B1" "$RAW/R1.fastq.gz"
  dl_prefix "$R2F" "$B2" "$RAW/R2.fastq.gz"
  # 截断 gzip: zcat 解出全部完整 record, 忽略尾部截断错误; 同 seed 抽样两端
  ( zcat "$RAW/R1.fastq.gz" 2>/dev/null || true ) | $SEQKIT sample -n $N -s $SEED -o "$OUT/R1.fq.gz"
  ( zcat "$RAW/R2.fastq.gz" 2>/dev/null || true ) | $SEQKIT sample -n $N -s $SEED -o "$OUT/R2.fq.gz"
fi

# ---- 2. Ensembl Plants 参考 ----
fetch_ref() { # $1=species-dir $2=pattern $3=outfile
  local F
  F=$(curl -sS "$EG/$1/" | grep -o "$2" | head -1)
  echo "Ensembl $1 file: $F"
  [ -n "$F" ] || { echo "目录里没找到 $2"; return 1; }
  curl -sS --retry 5 --retry-all-errors -C - -o "$3.gz" "$EG/$1/$F"
  gunzip -kf "$3.gz"
}
[ -f "$OUT/atha_cdna.fa" ] || fetch_ref arabidopsis_thaliana/cdna 'Arabidopsis_thaliana[^"]*cdna.all.fa.gz' "$OUT/atha_cdna.fa"
[ -f "$OUT/atha_pep.fa" ]  || fetch_ref arabidopsis_thaliana/pep  'Arabidopsis_thaliana[^"]*pep.all.fa.gz'  "$OUT/atha_pep.fa"
[ -f "$OUT/alyr_pep.fa" ]  || fetch_ref arabidopsis_lyrata/pep    'Arabidopsis_lyrata[^"]*pep.all.fa.gz'      "$OUT/alyr_pep.fa"

# ---- 3. 断言 ----
N1=$($SEQKIT stats -n "$OUT/R1.fq.gz" 2>/dev/null | tail -1 | awk '{print $4}')
N2=$($SEQKIT stats -n "$OUT/R2.fq.gz" 2>/dev/null | tail -1 | awk '{print $4}')
echo "R1=$N1 R2=$N2 (期望 $N)"
[ "$N1" = "$N2" ] || { echo "两端记录数不一致"; exit 1; }
if ! diff <(zcat "$OUT/R1.fq.gz" | head -4000 | awk 'NR%4==1{sub(/\/1$/,"");print}') \
          <(zcat "$OUT/R2.fq.gz" | head -4000 | awk 'NR%4==1{sub(/\/2$/,"");print}') >/dev/null; then
  # 同 seed 抽样索引不一致时回退: 以 R1 的 ID 为准过滤 R2
  echo "配对 ID 不一致, 回退 seqkit grep 配对"
  zcat "$OUT/R1.fq.gz" | awk 'NR%4==1{print substr($1,1,length($1)-2)}' > "$RAW/ids.txt"
  zcat "$OUT/R2.fastq.gz" 2>/dev/null | $SEQKIT grep -f "$RAW/ids.txt" -o "$OUT/R2.fq.gz"
  N2=$($SEQKIT stats -n "$OUT/R2.fq.gz" 2>/dev/null | tail -1 | awk '{print $4}')
  [ "$N1" = "$N2" ] || { echo "回退后仍不一致 R1=$N1 R2=$N2"; exit 1; }
fi
echo "OK: $(grep -c '>' $OUT/atha_cdna.fa) cdna / $(grep -c '>' $OUT/atha_pep.fa) atha pep / $(grep -c '>' $OUT/alyr_pep.fa) alyr pep"
