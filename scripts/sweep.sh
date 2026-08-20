#!/bin/bash
# sweep.sh <reads_prefix> <ref> <engine> —— 网格: max_cov{off,50,100,200} × K{19,23,25,29}
# 例: ./scripts/sweep.sh dev1sub dev1 classic
set -e
PRE=$1; REF=$2; ENGINE=${3:-classic}
cd "$(dirname "$0")/.."
ROOT=${ATHAL:-/storage/home/senior007/test/athal}/bench
echo "| max_cov | K | n | R | P | F1 | N50 |"
echo "|---|---|---|---|---|---|---|"
for cov in off 50 100 200; do
  for k in 19 23 25 29; do
    EXTRA="--KMER_SIZE $k"
    if [ "$cov" = off ]; then EXTRA="$EXTRA --no_normalize_reads"; else EXTRA="$EXTRA --normalize_max_read_cov $cov"; fi
    [ "$ENGINE" = recursive ] && EXTRA="$EXTRA --recursive"
    OUT=$(cargo run -q --release -p xtask -- bench --tag sw-${ENGINE}-c${cov}-k${k} \
      --reads $ROOT/reads/${PRE}_1.fq.gz,$ROOT/reads/${PRE}_2.fq.gz --ref $ROOT/refs/$REF.fa \
      -- $EXTRA --max_memory 8G 2>&1 | grep '^  转录本' || true)
    # 行格式: "  转录本 n/total | R=8.2% R70=10.0% P=80.5% F1=14.9% | N50=2789"
    ROW=$(printf '%s' "$OUT" | sed -e 's/^  转录本 //' -e 's/ R70=[0-9.]*%//' \
      -e 's/ | R=\([0-9.]*\)% P=\([0-9.]*\)% F1=\([0-9.]*\)% | N50=\([0-9]*\)$/ | \1 | \2 | \3 | \4/')
    echo "| $cov | $k | $ROW |"
  done
done
