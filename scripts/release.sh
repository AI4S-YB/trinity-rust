#!/bin/bash
# trinity-rust 发布打包脚本
# 产物: dist/trinity-rust-<ver>-<git短哈希>-<目标三元组>.tar.gz
set -euo pipefail
cd "$(dirname "$0")/.."

VER=$(grep -m1 '^version' Cargo.toml | sed 's/[^0-9.]//g')
HASH=$(git rev-parse --short HEAD)
TRIPLE=$(rustc -vV | grep '^host:' | awk '{print $2}')
NAME="trinity-rust-${VER}-${HASH}-${TRIPLE}"
DIST="dist/${NAME}"

echo "== 打包 ${NAME} =="
rm -rf "$DIST" && mkdir -p "$DIST/bin" "$DIST/docs"

cargo build --release --workspace

# 主二进制 + 阶段工具（交叉验证/手工调试用）
for b in trinity-cli inchworm trinity-kmer trinity-chrysalis butterfly; do
  cp "target/release/${b}" "$DIST/bin/"
  strip "$DIST/bin/${b}"
done

cp README.md "$DIST/"
cp docs/porting-map.md docs/benchmarks.md docs/backlog.md \
   docs/architecture-discovery.md \
   docs/athaliana-benchmark.md docs/athaliana-benchmark-hd.md \
   docs/rep10-report.md "$DIST/docs/" 2>/dev/null || true

cat > "$DIST/VERSION" <<EOF
trinity-rust ${VER} (commit ${HASH})
target: ${TRIPLE}
built: $(date -u +%Y-%m-%dT%H:%M:%SZ)
rustc: $(rustc --version | head -1)

主命令: bin/trinity-cli --seqType fq --left R1.fq.gz --right R2.fq.gz \
        --CPU N --max_memory XG --output outdir
与原版 Trinity v2.15.2 的差异与验证等级: 见 docs/（architecture-discovery.md,
athaliana-benchmark*.md, rep10-report.md）与仓库 spec/plans。
阶段工具为兼容验证保留, 主流程无需手动调用。
EOF

cd dist
tar -czf "${NAME}.tar.gz" "${NAME}"
sha256sum "${NAME}.tar.gz" | tee "${NAME}.tar.gz.sha256"
cd ..
echo "== 完成: dist/${NAME}.tar.gz ($(du -h dist/${NAME}.tar.gz | cut -f1)) =="
