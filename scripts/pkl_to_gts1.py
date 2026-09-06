#!/usr/bin/env python
"""edge_support_raw.pkl → GTS1 二进制（trinity-butterfly gtdecode::EdgeSupport::load 的格式）。
仅开发期对拍用: 把 python 原型算好的原始 reads 边支持喂给 Rust 解码器,
在**同一组图**（68349cb4 缓存）上验证 Rust 解码 ≈ 原型 C 变体。

GTS1（小端）: magic "GTS1" | u64 n | n × (u64 cid, i32 u, i32 v, u32 cnt),
按 (cid,u,v) 升序。
"""
import pickle
import struct
import sys

IN_PKL = "/storage/home/senior007/test/athal/bench/export/edge_support_raw.pkl"
OUT_BIN = "/tmp/edge_support_raw.gts1"


def main():
    sup = pickle.load(open(IN_PKL, "rb"))
    print(f"pkl 条目 {len(sup)}")
    recs = sorted((int(cid), int(u), int(v), int(c)) for (cid, u, v), c in sup.items())
    with open(OUT_BIN, "wb") as f:
        f.write(b"GTS1")
        f.write(struct.pack("<Q", len(recs)))
        for cid, u, v, c in recs:
            f.write(struct.pack("<QiiI", cid, u, v, c))
    print(f"写出 {OUT_BIN}")


if __name__ == "__main__":
    main()
