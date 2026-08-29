#!/usr/bin/env python
"""gtasm 数值对拍（验收 #2）: torch 重载 gtw 在同一 graph.out 上重算边概率,
与 Rust `gtasm_score_dump` 的输出逐边对比。

用法: gt_score_compare.py <graph.out> <weights.gtw> <rust_dump.txt>
容差: 逐边绝对差 ≤ 1e-4（f32 累加顺序差异量级）。
"""
import math
import struct
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, "/storage/home/senior007/test/athal/bench/export")
import torch

from gt_train3 import GTNet


def load_weights(path):
    """GTW1 纯 python 重载（与 export_gt_weights.load_weights 同）。"""
    out = {}
    with open(path, "rb") as f:
        assert f.read(4) == b"GTW1"
        (n,) = struct.unpack("<I", f.read(4))
        for _ in range(n):
            (nl,) = struct.unpack("<I", f.read(4))
            name = f.read(nl).decode()
            (nd,) = struct.unpack("<I", f.read(4))
            shape = struct.unpack(f"<{nd}I", f.read(4 * nd))
            cnt = int(np.prod(shape))
            arr = np.frombuffer(f.read(4 * cnt), dtype="<f4").reshape(shape)
            out[name] = arr.copy()
    return out


def load_comp(path):
    """单文件 load_comps 语义（gt_train3.load_comps 的逐行逻辑）。"""
    names, cnt, out = {}, {}, defaultdict(list)
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            to, frm, c = int(p[0]), int(p[1]), int(p[2])
            names[to] = p[3]
            cnt[to] = max(cnt.get(to, 0), c)
            out[frm].append(to)
    return names, cnt, out


def build_feats(names, cnt, out):
    """probe_e2e.build_all 的纯 numpy 复刻（不依赖 dgl——度从边表数）。"""
    vids = sorted(set(list(cnt) + list(out)))
    remap = {v: i for i, v in enumerate(vids)}
    src, dst = [], []
    for f, ts in out.items():
        for t in ts:
            if f in remap and t in remap:
                src.append(remap[f])
                dst.append(remap[t])
    n = len(vids)
    ind = np.zeros(n, dtype=np.float32)
    outd = np.zeros(n, dtype=np.float32)
    for a, b in zip(src, dst):
        outd[a] += 1
        ind[b] += 1
    logc = np.array([math.log1p(cnt.get(v, 0)) for v in vids], dtype=np.float32)
    x = torch.tensor(np.stack([logc, np.log1p(ind), np.log1p(outd)], 1), dtype=torch.float32)
    es, ed = logc[np.array(src)], logc[np.array(dst)]
    ev = np.array([math.log1p(cnt.get(vids[b], 1)) for a, b in zip(src, dst)], dtype=np.float32)
    e = torch.tensor(np.stack([es, ed, es - ed, ev], 1), dtype=torch.float32)
    return vids, src, dst, x, e


def main():
    graph_path, gtw_path, rust_dump = sys.argv[1], sys.argv[2], sys.argv[3]
    names, cnt, out = load_comp(graph_path)
    vids, src, dst, x, e = build_feats(names, cnt, out)

    w = load_weights(gtw_path)
    model = GTNet(3, 4, 64, 2)
    model.load_state_dict({k: torch.tensor(v) for k, v in w.items()})
    model.eval()

    # 不用 dgl: 手工 MHA 消息传递（与 gt_train3.MHA 同语义, 逐 head 权重）
    import torch.nn.functional as F

    with torch.no_grad():
        h = model.lin_x(x)
        ef = model.lin_e(e)
        src_t = torch.tensor(src)
        dst_t = torch.tensor(dst)
        for lay in model.layers:
            n = h.shape[0]
            heads, od = 8, 8
            Q = lay.attn.Q(h).view(-1, heads, od)
            K = lay.attn.K(h).view(-1, heads, od)
            V = lay.attn.V(h).view(-1, heads, od)
            pe = lay.attn.proj_e(ef).view(-1, heads, od)
            score = K[src_t] * Q[dst_t] / math.sqrt(od) * pe  # (nE, 8, 8)
            e_a = score
            wgt = torch.exp(score.sum(-1, keepdim=True).clamp(-5, 5))  # (nE, 8, 1)
            msg = V[src_t] * wgt  # (nE, 8, 8)
            wV = torch.zeros(n, heads, od).index_add_(0, dst_t, msg)
            z = torch.zeros(n, heads, 1).index_add_(0, dst_t, wgt)
            h_a = (wV / (z + 1e-6)).view(-1, 64)
            e_a = e_a.view(-1, 64)
            h = lay.ln1_h(h + lay.O_h(h_a))
            ef = lay.ln1_e(ef + lay.O_e(e_a))
            h = lay.ln2_h(h + lay.f_h(h))
            ef = lay.ln2_e(ef + lay.f_e(ef))
        cat = torch.cat([h[src_t], h[dst_t], ef], 1)
        prob = torch.sigmoid(model.pred(cat).squeeze(-1)).numpy()

    rust = {}
    with open(rust_dump) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rust[int(f[0])] = (int(f[1]), int(f[2]), float(f[3]))
    assert len(rust) == len(src), f"边数不符: rust {len(rust)} vs python {len(src)}"
    max_diff = 0.0
    n_bad = 0
    for ei in range(len(src)):
        rv, rw, rp = rust[ei]
        assert (rv, rw) == (vids[src[ei]], vids[dst[ei]]), f"边 {ei} 端点不符"
        d = abs(rp - float(prob[ei]))
        max_diff = max(max_diff, d)
        if d > 1e-4:
            n_bad += 1
            if n_bad <= 5:
                print(f"  边 {ei} ({rv}->{rw}): rust {rp:.8e} vs torch {prob[ei]:.8e} 差 {d:.2e}")
    print(f"边数 {len(src)}, 最大逐边差 {max_diff:.3e}, 超容差(1e-4)边数 {n_bad}")
    sys.exit(1 if n_bad else 0)


if __name__ == "__main__":
    main()
