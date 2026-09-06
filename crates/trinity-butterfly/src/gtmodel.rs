//! GTNet 权重加载（GTW1 二进制）与纯 f32 前向——`--path-engine gtasm` 的边打分器。
//!
//! 参考实现（python 原型，逐项核对）:
//! - 架构: `export/gt_train3.py` 的 GTNet（lin_x/lin_e → N×GTLayer → pred）;
//! - 权重导出: `export/export_gt_weights.py`。
//!
//! GTW1 格式（全小端）: magic "GTW1" | u32 n_tensors | 每 tensor:
//! u32 name_len | name bytes | u32 ndim | ndim×u32 shape | prod(shape)×f32（C 序）。
//! torch Linear 的 weight 形状为 (out, in)，前向 y[o] = b[o] + Σ_i W[o][i]·x[i]。
//!
//! 前向数值语义（与 torch 逐项对应，容差目标 1e-4）:
//! - 每边每 head: `score[hd][dd] = K[src][hd][dd] * Q[dst][hd][dd] / sqrt(head_dim)
//!   * proj_e(e)[hd][dd]`; 边激活 e_a = score 原样（64 维）。
//! - 边权重**逐 head**: `w[hd] = exp(clamp(Σ_dd score[hd][dd], -5, 5))`
//!   （python 是 `score.sum(-1, keepdim=True)`——只对末维 head_dim 求和,
//!   每边得到 heads 个权重, 不是单标量）。
//! - 聚合: `wV[dst][hd][dd] = Σ_{e∈in(dst)} w_e[hd]·V[src][hd][dd]`,
//!   `z[dst][hd] = Σ w_e[hd]`; `h_a = wV / (z + 1e-6)`（无入边节点 → 0）。
//! - `h = LN1_h(h + O_h(h_a))`; `e = LN1_e(e + O_e(e_a))`;
//!   `h = LN2_h(h + FFN_h(h))`; `e = LN2_e(e + FFN_e(e))`;
//!   FFN = Linear(d,2d) → ReLU → Linear(2d,d); LayerNorm eps=1e-5（γ/β, 有偏方差）。
//! - pred: 每边 `cat[h_src, h_dst, e]` (3d) → Linear(3d,32) → ReLU → Linear(32,1)
//!   → sigmoid 得边概率。

use std::path::Path;

use rustc_hash::FxHashMap;

use crate::{BflyError, Result};

/// GTW1 文件魔数。
const GTW_MAGIC: &[u8; 4] = b"GTW1";
/// 注意力头数（python MHA 固定 heads=8）。
pub const GT_HEADS: usize = 8;
/// torch LayerNorm 默认 eps。
const LN_EPS: f32 = 1e-5;
/// 注意力权重指数前的 clamp 界（python clamp(-5, 5)）。
const ATTN_CLAMP: f32 = 5.0;
/// 归一化分母修正（python wV / (z + 1e-6)）。
const Z_EPS: f32 = 1e-6;

/// torch `nn.Linear` 的 f32 镜像（weight (out,in) C 序 + bias）。
#[derive(Debug, Clone, PartialEq)]
pub struct Linear {
    pub out_features: usize,
    pub in_features: usize,
    weight: Vec<f32>,
    bias: Vec<f32>,
}

impl Linear {
    /// y = x @ W^T + b（单行）。
    #[inline]
    pub fn forward_row(&self, x: &[f32], out: &mut [f32]) {
        debug_assert_eq!(x.len(), self.in_features);
        debug_assert_eq!(out.len(), self.out_features);
        for (o, yo) in out.iter_mut().enumerate() {
            let row = &self.weight[o * self.in_features..(o + 1) * self.in_features];
            let mut acc = self.bias[o];
            for (w, &xv) in row.iter().zip(x.iter()) {
                acc += w * xv;
            }
            *yo = acc;
        }
    }
}

/// torch `nn.LayerNorm`（有 γ/β, 有偏方差, eps=1e-5）。
#[derive(Debug, Clone, PartialEq)]
pub struct LayerNorm {
    gamma: Vec<f32>,
    beta: Vec<f32>,
}

impl LayerNorm {
    fn forward_row(&self, x: &[f32], out: &mut [f32]) {
        let n = x.len() as f32;
        let mean = x.iter().sum::<f32>() / n;
        let var = x.iter().map(|v| (v - mean) * (v - mean)).sum::<f32>() / n;
        let inv = 1.0 / (var + LN_EPS).sqrt();
        for (i, yo) in out.iter_mut().enumerate() {
            *yo = (x[i] - mean) * inv * self.gamma[i] + self.beta[i];
        }
    }
}

/// 一层 GTLayer（MHA + O_h/O_e + 双 LayerNorm + 双 FFN）。
#[derive(Debug, Clone)]
struct GtLayer {
    q: Linear,
    k: Linear,
    v: Linear,
    proj_e: Linear,
    o_h: Linear,
    o_e: Linear,
    ln1_h: LayerNorm,
    ln1_e: LayerNorm,
    ln2_h: LayerNorm,
    ln2_e: LayerNorm,
    f_h0: Linear,
    f_h2: Linear,
    f_e0: Linear,
    f_e2: Linear,
}

/// GTNet 全体权重（从 GTW1 装载; 层数由张量名 `layers.{i}.` 推导）。
#[derive(Debug, Clone)]
pub struct GtModel {
    /// 隐维 d（lin_x 输出维; 要求 d % GT_HEADS == 0）。
    pub hidden: usize,
    /// 节点特征维（lin_x 输入维, 原型 = 3）。
    pub xdim: usize,
    /// 边特征维（lin_e 输入维, 原型 = 4）。
    pub edim: usize,
    lin_x: Linear,
    lin_e: Linear,
    layers: Vec<GtLayer>,
    pred0: Linear,
    pred2: Linear,
}

/// sigmoid（torch 同式: 1/(1+exp(-x))）。
#[inline]
fn sigmoid(x: f32) -> f32 {
    1.0 / (1.0 + (-x).exp())
}

/// 对 n 行输入逐行过 Linear → n×out_features 平坦缓冲。
fn matmul_rows(l: &Linear, x: &[f32], n: usize) -> Vec<f32> {
    let mut out = vec![0.0; n * l.out_features];
    for i in 0..n {
        l.forward_row(
            &x[i * l.in_features..(i + 1) * l.in_features],
            &mut out[i * l.out_features..(i + 1) * l.out_features],
        );
    }
    out
}

impl GtLayer {
    /// 单层前向（就地更新 h (n×d) 与 e (ne×d)）。
    fn forward(
        &self,
        h: &mut [f32],
        e: &mut [f32],
        src: &[usize],
        dst: &[usize],
        n: usize,
        d: usize,
    ) {
        let hd = d / GT_HEADS;
        let ne = src.len();
        let q = matmul_rows(&self.q, h, n);
        let k = matmul_rows(&self.k, h, n);
        let v = matmul_rows(&self.v, h, n);
        let pe = matmul_rows(&self.proj_e, e, ne);
        let sqrt_hd = (hd as f32).sqrt();

        // 每边激活 e_a 与逐 head 权重 w
        let mut e_a = vec![0f32; ne * d];
        let mut w = vec![0f32; ne * GT_HEADS];
        for (ei, (&s, &dd)) in src.iter().zip(dst.iter()).enumerate() {
            let kb = &k[s * d..(s + 1) * d];
            let qb = &q[dd * d..(dd + 1) * d];
            let pb = &pe[ei * d..(ei + 1) * d];
            let ea = &mut e_a[ei * d..(ei + 1) * d];
            for hd_i in 0..GT_HEADS {
                let mut sum = 0f32;
                for j in 0..hd {
                    let idx = hd_i * hd + j;
                    // torch: (K*Q)/sqrt(od) 再 *proj_e——保持同一运算次序
                    let sc = (kb[idx] * qb[idx]) / sqrt_hd * pb[idx];
                    ea[idx] = sc;
                    sum += sc;
                }
                w[ei * GT_HEADS + hd_i] = sum.clamp(-ATTN_CLAMP, ATTN_CLAMP).exp();
            }
        }

        // 聚合到 dst: wV/z（逐 head）
        let mut h_a = vec![0f32; n * d];
        let mut z = vec![0f32; n * GT_HEADS];
        for (ei, (&s, &dd)) in src.iter().zip(dst.iter()).enumerate() {
            for hd_i in 0..GT_HEADS {
                let we = w[ei * GT_HEADS + hd_i];
                z[dd * GT_HEADS + hd_i] += we;
                for j in 0..hd {
                    let idx = hd_i * hd + j;
                    h_a[dd * d + idx] += we * v[s * d + idx];
                }
            }
        }
        for i in 0..n {
            for hd_i in 0..GT_HEADS {
                let denom = z[i * GT_HEADS + hd_i] + Z_EPS;
                for j in 0..hd {
                    h_a[i * d + hd_i * hd + j] /= denom;
                }
            }
        }

        // 残差 + LN1（h 侧 / e 侧）
        add_linear_ln(h, &h_a, &self.o_h, &self.ln1_h, n, d);
        add_linear_ln(e, &e_a, &self.o_e, &self.ln1_e, ne, d);
        // FFN 残差 + LN2
        ffn_ln(h, &self.f_h0, &self.f_h2, &self.ln2_h, n, d);
        ffn_ln(e, &self.f_e0, &self.f_e2, &self.ln2_e, ne, d);
    }
}

/// base[i] = LN(base[i] + Linear(inp[i]))（逐行）。
fn add_linear_ln(base: &mut [f32], inp: &[f32], l: &Linear, ln: &LayerNorm, n: usize, d: usize) {
    let mut tmp = vec![0f32; d];
    let mut normed = vec![0f32; d];
    for i in 0..n {
        l.forward_row(&inp[i * d..(i + 1) * d], &mut tmp);
        for j in 0..d {
            tmp[j] += base[i * d + j];
        }
        ln.forward_row(&tmp, &mut normed);
        base[i * d..(i + 1) * d].copy_from_slice(&normed);
    }
}

/// base[i] = LN(base[i] + FFN(base[i])), FFN = Linear→ReLU→Linear（逐行）。
fn ffn_ln(base: &mut [f32], f0: &Linear, f2: &Linear, ln: &LayerNorm, n: usize, d: usize) {
    let mut a = vec![0f32; f0.out_features];
    let mut b = vec![0f32; d];
    let mut normed = vec![0f32; d];
    for i in 0..n {
        f0.forward_row(&base[i * d..(i + 1) * d], &mut a);
        for v in a.iter_mut() {
            if *v < 0.0 {
                *v = 0.0;
            }
        }
        f2.forward_row(&a, &mut b);
        for j in 0..d {
            b[j] += base[i * d + j];
        }
        ln.forward_row(&b, &mut normed);
        base[i * d..(i + 1) * d].copy_from_slice(&normed);
    }
}

impl GtModel {
    /// 层数（GTLayer 个数）。
    pub fn n_layers(&self) -> usize {
        self.layers.len()
    }

    /// 从 GTW1 文件加载。
    pub fn load(path: &Path) -> Result<GtModel> {
        let data = std::fs::read(path)
            .map_err(|e| BflyError::GtWeights(format!("cannot read {}: {e}", path.display())))?;
        Self::from_bytes(&data)
    }

    /// 从 GTW1 字节加载（解析 + 形状校验; 多余/缺失张量均报错）。
    pub fn from_bytes(data: &[u8]) -> Result<GtModel> {
        let tensors = parse_gtw(data)?;
        let lin_x = take_linear(&tensors, "lin_x")?;
        let lin_e = take_linear(&tensors, "lin_e")?;
        let hidden = lin_x.out_features;
        if hidden == 0 || hidden % GT_HEADS != 0 {
            return Err(BflyError::GtWeights(format!(
                "hidden ({hidden}) 必须为正且 % {GT_HEADS} == 0"
            )));
        }
        if lin_e.out_features != hidden {
            return Err(BflyError::GtWeights(format!(
                "lin_e 输出维 {} != hidden {hidden}",
                lin_e.out_features
            )));
        }
        // 层数 = `layers.{i}.` 前缀的最大 i + 1（要求 0..L 连续）
        let mut layer_idx: Vec<usize> = tensors
            .keys()
            .filter_map(|name| {
                name.strip_prefix("layers.")
                    .and_then(|rest| rest.split('.').next())
                    .and_then(|i| i.parse::<usize>().ok())
            })
            .collect();
        layer_idx.sort_unstable();
        layer_idx.dedup();
        let n_layers = layer_idx.last().map(|i| i + 1).unwrap_or(0);
        if layer_idx != (0..n_layers).collect::<Vec<_>>() {
            return Err(BflyError::GtWeights(format!(
                "layers.* 下标不连续: {layer_idx:?}"
            )));
        }
        let mut layers = Vec::with_capacity(n_layers);
        for i in 0..n_layers {
            let p = format!("layers.{i}");
            let layer = GtLayer {
                q: take_linear(&tensors, &format!("{p}.attn.Q"))?,
                k: take_linear(&tensors, &format!("{p}.attn.K"))?,
                v: take_linear(&tensors, &format!("{p}.attn.V"))?,
                proj_e: take_linear(&tensors, &format!("{p}.attn.proj_e"))?,
                o_h: take_linear(&tensors, &format!("{p}.O_h"))?,
                o_e: take_linear(&tensors, &format!("{p}.O_e"))?,
                ln1_h: take_ln(&tensors, &format!("{p}.ln1_h"))?,
                ln1_e: take_ln(&tensors, &format!("{p}.ln1_e"))?,
                ln2_h: take_ln(&tensors, &format!("{p}.ln2_h"))?,
                ln2_e: take_ln(&tensors, &format!("{p}.ln2_e"))?,
                f_h0: take_linear(&tensors, &format!("{p}.f_h.0"))?,
                f_h2: take_linear(&tensors, &format!("{p}.f_h.2"))?,
                f_e0: take_linear(&tensors, &format!("{p}.f_e.0"))?,
                f_e2: take_linear(&tensors, &format!("{p}.f_e.2"))?,
            };
            for (tag, l) in [
                ("attn.Q", &layer.q),
                ("attn.K", &layer.k),
                ("attn.V", &layer.v),
                ("attn.proj_e", &layer.proj_e),
                ("O_h", &layer.o_h),
                ("O_e", &layer.o_e),
            ] {
                if l.out_features != hidden || l.in_features != hidden {
                    return Err(BflyError::GtWeights(format!(
                        "{p}.{tag} 形状 ({},{}) 与 hidden {hidden} 不符",
                        l.out_features, l.in_features
                    )));
                }
            }
            layers.push(layer);
        }
        let pred0 = take_linear(&tensors, "pred.0")?;
        let pred2 = take_linear(&tensors, "pred.2")?;
        if pred0.in_features != 3 * hidden {
            return Err(BflyError::GtWeights(format!(
                "pred.0 输入维 {} != 3*hidden {}",
                pred0.in_features,
                3 * hidden
            )));
        }
        if pred2.in_features != pred0.out_features || pred2.out_features != 1 {
            return Err(BflyError::GtWeights(format!(
                "pred.2 形状 ({},{}) 与 pred.0 输出维 {} 不符",
                pred2.out_features, pred2.in_features, pred0.out_features
            )));
        }
        // 已知张量全部消费完才算格式吻合
        let mut known: Vec<String> = vec![
            "lin_x".into(),
            "lin_e".into(),
            "pred.0".into(),
            "pred.2".into(),
        ];
        for i in 0..n_layers {
            let p = format!("layers.{i}");
            for s in [
                "attn.Q", "attn.K", "attn.V", "attn.proj_e", "O_h", "O_e", "ln1_h", "ln1_e",
                "ln2_h", "ln2_e", "f_h.0", "f_h.2", "f_e.0", "f_e.2",
            ] {
                known.push(format!("{p}.{s}"));
            }
        }
        for name in tensors.keys() {
            let base = name.strip_suffix(".weight").or(name.strip_suffix(".bias"));
            match base {
                Some(b) if known.iter().any(|k| k == b) => {}
                _ => {
                    return Err(BflyError::GtWeights(format!(
                        "未知张量名 {name}（GTW 格式漂移?）"
                    )))
                }
            }
        }
        Ok(GtModel {
            hidden,
            xdim: lin_x.in_features,
            edim: lin_e.in_features,
            lin_x,
            lin_e,
            layers,
            pred0,
            pred2,
        })
    }

    /// 全边前向: x (n×xdim) / e (ne×edim) / 边表 (src,dst) → 每边概率（f32）。
    ///
    /// panic-free 前置由调用方保证: src/dst 等长且节点下标 < n。
    pub fn forward(&self, x: &[f32], e: &[f32], src: &[usize], dst: &[usize]) -> Vec<f32> {
        let d = self.hidden;
        let n = x.len() / self.xdim;
        let ne = src.len();
        debug_assert_eq!(x.len(), n * self.xdim);
        debug_assert_eq!(e.len(), ne * self.edim);
        debug_assert_eq!(dst.len(), ne);

        let mut h = matmul_rows(&self.lin_x, x, n);
        let mut ef = matmul_rows(&self.lin_e, e, ne);
        for layer in &self.layers {
            layer.forward(&mut h, &mut ef, src, dst, n, d);
        }

        let mut buf = vec![0f32; 3 * d];
        let mut h1 = vec![0f32; self.pred0.out_features];
        let mut logit = [0f32; 1];
        let mut probs = Vec::with_capacity(ne);
        for (ei, (&s, &dd)) in src.iter().zip(dst.iter()).enumerate() {
            buf[..d].copy_from_slice(&h[s * d..(s + 1) * d]);
            buf[d..2 * d].copy_from_slice(&h[dd * d..(dd + 1) * d]);
            buf[2 * d..].copy_from_slice(&ef[ei * d..(ei + 1) * d]);
            self.pred0.forward_row(&buf, &mut h1);
            for v in h1.iter_mut() {
                if *v < 0.0 {
                    *v = 0.0;
                }
            }
            self.pred2.forward_row(&h1, &mut logit);
            probs.push(sigmoid(logit[0]));
        }
        probs
    }
}

// ---------------------------------------------------------------------------
// GTW1 解析
// ---------------------------------------------------------------------------

/// 测试共用: 0 层全零 hidden=8 模型（全边概率恒 0.5）。供 gtdecode 测试复用。
#[cfg(test)]
pub(crate) fn zero_model_for_test() -> GtModel {
    let tensors: Vec<(&str, Vec<usize>, Vec<f32>)> = vec![
        ("lin_x.weight", vec![8, 3], vec![0.0; 24]),
        ("lin_x.bias", vec![8], vec![0.0; 8]),
        ("lin_e.weight", vec![8, 4], vec![0.0; 32]),
        ("lin_e.bias", vec![8], vec![0.0; 8]),
        ("pred.0.weight", vec![32, 24], vec![0.0; 768]),
        ("pred.0.bias", vec![32], vec![0.0; 32]),
        ("pred.2.weight", vec![1, 32], vec![0.0; 32]),
        ("pred.2.bias", vec![1], vec![0.0; 1]),
    ];
    let mut out = Vec::new();
    out.extend_from_slice(b"GTW1");
    out.extend_from_slice(&(tensors.len() as u32).to_le_bytes());
    for (name, shape, vals) in &tensors {
        out.extend_from_slice(&(name.len() as u32).to_le_bytes());
        out.extend_from_slice(name.as_bytes());
        out.extend_from_slice(&(shape.len() as u32).to_le_bytes());
        for &d in shape {
            out.extend_from_slice(&(d as u32).to_le_bytes());
        }
        for &v in vals {
            out.extend_from_slice(&v.to_le_bytes());
        }
    }
    GtModel::from_bytes(&out).unwrap()
}

/// name → (shape, data f32 C 序)。
type TensorMap = FxHashMap<String, (Vec<usize>, Vec<f32>)>;

struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        if self.pos + n > self.data.len() {
            return Err(BflyError::GtWeights(format!(
                "GTW 截断: 需要 {n} 字节 @{}, 余 {}",
                self.pos,
                self.data.len() - self.pos
            )));
        }
        let s = &self.data[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }

    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
}

fn parse_gtw(data: &[u8]) -> Result<TensorMap> {
    let mut r = Reader { data, pos: 0 };
    if r.take(4)? != GTW_MAGIC {
        return Err(BflyError::GtWeights("magic 不是 GTW1".to_string()));
    }
    let n_tensors = r.u32()? as usize;
    let mut out: TensorMap = FxHashMap::default();
    for _ in 0..n_tensors {
        let name_len = r.u32()? as usize;
        let name = String::from_utf8(r.take(name_len)?.to_vec())
            .map_err(|e| BflyError::GtWeights(format!("张量名非 UTF-8: {e}")))?;
        let ndim = r.u32()? as usize;
        if ndim > 4 {
            return Err(BflyError::GtWeights(format!("{name}: ndim {ndim} > 4")));
        }
        let mut shape = Vec::with_capacity(ndim);
        let mut count = 1usize;
        for _ in 0..ndim {
            let d = r.u32()? as usize;
            count = count
                .checked_mul(d)
                .ok_or_else(|| BflyError::GtWeights(format!("{name}: 形状溢出")))?;
            shape.push(d);
        }
        let bytes = r.take(4 * count)?;
        let mut vals = Vec::with_capacity(count);
        for c in bytes.chunks_exact(4) {
            vals.push(f32::from_le_bytes(c.try_into().unwrap()));
        }
        if out.insert(name.clone(), (shape, vals)).is_some() {
            return Err(BflyError::GtWeights(format!("重复张量名 {name}")));
        }
    }
    if r.pos != data.len() {
        return Err(BflyError::GtWeights(format!(
            "GTW 尾部多余 {} 字节",
            data.len() - r.pos
        )));
    }
    Ok(out)
}

/// 取一个 Linear（`<name>.weight` (out,in) + `<name>.bias` (out,)）。
fn take_linear(tensors: &TensorMap, name: &str) -> Result<Linear> {
    let (shape, w) = tensors
        .get(&format!("{name}.weight"))
        .ok_or_else(|| BflyError::GtWeights(format!("缺少 {name}.weight")))?;
    if shape.len() != 2 {
        return Err(BflyError::GtWeights(format!(
            "{name}.weight 应为 2 维, 实际 {shape:?}"
        )));
    }
    let (bshape, b) = tensors
        .get(&format!("{name}.bias"))
        .ok_or_else(|| BflyError::GtWeights(format!("缺少 {name}.bias")))?;
    if bshape.len() != 1 || bshape[0] != shape[0] {
        return Err(BflyError::GtWeights(format!(
            "{name}.bias 形状 {bshape:?} 与 weight {shape:?} 不符"
        )));
    }
    Ok(Linear {
        out_features: shape[0],
        in_features: shape[1],
        weight: w.clone(),
        bias: b.clone(),
    })
}

/// 取一个 LayerNorm（`<name>.weight` / `<name>.bias` 同长 1 维）。
fn take_ln(tensors: &TensorMap, name: &str) -> Result<LayerNorm> {
    let (gshape, g) = tensors
        .get(&format!("{name}.weight"))
        .ok_or_else(|| BflyError::GtWeights(format!("缺少 {name}.weight")))?;
    let (bshape, b) = tensors
        .get(&format!("{name}.bias"))
        .ok_or_else(|| BflyError::GtWeights(format!("缺少 {name}.bias")))?;
    if gshape.len() != 1 || bshape != gshape {
        return Err(BflyError::GtWeights(format!(
            "{name} LayerNorm 形状不符: {gshape:?} / {bshape:?}"
        )));
    }
    Ok(LayerNorm {
        gamma: g.clone(),
        beta: b.clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 测试内 GTW1 构造器（与 export_gt_weights.py 同格式）。
    fn build_gtw(tensors: &[(&str, Vec<usize>, Vec<f32>)]) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(b"GTW1");
        out.extend_from_slice(&(tensors.len() as u32).to_le_bytes());
        for (name, shape, vals) in tensors {
            out.extend_from_slice(&(name.len() as u32).to_le_bytes());
            out.extend_from_slice(name.as_bytes());
            out.extend_from_slice(&(shape.len() as u32).to_le_bytes());
            for &d in shape {
                out.extend_from_slice(&(d as u32).to_le_bytes());
            }
            for &v in vals {
                out.extend_from_slice(&v.to_le_bytes());
            }
        }
        out
    }

    /// 0 层最小模型（hidden=8）: lin_x/lin_e/pred 全零 → 所有边概率恒 0.5。
    fn zero_model_bytes() -> Vec<u8> {
        build_gtw(&[
            ("lin_x.weight", vec![8, 3], vec![0.0; 24]),
            ("lin_x.bias", vec![8], vec![0.0; 8]),
            ("lin_e.weight", vec![8, 4], vec![0.0; 32]),
            ("lin_e.bias", vec![8], vec![0.0; 8]),
            ("pred.0.weight", vec![32, 24], vec![0.0; 768]),
            ("pred.0.bias", vec![32], vec![0.0; 32]),
            ("pred.2.weight", vec![1, 32], vec![0.0; 32]),
            ("pred.2.bias", vec![1], vec![0.0; 1]),
        ])
    }

    #[test]
    fn roundtrip_zero_model_all_half() {
        let bytes = zero_model_bytes();
        let m = GtModel::from_bytes(&bytes).unwrap();
        assert_eq!(m.hidden, 8);
        assert_eq!(m.xdim, 3);
        assert_eq!(m.edim, 4);
        assert_eq!(m.layers.len(), 0);
        // 2 节点 1 边
        let x = vec![1.0f32; 2 * 3];
        let e = vec![1.0f32; 4];
        let p = m.forward(&x, &e, &[0], &[1]);
        assert_eq!(p, vec![0.5]);
    }

    #[test]
    fn pred_bias_flows_to_prob() {
        // pred.2.bias = 2 → logit 2 → sigmoid(2)
        let bytes = build_gtw(&[
            ("lin_x.weight", vec![8, 3], vec![0.0; 24]),
            ("lin_x.bias", vec![8], vec![0.0; 8]),
            ("lin_e.weight", vec![8, 4], vec![0.0; 32]),
            ("lin_e.bias", vec![8], vec![0.0; 8]),
            ("pred.0.weight", vec![32, 24], vec![0.0; 768]),
            ("pred.0.bias", vec![32], vec![0.0; 32]),
            ("pred.2.weight", vec![1, 32], vec![0.0; 32]),
            ("pred.2.bias", vec![1], vec![2.0]),
        ]);
        let m = GtModel::from_bytes(&bytes).unwrap();
        let p = m.forward(&[0.0; 6], &[0.0; 4], &[0], &[1]);
        let expect = sigmoid(2.0);
        assert!((p[0] - expect).abs() < 1e-7, "{} vs {expect}", p[0]);
    }

    #[test]
    fn zero_layer_wiring_hand_computed() {
        // 0 层: h = lin_x(x), e = lin_e(e), pred 直通——全程可手算。
        // hidden=8。lin_x.weight[0][0]=1 → h[i][0] = x[i][0]; 其余全零。
        // lin_e.weight[0][0]=2 → e[0][0] = 2*e_in[0]。
        // pred.0.weight[0][0]=1（取 h_src[0]）, pred.0.weight[1][8]=1（取 h_dst[0]）,
        // pred.0.weight[2][16]=1（取 e[0]）; pred.2.weight = [1,1,1,0...], bias 0.25。
        let mut lx_w = vec![0.0f32; 24];
        lx_w[0] = 1.0;
        let mut le_w = vec![0.0f32; 32];
        le_w[0] = 2.0;
        let mut p0_w = vec![0.0f32; 32 * 24];
        p0_w[0] = 1.0; // row0: h_src[0]
        p0_w[24 + 8] = 1.0; // row1: h_dst[0]
        p0_w[2 * 24 + 16] = 1.0; // row2: e[0]
        let mut p2_w = vec![0.0f32; 32];
        p2_w[0] = 1.0;
        p2_w[1] = 1.0;
        p2_w[2] = 1.0;
        let bytes = build_gtw(&[
            ("lin_x.weight", vec![8, 3], lx_w),
            ("lin_x.bias", vec![8], vec![0.0; 8]),
            ("lin_e.weight", vec![8, 4], le_w),
            ("lin_e.bias", vec![8], vec![0.0; 8]),
            ("pred.0.weight", vec![32, 24], p0_w),
            ("pred.0.bias", vec![32], vec![0.0; 32]),
            ("pred.2.weight", vec![1, 32], p2_w),
            ("pred.2.bias", vec![1], vec![0.25]),
        ]);
        let m = GtModel::from_bytes(&bytes).unwrap();
        // 节点 0: x=[0.5,0,0] → h0[0]=0.5; 节点 1: x=[1.5,0,0] → h1[0]=1.5
        // 边 0→1: e=[0.25,...] → ef[0]=0.5
        let x = vec![0.5, 0.0, 0.0, 1.5, 0.0, 0.0];
        let e = vec![0.25, 0.0, 0.0, 0.0];
        let p = m.forward(&x, &e, &[0], &[1]);
        // logit = relu(0.5) + relu(1.5) + relu(0.5) + 0.25 = 2.75
        let expect = sigmoid(2.75);
        assert!((p[0] - expect).abs() < 1e-6, "{} vs {expect}", p[0]);
    }

    #[test]
    fn layernorm_hand_computed() {
        // LN([1,0,...,0]) (dim 8): mean=0.125, var=0.109375
        let ln = LayerNorm {
            gamma: vec![1.0; 8],
            beta: vec![0.0; 8],
        };
        let x = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut out = [0.0f32; 8];
        ln.forward_row(&x, &mut out);
        let mean = 0.125f32;
        let var = 0.109375f32;
        let inv = 1.0 / (var + LN_EPS).sqrt();
        assert!((out[0] - (1.0 - mean) * inv).abs() < 1e-6);
        assert!((out[1] - (0.0 - mean) * inv).abs() < 1e-6);
    }

    #[test]
    fn loader_rejects_bad_magic_and_missing() {
        assert!(GtModel::from_bytes(b"XXXX").is_err());
        // 只有 lin_x → 缺 lin_e
        let bytes = build_gtw(&[
            ("lin_x.weight", vec![8, 3], vec![0.0; 24]),
            ("lin_x.bias", vec![8], vec![0.0; 8]),
        ]);
        assert!(GtModel::from_bytes(&bytes).is_err());
        // 未知张量 → 报错
        let mut t = vec![
            ("lin_x.weight", vec![8, 3], vec![0.0; 24]),
            ("lin_x.bias", vec![8], vec![0.0; 8]),
            ("lin_e.weight", vec![8, 4], vec![0.0; 32]),
            ("lin_e.bias", vec![8], vec![0.0; 8]),
            ("pred.0.weight", vec![32, 24], vec![0.0; 768]),
            ("pred.0.bias", vec![32], vec![0.0; 32]),
            ("pred.2.weight", vec![1, 32], vec![0.0; 32]),
            ("pred.2.bias", vec![1], vec![0.0]),
        ];
        t.push(("mystery.weight", vec![1], vec![0.0]));
        assert!(GtModel::from_bytes(&build_gtw(&t)).is_err());
    }

    #[test]
    fn single_layer_runs_and_is_deterministic() {
        // 单层 hidden=8（head_dim=1）: 只验证可跑/确定性/输出域（手算全链过长,
        // 数值对拍由 gtasm_score_dump + python 脚本完成）。
        let mut t: Vec<(&str, Vec<usize>, Vec<f32>)> = vec![
            ("lin_x.weight", vec![8, 3], vec![0.0; 24]),
            ("lin_x.bias", vec![8], vec![0.0; 8]),
            ("lin_e.weight", vec![8, 4], vec![0.0; 32]),
            ("lin_e.bias", vec![8], vec![0.0; 8]),
        ];
        for name in ["attn.Q", "attn.K", "attn.V", "attn.proj_e", "O_h", "O_e"] {
            t.push((Box::leak(format!("layers.0.{name}.weight").into_boxed_str()), vec![8, 8], vec![0.0; 64]));
            t.push((Box::leak(format!("layers.0.{name}.bias").into_boxed_str()), vec![8], vec![0.0; 8]));
        }
        for name in ["ln1_h", "ln1_e", "ln2_h", "ln2_e"] {
            t.push((Box::leak(format!("layers.0.{name}.weight").into_boxed_str()), vec![8], vec![1.0; 8]));
            t.push((Box::leak(format!("layers.0.{name}.bias").into_boxed_str()), vec![8], vec![0.0; 8]));
        }
        for name in ["f_h.0", "f_e.0"] {
            t.push((Box::leak(format!("layers.0.{name}.weight").into_boxed_str()), vec![16, 8], vec![0.0; 128]));
            t.push((Box::leak(format!("layers.0.{name}.bias").into_boxed_str()), vec![16], vec![0.0; 16]));
        }
        for name in ["f_h.2", "f_e.2"] {
            t.push((Box::leak(format!("layers.0.{name}.weight").into_boxed_str()), vec![8, 16], vec![0.0; 128]));
            t.push((Box::leak(format!("layers.0.{name}.bias").into_boxed_str()), vec![8], vec![0.0; 8]));
        }
        t.push(("pred.0.weight", vec![32, 24], vec![0.0; 768]));
        t.push(("pred.0.bias", vec![32], vec![0.0; 32]));
        t.push(("pred.2.weight", vec![1, 32], vec![0.0; 32]));
        t.push(("pred.2.bias", vec![1], vec![0.0]));
        let m = GtModel::from_bytes(&build_gtw(&t)).unwrap();
        assert_eq!(m.layers.len(), 1);
        let x = vec![0.1f32; 3 * 3];
        let e = vec![0.2f32; 2 * 4];
        let src = [0usize, 1];
        let dst = [1usize, 2];
        let p1 = m.forward(&x, &e, &src, &dst);
        let p2 = m.forward(&x, &e, &src, &dst);
        assert_eq!(p1, p2);
        assert!(p1.iter().all(|&p| (0.0..=1.0).contains(&p)));
        // 零权重 + 零偏置: 全链 0（LN(0)=0）→ logit 0 → 0.5
        assert!(p1.iter().all(|&p| (p - 0.5).abs() < 1e-7));
    }
}
