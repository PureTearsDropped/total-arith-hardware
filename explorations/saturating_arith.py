#!/usr/bin/env python3
"""**本物の飽和演算** — 除算器だけでなく、乗算器も飽和させる。

overflow_convention.py の実測でバグが本質を見せた:
    8.99 × 5 × 5 × **MAX** = **オーバーフローして Inf**
  ⟹ **MAX は除算器だけでは足りない。次の掛け算で即座に範囲を出て、NaN が戻る。**
  ⟹ **整数の飽和演算が「全演算で飽和する」理由が、これ。**

本物の飽和（整数 SIMD の paddsb / ARM の VQADD と同じ思想）:
    sat(a op b) = **必ず [−MAX, +MAX] に留める。Inf も NaN も、決して生まれない**

利用者の規約:
    a/0 = **+MAX** (a>0) ／ **−MAX** (a<0) ／ **0/0 = 0**（不定 ⟹ 中立）
    そして **全ての演算が飽和する**

4つ比べる: Inf（IEEE） ／ **飽和 ±MAX** ／ 0 ／ **飽和 ±MAX ＋ 頑健な損失**
"""
import warnings
import numpy as np
from scipy.optimize import least_squares
warnings.filterwarnings("ignore")
from zero_total import make

MAX = np.finfo(np.float64).max

def sat(v):
    """**Inf も NaN も、決して外に出さない。** 必ず [−MAX, +MAX]。"""
    v = np.nan_to_num(np.asarray(v, float), nan=0.0, posinf=MAX, neginf=-MAX)
    return np.clip(v, -MAX, MAX)

def sat_mul(a, b):
    with np.errstate(all="ignore"):
        return sat(np.asarray(a, float)*np.asarray(b, float))

def sat_pow(x, w, mode):
    """x^w。mode で x=0, w<0 の扱いを変える。**飽和版は結果も飽和させる。**"""
    x = np.asarray(x, float)
    if w == 0.0: return np.ones_like(x)
    with np.errstate(all="ignore"):
        v = np.power(np.abs(x), w)
    if w > 0:
        out = np.where(np.abs(x) == 0.0, 0.0, v)
    else:
        fill = {"inf": np.inf, "sat": MAX, "zero": 0.0}[mode]
        out = np.where(np.abs(x) == 0.0, fill, v)
    return sat(out) if mode == "sat" else out

def basis(W, X, mode):
    W = np.asarray(W).reshape(-1, X.shape[1])
    cols = []
    for k in range(len(W)):
        c = np.ones(len(X))
        for j in range(X.shape[1]):
            c = sat_mul(c, sat_pow(X[:, j], W[k, j], mode)) if mode == "sat" \
                else c*sat_pow(X[:, j], W[k, j], mode)
        cols.append(c)
    return np.c_[tuple(cols) + (np.ones(len(X)),)]

print("=" * 88)
print("① 飽和なら、Inf も NaN も生まれないか")
print("=" * 88)
print(f"  {'式':<28}{'素の float64':>18}{'**飽和**':>16}")
for tag, f_raw, f_sat in (("MAX × 25", lambda: MAX*25.0, lambda: sat_mul(MAX, 25.0)),
                          ("MAX − MAX", lambda: MAX-MAX, lambda: sat(MAX-MAX)),
                          ("MAX / MAX", lambda: MAX/MAX, lambda: sat(MAX/MAX)),
                          ("0 × MAX", lambda: 0.0*MAX, lambda: sat_mul(0.0, MAX)),
                          ("(1/0) × 25 【Inf】", lambda: np.inf*25.0, lambda: sat_mul(MAX, 25.0))):
    with np.errstate(all="ignore"):
        a = float(np.atleast_1d(f_raw())[0]); b = float(np.atleast_1d(f_sat())[0])
    sa = "**Inf**" if np.isinf(a) else ("**NaN**" if np.isnan(a) else f"{a:.4g}")
    sb = "**Inf**" if np.isinf(b) else ("**NaN**" if np.isnan(b) else f"{b:.4g}")
    print(f"  {tag:<28}{sa:>18}{sb:>16}")

def fit(X, y, mode, loss="linear", seeds=25):
    F = X.shape[1]; wt = 1.0/np.maximum(np.abs(y), 1e-300)
    def resid(w):
        Phi = basis(w, X, mode)
        if not np.all(np.isfinite(Phi)): return np.full(len(y), 1e12)
        c = np.linalg.lstsq(Phi*wt[:,None], y*wt, rcond=None)[0]
        r = (Phi @ c - y)*wt
        return sat(r) if mode == "sat" else np.where(np.isfinite(r), r, 1e12)
    best = None
    for s in range(seeds):
        w0 = np.random.default_rng(s+7).normal(0, 1, F)
        try:
            o = least_squares(resid, w0, loss=loss, f_scale=1.0, max_nfev=5000,
                              bounds=(-6*np.ones(F), 6*np.ones(F)))
            v = float(np.sum(o.fun**2))
            if best is None or v < best[0]: best = (v, o.x)
        except Exception: pass
    return best[1] if best else np.full(F, np.nan)

print()
print("=" * 88)
print("② 当てはめ — 分解能で r=0 と記録された点がある（真 w(r) = −2）")
print("=" * 88)
print(f"  {'r=0 の割合':>11}{'Inf（IEEE）':>13}{'**飽和 MAX**':>14}{'**0**':>10}{'飽和＋Cauchy':>15}")
for frac in (0.0, 0.01, 0.05, 0.15):
    X, Xo, yv, nb = make(2000, frac, seed=1)
    r = []
    for mode, loss in (("inf","linear"), ("sat","linear"), ("zero","linear"), ("sat","cauchy")):
        w = fit(Xo, yv, mode, loss)
        r.append(w[2] if np.isfinite(w[2]) else np.nan)
    print(f"  {frac:>10.0%}{r[0]:>13.4f}{r[1]:>14.4f}{r[2]:>10.4f}{r[3]:>15.4f}")
print(f"  {'真値':>11}{-2.0:>13.4f}{-2.0:>14.4f}{-2.0:>10.4f}{-2.0:>15.4f}")

print()
print("=" * 88)
print("③ 残差の大きさ — r=0 の点 1つ（真の y ~ 1e6）")
print("=" * 88)
print(f"  {'規約':<16}{'模型の値':>16}{'相対残差':>16}{'有界か':>10}")
for mode in ("inf", "sat", "zero"):
    v = sat_pow(np.array([0.0]), -2.0, mode)[0]
    if mode == "sat": v = float(sat_mul(sat_mul(8.99*5.0, 5.0), v)[0] if np.ndim(v)==0 else v)
    else: v = 8.99*25.0*float(v)
    rr = (v - 1e6)/1e6
    print(f"  {mode:<16}{v:>16.4g}{rr:>16.4g}{('**○**' if abs(rr) < 1e3 else '×'):>10}")
