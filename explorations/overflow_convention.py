#!/usr/bin/env python3
"""NaN でなく **正/負のオーバーフロー値（±MAX）** にしたらどうか — 利用者の案。

IEEE 754:  0/0 = **NaN**、 1/0 = **+Inf**
  Inf − Inf = **NaN**、 Inf/Inf = **NaN**、 0×Inf = **NaN**   ← **Inf は NaN を製造する**

±MAX（飽和除算）:
  MAX − MAX = **0** ✓、 MAX/MAX = **1** ✓、 0×MAX = **0** ✓   ← **NaN が どこからも 生まれない**
  そして **符号が残る**: +MAX は「正に発散」、−MAX は「負に発散」。0 は何も言わない

  a/0 = **+MAX** (a>0) ／ **−MAX** (a<0) ／ **0/0 = 0**（不定 ⟹ 中立）

だが当てはめでは: r^{−2} が r=0 で +MAX ≈ 1.8e308 → y~1e6 に対し **相対残差 1.8e302。有界でない**。
zero_total.py の実測では、0 の規約は残差が **1 で頭打ち**になるから勝っていた。

4つ比べる: Inf（IEEE） / **±MAX** / 0 / **±MAX ＋ 頑健な損失**
"""
import warnings
import numpy as np
from scipy.optimize import least_squares
warnings.filterwarnings("ignore")
from zero_total import make, R_RES

MAX = np.finfo(np.float64).max

print("=" * 84)
print("① 素の算術 — NaN は生まれるか")
print("=" * 84)
print(f"  {'式':<20}{'Inf 規約':>16}{'**±MAX 規約**':>18}")
for expr, a, b, op in (("x − x", np.inf, np.inf, "-"), ("x / x", np.inf, np.inf, "/"),
                       ("0 × x", 0.0, np.inf, "*")):
    with np.errstate(all="ignore"):
        vi = eval(f"a {op} b")
        am, bm = (MAX if np.isinf(a) else a), (MAX if np.isinf(b) else b)
        vm = eval(f"am {op} bm")
    si = "**NaN**" if np.isnan(vi) else f"{vi:.4g}"
    sm = "**NaN**" if np.isnan(vm) else f"{vm:.4g}"
    print(f"  {expr:<20}{si:>16}{sm:>18}")
print(f"\n  ⟹ **Inf は NaN を作る。MAX は作らない。**")

def pu(x, w, mode):
    x = np.asarray(x, float)
    if w == 0.0: return np.ones_like(x)
    with np.errstate(all="ignore"):
        v = np.power(np.abs(x), w)
    if w > 0: return np.where(np.abs(x) == 0.0, 0.0, v)
    fill = {"inf": np.inf, "max": MAX, "zero": 0.0}[mode]
    return np.where(np.abs(x) == 0.0, fill, v)

def basis(W, X, mode):
    W = np.asarray(W).reshape(-1, X.shape[1])
    cols = []
    for k in range(len(W)):
        c = np.ones(len(X))
        for j in range(X.shape[1]): c = c*pu(X[:, j], W[k, j], mode)
        cols.append(c)
    return np.c_[tuple(cols) + (np.ones(len(X)),)]

def fit(X, y, mode, loss="linear", seeds=25):
    F = X.shape[1]; wt = 1.0/np.maximum(np.abs(y), 1e-300)
    def resid(w):
        Phi = basis(w, X, mode)
        if not np.all(np.isfinite(Phi)): return np.full(len(y), 1e12)
        c = np.linalg.lstsq(Phi*wt[:,None], y*wt, rcond=None)[0]
        r = (Phi @ c - y)*wt
        return np.where(np.isfinite(r), r, 1e12)
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
print("=" * 84)
print("② 当てはめ — 分解能で r=0 と記録された点がある（真 w(r) = −2）")
print("=" * 84)
print(f"  {'r=0 の割合':>11}{'Inf（IEEE）':>13}{'**±MAX**':>12}{'**0**':>10}{'MAX＋Cauchy':>14}")
for frac in (0.0, 0.01, 0.05, 0.15):
    X, Xo, yv, nb = make(2000, frac, seed=1)
    r = []
    for mode, loss in (("inf","linear"), ("max","linear"), ("zero","linear"), ("max","cauchy")):
        w = fit(Xo, yv, mode, loss)
        r.append(w[2] if np.isfinite(w[2]) else np.nan)
    print(f"  {frac:>10.0%}{r[0]:>13.4f}{r[1]:>12.4f}{r[2]:>10.4f}{r[3]:>14.4f}")
print(f"  {'真値':>11}{-2.0:>13.4f}{-2.0:>12.4f}{-2.0:>10.4f}{-2.0:>14.4f}")

print()
print("=" * 84)
print("③ 残差は有界か — r=0 の点が1つあるときの相対残差")
print("=" * 84)
Xq = np.array([[5.0, 5.0, 0.0]])
yq = np.array([1e6])
print(f"  {'規約':<14}{'模型の値':>16}{'相対残差 (v−y)/y':>20}")
for mode in ("inf", "max", "zero"):
    v = float(pu(Xq[:,0],1.,mode)[0]*pu(Xq[:,1],1.,mode)[0]*pu(Xq[:,2],-2.,mode)[0])*8.99
    rr = (v - yq[0])/yq[0]
    print(f"  {mode:<14}{v:>16.4g}{rr:>20.4g}")
print(f"\n  ⟹ **MAX は有限だが、残差は 1e302。0 だけが 残差を −1 で頭打ちにする。**")
