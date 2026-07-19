#!/usr/bin/env python3
"""二つの説明は、違う予言をする。**壊れた点を極端に増やして、割る。**

  「**有界な嘘**」（本文書が書いていた話）: 壊れた点も小声で押している ⟹ **増えればいつか勝つ**
  「**てこがゼロ**」（実際）              : 壊れた点は白票 ⟹ **何枚あっても結果は変わらない**
"""
import warnings, numpy as np
from scipy.optimize import minimize_scalar
warnings.filterwarnings("ignore")

def pu_zero(x, w):
    with np.errstate(all="ignore"):
        v = np.power(np.abs(x), w)
    return np.where(np.abs(x) == 0.0, 0.0, v)          # **0^w = 0**

def fit_w(r, y):
    """**壊れた点を一つも除外しない。** 相対重み。§7.5.7 と同じ目的関数。"""
    def loss(w):
        with np.errstate(all="ignore"):
            A = pu_zero(r, w)[:, None]
            wt = 1.0/np.maximum(np.abs(y), 1e-300)
            Aw, yw = A*wt[:, None], y*wt
            c = np.linalg.lstsq(Aw, yw, rcond=None)[0]
            v = float(np.sum((Aw @ c - yw)**2))
            return v if np.isfinite(v) else 1e18
    return float(minimize_scalar(loss, bounds=(-6, 6), method="bounded",
                                 options={"xatol": 1e-10}).x), loss

print("="*94)
print("**白票は何枚あっても結果を変えないか** — 壊れた点を 1% から 99% まで（真値 −2.0000）")
print("="*94)
print(f"  {'壊れた点':>10}{'良い点の数':>12}{'**推定した w(r)**':>22}{'損失に居座る定数':>18}")
for frac in (0.01, 0.15, 0.50, 0.80, 0.90, 0.95, 0.99):
    got, konst = [], []
    for seed in range(8):
        rng = np.random.default_rng(seed)
        n = 4000
        r_true = np.abs(rng.normal(1.0, 0.3, n)) + 0.05
        y = 8.99e9*5.0/r_true**2 * np.exp(rng.normal(0, 0.01, n))
        k = rng.choice(n, int(frac*n), replace=False)
        r_rec = r_true.copy(); r_rec[k] = 0.0
        w, loss = fit_w(r_rec, y)
        got.append(w); konst.append(len(k) * 1.0)     # 壊れた点は 1 点あたり残差² = 1
    ngood = n - int(frac*n)
    print(f"  {f'{100*frac:.0f}%':>10}{ngood:>12,}{f'**{np.median(got):+.4f}**':>22}{f'{np.median(konst):,.0f}':>18}")
print()
print("="*94); print("損失の形を直接見る — 定数を引いたら、同じ関数か"); print("="*94)
rng = np.random.default_rng(0); n = 4000
r_true = np.abs(rng.normal(1.0, 0.3, n)) + 0.05
y = 8.99e9*5.0/r_true**2 * np.exp(rng.normal(0, 0.01, n))
k = rng.choice(n, int(0.5*n), replace=False)
r_rec = r_true.copy(); r_rec[k] = 0.0
_, loss_bad  = fit_w(r_rec, y)
_, loss_good = fit_w(r_true[~np.isin(np.arange(n), k)], y[~np.isin(np.arange(n), k)])
print(f"  {'w':>8}{'壊れ 50% の損失':>20}{'良い点だけの損失':>20}{'**差**':>16}")
for w in (-2.5, -2.0, -1.5, -1.0, 0.0, 1.0):
    a, b = loss_bad(w), loss_good(w)
    print(f"  {w:>8.1f}{a:>20.6f}{b:>20.6f}{f'**{a-b:.6f}**':>16}")
