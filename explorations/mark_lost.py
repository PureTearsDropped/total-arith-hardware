#!/usr/bin/env python3
"""**印が失われたとき、どの値が生き残るか。** ＝ `a/0 = 0` が、印があってもなお意味を持つ理由。

印は **読まれるとは限らない**（古いコード／印を知らないライブラリ／ファイルに書いて読み直す／他言語へ渡す）。
値は **必ず読まれる**。だから値は「**印が失われても最悪にならない**」ものを選ぶべきである。
"""
import warnings, numpy as np
from scipy.optimize import minimize_scalar
warnings.filterwarnings("ignore")
MAX = np.finfo(np.float64).max

def fit_w(r, y, k, mark_frac, rng):
    """mark_frac: 印が生き残った割合（1.0=全部届く、0.0=全部失われた）"""
    is_fab = (np.abs(r) == 0.0)
    survived = is_fab & (rng.random(len(r)) < mark_frac)      # **印が届いた作り物だけ落とせる**
    def loss(w):
        with np.errstate(all="ignore"):
            v = np.power(np.abs(r), w)
            A = np.where(is_fab, k, v)[:, None]
            wt = 1.0/np.maximum(np.abs(y), 1e-300)
            Aw, yw = A*wt[:, None], y*wt
            keep = ~survived
            Aw, yw = Aw[keep], yw[keep]
            if not np.all(np.isfinite(Aw)) or len(yw) < 5: return 1e18
            c = np.linalg.lstsq(Aw, yw, rcond=None)[0]
            vv = float(np.sum((Aw @ c - yw)**2))
            return vv if np.isfinite(vv) else 1e18
    return float(minimize_scalar(loss, bounds=(-6,6), method="bounded", options={"xatol":1e-10}).x)

def run(k, mark_frac, rscale=1.0):
    got = []
    for seed in range(6):
        rng = np.random.default_rng(seed); n = 2000
        r_true = (np.abs(rng.normal(1.0, 0.3, n)) + 0.05) * rscale
        y = 8.99e9*5.0/r_true**2 * np.exp(rng.normal(0, 0.01, n))
        idx = rng.choice(n, int(0.15*n), replace=False)
        r_rec = r_true.copy(); r_rec[idx] = 0.0
        got.append(fit_w(r_rec, y, k, mark_frac, np.random.default_rng(seed+100)))
    return np.median(got)

f = lambda v: f"**{v:+.4f}**" if abs(v+2) < 0.02 else f"{v:+.4f} ✗"
print("="*94)
print("**印が失われていく** — どの値が最後まで生き残るか（真値 −2.0000）")
print("="*94)
print(f"  {'印が届く割合':>14}{'k = **0**':>16}{'k = 1e−9':>16}{'k = 42':>14}{'k = **+MAX**':>18}")
for mf in (1.0, 0.9, 0.5, 0.1, 0.0):
    print(f"  {f'{100*mf:.0f}%':>14}{f(run(0.0,mf)):>16}{f(run(1e-9,mf)):>16}"
          f"{f(run(42.0,mf)):>14}{f(run(MAX,mf)):>18}")
print("""
  ⟹ **印が完全に届く（100%）なら、どの値でも同じ。** 値は自由。
     **印が一枚でも落ちると、値の質が効き始める。**
     **印が全部落ちた（0%）とき、生き残るのは 0 だけ。**""")
print()
print("="*94)
print("そして尺度を知らない場合（印 0%）")
print("="*94)
print(f"  {'r の典型値':>12}{'k = **0**':>16}{'k = 1e−9':>16}{'k = **+MAX**':>18}")
for rs in (1.0, 1e3, 1e5):
    print(f"  {rs:>12.0e}{f(run(0.0,0.0,rs)):>16}{f(run(1e-9,0.0,rs)):>16}{f(run(MAX,0.0,rs)):>18}")
