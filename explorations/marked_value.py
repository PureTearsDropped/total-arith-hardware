#!/usr/bin/env python3
"""**印があれば、値は何でもよくなるか。** §7.6.7⑤ の続き。

`why_zero.py` の結論「**0 は尺度を持たない唯一の数だから最良**」は、**印が無い世界での結論**だった。
3ビット符号（§7.6.7④）があれば `大きさ不明` の印が情報を運ぶので、**値が運ぶ必要が無い**。

そして、それが効くなら §7.5.7 の対立（**正直な ∞ は死ぬ／嘘の 0 は生きる**）が消える:
`r=0` での真値は `+∞` ⟹ **`+MAX` を返して印を立てれば、向きは正直なまま、行は無視できる。**
"""
import warnings, numpy as np
from scipy.optimize import minimize_scalar
warnings.filterwarnings("ignore")
MAX = np.finfo(np.float64).max

def fit_w(r, y, k, use_mark):
    """`0^w := k`。use_mark=True なら **印の立った行を損失から外す**（値は見ない）。"""
    mark = (np.abs(r) == 0.0)                    # **大きさ不明の印**（作り物の値である）
    def loss(w):
        with np.errstate(all="ignore"):
            v = np.power(np.abs(r), w)
            A = np.where(mark, k, v)[:, None]
            wt = 1.0/np.maximum(np.abs(y), 1e-300)
            Aw, yw = A*wt[:, None], y*wt
            keep = ~mark if use_mark else np.ones(len(y), bool)
            Aw, yw = Aw[keep], yw[keep]
            if not np.all(np.isfinite(Aw)) or len(yw) < 5: return 1e18
            c = np.linalg.lstsq(Aw, yw, rcond=None)[0]
            vv = float(np.sum((Aw @ c - yw)**2))
            return vv if np.isfinite(vv) else 1e18
    return float(minimize_scalar(loss, bounds=(-6,6), method="bounded",
                                 options={"xatol":1e-10}).x)

def run(k, use_mark, rscale):
    got = []
    for seed in range(6):
        rng = np.random.default_rng(seed)
        n = 2000
        r_true = (np.abs(rng.normal(1.0, 0.3, n)) + 0.05) * rscale
        y = 8.99e9*5.0/r_true**2 * np.exp(rng.normal(0, 0.01, n))
        idx = rng.choice(n, int(0.15*n), replace=False)
        r_rec = r_true.copy(); r_rec[idx] = 0.0
        got.append(fit_w(r_rec, y, k, use_mark))
    return np.median(got)

print("="*94)
print("**印を立てると、値の選択は自由になるか**  （壊れた点 15%、真値 −2.0000）")
print("="*94)
print(f"  {'0^w に入れる値':<26}{'印なし（値だけ）':>22}{'**印あり（3ビット）**':>26}")
for tag, k in [("**0**（尺度を持たない）", 0.0), ("1e−9", 1e-9), ("1", 1.0), ("42", 42.0),
               ("1e+6", 1e6), ("**+MAX**（真値の向き）", MAX)]:
    a, b = run(k, False, 1.0), run(k, True, 1.0)
    f = lambda v: f"**{v:+.4f}**" if abs(v+2) < 0.01 else f"{v:+.4f} ✗"
    print(f"  {tag:<26}{f(a):>22}{f(b):>26}")

print()
print("="*94)
print("**尺度を動かす** — 印があれば、尺度を知らなくてよくなるか")
print("="*94)
print(f"  {'r の典型値':>12}{'k=0 印なし':>16}{'k=1e−9 印なし':>18}{'**k=+MAX 印あり**':>24}")
for rs in (1.0, 1e3, 1e5):
    a, b, c = run(0.0, False, rs), run(1e-9, False, rs), run(MAX, True, rs)
    f = lambda v: f"**{v:+.4f}**" if abs(v+2) < 0.02 else f"{v:+.4f} ✗"
    print(f"  {rs:>12.0e}{f(a):>16}{f(b):>18}{f(c):>24}")
