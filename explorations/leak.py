#!/usr/bin/env python3
"""**白票は、浮動小数でも本当に白票か。** §7.6.7⑤ の「厳密に不可視」を疑う。

数学的には `0·c = y` の行は擬似逆行列に何も寄与しない（設計が 0 ⟹ ΦᵀΦ にも Φᵀy にも入らない）。
**だが `np.linalg.lstsq` は 2.60003 を返し、厳密解 13/5 = 2.6 と違った。** どこから漏れたか。
"""
import numpy as np
from fractions import Fraction as Fr

print("="*92); print("① 壊れた点の y を大きくしていく — 答えは動くか（動いてはならない）"); print("="*92)
print(f"  {'壊れた点の y':>14}{'厳密解':>12}{'**lstsq**':>18}{'**pinv**':>18}{'normal eq':>14}")
for ybad in (1.0, 1e6, 1e12, 1e18, 1e30, 1e100, 1e300):
    Phi = np.array([[0.0],[1.0],[2.0]])
    y = np.array([ybad, 3.0, 5.0])
    exact = Fr(13, 5)                                    # ΦᵀΦ c = Φᵀy → 5c = 13
    a = np.linalg.lstsq(Phi, y, rcond=None)[0][0]
    b = (np.linalg.pinv(Phi) @ y)[0]
    n = float(np.linalg.solve(Phi.T@Phi, Phi.T@y)[0])    # 正規方程式
    print(f"  {ybad:>14.0e}{float(exact):>12.4f}{a:>18.6f}{b:>18.6f}{n:>14.6f}")
print("""
  ⟹ **lstsq だけが動く。** pinv と正規方程式は動かない。
     ⟹ **漏れは数学ではなく、lstsq の実装（SVD ドライバ gelsd）にある。**
        設計が 0 でも、SVD は y 全体を回転させるので **巨大な成分が丸め誤差として混ざる。**""")

print(); print("="*92); print("② 我々の当てはめは lstsq を使っている（`separable_fit.py`）"); print("="*92)
import subprocess
r = subprocess.run(["grep","-c","lstsq","separable_fit.py"], capture_output=True, text=True)
print(f"  `separable_fit.py` の lstsq 呼び出し: **{r.stdout.strip()} 箇所**")
print("""
  ⟹ **§7.6.7⑤ の「厳密に不可視」は、数学としては正しいが、我々の実装では正しくない。**
     `blank_ballot.py` が 99% 壊れで **−1.9983**（−2.0000 でない）だったのは、
     **良い点が 40 個しかないからではなく、この漏れかもしれない。** 切り分ける。""")

print(); print("="*92); print("③ 切り分け — lstsq を pinv に替えたら、漏れは消えるか"); print("="*92)
from scipy.optimize import minimize_scalar
def fit_w(r, y, solver):
    def loss(w):
        with np.errstate(all="ignore"):
            v = np.power(np.abs(r), w)
            A = np.where(np.abs(r)==0.0, 0.0, v)[:, None]
            wt = 1.0/np.maximum(np.abs(y), 1e-300)
            Aw, yw = A*wt[:,None], y*wt
            if solver == "lstsq": c = np.linalg.lstsq(Aw, yw, rcond=None)[0]
            else:                 c = np.linalg.pinv(Aw) @ yw
            vv = float(np.sum((Aw @ c - yw)**2))
            return vv if np.isfinite(vv) else 1e18
    return float(minimize_scalar(loss, bounds=(-6,6), method="bounded", options={"xatol":1e-10}).x)
print(f"  {'壊れた点':>10}{'良い点':>10}{'**lstsq**':>16}{'**pinv**':>16}{'差':>12}")
for frac in (0.15, 0.90, 0.99):
    a, b = [], []
    for seed in range(6):
        rng = np.random.default_rng(seed); n = 4000
        rt = np.abs(rng.normal(1.0, 0.3, n)) + 0.05
        y = 8.99e9*5.0/rt**2 * np.exp(rng.normal(0, 0.01, n))
        idx = rng.choice(n, int(frac*n), replace=False)
        rr = rt.copy(); rr[idx] = 0.0
        a.append(fit_w(rr, y, "lstsq")); b.append(fit_w(rr, y, "pinv"))
    ma, mb = np.median(a), np.median(b)
    print(f"  {f'{100*frac:.0f}%':>10}{n-int(frac*n):>10,}{f'{ma:+.4f}':>16}{f'{mb:+.4f}':>16}{abs(ma-mb):>12.4f}")
