#!/usr/bin/env python3
"""**0 だけが符号を持たない。他は「分からない」のではなく「持っている」。**（利用者の指摘）

  `0` の特別さは、代数的には一つの事実の二つの顔である:
      **0 は x ↦ −x の唯一の不動点**（0 = −0）  ⟹ **符号を持たない**
      **0 は x ↦ λx の唯一の不動点**（0 = λ·0） ⟹ **尺度を持たない**
  他のどんな k も、符号と尺度の**両方を主張する**。**0 だけが何も主張しない。**

  検証: `0^w := k` の **k の符号を反転**して、答えが変わるか。
        変われば「k は符号を主張しており、主張には帰結がある」＝ **0 以外は嘘をつく**。
"""
import warnings, numpy as np
from scipy.optimize import minimize_scalar
warnings.filterwarnings("ignore")
MAX = np.finfo(np.float64).max

def fit_w(r, y, k):
    fab = (np.abs(r) == 0.0)
    def loss(w):
        with np.errstate(all="ignore"):
            v = np.power(np.abs(r), w)
            A = np.where(fab, k, v)[:, None]
            wt = 1.0/np.maximum(np.abs(y), 1e-300)
            Aw, yw = A*wt[:, None], y*wt
            if not np.all(np.isfinite(Aw)): return 1e18
            c = np.linalg.lstsq(Aw, yw, rcond=None)[0]
            vv = float(np.sum((Aw @ c - yw)**2))
            return vv if np.isfinite(vv) else 1e18
    return float(minimize_scalar(loss, bounds=(-6,6), method="bounded", options={"xatol":1e-10}).x)

def run(k):
    got = []
    for seed in range(6):
        rng = np.random.default_rng(seed); n = 2000
        r_true = np.abs(rng.normal(1.0, 0.3, n)) + 0.05
        y = 8.99e9*5.0/r_true**2 * np.exp(rng.normal(0, 0.01, n))
        idx = rng.choice(n, int(0.15*n), replace=False)
        r_rec = r_true.copy(); r_rec[idx] = 0.0
        got.append(fit_w(r_rec, y, k))
    return np.median(got)

print("="*92)
print("**k の符号を反転すると、答えは変わるか** — 変われば「k は符号を主張している」")
print("="*92)
print(f"  {'|k|':<16}{'k = **+|k|**':>18}{'k = **−|k|**':>18}{'**差**':>16}   判定")
for mag in (1e-9, 1e-3, 1.0, 42.0, 1e6, MAX):
    a, b = run(+mag), run(-mag)
    d = abs(a-b)
    tag = "**主張していない**" if d < 1e-3 else "**符号を主張している**"
    lab = "MAX" if mag == MAX else f"{mag:g}"
    print(f"  {lab:<16}{f'{a:+.4f}':>18}{f'{b:+.4f}':>18}{f'**{d:.4f}**':>16}   {tag}")
z = run(0.0)
print(f"  {'**0**':<16}{f'**{z:+.4f}**':>18}{f'**{run(-0.0):+.4f}**':>18}{'**0.0000**':>16}   **主張しようがない（0 = −0）**")

print()
print("="*92); print("代数: 不動点で数え上げる"); print("="*92)
xs = np.array([-1e300, -42.0, -1e-9, 0.0, 1e-9, 42.0, 1e300])
neg_fix   = xs[xs == -xs]
scal_fix  = xs[np.all([xs == l*xs for l in (0.5, 2.0, 7.0)], axis=0)]
print(f"  x ↦ −x の不動点  : **{list(neg_fix)}**   ⟹ **符号を持たない値は 0 だけ**")
print(f"  x ↦ λx の不動点  : **{list(scal_fix)}**   ⟹ **尺度を持たない値は 0 だけ**")
print(f"""
  ⟹ **二つの性質は別々の発見ではなく、同じ一点の二つの顔である。**
     `0` は **符号の主張** と **尺度の主張** の両方から降りている唯一の数。

  ⟹ そして利用者の指摘の要点:
     **`+MAX` は「符号が分からない」のではない。`+` を主張している。**
     真値は `a/0` ＝ 右から近づけば `+∞`、左から近づけば `−∞` ⟹ **符号は本当に定まらない**。
     **`+MAX` はそこに `+` と書いてしまう。`0` は書かない。**""")
