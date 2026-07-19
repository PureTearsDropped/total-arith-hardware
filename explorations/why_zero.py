#!/usr/bin/env python3
"""**`a/0 = 0` の 0 は、特別か。** 他の値を入れて、割る。

  `x^0 = 1`  は強制（空の積）／ `a×0 = 0` は定理（分配律）⟹ **評価が要るのは `a/0` だけ**。
  ゼロ除算は未定義なので **どんな値でもよい**。我々は 0 を選んでいる。**なぜ 0 か。**

  仮説（「白票」の話が正しいなら）: **0 は、設計行列の成分を消せる唯一の値**。
  設計が 0 ⟹ その行は `0·c = y` ⟹ **てこゼロ ⟹ 白票**。
  他のどんな値 k でも設計は `k` ⟹ **その行は本物の方程式 `k·c = y` になり、c を引っぱる**。
"""
import warnings, numpy as np
from scipy.optimize import minimize_scalar
warnings.filterwarnings("ignore")
MAX = np.finfo(np.float64).max

def make_pu(k):
    """`0^w := k`（w≠0）とする規約。k=0 が我々の選択。"""
    def pu(x, w):
        with np.errstate(all="ignore"):
            v = np.power(np.abs(x), w)
        return np.where(np.abs(x) == 0.0, k, v)
    return pu

def fit_w(r, y, k):
    pu = make_pu(k)
    def loss(w):
        with np.errstate(all="ignore"):
            A = pu(r, w)[:, None]
            wt = 1.0/np.maximum(np.abs(y), 1e-300)
            Aw, yw = A*wt[:, None], y*wt
            if not np.all(np.isfinite(Aw)): return 1e18
            c = np.linalg.lstsq(Aw, yw, rcond=None)[0]
            v = float(np.sum((Aw @ c - yw)**2))
            return v if np.isfinite(v) else 1e18
    return float(minimize_scalar(loss, bounds=(-6,6), method="bounded",
                                 options={"xatol":1e-10}).x)

print("="*92)
print("**`0^w := k` の k を変えて、w(r) を当てる**  （壊れた点 15%、真値 −2.0000、8 seed）")
print("="*92)
print(f"  {'k = 0^w の値':<24}{'**推定した w(r)**':>22}{'真値からのずれ':>18}")
CAND = [("**0**（我々の選択）", 0.0), ("1e−300", 1e-300), ("1e−30", 1e-30), ("1e−9", 1e-9),
        ("1（乗法の単位元）", 1.0), ("42", 42.0), ("1e+6", 1e6),
        ("**MAX**（飽和）", MAX), ("**inf**（正直）", np.inf)]
for tag, k in CAND:
    got = []
    for seed in range(8):
        rng = np.random.default_rng(seed)
        n = 2000
        r_true = np.abs(rng.normal(1.0, 0.3, n)) + 0.05
        y = 8.99e9*5.0/r_true**2 * np.exp(rng.normal(0, 0.01, n))
        idx = rng.choice(n, int(0.15*n), replace=False)
        r_rec = r_true.copy(); r_rec[idx] = 0.0
        got.append(fit_w(r_rec, y, k))
    m = np.median(got)
    print(f"  {tag:<24}{f'**{m:+.4f}**':>22}{f'{abs(m+2.0):.4f}':>18}")

print()
print("="*92); print("なぜ 0 だけか — 壊れた点が作る行を、そのまま見る"); print("="*92)
y_bad = 1.0e12; wt = 1.0/y_bad
print(f"  {'k':<16}{'行（設計 · c = 目標）':<30}{'この行が要求する c':>22}")
for tag, k in [("**0**", 0.0), ("1e−30", 1e-30), ("1", 1.0), ("42", 42.0)]:
    d = k*wt; t = y_bad*wt
    req = "**（要求できない）**" if d == 0 else f"{t/d:.4g}"
    print(f"  {tag:<16}{f'({d:.4g})·c = {t:.4g}':<30}{req:>22}")
print(f"""
  真の c = 8.99e9 × 5 = **4.5e+10**

  ⟹ **k = 0 のときだけ、設計が消えて「c への要求」が消える。** 白票。
     **k ≠ 0 なら、その行は本物の方程式になり、c を y/k へ引っぱる。**
     引っぱる力は k が小さいほど弱いが、**ゼロにはならない**。0 は極限ではなく、**別の状態**である。""")

# ================ 0 の取り柄は「尺度を知らなくてよい」ことか ================
print()
print("="*92)
print("**尺度を変える** — 本物の基底値 r^{−2} の桁を動かし、k を固定したまま測る")
print("="*92)
print(f"  {'r の典型値':>12}{'本物の基底 r^−2':>18}{'k=0':>14}{'k=1e−9':>14}{'k=1':>14}")
for rscale in (1e-3, 1.0, 1e3, 1e5):
    row = {}
    for k in (0.0, 1e-9, 1.0):
        got = []
        for seed in range(6):
            rng = np.random.default_rng(seed)
            n = 2000
            r_true = (np.abs(rng.normal(1.0, 0.3, n)) + 0.05) * rscale
            y = 8.99e9*5.0/r_true**2 * np.exp(rng.normal(0, 0.01, n))
            idx = rng.choice(n, int(0.15*n), replace=False)
            r_rec = r_true.copy(); r_rec[idx] = 0.0
            got.append(fit_w(r_rec, y, k))
        row[k] = np.median(got)
    base = 1.0/rscale**2
    f = lambda v: f"**{v:+.4f}**" if abs(v+2) < 0.01 else f"{v:+.4f} ✗"
    print(f"  {rscale:>12.0e}{base:>18.1e}{f(row[0.0]):>14}{f(row[1e-9]):>14}{f(row[1.0]):>14}")
print("""
  ⟹ **`0` は、データの尺度を知らなくても正しい唯一の選択である。**
     `k = 1e−9` は「小さい」のではなく「**この課題の基底値より小さい**」だけで、
     基底値が下がれば **同じ 1e−9 が崩壊する**。
     **`0` の取り柄はゼロであることではなく、`0` が どんな尺度より小さい ことである。**""")
