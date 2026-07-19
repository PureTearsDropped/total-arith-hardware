#!/usr/bin/env python3
"""**8番目の符号 = 欠測**。§7.5.7 の状況を正しく作って測る（3度目）。

    **欠けているのは `r` だけ**。`r` が分解能以下で **0 と記録されて**も、`y` は巨大な値として
    **ちゃんと測れている**。だから:
        規約 0  : 模型は c·0^w = **0** を予測。実測 y は巨大 ⟹ 相対残差 **(0−y)/y = −1** が損失に入る
        **欠測**: `r` に印がある ⟹ その点を損失に **入れない**

    前2版の誤り: 監査人が `y` にも 0 を立てたので **両辺 0 の空の式**になり、
                 最小二乗の行が自動的に消えていた（＝壊れた点が最初から損失に居なかった）。
"""
import warnings, numpy as np
from scipy.optimize import minimize_scalar
warnings.filterwarnings("ignore")

def pu_zero(x, w):
    with np.errstate(all="ignore"):
        v = np.power(np.abs(x), w)
    return np.where(np.abs(x) == 0.0, 0.0, v)      # **0^w = 0（w<0 は有界な嘘）**
def pu_inf(x, w):
    with np.errstate(all="ignore"):
        return np.power(np.abs(x), w)              # 正直: 0^{-2} = inf

def fit_w(r, y, mask, mode):
    pu = pu_inf if mode == "inf" else pu_zero
    use = mask if mode == "drop" else np.ones(len(y), bool)
    def loss(w):
        with np.errstate(all="ignore"):
            A = pu(r, w)[:, None]
            ok = np.isfinite(A).all(1) & np.isfinite(y) & use
            if ok.sum() < 5: return 1e18
            wt = 1.0/np.maximum(np.abs(y[ok]), 1e-300)
            Aw, yw = A[ok]*wt[:,None], y[ok]*wt
            c = np.linalg.lstsq(Aw, yw, rcond=None)[0]
            v = float(np.sum((Aw @ c - yw)**2))
            return v if np.isfinite(v) else 1e18
    return float(minimize_scalar(loss, bounds=(-6,6), method="bounded", options={"xatol":1e-8}).x)

print("="*96)
print("**8番目の符号 = 欠測**  — `r` だけが欠け、`y` は測れている（真値 w(r) = −2.0000）")
print("="*96)
print(f"  {'壊れた点':>10}{'規約 ∞（正直）':>18}{'規約 0（有界な嘘）':>22}{'**欠測の印（8番目）**':>26}")
for frac in (0.01, 0.05, 0.15, 0.30, 0.50):
    g = {m: [] for m in ("inf","zero","drop")}
    for seed in range(8):
        rng = np.random.default_rng(seed)
        n = 400
        r_true = np.abs(rng.normal(1.0, 0.3, n)) + 0.05
        y = 8.99e9*5.0/r_true**2 * np.exp(rng.normal(0, 0.01, n))   # **y は真の r で測れている**
        k = rng.choice(n, int(frac*n), replace=False)
        r_rec = r_true.copy(); r_rec[k] = 0.0                        # **r だけが 0 と記録された**
        mask = np.ones(n, bool); mask[k] = False                     # **8番目の符号**
        for m in ("inf","zero","drop"): g[m].append(fit_w(r_rec, y, mask, m))
    f = lambda v: f"{np.median(v):+.4f}"
    print(f"  {f'{100*frac:.0f}%':>10}{f(g['inf']):>18}{f(g['zero']):>22}{f'**{f(g[chr(100)+chr(114)+chr(111)+chr(112)])}**':>26}")
