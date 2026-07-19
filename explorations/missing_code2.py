#!/usr/bin/env python3
"""**8番目の符号 = 欠測**。§7.5.7 と同じ土俵（相対残差・非線形）で測り直す。

前版 `missing_code.py` は **測れていなかった**（対数線形の `y>0` マスクが、どの規約でも壊れた点を
等しく落としていた ⟹ 3規約が同じ数字）。**対数線形は y=0 を構造的に扱えない。**

ここでは §7.5.7 の目的関数で当てる:
    規約 0  : y_bad = 0 ⟹ 相対残差 (0 − y)/y = **−1** が損失に **入る**（有界な嘘）
    **欠測**: 印を読んで損失に **入れない**
"""
import warnings, numpy as np
from scipy.optimize import minimize_scalar
warnings.filterwarnings("ignore")

def power_unit_zero(x, w):
    """§7.5.7 の規約: 0^w = 0（w≠0）"""
    with np.errstate(all="ignore"):
        v = np.power(np.abs(x), w)
    return np.where(np.abs(x) == 0.0, 0.0, v)

def fit_w(r, y, mask, mode):
    """相対残差で w を当てる。mask=False の点は **欠測** として損失から外す（mode='drop' のとき）。"""
    use = np.ones(len(y), bool) if mode != "drop" else mask
    def loss(w):
        with np.errstate(all="ignore"):
            b = power_unit_zero(r, w)
            A = np.c_[b]
            ok = np.isfinite(A).all(1) & np.isfinite(y) & use
            if ok.sum() < 5: return 1e18
            # VarPro: 係数 c は閉形式（相対重み 1/|y|）
            wt = 1.0 / np.maximum(np.abs(y[ok]), 1e-300)
            Aw, yw = A[ok]*wt[:,None], y[ok]*wt
            c = np.linalg.lstsq(Aw, yw, rcond=None)[0]
            rr = Aw @ c - yw
            v = float(np.sum(rr**2))
            return v if np.isfinite(v) else 1e18
    return float(minimize_scalar(loss, bounds=(-6, 6), method="bounded",
                                 options={"xatol": 1e-8}).x)

print("="*92)
print("**8番目の符号 = 欠測**  — §7.5.7 と同じ相対残差の土俵（真値 w(r) = −2.0000）")
print("="*92)
print(f"  {'壊れた点':>10}{'規約 0（有界な嘘・損失に入る）':>30}{'**欠測の印（損失に入れない）**':>32}")
for frac in (0.01, 0.05, 0.15, 0.30, 0.50, 0.70):
    gz, gd = [], []
    for seed in range(8):
        rng = np.random.default_rng(seed)
        n = 400
        r = np.abs(rng.normal(1.0, 0.3, n)) + 0.05
        y = 8.99e9 * 5.0 / r**2 * np.exp(rng.normal(0, 0.01, n))
        k = rng.choice(n, int(frac*n), replace=False)
        r2 = r.copy(); r2[k] = 0.0
        y2 = y.copy(); y2[k] = 0.0            # **分解能以下 ⟹ 0 と記録された**
        mask = np.ones(n, bool); mask[k] = False   # **8番目の符号: この点は測れなかった**
        gz.append(fit_w(r2, y2, mask, "zero"))
        gd.append(fit_w(r2, y2, mask, "drop"))
    print(f"  {f'{100*frac:.0f}%':>10}{f'{np.median(gz):+.4f}':>30}{f'**{np.median(gd):+.4f}**':>32}")
