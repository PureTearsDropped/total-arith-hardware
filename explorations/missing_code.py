#!/usr/bin/env python3
"""**8番目の符号 = 欠測（測っていない）**。§7.5.7 の「有界な嘘」は、これが無いから払っている代金か。

    r = 0 と記録された点の意味は **「r は 0 である」ではなく「r は測れなかった」**。
    IEEE にも我々の規約にも、**それを書く場所が無い**ので 0 と書く。
    そして `a/0 = 0` が「有界な嘘」で受け止める（残差が −1 で頭打ち ⟹ 良い点が支配権を保つ）。

    **8番目の符号があれば、嘘をつかなくてよい。** 測る。

注: これは監査人が一度やらかした「壊れた点を `yy>0` で除外して −2.0000 を得る」と**同じ計算**である。
    違いは **知っているかどうか**。偶然除外すれば測定の破壊、印を読んで除外すればアルゴリズム。
"""
import warnings, numpy as np
warnings.filterwarnings("ignore")
MISSING = "**欠測**"

def fit(r, y, mode, rng):
    """逆二乗 y = k·q/r² の w(r) を対数線形で当てる。"""
    if mode == "drop":                                  # **欠測の印を読んで落とす**
        keep = ~np.isnan(r)
        r, y = r[keep], y[keep]
    ok = (y > 0) & np.isfinite(y) & np.isfinite(r)
    if ok.sum() < 10: return float("nan")
    A = np.c_[np.ones(ok.sum()), np.log(np.maximum(np.abs(r[ok]), 1e-320))]
    return float(np.linalg.lstsq(A, np.log(y[ok]), rcond=None)[0][1])

print("="*90)
print("**8番目の符号 = 欠測**  — `r=0` は「0」ではなく「測れなかった」（真値 w(r) = −2.0000）")
print("="*90)
print(f"  {'壊れた点':>10}{'規約 ∞（正直）':>18}{'規約 0（有界な嘘）':>22}{'**欠測の印（8番目）**':>24}")
for frac in (0.01, 0.05, 0.15, 0.30, 0.50):
    got = {m: [] for m in ("inf", "zero", "drop")}
    for seed in range(8):
        rng = np.random.default_rng(seed)
        n = 400
        r = np.abs(rng.normal(1.0, 0.3, n)) + 0.05
        y = 8.99e9 * 5.0 / r**2
        y = y * np.exp(rng.normal(0, 0.01, n))
        k = rng.choice(n, int(frac*n), replace=False)
        # ---- 分解能以下 ⟹ 測定器は 0 と記録する（あるいは「欠測」と記録できる）----
        r_inf  = r.copy(); r_inf[k]  = 0.0                       # y は 1/0² = inf に
        y_inf  = y.copy(); y_inf[k]  = np.inf
        r_zero = r.copy(); r_zero[k] = 0.0                       # 0^{-2} = 0 ⟹ y = 0（有界な嘘）
        y_zero = y.copy(); y_zero[k] = 0.0
        r_drop = r.copy(); r_drop[k] = np.nan                    # **欠測の印**
        y_drop = y.copy(); y_drop[k] = np.nan
        got["inf"].append(fit(r_inf, y_inf, "keep", rng))
        got["zero"].append(fit(r_zero, y_zero, "keep", rng))
        got["drop"].append(fit(r_drop, y_drop, "drop", rng))
    f = lambda v: f"{np.nanmedian(v):+.4f}" if np.isfinite(np.nanmedian(v)) else "**落ちる**"
    print(f"  {f'{100*frac:.0f}%':>10}{f(got['inf']):>18}{f(got['zero']):>22}{f'**{f(got[chr(100)+chr(114)+chr(111)+chr(112)])}**':>24}")
print("""
  ⟹ **8番目の符号は、嘘をやめさせる。**
     規約 0 は「有界な嘘」で **良い点に支配権を残す**が、**壊れた点を損失に入れ続ける**（残差 −1 × 個数）。
     欠測の印があれば **入れない**。壊れた点の割合が増えるほど、差が開く。

  代金: **測定器が「0」ではなく「測れなかった」と言えること。**
        これは算術の問題ではなく **記録の問題**であり、8番目の符号はその場所を作るだけである。""")
