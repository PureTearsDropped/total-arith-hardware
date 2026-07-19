#!/usr/bin/env python3
"""**符号は生き残るか** — 演算回路としての測定。物理は関係ない。

真値は **Fraction（正確な有理数）**。+ − × ÷ は有理数体で厳密に閉じているので、
乱数で作った式木の **真の符号** が誤差ゼロで分かる。ゼロ除算が起きた木は除外
（真値が定義されないので、どちらの規約が正しいとも言えない）。

測るのはただ一つ: **答えの符号が、真の符号と一致するか。**
"""
import warnings, numpy as np
from fractions import Fraction as Fr
warnings.filterwarnings("ignore")
from total_arith import MAX, MIN, saturate

rng = np.random.default_rng(7)
OPS = "+-*/"

def leaf():
    """桁範囲の広い葉。物理の実測はこうなる（電荷 1e−19、質量 1e30 が同じ式に入る）。"""
    e = int(rng.integers(-180, 181))
    m = int(rng.integers(1, 999))
    s = -1 if rng.random() < 0.5 else 1
    return s * m * Fr(10) ** e

def tree(d):
    if d == 0: return ("leaf", leaf())
    return (OPS[rng.integers(4)], tree(d - 1), tree(d - 1))

def ev_exact(t):
    if t[0] == "leaf": return t[1]
    a, b = ev_exact(t[1]), ev_exact(t[2])
    if a is None or b is None: return None
    if t[0] == "/" and b == 0: return None           # 真値が定義されない → 木ごと除外
    return {"+": a + b, "-": a - b, "*": a * b, "/": (a / b if b else None)}[t[0]]

def ev_ieee(t):
    if t[0] == "leaf": return float(t[1])
    a, b = ev_ieee(t[1]), ev_ieee(t[2])
    with np.errstate(all="ignore"):
        a, b = np.float64(a), np.float64(b)
        return float({"+": a + b, "-": a - b, "*": a * b, "/": a / b}[t[0]])

def ev_total(t):
    """**溢れ→±MAX ／ 潰れ→±MIN ／ a/0=0**。0 は本物の 0 からしか出ない。"""
    if t[0] == "leaf": return float(saturate(float(t[1]), zero_ok=(t[1] == 0)))
    a, b = ev_total(t[1]), ev_total(t[2])
    with np.errstate(all="ignore"):
        if t[0] == "/":
            if b == 0.0: return 0.0
            v, zo = np.float64(a) / np.float64(b), (a == 0.0)
        elif t[0] == "*":
            v, zo = np.float64(a) * np.float64(b), (a == 0.0 or b == 0.0)
        else:
            v = np.float64(a) + np.float64(b) if t[0] == "+" else np.float64(a) - np.float64(b)
            zo = (float(v) == 0.0 and abs(a) == abs(b))     # 真に打ち消した和差だけ 0 を許す
        return float(saturate(v, zero_ok=zo))

def sgn_f(x):  return 0 if x == 0 else (1 if x > 0 else -1)
def sgn_ie(x):
    if np.isnan(x): return None                    # **NaN は符号を持たない ⟹ 失われた**
    return 0 if x == 0 else (1 if x > 0 else -1)

print("=" * 90)
print("**符号は生き残るか** — 乱数の式木 × 正確な有理数の真値")
print("=" * 90)
print(f"  {'木の深さ':>8}{'式の数':>8}{'IEEE 符号が正':>16}{'**全域 符号が正**':>20}{'IEEE が NaN':>14}")
for depth in (1, 2, 3, 4, 5):
    n = ok_i = ok_t = nan_i = 0
    tries = 0
    while n < 4000 and tries < 60000:
        tries += 1
        t = tree(depth)
        ex = ev_exact(t)
        if ex is None: continue                     # ゼロ除算を含む木は除外
        n += 1
        s_true = sgn_f(ex)
        vi, vt = ev_ieee(t), ev_total(t)
        si, st = sgn_ie(vi), sgn_ie(vt)
        nan_i += (si is None)
        ok_i += (si == s_true); ok_t += (st == s_true)
    print(f"  {depth:>8}{n:>8,}{f'{100*ok_i/n:.2f}%':>16}{f'**{100*ok_t/n:.2f}%**':>20}{f'{100*nan_i/n:.2f}%':>14}")
