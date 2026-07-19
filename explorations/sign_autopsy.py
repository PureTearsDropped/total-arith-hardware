#!/usr/bin/env python3
"""全域規約が符号を落とす **1.6%** は、どの演算のせいか。＝「MIN と MAX をどうするか」の問題。"""
import warnings, numpy as np
from fractions import Fraction as Fr
from collections import Counter
warnings.filterwarnings("ignore")
from total_arith import MAX, MIN, saturate
rng = np.random.default_rng(7); OPS = "+-*/"

def leaf():
    return (-1 if rng.random() < .5 else 1) * int(rng.integers(1, 999)) * Fr(10) ** int(rng.integers(-180, 181))
def tree(d): return ("leaf", leaf()) if d == 0 else (OPS[rng.integers(4)], tree(d-1), tree(d-1))
def ev_exact(t):
    if t[0] == "leaf": return t[1]
    a, b = ev_exact(t[1]), ev_exact(t[2])
    if a is None or b is None or (t[0] == "/" and b == 0): return None
    return {"+": a+b, "-": a-b, "*": a*b, "/": (a/b if b else None)}[t[0]]

BLAME = Counter()
def ev_total(t, exact):
    """評価しながら、**真値は 0 でないのに 0 を返した**演算を記録する。"""
    if t[0] == "leaf": return float(t[1])
    a = ev_total(t[1], ev_exact(t[1])); b = ev_total(t[2], ev_exact(t[2]))
    with np.errstate(all="ignore"):
        if t[0] == "/":
            if b == 0.0: v, zo = 0.0, True
            else: v, zo = np.float64(a)/np.float64(b), (a == 0.0)
        elif t[0] == "*": v, zo = np.float64(a)*np.float64(b), (a == 0.0 or b == 0.0)
        else:
            v = np.float64(a)+np.float64(b) if t[0]=="+" else np.float64(a)-np.float64(b)
            zo = (float(v) == 0.0 and abs(a) == abs(b))
        r = float(saturate(v, zero_ok=zo))
    if exact is not None and exact != 0 and r == 0.0:
        big = abs(a) >= MAX*0.999 or abs(b) >= MAX*0.999
        BLAME[f"{t[0]}  " + ("**MAX 同士が打ち消した**" if big and t[0] in "+-"
              else "**MAX で割った/掛けた**" if big else "本物の 0 で割った" if t[0]=="/" and b==0.0
              else "その他")] += 1
    return r

n = ok = fail = 0; tries = 0
while n < 6000 and tries < 200000:
    tries += 1; t = tree(5); ex = ev_exact(t)
    if ex is None: continue
    n += 1
    r = ev_total(t, ex)
    s_true = 0 if ex == 0 else (1 if ex > 0 else -1)
    s_got = 0 if r == 0 else (1 if r > 0 else -1)
    ok += (s_true == s_got); fail += (s_true != s_got)

print("=" * 88); print(f"深さ5・{n:,} 式。全域規約が符号を落とした **{fail} 式（{100*fail/n:.2f}%）** の原因"); print("=" * 88)
for k, v in BLAME.most_common(6):
    print(f"  {v:>6,} 回   {k}")
print(f"""
  ⟹ **原因は一つ**: **MAX 同士の引き算・足し算**。1e+308 級の巨大な値が2つ来て、
     引くと **0** になる。真値は巨大なのに、**0 と答えてしまう**。

     これは Inf の `Inf − Inf = NaN` と **同じ穴**である。IEEE は「分からない」と言い、
     全域規約は「**0 だ**」と言う。どちらも真値（巨大）ではない。
     **飽和が情報を捨てた後では、引き算は復元できない。**

  ⟹ だからあなたの問い「**MIN と MAX の演算をどうするか**」の答えは:
     **MAX − MAX を 0 にしてはいけない。** ここだけは 0 が本物の 0 でない。
     候補は 3 つ:
       ① **MAX − MAX = MAX**（飽和を吸収的にする。符号は左辺から）
       ② **MAX に「飽和した」の 1 ビットを付ける**（= 汚染フラグ。NaN の軽量版）
       ③ **そのまま 0**（今の実装。1.6% で符号が死ぬ）""")
