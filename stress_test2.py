#!/usr/bin/env python3
"""さらに色んな計算 — 混在演算木・代数法則・極端な桁範囲・深い混在連鎖。

厳密な有理数を真値に、状態が 嘘をつかないか（健全性）と、
成り立つべき代数法則（分配律など）が 成り立つかを 見る。
"""
import sys, os
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
import numpy as np
from bfp_sed import BF, mul, div, add, EQ, GE, LE, NB, M
from sedenion_tensor_logic import ref_mult

def tvals(bf): return [Fr(m) * Fr(2) ** bf.E for m in bf.mant]
def sconjm(m): return [m[0]] + [-x for x in m[1:]]
def rmul(x, y): return ref_mult([Fr(a) for a in x], [Fr(b) for b in y])
def radd(x, y): return [a + b for a, b in zip(x, y)]

def sound(bf, truth):
    """健全性違反の件数（状態が 真値について 嘘をついた成分数）。"""
    rep = tvals(bf); v = 0
    for k in range(M):
        t, r, b = truth[k], rep[k], bf.bound[k]
        if b == EQ and t != r: v += 1
        elif b == GE and abs(t) < abs(r): v += 1
        elif b == LE and abs(t) > abs(r): v += 1
        if not bf.sunk[k] and t != 0 and r != 0 and (t > 0) != (r > 0): v += 1
    return v

rng = np.random.default_rng(4242)
report = []

# ---------------------------------------------- I. 混在演算木（mul/add/div をランダムに合成）
viol = 0
for trial in range(300):
    leaves = [BF([int(v) for v in rng.integers(-6, 6, M)], E=int(rng.integers(-3, 3)))
              for _ in range(4)]
    tru = [tvals(l) for l in leaves]
    cur, ctru = leaves[0], tru[0]
    for step in range(4):
        op = rng.integers(0, 3); nxt = leaves[rng.integers(0, 4)]; nt = tru[leaves.index(nxt)] if nxt in leaves else None
        j = int(rng.integers(0, 4)); nxt = leaves[j]; nt = tru[j]
        if op == 0:
            cur = mul(cur, nxt, W=28); ctru = rmul(ctru, nt)
        elif op == 1:
            cur = add(cur, nxt, W=28); ctru = radd(ctru, nt)
        else:
            N = sum(x * x for x in [int(round(float(t))) for t in nt])  # 近似 N（真値は 下で厳密に）
            Nf = sum(t * t for t in nt)
            if Nf == 0:
                cur = div(cur, nxt, W=28); ctru = [Fr(0)] * M
            else:
                cur = div(cur, nxt, W=28)
                ac = rmul(ctru, sconjm(nt)); ctru = [c / Nf for c in ac]
    viol += sound(cur, ctru)
report.append(("I. 混在演算木（mul/add/div ×4段, 300木）", viol))

# ---------------------------------------------- J. 分配律 (a+b)*c = a*c + b*c
viol = 0
for trial in range(300):
    a = [int(v) for v in rng.integers(-9, 9, M)]
    b = [int(v) for v in rng.integers(-9, 9, M)]
    c = [int(v) for v in rng.integers(-9, 9, M)]
    left = mul(add(BF(a), BF(b), W=40), BF(c), W=40)
    right = add(mul(BF(a), BF(c), W=40), mul(BF(b), BF(c), W=40), W=40)
    lv = [round(v) for v in left.values()]; rv = [round(v) for v in right.values()]
    ref = [int(x) for x in radd(rmul(a, c), rmul(b, c))]
    if lv != ref or rv != ref: viol += 1
report.append(("J. 分配律 (a+b)c = ac+bc（値一致すべき, 300）", viol))

# ---------------------------------------------- K. 非結合・非可換は ref と一致するか
viol = 0
for trial in range(200):
    a = [int(v) for v in rng.integers(-7, 7, M)]
    b = [int(v) for v in rng.integers(-7, 7, M)]
    c = [int(v) for v in rng.integers(-7, 7, M)]
    ab_c = mul(mul(BF(a), BF(b), W=48), BF(c), W=48)
    a_bc = mul(BF(a), mul(BF(b), BF(c), W=48), W=48)
    ref_l = [int(x) for x in rmul(rmul(a, b), c)]
    ref_r = [int(x) for x in rmul(a, rmul(b, c))]
    if [round(v) for v in ab_c.values()] != ref_l: viol += 1
    if [round(v) for v in a_bc.values()] != ref_r: viol += 1
    # 非結合ゆえ ref_l != ref_r が 普通（確認だけ）
report.append(("K. (ab)c と a(bc) が それぞれ ref と一致（非結合, 200）", viol))

# ---------------------------------------------- L. a−a=0, a+0=a, a*1=a, a/1=a, a/a=1
viol = 0
one = BF([1] + [0] * 15); zero = BF([0] * M)
for trial in range(300):
    a = [int(v) for v in rng.integers(-30, 30, M)]; A = BF(a, E=int(rng.integers(-4, 4)))
    if [round(v) for v in add(A, mul(A, BF([-1] + [0] * 15), W=40), W=40).values()] != [0] * M: viol += 1  # a+(−1)a=0
    if [round(v) for v in mul(A, one, W=40).values()] != a and A.E == 0: viol += 1                          # a*1=a
    if [round(v) for v in div(A, one, W=40).values()] != a and A.E == 0: viol += 1                          # a/1=a
    if sum(x * x for x in a) != 0:
        q = div(A, A, W=48)
        if [round(v) for v in q.values()] != [1] + [0] * 15: viol += 1                                      # a/a=1
report.append(("L. 恒等式 a−a=0, a*1=a, a/1=a, a/a=1（300）", viol))

# ---------------------------------------------- M. 極端な桁範囲（1e−19 と 1e30 を同居）
viol = 0
for trial in range(200):
    a = [0] * M
    a[0] = int(rng.integers(1, 100))                       # ~1
    a[1] = int(rng.integers(1, 100)) * 10 ** 9             # 巨大
    a[2] = 1                                               # 微小（相対的に）
    A = BF(a, E=int(rng.integers(-40, 40)))
    b = [int(v) for v in rng.integers(-5, 5, M)]
    c = mul(A, BF(b), W=16)
    viol += sound(c, rmul(tvals(A), [Fr(x) for x in b]))
report.append(("M. 極端な桁範囲（1〜1e11 同居, W=16, 200）", viol))

# ---------------------------------------------- N. 深い混在連鎖（50段）で 状態が 崩れないか
viol = 0
for trial in range(40):
    cur = BF([int(v) for v in rng.integers(-5, 5, M)], E=0); ctru = tvals(cur)
    for step in range(50):
        y = BF([int(v) for v in rng.integers(-4, 4, M)]); ny = tvals(y)
        if step % 3 == 2 and sum(t * t for t in ny) != 0:
            Nf = sum(t * t for t in ny)
            cur = div(cur, y, W=20, Emax=400); ctru = [c / Nf for c in rmul(ctru, sconjm(ny))]
        elif step % 3 == 1:
            cur = add(cur, y, W=20, Emax=400); ctru = radd(ctru, ny)
        else:
            cur = mul(cur, y, W=20, Emax=400); ctru = rmul(ctru, ny)
    viol += sound(cur, ctru)
report.append(("N. 深い混在連鎖 mul/add/div ×50段（40本）", viol))

# ---------------------------------------------- O. ノルム乗法性の破れ（零因子の証拠・法則でない）
# N(ab) == N(a)N(b) か（八元数までは真、セデニオンで破れる）
nfail = 0
for trial in range(200):
    a = [int(v) for v in rng.integers(-5, 5, M)]; b = [int(v) for v in rng.integers(-5, 5, M)]
    Nab = sum(x * x for x in [int(v) for v in rmul(a, b)])
    NaNb = sum(x * x for x in a) * sum(x * x for x in b)
    if Nab != NaNb: nfail += 1
report.append(("O. ノルム乗法性 N(ab)=N(a)N(b) が 破れる件数（法則でなく 代数の性質）", nfail))

# ---------------------------------------------- 結果
print("=" * 72)
print("さらに色んな計算 — 健全性違反 / 法則違反")
print("=" * 72)
for tag, v in report:
    if tag.startswith("O"):
        print(f"  {tag}: **{v}/200**（破れて 当然。セデニオンは 合成代数でない）")
    elif tag.startswith("K"):
        print(f"  {tag}: {v}（0 なら ref と一致・非結合を 正しく再現）")
    else:
        print(f"  {tag}: 違反 **{v}**")
print("=" * 72)
tot = sum(v for tag, v in report if not tag.startswith("O"))
print(f"\n代数法則・健全性の 違反 合計（O を除く）: **{tot}**")
