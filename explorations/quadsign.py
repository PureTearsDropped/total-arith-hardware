#!/usr/bin/env python3
"""**2ビット・4値**（利用者の提案）。

    00  **確定**              符号も大きさも本物
    01  **大きさ不明**        符号は確定。値は飽和した作り物（MAX/MIN、および a/0=0 の 0）
    10  **符号不明**          加減で大きさの比較に失敗した
    11  **型が違う**          √(−1)。実数の欄に置けない ⟹ ℂ へ行け（§7.6.6 の「NaN の兼務」の片割れ）

要点: **MAX と MIN は向きが逆**。MAX = 「真の大きさは **MAX 以上**」、MIN = 「**MIN 以下**」。
だから `MAX + (−MIN)` は **確実に正**。3値実装は「両方飽和 → ?」で諦めていた。
"""
import warnings, numpy as np
from fractions import Fraction as Fr
warnings.filterwarnings("ignore")
from total_arith import MAX, MIN, saturate
OPS = "+-*/"
BIG  = lambda v: abs(v) >= MAX*0.999          # 真の大きさは **これ以上**
TINY = lambda v: v != 0.0 and abs(v) <= MIN*1.001   # 真の大きさは **これ以下**
SATD = lambda v: BIG(v) or TINY(v)

def make(seed):
    rng = np.random.default_rng(seed)
    def leaf(): return (-1 if rng.random()<.5 else 1)*int(rng.integers(1,999))*Fr(10)**int(rng.integers(-180,181))
    def tree(d): return ("leaf", leaf()) if d==0 else (OPS[rng.integers(4)], tree(d-1), tree(d-1))
    return tree
def ev_exact(t):
    if t[0]=="leaf": return t[1]
    a,b = ev_exact(t[1]), ev_exact(t[2])
    if a is None or b is None or (t[0]=="/" and b==0): return None
    return {"+":a+b,"-":a-b,"*":a*b,"/":(a/b if b else None)}[t[0]]
raw = lambda v, zo: float(saturate(v, zero_ok=zo))

def ev3(t):
    """3値（前回）: 両方飽和した異符号の加減 → ?"""
    if t[0]=="leaf":
        v=float(t[1]); return (v, 0 if v==0 else (1 if v>0 else -1))
    (a,sa),(b,sb) = ev3(t[1]), ev3(t[2])
    with np.errstate(all="ignore"):
        if t[0] in "*/":
            if t[0]=="/": v,zo = (0.0,True) if b==0.0 else (np.float64(a)/np.float64(b), a==0.0)
            else: v,zo = np.float64(a)*np.float64(b), (a==0.0 or b==0.0)
            r = raw(v,zo)
            s = 0 if (t[0]=="/" and b==0.0) or sa==0 or sb==0 else (None if sa is None or sb is None else sa*sb)
            return (r,s)
        bb,sbb = (b,sb) if t[0]=="+" else (-b,(None if sb is None else -sb))
        v = np.float64(a)+np.float64(bb); r = raw(v,(float(v)==0.0 and abs(a)==abs(bb)))
        if sa is None or sbb is None: s=None
        elif sa==0: s=sbb
        elif sbb==0: s=sa
        elif sa==sbb: s=sa
        else: s = None if (SATD(a) and SATD(bb)) else (0 if r==0 else (1 if r>0 else -1))
        return (r,s)

def ev4(t):
    """**4値**: (値, 符号, 大きさ不明か)。**MAX（以上）と MIN（以下）の向きを使う。**"""
    if t[0]=="leaf":
        v=float(t[1]); return (v, 0 if v==0 else (1 if v>0 else -1), False)
    (a,sa,ma),(b,sb,mb) = ev4(t[1]), ev4(t[2])
    with np.errstate(all="ignore"):
        if t[0] in "*/":
            if t[0]=="/":
                if b==0.0: return (0.0, sa if sa else 0, True)     # **a/0=0 は「大きさ不明」の作り物**
                v,zo = np.float64(a)/np.float64(b), (a==0.0)
            else: v,zo = np.float64(a)*np.float64(b), (a==0.0 or b==0.0)
            r = raw(v,zo)
            s = 0 if sa==0 or sb==0 else (None if sa is None or sb is None else sa*sb)
            return (r, s, ma or mb or SATD(r))                     # 符号は積で残る。大きさだけ不明
        bb,sbb,mbb = (b,sb,mb) if t[0]=="+" else (-b,(None if sb is None else -sb),mb)
        v = np.float64(a)+np.float64(bb); r = raw(v,(float(v)==0.0 and abs(a)==abs(bb)))
        mu = ma or mbb or SATD(r)
        if sa is None or sbb is None: return (r, None, mu)
        if sa==0: return (r, sbb, mu)
        if sbb==0: return (r, sa, mu)
        if sa==sbb: return (r, sa, mu)                             # 同符号 ⟹ 符号は確実
        # ---- 異符号。**どちらが大きいか言えるか** ----
        if BIG(a) and not BIG(bb):  return (r, sa,  True)          # **MAX 以上 vs それ未満 ⟹ MAX 側が勝つ**
        if BIG(bb) and not BIG(a):  return (r, sbb, True)
        if TINY(a) and not TINY(bb):return (r, sbb, mu)            # **MIN 以下は無視できる ⟹ 相手の符号**
        if TINY(bb) and not TINY(a):return (r, sa,  mu)
        if SATD(a) and SATD(bb):    return (r, None, True)         # **同じ端で飽和 ⟹ 本当に分からない**
        return (r, (0 if r==0 else (1 if r>0 else -1)), mu)

sg = lambda x: 0 if x==0 else (1 if x>0 else -1)
N=0; acc={"3":np.zeros(3), "4":np.zeros(3)}; mflag=0
for seed in range(8):
    tree=make(seed); n=0; tries=0
    while n<2000 and tries<100000:
        tries+=1; t=tree(5); ex=ev_exact(t)
        if ex is None: continue
        n+=1; N+=1; s=sg(ex)
        r,s3 = ev3(t); u=(s3 is None)
        acc["3"] += [(not u) and s3==s, u, (not u) and s3!=s]
        r4,s4,m4 = ev4(t); u4=(s4 is None)
        acc["4"] += [(not u4) and s4==s, u4, (not u4) and s4!=s]
        mflag += m4
print("="*94)
print(f"**2ビット4値**（利用者）vs 3値  — 8 seed・{N:,} 式・深さ5")
print("="*94)
print(f"  {'':<30}{'符号が正しい':>14}{'符号不明と言った':>18}{'**黙って間違えた**':>20}")
for k,lab in (("3","3値（+ / − / ?）"), ("4","**4値（+ 大きさ不明）**")):
    c,u,w = acc[k]/N*100
    print(f"  {lab:<30}{f'{c:.2f}%':>14}{f'{u:.2f}%':>18}{f'**{w:.2f}%**':>20}")
print(f"\n  4値が「**大きさ不明**」を立てた式: **{100*mflag/N:.2f}%**"
      f"   ← **符号は使えるが、値は信じるなと言えている**")

# ============ 「大きさ不明」は、本当に大きさが違うときに立っているか ============
def ev4_tight(t):
    """**伝播を締める**: 一度立った「大きさ不明」を、消せるときは消す。
       `a/a = 1` は a が飽和していても **1 は正しい**（大きさは確定）。"""
    if t[0]=="leaf":
        v=float(t[1]); return (v, 0 if v==0 else (1 if v>0 else -1), False)
    (a,sa,ma),(b,sb,mb) = ev4_tight(t[1]), ev4_tight(t[2])
    with np.errstate(all="ignore"):
        if t[0] in "*/":
            if t[0]=="/":
                if b==0.0: return (0.0, sa if sa else 0, True)
                v,zo = np.float64(a)/np.float64(b), (a==0.0)
            else: v,zo = np.float64(a)*np.float64(b), (a==0.0 or b==0.0)
            r = raw(v,zo)
            s = 0 if sa==0 or sb==0 else (None if sa is None or sb is None else sa*sb)
            # **飽和した値同士の商・積は、飽和が打ち消し合えば大きさが戻る場合がある**が、
            # 一般には戻らない。**戻ると言えるのは、飽和が一度も起きていないときだけ。**
            return (r, s, ma or mb or SATD(r) or SATD(a) or SATD(b))
        bb,sbb,mbb = (b,sb,mb) if t[0]=="+" else (-b,(None if sb is None else -sb),mb)
        v = np.float64(a)+np.float64(bb); r = raw(v,(float(v)==0.0 and abs(a)==abs(bb)))
        mu = ma or mbb or SATD(r) or SATD(a) or SATD(bb)
        if sa is None or sbb is None: return (r, None, mu)
        if sa==0: return (r, sbb, mu)
        if sbb==0: return (r, sa, mu)
        if sa==sbb: return (r, sa, mu)
        if BIG(a) and not BIG(bb):  return (r, sa,  True)
        if BIG(bb) and not BIG(a):  return (r, sbb, True)
        if TINY(a) and not TINY(bb):return (r, sbb, mu)
        if TINY(bb) and not TINY(a):return (r, sa,  mu)
        if SATD(a) and SATD(bb):    return (r, None, True)
        return (r, (0 if r==0 else (1 if r>0 else -1)), mu)

print()
print("="*94)
print("「**大きさ不明**」の精度 — 立ったとき、本当に大きさが違うか（真値＝正確な有理数）")
print("="*94)
N2=0; tp=fp=tn=fn=0
for seed in range(8):
    tree=make(seed); n=0; tries=0
    while n<2000 and tries<100000:
        tries+=1; t=tree(5); ex=ev_exact(t)
        if ex is None: continue
        n+=1; N2+=1
        r,s,m = ev4(t)
        # **本当に大きさが違うか**: 相対誤差 > 1%（真値 0 なら値も 0 か）
        if ex == 0: wrong = (r != 0.0)
        else:
            try: wrong = abs(Fr(r).limit_denominator(10**30) - ex) > abs(ex)/100
            except (OverflowError, ValueError): wrong = True
        tp += (m and wrong); fp += (m and not wrong); fn += ((not m) and wrong); tn += ((not m) and not wrong)
prec = 100*tp/max(tp+fp,1); rec = 100*tp/max(tp+fn,1)
print(f"  立てた式: **{100*(tp+fp)/N2:.2f}%**   本当に大きさが違う式: **{100*(tp+fn)/N2:.2f}%**")
print(f"  **適合率（立ったとき当たっている）: {prec:.1f}%**")
print(f"  **再現率（違うとき立っている）    : {rec:.1f}%**")
print(f"  **黙って間違った大きさを返した   : {100*fn/N2:.2f}%**")
print(f"\n  ⟹ " + ("**言い過ぎではない。** 立った式の大半は、本当に大きさが違う" if prec > 80
      else f"**言い過ぎ。**立った式の {100-prec:.0f}% は大きさが合っている ⟹ **付箋として使えない**"))
