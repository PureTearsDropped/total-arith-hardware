#!/usr/bin/env python3
"""**2ビット = 独立した2つの事実**（利用者の修正）。型エラーは値の欄に入れない。

    00  確定       符号○ 大きさ○
    01  大きさ不明  符号○ 大きさ×      MAX / MIN / a/0=0 の 0
    10  **符号不明** 符号× **大きさ○**  **MIN − MIN**（真値は |r| ≤ 2·MIN ＝ **無視できると分かっている**）
    11  **両方不明** 符号× 大きさ×      **MAX − MAX**（真値は [0, ∞) のどこか）

`10` が空席でないのが要点。**符号が分からなくても、大きさが「無視できる」と分かっていれば値は使える。**
`(MIN − MIN) + 5.0` の符号は **+ で確定**。監査人の 3値/4値 はここを `?` にして捨てていた。
"""
import warnings, numpy as np
from fractions import Fraction as Fr
warnings.filterwarnings("ignore")
from total_arith import MAX, MIN, saturate
OPS = "+-*/"
BIG  = lambda v: abs(v) >= MAX*0.999
TINY = lambda v: v != 0.0 and abs(v) <= MIN*1.001
SATD = lambda v: BIG(v) or TINY(v)
NEG  = 4*MIN          # 「無視できる」の上限（MIN−MIN の真値は |r| ≤ 2·MIN）

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

def ev(t, use10):
    """(値, 符号 s∈{+1,−1,0,None}, 大きさ不明 m, **大きさの上界 ub**（None=不明）)。
       use10=False なら `10` の席を使わない（＝前回の4値。MIN−MIN も 11 扱い）。"""
    if t[0]=="leaf":
        v=float(t[1]); return (v, 0 if v==0 else (1 if v>0 else -1), False, abs(v))
    (a,sa,ma,ua),(b,sb,mb,ub) = ev(t[1],use10), ev(t[2],use10)
    with np.errstate(all="ignore"):
        if t[0] in "*/":
            if t[0]=="/":
                if b==0.0: return (0.0, sa if sa else 0, True, None)
                v,zo = np.float64(a)/np.float64(b), (a==0.0)
            else: v,zo = np.float64(a)*np.float64(b), (a==0.0 or b==0.0)
            r = raw(v,zo)
            s = 0 if sa==0 or sb==0 else (None if sa is None or sb is None else sa*sb)
            m = ma or mb or SATD(r) or SATD(a) or SATD(b)
            # 上界の伝播: 積は上界の積。商は分母の下界が要るので諦める
            u = (ua*ub if (t[0]=="*" and ua is not None and ub is not None and ua*ub < np.inf) else None)
            return (r, s, m, u)
        bb,sbb,mbb,ubb = (b,sb,mb,ub) if t[0]=="+" else (-b,(None if sb is None else -sb),mb,ub)
        v = np.float64(a)+np.float64(bb); r = raw(v,(float(v)==0.0 and abs(a)==abs(bb)))
        m  = ma or mbb or SATD(r) or SATD(a) or SATD(bb)
        u  = (ua+ubb if (ua is not None and ubb is not None and ua+ubb < np.inf) else None)
        # ---- **`10` の席: 符号は不明だが、大きさが無視できると分かっている相手** ----
        if use10:
            a_neg   = (sa  is None) and (ua  is not None) and (ua  <= NEG)
            bb_neg  = (sbb is None) and (ubb is not None) and (ubb <= NEG)
            if a_neg and (sbb is not None) and (ubb is None or ubb > NEG):
                return (r, sbb, m, u)          # **無視できる方は符号を持たなくてよい**
            if bb_neg and (sa is not None) and (ua is None or ua > NEG):
                return (r, sa, m, u)
        if sa is None or sbb is None: return (r, None, m, u)
        if sa==0: return (r, sbb, m, u)
        if sbb==0: return (r, sa, m, u)
        if sa==sbb: return (r, sa, m, u)
        if BIG(a) and not BIG(bb):   return (r, sa,  True, None)
        if BIG(bb) and not BIG(a):   return (r, sbb, True, None)
        if TINY(a) and not TINY(bb): return (r, sbb, m, u)
        if TINY(bb) and not TINY(a): return (r, sa,  m, u)
        if BIG(a) and BIG(bb):       return (r, None, True, None)    # **11: 両方不明**
        if TINY(a) and TINY(bb):     return (r, None, m, 2*MIN)      # **10: 符号だけ不明。大きさ ≤ 2·MIN**
        return (r, (0 if r==0 else (1 if r>0 else -1)), m, u)

sg = lambda x: 0 if x==0 else (1 if x>0 else -1)
N=0; acc={False: np.zeros(3), True: np.zeros(3)}
for seed in range(8):
    tree=make(seed); n=0; tries=0
    while n<2000 and tries<100000:
        tries+=1; t=tree(5); ex=ev_exact(t)
        if ex is None: continue
        n+=1; N+=1; s=sg(ex)
        for u10 in (False, True):
            r,ss,mm,uu = ev(t,u10); unk=(ss is None)
            acc[u10] += [(not unk) and ss==s, unk, (not unk) and ss!=s]
print("="*94)
print(f"**`10`（符号不明・大きさは無視できると確定）の席を使うか**  — 8 seed・{N:,} 式・深さ5")
print("="*94)
print(f"  {'':<40}{'符号が正しい':>14}{'符号不明と言った':>18}{'**黙って間違えた**':>20}")
for u10,lab in ((False,"3値（MIN−MIN も「分からない」）"), (True,"**4値（MIN−MIN は `10`）**")):
    c,u,w = acc[u10]/N*100
    print(f"  {lab:<40}{f'{c:.2f}%':>14}{f'{u:.2f}%':>18}{f'**{w:.2f}%**':>20}")
d = (acc[True]-acc[False])/N*100
print(f"\n  差: 符号が正 **{d[0]:+.3f}pt** ／ 不明 **{d[1]:+.3f}pt** ／ 黙って間違え **{d[2]:+.3f}pt**")
