#!/usr/bin/env python3
"""**符号を3値にする**（+ / − / **?**）。利用者の提案。監査人の飽和ビットとの比較。

飽和ビットは「この式のどこかで情報が失われた」しか言えないので、**言い過ぎる**
（`MAX + MAX` は両方正なら答えは確実に正なのに鳴る）。
3値符号は「**符号が分からない**」だけを言うので、**大きさの比較が要るときにしか ? にならない**。
"""
import warnings, numpy as np
from fractions import Fraction as Fr
warnings.filterwarnings("ignore")
from total_arith import MAX, MIN, saturate
OPS = "+-*/"
BIG  = lambda v: abs(v) >= MAX*0.999
TINY = lambda v: v != 0.0 and abs(v) <= MIN*1.001
SATD = lambda v: BIG(v) or TINY(v)            # **飽和した = 真の大きさが分からない**

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
def ev_ieee(t):
    if t[0]=="leaf": return float(t[1])
    a,b = np.float64(ev_ieee(t[1])), np.float64(ev_ieee(t[2]))
    with np.errstate(all="ignore"): return float({"+":a+b,"-":a-b,"*":a*b,"/":a/b}[t[0]])

def raw(v, zo):
    return float(saturate(v, zero_ok=zo))

def ev_flag(t):
    """**監査人の版**: 飽和ビット1個。両端が飽和した加減で鳴る（符号が同じでも鳴る）。"""
    if t[0]=="leaf": return (float(t[1]), False)
    (a,pa),(b,pb) = ev_flag(t[1]), ev_flag(t[2])
    with np.errstate(all="ignore"):
        if t[0]=="/": v,zo = (0.0,True) if b==0.0 else (np.float64(a)/np.float64(b), a==0.0)
        elif t[0]=="*": v,zo = np.float64(a)*np.float64(b), (a==0.0 or b==0.0)
        else:
            v = np.float64(a)+np.float64(b) if t[0]=="+" else np.float64(a)-np.float64(b)
            zo = (float(v)==0.0 and abs(a)==abs(b))
        r = raw(v, zo)
    ind = pa or pb or (SATD(a) and SATD(b) and t[0] in "+-")
    return (r, ind)

def ev_tri(t):
    """**利用者の版**: 値と、3値の符号 s ∈ {+1, −1, 0（本物の0）, None（**不明**）}。"""
    if t[0]=="leaf":
        v = float(t[1]); return (v, 0 if v==0 else (1 if v>0 else -1))
    (a,sa),(b,sb) = ev_tri(t[1]), ev_tri(t[2])
    with np.errstate(all="ignore"):
        if t[0] in "*/":
            if t[0]=="/":
                if b==0.0: v,zo = 0.0, True
                else: v,zo = np.float64(a)/np.float64(b), (a==0.0)
            else: v,zo = np.float64(a)*np.float64(b), (a==0.0 or b==0.0)
            r = raw(v, zo)
            # **積・商の符号は「符号の積」で決まる。大きさが分からなくても符号は残る**
            if t[0]=="/" and b==0.0: s = 0                       # a/0 = 0（本物の0）
            elif sa == 0 or sb == 0: s = 0
            elif sa is None or sb is None: s = None              # ? は吸収的
            else: s = sa*sb
            return (r, s)
        # ---- 加減 ----
        bb, sbb = (b, sb) if t[0]=="+" else (-b, (None if sb is None else -sb))
        v = np.float64(a)+np.float64(bb)
        r = raw(v, (float(v)==0.0 and abs(a)==abs(bb)))
        if sa is None or sbb is None: s = None
        elif sa == 0: s = sbb
        elif sbb == 0: s = sa
        elif sa == sbb: s = sa                                   # **同符号なら、大きさ不明でも符号は確実**
        else:
            # **異符号 ⟹ 大きさの比較が要る。両方とも真の大きさが不明なら ?**
            s = None if (SATD(a) and SATD(bb)) else (0 if r==0 else (1 if r>0 else -1))
        return (r, s)

sg = lambda x: 0 if x==0 else (1 if x>0 else -1)
N=0; acc = {"ieee": np.zeros(3), "flag": np.zeros(3), "tri": np.zeros(3)}
for seed in range(8):
    tree = make(seed); n=0; tries=0
    while n < 2000 and tries < 100000:
        tries+=1; t=tree(5); ex=ev_exact(t)
        if ex is None: continue
        n+=1; N+=1; s=sg(ex)
        vi = ev_ieee(t); isn = bool(np.isnan(vi))
        acc["ieee"] += [0 if isn else (sg(vi)==s), isn, 0 if isn else (sg(vi)!=s)]
        r,po = ev_flag(t)
        acc["flag"] += [(sg(r)==s) and not po, po, (sg(r)!=s) and not po]
        r2,s2 = ev_tri(t)
        unk = (s2 is None)
        acc["tri"] += [(not unk) and (s2==s), unk, (not unk) and (s2!=s)]

print("="*94)
print(f"**符号を3値にする**（+ / − / ?）  — 8 seed・{N:,} 式・深さ5・真値は正確な有理数")
print("="*94)
print(f"  {'':<34}{'符号が正しい':>14}{'正直に不明と言った':>18}{'**黙って間違えた**':>20}")
for k,lab in (("ieee","IEEE 754（NaN が警報）"), ("flag","飽和ビット1個（監査人）"), ("tri","**3値符号（利用者）**")):
    c,u,w = acc[k]/N*100
    print(f"  {lab:<34}{f'{c:.2f}%':>14}{f'{u:.2f}%':>18}{f'**{w:.2f}%**':>20}")
fc,fu,fw = acc["flag"]/N*100; tc,tu,tw = acc["tri"]/N*100
print(f"""
  ⟹ 3値符号は、飽和ビットより **警報が {fu/max(tu,1e-9):.1f} 倍少なく**、
     **符号を当てる率が {tc-fc:+.2f} ポイント**、黙って間違える率が **{fw:.2f}% → {tw:.2f}%**。""")
