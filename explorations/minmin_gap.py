#!/usr/bin/env python3
"""**MIN − MIN も不定である。** 監査人が MAX しか見ておらず数え落としていた分を数える。"""
import warnings, numpy as np
from fractions import Fraction as Fr
from collections import Counter
warnings.filterwarnings("ignore")
from total_arith import MAX, MIN, saturate
OPS = "+-*/"
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

BIG  = lambda v: abs(v) >= MAX*0.999          # 溢れて飽和した
TINY = lambda v: v != 0.0 and abs(v) <= MIN*1.001   # **潰れて飽和した** ← 見ていなかった

def ev(t, flag_min):
    """(値, 不定フラグ)。flag_min=True なら MIN 同士の加減にもフラグを立てる。"""
    if t[0]=="leaf": return (float(t[1]), False)
    (a,pa),(b,pb) = ev(t[1],flag_min), ev(t[2],flag_min)
    with np.errstate(all="ignore"):
        if t[0]=="/":
            if b==0.0: v,zo = 0.0, True
            else: v,zo = np.float64(a)/np.float64(b), (a==0.0)
        elif t[0]=="*": v,zo = np.float64(a)*np.float64(b), (a==0.0 or b==0.0)
        else:
            v = np.float64(a)+np.float64(b) if t[0]=="+" else np.float64(a)-np.float64(b)
            zo = (float(v)==0.0 and abs(a)==abs(b))
        r = float(saturate(v, zero_ok=zo))
    ind = (BIG(a) and BIG(b) and t[0] in "+-")
    if flag_min: ind = ind or (TINY(a) and TINY(b) and t[0] in "+-")
    return (r, pa or pb or ind)

sg = lambda x: 0 if x==0 else (1 if x>0 else -1)
BLAME = Counter(); N=0; res={False: np.zeros(3), True: np.zeros(3)}
for seed in range(8):
    tree = make(seed); n=0; tries=0
    while n < 2000 and tries < 100000:
        tries+=1; t=tree(5); ex=ev_exact(t)
        if ex is None: continue
        n+=1; N+=1; s=sg(ex)
        for fm in (False, True):
            r,po = ev(t,fm)
            res[fm] += [(sg(r)==s) and not po, po, (sg(r)!=s) and not po]
        # 落とした符号の原因を、MIN も見て数え直す
        def blame(t):
            if t[0]=="leaf": return float(t[1])
            a,b = blame(t[1]), blame(t[2])
            with np.errstate(all="ignore"):
                if t[0]=="/":
                    v,zo = (0.0,True) if b==0.0 else (np.float64(a)/np.float64(b), a==0.0)
                elif t[0]=="*": v,zo = np.float64(a)*np.float64(b), (a==0.0 or b==0.0)
                else:
                    v = np.float64(a)+np.float64(b) if t[0]=="+" else np.float64(a)-np.float64(b)
                    zo = (float(v)==0.0 and abs(a)==abs(b))
                r = float(saturate(v, zero_ok=zo))
            ex2 = ev_exact(t)
            if ex2 is not None and ex2 != 0 and r == 0.0:
                if BIG(a) and BIG(b) and t[0] in "+-":   BLAME["**MAX 同士の加減**"] += 1
                elif TINY(a) and TINY(b) and t[0] in "+-": BLAME["**MIN 同士の加減** ← 数え落としていた"] += 1
                elif t[0]=="/" and b==0.0:                BLAME["本物の 0 で割った"] += 1
                else:                                     BLAME["その他"] += 1
            return r
        if sg(ev(t,True)[0]) != s: blame(t)

print("="*88); print(f"符号を落とした原因を、**MIN 側も数えて**（8 seed・{N:,} 式・深さ5）"); print("="*88)
for k,v in BLAME.most_common():
    print(f"  {v:>6,} 回   {k}")
print(); print("="*88); print("**MIN 同士にもフラグを立てると**"); print("="*88)
print(f"  {'':<30}{'符号が正しい':>14}{'正直に不明':>14}{'**黙って間違えた**':>20}")
for fm, lab in ((False,"MAX だけにフラグ（前回）"), (True,"**MAX と MIN 両方**")):
    c,u,w = res[fm]/N*100
    print(f"  {lab:<30}{f'{c:.2f}%':>14}{f'{u:.2f}%':>14}{f'**{w:.2f}%**':>20}")
