#!/usr/bin/env python3
"""**同じ seed で横並び。** 「正しい符号を出す」か「正直に分からないと言う」か、それ以外は負け。"""
import warnings, numpy as np
from fractions import Fraction as Fr
warnings.filterwarnings("ignore")
from total_arith import MAX, saturate
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
def ev_ieee(t):
    if t[0]=="leaf": return float(t[1])
    a,b = np.float64(ev_ieee(t[1])), np.float64(ev_ieee(t[2]))
    with np.errstate(all="ignore"): return float({"+":a+b,"-":a-b,"*":a*b,"/":a/b}[t[0]])
SAT = lambda v: abs(v) >= MAX*0.999
def ev_sticky(t):
    if t[0]=="leaf": return (float(t[1]), False)
    (a,pa),(b,pb) = ev_sticky(t[1]), ev_sticky(t[2])
    with np.errstate(all="ignore"):
        if t[0]=="/":
            if b==0.0: v,zo = 0.0, True
            else: v,zo = np.float64(a)/np.float64(b), (a==0.0)
        elif t[0]=="*": v,zo = np.float64(a)*np.float64(b), (a==0.0 or b==0.0)
        else:
            v = np.float64(a)+np.float64(b) if t[0]=="+" else np.float64(a)-np.float64(b)
            zo = (float(v)==0.0 and abs(a)==abs(b))
        r = float(saturate(v, zero_ok=zo))
    return (r, pa or pb or (SAT(a) and SAT(b) and t[0] in "+-"))
sg = lambda x: 0 if x==0 else (1 if x>0 else -1)

print("="*94)
print("**演算回路としての判定** — 深さ5、8 seed、正確な有理数を真値に")
print("="*94)
print(f"  {'':<34}{'符号が正しい':>14}{'正直に不明と言った':>18}{'**黙って間違えた**':>20}{'警報の頻度':>12}")
tot = {k: np.zeros(3) for k in ("ieee","total")}; N = 0
for seed in range(8):
    tree = make(seed); n=0; tries=0
    while n < 2000 and tries < 100000:
        tries += 1; t = tree(5); ex = ev_exact(t)
        if ex is None: continue
        n += 1; s = sg(ex)
        vi = ev_ieee(t); isn = bool(np.isnan(vi))
        tot["ieee"] += [0 if isn else (sg(vi)==s), isn, 0 if isn else (sg(vi)!=s)]
        r, po = ev_sticky(t)
        tot["total"] += [(sg(r)==s) and not po, po, (sg(r)!=s) and not po]
    N += n
for k, lab in (("ieee","IEEE 754（NaN が警報）"), ("total","**この規約 + 飽和ビット**")):
    c, u, w = tot[k]/N*100
    print(f"  {lab:<34}{f'{c:.2f}%':>14}{f'{u:.2f}%':>18}{f'**{w:.2f}%**':>20}{f'{u:.2f}%':>12}")
print(f"""
  ⟹ **黙って間違える率が {tot['ieee'][2]/max(tot['total'][2],1):.1f} 倍ちがう。** しかも警報は **{tot['ieee'][1]/max(tot['total'][1],1):.1f} 倍少ない。**
     NaN は「分からない」と言うが、**言い過ぎ**（75% は鳴らなくてよかった）。
     この規約は **鳴るべき所（MAX 同士の打ち消し）でだけ鳴る**。""")
