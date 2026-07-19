#!/usr/bin/env python3
"""「MIN と MAX の演算をどうするか」の 3 候補を測り、IEEE の NaN が正当かも見る。"""
import warnings, numpy as np
from fractions import Fraction as Fr
warnings.filterwarnings("ignore")
from total_arith import MAX, MIN, saturate
rng = np.random.default_rng(11); OPS = "+-*/"
def leaf(): return (-1 if rng.random()<.5 else 1)*int(rng.integers(1,999))*Fr(10)**int(rng.integers(-180,181))
def tree(d): return ("leaf", leaf()) if d==0 else (OPS[rng.integers(4)], tree(d-1), tree(d-1))
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

def ev(t, policy):
    """policy: '0'=そのまま / 'absorb'=MAX−MAX→MAX / 'sticky'=飽和ビットを伝播（(値, 汚染)）"""
    if t[0]=="leaf": return (float(t[1]), False)
    (a,pa), (b,pb) = ev(t[1],policy), ev(t[2],policy)
    with np.errstate(all="ignore"):
        if t[0]=="/":
            if b==0.0: v,zo = 0.0, True
            else: v,zo = np.float64(a)/np.float64(b), (a==0.0)
        elif t[0]=="*": v,zo = np.float64(a)*np.float64(b), (a==0.0 or b==0.0)
        else:
            v = np.float64(a)+np.float64(b) if t[0]=="+" else np.float64(a)-np.float64(b)
            zo = (float(v)==0.0 and abs(a)==abs(b))
        r = float(saturate(v, zero_ok=zo))
    both_sat = SAT(a) and SAT(b) and t[0] in "+-"
    poison = pa or pb or both_sat        # **狭い汚染: MAX 同士の加減 だけ**。溢れただけでは立てない
    if policy=="absorb" and both_sat and r==0.0: r = MAX if a>0 else -MAX   # ① 吸収的
    return (r, poison)

sg = lambda x: 0 if x==0 else (1 if x>0 else -1)
n=0; res={p:0 for p in ("0","absorb","sticky")}; flag=0; nan=0; nan_just=0; genuine=0; tries=0
while n<6000 and tries<300000:
    tries+=1; t=tree(5); ex=ev_exact(t)
    if ex is None: continue
    n+=1; s=sg(ex)
    r0,_    = ev(t,"0")
    ra,_    = ev(t,"absorb")
    rs,po   = ev(t,"sticky")
    res["0"]      += (sg(r0)==s)
    res["absorb"] += (sg(ra)==s)
    res["sticky"] += (sg(rs)==s or po)          # 汚染フラグが立てば「間違えた」ではなく「分からないと言った」
    flag += po
    vi = ev_ieee(t); isn = bool(np.isnan(vi)); nan += isn
    bad0 = (sg(r0)!=s)
    genuine += bad0
    if isn and bad0: nan_just += 1

print("="*88); print(f"深さ5・{n:,} 式。「MIN/MAX の演算をどうするか」3 候補"); print("="*88)
print(f"  {'候補':<40}{'符号が正（か、正直に不明と言った）':>32}")
print(f"  {'③ MAX−MAX = 0（今の実装）':<40}{f'**{100*res[chr(48)]/n:.2f}%**':>32}")
print(f"  {'① MAX−MAX = MAX（吸収的）':<40}{f'**{100*res[chr(97)+chr(98)+chr(115)+chr(111)+chr(114)+chr(98)]/n:.2f}%**':>32}")
print(f"  {'② 飽和ビットを伝播（NaN の軽量版）':<40}{f'**{100*res[chr(115)+chr(116)+chr(105)+chr(99)+chr(107)+chr(121)]/n:.2f}%**':>32}")
print(f"\n  ② が立てたフラグ: **{100*flag/n:.2f}%** の式")
print("="*88); print("IEEE の NaN は、正当に鳴っているか"); print("="*88)
print(f"  NaN が出た式            : **{100*nan/n:.2f}%**")
print(f"  本当に符号が復元不能な式 : **{100*genuine/n:.2f}%**")
print(f"  NaN のうち、正当だった割合: **{100*nan_just/max(nan,1):.1f}%**"
      f"  ⟹ 残り **{100-100*nan_just/max(nan,1):.1f}%** は **鳴らなくてよかった**")
