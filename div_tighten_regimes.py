#!/usr/bin/env python3
"""締める余地を **実用域** で 測る:
  R1. 近スカラー除数（支配的 b_0 + 小さな LE 摂動）= 大きさが 決まる 実用ケース
  R2. 有界除数のみ（LE/EQ・非零成分）= affine が 原理的に 効きうる 唯一の域
     （GE/NB は 非有界 ⟹ affine の 雑音記号で 表現不可 ⟹ そもそも 対象外）
各域で 境界なしのうち「見かけ（締められる）」割合 = affine の 上限。
"""
import sys, os
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
import numpy as np
from bfp_sed import BF, div, EQ, GE, LE, NB, M
from sedenion_tensor_logic import ref_mult

rng = np.random.default_rng(20260723)

def true_samples(m, flag, sunk):
    if flag == EQ:   base = [Fr(m)]
    elif flag == GE: base = [Fr(m)*r for r in (1,3,10)] if m != 0 else [Fr(0)]
    elif flag == LE: base = [Fr(m), Fr(m,2), Fr(m,5), Fr(0)] if m != 0 else [Fr(0)]
    else:            base = [Fr(m)*r for r in (1,3,10)] + [Fr(0)]
    if sunk: base = base + [-x for x in base]
    return base

def true_range(am, bm, abound, bbound, asunk, bsunk, n_mc):
    lo=[None]*M; hi=[None]*M
    ar=[true_samples(am[i],abound[i],asunk[i]) for i in range(M)]
    br=[true_samples(bm[i],bbound[i],bsunk[i]) for i in range(M)]
    for _ in range(n_mc):
        at=[ar[i][int(rng.integers(0,len(ar[i])))] for i in range(M)]
        bt=[br[i][int(rng.integers(0,len(br[i])))] for i in range(M)]
        Nb=sum(x*x for x in bt)
        q=[Fr(0)]*M if Nb==0 else [c/Nb for c in ref_mult(at,[bt[0]]+[-x for x in bt[1:]])]
        for k in range(M):
            if lo[k] is None or q[k]<lo[k]: lo[k]=q[k]
            if hi[k] is None or q[k]>hi[k]: hi[k]=q[k]
    return lo,hi

def classify(gen_ab, n_trials, tag, n_mc=120):
    A=B=C=D=tot=0
    for _ in range(n_trials):
        am,bm,abd,bbd,asu,bsu = gen_ab()
        A_=BF(am,bound=abd,sunk=asu); B_=BF(bm,bound=bbd,sunk=bsu)
        try: R=div(A_,B_,W=18)
        except Exception: continue
        lo,hi=true_range(am,bm,abd,bbd,asu,bsu,n_mc)
        for k in range(M):
            tot+=1
            if R.bound[k]==NB:
                if lo[k] is not None and lo[k]<0<hi[k]: A+=1
                else: B+=1
            elif R.bound[k] in (GE,LE): C+=1
            else: D+=1
    nb=A+B
    print(f"[{tag}]  成分{tot}  境界なし本物 {100*A/tot:4.1f}%  見かけ {100*B/tot:4.1f}%  量的 {100*C/tot:4.1f}%  厳密 {100*D/tot:4.1f}%")
    if nb: print(f"          ⟹ 境界なしのうち 見かけ = **{100*B/nb:.1f}%**（affine の上限）")
    return B, nb

def gen_near_scalar():
    """b_0 が 支配的（大）+ 1..M に 小さな LE 摂動。a も 支配的 + LE。"""
    bm=[0]*M; bm[0]=int(rng.choice([-100,100]))
    for i in range(1,M):
        if rng.random()<0.4: bm[i]=int(rng.integers(-5,6))
    am=[0]*M; am[0]=int(rng.choice([-100,100]))
    for i in range(1,M):
        if rng.random()<0.4: am[i]=int(rng.integers(-5,6))
    abd=[EQ]*M; bbd=[EQ]*M
    abd[0]=int(rng.choice([EQ,GE,LE]))
    for i in range(1,M):
        if bm[i]!=0: bbd[i]=LE           # 小摂動は |·|≤ の 有界
        if am[i]!=0 and rng.random()<0.5: abd[i]=LE
    return am,bm,abd,bbd,[False]*M,[False]*M

def gen_bounded_only():
    """除数の 全非零成分が LE/EQ（有界）= affine が 表現できる 域。"""
    bm=[int(v) for v in rng.integers(-6,6,M)]
    while all(bm[i]==0 for i in range(1,M)): bm=[int(v) for v in rng.integers(-6,6,M)]
    am=[int(v) for v in rng.integers(-6,6,M)]
    bbd=[LE if bm[i]!=0 else EQ for i in range(M)]     # 除数 全有界
    abd=[int(rng.choice([EQ,LE])) for _ in range(M)]   # 被除数も 有界
    return am,bm,abd,bbd,[False]*M,[False]*M

if __name__ == "__main__":
    print("="*80)
    classify(gen_near_scalar, 800, "R1 近スカラー除数（実用・大きさが決まる域）")
    print("-"*80)
    classify(gen_bounded_only, 800, "R2 有界除数のみ（affine が 効きうる 唯一の域）")
    print("="*80)
