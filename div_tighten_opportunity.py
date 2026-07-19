#!/usr/bin/env python3
"""セデニオン除数の除算で「締める余地」を測る（affine を作る前の 上限測定）。

境界なしに なっている 商成分を 三分類:
  A. 真の範囲が **本当に 0 を跨ぐ** → 境界なしが 正しい・どんな手法でも 締められない
  B. 真の範囲が 0 を跨がない のに 境界なし → **依存性ゆえの 見かけ・締められる**（affine の的）
  C. 既に 量的境界（≥/≤）

B の割合が affine で 取れる 上限。B≈0 なら affine は 無益 → 正直に そう報告する。
"""
import sys, os
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
import numpy as np
from bfp_sed import BF, div, EQ, GE, LE, NB, M
from sedenion_tensor_logic import ref_mult

rng = np.random.default_rng(20260722)

def true_samples(m, flag, sunk, reps):
    """フラグ区間内の 真値候補（Fraction）。GE/NB は 有界近似（reps 倍まで）で 範囲推定。"""
    if flag == EQ:   base = [Fr(m)]
    elif flag == GE: base = [Fr(m)*r for r in reps] if m != 0 else [Fr(0)]
    elif flag == LE: base = [Fr(m), Fr(m,2), Fr(m,5), Fr(0)] if m != 0 else [Fr(0)]
    else:            base = [Fr(m)*r for r in reps] + [Fr(0)]      # NB
    if sunk: base = base + [-x for x in base]
    return base

def true_range_qk(am, bm, abound, bbound, asunk, bsunk, n_mc):
    """モンテカルロで 各商成分 q_k の 真の [min,max] を 推定（Fraction）。"""
    lo = [None]*M; hi = [None]*M
    ar = [true_samples(am[i], abound[i], asunk[i], (1,3,10)) for i in range(M)]
    br = [true_samples(bm[i], bbound[i], bsunk[i], (1,3,10)) for i in range(M)]
    for _ in range(n_mc):
        at = [ar[i][int(rng.integers(0,len(ar[i])))] for i in range(M)]
        bt = [br[i][int(rng.integers(0,len(br[i])))] for i in range(M)]
        Nb = sum(x*x for x in bt)
        if Nb == 0:
            q = [Fr(0)]*M
        else:
            ac = ref_mult(at, [bt[0]] + [-x for x in bt[1:]])
            q = [c/Nb for c in ac]
        for k in range(M):
            if lo[k] is None or q[k] < lo[k]: lo[k] = q[k]
            if hi[k] is None or q[k] > hi[k]: hi[k] = q[k]
    return lo, hi

def main():
    print("="*78)
    print("セデニオン除数の除算: 境界なし成分の 締める余地（affine の上限）")
    print("="*78)
    A_cross = 0   # 真に 0 跨ぎ（締められない）
    B_spur  = 0   # 見かけの 境界なし（締められる = affine の的）
    C_quant = 0   # 既に 量的境界
    D_exact = 0
    tot = 0
    n_trials = 800
    for _ in range(n_trials):
        # 除数を セデニオン（1..M に 非零）にして 一般分岐へ
        am = [int(v) for v in rng.integers(-6, 6, M)]
        bm = [int(v) for v in rng.integers(-6, 6, M)]
        while all(bm[i]==0 for i in range(1,M)):
            bm = [int(v) for v in rng.integers(-6, 6, M)]
        abound = [int(rng.choice([EQ, GE, LE])) for _ in range(M)]
        bbound = [int(rng.choice([EQ, GE, LE])) for _ in range(M)]
        asunk = [bool(rng.random()<0.1) for _ in range(M)]
        bsunk = [bool(rng.random()<0.1) for _ in range(M)]
        A = BF(am, bound=abound, sunk=asunk); B = BF(bm, bound=bbound, sunk=bsunk)
        try: R = div(A, B, W=18)
        except Exception: continue
        lo, hi = true_range_qk(am, bm, abound, bbound, asunk, bsunk, n_mc=120)
        for k in range(M):
            tot += 1
            if R.bound[k] == NB:
                if lo[k] is not None and lo[k] < 0 < hi[k]:
                    A_cross += 1                 # 真に跨ぐ
                else:
                    B_spur += 1                  # 見かけ（締められる）
            elif R.bound[k] in (GE, LE):
                C_quant += 1
            else:
                D_exact += 1
    print(f"  試行 {n_trials}・成分 {tot}")
    print(f"  A. 境界なし・真に0跨ぎ（締められない）   : {A_cross:6d}  {100*A_cross/tot:5.1f}%")
    print(f"  B. 境界なし・見かけ（**affine の的**）    : {B_spur:6d}  {100*B_spur/tot:5.1f}%")
    print(f"  C. 既に 量的境界（≥/≤）                   : {C_quant:6d}  {100*C_quant/tot:5.1f}%")
    print(f"  D. 厳密（=/0）                            : {D_exact:6d}  {100*D_exact/tot:5.1f}%")
    print("="*78)
    nb_tot = A_cross + B_spur
    if nb_tot:
        print(f"  ⟹ 境界なしのうち **{100*B_spur/nb_tot:.1f}%** が 見かけ（affine で 締めうる 上限）")
    print("  （B が 小さければ affine は 無益・境界なしは 本物の 符号曖昧さ）")

if __name__ == "__main__":
    main()
