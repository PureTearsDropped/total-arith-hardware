#!/usr/bin/env python3
"""中心化形式の 効果測定: 有界近スカラー除数で
  (a) 健全性（真値サンプルで 状態が 嘘つかないか・0 違反 必須）
  (b) 締まり（境界なし → 量的境界に なった 割合）中心化 ON/OFF 比較。
"""
import sys, os
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
import numpy as np
import bfp_sed
from bfp_sed import BF, div, EQ, GE, LE, NB, M
from sedenion_tensor_logic import ref_mult

rng = np.random.default_rng(20260724)

def true_samples(m, flag, sunk):
    if flag == EQ:   base = [Fr(m)]
    elif flag == LE: base = [Fr(m), Fr(m,2), Fr(m,3), Fr(m,7), Fr(0)] if m != 0 else [Fr(0)]
    elif flag == GE: base = [Fr(m)*r for r in (1,2,5)] if m != 0 else [Fr(0)]
    else:            base = [Fr(m)*r for r in (1,3)] + [Fr(0)]
    if sunk: base = base + [-x for x in base]
    return base

def gen_bounded_near_scalar():
    """b_0 支配的（EQ/LE・非零大）+ 小 LE 摂動、a も 有界（EQ/LE）。箱 有界 ⟹ 中心化 発動。"""
    bm = [0]*M; bm[0] = int(rng.choice([-64, 64, -100, 100]))
    for i in range(1, M):
        if rng.random() < 0.5: bm[i] = int(rng.integers(-6, 7))
    am = [int(rng.integers(-40, 41)) for _ in range(M)]
    abd = [int(rng.choice([EQ, LE])) for _ in range(M)]       # 全有界
    bbd = [EQ]*M
    bbd[0] = int(rng.choice([EQ, LE]))
    for i in range(1, M):
        if bm[i] != 0 and rng.random() < 0.7: bbd[i] = LE     # 摂動は 有界 LE
    return am, bm, abd, bbd

def measure(use_centered, n_trials=1500, n_samples=10):
    bfp_sed.USE_CENTERED = use_centered
    viol = checks = 0
    nb = quant = exact = tot = 0
    for _ in range(n_trials):
        am, bm, abd, bbd = gen_bounded_near_scalar()
        if all(bm[i]==0 for i in range(1,M)): continue        # スカラー除数は 別分岐
        A = BF(am, bound=abd); B = BF(bm, bound=bbd)
        try: R = div(A, B, W=18)
        except Exception: continue
        Ep = R.E
        for k in range(M):
            tot += 1
            if R.bound[k] == NB: nb += 1
            elif R.bound[k] in (GE, LE): quant += 1
            else: exact += 1
        ar = [true_samples(am[i], abd[i], False) for i in range(M)]
        br = [true_samples(bm[i], bbd[i], False) for i in range(M)]
        for _s in range(n_samples):
            at = [ar[i][int(rng.integers(0,len(ar[i])))] for i in range(M)]
            bt = [br[i][int(rng.integers(0,len(br[i])))] for i in range(M)]
            Nb = sum(x*x for x in bt)
            actual = [Fr(0)]*M if Nb==0 else [c/Nb for c in ref_mult(at,[bt[0]]+[-x for x in bt[1:]])]
            for k in range(M):
                checks += 1
                rep = R.mant[k]*(2**Ep) if Ep>=0 else Fr(R.mant[k], 2**(-Ep))
                r = Fr(rep); t = actual[k]; bd = R.bound[k]
                lie = ((bd==EQ and t!=r) or (bd==GE and abs(t)<abs(r)) or
                       (bd==LE and abs(t)>abs(r)) or
                       (not R.sunk[k] and t!=0 and r!=0 and (t>0)!=(r>0)))
                if lie: viol += 1
    tag = "中心化 ON " if use_centered else "中心化 OFF"
    print(f"  [{tag}] 違反 **{viol}**/{checks}  ｜ 量的 {100*quant/tot:4.1f}%  境界なし {100*nb/tot:4.1f}%  厳密 {100*exact/tot:4.1f}%")
    return viol, quant, nb, tot

if __name__ == "__main__":
    print("="*82)
    print("中心化形式の効果: 有界近スカラー除数（大きさが決まる 実用域）")
    print("="*82)
    v_off, q_off, nb_off, tot = measure(False)
    v_on,  q_on,  nb_on,  _   = measure(True)
    print("="*82)
    print(f"  健全性: OFF 違反 {v_off}・ON 違反 {v_on}（両方 0 必須）")
    if nb_off:
        recovered = (nb_off - nb_on)
        print(f"  締まり: 境界なし {100*nb_off/tot:.1f}% → {100*nb_on/tot:.1f}%  "
              f"（**{100*recovered/nb_off:.1f}%** の 境界なしを 量的境界へ 回収）")
    print("="*82)
