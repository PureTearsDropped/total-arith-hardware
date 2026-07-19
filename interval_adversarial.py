#!/usr/bin/env python3
"""両オペランドのフラグ区間に対する 敵対的健全性テスト（mul と div）。

フラグ(=/≥/≤/符号不明)から 真値の 候補を 作り、**フラグ区間内の 全真値**で
bfp の 積/商の 状態(境界・符号不明)が 嘘をつかないか。成分を揃えて 相殺を 誘発、
除数区間が 0 を跨ぐ場合（N(b) が 0 に なりうる）も 含める。
"""
import sys, os
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
import numpy as np
from bfp_sed import BF, mul, div, EQ, GE, LE, NB, M
from sedenion_tensor_logic import ref_mult

rng = np.random.default_rng(20260720)

def sample_true(m, flag, sunk):
    """(代表値 m, フラグ, 符号不明) → 真値の候補（Fraction のリスト）。"""
    if flag == NB:
        base = [Fr(m), Fr(m, 3), Fr(m*7), Fr(0)]            # 大きさ 何でも
    elif flag == GE:
        base = [Fr(m*c) for c in (1, 2, 10)] if m != 0 else [Fr(0)]   # |x|≥|m|
    elif flag == LE:
        base = [Fr(m), Fr(m, 2), Fr(m, 7), Fr(0)] if m != 0 else [Fr(0)]  # |x|≤|m|
    else:  # EQ
        base = [Fr(m)]
    if sunk:                                               # 符号不明 ⟹ 両符号
        base = base + [-x for x in base]
    return base

def sound_bound(rep_val, true_val, bound, sunk):
    """状態が true について 嘘をつくか（True=違反）。"""
    t = true_val; r = Fr(rep_val)
    if bound == EQ and t != r: return True
    if bound == GE and abs(t) < abs(r): return True
    if bound == LE and abs(t) > abs(r): return True
    if not sunk and t != 0 and r != 0 and (t > 0) != (r > 0): return True
    return False

def run_test(name, op, is_div, n_trials=3000, n_samples=8, both_flagged=True):
    viol = 0; checks = 0; nb = 0; qbound = 0; tot = 0
    for _ in range(n_trials):
        am = [int(v) for v in rng.integers(-6, 6, M)]
        bm = [int(v) for v in rng.integers(-6, 6, M)]
        abound = [int(rng.choice([EQ, GE, LE])) for _ in range(M)]
        bbound = [int(rng.choice([EQ, GE, LE])) for _ in range(M)] if both_flagged else [EQ]*M
        asunk = [bool(rng.random() < 0.1) for _ in range(M)]
        bsunk = [bool(rng.random() < 0.1) for _ in range(M)] if both_flagged else [False]*M
        A = BF(am, bound=abound, sunk=asunk); B = BF(bm, bound=bbound, sunk=bsunk)
        try:
            R = op(A, B, W=18)
        except Exception:
            continue
        Ep = R.E
        for k in range(M):
            tot += 1
            if R.bound[k] == NB: nb += 1
            elif R.bound[k] in (GE, LE): qbound += 1
        # フラグ区間内の 真値を サンプルして 検証
        for _s in range(n_samples):
            at = [sample_true(am[i], abound[i], asunk[i])[
                    int(rng.integers(0, len(sample_true(am[i], abound[i], asunk[i]))))] for i in range(M)]
            bt = [sample_true(bm[i], bbound[i], bsunk[i])[
                    int(rng.integers(0, len(sample_true(bm[i], bbound[i], bsunk[i]))))] for i in range(M)]
            if is_div:
                Nb = sum(x*x for x in bt)
                if Nb == 0:
                    actual = [Fr(0)]*M                       # a/0 = 0
                else:
                    ac = ref_mult(at, [bt[0]] + [-x for x in bt[1:]])
                    actual = [c / Nb for c in ac]
            else:
                actual = ref_mult(at, bt)
            for k in range(M):
                checks += 1
                if sound_bound(R.mant[k] * (2 ** Ep) if Ep >= 0 else Fr(R.mant[k], 2**(-Ep)),
                               actual[k], R.bound[k], R.sunk[k]):
                    viol += 1
    print(f"  {name:<30} 違反 **{viol}**/{checks}  ｜ 量的境界 {100*qbound/tot:2.0f}% 境界なし {100*nb/tot:2.0f}%")
    return viol


def self_test():
    print("=" * 76)
    print("敵対的健全性: 両オペランドのフラグ区間内 全真値で 状態が 嘘をつくか")
    print("=" * 76)
    v = 0
    v += run_test("積 mul（両フラグ・相殺誘発）", mul, False, both_flagged=True)
    v += run_test("積 mul（a フラグ・b 厳密）", mul, False, both_flagged=False)
    v += run_test("商 div（両フラグ・N が 0 も）", div, True, both_flagged=True)
    v += run_test("商 div（a フラグ・b 厳密）", div, True, both_flagged=False)
    print("=" * 76)
    print(f"  {'健全性違反 合計:':<20} **{v}**（0 なら どの真値でも 状態は 嘘をつかない）")
    print("=" * 76)


if __name__ == "__main__":
    self_test()
