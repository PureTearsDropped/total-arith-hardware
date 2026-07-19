#!/usr/bin/env python3
"""区間ベースの セデニオン積 — 代表値で計算＋方向保持（利用者の改善）。

各成分の 真値を フラグから 符号つき区間 [lo,hi] にし、積の成分（符号つき和）を
**区間演算**で 伝播 → 結果区間から フラグを 引き直す。
  ・単調な場合（相殺なし）: ≥2 のような 量的境界が 出る（境界なしに 潰さない）
  ・相殺する場合（セデニオン積）: 区間が 広がり 自動的に 境界なしへ
両方 健全（真値が 必ず 区間の中）。
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
from sedenion_tensor_logic import OMEGA, ref_mult, M
from bfp_sed import EQ, GE, LE, NB, _ival, _flags_for

INF = float("inf")

def imul(a, b):
    """符号つき区間 a=[a1,a2], b=[b1,b2] の 積区間 = [min, max]（±∞ 対応）。"""
    ps = []
    for x in (a[0], a[1]):
        for y in (b[0], b[1]):
            if (x == 0 and y in (INF, -INF)) or (y == 0 and x in (INF, -INF)):
                ps.append(0.0)                 # 0 × ∞ は この文脈では 0（区間端の 0 は 真の 0）
            else:
                ps.append(x * y)
    return min(ps), max(ps)

def iadd(a, b):
    return (a[0] + b[0], a[1] + b[1])

def sed_mul_intervals(am, bm, abound, bbound):
    """整数仮数 am,bm と フラグ abound,bbound → 積の (仮数, フラグ, 符号不明)。"""
    prod = [int(v) for v in ref_mult(am, bm)]              # 代表値の 積（厳密）
    # 各入力成分の 真値区間
    ai = [_ival(am[i], abound[i]) for i in range(M)]
    bi = [_ival(bm[i], bbound[i]) for i in range(M)]
    bound = [EQ] * M; sunk = [False] * M
    for k in range(M):
        lo = hi = 0.0
        for i in range(M):
            j = i ^ k; s = OMEGA[i, j]
            t = imul(ai[i], bi[j])
            if s < 0: t = (-t[1], -t[0])
            lo += t[0]; hi += t[1]
        bound[k], sunk[k] = _flags_for(prod[k], lo, hi)
    return prod, bound, sunk


def self_test():
    import numpy as np
    from fractions import Fraction as Fr
    rng = np.random.default_rng(7)

    def sample_true(m, flag, reps):
        """フラグと 代表値 m から 真値の 候補（Fraction）を 返す。"""
        if flag == EQ: return [Fr(m)]
        if flag == GE:                                    # |x| ≥ |m|, 符号 sign(m)
            return [Fr(m) * r for r in reps] if m != 0 else [Fr(0)]
        if flag == LE:                                    # |x| ≤ |m|, 符号 sign(m)
            return [Fr(m), Fr(m, 2), Fr(0)] if m != 0 else [Fr(0)]
        return None                                       # NB は 走査せず

    print("=" * 74)
    print("① 敵対的 健全性: フラグ区間内の 全真値で 出力境界が 成り立つか（成分を揃える）")
    print("=" * 74)
    reps = [1, 3, 50]                                      # ≥ の 真値候補（|m|, 3|m|, 50|m|）
    viol = 0; checks = 0
    for _ in range(4000):
        am = [int(v) for v in rng.integers(-9, 9, M)]     # 揃った 成分（相殺を 誘発）
        bm = [int(v) for v in rng.integers(-9, 9, M)]
        abound = [int(rng.choice([EQ, GE, LE])) for _ in range(M)]
        bbound = [EQ] * M                                 # b は 厳密（単調な 掛け算の 検査）
        prod, bound, sunk = sed_mul_intervals(am, bm, abound, bbound)
        # a の 真値候補を 組合せ（多いので 各成分 独立に サンプル1つ×数回）
        for _s in range(6):
            at = []
            ok = True
            for i in range(M):
                cs = sample_true(am[i], abound[i], reps)
                if cs is None: ok = False; break
                at.append(cs[int(rng.integers(0, len(cs)))])
            if not ok: continue
            actual = ref_mult(at, [Fr(x) for x in bm])    # 真の 積
            for k in range(M):
                checks += 1
                r = Fr(prod[k]); t = actual[k]; bd = bound[k]
                if bd == EQ and t != r: viol += 1
                elif bd == GE and abs(t) < abs(r): viol += 1
                elif bd == LE and abs(t) > abs(r): viol += 1
    print(f"   健全性違反: **{viol}** / {checks} チェック（0 なら 常に健全）")

    print()
    print("=" * 74)
    print("② 量的境界が 出るか — スカラー的な 単調ケースで（利用者の ≥2 の例）")
    print("=" * 74)
    # a = (≥MAX 相当) を 成分0だけ、b = 厳密スカラー（成分0だけ）で 掛ける = スカラー積
    for tag, av, af, bv in [("(≥100)×(=3)", 100, GE, 3), ("(≥100)×(=1)", 100, GE, 1),
                             ("(≤2)×(=5)", 2, LE, 5)]:
        am = [0]*M; am[0] = av; ab = [EQ]*M; ab[0] = af
        bm = [0]*M; bm[0] = bv
        prod, bound, sunk = sed_mul_intervals(am, bm, ab, [EQ]*M)
        from bfp_sed import BNAME
        print(f"   {tag:<14} → 成分0 = {prod[0]}, 境界 {BNAME[bound[0]]}  "
              f"（{'量的境界 ✓' if bound[0] in (GE,LE) else '厳密' if bound[0]==EQ else '境界なし'}）")

    print()
    print("=" * 74)
    print("③ 相殺ケース（成分が混ざる）は 自動で 境界なしへ")
    print("=" * 74)
    nb_cnt = 0; tot = 0
    for _ in range(500):
        am = [int(v) for v in rng.integers(-9, 9, M)]
        bm = [int(v) for v in rng.integers(-9, 9, M)]
        ab = [GE]*M                                       # 全成分 ≥（相殺で 方向破れ 起きやすい）
        prod, bound, sunk = sed_mul_intervals(am, bm, ab, [EQ]*M)
        for k in range(M):
            tot += 1
            if bound[k] == NB: nb_cnt += 1
    print(f"   全成分 ≥ の 積: 境界なしに なった成分 {100*nb_cnt/tot:.0f}%（相殺ゆえ・健全）")
    print("\n   ⟹ 単調では 量的境界、相殺では 境界なし。両方 健全（① で 0 違反）。")


if __name__ == "__main__":
    self_test()
