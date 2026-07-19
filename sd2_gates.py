#!/usr/bin/env python3
"""底2 符号つき桁の 3:2 圧縮器を、AND/OR/NOT/XOR の 論理回路として実装。

**構造（最小に近い）**: 符号つき桁 t = p − n を (p,n) で持つと、
    s = x+y+c = (xp+yp+cp) − (xn+yn+cn)
      = (2·Pc + Ps0) − (2·Nc + Ns0)              ← p レール と n レールに 全加算器 1 個ずつ
      = 2·(Pc − Nc) + (Ps0 − Ns0)
    ⟹ **low = Ps0 − Ns0、high = Pc − Nc**（どちらも {−1,0,1}）

  Pc, Nc = 多数決（全加算器の桁上げ）／ Ps0, Ns0 = パリティ（全加算器の和）
  ＝ **正レールの全加算器 + 負レールの全加算器 + レールごとの引き算**。
  桁上げ伝播のない carry-free 圧縮（Avizienis）。redundant（複数表現可）だが 値は厳密に保たれる。

  ゲートは (p,n) 上の AND/OR/NOT/XOR のみ（現行 sed の enc/dec/bnot と同じ土俵）。
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
from sedenion_tensor_logic import enc, dec, OMEGA, ref_mult, M, gate9
from sd2_core import to_sd2, from_sd2

# ------------------------------------------------- 論理ゲート（呼び出し回数を数える）
def AND(a, b, st): st['and'] += 1; return a & b
def OR (a, b, st): st['or']  += 1; return a | b
def NOT(a,    st): st['not'] += 1; return 1 - a
def XOR(a, b, st): st['xor'] += 1; return a ^ b

def full_adder(a, b, c, st):
    """教科書の最小全加算器（a⊕b を共有）: (sum=パリティ, carry=多数決)。5 ゲート。"""
    ab = XOR(a, b, st)
    s  = XOR(ab, c, st)                       # sum = a⊕b⊕c
    cy = OR(AND(a, b, st), AND(c, ab, st), st)  # carry = ab ∨ c·(a⊕b)
    return s, cy

# ------------------------------------------------- 3:2 圧縮器（論理回路）
def sd2_compress3_gates(x, y, c, st):
    """x+y+c = low + 2·high（low, high ∈ {−1,0,1}）を 論理ゲートで。

    正レール (xp,yp,cp) と 負レール (xn,yn,cn) に 全加算器を 1 個ずつ。
    low = Ps0 − Ns0（和どうしの差）、high = Pc − Nc（桁上げどうしの差）。
    """
    xp, xn = enc(x); yp, yn = enc(y); cp, cn = enc(c)
    Ps0, Pc = full_adder(xp, yp, cp, st)      # 正レール全加算器
    Ns0, Nc = full_adder(xn, yn, cn, st)      # 負レール全加算器
    nNs0 = NOT(Ns0, st); nPs0 = NOT(Ps0, st)
    nNc  = NOT(Nc,  st); nPc  = NOT(Pc,  st)
    lp = AND(Ps0, nNs0, st); ln = AND(Ns0, nPs0, st)   # low  = Ps0 − Ns0（両立時 0）
    hp = AND(Pc,  nNc,  st); hn = AND(Nc,  nPc,  st)   # high = Pc − Nc
    return dec(lp, ln), dec(hp, hn)

# ------------------------------------------------- 圧縮木・リップル（ゲート版）
def static_tree_g(columns, st):
    cols = [list(c) for c in columns]
    k = 0
    while k < len(cols):
        while len(cols[k]) > 2:
            x = cols[k].pop(); y = cols[k].pop(); z = cols[k].pop()
            low, high = sd2_compress3_gates(x, y, z, st)
            cols[k].append(low)
            if k + 1 >= len(cols): cols.append([])
            cols[k + 1].append(high)
        k += 1
    return cols

def ripple_g(cols, st):
    out, carry = [], 0
    for c in cols:
        r0 = c[0] if len(c) > 0 else 0
        r1 = c[1] if len(c) > 1 else 0
        low, carry = sd2_compress3_gates(r0, r1, carry, st)
        out.append(low)
    out.append(carry)
    return out

def sedenion_mult_gates(xw, yw, st):
    """セデニオン乗算（底2）を、圧縮器まで 論理ゲートで。gate9 は桁の積。"""
    K1, K2 = len(xw[0]), len(yw[0])
    cols = [[[] for _ in range(K1 + K2)] for _ in range(M)]
    for i in range(M):
        for j in range(M):
            k, s = i ^ j, OMEGA[i, j]
            for p in range(K1):
                for q in range(K2):
                    t = gate9(xw[i][p], yw[j][q]); st['g9'] += 1
                    cols[k][p + q].append(s * t)
    return [ripple_g(static_tree_g(cols[k], st), st) for k in range(M)]


def self_test(seed=20260718):
    import numpy as np
    rng = np.random.default_rng(seed)

    print("=" * 74)
    print("① 圧縮器は 27 通り すべて 値を保つ（low + 2·high = x+y+c、各 ∈ {−1,0,1}）")
    print("=" * 74)
    st = dict(and_=0)
    for x in (-1, 0, 1):
        for y in (-1, 0, 1):
            for c in (-1, 0, 1):
                low, high = sd2_compress3_gates(x, y, c, {'and':0,'or':0,'not':0,'xor':0})
                assert low + 2 * high == x + y + c, (x, y, c, low, high)
                assert low in (-1, 0, 1) and high in (-1, 0, 1)
    print("  27/27 値保存・桁は {−1,0,1} 内  ✓（redundant だが値は厳密）")

    print()
    print("=" * 74)
    print("② ゲート数（1 圧縮器あたり）")
    print("=" * 74)
    st = {'and': 0, 'or': 0, 'not': 0, 'xor': 0}
    sd2_compress3_gates(1, -1, 1, st)
    total = sum(st.values())
    print(f"    AND {st['and']}  OR {st['or']}  NOT {st['not']}  XOR {st['xor']}  ＝ **{total} ゲート**")
    print(f"    内訳: 全加算器×2（各5：XOR2+AND2+OR1）＋ レール引き算（AND4+NOT4）")
    xor_expanded = total - st['xor'] + st['xor'] * 3     # XOR = (a∧¬b)∨(¬a∧b) ≈ 3 段
    print(f"    （XOR を AND/OR/NOT に展開すると 約 **{xor_expanded} ゲート**。参考: 底3 gate27 は約30式項）")

    print()
    print("=" * 74)
    print("③ 論理ゲートだけで セデニオン乗算 == 参照（底2、端から端まで）")
    print("=" * 74)
    K = 6
    st = {'and':0,'or':0,'not':0,'xor':0,'g9':0}
    for _ in range(10):
        x = [int(v) for v in rng.integers(-2**K // 2, 2**K // 2, M)]
        y = [int(v) for v in rng.integers(-2**K // 2, 2**K // 2, M)]
        got = [from_sd2(w) for w in
               sedenion_mult_gates([to_sd2(c, K) for c in x], [to_sd2(c, K) for c in y], st)]
        assert got == ref_mult(x, y), (got, ref_mult(x, y))
    print("  10 件 厳密一致（圧縮器を 論理ゲートに落としても 厳密性は保たれる）  ✓")
    print(f"  1 積あたり ゲート: gate9 {st['g9']//10}, 圧縮器 AND/OR/NOT/XOR "
          f"{(st['and']+st['or']+st['not']+st['xor'])//10}")

    print()
    print("=" * 74)
    print("④ 零因子 (e1+e10)(e4−e15) = 0 も 論理ゲートで 厳密 0")
    print("=" * 74)
    a = [0]*M; a[1] = 1; a[10] = 1
    b = [0]*M; b[4] = 1; b[15] = -1
    got = [from_sd2(w) for w in
           sedenion_mult_gates([to_sd2(c, 4) for c in a], [to_sd2(c, 4) for c in b],
                               {'and':0,'or':0,'not':0,'xor':0,'g9':0})]
    assert all(g == 0 for g in got), got
    print("  全成分 0  ✓")
    print()
    print("すべて通過 — 3:2 圧縮器は AND/OR/NOT/XOR の組み合わせで動く。")
    print("構造 = 正レール全加算器 + 負レール全加算器 + レール引き算（最小に近い）。")


if __name__ == "__main__":
    self_test()
