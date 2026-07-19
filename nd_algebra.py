#!/usr/bin/env python3
"""任意の 2^k 次元（Cayley–Dickson の塔）へ 一般化 — 実/複素/四元/八元/セデニオン/32/64…

このシステムで 次元に依存するのは **符号表 OMEGA だけ**。それは Cayley–Dickson 構成で
再帰的に 定義され、**任意の 2^k で XOR routing（e_i·e_j = ±e_{i⊕j}）が 成り立つ**。
桁の積 gate9・sd2 圧縮器・ブロック浮動・状態は 全て 成分ごと ⟹ 次元に依らず そのまま動く。
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
import numpy as np
from sedenion_tensor_logic import _cd, _conj, gate9
from sd2_core import to_sd2, from_sd2, static_tree_sd2, ripple_sd2

def cd_omega(M):
    """次元 M=2^k の Cayley–Dickson 符号表。OMEGA[i,j] ∈ {−1,+1}, 経路 = i⊕j。"""
    E = np.eye(M)
    OM = np.zeros((M, M), dtype=int)
    for i in range(M):
        for j in range(M):
            v = _cd(E[i], E[j])
            k = int(np.argmax(np.abs(v)))
            assert k == (i ^ j), f"XOR routing 破れ M={M} ({i},{j})"
            OM[i, j] = int(np.sign(v[k]))
    return OM

def ref_mult_M(x, y, OM, M):
    r = [0] * M
    for i in range(M):
        for j in range(M):
            r[i ^ j] += OM[i, j] * x[i] * y[j]
    return r

def sd2_mult_M(xw, yw, OM, M):
    """底2 符号つき桁で 次元 M のセデニオン様積（成分ごと・次元非依存の骨格）。"""
    K1, K2 = len(xw[0]), len(yw[0])
    cols = [[[] for _ in range(K1 + K2)] for _ in range(M)]
    for i in range(M):
        for j in range(M):
            k, s = i ^ j, OM[i, j]
            for p in range(K1):
                for q in range(K2):
                    cols[k][p + q].append(s * gate9(xw[i][p], yw[j][q]))
    return [ripple_sd2(static_tree_sd2(cols[k])) for k in range(M)]

def properties(OM, M):
    """可換・結合・零因子（Cayley–Dickson で 段階的に 失われる性質）。"""
    rng = np.random.default_rng(0)
    comm = assoc = True
    for _ in range(60):
        a = [int(v) for v in rng.integers(-4, 4, M)]
        b = [int(v) for v in rng.integers(-4, 4, M)]
        c = [int(v) for v in rng.integers(-4, 4, M)]
        if ref_mult_M(a, b, OM, M) != ref_mult_M(b, a, OM, M): comm = False
        if ref_mult_M(ref_mult_M(a, b, OM, M), c, OM, M) != \
           ref_mult_M(a, ref_mult_M(b, c, OM, M), OM, M): assoc = False
    # 零因子: a≠0, b≠0 で a·b=0 が 在るか（既知の構成 e_i+e_j 型を探索）
    zdiv = False
    for i in range(1, M):
        for j in range(i + 1, M):
            for p in range(1, M):
                for q in range(p + 1, M):
                    a = [0]*M; a[i]=1; a[j]=1
                    b = [0]*M; b[p]=1; b[q]=-1
                    if any(a) and any(b) and all(v == 0 for v in ref_mult_M(a, b, OM, M)):
                        zdiv = True; break
                if zdiv: break
            if zdiv: break
        if zdiv: break
    return comm, assoc, zdiv


def self_test():
    NAMES = {1:"実数", 2:"複素", 4:"四元数", 8:"八元数", 16:"セデニオン", 32:"32次元", 64:"64次元"}
    print("=" * 78)
    print("Cayley–Dickson の塔 — 任意の 2^k 次元で このシステムが 動くか")
    print("=" * 78)
    print(f"  {'M':>4} {'名前':<10}{'XOR routing':>12}{'可換':>6}{'結合':>6}{'零因子':>8}{'底2積==ref':>12}")
    for k in range(7):
        M = 1 << k
        OM = cd_omega(M)                          # XOR routing は assert で 内部確認済み
        comm, assoc, zdiv = properties(OM, M)
        # 底2 符号つき桁の積が 参照と一致するか
        rng = np.random.default_rng(7)
        ok = True
        for _ in range(5):
            x = [int(v) for v in rng.integers(-8, 8, M)]
            y = [int(v) for v in rng.integers(-8, 8, M)]
            got = [from_sd2(w) for w in sd2_mult_M([to_sd2(c, 6) for c in x],
                                                   [to_sd2(c, 6) for c in y], OM, M)]
            if got != ref_mult_M(x, y, OM, M): ok = False
        print(f"  {M:>4} {NAMES[M]:<10}{'✓':>12}{('✓' if comm else '×'):>6}"
              f"{('✓' if assoc else '×'):>6}{('有' if zdiv else '無'):>8}{('✓' if ok else '×'):>12}")
    print("""
  ⟹ **XOR routing と 底2積 == ref は 全次元で ✓**（符号表 OMEGA だけ 次元依存・自動生成）。
     Cayley–Dickson の 性質喪失が そのまま 現れる:
       M=2 まで 可換 → M=4（四元数）で 非可換 → M=8（八元数）で 非結合 → M=16 で 零因子。
     **八元数(8)まで 零因子なし（除算 綺麗）／セデニオン(16)以上 零因子あり（境界トリットで処理）。**

  ⟹ **一般化は ほぼ無料。** 次元に依存するのは OMEGA だけで、それは 再帰的に 定義される。
     桁の演算・圧縮器・ブロック浮動・状態（符号不明/境界）は 全て 成分ごと ⟹ そのまま 任意次元。""")
    print("  注: 真に任意の n（2の冪でない）は XOR routing を 持たない ⟹ 一般の構造定数テンソルが要り、")
    print("      乗算器フリーの利点は 失う。**2^k の塔が この設計の 自然な一般化。**")


if __name__ == "__main__":
    self_test()
