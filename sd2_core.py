#!/usr/bin/env python3
"""底2・符号つき桁（signed-digit binary）の芯 — SPEC §1,§2,§5 の実装。

桁 d ∈ {−1, 0, +1}、位置の重み 2^k（← 現行 sed の底3 を底2 にした版）。
底が入っていた 2 箇所だけを差し替える（SPEC §1.2 で切り分け済み）:
    gate27（3:2 圧縮器、桁上げ重み 3）  →  **sd2_compress3**（桁上げ重み 2）
    to_trits / from_trits（3^k）       →  **to_sd2 / from_sd2**（2^k）
底に依存しない gate9（桁の積）と 重み配置 p+q は そのまま流用する。
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
from sedenion_tensor_logic import OMEGA, ref_mult, M, gate9   # gate9 は底に依存しない

# ------------------------------------------------- 符号化（底2、符号つき桁）
def to_sd2(n, width):
    """整数 → {−1,0,+1} 桁の列（重み 2^k、下位から）。符号は全桁に配る（sign-magnitude）。"""
    n = int(n); neg = n < 0; a = abs(n)
    out = [((a >> k) & 1) for k in range(width)]
    if a >> width:
        raise OverflowError(f"width {width} too small for {n}")
    return [-d for d in out] if neg else out

def from_sd2(word):
    """{−1,0,+1} 桁の列 → 整数（値 = Σ d_k · 2^k）。"""
    return int(sum(d * (1 << k) for k, d in enumerate(word)))

# ------------------------------------------------- 底2 の 3:2 圧縮器
def sd2_compress3(x, y, c):
    """3 桁 {−1,0,+1} を x+y+c = low + 2·high に圧縮（low, high ∈ {−1,0,+1}）。

    底3 の gate27（x+y+c = same + 3·high）の 底2 版。桁上げ high は 重み 2 へ。
    x+y+c ∈ [−3,3] は low + 2·high（各 ∈ {−1,0,1}）で必ず表せる:
        3 = 1+2·1   2 = 0+2·1   1 = 1   0 = 0   −1 = −1   −2 = 0+2·(−1)   −3 = −1+2·(−1)
    """
    s = x + y + c
    high = 0
    if s > 1:
        high = 1; s -= 2
    elif s < -1:
        high = -1; s += 2
    return s, high            # low(=s), high

# ------------------------------------------------- 圧縮木 + 最終桁上げ（底2版）
def static_tree_sd2(columns):
    """各重み位置を 高々 2 桁に落とす（静的スケジュール、桁上げは w→w+1）。"""
    cols = [list(c) for c in columns]
    k = 0
    while k < len(cols):
        while len(cols[k]) > 2:
            x = cols[k].pop(); y = cols[k].pop(); z = cols[k].pop()
            low, high = sd2_compress3(x, y, z)
            cols[k].append(low)
            if k + 1 >= len(cols):
                cols.append([])
            cols[k + 1].append(high)
        k += 1
    return cols

def ripple_sd2(cols):
    """(row0, row1) の 2 桁/重み を 単一桁へ（唯一のリップル）。"""
    out, carry = [], 0
    for c in cols:
        r0 = c[0] if len(c) > 0 else 0
        r1 = c[1] if len(c) > 1 else 0
        low, carry = sd2_compress3(r0, r1, carry)
        out.append(low)
    out.append(carry)
    return out

def sd2_add(a, b):
    """2 数（{−1,0,1} 桁の列）を足す — carry-free 圧縮 + 単一リップル。"""
    W = max(len(a), len(b))
    cols = [[] for _ in range(W)]
    for k in range(W):
        if k < len(a): cols[k].append(a[k])
        if k < len(b): cols[k].append(b[k])
    return ripple_sd2(static_tree_sd2(cols))

def sd2_mul(a, b):
    """2 数の積 — 部分積 a_i·b_j を 重み i+j に置く（p+q 配置）。"""
    cols = [[] for _ in range(len(a) + len(b))]
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            cols[i + j].append(gate9(ai, bj))     # 桁の積（底に依存しない）
    return ripple_sd2(static_tree_sd2(cols))

# ------------------------------------------------- セデニオン乗算（底2、厳密整数）
def sedenion_mult_sd2(xw, yw):
    """xw, yw: 16 個の {−1,0,1} 桁ワード → 16 個。SPEC §5 の 底2 版。
    符号 OMEGA = 配線の入れ替え、経路 k = i^j、桁上げは底2 圧縮器。"""
    K1, K2 = len(xw[0]), len(yw[0])
    cols = [[[] for _ in range(K1 + K2)] for _ in range(M)]
    for i in range(M):
        for j in range(M):
            k, s = i ^ j, OMEGA[i, j]
            for p in range(K1):
                for q in range(K2):
                    t = gate9(xw[i][p], yw[j][q])
                    cols[k][p + q].append(s * t)          # s·t ∈ {−1,0,1}
    return [ripple_sd2(static_tree_sd2(cols[k])) for k in range(M)]


def self_test(seed=20260718):
    import numpy as np
    rng = np.random.default_rng(seed)

    # 符号化の往復
    for _ in range(200):
        n = int(rng.integers(-10**6, 10**6))
        assert from_sd2(to_sd2(n, 24)) == n
    print("to_sd2/from_sd2 往復: 200 件 厳密")

    # スカラーの加算・乗算
    for _ in range(200):
        a = int(rng.integers(-1000, 1000)); b = int(rng.integers(-1000, 1000))
        assert from_sd2(sd2_add(to_sd2(a, 16), to_sd2(b, 16))) == a + b
        assert from_sd2(sd2_mul(to_sd2(a, 16), to_sd2(b, 16))) == a * b
    print("sd2_add / sd2_mul: 200 件 厳密（底2 の carry-save 圧縮）")

    # 圧縮器の網羅チェック
    for x in (-1, 0, 1):
        for y in (-1, 0, 1):
            for c in (-1, 0, 1):
                low, high = sd2_compress3(x, y, c)
                assert low + 2 * high == x + y + c
                assert low in (-1, 0, 1) and high in (-1, 0, 1)
    print("sd2_compress3: 27 通り 全て x+y+c = low + 2·high、各桁 ∈ {−1,0,1}")

    # セデニオン乗算 == 参照
    K = 6
    for _ in range(10):
        x = [int(v) for v in rng.integers(-2**K // 2, 2**K // 2, M)]
        y = [int(v) for v in rng.integers(-2**K // 2, 2**K // 2, M)]
        got = [from_sd2(w) for w in
               sedenion_mult_sd2([to_sd2(c, K) for c in x], [to_sd2(c, K) for c in y])]
        assert got == ref_mult(x, y), (got, ref_mult(x, y))
    print("sedenion_mult_sd2 == ref_mult: 10 件 厳密（底2 セデニオン乗算）")

    # 零因子 (e1+e10)(e4−e15) = 0 が 厳密に 0
    a = [0]*M; a[1] = 1; a[10] = 1
    b = [0]*M; b[4] = 1; b[15] = -1
    got = [from_sd2(w) for w in
           sedenion_mult_sd2([to_sd2(c, 4) for c in a], [to_sd2(c, 4) for c in b])]
    assert all(g == 0 for g in got), got
    print("零因子 (e1+e10)(e4−e15) = 0: 底2 でも 厳密に 0")
    print("すべて通過 — 底2 符号つき桁の芯は 底3 と同じ厳密性を持つ。")


if __name__ == "__main__":
    self_test()
