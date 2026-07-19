#!/usr/bin/env python3
# ⚠️ 生成AI使用・要検証
"""表現論のレンズ — 非可換群の配線 = 行列積の直和（Wedderburn 分解）を 実測で。

  群のフーリエ変換 â(ρ) = Σ_g a_g ρ(g)（既約表現 ρ ごと）は 畳み込みを 積に:
      (a ∗ b)^(ρ) = â(ρ) · b̂(ρ)      ← 1次元 ρ なら スカラー積、2次元なら 2×2 行列積
  ⟹ **一つの 群配線パス = 複数の 独立な 積を 同時に 実行**（表現論の SIMD）。

  検証:
   ① D4（位数8）: R[D4] ≅ R⁴ ⊕ M₂(R) — 配線1回 = スカラー4 + 2×2行列積1 同時
   ② Q8（位数8）: R[Q8] ≅ R⁴ ⊕ H — **四元数の配線は Q8 群配線の 1ブロック**
   ③ S3（位数6）: R[S3] ≅ R² ⊕ M₂(R) — 置換の群（対称群）でも 同じ
  （Cohn–Umans「群代数で 行列積」の 入り口が この分解）
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from fractions import Fraction as Fr
from nd_algebra import cd_omega, ref_mult_M
from wiring_registry import _d4_mul, _d4_elems

rng = np.random.default_rng(20260801)


def group_conv(mul, a, b):
    """群配線の 畳み込み c[g·h] += a[g]b[h]（レジストリの bilinear と 同じもの）。"""
    n = len(mul); c = [0]*n
    for i in range(n):
        for j in range(n):
            c[mul[i][j]] += a[i]*b[j]
    return c


# ---------------------------------------------------------------- ① D4
def d4_irreps():
    """D4 の 既約表現: 1次元×4 ＋ 2次元×1（元 (f,k) = s^f r^k）。"""
    R = np.array([[0, -1], [1, 0]])                      # r → 90° 回転
    S = np.array([[1, 0], [0, -1]])                      # s → 鏡映
    reps = []
    for alpha in (1, -1):                                # 1次元: s→α, r→β
        for beta in (1, -1):
            reps.append(('1d', lambda f, k, a=alpha, b=beta: np.array([[a**f * b**k]])))
    reps.append(('2d', lambda f, k: np.linalg.matrix_power(S, f) @ np.linalg.matrix_power(R, k)))
    return reps

def fourier(a, elems, rep):
    return sum(a[i] * rep(f, k) for i, (f, k) in enumerate(elems))


def test_d4():
    print("=" * 78)
    print("① D4: 群配線 1 回 = スカラー積 4 個 ＋ 2×2 行列積 1 個（同時実行）")
    print("=" * 78)
    reps = d4_irreps()
    bad = 0; trials = 2000
    for _ in range(trials):
        a = [int(v) for v in rng.integers(-9, 10, 8)]
        b = [int(v) for v in rng.integers(-9, 10, 8)]
        c = group_conv(_d4_mul, a, b)                    # 群配線（1 パス）
        for _, rep in reps:                              # 各既約で ĉ = â·b̂ か
            ca = fourier(a, _d4_elems, rep); cb = fourier(b, _d4_elems, rep)
            cc = fourier(c, _d4_elems, rep)
            if not np.allclose(cc, ca @ cb): bad += 1
    dims = [1, 1, 1, 1, 2]
    print(f"  既約表現の次元: {dims} → Σd² = {sum(d*d for d in dims)} = 群位数 8 ✓（Wedderburn）")
    print(f"  ĉ(ρ) == â(ρ)·b̂(ρ)（全 5 既約 × {trials} 試行）: 違反 **{bad}**")
    print(f"  ⟹ D4 配線 1 回 ＝ スカラー×4 ＋ 2×2 行列積×1 を **同時に**（基底変換は 固定＝配線）")


# ---------------------------------------------------------------- ② Q8
def q8_table():
    """四元数群 Q8 = {±1, ±i, ±j, ±k}。積表と、H への 埋め込み φ(g) ∈ R⁴。"""
    units = [(1, 0), (-1, 0), (1, 1), (-1, 1), (1, 2), (-1, 2), (1, 3), (-1, 3)]  # (符号, 軸0=1,1=i,2=j,3=k)
    OM = cd_omega(4)                                     # 四元数の 符号表（既検証）
    def mulq(u, v):
        (sa, ea), (sb, eb) = u, v
        k = ea ^ eb; s = sa * sb * int(OM[ea, eb])
        return (s, k)
    idx = {u: i for i, u in enumerate(units)}
    mul = [[idx[mulq(units[i], units[j])] for j in range(8)] for i in range(8)]
    def phi(i):
        s, e = units[i]; v = [0]*4; v[e] = s
        return v                                          # H の 基底ベクトルとして
    return mul, units, phi, idx

def test_q8():
    print()
    print("=" * 78)
    print("② Q8: 群配線の 1 ブロックが **四元数の配線 そのもの**（R[Q8] ≅ R⁴ ⊕ H）")
    print("=" * 78)
    mul, units, phi, idx = q8_table()
    OM = cd_omega(4)
    # 群の性質（独立検証）
    assoc = all(mul[mul[i][j]][k] == mul[i][mul[j][k]] for i in range(8) for j in range(8) for k in range(8))
    noncomm = any(mul[i][j] != mul[j][i] for i in range(8) for j in range(8))
    print(f"  Q8 群表: 結合律 {assoc}・非可換 {noncomm}・i·j=k {mul[idx[(1,1)]][idx[(1,2)]] == idx[(1,3)]}"
          f"・j·i=−k {mul[idx[(1,2)]][idx[(1,1)]] == idx[(-1,3)]}")
    # H ブロック: Φ(a) = Σ a_g φ(g)、Φ(a∗b) == Φ(a)·Φ(b)（四元数積 = レジストリの quaternion）
    bad = 0; trials = 2000
    for _ in range(trials):
        a = [int(v) for v in rng.integers(-9, 10, 8)]
        b = [int(v) for v in rng.integers(-9, 10, 8)]
        c = group_conv(mul, a, b)
        Phi = lambda x: [sum(x[i]*phi(i)[t] for i in range(8)) for t in range(4)]
        lhs = Phi(c)
        rhs = ref_mult_M(Phi(a), Phi(b), OM, 4)          # 四元数配線（既検証）で 掛ける
        if lhs != rhs: bad += 1
    # 1次元×4（Q8/{±1} ≅ Z2×Z2 の 指標）
    bad1 = 0
    for _ in range(500):
        a = [int(v) for v in rng.integers(-9, 10, 8)]
        b = [int(v) for v in rng.integers(-9, 10, 8)]
        c = group_conv(mul, a, b)
        for ci in (1, -1):
            for cj in (1, -1):
                def chi(i):
                    s, e = units[i]                       # χ(±1)=1, χ(i)=ci, χ(j)=cj, χ(k)=ci·cj
                    return [1, ci, cj, ci*cj][e]
                fa = sum(a[i]*chi(i) for i in range(8)); fb = sum(b[i]*chi(i) for i in range(8))
                fc = sum(c[i]*chi(i) for i in range(8))
                if fc != fa*fb: bad1 += 1
    print(f"  Φ(a∗b) == Φ(a)·Φ(b)（H=四元数ブロック・{trials} 試行）: 違反 **{bad}**")
    print(f"  1次元指標 ×4: 違反 **{bad1}**   次元 4·1 + 4(H) = 8 ✓")
    print(f"  ⟹ **四元数という『数』は、群 Q8 の配線の 1 ブロック**（数の配線 ⊂ 群の配線）")


# ---------------------------------------------------------------- ③ S3
def s3_table():
    """対称群 S3（3 文字の置換・位数 6）。経路 = 置換の合成。"""
    import itertools
    perms = list(itertools.permutations(range(3)))
    idx = {p: i for i, p in enumerate(perms)}
    comp = lambda p, q: tuple(p[q[t]] for t in range(3))          # (p∘q)(t) = p[q[t]]
    mul = [[idx[comp(perms[i], perms[j])] for j in range(6)] for i in range(6)]
    return mul, perms

def s3_std_rep(p):
    """標準 2 次元表現: {x∈R³: Σx=0} に 置換行列を 制限（基底 e0−e1, e1−e2・整数行列）。"""
    P = np.zeros((3, 3), dtype=int)
    for t in range(3): P[p[t], t] = 1
    B = np.array([[1, 0], [-1, 1], [0, -1]])                      # 基底: e0−e1, e1−e2
    # P·B の 各列を B の 列の 整数結合で 表す（(x0,x1,x2), Σ=0 ⟹ 係数 = (x0, −x2)）
    PB = P @ B
    return np.array([[PB[0, 0], PB[0, 1]], [-PB[2, 0], -PB[2, 1]]])

def test_s3():
    print()
    print("=" * 78)
    print("③ S3（対称群・置換の合成が 経路）: R[S3] ≅ R² ⊕ M₂(R)")
    print("=" * 78)
    mul, perms = s3_table()
    sgn = lambda p: int(np.sign(np.linalg.det(np.eye(3)[list(p)])))
    bad = 0; trials = 2000
    for _ in range(trials):
        a = [int(v) for v in rng.integers(-9, 10, 6)]
        b = [int(v) for v in rng.integers(-9, 10, 6)]
        c = group_conv(mul, a, b)
        # 1次元: 自明 と 符号
        for ch in (lambda p: 1, sgn):
            fa = sum(a[i]*ch(perms[i]) for i in range(6)); fb = sum(b[i]*ch(perms[i]) for i in range(6))
            if sum(c[i]*ch(perms[i]) for i in range(6)) != fa*fb: bad += 1
        # 2次元 標準表現
        F = lambda x: sum(x[i] * s3_std_rep(perms[i]) for i in range(6))
        if not np.array_equal(F(c), F(a) @ F(b)): bad += 1
    print(f"  自明・符号・標準2次元 の 3 既約で ĉ == â·b̂（{trials} 試行）: 違反 **{bad}**")
    print(f"  次元 1 + 1 + 4 = 6 ✓   ⟹ 置換の群でも 配線 1 回 = スカラー2 + 2×2行列積1")


if __name__ == "__main__":
    test_d4()
    test_q8()
    test_s3()
    print()
    print("=" * 78)
    print("まとめ — 表現論のレンズ")
    print("=" * 78)
    print("""  非可換群の 配線は、固定の 基底変換（フーリエ＝これも配線）で **行列積の直和** に 割れる:
    D4 → R⁴ ⊕ M₂(R)   ／   Q8 → R⁴ ⊕ H（四元数！）   ／   S3 → R² ⊕ M₂(R)
  ⟹ ① 一つの 群配線パス = 複数の 独立な 積の 同時実行（表現論の SIMD）
    ② 『数』（四元数）は 群配線の 1 ブロック — 数の配線 ⊂ 群の配線
    ③ 大きい非可換群ほど 大きい行列ブロックを 含む ＝ Cohn–Umans（群代数で 高速行列積）の 入り口""")
