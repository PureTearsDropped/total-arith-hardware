#!/usr/bin/env python3
"""並列アレイ: 同じ XOR データパスのユニットを **並べて(SIMT)・繋いで(タイル化GEMM)** 使う。

  原子演算 = 1 ユニット = 1 パスで b×b ブロック行列積（b=2^q）。＝ テンソルコアの 1 命令に相当。
    ・符号表 OMEGA を 差し替えるだけで: symplectic→行列積 / Cayley–Dickson→複素・四元・セデニオン
    ・乗算器なし（gate9 の 符号つき桁セル）・経路は 全部 XOR（k=i⊕j）

  軸1【SIMT・並列】  P 個のユニットを 並べ、各々 独立データを 1 クロックで（スループット ∝ P）
  軸2【接続・繋ぐ】  大 N×N 行列積を (N/b)³ 個の ブロック積に タイル化、accumulate 網で 繋ぐ
                    ＝ GPU の GEMM データフロー そのもの。ここでは **厳密に一致するか** を検証。
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
import numpy as np
from nd_algebra import ref_mult_M, sd2_mult_M, cd_omega
from matrix_algebra import symplectic_omega, coeff_to_mat, mat_to_coeff
from sd2_core import to_sd2, from_sd2


def unit_block_product(ca, cb, OM, M):
    """1 ユニット = 係数領域での twisted 積。symplectic OM なら b×b 行列積 A·B の係数。"""
    return ref_mult_M(ca, cb, OM, M)

def unit_block_product_sd2(ca, cb, OM, M, K=14):
    """同じ原子演算を 底2 符号つき桁の 乗算器フリー回路で（ハード相当）。"""
    return [from_sd2(w) for w in sd2_mult_M([to_sd2(v, K) for v in ca],
                                            [to_sd2(v, K) for v in cb], OM, M)]


def tiled_gemm(Acoef, Bcoef, T, OM, M, use_sd2=False):
    """T×T ブロックの 大行列積 C=A·B を ユニット群で。各 C_ij = Σ_k unit(A_ik, B_kj)。
       戻り: C の 係数ブロック[T][T], 使った ユニット積回数。"""
    unit = unit_block_product_sd2 if use_sd2 else unit_block_product
    Ccoef = [[None]*T for _ in range(T)]
    n_mults = 0
    for i in range(T):
        for j in range(T):
            acc = [0]*M
            for k in range(T):
                p = unit(Acoef[i][k], Bcoef[k][j], OM, M); n_mults += 1
                acc = [x+y for x, y in zip(acc, p)]      # accumulate 網（線形・安い）
            Ccoef[i][j] = acc
    return Ccoef, n_mults


def self_test():
    rng = np.random.default_rng(20260725)

    print("=" * 80)
    print("原子演算: 1 ユニット = 1 パスで b×b 行列積（テンソルコア相当）。まず 厳密性。")
    print("=" * 80)
    for q in (1, 2, 3):
        M = 1 << (2*q); b = 1 << q
        OM, elems, Bbasis = symplectic_omega(q)
        ok = ok_sd2 = True
        for _ in range(20):
            ca = [int(v) for v in rng.integers(-3, 4, M)]
            cb = [int(v) for v in rng.integers(-3, 4, M)]
            A = coeff_to_mat(ca, Bbasis); Bm = coeff_to_mat(cb, Bbasis)   # ユニットが表す b×b 行列
            cprod = unit_block_product(ca, cb, OM, M)
            if not np.array_equal(coeff_to_mat(cprod, Bbasis), A @ Bm): ok = False
            if unit_block_product_sd2(ca, cb, OM, M) != cprod: ok_sd2 = False
        print(f"  q={q}  M={M:>3}  b×b={b}×{b}  ユニット積==A·B: {'✓' if ok else '×'}   "
              f"乗算器フリー回路==代数: {'✓' if ok_sd2 else '×'}")

    print()
    print("=" * 80)
    print("軸2【繋ぐ】 大 N×N 行列積を ブロックに タイル化 → ユニット群 + accumulate で 厳密一致か")
    print("=" * 80)
    q = 2; M = 1 << (2*q); b = 1 << q                    # b=4 のユニット
    OM, elems, Bbasis = symplectic_omega(q)
    for T in (2, 3, 4):
        N = T * b
        # 大行列を T×T の 係数ブロックで（各ブロック = ユニットの native 型）
        Acoef = [[[int(v) for v in rng.integers(-2, 3, M)] for _ in range(T)] for _ in range(T)]
        Bcoef = [[[int(v) for v in rng.integers(-2, 3, M)] for _ in range(T)] for _ in range(T)]
        # 参照: 実際の N×N 行列を 組んで numpy で 直接積
        Afull = np.block([[coeff_to_mat(Acoef[i][k], Bbasis) for k in range(T)] for i in range(T)])
        Bfull = np.block([[coeff_to_mat(Bcoef[i][k], Bbasis) for k in range(T)] for i in range(T)])
        Cfull = Afull @ Bfull
        # 繋いだユニット群で
        Ccoef, n_mults = tiled_gemm(Acoef, Bcoef, T, OM, M)
        match = all(np.array_equal(coeff_to_mat(Ccoef[i][j], Bbasis),
                                   Cfull[i*b:(i+1)*b, j*b:(j+1)*b]) for i in range(T) for j in range(T))
        print(f"  N={N:>2}×{N:<2} = {T}×{T} ブロック(b={b})   ユニット積 {n_mults:>3} 回(=T³={T**3})   "
              f"大行列積と一致: {'✓' if match else '×'}")

    print()
    print("=" * 80)
    print("設計数字（1 ユニット = b×b 行列積 / 1 XOR パス）")
    print("=" * 80)
    print(f"  {'q':>2}{'M=2^2q':>8}{'b×b':>7}{'ユニット内 積セル(M²)':>18}{'gate9/パス(≈M²K²)':>18}")
    K = 8
    for q in (1, 2, 3, 4):
        M = 1 << (2*q); b = 1 << q
        print(f"  {q:>2}{M:>8}{f'{b}×{b}':>7}{M*M:>18}{M*M*K*K:>18}")
    print(f"""
  ⟹ **繋いだユニット群が 大行列積を 厳密に 出す**（T³ 個の ブロック積 + 線形 accumulate）。
     並列(SIMT): P ユニットで P 原子演算/クロック。接続(GEMM): 出力 T² ブロックは 独立 ⟹ 並列、
     内側 Σ_k は accumulate 網。＝ GPU の テンソルコア + タイル化 と 同じ データフロー。
  ⟹ 差別化: **同じ配線で OMEGA を 差し替えるだけ** で b×b 行列積 ↔ 複素/四元/セデニオン。
     乗算器なし（符号つき桁セル）。原子演算コストは M²=b⁴（＝直接 b×b 積 b³ の b 倍）＝
     速さでなく **一様な再構成性・乗算器フリー** が 買うもの（正直な トレードオフ）。""")


if __name__ == "__main__":
    self_test()
