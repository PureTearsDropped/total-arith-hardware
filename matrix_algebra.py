#!/usr/bin/env python3
"""同じ XOR 配線のハードで、**符号表（twist）を変えると 行列積になる**（利用者の洞察）。

  Cayley–Dickson twist → セデニオン等（非結合・零因子）        ← nd_algebra.py
  **symplectic (Pauli/clock-shift) twist → 2^q × 2^q 行列代数**  ← ここ

  基底 e_g（g ∈ (Z/2)^{2q}）を パウリ演算子 ⊗_i (X^{a_i} Z^{b_i}) に対応させると、
  twisted 積 e_g·e_h = ±e_{g⊕h} が **行列の積** になる（符号 = symplectic cocycle）。
  ⟹ M=2^{2q} 成分 = 2^q × 2^q 行列。M を変えると 行列サイズが変わる。
  そして 経路は 同じ XOR、変わるのは 符号表だけ ⟹ **同じ乗算器で 行列積**。
"""
import sys, os
from itertools import product as iproduct
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
import numpy as np
from sd2_core import to_sd2, from_sd2
from nd_algebra import sd2_mult_M

X = np.array([[0, 1], [1, 0]]); Z = np.array([[1, 0], [0, -1]])
def pauli(a, b):
    return np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b)   # X^a Z^b

def basis_matrix(g, q):
    """g=(a0,b0,...,a_{q-1},b_{q-1}) → ⊗_i X^{a_i} Z^{b_i}（2^q × 2^q）。"""
    Mm = np.array([[1]])
    for i in range(q):
        Mm = np.kron(Mm, pauli(g[2*i], g[2*i+1]))
    return Mm

def symplectic_omega(q):
    """基底行列から 直接 OMEGA を導く: e_g·e_h = OMEGA[g,h]·e_{g⊕h}。"""
    d = 1 << (2*q); D = 1 << q
    elems = [tuple(t) for t in iproduct((0, 1), repeat=2*q)]
    idx = {g: i for i, g in enumerate(elems)}
    B = [basis_matrix(g, q) for g in elems]
    OM = np.zeros((d, d), dtype=int)
    for gi, g in enumerate(elems):
        for hi, h in enumerate(elems):
            prod = B[gi] @ B[hi]
            k = idx[tuple((a ^ b) for a, b in zip(g, h))]      # 経路 = XOR
            # prod = s · B[k] の s を 取り出す（B[k] の 非零成分で 比較）
            r, c = np.unravel_index(np.argmax(np.abs(B[k])), B[k].shape)
            OM[gi, hi] = int(round(prod[r, c] / B[k][r, c]))
            assert np.array_equal(prod, OM[gi, hi] * B[k]), "XOR/符号 で 表せない"
    return OM, elems, B

def mat_to_coeff(A, elems, B):
    """行列 A → 基底係数（trace 直交性: <e_g,e_h> = 2^q δ）。"""
    D = A.shape[0]
    return [int(round(np.trace(B[i].T @ A) / D)) for i in range(len(elems))]

def coeff_to_mat(c, B):
    return sum(c[i] * B[i] for i in range(len(c)))


def self_test():
    print("=" * 76)
    print("符号表を symplectic にすると 同じ XOR ハードが 2^q×2^q 行列積 になる")
    print("=" * 76)
    print(f"  {'q':>3}{'M=2^2q':>8}{'行列':>10}{'結合的':>8}{'e_g·e_h=行列積':>16}{'sd2積で行列積':>16}")
    rng = np.random.default_rng(0)
    for q in (1, 2, 3):
        M = 1 << (2*q); D = 1 << q
        OM, elems, B = symplectic_omega(q)          # XOR/符号で 表せることは assert 済み
        # 結合性（cocycle ⟹ 結合的、セデニオンと違う）
        from nd_algebra import ref_mult_M
        assoc = True
        for _ in range(30):
            a=[int(v) for v in rng.integers(-3,3,M)]; b=[int(v) for v in rng.integers(-3,3,M)]
            c=[int(v) for v in rng.integers(-3,3,M)]
            if ref_mult_M(ref_mult_M(a,b,OM,M),c,OM,M)!=ref_mult_M(a,ref_mult_M(b,c,OM,M),OM,M):
                assoc=False
        # 代数側（整数の係数ベクトル）から始める ⟹ 丸め不要。
        # A=Σc_g e_g, B=Σd_h e_h の 行列積 A·B が、係数の twisted 積 で 出るか。
        ok_alg = ok_sd2 = True
        for _ in range(20):
            ca=[int(v) for v in rng.integers(-4,4,M)]; cb=[int(v) for v in rng.integers(-4,4,M)]
            A=coeff_to_mat(ca,B); Bm=coeff_to_mat(cb,B)
            prod_coeff=ref_mult_M(ca,cb,OM,M)
            if not np.array_equal(coeff_to_mat(prod_coeff,B), A@Bm): ok_alg=False
            # 同じことを 底2 符号つき桁の 乗算器で
            got=[from_sd2(w) for w in sd2_mult_M([to_sd2(v,12) for v in ca],
                                                 [to_sd2(v,12) for v in cb], OM, M)]
            if got!=prod_coeff: ok_sd2=False
        print(f"  {q:>3}{M:>8}{f'{D}×{D}':>10}{('✓' if assoc else '×'):>8}"
              f"{('✓' if ok_alg else '×'):>16}{('✓' if ok_sd2 else '×'):>16}")
    print("""
  ⟹ **利用者の通り。M=2^{2q} + symplectic 符号表 で、2^q × 2^q 行列積が そのまま出る。**
     q=1 → 2×2、q=2 → 4×4、q=3 → 8×8。M を変えると 行列サイズが 変わる。
     経路は 全部 同じ XOR（k=i⊕j）、変わるのは 符号表 OMEGA だけ ⟹ **同じ乗算器**。
     結合的（cocycle）＝ セデニオン(非結合)と 違う代数が、**同じハードから twist で 出る**。

  ⟹ **一つの XOR 配線ハードが:**
       Cayley–Dickson twist → 実/複素/四元/八元/セデニオン（除算代数の塔）
       symplectic twist     → 2×2, 4×4, 8×8, … 行列積
     **twist（小さい符号表）を選ぶだけで 用途が変わる。** ＝ 色んな行列積・色んな代数。""")


if __name__ == "__main__":
    self_test()
