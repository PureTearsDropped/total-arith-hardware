#!/usr/bin/env python3
"""セデニオンを 16×16 行列にする。零因子 = 逆行列が無いもの。擬似逆行列は在る。"""
import numpy as np
from sedenion_tensor_logic import OMEGA, ref_mult, M
np.set_printoptions(precision=3, suppress=True)

def L(x):
    """x を左から掛ける 16×16 行列。 (L(x) @ y) == x*y"""
    A = np.zeros((M, M))
    for i in range(M):
        for j in range(M):
            A[i ^ j, j] += OMEGA[i, j] * x[i]
    return A

print("=" * 80)
print("① 定義の確認: L(x) @ y == x*y")
print("=" * 80)
rng = np.random.default_rng(0)
x = rng.integers(-5, 5, M).astype(float)
y = rng.integers(-5, 5, M).astype(float)
print(f"  L(x) @ y == ref_mult(x,y) : **{np.allclose(L(x) @ y, ref_mult(list(x), list(y)))}**")

print()
print("=" * 80)
print("② 零因子の行列は、特異か")
print("=" * 80)
a = np.zeros(M); a[1] = 1; a[10] = 1        # e1 + e10
b = np.zeros(M); b[4] = 1; b[15] = -1       # e4 - e15
cases = [("**e1 + e10**（零因子）", a), ("**e4 − e15**（零因子）", b),
         ("e1（単位）", np.eye(M)[1]), ("ランダム", rng.integers(-5,5,M).astype(float))]
print(f"  {'x':<24}{'rank L(x)':>12}{'det':>16}{'逆行列':>12}")
for tag, v in cases:
    A = L(v)
    r = np.linalg.matrix_rank(A)
    d = np.linalg.det(A)
    print(f"  {tag:<24}{r:>12}{d:>16.4g}{'**無い**' if r < M else 'ある':>12}")

print()
print("=" * 80)
print("③ 逆行列が無い ⟺ 零因子。何を掛けると 0 になるか = 零空間")
print("=" * 80)
from scipy.linalg import null_space
A = L(a)
N = null_space(A)
print(f"  L(e1+e10) の零空間の次元: **{N.shape[1]}**  （16 − rank = {M - np.linalg.matrix_rank(A)}）")
print(f"  b = e4 − e15 は、その零空間に居るか: **{np.allclose(A @ b, 0)}**")
print(f"  ⟹ **(e1+e10)(e4−e15) = 0 は「b が L(a) の零空間に居る」ということ。**")

print()
print("=" * 80)
print("④ 擬似逆行列は、在る")
print("=" * 80)
P = np.linalg.pinv(A)
print(f"  pinv(L(e1+e10)) は存在するか      : **True**（擬似逆行列は必ず在る）")
print(f"  それは 0 か                      : **{np.allclose(P, 0)}**")
print(f"  Penrose (1) A P A == A          : **{np.allclose(A@P@A, A)}**")
print(f"  Penrose (2) P A P == P          : **{np.allclose(P@A@P, P)}**")
print(f"  rank(pinv) == rank(A)           : **{np.linalg.matrix_rank(P) == np.linalg.matrix_rank(A)}**")

print()
print("=" * 80)
print("⑤ そして 実数の a/0 と、同じ形")
print("=" * 80)
print(f"""  実数で 0 を 1×1 行列にすると:  L(0) = [[0]]、rank **0**、逆行列 **無い**
  セデニオンの零因子 e1+e10:      L(a)  = 16×16、rank **{np.linalg.matrix_rank(L(a))}**、逆行列 **無い**

  ⟹ **「逆元が無い」は、どちらも「L が特異」である。同じ現象。**
     実数の 0 は、**rank 0 という極端な場合**にすぎない。

  ⟹ そして **擬似逆行列は両方に在る**:
       pinv([[0]])   = [[0]]        ⟹ **x/0 = 0**（§7.6.7⑪）
       pinv(L(a))    = 0 でない     ⟹ **零因子で「割る」ことに、答えが在る**
""")
v = rng.integers(-5, 5, M).astype(float)
sol = P @ v
print(f"  例: 「(e1+e10) · z = v を解け」の最小ノルム最小二乗解 z = pinv(L(a)) @ v")
print(f"      残差 ‖L(a)z − v‖ = **{np.linalg.norm(A@sol - v):.4f}**  （0 でない = 厳密解は無い）")
print(f"      ‖z‖              = **{np.linalg.norm(sol):.4f}**  （その中で最小ノルム）")
