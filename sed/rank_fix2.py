#!/usr/bin/env python3
"""直す: 積の値を見ず、L(a) の **構造**（正規化してから rank / 零空間）だけを見る。

前版の誤り: `L(a) @ b` を浮動小数で評価したので、1e−400 が 0 に潰れて
            「零空間に居る」と誤診した。**判定そのものが同じ病気にかかっていた。**
直し:       rank は **尺度不変**（rank(λA) = rank(A)）。**a と b を正規化してから**
            零空間の所属を見れば、桁の大小は判定に入らない。
"""
import numpy as np
from scipy.linalg import null_space
from sedenion_tensor_logic import OMEGA, M

def L(x):
    A = np.zeros((M, M))
    for i in range(M):
        for j in range(M):
            A[i ^ j, j] += OMEGA[i, j] * x[i]
    return A

def zero_ok_scalar(a, b):
    return not np.any(a) or not np.any(b)

def zero_ok_struct(a, b):
    """**構造だけを見る。** a, b を正規化 ⟹ 桁の大小は判定に入らない。"""
    a, b = np.asarray(a, float), np.asarray(b, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0: return True                # 本物の 0
    A = L(a / na)                                     # **正規化: rank は不変**
    bh = b / nb
    return bool(np.linalg.norm(A @ bh) < 1e-10)       # bh が零空間に居るか

rng = np.random.default_rng(0)
a1 = np.zeros(M); a1[1] = 1; a1[10] = 1
b1 = np.zeros(M); b1[4] = 1; b1[15] = -1
a2 = np.zeros(M); a2[1] = 1e-200
b2 = np.zeros(M); b2[2] = 1e-200
a3 = np.zeros(M)
b3 = rng.integers(-5, 5, M).astype(float)
a4 = np.zeros(M); a4[1] = 1e-200; a4[10] = 1e-200    # **零因子だが、桁も小さい**
b4 = np.zeros(M); b4[4] = 1e-200; b4[15] = -1e-200

print("=" * 88)
print("**0 が本物か** — 正規化してから構造を見る")
print("=" * 88)
print(f"  {'場合':<34}{'真の積':>14}{'スカラー':>12}{'**構造**':>14}")
for tag, a, b, want in [
    ("零因子 (e1+e10)(e4−e15)",          a1, b1, True),
    ("潰れ 1e−200 × 1e−200",             a2, b2, False),
    ("本物の 0  0 × ランダム",              a3, b3, True),
    ("**零因子 かつ 桁も小さい**",            a4, b4, True),
]:
    s, r = zero_ok_scalar(a, b), zero_ok_struct(a, b)
    truth = "**厳密に 0**" if want else "**0 でない**"
    print(f"  {tag:<34}{truth:>14}"
          f"{f'{chr(26412)+chr(29289) if s else chr(28291)+chr(12428)} {chr(10003) if s==want else chr(10007)}':>12}"
          f"{f'**{chr(26412)+chr(29289) if r else chr(28291)+chr(12428)}** {chr(10003) if r==want else chr(10007)}':>14}")

print()
print("=" * 88)
print("実数での後方互換（1×1 行列）")
print("=" * 88)
print(f"  {'a':>10}{'b':>10}{'真の積':>14}{'スカラー':>12}{'**構造**':>14}")
def zero_ok_struct_1d(a, b):
    if a == 0.0 or b == 0.0: return True
    return bool(abs((a/abs(a)) * (b/abs(b))) < 1e-10)     # 正規化 ⟹ |±1 × ±1| = 1 ≠ 0
for a, b in [(0.0, 5.0), (5.0, 0.0), (0.0, 0.0), (1e-200, 1e-200), (2.0, 3.0)]:
    s = (a == 0.0 or b == 0.0)
    r = zero_ok_struct_1d(a, b)
    truth = "**厳密に 0**" if (a == 0 or b == 0) else ("**潰れ**" if a*b == 0 else f"{a*b:.3g}")
    print(f"  {a:>10.3g}{b:>10.3g}{truth:>14}{str(s):>12}{f'**{r}**':>14}")
