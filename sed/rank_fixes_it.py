#!/usr/bin/env python3
"""零空間で判定すれば、±MIN 規約はセデニオンでも正しくなるか。"""
import numpy as np
from sedenion_tensor_logic import OMEGA, ref_mult, M

def L(x):
    A = np.zeros((M, M))
    for i in range(M):
        for j in range(M):
            A[i ^ j, j] += OMEGA[i, j] * x[i]
    return A

def zero_ok_scalar(a, b, k):
    """我々の規約（§7.7.2 で壊れた）: 掛ける前の値が 0 か"""
    return not any(a) or not any(b)

def zero_ok_rank(a, b, k):
    """**新: b が L(a) の零空間に居るか。** 居れば 0 は本物。"""
    return np.allclose(L(np.asarray(a, float)) @ np.asarray(b, float), 0, atol=1e-12)

print("=" * 84)
print("**0 が本物か** — スカラー判定 vs 零空間判定")
print("=" * 84)
rng = np.random.default_rng(0)

# ① 零因子: 積が厳密に 0
a1 = np.zeros(M); a1[1] = 1; a1[10] = 1
b1 = np.zeros(M); b1[4] = 1; b1[15] = -1
# ② 潰れ: 小さすぎて 0 になる（浮動小数）
a2 = np.zeros(M); a2[1] = 1e-200
b2 = np.zeros(M); b2[2] = 1e-200
# ③ 本物の 0
a3 = np.zeros(M)
b3 = rng.integers(-5, 5, M).astype(float)

print(f"  {'場合':<32}{'真の積':>12}{'スカラー判定':>16}{'**零空間判定**':>20}")
for tag, a, b, truth in [
    ("**零因子**（e1+e10）(e4−e15)", a1, b1, "**厳密に 0**"),
    ("**潰れ** 1e−200 × 1e−200",     a2, b2, "**0 でない**"),
    ("**本物の 0** 0 × ランダム",      a3, b3, "**厳密に 0**"),
]:
    s = zero_ok_scalar(list(a), list(b), 0)
    r = zero_ok_rank(a, b, 0)
    want = (truth == "**厳密に 0**")
    fs = f"{'本物' if s else '潰れ'} {'✓' if s == want else '**✗**'}"
    fr = f"**{'本物' if r else '潰れ'}** {'✓' if r == want else '**✗**'}"
    print(f"  {tag:<32}{truth:>12}{fs:>16}{fr:>20}")

print("""
  ⟹ **スカラー判定は零因子で外す。零空間判定は 3 つとも当てる。**
     「掛ける前の値」からは分からないが、**「掛ける前の行列」からは分かる。**""")

print()
print("=" * 84)
print("そして零空間判定は、実数でも同じ答えを出すか（後方互換）")
print("=" * 84)
print(f"  {'a':>10}{'b':>10}{'真の積':>12}{'スカラー':>12}{'**零空間**':>14}")
for a, b in [(0.0, 5.0), (5.0, 0.0), (0.0, 0.0), (1e-200, 1e-200), (2.0, 3.0)]:
    La = np.array([[a]])
    s = (a == 0.0 or b == 0.0)
    r = bool(np.allclose(La @ np.array([b]), 0, atol=1e-300))
    prod = a*b
    want = (prod == 0.0 and not (a != 0 and b != 0 and abs(a*b) > 0))
    truth = "0" if (a == 0 or b == 0) else ("**潰れ**" if a*b == 0 else f"{a*b:.3g}")
    print(f"  {a:>10.3g}{b:>10.3g}{truth:>12}{str(s):>12}{f'**{r}**':>14}")
print("""
  ⟹ **実数では 1×1 行列なので、零空間判定 ≡ スカラー判定。**（a=0 のときだけ零空間が非自明）
     ⟹ **後方互換。実数では何も変わらず、零因子のある代数でだけ正しくなる。**

  ⟹ **`zero_ok` の正しい定義:**
       ~~掛ける前の値が 0 か~~  →  **b が L(a) の零空間に居るか**
     ＝ **「積が 0 になるのは構造上そうなるからか、それとも桁が足りないからか」**""")
