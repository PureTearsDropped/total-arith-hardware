#!/usr/bin/env python3
"""直す(2): ノルムが二乗して潰れていた。**L∞（max|x|）は二乗しない。**"""
import numpy as np
from sedenion_tensor_logic import OMEGA, M

def L(x):
    A = np.zeros((M, M))
    for i in range(M):
        for j in range(M):
            A[i ^ j, j] += OMEGA[i, j] * x[i]
    return A

linf = lambda v: float(np.max(np.abs(v)))          # **二乗しない**

def zero_ok_scalar(a, b):
    return not np.any(a) or not np.any(b)

def zero_ok_struct(a, b):
    """構造だけを見る。正規化は **L∞**（二乗しないので潰れない）。"""
    a, b = np.asarray(a, float), np.asarray(b, float)
    na, nb = linf(a), linf(b)
    if na == 0.0 or nb == 0.0: return True                # 本物の 0
    A = L(a / na)                                          # rank は尺度不変
    return bool(linf(A @ (b / nb)) < 1e-10)

rng = np.random.default_rng(0)
cases = []
a = np.zeros(M); a[1] = 1; a[10] = 1
b = np.zeros(M); b[4] = 1; b[15] = -1
cases.append(("零因子 (e1+e10)(e4−e15)", a, b, True))
a = np.zeros(M); a[1] = 1e-200
b = np.zeros(M); b[2] = 1e-200
cases.append(("潰れ 1e−200·e1 × 1e−200·e2", a, b, False))
cases.append(("本物の 0  0 × ランダム", np.zeros(M), rng.integers(-5,5,M).astype(float), True))
a = np.zeros(M); a[1] = 1e-200; a[10] = 1e-200
b = np.zeros(M); b[4] = 1e-200; b[15] = -1e-200
cases.append(("**零因子 かつ 桁も小さい**", a, b, True))
a = np.zeros(M); a[1] = 1e200; a[10] = 1e200
b = np.zeros(M); b[4] = 1e200; b[15] = -1e200
cases.append(("**零因子 かつ 桁が大きい**", a, b, True))
cases.append(("普通の積", rng.integers(-5,5,M).astype(float), rng.integers(-5,5,M).astype(float), False))

print("=" * 92)
print("**0 が本物か** — L∞ で正規化してから構造を見る")
print("=" * 92)
print(f"  {'場合':<34}{'真の積':>14}{'スカラー':>14}{'**構造(L∞)**':>18}")
ns = nr = 0
for tag, a, b, want in cases:
    s, r = zero_ok_scalar(a, b), zero_ok_struct(a, b)
    ns += (s == want); nr += (r == want)
    truth = "**厳密に 0**" if want else "**0 でない**"
    fs = ("本物" if s else "潰れ") + (" ✓" if s == want else " **✗**")
    fr = "**" + ("本物" if r else "潰れ") + "**" + (" ✓" if r == want else " **✗**")
    print(f"  {tag:<34}{truth:>14}{fs:>14}{fr:>18}")
print(f"\n  正解数:  スカラー **{ns}/{len(cases)}**   ／   **構造(L∞) {nr}/{len(cases)}**")

print()
print("=" * 92)
print("実数での後方互換（1×1）")
print("=" * 92)
print(f"  {'a':>10}{'b':>10}{'真の積':>14}{'スカラー':>12}{'**構造**':>14}")
for a, b in [(0.0, 5.0), (5.0, 0.0), (0.0, 0.0), (1e-200, 1e-200), (2.0, 3.0), (1e200, 1e200)]:
    s = (a == 0.0 or b == 0.0)
    r = zero_ok_struct(np.array([a]+[0.0]*15), np.array([b]+[0.0]*15))
    truth = "**厳密に 0**" if (a == 0 or b == 0) else ("**潰れ**" if a*b == 0 else
             "**溢れ**" if not np.isfinite(a*b) else f"{a*b:.3g}")
    print(f"  {a:>10.3g}{b:>10.3g}{truth:>14}{str(s):>12}{f'**{r}**':>14}")
