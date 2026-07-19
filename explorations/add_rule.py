#!/usr/bin/env python3
"""**加減の境界は、符号で変わる**（利用者）。§7.7.3.4 の「MIN−MIN → ≤2·MIN」を疑う。

  |a| ≤ p, |b| ≤ q のとき:
      **同符号**: a, b ∈ (0,p]×(0,q]  ⟹ a+b ∈ (0, p+q]     ⟹ **≤ p+q**、符号は確定
      **異符号**: a ∈ (0,p], −b ∈ [−q,0) ⟹ a−b ∈ [−q, p)   ⟹ **≤ max(p,q)** ← **締まる**、符号は不明

  ⟹ **三角不等式 |a+b| ≤ |a|+|b| が締まるのは、同符号のときだけ。**
  ⟹ **符号トリットが、境界の式を選ぶ。二つは 加減では 独立でない。**
"""
from fractions import Fraction as Fr
import numpy as np
MAX = np.finfo(np.float64).max
MIN = np.nextafter(0.0, 1.0)

def scan(p, q, same_sign, is_upper):
    """真値を走査して、|a±b| の実際の上界（or 下界）と符号を求める。"""
    P, Q = Fr(p), Fr(q)
    if is_upper:  # |a| ≤ p, |b| ≤ q
        A = [P, P/2, P/1000, P/10**60]
        B = [Q, Q/2, Q/1000, Q/10**60]
    else:         # |a| ≥ p, |b| ≥ q
        A = [P, P*2, P*1000, P*10**60]
        B = [Q, Q*2, Q*1000, Q*10**60]
    mags, signs = [], set()
    for a in A:
        for b in B:
            r = a + b if same_sign else a - b
            mags.append(abs(r)); signs.add(0 if r == 0 else (1 if r > 0 else -1))
    return (max(mags) if is_upper else min(mags)), (len(signs) == 1)

print("=" * 92)
print("**加減の境界は、符号で変わるか** — 真値を走査して実測")
print("=" * 92)
print(f"  {'式':<30}{'実際の境界':>24}{'監査人が書いた':>20}{'符号':>12}")
rows = []
def ratio(m, unit):
    """**MIN / MAX の何倍か**で表す（1e400 は float にできない）。"""
    r = Fr(m) / Fr(unit)
    return f"**{r}·{'MIN' if unit == MIN else 'MAX'}**" if r.denominator == 1 else f"**{float(r):.4g}×**"
# --- 上端（MIN、|x| ≤ MIN）---
m, s = scan(MIN, MIN, True, True)
rows.append(("**(+MIN) + (+MIN)**", f"≤ {ratio(m, MIN)}", "≤ 2·MIN", s))
m, s = scan(MIN, MIN, False, True)
rows.append(("**(+MIN) − (+MIN)**", f"≤ {ratio(m, MIN)}", "**≤ 2·MIN ← 2倍緩い**", s))
# --- 下端（MAX、|x| ≥ MAX）---
m, s = scan(MAX, MAX, True, False)
rows.append(("**(+MAX) + (+MAX)**", f"≥ {ratio(m, MAX)}", "≥ 2·MAX", s))
m, s = scan(MAX, MAX, False, False)
mn, s2 = scan(MAX, MAX, False, False)
rows.append(("**(+MAX) − (+MAX)**", f"≥ {ratio(mn, MAX)} ← **0。境界なし**", "境界なし", s2))
for tag, got, said, s in rows:
    print(f"  {tag:<30}{got:>24}{said:>20}{('確定' if s else '**不明**'):>12}")

print("""
  ⟹ **利用者の通り。`(+MIN) − (+MIN)` の境界は MIN であって 2·MIN ではない。**
     監査人の §7.7.3.4 は **足し算の境界を引き算に当てていた。2 倍緩い。**

  ⟹ そして **MAX 側では逆**: 引き算は境界を失う（`境界なし`）。**締まるどころか消える。**""")

print()
print("=" * 92)
print("**加減の規則（訂正版）** — 符号トリットが、境界の式を選ぶ")
print("=" * 92)
print("""  |a| ≤ p, |b| ≤ q（**上端どうし**）:
      **同符号** → **≤ p + q**       符号は確定       （三角不等式。締まるのはここだけ）
      **異符号** → **≤ max(p, q)**   **符号は不明**   ← **足すのでなく、大きい方**

  |a| ≥ p, |b| ≥ q（**下端どうし**）:
      **同符号** → **≥ p + q**       符号は確定
      **異符号** → **境界なし**       **符号は不明**   ← 真値が 0 を跨ぐ

  片方が確定（|a| = v）:
      **≤ q と足す** → |a ± b| ≤ v + q、かつ ≥ v − q  ⟹ **確定に近い**（両側の境界）
      **≥ q と足す** → **≥ q − v**（同符号なら ≥ q + v）

  ⟹ **符号トリットと境界トリットは、加減では独立でない。**
     **符号が、どの式を使うかを決める。**（積商では独立に流れた。§7.7.4.3）
     **⟹ 加減のゲートは、二つのトリットを両方読む。ここが唯一の交差点。**""")
