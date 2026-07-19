#!/usr/bin/env python3
"""**大きさのビットを二つに割る**（利用者）: le = 「|x| ≤ v」／ ge = 「|x| ≥ v」。

    le  ge
     1   1  → **確定**（|x| = v）
     1   0  → **|x| ≤ v**   ← **MIN は v=MIN の特別な場合にすぎない**
     0   1  → **|x| ≥ v**   ← **MAX は v=MAX の場合**
     0   0  → **境界なし**

問い: ① 健全か（真値が必ず境界の中に居るか） ② 今より情報を持つか
"""
import numpy as np
from fractions import Fraction as Fr
from total_arith import MIN, MAX, saturate

class B:
    """(値, le, ge)。|x| についての主張。符号は別に持つ（ここでは大きさだけ見る）。"""
    __slots__ = ("v", "le", "ge")
    def __init__(self, v, le=True, ge=True): self.v, self.le, self.ge = float(v), le, ge
    def __repr__(self):
        if self.le and self.ge: return f"={self.v:.4g}"
        if self.le: return f"**≤{self.v:.4g}**"
        if self.ge: return f"**≥{self.v:.4g}**"
        return "**境界なし**"

def sat_mag(v):
    """溢れ→(MAX, ge) ／ 潰れ→(MIN, le) ／ それ以外→確定。**大きさだけ。**"""
    with np.errstate(all="ignore"):
        v = float(np.abs(v))
        if np.isinf(v) or v > MAX: return B(MAX, le=False, ge=True)     # **|x| ≥ MAX**
        if v == 0.0:               return B(MIN, le=True,  ge=False)     # **|x| ≤ MIN**
        return B(v)

def mul(x, y):
    with np.errstate(all="ignore"):
        v = float(np.float64(x.v) * np.float64(y.v))
    s = sat_mag(v)
    # |ab| ≤ ... : 両方に上界があれば積が上界。|ab| ≥ ... : 両方に下界があれば積が下界
    le = x.le and y.le and s.le
    ge = x.ge and y.ge and s.ge
    if not (le or ge): return B(s.v, False, False)
    return B(s.v, le, ge)

def add(x, y, same_sign):
    with np.errstate(all="ignore"):
        v = float(np.float64(x.v) + np.float64(y.v))
    s = sat_mag(v)
    le = x.le and y.le and s.le                      # |a+b| ≤ |a|+|b|（三角不等式。常に成り立つ）
    ge = same_sign and x.ge and y.ge and s.ge        # **同符号のときだけ** |a+b| ≥ |a|+|b|
    if not (le or ge): return B(s.v, False, False)
    return B(s.v, le, ge)

def sqrt(x):
    with np.errstate(all="ignore"):
        v = float(np.sqrt(np.float64(x.v)))
    s = sat_mag(v)
    return B(s.v, x.le and s.le, x.ge and s.ge)      # **単調 ⟹ 境界の向きが保たれる**

print("=" * 92)
print("① **健全か** — 真値が、主張した境界の中に居るか")
print("=" * 92)
def holds(b, truth):
    """**有理数のまま比べる**（1e600 は float にできない）。"""
    t = abs(Fr(truth))
    bv = Fr(b.v)
    if b.le and b.ge: return abs(t - bv) <= abs(bv)/10**9
    if b.le: return t <= bv * Fr(1000000001, 1000000000)
    if b.ge: return t >= bv * Fr(999999999, 1000000000)
    return True
print(f"  {'式':<32}{'真値':>14}{'主張':>22}{'健全':>8}")
ok = n = 0
tests = [
    ("1e−200 × 1e−200",  mul(sat_mag(1e-200), sat_mag(1e-200)),  Fr(10)**-400),
    ("1e300 × 1e300",    mul(sat_mag(1e300),  sat_mag(1e300)),   Fr(10)**600),
    ("2.0 × 3.0",        mul(sat_mag(2.0),    sat_mag(3.0)),     Fr(6)),
    ("**√(1e−200 × 1e−200)**", sqrt(mul(sat_mag(1e-200), sat_mag(1e-200))), Fr(10)**-200),
    ("**(≥MAX) × (≤MIN)**", mul(sat_mag(1e400), sat_mag(1e-400)), Fr(1)),
    ("(≤MIN) + (≤MIN)",  add(sat_mag(1e-400), sat_mag(1e-400), True), Fr(2)*Fr(10)**-400),
]
for tag, b, truth in tests:
    h = holds(b, truth); ok += h; n += 1
    e = len(str(Fr(truth).denominator)) - len(str(Fr(truth).numerator))
    ft = f"~1e{-e:+d}" if abs(e) > 20 else f"{float(truth):.4g}"
    print(f"  {tag:<32}{ft:>14}{repr(b):>22}{('✓' if h else '**✗**'):>8}")
print(f"\n  健全: **{ok}/{n}**")

print()
print("=" * 92)
print("② **今より情報を持つか** — 監査人の1ビット版と並べる")
print("=" * 92)
print(f"  {'式':<32}{'1ビット版（監査人）':>26}{'**2ビット版（利用者）**':>28}")
rows = [
    ("1e−200 × 1e−200",           "MIN（= 大きさ不明）",  mul(sat_mag(1e-200), sat_mag(1e-200))),
    ("**√(1e−200 × 1e−200)**",    "**2.22e−162 + 不明**", sqrt(mul(sat_mag(1e-200), sat_mag(1e-200)))),
    ("1e300 × 1e300",             "MAX（= 大きさ不明）",  mul(sat_mag(1e300), sat_mag(1e300))),
    ("**(≥MAX) × (≤MIN)**",       "**MIN or MAX（嘘）**", mul(sat_mag(1e400), sat_mag(1e-400))),
]
for tag, old, b in rows:
    print(f"  {tag:<32}{old:>26}{repr(b):>28}")
print("""
  ⟹ **√ の行が要点**: 1ビット版は「2.22e−162 だが大きさ不明」＝ **何も分からない**と言う。
     2ビット版は「**|x| ≤ 2.22e−162**」と言う。**真値 1e−200 は、その中に居る。**
     **1ビット版が捨てていた情報を、2ビット版は持っている。**

  ⟹ そして **(≥MAX) × (≤MIN) は 2ビット版でだけ正しく「境界なし」と言える。**
     1ビット版は MIN か MAX のどちらかの値を返すので、**在りもしない境界を主張してしまう。**""")
