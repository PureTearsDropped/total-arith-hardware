#!/usr/bin/env python3
"""**符号不明になれるのは、和と差だけか**（利用者）。総当たりで数える。

  §7.6.4.1 の証明: 積・商は符号を「符号の積」で決める ⟹ 大きさを失っても符号は残る。
                   加減だけが符号を「大きさの比較」で決める ⟹ 大きさを失うと符号も失う。

  ここでは **証明ではなく、数える**。飽和した入力で、各演算が符号を決められない回数。
"""
import warnings, numpy as np
from fractions import Fraction as Fr
warnings.filterwarnings("ignore")
MAX = np.finfo(np.float64).max
MIN = np.nextafter(0.0, 1.0)

def sgn(x): return 0 if x == 0 else (1 if x > 0 else -1)

# 飽和した値 = 「真値は範囲のどこか」。**その範囲を実際に走査して、符号が一意に決まるか見る。**
def candidates(v):
    """値 v が表す真値の候補（有理数）。MAX/MIN は境界なので複数。"""
    if abs(v) == MAX:   return [Fr(v)*k for k in (1, 10, 10**6, 10**60)]      # **|x| ≥ MAX**
    if abs(v) == MIN:   return [Fr(v)/k for k in (1, 10, 10**6, 10**60)]      # **|x| ≤ MIN**
    return [Fr(v)]                                                             # 確定

def sign_determined(op, a, b):
    """a, b が表す真値の**全ての組**で、結果の符号が一つに定まるか。"""
    signs = set()
    for A in candidates(a):
        for B in candidates(b):
            if op == "/" and B == 0: return None
            r = {"+": A+B, "-": A-B, "*": A*B, "/": (A/B if B else None)}[op]
            signs.add(0 if r == 0 else (1 if r > 0 else -1))
    return len(signs) == 1

print("=" * 84)
print("**符号が決まらないのは、どの演算か** — 飽和した入力の真値を走査して数える")
print("=" * 84)
vals = [MAX, -MAX, MIN, -MIN, 1.0, -1.0, 5.0, -5.0]
tab = {op: {"unk": 0, "n": 0} for op in "+-*/"}
for op in "+-*/":
    for a in vals:
        for b in vals:
            d = sign_determined(op, a, b)
            if d is None: continue
            tab[op]["n"] += 1
            tab[op]["unk"] += (not d)
print(f"  {'演算':<8}{'組':>8}{'**符号が決まらない**':>24}")
for op, name in (("+", "**和**"), ("-", "**差**"), ("*", "**積**"), ("/", "**商**")):
    t = tab[op]
    fu = f"**{t['unk']}**  ({100*t['unk']/t['n']:.0f}%)" if t['unk'] else "**0 ← 起きない**"
    print(f"  {name:<8}{t['n']:>8}{fu:>24}")
print("""
  ⟹ **利用者の通り: 符号不明になれるのは 和と差だけ。積と商では 0 回。**
     （証明は §7.6.4.1。ここでは 飽和した値の真値を実際に走査して、数で確認した）""")

print()
print("=" * 84)
print("**表が閉じる** — それぞれの演算が、ちょうど一つずつ、別のものを失う")
print("=" * 84)
print(f"""  {'':16}{'**符号を失う**':>18}{'**大きさを下へ失う**':>22}
  {'**和・差**':16}{'**できる**':>18}{'できない':>22}
  {'**積・商**':16}{'できない':>18}{'**できる**':>22}

  ⟹ **完全に相補的。** 和差は符号を、積商は大きさを失う。**両方を失う演算は無い。**

  ⟹ **回路として、二つの経路が交差しない:**
       **符号トリットが 0 になる** ← **加減の側だけ**
       **境界トリットが −1（≤）になる** ← **積商の側だけ**
       **境界トリットが +1（≥）になる** ← 四つとも（溢れは全演算で起きる）

  ⟹ そして 二つは **既にハードウェアにある**:
       **≥ / ≤** = MXCSR の **OE / UE**（1985 年から）
       **符号不明** = ...**無い**。IEEE はここだけ NaN で潰している（§7.6.6）""")
