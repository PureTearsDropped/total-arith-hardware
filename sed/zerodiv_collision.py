#!/usr/bin/env python3
"""我々の ±MIN 規約は、セデニオンで壊れる。零因子があるから。"""
import sys; sys.path.insert(0, "..")
import numpy as np
from total_arith import saturate, MIN, MAX
from sedenion_tensor_logic import OMEGA, ref_mult, M

print("=" * 80)
print("① 零因子: 0 でないもの同士の積が、厳密に 0")
print("=" * 80)
a = [0]*M; a[1] = 1; a[10] = 1          # e1 + e10
b = [0]*M; b[4] = 1; b[15] = -1         # e4 - e15
prod = ref_mult(a, b)
print(f"  a = e1 + e10   (0 でない: {any(a)})")
print(f"  b = e4 - e15   (0 でない: {any(b)})")
print(f"  a * b = {prod}")
print(f"  全部 0 か: **{all(v == 0 for v in prod)}**")

print()
print("=" * 80)
print("② 我々の ±MIN 規約を、この積にかけると")
print("=" * 80)
print("""  total_arith.py の規約:
      saturate(v, zero_ok = (a == 0 or b == 0))
  = 「**両方 0 でないのに結果が 0 なら、それは潰れだ**」 ⟹ ±MIN を返す

  セデニオンでは a も b も 0 でない。だから zero_ok = False。""")
for k in (0, 1, 4, 5, 10, 11, 14, 15):
    v = saturate(float(prod[k]), zero_ok=False)
    print(f"  成分 {k:>2}: 真値 **{prod[k]}** → 我々の規約は **{'+MIN' if v == MIN else v}**  "
          f"{'← **厳密な 0 を、潰れと誤診した**' if prod[k] == 0 and v == MIN else ''}")

print()
print("=" * 80)
print("③ なぜ壊れるか")
print("=" * 80)
print("""  我々の規約の隠れた仮定:

      **a·b = 0  ⟹  a = 0 または b = 0**       （零因子が無い）

  実数では真。**セデニオンでは偽。**（八元数までは真。16次元で初めて壊れる）

  ⟹ **±MIN 規約は「零因子の無い代数」でしか使えない。**
     セデニオンでは、**0 が本物かどうかを、掛ける前の値からは判定できない。**""")

print()
print("=" * 80)
print("④ ではセデニオン回路は、どうやって 0 を見分けているか")
print("=" * 80)
print("""  **見分ける必要が無い。厳密だから。**

  gate9 / gate27 / ripple は **整数の厳密演算**（浮動小数を一切使わない）。
  溢れも潰れも起きないので、**0 は必ず本物の 0**。

  ⟹ **状態符号（3ビット）は、この回路では出番が無い。**
     出番があるのは **rescale**（低い桁を捨てる = 唯一の非可逆な操作）だけで、
     それも実測 **1024 回中 0 回**しか偽ゼロを作らなかった。

  ⟹ **厳密な演算に、不定は生まれない。不定は、捨てたときにだけ生まれる。**""")
