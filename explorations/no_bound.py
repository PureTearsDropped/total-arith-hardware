#!/usr/bin/env python3
"""**「境界がある」と「境界が無い」は別**（利用者）。

    **MIN**      「|x| ≤ MIN」        ← 知っている（上界）
    **MAX**      「|x| ≥ MAX」        ← 知っている（下界）
    **大きさ不明** 「MAX 以上かも、MAX と MIN の間かも、MIN 以下かも」← **何も知らない**

§7.6.7④ の表は MAX/MIN を「大きさ不明」に入れていた。**誤り。** どの演算が境界を失うか数える。
"""
import numpy as np
from total_arith import MIN, MAX

# 大きさの状態: "exact" / "le_min"(|x|≤MIN) / "ge_max"(|x|≥MAX) / **"none"(境界なし)**
def mag_mul(x, y):
    if x == "exact" and y == "exact": return "exact"
    if "none" in (x, y): return "none"
    if x == "le_min" and y == "le_min": return "le_min"     # ≤ MIN² ≤ MIN
    if x == "ge_max" and y == "ge_max": return "ge_max"     # ≥ MAX² ≥ MAX
    if {x, y} == {"le_min", "ge_max"}: return "none"        # **≥MAX × ≤MIN ⟹ 何でもありうる**
    if "le_min" in (x, y): return "le_min"                  # 有限 × ≤MIN ⟹ 小さいまま（有限が有界なら）
    if "ge_max" in (x, y): return "ge_max"
    return "none"

def mag_add(x, y, same_sign):
    if x == "exact" and y == "exact": return "exact"
    if "none" in (x, y): return "none"
    if x == "le_min" and y == "le_min": return "le_min"     # |a+b| ≤ 2·MIN
    if x == "ge_max" and y == "ge_max":
        return "ge_max" if same_sign else "none"            # **同符号なら ≥2MAX、異符号なら何でも**
    if {x, y} == {"le_min", "ge_max"}: return "ge_max"      # 巨大 ± 微小 ⟹ 巨大のまま
    if "ge_max" in (x, y): return "ge_max"
    if "le_min" in (x, y): return "exact"                   # 確定 ± 微小 ⟹ ほぼ確定（丸めは別）
    return "none"

print("=" * 90)
print("**どの演算が、境界を失うか**")
print("=" * 90)
names = {"exact": "確定", "le_min": "**|x|≤MIN**", "ge_max": "**|x|≥MAX**", "none": "**境界なし**"}
print(f"  {'式':<28}{'結果の大きさ':>18}   真値は")
rows = [
    ("MIN × MIN",              mag_mul("le_min", "le_min"),  "≤ MIN² ≤ MIN"),
    ("MAX × MAX",              mag_mul("ge_max", "ge_max"),  "≥ MAX² ≥ MAX"),
    ("**MAX × MIN**",          mag_mul("ge_max", "le_min"),  "**≥MAX と ≤MIN の積 ⟹ 何でもありうる**"),
    ("MIN + MIN",              mag_add("le_min", "le_min", True),  "≤ 2·MIN"),
    ("MAX + MAX（同符号）",      mag_add("ge_max", "ge_max", True),  "≥ 2·MAX"),
    ("**MAX − MAX（異符号）**",  mag_add("ge_max", "ge_max", False), "**[0, ∞) のどこか**"),
    ("**MIN − MIN（異符号）**",  mag_add("le_min", "le_min", False), "≤ 2·MIN  ← **境界は在る**"),
    ("MAX + 5.0",              mag_add("ge_max", "exact", True),   "≥ MAX − 5 ≈ MAX"),
    ("5.0 × 3.0",              mag_mul("exact", "exact"),          "= 15"),
]
for tag, st, truth in rows:
    print(f"  {tag:<28}{names[st]:>18}   {truth}")

n_none = sum(1 for _, st, _ in rows if st == "none")
print(f"""
  ⟹ **境界を失うのは {n_none}/{len(rows)} だけ:**
       **MAX × MIN**（下界 × 上界 ⟹ 打ち消し合う）
       **MAX − MAX**（同じ下界どうしの引き算）

  ⟹ **MIN − MIN は境界を失わない**（≤ 2·MIN）。§7.6.7④ が `01-` に置いたのは正しかった。
     **だが「大きさ確定」と書いたのは誤り。正しくは「大きさに境界が在る」。**

  ⟹ **§7.6.7④ が MAX/MIN を「大きさ不明」に入れていたのは誤り。** 両方とも境界を知っている。""")

print()
print("=" * 90)
print("**そして 4 番目の状態が要る** — 大きさの欄")
print("=" * 90)
print("""  §7.6.7④ の bit2 は「飽和した」の意味だった。**それは大きさの状態ではない。**

  正しい大きさの欄（4 値）:
      **確定**       値がそのまま答え
      **|x| ≤ MIN**  上界。**値が MIN であること自体が、それを言っている**
      **|x| ≥ MAX**  下界。**値が MAX であること自体が、それを言っている**
      **境界なし**    MAX×MIN と MAX−MAX からだけ生まれる

  ⟹ **上の 2 つは、値を見れば分かる**（MIN か MAX か）。**ビットが要るのは「境界なし」だけ。**

  ⟹ 符号化は変わらず 3 ビットのまま。**bit2 の意味だけが変わる:**
       ~~bit2 = 飽和した（MAX か MIN）~~
       **bit2 = 大きさに境界が無い**（MAX×MIN / MAX−MAX からのみ）

  ⟹ **MAX も MIN も bit2 = 0 になる。** 境界を知っているので。""")
