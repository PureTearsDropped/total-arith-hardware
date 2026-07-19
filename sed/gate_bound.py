#!/usr/bin/env python3
"""**境界トリットのゲート**。AND/OR/NOT だけで作れるか。ゲート数は。

  境界トリット b ∈ {−1, 0, +1} = **|x| ≤ v** / **|x| = v** / **|x| ≥ v**
  「境界なし」 = （値 0, b = +1）＝ **|x| ≥ 0** ＝ 常に真 ＝ 何も言っていない
"""
from sedenion_tensor_logic import enc, dec, bnot, gate9

def gate_bmul(a, b):
    """**境界の積**。+1 が支配、次に −1、なければ 0。
       (p,n) 符号化で **OR 2 個 + AND 1 個 + NOT 1 個 = 4 ゲート**（gate9 より安い）。"""
    ap, an = enc(a); bp, bn = enc(b)
    op = ap | bp                                   # どちらかが ≥ なら ≥
    on = (an | bn) & bnot(op)                      # そうでなく、どちらかが ≤ なら ≤
    return dec(op, on)

def vacuous(a, b):
    """**値を 0 に潰すか**（≤ と ≥ がぶつかった ⟹ 境界なし）。AND 2 + OR 1 = 3 ゲート。"""
    ap, an = enc(a); bp, bn = enc(b)
    return (an & bp) | (ap & bn)

NB = {-1: "**≤v**", 0: "**=v**", 1: "**≥v**"}
def want(a, b):
    if a == 0 and b == 0: return 0, False
    if a == -b and a != 0: return 1, True          # **境界なし = ≥0**
    if a == 0: return b, False
    if b == 0: return a, False
    return a, False

print("=" * 84)
print("**境界の積のゲート** — AND/OR/NOT だけで作れるか")
print("=" * 84)
print(f"  {'a':>8}{'b':>10}{'欲しい (b, 値を0に)':>24}{'**gate_bmul**':>16}{'一致':>8}")
ok = n = 0
for a in (-1, 0, 1):
    for b in (-1, 0, 1):
        wb, wv = want(a, b)
        gb, gv = gate_bmul(a, b), bool(vacuous(a, b))
        c = (gb == wb and gv == wv); ok += c; n += 1
        print(f"  {NB[a]:>8}{NB[b]:>10}{f'{NB[wb]}, {wv}':>24}{f'{NB[gb]}, {gv}':>16}"
              f"{('✓' if c else '**✗**'):>8}")
print(f"\n  一致: **{ok}/{n}**")

print()
print("=" * 84)
print("**ゲート数**")
print("=" * 84)
print("""  符号トリット : **gate9**（既存。4 AND + 2 OR = 6 ゲート）  ← **新規ゼロ**
  境界トリット : **gate_bmul**（2 OR + 1 AND + 1 NOT = **4 ゲート**）
                 **vacuous**  （2 AND + 1 OR = **3 ゲート**）
                 値を 0 に潰すマスク: 16 本の AND（1 トリットあたり 1）

  ⟹ **状態の伝播 = 1 セデニオン成分あたり 7 ゲート**（符号 6 + 境界 7 のうち新規 7）。
     比較: セデニオン積 1 個は **gate9 が 9,216 個**（K=6 のとき）。
     ⟹ **状態符号の費用は、演算の 0.1% 未満。**""")

print()
print("=" * 84)
print("**符号と境界は独立に流れるか** — 干渉しないことの確認")
print("=" * 84)
print(f"  {'符号(a,b)':>14}{'境界(a,b)':>16}{'符号の出力':>12}{'境界の出力':>14}")
NS = {1: "正", 0: "**不明**", -1: "負"}
for sa, sb, ba, bb in [(1, 1, -1, -1), (1, -1, -1, 1), (0, 1, 0, 0), (-1, 0, 1, 1)]:
    print(f"  {f'{NS[sa]},{NS[sb]}':>14}{f'{NB[ba]},{NB[bb]}':>16}"
          f"{NS[gate9(sa, sb)]:>12}{NB[gate_bmul(ba, bb)]:>14}")
print("""
  ⟹ **符号は gate9、境界は gate_bmul。互いの入力を見ない ⟹ 並列に流れる。**
     （§7.6.4.1 の「積は符号を保つ」が、**配線として** そうなっている）""")
