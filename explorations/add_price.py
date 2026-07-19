#!/usr/bin/env python3
"""**加算は必ず何かを払う。符号か、締まりか**（利用者）。

  利用者: 「でも符号不明になる。**同じ符号同士の min の加算は範囲が広がる**」

      **同符号**: 符号は残る      **範囲が広がる**（p, q → p+q）
      **異符号**: **符号を失う**   範囲は広がらない（→ max(p,q)。**p+q より締まる**）

  ⟹ タダの加算は在るか。累積するとどうなるか。
"""
from fractions import Fraction as Fr
import numpy as np
MIN = np.nextafter(0.0, 1.0)
MAX = np.finfo(np.float64).max

def scan(p, q, same, upper=True):
    P, Q = Fr(p), Fr(q)
    A = [P, P/2, P/1000] if upper else [P, P*2, P*1000]
    B = [Q, Q/2, Q/1000] if upper else [Q, Q*2, Q*1000]
    mags, signs = [], set()
    for a in A:
        for b in B:
            r = a + b if same else a - b
            mags.append(abs(r)); signs.add(0 if r == 0 else (1 if r > 0 else -1))
    return (max(mags) if upper else min(mags)), (len(signs) == 1)

print("=" * 86)
print("① **加算の代金** — 何を払うか")
print("=" * 86)
print(f"  {'場合':<22}{'境界':>18}{'p+q と比べて':>16}{'符号':>12}{'払ったもの':>16}")
p = q = Fr(MIN)
for tag, same in (("**同符号**（+MIN + +MIN）", True), ("**異符号**（+MIN − +MIN）", False)):
    m, s = scan(MIN, MIN, same, True)
    r = Fr(m) / (p + q)
    cmp = f"**{r}×**" if r != 1 else "**同じ**"
    paid = "**範囲**（2倍に広がる）" if s else "**符号**"
    print(f"  {tag:<22}{f'≤ {Fr(m)/Fr(MIN)}·MIN':>18}{cmp:>16}"
          f"{('確定' if s else '**不明**'):>12}{paid:>16}")
print("""
  ⟹ **同符号は 符号を残して 範囲を 2 倍に広げる。異符号は 範囲を保って 符号を捨てる。**
     **ちょうど一つずつ。タダの加算は無い。**""")

print()
print("=" * 86)
print("② **タダの加算は在るか**")
print("=" * 86)
print(f"  {'場合':<30}{'結果':>24}{'払ったもの':>18}")
rows = [
    ("**確定 + 確定**（溢れなし）",  "**確定**",                    "**何も払わない**"),
    ("確定 + (≤q)",                "**v−q ≤ |x| ≤ v+q**",        "**締まり**（両側の境界へ）"),
    ("(≤p) + (≤q) 同符号",          "≤ p+q",                     "**締まり**（2倍）"),
    ("(≤p) + (≤q) 異符号",          "≤ max(p,q)",                "**符号**"),
    ("(≥p) + (≥q) 同符号",          "≥ p+q",                     "**何も払わない**（下界は締まる）"),
    ("(≥p) + (≥q) 異符号",          "**境界なし**",                "**符号と締まり 両方**"),
]
for tag, res, paid in rows:
    print(f"  {tag:<30}{res:>24}{paid:>18}")
print("""
  ⟹ **タダなのは 2 つだけ:**
       **確定 + 確定**（溢れなければ。これは普通の演算）
       **(≥p) + (≥q) 同符号** — **下界どうしは、足すと締まる**（p+q > max(p,q)）

  ⟹ **下界（MAX 側）は、同符号の加算で 情報が増える。上界（MIN 側）は 減る。**
     **ここでも 向きが逆。**（§7.7.3.5）""")

print()
print("=" * 86)
print("③ **累積すると** — 同符号の MIN を N 個足す")
print("=" * 86)
print(f"  {'N':>10}{'境界':>18}{'まだ MIN 以下か':>18}")
for N in (1, 2, 10, 1000, 10**6, 10**16):
    b = Fr(N) * Fr(MIN)
    print(f"  {N:>10,}{f'≤ {N}·MIN':>18}{str(b <= Fr(MIN)):>18}")
TINY = np.finfo(np.float64).tiny
print(f"""
  ⟹ **N 個足すと ≤ N·MIN。線形に緩む。**
     N=1e16 で {float(Fr(10**16)*Fr(MIN)):.3g}。正規化数の底は **{TINY:.4g}**
     ⟹ **{float(Fr(10**16)*Fr(MIN)) > TINY and '正規化の範囲に **出る**' or '出ない'}**

  ⟹ **そしてそれでよい。** §7.7.3.5 の 2 ビット方式は **どの値の上にも境界を立てられる**ので、
     `≤ 4.94e−308`（**普通の数の上の上界**）は正しい状態である。
     **「出られない」なら、むしろ設計が壊れていた**（MIN の上にしか境界を置けないことになる）。

  ⟹ **「≤MIN」は 加算で壊れない。値が動いて、緩むだけ。**""")
