#!/usr/bin/env python3
"""**どの演算が、どのフラグを立てるか**（利用者）。

  利用者: 「**積や和でオーバーフローしたら〜以上フラグを、
             積や商で min 以下になったら〜以下フラグをつける？**」

  ＝ **溢れ → ≥ フラグ ／ 潰れ → ≤ フラグ**。そして **和は潰れの側に入っていない。**
  それは正しいか。総当たりで数える。
"""
import warnings, numpy as np
from fractions import Fraction as Fr
warnings.filterwarnings("ignore")
MAX = np.finfo(np.float64).max
MIN = np.nextafter(0.0, 1.0)
TINY = np.finfo(np.float64).tiny        # 最小の正規化数 2.2e-308
f = np.float64

def exact(op, a, b):
    A, B = Fr(a), Fr(b)
    if op == "/" and B == 0: return None
    return {"+": A+B, "-": A-B, "*": A*B, "/": (A/B if B else None)}[op]

def flags(op, a, b):
    """真値と float の結果を比べ、**溢れたか／潰れたか**を返す。"""
    ex = exact(op, a, b)
    if ex is None: return None
    with np.errstate(all="ignore"):
        got = float({"+": f(a)+f(b), "-": f(a)-f(b), "*": f(a)*f(b), "/": f(a)/f(b)}[op])
    over  = not np.isfinite(got) or (abs(ex) > Fr(MAX))          # 表現できるより大きい
    under = (ex != 0) and (abs(ex) < Fr(MIN))                    # 表現できるより小さい
    return over, under

print("=" * 84)
print("**どの演算が、どのフラグを立てられるか** — 総当たり")
print("=" * 84)
vals = [0.0, 1.0, -1.0, 1e-320, -1e-320, TINY, 1e-200, 1e200, MAX, -MAX, MIN, 1e-308, 1e308]
tab = {op: {"over": 0, "under": 0, "n": 0} for op in "+-*/"}
for op in "+-*/":
    for a in vals:
        for b in vals:
            r = flags(op, a, b)
            if r is None: continue
            o, u = r
            tab[op]["n"] += 1; tab[op]["over"] += o; tab[op]["under"] += u
print(f"  {'演算':<8}{'組':>8}{'**溢れた（→ ≥ フラグ）**':>26}{'**潰れた（→ ≤ フラグ）**':>26}")
for op, name in (("+", "**和**"), ("-", "**差**"), ("*", "**積**"), ("/", "**商**")):
    t = tab[op]
    fo = f"**{t['over']}**" if t['over'] else "**0 ← 起きない**"
    fu = f"**{t['under']}**" if t['under'] else "**0 ← 起きない**"
    print(f"  {name:<8}{t['n']:>8}{fo:>26}{fu:>26}")

print("""
  ⟹ **利用者の規則が、そのまま出た:**
       **溢れ → ≥ フラグ**  : 和・差・積・商 すべて
       **潰れ → ≤ フラグ**  : **積と商だけ。和と差では 一度も起きない**

  ⟹ 理由（Sterbenz の補題）: `a/2 ≤ b ≤ 2a` なら `a − b` は **厳密**（誤差ゼロ）。
     近い数の引き算は表現できる。**加減は、大きさを表現の下へ押し出せない。**
     **押し出せるのは、掛け算と割り算だけ。**""")

print()
print("=" * 84)
print("**そして そのフラグは、既にハードウェアにある**")
print("=" * 84)
print("""  x86 の MXCSR レジスタ:
      **OE（Overflow）**   ビット 3   ← **我々の ≥ フラグ**
      **UE（Underflow）**  ビット 4   ← **我々の ≤ フラグ**
  IEEE 754 が 1985 年から要求している 5 つの例外のうちの 2 つ。

  ⟹ **FPU は既に、毎回これを計算している。**
     違いは **置き場所** だけ:
       IEEE   : **粘着ビットとしてレジスタへ**（どの値が溢れたかは分からない。誰も読まない）
       この規約: **値そのものに付ける**（どの値が、どちらへ溢れたかが分かる）""")
