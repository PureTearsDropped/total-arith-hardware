#!/usr/bin/env python3
"""**MIN/MAX は点ではなく、境界である**（利用者）。

    **MIN = 「絶対値 ≤ MIN。値は不明」**   ← 上界
    **MAX = 「絶対値 ≥ MAX。値は不明」**   ← 下界（向きが逆）

この読みなら、飽和した値は **嘘をついていない。緩い境界を言っているだけ**である。
そして「大きさ不明」の印が、それを言っている。

問い: **すべての演算が、健全な境界を返すか。**（＝ 真値が必ず境界の中に居るか）
"""
import numpy as np
from total_arith import MIN, MAX, saturate, total_mul

# 状態: (値, 符号確定か, 大きさ確定か)。MIN/MAX は 大きさ不確定。
def wrap(v):
    v = float(v)
    mag_ok = not (abs(v) == MIN or abs(v) == MAX)
    return (v, True, mag_ok)

def tmul(x, y):
    (a, sa, ma), (b, sb, mb) = x, y
    v = float(saturate(np.float64(a)*np.float64(b), zero_ok=(a == 0.0 or b == 0.0)))
    mag_ok = ma and mb and abs(v) != MIN and abs(v) != MAX
    return (v, sa and sb, mag_ok)

def tadd(x, y):
    (a, sa, ma), (b, sb, mb) = x, y
    v = float(saturate(np.float64(a)+np.float64(b),
                       zero_ok=(float(np.float64(a)+np.float64(b)) == 0.0 and abs(a) == abs(b))))
    BIG  = lambda t: abs(t) >= MAX*0.999
    TINY = lambda t: t != 0.0 and abs(t) <= MIN*1.001
    # **符号不明: 同じ端で飽和した異符号の加減**（§7.6.4.1）
    sign_ok = sa and sb
    if np.sign(a) != np.sign(b) and ((BIG(a) and BIG(b)) or (TINY(a) and TINY(b))):
        sign_ok = False
    mag_ok = ma and mb and abs(v) != MIN and abs(v) != MAX
    return (v, sign_ok, mag_ok)

print("=" * 88)
print("① **境界は健全か** — 真値が、返された境界の中に居るか")
print("=" * 88)
print("""  読み:  値が **MIN** なら「真の |x| ≤ MIN」
         値が **MAX** なら「真の |x| ≥ MAX」
         それ以外（大きさ確定）なら「真の x = その値」""")
def contains(state, true_val):
    v, so, mo = state
    if mo: return abs(true_val - v) <= abs(v)*1e-9 + 1e-320   # 大きさ確定 ⟹ 一致すべき
    if abs(v) == MIN: return abs(true_val) <= MIN              # **上界**
    if abs(v) == MAX: return abs(true_val) >= MAX              # **下界**
    return True
print(f"\n  {'式':<30}{'真値':>14}{'返る値':>14}{'大きさ':>10}{'健全':>8}")
tests = [
    ("1e−200 × 1e−200", 1e-200, 1e-200, "*", 1e-400),
    ("1e−300 × 1e−300", 1e-300, 1e-300, "*", 1e-600),
    ("MIN × MIN",       MIN,    MIN,    "*", MIN*MIN if MIN*MIN > 0 else 0.0),
    ("1e300 × 1e300",   1e300,  1e300,  "*", float('inf')),
    ("2.0 × 3.0",       2.0,    3.0,    "*", 6.0),
]
ok = 0
for tag, a, b, _, truth in tests:
    st = tmul(wrap(a), wrap(b))
    c = contains(st, truth)
    ok += c
    fv = "**MIN**" if abs(st[0]) == MIN else ("**MAX**" if abs(st[0]) == MAX else f"{st[0]:.4g}")
    ft = f"{truth:.4g}" if np.isfinite(truth) else "**∞**"
    print(f"  {tag:<30}{ft:>14}{fv:>14}{('確定' if st[2] else '**不明**'):>10}{('✓' if c else '**✗**'):>8}")
print(f"\n  健全: **{ok}/{len(tests)}**")

print()
print("=" * 88)
print("② **加減の境界** — MIN+MIN は MIN か、2·MIN か")
print("=" * 88)
print(f"""  読みに従うと:
    |a| ≤ MIN, |b| ≤ MIN  ⟹  **|a+b| ≤ 2·MIN**    ← MIN より **大きい**。表現できる
    |a| ≥ MAX, |b| ≥ MAX  ⟹  **|a+b| ≥ 2·MAX**    ← 溢れる ⟹ MAX（下界としては正しい）

  実装が返す値:
    MIN + MIN = **{'MIN' if float(saturate(np.float64(MIN)+np.float64(MIN), zero_ok=False)) == MIN else f'{float(np.float64(MIN)+np.float64(MIN)):.4g} = 2·MIN'}**
    MAX + MAX = **{'MAX' if float(saturate(np.float64(MAX)+np.float64(MAX), zero_ok=False)) == MAX else '?'}**

  ⟹ 2·MIN は表現できるので、**MIN+MIN = 2·MIN が正しい上界**（より締まる）。
     MAX+MAX は表現できないので、**MAX（緩いが健全な下界）**。""")

print()
print("=" * 88)
print("③ **印を読めば、診断器は直るか** — L2 が 大きさ不明 なら L∞ へ落ちる")
print("=" * 88)
def norm_marked(v):
    """L2 を試す。**大きさ不明の印が立ったら L∞ へ落ちる。**"""
    sq = [tmul(wrap(x), wrap(x)) for x in v]
    acc = (0.0, True, True)
    for s in sq: acc = tadd(acc, s)
    if not acc[2]:                                  # **大きさ不明 ⟹ この道は使えない**
        return float(np.max(np.abs(v))), "L∞ へ落ちた"
    return float(np.sqrt(acc[0])), "L2 のまま"
print(f"  {'v':<24}{'返る値':>16}{'どうしたか':>16}{'正しいか':>12}")
for tag, v, want in [("[1e−200, 0, ...]", [1e-200]+[0.0]*15, 1e-200),
                     ("[1e−300, 0, ...]", [1e-300]+[0.0]*15, 1e-300),
                     ("[3.0, 4.0, 0...]", [3.0, 4.0]+[0.0]*14, 5.0),
                     ("[1.0, 0, ...]",    [1.0]+[0.0]*15,     1.0)]:
    got, how = norm_marked(v)
    # L∞ に落ちた場合は max|x| が正解（正規化に使える）／L2 のままなら真のノルム
    good = np.isclose(got, want, rtol=1e-9) or (how.startswith("L∞") and np.isclose(got, max(abs(x) for x in v)))
    print(f"  {tag:<24}{got:>16.4g}{how:>16}{('**✓**' if good else '**✗**'):>12}")
