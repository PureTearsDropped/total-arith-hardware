#!/usr/bin/env python3
# ⚠️ 生成AI使用・要検証
"""配線パターンのカタログ — (U,V,W) 双線形アルゴリズム を 育てる。ゲートで 厳密検証。

  各パターン = 一つの (U,V,W) ＝ 構造テンソルの rank-R 分解 ＝ 「R 個の積で 計算する 配線」。
  ランク削減（R < 素朴）＝ 非群の 速いアルゴリズム。全部 同じ ゲートユニットで 動く。

  収録:
    Strassen 2×2 行列積   R=7 (素朴8)   ← gate_bilinear に既出
    Karatsuba 多項式積    R=3 (素朴4)   ← ここ
    Gauss 複素積          R=3 (素朴4)   ← ここ（群の複素は R=M²=4、これは 3 に削減）
  （多層＝算術回路の例は newton_recip.py: Newton 逆数）
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from gate_bilinear import new_counter, bilinear_unit_gates, from_sd

# ---- パターン定義: (名前, U, V, W, 入力m, 入力n, 参照関数, R, 素朴R) ----
def karatsuba_ref(a, b):
    """2 桁 (a0+a1 t)(b0+b1 t) の 係数 (c0,c1,c2)。"""
    a0, a1 = a; b0, b1 = b
    return [a0*b0, a0*b1 + a1*b0, a1*b1]

def gauss_ref(a, b):
    """複素 (a0+a1 i)(b0+b1 i) → (実, 虚)。"""
    a0, a1 = a; b0, b1 = b
    return [a0*b0 - a1*b1, a0*b1 + a1*b0]

PATTERNS = [
    dict(name="Karatsuba 多項式積", R=3, naive=4, m=2, n=2, ref=karatsuba_ref,
         U=[[1,0],[0,1],[1,1]],
         V=[[1,0],[0,1],[1,1]],
         W=[[1,0,0],[-1,-1,1],[0,1,0]]),
    dict(name="Gauss 複素積",       R=3, naive=4, m=2, n=2, ref=gauss_ref,
         U=[[1,1],[1,0],[0,1]],           # (a0+a1, a0, a1)
         V=[[1,0],[-1,1],[1,1]],          # (b0, b1−b0, b0+b1)
         W=[[1,0,-1],[1,1,0]]),           # 実=p1−p3, 虚=p1+p2
]


def self_test():
    import numpy as np
    rng = np.random.default_rng(20260730)
    print("="*76)
    print("配線パターン・カタログ — (U,V,W)=双線形アルゴリズム、ゲートで 厳密検証")
    print("="*76)
    print(f"  {'パターン':<20}{'積の数R':>8}{'素朴':>6}{'削減':>6}   ゲート版==参照")
    for p in PATTERNS:
        st = new_counter(); bad = 0
        for _ in range(2000):
            a = [int(v) for v in rng.integers(-40, 41, p['m'])]
            b = [int(v) for v in rng.integers(-40, 41, p['n'])]
            out = bilinear_unit_gates(p['U'], p['V'], p['W'], a, b, 14, st)
            c = [from_sd(o) for o in out]
            if c != p['ref'](a, b): bad += 1
        mark = "✓" if bad == 0 else f"× ({bad})"
        print(f"  {p['name']:<20}{p['R']:>8}{p['naive']:>6}{p['naive']-p['R']:>5}積   {mark}")

    print()
    print("  ⟹ どれも 同じ ゲートユニット（bilinear_unit_gates）で 厳密。配線(U,V,W)を 変えるだけ。")
    print("     ランク削減 R<素朴 ＝ 非群の 速いアルゴリズム ＝ 別の 配線パターン。")
    print("     多層（算術回路）の 例は newton_recip.py（Newton 逆数・近似商＋ge境界）。")
    print("     次に 足せる候補: 3×3行列(Laderman 23積) / 逆平方根rsqrt / exp / DFTバタフライ / 二面体群…")


if __name__ == "__main__":
    self_test()
