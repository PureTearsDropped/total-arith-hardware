#!/usr/bin/env python3
# ⚠️ 生成AI使用・要検証
"""複数値 同時加算器 — 3 つの 形を ゲートで 実装・実測。

  ① 木型   : m 入力を 層別 3:2 圧縮で 同時に（sd_sum_fast・O(log m)＝fan-in 2 の 下限）
  ② 蓄積型 : 冗長 2 行の アキュムレータに 1 個ずつ 吐き込む — **1 加算 = 圧縮 1 層 = 定数深さ**
             （キャリー解決なし・クロックごとに 1 値 = ハードの MAC の 正体）
  ③ BFP 多数同時加算: 共通指数へ 整列 → 冗長木で 厳密に 足す → **丸め（正規化）は 最後に 1 回**
             （IEEE の「足すたび丸める」誤差蓄積が ない = fused dot product）
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from fractions import Fraction as Fr
from gate_bilinear import (new_counter, ZERO, to_sd, from_sd, compress3,
                           canonicalize, nz)
from gate_exponent import bus_const, bus_val, mux_bus, barrel_shift_left_digits
from gate_fast import (B, wrap, wrapb, depth_of, sd_sum_fast, sd_add2,
                       canonicalize_fast, bus_sub_fast, bus_lt_fast,
                       block_normalize_g_fast, or_tree)


# ---------------------------------------------------------------- ② 蓄積型（carry-save）
def accu_new(width):
    """アキュムレータ = 冗長 2 行（キャリー未解決のまま 保持）。"""
    return [[ZERO] * width, [ZERO] * width]

def accu_add_clean(acc, x, st):
    """1 値を 吐き込む: 2 行 + 新値 = 3 行 → 桁ごと 3:2 圧縮 1 層 → 2 行。
       low は 同位置・high は 1 桁 上へ。**キャリー解決なし・深さ 定数/加算**（クロック 1 回 相当）。"""
    A, Bo = acc
    n = max(len(A), len(Bo), len(x)) + 1
    g = lambda R, i: R[i] if i < len(R) else ZERO
    lo = []; hi = [ZERO]
    for i in range(n):
        l, h = compress3(g(A, i), g(Bo, i), g(x, i), st)
        lo.append(l); hi.append(h)
    return [lo + [ZERO], hi]

def accu_value(acc):
    """読み出し（検証用）: 2 行の 値の和。"""
    return from_sd(acc[0]) + from_sd(acc[1])

def accu_resolve(acc, st):
    """最後に 1 回だけ 解決: sd_add2（定数深さ）→ 1 行。"""
    return sd_add2(acc[0], acc[1], st)


# ---------------------------------------------------------------- ③ BFP 多数同時加算
def bus_min_fast(A, Bb, st):
    lt = bus_lt_fast(A, Bb, st)
    return mux_bus(lt, A, Bb, st)

def bfg_sum_many(xs, Dmax_bits, st):
    """N 個の BFP を 同時に: 共通 Emin へ 全整列（バレル）→ sd_sum_fast（木）→ E=Emin。
       戻りは 未正規化（丸めゼロ・厳密）。正規化 = 最後に 1 回だけ（呼び出し側）。"""
    Es = [x.Ebus for x in xs]
    while len(Es) > 1:                                    # Emin トーナメント木
        Es = [bus_min_fast(Es[i], Es[i + 1], st) if i + 1 < len(Es) else Es[i]
              for i in range(0, len(Es), 2)]
    Emin = Es[0]
    ow = max(len(x.mant) for x in xs) + (1 << Dmax_bits)
    rows = []
    for x in xs:
        d = bus_sub_fast(x.Ebus, Emin, st)[:Dmax_bits]    # ≥0（Emin が 最小）
        rows.append(barrel_shift_left_digits(x.mant, d, ow, st))
    from gate_bfp2 import BFg
    return BFg(sd_sum_fast(rows, st), Emin)


# ---------------------------------------------------------------- 自己テスト
def self_test():
    import numpy as np
    rng = np.random.default_rng(20260806)

    print("=" * 78)
    print("② 蓄積型 — 1 加算 = 圧縮 1 層（定数）。キャリー解決は 最後に 1 回")
    print("=" * 78)
    bad = 0
    for _ in range(800):
        m = int(rng.integers(2, 20))
        xs = [int(v) for v in rng.integers(-2000, 2000, m)]
        acc = accu_new(14)
        st = new_counter()
        for v in xs:
            acc = accu_add_clean(acc, to_sd(v, 14), st)
        if accu_value(acc) != sum(xs): bad += 1               # 冗長のままでも 値は 厳密
        if from_sd(accu_resolve(acc, st)) != sum(xs): bad += 1
    print(f"  値保存（冗長のまま・解決後とも）: 違反 {bad}/1600 ✓")
    # 深さ実測: 1 加算あたりの 増分が 一定か（= クロック 1 回分の 組合せ深さ）
    acc = accu_new(14); ds = []
    st = new_counter()
    for v in (123, -456, 789, -321, 654, -987, 111, -222):
        acc = accu_add_clean(acc, wrap(to_sd(v, 14)), st)
        ds.append(depth_of(acc))
    inc = [ds[i + 1] - ds[i] for i in range(len(ds) - 1)]
    print(f"  深さの 伸び（8 回 吐き込み）: {ds}  増分 {inc}")
    print(f"  ⟹ **1 加算あたり 定数 {inc[0]} 段**（レジスタで 切れば クロックごとに 1 値・MAC の 正体）")

    print()
    print("=" * 78)
    print("①vs② 木型と 直列鎖 — 同じ 16 値の 和")
    print("=" * 78)
    xs = [int(v) for v in rng.integers(-2000, 2000, 16)]
    tree = depth_of(sd_sum_fast([wrap(to_sd(v, 14)) for v in xs], new_counter()))
    acc = accu_new(14); st = new_counter()
    for v in xs: acc = accu_add_clean(acc, wrap(to_sd(v, 14)), st)
    serial = depth_of(acc)
    print(f"  木型（全部 同時・組合せ）: 深さ {tree}   蓄積型（直列 組合せなら）: {serial}")
    print(f"  ⟹ 一度に 全部 → 木 O(log m)（fan-in2 の 下限）／ 流れてくる → 蓄積（定数/クロック）")

    print()
    print("=" * 78)
    print("③ BFP 多数同時加算 — 指数バラバラの N 値を 整列→厳密木→丸めは 最後に 1 回")
    print("=" * 78)
    from gate_bfp2 import to_bfg, from_bfg, bfg_add, EW
    bad = 0
    for _ in range(300):
        N = int(rng.integers(3, 9))
        vals = [int(rng.integers(-300, 300)) for _ in range(N)]
        Es = [int(rng.integers(0, 8)) for _ in range(N)]
        xs = [to_bfg(v, 12, E) for v, E in zip(vals, Es)]
        st = new_counter()
        s = bfg_sum_many(xs, 3, st)
        exact = sum(Fr(v) * Fr(2) ** E for v, E in zip(vals, Es))
        if from_bfg(s) != exact: bad += 1                     # 丸めゼロ = 厳密
    print(f"  bfg_sum_many == 厳密和（丸めゼロ・Fraction 照合）: 違反 {bad}/300 ✓")
    # 深さ: 直列 bfg_add 鎖 と 同時加算器
    N = 8
    vals = [int(v) for v in rng.integers(-300, 300, N)]
    Es = [int(rng.integers(0, 8)) for _ in range(N)]
    def wrapbfg(x):
        x.mant = wrap(x.mant); x.Ebus = wrapb(x.Ebus); return x
    xs = [wrapbfg(to_bfg(v, 12, E)) for v, E in zip(vals, Es)]
    st = new_counter()
    d_many = depth_of(bfg_sum_many(xs, 3, st).mant)
    xs2 = [wrapbfg(to_bfg(v, 12, E)) for v, E in zip(vals, Es)]
    accx = xs2[0]; st = new_counter()
    from gate_fast import bfg_add_fast
    for x in xs2[1:]:
        accx = bfg_add_fast(accx, x, 3, st)
    d_chain = depth_of(accx.mant)
    print(f"  深さ実測（N=8・指数バラバラ）: 直列 bfg_add 鎖 {d_chain} → 同時加算器 **{d_many}**")
    print(f"  ⟹ 整列は 全値 並列・和は 木・**正規化（丸め）は 最後に 1 回だけ**")
    print(f"     = IEEE の「足すたび 丸める」誤差蓄積が 構造的に ない（fused dot product の 形）")

    print()
    print("複数値 同時加算器: 木型 O(log m)（下限）・蓄積型 定数/クロック・BFP 版 丸め 1 回。全部 ゲート。")


if __name__ == "__main__":
    self_test()
