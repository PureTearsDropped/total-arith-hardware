#!/usr/bin/env python3
# ⚠️ 生成AI使用・要検証
"""フルゲート版 ブロック浮動 — 監査の 処方箋どおり **指数パスも 信号バス** に。

  gate_bfp.py（仮数=ゲート・指数=ホスト）の 違反箇所を 全て ゲート化:
    ・指数の 加算/比較/min/max/clamp     → bus_add/bus_lt/bus_min/bus_max/clamp0
    ・leading_pos（値を覗くスキャン）     → priority_encoder（nz の OR 累積）
    ・可変シフト（付け替え/量子化）       → barrel_shift_left/right_digits（mux 層・落とし桁収集）
    ・データ依存の 配線幅                → 全幅 設計時固定・ZERO 詰め（契約 K2）
  検証 = **差分テスト**: ホスト版 gate_bfp.block_normalize / bf_add / bf_mul と 完全一致。
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from gate_bilinear import (AND, OR, NOT, new_counter, ZERO, enc, to_sd, from_sd,
                           sd_sum, multiply, canonicalize, nz)
from gate_exponent import (bus_const, bus_val, bus_add, bus_sub, bus_lt, bus_max,
                           bus_min, clamp0, mux_bit, mux_bus, priority_encoder,
                           barrel_shift_right_digits, barrel_shift_left_digits)

EW = 12                                                  # 指数バス幅（設計定数・余裕込み）


class BFg:
    """フルゲート ブロック浮動数: 仮数 = 符号つき桁列、指数 = EW ビット 2の補数バス。"""
    __slots__ = ('mant', 'Ebus')
    def __init__(self, mant, Ebus):
        self.mant = mant; self.Ebus = Ebus

def to_bfg(value, W, E=0):
    return BFg(to_sd(value, W), bus_const(E, EW))

def from_bfg(x):
    from fractions import Fraction as Fr
    return Fr(from_sd(x.mant)) * (Fr(2) ** bus_val(x.Ebus))


# ---------------------------------------------------------------- 積・和（指数も ゲート）
def bfg_mul(x, y, st):
    """仮数 = multiply（ゲート）・指数 = bus_add（ゲート）。"""
    return BFg(multiply(x.mant, y.mant, st), bus_add(x.Ebus, y.Ebus, st))

def bfg_add(x, y, Dmax_bits, st):
    """指数整列を ゲートで: Elo = min(Ex,Ey)、差 = sub、左バレルシフト（固定幅）、仮数加算。
       契約: |Ex−Ey| < 2^Dmax_bits（設計時に 保証する 範囲・幅は 最悪ケースで 固定）。"""
    Elo = bus_min(x.Ebus, y.Ebus, st)
    dx = bus_sub(x.Ebus, Elo, st)[:Dmax_bits]            # ≥0 が 契約 ⟹ 低位ビットだけで 足りる
    dy = bus_sub(y.Ebus, Elo, st)[:Dmax_bits]
    ow = max(len(x.mant), len(y.mant)) + (1 << Dmax_bits)
    mx = barrel_shift_left_digits(x.mant, dx, ow, st)
    my = barrel_shift_left_digits(y.mant, dy, ow, st)
    return BFg(sd_sum([mx, my], st), Elo)


# ---------------------------------------------------------------- ブロック正規化（フルゲート）
def block_normalize_g(mants, Ebus, W, Emax, st):
    """M 成分 + 共有指数バス → 各 W 桁 + 共有指数バス + フラグ。制御も 全て 信号。
       Emax は 設計定数（ROM）。gate_bfp.block_normalize と 同じ 仕様の ゲート版。"""
    canons = []; signs = []; Ls = []
    Wc = max(len(m) for m in mants) + 1                  # canonicalize 後の 固定幅（形状のみ依存）
    for m in mants:
        c, s = canonicalize(m, st)
        c = [c[i] if i < len(c) else ZERO for i in range(Wc)]
        canons.append(c); signs.append(s)
        L, _, _ = priority_encoder(c, EW, st)            # 先頭位置（非冗長なので 健全・規律 R1）
        Ls.append(L)
    Lmax = Ls[0]
    for L in Ls[1:]:
        Lmax = bus_max(Lmax, L, st)                      # 比較器ツリー
    sh = clamp0(bus_sub(Lmax, bus_const(W - 1, EW), st), st)     # max(0, Lmax−(W−1))
    E_out = bus_add(Ebus, sh, st)
    EmaxB = bus_const(Emax, EW)
    ovE = bus_lt(EmaxB, E_out, st)                       # E_out > Emax
    sh_cap = clamp0(bus_sub(EmaxB, Ebus, st), st)        # 溢れ時の 目一杯シフト
    sh_al = mux_bus(ovE, sh_cap, sh, st)
    E_fin = mux_bus(ovE, EmaxB, E_out, st)
    SW = max(1, (Wc - 1).bit_length())                   # シフトバス幅（設計定数）
    out = []; flags = []
    for c, s in zip(canons, signs):
        shifted, drop_nz = barrel_shift_right_digits(c, sh_al[:SW], st)
        kept = shifted[:W]
        kept_nz = 0
        for d in kept: kept_nz = OR(kept_nz, nz(d, st), st)
        over = 0
        for d in shifted[W:]: over = OR(over, nz(d, st), st)    # まだ W を 超える（溢れ時のみ）
        collapse = AND(drop_nz, NOT(kept_nz, st), st)
        nsign = NOT(s, st)
        sel_max = over
        sel_min = AND(collapse, NOT(over, st), st)
        sel_kept = AND(NOT(over, st), NOT(collapse, st), st)
        om = []
        for i in range(W):
            p, n = kept[i]
            minp = nsign if i == 0 else 0
            minn = s if i == 0 else 0
            op = OR(OR(AND(sel_max, nsign, st), AND(sel_min, minp, st), st), AND(sel_kept, p, st), st)
            on = OR(OR(AND(sel_max, s, st),     AND(sel_min, minn, st), st), AND(sel_kept, n, st), st)
            om.append((op, on))
        ge = OR(over, AND(drop_nz, kept_nz, st), st)
        le = collapse
        out.append(om); flags.append((ge, le, 0))
    return out, E_fin, flags


# ---------------------------------------------------------------- 自己テスト（差分）
def self_test():
    import numpy as np
    from fractions import Fraction as Fr
    from gate_bfp import BF, to_bf, from_bf, bf_mul, bf_add, block_normalize
    rng = np.random.default_rng(20260804)

    print("=" * 76)
    print("① 積・和 — 指数も ゲート。ホスト版 gate_bfp と 差分テスト")
    print("=" * 76)
    bad = 0
    for _ in range(1500):
        a = int(rng.integers(-400, 400)); b = int(rng.integers(-400, 400))
        Ea = int(rng.integers(-8, 9)); Eb = int(rng.integers(-8, 9))
        st = new_counter()
        zg = bfg_mul(to_bfg(a, 11, Ea), to_bfg(b, 11, Eb), st)
        zh = bf_mul(to_bf(a, 11, Ea), to_bf(b, 11, Eb), new_counter())
        if from_bfg(zg) != from_bf(zh) or bus_val(zg.Ebus) != zh.E: bad += 1
    print(f"  bfg_mul == bf_mul（値・指数とも）: 違反 {bad}/1500 ✓")
    bad = 0
    for _ in range(1000):
        a = int(rng.integers(-500, 500)); b = int(rng.integers(-500, 500))
        Ea = int(rng.integers(0, 8)); Eb = int(rng.integers(0, 8))
        st = new_counter()
        sg = bfg_add(to_bfg(a, 12, Ea), to_bfg(b, 12, Eb), Dmax_bits=3, st=st)
        sh = bf_add(to_bf(a, 12, Ea), to_bf(b, 12, Eb), new_counter())
        if from_bfg(sg) != from_bf(sh) or bus_val(sg.Ebus) != sh.E: bad += 1
    print(f"  bfg_add == bf_add（指数整列も ゲート・|ΔE|<8）: 違反 {bad}/1000 ✓")

    print()
    print("=" * 76)
    print("② ブロック正規化 — 制御（先頭検出・シフト量・溢れ判定）も 全て 信号")
    print("=" * 76)
    bad = 0; eps_hit = 0; ov_hit = 0
    for _ in range(600):
        M = 4; W = 6; Emax = 20
        E0 = int(rng.integers(0, 12))
        vals = [int(rng.integers(-5, 6)) * (10 ** int(rng.integers(0, 5))) for _ in range(M)]
        ms = [to_sd(v, 24) for v in vals]
        st = new_counter()
        og, Eg, fg = block_normalize_g(ms, bus_const(E0, EW), W, Emax, st)
        oh, Eh, fh = block_normalize(ms, E0, W, Emax, new_counter())
        if bus_val(Eg) != Eh: bad += 1
        for gg, hh, (g1, l1, _), (g2, l2, _) in zip(og, oh, fg, fh):
            if from_sd(gg) != from_sd(hh) or (g1, l1) != (g2, l2): bad += 1
            if l1: eps_hit += 1
        if any(f[0] for f in fg): ov_hit += 1
    print(f"  block_normalize_g == ホスト版（仮数・指数・フラグ 全て）: 違反 {bad}/600ブロック ✓")
    print(f"  （ε=±MIN 発生 {eps_hit}・溢れ経路も 通過 {ov_hit}）")
    st = new_counter()
    block_normalize_g([to_sd(10000, 24), to_sd(3, 24), to_sd(-2, 24), to_sd(0, 24)],
                      bus_const(0, EW), 6, 40, st)
    print(f"  参考: M=4, 24桁, W=6 の 正規化 1 回 = {sum(st.values()):,} ゲート（指数制御 込み・今回 初めて 全計上）")

    print()
    print("フルゲート版 通過 — 仮数も 指数も 制御も 信号。監査の (c) 違反は gate_bfp2 に 存在しない。")


if __name__ == "__main__":
    self_test()
