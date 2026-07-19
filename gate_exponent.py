#!/usr/bin/env python3
# ⚠️ 生成AI使用・要検証
"""指数パスの 部品 — 監査で 指摘された「ホスト計算」を 全て ゲート化する 部品ライブラリ。

  指数 = 固定幅 EW ビットの 2の補数バス（ただのビット列・仮数の (p,n) 桁とは 別種の 線）。
  部品（全て AND/OR/NOT/XOR から・st 計上）:
    bus_const / bus_val          … 定数の ROM 化・読み出し（境界）
    bus_add / bus_sub            … リップル加算器・減算器（a + ~b + 1）
    bus_lt / bus_max / bus_min   … 比較（差の符号）＋ mux
    clamp0                       … max(0, x)（符号ビットで マスク）
    mux_bit / mux_bus / mux_digits … 2:1 選択（マスクゲート）
    priority_encoder             … 先頭非零位置 → one-hot ＋ 二進（leading_pos の ゲート版）
    barrel_shift_right_digits    … 可変右シフト（log 段 mux・落とし桁の 非零も 段ごとに 収集）
  契約 K1/K2: シフト量・選択は 全て 信号バスで 受ける。幅は 固定。
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from gate_bilinear import AND, OR, NOT, XOR, full_adder, ZERO, nz, new_counter


# ---------------------------------------------------------------- バスの 境界I/O
def bus_const(v, EW):
    """定数 v → EW ビット 2の補数（合成時定数の ROM 化・ゲート0）。"""
    return [(v >> i) & 1 for i in range(EW)]

def bus_val(B):
    """バス → 符号つき整数（読み出しのみ）。"""
    v = sum(b << i for i, b in enumerate(B[:-1]))
    return v - (B[-1] << (len(B) - 1))


# ---------------------------------------------------------------- 加減算・比較
def bus_add(A, B, st, cin=0):
    """リップル加算（同幅・桁上げ捨て = 2の補数の 環）。"""
    out = []; c = cin
    for a, b in zip(A, B):
        s, c = full_adder(a, b, c, st); out.append(s)
    return out

def bus_sub(A, B, st):
    """A − B = A + ~B + 1。"""
    return bus_add(A, [NOT(b, st) for b in B], st, cin=1)

def bus_lt(A, B, st):
    """A < B（符号つき）= (A−B) の 符号ビット。幅に 1 ビットの 余裕を 前提（契約 K2）。"""
    return bus_sub(A, B, st)[-1]

def mux_bit(s, a, b, st):
    """s ? a : b（4 ゲート）。"""
    return OR(AND(s, a, st), AND(NOT(s, st), b, st), st)

def mux_bus(s, A, B, st):
    return [mux_bit(s, a, b, st) for a, b in zip(A, B)]

def mux_digits(s, DA, DB, st):
    """符号つき桁列の 2:1 選択（両レールを mux）。"""
    return [(mux_bit(s, pa, pb, st), mux_bit(s, na, nb, st))
            for (pa, na), (pb, nb) in zip(DA, DB)]

def bus_max(A, B, st):
    lt = bus_lt(A, B, st)
    return mux_bus(lt, B, A, st)

def bus_min(A, B, st):
    lt = bus_lt(A, B, st)
    return mux_bus(lt, A, B, st)

def clamp0(A, st):
    """max(0, A) = 符号ビットが 1 なら 全ビット 0（マスク）。"""
    ns = NOT(A[-1], st)
    return [AND(b, ns, st) for b in A]


# ---------------------------------------------------------------- 優先エンコーダ
def priority_encoder(digits, EW, st):
    """符号つき桁列の 最上位 非零位置 L を 二進バスで。戻り (L_bus, none_flag, onehot)。
       leading_pos の ゲート版: nz を 高位から OR 累積 → 「自分 非零 かつ 上位 全零」= one-hot
       → one-hot を 二進へ（各ビット = 該当位置の OR）。"""
    n = len(digits)
    nzs = [nz(d, st) for d in digits]
    above = [0] * n                          # above[i] = OR of nzs[j], j>i
    acc = 0
    for i in range(n - 1, -1, -1):
        above[i] = acc
        acc = OR(acc, nzs[i], st)
    none = NOT(acc, st)                      # 全零
    onehot = [AND(nzs[i], NOT(above[i], st), st) for i in range(n)]
    L = []
    for k in range(EW):
        bit = 0
        for i in range(n):
            if (i >> k) & 1:
                bit = OR(bit, onehot[i], st)
        L.append(bit)
    return L, none, onehot


# ---------------------------------------------------------------- バレルシフタ
def barrel_shift_right_digits(digits, S, st):
    """符号つき桁列を 信号バス S だけ 右へ（低位を 落とす）。log 段の mux 層。
       戻り (シフト結果[固定幅], dropped_nz) — 落とした桁の 非零を 段ごとに 収集（ge/le 判定用）。
       契約: 0 ≤ S < 2^len(S)。幅 固定・ZERO 詰め。"""
    n = len(digits)
    cur = list(digits)
    dropped = 0
    for j, sbit in enumerate(S):
        k = 1 << j
        # この段で 落ちる 低位 k 桁の 非零（s=1 のときだけ）
        dnz = 0
        for i in range(min(k, n)):
            dnz = OR(dnz, nz(cur[i], st), st)
        dropped = OR(dropped, AND(sbit, dnz, st), st)
        cur = [mux_digits_one(sbit, cur[i + k] if i + k < n else ZERO, cur[i], st)
               for i in range(n)]
    return cur, dropped

def mux_digits_one(s, da, db, st):
    (pa, na), (pb, nb) = da, db
    return (mux_bit(s, pa, pb, st), mux_bit(s, na, nb, st))

def barrel_shift_left_digits(digits, S, out_width, st):
    """符号つき桁列を 信号バス S だけ 左へ（×2^S・低位 ZERO 詰め）。固定 out_width（契約 K2:
       設計時に 最大シフトを 見込んだ 幅。上位に はみ出た 桁は 設計違反＝呼び出し側の 幅選定責務）。"""
    cur = [digits[i] if i < len(digits) else ZERO for i in range(out_width)]
    for j, sbit in enumerate(S):
        k = 1 << j
        cur = [mux_digits_one(sbit, cur[i - k] if i - k >= 0 else ZERO, cur[i], st)
               for i in range(out_width)]
    return cur


# ---------------------------------------------------------------- 自己テスト
def self_test():
    import numpy as np
    from gate_bilinear import to_sd, from_sd, enc
    rng = np.random.default_rng(20260803)
    EW = 10

    print("=" * 74)
    print("指数バス部品 — 加減算・比較・max/min/clamp（全 2の補数・ゲート）")
    print("=" * 74)
    bad = 0
    for _ in range(4000):
        a = int(rng.integers(-200, 200)); b = int(rng.integers(-200, 200))
        A = bus_const(a, EW); B = bus_const(b, EW); st = new_counter()
        if bus_val(bus_add(A, B, st)) != a + b: bad += 1
        if bus_val(bus_sub(A, B, st)) != a - b: bad += 1
        if bus_lt(A, B, st) != (1 if a < b else 0): bad += 1
        if bus_val(bus_max(A, B, st)) != max(a, b): bad += 1
        if bus_val(bus_min(A, B, st)) != min(a, b): bad += 1
        if bus_val(clamp0(A, st)) != max(0, a): bad += 1
    print(f"  add/sub/lt/max/min/clamp0: 違反 {bad}/24000 ✓")

    print()
    print("=" * 74)
    print("優先エンコーダ — leading_pos の ゲート版（one-hot → 二進）")
    print("=" * 74)
    bad = 0
    for _ in range(3000):
        v = int(rng.integers(-(1 << 14), 1 << 14))
        D = to_sd(v, 16)
        st = new_counter()
        L, none, onehot = priority_encoder(D, EW, st)
        if v == 0:
            if none != 1: bad += 1
        else:
            true_L = abs(v).bit_length() - 1
            if none != 0 or bus_val(L + [0]) != true_L: bad += 1
    # 冗長表現でも 正しいか？ → 正しくない可能性（冗長の罠・canonicalize 前提）を 明示テスト
    red = [enc(-1)] * 10 + [enc(1)]                     # 値 1 の 冗長形（先頭は 位置10）
    L, none, _ = priority_encoder(red, EW, new_counter())
    trap = bus_val(L + [0])
    print(f"  非冗長入力 3000 件: 違反 {bad} ✓")
    print(f"  冗長入力の 罠: 値1 の 冗長形 → 位置 {trap}（≠0）⟹ **canonicalize 後にのみ 使う**（規律 R1）")

    print()
    print("=" * 74)
    print("バレルシフタ — 可変右シフト（mux 層）＋ 落とし桁 非零の 収集")
    print("=" * 74)
    bad = 0
    for _ in range(3000):
        v = int(rng.integers(-(1 << 14), 1 << 14))
        sh = int(rng.integers(0, 16))
        D = to_sd(v, 18)
        st = new_counter()
        out, dropped = barrel_shift_right_digits(D, bus_const(sh, 5), st)
        want = abs(v) >> sh
        want = want if v >= 0 else -want
        got = from_sd(out)
        want_drop = 1 if (abs(v) & ((1 << sh) - 1)) != 0 else 0
        if got != want or dropped != want_drop: bad += 1
    st = new_counter()
    barrel_shift_right_digits(to_sd(12345, 18), bus_const(7, 5), st)
    print(f"  シフト値・落とし桁非零: 違反 {bad}/3000 ✓（1 回の 18桁シフトで {sum(st.values())} ゲート）")

    print()
    print("指数部品 揃った — 全て AND/OR/NOT/XOR・固定幅・信号バス駆動（契約 K1/K2 準拠）。")


if __name__ == "__main__":
    self_test()
