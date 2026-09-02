#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証（research ブランチ・未監査）
"""gate_funcs — Level 2: exp / expm1 を **端から端までゲート**で（指数バスも含めて ホスト計算ゼロ）。
   funcs_spec.Spec（Level 1）と bit 一致することを差分テストで確かめ、ゲート数と深さを実測する。

  配線（f32 構成: Win=Wout=24, Wk=40, Pc=56, Q=40, XC=7; f64 も同じ図で幅だけ違う）:
    入力  x = (Win 桁 符号つき桁, E バス)   ── canonicalize(R1) → (大きさビット, 符号) → 先頭位置 L（優先エンコーダ木）
    ld = L + E（バス加算）: small = ld ≤ −3, big = ld ≥ XC, zero = 全桁 0
    中経路: xf = N·2^(E+Pc)（左バレル; 幅 XC+Pc+1。ld∈[−2,XC) では E+Pc ≥ Pc−Win−1 > 0 で 切り捨ては起きない）
            big なら xf = ±2^(XC+Pc) に mux（clamp）
            k = floor((xf·INV_LN2 + 2^(Pc+Q−1)) / 2^(Pc+Q)) …… |xf|·INV_LN2 は 二進 Wallace（定数の立っているビット = 配線）、
                ±M + 2^(Pc+Q−1) は 条件反転 + Kogge–Stone、floor は 上位ビットを取るだけ（配線）
            r = xf − k·LN2（k は 2の補数バス → 桁、CSD 配線 + SD 圧縮; 2^-Pc 格子で厳密）
            u = trunc_Wk(r)（canonicalize → 先頭位置 → 枠つき右バレル → 正規化 Wk 桁, Eu）
    小経路: u = x（|N| < 2^Win ≤ 2^Wk なので 切り捨て無し）。u と Eu を small で mux。
    Estrin: 値は全て **符号-大きさ**（canonical）で持つ。信号×信号 = 大きさの二進 Wallace（AND + 全加算器）+ 符号 XOR、
            定数×信号 = 定数の立っているビットの行を 二進 Wallace、和 = 窓つき BF 加算（sticky 桁; ここだけ SD）
            → canonicalize → 正規化 Wk 桁（値は spec の trunc_Wk と一致）。
    出口:   v = 1 + u·q（BF 加算）→ exp: 2^k は 指数バスへの加算（配線）、expm1: (v·2^k) − 1 は 位置 0 を底にした窓
            → finish: v ∓ v·2^-e（定数シフト = 配線 + 二進加減算）→ 丸め（trunc/away/near-even, 桁上げ再正規化）
            → 溢れ ±MAX / 潰れ ±MIN の mux とフラグ（mode は合成時定数 → 三つのネットリスト）。

  ゲート数は 2 通り数える:
    raw  = ライブラリの呼び出し全部（K3 の「最適化していない上界」）
    fold = 定数入力の折り畳み後（AND(1,b)=b, AND(0,b)=0, OR(0,b)=b, OR(1,b)=1, XOR(0,b)=b, XOR(1,b)=¬b,
           定数どうし = 定数）。合成器なら必ずやる範囲だけ。深さは fold の側で 影の評価器 B が測る。

  窓つき BF 加算の正しさ（spec の「厳密和を trunc_Wk」との一致条件）:
    大きい指数側を 上位に置き、小さい指数側を d だけ右へ（落ちた桁の非零 → sticky 桁 ±1 を窓の 1 つ下に置く。
    符号は落とした側の符号 = canonical だから 桁の符号は 数の符号）。整数 Z と 端数 δ（|δ|<1, 符号 σ）について
    Z+δ と Z+σ/2 は 同じ開区間 (Z, Z+1)（σ>0）に居るので、位置 ≥ 0 の切り捨てと 先頭位置は一致する。
    ⟹ 和の先頭位置が 窓の底から Wk−1 以上 上にあれば厳密。G 桁の保護桁で 相殺 G 桁まで許す。
    条件は ホスト側 assert で監視する（破れたら bit 一致が壊れるので テストで見える）。

  Wallace 木の終わり方: 高さ 3 の列を圧縮すると 隣へ桁上げが出て 隣が 3 になる → 層が 幅だけ続く（波及）。
    高さ ≤ 3 で層を止め、最後の 1 層は 高さ 2 の列にも 半圧縮器（HA）を当てる → 全列 ≤ 2 → 加算器 1 個。
"""
import os, sys, math
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)

import gate_bilinear as GB
import gate_exponent as GE
import gate_fast as GF
from gate_bilinear import ZERO, to_sd, from_sd, neg, enc, dec, new_counter
from gate_fast import B, depth_of, or_tree, sd_add2, bin_add_fast
from gate_exponent import bus_const, bus_val, mux_bit, mux_bus, clamp0
from tape_eval import estrin_plan, dy_norm, dy_to_fr
from funcs_spec import Spec


# ============================================================ 折り畳み計数（定数は数えない・深さは B が運ぶ）
def _isc(a): return not isinstance(a, B)          # 定数（ホスト int 0/1）か

class FoldCounter(dict):
    """st: raw（呼び出し全部）と fold（定数折り畳み後）を同時に数える。"""
    def __init__(self):
        super().__init__(); self.raw = new_counter(); self.fold = new_counter()

def _AND(a, b, st):
    st.raw['AND'] += 1
    if _isc(a) and _isc(b): return a & b
    if _isc(a): a, b = b, a
    if _isc(b): return a if b else 0
    st.fold['AND'] += 1; return a & b

def _OR(a, b, st):
    st.raw['OR'] += 1
    if _isc(a) and _isc(b): return a | b
    if _isc(a): a, b = b, a
    if _isc(b): return 1 if b else a
    st.fold['OR'] += 1; return a | b

def _XOR(a, b, st):
    st.raw['XOR'] += 1
    if _isc(a) and _isc(b): return a ^ b
    if _isc(a): a, b = b, a
    if _isc(b):
        if not b: return a
        st.fold['NOT'] += 1; return a ^ 1
    st.fold['XOR'] += 1; return a ^ b

def _NOT(a, st):
    st.raw['NOT'] += 1
    if _isc(a): return a ^ 1
    st.fold['NOT'] += 1; return a ^ 1

for _m in (GB, GE, GF):
    _m.AND, _m.OR, _m.XOR, _m.NOT = _AND, _OR, _XOR, _NOT
AND, OR, XOR, NOT = _AND, _OR, _XOR, _NOT
from gate_bilinear import full_adder, compress3


def sig_digits(v, K):
    """ホスト整数 → K 桁の 信号桁（B 化: 深さ 0）。"""
    return [(B(p), B(n)) for (p, n) in to_sd(v, K)]

def sig_bus(v, EW):
    return [B(b) for b in bus_const(v, EW)]

def bit_v(b): return int(b.v if isinstance(b, B) else b)

def val_digits(digits):
    return sum((bit_v(p) - bit_v(n)) << i for i, (p, n) in enumerate(digits))

def val_bits(bits):
    return sum(bit_v(b) << i for i, b in enumerate(bits))

def val_bus(bus):
    return bus_val([bit_v(b) for b in bus])

def total(c): return sum(c.values())


# ============================================================ 配線だけの部品（ゲート 0）
def to_csd(m):
    """整数 m → CSD 桁列 [(p,n)] 低位から（非零桁が隣接しない・非零数 ≤ ⌈bits/2⌉+1）。合成時定数。"""
    out = []
    while m != 0:
        if m & 1:
            d = 2 - (m & 3)                        # m ≡ 1 (mod 4) → +1, ≡ 3 → −1
            out.append(enc(d)); m -= d
        else:
            out.append(ZERO)
        m >>= 1
    return out or [ZERO]

def zext(bus, EW):
    return list(bus) + [0] * (EW - len(bus))

def sext(bus, EW):
    return list(bus) + [bus[-1]] * (EW - len(bus))

def bus_to_digits(bus):
    """2の補数バス → 符号つき桁列（配線: 上位ビットだけ 重み −2^(w−1) = 負レール）。"""
    w = len(bus)
    return [(bus[i], 0) for i in range(w - 1)] + [(0, bus[w - 1])]

def digits_of_bits(bits, s, st):
    """大きさビット列 + 符号 s → 符号つき桁（正: (b,0) 負: (0,b)）。桁あたり AND 2。"""
    ns = NOT(s, st)
    return [(AND(b, ns, st), AND(b, s, st)) for b in bits]

def trim_zeros(digits):
    """上位の 定数 0 桁（ホスト (0,0)）を落とす（配線）。"""
    n = len(digits)
    while n > 1 and digits[n - 1] == ZERO: n -= 1
    return list(digits[:n])

def trim_bits(bits):
    n = len(bits)
    while n > 1 and _isc(bits[n - 1]) and bits[n - 1] == 0: n -= 1
    return list(bits[:n])

def zero_E(EW):
    """値 0 の指数の約束: −2^(EW−2)。実在する指数は |E| < 2^(EW−2) なので 差が EW ビットに収まり、
       零は加算で必ず「小さい側」になる（全部落ちて sticky 0）。"""
    return bus_const(-(1 << (EW - 2)), EW)


# ============================================================ 圧縮木（定数 0 は列に積まない = 配線）
def _layers(cols, comp3, comp2, st):
    """列ごとの 3:2 圧縮を 全列同時に 1 層ずつ。高さ ≤ 3 で止め、最後の層は 高さ 2 にも comp2 を当てて 全列 ≤ 2 に。"""
    while max(len(c) for c in cols) > 3:
        nxt = [[] for _ in range(len(cols) + 1)]
        for k, c in enumerate(cols):
            i = 0
            while len(c) - i >= 3:
                low, high = comp3(c[i], c[i + 1], c[i + 2], st)
                nxt[k].append(low); nxt[k + 1].append(high)
                i += 3
            nxt[k].extend(c[i:])
        cols = nxt
    if max(len(c) for c in cols) == 3:
        nxt = [[] for _ in range(len(cols) + 1)]
        for k, c in enumerate(cols):
            if len(c) == 3:
                low, high = comp3(c[0], c[1], c[2], st); nxt[k].append(low); nxt[k + 1].append(high)
            elif len(c) == 2:
                low, high = comp2(c[0], c[1], st); nxt[k].append(low); nxt[k + 1].append(high)
            else:
                nxt[k].extend(c)
        cols = nxt
    assert max(len(c) for c in cols) <= 2
    return cols

def sd_sum_cols(rows, st):
    """符号つき桁の 行の和（層別 Wallace + sd_add2）。ホスト定数 ZERO は列に積まない。"""
    width = max(len(x) for x in rows)
    cols = [[] for _ in range(width + 2)]
    for x in rows:
        for i, dg in enumerate(x):
            if dg != ZERO: cols[i].append(dg)
    if max(len(c) for c in cols) == 0: return [ZERO]
    cols = _layers(cols, compress3, lambda a, b, st: compress3(a, b, ZERO, st), st)
    X = [c[0] if len(c) > 0 else ZERO for c in cols]
    Y = [c[1] if len(c) > 1 else ZERO for c in cols]
    if all(y == ZERO for y in Y): return X
    return sd_add2(X, Y, st)

def _ha(a, b, st):
    return XOR(a, b, st), AND(a, b, st)

def bin_sum_cols(rows, st):
    """二進（非負）行の和（層別 Wallace: 全加算器 5 ゲート + 半加算器 2 ゲート → Kogge–Stone 1 個）。"""
    width = max(len(x) for x in rows)
    cols = [[] for _ in range(width + 2)]
    for x in rows:
        for i, b in enumerate(x):
            if not (_isc(b) and b == 0): cols[i].append(b)
    if max(len(c) for c in cols) == 0: return [0]
    cols = _layers(cols, full_adder, _ha, st)
    n = len(cols) + 1
    X = [c[0] if len(c) > 0 else 0 for c in cols] + [0]
    Y = [c[1] if len(c) > 1 else 0 for c in cols] + [0]
    if all(_isc(y) and y == 0 for y in Y): return X
    return bin_add_fast(X, Y, st)

def mul_mag(xm, ym, st):
    """大きさ × 大きさ（二進）: 部分積 AND → 二進 Wallace。戻り ビット列（幅 |x|+|y|+α、上位定数 0 は落とす）。"""
    parts = [[0] * (i + j) + [AND(x, y, st)] for i, x in enumerate(xm) for j, y in enumerate(ym)]
    return trim_bits(bin_sum_cols(parts, st))

def mul_mag_const(xm, m, st):
    """大きさ × 合成時定数 |m|（立っているビットの行 = 配線）→ 二進 Wallace。"""
    a = abs(m); rows = [[0] * j + list(xm) for j in range(a.bit_length()) if (a >> j) & 1]
    if len(rows) == 1: return rows[0]
    return trim_bits(bin_sum_cols(rows, st))


# ============================================================ canonicalize（SD → 大きさ + 符号）・優先エンコーダ・バレル（ビット）
def canon_mag(digits, st):
    """符号つき桁列 → (大きさビット列, 符号)。P−N と N−P を 並列 Kogge–Stone で作り 符号で選ぶ（canonicalize_fast と同値）。
       最上位ビット（|値| < 2^len なので 恒等的 0）は落とす。"""
    w = len(digits) + 1
    P = [(digits[i][0] if i < len(digits) else 0) for i in range(w)]
    N = [(digits[i][1] if i < len(digits) else 0) for i in range(w)]
    D = bin_add_fast(P, [NOT(b, st) for b in N], st, cin=1)                 # P − N
    D2 = bin_add_fast(N, [NOT(b, st) for b in P], st, cin=1)                # N − P
    sign = D[-1]
    return mux_bus(sign, D2, D, st)[:-1], sign

def pe_bits(bits, LW, st):
    """ビット列の 最上位 1 の位置 L（LW ビット）。戻り (L, none, onehot)。suffix-OR 倍化木 + OR 木。"""
    n = len(bits)
    suf = list(bits); k = 1
    while k < n:
        suf = [OR(suf[i], suf[i + k], st) if i + k < n else suf[i] for i in range(n)]
        k <<= 1
    none = NOT(suf[0], st)
    onehot = [AND(bits[i], NOT(suf[i + 1], st), st) if i + 1 < n else bits[i] for i in range(n)]
    L = [or_tree([onehot[i] for i in range(n) if (i >> b) & 1], st) for b in range(LW)]
    return L, none, onehot

def barrel_bits_window(bits, S, Wout, st):
    """bits を S だけ右へ（0 ≤ S < 2^len(S)）、出力は 下位 Wout ビットだけ。上位段（大きい k）から順に施すと
       段 j の後に必要な幅は Wout + (2^j − 1) で済む。戻り (Wout ビット, dropped_nz)。"""
    cur = list(bits); dropped = 0
    for j in range(len(S) - 1, -1, -1):
        sbit = S[j]; k = 1 << j
        need = Wout + (k - 1); n = len(cur)
        dnz = or_tree(cur[:min(k, n)], st)
        dropped = OR(dropped, AND(sbit, dnz, st), st)
        cur = [mux_bit(sbit, cur[i + k] if i + k < n else 0, cur[i], st) for i in range(min(need, n))]
    return cur[:Wout] + [0] * (Wout - len(cur)), dropped

def barrel_left_bits(bits, S, out_width, st):
    """bits を S だけ左へ（×2^S・低位 0 詰め）。固定 out_width（K2: 最大シフトを見込んだ幅）。"""
    cur = [bits[i] if i < len(bits) else 0 for i in range(out_width)]
    for j, sbit in enumerate(S):
        k = 1 << j
        cur = [mux_bit(sbit, cur[i - k] if i - k >= 0 else 0, cur[i], st) for i in range(out_width)]
    return cur

def sat_amount(D, SB, st):
    """符号つきバス D → 右シフト量 min(max(D,0), 2^SB−1) の SB ビット（上位に非零があれば 全 1）。"""
    Dc = clamp0(D, st)
    hi = or_tree(Dc[SB:-1], st) if len(Dc) - 1 > SB else 0
    return [OR(b, hi, st) for b in Dc[:SB]]

def and_tree(bits, st):
    xs = list(bits)
    if not xs: return 1
    while len(xs) > 1:
        nxt = [AND(xs[i], xs[i + 1], st) for i in range(0, len(xs) - 1, 2)]
        if len(xs) % 2: nxt.append(xs[-1])
        xs = nxt
    return xs[0]


# ============================================================ BF（符号-大きさ + 指数バス）
class BF:
    """mag: 大きさビット列（正規化なら 先頭 1 が 位置 W−1）, sign, E: EW ビットバス（LSB の指数）, zero: 値 0 の信号。
       digits(st) で 符号つき桁列を 作る（AND 2/桁、キャッシュ）。"""
    __slots__ = ("mag", "sign", "E", "zero", "_mant")
    def __init__(self, mag, sign, E, zero=0):
        self.mag, self.sign, self.E, self.zero, self._mant = list(mag), sign, list(E), zero, None
    def digits(self, st):
        if self._mant is None: self._mant = digits_of_bits(self.mag, self.sign, st)
        return self._mant
    def value(self):
        m = val_bits(self.mag)
        return dy_norm(-m if bit_v(self.sign) else m, val_bus(self.E))

def const_bf(v, E, W, EW):
    """合成時定数 v·2^E → 正規化 W ビットの BF（ホストビット。|v| < 2^W）。"""
    assert v != 0 and abs(v) < (1 << W)
    a = abs(v); L = a.bit_length() - 1; sh = (W - 1) - L
    a <<= sh
    return BF([(a >> i) & 1 for i in range(W)], 1 if v < 0 else 0, bus_const(E - sh, EW))

def norm_mag(mag, sign, Ebus, W, st, EW=None):
    """大きさ·2^E → 正規化 W ビット BF（値 = spec の trunc_W）。戻り (BF, dropped_nz)。
       先頭位置 L → 枠: 下に W−1 ビットの 0 を足してから L だけ右へ（先頭が 位置 W−1 に来る）。E' = E + L − (W−1)。
       零なら E' = zero_E（加算で必ず「小さい側」になり 全部落ちて sticky 0）。"""
    EW = EW or len(Ebus)
    mag = trim_bits(mag); Wc = len(mag)
    LW = max(1, (Wc - 1).bit_length())
    L, none, _ = pe_bits(mag, LW, st)
    kept, dropped = barrel_bits_window([0] * (W - 1) + mag, L, W, st)
    E1 = bus_add_fast(Ebus, zext(L, EW), st)
    E2 = bus_add_fast(E1, bus_const(-(W - 1), EW), st)
    Eo = mux_bus(none, zero_E(EW), E2, st)
    return BF(kept, AND(sign, NOT(none, st), st), Eo, none), dropped

def bf_normalize(digits, Ebus, W, st, EW=None):
    """任意の（冗長でもよい）桁列·2^E → 正規化 W ビット BF。"""
    mag, s = canon_mag(trim_zeros(digits), st)
    return norm_mag(mag, s, Ebus, W, st, EW)

def bus_add_fast(A, Bb, st, cin=0): return bin_add_fast(A, Bb, st, cin)
def bus_sub_fast(A, Bb, st): return bin_add_fast(A, [NOT(b, st) for b in Bb], st, cin=1)
def bus_lt_fast(A, Bb, st): return bus_sub_fast(A, Bb, st)[-1]
def bus_max_fast(A, Bb, st):
    lt = bus_lt_fast(A, Bb, st); return mux_bus(lt, Bb, A, st)


def bf_add_win(X, Y, W, G, st, EW=None, Dext=0):
    """厳密和の trunc_W と一致する 窓つき加算。X, Y は BF（canonical; 正規化でなくてもよい）。戻り (BF, S, Es)。
       両者を G ビット上げて置き、指数の小さい側を d=|Ex−Ey| だけ右へ（飽和）。落ちた桁 → sticky 桁（index 0）。
       正しさの条件（和の先頭 index ≥ W）は ホスト側で assert。
       Dext > 0: 窓を 下へ Dext 桁伸ばす（相殺が G 桁で収まらない木のため）。d ≤ G+Dext なら 何も落ちず厳密、
       d > G+Dext なら 小さい側の先頭は 大きい側の先頭より 2 桁以上下（G+Dext ≥ len_小 − len_大 + 3 のとき）なので
       相殺は 1 桁以内 → どちらでも 和の先頭 ≥ W。"""
    EW = EW or len(X.E)
    dxy = bus_sub_fast(X.E, Y.E, st); dyx = bus_sub_fast(Y.E, X.E, st)
    lt = dxy[-1]                                               # Ex < Ey
    Ehi = mux_bus(lt, Y.E, X.E, st)
    Wwin = max(len(X.mag), len(Y.mag)) + G + Dext
    SB = max(1, Wwin.bit_length())                             # 2^SB − 1 ≥ Wwin → 全部落とせる
    rows = []; dropped_nz = 0
    for Z, D in ((X, dyx), (Y, dxy)):                          # X は Ey−Ex だけ, Y は Ex−Ey だけ 右へ
        amt = sat_amount(D, SB, st)
        sh, dnz = barrel_bits_window([0] * (G + Dext) + Z.mag, amt, Wwin, st)
        sticky = (AND(dnz, NOT(Z.sign, st), st), AND(dnz, Z.sign, st))
        rows.append([sticky] + digits_of_bits(sh, Z.sign, st))
        dropped_nz |= bit_v(dnz)
    S = sd_add2(rows[0], rows[1], st)
    Es = bus_add_fast(Ehi, bus_const(-G - Dext - 1, EW), st)
    out, _ = bf_normalize(S, Es, W, st, EW)
    v = val_digits(S)                                          # 監視（ゲートではない）: 何か落としたときだけ 先頭 ≥ W が要る
    if dropped_nz and v != 0 and abs(v).bit_length() - 1 < W:
        raise AssertionError(f"bf_add_win: 窓の相殺が G={G}+{Dext} を超えた（先頭 index {abs(v).bit_length()-1} < {W}）")
    return out, S, Es

def bf_mul_sig(X, Y, st, EW=None):
    """信号×信号 → BF（大きさは二進積、符号は XOR）。零なら 指数を 約束値に。"""
    EW = EW or len(X.E)
    z = OR(X.zero, Y.zero, st)
    E = mux_bus(z, zero_E(EW), bus_add_fast(X.E, Y.E, st), st)
    return BF(mul_mag(X.mag, Y.mag, st), XOR(X.sign, Y.sign, st), E, z)

def bf_mul_const(X, m, Em, st, EW=None):
    """信号 × 合成時定数 m·2^Em → BF。"""
    EW = EW or len(X.E)
    assert m != 0
    E = mux_bus(X.zero, zero_E(EW), bus_add_fast(X.E, bus_const(Em, EW), st), st)
    sign = NOT(X.sign, st) if m < 0 else X.sign
    return BF(mul_mag_const(X.mag, m, st), sign, E, X.zero)


def estrin_gate_full(tape_coeffs, u, W, G, st, EW=None):
    """estrin_spec と同じ木を BF で歩く（各節点 正規化 W ビット = trunc_W）。戻り BF。"""
    EW = EW or len(u.E)
    n = max(k for k, m, E in tape_coeffs)
    cd = {k: dy_norm(m, E) for k, m, E in tape_coeffs}
    L, nodes, root = estrin_plan(n)
    pw = [u]
    for l in range(1, L):
        sq = bf_mul_sig(pw[-1], pw[-1], st, EW)
        pw.append(norm_mag(sq.mag, sq.sign, sq.E, W, st, EW)[0])
    val = {}
    for nid, kind, a, b, lvl in nodes:
        if kind == 'leaf':
            m, E = cd.get(a, (0, 0))
            val[nid] = const_bf(m, E, W, EW) if m != 0 else None
            continue
        if nodes[b][1] == 'leaf':
            m, E = cd[nodes[b][2]]
            prod = bf_mul_const(pw[lvl], m, E, st, EW)
        else:
            prod = bf_mul_sig(val[b], pw[lvl], st, EW)
        if a is not None and val[a] is not None:
            val[nid] = bf_add_win(val[a], prod, W, G, st, EW)[0]
        else:
            val[nid] = norm_mag(prod.mag, prod.sign, prod.E, W, st, EW)[0]
    return val[root]


# ============================================================ 出口ユニット（finish_spec / round_spec のゲート版）
def finish_gate(P, e, mode, Wout, Emin, Emax, st, EW):
    """P: BF（W ビット正規化）。mode ∈ lo/hi/near は合成時定数。戻り ((Wout 桁, E バス), (ge, le, sunk))。
       lo: trunc_0(P − P·2^-e) → ge / hi: away_0(P + P·2^(1−e)) → le / near: 最近接（偶数）→ 11。"""
    W = len(P.mag); s = P.sign
    if mode == 'near':
        mag = list(P.mag); Ev = P.E
    else:
        A = [0] * e + list(P.mag)                              # |P|·2^e（単位 2^(E−e)）
        if mode == 'lo':                                       # |P|(2^e − 1)
            Bm = list(P.mag) + [0] * e
            mag = bin_add_fast(A, [NOT(b, st) for b in Bm], st, cin=1)
        else:                                                  # |P|(2^e + 2)
            A = A + [0]; Bm = [0] + list(P.mag) + [0] * e
            mag = bin_add_fast(A, Bm, st)
        Ev = bus_add_fast(P.E, bus_const(-e, EW), st)
    mag = trim_bits(mag)
    Wc = len(mag); LW = max(1, (Wc - 1).bit_length())
    T, none, _ = pe_bits(mag, LW, st)
    shA = bus_add_fast(zext(T, EW), bus_const(-(Wout - 1), EW), st)       # T − (Wout−1)  (> 0: W > Wout)
    shB = bus_sub_fast(bus_const(Emin, EW), Ev, st)                        # Emin − Ev
    sh = bus_max_fast(shA, shB, st)
    SB = max(1, (Wc + 1).bit_length())
    amt = sat_amount(sh, SB, st)
    kept, stk = barrel_bits_window([0] + mag, amt, Wout + 1, st)          # index 0 = 保護ビット g
    gb = kept[0]; qb = kept[1:]
    if mode == 'lo':      inc = 0
    elif mode == 'hi':    inc = OR(gb, stk, st)
    else:                 inc = AND(gb, OR(stk, qb[0], st), st)
    q = bin_add_fast(qb, [0] * Wout, st, cin=inc)
    cout = AND(inc, and_tree(qb, st), st)                                  # 全 1 + 1 → 2^Wout
    qf = mux_bus(cout, [0] * (Wout - 1) + [1], q, st)
    Eo = bus_add_fast(Ev, sh, st, cin=cout)
    collapse = NOT(or_tree(qf, st), st)
    overflow = bus_lt_fast(bus_const(Emax, EW), Eo, st)
    ns = NOT(s, st)
    normal = AND(NOT(collapse, st), NOT(overflow, st), st)
    col_digit = ZERO if mode == 'lo' else (ns, s)                          # 潰れ: lo→0, hi/near→±MIN
    digs = []
    for i in range(Wout):
        p = OR(OR(AND(normal, AND(qf[i], ns, st), st), AND(overflow, ns, st), st),
               AND(collapse, col_digit[0], st) if i == 0 else 0, st)
        n = OR(OR(AND(normal, AND(qf[i], s, st), st), AND(overflow, s, st), st),
               AND(collapse, col_digit[1], st) if i == 0 else 0, st)
        digs.append((p, n))
    Eout = mux_bus(overflow, bus_const(Emax, EW), mux_bus(collapse, bus_const(Emin, EW), Eo, st), st)
    base = {'lo': (1, 0), 'hi': (0, 1), 'near': (1, 1)}[mode]
    col_f = (1, 0) if mode == 'lo' else (0, 1)
    ov_f = (1, 1) if mode == 'hi' else (1, 0)
    ge = OR(OR(AND(normal, base[0], st), AND(collapse, col_f[0], st), st), AND(overflow, ov_f[0], st), st)
    le = OR(OR(AND(normal, base[1], st), AND(collapse, col_f[1], st), st), AND(overflow, ov_f[1], st), st)
    nz_ = NOT(P.zero, st)                                                  # P = 0 → (0, Emin), 000
    digs = [(AND(p, nz_, st), AND(n, nz_, st)) for (p, n) in digs]
    Eout = mux_bus(P.zero, bus_const(Emin, EW), Eout, st)
    return (digs, Eout), (AND(ge, nz_, st), AND(le, nz_, st), 0)


def propagate_gate(fn, fin, s, zero, fc, st):
    """Spec._propagate の exp / expm1 部分（入力フラグ × 符号 × 計算フラグ）。"""
    gi, li, si = fin; cg, cl, cs = fc
    none_in = NOT(OR(OR(gi, li, st), si, st), st)              # 入力フラグ 000 → 計算フラグそのまま
    both = AND(gi, li, st)
    xpos = AND(NOT(s, st), NOT(zero, st), st); xneg = AND(s, NOT(zero, st), st)
    if fn == 'exp':
        dg = OR(AND(gi, xpos, st), AND(li, xneg, st), st); dl = OR(AND(li, xpos, st), AND(gi, xneg, st), st)
        sunk_case = 0
    else:
        dg, dl = gi, li
        sunk_case = OR(si, zero, st)                            # sunk / x=0 → 111
    og = AND(NOT(cl, st), dg, st); ol = AND(NOT(cg, st), dl, st)
    neither = NOT(OR(og, ol, st), st)
    eleven = OR(OR(OR(si, zero, st), both, st), neither, st)    # 11 を返す場合の集合
    ge = mux_bit(none_in, cg, OR(eleven, og, st), st)
    le = mux_bit(none_in, cl, OR(eleven, ol, st), st)
    sk = mux_bit(none_in, cs, OR(sunk_case, AND(neither, cs, st), st), st)   # (1,1,cs) の cs は 計算側（exp では 0）
    return ge, le, sk


# ============================================================ exp / expm1 の本体
class GateExp:
    """f32/f64 構成の exp・expm1 ネットリスト（Spec と同じ定数・同じ木）。"""
    def __init__(self, spec: Spec, G=4, EW=None):
        self.s = spec; self.G = G
        c = spec
        # 指数バス幅（K2: 現れる指数の最悪値から）。小経路の u = x ~ 2^Emin を Estrin の冪 u^(2^(L−1)) まで
        # 上げると 指数は 2^(L−1)·Emin 程度まで下がる（値は捨てられるが バスの比較が壊れてはいけない）。
        # 差が EW ビットに収まるよう |E| < 2^(EW−2) を取る（zero_E も その外側）。
        n = max(k for k, m, E in c.tapes["exp"].coeffs)
        L = estrin_plan(n)[0]
        lo = (1 << (L - 1)) * (-c.Emin + c.Wk) + 2 * c.Wk + 64
        hi = c.Emax + c.Win + c.Pc + c.Q + c.kmax + 8
        self.EW = EW or (max(lo, hi).bit_length() + 2)
        self.KB = c.kmax.bit_length() + 1                       # k の 2の補数幅
        self.AB = (c.XC + c.Win).bit_length()                   # 左バレルの量 E + Win + 1 ∈ [0, XC+Win]
        self.XW = c.XC + c.Pc + 1                               # xf の幅（±2^(XC+Pc) を含む）

    # ---------------------------------------------------------------- 縮小: x → (k バス, u, small, big, zero, s)
    def reduce(self, xd, Eb, st):
        c, EW = self.s, self.EW
        Win, Wk, Pc, Q, XC = c.Win, c.Wk, c.Pc, c.Q, c.XC
        xm, s = canon_mag(xd, st); xm = trim_bits(xm)
        LW = max(1, (len(xm) - 1).bit_length())
        L, none, _ = pe_bits(xm, LW, st)
        ld = bus_add_fast(Eb, zext(L, EW), st)
        small = bus_lt_fast(ld, bus_const(-2, EW), st)                        # ld ≤ −3
        big = NOT(bus_lt_fast(ld, bus_const(XC, EW), st), st)                # ld ≥ XC
        # 中経路: |xf| = |N|·2^(E+Pc) = (|N| << (Pc−Win−1)) << (E+Win+1)
        amt = bus_add_fast(Eb, bus_const(Win + 1, EW), st)[:self.AB]
        xfm = barrel_left_bits([0] * (Pc - Win - 1) + xm, amt, self.XW, st)
        xfm = mux_bus(big, [0] * (XC + Pc) + [1], xfm, st)                    # clamp: |xf| = 2^(XC+Pc)
        # k = floor((xf·INV_LN2 + 2^(Pc+Q−1)) / 2^(Pc+Q)): |xf|·INV_LN2（二進 Wallace）→ ±M + h（条件反転 + KS）→ 上位ビット
        M = mul_mag_const(xfm, c.INV_LN2, st)
        w = max(len(M) + 2, Pc + Q + self.KB + 1)
        Mx = [XOR(M[i] if i < len(M) else 0, s, st) for i in range(w)]        # s なら ~M（+1 は cin）
        D = bin_add_fast(Mx, bus_const(1 << (Pc + Q - 1), w), st, cin=s)
        kbus = D[Pc + Q: Pc + Q + self.KB]
        # r = xf − k·LN2（2^-Pc 格子で厳密）: SD 行（xf の桁 と ∓k·2^j）
        xfd = digits_of_bits(xfm, s, st)
        kd = bus_to_digits(kbus)
        rows = [xfd]
        for j, (p, n) in enumerate(to_csd(c.LN2)):
            if p:   rows.append([ZERO] * j + neg(kd))
            elif n: rows.append([ZERO] * j + list(kd))
        r = trim_zeros(sd_sum_cols(rows, st))
        # u = trunc_Wk(r)（中経路）／ u = x（小経路）: 正規化器の入力を mux
        xdig = digits_of_bits(xm, s, st)
        width = max(len(r), len(xdig))
        rr = r + [ZERO] * (width - len(r)); cc = xdig + [ZERO] * (width - len(xdig))
        nin = GE.mux_digits(small, cc, rr, st)
        Ein = mux_bus(small, Eb, bus_const(-Pc, EW), st)
        u, _ = bf_normalize(nin, Ein, Wk, st, EW)
        kbus_ext = [AND(b, NOT(small, st), st) for b in sext(kbus, EW)]      # 小経路は k = 0
        return kbus_ext, u, small, big, none, s

    # ---------------------------------------------------------------- 前半: 共通（テープ 1 回）
    def core(self, xd, Eb, st, marks=None):
        c, EW, G = self.s, self.EW, self.G
        Wk = c.Wk
        def mark(name):
            if marks is not None: marks[name] = total(st.fold)
        kbus, u, small, big, zero, s = self.reduce(xd, Eb, st)
        mark('reduce')
        q = estrin_gate_full(c.tapes["exp"].coeffs, u, Wk, G, st, EW)
        mark('estrin')
        P = bf_mul_sig(u, q, st, EW)                                          # u·q（厳密）
        one = const_bf(1, 0, Wk, EW)
        v, _, _ = bf_add_win(one, P, Wk, G, st, EW)                           # trunc_Wk(1 + u·q)
        Pe = BF(v.mag, v.sign, bus_add_fast(v.E, kbus, st), v.zero)          # ·2^k は 指数バス
        mark('exp_tail')
        Pm = self.expm1_tail(P, kbus, st)
        mark('expm1_tail')
        return Pe, Pm, zero, s

    # ---------------------------------------------------------------- expm1: trunc_Wk((1 + u·q)·2^k − 1) を厳密に
    def expm1_tail(self, P, kbus, st):
        """k=0: trunc(P)。k≠0: 錨 a = max(k,0) の周りに Wk+6 桁の窓（index i ↔ 絶対位置 a + i − (Wk+4)）を張り
           三行を足す: 2^k（k>0: index Wk+4 の定数 / k<0: index Wk+4−|k| の one-hot）、−1（k>0: index Wk+4−k の one-hot /
           k<0: index Wk+4 の定数）、P·2^k（枠バレル）。結果の先頭は k>0 で [k−3, k]、k<0 で {−1,−2} なので窓に収まる。
           sticky（index 0 の下）: 落ちた合計の符号 =
             k ≤ −(Wk+5): 2^k(1+P) ごと窓の下 → +
             k ≥ Wk+5:   −1 も窓の下。P·2^k の絶対位置 [0, 窓底) に非零ビットがあれば +、無ければ −（−1 が勝つ）
                          → バレルを 2 段に分け（1 段目: 絶対 0 未満を落とす / 2 段目: 窓底まで）2 段目の dropped で判定
             それ以外:     P の落ちた桁の符号 = P の符号"""
        c, EW = self.s, self.EW
        Wk = c.Wk
        Wwin = Wk + 6                                                         # index Wk+4 ↔ 相対位置 0
        A0 = Wk + 4
        ks = kbus[-1]                                                         # k < 0
        kzero = NOT(or_tree(kbus, st), st)
        kpos = AND(NOT(ks, st), NOT(kzero, st), st)
        # 行 1: 錨の位置の定数桁（k>0: +2^k, k<0: −1）
        row1 = [ZERO] * Wwin; row1[A0] = (kpos, ks)
        # 行 2: one-hot（k>0: −1 @ A0−k, k<0: +2^k @ A0+k）: 値 v ∈ [−(Wk+4), −1] ∪ [1, Wk+4] のデコーダ
        vals = [v for v in range(-A0, A0 + 1) if v != 0]
        dec_ = _decoder(kbus, vals, st)
        row2 = [ZERO] * Wwin
        for j in range(A0):
            m = A0 - j                                                        # |k| = m ⟹ index j
            row2[j] = (dec_[-m], dec_[m])
        below = bus_lt_fast(kbus, bus_const(-A0, EW), st)                     # k ≤ −(Wk+5): 全部 窓の下
        onebelow = NOT(bus_lt_fast(kbus, bus_const(A0 + 1, EW), st), st)      # k ≥ Wk+5: −1 が窓の下
        split = NOT(bus_lt_fast(kbus, bus_const(A0, EW), st), st)            # k ≥ Wk+4: 絶対 0 が窓底以下 → 2 段に分ける
        # 行 3: P·2^k を 相対位置 EP + min(k,0) に置く: 枠 = 左に Wwin 前置 → 右へ amt = 2 − EP − min(k,0)
        pre = [0] * Wwin + list(P.mag)
        kneg = mux_bus(ks, kbus, bus_const(0, EW), st)                        # min(k, 0)
        amt_full = bus_sub_fast(bus_sub_fast(bus_const(2, EW), P.E, st), kneg, st)
        SB = max(1, len(pre).bit_length())
        sh, dnz = barrel_bits_window(pre, sat_amount(amt_full, SB, st), Wwin, st)
        stk_p = AND(dnz, NOT(P.sign, st), st); stk_n = AND(dnz, P.sign, st)  # |k| ≤ Wk+4: 落ちるのは P の尾だけ
        # k ≥ Wk+5: −1 も窓の下。落ちた部分 = ±d − 1（d = P·2^k の 窓底未満の大きさ < unit = 2^(k−Wk−4)）は unit を超えうるので
        #   P ≥ 0: ρ = d − 1 ∈ [−1, unit)      → 符号は d vs 1
        #   P < 0: 窓底に −1 桁を余分に置き ρ = (unit−1) − d ∈ (−1, unit) → 符号は d vs unit−1
        #   d の整数ビット（絶対位置 [0, b)）と小数部（絶対位置 < 0）を P のビットに 温度計マスクを掛けて読む:
        #   P のビット i の絶対位置 = i − t（t = −(EP+k)）、窓底 b ↔ i = hi = −EP − (Wk+4)（k に依らない）
        t = bus_sub_fast(bus_const(0, EW), bus_add_fast(P.E, kbus, st), st)
        hi = bus_sub_fast(bus_const(-A0, EW), P.E, st)
        lt, eq = _thermo(t, Wk, st)                                           # i < t, i == t
        lth, _ = _thermo(hi, Wk, st)                                          # i < hi
        f_nz = or_tree([AND(b, l, st) for b, l in zip(P.mag, lt)], st)       # 小数部 ≠ 0
        bit0 = or_tree([AND(b, e, st) for b, e in zip(P.mag, eq)], st)       # 絶対位置 0 のビット
        or1 = or_tree([AND(AND(b, NOT(OR(l, e, st), st), st), h, st) for b, l, e, h in zip(P.mag, lt, eq, lth)], st)
        all0 = AND(and_tree([OR(OR(b, l, st), NOT(h, st), st) for b, l, h in zip(P.mag, lt, lth)], st),
                   AND(NOT(t[-1], st), NOT(bus_lt_fast(bus_const(Wk, EW), hi, st), st), st), st)   # [0,b) が P の範囲内で全 1
        pos_p = OR(or1, AND(bit0, f_nz, st), st); neg_p = AND(NOT(bit0, st), NOT(or1, st), st)     # P ≥ 0: d > 1 / d < 1
        pos_n = NOT(all0, st); neg_n = AND(all0, f_nz, st)                                          # P < 0: d < unit−1 / d > unit−1
        ob_p = mux_bit(P.sign, pos_n, pos_p, st); ob_n = mux_bit(P.sign, neg_n, neg_p, st)
        stk_p = mux_bit(onebelow, ob_p, stk_p, st); stk_n = mux_bit(onebelow, ob_n, stk_n, st)
        stk = (OR(stk_p, below, st), AND(stk_n, NOT(below, st), st))
        row1[0] = (0, AND(onebelow, P.sign, st))                              # P < 0 ∧ k ≥ Wk+5: 窓底の −1 桁
        row3 = [stk] + digits_of_bits(sh, P.sign, st); row1 = [ZERO] + row1; row2 = [ZERO] + row2
        S = sd_sum_cols([row1, row2, row3], st)
        Es = bus_add_fast(mux_bus(ks, bus_const(0, EW), kbus, st), bus_const(-A0 - 1, EW), st)   # max(k,0) − (Wk+5)
        wk, _ = bf_normalize(S, Es, Wk, st, EW)
        w0, _ = norm_mag(P.mag, P.sign, P.E, Wk, st, EW)                      # k = 0: trunc(P)
        out = BF(mux_bus(kzero, w0.mag, wk.mag, st), mux_bit(kzero, w0.sign, wk.sign, st),
                 mux_bus(kzero, w0.E, wk.E, st), mux_bit(kzero, w0.zero, wk.zero, st))
        if not bit_v(kzero):                                                  # 監視（ゲートではない）
            v = val_digits(S)
            assert v != 0 and abs(v).bit_length() - 1 >= Wk, "expm1 窓: 先頭が低すぎる"
        return out

    # ---------------------------------------------------------------- 全体（core 1 回 → 出口 6 個: exp/expm1 × lo/hi/near）
    def run_all(self, xd, Eb, st, fin=(0, 0, 0), fns=('exp', 'expm1'), modes=('lo', 'hi', 'near'), marks=None):
        """戻り {(fn, mode): ((Wout 桁, E バス), (ge, le, sunk))}。marks に dict を渡すと ブロック別 fold ゲート数を書く。"""
        c, EW = self.s, self.EW
        def mark(name):
            if marks is not None: marks[name] = total(st.fold)
        mark('start')
        Pe, Pm, zero, s = self.core(xd, Eb, st, marks)
        out = {}
        for fn in fns:
            P = Pe if fn == 'exp' else Pm
            for mode in modes:
                (digs, Eo), fc = finish_gate(P, c.e[fn], mode, c.Wout, c.Emin, c.Emax, st, EW)
                if fn == 'exp':                                               # x = 0: exp → 1（正規形）
                    one_d = to_sd(1 << (c.Wout - 1), c.Wout); one_E = bus_const(-(c.Wout - 1), EW)
                    digs = GE.mux_digits(zero, one_d, digs, st); Eo = mux_bus(zero, one_E, Eo, st)
                else:                                                         # x = 0: expm1 → (0, Emin)
                    nz_ = NOT(zero, st)
                    digs = [(AND(p, nz_, st), AND(n, nz_, st)) for (p, n) in digs]
                    Eo = mux_bus(zero, bus_const(c.Emin, EW), Eo, st)
                fc = tuple(AND(f, NOT(zero, st), st) for f in fc)
                flags = propagate_gate(fn, fin, s, zero, fc, st)
                out[(fn, mode)] = ((digs, Eo), flags)
                mark(f'finish:{fn}:{mode}')
        return out


def _decoder(bus, values, st):
    """2の補数バス → {v: (bus == v)} の one-hot（2 段: 下位 h ビットと 上位 w−h ビットを 前段デコードして AND）。"""
    w = len(bus); h = w // 2
    nb = [NOT(b, st) for b in bus]
    def lits(lo, hi_, pat):
        return [bus[i] if (pat >> (i - lo)) & 1 else nb[i] for i in range(lo, hi_)]
    lowpat = {}; highpat = {}; out = {}
    for v in values:
        bits = bus_const(v, w)
        lp = sum(bits[i] << i for i in range(h)); hp = sum(bits[i] << (i - h) for i in range(h, w))
        if lp not in lowpat: lowpat[lp] = and_tree(lits(0, h, lp), st)
        if hp not in highpat: highpat[hp] = and_tree(lits(h, w, hp), st)
        out[v] = AND(lowpat[lp], highpat[hp], st)
    return out


def _thermo(t, n, st):
    """符号つきバス t と 位置 i ∈ [0, n) について (i < t, i == t) のビット列。t < 0 → 全 0 / t ≥ n → lt 全 1・eq 全 0。
       t を [0, n] に飽和 → one-hot デコード → suffix-OR（倍化木）。"""
    EW = len(t)
    neg = t[-1]
    ge_n = NOT(bus_lt_fast(t, bus_const(n, EW), st), st)
    LW = n.bit_length() + 1
    tc = mux_bus(neg, bus_const(0, LW), mux_bus(ge_n, bus_const(n, LW), t[:LW], st), st)
    oh = _decoder(tc, list(range(n + 1)), st)
    eq = [oh[i] for i in range(n)]
    suf = [oh[i] for i in range(n + 1)]; k = 1                                # suf[i] = OR_{j ≥ i} oh[j]
    while k < n + 1:
        suf = [OR(suf[i], suf[i + k], st) if i + k < n + 1 else suf[i] for i in range(n + 1)]
        k <<= 1
    lt = [suf[i + 1] for i in range(n)]                                       # i < t ⟺ ∃ j > i: t == j
    return lt, eq


# ============================================================ 差分テスト（Spec と bit 一致）
def eval_point(ge, N, E, fin=(0, 0, 0), digits=None, want_depth=False):
    """1 点を評価。digits を渡せば 冗長表現の入力（値 N と一致していること）。戻り (結果 dict, st, marks, depth)。"""
    c = ge.s
    st = FoldCounter(); marks = {}
    if digits is None:
        xd = sig_digits(N, c.Win)
    else:
        xd = [(B(p), B(n)) for (p, n) in digits]
        assert val_digits(xd) == N
    Eb = sig_bus(E, ge.EW)
    finb = tuple(B(f) for f in fin)
    res = ge.run_all(xd, Eb, st, fin=finb, marks=marks)
    out = {}; dmax = 0
    for key, ((digs, Eo), flags) in res.items():
        out[key] = ((val_digits(digs), val_bus(Eo)), tuple(bit_v(f) for f in flags))
        if want_depth:
            dmax = max(dmax, depth_of(digs), depth_of(Eo), depth_of(list(flags)))
    return out, st, marks, dmax

def spec_point(spec, N, E, fin=(0, 0, 0)):
    out = {}
    for fn in ('exp', 'expm1'):
        for mode in ('lo', 'hi', 'near'):
            out[(fn, mode)] = getattr(spec, fn)((N, E), fin, mode)
    return out

def compare_point(ge, N, E, fin=(0, 0, 0), digits=None):
    got, st, marks, _ = eval_point(ge, N, E, fin, digits)
    want = spec_point(ge.s, N, E, fin)
    bad = {k: (got[k], want[k]) for k in want if got[k] != want[k]}
    return bad, st, marks


_GE = {}

def _worker(args):
    cfg, N, E, fin, digits = args
    ge = _GE[cfg]
    try:
        bad, st, marks = compare_point(ge, N, E, fin, digits)
        return (N, E, fin, digits is not None, bad, None)
    except AssertionError as ex:
        return (N, E, fin, digits is not None, None, str(ex))

def _adversarial_inputs(spec, rng, n_rand=200):
    """小・中・大の経路の境界、k·ln2 の近く、極小・極大、潰れ/溢れの境界、対数一様ランダム。"""
    c = spec
    Win, Wout = c.Win, c.Wout
    pts = []
    def add(N, E):
        if N != 0 and abs(N) < (1 << Win): pts.append((int(N), int(E)))
    for ld in (-4, -3, -2, -1, 0, 1, c.XC - 1, c.XC, c.XC + 1):
        for sgn in (1, -1):
            add(sgn * ((1 << (Win - 1)) | 1), ld - (Win - 1))
            add(sgn * ((1 << Win) - 1), ld - (Win - 1))
            add(sgn * 1, ld)
    from mpmath import mp, mpf
    mp.dps = 60
    for k in list(range(1, 9)) + [c.kmax - 2, c.kmax - 1, c.kmax]:
        for sgn in (1, -1):
            t = mpf(k) * mp.log(2) * sgn
            L = int(mp.floor(mp.log(abs(t), 2))); E = L - (Win - 1)
            M = int(mp.nint(t / mpf(2) ** E))
            for d in (-1, 0, 1): add(M + d, E)
    add(1, c.Emin); add(-1, c.Emin); add((1 << Wout) - 1, c.Emax); add(-((1 << Wout) - 1), c.Emax)
    add(1 << (Win - 1), c.Emin); add(3, c.Emin)
    for t in (mpf(c.Emax + Wout) * mp.log(2), mpf(c.Emax + Wout + 1) * mp.log(2), mpf(-c.Emin) * mp.log(2),
              mpf(-c.Emin + 1) * mp.log(2), mpf(-c.Emin - 1) * mp.log(2)):
        for sgn in (1, -1):
            L = int(mp.floor(mp.log(abs(t), 2))); E = L - (Win - 1); M = int(mp.nint(t / mpf(2) ** E)) * sgn
            for d in (-2, -1, 0, 1, 2): add(M + d, E)
    for _ in range(n_rand):
        M = int(rng.integers(1, 1 << Win)); s = 1 if rng.random() < 0.5 else -1
        ld = int(rng.integers(c.Emin, c.XC + 3))
        add(s * M, ld - (M.bit_length() - 1))
    return pts

def _redundant_digits(rng, N, Win):
    """値 N の 冗長な SD 表現（(1,1) の冗長零・借り/桁上げの混在）。"""
    for _ in range(50):
        d = [0] * Win; rem = N
        for i in range(Win):
            choices = [t for t in (-1, 0, 1) if (rem - t) % 2 == 0]
            t = int(rng.choice(choices)); d[i] = t; rem = (rem - t) // 2
        if rem != 0: continue
        return [(1, 1) if (t == 0 and rng.random() < 0.2) else enc(t) for t in d]
    return None

def self_test(cfg="f32", n_rand=200, seed=20260902, procs=None):
    import time, numpy as np, multiprocessing
    spec = Spec(cfg); ge = GateExp(spec)
    _GE[cfg] = ge
    rng = np.random.default_rng(seed)
    print(f"[{cfg}] EW={ge.EW} KB={ge.KB} AB={ge.AB} XW={ge.XW} G={ge.G} Wk={spec.Wk} Pc={spec.Pc} Q={spec.Q} XC={spec.XC} kmax={spec.kmax}")
    t0 = time.time()
    N0, E0 = (1 << (spec.Win - 1)) | 12345, -(spec.Win - 1)     # x ≈ 1.0015
    got, st, marks, dmax = eval_point(ge, N0, E0, want_depth=True)
    want = spec_point(spec, N0, E0)
    print(f"  1 点 {time.time()-t0:.1f}s  raw={total(st.raw)}  fold={total(st.fold)}  depth={dmax}")
    prev = 0
    for k in marks:
        print(f"    {k:18s} {marks[k]-prev:8d}"); prev = marks[k]
    for k in want:
        print(f"    {k}: gate={got[k]} spec={want[k]} {'OK' if got[k]==want[k] else 'MISMATCH'}")
    pts = _adversarial_inputs(spec, rng, n_rand)
    jobs = [(cfg, N, E, (0, 0, 0), None) for (N, E) in pts]
    for (N, E) in pts[:12]:
        for fin in ((1,0,0),(0,1,0),(1,1,0),(0,0,1),(1,0,1),(0,1,1),(1,1,1)):
            jobs.append((cfg, N, E, fin, None))
    for (N, E) in pts[::5]:
        d = _redundant_digits(rng, N, spec.Win)
        if d is not None: jobs.append((cfg, N, E, (0, 0, 0), d))
    jobs.append((cfg, 0, 0, (0, 0, 0), None)); jobs.append((cfg, 0, -5, (1, 0, 0), None))
    jobs.append((cfg, 0, 3, (0, 0, 1), [(1, 1)] * spec.Win))
    t0 = time.time()
    ctx = multiprocessing.get_context("fork")
    with ctx.Pool(procs or min(24, len(jobs))) as pool:
        res = pool.map(_worker, jobs, chunksize=1)
    nbad = nerr = 0
    for (N, E, fin, red, bad, err) in res:
        if err is not None:
            nerr += 1
            if nerr <= 10: print(f"  ASSERT x=({N},{E}) fin={fin} red={red}: {err}")
        elif bad:
            nbad += 1
            if nbad <= 10: print(f"  MISMATCH x=({N},{E}) fin={fin} red={red}: {bad}")
    print(f"  {len(jobs)} 点 {time.time()-t0:.1f}s  不一致 {nbad}  assert {nerr}")
    return nbad == 0 and nerr == 0


def test_expm1_tail(cfg="f32", n_rand=300, seed=1):
    """expm1_tail 単体: P（正規化 Wk ビット, 任意の指数）と k を直接注入し、Fraction で厳密な trunc_Wk((1+P)·2^k − 1) と比べる。
       狙って踏む経路: k = 0 / |k| ≤ Wk+4 / k ≤ −(Wk+5) / k ≥ Wk+5 で P·2^k の下位整数ビットが 全零（−1 の借り）と 非零。"""
    import numpy as np
    from fractions import Fraction as Fr
    spec = Spec(cfg); ge = GateExp(spec); Wk, EW = spec.Wk, ge.EW
    rng = np.random.default_rng(seed)
    cases = []
    ks = sorted(set([0, 1, 2, 3, -1, -2, -3, Wk + 2, Wk + 3, Wk + 4, Wk + 5, Wk + 6, Wk + 7, Wk + 9, Wk + 20, Wk + 60,
                     -(Wk + 3), -(Wk + 4), -(Wk + 5), -(Wk + 6), -(Wk + 30), spec.kmax, -spec.kmax]))
    for k in ks:
        for _ in range(n_rand // len(ks) + 1):
            # |P| < 0.42, 先頭位置 ≤ −2。指数は 深いものも（u が小さい場合）
            lead = int(rng.integers(-1 - Wk * 2, -1)) if rng.random() < 0.5 else int(rng.integers(-6, -1))
            mag = (1 << (Wk - 1)) | int(rng.integers(0, 1 << (Wk - 1)))
            if lead == -2 and mag >= int(0.84 * (1 << Wk)): mag = (1 << (Wk - 1)) | int(rng.integers(0, 1 << (Wk - 2)))
            E = lead - (Wk - 1)
            if rng.random() < 0.4:                                      # 下位を 0 で埋める（P·2^k の整数部下位が全零になりやすい）
                z = int(rng.integers(1, Wk - 2)); mag = (mag >> z) << z
            sgn = -1 if rng.random() < 0.5 else 1
            cases.append((sgn * mag, E, k))
        cases.append(((1 << (Wk - 1)), -(Wk - 1) - 3, k))                 # P = +2^-3（P·2^k は 2 のべき）
        cases.append((-(1 << (Wk - 1)), -(Wk - 1) - 3, k))
    bad = 0
    for (m, E, k) in cases:
        st = FoldCounter()
        P = BF([B((abs(m) >> i) & 1) for i in range(Wk)], B(1 if m < 0 else 0), sig_bus(E, EW), B(0))
        out = ge.expm1_tail(P, sig_bus(k, EW), st)
        got = (val_bits(out.mag) * (-1 if bit_v(out.sign) else 1), val_bus(out.E))
        R = (1 + Fr(m) * Fr(2) ** E) * Fr(2) ** k - 1
        if R == 0:
            want = (0, -(1 << (EW - 2)))
        else:
            a = abs(R); L = a.numerator.bit_length() - a.denominator.bit_length()
            if Fr(2) ** L > a: L -= 1                                   # 2^L ≤ |R| < 2^(L+1)
            Eo = L - (Wk - 1); q = int(a / Fr(2) ** Eo)                 # trunc
            want = ((-q if R < 0 else q), Eo)
        if got != want:
            bad += 1
            if bad <= 10: print(f"  MISMATCH P=({m},{E}) k={k}: gate={got} exact={want}")
    print(f"  expm1_tail[{cfg}] {len(cases)} 例  不一致 {bad}")
    return bad == 0


if __name__ == "__main__":
    cfg = sys.argv[1] if len(sys.argv) > 1 else "f32"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    ok = test_expm1_tail(cfg) and self_test(cfg, n)
    print("PASS" if ok else "FAIL")
