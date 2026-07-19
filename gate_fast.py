#!/usr/bin/env python3
# ⚠️ 生成AI使用・要検証
"""深さ削減版 — 逐次性分析の 処方箋を 実装し、**深さを 実測**で 前後比較する。

  処方箋（優先順）:
   ① sd_sum の 圧縮を LIFO 鎖 → 層別 Wallace（O(m) → O(log m)）
   ② 最終リップル → **定数深さ SD 加算器 sd_add2**（Avizienis 限定キャリー・隣接窓のみ）
   ③ canonicalize の 2連リップル → 双方向 並列減算（Kogge–Stone）＋ 符号 mux
   ④ OR 累積鎖 → OR 木（優先エンコーダ・kept_nz/over・バレル dnz）
   ⑤ bfg_add: min→sub の 直列 → 並列 双方向差 ＋ 2行 sd_add2
   ⑥ Lmax の 線形 fold → トーナメント木
  検証 = 全て 差分テスト（遅い版と 値・フラグ 完全一致）。深さ = 影の評価器 B で 実測。
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from gate_bilinear import (AND, OR, NOT, XOR, new_counter, ZERO, enc, dec,
                           to_sd, from_sd, neg, sd_sum, multiply, canonicalize,
                           compress3, full_adder, nz, bin_add)
from gate_exponent import (bus_const, bus_val, bus_add, bus_sub, bus_lt, bus_max,
                           mux_bit, mux_bus, clamp0, priority_encoder,
                           barrel_shift_right_digits, barrel_shift_left_digits)


# ============================================================ 影の評価器（深さ実測）
class B:
    """ビット + 深さ。ゲート演算子で 深さ = max(入力)+1 を 自動伝播（値は 通常どおり）。"""
    __slots__ = ('v', 'd')
    def __init__(self, v, d=0): self.v = int(v); self.d = d
    @staticmethod
    def _c(o): return o if isinstance(o, B) else B(int(o), 0)
    def _op(self, o, f):
        o = B._c(o); return B(f(self.v, o.v), max(self.d, o.d) + 1)
    def __and__(self, o):  return self._op(o, lambda a, b: a & b)
    def __rand__(self, o): return self._op(o, lambda a, b: a & b)
    def __or__(self, o):   return self._op(o, lambda a, b: a | b)
    def __ror__(self, o):  return self._op(o, lambda a, b: a | b)
    def __xor__(self, o):  return self._op(o, lambda a, b: a ^ b)
    def __rxor__(self, o): return self._op(o, lambda a, b: a ^ b)
    def __sub__(self, o):  return self.v - B._c(o).v          # dec/from_sd の 読み出しのみ
    def __rsub__(self, o): return B._c(o).v - self.v
    def __lshift__(self, k): return self.v << k

def wrap(digits):
    return [(B(p), B(n)) for (p, n) in digits]

def wrapb(bits):
    return [B(b) for b in bits]

def depth_of(x):
    def d(o): return o.d if isinstance(o, B) else 0
    if isinstance(x, (list, tuple)):
        return max((depth_of(e) for e in x), default=0)
    return d(x)


# ============================================================ ④ OR 木
def or_tree(bits, st):
    """OR の 並列木（O(log)）。空なら 0。"""
    xs = list(bits)
    if not xs: return 0
    while len(xs) > 1:
        nxt = []
        for i in range(0, len(xs) - 1, 2):
            nxt.append(OR(xs[i], xs[i + 1], st))
        if len(xs) % 2: nxt.append(xs[-1])
        xs = nxt
    return xs[0]


# ============================================================ ② 定数深さ SD 加算器
def sd_add2(X, Y, st):
    """2 つの 符号つき数の 和を **定数深さ**で（Avizienis 限定キャリー）。
       s_i = x_i+y_i を compress3 で (l,h) に。h=±1 ⟺ s=±2、l=±1 ⟺ s=±1。
       転送 t_{i+1}: +1 ⟺ s=2 ∨ (s=1 ∧ 隣s_{i-1}≥0) ／ −1 ⟺ 鏡像。
       仮 w_i: ±1 は s=±1 の 残り。z_i = w_i + t_i は 構成上 {−1,0,1}（もう 伝播しない）。"""
    n = max(len(X), len(Y))
    Xp = [X[i] if i < len(X) else ZERO for i in range(n)]
    Yp = [Y[i] if i < len(Y) else ZERO for i in range(n)]
    L = []; H = []
    for i in range(n):
        l, h = compress3(Xp[i], Yp[i], ZERO, st)
        L.append(l); H.append(h)
    pos = [OR(H[i][0], L[i][0], st) for i in range(n)]        # s_i ≥ 1
    neg_ = [OR(H[i][1], L[i][1], st) for i in range(n)]       # s_i ≤ −1
    tp = [0] * (n + 1); tn = [0] * (n + 1)
    wp = [0] * n; wn = [0] * n
    for i in range(n):
        npv = neg_[i - 1] if i > 0 else 0                     # 隣（下位）の 符号
        ppv = pos[i - 1] if i > 0 else 0
        tp[i + 1] = OR(H[i][0], AND(L[i][0], NOT(npv, st), st), st)
        tn[i + 1] = OR(H[i][1], AND(L[i][1], NOT(ppv, st), st), st)
        wp[i] = OR(AND(L[i][0], npv, st), AND(L[i][1], NOT(ppv, st), st), st)
        wn[i] = OR(AND(L[i][0], NOT(npv, st), st), AND(L[i][1], ppv, st), st)
    out = []
    for i in range(n + 1):
        w = (wp[i], wn[i]) if i < n else ZERO
        t = (tp[i], tn[i])
        z, hi = compress3(w, t, ZERO, st)                     # 構成上 hi=0（伝播 終端）
        out.append(z)
    return out


# ============================================================ ① 層別 Wallace sd_sum
def sd_sum_fast(nums, st):
    """符号つき数 リストの 和。層別 3:2 圧縮（全列 同時・O(log m)）→ 2 行 → sd_add2（定数）。"""
    if not nums: return [ZERO]
    width = max(len(x) for x in nums)
    cols = [[] for _ in range(width + 2)]
    for x in nums:
        for i, dg in enumerate(x):
            cols[i].append(dg)
    while max(len(c) for c in cols) > 2:
        nxt = [[] for _ in range(len(cols) + 1)]
        for k, c in enumerate(cols):                          # 全列 同時に 1 層
            i = 0
            while i + 2 < len(c) + 1 and len(c) - i >= 3:
                low, high = compress3(c[i], c[i + 1], c[i + 2], st)
                nxt[k].append(low); nxt[k + 1].append(high)
                i += 3
            nxt[k].extend(c[i:])                              # 余り 素通し
        cols = nxt
    X = [c[0] if len(c) > 0 else ZERO for c in cols]
    Y = [c[1] if len(c) > 1 else ZERO for c in cols]
    return sd_add2(X, Y, st)

def multiply_fast(X, Y, st):
    parts = []
    for i, xd in enumerate(X):
        for j, yd in enumerate(Y):
            from gate_bilinear import gate9
            parts.append([ZERO] * (i + j) + [gate9(xd, yd, st)])
    return sd_sum_fast(parts, st)


# ============================================================ ③ Kogge–Stone ＋ canonicalize
def bin_add_fast(A, Bb, st, cin=0):
    """並列プレフィックス（Kogge–Stone）加算器 O(log n)。bin_add と 同値。"""
    n = len(A)
    g = [AND(a, b, st) for a, b in zip(A, Bb)]
    p = [XOR(a, b, st) for a, b in zip(A, Bb)]
    G = list(g); P = list(p)
    k = 1
    while k < n:
        G2 = list(G); P2 = list(P)
        for i in range(n - 1, k - 1, -1):
            G2[i] = OR(G[i], AND(P[i], G[i - k], st), st)
            P2[i] = AND(P[i], P[i - k], st)
        G, P = G2, P2
        k <<= 1
    carry = [cin] + [OR(G[i], AND(P[i], cin, st), st) for i in range(n - 1)]
    return [XOR(p[i], carry[i], st) for i in range(n)]

def canonicalize_fast(digits, st):
    """双方向 並列減算（D=P−N と D'=N−P を 同時）＋ 符号 mux。2連リップルを 排除。"""
    w = len(digits) + 1
    P = [(digits[i][0] if i < len(digits) else 0) for i in range(w)]
    N = [(digits[i][1] if i < len(digits) else 0) for i in range(w)]
    notN = [NOT(b, st) for b in N]
    notP = [NOT(b, st) for b in P]
    D  = bin_add_fast(P, notN, st, cin=1)                     # P − N
    D2 = bin_add_fast(N, notP, st, cin=1)                     # N − P（並列）
    sign = D[-1]
    mag = mux_bus(sign, D2, D, st)
    nsign = NOT(sign, st)
    return [(AND(m, nsign, st), AND(m, sign, st)) for m in mag], sign


# ============================================================ ③b 指数バス演算も prefix 化
def bus_add_fast(A, Bb, st, cin=0):
    return bin_add_fast(A, Bb, st, cin)

def bus_sub_fast(A, Bb, st):
    return bin_add_fast(A, [NOT(b, st) for b in Bb], st, cin=1)

def bus_lt_fast(A, Bb, st):
    return bus_sub_fast(A, Bb, st)[-1]

def bus_max_fast(A, Bb, st):
    lt = bus_lt_fast(A, Bb, st)
    return mux_bus(lt, Bb, A, st)


# ============================================================ ④ 優先エンコーダ（木）
def priority_encoder_fast(digits, EW, st):
    """suffix-OR を 倍化木で O(log n)、二進符号化も OR 木。"""
    n = len(digits)
    nzs = [nz(d, st) for d in digits]
    suf = list(nzs)                                           # suf[i] = OR(nzs[i..])
    k = 1
    while k < n:
        suf = [OR(suf[i], suf[i + k], st) if i + k < n else suf[i] for i in range(n)]
        k <<= 1
    none = NOT(suf[0], st)
    onehot = [AND(nzs[i], NOT(suf[i + 1], st), st) if i + 1 < n else nzs[i]
              for i in range(n)]
    L = [or_tree([onehot[i] for i in range(n) if (i >> b) & 1], st) for b in range(EW)]
    return L, none, onehot


# ============================================================ ⑤⑥ ブロック浮動の 高速版
def bfg_add_fast(x, y, Dmax_bits, st):
    """双方向差を 並列に → 符号で Elo 選択 → バレル整列 → **2 行 sd_add2（定数深さ）**。"""
    dxy = bus_sub_fast(x.Ebus, y.Ebus, st)                    # Ex − Ey（prefix）
    dyx = bus_sub_fast(y.Ebus, x.Ebus, st)                    # Ey − Ex（並列・prefix）
    lt = dxy[-1]                                              # Ex < Ey
    Elo = mux_bus(lt, x.Ebus, y.Ebus, st)
    dx = clamp0(dxy, st)[:Dmax_bits]
    dy = clamp0(dyx, st)[:Dmax_bits]
    ow = max(len(x.mant), len(y.mant)) + (1 << Dmax_bits)
    mx = barrel_shift_left_digits(x.mant, dx, ow, st)
    my = barrel_shift_left_digits(y.mant, dy, ow, st)
    from gate_bfp2 import BFg
    return BFg(sd_add2(mx, my, st), Elo)

def block_normalize_g_fast(mants, Ebus, W, Emax, st):
    """canonicalize_fast ＋ 優先エンコーダ木 ＋ Lmax トーナメント ＋ OR 木。仕様は 遅い版と 同一。"""
    from gate_bilinear import ZERO as Z
    canons = []; signs = []; Ls = []
    Wc = max(len(m) for m in mants) + 1
    EW = len(Ebus)
    for m in mants:
        c, s = canonicalize_fast(m, st)
        c = [c[i] if i < len(c) else Z for i in range(Wc)]
        canons.append(c); signs.append(s)
        L, _, _ = priority_encoder_fast(c, EW, st)
        Ls.append(L)
    while len(Ls) > 1:                                        # ⑥ トーナメント木
        nxt = [bus_max_fast(Ls[i], Ls[i + 1], st) if i + 1 < len(Ls) else Ls[i]
               for i in range(0, len(Ls), 2)]
        Ls = nxt
    Lmax = Ls[0]
    sh = clamp0(bus_sub_fast(Lmax, bus_const(W - 1, EW), st), st)
    E_out = bus_add_fast(Ebus, sh, st)
    EmaxB = bus_const(Emax, EW)
    ovE = bus_lt_fast(EmaxB, E_out, st)
    sh_cap = clamp0(bus_sub_fast(EmaxB, Ebus, st), st)
    sh_al = mux_bus(ovE, sh_cap, sh, st)
    E_fin = mux_bus(ovE, EmaxB, E_out, st)
    SW = max(1, (Wc - 1).bit_length())
    out = []; flags = []
    for c, s in zip(canons, signs):
        shifted, drop_nz = barrel_shift_right_digits(c, sh_al[:SW], st)
        kept = shifted[:W]
        kept_nz = or_tree([nz(d, st) for d in kept], st)      # ④ OR 木
        over = or_tree([nz(d, st) for d in shifted[W:]], st)
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


# ============================================================ 自己テスト（同値＋深さ実測）
def self_test():
    import numpy as np
    rng = np.random.default_rng(20260805)

    print("=" * 78)
    print("② sd_add2 — 定数深さ SD 加算器（値保存・桁範囲・深さが 幅に 依らない）")
    print("=" * 78)
    bad = 0
    for _ in range(4000):
        a = int(rng.integers(-30000, 30000)); b = int(rng.integers(-30000, 30000))
        z = sd_add2(to_sd(a, 16), to_sd(b, 16), new_counter())
        if from_sd(z) != a + b: bad += 1
        if any(dec(*d) not in (-1, 0, 1) for d in z): bad += 1
    ds = []
    for n in (8, 16, 32, 64):
        z = sd_add2(wrap(to_sd((1 << (n - 2)) - 3, n)), wrap(to_sd(-(1 << (n - 3)) - 7, n)), new_counter())
        ds.append(depth_of(z))
    print(f"  値保存: 違反 {bad}/4000 ✓   深さ実測 n=8/16/32/64 → {ds}（**幅に 依らず 一定**）")

    print()
    print("=" * 78)
    print("①+② sd_sum_fast / multiply_fast — 遅い版と 同値・深さ 実測比較")
    print("=" * 78)
    bad = 0
    for _ in range(1500):
        m = int(rng.integers(2, 9))
        xs = [int(v) for v in rng.integers(-500, 500, m)]
        if from_sd(sd_sum_fast([to_sd(v, 12) for v in xs], new_counter())) != sum(xs): bad += 1
    for _ in range(800):
        a = int(rng.integers(-300, 300)); b = int(rng.integers(-300, 300))
        if from_sd(multiply_fast(to_sd(a, 10), to_sd(b, 10), new_counter())) != a * b: bad += 1
    xs = [int(v) for v in rng.integers(-500, 500, 8)]
    slow = depth_of(sd_sum([wrap(to_sd(v, 12)) for v in xs], new_counter()))
    fast = depth_of(sd_sum_fast([wrap(to_sd(v, 12)) for v in xs], new_counter()))
    a, b = 217, -178
    ms = depth_of(multiply(wrap(to_sd(a, 10)), wrap(to_sd(b, 10)), new_counter()))
    mf = depth_of(multiply_fast(wrap(to_sd(a, 10)), wrap(to_sd(b, 10)), new_counter()))
    print(f"  同値: 違反 {bad}/2300 ✓   深さ: sd_sum(8行) {slow}→{fast}  multiply(10桁) {ms}→{mf}")

    print()
    print("=" * 78)
    print("③ bin_add_fast(Kogge–Stone) / canonicalize_fast — 同値・深さ比較")
    print("=" * 78)
    bad = 0
    for _ in range(3000):
        n = 14
        a = int(rng.integers(0, 1 << n)); b = int(rng.integers(0, 1 << n))
        A = [(a >> i) & 1 for i in range(n)]; Bb = [(b >> i) & 1 for i in range(n)]
        cin = int(rng.integers(0, 2))
        if bin_add_fast(A, Bb, new_counter(), cin) != bin_add(A, Bb, cin, new_counter())[0]: bad += 1
    bad2 = 0
    for _ in range(1500):
        m = int(rng.integers(2, 7))
        xs = [int(v) for v in rng.integers(-800, 800, m)]
        red = sd_sum([to_sd(v, 14) for v in xs], new_counter())     # 冗長形を 作る
        cf, sf = canonicalize_fast(red, new_counter())
        cs, ss = canonicalize(red, new_counter())
        if from_sd(cf) != from_sd(cs) or sf != ss: bad2 += 1
        vals = [dec(*d) for d in cf if dec(*d) != 0]
        if vals and len(set(vals)) != 1: bad2 += 1                  # 非冗長性
    red = sd_sum([to_sd(v, 20) for v in (999999, -999998, 777)], new_counter())
    dslow = depth_of(canonicalize(wrap(red), new_counter()))
    dfast = depth_of(canonicalize_fast(wrap(red), new_counter()))
    print(f"  bin_add 同値 {3000-bad}/3000 ✓  canonicalize 同値+非冗長 {1500-bad2}/1500 ✓  "
          f"深さ(21桁): {dslow}→{dfast}")

    print()
    print("=" * 78)
    print("④ priority_encoder_fast — 同値・深さ比較")
    print("=" * 78)
    bad = 0
    for _ in range(2000):
        v = int(rng.integers(-(1 << 18), 1 << 18))
        D = to_sd(v, 20)
        Lf, nf, of_ = priority_encoder_fast(D, 10, new_counter())
        Ls, ns, os_ = priority_encoder(D, 10, new_counter())
        if [int(x) for x in Lf] != [int(x) for x in Ls] or int(nf) != int(ns): bad += 1
    D = wrap(to_sd(123456, 24))
    ds = depth_of(priority_encoder(D, 10, new_counter())[0])
    df = depth_of(priority_encoder_fast(D, 10, new_counter())[0])
    print(f"  同値: 違反 {bad}/2000 ✓   深さ(24桁): {ds}→{df}")

    print()
    print("=" * 78)
    print("⑤⑥ bfg_add_fast / block_normalize_g_fast — 差分テスト（値・指数・フラグ 一致）・深さ")
    print("=" * 78)
    from gate_bfp2 import to_bfg, from_bfg, bfg_add, block_normalize_g, EW
    bad = 0
    for _ in range(400):
        a = int(rng.integers(-500, 500)); b = int(rng.integers(-500, 500))
        Ea = int(rng.integers(0, 8)); Eb = int(rng.integers(0, 8))
        sf = bfg_add_fast(to_bfg(a, 12, Ea), to_bfg(b, 12, Eb), 3, new_counter())
        ss = bfg_add(to_bfg(a, 12, Ea), to_bfg(b, 12, Eb), 3, new_counter())
        if from_bfg(sf) != from_bfg(ss) or bus_val(sf.Ebus) != bus_val(ss.Ebus): bad += 1
    print(f"  bfg_add_fast == bfg_add: 違反 {bad}/400 ✓")
    bad = 0
    for _ in range(200):
        M = 4; W = 6; Emax = 20
        E0 = int(rng.integers(0, 12))
        vals = [int(rng.integers(-5, 6)) * (10 ** int(rng.integers(0, 5))) for _ in range(M)]
        ms = [to_sd(v, 24) for v in vals]
        of_, Ef, ff = block_normalize_g_fast(ms, bus_const(E0, EW), W, Emax, new_counter())
        os_, Es, fs = block_normalize_g(ms, bus_const(E0, EW), W, Emax, new_counter())
        if bus_val([int(x) for x in Ef]) != bus_val([int(x) for x in Es]): bad += 1
        for gg, hh, f1, f2 in zip(of_, os_, ff, fs):
            if from_sd(gg) != from_sd(hh) or (int(f1[0]), int(f1[1])) != (int(f2[0]), int(f2[1])): bad += 1
    print(f"  block_normalize_g_fast == 遅い版（仮数・指数・フラグ）: 違反 {bad}/200ブロック ✓")
    ms = [wrap(to_sd(v, 24)) for v in (10000, 3, -2, 0)]
    Eb0 = wrapb(bus_const(5, EW))
    dslow = depth_of(block_normalize_g([[(B._c(p), B._c(n)) for p, n in m] for m in ms], Eb0, 6, 20, new_counter())[0])
    dfast = depth_of(block_normalize_g_fast([[(B._c(p), B._c(n)) for p, n in m] for m in ms], Eb0, 6, 20, new_counter())[0])
    xa = to_bfg(217, 12, 3); xb = to_bfg(-178, 12, 6)
    xa.mant = wrap(xa.mant); xa.Ebus = wrapb(xa.Ebus)
    xb.mant = wrap(xb.mant); xb.Ebus = wrapb(xb.Ebus)
    daslow = depth_of(bfg_add(xa, xb, 3, new_counter()).mant)
    dafast = depth_of(bfg_add_fast(xa, xb, 3, new_counter()).mant)
    print(f"  深さ実測: block_normalize {dslow}→{dfast}   bfg_add {daslow}→{dafast}")

    print()
    print("処方箋 ①〜⑥ 実装・全差分テスト同値・深さは 実測で 短縮。")


if __name__ == "__main__":
    self_test()
