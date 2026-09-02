#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証（research ブランチ・未監査）
"""gate_logroot — Level 2: log / sqrt / rsqrt を **端から端までゲート**で（gate_funcs の部品を使い回す）。
   funcs_spec.Spec（Level 1）と bit 一致することを差分テストで確かめ、ゲート数と深さを実測する。

  配線（f32: Win=Wout=24, Wk=40, Pc=56 / f64: 53, 64, 80）:
    入力 x = (Win 桁 符号つき桁, E バス) ── canonicalize(R1) → 正規化 Win ビット m（先頭 1 が位置 Win−1）, Ex = E + L
    h:  log  = (m の上位 8 ビット ≥ 181)（9 ビット比較器 1 個）   … x = m·2^Ex, m ∈ [1,2)
        根   = Ex の最下位ビット（偶奇。2の補数バスなのでそのまま）
    t = m·2^-h − 1: m の桁列に 位置 Win−1+h の −1 桁を足して canonicalize（厳密: Win+1 ビット ≤ Wk）
        → 正規化 Wk ビット u。指数は −(Win−1+h)（h で mux）。x = 0 なら t := 0（Estrin を定義域の外で走らせない）。
    Ex' = Ex + h
    log:  q = Estrin(log テープ, u) → P = trunc_Wk(Ex'·LN2·2^-Pc + u·q)
          Ex'·LN2 は |Ex'|（条件反転 + 1）× 定数 LN2（立っているビットの行 = 配線 → 二進 Wallace）、符号 = Ex' の符号、
          Ex' = 0 なら 零 BF（指数 zero_E）。和は 窓つき BF 加算（|Ex' ln2| ≥ 0.69 > 0.35 ≥ |u q| なので 相殺は 2 桁以内）。
    sqrt / rsqrt: 2 本のテープ p0/p1（同じ次数）は **係数を h で選ぶ**（定数ビットの mux は 配線 + NOT）。
          係数が信号になるので 葉の積は 定数積でなく 信号積（Wallace）になる … 木を 2 本並べる 2 倍より安い。
          p1 テープ（t ∈ [−1/2, 0)）の係数は |c_k| > 1 で Wk より広い（f32 rsqrt 45 / f64 sqrt 69・rsqrt 76 ビット）:
          spec が厳密に使うので ゲートも そのままの幅で持つ（K1）。
          指数は v.E ± (Ex' >> 1)（Ex' は偶数なので算術シフト = 配線）。
    出口: finish_gate（gate_funcs と同じ）× lo/hi/near。x = 0: log → −MAX (1,0,0)、根 → (0, Emin) 000。
          x < 0: 値は |x| で計算し フラグ 111（spec と同じ: 実数では何も主張しない）。
    フラグ伝播: Spec._propagate の log / sqrt / rsqrt 部分をそのままゲートに（x ≷ 1 は Ex の符号と m の下位ビットで）。

  根のテープの 窓（exp との違い）: exp / log / sin / cos の Estrin は 節点多項式の根が定義域の外にあり 相殺 < 1 bit
    （cancel_report: 格子上の最大 0.3 / 0.8 / 0.2 / 0.5 bit）なので 固定 G=4 の窓で足りる。sqrt / rsqrt は p1 で
    節点の根が定義域の内側（f32 1+2 個、f64 4+4 個）にあり 相殺が非有界、p0 も根は無いが f64 rsqrt_p0 で 4.4 bit > G。
    → estrin_gate_mux は 窓を 下へ Dext = W + extra + 3 − G 桁伸ばす（bf_add_win の Dext）。落ちるのは 先頭位置の差が
    ≥ 2 のときだけになり（相殺 ≤ 1 bit）、落ちなければ厳密。代金は 節点の加算 1 回が f32 13.9k → 22.1k・f64 23.3k → 36.5k
    （W + 2W ビット）、関数全体で f32 sqrt +25% / rsqrt +33%、f64 +23% / +25%。
  ゲート数の数え方は gate_funcs.py 冒頭と同じ。R4: 係数 mux の Estrin は 単体テスト（test_estrin_mux）で u と h を直接
    注入し、節点多項式の根 ± 1 ulp（相殺が最も深い入力）も踏む。
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from gate_funcs import (FoldCounter, AND, OR, XOR, NOT, B, ZERO, to_sd, sig_digits, sig_bus, bit_v,
                        val_digits, val_bits, val_bus, total, zext, sext, bus_to_digits, trim_bits, zero_E,
                        mul_mag_const, canon_mag, or_tree, BF, const_bf, norm_mag, bus_add_fast, bus_sub_fast,
                        bus_lt_fast, bf_add_win, bf_mul_sig, estrin_gate_full, finish_gate, GateExp,
                        _redundant_digits, bus_const, mux_bit, mux_bus, bin_add_fast, depth_of)
import gate_exponent as GE
from tape_eval import estrin_plan, dy_norm, dy_to_fr
from fractions import Fraction as Fr
from funcs_spec import Spec


# ============================================================ 係数を h で選ぶ Estrin
def mux_const_bits(h, A, Bc, st):
    """合成時定数ビット列 A（h=1）/ Bc（h=0）の選択 = 配線（等しいビットは定数、違うビットは h か ¬h）。"""
    out = []
    for a, b in zip(A, Bc):
        if a == b: out.append(a)
        elif a: out.append(h)
        else: out.append(NOT(h, st))
    return out

def const_bf_wide(v, E, W, EW):
    """合成時定数 v·2^E → BF。|v| < 2^W なら 正規化 W ビット（const_bf）、それより広い係数（p1 テープ: t ∈ [−1/2, 0) の
       単項式基底では |c_k| > 1 になる）は そのままの幅で持つ（定数は厳密・K1。窓つき加算は 幅の違う mag を受ける）。"""
    if abs(v) < (1 << W): return const_bf(v, E, W, EW)
    a = abs(v); n = a.bit_length()
    return BF([(a >> i) & 1 for i in range(n)], 1 if v < 0 else 0, bus_const(E, EW))

DEXT_OVERRIDE = None        # 計測用（None = 自動）。0 にすると exp と同じ固定 G の窓（根のテープでは監視が破れる）

def estrin_gate_mux(coeffs0, coeffs1, h, u, W, G, st, EW=None):
    """estrin_gate_full と同じ木。葉の係数は h ? coeffs1 : coeffs0（同じ次数）。係数が信号なので 葉の積は 信号積。
       節点の和は 窓を Dext = W + extra + 3 − G 桁下へ伸ばす（extra = 係数の W を超える幅）: 根のテープは 節点
       （葉の対 c_k + c_{k+1}·t など）の根が定義域の内側にあり 相殺が無限に深くなり得るので、固定 G では救えない。
       正しさ（bf_add_win 参照）: 小さい側の LSB が 窓に入る（d ≤ G+Dext）なら 何も落ちず厳密（相殺の深さは無関係）。
       落ちる（d > G+Dext）のは 先頭位置の差が ≥ 2 のときだけ ⇐ G+Dext ≥ li_S − li_U + 1（li = 先頭位置 = 幅−1;
       積 2W+extra ビット vs 正規化 W ビットで 最大 W + extra）。+3 は余裕。"""
    EW = EW or len(u.E)
    n0 = max(k for k, m, E in coeffs0); n1 = max(k for k, m, E in coeffs1)
    assert n0 == n1, "2 本のテープは同じ次数のこと"
    cd0 = {k: dy_norm(m, E) for k, m, E in coeffs0}; cd1 = {k: dy_norm(m, E) for k, m, E in coeffs1}
    extra = max(0, max(abs(v[0]).bit_length() for v in list(cd0.values()) + list(cd1.values())) - W)
    Dext = max(0, W + extra + 3 - G) if DEXT_OVERRIDE is None else DEXT_OVERRIDE
    L, nodes, root = estrin_plan(n0)
    pw = [u]
    for l in range(1, L):
        sq = bf_mul_sig(pw[-1], pw[-1], st, EW)
        pw.append(norm_mag(sq.mag, sq.sign, sq.E, W, st, EW)[0])
    def leaf(k):
        a, b = cd0.get(k, (0, 0)), cd1.get(k, (0, 0))
        if a[0] == 0 and b[0] == 0: return None
        assert a[0] != 0 and b[0] != 0, "片方だけ零の係数は未対応"
        A = const_bf_wide(a[0], a[1], W, EW); Bc = const_bf_wide(b[0], b[1], W, EW)
        n = max(len(A.mag), len(Bc.mag))                                     # 幅を揃える（上位 0 詰め = 配線）
        return BF(mux_const_bits(h, Bc.mag + [0] * (n - len(Bc.mag)), A.mag + [0] * (n - len(A.mag)), st),
                  mux_bit(h, Bc.sign, A.sign, st), mux_const_bits(h, Bc.E, A.E, st), 0)
    val = {}
    for nid, kind, a, b, lvl in nodes:
        if kind == 'leaf':
            val[nid] = leaf(a); continue
        prod = bf_mul_sig(val[b], pw[lvl], st, EW)
        if a is not None and val[a] is not None:
            val[nid] = bf_add_win(val[a], prod, W, G, st, EW, Dext)[0]
        else:
            val[nid] = norm_mag(prod.mag, prod.sign, prod.E, W, st, EW)[0]
    return val[root]


# ============================================================ フラグ伝播（Spec._propagate の log / sqrt / rsqrt）
def propagate_gate_lr(fn, fin, s, zero, fc, cmp1, st):
    """fin=(gi,li,si) 入力, fc=(cg,cl,cs) 計算, cmp1=(gt1, lt1, eq1)（log だけ使う）。x<0 の 111 は呼び側で。"""
    gi, li, si = fin; cg, cl, cs = fc
    none_in = NOT(OR(OR(gi, li, st), si, st), st)
    both = AND(gi, li, st)
    if fn == 'log':
        gt1, lt1, eq1 = cmp1
        sunk_case = OR(OR(OR(si, zero, st), both, st), OR(AND(gt1, li, st), AND(lt1, gi, st), st), st)
        dg = OR(OR(AND(gt1, gi, st), AND(lt1, li, st), st), eq1, st); dl = 0
    elif fn == 'sqrt':
        sunk_case = OR(si, zero, st); dg, dl = gi, li
    else:
        sunk_case = OR(si, zero, st); dg, dl = li, gi
    og = AND(NOT(cl, st), dg, st); ol = AND(NOT(cg, st), dl, st)
    neither = NOT(OR(og, ol, st), st)
    eleven = OR(OR(sunk_case, both, st), neither, st)
    ge = mux_bit(none_in, cg, OR(eleven, og, st), st)
    le = mux_bit(none_in, cl, OR(eleven, ol, st), st)
    sk = mux_bit(none_in, cs, OR(sunk_case, AND(AND(neither, NOT(both, st), st), cs, st), st), st)
    return ge, le, sk


# ============================================================ 本体
class GateLogRoot:
    """f32/f64 構成の log・sqrt・rsqrt ネットリスト（Spec と同じ定数・同じ木）。"""
    FNS = ('log', 'sqrt', 'rsqrt')

    def __init__(self, spec: Spec, G=4, EW=None):
        self.s = spec; self.G = G
        c = spec
        n = max(max(k for k, m, E in c.tapes[t].coeffs) for t in ("log", "sqrt_p0", "sqrt_p1", "rsqrt_p0", "rsqrt_p1"))
        L = estrin_plan(n)[0]                                    # GateExp と同じ式（K2: 現れる指数の最悪値から）
        lo = (1 << (L - 1)) * (-c.Emin + c.Wk) + 2 * c.Wk + 64
        hi = c.Emax + c.Win + c.Pc + c.Q + c.kmax + 8
        self.EW = EW or max(max(lo, hi).bit_length() + 2, GateExp(spec, G).EW)
        self.XB = (max(c.Emax + c.Win, -c.Emin) + 1).bit_length()   # |Ex'| の幅（K2）

    # ---------------------------------------------------------------- 縮小: x → m（正規化 Win ビット）, Ex = E + L, s, zero
    def reduce(self, xd, Eb, st):
        c, EW = self.s, self.EW
        xm, s = canon_mag(xd, st); xm = trim_bits(xm)
        m, dropped = norm_mag(xm, s, Eb, c.Win, st, EW)          # E' = E + L − (Win−1)
        assert bit_v(dropped) == 0                                # 監視: |N| < 2^Win なので落ちない
        Ex = bus_add_fast(m.E, bus_const(c.Win - 1, EW), st)
        return m, Ex, s, m.zero

    def _t_of(self, m, h, zero, st):
        """t = m·2^-h − 1 → 正規化 Wk ビット BF（厳密）。x = 0 なら t = 0。"""
        c, EW = self.s, self.EW
        Win = c.Win
        nz = NOT(zero, st); nh = NOT(h, st)
        digs = [(AND(m.mag[i], nz, st), 0) for i in range(Win - 1)] \
             + [(AND(m.mag[Win - 1], nz, st), AND(nh, nz, st))] + [(0, AND(h, nz, st))]
        tm, ts = canon_mag(digs, st)
        Et = mux_bus(h, bus_const(-Win, EW), bus_const(-(Win - 1), EW), st)
        u, dropped = norm_mag(tm, ts, Et, c.Wk, st, EW)
        assert bit_v(dropped) == 0                                # 監視: t は Win+1 ビット ≤ Wk
        return u

    def _cmp1(self, m, Ex, st):
        """(x > 1, x < 1, x == 1)（x = m·2^Ex, m ∈ [1,2)）。"""
        exneg = Ex[-1]; exzero = NOT(or_tree(Ex, st), st)
        restnz = or_tree(m.mag[:self.s.Win - 1], st)
        gt1 = AND(NOT(exneg, st), OR(NOT(exzero, st), restnz, st), st)
        eq1 = AND(exzero, NOT(restnz, st), st)
        return gt1, exneg, eq1

    # ---------------------------------------------------------------- log: P = trunc_Wk(Ex'·ln2 + u·q(u))
    def log_P(self, m, Ex, zero, st, marks=None):
        c, EW, G = self.s, self.EW, self.G
        Win, Wk, Pc = c.Win, c.Wk, c.Pc
        def mark(name):
            if marks is not None: marks[name] = total(st.fold)
        top8 = m.mag[Win - 8:]
        h = NOT(bus_lt_fast(zext(top8, 9), bus_const(181, 9), st), st)       # top8 ≥ 181
        u = self._t_of(m, h, zero, st)
        Ex1 = bus_add_fast(Ex, zext([h], EW), st)
        mark('reduce:log')
        q = estrin_gate_full(c.tapes["log"].coeffs, u, Wk, G, st, EW)        # テープ 1 本（定数積）
        mark('estrin:log')
        P0 = bf_mul_sig(u, q, st, EW)                                         # u·q（厳密 2Wk ビット）
        exs = Ex1[-1]
        exm = bin_add_fast([XOR(b, exs, st) for b in Ex1], [0] * EW, st, cin=exs)[:self.XB]   # |Ex'|（K2: XB ビット）
        ez = NOT(or_tree(exm, st), st)
        Cm = mul_mag_const(exm, c.LN2, st)                                    # |Ex'|·LN2（二進 Wallace）
        C = BF(Cm, AND(exs, NOT(ez, st), st), mux_bus(ez, zero_E(EW), bus_const(-Pc, EW), st), ez)
        P, _, _ = bf_add_win(C, P0, Wk, G, st, EW)
        mark('tail:log')
        return P

    # ---------------------------------------------------------------- sqrt / rsqrt: P = Estrin(p_h, u)·2^(±Ex'/2)
    def root_P(self, fn, m, Ex, zero, st, marks=None):
        c, EW, G = self.s, self.EW, self.G
        Wk = c.Wk
        def mark(name):
            if marks is not None: marks[name] = total(st.fold)
        h = Ex[0]
        u = self._t_of(m, h, zero, st)
        Ex1 = bus_add_fast(Ex, zext([h], EW), st)                             # 偶数
        half = Ex1[1:] + [Ex1[-1]]                                            # Ex' >> 1（算術・配線）
        mark(f'reduce:{fn}')
        v = estrin_gate_mux(c.tapes[f"{fn}_p0"].coeffs, c.tapes[f"{fn}_p1"].coeffs, h, u, Wk, G, st, EW)
        mark(f'estrin:{fn}')
        E = bus_add_fast(v.E, half, st) if fn == 'sqrt' else bus_sub_fast(v.E, half, st)
        E = mux_bus(v.zero, zero_E(EW), E, st)
        mark(f'tail:{fn}')
        return BF(v.mag, v.sign, E, v.zero)

    # ---------------------------------------------------------------- 全体（縮小 1 回 → 関数ごとの前半 → 出口 3 個）
    def run_all(self, xd, Eb, st, fin=(0, 0, 0), fns=FNS, modes=('lo', 'hi', 'near'), marks=None):
        """戻り {(fn, mode): ((Wout 桁, E バス), (ge, le, sunk))}。marks に dict を渡すと ブロック別 fold ゲート数を書く。"""
        c, EW = self.s, self.EW
        Wout = c.Wout
        def mark(name):
            if marks is not None: marks[name] = total(st.fold)
        mark('start')
        m, Ex, s, zero = self.reduce(xd, Eb, st)
        cmp1 = self._cmp1(m, Ex, st)
        neg = AND(s, NOT(zero, st), st); nz_ = NOT(zero, st)
        mark('reduce')
        out = {}
        for fn in fns:
            P = self.log_P(m, Ex, zero, st, marks) if fn == 'log' else self.root_P(fn, m, Ex, zero, st, marks)
            for mode in modes:
                (digs, Eo), fc = finish_gate(P, c.e[fn], mode, Wout, c.Emin, c.Emax, st, EW)
                if fn == 'log':                                               # x = 0: −MAX, (1,0,0)
                    maxd = to_sd(-((1 << Wout) - 1), Wout)
                    digs = GE.mux_digits(zero, maxd, digs, st); Eo = mux_bus(zero, bus_const(c.Emax, EW), Eo, st)
                    fc = (OR(fc[0], zero, st), AND(fc[1], nz_, st), AND(fc[2], nz_, st))
                else:                                                         # x = 0: (0, Emin), 000
                    digs = [(AND(p, nz_, st), AND(n, nz_, st)) for (p, n) in digs]
                    Eo = mux_bus(zero, bus_const(c.Emin, EW), Eo, st)
                    fc = tuple(AND(f, nz_, st) for f in fc)
                flags = propagate_gate_lr(fn, fin, s, zero, fc, cmp1, st)
                flags = tuple(OR(f, neg, st) for f in flags)                 # x < 0 → 111
                out[(fn, mode)] = ((digs, Eo), flags)
                mark(f'finish:{fn}:{mode}')
        return out


# ============================================================ 差分テスト（Spec と bit 一致）
def eval_point(g, N, E, fin=(0, 0, 0), digits=None, want_depth=False):
    c = g.s
    st = FoldCounter(); marks = {}
    if digits is None:
        xd = sig_digits(N, c.Win)
    else:
        xd = [(B(p), B(n)) for (p, n) in digits]
        assert val_digits(xd) == N
    Eb = sig_bus(E, g.EW)
    finb = tuple(B(f) for f in fin)
    res = g.run_all(xd, Eb, st, fin=finb, marks=marks)
    out = {}; depth = {}
    for key, ((digs, Eo), flags) in res.items():
        out[key] = ((val_digits(digs), val_bus(Eo)), tuple(bit_v(f) for f in flags))
        if want_depth:
            depth[key] = max(depth_of(digs), depth_of(Eo), depth_of(list(flags)))
    return out, st, marks, depth

def spec_point(spec, N, E, fin=(0, 0, 0)):
    out = {}
    for fn in GateLogRoot.FNS:
        for mode in ('lo', 'hi', 'near'):
            out[(fn, mode)] = getattr(spec, fn)((N, E), fin, mode)
    return out

def compare_point(g, N, E, fin=(0, 0, 0), digits=None):
    got, st, marks, _ = eval_point(g, N, E, fin, digits)
    want = spec_point(g.s, N, E, fin)
    bad = {k: (got[k], want[k]) for k in want if got[k] != want[k]}
    return bad, st, marks


_GL = {}

def _worker(args):
    cfg, N, E, fin, digits = args
    g = _GL[cfg]
    try:
        bad, st, marks = compare_point(g, N, E, fin, digits)
        return (N, E, fin, digits is not None, bad, None)
    except AssertionError as ex:
        return (N, E, fin, digits is not None, None, str(ex))

def _adversarial_inputs(spec, rng, n_rand=200):
    """2 のべき（t=0）、1 の近く、top8=181 の境界（log の h）、Ex の偶奇（根の h）、完全平方、MIN/MAX、対数一様ランダム。
       同じ値でも L の違う表現（N=1 と N=2^(Win−1)）を入れて 正規化器を踏む。"""
    c = spec; Win = c.Win
    pts = []
    def add(N, E):
        if N != 0 and abs(N) < (1 << Win): pts.append((int(N), int(E)))
    top = 1 << (Win - 1)
    Es = [0, -(Win - 1), -Win, 1, -1, 2, -2, 3, -3, c.Emin, c.Emin + 1, c.Emax, c.Emax - 1, c.Emax - (Win - 1),
          -(Win - 1) + 1, -(Win - 1) - 1]
    for E in Es + list(range(-6, 7)):
        for sgn in (1, -1):
            add(sgn, E); add(sgn * top, E); add(sgn * top, E - (Win - 1))
    for d in (-3, -2, -1, 1, 2, 3):
        for sgn in (1, -1):
            add(sgn * (top + d), -(Win - 1)); add(sgn * (top + d), -(Win - 1) + 1); add(sgn * (top + d), -Win)
    for sgn in (1, -1):
        add(sgn * ((1 << Win) - 1), -Win); add(sgn * ((1 << Win) - 3), -Win); add(sgn * ((1 << Win) - 1), -(Win - 1))
    for t8 in (180, 181, 182, 128, 129, 255, 254):
        for low in (0, (1 << (Win - 8)) - 1, int(rng.integers(0, 1 << (Win - 8))), 1):
            N = (t8 << (Win - 8)) | low
            for E in (0, -(Win - 1), -Win, 5, -7, c.Emin, c.Emax - (Win - 1)):
                for sgn in (1, -1): add(sgn * N, E)
    for _ in range(20):
        N = int(rng.integers(1, 1 << Win)); E = int(rng.integers(c.Emin, c.Emax - Win))
        for d in (0, 1, -1): add(N, E + d); add(-N, E + d)
    for k in (3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 4095, 4097):
        for E in (0, 2, -2, 4, -4, -10, 10, c.Emin, c.Emin + 1):
            add(k * k, E); add(-k * k, E)
    add(1, c.Emin); add(-1, c.Emin); add(3, c.Emin); add(top | 1, c.Emin); add(top, c.Emin + 1)
    add((1 << Win) - 1, c.Emax); add(-((1 << Win) - 1), c.Emax); add(top, c.Emax); add(top | 1, c.Emax)
    for _ in range(n_rand):
        M = int(rng.integers(1, 1 << Win)); s = 1 if rng.random() < 0.5 else -1
        ld = int(rng.integers(c.Emin, c.Emax + Win))
        E = min(max(ld - (M.bit_length() - 1), c.Emin), c.Emax)
        add(s * M, E)
    return pts

def self_test(cfg="f32", n_rand=200, seed=20260902, procs=None):
    import time, numpy as np, multiprocessing
    spec = Spec(cfg); g = GateLogRoot(spec)
    _GL[cfg] = g
    rng = np.random.default_rng(seed)
    print(f"[{cfg}] EW={g.EW} XB={g.XB} G={g.G} Win={spec.Win} Wk={spec.Wk} Pc={spec.Pc} "
          f"deg log={max(k for k,_,_ in spec.tapes['log'].coeffs)} sqrt={max(k for k,_,_ in spec.tapes['sqrt_p0'].coeffs)} "
          f"rsqrt={max(k for k,_,_ in spec.tapes['rsqrt_p0'].coeffs)}")
    t0 = time.time()
    N0, E0 = (1 << (spec.Win - 1)) | 12345, -(spec.Win - 1)     # x ≈ 1.0015
    got, st, marks, depth = eval_point(g, N0, E0, want_depth=True)
    want = spec_point(spec, N0, E0)
    print(f"  1 点 {time.time()-t0:.1f}s  raw={total(st.raw)}  fold={total(st.fold)}")
    prev = 0; per_fn = {}
    for k in marks:
        d = marks[k] - prev; prev = marks[k]
        print(f"    {k:18s} {d:8d}")
        head = k.split(':')
        fn = head[1] if len(head) > 1 else 'common'
        per_fn[fn] = per_fn.get(fn, 0) + d
    for fn in GateLogRoot.FNS:
        dm = max(depth[(fn, mode)] for mode in ('lo', 'hi', 'near'))
        print(f"    {fn:6s}: fold {per_fn.get('common',0) + per_fn.get(fn,0):8d} (共通 {per_fn.get('common',0)} + 固有 {per_fn.get(fn,0)})  depth {dm}")
    for k in want:
        print(f"    {k}: gate={got[k]} spec={want[k]} {'OK' if got[k]==want[k] else 'MISMATCH'}")
    pts = _adversarial_inputs(spec, rng, n_rand)
    jobs = [(cfg, N, E, (0, 0, 0), None) for (N, E) in pts]
    flagpts = pts[:6] + [(1, 0), ((1 << (spec.Win - 1)) + 1, -(spec.Win - 1)), ((1 << spec.Win) - 1, -spec.Win),
                         (-5, 0), (5, 0), (5, -3), (-3, -2)]
    for (N, E) in flagpts:
        for fin in ((1,0,0),(0,1,0),(1,1,0),(0,0,1),(1,0,1),(0,1,1),(1,1,1)):
            jobs.append((cfg, N, E, fin, None))
    for (N, E) in pts[::5]:
        d = _redundant_digits(rng, N, spec.Win)
        if d is not None: jobs.append((cfg, N, E, (0, 0, 0), d))
    for E in (0, -5, 3, spec.Emin, spec.Emax):
        jobs.append((cfg, 0, E, (0, 0, 0), None)); jobs.append((cfg, 0, E, (1, 0, 0), None))
    jobs.append((cfg, 0, 3, (0, 0, 1), [(1, 1)] * spec.Win)); jobs.append((cfg, 0, -2, (0, 1, 0), [(1, 1)] * spec.Win))
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


def node_polys(coeffs):
    """Estrin の各 mac 節点（a ≠ None）を t の多項式で: [(nid, lvl, pa, pb_shifted, p)]（係数は昇冪の float 配列）。
       節点の値 = pa(t) + pb_shifted(t) で、相殺の深さ = log2(max(|pa|,|pb|)/|p|)。"""
    import numpy as np
    cd = {k: float(dy_to_fr(dy_norm(m, E))) for k, m, E in coeffs}
    nmax = max(cd); L, nodes, root = estrin_plan(nmax)
    poly = {}; out = []
    for nid, kind, a, b, lvl in nodes:
        if kind == 'leaf':
            poly[nid] = np.array([cd.get(a, 0.0)]); continue
        pb = np.concatenate([np.zeros(1 << lvl), poly[b]])               # val[b]·t^(2^lvl)
        pa = poly[a] if a is not None else np.zeros(1)
        p = np.zeros(max(len(pa), len(pb))); p[:len(pa)] += pa; p[:len(pb)] += pb
        poly[nid] = p
        if a is not None: out.append((nid, lvl, pa, pb, p))
    return out

def node_roots(coeffs, lo, hi):
    """定義域 [lo,hi] 内にある 節点多項式の実根（相殺が無限に深くなる点）。"""
    import numpy as np
    out = []
    for nid, lvl, pa, pb, p in node_polys(coeffs):
        if len(p) > 1:
            for r in np.roots(p[::-1]):
                if abs(r.imag) < 1e-12 and lo <= r.real <= hi: out.append(float(r.real))
    return out

def cancel_report(cfg="f32", npts=100001):
    """全テープの Estrin 節点の相殺の深さ（ビット）: 定義域内の根の個数 と 格子上の最大。根が無ければ 深さは有界
       （固定 G で足りるかの根拠）、根があれば 非有界（Dext の窓が要る）。"""
    import numpy as np
    spec = Spec(cfg); rows = []
    for name, tp in spec.tapes.items():
        lo, hi = (float(Fr(v)) for v in tp.d["interval"])
        ts = np.linspace(lo, hi, npts)
        worst = 0.0; nroot = len(node_roots(tp.coeffs, lo, hi))
        for nid, lvl, pa, pb, p in node_polys(tp.coeffs):
            A = np.polyval(pa[::-1], ts); Bv = np.polyval(pb[::-1], ts); S = A + Bv
            big = np.maximum(np.abs(A), np.abs(Bv)); ok = S != 0
            worst = max(worst, float(np.max(np.log2(big[ok] / np.abs(S[ok])))))
        rows.append((name, nroot, worst))
        print(f"  {cfg} {name:10s} 定義域内の節点根 {nroot}  格子上の最大相殺 {worst:5.1f} bit")
    return rows

def test_estrin_mux(cfg="f32", n_rand=24, seed=3):
    """R4: 係数 mux の Estrin 単体。u（Wk ビット, [−1/2, 1)）と h を直接注入し、estrin_spec(テープ h, u) と bit 一致を見る。
       広い係数（|c_k| > 1 の p1 テープ）と 木の全節点の窓つき加算を、入力経路を通さずに踏む。"""
    import numpy as np
    from tape_eval import estrin_spec
    spec = Spec(cfg); g = GateLogRoot(spec); Wk, EW = spec.Wk, g.EW
    rng = np.random.default_rng(seed)
    bad = n = 0
    for fn in ("sqrt", "rsqrt"):
        c0 = spec.tapes[f"{fn}_p0"].coeffs; c1 = spec.tapes[f"{fn}_p1"].coeffs
        us = [(0, 0), (1, -1), (-1, -1), (1, -Wk), (-1, -Wk), ((1 << Wk) - 1, -Wk), (-((1 << Wk) - 1), -(Wk + 1))]
        for _ in range(n_rand):
            m = int.from_bytes(rng.bytes(16), 'little') % ((1 << Wk) - 1) + 1   # Wk=64 は int64 に収まらない
            lead = int(rng.integers(-Wk, 0))
            us.append(((-1 if rng.random() < 0.5 else 1) * m, lead - (m.bit_length() - 1)))
        nroots = 0
        for coeffs, lo, hi in ((c0, 0.0, 1.0), (c1, -0.5, 0.0)):
            for r in node_roots(coeffs, lo, hi):                          # 節点の根 ± 1 ulp: 相殺が最も深い入力
                q = int(round(r * 2 ** Wk))
                for d in (-1, 0, 1):
                    if q + d != 0: us.append((q + d, -Wk)); nroots += 1
        print(f"  {fn}: 節点の根（定義域内） {nroots // 3} 個")
        for h in (0, 1):
            for (N, E) in us:
                uf = dy_to_fr((N, E))
                if h == 0 and not (0 <= uf < 1): continue                    # p0: t ∈ [0, 1)
                if h == 1 and not (-Fr(1, 2) <= uf <= 0): continue           # p1: t ∈ [−1/2, 0]
                st = FoldCounter()
                if N == 0:
                    u = BF([B(0)] * Wk, B(0), sig_bus(-(1 << (EW - 2)), EW), B(1))
                else:
                    a = abs(N); L = a.bit_length() - 1; sh = (Wk - 1) - L
                    u = BF([B((a << sh >> i) & 1) for i in range(Wk)], B(1 if N < 0 else 0), sig_bus(E - sh, EW), B(0))
                v = estrin_gate_mux(c0, c1, B(h), u, Wk, g.G, st, EW)
                got = dy_norm(val_bits(v.mag) * (-1 if bit_v(v.sign) else 1), val_bus(v.E))
                want = estrin_spec(c1 if h else c0, dy_norm(N, E), Wk)
                n += 1
                if got != want:
                    bad += 1
                    if bad <= 10: print(f"  MISMATCH {fn} h={h} u=({N},{E}): gate={got} spec={want}")
    print(f"  estrin_mux[{cfg}] {n} 例  不一致 {bad}")
    return bad == 0


if __name__ == "__main__":
    cfg = sys.argv[1] if len(sys.argv) > 1 else "f32"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    ok = test_estrin_mux(cfg) and self_test(cfg, n)
    print("PASS" if ok else "FAIL")
