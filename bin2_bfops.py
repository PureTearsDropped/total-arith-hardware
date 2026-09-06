#!/usr/bin/env python3
"""Binary counterpart of gate_bfp.py: block floating point on carry-save mantissas.
number = (rows r0, r1 of Wc bits, E host int); value = ((r0 + r1) mod 2^Wc read as two's complement) · 2^E.

Contract that differs from the signed-digit layer: the block width Wc is fixed up front and every row lives at
Wc (all arithmetic is exact mod 2^Wc; the unit guarantees |true value| < 2^(Wc−1) by construction, tracked
statically as Wv = bound on the value width).  A carry-save number cannot be widened later (its rows are not
sign-extendable: the top bits of a row are not sign bits), so the width has to be declared where the signed-digit
layer just grows lists.  Exponent alignment (shift up = ×2^k) is wiring at fixed Wc: prepend k zeros, drop k top bits.

Multiplication of two mantissas, both possibly redundant:
  mode 'resolve'   : resolve the redundant operand(s) with one Kogge–Stone each, then a Baugh–Wooley tree of
                     Wv_x × Wv_y partial products with the exact integer constant folded into the accumulator
                     (no sign extension of rows: the constant is exact, not reduced mod 2^W).
  mode 'redundant' : carry-free — (r0 + r1)(s0 + s1) mod 2^Wc = the four cross products of the *unsigned* row
                     patterns mod 2^Wc, all into one Dadda tree; no carry chain anywhere, but Σ_{i+j<Wc} ≈ Wc²/2
                     AND gates per cross product (4 sets) instead of Wv_x·Wv_y (1 set).
Constant coefficients (bf_lincomb): CSD digits of c on the rows themselves (shift = wiring, negative digit =
complemented rows + constant); no resolve, no partial products.
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from gate_bilinear import AND, OR, NOT, XOR, new_counter
from bin2_gates import to_bits, from_bits, cs_val, dadda_rows, resolve, neg_cs, add_cs

def _z(b): return isinstance(b, int) and b == 0
def _one(b): return isinstance(b, int) and b == 1
def _and(a, b, st):
    if _z(a) or _z(b): return 0
    if _one(a): return b
    if _one(b): return a
    return AND(a, b, st)
def _not(a, st): return (1 - a) if isinstance(a, int) else NOT(a, st)

class BFc:
    """carry-save block-floating number: rows (r0, r1) at width Wc, exponent E (host int), Wv = static bound on
    the number of bits the value needs (incl. sign); nr = True when r1 is all-zero (non-redundant)."""
    __slots__ = ('rows', 'E', 'Wv')
    def __init__(self, rows, E=0, Wv=None):
        self.rows = (list(rows[0]), list(rows[1])); self.E = E
        self.Wv = Wv if Wv is not None else len(rows[0])
    @property
    def Wc(self): return len(self.rows[0])
    @property
    def nr(self): return all(_z(b) for b in self.rows[1])

def to_bfc(value, Wc, E=0, Wv=None):
    """integer → non-redundant BFc at width Wc (row1 = 0)"""
    return BFc((to_bits(value, Wc), [0] * Wc), E, Wv or Wc)

def from_bfc(x): return Fr(cs_val(x.rows, x.Wc)) * Fr(2) ** x.E

def bfc_from_bits(bits, Wc, E=0):
    """non-redundant two's complement input of W bits → BFc at width Wc (sign extension = wiring)"""
    W = len(bits); return BFc((list(bits) + [bits[-1]] * (Wc - W), [0] * Wc), E, W)

def bfc_from_sm(sign, mag, Wc, E=0, st=None):
    """normaliser output (sign, W-bit magnitude) → BFc: −m = ~m + 1, i.e. r0 = m ⊕ s (sign-extended with s), r1 = s at bit 0"""
    W = len(mag)
    r0 = [XOR(m, sign, st) if not isinstance(m, int) else (_not(sign, st) if m else sign) for m in mag] + [sign] * (Wc - W)
    return BFc((r0, [sign] + [0] * (Wc - 1)), E, W + 1)

# ------------------------------------------------------------------ partial products into a fixed-width column set
def pp_signed_into(cols, X, Y, Wc, st, neg=False):
    """Baugh–Wooley partial products of two's complement X×Y (exact integer identity, constant returned, not
    reduced mod 2^len) into cols (weights ≥ Wc dropped).  neg=True negates the product: complement every bit
    (folded into the inversion parity) and return the negated constant minus the sum of all pp weights."""
    Wx, Wy = len(X), len(Y); W = Wx + Wy
    for p, x in enumerate(X):
        for q, y in enumerate(Y):
            if p + q >= Wc: continue
            inv = ((p == Wx - 1) != (q == Wy - 1)) != neg
            t = _and(x, y, st)
            if inv: t = _not(t, st)
            if not _z(t): cols[p + q].append(t)
    c = -((1 << (W - 1)) - (1 << (Wx - 1)) - (1 << (Wy - 1)))
    if neg: c = -c - sum((1 << (p + q)) for p in range(Wx) for q in range(Wy) if p + q < Wc)
    return c

def pp_unsigned_into(cols, R, S, Wc, st):
    """plain partial products of two unsigned bit patterns mod 2^Wc"""
    for p, x in enumerate(R):
        if _z(x): continue
        for q, y in enumerate(S):
            if p + q >= Wc or _z(y): continue
            cols[p + q].append(_and(x, y, st))

def _const_into(cols, c, Wc):
    c &= (1 << Wc) - 1
    for w in range(Wc):
        if (c >> w) & 1: cols[w].append(1)

def _finish(cols, Wc, st):
    cols = [[b for b in c if not _z(b)] for c in cols]
    return dadda_rows(cols, Wc, st)

def _operand_bits(x, st):
    """non-redundant Wv-bit two's complement operand of x (resolve with one Kogge–Stone if redundant)"""
    if x.nr: return x.rows[0][:x.Wv]
    return resolve(x.rows, st)[:x.Wv]

# ------------------------------------------------------------------ operations (mirror of gate_bfp)
def shift_up(x, k):
    """×2^k at fixed width (wiring): value·2^k mod 2^Wc"""
    if k == 0: return x
    Wc = x.Wc
    return BFc(([0] * k + x.rows[0][:Wc - k], [0] * k + x.rows[1][:Wc - k]), x.E - k, x.Wv + k)

def bfc_mul(x, y, st, mode='resolve'):
    assert x.Wc == y.Wc, "operands must share the block width"
    Wc = x.Wc; cols = [[] for _ in range(Wc)]
    assert x.Wv + y.Wv <= Wc, f"product needs {x.Wv + y.Wv} bits > Wc={Wc}"
    if mode == 'resolve':
        c = pp_signed_into(cols, _operand_bits(x, st), _operand_bits(y, st), Wc, st); _const_into(cols, c, Wc)
    else:
        for R in x.rows:
            for S in y.rows: pp_unsigned_into(cols, R, S, Wc, st)
    return BFc(_finish(cols, Wc, st), x.E + y.E, x.Wv + y.Wv)

def bfc_neg(x, st):
    return BFc(neg_cs(x.rows, x.Wc, st), x.E, x.Wv + 1)

def bfc_add(x, y, st): return bfc_sum([x, y], st)

def bfc_sub(x, y, st): return bfc_sum([x, bfc_neg(y, st)], st)

def bfc_sum(terms, st):
    """Σ terms: align to Elo (shift up = wiring), one Dadda tree over all rows"""
    Wc = terms[0].Wc; assert all(t.Wc == Wc for t in terms)
    Elo = min(t.E for t in terms); rows = []; Wv = 0
    for t in terms:
        s = shift_up(t, t.E - Elo); rows += [s.rows[0], s.rows[1]]; Wv = max(Wv, s.Wv)
    Wv += max(1, (len(terms) - 1).bit_length())
    assert Wv <= Wc, f"sum needs {Wv} bits > Wc={Wc}"
    return BFc(add_cs(rows, Wc, st), Elo, Wv)

def csd(c):
    """canonical signed-digit form of an integer: list of (k, ±1) with c = Σ ±2^k, no two adjacent"""
    out = []; k = 0
    while c != 0:
        if c & 1:
            d = 2 - (c & 3)          # 1 → +1, 3 → −1
            out.append((k, d)); c -= d
        c >>= 1; k += 1
    return out

def bfc_scale(x, c, st):
    """constant c × carry-save number, carry-free: CSD digits of c applied to the rows themselves
    (2^k = wiring; −2^k·(r0+r1) = 2^k·(~r0 + ~r1 + 2) mod 2^Wc); one Dadda tree.  c = ±1 costs nothing / 2 inverter rows."""
    if c == 0: return None
    Wc = x.Wc; cols = [[] for _ in range(Wc)]; const = 0
    for k, d in csd(c):
        for r in x.rows:
            for w, b in enumerate(r):
                if w + k >= Wc: continue
                t = b if d > 0 else _not(b, st)
                if not _z(t): cols[w + k].append(t)
        if d < 0: const += 2 << k
    _const_into(cols, const, Wc)
    return BFc(_finish(cols, Wc, st), x.E, x.Wv + max(1, abs(c).bit_length()))

def bfc_lincomb(coeffs, xs, st):
    """Σ coeffs·xs on BFc: ±1 as rows (negation = complemented rows + constant 2), other constants by CSD on the
    rows, one alignment and one Dadda tree for the whole combination."""
    Wc = xs[0].Wc; cols = [[] for _ in range(Wc)]; const = 0
    terms = [(c, x) for c, x in zip(coeffs, xs) if c != 0]
    if not terms: return BFc(([0] * Wc, [0] * Wc), 0, 1)
    Elo = min(x.E for _, x in terms); Wv = 0
    for c, x in terms:
        s = shift_up(x, x.E - Elo); Wv = max(Wv, s.Wv + max(1, abs(c).bit_length()))
        for k, d in csd(c):
            for r in s.rows:
                for w, b in enumerate(r):
                    if w + k >= Wc: continue
                    t = b if d > 0 else _not(b, st)          # a constant 0 bit complements to 1: keep it
                    if not _z(t): cols[w + k].append(t)
            if d < 0: const += 2 << k
    _const_into(cols, const, Wc)
    Wv += max(1, (len(terms) - 1).bit_length())
    assert Wv <= Wc, f"lincomb needs {Wv} bits > Wc={Wc}"
    return BFc(_finish(cols, Wc, st), Elo, Wv)

def bfc_bilinear_unit(U, V, W, a, b, st, mode='resolve'):
    """c = W·((U·a) ⊙ (V·b)) on BFc (unnormalised carry-save outputs), mirror of bf_bilinear_unit"""
    R = len(U)
    left = [bfc_lincomb(U[r], a, st) for r in range(R)]
    right = [bfc_lincomb(V[r], b, st) for r in range(R)]
    prod = [bfc_mul(left[r], right[r], st, mode) for r in range(R)]
    return [bfc_lincomb(W[k], prod, st) for k in range(len(W))]

def _wv_lincomb(coeffs, Wv_in):
    terms = [c for c in coeffs if c != 0]
    if not terms: return 1
    return max(Wv_in + max(1, abs(c).bit_length()) for c in terms) + max(1, (len(terms) - 1).bit_length())

def unit_width(Win, U, V, W):
    """block width Wc that makes the unit exact (same static bookkeeping as the operations)"""
    wl = max(_wv_lincomb(r, Win) for r in U); wr = max(_wv_lincomb(r, Win) for r in V)
    return max(_wv_lincomb(r, wl + wr) for r in W)

# ------------------------------------------------------------------ self-test: golden equivalence with gate_bfp (signed digits)
def self_test(seed=20260906):
    from gate_bfp import to_bf, from_bf, bf_mul, bf_add, bf_sub, bf_lincomb, bf_bilinear_unit, block_normalize
    from gate_bilinear import to_sd, from_sd
    from gate_exponent import bus_val
    from bin2_bfp import block_normalize_bin, sm_val
    rng = random.Random(seed); ok = True
    def report(name, bad, n):
        nonlocal ok; ok = ok and bad == 0
        print(f"  {name}: {'all equal' if bad == 0 else f'{bad}/{n} DIFFER'}")
    print("bin2_bfops vs gate_bfp (same values, binary carry-save vs signed digits)")
    # B1: framing and free shift
    bad = 0
    for v in (-2033, -1, 0, 1, 2033):
        x = to_bfc(v, 16, 5, 12); s = shift_up(x, 3)
        if from_bfc(x) != Fr(v) * 32 or from_bfc(s) != from_bfc(x): bad += 1
    for _ in range(300):                       # normaliser output → BFc round trip
        v = rng.randint(-63, 63); s = 1 if v < 0 else 0; m = to_bits(abs(v), 7)
        x = bfc_from_sm(s, m, 12, 0, new_counter())
        if from_bfc(x) != v: bad += 1
    report("B1 framing / shift_up / from (sign, magnitude)", bad, 305)
    # B2: products (both modes) == SD product
    bad = 0; n = 2000
    for _ in range(n):
        a, b = rng.randint(-400, 400), rng.randint(-400, 400); Ea, Eb = rng.randint(-4, 4), rng.randint(-4, 4)
        g = bf_mul(to_bf(a, 11, Ea), to_bf(b, 11, Eb), new_counter())
        for mode in ('resolve', 'redundant'):
            z = bfc_mul(to_bfc(a, 24, Ea, 11), to_bfc(b, 24, Eb, 11), new_counter(), mode)
            if from_bfc(z) != from_bf(g) or z.E != g.E: bad += 1
    report("B2 product (resolve and redundant modes)", bad, 2 * n)
    # B3: sums / differences with unequal exponents, including redundant operands
    bad = 0; n = 2000
    for _ in range(n):
        a, b, c = (rng.randint(-500, 500) for _ in range(3)); Ea, Eb, Ec = (rng.randint(-3, 7) for _ in range(3))
        X, Y, Z = to_bf(a, 12, Ea), to_bf(b, 12, Eb), to_bf(c, 12, Ec)
        x, y, z = to_bfc(a, 48, Ea, 12), to_bfc(b, 48, Eb, 12), to_bfc(c, 48, Ec, 12)
        st = new_counter()
        if from_bfc(bfc_add(x, y, st)) != from_bf(bf_add(X, Y, st)): bad += 1
        if from_bfc(bfc_sub(x, y, st)) != from_bf(bf_sub(X, Y, st)): bad += 1
        s = bfc_add(x, y, st)                                      # redundant + non-redundant, then × redundant
        if from_bfc(bfc_sub(s, z, st)) != from_bf(bf_sub(bf_add(X, Y, st), Z, st)): bad += 1
        if from_bfc(bfc_mul(s, s, st)) != from_bf(bf_add(X, Y, st)) ** 2: bad += 1
        if from_bfc(bfc_mul(s, s, st, 'redundant')) != from_bf(bf_add(X, Y, st)) ** 2: bad += 1
    report("B3 sum / difference / redundant chains", bad, 5 * n)
    # B3': linear combinations with constant coefficients (CSD on rows) == bf_lincomb
    bad = 0; n = 1000
    for _ in range(n):
        vs = [rng.randint(-200, 200) for _ in range(4)]; Es = [rng.randint(0, 3) for _ in range(4)]
        cs = [rng.choice([0, 1, -1, 3, -5, 7, 11, -13]) for _ in range(4)]
        g = bf_lincomb(cs, [to_bf(v, 10, E) for v, E in zip(vs, Es)], new_counter())
        z = bfc_lincomb(cs, [to_bfc(v, 24, E, 10) for v, E in zip(vs, Es)], new_counter())
        if from_bfc(z) != from_bf(g): bad += 1
    report("B3' lincomb with constants (CSD on rows, no partial products)", bad, n)
    # B4/B5: Strassen and groups end to end, then the two normalisers: values, exponent, flags must agree
    U_STR = [[1,0,0,1],[0,0,1,1],[1,0,0,0],[0,0,0,1],[1,1,0,0],[-1,0,1,0],[0,1,0,-1]]
    V_STR = [[1,0,0,1],[1,0,0,0],[0,1,0,-1],[-1,0,1,0],[0,0,0,1],[1,1,0,0],[0,0,1,1]]
    W_STR = [[1,0,0,1,-1,0,1],[0,0,1,0,1,0,0],[0,1,0,1,0,0,0],[1,-1,1,0,0,1,0]]
    from nd_algebra import cd_omega, ref_mult_M
    def group_uvw(OM, M):
        U = [[1 if c == i else 0 for c in range(M)] for i in range(M) for j in range(M)]
        V = [[1 if c == j else 0 for c in range(M)] for i in range(M) for j in range(M)]
        Wm = [[int(OM[i][j]) if (i ^ j) == k else 0 for i in range(M) for j in range(M)] for k in range(M)]
        return U, V, Wm
    from gate_exponent import bus_const
    EW = 12
    for name, (U, V, Wm), Win, Wn, n in (("Strassen 2×2", (U_STR, V_STR, W_STR), 9, 6, 400),
                                          ("complex", group_uvw(cd_omega(2), 2), 8, 5, 300),
                                          ("quaternion", group_uvw(cd_omega(4), 4), 8, 6, 300)):
        Wc = unit_width(Win, U, V, Wm); bad = 0; lies = 0; rounded = 0
        for _ in range(n):
            lim = 1 << (Win - 1)
            a = [rng.randint(-lim, lim - 1) for _ in range(len(U[0]))]; b = [rng.randint(-lim, lim - 1) for _ in range(len(V[0]))]
            E0 = rng.randint(0, 3); Emax = rng.choice([40, E0 + 4])
            st_sd = new_counter(); st_b = new_counter()
            outs_sd = bf_bilinear_unit(U, V, Wm, [to_bf(v, Win, E0) for v in a], [to_bf(v, Win, E0) for v in b], st_sd)
            om_sd, E_sd, fl_sd = block_normalize([o.mant for o in outs_sd], outs_sd[0].E, Wn, Emax, st_sd)
            outs_b = bfc_bilinear_unit(U, V, Wm, [to_bfc(v, Wc, E0, Win) for v in a], [to_bfc(v, Wc, E0, Win) for v in b], st_b)
            assert all(o.E == outs_sd[0].E for o in outs_b)
            om_b, E_b, fl_b = block_normalize_bin([o.rows for o in outs_b], bus_const(outs_b[0].E, EW), Wn, Emax, st_b)
            if any(from_bfc(ob) != from_bf(os) for ob, os in zip(outs_b, outs_sd)): bad += 1; continue
            if bus_val(E_b) % (1 << EW) != E_sd % (1 << EW): bad += 1; continue
            for ob, os, fb, fs in zip(om_b, om_sd, fl_b, fl_sd):
                if sm_val(ob) != from_sd(os) or (int(fb[0]), int(fb[1])) != (int(fs[0]), int(fs[1])): bad += 1; break
                if int(fb[0]) or int(fb[1]): rounded += 1
        report(f"B5 {name} unit + normaliser (values, E, ge/le; {rounded} rounded components)", bad, n)
    # gate counts (symbolic run).  gate_bfp itself sits on the unoptimised gate_bilinear layer, so the comparable
    # signed-digit number is the fast fused group unit (mul_fused.group_component, the sed_comp generator) per component.
    from gate_fast import B, depth_of
    from mul_fused import group_component
    print("  gates per unit (all components, unnormalised carry-save / signed-digit outputs):")
    for name, (U, V, Wm), Win, M in (("Strassen 2×2", (U_STR, V_STR, W_STR), 9, None), ("quaternion", group_uvw(cd_omega(4), 4), 8, 4)):
        Wc = unit_width(Win, U, V, Wm); res = {}
        for mode in ('resolve', 'redundant'):
            st = new_counter()
            ba = [bfc_from_bits([B(0) for _ in range(Win)], Wc) for _ in range(len(U[0]))]
            bb = [bfc_from_bits([B(0) for _ in range(Win)], Wc) for _ in range(len(V[0]))]
            outs = bfc_bilinear_unit(U, V, Wm, ba, bb, st, mode)
            res[mode] = (sum(st.values()), max(depth_of(o.rows[0] + o.rows[1]) for o in outs))
        sd = ""
        if M:
            OM = cd_omega(M); OMl = [[int(OM[i, j]) for j in range(M)] for i in range(M)]; st = new_counter(); dmax = 0
            for k in range(M):
                a = [[(B(0), B(0)) for _ in range(Win)] for _ in range(M)]; b = [[(B(0), B(0)) for _ in range(Win)] for _ in range(M)]
                Z = group_component(a, b, OMl, M, k, st); dmax = max(dmax, depth_of([p for d in Z for p in d]))
            sd = f"SD fast fused (group_component ×{M}) {sum(st.values())} gates / depth {dmax} | "
        print(f"    {name} (Win={Win}, Wc={Wc}): {sd}binary resolve {res['resolve'][0]} / {res['resolve'][1]} | binary redundant {res['redundant'][0]} / {res['redundant'][1]}")
    return ok

if __name__ == "__main__":
    sys.exit(0 if self_test() else 1)
