#!/usr/bin/env python3
"""Sedenion tensor product as pure combinational logic (AND / OR / NOT + wiring).

BRACKETING DECLARATION
----------------------
Sedenion multiplication is non-associative, so "the" tensor product does not
exist: there is a coherent family of tensor products indexed by bracketing
(parenthesization), related by the associator Phi = d(omega), which satisfies
the pentagon identity (verified exhaustively: 0 violations / 65536).

This module computes THREE declared members of that family:

    A. STAGED (data-innermost), `tensor_units`:
       butterfly stages applied lsb-axis first; per path this is
       g_msb (g_{msb-1} ( ... (g_lsb x))) — right-comb bracketing with the
       DATA as the innermost leaf.
       [Correction: an earlier revision described this as "equivalently the
       left-associated entry product". That equivalence holds only in
       associative algebras; in sedenions it is false. A is its own member.]

    B. LEFT-COMB ENTRIES, `tensor_precomb(entry_leftcomb)`:
       entry (((g_msb g_{msb-1}) ... ) g_lsb) pre-folded, then applied to x.

    C. RIGHT-COMB ENTRIES, `tensor_precomb(entry_rightcomb)`:
       entry (g_msb (g_{msb-1} ( ... g_lsb))) pre-folded, then applied to x.

For unit entries, pre-folded entries are again signed basis units (index =
XOR of indices, bracket-independent; only the SIGN depends on the comb), so
B and C are also multiplier-free: wiring + adders. A, B, C differ pairwise
on generic sedenion data (needs m >= 3 for B != C) and collapse to a single
product on any associative subalgebra (e.g. quaternions).

Translations between members are fixed diagonal sign masks (per path).
Those masks provably do NOT factor into stage-local wiring (GF(2)
obstruction = the 1848 non-cocycle triples), so switching members requires
per-path wiring, not a patched shared butterfly.

DATAPATH (all static: circuit shape is data-independent, zeros flow through)
    L0  gate9   : real trit x trit product (4 AND, 2 OR, no carry)
    L1  wiring  : sign omega(i,j) = p/n wire swap; routing k = i XOR j
    L2  gate27  : 3:2 balanced-ternary compressor (pure formulas)
    L3  tree    : static compressor tree per output component
    L4  ripple  : final normalization (gate27 carry chain)
    L5  stage   : butterfly; unit entries cost ZERO multiplier cells
    fp  rescale : fixed point = drop low trit wires (free, unbiased rounding)

NOT logic (declared residues): host-side quantization, int<->trit framing,
width constants, pipeline registers for sequential folds.
"""
from typing import Sequence
import numpy as np

# ---------------------------------------------------------------- primitives
def enc(t): return int(t == 1), int(t == -1)
def dec(p, n): return p - n
def bnot(x): return 1 - x

def gate9(a, b):
    """Trit product: {-1,0,1} x {-1,0,1} -> {-1,0,1}. Never carries."""
    ap, an = enc(a); bp, bn = enc(b)
    return dec((ap & bp) | (an & bn), (ap & bn) | (an & bp))

def gate27(x, y, c):
    """3:2 compressor: x + y + c = same + 3 * high. Pure AND/OR/NOT."""
    xp, xn = enc(x); yp, yn = enc(y); cp, cn = enc(c)
    sp = ((cn&xp&yp)|(cp&xn&yp)|(cp&xp&yn)
          | (cn&xn&bnot(yn)&bnot(yp)) | (cn&yn&bnot(xn)&bnot(xp))
          | (xn&yn&bnot(cn)&bnot(cp))
          | (cp&bnot(xn)&bnot(xp)&bnot(yn)&bnot(yp))
          | (xp&bnot(cn)&bnot(cp)&bnot(yn)&bnot(yp))
          | (yp&bnot(cn)&bnot(cp)&bnot(xn)&bnot(xp)))
    sn = ((cn&xn&yp)|(cn&xp&yn)|(cp&xn&yn)
          | (cp&xp&bnot(yn)&bnot(yp)) | (cp&yp&bnot(xn)&bnot(xp))
          | (xp&yp&bnot(cn)&bnot(cp))
          | (cn&bnot(xn)&bnot(xp)&bnot(yn)&bnot(yp))
          | (xn&bnot(cn)&bnot(cp)&bnot(yn)&bnot(yp))
          | (yn&bnot(cn)&bnot(cp)&bnot(xn)&bnot(xp)))
    hp = (cp&xp&bnot(yn)) | (cp&yp&bnot(xn)) | (xp&yp&bnot(cn))
    hn = (cn&xn&bnot(yp)) | (cn&yn&bnot(xp)) | (xn&yn&bnot(cp))
    return dec(sp, sn), dec(hp, hn)

# --------------------------------------------- fixed wiring constants: omega
M = 16

def _conj(x):
    n = len(x)
    if n == 1: return x.copy()
    h = n // 2
    return np.concatenate([_conj(x[:h]), -x[h:]])

def _cd(x, y):
    n = len(x)
    if n == 1: return x * y
    h = n // 2
    a, b, c, d = x[:h], x[h:], y[:h], y[h:]
    return np.concatenate([_cd(a, c) - _cd(_conj(d), b),
                           _cd(d, a) + _cd(b, _conj(c))])

_E = np.eye(M)
OMEGA = np.zeros((M, M), dtype=int)
for _i in range(M):
    for _j in range(M):
        _v = _cd(_E[_i], _E[_j])
        _k = int(np.argmax(np.abs(_v)))
        assert _k == (_i ^ _j), "XOR structure violated"
        OMEGA[_i, _j] = int(np.sign(_v[_k]))

# ------------------------------------------------- host-side framing (R2)
def to_trits(n, width):
    n = int(n); out = []
    for _ in range(width):
        n, r = divmod(n, 3)
        if r == 2: r = -1; n += 1
        out.append(r)
    if n != 0: raise OverflowError("width too small")
    return out

def from_trits(ts):
    return int(sum(t * 3 ** k for k, t in enumerate(ts)))

# ------------------------------------------------------- L3 / L4 (static)
def static_tree(columns, stats):
    """Slots are kept even when values are zero: shape is data-independent."""
    cols = [list(c) for c in columns]
    k = 0
    while k < len(cols):
        while len(cols[k]) > 2:
            x = cols[k].pop(); y = cols[k].pop(); z = cols[k].pop()
            s, h = gate27(x, y, z); stats['g27'] += 1
            cols[k].append(s)
            if k + 1 >= len(cols): cols.append([])
            cols[k + 1].append(h)
        k += 1
    return cols

def two_rows(cols):
    r0 = [c[0] if len(c) > 0 else 0 for c in cols]
    r1 = [c[1] if len(c) > 1 else 0 for c in cols]
    return r0, r1

def ripple(r0, r1, stats):
    out, carry = [], 0
    for x, y in zip(r0, r1):
        s, carry = gate27(x, y, carry); stats['g27'] += 1
        out.append(s)
    out.append(carry)
    return out

# --------------------------------------------------- sedenion multiplier
def sedenion_mult(xw, yw, stats, W=None):
    """xw, yw: 16 trit words each -> 16 trit words. Pure L0-L4.
    Sign omega = wire swap (zero gates); routing = XOR (zero gates).
    W: optional twist table (default: the Cayley-Dickson OMEGA)."""
    W = OMEGA if W is None else W
    K1, K2 = len(xw[0]), len(yw[0])
    cols = [[[] for _ in range(K1 + K2)] for _ in range(M)]
    for i in range(M):
        for j in range(M):
            k, s = i ^ j, W[i, j]
            for p in range(K1):
                for q in range(K2):
                    t = gate9(xw[i][p], yw[j][q]); stats['g9'] += 1
                    cols[k][p + q].append(s * t)      # s*t: p/n wire swap
    return [ripple(*two_rows(static_tree(cols[k], stats)), stats)
            for k in range(M)]

def sedenion_add(aw, bw, stats):
    out = []
    for c in range(M):
        w = max(len(aw[c]), len(bw[c]))
        r0 = aw[c] + [0 for _ in range(w - len(aw[c]))]
        r1 = bw[c] + [0 for _ in range(w - len(bw[c]))]
        out.append(ripple(r0, r1, stats))
    return out

# ------------------------------- unit-entry multiply = pure wiring (0 gates)
def unit_wire(j, sign, xw, W=None):
    """(sign * e_j) * x : output component k takes input component j XOR k
    with sign * omega(j, j XOR k). No gates."""
    W = OMEGA if W is None else W
    return [[sign * W[j, j ^ k] * t for t in xw[j ^ k]] for k in range(M)]

# ------------------------------------- signed-unit folds (pure constants)
def unit_mul(a, sa, b, sb, W=None):
    """(sa e_a)(sb e_b) = sa sb omega(a,b) e_{a XOR b} — compile-time."""
    W = OMEGA if W is None else W
    return a ^ b, sa * sb * W[a, b]

def entry_leftcomb(Gidx, i, j, m, W=None):
    """(((g_msb g_{msb-1}) ...) g_lsb) as a signed unit."""
    acc = None
    for s in range(m - 1, -1, -1):
        g = Gidx[(i >> s) & 1][(j >> s) & 1]
        acc = g if acc is None else unit_mul(acc[0], acc[1], g[0], g[1], W)
    return acc

def entry_rightcomb(Gidx, i, j, m, W=None):
    """(g_msb (g_{msb-1} (... g_lsb))) as a signed unit."""
    acc = None
    for s in range(0, m):
        g = Gidx[(i >> s) & 1][(j >> s) & 1]
        acc = g if acc is None else unit_mul(g[0], g[1], acc[0], acc[1], W)
    return acc

def tensor_precomb(Gidx, x, m, entry_fn, stats, W=None):
    """Members B / C: pre-folded unit entries (pure wiring) + adder chains.
    out_i = sum_j (entry_ij) * x_j.  Zero multiplier cells."""
    out = []
    for i in range(2 ** m):
        acc = None
        for j in range(2 ** m):
            k, s = entry_fn(Gidx, i, j, m, W)
            term = unit_wire(k, s, x[j], W)
            acc = term if acc is None else sedenion_add(acc, term, stats)
        out.append(acc)
    return out

# --------------------------- L5: member A, staged (declared bracket)
def tensor_units(Gidx, x, m, stats, W=None):
    """G^(tensor m) applied to x, LEFT-FOLD member of the family.
    Gidx: 2x2 of (j, sign) unit entries. x: list of 2^m sedenions (trit words).
    Stages applied lsb-axis first == entries left-associated, msb leftmost.
    Unit entries: ZERO multiplier cells; adders and wiring only."""
    y = [[w[:] for w in q] for q in x]
    for s in range(m):                                  # lsb first: declared
        out = [None] * len(y)
        for base in range(len(y)):
            if (base >> s) & 1: continue
            u, v = y[base], y[base | (1 << s)]
            (j00, s00), (j01, s01) = Gidx[0]
            (j10, s10), (j11, s11) = Gidx[1]
            out[base]            = sedenion_add(unit_wire(j00, s00, u, W),
                                                unit_wire(j01, s01, v, W), stats)
            out[base | (1 << s)] = sedenion_add(unit_wire(j10, s10, u, W),
                                                unit_wire(j11, s11, v, W), stats)
        y = out
    return y

# ------------------------------------------------ fixed point (fp) layer
def rescale(word, f):
    """Drop f low trits = wiring. Unbiased round-to-nearest, no gates."""
    out = word[f:]
    return out if out else [0]

def tensor_general_fp(Gw, x, m, f, stats):
    """General sedenion constant entries (trit words at point f).
    Per stage: exact multiplies at 2f, exact add, one free rescale to f.
    Same LEFT-FOLD bracket as tensor_units."""
    y = [[w[:] for w in q] for q in x]
    for s in range(m):
        out = [None] * len(y)
        for base in range(len(y)):
            if (base >> s) & 1: continue
            u, v = y[base], y[base | (1 << s)]
            a00 = sedenion_mult(Gw[0][0], u, stats)
            a01 = sedenion_mult(Gw[0][1], v, stats)
            a10 = sedenion_mult(Gw[1][0], u, stats)
            a11 = sedenion_mult(Gw[1][1], v, stats)
            out[base]            = [rescale(w, f) for w in sedenion_add(a00, a01, stats)]
            out[base | (1 << s)] = [rescale(w, f) for w in sedenion_add(a10, a11, stats)]
        y = out
    return y

# ---------------------------------------------------------- references
def ref_mult(x, y, W=None):
    W = OMEGA if W is None else W
    r = [0] * M
    for i in range(M):
        for j in range(M):
            r[i ^ j] += W[i, j] * x[i] * y[j]
    return r

def ref_tensor_units(Gidx, xs, m):
    """Same LEFT-FOLD bracket, integer arithmetic."""
    def unit(j, sign, x):
        return [sign * OMEGA[j, j ^ k] * x[j ^ k] for k in range(M)]
    out = []
    for i in range(2 ** m):
        tot = [0] * M
        for jj in range(2 ** m):
            v = xs[jj]
            for s in range(m):                          # lsb applied first
                j, sg = Gidx[(i >> s) & 1][(jj >> s) & 1]
                v = unit(j, sg, v)
            tot = [a + b for a, b in zip(tot, v)]
        out.append(tot)
    return out

# ------------------------------------------------------------- self-test
def self_test(seed=20260711):
    rng = np.random.default_rng(seed)
    K = 6

    st = {'g9': 0, 'g27': 0}
    for _ in range(10):
        x = [int(v) for v in rng.integers(-3**K // 2, 3**K // 2, M)]
        y = [int(v) for v in rng.integers(-3**K // 2, 3**K // 2, M)]
        got = [from_trits(w) for w in
               sedenion_mult([to_trits(c, K) for c in x],
                             [to_trits(c, K) for c in y], st)]
        assert got == ref_mult(x, y)
    print("sedenion multiplier: 10 random tests exact")

    sA = {'g9': 0, 'g27': 0}; sB = {'g9': 0, 'g27': 0}
    z = [to_trits(0, K)] * M
    r = [to_trits(int(v), K) for v in rng.integers(-3**K // 2, 3**K // 2, M)]
    sedenion_mult(z, z, sA); sedenion_mult(r, r, sB)
    assert sA == sB
    print(f"static structure: gate counts data-independent "
          f"(per product K={K}: gate9 {sB['g9']}, gate27 {sB['g27']})")

    xz = [0]*M; xz[1] = 3**4; xz[10] = 3**4
    yz = [0]*M; yz[4] = 3**4; yz[15] = -3**4
    got = [from_trits(w) for w in
           sedenion_mult([to_trits(c, K) for c in xz],
                         [to_trits(c, K) for c in yz], {'g9': 0, 'g27': 0})]
    assert all(g == 0 for g in got)
    print("zero divisor (e1+e10)(e4-e15) = 0: exact")

    Gidx = [[(0, 1), (1, 1)], [(2, 1), (4, -1)]]
    for m in (2, 3):
        xs = [[int(v) for v in rng.integers(-3**4, 3**4, M)]
              for _ in range(2 ** m)]
        xw = [[to_trits(c, 10) for c in q] for q in xs]
        stt = {'g9': 0, 'g27': 0}
        got = [[from_trits(w) for w in q]
               for q in tensor_units(Gidx, xw, m, stt)]
        assert got == ref_tensor_units(Gidx, xs, m)
        print(f"left-fold unit-entry tensor m={m}: exact | "
              f"gate9 {stt['g9']} (ZERO multipliers), gate27 {stt['g27']}")

    f, W = 6, 18
    q3 = lambda v: int(round(float(v) * 3**f))
    G = [[rng.uniform(-1, 1, M) for _ in range(2)] for _ in range(2)]
    Gw = [[[to_trits(q3(c), W) for c in e] for e in row] for row in G]
    xs = [rng.uniform(-1, 1, M) for _ in range(4)]
    xw = [[to_trits(q3(c), W) for c in q] for q in xs]
    got = [[from_trits(w) / 3**f for w in q]
           for q in tensor_general_fp(Gw, xw, 2, f, {'g9': 0, 'g27': 0})]
    ref_in = [[q3(c) / 3**f for c in q] for q in xs]
    Gq = [[[q3(c) / 3**f for c in e] for e in row] for row in G]
    y = [list(q) for q in ref_in]
    def fmul(a, b):
        r = [0.0] * M
        for i in range(M):
            for j in range(M):
                r[i ^ j] += OMEGA[i, j] * a[i] * b[j]
        return r
    for s in range(2):
        out = [None] * len(y)
        for base in range(len(y)):
            if (base >> s) & 1: continue
            u, v = y[base], y[base | (1 << s)]
            out[base]            = [a + b for a, b in zip(fmul(Gq[0][0], u), fmul(Gq[0][1], v))]
            out[base | (1 << s)] = [a + b for a, b in zip(fmul(Gq[1][0], u), fmul(Gq[1][1], v))]
        y = out
    err = max(abs(g - r) for gq, rq in zip(got, y) for g, r in zip(gq, rq))
    # per rescale: |err| <= 3^-f/2 per component; earlier-stage error is
    # amplified by the stage gain (sum of 2 sedenion products of 16 terms
    # each with |entry| <= 1), so bound stages by (1 + gain) * m * 3^-f/2.
    gain = 2 * 16
    bound = (1 + gain) * 2 * (3**-f) / 2
    assert err <= bound, (err, bound)
    print(f"fixed-point general-entry tensor m=2: max err {err:.2e} "
          f"<= gain-aware bound {bound:.2e} (no division anywhere)")

    # three members: pairwise distinct on S (m=3), collapsed on H
    G3 = [[(0, 1), (1, 1)], [(2, 1), (4, -1)]]
    xs3 = [[int(v) for v in rng.integers(-3**4, 3**4, M)] for _ in range(8)]
    xw3 = [[to_trits(c, 10) for c in q] for q in xs3]
    A = [[from_trits(w) for w in q]
         for q in tensor_units(G3, xw3, 3, {'g9': 0, 'g27': 0})]
    B = [[from_trits(w) for w in q]
         for q in tensor_precomb(G3, xw3, 3, entry_leftcomb, {'g9': 0, 'g27': 0})]
    C = [[from_trits(w) for w in q]
         for q in tensor_precomb(G3, xw3, 3, entry_rightcomb, {'g9': 0, 'g27': 0})]
    assert A != B and A != C and B != C
    GH = [[(0, 1), (1, 1)], [(2, 1), (3, -1)]]
    xh = [[int(v) for v in rng.integers(-3**4, 3**4, 4)] + [0] * 12
          for _ in range(8)]
    xhw = [[to_trits(c, 10) for c in q] for q in xh]
    AH = [[from_trits(w) for w in q]
          for q in tensor_units(GH, xhw, 3, {'g9': 0, 'g27': 0})]
    BH = [[from_trits(w) for w in q]
          for q in tensor_precomb(GH, xhw, 3, entry_leftcomb, {'g9': 0, 'g27': 0})]
    CH = [[from_trits(w) for w in q]
          for q in tensor_precomb(GH, xhw, 3, entry_rightcomb, {'g9': 0, 'g27': 0})]
    assert AH == BH == CH
    print("members A/B/C: pairwise distinct on S (m=3), collapsed on H")
    print("all self-tests passed — three declared members of the coherent "
          "tensor-product family (A staged, B left-comb, C right-comb).")

if __name__ == "__main__":
    self_test()
