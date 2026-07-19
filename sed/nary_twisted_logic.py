#!/usr/bin/env python3
"""Twisted group algebra over (Z/n)^k with phase twists — exact integer core.

This generalizes the balanced-ternary sedenion line from (Z/2)^4 with sign
twists to (Z/n)^k with n-th-root-of-unity phase twists. Focus: n = 3.

WHAT GENERALIZES (verified in self_test):
  - routing: XOR -> digit-wise addition mod n (fixed wiring per index pair)
  - twist:   signs {+-1} -> phases zeta^t (zeta = primitive n-th root)
  - associativity  <=>  the twist is a 2-cocycle (d omega = 0)
  - coherence: pentagon identity holds for EVERY twist (d(d omega) = 0),
    so any twist yields a coherent family of bracket-indexed products
  - the three declared members A (staged, data-innermost), B (left-comb
    entries), C (right-comb entries): distinct for non-cocycle twists,
    collapsed to one product for cocycle (associative) twists

EXACTNESS / HARDWARE NOTES:
  - values are Eisenstein integers a + b*zeta (zeta^2 = -1 - zeta) stored as
    integer pairs; no floats anywhere in the algebra
  - multiplying by zeta^t is the integer matrix [[0,-1],[1,-1]]^t:
    t=1: (a,b) -> (-b, a-b); t=2: (a,b) -> (b-a, -a) — negations and one
    addition, i.e. adder-level cost, NO multiplier; hence unit-entry tensor
    transforms remain multiplier-free, as in the n = 2 (sedenion) case
  - indices over (Z/3)^k are trit vectors: native to balanced-ternary
    datapaths; the balanced digit set {-1,0,+1} is Z/3 written so that
    every element has a pair (-x), i.e. the pair principle survives odd n
  - the standard symplectic cocycle on (Z/3)^(2q) generates the q-qutrit
    Pauli (clock & shift) algebra = full matrix algebra M_{3^q}

BRACKETING: products of >= 3 factors require a declared bracket, exactly as
in the sedenion module. Members A/B/C below carry their declaration.
"""
from itertools import product as iproduct
import numpy as np

# ------------------------------------------------------------ Eisenstein Z[zeta]
def eis_add(u, v):
    return (u[0] + v[0], u[1] + v[1])

def eis_neg(u):
    return (-u[0], -u[1])

def eis_mul(u, v):
    """general product (a+b z)(c+d z), z^2 = -1-z. Integer multiplies:
    this is the 'general entry' path (analogue of gate9 arrays)."""
    a, b = u; c, d = v
    return (a * c - b * d, a * d + b * c - b * d)

def phase(t, u):
    """multiply by zeta^t — ADDER-LEVEL ONLY (no multiplier), t in {0,1,2}."""
    a, b = u
    t %= 3
    if t == 0: return (a, b)
    if t == 1: return (-b, a - b)
    return (b - a, -a)

ZERO, ONE = (0, 0), (1, 0)

# ------------------------------------------------------------ group (Z/n)^k
class Group:
    def __init__(self, n, k):
        self.n, self.k = n, k
        self.elems = [tuple(t) for t in iproduct(range(n), repeat=k)]
        self.idx = {g: i for i, g in enumerate(self.elems)}
        self.size = n ** k
    def add(self, g, h):
        return tuple((x + y) % self.n for x, y in zip(g, h))
    def neg(self, g):
        return tuple((-x) % self.n for x in g)

def omega_symplectic(G):
    """standard 2-cocycle on (Z/n)^(2q): pairs axes (0,1),(2,3),...
    omega(g,h) = zeta^{sum g[2i] * h[2i+1]}  -> associative (clock&shift)."""
    q = G.k // 2
    W = np.zeros((G.size, G.size), dtype=int)
    for g in G.elems:
        for h in G.elems:
            W[G.idx[g], G.idx[h]] = sum(g[2*i] * h[2*i+1] for i in range(q)) % G.n
    return W

def omega_random(G, seed=0):
    rng = np.random.default_rng(seed)
    W = rng.integers(0, G.n, (G.size, G.size))
    W[0, :] = 0; W[:, 0] = 0          # e_0 stays the two-sided unit
    return W

# --------------------------------------------- algebra element multiplication
def alg_mult(G, W, x, y):
    """x, y: lists of Eisenstein pairs indexed by G. General product."""
    out = [ZERO] * G.size
    for gi, g in enumerate(G.elems):
        if x[gi] == ZERO: pass
        for hi, h in enumerate(G.elems):
            t = W[gi, hi]
            p = eis_mul(x[gi], y[hi])
            ki = G.idx[G.add(g, h)]
            out[ki] = eis_add(out[ki], phase(t, p))
    return out

# ------------------------------------------------- structure diagnostics
def assoc_defects(G, W):
    d = 0
    for a in range(G.size):
        for b in range(G.size):
            ab = G.idx[G.add(G.elems[a], G.elems[b])]
            for c in range(G.size):
                bc = G.idx[G.add(G.elems[b], G.elems[c])]
                if (W[a, b] + W[ab, c] - W[b, c] - W[a, bc]) % G.n:
                    d += 1
    return d

def pentagon_defects(G, W):
    n = G.n
    def phi(a, b, c):
        ab = G.idx[G.add(G.elems[a], G.elems[b])]
        bc = G.idx[G.add(G.elems[b], G.elems[c])]
        return (W[a, b] + W[ab, c] - W[b, c] - W[a, bc]) % n
    d = 0
    for a in range(G.size):
        for b in range(G.size):
            ab = G.idx[G.add(G.elems[a], G.elems[b])]
            for c in range(G.size):
                bc = G.idx[G.add(G.elems[b], G.elems[c])]
                for e in range(G.size):
                    ce = G.idx[G.add(G.elems[c], G.elems[e])]
                    if (phi(b, c, e) + phi(a, bc, e) + phi(a, b, c)
                            - phi(ab, c, e) - phi(a, b, ce)) % n:
                        d += 1
    return d

# --------------------------------- signed(phased)-unit folds (compile-time)
def unit_mul(G, W, g, t, h, s):
    """(zeta^t e_g)(zeta^s e_h) = zeta^{t+s+omega(g,h)} e_{g+h}."""
    return G.add(g, h), (t + s + W[G.idx[g], G.idx[h]]) % G.n

def entry_leftcomb(G, W, Gidx, i, j, m):
    acc = None
    for s in range(m - 1, -1, -1):
        g = Gidx[(i >> s) & 1][(j >> s) & 1]
        acc = g if acc is None else unit_mul(G, W, acc[0], acc[1], g[0], g[1])
    return acc

def entry_rightcomb(G, W, Gidx, i, j, m):
    acc = None
    for s in range(0, m):
        g = Gidx[(i >> s) & 1][(j >> s) & 1]
        acc = g if acc is None else unit_mul(G, W, g[0], g[1], acc[0], acc[1])
    return acc

# ---------------------------------------- phased-unit action = wiring + phase
def unit_wire(G, W, g, t, x):
    """(zeta^t e_g) * x : component (g+h) receives zeta^{t+omega(g,h)} x[h].
    Permutation wiring + adder-level phase. NO multiplier."""
    out = [ZERO] * G.size
    for hi, h in enumerate(G.elems):
        ki = G.idx[G.add(g, h)]
        out[ki] = phase(t + W[G.idx[g], hi], x[hi])
    return out

def vec_add(x, y):
    return [eis_add(a, b) for a, b in zip(x, y)]

# ------------------------------ tensor-product members (declared brackets)
def tensor_units_staged(G, W, Gidx, x, m):
    """member A: butterfly stages lsb-first, data innermost:
    g_msb ( g_{msb-1} ( ... (g_lsb x))). Multiplier-free."""
    y = [list(q) for q in x]
    for s in range(m):
        out = [None] * len(y)
        for base in range(len(y)):
            if (base >> s) & 1: continue
            u, v = y[base], y[base | (1 << s)]
            (g00, t00), (g01, t01) = Gidx[0]
            (g10, t10), (g11, t11) = Gidx[1]
            out[base]            = vec_add(unit_wire(G, W, g00, t00, u),
                                           unit_wire(G, W, g01, t01, v))
            out[base | (1 << s)] = vec_add(unit_wire(G, W, g10, t10, u),
                                           unit_wire(G, W, g11, t11, v))
        y = out
    return y

def tensor_units_precomb(G, W, Gidx, x, m, entry_fn):
    """members B / C: entries pre-folded (still phased units -> wiring)."""
    out = []
    for i in range(2 ** m):
        acc = None
        for j in range(2 ** m):
            g, t = entry_fn(G, W, Gidx, i, j, m)
            term = unit_wire(G, W, g, t, x[j])
            acc = term if acc is None else vec_add(acc, term)
        out.append(acc)
    return out

# ------------------------------------------------------------- self-test
def self_test():
    n = 3
    G = Group(n, 2)                              # (Z/3)^2 : one qutrit's worth
    Ws = omega_symplectic(G)

    print(f"(Z/{n})^2, standard symplectic cocycle:")
    print("  associativity defects:", assoc_defects(G, Ws), "/", G.size ** 3)
    assert assoc_defects(G, Ws) == 0

    # clock & shift: U = e_(1,0), V = e_(0,1);  U V = zeta^? V U
    U = [ZERO] * 9; U[G.idx[(1, 0)]] = ONE
    V = [ZERO] * 9; V[G.idx[(0, 1)]] = ONE
    UV = alg_mult(G, Ws, U, V)
    VU = alg_mult(G, Ws, V, U)
    zVU = [phase(1, c) for c in VU]
    assert UV == zVU
    print("  clock-shift relation U V = zeta V U: True (qutrit Pauli algebra)")

    # left-multiplication operators span dim 9 -> full matrix algebra M_3
    def eis_to_c(u):
        a, b = u
        return a + b * np.exp(2j * np.pi / 3)
    Ls = []
    for gi in range(9):
        e = [ZERO] * 9; e[gi] = ONE
        cols = []
        for hi in range(9):
            b = [ZERO] * 9; b[hi] = ONE
            cols.append([eis_to_c(u) for u in alg_mult(G, Ws, e, b)])
        Ls.append(np.array(cols).T.ravel())
    rank = np.linalg.matrix_rank(np.stack(Ls))
    print("  algebra dimension (should be 9 = M_3):", rank)
    assert rank == 9

    Wr = omega_random(G, seed=0)
    ad = assoc_defects(G, Wr)
    pd = pentagon_defects(G, Wr)
    print(f"random phase twist: assoc defects {ad} (>0), pentagon defects {pd} / {G.size**4}")
    assert ad > 0 and pd == 0

    # zeta multiplication is adder-level: matches complex reference
    z = np.exp(2j * np.pi / 3)
    for t in range(3):
        a, b = 7, -4
        got = phase(t, (a, b))
        ref = (a + b * z) * z ** t
        assert np.isclose(complex(got[0] + got[1] * np.cos(2*np.pi/3),
                                  got[1] * np.sin(2*np.pi/3)), ref)
    print("phase(t) integer closed forms == zeta^t (adder-level, no multiplier)")

    # tensor members A/B/C: collapse under cocycle, split under random twist
    rng = np.random.default_rng(1)
    Gidx = [[((0, 0), 0), ((1, 0), 0)], [((0, 1), 0), ((1, 1), 1)]]
    m = 3
    xs = [[(int(rng.integers(-9, 9)), int(rng.integers(-9, 9)))
           for _ in range(9)] for _ in range(2 ** m)]
    for name, W in (("cocycle (assoc)", Ws), ("random twist   ", Wr)):
        A = tensor_units_staged(G, W, Gidx, xs, m)
        B = tensor_units_precomb(G, W, Gidx, xs, m,
                                 lambda G_, W_, Gi, i, j, mm: entry_leftcomb(G_, W_, Gi, i, j, mm))
        C = tensor_units_precomb(G, W, Gidx, xs, m,
                                 lambda G_, W_, Gi, i, j, mm: entry_rightcomb(G_, W_, Gi, i, j, mm))
        rel = "A == B == C (family collapsed)" if A == B == C else \
              f"pairwise: A==B {A==B}, A==C {A==C}, B==C {B==C} (family split)"
        print(f"tensor members under {name}: {rel}")
    print("all self-tests passed.")

if __name__ == "__main__":
    self_test()
