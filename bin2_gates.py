#!/usr/bin/env python3
"""Binary {0,1} carry-save arithmetic in the same gate style as sd2_gates.py — the two's-complement counterpart of
the signed-digit rails.

NUMBER FORMAT
  A value is W bits of two's complement.  The *redundant* form is two rows (r0, r1) of W bits each with
  value (r0 + r1) mod 2^W, read as two's complement — "carry-save".  Nothing carries across bit positions
  until `resolve` (one Kogge–Stone adder at the boundary), exactly as the signed-digit datapath only
  canonicalises at its boundary.

OPERATIONS (each a gate graph; works on ints 0/1 for the golden model, on gate_fast.B for depth counting,
and on rtl/emit_sv.T for tracing to SystemVerilog)
  mul_cs(X, Y, st)          W×W two's complement → (r0, r1) of 2W bits.  Baugh–Wooley partial products
                            (the two sign-weighted rows are complemented, +2^W + 2^(2W-1) as constants),
                            Dadda reduction with full_adder / half adder to two rows.  No carry-propagate.
  add_cs(rows, W, st)       any number of W-bit rows → two rows (Dadda).  This is the accumulator: a
                            multiply-accumulate chain never leaves carry-save form.
  neg_cs((r0, r1), W, st)   two's complement negation of a carry-save number: (~r0, ~r1, +2).
  resolve((r0, r1), st)     Kogge–Stone (gate_fast.bin_add_fast) → W bits.  The only O(log W)-depth carry.

WHY (measured in the equivalent-logic-physical-cost repository, sky130, cell-only STA, equal ±1024 range):
  binary 11×11 Dadda+KS 770 cells / 5485 µm² / 3.0 ns  vs  signed-digit 10-digit two-rail KS 1610 / 11228 / 3.1 ns;
  tree-only (carry-save out) 484 / 3505 / 1.8 ns  vs  signed-digit tree-only 1140 / 7965 / 1.9 ns.
  Same carry-free accumulation, half the area.  What signed digits buy is representation (negation by
  wiring, symmetric digits, flag semantics), not speed or area.
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate_bilinear import AND, OR, NOT, XOR, new_counter, full_adder

# ------------------------------------------------------------------ host-side framing
def to_bits(v, W):
    """integer (two's complement range) → W bits, LSB first, as ints 0/1"""
    v = int(v); assert -(1 << (W - 1)) <= v < (1 << (W - 1)), f"{v} does not fit {W} bits"
    return [(v >> i) & 1 for i in range(W)]

def from_bits(bits):
    """W bits (LSB first; ints or objects with .v) → signed integer"""
    W = len(bits); u = sum(_bv(b) << i for i, b in enumerate(bits))
    return u - (1 << W) if u >> (W - 1) else u

def cs_val(rows, W):
    """carry-save rows → signed integer, (Σ rows) mod 2^W read as two's complement"""
    u = sum(sum(_bv(b) << i for i, b in enumerate(r)) for r in rows) % (1 << W)
    return u - (1 << W) if u >> (W - 1) else u

def _bv(b): return int(b.v) if hasattr(b, 'v') else int(b)

# ------------------------------------------------------------------ half adder (2 gates)
def half_adder(a, b, st):
    return XOR(a, b, st), AND(a, b, st)

# ------------------------------------------------------------------ Dadda reduction to two rows
def dadda_rows(cols, W, st):
    """cols: list (length W) of lists of bits at weight 2^k.  Reduce every column to ≤ 2 bits with the
    Dadda height sequence (2, 3, 4, 6, 9, …), carries going one column up; bits beyond weight W-1 are
    dropped (arithmetic mod 2^W).  Returns (row0, row1), each W bits (constant 0 where empty)."""
    cols = [list(c) for c in cols] + [[]]
    height = max(len(c) for c in cols)
    targets = []; t = 2
    while t < height: targets.append(t); t = t * 3 // 2
    for t in reversed(targets):
        for k in range(W):
            while len(cols[k]) > t:
                if len(cols[k]) - t >= 2:
                    a, b, c = cols[k].pop(0), cols[k].pop(0), cols[k].pop(0)
                    s, cy = full_adder(a, b, c, st)
                else:
                    a, b = cols[k].pop(0), cols[k].pop(0)
                    s, cy = half_adder(a, b, st)
                cols[k].append(s)
                if k + 1 < W: cols[k + 1].append(cy)
    row0 = [c[0] if len(c) > 0 else 0 for c in cols[:W]]
    row1 = [c[1] if len(c) > 1 else 0 for c in cols[:W]]
    return row0, row1

# ------------------------------------------------------------------ operations
def mul_cs(X, Y, st):
    """two's complement X (Wx bits) × Y (Wy bits) → carry-save rows of Wx+Wy bits (Baugh–Wooley)."""
    Wx, Wy = len(X), len(Y); W = Wx + Wy
    cols = [[] for _ in range(W)]
    for i, x in enumerate(X):
        for j, y in enumerate(Y):
            t = AND(x, y, st)
            if (i == Wx - 1) != (j == Wy - 1):      # exactly one sign bit: negative weight → complement
                t = NOT(t, st)
            cols[i + j].append(t)
    cols[Wx].append(1); cols[W - 1].append(1)       # Baugh–Wooley correction: +2^Wx + 2^(W-1)  (mod 2^W)
    return dadda_rows(cols, W, st)

def add_cs(rows, W, st):
    """merge any number of W-bit rows (two's complement, sign bits already extended to W) → two rows."""
    cols = [[] for _ in range(W)]
    for r in rows:
        for k, b in enumerate(r[:W]):
            if not (isinstance(b, int) and b == 0): cols[k].append(b)
    return dadda_rows(cols, W, st)

def sign_extend(row, W, st=None):
    """extend a shorter two's complement row to W bits (copies the sign bit)"""
    return list(row) + [row[-1]] * (W - len(row))

def neg_cs(rows, W, st):
    """−(r0 + r1) = ~r0 + ~r1 + 2  (mod 2^W): two complemented rows plus the constant 2 as a third row."""
    r0, r1 = rows
    n0 = [NOT(b, st) if not isinstance(b, int) else (1 - b) for b in r0]
    n1 = [NOT(b, st) if not isinstance(b, int) else (1 - b) for b in r1]
    two = [0, 1] + [0] * (W - 2)
    return add_cs([n0, n1, two], W, st)

def resolve(rows, st):
    """carry-save → binary: one Kogge–Stone adder (the only carry chain, O(log W) depth)."""
    from gate_fast import bin_add_fast
    r0, r1 = rows
    return bin_add_fast(list(r0), list(r1), st)

def mac_cs(X, Y, acc, st):
    """multiply-accumulate without a carry chain: acc (two rows, W bits) += X·Y (sign-extended)."""
    W = len(acc[0]); p0, p1 = mul_cs(X, Y, st)
    return add_cs([acc[0], acc[1], sign_extend(p0, W), sign_extend(p1, W)], W, st)

# ------------------------------------------------------------------ self-test (golden ints) + gate/depth report
def self_test(seed=20260906):
    rng = random.Random(seed)
    for W in (8, 11):
        lim = 1 << (W - 1)
        for _ in range(500):
            a = rng.randrange(-lim, lim); b = rng.randrange(-lim, lim)
            st = new_counter()
            rows = mul_cs(to_bits(a, W), to_bits(b, W), st)
            assert cs_val(rows, 2 * W) == a * b, (a, b, cs_val(rows, 2 * W))
            assert from_bits(resolve(rows, st)) == a * b
            n = neg_cs(rows, 2 * W, st); assert cs_val(n, 2 * W) == -a * b
        # dot products of K terms stay in carry-save form, one resolve at the end
        for K in (2, 4, 8):
            for _ in range(100):
                xs = [rng.randrange(-lim, lim) for _ in range(K)]; ys = [rng.randrange(-lim, lim) for _ in range(K)]
                Wa = 2 * W + 4; st = new_counter()
                acc = ([0] * Wa, [0] * Wa)
                for x, y in zip(xs, ys): acc = mac_cs(to_bits(x, W), to_bits(y, W), acc, st)
                want = sum(x * y for x, y in zip(xs, ys))
                assert cs_val(acc, Wa) == want and from_bits(resolve(acc, st)) == want
    print("bin2_gates self_test: mul_cs / resolve / neg_cs / mac_cs exact on random two's complement inputs ✓")
    # gate counts and depth with gate_fast.B (depth = longest gate path)
    from gate_fast import B, depth_of
    for W in (8, 11, 16, 24):
        st = new_counter(); X = [B(0) for _ in range(W)]; Y = [B(0) for _ in range(W)]
        rows = mul_cs(X, Y, st); d_tree = depth_of(rows); g_tree = sum(st.values())
        z = resolve(rows, st); d_all = depth_of(z); g_all = sum(st.values())
        print(f"  {W:2d}×{W:<2d}: tree only {g_tree:6d} gates depth {d_tree:3d} | + Kogge–Stone {g_all:6d} gates depth {d_all:3d}")

if __name__ == "__main__":
    self_test()
