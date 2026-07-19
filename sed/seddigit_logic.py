#!/usr/bin/env python3
"""Digit-major sedenion arithmetic: per-trit sedenion-digit cells with
explicit carry/borrow wiring.

VIEW CHANGE (requested reorganization):
  component-major (sedenion_tensor_logic): a number is 16 trit WORDS.
  digit-major (this module):               a number is a word of SEDENION
                                           DIGITS — each digit is a 16-trit
                                           vector (one trit per component).
The two are transposes of the same trit array; arithmetic must agree
exactly (verified in self_test against the component-major module).

CELLS AND THEIR CARRY WIRING (all fixed wiring, data-independent):
  sed_digit_add : 16 gate27 in parallel.
                  carry wiring: (component c, weight w) -> (c, w+1).
                  In balanced ternary the high trit is SIGNED, so carry
                  (+1) and borrow (-1) are the SAME wire pair — 桁上げ and
                  桁下げ are unified by the (p,n) encoding.
  unit action   : per-digit wiring only (component permutation + sign).
                  The algebra acts identically on every digit slice; it
                  never crosses weights. Only the NUMBER SYSTEM (carries)
                  connects slices.
  sed_digit_mult_cell : one cell per weight pair (p,q).
                  256 gate9 products, XOR routing, omega sign wiring,
                  per-component compressor tree + ripple INSIDE the cell.
                  carry wiring: cell (p,q) emits digits to global weights
                  p+q, p+q+1, ..., p+q+R-1. A single digit-product sums up
                  to 16 trits per component; the VALUE fits in 4 trits
                  (|sum| <= 16 <= 40), and the ripple's final carry slot
                  adds one structurally reserved wire, so R = 5 with the
                  top wire frequently constant zero (synthesizable away).
"""
import numpy as np
from sedenion_tensor_logic import (gate9, gate27, OMEGA, to_trits, from_trits,
                                   static_tree, two_rows, ripple,
                                   sedenion_mult, ref_mult, tensor_units,
                                   ref_tensor_units)

M = 16

# ------------------------------------------------ view transposition (framing)
def word_to_digits(xw):
    """component-major trit words -> weight-major sedenion digits."""
    K = max(len(w) for w in xw)
    return [[(xw[c][k] if k < len(xw[c]) else 0) for c in range(M)]
            for k in range(K)]

def digits_to_ints(digits):
    """decode: per component integer (host-side reading)."""
    return [from_trits([d[c] for d in digits]) for c in range(M)]

# ------------------------------------------------------- sedenion digit adder
def sed_digit_add(a, b, cin, stats):
    """digit + digit + carry-digit -> (sum digit, carry digit).
    16 parallel gate27; carry wiring: (c, w) -> (c, w+1), signed
    (carry = +1, borrow = -1 on the same wires)."""
    s = [0] * M
    cout = [0] * M
    for c in range(M):
        s[c], cout[c] = gate27(a[c], b[c], cin[c])
        stats['g27'] += 1
    return s, cout

def sed_word_add(A, B, stats):
    """add two digit-major words with a per-component carry chain."""
    W = max(len(A), len(B))
    A = A + [[0]*M for _ in range(W - len(A))]
    B = B + [[0]*M for _ in range(W - len(B))]
    out, carry = [], [0]*M
    for w in range(W):
        s, carry = sed_digit_add(A[w], B[w], carry, stats)
        out.append(s)
    out.append(carry)
    return out

# ---------------------------------------------------- unit action (per digit)
def sed_digit_unit(j, sign, d):
    """(sign e_j) acting on ONE digit: wiring only, no weight crossing."""
    return [sign * OMEGA[j, j ^ k] * d[j ^ k] for k in range(M)]

def sed_word_unit(j, sign, digits):
    return [sed_digit_unit(j, sign, d) for d in digits]

# --------------------------------------------------- sedenion digit mult cell
def sed_digit_mult_cell(xd, yd, stats):
    """digit x digit -> list of digits at relative weights 0..R-1.
    Internal: 256 gate9, XOR routing + omega sign wiring, per-component
    static tree + ripple. Static shape (zeros flow through gates)."""
    piles = [[[]] for _ in range(M)]          # per component, rel weight 0
    for i in range(M):
        for j in range(M):
            t = gate9(xd[i], yd[j])
            stats['g9'] += 1
            piles[i ^ j][0].append(OMEGA[i, j] * t)   # sign = wire swap
    comp_words = []
    for c in range(M):
        cols = static_tree(piles[c], stats)
        word = ripple(*two_rows(cols), stats)
        comp_words.append(word)
    R = max(len(w) for w in comp_words)
    return [[(comp_words[c][r] if r < len(comp_words[c]) else 0)
             for c in range(M)] for r in range(R)]

# ------------------------------------------------- full digit-major multiplier
def sed_mult_digits(X, Y, stats):
    """X, Y: digit-major words. One mult cell per weight pair (p,q);
    cell (p,q) wires its output digits to global weights p+q+r."""
    piles = [[[] for _ in range(len(X) + len(Y) + 4)] for _ in range(M)]
    for p, xd in enumerate(X):
        for q, yd in enumerate(Y):
            for r, dig in enumerate(sed_digit_mult_cell(xd, yd, stats)):
                for c in range(M):
                    piles[c][p + q + r].append(dig[c])   # fixed carry wiring
    out_words = []
    for c in range(M):
        cols = static_tree(piles[c], stats)
        out_words.append(ripple(*two_rows(cols), stats))
    W = max(len(w) for w in out_words)
    return [[(out_words[c][w] if w < len(out_words[c]) else 0)
             for c in range(M)] for w in range(W)]

# --------------------------------- tensor member A, digit-sliced organization
def tensor_units_digits(Gidx, Xs, m, stats):
    """member A (staged, lsb-first) on digit-major data.
    The unit wiring is applied per digit slice (weight-independent);
    carries arise ONLY from the additions."""
    Y = [ [d[:] for d in X] for X in Xs ]
    for s in range(m):
        out = [None] * len(Y)
        for base in range(len(Y)):
            if (base >> s) & 1: continue
            u, v = Y[base], Y[base | (1 << s)]
            (j00, s00), (j01, s01) = Gidx[0]
            (j10, s10), (j11, s11) = Gidx[1]
            out[base]            = sed_word_add(sed_word_unit(j00, s00, u),
                                                sed_word_unit(j01, s01, v), stats)
            out[base | (1 << s)] = sed_word_add(sed_word_unit(j10, s10, u),
                                                sed_word_unit(j11, s11, v), stats)
        Y = out
    return Y

# ------------------------------------------------------------------ self-test
def self_test(seed=20260711):
    rng = np.random.default_rng(seed)
    K = 6

    st = {'g9': 0, 'g27': 0}
    for _ in range(8):
        x = [int(v) for v in rng.integers(-3**K // 2, 3**K // 2, M)]
        y = [int(v) for v in rng.integers(-3**K // 2, 3**K // 2, M)]
        X = word_to_digits([to_trits(c, K) for c in x])
        Y = word_to_digits([to_trits(c, K) for c in y])
        got = digits_to_ints(sed_mult_digits(X, Y, st))
        assert got == ref_mult(x, y)
    print("digit-major multiplier == reference: 8 random tests exact")

    # agreement with the component-major implementation (same arithmetic)
    x = [int(v) for v in rng.integers(-3**K // 2, 3**K // 2, M)]
    y = [int(v) for v in rng.integers(-3**K // 2, 3**K // 2, M)]
    a = digits_to_ints(sed_mult_digits(
        word_to_digits([to_trits(c, K) for c in x]),
        word_to_digits([to_trits(c, K) for c in y]), {'g9': 0, 'g27': 0}))
    b = [from_trits(w) for w in sedenion_mult(
        [to_trits(c, K) for c in x], [to_trits(c, K) for c in y],
        {'g9': 0, 'g27': 0})]
    assert a == b
    print("digit-major == component-major (transposed views agree)")

    # static shape: gate counts data-independent
    sA, sB = {'g9': 0, 'g27': 0}, {'g9': 0, 'g27': 0}
    Z = word_to_digits([to_trits(0, K)] * M)
    Rn = word_to_digits([to_trits(int(v), K)
                         for v in rng.integers(-3**K//2, 3**K//2, M)])
    sed_mult_digits(Z, Z, sA); sed_mult_digits(Rn, Rn, sB)
    assert sA == sB
    print(f"static structure OK (per K={K} product: gate9 {sB['g9']}, "
          f"gate27 {sB['g27']})")

    # cell carry span: how far one digit-product cell reaches
    st2 = {'g9': 0, 'g27': 0}
    ones = [1] * M
    span = len(sed_digit_mult_cell(ones, ones, st2))
    print(f"mult cell carry wiring span: rel weights 0..{span-1} "
          f"(per-component |sum| <= 16 -> {span} trits)")

    # tensor member A, digit-sliced == component-major tensor
    Gidx = [[(0, 1), (1, 1)], [(2, 1), (4, -1)]]
    m = 3
    xs = [[int(v) for v in rng.integers(-3**4, 3**4, M)] for _ in range(2**m)]
    Xd = [word_to_digits([to_trits(c, 10) for c in q]) for q in xs]
    st3 = {'g9': 0, 'g27': 0}
    A_dig = [digits_to_ints(q) for q in tensor_units_digits(Gidx, Xd, m, st3)]
    A_ref = ref_tensor_units(Gidx, xs, m)
    assert A_dig == A_ref
    print(f"digit-sliced tensor member A == reference (m={m}) | "
          f"gate9 {st3['g9']} (unit wiring: ZERO multipliers), "
          f"gate27 {st3['g27']} (adders only; carries are the only "
          f"inter-slice wires)")
    print("all self-tests passed.")

if __name__ == "__main__":
    self_test()
