#!/usr/bin/env python3
"""Carry-save arithmetic at SEDENION-DIGIT granularity.

NUMBER FORMAT (redundant, as requested):
  at every weight position w, hold up to TWO sedenion digits:
      S[w]  = the no-carry answer digit   (16 trits)
      C[w]  = the carry digit             (16 trits)
  i.e. value = sum_w 3^w * (S[w] + C[w]), componentwise.
  Carries (+1) and borrows (-1) share the same wires (signed trits).

GATES — everything is a gate between 16-trit vectors:
  sed_gate27(a, b, c) -> (same, high):
      three sedenion digits in, one same-weight digit + one carry digit out.
      Internally 16 parallel gate27; NO wiring between components.
  sed_partial(xd, yd, j) -> one sedenion digit:
      the j-th partial product of two digits:
      out[k] = omega(k XOR j, j) * gate9(xd[k XOR j], yd[j]).
      16 gate9 + fixed wiring. A digit product is EXACTLY 16 such digits
      summed at one weight — so multiplication never leaves digit-vector
      land.

WHY THIS FORM: accumulating another product into the redundant number is
O(1) depth (a short sed_gate27 tree) — no carry ripple until ONE final
normalization at the very end. This is the multiply-accumulate shape that
makes low-clock-count matrix products possible, expressed at sedenion
granularity.
"""
import numpy as np
from sedenion_tensor_logic import (gate9, gate27, OMEGA,
                                   to_trits, from_trits, ref_mult)

M = 16
ZD = [0] * M                      # zero digit

# ------------------------------------------------------------- digit gates
def sed_gate27(a, b, c, stats):
    """3 sedenion digits -> (same-weight digit, carry digit)."""
    same = [0] * M
    high = [0] * M
    for k in range(M):
        same[k], high[k] = gate27(a[k], b[k], c[k])
        stats['g27'] += 1
    return same, high

def sed_partial(xd, yd, j, stats):
    """j-th partial-product digit of two sedenion digits (16 gate9 + wiring)."""
    out = [0] * M
    for k in range(M):
        t = gate9(xd[k ^ j], yd[j])
        stats['g9'] += 1
        out[k] = OMEGA[k ^ j, j] * t
    return out

# --------------------------------------------- redundant (carry-save) numbers
def cs_zero(width):
    return [[] for _ in range(width)]

def cs_inject(acc, w, digit):
    while len(acc) <= w:
        acc.append([])
    acc[w].append(digit)

def cs_compress(acc, stats):
    """reduce every weight position to at most 2 digits (static schedule)."""
    w = 0
    while w < len(acc):
        while len(acc[w]) > 2:
            a = acc[w].pop(); b = acc[w].pop(); c = acc[w].pop()
            same, high = sed_gate27(a, b, c, stats)
            acc[w].append(same)
            if w + 1 >= len(acc):
                acc.append([])
            acc[w + 1].append(high)
        w += 1
    return acc

def cs_finalize(acc, stats):
    """one final normalization: (S, C) rows -> single digit per weight.
    sed_gate27 chain with a carry digit; the ONLY ripple in the design."""
    out = []
    carry = ZD[:]
    for w in range(len(acc)):
        a = acc[w][0] if len(acc[w]) > 0 else ZD[:]
        b = acc[w][1] if len(acc[w]) > 1 else ZD[:]
        s, carry = sed_gate27(a, b, carry, stats)
        out.append(s)
    out.append(carry)
    return out

def digits_to_ints(digits):
    return [from_trits([d[k] for d in digits]) for k in range(M)]

# --------------------------------------------------- multiply / accumulate
def cs_mac(acc, X, Y, stats):
    """accumulate X*Y into the redundant accumulator.
    X, Y digit-major words. Per digit pair: 16 partial digits injected,
    then local compression — constant depth, NO ripple."""
    for p, xd in enumerate(X):
        for q, yd in enumerate(Y):
            for j in range(M):
                cs_inject(acc, p + q, sed_partial(xd, yd, j, stats))
    return cs_compress(acc, stats)

def word_to_digits(xw):
    K = max(len(w) for w in xw)
    return [[(xw[c][k] if k < len(xw[c]) else 0) for c in range(M)]
            for k in range(K)]

# ------------------------------------------------------------------ self-test
def self_test(seed=20260711):
    rng = np.random.default_rng(seed)
    K = 6

    # single product: CS multiply == reference
    st = {'g9': 0, 'g27': 0}
    for _ in range(6):
        x = [int(v) for v in rng.integers(-3**K // 2, 3**K // 2, M)]
        y = [int(v) for v in rng.integers(-3**K // 2, 3**K // 2, M)]
        acc = cs_mac(cs_zero(2 * K), word_to_digits([to_trits(c, K) for c in x]),
                     word_to_digits([to_trits(c, K) for c in y]), st)
        got = digits_to_ints(cs_finalize(acc, st))
        assert got == ref_mult(x, y)
    print("carry-save sedenion multiply == reference: 6 random tests exact")

    # MAC chain: sum of 8 products, ONE final normalization at the end
    st2 = {'g9': 0, 'g27': 0}
    acc = cs_zero(2 * K)
    ref = [0] * M
    for _ in range(8):
        x = [int(v) for v in rng.integers(-3**K // 2, 3**K // 2, M)]
        y = [int(v) for v in rng.integers(-3**K // 2, 3**K // 2, M)]
        acc = cs_mac(acc, word_to_digits([to_trits(c, K) for c in x]),
                     word_to_digits([to_trits(c, K) for c in y]), st2)
        r = ref_mult(x, y)
        ref = [a + b for a, b in zip(ref, r)]
        # accumulator stays flat: at most 2 digits per weight, no ripple yet
        assert all(len(pile) <= 2 for pile in acc)
    got = digits_to_ints(cs_finalize(acc, st2))
    assert got == ref
    print("8-term MAC in redundant form, single final normalize: exact "
          "(accumulator held <= 2 digits per weight throughout)")

    # static shape
    sA, sB = {'g9': 0, 'g27': 0}, {'g9': 0, 'g27': 0}
    Z = word_to_digits([to_trits(0, K)] * M)
    Rn = word_to_digits([to_trits(int(v), K)
                         for v in rng.integers(-3**K // 2, 3**K // 2, M)])
    cs_finalize(cs_mac(cs_zero(2 * K), Z, Z, sA), sA)
    cs_finalize(cs_mac(cs_zero(2 * K), Rn, Rn, sB), sB)
    assert sA == sB
    print(f"static structure OK (K={K} product+finalize: gate9 {sB['g9']}, "
          f"gate27 {sB['g27']})")
    print("all self-tests passed — every operation is a gate between "
          "16-trit sedenion digits; carries are digit vectors on fixed "
          "(w -> w+1) wiring.")

if __name__ == "__main__":
    self_test()
