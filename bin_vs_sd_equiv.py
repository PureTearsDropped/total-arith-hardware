#!/usr/bin/env python3
"""Do the binary carry-save circuits (bin2_gates.py) produce the same VALUES as the signed-digit ones (sd2_core.py,
sd2_gates.py)?  Same integers in, compare: scalar products (edge cases exhaustively + random), multiply-accumulate
chains, and the full 16-component sedenion product (signs from the OMEGA table folded with neg_cs)."""
import sys, os, random, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sed'))
from bin2_gates import to_bits, from_bits, cs_val, mul_cs, resolve, mac_cs, neg_cs, add_cs, sign_extend, is_zero_cs, sign_cs
from gate_bilinear import new_counter
from sd2_core import to_sd2, from_sd2, sd2_mul, sedenion_mult_sd2
from sedenion_tensor_logic import OMEGA, ref_mult, M
from gate_fast import multiply_fast, wrap
rng = random.Random(20260906)
K = 10                      # signed-digit width; binary uses K+1 bits (range ±1024 ⊇ ±1023)
lim = (1 << K) - 1          # 1023
def sd_mul_val(a, b):       # signed-digit golden (integer level, sd2_core)
    return from_sd2(sd2_mul(to_sd2(a, K), to_sd2(b, K)))
def sd_mul_gates(a, b):     # signed-digit gate-level golden (the code emit_sv traces for sd_mult10)
    st = new_counter(); Z = multiply_fast(wrap(to_sd2_pairs(a)), wrap(to_sd2_pairs(b)), st)
    bv = lambda b: int(b.v) if hasattr(b, 'v') else int(b)
    return sum((bv(p) - bv(n)) * (1 << i) for i, (p, n) in enumerate(Z))
def to_sd2_pairs(v):
    return [(1 if d > 0 else 0, 1 if d < 0 else 0) for d in to_sd2(v, K)]
def bin_mul_val(a, b, resolved=True):
    st = new_counter(); rows = mul_cs(to_bits(a, K + 1), to_bits(b, K + 1), st)
    return from_bits(resolve(rows, st)) if resolved else cs_val(rows, 2 * (K + 1))

edges = [0, 1, -1, 2, -2, 3, -3, 511, -511, 512, -512, 1022, -1022, lim, -lim]
cases = list(itertools.product(edges, edges)) + [(rng.randint(-lim, lim), rng.randint(-lim, lim)) for _ in range(20000)]
bad = 0
for a, b in cases:
    v_sd = sd_mul_val(a, b); v_sdg = sd_mul_gates(a, b); v_b = bin_mul_val(a, b); v_bt = bin_mul_val(a, b, False)
    if not (v_sd == v_sdg == v_b == v_bt == a * b): bad += 1
print(f"scalar product: {len(cases)} cases (225 edge combinations + 20000 random), signed-digit golden == signed-digit gates == binary resolved == binary rows == a·b: {'all equal' if bad == 0 else f'{bad} DIFFER'}")

# multiply-accumulate chains: signed-digit accumulate (sd2_add of products) vs binary carry-save chain
from sd2_core import sd2_add
bad = 0
for Kt in (2, 4, 8, 16):
    for _ in range(300):
        xs = [rng.randint(-lim, lim) for _ in range(Kt)]; ys = [rng.randint(-lim, lim) for _ in range(Kt)]
        acc_sd = to_sd2(0, 2 * K + 6)
        for x, y in zip(xs, ys): acc_sd = sd2_add(acc_sd, sd2_mul(to_sd2(x, K), to_sd2(y, K)))
        Wa = 2 * (K + 1) + 6; st = new_counter(); acc = ([0] * Wa, [0] * Wa)
        for x, y in zip(xs, ys): acc = mac_cs(to_bits(x, K + 1), to_bits(y, K + 1), acc, st)
        want = sum(x * y for x, y in zip(xs, ys))
        if not (from_sd2(acc_sd) == cs_val(acc, Wa) == want): bad += 1
print(f"multiply-accumulate chains (K = 2, 4, 8, 16 terms, 300 each): {'all equal' if bad == 0 else f'{bad} DIFFER'}")

# sedenion product: sedenion_mult_sd2 (signed digits) vs binary carry-save with OMEGA signs
def sed_mult_cs(xi, yi, Wd, st):
    """xi, yi: 16 integers; returns 16 carry-save (row0,row1) accumulators, width 2Wd+6"""
    Wa = 2 * Wd + 6; out = []
    X = [to_bits(v, Wd) for v in xi]; Y = [to_bits(v, Wd) for v in yi]
    for k in range(M):
        rows = []
        for i in range(M):
            j = i ^ k; s = OMEGA[i, j]
            p = mul_cs(X[i], Y[j], st); p = (sign_extend(p[0], Wa), sign_extend(p[1], Wa))
            if s < 0: p = neg_cs(p, Wa, st)
            rows += [p[0], p[1]]
        out.append(add_cs(rows, Wa, st))
    return out
bad = 0; N = 200; Wd = 6
for _ in range(N):
    xi = [rng.randint(-(1 << (Wd - 1)) + 1, (1 << (Wd - 1)) - 1) for _ in range(M)]; yi = [rng.randint(-(1 << (Wd - 1)) + 1, (1 << (Wd - 1)) - 1) for _ in range(M)]
    sd = [from_sd2(w) for w in sedenion_mult_sd2([to_sd2(v, Wd) for v in xi], [to_sd2(v, Wd) for v in yi])]
    st = new_counter(); cs = [cs_val(r, 2 * Wd + 6) for r in sed_mult_cs(xi, yi, Wd, st)]
    ref = ref_mult(xi, yi)
    if not (sd == cs == ref): bad += 1
print(f"sedenion product (16 components, {Wd}-bit digits, {N} random pairs): signed-digit == binary carry-save == reference: {'all equal' if bad == 0 else f'{bad} DIFFER'}")
# zero divisors: a·b = 0 with a, b ≠ 0 — the binary rows must report structural zero via is_zero_cs
zd = 0; good = 0
for a, sb, b, c, sd_, d in [(1, 1, 10, 4, -1, 15), (1, -1, 10, 4, 1, 15), (1, 1, 10, 5, 1, 14)]:
    xi = [0] * M; yi = [0] * M; xi[a] = 3; xi[b] = 3 * sb; yi[c] = 5; yi[d] = 5 * sd_
    st = new_counter(); rows = sed_mult_cs(xi, yi, Wd, st)
    zero_flags = [is_zero_cs(r, st) for r in rows]; vals = [cs_val(r, 2 * Wd + 6) for r in rows]
    zd += 1; good += (all(zero_flags) and all(v == 0 for v in vals) and ref_mult(xi, yi) == [0] * M)
print(f"zero-divisor pairs: {good}/{zd} give all-zero components with is_zero_cs = 1 on every carry-save row pair (structural zero detected without resolving)")
# sign of products: sign_cs on carry-save rows == sign of the exact product
bad = 0
for _ in range(5000):
    a = rng.randint(-lim, lim); b = rng.randint(-lim, lim); st = new_counter(); rows = mul_cs(to_bits(a, K + 1), to_bits(b, K + 1), st)
    if sign_cs(rows, st) != (1 if a * b < 0 else 0): bad += 1
    if is_zero_cs(rows, st) != (1 if a * b == 0 else 0): bad += 1
print(f"sign_cs / is_zero_cs on product rows (5000 random): {'all correct' if bad == 0 else f'{bad} wrong'}")
