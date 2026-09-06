#!/usr/bin/env python3
"""Fused sedenion component in binary carry-save form — the counterpart of mul_fused.group_component.
component k = Σᵢ σ(i, i⊕k) · aᵢ · b_{i⊕k}: 16 Baugh–Wooley partial-product sets, the σ = −1 ones negated as rows
(complement + 2: inverters and one constant, no adder), all reduced by one Dadda tree to two rows.
`sed_comp_cs` returns the two rows; `sed_comp_bin` resolves them with one Kogge–Stone."""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate_bilinear import AND, OR, NOT, XOR, new_counter, full_adder
from bin2_gates import mul_cs, neg_cs, add_cs, resolve, sign_extend, to_bits, from_bits, cs_val, half_adder

def sed_comp_cs(a, b, OM, M, k, Wa, st):
    """a, b: M two's complement digit-words (lists of bits); OM: sign table.  One Dadda tree over all
    16 Baugh–Wooley partial-product sets.  A σ = −1 product is negated *as partial products*:
    −Σ 2^w t_w = Σ 2^w (1−t_w) − Σ 2^w, so its bits are complemented (folded into the Baugh–Wooley
    inversion parity: no extra gate) and one constant per product is accumulated; all constants
    fold into a single constant word mod 2^Wa.  Returns (row0, row1) of Wa bits."""
    from bin2_gates import dadda_rows
    cols = [[] for _ in range(Wa)]; const = 0
    for i in range(M):
        j = i ^ k; X, Y = a[i], b[j]; Wx, Wy = len(X), len(Y); W = Wx + Wy
        neg = OM[i][j] < 0
        for p, x in enumerate(X):
            for q, y in enumerate(Y):
                if p + q >= Wa: continue
                t = AND(x, y, st)
                inv = ((p == Wx - 1) != (q == Wy - 1)) != neg
                if inv: t = NOT(t, st)
                cols[p + q].append(t)
        # exact integer identity −t·2^w = (1−t)·2^w − 2^w on the two cross rows: constant of +P as an
        # integer (not reduced mod 2^W), so no sign extension is needed inside the wider accumulator
        c = -((1 << (W - 1)) - (1 << (Wx - 1)) - (1 << (Wy - 1)))
        if neg: c = -c - sum((1 << (p + q)) for p in range(Wx) for q in range(Wy))
        const += c
    const &= (1 << Wa) - 1
    for w in range(Wa):
        if (const >> w) & 1: cols[w].append(1)
    return dadda_rows(cols, Wa, st)

def sed_comp_bin(a, b, OM, M, k, Wa, st):
    return resolve(sed_comp_cs(a, b, OM, M, k, Wa, st), st)

def self_test(seed=20260906, K=6, M=16, k=1, N=200):
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sed'))
    from sedenion_tensor_logic import OMEGA, ref_mult
    OM = [[int(OMEGA[i, j]) for j in range(M)] for i in range(M)]
    rng = random.Random(seed); lim = 1 << (K - 1); Wa = 2 * K + 5; bad = 0
    for _ in range(N):
        xi = [rng.randrange(-lim, lim) for _ in range(M)]; yi = [rng.randrange(-lim, lim) for _ in range(M)]
        st = new_counter()
        rows = sed_comp_cs([to_bits(v, K) for v in xi], [to_bits(v, K) for v in yi], OM, M, k, Wa, st)
        z = resolve(rows, st)
        want = ref_mult(xi, yi)[k]
        if cs_val(rows, Wa) != want or from_bits(z) != want: bad += 1
    print(f"bin2_sed component {k}: {N} random sedenion pairs, carry-save rows and resolved value == ref_mult: {'all equal' if bad == 0 else f'{bad} DIFFER'}")
    st = new_counter(); from gate_fast import B, depth_of
    a = [[B(0) for _ in range(K)] for _ in range(M)]; b = [[B(0) for _ in range(K)] for _ in range(M)]
    rows = sed_comp_cs(a, b, OM, M, k, Wa, st); g1 = sum(st.values()); d1 = depth_of(rows)
    z = resolve(rows, st); print(f"  gates: rows {g1} (depth {d1}) | resolved {sum(st.values())} (depth {depth_of(z)})")
    return bad == 0

if __name__ == "__main__":
    sys.exit(0 if self_test() else 1)
