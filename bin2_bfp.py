#!/usr/bin/env python3
"""Binary block normaliser — the two's-complement / carry-save counterpart of gate_fast.block_normalize_g_fast,
same contract, same flags.

INPUT   M carry-save numbers (r0, r1) of Wc bits (two's complement, value (r0+r1) mod 2^Wc), one shared exponent bus.
OUTPUT  M mantissas as (sign, magnitude[W]) — the binary analogue of W canonical signed digits — the final exponent
        bus, and per component the flags (ge, le, 0) with exactly the signed-digit rules:
          truncation of nonzero low bits with something kept → ge ; overflow (exponent capped, bits above W) → ±MAX, ge ;
          collapse (everything nonzero dropped, nothing kept) → ±MIN, le ; true zero stays 0 with no flag.
STEPS   resolve D = r0+r1 and, in parallel, D2 = −(r0+r1) (rows complemented + 2, then Kogge–Stone), sign = MSB of D,
        magnitude = sign ? D2 : D   (this is canonicalize_fast's bidirectional trick on rows)
        → leading-one position per component (suffix-OR doubling tree) → Lmax tournament → shift / exponent / saturation
        exactly as block_normalize_g_fast → barrel shift of the magnitude collecting dropped-nonzero → select
        MAX / MIN / kept → flags.
Magnitude truncation (not two's-complement arithmetic shift) keeps "dropped bits ⇒ |true| ≥ |shown|" i.e. the ≥ flag.
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate_bilinear import AND, OR, NOT, XOR, new_counter
from gate_fast import bin_add_fast, or_tree, bus_max_fast, bus_sub_fast, bus_add_fast, bus_lt_fast
from gate_exponent import bus_const, mux_bus, clamp0, mux_bit
from bin2_gates import neg_cs, to_bits, from_bits, cs_val

def _ext(row, W):                      # sign-extend a two's complement row to W bits
    return list(row) + [row[-1]] * (W - len(row))

def sign_magnitude_cs(rows, st):
    """(sign, magnitude bits) of a carry-save two's complement number; the magnitude is read as an unsigned Wc-bit
    number, so even −2^(Wc−1) fits (|−2^(Wc−1)| = 2^(Wc−1) < 2^Wc)."""
    r0, r1 = rows; Wc = len(r0)
    e0, e1 = list(r0), list(r1)                                # rows are arbitrary bit patterns: no per-row sign extension
    D = bin_add_fast(e0, e1, st)                              # (r0 + r1) mod 2^Wc = the value in two's complement
    n0, n1 = neg_cs((e0, e1), Wc, st)                          # −(r0 + r1) as two rows (parallel); read unsigned = |value|
    D2 = bin_add_fast(list(n0), list(n1), st)
    sign = D[-1]
    mag = mux_bus(sign, D2, D, st)                             # mux_bus(s, A, B): s ? A : B
    return sign, mag

def leading_one_bits(bits, EW, st):
    """position L of the highest 1 (EW-bit bus), and `none` (all zero) — priority_encoder_fast on plain bits"""
    n = len(bits); suf = list(bits); k = 1
    while k < n:
        suf = [OR(suf[i], suf[i + k], st) if i + k < n else suf[i] for i in range(n)]; k <<= 1
    none = NOT(suf[0], st)
    onehot = [AND(bits[i], NOT(suf[i + 1], st), st) if i + 1 < n else bits[i] for i in range(n)]
    L = [or_tree([onehot[i] for i in range(n) if (i >> b) & 1], st) for b in range(EW)]
    return L, none

def barrel_right_bits(bits, S, st):
    """shift right by bus S (low bits fall off); returns (shifted[fixed width], dropped_nonzero)"""
    n = len(bits); cur = list(bits); dropped = 0
    for j, sbit in enumerate(S):
        k = 1 << j; dnz = 0
        for i in range(min(k, n)): dnz = OR(dnz, cur[i], st)
        dropped = OR(dropped, AND(sbit, dnz, st), st)
        cur = [mux_bit(sbit, cur[i + k] if i + k < n else 0, cur[i], st) for i in range(n)]
    return cur, dropped

def block_normalize_bin(rows_list, Ebus, W, Emax, st):
    """rows_list: M carry-save numbers; returns ([(sign, mag[W])]*M, E_fin, [(ge, le, 0)]*M)"""
    EW = len(Ebus); signs = []; mags = []; Ls = []
    for rows in rows_list:
        s, m = sign_magnitude_cs(rows, st); signs.append(s); mags.append(m)
        L, _ = leading_one_bits(m, EW, st); Ls.append(L)
    Wc = len(mags[0])
    while len(Ls) > 1:                                           # tournament
        Ls = [bus_max_fast(Ls[i], Ls[i + 1], st) if i + 1 < len(Ls) else Ls[i] for i in range(0, len(Ls), 2)]
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
    for s, m in zip(signs, mags):
        shifted, drop_nz = barrel_right_bits(m, sh_al[:SW], st)
        kept = shifted[:W]
        kept_nz = or_tree(list(kept), st)
        over = or_tree(list(shifted[W:]), st)
        collapse = AND(drop_nz, NOT(kept_nz, st), st)
        sel_max = over
        sel_min = AND(collapse, NOT(over, st), st)
        sel_kept = AND(NOT(over, st), NOT(collapse, st), st)
        om = []
        for i in range(W):
            t = AND(sel_kept, kept[i], st)
            if i == 0: t = OR(t, sel_min, st)                 # ±MIN: magnitude 1
            om.append(OR(t, sel_max, st))                     # ±MAX: all ones
        ge = OR(over, AND(drop_nz, kept_nz, st), st)
        le = collapse
        out.append((s, om)); flags.append((ge, le, 0))
    return out, E_fin, flags

# ------------------------------------------------------------------ golden helpers
def sm_val(sm):
    s, m = sm; v = sum(int(b.v if hasattr(b, 'v') else b) << i for i, b in enumerate(m)); s = int(s.v if hasattr(s, 'v') else s)
    return -v if s else v

def self_test(seed=20260906):
    """same values into the signed-digit normaliser and the binary one: mantissas, exponent, flags must agree"""
    from gate_fast import block_normalize_g_fast
    from gate_bilinear import to_sd, from_sd
    from gate_exponent import bus_val
    rng = random.Random(seed); M, Win, W, Emax, EW = 4, 24, 6, 20, 12; Wc = Win + 2
    cases = 0; bad = 0
    for _ in range(600):
        E0 = rng.randint(0, 24)
        vals = [rng.choice([0, rng.randint(-5, 5) * (10 ** rng.randint(0, 6)), rng.randint(-(1 << Win) + 1, (1 << Win) - 1), rng.randint(-40, 40)]) for _ in range(M)]
        # binary inputs as random carry-save splits of each value (mod 2^Wc)
        rows = []
        for v in vals:
            a = rng.getrandbits(Wc); b = (v - a) % (1 << Wc)
            rows.append(([(a >> i) & 1 for i in range(Wc)], [(b >> i) & 1 for i in range(Wc)]))
        st = new_counter()
        ob, Eb, fb = block_normalize_bin(rows, bus_const(E0, EW), W, Emax, st)
        og, Eg, fg = block_normalize_g_fast([to_sd(v, Win) for v in vals], bus_const(E0, EW), W, Emax, new_counter())
        cases += 1
        if bus_val(Eb) % (1 << EW) != bus_val(Eg) % (1 << EW): bad += 1; continue
        for i in range(M):
            if sm_val(ob[i]) != from_sd(og[i]) or (int(fb[i][0]), int(fb[i][1])) != (int(fg[i][0]), int(fg[i][1])):
                bad += 1
                if bad <= 4: print(f"  diff: v={vals[i]} E0={E0} bin=({sm_val(ob[i])}, ge={int(fb[i][0])}, le={int(fb[i][1])}) sd=({from_sd(og[i])}, ge={int(fg[i][0])}, le={int(fg[i][1])}) E bin={bus_val(Eb)} sd={bus_val(Eg)}")
                break
    print(f"bin2_bfp vs block_normalize_g_fast: {cases} blocks × {M} components (values, exponent, ge/le flags): {'all equal' if bad == 0 else f'{bad} DIFFER'}")
    return bad == 0

if __name__ == "__main__":
    ok = self_test(); sys.exit(0 if ok else 1)
