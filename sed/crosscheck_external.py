#!/usr/bin/env python3
"""Independent cross-check of the Cayley–Dickson tables (quaternion / octonion / sedenion) against things this
repository did not write:

  * an exact rational Cayley–Dickson product written from the textbook doubling
        (a,b)(c,d) = (a c − d̄ b,  d a + b c̄),   conj(a,b) = (ā, −b)
  * numpy-quaternion            (pip install numpy-quaternion)            — optional
  * Quaternions.jl / Octonions.jl via crosscheck_ref.jl                 — optional (needs julia + packages)
  * convention-free invariants: associativity (ℍ), alternativity / Moufang / norm multiplicativity (𝕆),
    flexibility / power-associativity / zero divisors / norm failure (sedenions)

Basis convention of this repository:  e_i · e_j = OMEGA[i,j] · e_{i XOR j}   (Cayley–Dickson "XOR" labelling).
Octonions.jl uses a different labelling of the seven imaginary units; the two tables are related by a signed
permutation (1344 of them exist — the automorphisms of the Fano plane times signs), e.g.
    e_i -> s_i f_{p(i)},  p = (0,1,2,3,4,7,6,5),  s = (+,−,+,−,−,−,−,−)
so element-wise numbers are NOT interchangeable between the two without this change of basis.

Run:  python3 sed/crosscheck_external.py            (exact reference + invariants always; libraries if present)
      JULIA=/path/to/julia python3 sed/crosscheck_external.py   to include the Julia references
"""
import os, sys, random, subprocess, itertools
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sedenion_tensor_logic import OMEGA, ref_mult
rng = random.Random(2026)

def theirs(x, y, d=16):
    xx = list(x) + [0] * (16 - len(x)); yy = list(y) + [0] * (16 - len(y)); return ref_mult(xx, yy)[:d]
def conj(x):
    n = len(x); h = n // 2
    return [x[0]] if n == 1 else conj(x[:h]) + [-v for v in x[h:]]
def cd(x, y):
    n = len(x)
    if n == 1: return [x[0] * y[0]]
    h = n // 2; a, b, c, d = x[:h], x[h:], y[:h], y[h:]
    ac, db, da, bc = cd(a, c), cd(conj(d), b), cd(d, a), cd(b, conj(c))
    return [p - q for p, q in zip(ac, db)] + [p + q for p, q in zip(da, bc)]
def rvec(d): return [Fr(rng.randint(-5, 5)) for _ in range(d)]
def T(x, y, d): return [Fr(v) for v in theirs(x, y, d)]
def norm2(x): return sum(v * v for v in x)
fails = 0
def report(name, ok, detail=""):
    global fails
    print(("  PASS  " if ok else "  FAIL  ") + name + (f"  ({detail})" if detail else "")); fails += (not ok)

print("exact from-definition Cayley–Dickson vs OMEGA / XOR table")
bad = 0
for i in range(16):
    for j in range(16):
        e = [Fr(0)] * 16; f = [Fr(0)] * 16; e[i] = Fr(1); f[j] = Fr(1); v = cd(e, f)
        k = [t for t in range(16) if v[t] != 0]
        if len(k) != 1 or k[0] != (i ^ j) or v[k[0]] != OMEGA[i, j]: bad += 1
report("basis table (256 products)", bad == 0, f"{256-bad}/256")
bad = sum(1 for _ in range(1000) if (lambda x, y: cd(x, y) != T(x, y, 16))(rvec(16), rvec(16)))
report("random sedenion products, exact", bad == 0, f"{1000-bad}/1000")

print("invariants (exact rationals)")
report("quaternions associative", all(T(T(x, y, 4), z, 4) == T(x, T(y, z, 4), 4) for x, y, z in (tuple(rvec(4) for _ in range(3)) for _ in range(200))))
report("octonions alternative", all(T(x, T(x, y, 8), 8) == T(T(x, x, 8), y, 8) and T(T(y, x, 8), x, 8) == T(y, T(x, x, 8), 8) for x, y in ((rvec(8), rvec(8)) for _ in range(200))))
report("octonions Moufang", all(T(T(T(z, x, 8), z, 8), y, 8) == T(z, T(x, T(z, y, 8), 8), 8) for x, y, z in (tuple(rvec(8) for _ in range(3)) for _ in range(200))))
report("octonions norm multiplicative", all(norm2(T(x, y, 8)) == norm2(x) * norm2(y) for x, y in ((rvec(8), rvec(8)) for _ in range(200))))
report("octonions not associative", any(T(T(x, y, 8), z, 8) != T(x, T(y, z, 8), 8) for x, y, z in (tuple(rvec(8) for _ in range(3)) for _ in range(100))))
report("sedenions flexible", all(T(x, T(y, x, 16), 16) == T(T(x, y, 16), x, 16) for x, y in ((rvec(16), rvec(16)) for _ in range(200))))
report("sedenions power-associative", all(T(T(x, x, 16), x, 16) == T(x, T(x, x, 16), 16) and T(T(T(x, x, 16), x, 16), x, 16) == T(T(x, x, 16), T(x, x, 16), 16) for x in (rvec(16) for _ in range(200))))
report("sedenions norm NOT multiplicative", any(norm2(T(x, y, 16)) != norm2(x) * norm2(y) for x, y in ((rvec(16), rvec(16)) for _ in range(100))))
zd = 0
for a in range(1, 16):
    for b in range(a + 1, 16):
        for c in range(1, 16):
            for d in range(c + 1, 16):
                for sb in (1, -1):
                    for sd in (1, -1):
                        x = [Fr(0)] * 16; y = [Fr(0)] * 16; x[a] = Fr(1); x[b] = Fr(sb); y[c] = Fr(1); y[d] = Fr(sd)
                        if all(v == 0 for v in T(x, y, 16)): zd += 1
report("sedenion zero divisors (e_a±e_b)(e_c±e_d)=0: 336 = 84×4", zd == 336, str(zd))

try:
    import numpy as np, quaternion
    bad = 0
    for _ in range(500):
        x = [rng.uniform(-3, 3) for _ in range(4)]; y = [rng.uniform(-3, 3) for _ in range(4)]
        p = np.quaternion(*x) * np.quaternion(*y)
        if max(abs(a - b) for a, b in zip([p.w, p.x, p.y, p.z], theirs(x, y, 4))) > 1e-9: bad += 1
    report("numpy-quaternion agrees (e1=i,e2=j,e3=k)", bad == 0, f"{500-bad}/500")
except ImportError: print("  skip  numpy-quaternion not installed")

JULIA = os.environ.get("JULIA")
if JULIA:
    ref = os.path.join(os.path.dirname(os.path.abspath(__file__)), "crosscheck_ref.jl")
    def jl(lines):
        import tempfile; fd, p = tempfile.mkstemp(suffix=".txt"); os.write(fd, ("\n".join(lines) + "\n").encode()); os.close(fd)
        return subprocess.run([JULIA, ref, p], capture_output=True, text=True).stdout.strip().splitlines()
    qs = [([rng.uniform(-3, 3) for _ in range(4)], [rng.uniform(-3, 3) for _ in range(4)]) for _ in range(300)]
    out = jl([" ".join(map(repr, x + y)) for x, y in qs])
    bad = sum(1 for (x, y), l in zip(qs, out) if max(abs(a - b) for a, b in zip(map(float, l.split()), theirs(x, y, 4))) > 1e-9)
    report("Quaternions.jl agrees", bad == 0 and len(out) == 300, f"{len(out)-bad}/{len(out)}")
    lines = []
    for i in range(8):
        for j in range(8):
            e = [0.0] * 8; f = [0.0] * 8; e[i] = 1; f[j] = 1; lines.append(" ".join(map(repr, e + f)))
    out = jl(lines); import numpy as np
    K = [[0] * 8 for _ in range(8)]; S = [[0] * 8 for _ in range(8)]
    for n, l in enumerate(out):
        i, j = divmod(n, 8); v = list(map(float, l.split())); k = max(range(8), key=lambda t: abs(v[t])); K[i][j] = k; S[i][j] = 1 if v[k] > 0 else -1
    iso = 0; example = None
    for perm in itertools.permutations(range(1, 8)):
        p = [0] + list(perm)
        if any(K[p[i]][p[j]] != p[i ^ j] for i in range(1, 8) for j in range(1, 8) if i != j): continue
        for bits in range(128):
            s = [1] + [1 if (bits >> t) & 1 else -1 for t in range(7)]
            if all(s[i] * s[j] * S[p[i]][p[j]] == OMEGA[i, j] * s[i ^ j] for i in range(8) for j in range(8)):
                iso += 1; example = example or (p, s)
    report("Octonions.jl isomorphic via signed permutation (1344 expected)", iso == 1344, f"{iso}, e.g. p={example[0]} s={example[1]}" if example else str(iso))
else: print("  skip  set JULIA=/path/to/julia (with Quaternions.jl, Octonions.jl) for the Julia references")
print("FAILURES:", fails); sys.exit(1 if fails else 0)
