#!/usr/bin/env python3
# ⚠️ 生成AI使用・要検証
"""配線パターンのレジストリ — 全部を 一つに束ね「名前で 配線を 選んで 実行」。

  すべての 双線形パターンを (U,V,W) に 正規化（群/畳み込み/行列/Strassen/Karatsuba/Gauss）。
  多層（算術回路）は callable で 登録（Newton 逆数）。統一インターフェース:
     run(name, a, b, backend='gate'|'algebra')   ／   catalog()   ／   REGISTRY

  各パターン ＝ 一つの 配線。ファブリック（18ゲット圧縮器＋accumulate）は 固定、配線だけ 差し替え。
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from fractions import Fraction as Fr
from nd_algebra import cd_omega, ref_mult_M
from matrix_algebra import symplectic_omega
from gate_bilinear import new_counter, bilinear_unit_gates, from_sd
from newton_recip import newton_recip_int

REGISTRY = {}

def _reg_bilinear(name, U, V, W, m, n, ref, computes, R):
    REGISTRY[name] = dict(kind='bilinear', U=U, V=V, W=W, m=m, n=n, ref=ref, computes=computes, R=R)

def _reg_circuit(name, run, computes, demo):
    REGISTRY[name] = dict(kind='circuit', run=run, computes=computes, R=None, demo=demo)

# ---- 配線 → (U,V,W): 捻れ群代数（経路 mul・符号 sig）を rank=M² 分解に ----
def wiring_to_uvw(mul, sig, M):
    R = M * M
    U = [[0]*M for _ in range(R)]; V = [[0]*M for _ in range(R)]; W = [[0]*R for _ in range(M)]
    for i in range(M):
        for j in range(M):
            r = i*M + j; U[r][i] = 1; V[r][j] = 1; W[mul[i][j]][r] = int(sig[i][j])
    return U, V, W

def _omega_uvw(OM, M):                                   # 経路 = XOR
    mul = [[i ^ j for j in range(M)] for i in range(M)]
    return wiring_to_uvw(mul, [[int(OM[i][j]) for j in range(M)] for i in range(M)], M)

def _cyclic_uvw(M):                                      # 経路 = mod 足し算・符号 +1
    mul = [[(i+j) % M for j in range(M)] for i in range(M)]
    return wiring_to_uvw(mul, [[1]*M for _ in range(M)], M)

def _matmul_uvw(n):                                      # n×n 素朴行列積の (U,V,W)・R=n³
    M = n*n; U = []; V = []; Wrows = [[] for _ in range(M)]; R = 0
    Ur = []; Vr = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                u = [0]*M; u[i*n+k] = 1; Ur.append(u)     # a_ik
                v = [0]*M; v[k*n+j] = 1; Vr.append(v)     # b_kj
                for c in range(M): Wrows[c].append(1 if c == i*n+j else 0)
                R += 1
    return Ur, Vr, Wrows

def _clifford_omega(n):                                  # 幾何積の符号（Cl(n,0)・XOR経路）
    M = 1 << n
    def sign(A, B):
        s = 0; a = A >> 1
        while a:
            s += bin(a & B).count('1'); a >>= 1          # 並べ替えの 転倒数（e_i²=+1）
        return -1 if (s & 1) else 1
    return np.array([[sign(i, j) for j in range(M)] for i in range(M)], dtype=int)

def _rsqrt_int(x, F, k):                                 # Newton 逆平方根 1/√x（下界=ge健全）
    S = 1 << F; y = S >> ((x.bit_length() + 1) // 2)
    if y == 0: y = 1
    for _ in range(k):                                   # y ← y(3 − x y²)/2（固定小数）
        xyy = x * ((y*y + S - 1) >> F)                   # x·⌈y²/S⌉（上丸め）⟹ 3−xy² 過小 ⟹ 下界保持
        w = 3*S - xyy
        y = (y * w) >> (F + 1) if w > 0 else y           # ⌊·⌋（下丸め）
    return y, S


# ============================================================ 登録
# 群（Cayley–Dickson・XOR経路）: 複素/四元/セデニオン
for M, nm in [(2, "complex"), (4, "quaternion"), (16, "sedenion")]:
    OM = cd_omega(M); U, V, W = _omega_uvw(OM, M)
    _reg_bilinear(nm, U, V, W, M, M,
                  (lambda OM, M: (lambda a, b: ref_mult_M(a, b, OM, M)))(OM, M),
                  f"除算代数({nm})", M*M)
# 行列積（symplectic・XOR経路）: 2×2, 4×4
for q in (1, 2):
    M = 1 << (2*q); OM, _, _ = symplectic_omega(q); U, V, W = _omega_uvw(OM, M)
    _reg_bilinear(f"matmul{1<<q}x{1<<q}", U, V, W, M, M,
                  (lambda OM, M: (lambda a, b: ref_mult_M(a, b, OM, M)))(OM, M),
                  f"{1<<q}×{1<<q} 行列積", M*M)
# 巡回畳み込み（cyclic）
for M in (4, 8):
    U, V, W = _cyclic_uvw(M)
    _reg_bilinear(f"conv_cyclic{M}", U, V, W, M, M,
                  (lambda M: (lambda a, b: [sum(a[i]*b[(k-i) % M] for i in range(M)) for k in range(M)]))(M),
                  f"巡回畳み込み Z/{M}", M*M)
# ランク削減 双線形（明示 U,V,W）
_reg_bilinear("strassen2x2",
    [[1,0,0,1],[0,0,1,1],[1,0,0,0],[0,0,0,1],[1,1,0,0],[-1,0,1,0],[0,1,0,-1]],
    [[1,0,0,1],[1,0,0,0],[0,1,0,-1],[-1,0,1,0],[0,0,0,1],[1,1,0,0],[0,0,1,1]],
    [[1,0,0,1,-1,0,1],[0,0,1,0,1,0,0],[0,1,0,1,0,0,0],[1,-1,1,0,0,1,0]], 4, 4,
    lambda a, b: list((np.array(a).reshape(2,2) @ np.array(b).reshape(2,2)).flatten()),
    "2×2 行列積(Strassen)", 7)
_reg_bilinear("karatsuba",
    [[1,0],[0,1],[1,1]], [[1,0],[0,1],[1,1]], [[1,0,0],[-1,-1,1],[0,1,0]], 2, 2,
    lambda a, b: [a[0]*b[0], a[0]*b[1]+a[1]*b[0], a[1]*b[1]], "多項式積(Karatsuba)", 3)
_reg_bilinear("gauss_complex",
    [[1,1],[1,0],[0,1]], [[1,0],[-1,1],[1,1]], [[1,0,-1],[1,1,0]], 2, 2,
    lambda a, b: [a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]], "複素積(Gauss 3積)", 3)
# 3×3 行列積（素朴 R=27・任意サイズを 配線に）
_U33, _V33, _W33 = _matmul_uvw(3)
_reg_bilinear("matmul3x3", _U33, _V33, _W33, 9, 9,
    lambda a, b: list((np.array(a).reshape(3,3) @ np.array(b).reshape(3,3)).flatten()),
    "3×3 行列積(素朴)", 27)
# Clifford 幾何代数（XOR経路・幾何積の符号）: Cl(2), Cl(3)
for n, nm in [(2, "clifford2"), (3, "clifford3")]:
    M = 1 << n; OM = _clifford_omega(n); U, V, W = _omega_uvw(OM, M)
    _reg_bilinear(nm, U, V, W, M, M,
        (lambda OM, M: (lambda a, b: ref_mult_M(a, b, OM, M)))(OM, M),
        f"幾何代数 Cl({n},0)", M*M)
# 負巡回畳み込み（x^M+1 の 多項式積）: 巻き付きに −1。M=2 は 複素数 そのもの（一般化）
for M in (2, 4, 8):
    mul = [[(i+j) % M for j in range(M)] for i in range(M)]
    sig = [[(-1 if i+j >= M else 1) for j in range(M)] for i in range(M)]
    U, V, W = wiring_to_uvw(mul, sig, M)
    _reg_bilinear(f"negacyclic{M}", U, V, W, M, M,
        (lambda M: (lambda a, b: [sum(a[i]*b[(k-i) % M]*(-1 if i+((k-i) % M) >= M else 1)
                                      for i in range(M)) for k in range(M)]))(M),
        f"負巡回畳み込み x^{M}+1" + ("（=複素）" if M == 2 else ""), M*M)
# 二面体群 D4（位数8・非可換・非XOR経路）: 回転 r^i と 鏡映 s r^i。初の 一般群経路
def _d4_table():
    # 元 = (f, r): f∈{0,1}(鏡映), r∈Z/4。積 (f1,r1)(f2,r2) = (f1^f2, r2 + (−1)^{f2}·r1 …
    # 正: s·r = r^{-1}·s ⟹ (f1,r1)(f2,r2) = (f1⊕f2, (r1·(-1)^{f2}? ) — 具体で 定義して 検証で 固める
    def mulg(g, h):
        f1, r1 = g; f2, r2 = h
        return (f1 ^ f2, (r2 + (r1 if f2 == 0 else -r1)) % 4)
    elems = [(f, r) for f in (0, 1) for r in range(4)]
    idx = {g: i for i, g in enumerate(elems)}
    return [[idx[mulg(elems[i], elems[j])] for j in range(8)] for i in range(8)], elems
_d4_mul, _d4_elems = _d4_table()
_U_d4, _V_d4, _W_d4 = wiring_to_uvw(_d4_mul, [[1]*8 for _ in range(8)], 8)
def _d4_ref(a, b):
    c = [0]*8
    for i in range(8):
        for j in range(8):
            c[_d4_mul[i][j]] += a[i]*b[j]
    return c
_reg_bilinear("dihedral8", _U_d4, _V_d4, _W_d4, 8, 8, _d4_ref, "二面体群 D4（非可換・非XOR）", 64)
# Strassen 4×4 = Strassen ⊗ Strassen（クロネッカー合成・R=49 vs 素朴64）
def _strassen_kron():
    e = REGISTRY["strassen2x2"]; U1, V1, W1 = np.array(e['U']), np.array(e['V']), np.array(e['W'])
    def colmap(n):                                       # 2 段 flatten の 添字合わせ
        # 入力 4×4 の flatten (row*4+col) → テンソル積の 列 ((i1*2+k1)*4 + (i2*2+k2))
        P = np.zeros((16, 16), dtype=int)
        for i1 in range(2):
            for i2 in range(2):
                for k1 in range(2):
                    for k2 in range(2):
                        src = (i1*2+i2)*4 + (k1*2+k2)     # 4×4 の 素直な flatten
                        dst = (i1*2+k1)*4 + (i2*2+k2)     # kron の 列
                        P[dst, src] = 1
        return P
    P = colmap(4)
    U4 = np.kron(U1, U1) @ P; V4 = np.kron(V1, V1) @ P
    W4 = P.T @ np.kron(W1, W1)                           # 出力も 同じ 置換で 戻す
    return U4.tolist(), V4.tolist(), W4.tolist()
_U44, _V44, _W44 = _strassen_kron()
_reg_bilinear("strassen4x4", _U44, _V44, _W44, 16, 16,
    lambda a, b: list((np.array(a).reshape(4,4) @ np.array(b).reshape(4,4)).flatten()),
    "4×4 行列積(Strassen⊗Strassen)", 49)
# 多層 算術回路: Newton 逆数・逆平方根（近似＋ge境界）
def _newton_div(a, b, F=24, k=6):
    y = newton_recip_int(abs(b), F, k)                   # 1/|b| の 下界
    return Fr(a * y, 1 << F) * (1 if b > 0 else -1), "ge"
def _newton_rsqrt(x, _b=None, F=24, k=8):
    y, S = _rsqrt_int(x, F, k)                            # 1/√x の 下界
    return Fr(y, S), "ge"
_reg_circuit("newton_recip", _newton_div, "近似除算 a/b＋ge境界(多層)",
             dict(args=(22, 7), true=22/7))
_reg_circuit("newton_rsqrt", _newton_rsqrt, "近似逆平方根 1/√x＋ge境界(多層)",
             dict(args=(22, None), true=1/22**0.5))


# ============================================================ 統一 実行 ＆ カタログ
def run(name, a, b, backend='gate', K=14):
    """名前で 配線を 選んで 実行。bilinear は gate/algebra、circuit は callable。"""
    e = REGISTRY[name]
    if e['kind'] == 'circuit':
        return e['run'](a, b)
    U, V, W = e['U'], e['V'], e['W']
    if backend == 'algebra':
        r = np.array(W) @ ((np.array(U) @ np.array(a)) * (np.array(V) @ np.array(b)))
        return [int(round(v)) for v in r]
    st = new_counter()
    return [from_sd(o) for o in bilinear_unit_gates(U, V, W, list(a), list(b), K, st)]

def catalog():
    print(f"  {'名前':<16}{'種類':<10}{'R':>5}   計算内容")
    for nm, e in REGISTRY.items():
        print(f"  {nm:<16}{e['kind']:<10}{str(e['R'] or '-'):>5}   {e['computes']}")


def self_test():
    rng = np.random.default_rng(20260731)
    print("="*78)
    print("配線レジストリ — 名前で 配線を 選んで 実行。全パターン 一括検証")
    print("="*78)
    catalog()
    print()
    # 配線正規形の門番: bilinear の (U,V,W) は 三値 {−1,0,+1} が 契約。
    # 三値なら 配線段は 経路+符号+加算だけ（係数の 丸めが どのバックエンドにも 無い）
    # ＝ R が「真の乗算回数」の 正直な請求書になる。二進有理は べき2対角に 括り出し、
    # 無理数係数は CSD 三値化+誤差票つきでのみ 入場（TBM_SPEC「配線正規形」参照）。
    for nm, e in REGISTRY.items():
        if e['kind'] != 'bilinear':
            continue
        vals = set()
        for kk in ('U', 'V', 'W'):
            vals |= set(np.unique(np.array(e[kk], dtype=float)).tolist())
        assert vals <= {-1.0, 0.0, 1.0}, f"三値正規形 破れ: {nm} 係数 {sorted(vals)}"
    print("  三値門番: bilinear 全件 (U,V,W) ⊆ {−1,0,+1} ✓（配線正規形）")
    print()
    print("  検証（algebra=(U,V,W)代数、gate=ゲート版）:")
    small_gate = {"complex","quaternion","matmul2x2","conv_cyclic4","strassen2x2",
                  "karatsuba","gauss_complex","clifford2"}
    for nm, e in REGISTRY.items():
        if e['kind'] == 'circuit':
            d = e['demo']; q, fl = run(nm, *d['args'])
            approx = abs(float(q) - d['true']) < abs(d['true'])*0.02 + 1e-9
            sound = float(q) <= d['true'] + 1e-9          # ge: 近似 ≤ 真（下界）
            print(f"  {nm:<16} circuit: ≈{float(q):.5f} (真{d['true']:.5f}) {fl}  健全&近似 "
                  f"{'✓' if (approx and sound) else '×'}")
            continue
        m, n = e['m'], e['n']
        bad_a = bad_g = 0
        for _ in range(300):
            a = [int(v) for v in rng.integers(-9, 10, m)]
            b = [int(v) for v in rng.integers(-9, 10, n)]
            ref = e['ref'](a, b)
            if run(nm, a, b, backend='algebra') != ref: bad_a += 1
            if nm in small_gate and [from_sd(o) for o in
                bilinear_unit_gates(e['U'], e['V'], e['W'], a, b, 13, new_counter())] != ref: bad_g += 1
        g = "gate✓" if (nm in small_gate and bad_g == 0) else ("gate× " if nm in small_gate else "gate—")
        print(f"  {nm:<16} algebra {'✓' if bad_a==0 else f'×{bad_a}':<4} {g}")
    print()
    print("  ⟹ 一つのファブリック＋レジストリ。run(名前, a, b) で 配線を 選んで 実行。")
    print("     bilinear は (U,V,W) 一枚＝群/行列/Strassen/Karatsuba/Gauss、circuit は 多層(Newton)。")


if __name__ == "__main__":
    self_test()
