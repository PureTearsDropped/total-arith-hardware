#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証（research ブランチ・未監査）
"""tape_eval — 係数テープを食う評価器（Estrin・固定配線・深さ O(log n)）。三つの層が **同じ木** を歩く:

  Level 1  spec  : 二進小数 (N, E) = N·2^E の整数演算。積・和は厳密、各節点の後で Wk 桁に **切り捨て**
                   (canonical な大きさの下位桁を落とす = ゲートの block_normalize と同じ向き)。
  Level 0  numpy : float64 で同じ木・同じ切り捨てを模倣（Wk ≤ 52 の構成でだけ bit 一致; 総当たり用）。
  Level 2  gate  : gate_bfp の BF（符号つき桁列 + ホスト指数）で同じ木。K1: 係数は合成時定数 ⟹
                   定数×信号 の部分積は配線（桁 ±1 は選択/neg、0 は無し）。信号×信号だけ multiply。

  Estrin 木 (係数 c_0..c_n, 引数 u):
     節点(lo,hi,ℓ) = 節点(lo, mid) + 節点(mid, hi)·u^(2^ℓ)    半分ずつ、深さ ⌈log2(n+1)⌉
     葉 = 係数(定数)。 u^(2^ℓ) は 前の冪の二乗（切り捨てつき）。
     零係数の右半分は枝刈り(積 0 は配線)。

  誤差の勘定 (estrin_error_bound): 各切り捨て = 相対 τ = 2^-(Wk−1) 以下。冪 u^(2^ℓ) の相対誤差
  ρ_ℓ = (1+τ)^(2^ℓ−1) − 1。節点の絶対誤差を |c|・R=max|u| で上から評価して根まで運ぶ。全部 Fraction。
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from fractions import Fraction as Fr
import math


# ============================================================ 木の構造（三層で共有）
def estrin_plan(n):
    """次数 n → (L, nodes)。L = 冪の段数 (u, u², …, u^(2^(L−1)))。
       nodes は 後順(子→親) のリスト: (id, kind, lo_id, hi_id, pow_level) 。kind: 'leaf' k / 'mac'。"""
    m = n + 1
    L = max(1, math.ceil(math.log2(m))) if m > 1 else 0
    size = 1 << L
    nodes = []
    def build(lo, hi, lvl):
        # 区間 [lo,hi) の係数 (幅 2^lvl)。 戻り: node id か None(全零)
        if hi - lo == 1:
            if lo > n: return None
            nid = len(nodes); nodes.append((nid, 'leaf', lo, None, None)); return nid
        mid = (lo + hi) // 2
        a = build(lo, mid, lvl - 1); b = build(mid, hi, lvl - 1)
        if b is None: return a                          # 右半分が全零 → 積は配線 (0)
        nid = len(nodes); nodes.append((nid, 'mac', a, b, lvl - 1)); return nid
    root = build(0, size, L)
    return L, nodes, root

def plan_stats(n):
    """(信号×信号 の積の数, 定数×信号 の積の数, MAC 深さ)。冪の二乗も 信号×信号。"""
    L, nodes, root = estrin_plan(n)
    sig = max(0, L - 1)                                   # u², u⁴, … (u 自体は入力)
    const = 0
    depth = {}
    for nid, kind, a, b, lvl in nodes:
        if kind == 'leaf': depth[nid] = 0; continue
        if nodes[b][1] == 'leaf': const += 1              # 係数×u^(2^lvl): 定数積
        else: sig += 1
        depth[nid] = 1 + max(depth[a] if a is not None else 0, depth[b])
    return sig, const, depth[root] if root is not None else 0


# ============================================================ Level 1: 二進小数 spec
def dy_norm(N, E):
    """N·2^E を 末尾の 0 を落として正規化（表現の一意化・値不変）。"""
    if N == 0: return (0, 0)
    t = (N & -N).bit_length() - 1
    return (N >> t, E + t)

def dy_mul(x, y): return dy_norm(x[0] * y[0], x[1] + y[1])

def dy_add(x, y):
    E = min(x[1], y[1])
    return dy_norm((x[0] << (x[1] - E)) + (y[0] << (y[1] - E)), E)

def dy_trunc(x, W):
    """|N| を 上位 W 桁に切り捨て（0 方向）。戻り ((N',E'), dropped_nonzero)。
       = canonicalize → 先頭位置 → 低位 sh 桁を落とす（block_normalize の 単成分版）。"""
    N, E = x
    a = abs(N); L = a.bit_length() - 1
    sh = max(0, L - (W - 1))
    if sh == 0: return x, 0
    q = a >> sh; dropped = int((a & ((1 << sh) - 1)) != 0)
    return dy_norm(q if N > 0 else -q, E + sh), dropped

def dy_from_fr(f):
    """Fraction (二進小数) → (N,E)。分母が 2 の冪でなければ例外。"""
    d = f.denominator; assert d & (d - 1) == 0, "not dyadic"
    return dy_norm(f.numerator, -(d.bit_length() - 1))

def dy_to_fr(x): return Fr(x[0]) * Fr(2) ** x[1]

def estrin_spec(tape_coeffs, u, Wk):
    """tape_coeffs: [(k, m, E)] 昇順、u: (N,E)。戻り (N,E) = 木の値（各節点で Wk 桁切り捨て）。"""
    n = max(k for k, m, E in tape_coeffs)
    cd = {k: dy_norm(m, E) for k, m, E in tape_coeffs}
    L, nodes, root = estrin_plan(n)
    pw = [u]
    for l in range(1, L):
        pw.append(dy_trunc(dy_mul(pw[-1], pw[-1]), Wk)[0])
    val = {}
    for nid, kind, a, b, lvl in nodes:
        if kind == 'leaf': val[nid] = cd.get(a, (0, 0)); continue
        acc = dy_mul(val[b], pw[lvl])
        if a is not None: acc = dy_add(val[a], acc)
        val[nid] = dy_trunc(acc, Wk)[0]
    return val[root]


# ============================================================ 誤差の上界（Fraction）
def estrin_error_bound(tape_coeffs, R, Wk):
    """|u| ≤ R での 木の絶対誤差の上界（切り捨てだけ; 係数は厳密）と |値| の上界。"""
    n = max(k for k, m, E in tape_coeffs)
    cd = {k: abs(Fr(m) * Fr(2) ** E) for k, m, E in tape_coeffs}
    R = Fr(R); tau = Fr(1, 1 << (Wk - 1))
    L, nodes, root = estrin_plan(n)
    rho = [Fr(0)] + [(1 + tau) ** ((1 << l) - 1) - 1 for l in range(1, L)]   # u^(2^l) の相対誤差
    V = {}; Eb = {}
    for nid, kind, a, b, lvl in nodes:
        if kind == 'leaf': V[nid] = cd.get(a, Fr(0)); Eb[nid] = Fr(0); continue
        Rp = R ** (1 << lvl)
        v = V[b] * Rp; e = Eb[b] * Rp * (1 + rho[lvl]) + V[b] * Rp * rho[lvl]
        if a is not None: v += V[a]; e += Eb[a]
        e += tau * v                                        # 自身の切り捨て (|値| ≤ v)
        V[nid] = v; Eb[nid] = e
    return Eb[root], V[root]


# ============================================================ Level 0: numpy 模倣（Wk ≤ 52）
def np_trunc(v, W):
    """float64 配列を 上位 W 桁に 0 方向へ切り捨て（W ≤ 52 で厳密）。"""
    import numpy as np
    m, e = np.frexp(v)                                    # v = m·2^e, 0.5 ≤ |m| < 1
    s = np.ldexp(1.0, W)                                  # m·2^W は 整数部 W 桁
    return np.ldexp(np.trunc(m * s), e - W)

def estrin_np(tape_coeffs, u, Wk):
    import numpy as np
    n = max(k for k, m, E in tape_coeffs)
    cd = {k: math.ldexp(m, E) for k, m, E in tape_coeffs}
    L, nodes, root = estrin_plan(n)
    pw = [u]
    for l in range(1, L):
        pw.append(np_trunc(pw[-1] * pw[-1], Wk))
    val = {}
    for nid, kind, a, b, lvl in nodes:
        if kind == 'leaf': val[nid] = cd.get(a, 0.0); continue
        acc = val[b] * pw[lvl]
        if a is not None: acc = val[a] + acc
        val[nid] = np_trunc(acc, Wk)
    return val[root]


# ============================================================ Level 2: ゲート（BF, ホスト指数）
from gate_bfp import BF, from_bf, shift_up, bf_mul, bf_add
from gate_bilinear import to_sd, from_sd, neg, sd_sum, canonicalize, ZERO, nz, OR, AND, NOT, dec

def const_digits(m):
    """合成時定数 m → 符号つき桁（ホスト int 0/1、ゲートを数えない）。"""
    return to_sd(m, max(1, abs(m).bit_length()))

def mul_const(X, cdig, st):
    """信号 X × 定数桁列（K1）: 桁 ±1 は 選択/neg の配線、0 は無し ⟹ 部分積のゲート 0・集約だけ数える。"""
    rows = []
    for j, (p, n) in enumerate(cdig):
        if p == 1:   rows.append([ZERO] * j + list(X))
        elif n == 1: rows.append([ZERO] * j + neg(X))
    if not rows: return [ZERO]
    if len(rows) == 1: return rows[0]
    return sd_sum(rows, st)

def bf_trunc(x, W, st):
    """BF を canonicalize → 先頭位置(ホスト) → 低位を落として W 桁。戻り (BF, dropped_nz)。"""
    c, s = canonicalize(x.mant, st)
    L = -1
    for i in range(len(c) - 1, -1, -1):
        if dec(*c[i]) != 0: L = i; break
    sh = max(0, L - (W - 1))
    dropped = 0
    for i in range(sh): dropped = OR(dropped, nz(c[i], st), st)
    kept = [c[i] if i < len(c) else ZERO for i in range(sh, sh + W)]
    return BF(kept, x.E + sh), dropped

def estrin_gate(tape_coeffs, u_bf, Wk, st):
    """u_bf: BF(Wk 桁)。戻り BF(Wk 桁)。"""
    n = max(k for k, m, E in tape_coeffs)
    cd = {k: (m, E) for k, m, E in tape_coeffs}
    L, nodes, root = estrin_plan(n)
    pw = [u_bf]
    for l in range(1, L):
        pw.append(bf_trunc(bf_mul(pw[-1], pw[-1], st), Wk, st)[0])
    val = {}
    for nid, kind, a, b, lvl in nodes:
        if kind == 'leaf':
            m, E = cd.get(a, (0, 0)); val[nid] = BF(const_digits(m), E); continue
        if nodes[b][1] == 'leaf':                         # 定数 × 信号: 配線 + 集約
            m, E = cd.get(nodes[b][2], (0, 0))
            prod = BF(mul_const(pw[lvl].mant, const_digits(m), st), pw[lvl].E + E)
        else:
            prod = bf_mul(val[b], pw[lvl], st)
        acc = bf_add(val[a], prod, st) if a is not None else prod
        val[nid] = bf_trunc(acc, Wk, st)[0]
    return val[root]


# ============================================================ 出口ユニット: 片側境界と丸め
def finish_spec(P, e, mode, Wout, Emin, Emax, sticky=True):
    """P=(N,E) 計算値, 2^-e = 全相対誤差の上界 (|P−f| ≤ 2^-e·|f|)。
         mode 'lo' : trunc_0(P − P·2^-e)      → |shown| ≤ |f|  → ge
         mode 'hi' : away_0(P + P·2^-(e−1))   → |shown| ≥ |f|  → le
         mode 'near': nearest(P)              → 界なし (11)   (P=0 なら厳密 0)
       sticky=True: P は近似（多項式経路を通った）なので、丸めが厳密でも 00 を主張しない。
       （exp/log/sin/cos の値は x≠0/1 で無理数、sqrt は完全平方でだけ厳密 — 検出器なしでは主張できない）
       戻り ((N,E) 正規化 W_out 桁, (ge, le, sunk))。指数域 [Emin,Emax] は 出力 LSB の指数の範囲:
       上に溢れ → ±MAX(ge)、下に潰れ → ±MIN(le)、下で桁落ち → 切り捨て(ge)。"""
    N, E = P
    assert e >= Wout + 1, "near の溢れ/潰れ判定は e ≥ Wout+1 を前提にする"
    if N == 0: return (0, Emin), (0, 0, 0)
    if mode == 'lo':
        v = dy_add(P, (-N, E - e)); rnd = 'trunc'; base = (1, 0)
    elif mode == 'hi':
        v = dy_add(P, (N, E - e + 1)); rnd = 'away'; base = (0, 1)
    else:
        v = P; rnd = 'near'; base = (1, 1)
    return round_spec(v, rnd, Wout, Emin, Emax, base, sticky)

def round_spec(v, rnd, Wout, Emin, Emax, base=None, sticky=False):
    """(N,E) を W_out 桁・LSB 指数 ∈ [Emin,Emax] に丸める。rnd ∈ trunc/away/near。
       base = 非厳密のときの (ge,le)（既定: trunc→ge, away→le, near→11）。sticky: 入力が既に非厳密。"""
    if base is None: base = {'trunc': (1, 0), 'away': (0, 1), 'near': (1, 1)}[rnd]
    N, E = v
    if N == 0: return (0, Emin), (0, 0, 0)
    sgn = 1 if N > 0 else -1; a = abs(N)
    L = a.bit_length() - 1
    Eo = E + L - (Wout - 1)                               # 正規化したときの LSB 指数
    if Eo < Emin: Eo = Emin                               # 非正規化域: 桁を落として Emin に留める
    sh = Eo - E
    if sh <= 0:
        q = a << (-sh); inexact = 0
    else:
        q = a >> sh; rem = a & ((1 << sh) - 1); inexact = int(rem != 0)
        if inexact:
            if rnd == 'away' or (rnd == 'near' and (rem << 1) >= (1 << sh) and not (rem << 1 == (1 << sh) and q % 2 == 0)):
                q += 1
    if q.bit_length() > Wout:                             # 繰り上がりで桁が増えた
        q >>= 1; Eo += 1                                  # 2^Wout → 2^(Wout−1)·2 (下位は 0 で厳密)
    inexact = inexact or int(bool(sticky))
    ge, le = base if inexact else (0, 0)
    if rnd == 'near' and inexact: ge, le = 1, 1
    if q == 0:                                            # 潰れ
        # away/near: v ≥ f（away）または |f−v| ≤ 2^-e|v| かつ v < MIN/2（near, e ≥ Wout+1）→ f < MIN ⟹ ±MIN・le
        # trunc (lo): v ≤ f しか知らず f ≥ MIN かもしれない → 下界は 0（0・ge = 「|true| ≥ 0」）
        if rnd == 'trunc': return (0, Emin), (1, 0, 0)
        return (sgn, Emin), (0, 1, 0)
    if Eo > Emax:                                         # 溢れ
        # trunc/near: f ≥ v(1−2^-e) ≥ MAX（near は e ≥ Wout+1 で保証）→ ±MAX・ge
        # away (hi): v ≥ f だけで f < MAX かもしれない → 界なし 11
        if rnd == 'away': return (sgn * ((1 << Wout) - 1), Emax), (1, 1, 0)
        return (sgn * ((1 << Wout) - 1), Emax), (1, 0, 0)
    # trunc で inexact → 大きさは小さい側 ⇒ ge; away → le; 厳密 → 00
    return (sgn * q, Eo), (ge, le, 0)


# ============================================================ self-test
def self_test():
    import numpy as np, random
    print("=" * 80); print("tape_eval self-test"); print("=" * 80)
    # 木の統計
    for n in (3, 7, 8, 11, 12, 16, 23):
        sig, const, depth = plan_stats(n)
        print(f"  次数 {n:2d}: 信号積 {sig:2d} 定数積 {const:2d} MAC深さ {depth}")
    # (1) spec == 厳密 Horner (Wk 無限大) を Fraction で
    rng = random.Random(1)
    tape = [(k, rng.randrange(-(1 << 40), 1 << 40), -40) for k in range(0, 12)]
    for _ in range(50):
        u = (rng.randrange(-(1 << 40), 1 << 40), -41)
        v = estrin_spec(tape, u, 10 ** 6)
        uf = dy_to_fr(u); exact = sum(Fr(m) * Fr(2) ** E * uf ** k for k, m, E in tape)
        assert dy_to_fr(v) == exact
    print("  (1) spec(Wk=∞) == 厳密多項式 ✓ (50 乱数)")
    # (2) 切り捨て誤差が 上界の内側
    worst = 0
    for _ in range(200):
        u = (rng.randrange(-(1 << 40), 1 << 40), -41)
        v = estrin_spec(tape, u, 30)
        uf = dy_to_fr(u); exact = sum(Fr(m) * Fr(2) ** E * uf ** k for k, m, E in tape)
        eb, vb = estrin_error_bound(tape, Fr(1, 2), 30)
        err = abs(dy_to_fr(v) - exact); worst = max(worst, err / eb)
        assert err <= eb
    print(f"  (2) 切り捨て誤差 ≤ 上界 ✓ (200 乱数, 最悪 誤差/上界 = {float(worst):.3f})")
    # (3) numpy == spec (Wk=30, 係数 40 ビット, u 41 ビット)
    us = rng.choices(range(-(1 << 40), 1 << 40), k=1000)
    un = np.array([math.ldexp(x, -41) for x in us])
    vn = estrin_np(tape, un, 30)
    for i in range(len(us)):
        v = estrin_spec(tape, (us[i], -41), 30)
        assert vn[i] == float(Fr(v[0]) * Fr(2) ** v[1]), (i, vn[i], v)
    print("  (3) numpy 模倣 == spec ✓ (1000 乱数, bit 一致)")
    # (4) ゲート == spec
    from gate_bilinear import new_counter
    st = new_counter()
    tape8 = [(k, rng.randrange(-(1 << 20), 1 << 20), -20) for k in range(0, 8)]
    for i in range(20):
        ui = rng.randrange(-(1 << 24), 1 << 24)
        ub = BF(to_sd(ui, 25), -25)
        vg = estrin_gate(tape8, ub, 24, st)
        vs = estrin_spec(tape8, dy_norm(ui, -25), 24)
        assert from_bf(vg) == dy_to_fr(vs), (i, from_bf(vg), dy_to_fr(vs))
    tot = sum(st.values())
    print(f"  (4) ゲート == spec ✓ (20 乱数, 次数7 Wk=24; 1 回あたり {tot // 20:,} ゲート)")
    # (5) 出口ユニット: lo ≤ 真 ≤ hi、near は最近接
    for _ in range(2000):
        N = rng.randrange(1, 1 << 60) * rng.choice([1, -1]); E = rng.randrange(-80, 20)
        e = rng.randrange(25, 58)
        # 真値 f は P·(1+δ), |δ| ≤ 2^-e の 任意の点
        d = Fr(rng.randrange(-(1 << 30), (1 << 30) + 1), 1 << 30) * Fr(1, 1 << e)
        f = dy_to_fr((N, E)) * (1 + d)
        (lo, fl), (hi, fh), (nr, fn) = (finish_spec((N, E), e, m, 24, -149, 104) for m in ('lo', 'hi', 'near'))
        lo_v, hi_v, nr_v = (dy_to_fr(v) for v in (lo, hi, nr))
        assert abs(lo_v) <= abs(f) <= abs(hi_v), (f, lo_v, hi_v)
        assert fl[0] == 1 or lo_v == f; assert fh[1] == 1 or hi_v == f
        # near: |nr − P| ≤ 半 ulp
        L = abs(N).bit_length() - 1 + E; ulp = Fr(2) ** max(L - 23, -149)
        assert abs(nr_v - dy_to_fr((N, E))) <= ulp / 2
    print("  (5) 出口ユニット: lo ≤ |f| ≤ hi・near は半ulp ✓ (2000 乱数)")
    # (6) 溢れ・潰れ
    v, f = round_spec((1, 200), 'trunc', 24, -149, 104); assert v == ((1 << 24) - 1, 104) and f == (1, 0, 0)
    v, f = round_spec((1, 200), 'away', 24, -149, 104); assert v == ((1 << 24) - 1, 104) and f == (1, 1, 0)
    v, f = round_spec((1, -300), 'trunc', 24, -149, 104); assert v == (0, -149) and f == (1, 0, 0)
    v, f = round_spec((1, -300), 'away', 24, -149, 104); assert v == (1, -149) and f == (0, 1, 0)
    v, f = round_spec((1, -300), 'near', 24, -149, 104); assert v == (1, -149) and f == (0, 1, 0)
    v, f = round_spec((-3, -150), 'trunc', 24, -149, 104); assert v == (-1, -149) and f == (1, 0, 0)
    print("  (6) 溢れ→MAX(ge; hi は 11)・潰れ→MIN(le; lo は 0·ge)・非正規化 切り捨て(ge) ✓")
    print("self-test 通過")

if __name__ == "__main__":
    self_test()
