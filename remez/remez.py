#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証（research ブランチ・未監査）
"""remez — 係数テープの生成器（オフライン・mpmath）。機械に載るのは テープと縮小の配線だけ。

  目的: 初等関数を「指数分解 + 固定配線の積和」で計算するための 最良近似多項式(Remez)を
  作り、係数を 二進格子 2^-P に量子化し、量子化後の誤差を **厳密有理数で証明つきで** 上界する。

  出力（テープ）= gate_series の K1 形式: 係数は (m, E) = m·2^E の合成時定数のリスト。

  三つの部品:
    ① remez()        — 等リップル定理に基づく交換アルゴリズム(重みつき: 相対/絶対)。
                       参照点 n+2 → 線形系で (c, E) → 誤差の極値へ参照点を移動 → 収束まで。
    ② quantize_refit() — 係数を 2^-P 格子に丸める。上位係数から順に丸めて残りを再最適化
                       (Sollya fpminimax の簡易版) + 下位 2 係数を格子上で局所探索。
    ③ supnorm_exact() — **証明**: 量子化済み多項式 p と 目標 g の相対誤差 |p−g|/|g| の
                       上界を、部分区間ごとの Taylor 模型で **Fraction だけで** 計算する
                       (浮動小数点を一切使わない・区間演算ライブラリも使わない)。
                       g は有理係数の冪級数として与え、尾は比判定で有理数の上界に閉じる。
                       ②の「密格子上の実測 ε」は証明でなく、③の ε_rig ≥ 実測 が証明。

  片側版は p ± ε で作る: |p−g| ≤ ε|g| なら g ∈ [p(1−ε), p(1+2ε)] (ε≤1/2)。ε を 2 の冪 2^-e に
  切り上げれば ±p·2^-e は指数の付け替え(ゲート 0 個)+加算 1 回で作れる。テープは e も運ぶ。

  目標関数は「縮小後の核」を 冪級数で持つ（全部 原点で正則・有理係数）:
    expm1(r)/r = Σ r^k/(k+1)!            sin(√y)/√y = Σ (−1)^k y^k/(2k+1)!
    cos(√y)    = Σ (−1)^k y^k/(2k)!      log1p(t)/t = Σ (−1)^k t^k/(k+1)
    √(1+t)     = Σ C(½,k) t^k            1/√(1+t)   = Σ C(−½,k) t^k
"""
from fractions import Fraction as Fr
import math, json, sys, os
from mpmath import mp, mpf, matrix, lu_solve, cos as mcos, pi as mpi

mp.dps = 50

def _mp(x):
    """Fraction/int/str/mpf → mpf（Fraction は分子/分母で厳密に）。"""
    if isinstance(x, Fr): return mpf(x.numerator) / x.denominator
    return mpf(x)


# ============================================================ 目標関数 = 有理係数の冪級数
class Series:
    """g(x) = Σ_k a_k x^k（a_k 有理）。 |a_{k+1}/a_k| ≤ q(k) を与えて尾を閉じる。
       名前と mpmath 版（Remez の評価用・高精度）も持つ。"""
    def __init__(self, name, coef, ratio, mpf_fn, note=""):
        self.name, self.coef, self.ratio, self.mpf_fn, self.note = name, coef, ratio, mpf_fn, note

    # --- 厳密: j 階導関数を 点 c で（有理数の区間 [lo,hi]）。 x^k の j 階導関数 = k!/(k−j)! x^{k−j}
    def deriv_at(self, c, j, K=None):
        c = Fr(c); K = K or 80
        s = Fr(0)
        for k in range(j, K):
            s += self.coef(k) * Fr(math.perm(k, j)) * c ** (k - j)
        tail = self._tail(abs(c), j, K)
        return s - tail, s + tail

    # --- 厳密: sup_{|x|≤R} |g^{(j)}(x)| の上界（粗いが正しい）
    def deriv_sup(self, R, j, K=None):
        R = Fr(R); K = K or 80
        s = Fr(0)
        for k in range(j, K):
            s += abs(self.coef(k)) * Fr(math.perm(k, j)) * R ** (k - j)
        return s + self._tail(R, j, K)

    def _tail(self, R, j, K):
        """Σ_{k≥K} |a_k| k!/(k−j)! R^{k−j} の上界。項比 ≤ q(k)·(k+1)/(k+1−j)·R ≤ q(K)·(K+1)/(K+1−j)·R =: ρ < 1 を要求。"""
        R = Fr(R)
        first = abs(self.coef(K)) * Fr(math.perm(K, j)) * R ** (K - j)
        rho = self.ratio(K) * Fr(K + 1, K + 1 - j) * R
        assert rho < 1, ("tail does not close", self.name, float(rho))
        return first / (1 - rho)

    # --- 厳密: sup_{x∈[a,b]} |g^{(j)}(x)| の上界。既定は |x| ≤ max(|a|,|b|) の粗い版。
    def deriv_sup_on(self, a, b, j, K=None):
        return self.deriv_sup(max(abs(Fr(a)), abs(Fr(b))), j, K)


def _sqrt_interval(q, bits=160):
    """有理数 q>0 の平方根を 有理区間 [lo,hi] で（幅 < 2^-bits・整数 isqrt のみ）。"""
    q = Fr(q); num, den = q.numerator, q.denominator
    n = math.isqrt(num * den << (2 * bits))
    lo = Fr(n, den << bits)
    hi = Fr(n + 1, den << bits)
    return lo, hi


class BinomialSeries(Series):
    """g(x) = (1+x)^α（α = ±1/2）。級数は |x|<1 でしか閉じないので 導関数は閉形式で:
         g^{(j)}(x) = α(α−1)…(α−j+1) · (1+x)^{α−j}
       (1+x)^{±1/2} は isqrt の有理区間で挟む（浮動小数点なし）。"""
    def __init__(self, name, alpha2, mpf_fn):
        # alpha2 = 2α ∈ {+1, −1}
        super().__init__(name, _binom_half(alpha2), lambda k: Fr(1), mpf_fn)
        self.alpha = Fr(alpha2, 2)

    def _fall(self, j):
        f = Fr(1)
        for i in range(j): f *= (self.alpha - i)
        return f

    def _pow_alpha_interval(self, x):
        """(1+x)^α の有理区間。"""
        lo, hi = _sqrt_interval(1 + Fr(x))
        if self.alpha > 0: return lo, hi
        return 1 / hi, 1 / lo

    def deriv_at(self, c, j, K=None):
        c = Fr(c)
        assert c > -1
        f = self._fall(j)
        plo, phi = self._pow_alpha_interval(c)
        d = (1 + c) ** j
        vals = (f * plo / d, f * phi / d)
        return min(vals), max(vals)

    def deriv_sup_on(self, a, b, j, K=None):
        # (1+x)^{α−j} は (−1,∞) で単調 → 端点の大きい方（上界側）
        a, b = Fr(a), Fr(b)
        assert a > -1
        f = abs(self._fall(j))
        out = Fr(0)
        for x in (a, b):
            _, phi = self._pow_alpha_interval(x)
            out = max(out, f * phi / (1 + x) ** j)
        return out

    def deriv_sup(self, R, j, K=None):
        return self.deriv_sup_on(-Fr(R), Fr(R), j, K)


def _binom_half(sgn):
    # C(±1/2, k) を逐次生成（有理）
    cache = {0: Fr(1)}
    def coef(k):
        if k not in cache:
            for i in range(max(cache) + 1, k + 1):
                cache[i] = cache[i - 1] * (Fr(sgn, 2) - (i - 1)) / i
        return cache[k]
    return coef

SERIES = {
    # expm1(r)/r
    "expm1_over_x": Series("expm1_over_x", lambda k: Fr(1, math.factorial(k + 1)),
                           lambda k: Fr(1, k + 2), lambda r: (mp.expm1(r) / r if r != 0 else mpf(1))),
    # sin(√y)/√y   (y = r²)
    "sinc_sqrt": Series("sinc_sqrt", lambda k: Fr((-1) ** k, math.factorial(2 * k + 1)),
                        lambda k: Fr(1, (2 * k + 2) * (2 * k + 3)),
                        lambda y: (mp.sin(mp.sqrt(y)) / mp.sqrt(y) if y != 0 else mpf(1))),
    # cos(√y)
    "cos_sqrt": Series("cos_sqrt", lambda k: Fr((-1) ** k, math.factorial(2 * k)),
                       lambda k: Fr(1, (2 * k + 1) * (2 * k + 2)), lambda y: mp.cos(mp.sqrt(y))),
    # log1p(t)/t
    "log1p_over_x": Series("log1p_over_x", lambda k: Fr((-1) ** k, k + 1),
                           lambda k: Fr(k + 1, k + 2), lambda t: (mp.log1p(t) / t if t != 0 else mpf(1))),
    # √(1+t)   （級数は |t|<1 のみ → 導関数は閉形式、BinomialSeries）
    "sqrt1p": BinomialSeries("sqrt1p", +1, lambda t: mp.sqrt(1 + t)),
    # 1/√(1+t)
    "rsqrt1p": BinomialSeries("rsqrt1p", -1, lambda t: 1 / mp.sqrt(1 + t)),
}


# ============================================================ ① Remez 交換アルゴリズム（重みつき）
def _grid(a, b, N):
    """[a,b] の密格子: Chebyshev 分布（端に密）+ 一様 を合わせて昇順に。"""
    a, b = _mp(a), _mp(b)
    pts = set()
    for i in range(N + 1):
        pts.add((a + b) / 2 + (b - a) / 2 * mcos(mpi * i / N))
        pts.add(a + (b - a) * i / N)
    return sorted(pts)

def _poly(coeffs, exps, x):
    return sum(c * x ** e for c, e in zip(coeffs, exps))

def remez(g, a, b, exps, weight="rel", fixed=None, N=None, maxit=40, tol=mpf("1e-12"), verbose=False):
    """min max_{[a,b]} w(x)|p(x)−g(x)|, p = Σ_{e∈exps} c_e x^e + Σ_{fixed}。
       weight: "rel" → w=1/|g|、"abs" → w=1、callable → w(x)。
       戻り: (coeffs[list of mpf, exps 順], E=水平化誤差 |E|, 参照点)"""
    fixed = fixed or {}
    a, b = _mp(a), _mp(b)
    m = len(exps)                                        # 自由係数の数 → 参照点 m+1
    if weight == "rel":   w = lambda x: 1 / abs(g(x))
    elif weight == "abs": w = lambda x: mpf(1)
    else:                 w = weight
    def geff(x): return g(x) - sum(_mp(c) * x ** e for e, c in fixed.items())
    N = N or max(4000, 60 * (m + 1))
    grid = _grid(a, b, N)
    gv = [geff(x) for x in grid]; wv = [w(x) for x in grid]
    # 初期参照点: Chebyshev 極点
    ref = [(a + b) / 2 + (b - a) / 2 * mcos(mpi * i / m) for i in range(m + 1)][::-1]
    coeffs = None; E = mpf(0)
    for it in range(maxit):
        # 線形系: Σ c_e x_i^e − (−1)^i E / w(x_i) = geff(x_i)
        A = matrix(m + 1, m + 1); rhs = matrix(m + 1, 1)
        for i, x in enumerate(ref):
            for j, e in enumerate(exps): A[i, j] = x ** e
            A[i, m] = -((-1) ** i) / w(x)
            rhs[i] = geff(x)
        sol = lu_solve(A, rhs)
        coeffs = [sol[j] for j in range(m)]; E = abs(sol[m])
        # 誤差を密格子で
        err = [wv[i] * (_poly(coeffs, exps, x) - gv[i]) for i, x in enumerate(grid)]
        emax = max(abs(v) for v in err)
        if verbose: print(f"    it{it}: E={mp.nstr(E,6)} max={mp.nstr(emax,6)}")
        if emax - E <= tol * emax: break
        # 極値抽出 → 交互符号の m+1 点へ
        ext = []
        for i in range(len(grid)):
            l = abs(err[i - 1]) if i > 0 else -1
            r = abs(err[i + 1]) if i + 1 < len(grid) else -1
            if abs(err[i]) >= l and abs(err[i]) >= r and err[i] != 0:
                ext.append(i)
        alt = []                                          # 同符号の連続は 大きい方だけ
        for i in ext:
            if alt and (err[alt[-1]] > 0) == (err[i] > 0):
                if abs(err[i]) > abs(err[alt[-1]]): alt[-1] = i
            else: alt.append(i)
        while len(alt) > m + 1:                           # 余分は 小さい隣接対 or 端を落とす
            if len(alt) - (m + 1) >= 2:
                k = min(range(len(alt) - 1), key=lambda j: max(abs(err[alt[j]]), abs(err[alt[j + 1]])))
                del alt[k:k + 2]
            else:
                if abs(err[alt[0]]) < abs(err[alt[-1]]): alt.pop(0)
                else: alt.pop()
        if len(alt) < m + 1:                              # 極値不足(退化) → そのまま返す
            break
        ref = [grid[i] for i in alt]
    return coeffs, E, ref

def max_err(g, a, b, exps, coeffs, weight="rel", N=6000):
    """密格子上の 重みつき誤差の最大値（実測・証明でない）。"""
    if weight == "rel":   w = lambda x: 1 / abs(g(x))
    elif weight == "abs": w = lambda x: mpf(1)
    else:                 w = weight
    grid = _grid(a, b, N)
    return max(abs(w(x) * (_poly(coeffs, exps, x) - g(x))) for x in grid)


# ============================================================ ② 量子化 + 再最適化
def quantize_refit(g, a, b, exps, P, weight="rel", fixed=None, search=2, span=2, verbose=False):
    """係数を 2^-P 格子へ。上位から順に丸め → 残りを Remez で再最適化 → 下位 search 個を
       ±span 格子ステップで総当たり。戻り: ({e: int m}, 実測ε) 係数は m·2^-P。"""
    fixed = dict(fixed or {})
    free = list(exps)
    ints = {}
    coeffs, E, _ = remez(g, a, b, free, weight, fixed)
    eps0 = max_err(g, a, b, free, coeffs, weight)
    for e in sorted(free, reverse=True):
        c = coeffs[free.index(e)]
        m_int = int(mp.nint(c * (1 << P)))
        ints[e] = m_int; fixed[e] = Fr(m_int, 1 << P); free.remove(e)
        if free:
            coeffs, E, _ = remez(g, a, b, free, weight, fixed)
    def err_of(ints_):
        fx = {e: Fr(m, 1 << P) for e, m in ints_.items()}
        es = sorted(fx); cs = [_mp(fx[e]) for e in es]
        return max_err(g, a, b, es, cs, weight)
    eps_q = err_of(ints)
    # 局所探索: 下位 search 個
    low = sorted(exps)[:search]
    import itertools
    best = (eps_q, dict(ints))
    for deltas in itertools.product(range(-span, span + 1), repeat=len(low)):
        if all(d == 0 for d in deltas): continue
        trial = dict(ints)
        for e, d in zip(low, deltas): trial[e] += d
        v = err_of(trial)
        if v < best[0]: best = (v, trial)
    eps_s, ints = best
    if verbose:
        print(f"    Remez ε={mp.nstr(eps0,4)} → 量子化(P={P}) ε={mp.nstr(eps_q,4)} → 局所探索 ε={mp.nstr(eps_s,4)}")
    return ints, eps_s, eps0


# ============================================================ ③ 厳密な上界（Fraction のみ）
def _poly_deriv_at(cf, c, j):
    """cf: {e: Fraction}。p^{(j)}(c) 厳密。"""
    return sum(v * Fr(math.perm(e, j)) * Fr(c) ** (e - j) for e, v in cf.items() if e >= j)

def _poly_deriv_sup(cf, R, j):
    return sum(abs(v) * Fr(math.perm(e, j)) * Fr(R) ** (e - j) for e, v in cf.items() if e >= j)

def supnorm_exact(cf, ser, a, b, J=8, nsub=64, rel=True, K=80, max_depth=12):
    """|p−g| / |g| （rel=False なら |p−g|）の [a,b] 上の上界を Fraction で。
       部分区間 [lo,hi] (中心 c, 半径 ρ) で Taylor 模型:
         |e(x)| ≤ Σ_{j≤J} |e^{(j)}(c)| ρ^j/j! + (sup|p^{(J+1)}| + sup|g^{(J+1)}|) ρ^{J+1}/(J+1)!
       e^{(j)}(c) = p^{(j)}(c) − g^{(j)}(c) は 点評価 → 厳密（g は級数 + 尾の区間）。
       |g| の下界: |g(c)| − sup|g'|·ρ。項が減衰しなければ 区間を二分（適応・常に正しい）。"""
    a, b = Fr(a), Fr(b)
    Rmax = max(abs(a), abs(b))
    supP = _poly_deriv_sup(cf, Rmax, J + 1)
    supG = ser.deriv_sup_on(a, b, J + 1, K)
    supG1 = ser.deriv_sup_on(a, b, 1, K)
    fact = [Fr(math.factorial(j)) for j in range(J + 2)]
    worst = Fr(0); pieces = 0; gmin_all = None
    stack = [(a + (b - a) * i / nsub, a + (b - a) * (i + 1) / nsub, 0) for i in range(nsub)]
    while stack:
        lo, hi, depth = stack.pop()
        c = (lo + hi) / 2; rho = (hi - lo) / 2
        terms = []
        for j in range(J + 1):
            glo, ghi = ser.deriv_at(c, j, K)
            pj = _poly_deriv_at(cf, c, j)
            ej = max(abs(pj - glo), abs(pj - ghi))           # |e^{(j)}(c)| の上界
            terms.append(ej * rho ** j / fact[j])
        rem = (supP + supG) * rho ** (J + 1) / fact[J + 1]
        bound = sum(terms) + rem
        # 適応: 末尾項 or 剰余 が 主項の 1/64 を超えるなら二分（正しさには影響しない）
        if depth < max_depth and (terms[-1] + rem) > terms[0] / 64 and terms[0] > 0:
            stack.append((lo, c, depth + 1)); stack.append((c, hi, depth + 1)); continue
        glo, ghi = ser.deriv_at(c, 0, K)
        gmin = min(abs(glo), abs(ghi)) - supG1 * rho              # 部分区間での |g| の下界
        assert gmin > 0, "g が 0 に近すぎて相対誤差を閉じられない"
        gmin_all = gmin if gmin_all is None else min(gmin_all, gmin)
        if rel: bound = bound / gmin
        worst = max(worst, bound); pieces += 1
    return worst, pieces, gmin_all


# ============================================================ テープ（K1 形式）の生成と保存
def dyadic_ceil_pow2(x):
    """x (Fraction>0) 以上の最小の 2 の冪の指数 e: 2^e ≥ x。"""
    e = 0
    while Fr(2) ** e < x: e += 1
    while Fr(2) ** (e - 1) >= x: e -= 1
    return e

def make_tape(name, ser_key, a, b, exps, P, weight="rel", form="poly", arg="x", J=8, verbose=True, **meta):
    """テープ = 係数 (m, −P) のリスト + 証明つき相対誤差。form/arg は評価器への指示:
         form: "poly" p(u) / "x_times" x·p(u) / "one_plus_x_times" 1 + x·p(u)
         arg : "x" u=x / "x2" u=x²"""
    ser = SERIES[ser_key]
    if verbose: print(f"  [{name}] {ser_key} on [{float(a):.6g},{float(b):.6g}] 次数 exps={exps} P={P} 重み={weight}")
    ints, eps_meas, eps0 = quantize_refit(ser.mpf_fn, a, b, exps, P, weight, verbose=verbose)
    cf = {e: Fr(m, 1 << P) for e, m in ints.items()}
    eps_rig, pieces, gmin = supnorm_exact(cf, ser, a, b, J=J, rel=(weight == "rel"))
    assert eps_rig >= Fr(mp.nstr(eps_meas, 20)) * Fr(999, 1000), "証明上界が実測を下回った(バグ)"
    e_pow2 = dyadic_ceil_pow2(eps_rig)
    if verbose:
        print(f"    ε 実測 {mp.nstr(eps_meas,4)} / 証明 {float(eps_rig):.4g} ({pieces} 区間) / 2^{e_pow2} に切上げ"
              f"  = {-e_pow2} ビット")
    tape = {"name": name, "series": ser_key, "interval": [str(Fr(a)), str(Fr(b))],
            "form": form, "arg": arg, "P": P, "weight": weight,
            "coeffs": [[e, m, -P] for e, m in sorted(ints.items())],
            "eps_measured": mp.nstr(eps_meas, 12), "eps_remez_unquantized": mp.nstr(eps0, 12),
            "eps_rigorous": str(_ceil_to_bits(eps_rig, 120)),          # 2^-120 粒度で切上げ(上界を保つ)
            "eps_rigorous_float": float(eps_rig), "eps_pow2_exponent": e_pow2,
            "gmin_rigorous": str(_floor_to_bits(gmin, 120)), "gmin_float": float(gmin), **meta}
    return tape

def _floor_to_bits(x, bits):
    """x を 2^-bits 粒度で切り下げ（下界を保つ）。"""
    q = Fr(1, 1 << bits)
    return Fr(math.floor(x / q)) * q

def _ceil_to_bits(x, bits):
    """x を 2^-bits 粒度で切り上げ（上界を保つ）。"""
    q = Fr(1, 1 << bits)
    n = x / q
    return Fr(math.ceil(n)) * q

def save_tape(tape, path):
    with open(path, "w") as f: json.dump(tape, f, indent=1, ensure_ascii=False)

def load_tape(path):
    with open(path) as f: return json.load(f)

def tape_coeffs_fr(tape):
    """{e: Fraction}"""
    return {e: Fr(m) * Fr(2) ** E for e, m, E in tape["coeffs"]}


# ============================================================ self-test
def self_test():
    print("=" * 80)
    print("remez.py self-test — 交換アルゴリズム / 量子化再最適化 / 厳密上界")
    print("=" * 80)
    ser = SERIES["expm1_over_x"]
    a, b = -Fr(35, 100), Fr(35, 100)
    # ① 等リップル: 参照点での誤差が交互符号・等大きさ
    exps = list(range(0, 8))
    coeffs, E, ref = remez(ser.mpf_fn, a, b, exps, "rel")
    errs = [(_poly(coeffs, exps, x) - ser.mpf_fn(x)) / abs(ser.mpf_fn(x)) for x in ref]
    alt = all((errs[i] > 0) != (errs[i + 1] > 0) for i in range(len(errs) - 1))
    lev = max(abs(abs(v) - E) for v in errs) / E
    em = max_err(ser.mpf_fn, a, b, exps, coeffs, "rel")
    print(f"① expm1(r)/r 次数7 on [-.35,.35]: E={mp.nstr(E,4)} 交互符号={alt} 水平化ずれ={mp.nstr(lev,2)}"
          f" 密格子max={mp.nstr(em,4)} (max/E={mp.nstr(em/E,6)})")
    assert alt and lev < 1e-8 and em / E < 1 + 1e-6
    # ② 量子化 + 再最適化: 素の丸めより良く、P を増やせば Remez に近づく
    for P in (30, 40):
        ints, eps_s, eps0 = quantize_refit(ser.mpf_fn, a, b, exps, P, "rel")
        naive = {e: int(mp.nint(c * (1 << P))) for e, c in zip(exps, coeffs)}
        cs = [mpf(naive[e]) / (1 << P) for e in exps]
        eps_naive = max_err(ser.mpf_fn, a, b, exps, cs, "rel")
        print(f"② P={P}: Remez ε={mp.nstr(eps0,3)}  素の丸め ε={mp.nstr(eps_naive,3)}  再最適化 ε={mp.nstr(eps_s,3)}"
              f"  (再最適化/Remez = {mp.nstr(eps_s/eps0,3)})")
        assert eps_s <= eps_naive * mpf("1.000001")
    # ③ 厳密上界 ≥ 実測、かつ 実測の数倍以内（緩すぎない）
    ints, eps_s, _ = quantize_refit(ser.mpf_fn, a, b, exps, 40, "rel")
    cf = {e: Fr(m, 1 << 40) for e, m in ints.items()}
    rig, pieces, gmin = supnorm_exact(cf, ser, a, b, J=8)
    ratio = float(rig) / float(eps_s)
    print(f"③ 厳密上界 ε_rig={float(rig):.4g} ({pieces}区間) / 実測 {mp.nstr(eps_s,4)} = {ratio:.4f}")
    assert 1.0 <= ratio < 1.5
    # ③' 陰性対照: 係数を 1 格子だけ壊すと 上界も上がる（検出できる）
    bad = dict(cf); bad[0] += Fr(1, 1 << 40)
    rig_bad, _, _ = supnorm_exact(bad, ser, a, b, J=8)
    print(f"③' 陰性対照(c0 を 2^-40 ずらす): ε_rig={float(rig_bad):.4g} (> {float(rig):.4g}) ✓" if rig_bad > rig else "✗")
    assert rig_bad > rig
    # ④ 全級数の尾が閉じ、mpmath 値と級数値が一致する（二経路の突き合わせ）
    for key, pt in [("expm1_over_x", Fr(3, 10)), ("sinc_sqrt", Fr(6, 10)), ("cos_sqrt", Fr(6, 10)),
                    ("log1p_over_x", Fr(1, 2)), ("sqrt1p", Fr(1, 2)), ("rsqrt1p", Fr(1, 2)),
                    ("log1p_over_x", -Fr(1, 4)), ("sqrt1p", -Fr(1, 2)), ("rsqrt1p", -Fr(1, 2))]:
        s = SERIES[key]; lo, hi = s.deriv_at(pt, 0, 160)
        v = s.mpf_fn(_mp(pt)); tol = mpf(10) ** (-(mp.dps - 5))     # mpmath 側の丸めだけ許す
        ok = _mp(lo) - tol <= v <= _mp(hi) + tol
        print(f"④ {key:<14} x={float(pt):+.3f}: 級数区間幅 {float(hi-lo):.1e} ∋ mpmath値(±1e-{mp.dps-5}) {'✓' if ok else '✗'}")
        assert ok
    print("self-test 通過")

if __name__ == "__main__":
    self_test()
