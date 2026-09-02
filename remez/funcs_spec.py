#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証（research ブランチ・未監査）
"""funcs_spec — Level 1: 初等関数 exp / expm1 / log / sqrt / rsqrt / sin / cos の **二進小数 spec**。

  機械に載る配線 = 「縮小（指数バス操作）+ テープ（Estrin 積和）+ 出口（片側境界と丸め）」を、
  Python の整数だけで（浮動小数点なしで）定義する。ゲート実装（Level 2）はこれと bit 一致、
  numpy 模倣（Level 0, funcs_np.py）は float32 全入力総当たり用の近似双子。

  入力: x = (N, E) = N·2^E（|N| < 2^W_in、W_in = W_out）と 順序層フラグ (ge, le, sunk)。
  出力: ((N, E) W_out 桁, (ge, le, sunk))。mode ∈ {'lo','hi','near'}:
      'lo'  : |shown| ≤ |f| を保証する下側境界（ge を立てる）
      'hi'  : |shown| ≥ |f| を保証する上側境界（le を立てる）
      'near': 最近接丸め（界なし 11、厳密なら 00）
  出口で使う 2^-e（全相対誤差の証明つき上界）は 関数ごとに `Spec.e[name]` にあり、
  その導出は `_budgets()` に書いてある（テープの ε_rig・Estrin 切り捨て上界・定数の丸め・
  縮小の切り捨て を Fraction で足す）。

  縮小（全部 指数バス + 定数 + mux。値で分岐しない）:
    exp/expm1 : |x| < 1/4（先頭位置 ≤ −3）→ k=0, r=x（厳密）。それ以外 → x を 2^-Pc 格子に切り捨て、
                k = round(x·(1/ln2)_Q), r = x − k·LN2_Pc（固定小数点で厳密）, u = trunc_Wk(r)。
                |x| ≥ 2^XC は ±2^XC に clamp（結果は どのみち 飽和/潰れ）。
                exp = 2^k·trunc_Wk(1 + u·q(u)),  expm1 = trunc_Wk((1 + u·q(u))·2^k − 1)（厳密に引いてから切り捨て）
    log       : x = m·2^Ex, m ∈ [1,2)。上位 8 桁 ≥ 181 (m ≥ 181/128 ≈ √2) なら m/2, Ex+1。
                t = m − 1（厳密）, log = trunc_Wk(Ex·LN2_Pc + t·q(t))。x<0 は log|x| を 11+sunk で返す。x=0 → −MAX (ge)。
    sqrt/rsqrt: Ex 奇数 → m/2, Ex+1, t ∈ [−1/2,0) (片 p1); 偶数 → t = m−1 ∈ [0,1) (片 p0)。
                sqrt = p(t)·2^(Ex/2), rsqrt = p(t)·2^(−Ex/2)（指数の算術シフト）。rsqrt(0) = 0（a/0 = 0 の規約）。
    sin/cos   : |x| < 1/2（先頭位置 ≤ −2）→ k=0, r=x。それ以外 → 窓つき 2/π（Payne–Hanek）:
                T = (M·window(E)) mod 2^(F+2), k = round(T·2^-F), frac = T − k·2^F,
                r = trunc_Wk(trunc_{Wk+4}(frac·2^-F)·(π/2)_Pc)。y = trunc_Wk(r²), s=s(y), c=c(y) を **両方** 評価し、
                象限 k mod 4 と x の符号で mux。F = W_in + fbits + Wk + 2（2^-fbits ≤ min|frac|: f32 は総当たり、
                f64 は文献値 2^-62 を **引用**、fracmin.py）。

  入力フラグの伝播（順序層。単調性から言えるものだけ、言えなければ 11 / 符号が飛ぶなら sunk）:
    exp  : ge∧x>0 or le∧x<0 → f(true) ≥ f(shown);  le∧x>0 or ge∧x<0 → ≤。sunk → 11（exp>0 なので sunk は立てない）
    expm1: 奇関数で |·| 単調 → ge→ge, le→le。sunk → 11+sunk
    log  : x≥1∧ge → ge, x≤1∧le → ge（|log| が増える向き）, x>1∧le / x<1∧ge → 0 を跨ぐので 11+sunk
    sqrt : ge→ge, le→le。 rsqrt: ge→le, le→ge。 sin/cos: 非厳密入力 → 11+sunk
    計算側のフラグと合成: out_ge = (計算が上側を否定しない) ∧ (入力の向きが ge を許す)、le も同様。
    両方立たなければ 11。
"""
import os, sys, json, math
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf
from tape_eval import dy_norm, dy_mul, dy_add, dy_trunc, dy_to_fr, estrin_spec, finish_spec
from tapes import CONFIGS, TAPE_DIR
from remez import dyadic_ceil_pow2, SERIES

HERE = os.path.dirname(os.path.abspath(__file__))
FNS = ("exp", "expm1", "log", "sqrt", "rsqrt", "sin", "cos")


def _round_const(v, bits):
    """mpf → round(v·2^bits)（合成時定数）。"""
    return int(mp.nint(v * mpf(2) ** bits))


def _lead(N):
    return abs(N).bit_length() - 1


def fix_trunc(x, P):
    """(N,E) → trunc_0(N·2^E·2^P) 整数（固定小数点 2^-P 格子・0 方向）。"""
    N, E = x
    if E + P >= 0: return N << (E + P)
    a = abs(N) >> (-(E + P))
    return a if N >= 0 else -a


def window(TOP, NB, E, F):
    """floor(2^(E+F)·2/π) mod 2^(F+2)（2/π の bit 列 TOP = floor(2/π·2^NB) から切り出す）。"""
    s = NB - E - F
    assert s >= 0, "2/π の bit 列が短い"
    return (TOP >> s) & ((1 << (F + 2)) - 1)


class Tape:
    def __init__(self, d):
        self.d = d
        self.coeffs = [tuple(c) for c in d["coeffs"]]
        self.eps = Fr(d["eps_rigorous"])
        self.e_est = Fr(d["estrin_abs_error_bound"])
        self.gmin = Fr(d["gmin_rigorous"])
        self.degree = d["degree"]


class Spec:
    def __init__(self, cfgname="f64", fbits_f64=62):
        c = CONFIGS[cfgname]
        self.cfg = cfgname
        self.Wout = c["Wout"]; self.Win = self.Wout; self.Wk = c["Wk"]
        self.Emin = c["Emin"]; self.Emax = c["Emax"]
        self.tau = Fr(1, 1 << (self.Wk - 1))
        self.Pc = self.Wk + 16                 # 定数 (ln2, π/2) の桁数 = exp の固定小数点格子
        self.Q = self.Wk                       # 1/ln2 の桁数
        mp.dps = 80
        self.LN2 = _round_const(mp.log(2), self.Pc)
        self.INV_LN2 = _round_const(1 / mp.log(2), self.Q)
        self.PIO2 = _round_const(mp.pi / 2, self.Pc)
        # exp の clamp: 2^XC ≥ (max(Emax+Wout+1, −Emin))·ln2 なら |x| ≥ 2^XC は 飽和/潰れ確定
        lim = max(self.Emax + self.Wout + 1, -self.Emin) * float(mp.log(2))
        self.XC = math.ceil(math.log2(lim))
        self.kmax = ((1 << self.XC) * self.INV_LN2 >> self.Q) + 2
        # sin/cos の窓
        if cfgname == "f64":
            self.fbits = fbits_f64; self.fbits_source = "cited (binary64 worst case, Muller; not exhaustive)"
        else:
            fm = json.load(open(os.path.join(TAPE_DIR, f"fracmin_{cfgname}.json")))
            self.fbits = fm["fbits"]; self.fbits_source = "exhaustive (fracmin.py)"
        self.F = self.Win + self.fbits + self.Wk + 2
        self.NB = self.Emax + self.Win + self.F + 8       # 正規化で E は Emax+Win−1 まで上がる
        mp.dps = int(self.NB * 0.302) + 40
        self.TOP = int(mp.floor(2 / mp.pi * mpf(2) ** self.NB))
        mp.dps = 80
        # テープ
        self.tapes = {}
        for k in ("exp", "log", "sqrt_p0", "sqrt_p1", "rsqrt_p0", "rsqrt_p1", "sin", "cos"):
            self.tapes[k] = Tape(json.load(open(os.path.join(TAPE_DIR, f"{k}_{cfgname}.json"))))
        self.MAX = ((1 << self.Wout) - 1, self.Emax)
        self.MIN = (1, self.Emin)
        self.e, self.budget = self._budgets()

    # ------------------------------------------------------------------ 誤差予算（Fraction、証明つきの上界の和）
    def _budgets(self):
        tau, Pc, T = self.tau, self.Pc, self.tapes
        B = {}
        # --- exp:  計算する (1 + u·q̂) と exp(r_true) の絶対差 A（exp(r) ≤ e^0.35 < 1.42）:
        #   x の格子切り捨て 2^-Pc, LN2 の丸め k·2^-(Pc+1), u の切り捨て |r|τ ≤ 0.35τ  → exp の変化 ≤ 1.42·(…)
        #   テープ ε|expm1(u)| ≤ 0.42ε, Estrin |u|e_est ≤ 0.35 e_est
        A = Fr(142, 100) * (Fr(1, 1 << Pc) + self.kmax * Fr(1, 1 << (Pc + 1)) + Fr(35, 100) * tau) \
            + Fr(42, 100) * T["exp"].eps + Fr(35, 100) * T["exp"].e_est
        # 相対: exp(r) ≥ e^-0.35 > 0.704、最後の切り捨て τ
        B["exp"] = A / Fr(704, 1000) + tau
        # --- expm1: 小さい経路 (|x|<1/4, u=x 厳密): ε + e_est/gmin + τ。
        #            大きい経路: |expm1(x)| ≥ 1−e^-0.25 > 0.221, 2^k ≤ 1.42·e^x → 相対 ≤ 6.5A + τ
        B["expm1"] = max(T["exp"].eps + T["exp"].e_est / T["exp"].gmin + tau, Fr(65, 10) * A + tau)
        # --- log: Ex=0: t·q̂ の相対 = ε + e_est/gmin + τ
        #          Ex≠0: |log x| ≥ ln2 − log√2 > 0.34, |Ex| ≤ Exmax, |log1p(t)| ≤ 0.347, |t| ≤ 0.414
        Exmax = max(self.Emax + self.Wout, -self.Emin) + 1
        B["log"] = max(T["log"].eps + T["log"].e_est / T["log"].gmin + tau,
                       (Exmax * Fr(1, 1 << (Pc + 1)) + Fr(347, 1000) * T["log"].eps + Fr(414, 1000) * T["log"].e_est)
                       / Fr(34, 100) + tau)
        # --- sqrt / rsqrt: 根の切り捨ては e_est に含まれる
        for g in ("sqrt", "rsqrt"):
            B[g] = max(T[f"{g}_p{i}"].eps + T[f"{g}_p{i}"].e_est / T[f"{g}_p{i}"].gmin for i in (0, 1))
        # --- sin / cos: 大きい経路の u=r の相対誤差:
        #   frac の絶対誤差 ≤ 2^(Win−F), |frac| ≥ 2^-fbits → 2^(Win−F+fbits) = 2^-(Wk+2);  frac の切り捨て 2^-(Wk+3);
        #   (π/2)_Pc の丸め 2^-(Pc+1);  最後の切り捨て τ
        ru = (1 + Fr(1, 1 << (self.F - self.Win - self.fbits))) * (1 + Fr(1, 1 << (self.Wk + 3))) \
             * (1 + Fr(1, 1 << (Pc + 1))) * (1 + tau) - 1
        ry = (1 + ru) ** 2 * (1 + tau) - 1                    # y = trunc(u²)
        ymax = Fr(31, 50)
        sup_s1 = SERIES["sinc_sqrt"].deriv_sup_on(0, ymax, 1)
        sup_c1 = SERIES["cos_sqrt"].deriv_sup_on(0, ymax, 1)
        def rel_poly(tp, sup1):
            dy = sup1 * ymax * ry                              # y の誤差による g の変化（絶対）
            return tp.eps * (1 + dy / tp.gmin) + (tp.e_est + dy) / tp.gmin
        rs, rc = rel_poly(T["sin"], sup_s1), rel_poly(T["cos"], sup_c1)
        rsr = (1 + ru) * (1 + rs) * (1 + tau) - 1              # sin_r = trunc(u·s)
        B["sin"] = B["cos"] = max(rsr, rc)                     # 象限で混ざるので両方の最大
        e = {k: -dyadic_ceil_pow2(v) for k, v in B.items()}
        return e, B

    def claims(self):
        """関数ごとの 主張（ulp 単位）。P は 真値 f から 相対 2^-e 以内（budget ≤ 2^-e）:
             near: v = P                → |v−f| ≤ 2^-e·f            → ≤ 0.5 + 2^(Wout−e)
             lo:   v = P(1 − 2^-e)      → 0 ≤ f−v ≤ (2·2^-e − 2^-2e)·f → ≤ 1 + 2^(Wout−e+1)
             hi:   v = P(1 + 2^(1−e))   → 0 ≤ v−f ≤ (3·2^-e + 2^(1−2e))·f → ≤ 1 + 3·2^(Wout−e)
           （hi は lo より 2^-e ぶん広い: 2^(1−e) の幅出しが P の誤差を打ち消すのに要るため。
             f32 全数掃引で expm1/hi が 1.00101 ulp に達し、旧主張 1 + 2^(Wout−e+1) = 1.000977 を超えて発覚。）
           side_ulp は 互換のため hi 側（大きい方）を返す。"""
        out = {}
        for k, e in self.e.items():
            lo = 1 + 2.0 ** (self.Wout - e + 1); hi = 1 + 3 * 2.0 ** (self.Wout - e)
            out[k] = dict(e=e, near_ulp=0.5 + 2.0 ** (self.Wout - e), side_ulp=hi, side_ulp_lo=lo, side_ulp_hi=hi,
                          rel_bound=float(self.budget[k]))
        return out

    # ------------------------------------------------------------------ 入力の検査とフラグ伝播
    def _check(self, x):
        N, E = x
        assert abs(N) < (1 << self.Win), "入力は W_in 桁まで"
        return dy_norm(N, E)

    def _propagate(self, fn, x, fin, fc):
        gi, li, si = fin
        if (gi, li, si) == (0, 0, 0): return fc
        cg, cl, cs = fc
        xs = (x[0] > 0) - (x[0] < 0)
        if si or fn in ("sin", "cos") or (fn == "log" and gi and li) or xs == 0:
            return (1, 1, 0) if fn == "exp" else (1, 1, 1)   # exp > 0: 符号は沈まない（値の界だけ失う）
        if gi and li:                                          # 入力に界なし → 向きが無い（符号は分かる）
            return (1, 1, 0)
        if fn == "exp":
            dg = (gi and xs > 0) or (li and xs < 0); dl = (li and xs > 0) or (gi and xs < 0)
        elif fn in ("expm1", "sqrt"):
            dg, dl = gi, li
        elif fn == "rsqrt":
            dg, dl = li, gi
        elif fn == "log":
            one = dy_to_fr(x) - 1
            if one > 0:
                if li: return (1, 1, 1)
                dg, dl = gi, 0
            elif one < 0:
                if gi: return (1, 1, 1)
                dg, dl = li, 0
            else:
                dg, dl = 1, 0
        og = int(bool(cl == 0 and dg)); ol = int(bool(cg == 0 and dl))
        if not (og or ol): return (1, 1, cs)
        return (og, ol, 0)

    def _finish(self, fn, P, mode):
        return finish_spec(P, self.e[fn], mode, self.Wout, self.Emin, self.Emax)

    # ------------------------------------------------------------------ exp / expm1
    def _exp_core(self, x):
        N, E = x
        if N == 0: return 0, (0, 0), False
        ld = _lead(N) + E
        if ld <= -3:
            return 0, dy_trunc(x, self.Wk)[0], False           # 小さい経路: r = x
        if ld >= self.XC:
            xf = (1 if N > 0 else -1) << (self.XC + self.Pc)   # clamp
        else:
            xf = fix_trunc(x, self.Pc)
        k = (xf * self.INV_LN2 + (1 << (self.Pc + self.Q - 1))) >> (self.Pc + self.Q)
        r = dy_norm(xf - k * self.LN2, -self.Pc)
        assert abs(dy_to_fr(r)) <= Fr(35, 100)
        return k, dy_trunc(r, self.Wk)[0], True

    def _exp_pair(self, x):
        """x ≠ 0 → (P_exp, P_expm1)（出口前の Wk 桁の値。テープは 1 回だけ評価する）。"""
        k, u, _ = self._exp_core(x)
        q = estrin_spec(self.tapes["exp"].coeffs, u, self.Wk)
        v = dy_add((1, 0), dy_mul(u, q))                       # 1 + u·q  厳密
        Pe = dy_trunc(v, self.Wk)[0]
        Pe = (Pe[0], Pe[1] + k)
        w = dy_add((v[0], v[1] + k), (-1, 0))                  # ·2^k − 1  厳密（ゲートでは sticky で同値）
        Pm = dy_trunc(w, self.Wk)[0]
        return Pe, Pm

    def exp(self, x, flags=(0, 0, 0), mode="near"):
        x = self._check(x)
        if x[0] == 0:                                          # exp(0) = 1 を W_out 桁の正規形で（ゲートと同じ表現）
            return (1 << (self.Wout - 1), -(self.Wout - 1)), self._propagate("exp", x, flags, (0, 0, 0))
        P, _ = self._exp_pair(x)
        out, fc = self._finish("exp", P, mode)
        return out, self._propagate("exp", x, flags, fc)

    def expm1(self, x, flags=(0, 0, 0), mode="near"):
        x = self._check(x)
        if x[0] == 0: return (0, self.Emin), self._propagate("expm1", x, flags, (0, 0, 0))
        _, P = self._exp_pair(x)
        out, fc = self._finish("expm1", P, mode)
        return out, self._propagate("expm1", x, flags, fc)

    # ------------------------------------------------------------------ log
    def _log_P(self, x):
        """x ≠ 0 → P = trunc_Wk(Ex·ln2 + t·q(t))（log|x|）。"""
        N, E = x
        a = abs(N); L = _lead(N); Ex = E + L
        top8 = (a >> (L - 7)) if L >= 7 else (a << (7 - L))
        h = 1 if top8 >= 181 else 0
        t = dy_norm(a - (1 << (L + h)), -(L + h))
        Ex += h
        u = dy_trunc(t, self.Wk)[0]
        assert u == t
        q = estrin_spec(self.tapes["log"].coeffs, u, self.Wk)
        v = dy_add(dy_norm(Ex * self.LN2, -self.Pc), dy_mul(u, q))
        return dy_trunc(v, self.Wk)[0]

    def log(self, x, flags=(0, 0, 0), mode="near"):
        x = self._check(x)
        N, E = x
        if N == 0:
            out = (-self.MAX[0], self.MAX[1])
            return out, self._propagate("log", x, flags, (1, 0, 0))
        P = self._log_P(x)
        out, fc = self._finish("log", P, mode)
        if N < 0: return out, (1, 1, 1)                        # 複素: 実数では何も主張しない
        return out, self._propagate("log", x, flags, fc)

    # ------------------------------------------------------------------ sqrt / rsqrt
    def _root_core(self, x):
        N, E = x
        a = abs(N); L = _lead(N); Ex = E + L
        h = Ex & 1
        t = dy_norm(a - (1 << (L + h)), -(L + h))
        Ex += h
        assert dy_trunc(t, self.Wk)[0] == t
        return t, Ex, h

    def _sqrt_P(self, x):
        t, Ex, h = self._root_core(x)
        v = estrin_spec(self.tapes[f"sqrt_p{h}"].coeffs, t, self.Wk)
        return (v[0], v[1] + Ex // 2)

    def _rsqrt_P(self, x):
        t, Ex, h = self._root_core(x)
        v = estrin_spec(self.tapes[f"rsqrt_p{h}"].coeffs, t, self.Wk)
        return (v[0], v[1] - Ex // 2)

    def sqrt(self, x, flags=(0, 0, 0), mode="near"):
        x = self._check(x)
        if x[0] == 0: return (0, self.Emin), self._propagate("sqrt", x, flags, (0, 0, 0))
        P = self._sqrt_P(x)
        out, fc = self._finish("sqrt", P, mode)
        if x[0] < 0: return out, (1, 1, 1)
        return out, self._propagate("sqrt", x, flags, fc)

    def rsqrt(self, x, flags=(0, 0, 0), mode="near"):
        x = self._check(x)
        if x[0] == 0: return (0, self.Emin), self._propagate("rsqrt", x, flags, (0, 0, 0))   # a/0 = 0
        P = self._rsqrt_P(x)
        out, fc = self._finish("rsqrt", P, mode)
        if x[0] < 0: return out, (1, 1, 1)
        return out, self._propagate("rsqrt", x, flags, fc)

    # ------------------------------------------------------------------ sin / cos
    def _trig_core(self, x):
        """→ (neg, k mod 4, r=(N,E) Wk 桁, large_path)"""
        N, E = x
        neg = N < 0; a = abs(N)
        ld = _lead(N) + E
        if ld <= -2:
            return neg, 0, dy_trunc((a, E), self.Wk)[0], False
        C = window(self.TOP, self.NB, E, self.F)
        T = (a * C) & ((1 << (self.F + 2)) - 1)
        k = (T + (1 << (self.F - 1))) >> self.F
        frac = T - (k << self.F)
        fd = dy_trunc(dy_norm(frac, -self.F), self.Wk + 4)[0]
        r = dy_trunc(dy_mul(fd, (self.PIO2, -self.Pc)), self.Wk)[0]
        return neg, k & 3, r, True

    def _trig(self, x):
        neg, k, r, _ = self._trig_core(x)
        y = dy_trunc(dy_mul(r, r), self.Wk)[0]
        s = estrin_spec(self.tapes["sin"].coeffs, y, self.Wk)
        c = estrin_spec(self.tapes["cos"].coeffs, y, self.Wk)
        sr = dy_trunc(dy_mul(r, s), self.Wk)[0]
        cr = c
        m = lambda v: (-v[0], v[1])
        sin_v = [sr, cr, m(sr), m(cr)][k]
        cos_v = [cr, m(sr), m(cr), sr][k]
        if neg: sin_v = m(sin_v)
        return sin_v, cos_v

    def sin(self, x, flags=(0, 0, 0), mode="near"):
        x = self._check(x)
        if x[0] == 0: return (0, self.Emin), self._propagate("sin", x, flags, (0, 0, 0))
        sv, _ = self._trig(x)
        out, fc = self._finish("sin", sv, mode)
        return out, self._propagate("sin", x, flags, fc)

    def cos(self, x, flags=(0, 0, 0), mode="near"):
        x = self._check(x)
        _, cv = self._trig(x)
        out, fc = self._finish("cos", cv, mode)
        return out, self._propagate("cos", x, flags, fc)

    def fn(self, name, x, flags=(0, 0, 0), mode="near"):
        return getattr(self, name)(x, flags, mode)

    def describe(self):
        print(f"[{self.cfg}] Wout={self.Wout} Wk={self.Wk} Pc={self.Pc} Q={self.Q} XC={self.XC} kmax={self.kmax} "
              f"F={self.F} fbits={self.fbits} ({self.fbits_source}) NB={self.NB}")
        for k in FNS:
            c = self.claims()[k]
            print(f"  {k:6s} e={c['e']:3d}  相対誤差上界 {c['rel_bound']:.3e}  near ≤ {c['near_ulp']:.4f} ulp  "
                  f"lo/hi ≤ {c['side_ulp']:.4f} ulp")


# ====================================================================== 真値（mpmath）と ulp
def truth(name, xf, dps=80):
    mp.dps = dps
    x = mpf(xf.numerator) / xf.denominator
    if name == "exp": return mp.exp(x)
    if name == "expm1": return mp.expm1(x)
    if name == "log": return mp.log(abs(x)) if x != 0 else None
    if name == "sqrt": return mp.sqrt(abs(x))
    if name == "rsqrt": return 1 / mp.sqrt(abs(x)) if x != 0 else mpf(0)
    if name == "sin": return mp.sin(x)
    if name == "cos": return mp.cos(x)


def ulp_of(f, Wout, Emin):
    """真値 f の ulp = 2^(⌊log2|f|⌋ − (Wout−1))、Emin で下限。"""
    if f == 0: return mpf(2) ** Emin
    e = int(mp.floor(mp.log(abs(f), 2)))
    return mpf(2) ** max(e - (Wout - 1), Emin)


def check_one(spec, name, x, mode, flags=(0, 0, 0), strict=True):
    """1 点の検査: 値の ulp 誤差、フラグの正直さ。戻り (ulp_err or None, note)"""
    out, fl = spec.fn(name, x, flags, mode)
    f = truth(name, dy_to_fr(x))
    if f is None: return None, "log0"
    shown = mpf(out[0]) * mpf(2) ** out[1]
    ge, le, sunk = fl
    if sunk: return None, "sunk"
    if (ge, le) == (0, 0):                                     # 厳密の主張
        assert shown == f or abs(shown - f) <= abs(f) * mpf(2) ** -70, (name, x, mode, "exact claim false", shown, f)
        return mpf(0), "exact"
    if ge and not le:
        assert abs(shown) <= abs(f), (name, x, mode, "ge lies", shown, f)
    if le and not ge:
        assert abs(shown) >= abs(f), (name, x, mode, "le lies", shown, f)
    if out == spec.MAX or out == (-spec.MAX[0], spec.MAX[1]) or (abs(out[0]) <= 1 and out[1] == spec.Emin):
        return None, "sat"                                     # 飽和/潰れ（lo の潰れは 0·ge）
    if f != 0 and (shown > 0) != (f > 0): raise AssertionError((name, x, mode, "sign", shown, f))
    err = abs(shown - f) / ulp_of(f, spec.Wout, spec.Emin)
    note = ""
    if strict:
        cl = spec.claims()[name]
        lim = cl["near_ulp"] if mode == "near" else cl["side_ulp"]
        # 'hi' が 2 の冪を上に跨ぐと shown の ulp は f の ulp の 2 倍（外向き丸めの必然: 上界 > 2^j なら shown ≥ 2^j + ulp）
        if mode == "hi" and ulp_of(shown, spec.Wout, spec.Emin) > ulp_of(f, spec.Wout, spec.Emin):
            lim = 2 * lim; note = "cross"
        assert err <= lim + mpf(10) ** -9, (name, x, mode, "ulp claim violated", float(err), lim)
    return err, note


def sample_inputs(spec, name, rng, n):
    """対数一様 + 境界 + 敵対例。"""
    W, Emin, Emax = spec.Wout, spec.Emin, spec.Emax
    xs = []
    for _ in range(n):
        E = rng.randint(Emin, Emax)
        N = rng.getrandbits(W) | (1 << (W - 1))
        if rng.random() < 0.5 and name in ("exp", "expm1", "sin", "cos", "log", "sqrt", "rsqrt"): N = -N
        xs.append((N, E))
    # 境界
    xs += [(1, 0), (-1, 0), ((1 << W) - 1, Emax), (1, Emin), (-(1 << W) + 1, Emax), (-1, Emin),
           (1, -1), (3, -1), (1, 10), (1, -10), ((1 << W) - 1, -W), ((1 << W) - 1, -W + 1), ((1 << (W - 1)) + 1, -(W - 1))]
    mp.dps = 80
    if name in ("exp", "expm1"):
        for k in (1, 2, 3, 7, 10, 100, 700, 1000):
            for sgn in (1, -1):
                v = sgn * k * mp.log(2)
                N = int(mp.nint(v * mpf(2) ** (W - 11)))
                xs += [(N, -(W - 11)), (N + 1, -(W - 11)), (N - 1, -(W - 11))]
        xs += [(N, E) for N in (1, 3, 5) for E in (-W, -W + 1, -2 * W, -3, -4)]
        xs += [(1, spec.XC), (1, spec.XC - 1), (3, spec.XC - 1), (-1, spec.XC), (-3, spec.XC - 1)]
    if name == "log":
        for j in range(1, W + 3):
            xs += [((1 << W) - 1, -W + 1), ((1 << (W - 1)) + (1 << max(0, W - 1 - j)), -(W - 1)),
                   ((1 << W) - (1 << max(0, W - 1 - j)), -W)]
        xs += [(181, -7), (181 * 2 + 1, -8), (181 * 2 - 1, -8), (1, 0), (1, 1), (3, 0), (1, Emax), (1, Emin)]
    if name in ("sqrt", "rsqrt"):
        xs += [(1, 2), (1, 3), (1, -2), (1, -3), (3, 0), (3, 1), (1, Emax), (1, Emin), (1, Emin + 1), ((1 << W) - 1, Emax - 1)]
    if name in ("sin", "cos"):
        for k in (1, 2, 3, 4, 5, 6, 100, 1001, 12345):
            v = k * mp.pi / 2
            for sh in (W - 5, W - 20, 5):
                N = int(mp.nint(v * mpf(2) ** sh))
                if N < (1 << W): xs += [(N, -sh), (N + 1, -sh), (N - 1, -sh)]
        if W == 53: xs += [(6381956970095103, 797), (-6381956970095103, 797)]
        xs += [(1, E) for E in (-1, -2, 0, 1, 2, 5, 10, 60, 100, 500, 900)] + [(1, E) for E in (-W, -3 * W)]
    return [x for x in xs if abs(x[0]) < (1 << W) and Emin <= x[1] <= Emax]


def self_test(cfgs=("f32", "f64"), n=400, seed=1):
    import random, time
    for cfg in cfgs:
        spec = Spec(cfg)
        spec.describe()
        rng = random.Random(seed)
        for name in FNS:
            t0 = time.time()
            worst = {m: mpf(0) for m in ("lo", "hi", "near", "hi-cross")}
            cnt = 0
            for x in sample_inputs(spec, name, rng, n):
                for mode in ("lo", "hi", "near"):
                    err, note = check_one(spec, name, x, mode)
                    if err is not None:
                        key = "hi-cross" if note == "cross" else mode
                        worst[key] = max(worst[key], err)
                    cnt += 1
            print(f"  {name:6s} {cnt:5d} 検査  最悪 ulp: lo {float(worst['lo']):.4f} hi {float(worst['hi']):.4f} "
                  f"(2^j 跨ぎ {float(worst['hi-cross']):.4f}) near {float(worst['near']):.4f}   ({time.time() - t0:.1f}s)")
        flag_test(spec)
    print("funcs_spec self-test OK")


def flag_test(spec, n=300, seed=2):
    """入力フラグの伝播: shown x と フラグ から許される true x の集合を作り、出力の主張を反証しにいく。"""
    import random
    rng = random.Random(seed)
    W = spec.Wout
    bad = 0; total = 0
    for name in FNS:
        for _ in range(n):
            E = rng.randint(-40, 40) if name != "exp" else rng.randint(-40, 6)
            N = rng.getrandbits(W) | (1 << (W - 1))
            if rng.random() < 0.5: N = -N
            x = (N, E)
            fin = rng.choice([(1, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1)])
            out, fo = spec.fn(name, x, fin, "near")
            shown = mpf(out[0]) * mpf(2) ** out[1]
            xf = dy_to_fr(x)
            # true x の候補: 大きさを ×(1±δ) 変えたもの、sunk なら符号も反転
            cands = []
            for d in (Fr(1, 1 << 20), Fr(1, 4), Fr(3, 2), Fr(1, 1)):
                if fin[0] and not fin[1]: cands.append(xf * (1 + d))
                if fin[1] and not fin[0]: cands.append(xf / (1 + d))
                if fin[0] and fin[1]: cands += [xf * (1 + d), xf / (1 + d)]
            if fin[2]: cands += [-c for c in cands]
            for c in cands:
                if name in ("log", "sqrt", "rsqrt") and c <= 0: continue
                ft = truth(name, c)
                if ft is None: continue
                total += 1
                ge, le, sunk = fo
                if not sunk and ft != 0 and (ft > 0) != (shown > 0) and shown != 0:
                    bad += 1; print("   sign lie", name, x, fin, fo, float(shown), ft)
                if ge and not le and abs(shown) > abs(ft) * (1 + mpf(2) ** -60):
                    bad += 1; print("   ge lie", name, x, fin, fo, float(shown), ft)
                if le and not ge and abs(shown) < abs(ft) * (1 - mpf(2) ** -60):
                    bad += 1; print("   le lie", name, x, fin, fo, float(shown), ft)
                if (ge, le) == (0, 0) and abs(shown - ft) > abs(ft) * mpf(2) ** -60:
                    bad += 1; print("   exact lie", name, x, fin, fo, float(shown), ft)
    assert bad == 0, f"フラグ伝播の嘘 {bad}/{total}"
    print(f"  flag propagation: {total} true-candidates, lies 0")


if __name__ == "__main__":
    cfgs = sys.argv[1:] or ("f32", "f64")
    self_test(cfgs)
