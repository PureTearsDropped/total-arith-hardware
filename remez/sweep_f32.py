#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証（research ブランチ・未監査）
"""sweep_f32 — Level 1 spec (funcs_spec.Spec, f32 構成) の **float32 幅 全入力総当たり**。

  双子（numpy 模倣）は作らない: spec 自体が 1 点 15 µs なので、24 コアで全 2^32 入力を数時間で回せる。
  検査対象はしたがって「機械に載る spec そのもの」で、近似の双子ではない。

  真値: numpy の float64 関数（glibc; 誤差 < 1 ulp_f64 ≈ 2^-52 相対）。f32 の ulp は 2^-23 相対なので
  ulp 誤差の分解能は 2^-29 ulp。片側境界（lo ≤ f ≤ hi）の判定は 2^-48 の許容をつけ、許容内で際どい例と
  違反候補は最後に mpmath（80 桁）で再検して確定する。

  入力集合（f32 と同じ語）: N ∈ [2^23, 2^24), E ∈ [−149, 104]（正規）+ N ∈ [1, 2^23), E = −149（非正規）, ±。
    exp/expm1 : 両符号（別経路）。 log/sqrt/rsqrt : 正のみ（負は |x| を 111 で返すだけ、値は同じ経路）。
    sin/cos   : 正のみ（負は符号の鏡: _trig_core は |N| で縮小し sin だけ符号を反転）。x=0 は別に検査。

  各点で検査するもの（mode ∈ lo/hi/near）:
    ・ ge だけ → |shown| ≤ |f|,  le だけ → |shown| ≥ |f|（片側境界の正直さ）
    ・ 00（厳密の主張）→ shown == f
    ・ 符号一致、飽和 ±MAX・潰れ ±MIN/0 の記録
    ・ ulp 誤差 ≤ 主張（near: 0.5+2^(Wout−e), lo/hi: 1+2^(Wout−e+1), hi が 2^j を跨ぐと 2 倍）
    ・ 統計: 最悪 ulp と その x、片側境界の最小余裕（f と shown の距離 / ulp）、near が最近接でない件数

  使い方: python -u sweep_f32.py [group ...] [--quick]     group ∈ exp log sqrt rsqrt trig
          結果は sweep/f32_<group>.json、進捗は stdout。
"""
import os, sys, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from funcs_spec import Spec, truth, ulp_of, FNS
from tape_eval import finish_spec, dy_to_fr
from mpmath import mp, mpf

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "sweep")
MODES = ("lo", "hi", "near")
GROUPS = {"exp": ("exp", "expm1"), "log": ("log",), "sqrt": ("sqrt",), "rsqrt": ("rsqrt",), "trig": ("sin", "cos")}
SIGNS = {"exp": (1, -1), "log": (1,), "sqrt": (1,), "rsqrt": (1,), "trig": (1,)}
TOL = 2.0 ** -48                       # 真値（float64）の誤差に対する許容
CHUNK = 1 << 21

SPEC = None                            # fork で継承（各ワーカーで再構築しない）


def _truth_np(fn, xv):
    with np.errstate(all="ignore"):
        if fn == "exp": return np.exp(xv)
        if fn == "expm1": return np.expm1(xv)
        if fn == "log": return np.log(xv)
        if fn == "sqrt": return np.sqrt(xv)
        if fn == "rsqrt": return 1.0 / np.sqrt(xv)
        if fn == "sin": return np.sin(xv)
        if fn == "cos": return np.cos(xv)


def _cores(group):
    s = SPEC
    if group == "exp": return lambda x: s._exp_pair(x)
    if group == "log": return lambda x: (s._log_P(x),)
    if group == "sqrt": return lambda x: (s._sqrt_P(x),)
    if group == "rsqrt": return lambda x: (s._rsqrt_P(x),)
    if group == "trig": return lambda x: s._trig(x)


def ulp32(f, Wout, Emin):
    if f == 0.0 or math.isinf(f): return 0.0
    _, ex = math.frexp(abs(f))           # |f| = m·2^ex, m ∈ [0.5,1)
    return math.ldexp(1.0, max(ex - 1 - (Wout - 1), Emin))


class Stat:
    __slots__ = ("n", "exact", "sat", "collapse", "max_err", "max_cross", "min_margin", "n_not_cr", "viol", "marg")

    def __init__(self):
        self.n = 0; self.exact = 0; self.sat = 0; self.collapse = 0
        self.max_err = (-1.0, None); self.max_cross = (-1.0, None); self.min_margin = (math.inf, None)
        self.n_not_cr = 0; self.viol = []; self.marg = []

    def dump(self):
        return dict(n=self.n, exact=self.exact, sat=self.sat, collapse=self.collapse,
                    max_err=self.max_err, max_cross=self.max_cross,
                    min_margin=[self.min_margin[0] if self.min_margin[0] != math.inf else None, self.min_margin[1]],
                    n_not_cr=self.n_not_cr, viol=self.viol, marg=self.marg)


def check(st, fn, mode, x, out, fc, f, lim_near, lim_side, Wout, Emin, MAXN, Emax):
    st.n += 1
    shown = math.ldexp(out[0], out[1])
    ge, le, sunk = fc
    if sunk:
        if len(st.viol) < 20: st.viol.append(("sunk", x, mode, out, fc, f)); return
    if (ge, le) == (0, 0):
        st.exact += 1
        if shown != f and len(st.viol) < 20: st.viol.append(("exact-claim", x, mode, out, fc, f))
        return
    af, ash = abs(f), abs(shown)
    if ge and not le and ash > af:
        if ash > af * (1 + TOL):
            if len(st.viol) < 20: st.viol.append(("ge-lies", x, mode, out, fc, f))
        elif len(st.marg) < 20: st.marg.append(("ge-marginal", x, mode, out, fc, f))
    if le and not ge and ash < af:
        if ash < af * (1 - TOL):
            if len(st.viol) < 20: st.viol.append(("le-lies", x, mode, out, fc, f))
        elif len(st.marg) < 20: st.marg.append(("le-marginal", x, mode, out, fc, f))
    if abs(out[0]) == MAXN and out[1] == Emax:
        st.sat += 1; return
    if abs(out[0]) <= 1 and out[1] == Emin:                    # ±MIN（le）または 0（lo の潰れ, ge）
        st.collapse += 1; return
    if f != 0.0 and shown != 0.0 and (shown > 0) != (f > 0):
        if len(st.viol) < 20: st.viol.append(("sign", x, mode, out, fc, f))
        return
    u = ulp32(f, Wout, Emin)
    if u == 0.0:                                               # f = 0 / inf なのに 飽和も潰れもしていない
        if len(st.viol) < 20: st.viol.append(("no-sat", x, mode, out, fc, f))
        return
    err = abs(shown - f) / u
    if mode == "near":
        lim = lim_near
        if err > 0.5 + 2.0 ** -24: st.n_not_cr += 1
    else:
        lim = lim_side
        m = (af - ash) / u if mode == "lo" else (ash - af) / u
        if m < st.min_margin[0]: st.min_margin = (m, x)
        if mode == "hi" and ulp32(shown, Wout, Emin) > u:
            lim = 2 * lim
            if err > st.max_cross[0]: st.max_cross = (err, x)
            if err > lim + 2.0 ** -20 and len(st.viol) < 20: st.viol.append(("ulp", x, mode, out, fc, f, err, lim))
            return
    if err > st.max_err[0]: st.max_err = (err, x)
    if err > lim + 2.0 ** -20 and len(st.viol) < 20: st.viol.append(("ulp", x, mode, out, fc, f, err, lim))


def run_job(job):
    group, sign, E, N0, N1 = job
    s = SPEC
    fns = GROUPS[group]
    core = _cores(group)
    Ns = np.arange(N0, N1, dtype=np.int64)
    xv = np.ldexp(Ns.astype(np.float64) * sign, E)
    tr = {fn: _truth_np(fn, xv).tolist() for fn in fns}
    cl = s.claims()
    lim = {fn: (cl[fn]["near_ulp"], cl[fn]["side_ulp_lo"], cl[fn]["side_ulp_hi"]) for fn in fns}
    Wout, Emin, Emax, MAXN = s.Wout, s.Emin, s.Emax, s.MAX[0]
    st = {(fn, mode): Stat() for fn in fns for mode in MODES}
    for i in range(N1 - N0):
        x = (sign * (N0 + i), E)
        Ps = core(x)
        for fn, P in zip(fns, Ps):
            f = tr[fn][i]
            ln, llo, lhi = lim[fn]
            for mode in MODES:
                out, fc = finish_spec(P, s.e[fn], mode, Wout, Emin, Emax)
                check(st[(fn, mode)], fn, mode, x, out, fc, f, ln, lhi if mode == "hi" else llo, Wout, Emin, MAXN, Emax)
    return job, {f"{fn}/{mode}": v.dump() for (fn, mode), v in st.items()}


def merge(acc, part):
    for k, d in part.items():
        a = acc.setdefault(k, dict(n=0, exact=0, sat=0, collapse=0, max_err=(-1.0, None), max_cross=(-1.0, None),
                                   min_margin=(math.inf, None), n_not_cr=0, viol=[], marg=[]))
        for f in ("n", "exact", "sat", "collapse", "n_not_cr"): a[f] += d[f]
        for f in ("max_err", "max_cross"):
            if d[f][0] > a[f][0]: a[f] = tuple(d[f])
        if d["min_margin"][0] is not None and d["min_margin"][0] < a["min_margin"][0]: a["min_margin"] = tuple(d["min_margin"])
        a["viol"] = (a["viol"] + d["viol"])[:200]
        a["marg"] = (a["marg"] + d["marg"])[:200]
    return acc


def jobs_for(group, quick=False):
    s = SPEC
    W, Emin, Emax = s.Win, s.Emin, s.Emax
    Es = list(range(Emin, Emax + 1)) if not quick else [-149, -30, -3, -2, -1, 0, 1, 2, 6, 7, 104]
    out = []
    for sign in SIGNS[group]:
        for E in Es:
            for N0 in range(1 << (W - 1), 1 << W, CHUNK):
                out.append((group, sign, E, N0, N0 + CHUNK))
        for N0 in range(1, 1 << (W - 1), CHUNK):                       # 非正規
            out.append((group, sign, Emin, N0, min(N0 + CHUNK, 1 << (W - 1))))
    return out


def recheck(acc, cfg="f32"):
    """違反候補と際どい例を mpmath で確定する。"""
    s = SPEC
    final = {}
    for k, a in acc.items():
        fn, mode = k.split("/")
        conf = []
        for item in a["viol"] + a["marg"]:
            kind, x, mode_, out, fc, f = item[:6]
            x = (int(x[0]), int(x[1])); out = (int(out[0]), int(out[1]))
            ft = truth(fn, dy_to_fr(x))
            shown = mpf(out[0]) * mpf(2) ** out[1]
            ge, le, _ = fc
            bad = None
            if kind == "exact-claim" and shown != ft: bad = "exact-claim"
            elif kind.startswith("ge") and ge and not le and abs(shown) > abs(ft): bad = "ge-lies"
            elif kind.startswith("le") and le and not ge and abs(shown) < abs(ft): bad = "le-lies"
            elif kind == "ulp":
                err = abs(shown - ft) / ulp_of(ft, s.Wout, s.Emin)
                if err > item[7] + mpf(10) ** -9: bad = f"ulp {float(err):.6f} > {item[7]}"
            elif kind in ("sign", "no-sat", "sunk"): bad = kind
            if bad: conf.append((bad, x, mode, out, fc, mp.nstr(ft, 20)))
        final[k] = conf
    return final


def main():
    global SPEC
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    quick = "--quick" in sys.argv
    groups = args or list(GROUPS)
    SPEC = Spec("f32")
    SPEC.describe()
    os.makedirs(OUT_DIR, exist_ok=True)
    import multiprocessing
    nproc = min(24, os.cpu_count() or 1)
    # x = 0 と 負の log/sqrt/rsqrt（111）を別に
    for fn in FNS:
        for mode in MODES:
            out, fc = SPEC.fn(fn, (0, 0), (0, 0, 0), mode)
            print(f"  x=0 {fn:6s} {mode:4s} → {out} {fc}")
    for fn in ("log", "sqrt", "rsqrt"):
        out, fc = SPEC.fn(fn, (-(1 << 23) - 12345, -20), (0, 0, 0), "near")
        assert fc == (1, 1, 1), (fn, fc)
    print(f"  負の log/sqrt/rsqrt: 111 ✓")
    for group in groups:
        jobs = jobs_for(group, quick)
        t0 = time.time(); acc = {}; done = 0
        print(f"[{group}] {len(jobs)} jobs × {CHUNK} 点, {nproc} workers")
        with multiprocessing.get_context("fork").Pool(nproc) as pool:
            for job, part in pool.imap_unordered(run_job, jobs, chunksize=1):
                acc = merge(acc, part); done += 1
                if done % 50 == 0 or done == len(jobs):
                    el = time.time() - t0
                    print(f"  {done}/{len(jobs)}  {el:.0f}s  残り≈{el / done * (len(jobs) - done):.0f}s", flush=True)
        final = recheck(acc)
        total_bad = sum(len(v) for v in final.values())
        print(f"[{group}] 完了 {time.time() - t0:.0f}s  確定違反 {total_bad}")
        for k in sorted(acc):
            a = acc[k]
            mm = a["min_margin"]
            print(f"  {k:12s} n={a['n']:>11d} exact={a['exact']} sat={a['sat']} collapse={a['collapse']} "
                  f"最悪ulp={a['max_err'][0]:.6f}@{a['max_err'][1]} 跨ぎ={a['max_cross'][0]:.6f}@{a['max_cross'][1]} "
                  f"最小余裕={mm[0] if mm[0] == math.inf else f'{mm[0]:.3e}'}@{mm[1]} "
                  f"非最近接={a['n_not_cr']} 違反候補={len(a['viol'])} 際どい={len(a['marg'])} 確定違反={len(final[k])}")
            for c in final[k][:10]: print("     ✗", c)
        res = dict(config="f32", group=group, quick=quick, seconds=round(time.time() - t0, 1), jobs=len(jobs),
                   claims={fn: SPEC.claims()[fn] for fn in GROUPS[group]},
                   stats={k: dict(a, viol=a["viol"][:20], marg=a["marg"][:20]) for k, a in acc.items()},
                   confirmed_violations=final, truth="numpy float64 (glibc), rechecked with mpmath dps=80")
        with open(os.path.join(OUT_DIR, f"f32_{group}{'_quick' if quick else ''}.json"), "w") as f:
            json.dump(res, f, indent=1, default=str)


if __name__ == "__main__":
    main()
