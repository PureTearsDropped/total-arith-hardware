#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証（research ブランチ・未監査）
"""tapes — 関数ごとの係数テープを設計して tapes/*.json に保存する（オフライン・一度きり）。

  縮小後の核 と 区間（funcs_spec.py の縮小配線と一対一）:
    exp    : exp(r) = 1 + r·q(r),  q ≈ expm1(r)/r,  r ∈ [−0.35, 0.35]  (|r| ≤ ln2/2 + 余裕)
    log    : log1p(t) = t·q(t),    q ≈ log1p(t)/t,  t ∈ [181/256−1, 181/128−1]  (m の閾値 τ=181/128)
    sqrt   : √(1+t)   2 片:  p0 t ∈ [0,1) (E 偶),  p1 t ∈ [−1/2,0) (E 奇 → m/2, E+1)
    rsqrt  : 1/√(1+t) 同上
    sin    : sin(r) = r·s(y),  s ≈ sin(√y)/√y,  y = r² ∈ [0, 0.62]  (|r| ≤ π/4 + 余裕)
    cos    : cos(r) = c(y),    c ≈ cos(√y)

  構成 (config): f64 = 出力 53 桁・核 Wk=64・格子 P=64・ε 目標 2^-60
                 f32 = 出力 24 桁・核 Wk=40・格子 P=40・ε 目標 2^-36
  次数は 目標を満たす最小を探索（同じ関数の 2 片は 次数を揃える = 配線が同じ）。
  各テープに 証明つき ε_rig, |g| の下界 gmin, 各 Wk での Estrin 切り捨て誤差の上界 を記録。
"""
import sys, os, json, time
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from remez import (SERIES, remez, max_err, quantize_refit, supnorm_exact, make_tape, save_tape,
                   dyadic_ceil_pow2, _mp, mp, mpf)
from tape_eval import estrin_error_bound, plan_stats

HERE = os.path.dirname(os.path.abspath(__file__))
TAPE_DIR = os.path.join(HERE, "tapes")

CONFIGS = {
    "f64": dict(Wout=53, Wk=64, P=64, eps_target=Fr(1, 1 << 60), Emin=-1074, Emax=971),
    "f32": dict(Wout=24, Wk=40, P=40, eps_target=Fr(1, 1 << 36), Emin=-149, Emax=104),
}

KERNELS = {
    "exp":      dict(series="expm1_over_x", a=Fr(-7, 20), b=Fr(7, 20), form="one_plus_x_times", arg="x"),
    "log":      dict(series="log1p_over_x", a=Fr(181, 256) - 1, b=Fr(181, 128) - 1, form="x_times", arg="x"),
    "sqrt_p0":  dict(series="sqrt1p",  a=Fr(0), b=Fr(1), form="poly", arg="x", group="sqrt"),
    "sqrt_p1":  dict(series="sqrt1p",  a=Fr(-1, 2), b=Fr(0), form="poly", arg="x", group="sqrt"),
    "rsqrt_p0": dict(series="rsqrt1p", a=Fr(0), b=Fr(1), form="poly", arg="x", group="rsqrt"),
    "rsqrt_p1": dict(series="rsqrt1p", a=Fr(-1, 2), b=Fr(0), form="poly", arg="x", group="rsqrt"),
    "sin":      dict(series="sinc_sqrt", a=Fr(0), b=Fr(31, 50), form="x_times", arg="x2"),
    "cos":      dict(series="cos_sqrt",  a=Fr(0), b=Fr(31, 50), form="poly", arg="x2"),
}


def find_degree(kname, cfg, n_min=2, n_max=40, verbose=True):
    """量子化前の Remez ε が 目標/4 以下になる最小次数（速い予備探索）。"""
    K = KERNELS[kname]; ser = SERIES[K["series"]]
    target = cfg["eps_target"]
    for n in range(n_min, n_max + 1):
        coeffs, E, _ = remez(ser.mpf_fn, K["a"], K["b"], list(range(n + 1)), "rel")
        em = max_err(ser.mpf_fn, K["a"], K["b"], list(range(n + 1)), coeffs, "rel")
        if verbose: print(f"    {kname}: n={n} ε_remez={mp.nstr(em, 3)}", flush=True)
        if em <= _mp(target) / 4: return n
    raise RuntimeError("degree search failed")


def build_tape(kname, cfgname, n, verbose=True):
    K = KERNELS[kname]; cfg = CONFIGS[cfgname]
    name = f"{kname}_{cfgname}"
    for attempt in range(3):
        tape = make_tape(name, K["series"], K["a"], K["b"], list(range(n + 1)), cfg["P"], "rel",
                         form=K["form"], arg=K["arg"], verbose=verbose, degree=n, config=cfgname,
                         Wk=cfg["Wk"], Wout=cfg["Wout"])
        if Fr(tape["eps_rigorous"]) <= cfg["eps_target"]: break
        if verbose: print(f"    {name}: ε_rig {tape['eps_rigorous_float']:.3g} > 目標 → 次数 +1", flush=True)
        n += 1
    else:
        raise RuntimeError(f"{name}: ε 目標に届かない")
    # Estrin 切り捨て誤差（絶対）と |値| の上界、|g| の下界 → 相対
    R = max(abs(Fr(K["a"])), abs(Fr(K["b"])))
    eb, vb = estrin_error_bound(tape["coeffs"], R, cfg["Wk"])
    gmin = Fr(tape["gmin_rigorous"])
    tape["estrin_abs_error_bound"] = str(eb); tape["estrin_value_bound"] = str(vb)
    tape["estrin_rel_error_bound_float"] = float(eb / gmin)
    sig, const, depth = plan_stats(n)
    tape["estrin_stats"] = dict(signal_mults=sig, const_mults=const, mac_depth=depth)
    return tape


def design_one(args):
    kname, cfgname = args
    t0 = time.time()
    cfg = CONFIGS[cfgname]
    n = find_degree(kname, cfg, verbose=False)
    tape = build_tape(kname, cfgname, n, verbose=False)
    tape["design_seconds"] = round(time.time() - t0, 1)
    return kname, cfgname, tape


def design_all(names=None, configs=None, workers=8):
    import multiprocessing
    os.makedirs(TAPE_DIR, exist_ok=True)
    jobs = [(k, c) for c in (configs or CONFIGS) for k in (names or KERNELS)]
    with multiprocessing.get_context("fork").Pool(workers) as pool:      # forkserver は sandbox で socket 不可
        results = pool.map(design_one, jobs)
    tapes = {}
    for kname, cfgname, tape in results:
        tapes[(kname, cfgname)] = tape
    # 同じグループ（sqrt_p0/p1 …）は 次数を揃える: 高い方に合わせて作り直し
    for cfgname in (configs or CONFIGS):
        groups = {}
        for kname in (names or KERNELS):
            g = KERNELS[kname].get("group")
            if g: groups.setdefault(g, []).append(kname)
        for g, members in groups.items():
            nmax = max(tapes[(k, cfgname)]["degree"] for k in members)
            for k in members:
                if tapes[(k, cfgname)]["degree"] < nmax:
                    print(f"  {k}_{cfgname}: 次数 {tapes[(k, cfgname)]['degree']} → {nmax} に揃える", flush=True)
                    tapes[(k, cfgname)] = build_tape(k, cfgname, nmax, verbose=False)
    for (kname, cfgname), tape in tapes.items():
        save_tape(tape, os.path.join(TAPE_DIR, f"{kname}_{cfgname}.json"))
    return tapes


def report(tapes=None):
    if tapes is None:
        tapes = {}
        for fn in sorted(os.listdir(TAPE_DIR)):
            if fn.endswith(".json"):
                t = json.load(open(os.path.join(TAPE_DIR, fn)))
                tapes[(t["name"].rsplit("_", 1)[0], t["config"])] = t
    print(f"{'tape':<14}{'次数':>4} {'ε_remez':>9} {'ε_量子化':>9} {'ε_証明':>9} {'2^-e':>5} {'gmin':>6} "
          f"{'切捨誤差':>9} {'信号積':>5} {'定数積':>5} {'深さ':>3} {'秒':>6}")
    for (k, c), t in sorted(tapes.items(), key=lambda kv: (kv[0][1], list(KERNELS).index(kv[0][0]))):
        s = t["estrin_stats"]
        print(f"{t['name']:<14}{t['degree']:>4} {float(t['eps_remez_unquantized']):>9.2e} "
              f"{float(t['eps_measured']):>9.2e} {t['eps_rigorous_float']:>9.2e} {t['eps_pow2_exponent']:>5} "
              f"{t['gmin_float']:>6.3f} {t['estrin_rel_error_bound_float']:>9.2e} "
              f"{s['signal_mults']:>5} {s['const_mults']:>5} {s['mac_depth']:>3} {t.get('design_seconds', 0):>6.0f}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--config", nargs="*", default=None)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    if a.report:
        report()
    else:
        t0 = time.time()
        tapes = design_all(a.only, a.config, a.workers)
        print(f"合計 {time.time() - t0:.0f} 秒")
        report(tapes)
