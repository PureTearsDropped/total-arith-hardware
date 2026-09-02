#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証（research ブランチ・未監査）
"""fracmin — sin/cos の引数縮小 x = k·(π/2) + r の最悪ケース（|x·2/π − k| の最小）を求める。

  縮小の精度要求は **相対**（sin(x) ≈ r が小さいとき r の相対精度が出力の相対精度）なので、
  窓つき定数（Payne–Hanek）の分数ビット数 F は
      F ≥ W_in + fbits + Wk + 2,   2^-fbits ≤ min_x |frac(x)|  (|x| ≥ 1/2 の全入力)
  で決める。f32 構成は **全入力総当たり**（本スクリプト）。f64 構成は総当たり不能なので
  文献値を引く（binary64 の最悪 x = 6381956970095103·2^797, |r| ≈ 2^-60.9 → |frac| ≈ 2^-61.5;
  Muller, *Elementary Functions* 3rd ed. §11 / Kahan の "Minimizing q*m−n"）。本スクリプトは
  その 1 点を mpmath で追検算するだけで、f64 の全域は **引用**であって検証ではない。

  使い方: python fracmin.py f32   → tapes/fracmin_f32.json（min|frac|, その x, 使った窓幅）
"""
import sys, os, json, time
from fractions import Fraction as Fr
from mpmath import mp, mpf

HERE = os.path.dirname(os.path.abspath(__file__))
TAPE_DIR = os.path.join(HERE, "tapes")

F0 = 200                       # 探索に使う窓の分数ビット（求める最小値 2^-40 程度に対し十分）


def two_over_pi_bits(NB):
    mp.dps = int(NB * 0.302) + 40
    return int(mp.floor(2 / mp.pi * mpf(2) ** NB))


def window(TOP, NB, E, F):
    """floor(2^(E+F)·2/π) mod 2^(F+2)"""
    s = NB - E - F
    v = TOP >> s if s >= 0 else TOP << (-s)
    return v & ((1 << (F + 2)) - 1)


def scan_bin(args):
    """指数ビン E: M ∈ [2^(W−1), 2^W) の全 M について frac の最小絶対値。"""
    W, E, TOP, NB = args
    C = window(TOP, NB, E, F0)
    mask = (1 << (F0 + 2)) - 1
    half = 1 << (F0 - 1)
    best = None
    for M in range(1 << (W - 1), 1 << W):
        T = (M * C) & mask
        k = (T + half) >> F0
        frac = T - (k << F0)
        a = -frac if frac < 0 else frac
        if best is None or a < best[0]:
            best = (a, M, k & 3)
    return E, best


def main(cfg="f32"):
    from tapes import CONFIGS
    c = CONFIGS[cfg]
    W = c["Wout"]                                          # 入力の桁数 = 出力の桁数
    E_lo = -W                                              # 先頭位置 ≥ −1 ⟺ |x| ≥ 1/2 (E + W − 1 ≥ −1)
    E_hi = c["Emax"]
    NB = E_hi + F0 + 8
    TOP = two_over_pi_bits(NB)
    import multiprocessing
    t0 = time.time()
    jobs = [(W, E, TOP, NB) for E in range(E_lo, E_hi + 1)]
    with multiprocessing.get_context("fork").Pool(min(24, len(jobs))) as pool:
        res = pool.map(scan_bin, jobs)
    best = min(res, key=lambda r: r[1][0])
    E, (a, M, q) = best
    frac = Fr(a, 1 << F0)
    fbits = 0
    while Fr(1, 1 << (fbits + 1)) >= frac: fbits += 1      # 2^-fbits ≤ frac < 2^-(fbits-1)… 保守側: 2^-(fbits+1) < frac
    fbits += 1                                             # 2^-fbits ≤ frac を保証
    assert Fr(1, 1 << fbits) <= frac
    # 追検算: mpmath で x·2/π の小数部
    mp.dps = 120
    x = mpf(M) * mpf(2) ** E
    t = x * 2 / mp.pi
    fr = t - mp.nint(t)
    out = dict(config=cfg, W_in=W, E_range=[E_lo, E_hi], F_search=F0,
               fracmin=str(frac), fracmin_float=float(frac), fbits=fbits,
               worst_M=M, worst_E=E, worst_x=mp.nstr(x, 12), worst_quadrant=q,
               mpmath_check=mp.nstr(fr, 12), seconds=round(time.time() - t0, 1),
               note="exhaustive over all inputs with |x| >= 1/2 (leading position >= -1); frac = x*2/pi - round(...)")
    os.makedirs(TAPE_DIR, exist_ok=True)
    with open(os.path.join(TAPE_DIR, f"fracmin_{cfg}.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


def check_f64_citation():
    """引用した binary64 の最悪ケースを 1 点だけ追検算する。"""
    mp.dps = 400                                            # x ≈ 2^850 なので小数部に 900 bit 以上要る
    M, E = 6381956970095103, 797
    x = mpf(M) * mpf(2) ** E
    t = x * 2 / mp.pi
    fr = t - mp.nint(t)
    print(f"binary64 cited worst case x = {M}·2^{E}: frac = {mp.nstr(fr, 8)} = 2^{mp.nstr(mp.log(abs(fr), 2), 5)}"
          f"  (|r| = {mp.nstr(abs(fr) * mp.pi / 2, 8)})")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "f64":
        check_f64_citation()
    else:
        main("f32")
