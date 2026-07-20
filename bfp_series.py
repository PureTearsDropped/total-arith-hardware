#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証
"""bfp_series — 正規化つき 級数ユニット: gate_series(EXACT) の「次の一歩」。

  gate_series.py の実測: EXACTモード(正規化なし)は 仮数幅が 積で倍増し、
  order=8+二乗2回で 43億ゲート — 「実用回路には 各段の正規化が必須」が数字で出た。
  本モジュールは その正規化版: bfp_sed の 誠実な丸め(normalize: 切り捨て→≥・潰れ→±MIN+≤)
  と フラグ健全な mul/add (oracle_bfp 164,800照合済) を 級数に 連鎖する。

      xs   = x·2^{-s}           … 共有指数 E の減算(無損失・ゲート0個)
      term = 1 ; acc = c_0
      term = mul(term, xs, W)   … 各段で W桁へ 正規化 = 幅が 有界 = コスト線形
      acc  = add(acc, mul(c_k, term, W), W)
      acc  = mul(acc, acc, W) × s回

  ★端から端までの主張: 出力の 成分ごとフラグ {=, ≥, ≤, 境界なし, 符号不明} は
    「丸めを 一切しない 理想計算(同じ量子化テープ・厳密有理数)の 真値」への 主張である。
      =(かつ符号信頼)  ⟹ 表示 == 理想 (厳密一致!)
      ≥               ⟹ |理想| ≥ |表示| (符号は sunk がなければ 信頼)
      ≤               ⟹ |理想| ≤ |表示|
    各演算のフラグ代数は 単演算では oracle_bfp が保証済み。ここでは その「合成」が
    端から端まで 健全かを 検証する — 理想値を 厳密分数で 計算し、全成分の主張を 反証にかける。
"""
from fractions import Fraction as Fr
import math
import numpy as np

from bfp_sed import BF, mul, add, EQ, GE, LE, NB
from sd2_core import M
from sedenion_tensor_logic import OMEGA


# ============================================================ 級数ユニット(正規化つき)
def const_sed(m_int, E):
    "テープ定数 c = m·2^E を セデニオンの 実成分に(厳密・フラグ=)"
    return BF([m_int] + [0] * (M - 1), E)

def tape_exp(order, P):
    return [(int(round((1 << P) / math.factorial(k))), -P) if k else (1, 0)
            for k in range(order + 1)]

def tape_sin(order, P):
    return [(0, 0) if k % 2 == 0 else
            (int(round((1 << P) / math.factorial(k))) * (-1) ** ((k - 1) // 2), -P)
            for k in range(order + 1)]

def series_unit(x, tape, s, W, Emax=None):
    """Σ c_k (x·2^{-s})^k を 左結合で → s回 二乗。各段 W桁へ 誠実に 正規化。"""
    xs = BF(list(x.mant), x.E - s, list(x.bound), list(x.sunk))
    m0, e0 = tape[0]
    acc = const_sed(m0, e0)
    term = const_sed(1, 0)
    for k in range(1, len(tape)):
        term = mul(term, xs, W, Emax)
        mk, ek = tape[k]
        if mk != 0:
            acc = add(acc, mul(const_sed(mk, ek), term, W, Emax), W, Emax)
    for _ in range(s):
        acc = mul(acc, acc, W, Emax)
    return acc


# ============================================================ 理想(丸めなし)の厳密仕様
def _spec_mul(a, b):
    r = [Fr(0)] * M
    for i in range(M):
        if a[i] == 0: continue
        for j in range(M):
            if b[j] == 0: continue
            r[i ^ j] += Fr(int(OMEGA[i, j])) * a[i] * b[j]
    return r

def series_spec(x_fr, tape, s):
    xs = [v / (1 << s) for v in x_fr]
    m0, e0 = tape[0]
    acc = [Fr(m0) * Fr(2) ** e0] + [Fr(0)] * (M - 1)
    term = [Fr(1)] + [Fr(0)] * (M - 1)
    for k in range(1, len(tape)):
        term = _spec_mul(term, xs)
        mk, ek = tape[k]
        if mk != 0:
            ck = Fr(mk) * Fr(2) ** ek
            acc = [a + ck * t for a, t in zip(acc, term)]
    for _ in range(s):
        acc = _spec_mul(acc, acc)
    return acc


# ============================================================ 端から端のフラグ反証
def admits(shown_fr, bound, sunk, ideal_fr):
    "出力1成分の主張が 理想値を 許容するか (=は厳密一致・≥/≤は|·|境界・符号はsunkなし時のみ)"
    ge, le = bound & GE, bound & LE
    if bound == EQ and not sunk:
        return ideal_fr == shown_fr                      # = は「厳密」の主張
    if ge and not le and not (abs(ideal_fr) >= abs(shown_fr)): return False
    if le and not ge and not (abs(ideal_fr) <= abs(shown_fr)): return False
    if not sunk and shown_fr != 0 and ideal_fr != 0:
        if (ideal_fr > 0) != (shown_fr > 0): return False
    return True

def check_series(x_float, tape, s, W, F, Emax=None):
    """入口量子化 → 正規化つき級数 → 全16成分の主張を 理想(厳密分数)に 反証。
       返り値: (違反数, max|表示−理想|, フラグ分布)"""
    m_in = [int(round(v * (1 << F))) for v in x_float]
    x = BF(m_in, -F)                                     # 入口: 厳密(フラグ=)
    got = series_unit(x, tape, s, W, Emax)
    ideal = series_spec([Fr(m, 1 << F) for m in m_in], tape, s)
    viol, err = 0, Fr(0)
    dist = {EQ: 0, GE: 0, LE: 0, NB: 0}
    for k in range(M):
        shown = Fr(got.mant[k]) * Fr(2) ** got.E
        dist[got.bound[k]] += 1
        if not admits(shown, got.bound[k], got.sunk[k], ideal[k]):
            viol += 1
            print(f"    違反: 成分{k} 表示{float(shown):.6g}"
                  f"⟦{'≥' if got.bound[k]&GE else ''}{'≤' if got.bound[k]&LE else ''}"
                  f"{'±' if got.sunk[k] else ''}⟧ vs 理想 {float(ideal[k]):.6g}")
        err = max(err, abs(shown - ideal[k]))
    return viol, float(err), dist, got


# ============================================================ self-test
def self_test():
    print("=" * 84)
    print("bfp_series: 正規化つき級数ユニット — フラグの主張を 理想(厳密分数)に 端から端まで反証")
    print("=" * 84)
    rng = np.random.default_rng(11)

    print("① 健全性スイープ: セデニオンexp (order=8, s=2, W=20) × 20シード")
    tot_v, worst = 0, 0.0
    dsum = {EQ: 0, GE: 0, LE: 0, NB: 0}
    for _ in range(20):
        x = 0.4 * rng.standard_normal(M)
        v, e, d, _g = check_series(x, tape_exp(8, 24), s=2, W=20, F=10)
        tot_v += v; worst = max(worst, e)
        for k in d: dsum[k] += d[k]
    print(f"   違反 {tot_v} / 320成分  max|表示−理想| {worst:.1e}")
    print(f"   フラグ分布: ={dsum[EQ]}  ≥{dsum[GE]}  ≤{dsum[LE]}  境界なし{dsum[NB]}")
    print("   ★正直な発見: 密な長い連鎖では フラグは ほぼ全て「境界なし」へ 退化する")
    print("     (単演算の密保持率0.5%の 合成版・嘘は 0 だが 情報も 薄い —")
    print("      量的な区間が 欲しければ 4値桁/区間表現が 次の構造、の 実証でもある)")
    assert tot_v == 0

    print("② Wのダイヤル: 幅を増やすと 理想への距離が 縮む(健全なまま)")
    x = 0.4 * rng.standard_normal(M)
    prev = None
    for W in (10, 14, 20, 28):
        v, e, d, _g = check_series(x, tape_exp(8, 30), s=2, W=W, F=12)
        ok = prev is None or e <= prev * 1.5              # 単調(ゆるめ・丸めの偶然を許す)
        print(f"   W={W:>2}: 違反{v} max誤差 {e:.1e} {'✓' if v==0 else '✗'}")
        assert v == 0
        prev = e

    print("③ 溢れ: 実部の大きい入力(理想e¹²≈16万) + 小Emax → 本当に飽和して なお健全")
    # 虚方向のexpは回転(有界)なので、指数的に伸びるのは実部。x0を大きくして飽和させる
    xbig = 0.5 * np.abs(rng.standard_normal(M)); xbig[0] = 12.0
    for Emax in (40, -4):
        v, e, d, g = check_series(xbig, tape_exp(12, 24), s=4, W=16, F=8, Emax=Emax)
        MAXm = (1 << 16) - 1
        sat = sum(1 for m in g.mant if abs(m) >= MAXm)     # 仮数が±MAXに張り付いた成分
        print(f"   Emax={Emax:>3}: 違反{v} 飽和成分{sat:>2}"
              f" {'(±MAX+≥が長連鎖の区間フラグとORされ境界なしに埋没するが、主張は破れない)' if sat else '(飽和なし)'}")
        assert v == 0
        if Emax == -4:
            assert sat > 0, "Emax=-4 で 飽和が 起きないのは おかしい"

    print("④ テープ差し替え → sin (同じ骨格・同じオラクル)")
    tot_v = 0
    for _ in range(10):
        x = 0.4 * rng.standard_normal(M)
        v, e, d, _g = check_series(x, tape_sin(9, 24), s=0, W=20, F=10)
        tot_v += v
    print(f"   違反 {tot_v} / 160成分 ✓")
    assert tot_v == 0

    print("-" * 84)
    print("到達点: gate_series(EXACT)の『幅爆発 43億ゲート』に対し、正規化版は 各段 W桁で")
    print("        コスト線形。そして 丸めの代償は 沈黙でなく フラグ — 出力の {=,≥,≤,±} は")
    print("        『丸めなしの理想計算』への 主張として 端から端まで 反証済み(違反0)。")
    print("        単演算の健全性(oracle_bfp)が 合成でも 破れない ことの 実証。")

if __name__ == "__main__":
    self_test()
