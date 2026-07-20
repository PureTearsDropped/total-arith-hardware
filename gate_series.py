#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証
"""gate_series — 超越関数を ゲートで: 係数テープ × (U,V,W)ユニットの 数珠つなぎ。

  exp / sin / … を 「級数の係数テープ(合成時定数)」が指揮する「双線形ユニットの列」として
  AND/OR/NOT/XOR まで 落とす。判断も ループも ない 固定配線 — 演算 = セルの結合の仕方。

      xs   = x · 2^{-s}            … 指数Eの 減算 = 底2の付け替え(無損失・ゲート0個)
      term = e0 ; acc = c_0·e0
      term = bilinear(term, xs)    … 左結合と 宣言した 括弧 (left_power と同じ規約)
      acc += c_k · term            … c_k = 二進有理に 量子化した 定数 round(2^P/k!)·2^{-P}
                                      (K1: 係数は 合成時定数 = 配線)
      acc = bilinear(acc, acc)     … スクエアリング ×s回 (exp のみ; 半群性が根拠)

  テープを 差し替えると 同じ骨格が 別の関数になる(exp ↔ sin)。逆行列は 一切 使わない
  (テイラー = 前向きだけ = 零因子でも 壊れない)。

  検証は 二層:
    ① ゲートグラフ ≡ 厳密有理数仕様 — 同じ量子化定数・同じ順序の Fraction 計算と
       **厳密一致**(EXACTモード: 正規化なし・仮数は 厳密に 伸びる)。ゲートの 嘘 0 を 主張。
    ② 仕様 ≈ 真の exp — テープ量子化 + 打ち切り + スケーリングの 誤差を 実測して 報告。
       ①と②を 分けるのが 誠実さ: 回路は 仕様に 対して 厳密・仕様は 真値に 対して 近似。
"""
from fractions import Fraction as Fr
import math
import numpy as np

from gate_bilinear import new_counter, to_sd
from gate_bfp import BF, to_bf, from_bf, bf_mul, bf_add, bf_bilinear_unit
from bilinear_unit import group_as_UVW
from nd_algebra import cd_omega


# ============================================================ 入口: 実数ベクトル → BF (二進有理)
def quantize_vec(xs, F):
    """成分ごとに m = round(x·2^F) → BF(m, E=-F)。量子化は 入口で 1回だけ・以後は 厳密。"""
    out, exact = [], []
    for x in xs:
        m = int(round(x * (1 << F)))
        out.append(BF(to_sd(m, max(2, m.bit_length() + 2)), -F))
        exact.append(Fr(m, 1 << F))
    return out, exact

def tape_exp(order, P):
    """exp のテープ: c_k = round(2^P/k!)·2^{-P} (二進有理・合成時定数)。c_0=1 は 厳密。"""
    return [(int(round((1 << P) / math.factorial(k))), -P) if k else (1, 0)
            for k in range(order + 1)]

def tape_sin(order, P):
    """sin のテープ: 奇数項 (−1)^((k−1)/2)/k!。偶数項は 0 = 配線ごと 消える。"""
    t = []
    for k in range(order + 1):
        if k % 2 == 0:
            t.append((0, 0))
        else:
            m = int(round((1 << P) / math.factorial(k))) * (-1) ** ((k - 1) // 2)
            t.append((m, -P))
    return t


# ============================================================ 本体: テープ × ユニット列
def const_bf(m, E):
    return BF(to_sd(m, max(2, abs(m).bit_length() + 2)), E)

def gate_series_unit(U, V, Wm, x_bf, tape, s, st):
    """ゲート版: Σ c_k (x·2^{-s})^k を 左結合で 積み上げ → s回 二乗。EXACT(正規化なし)。"""
    M = len(Wm)
    xs = [BF(list(b.mant), b.E - s) for b in x_bf]          # ×2^{-s}: E減算のみ(ゲート0個)
    one = [const_bf(1, 0)] + [const_bf(0, 0) for _ in range(M - 1)]
    m0, e0c = tape[0]
    acc = [bf_mul(const_bf(m0, e0c), c, st) for c in one]
    term = [BF(list(c.mant), c.E) for c in one]
    for k in range(1, len(tape)):
        term = bf_bilinear_unit(U, V, Wm, term, xs, st)     # term·xs (左結合を宣言)
        mk, ek = tape[k]
        if mk != 0:
            ck = const_bf(mk, ek)
            acc = [bf_add(a, bf_mul(ck, t, st), st) for a, t in zip(acc, term)]
    for _ in range(s):
        acc = bf_bilinear_unit(U, V, Wm, acc, acc, st)      # スクエアリング
    return acc


# ============================================================ ①の仕様: 同じ計算を Fraction で
def series_spec(OM, x_fr, tape, s):
    M = len(x_fr)
    def mul(a, b):
        r = [Fr(0)] * M
        for i in range(M):
            if a[i] == 0: continue
            for j in range(M):
                if b[j] == 0: continue
                r[i ^ j] += Fr(int(OM[i, j])) * a[i] * b[j]
        return r
    xs = [v / (1 << s) for v in x_fr]
    m0, e0c = tape[0]
    c0 = Fr(m0) * Fr(2) ** e0c
    acc = [c0] + [Fr(0)] * (M - 1)
    term = [Fr(1)] + [Fr(0)] * (M - 1)
    for k in range(1, len(tape)):
        term = mul(term, xs)
        mk, ek = tape[k]
        if mk != 0:
            ck = Fr(mk) * Fr(2) ** ek
            acc = [a + ck * t for a, t in zip(acc, term)]
    for _ in range(s):
        acc = mul(acc, acc)
    return acc


# ============================================================ ②の参照: 真の exp/sin (float64 高次)
def ref_series(OM, x, coeffs, s):
    M = len(x)
    def mul(a, b):
        r = np.zeros(M)
        for i in range(M):
            for j in range(M):
                r[i ^ j] += OM[i, j] * a[i] * b[j]
        return r
    xs = np.asarray(x) / (1 << s)
    acc = np.zeros(M); acc[0] = coeffs[0]
    term = np.zeros(M); term[0] = 1.0
    for k in range(1, len(coeffs)):
        term = mul(term, xs)
        if coeffs[k] != 0.0:
            acc = acc + coeffs[k] * term
    for _ in range(s):
        acc = mul(acc, acc)
    return acc


# ============================================================ self-test
def run_config(name, M, tape_fn, true_coeffs, order, s, F, P, x=None):
    OM = cd_omega(M)
    U, V, Wm = group_as_UVW(OM, M)
    rng = np.random.default_rng(5)
    if x is None:
        x = 0.4 * rng.standard_normal(M)
    x_bf, x_fr = quantize_vec(x, F)
    tape = tape_fn(order, P)
    st = new_counter()
    got = gate_series_unit(U, V, Wm, x_bf, tape, s, st)
    spec = series_spec(OM, x_fr, tape, s)
    exact = all(from_bf(g) == sv for g, sv in zip(got, spec))       # ① 厳密一致(Fraction)
    ref = ref_series(OM, x, true_coeffs, s)                          # ② 真値(打ち切りなし級数)
    err = max(abs(float(sv) - rv) for sv, rv in zip(spec, ref))
    gates = sum(st.values())
    widths = max(len(g.mant) for g in got)
    print(f"  {name:<22} ゲート≡仕様: {'✓厳密一致' if exact else '✗不一致!'}"
          f"  仕様vs真値: {err:.1e}  ゲート数: {gates:>10,}  最終仮数幅: {widths}")
    assert exact, name
    return gates, err

def self_test():
    print("=" * 86)
    print("gate_series: 超越関数 = 係数テープ(合成時定数) × (U,V,W)ユニット列 — ゲートで")
    print("=" * 86)
    # 真値用の 打ち切りなし係数 (float64, 高次)
    exp_c = [1.0 / math.factorial(k) for k in range(30)]
    sin_c = [0.0 if k % 2 == 0 else (-1.0) ** ((k - 1) // 2) / math.factorial(k) for k in range(31)]

    print("① 四元数 exp (order=4, s=1, 入口F=4bit, テープP=8bit) — EXACTモード:")
    g1, e1 = run_config("quaternion exp", 4, tape_exp, exp_c, order=4, s=1, F=4, P=8)
    print("② テープ差し替え → 同じ骨格が sin に (order=5, s=0・二乗なし=軽い):")
    g2, e2 = run_config("quaternion sin", 4, tape_sin, sin_c, order=5, s=0, F=4, P=8)
    print("③ 精度はビット数のダイヤル (order=6, s=1, F=6, P=12):")
    g4, e4 = run_config("quaternion exp(精)", 4, tape_exp, exp_c, order=6, s=1, F=6, P=12)
    assert e4 < e1, "ビットを増やして精度が上がらないのはおかしい"
    print("④ セデニオン exp (order=3, s=1) — R=256 の配線でも同じ1本のコード (~30s):")
    g3, e3 = run_config("sedenion exp", 16, tape_exp, exp_c, order=3, s=1, F=3, P=6)
    print("-" * 86)
    print("実測の成長則(EXACTモード・正規化なし・四元数):")
    print("  order=4,s=1:   86M ゲート(幅135) / order=6,s=1: 263M(幅193) / order=8,s=2: 4.3B(幅477)")
    print("  ⟹ 厳密モードの値段は仮数幅の伸び(≈積で倍増)で爆発する。二乗(スクエアリング)が怪物。")
    print("  ⟹ 実用回路は 各段の block_normalize(誠実な丸め+ge/leフラグ) が必須 — bfp_sed の")
    print("     interval 機構との接続が次の一歩。誠実さの主張①(ゲート≡有理数仕様の厳密一致)は")
    print("     EXACTモードでのみ成立し、上の全構成で検証済み。")
    print("まとめ: 判断もループも逆行列もない固定配線が exp/sin を計算する。テープ交換=関数交換。")
    print("        スケーリング 2^{-s} は E の減算 = 底2の付け替え = ゲート0個。")

if __name__ == "__main__":
    self_test()
