#!/usr/bin/env python3
# ⚠️ 生成AI使用・要検証
"""乗算器 × 複数値同時加算 — 積を 解決しない 融合（fused）乗算器。

  観察: 乗算器の 中身は「部分積の 多入力加算」そのもの（multiply_fast が 既に 木）。
  さらに 進める:
   ① multiply_red : 積を **冗長 2 行のまま** 返す（最後の sd_add2 すら しない）
   ② mac_dot      : 内積 Σ aᵢ·bᵢ — **全乗算の 全部分積を 1 本の 木**に 流し込む
                    （積ごとの 解決なし・解決は 全体で 1 回 = fused dot product）
   ③ 群積の 成分   : (a·b)_k = Σᵢ ±aᵢ·b_{i⊕k} — 符号は neg（配線）・16 積を 1 本の 木
                    ＝ セデニオン/四元数 乗算器の 各成分が 融合 MAC そのもの
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from gate_bilinear import new_counter, ZERO, to_sd, from_sd, neg, gate9
from gate_fast import (wrap, depth_of, sd_sum_fast, sd_add2, multiply_fast, or_tree)


def _pp_rows(X, Y, st):
    """部分積の 行リスト（gate9・位置 i+j）。乗算器の 素材。"""
    rows = []
    for i, xd in enumerate(X):
        for j, yd in enumerate(Y):
            rows.append([ZERO] * (i + j) + [gate9(xd, yd, st)])
    return rows

def _compress_to2(rows, st):
    """層別 3:2 圧縮で 2 行まで（解決しない）。sd_sum_fast の 前半だけ。"""
    if not rows: return [ZERO], [ZERO]
    width = max(len(x) for x in rows)
    cols = [[] for _ in range(width + 2)]
    for x in rows:
        for i, dg in enumerate(x):
            cols[i].append(dg)
    from gate_bilinear import compress3
    while max(len(c) for c in cols) > 2:
        nxt = [[] for _ in range(len(cols) + 1)]
        for k, c in enumerate(cols):
            i = 0
            while len(c) - i >= 3:
                low, high = compress3(c[i], c[i + 1], c[i + 2], st)
                nxt[k].append(low); nxt[k + 1].append(high)
                i += 3
            nxt[k].extend(c[i:])
        cols = nxt
    X = [c[0] if len(c) > 0 else ZERO for c in cols]
    Y = [c[1] if len(c) > 1 else ZERO for c in cols]
    return X, Y

def multiply_red(X, Y, st):
    """① 融合用 乗算器: 積を **冗長 2 行**で 返す（解決なし）。"""
    return _compress_to2(_pp_rows(X, Y, st), st)

def mac_dot(pairs, st):
    """② 内積 Σ aᵢ·bᵢ: 全乗算の 全部分積を 1 本の 木へ → 解決 1 回（sd_add2）。"""
    rows = []
    for X, Y in pairs:
        rows.extend(_pp_rows(X, Y, st))
    A, Bo = _compress_to2(rows, st)
    return sd_add2(A, Bo, st)

def group_component(a_digits, b_digits, OM, M, k, st):
    """③ 群積の 成分 k = Σᵢ σ(i,i⊕k)·aᵢ·b_{i⊕k}。符号は neg（配線）・M 積を 1 本の 木。"""
    rows = []
    for i in range(M):
        j = i ^ k
        pp = _pp_rows(a_digits[i], b_digits[j], st)
        if OM[i][j] < 0:
            pp = [neg(r) for r in pp]                      # 符号 = 配線（ゲート0）
        rows.extend(pp)
    A, Bo = _compress_to2(rows, st)
    return sd_add2(A, Bo, st)


def self_test():
    import numpy as np
    rng = np.random.default_rng(20260807)

    print("=" * 78)
    print("① multiply_red — 積を 冗長 2 行のまま（解決なし・次段へ 直行）")
    print("=" * 78)
    bad = 0
    for _ in range(1500):
        a = int(rng.integers(-300, 300)); b = int(rng.integers(-300, 300))
        A, Bo = multiply_red(to_sd(a, 10), to_sd(b, 10), new_counter())
        if from_sd(A) + from_sd(Bo) != a * b: bad += 1
    A, Bo = multiply_red(wrap(to_sd(217, 10)), wrap(to_sd(-178, 10)), new_counter())
    d_red = depth_of([A, Bo])
    d_full = depth_of(multiply_fast(wrap(to_sd(217, 10)), wrap(to_sd(-178, 10)), new_counter()))
    print(f"  2 行の 和 == 積: 違反 {bad}/1500 ✓   深さ: 解決つき {d_full} → 冗長のまま **{d_red}**")

    print()
    print("=" * 78)
    print("② mac_dot — 内積の 全部分積を 1 本の 木・解決も 丸めも 全体で 1 回")
    print("=" * 78)
    bad = 0
    for _ in range(600):
        N = int(rng.integers(2, 9))
        av = [int(v) for v in rng.integers(-200, 200, N)]
        bv = [int(v) for v in rng.integers(-200, 200, N)]
        got = from_sd(mac_dot([(to_sd(x, 10), to_sd(y, 10)) for x, y in zip(av, bv)],
                              new_counter()))
        if got != sum(x * y for x, y in zip(av, bv)): bad += 1
    N = 8
    av = [int(v) for v in rng.integers(-200, 200, N)]
    bv = [int(v) for v in rng.integers(-200, 200, N)]
    d_fused = depth_of(mac_dot([(wrap(to_sd(x, 10)), wrap(to_sd(y, 10)))
                                for x, y in zip(av, bv)], new_counter()))
    # 対照: 積を 個別に 解決 → 結果を 木で 加算
    st = new_counter()
    prods = [multiply_fast(wrap(to_sd(x, 10)), wrap(to_sd(y, 10)), st) for x, y in zip(av, bv)]
    d_sep = depth_of(sd_sum_fast(prods, st))
    print(f"  内積 == Σaᵢbᵢ: 違反 {bad}/600 ✓")
    print(f"  深さ（N=8）: 個別解決→加算 {d_sep} → **融合 {d_fused}**（途中解決 {N} 回が 消える）")

    print()
    print("=" * 78)
    print("③ 群積の 成分 — セデニオン/四元数の 各成分 = 融合 MAC（符号は 配線）")
    print("=" * 78)
    from nd_algebra import cd_omega, ref_mult_M
    for M, name in [(4, "四元数"), (16, "セデニオン")]:
        OM = cd_omega(M)
        OMl = [[int(OM[i, j]) for j in range(M)] for i in range(M)]
        bad = 0
        for _ in range(120 if M == 4 else 30):
            a = [int(v) for v in rng.integers(-9, 10, M)]
            b = [int(v) for v in rng.integers(-9, 10, M)]
            ref = ref_mult_M(a, b, OM, M)
            st = new_counter()
            for k in range(M):
                got = from_sd(group_component([to_sd(v, 6) for v in a],
                                              [to_sd(v, 6) for v in b], OMl, M, k, st))
                if got != ref[k]: bad += 1
        ad = [wrap(to_sd(v, 6)) for v in a]; bd = [wrap(to_sd(v, 6)) for v in b]
        d_comp = depth_of(group_component(ad, bd, OMl, M, 1, new_counter()))
        print(f"  {name:<6} 全成分 == 参照: 違反 {bad} ✓   成分 1 個の 融合深さ **{d_comp}**"
              f"（{M} 積 × 36 部分積を 1 本の 木）")

    print()
    print("乗算器 = 部分積の 多入力加算器。融合すれば 積の 解決は 消え、")
    print("内積・群積成分は「全部分積 → 1 本の 木 → 解決 1 回」= MAC が 原子演算になる。")


if __name__ == "__main__":
    self_test()
