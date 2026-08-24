#!/usr/bin/env python3
"""印は 下へ リップルする — 小数点とは 無関係の 一様規則。

ユーザの 気づき (2026-08-25): `−1 + 1` は 位置に よらず 下へ キャリーする。
小数部だけの 話ではない。⟹ **桁上げは 上へ、相殺の印は 下へ**。逆向きの 二本の リップル。

前の 実装は 「全位置の 相殺を OR で 集めて 予約桁に 置く」だった。それは この
下向きリップルの **結果を 先回りして 書いたもの**。局所規則として 書き直すと
特別扱いが ゼロに なる。

**検査**: 局所リップル版が OR 版と 一致するか / 値は 正しいままか。
"""
from fractions import Fraction as F
import random

D = 20            # 小数部の 桁数 (予約桁は その下に 1 つ)


def val(frac):
    return sum(F(d, 1 << (i + 1)) for i, d in enumerate(frac))


def plain_add(x, y):
    n = len(x); out = [0] * n; c = 0
    for i in range(n - 1, -1, -1):
        s = x[i] + y[i] + c; c = 0
        while s > 1:
            s -= 2; c += 1
        while s < -1:
            s += 2; c -= 1
        out[i] = s
    return out, c


def add_or(x, y):
    """OR 版: 全位置の 相殺を 集めて 予約桁に 置く。"""
    n = len(x)
    mark = 1 if any(x[i] + y[i] == 0 and (x[i] != 0 or y[i] != 0) for i in range(n)) else 0
    xs = list(x) + [0]; ys = list(y) + [0]
    out = [0] * (n + 1); c = 0
    for i in range(n, -1, -1):
        s = xs[i] + ys[i] + (mark if i == n else 0) + c; c = 0
        while s > 1:
            s -= 2; c += 1
        while s < -1:
            s += 2; c -= 1
        out[i] = s
    return out, c


def add_ripple(x, y):
    """局所リップル版: 各位置で 相殺を 検出したら 印を **一つ下へ 渡す**。

    印は 値を 持たず、位置を 下りるだけ。最下段 (予約桁) で 初めて 値に なる。
    ⟹ どの 位置でも 同じ 規則。小数点も 整数部も 区別しない。
    """
    n = len(x)
    xs = list(x) + [0]; ys = list(y) + [0]
    # 上から 下へ: 印を 受け取り、自分の 相殺と 合わせて 下へ 渡す
    passing = 0
    landed = 0
    for i in range(n + 1):
        cancel = 1 if (xs[i] + ys[i] == 0 and (xs[i] != 0 or ys[i] != 0)) else 0
        passing = 1 if (passing or cancel) else 0
        if i == n:                       # 予約桁に 到達 ⟹ ここで 値に なる
            landed = passing
    out = [0] * (n + 1); c = 0
    for i in range(n, -1, -1):
        s = xs[i] + ys[i] + (landed if i == n else 0) + c; c = 0
        while s > 1:
            s -= 2; c += 1
        while s < -1:
            s += 2; c -= 1
        out[i] = s
    return out, c


def main():
    rng = random.Random(0)
    print("=" * 74)
    print("下向きリップル版 vs OR 版 — 一致するか")
    print("=" * 74)
    diff = 0; n = 0; bad = 0
    for _ in range(6000):
        x = [rng.choice((-1, 0, 1)) for _ in range(D)]
        y = [rng.choice((-1, 0, 1)) for _ in range(D)]
        a, ca = add_or(x, y)
        b, cb = add_ripple(x, y)
        n += 1
        if a != b or ca != cb:
            diff += 1
        ref, c0 = plain_add(x, y)
        if (val(b) + cb) < (val(ref) + c0):
            bad += 1
    print(f"\n  {n} 例中 二つの版が 食い違った: {diff}")
    print(f"  {n} 例中 真の和を 下回った (健全性): {bad}")
    print()
    print("  ⟹ 一致するなら、OR は 下向きリップルの 別名にすぎない。")
    print("     回路は 「桁上げは 上へ、印は 下へ」の 二本の リップルだけ。")
    print("     小数点の 位置を 知る 必要が ない = **完全に 一様**")


if __name__ == "__main__":
    main()
