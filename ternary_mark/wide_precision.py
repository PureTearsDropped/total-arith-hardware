#!/usr/bin/env python3
"""ternary_mark の 符号化を float128 / float256 の 幅で 回す。

**主張**: 今日の 規則は 桁数に 依存しない。
  ・端のトリットを 無限に 繰り返す (padding)
  ・相殺は 下向きリップルで 予約桁へ
  ・桁上げは 上へ
どれも 隣の桁としか 話さないので、W を 変えても 規則は 同じ。
⟹ **W = 113 で float128 精度、W = 237 で float256 精度** に なるだけ のはず。

「はず」で 終わらせずに、その 幅で 実際に 回して 健全性と 費用を 測る。

註: IEEE binary128 は 仮数 113bit・指数 15bit、binary256 は 仮数 237bit・指数 19bit。
ここでは **仮数の 桁数**を それに 合わせる (指数は 共有なので 別勘定・広げるのは 安い)。
"""
from fractions import Fraction as F
import random
import time


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


def mark_add(x, y):
    """相殺の 印を 予約桁へ。規則は 幅に 依存しない。"""
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


def val(frac):
    return sum(F(d, 1 << (i + 1)) for i, d in enumerate(frac))


def padded(frac):
    v = val(frac)
    if frac:
        v += F(frac[-1], 1 << len(frac))
    return v


def claim(frac):
    if not frac or frac[-1] == 0:
        return 'exact'
    return 'lower' if frac[-1] < 0 else 'upper'


def sound(got, truth, cl):
    if cl == 'exact':
        return got == truth
    return got >= truth if cl == 'upper' else got <= truth


WIDTHS = [(24, 'float32 相当'), (53, 'float64 相当'),
          (113, 'float128 相当'), (237, 'float256 相当'), (512, '512bit')]


def main():
    rng = random.Random(0)
    print("=" * 78)
    print("ternary_mark を float128 / float256 の 幅で 回す")
    print("=" * 78)

    print("\n① 健全性 — 幅を 変えても 嘘は 出ないか (各 800 例)")
    print(f"   {'仮数桁 W':>10} {'相当':>16} {'厳密と主張':>11} {'嘘':>5} {'ずれ':>12}")
    for W, lab in WIDTHS:
        lie = nex = 0
        dmax = F(0)
        for _ in range(800):
            x = [rng.choice((-1, 0, 1)) for _ in range(W)]
            y = [rng.choice((-1, 0, 1)) for _ in range(W)]
            ref, c0 = plain_add(x, y); truth = val(ref) + c0
            out, c1 = mark_add(x, y); got = padded(out) + c1
            cl = claim(out)
            if cl == 'exact':
                nex += 1
            if not sound(got, truth, cl):
                lie += 1
            dmax = max(dmax, abs(got - truth))
        print(f"   {W:>10} {lab:>16} {nex:>11} {lie:>5} {float(dmax):>12.3e}")

    print("\n② 減算 (符号反転 + 加算) — 幅を 変えても 通るか")
    print(f"   {'仮数桁 W':>10} {'嘘':>5} {'ずれ':>12}")
    for W, lab in WIDTHS[:4]:
        lie = 0; dmax = F(0)
        for _ in range(500):
            x = [rng.choice((-1, 0, 1)) for _ in range(W)]
            y = [rng.choice((-1, 0, 1)) for _ in range(W)]
            ny = [-d for d in y]
            ref, c0 = plain_add(x, ny); truth = val(ref) + c0
            out, c1 = mark_add(x, ny); got = padded(out) + c1
            if not sound(got, truth, claim(out)):
                lie += 1
            dmax = max(dmax, abs(got - truth))
        print(f"   {W:>10} {lie:>5} {float(dmax):>12.3e}")

    print("\n③ ずれは 予約桁 1 つぶんで 一定か (幅に 依存しないか)")
    for W, lab in WIDTHS:
        print(f"   W={W:>4}: 予約桁の padding 込みの 重み = 2·2^-{W + 1} = {2 * 2.0 ** -(W + 1):.3e}")

    print("\n④ 費用 — 加算は 桁数に 線形か")
    print(f"   {'仮数桁 W':>10} {'1 回の 時間':>14} {'W=24 比':>10} {'桁あたり':>12}")
    base = None
    for W, lab in WIDTHS:
        x = [rng.choice((-1, 0, 1)) for _ in range(W)]
        y = [rng.choice((-1, 0, 1)) for _ in range(W)]
        t0 = time.perf_counter()
        for _ in range(2000):
            mark_add(x, y)
        dt = (time.perf_counter() - t0) / 2000
        if base is None:
            base = dt
        print(f"   {W:>10} {dt * 1e6:>12.2f}µs {dt / base:>9.1f}倍 {dt / W * 1e9:>10.1f}ns")

    print("\n   ⟹ 桁あたりが ほぼ 一定なら 線形。規則が 局所である ことの 帰結")


if __name__ == "__main__":
    main()
