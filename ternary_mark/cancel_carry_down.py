#!/usr/bin/env python3
"""相殺したら 下の桁に 印を 押し出す — 壁を 越えられるか。

ユーザの 案 (2026-08-25)。これまでの 四規約は 全部「相殺したら 情報が 消える」で
詰まった。この案は **消さずに 下へ 逃がす**:

    位置 k で −1 と +1 が 出会って 0 に なったら、0 を 書いた 上で
    **一つ下の 位置 (k+1) に +1 を 足す**

値は 2^-(k+1) だけ ずれる。しかし **ずれの 向きが 常に 上** なら、結果は
「真の 和 以上」= **LE の 片側境界**に なる。向き不明を、緩いが 確実な 境界に 変換する。

**測る**: ① ずれの 向きは 常に 一定か ② 健全か (真の和 ≤ 出力) ③ どれだけ 緩いか
"""
from fractions import Fraction as F
import random

D = 20                        # 小数部の 桁数


def val(frac):
    return sum(F(d, 1 << (i + 1)) for i, d in enumerate(frac))


def plain_add(x, y):
    """普通の 符号つき桁 加算 (低位→高位に 桁上げ)。厳密。"""
    n = len(x)
    out = [0] * n
    c = 0
    for i in range(n - 1, -1, -1):
        s = x[i] + y[i] + c
        c = 0
        while s > 1:
            s -= 2; c += 1
        while s < -1:
            s += 2; c -= 1
        out[i] = s
    return out, c


def cancel_add(x, y):
    """相殺したら 下へ 印を 押す 加算。

    位置 i で 入力が 打ち消し合って 0 に なったら、その 一つ下に +1 を 足す。
    押し出した 分は 値を 上へ ずらす ⟹ 出力 ≥ 真の和 が 期待される。
    """
    n = len(x)
    xs, ys = list(x), list(y)
    push = [0] * (n + 1)                       # 下へ 押し出す 印
    for i in range(n):
        if xs[i] + ys[i] == 0 and (xs[i] != 0 or ys[i] != 0):
            push[i + 1] += 1                   # **一つ下に +1**
    # 押し出しを 含めて 足す
    n2 = n
    out = [0] * n2
    c = 0
    for i in range(n2 - 1, -1, -1):
        s = xs[i] + ys[i] + push[i] + c
        c = 0
        while s > 1:
            s -= 2; c += 1
        while s < -1:
            s += 2; c -= 1
        out[i] = s
    return out, c, sum(push)


def main():
    rng = random.Random(0)
    print("=" * 74)
    print("相殺したら 下の桁に +1 を 押す — 健全性と 緩さ")
    print("=" * 74)

    print("\n① ずれの 向きは 常に 一定か (出力 − 真の和)")
    signs = {'+': 0, '0': 0, '-': 0}
    diffs = []
    for _ in range(4000):
        x = [rng.choice((-1, 0, 1)) for _ in range(D)]
        y = [rng.choice((-1, 0, 1)) for _ in range(D)]
        ref, c0 = plain_add(x, y)
        got, c1, npush = cancel_add(x, y)
        d = (val(got) + c1) - (val(ref) + c0)
        diffs.append(d)
        signs['+' if d > 0 else ('-' if d < 0 else '0')] += 1
    print(f"   出力 > 真の和: {signs['+']}   等しい: {signs['0']}   出力 < 真の和: {signs['-']}")

    print("\n② 健全性 (常に 出力 ≥ 真の和 なら LE の 境界として 使える)")
    bad = sum(1 for d in diffs if d < 0)
    print(f"   4000 例中 出力が 真の和を 下回った: {bad}")

    print("\n③ 緩さ (ずれの 大きさ)")
    ds = sorted(float(d) for d in diffs)
    print(f"   中央値 {ds[len(ds)//2]:.6f} / 最大 {ds[-1]:.6f} / 最小 {ds[0]:.6f}")
    print(f"   参考: 最下位の 重み 2^-{D} = {2.0**-D:.3e}")

    print("\n④ 相殺が 起きなかった 場合は 厳密の まま か")
    n = ex = 0
    for _ in range(4000):
        x = [rng.choice((-1, 0, 1)) for _ in range(D)]
        y = [0 if xi != 0 else rng.choice((-1, 0, 1)) for xi in x]   # 相殺が 起きない 組
        ref, c0 = plain_add(x, y)
        got, c1, npush = cancel_add(x, y)
        n += 1
        if npush == 0 and (val(got) + c1) == (val(ref) + c0):
            ex += 1
    print(f"   相殺なしの 組 {n} 例中 厳密の まま: {ex}")

    print("\n⑤ 相殺の 有無が 出力から 読めるか")
    print("   押し出しが あれば 出力は 真の和より 上 ⟹ 「LE」を 主張できる")
    print("   押し出しが 無ければ 厳密。**出力の 桁だけで 区別できるかは 別問題**")


if __name__ == "__main__":
    main()
