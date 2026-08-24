#!/usr/bin/env python3
"""4 値の 乗算 v2 — 導出した 3:2 圧縮器で 差し替える。

v1 は 桁上げの 分解を **全部 union** していたので、`?` の 無い 入力でも 63 倍 緩かった
（5×6=30 に 対し 63 通りを 主張）。ここでは 今日 導出・全列挙検証した 4 値 圧縮器
（64 通りで 健全かつ 最小）を 使う。

  nq ≥ 2 → (low, high) = (?, ?)
  nq = 1 → low = ?、他 2 つの 和 t2 が 偶数なら high = t2/2、奇数なら high = ?
  nq = 0 → 三値の carry-free 分解（表引き）

オラクルは v1 と 同じ **機構から独立**なもの: 入力を 全具体化 → 真の積の 集合 →
出力の 語が 表す 集合に 含まれるか。
"""
import itertools

QS = {-1: (-1,), 0: (0,), 1: (1,), '?': (-1, 0, 1)}
NAMES = (-1, 0, 1, '?')

# 三値 carry-free 分解: t = low + 2·high（low,high ∈ {−1,0,1}）
TERN = {-3: (-1, -1), -2: (0, -1), -1: (-1, 0), 0: (0, 0),
        1: (1, 0), 2: (0, 1), 3: (1, 1)}


def compress3_4(a, b, c):
    """4 値 3:2 圧縮器（今日 64 通りで 導出・検証した 規則）。"""
    unk = [d == '?' for d in (a, b, c)]
    nq = sum(unk)
    if nq >= 2:
        return '?', '?'
    if nq == 1:
        t2 = sum(d for d in (a, b, c) if d != '?')
        return ('?', t2 // 2) if t2 % 2 == 0 else ('?', '?')
    return TERN[a + b + c]


def gate9_4(a, b):
    s = frozenset(x * y for x in QS[a] for y in QS[b])
    for n in NAMES:
        if frozenset(QS[n]) == s:
            return n
    raise AssertionError


def reduce_columns(cols):
    """Wallace 木: 各列を 3:2 圧縮器で 減らす。high は 一つ上の 列へ。

    cols[k] = 位置 k（MSB 先頭）の 桁のリスト。
    """
    n = len(cols)
    changed = True
    while changed:
        changed = False
        for k in range(n - 1, -1, -1):          # 低位から
            while len(cols[k]) >= 3:
                a = cols[k].pop(); b = cols[k].pop(); c = cols[k].pop()
                lo, hi = compress3_4(a, b, c)
                cols[k].append(lo)
                if hi != 0:
                    if k - 1 >= 0:
                        cols[k - 1].append(hi)
                    # k-1 < 0 なら 桁溢れ（この試験では 起きない 幅を 取る）
                changed = True
        # 2 個 残った 列は 0 を 足して 3 個に して 潰す
        for k in range(n - 1, -1, -1):
            if len(cols[k]) == 2:
                a = cols[k].pop(); b = cols[k].pop()
                lo, hi = compress3_4(a, b, 0)
                cols[k].append(lo)
                if hi != 0 and k - 1 >= 0:
                    cols[k - 1].append(hi)
                changed = True
    return [c[0] if c else 0 for c in cols]


def mul4(a, b):
    n = len(a)
    width = 2 * n
    cols = [[] for _ in range(width)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            p = gate9_4(x, y)
            if p != 0:
                cols[i + j + 1].append(p)
    return reduce_columns(cols)


def val(digits):
    n = len(digits)
    return sum(d * (1 << (n - 1 - i)) for i, d in enumerate(digits))


def value_set(word):
    return frozenset(val(c) for c in itertools.product(*[QS[d] for d in word]))


def main():
    print("=" * 72)
    print("4 値 乗算 v2 — 導出した 圧縮器で 差し替え（オラクルは v1 と 同じ）")
    print("=" * 72)

    print("\n① まず 具体的な 入力（? なし）で 厳密に なるか — v1 は ここで 63 倍だった")
    for a, b in (([1, 0, 1], [1, 1, 0]), ([1, 1], [1, 1]), ([-1, 1, -1], [1, 0, 1])):
        out = mul4(list(a), list(b))
        got = value_set(out)
        truth = val(a) * val(b)
        mark = '✓ 厳密' if got == {truth} else f'✗ {len(got)} 通り'
        print(f"   {a} × {b} = {truth:>5}  出力={out}  {mark}")

    print("\n② 全列挙で 健全性と 締まり")
    for n in (2, 3):
        bad = tot = 0
        ratios = []
        from collections import defaultdict
        byq = defaultdict(list)
        for a in itertools.product(NAMES, repeat=n):
            for b in itertools.product(NAMES, repeat=n):
                truth = frozenset(x * y for x in value_set(a) for y in value_set(b))
                got = value_set(mul4(list(a), list(b)))
                tot += 1
                if not truth <= got:
                    bad += 1
                if truth:
                    r = len(got) / len(truth)
                    ratios.append(r)
                    byq[sum(1 for d in a + b if d == '?')].append(r)
        ratios.sort()
        print(f"\n   【桁数 {n}】{tot} 組")
        print(f"   健全性: 取りこぼし {bad}")
        print(f"   締まり: 中央値 {ratios[len(ratios)//2]:.2f} / 最大 {ratios[-1]:.2f}")
        print(f"   {'?の総数':>8} {'組数':>7} {'締まり中央値':>13} {'最大':>8}")
        for k in sorted(byq):
            v = sorted(byq[k])
            print(f"   {k:>8} {len(v):>7} {v[len(v)//2]:>13.2f} {v[-1]:>8.2f}")


if __name__ == "__main__":
    main()
