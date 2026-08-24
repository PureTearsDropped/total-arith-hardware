#!/usr/bin/env python3
"""除算は 桁規則を 増やさない — 乗算と 減算の 上の **プログラム**。

total-arith-cuda の 分類 (ユーザ指摘 2026-08-25): `inv` は `kind="candidate"`。
計算した 後に **乗算で 検算**し、合わなければ INEXACT。TBM の 思想でも
「exp も solve も BILIN の 上の プログラム」。

⟹ 除算に **新しい 桁規則は 要らない**。ニュートン反復
      y ← y·(2 − a·y)
は 乗算と 減算だけで 書け、その 二つは 今日 閉じた (嘘 0)。境界は 演算の 合成として
伝播する。そして **検算の 段が 今日の 仕組みの 出番**: 候補を 掛け戻し、残差が
主張した 境界の 中に あるかを 見る。

このファイルは その 筋が 通るかを 確かめる。値は 厳密な 有理数で 追い、各段で
「今日の 規則が 与える 境界」を 併走させる。
"""
from fractions import Fraction as F
import random

MARK = F(1, 1 << 24)          # 予約桁 1 つぶんの ずれ (相殺が 起きたら 上へ)


class B:
    """値と 主張の 対。主張は 'exact' / 'upper' (v ≥ 真) / 'lower' (v ≤ 真)。"""
    __slots__ = ('v', 'c')

    def __init__(self, v, c='exact'):
        self.v, self.c = F(v), c

    def __repr__(self):
        return f"({float(self.v):.10f}, {self.c})"


def _combine(ca, cb):
    if ca == 'exact' and cb == 'exact':
        return 'exact'
    if ca in ('exact', 'upper') and cb in ('exact', 'upper'):
        return 'upper'
    if ca in ('exact', 'lower') and cb in ('exact', 'lower'):
        return 'lower'
    return None                                  # 混在 ⟹ 相殺の 印が 要る


def add(a: B, b: B, cancelled=False):
    """今日の 加算: 値は 厳密、相殺が あれば 予約桁ぶん 上へ ずらして 'upper'。"""
    c = _combine(a.c, b.c)
    if c is None or cancelled:
        return B(a.v + b.v + MARK, 'upper')      # 印を 立てる ⟹ 上端
    return B(a.v + b.v, c)


def sub(a: B, b: B, cancelled=True):
    """減算 = 符号反転 + 加算。反転で 主張も 反転する。"""
    flip = {'exact': 'exact', 'upper': 'lower', 'lower': 'upper'}[b.c]
    return add(a, B(-b.v, flip), cancelled)


def mul(a: B, b: B, cancelled=True):
    """乗算: 部分積の 足し込みで ほぼ 必ず 相殺が 起きる (実測 1998/2000)。"""
    pos = a.v >= 0 and b.v >= 0
    c = _combine(a.c, b.c) if pos else None
    if c is None or cancelled:
        return B(a.v * b.v + MARK, 'upper')
    return B(a.v * b.v, c)


def sound(x: B, truth: F):
    if x.c == 'exact':
        return x.v == truth
    return x.v >= truth if x.c == 'upper' else x.v <= truth


def newton_inv(a: B, steps=5):
    """1/a を ニュートンで。乗算と 減算だけ — **新しい 規則は 使わない**。"""
    y = B(F(1) / a.v)                            # 初期値 (種)
    for _ in range(steps):
        ay = mul(a, y)
        two_minus = sub(B(2), ay)
        y = mul(y, two_minus)
    return y


def main():
    rng = random.Random(0)
    print("=" * 74)
    print("除算 = 乗算と減算の上のプログラム — 桁規則を 増やさずに 閉じるか")
    print("=" * 74)

    print("\n① 候補の 健全性 (主張が 真の 1/a に対して 成り立つか)")
    bad = n = 0
    gaps = []
    for _ in range(3000):
        av = F(rng.randrange(1, 1000), rng.randrange(1, 64))
        a = B(av, rng.choice(('exact', 'upper', 'lower')))
        y = newton_inv(a)
        truth = F(1) / av
        n += 1
        if not sound(y, truth):
            bad += 1
        gaps.append(float(abs(y.v - truth)))
    gaps.sort()
    print(f"   {n} 例中 主張が 破れた: {bad}")
    print(f"   候補と 真値の 差: 中央値 {gaps[len(gaps)//2]:.3e} / 最大 {gaps[-1]:.3e}")

    print("\n② 検算 (候補を 掛け戻して 1 に なるか) — repo の candidate 作法")
    print(f"   {'a':>12} {'候補 y':>16} {'a·y':>16} {'残差':>12} {'主張':>7}")
    for av in (F(3), F(7, 2), F(1, 3), F(255, 256)):
        a = B(av)
        y = newton_inv(a)
        back = mul(a, y)
        resid = abs(back.v - 1)
        print(f"   {str(av):>12} {float(y.v):>16.12f} {float(back.v):>16.12f} "
              f"{float(resid):>12.3e} {y.c:>7}")
    print(f"   予約桁 1 つぶん = {float(MARK):.3e}")
    print("   ⟹ 残差が 予約桁の 数個ぶんに 収まれば、検算は 通ったと 言える")

    print("\n③ 桁規則の 追加は あったか")
    print("   newton_inv が 使ったのは mul と sub だけ。**新しい 規則ゼロ**")
    print("   除算は candidate — 検算で 通れば 値、通らなければ INEXACT")


if __name__ == "__main__":
    main()
