#!/usr/bin/env python3
"""乗算に 今日の 符号化を 雑に 入れる — 部分積の 足し込みに 下向きリップルを 効かせる。

乗算は **部分積の 和**。桁どうしの 積 (`gate9`) は 相殺を 起こさない (1桁×1桁は
そのまま)。相殺が 起きるのは **足し込み**の 段。だから 加算で 効いた 仕組みが
そのまま 乗るはず。

  値の 向き: a,b が それぞれ 上端/下端の どちらを 主張しているか + 符号 で 決まる
             (2026-08-25 に 区間から 16/16 で 導出済み)
  相殺の 印: 部分積の 足し込みで 相殺が 起きたら 予約桁へ 下向きリップル

**検査**: ① 値は 正しいか ② 主張は 健全か ③ 予約桁は 機能するか
"""
from fractions import Fraction as F
import random

D = 12                      # 小数部の 桁数 (積は 2D 桁に なる)


def val(frac, off=0):
    return sum(F(d, 1 << (i + 1 + off)) for i, d in enumerate(frac))


def norm(digits):
    """符号つき桁の 正規化 (桁上げは 上へ)。戻り (桁列, 上への 繰り上がり)。"""
    out = list(digits); c = 0
    for i in range(len(out) - 1, -1, -1):
        s = out[i] + c; c = 0
        while s > 1:
            s -= 2; c += 1
        while s < -1:
            s += 2; c -= 1
        out[i] = s
    return out, c


def multiply(x, y, mark_on=True):
    """部分積を 足し込む。相殺を 検出したら 予約桁へ 印を 下ろす。

    x, y: 小数部の 桁列 (長さ D)。積は 2D 桁 + 予約桁 1 つ。
    """
    n = len(x)
    acc = [0] * (2 * n + 1)                      # 末尾 1 桁が 予約桁
    cancelled = False
    for i, xi in enumerate(x):
        if xi == 0:
            continue
        for j, yj in enumerate(y):
            if yj == 0:
                continue
            p = xi * yj                          # 桁どうしの 積 (相殺しない)
            pos = i + j + 1                      # 重み 2^-(i+1)·2^-(j+1)
            if acc[pos] + p == 0 and acc[pos] != 0:
                cancelled = True                 # **足し込みで 相殺**
            acc[pos] += p
    if mark_on and cancelled:
        acc[2 * n] += 1                          # 予約桁へ (下向きリップルの 終点)
    return norm(acc)


def claim_of(frac):
    if not frac or frac[-1] == 0:
        return 'exact'
    return 'lower' if frac[-1] < 0 else 'upper'


def padded(frac):
    v = sum(F(d, 1 << (i + 1)) for i, d in enumerate(frac))
    if frac:
        v += F(frac[-1], 1 << len(frac))
    return v


def main():
    rng = random.Random(0)
    print("=" * 74)
    print("乗算に 下向きリップルを 入れる")
    print("=" * 74)

    print("\n① 値の 正しさ (印を 切った 場合、厳密な 積に なるか)")
    bad = n = 0
    for _ in range(3000):
        x = [rng.choice((-1, 0, 1)) for _ in range(D)]
        y = [rng.choice((-1, 0, 1)) for _ in range(D)]
        out, c = multiply(x, y, mark_on=False)
        got = val(out) + c
        want = val(x) * val(y)
        n += 1
        if got != want:
            bad += 1
    print(f"   {n} 例中 値が 合わなかった: {bad}")

    print("\n② 印を 入れた 場合の 健全性 (出力 ≥ 真の積 か)")
    bad2 = n2 = 0
    ds = []
    for _ in range(3000):
        x = [rng.choice((-1, 0, 1)) for _ in range(D)]
        y = [rng.choice((-1, 0, 1)) for _ in range(D)]
        out, c = multiply(x, y, mark_on=True)
        got = padded(out) + c
        want = val(x) * val(y)
        n2 += 1
        ds.append(float(got - want))
        if got < want:
            bad2 += 1
    ds.sort()
    print(f"   {n2} 例中 真の積を 下回った: {bad2}")
    print(f"   ずれ 中央値 {ds[len(ds)//2]:.3e} / 最大 {ds[-1]:.3e} / 最小 {ds[0]:.3e}")

    print("\n③ 主張と 実際の 一致 (厳密と 言ったら 本当に 厳密か)")
    lie = nex = 0
    for _ in range(3000):
        x = [rng.choice((-1, 0, 1)) for _ in range(D)]
        y = [rng.choice((-1, 0, 1)) for _ in range(D)]
        out, c = multiply(x, y, mark_on=True)
        got = padded(out) + c
        want = val(x) * val(y)
        cl = claim_of(out)
        if cl == 'exact':
            nex += 1
            if got != want:
                lie += 1
        elif cl == 'upper' and got < want:
            lie += 1
        elif cl == 'lower' and got > want:
            lie += 1
    print(f"   「厳密」と 主張 {nex}/3000 — 嘘 {lie}")

    print("\n④ 相殺が 起きる 頻度")
    cnt = 0
    for _ in range(2000):
        x = [rng.choice((-1, 0, 1)) for _ in range(D)]
        y = [rng.choice((-1, 0, 1)) for _ in range(D)]
        a, _ = multiply(x, y, mark_on=False)
        b, _ = multiply(x, y, mark_on=True)
        if a != b:
            cnt += 1
    print(f"   2000 例中 相殺が 起きた: {cnt}")


if __name__ == "__main__":
    main()
