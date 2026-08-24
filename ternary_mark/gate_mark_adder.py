#!/usr/bin/env python3
"""ゲート段の 「印つき加算器」 — 今日の 符号化を total-arith-hardware の 部品で 組む。

research ブランチの 既存の 問いへの 答え:
  ・`?` を 4 値目に する 計画 → **要らない**。不確かさは padding と 予約桁に 住み、
    桁は 三値の まま (配線正規形を 変えない)
  ・`add_price.py` の 「加算は 符号か 締まりの どちらかを 払う」 → **予約桁 1 つを
    払えば 両方 残る**。しかも 払う量は 固定で 相殺が 何回 起きても 累積しない

**回路**:
  ・値の 経路: 既存の 符号つき桁 加算 (`sd_sum`) そのまま
  ・相殺の 検出: 位置 i で `x = +1 かつ y = −1` または `x = −1 かつ y = +1`
      (p,n) レールなら  `(xp & yn) | (xn & yp)`  = AND 2 + OR 1 = **3 ゲート/桁**
  ・印の 集約: OR 木 ⟹ 予約桁へ (下向きリップルの 終点)

**検査**: ① 値が 既存加算と 一致 ② 印が Python 模型と 一致 ③ ゲート数
"""
import sys, random
sys.path.insert(0, __import__('os').path.join(
    __import__('os').path.dirname(__file__), '..'))

import gate_bilinear as g


def cancel_bit(xd, yd, st):
    """位置 1 桁ぶんの 相殺検出。(p,n) レールで 3 ゲート。"""
    xp, xn = xd
    yp, yn = yd
    a = g.AND(xp, yn, st)          # x=+1 かつ y=−1
    b = g.AND(xn, yp, st)          # x=−1 かつ y=+1
    return g.OR(a, b, st)


def or_tree(bits, st):
    """印の 集約。段数 log。"""
    cur = list(bits)
    while len(cur) > 1:
        nxt = []
        for i in range(0, len(cur) - 1, 2):
            nxt.append(g.OR(cur[i], cur[i + 1], st))
        if len(cur) % 2:
            nxt.append(cur[-1])
        cur = nxt
    return cur[0] if cur else 0


def mark_add(X, Y, st):
    """印つき加算。X, Y は 符号つき桁の 列 (低位から)。

    戻り (和の 桁列, 印)。印は 予約桁に 立てる ビット。
    """
    marks = [cancel_bit(X[i], Y[i], st) for i in range(min(len(X), len(Y)))]
    mark = or_tree(marks, st)
    s = g.sd_sum([list(X), list(Y)], st)
    return s, mark


def py_model(xv, yv):
    """今日の Python 模型 (down_ripple と 同じ 規則)。"""
    mark = 1 if any(xv[i] + yv[i] == 0 and (xv[i] != 0 or yv[i] != 0)
                    for i in range(len(xv))) else 0
    return mark


def main():
    rng = random.Random(0)
    print("=" * 74)
    print("ゲート段の 印つき加算器")
    print("=" * 74)

    print("\n① 値が 既存の 符号つき桁 加算と 一致するか")
    bad = n = 0
    for _ in range(400):
        L = 12
        xv = [rng.choice((-1, 0, 1)) for _ in range(L)]
        yv = [rng.choice((-1, 0, 1)) for _ in range(L)]
        X = [g.enc(d) for d in xv]
        Y = [g.enc(d) for d in yv]
        st = g.new_counter()
        s, mk = mark_add(X, Y, st)
        got = g.from_sd(s)
        want = g.from_sd([g.enc(d) for d in xv]) + g.from_sd([g.enc(d) for d in yv])
        n += 1
        if got != want:
            bad += 1
    print(f"   {n} 例中 値が 合わなかった: {bad}")

    print("\n② 印が Python 模型と 一致するか")
    bad2 = n2 = 0
    for _ in range(2000):
        L = 12
        xv = [rng.choice((-1, 0, 1)) for _ in range(L)]
        yv = [rng.choice((-1, 0, 1)) for _ in range(L)]
        X = [g.enc(d) for d in xv]
        Y = [g.enc(d) for d in yv]
        st = g.new_counter()
        s, mk = mark_add(X, Y, st)
        n2 += 1
        if int(mk) != py_model(xv, yv):
            bad2 += 1
    print(f"   {n2} 例中 印が 食い違った: {bad2}")

    print("\n③ ゲート数 — 印の 費用は いくらか")
    print(f"   {'桁数':>6} {'値の経路':>10} {'印の経路':>10} {'合計':>10} {'印の割合':>10}")
    for L in (4, 8, 12, 16, 24):
        xv = [rng.choice((-1, 0, 1)) for _ in range(L)]
        yv = [rng.choice((-1, 0, 1)) for _ in range(L)]
        X = [g.enc(d) for d in xv]; Y = [g.enc(d) for d in yv]
        st1 = g.new_counter(); g.sd_sum([list(X), list(Y)], st1)
        v = sum(st1.values())
        st2 = g.new_counter()
        marks = [cancel_bit(X[i], Y[i], st2) for i in range(L)]
        or_tree(marks, st2)
        mk = sum(st2.values())
        print(f"   {L:>6} {v:>10,} {mk:>10,} {v + mk:>10,} {mk / (v + mk):>9.1%}")
    print("\n   ⟹ 相殺検出は 3 ゲート/桁 + OR 木。値の 経路に比べて 小さい")


if __name__ == "__main__":
    main()
