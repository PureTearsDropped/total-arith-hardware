#!/usr/bin/env python3
"""**素子を差し替えるだけで上の層は動くか** — ユーザの見立ての検証。

`sd_sum` は 列の 圧縮も リップルも `compress3` だけで 組まれ、`multiply` は
`sd_sum` + `gate9`。⟹ この 2 つを 4 値版に すれば 上は 無改造の はず。

ここでは `gate_bilinear` の `compress3` と `gate9` を 4 値版に **差し替え**、
`sd_sum` / `multiply` を **一行も 触らずに** 呼んで、独立オラクルで 検査する。

オラクル: 入力を 全具体化 → 真の 和/積の 集合 → 出力が 表す 集合に 含まれるか。
"""
import sys, itertools
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__file__), '..', '..'))
import gate_bilinear as g
from gate_compress3_four import compress3_4_gates, ENC4, DEC4, NAMES, QS

_orig_compress3 = g.compress3
_orig_gate9 = g.gate9


def gate9_4_gates(x, y, st):
    """4 値の 桁積。真の 0 が 吸収 → 0 / どちらかが ? → ? / 残りは 三値 gate9。"""
    xp, xn = x
    yp, yn = y
    ux = g.AND(xp, xn, st)
    uy = g.AND(yp, yn, st)
    zx = g.NOT(g.OR(xp, xn, st), st)        # x が 真の 0
    zy = g.NOT(g.OR(yp, yn, st), st)
    anyz = g.OR(zx, zy, st)
    anyu = g.OR(ux, uy, st)
    # 三値 gate9 は ? を 0 に マスクして から
    xm = (g.XOR(xp, ux, st), g.XOR(xn, ux, st))
    ym = (g.XOR(yp, uy, st), g.XOR(yn, uy, st))
    tp, tn = _orig_gate9(xm, ym, st)
    unk = g.AND(anyu, g.NOT(anyz, st), st)   # ? かつ 真の 0 でない ⟹ ?
    op = g.AND(g.OR(tp, unk, st), g.NOT(anyz, st), st)
    on = g.AND(g.OR(tn, unk, st), g.NOT(anyz, st), st)
    return (op, on)


def install():
    g.compress3 = compress3_4_gates
    g.gate9 = gate9_4_gates


def restore():
    g.compress3 = _orig_compress3
    g.gate9 = _orig_gate9


def dec(d):
    return DEC4[(int(d[0]), int(d[1]))]


def val_lsb(digits):
    """低位から の 桁列（gate_bilinear の 慣習）→ 値。"""
    return sum(d * (1 << i) for i, d in enumerate(digits))


def value_set_lsb(word):
    return frozenset(val_lsb(c) for c in itertools.product(*[QS[d] for d in word]))


def main():
    print("=" * 72)
    print("素子を差し替えるだけで上の層は動くか")
    print("=" * 72)

    install()
    try:
        print("\n① sd_sum（無改造）を 4 値入力で — 加算の 健全性")
        bad = tot = 0
        ratios = []
        for n in (2, 3):
            for a in itertools.product(NAMES, repeat=n):
                for b in itertools.product(NAMES, repeat=n):
                    st = g.new_counter()
                    A = [ENC4[d] for d in a]
                    B = [ENC4[d] for d in b]
                    Z = g.sd_sum([A, B], st)
                    out = [dec(d) for d in Z]
                    truth = frozenset(x + y for x in value_set_lsb(a)
                                      for y in value_set_lsb(b))
                    got = value_set_lsb(out)
                    tot += 1
                    if not truth <= got:
                        bad += 1
                    if truth:
                        ratios.append(len(got) / len(truth))
            ratios.sort()
            print(f"   桁数 {n}: {tot} 組  取りこぼし {bad}  "
                  f"締まり 中央値 {ratios[len(ratios)//2]:.2f} / 最大 {ratios[-1]:.2f}")
            bad = tot = 0
            ratios = []

        print("\n② multiply（無改造）を 4 値入力で — 乗算の 健全性")
        for n in (2, 3):
            bad = tot = 0
            ratios = []
            for a in itertools.product(NAMES, repeat=n):
                for b in itertools.product(NAMES, repeat=n):
                    st = g.new_counter()
                    A = [ENC4[d] for d in a]
                    B = [ENC4[d] for d in b]
                    Z = g.multiply(A, B, st)
                    out = [dec(d) for d in Z]
                    truth = frozenset(x * y for x in value_set_lsb(a)
                                      for y in value_set_lsb(b))
                    got = value_set_lsb(out)
                    tot += 1
                    if not truth <= got:
                        bad += 1
                    if truth:
                        ratios.append(len(got) / len(truth))
            ratios.sort()
            print(f"   桁数 {n}: {tot} 組  取りこぼし {bad}  "
                  f"締まり 中央値 {ratios[len(ratios)//2]:.2f} / 最大 {ratios[-1]:.2f}")

        print("\n③ 三値だけの 入力で 既存と 一致するか（後方互換）")
        bad3 = n3 = 0
        for a in itertools.product((-1, 0, 1), repeat=3):
            for b in itertools.product((-1, 0, 1), repeat=3):
                st = g.new_counter()
                got = val_lsb([dec(d) for d in g.multiply(
                    [ENC4[x] for x in a], [ENC4[y] for y in b], st)])
                want = val_lsb(list(a)) * val_lsb(list(b))
                n3 += 1
                if got != want:
                    bad3 += 1
        print(f"   {n3} 組中 値が 食い違い: {bad3}")
    finally:
        restore()


if __name__ == "__main__":
    main()
