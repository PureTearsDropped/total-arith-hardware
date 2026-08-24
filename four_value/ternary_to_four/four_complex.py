#!/usr/bin/env python3
"""複素数の 代数を 4 値桁で — 配線を 変えると `?` の 伝わり方は 変わるか。

`bilinear_unit` は `c = W·((U·a)⊙(V·b))`。同じ 複素代数を 2 通りの 配線で 計算する:

  素朴 R=4:  p = [a0b0, a1b1, a0b1, a1b0]
             Re = p0 − p1、Im = p2 + p3
  ガウス R=3: p = [a0b0, a1b1, (a0+a1)(b0+b1)]
             Re = p0 − p1、Im = p2 − p0 − p1

積の 個数は 4 → 3 に 減るが、**`(a0+a1)` で 不明が 早く 混ざる**。
締まりに 差が 出るか。オラクルは 独立（全具体化して 真の 複素積を 集める）。
"""
import sys, itertools
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__file__), '..', '..'))
import gate_bilinear as g
from gate_compress3_four import ENC4, DEC4, NAMES, QS
from swap_primitives import install, restore, dec, val_lsb, value_set_lsb

# (U, V, W) — 成分は 三値 {−1,0,+1} なので 前後は 配線だけ
NAIVE = ([[1, 0], [0, 1], [1, 0], [0, 1]],
         [[1, 0], [0, 1], [0, 1], [1, 0]],
         [[1, -1, 0, 0], [0, 0, 1, 1]])
GAUSS = ([[1, 0], [0, 1], [1, 1]],
         [[1, 0], [0, 1], [1, 1]],
         [[1, -1, 0], [-1, -1, 1]])


def lincomb_d(coeffs, nums, st):
    """digit 列に 対する 線形結合（±1 は 配線・0 は 捨てる）。"""
    terms = []
    for c, x in zip(coeffs, nums):
        if c == 0:
            continue
        terms.append(x if c == 1 else g.neg(x))
    return g.sd_sum(terms, st) if terms else [g.ZERO]


def bilinear_d(UVW, A, B, st):
    """digit 列 (4 値) を 直に 受ける 双線形ユニット。"""
    U, V, W = UVW
    R = len(U)
    left = [lincomb_d(U[r], A, st) for r in range(R)]
    right = [lincomb_d(V[r], B, st) for r in range(R)]
    prod = [g.multiply(left[r], right[r], st) for r in range(R)]
    return [lincomb_d(W[k], prod, st) for k in range(len(W))]


def complex_truth(a_words, b_words):
    """真の 複素積の 集合（成分ごと）。全具体化して 集める。"""
    re_set, im_set = set(), set()
    a0s, a1s = value_set_lsb(a_words[0]), value_set_lsb(a_words[1])
    b0s, b1s = value_set_lsb(b_words[0]), value_set_lsb(b_words[1])
    for a0 in a0s:
        for a1 in a1s:
            for b0 in b0s:
                for b1 in b1s:
                    re_set.add(a0 * b0 - a1 * b1)
                    im_set.add(a0 * b1 + a1 * b0)
    return frozenset(re_set), frozenset(im_set)


def run(name, UVW, n, patterns):
    bad = tot = 0
    rr, ri = [], []
    for a0 in patterns:
        for a1 in patterns:
            for b0 in patterns:
                for b1 in patterns:
                    st = g.new_counter()
                    A = [[ENC4[d] for d in a0], [ENC4[d] for d in a1]]
                    B = [[ENC4[d] for d in b0], [ENC4[d] for d in b1]]
                    out = bilinear_d(UVW, A, B, st)
                    got_re = value_set_lsb([dec(d) for d in out[0]])
                    got_im = value_set_lsb([dec(d) for d in out[1]])
                    t_re, t_im = complex_truth([a0, a1], [b0, b1])
                    tot += 1
                    if not (t_re <= got_re and t_im <= got_im):
                        bad += 1
                    if t_re:
                        rr.append(len(got_re) / len(t_re))
                    if t_im:
                        ri.append(len(got_im) / len(t_im))
    rr.sort(); ri.sort()
    print(f"   {name:>10} R={len(UVW[0])}  {tot:>5} 組  取りこぼし {bad:>3}  "
          f"締まり Re 中央 {rr[len(rr)//2]:>5.2f} / Im 中央 {ri[len(ri)//2]:>5.2f}  "
          f"最大 {max(rr[-1], ri[-1]):>6.2f}")
    return bad


def main():
    print("=" * 78)
    print("複素数の 代数を 4 値桁で — 配線 2 通り")
    print("=" * 78)
    install()
    try:
        print("\n① 三値だけの 入力（後方互換・厳密であるべき）")
        pats3 = list(itertools.product((-1, 0, 1), repeat=2))
        for name, UVW in (("素朴", NAIVE), ("ガウス", GAUSS)):
            run(name, UVW, 2, pats3)

        print("\n② `?` を 1 桁 含む 入力")
        pats = [(-1, 0), (1, 0), (0, 1), ('?', 0), (0, '?'), (1, '?')]
        for name, UVW in (("素朴", NAIVE), ("ガウス", GAUSS)):
            run(name, UVW, 2, pats)

        print("\n③ `?` を 多く 含む 入力")
        pats = [('?', '?'), ('?', 0), (1, '?'), (0, 0), (1, 1)]
        for name, UVW in (("素朴", NAIVE), ("ガウス", GAUSS)):
            run(name, UVW, 2, pats)

        print("\n④ ゲート数の 比較（同じ 入力で）")
        for name, UVW in (("素朴", NAIVE), ("ガウス", GAUSS)):
            st = g.new_counter()
            A = [[ENC4[d] for d in (1, '?')], [ENC4[d] for d in (0, 1)]]
            B = [[ENC4[d] for d in (1, 1)], [ENC4[d] for d in (-1, 0)]]
            bilinear_d(UVW, A, B, st)
            print(f"   {name:>10} R={len(UVW[0])}: {sum(st.values()):>6,} ゲート")
    finally:
        restore()


if __name__ == "__main__":
    main()
