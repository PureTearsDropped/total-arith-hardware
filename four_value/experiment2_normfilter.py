#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証
"""4値桁研究・実験2 — ノルム多項式の網（±2±5 の思考実験の検証と限界）。

構成(古典: レゾルベント/ノルム形式): x = Σ sᵢaᵢ (sᵢ∈{±1}不明) の取りうる値の集合は
P(x) = ∏(x − Σsᵢaᵢ) の根。倍化合成 P_{S∪{a}}(x) = P_S(x−a)·P_S(x+a) により
**符号を一度も列挙せずに** 係数が決まる(aᵢ²だけの関数)。?桁は x³−x の合成。

結果:
  ① ±2±5: 幻影(±1,±5)を残差±384で棄却・本物(±3,±7)のみ通過
  ② n=5(次数32)まで 偽陽性0 — 網に穴なし
  ③ ?×m桁(次数3^m)の根 = **ちょうど連続区間** ⟹ ?尻尾には多項式は冗長(区間で足りる)
  ④ ルーズな桁独立集合(可能性キャリーのスープ)を網で濾すと真集合と完全一致

分業法則: 密な切り捨て(?尻尾)=区間が厳密で最安 / **疎な符号不明(SUNK)=真集合に隙間
→ 多項式の網だけが厳密に運べる**(コスト2^n) / フラグ=両者の粗い影。
SUNK ≡ x²−v² の根・?桁 ≡ x³−x の根・倍化合成の実体=畳み込み・評価=powersum。
"""
import itertools

def polymul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i + j] += a * b
    return r

def shift(p, c):
    out = [0] * len(p)
    for k, pk in enumerate(p):
        if pk == 0: continue
        row = [1]
        for _ in range(k):
            row = polymul(row, [-c, 1])
        for i, r in enumerate(row):
            out[i] += pk * r
    return out

def compose_pm(terms):
    P = [0, 1]
    for a in terms:
        P = polymul(shift(P, a), shift(P, -a))
    return P

def compose_trit(P, w):
    return polymul(polymul(shift(P, w), P), shift(P, -w))

def peval(p, x):
    v = 0
    for c in reversed(p): v = v * x + c
    return v

def main():
    P = compose_pm([2, 5])
    assert P == [441, 0, -58, 0, 1]
    ts = sorted({s1 * 2 + s2 * 5 for s1 in (-1, 1) for s2 in (-1, 1)})
    assert all(peval(P, x) == 0 for x in ts)
    assert all(peval(P, x) != 0 for x in range(-8, 9) if x not in ts)
    print("① ±2±5: P=x⁴−58x²+441・幻影全棄却・本物全通過 ✓")

    terms = [2, 5, 11, 23, 47]
    for n in range(1, 6):
        Pn = compose_pm(terms[:n])
        tsn = sorted({sum(s * a for s, a in zip(ss, terms[:n]))
                      for ss in itertools.product((-1, 1), repeat=n)})
        assert all(peval(Pn, x) == 0 for x in tsn)
        assert not [x for x in range(min(tsn), max(tsn) + 1)
                    if x not in tsn and peval(Pn, x) == 0]
    print("② n=1..5(次数2..32): 全通過・偽陽性0 ✓")

    for m in (1, 2, 3):
        Pt = [0, 1]
        for i in range(m):
            Pt = compose_trit(Pt, 1 << i)
        tst = sorted({sum(d << i for i, d in enumerate(ds))
                      for ds in itertools.product((-1, 0, 1), repeat=m)})
        assert tst == list(range(-(2 ** m - 1), 2 ** m))       # 連続区間(実験0)
        assert all(peval(Pt, x) == 0 for x in tst)
        assert not [x for x in range(tst[0] - 2, tst[-1] + 3)
                    if x not in tst and peval(Pt, x) == 0]
    print("③ ?×m(次数3^m): 根=連続区間 ⟹ ?尻尾に多項式は冗長 ✓")

    kept = [x for x in range(-7, 8) if peval(P, x) == 0]
    assert kept == ts
    print("④ 幻影入りスープ(−7..7) → 網 → 真集合と一致 ✓")
    print("結論: 疎な符号不明の厳密台帳=ノルム多項式(2^n)・密な切り捨て=区間・フラグ=影。")

if __name__ == "__main__":
    main()
