#!/usr/bin/env python3
"""4 値 3:2 圧縮器を **ゲート段**で — 既存の compress3 を そのまま 包む。

`sd_sum` は 列の 圧縮も リップルも **compress3 だけ**で 組まれ、`multiply` は
`sd_sum` + `gate9`。⟹ **この 2 つを 4 値版に 差し替えれば 上の 層は 無改造**。

**包み方**（今日 導出した 規則を 既存部品で 書き直したもの）:
  unk_x = xp & xn                       … その桁が `?`（(1,1)）か
  マスク: xp' = xp ^ unk_x, xn' = xn ^ unk_x   … `?` を 0 に する（配線 2 ゲート）
  (lo3, hi3) = compress3(マスク後の 3 つ)      … **既存の 18 ゲートを そのまま**
  nq = unk_a + unk_b + unk_c
    nq = 0 → (lo3, hi3)
    nq = 1 → low = ?、high = hi3 if lo3 == 0 else ?
              （lo3 == 0 ⟺ 他 2 つの 和 t2 が 偶数。t2 = lo3 + 2·hi3 で lo3∈{−1,0,1}）
    nq ≥ 2 → (?, ?)

これが 導出表（64 通り・健全かつ最小）と 一致するかを 全列挙で 確かめる。
"""
import sys, itertools
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__file__), '..', '..'))
import gate_bilinear as g

# 差し替え（monkeypatch）で 自分自身を 呼ばないよう、**元の関数を 直接 掴んでおく**。
# g.compress3 を 参照すると 差し替え後に 無限再帰する（2026-08-25 に 踏んだ）。
_TERN_COMPRESS3 = g.compress3

QS = {-1: (-1,), 0: (0,), 1: (1,), '?': (-1, 0, 1)}
NAMES = (-1, 0, 1, '?')
ENC4 = {0: (0, 0), 1: (1, 0), -1: (0, 1), '?': (1, 1)}
DEC4 = {v: k for k, v in ENC4.items()}


def compress3_4_gates(x, y, z, st):
    """(p,n) レールで 4 値 3:2 圧縮。戻り ((lp,ln), (hp,hn))。"""
    def unk(d):
        return g.AND(d[0], d[1], st)

    ua, ub, uc = unk(x), unk(y), unk(z)

    def mask(d, u):
        return (g.XOR(d[0], u, st), g.XOR(d[1], u, st))

    xm, ym, zm = mask(x, ua), mask(y, ub), mask(z, uc)
    (lp, ln), (hp, hn) = _TERN_COMPRESS3(xm, ym, zm, st)  # ← 既存を そのまま

    nq1 = g.OR(g.OR(ua, ub, st), uc, st)                  # ? が 1 個以上
    m2 = g.OR(g.OR(g.AND(ua, ub, st), g.AND(ua, uc, st), st),
              g.AND(ub, uc, st), st)                      # ? が 2 個以上
    lo_nz = g.OR(lp, ln, st)                              # lo3 ≠ 0
    hi_unk = g.OR(m2, g.AND(nq1, lo_nz, st), st)          # high を ? に する 条件

    out_lp = g.OR(lp, nq1, st); out_ln = g.OR(ln, nq1, st)
    out_hp = g.OR(hp, hi_unk, st); out_hn = g.OR(hn, hi_unk, st)
    return (out_lp, out_ln), (out_hp, out_hn)


# ---------------------------------------------------------------- 導出表（対照）
# 三値の 分解は **冗長**（t=±1 は 2 通り）。参照は 既存 compress3 の 選択に 合わせる —
# 私の 勝手な 表と 比べても 意味が ない（2026-08-25: それで 6 件の 見かけの 食い違いが 出た）。
def _tern_from_compress3(a, b, c):
    st = g.new_counter()
    (lp, ln), (hp, hn) = _TERN_COMPRESS3(g.enc(a), g.enc(b), g.enc(c), st)
    return DEC4[(int(lp), int(ln))], DEC4[(int(hp), int(hn))]


def ref_rule(a, b, c):
    unk = [d == '?' for d in (a, b, c)]
    nq = sum(unk)
    if nq >= 2:
        return '?', '?'
    if nq == 1:
        t2 = sum(d for d in (a, b, c) if d != '?')
        return ('?', t2 // 2) if t2 % 2 == 0 else ('?', '?')
    return _tern_from_compress3(a, b, c)


def repr_set(lo, hi):
    return frozenset(l + 2 * h for l in QS[lo] for h in QS[hi])


def main():
    print("=" * 70)
    print("4 値 3:2 圧縮器（ゲート段）— 既存 compress3 を包む")
    print("=" * 70)

    print("\n① ゲート段が 導出表と 一致するか（64 通り 全列挙）")
    bad = 0
    for a, b, c in itertools.product(NAMES, repeat=3):
        st = g.new_counter()
        (lp, ln), (hp, hn) = compress3_4_gates(ENC4[a], ENC4[b], ENC4[c], st)
        got = (DEC4[(int(lp), int(ln))], DEC4[(int(hp), int(hn))])
        want = ref_rule(a, b, c)
        if got != want:
            bad += 1
            if bad <= 3:
                print(f"   ✗ {a},{b},{c}: ゲート={got} 表={want}")
    print(f"   64 通り中 食い違い: {bad}")

    print("\n② 健全性（真の 和 ⊆ 出力が 表す 集合）")
    bad2 = 0
    for a, b, c in itertools.product(NAMES, repeat=3):
        t = frozenset(x + y + z for x in QS[a] for y in QS[b] for z in QS[c])
        st = g.new_counter()
        (lp, ln), (hp, hn) = compress3_4_gates(ENC4[a], ENC4[b], ENC4[c], st)
        s = repr_set(DEC4[(int(lp), int(ln))], DEC4[(int(hp), int(hn))])
        if not t <= s:
            bad2 += 1
    print(f"   64 通り中 取りこぼし: {bad2}")

    print("\n③ ゲート数")
    st = g.new_counter()
    _TERN_COMPRESS3(g.enc(1), g.enc(-1), g.enc(0), st)
    tern = sum(st.values())
    st2 = g.new_counter()
    compress3_4_gates(ENC4[1], ENC4[-1], ENC4['?'], st2)
    four = sum(st2.values())
    print(f"   三値 compress3 : {tern:>3} ゲート")
    print(f"   4 値 compress3 : {four:>3} ゲート  = {four/tern:.2f} 倍")
    print(f"   内訳の 追加分  : {four - tern:>3} ゲート")


if __name__ == "__main__":
    main()
