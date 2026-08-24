#!/usr/bin/env python3
"""溢れた語の 演算 — 「MAX 以上」「−MAX 以下」は 半無限区間。区間演算が そのまま 効く。

**動機** (2026-08-25 ユーザ): 溢れは 「以上/以下」なのだから 演算できるはず。
現行の `tot_add` は 溢れが 絡むと ほぼ 一律 `SUNK` に 落とすが、区間で 見ると 分かれる:

    [MAX, ∞) + 有限        = [MAX+有限, ∞)    まだ 溢れている・**向きは 分かる**
    [MAX, ∞) + [MAX, ∞)    = [2MAX, ∞)        もっと 溢れている・向きは 分かる
    [MAX, ∞) + (−∞, −MAX]  = (−∞, ∞)          **本当に 何も 言えない** ← これだけが SUNK

さらに 今日 確認した 「連なりが 場の端まで 届けば 発散」を 使うと、**有限の 溢れ**
(上に 0 が ある = k 桁ぶん) と **無限の 溢れ** (端まで 届く) が 区別できる。
有限なら 区間は 有限幅で、加算しても 有限の まま。

このファイルは 溢れ込みの 加算を 区間で 行い、結果を 語に 戻す。
"""
from total_word import Word, read_tail, enc_head, dec_head, F_GE, F_LE, F_SUNK

INF = float('inf')


def head_is_saturated(head):
    """連なりが 場の 端まで 届いているか (= 上へ 拡張すると 発散 = 無限の 溢れ)。"""
    return len(head) > 0 and head[0] != 0


def word_interval(w: Word):
    """溢れを 含めた 真値の 区間。無限の 溢れは 半無限に なる。"""
    o, k, s = dec_head(w.head)
    if o == 'none':
        return w.interval()
    C = w.ceiling()
    if o == 'unknown':
        hi = INF if head_is_saturated(w.head) else C * (2 ** k)
        return -hi, hi
    lo_m = C * (2 ** (k - 1))
    hi_m = INF if head_is_saturated(w.head) else C * (2 ** k)
    return (lo_m, hi_m) if s > 0 else (-hi_m, -lo_m)


def classify(lo, hi):
    """区間 → 溢れの 状態 (語に 戻すため)。"""
    if lo == -INF and hi == INF:
        return 'unknown', 0
    if hi == INF:
        return 'over+', 0
    if lo == -INF:
        return 'over-', 0
    return 'finite', 0


def add_iv(a: Word, b: Word):
    """溢れ込みの 加算 — 区間で 足すだけ。戻り (区間, 状態名)。"""
    alo, ahi = word_interval(a)
    blo, bhi = word_interval(b)
    lo = -INF if (alo == -INF or blo == -INF) else alo + blo
    hi = INF if (ahi == INF or bhi == INF) else ahi + bhi
    return (lo, hi), classify(lo, hi)[0]


def legacy_says(a: Word, b: Word):
    """現行 `tot_add` が この 組で 何を 言うか (規則を 再現)。

    同符号なら 単純和が 健全 / 相殺しうるなら 境界を 落として SUNK。
    """
    fa, fb = a.to_legacy_flag(), b.to_legacy_flag()
    va, vb = a.nominal(), b.nominal()
    same = (va > 0) == (vb > 0) and va != 0 and vb != 0
    if (fa | fb) == 0:
        return 'exact'
    if same:
        if (fa & F_GE) and (fb & F_GE) and not (fa & F_LE) and not (fb & F_LE):
            return 'GE'
        if (fa & F_LE) and (fb & F_LE) and not (fa & F_GE) and not (fb & F_GE):
            return 'LE'
    return 'none+SUNK'


def make_over(mant, h, k, s, sat, E=0, tail=None):
    """溢れた語を 作る。sat=True なら 連なりが 端まで 届く (無限の 溢れ)。

    **不変条件**: 溢れの 向き s は **仮数の 符号と 一致**しなければならない
    (`_sat` は 符号を 保つので 正の 数は 正側にしか 溢れない)。初版の テストは
    仮数を 正の まま 先頭だけ 負に して 矛盾した 語を 作っていた — 2026-08-25。
    """
    assert mant, "仮数が 空"
    msign = 1 if mant[0] > 0 else -1
    assert msign == s, f"溢れの 向き {s} が 仮数の 符号 {msign} と 矛盾"
    head = [0] * (h - k) + [s] * k if not sat else [s] * h
    return Word(head, mant, tail or [], E)


def main():
    W, H = 8, 3
    full = [1] * W                      # 255
    print("=" * 76)
    print("溢れ込みの 加算 — 区間で 足すと 何が 区別できるか")
    print("=" * 76)

    A_inf = make_over(full, H, H, +1, sat=True)      # [MAX, ∞)
    B_inf = make_over([-1] * W, H, H, -1, sat=True)  # (−∞, −MAX]  仮数も 負
    A_fin = make_over(full, H, 2, +1, sat=False)     # 有限の 溢れ (2 桁ぶん)
    C_fin = Word(enc_head('none', H), [1, 0, 0, 0, 0, 0, 0, 0], [])   # 普通の 128

    cases = [
        ("無限の溢れ(+) + 普通の値", A_inf, C_fin),
        ("無限の溢れ(+) + 無限の溢れ(+)", A_inf, A_inf),
        ("無限の溢れ(+) + 無限の溢れ(−)", A_inf, B_inf),
        ("有限の溢れ(+) + 普通の値", A_fin, C_fin),
        ("有限の溢れ(+) + 有限の溢れ(+)", A_fin, A_fin),
    ]
    print(f"\n{'場面':>30} {'区間で足した結果':>26} {'判定':>9} | {'現行':>11}")
    for lab, x, y in cases:
        (lo, hi), st = add_iv(x, y)
        f = lambda v: ("∞" if v == INF else ("−∞" if v == -INF else f"{v:g}"))
        print(f"{lab:>30} [{f(lo):>10},{f(hi):>10}] {st:>9} | {legacy_says(x, y):>11}")

    print("\n⟹ 現行は どれも 境界なし+SUNK に 落とすが、区間なら")
    print("   ・上3件の うち **向きが 分かる のは 2 件**、本当に 何も 言えないのは 1 件だけ")
    print("   ・有限の 溢れ どうしは **有限の まま**（区間幅が 残る）")

    print("\n" + "=" * 76)
    print("健全性: 区間の 和は 取りこぼさないか (端点を 総当たり)")
    print("=" * 76)
    bad = n = 0
    words = [A_inf, B_inf, A_fin, C_fin,
             make_over([-1] * W, H, 1, -1, sat=False),
             Word(enc_head('none', H), [1, 1, 0, 0, 0, 0, 0, 1], [-1, -1])]
    for x in words:
        for y in words:
            (lo, hi), _ = add_iv(x, y)
            xl, xh = word_interval(x); yl, yh = word_interval(y)
            for u in (xl, xh):
                for v in (yl, yh):
                    if u in (INF, -INF) or v in (INF, -INF):
                        continue
                    n += 1
                    if not (lo - 1e-9 <= u + v <= hi + 1e-9):
                        bad += 1
    print(f"  端点の 組 {n} 通り中 取りこぼし {bad}")


if __name__ == "__main__":
    main()
