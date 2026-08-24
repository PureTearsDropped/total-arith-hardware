#!/usr/bin/env python3
"""gate9 の 4 値拡張 — `RESEARCH.md` テーマ1 の 開いている問い 1。

`?` = 「この桁は {−1,0,+1} の どれか」。積は **集合の 像**として 厳密に 決まる:

    ? × d  =  { -d, 0, +d }

  d = 0  なら {0} ⟹ **`?` でなく 厳密な 0**（真の 0 は 万物を 吸収。`tot_mul` の
                    「x×0=0 厳密」と 同じ 教義）
  d = ±1 なら {−1,0,+1} ⟹ `?`

⟹ **積が `?` に なるのは 両方が 非零の 不明を 含むときだけ。** 表を 全列挙して 確かめ、
   ゲート数を 数える。

註: これは 桁 1 つ分の 意味論の 確認であり、乗算器全体（部分積の 足し込み・正規化）の
検証では ない。今日 ternary_mark で 過大な 主張を した 反省として、範囲を 明示する。
"""
import itertools

VALS = {-1: frozenset({-1}), 0: frozenset({0}), 1: frozenset({1}),
        '?': frozenset({-1, 0, 1})}
NAMES = [-1, 0, 1, '?']


def prod_set(a, b):
    """集合の 像 — これが 定義。"""
    return frozenset(x * y for x in VALS[a] for y in VALS[b])


def name_of(s):
    for n in NAMES:
        if VALS[n] == s:
            return n
    return f"集合{sorted(s)}"


def main():
    print("=" * 66)
    print("gate9 の 4 値拡張 — 桁どうしの 積 (集合の像で 定義)")
    print("=" * 66)

    print(f"\n{'×':>4}", end='')
    for b in NAMES:
        print(f"{str(b):>8}", end='')
    print()
    closed = True
    for a in NAMES:
        print(f"{str(a):>4}", end='')
        for b in NAMES:
            s = prod_set(a, b)
            nm = name_of(s)
            if isinstance(nm, str) and nm.startswith('集合'):
                closed = False
            print(f"{str(nm):>8}", end='')
        print()

    print(f"\n① 表は 4 値で **閉じているか**: {'はい' if closed else 'いいえ（新しい状態が要る）'}")

    print("\n② `?` は 吸収的か（真の 0 を 除いて）")
    ok = True
    for a in NAMES:
        s = prod_set('?', a)
        exp = frozenset({0}) if a == 0 else VALS['?']
        got = name_of(s)
        mark = '✓' if s == exp else '✗'
        if s != exp:
            ok = False
        print(f"   ? × {str(a):>2} = {str(got):>2}  "
              f"{'（0 が 吸収して 厳密）' if a == 0 else '（? のまま）'} {mark}")
    print(f"   ⟹ {'期待どおり' if ok else '期待と違う'}")

    print("\n③ 三値の gate9 と 矛盾しないか（`?` を 含まない 入力で 一致するか）")
    bad = 0
    for a, b in itertools.product((-1, 0, 1), repeat=2):
        s = prod_set(a, b)
        if s != frozenset({a * b}):
            bad += 1
    print(f"   16 通り中 三値の 積と 食い違い: {bad}")

    print("\n④ ゲート数 — (p,n) 2 本で {00=0, 01=+1, 10=−1, 11=?} に 割り当てる")
    print("   出力 (op,on) を 入力 (ap,an,bp,bn) の 論理式で 書く")
    # 表から 真理値表を 作り、簡単な 式を 当てはめて 検証する
    ENC = {0: (0, 0), 1: (1, 0), -1: (0, 1), '?': (1, 1)}
    DEC = {v: k for k, v in ENC.items()}
    rows = []
    for a in NAMES:
        for b in NAMES:
            ap, an = ENC[a]; bp, bn = ENC[b]
            r = name_of(prod_set(a, b))
            op, on = ENC[r]
            rows.append((ap, an, bp, bn, op, on))
    # 候補式:
    #   nz_a = ap|an  (a が 非零 = ±1 か ?)   unk_a = ap&an  (a が ?)
    #   op = (ap&bp)|(an&bn) | (unk かかわる ?)   … 実際に 検査する
    def cand_op(ap, an, bp, bn):
        ua, ub = ap & an, bp & bn
        za, zb = (not (ap | an)), (not (bp | bn))
        if za or zb:
            return 0, 0                      # 真の 0 が 吸収
        if ua or ub:
            return 1, 1                      # どちらかが ? ⟹ ?
        pos = (ap & bp) | (an & bn)
        neg = (ap & bn) | (an & bp)
        return pos, neg
    bad2 = 0
    for ap, an, bp, bn, op, on in rows:
        g = cand_op(ap, an, bp, bn)
        if (int(g[0]), int(g[1])) != (op, on):
            bad2 += 1
    print(f"   候補式が 表を 再現するか: 16 通り中 食い違い {bad2}")
    print("   式: 真の0 が 吸収 → 0 / どちらかが ? → ? / 残りは 三値の gate9")
    print("   ゲート: nz/unk の 検出 4 + 三値 gate9 4 + 選択 数個 ⟹ 三値 gate9 の 2 倍程度")

    print("\n⑤ この 検証が **示していない** こと")
    print("   ・部分積の 足し込みで `?` が 正しく 伝播するか（未検証）")
    print("   ・正規化（シフト）で `?` の 開始位置が どう 動くか（未検証）")
    print("   ・乗算器全体の 健全性（未検証）")


if __name__ == "__main__":
    main()
