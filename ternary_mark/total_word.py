#!/usr/bin/env python3
"""全域語 (TotalWord) — フラグを 別ビットで 持たず、**桁の 並びだけ**で 状態を 表す 符号化。

2026-08-24〜25 の 設計 (発案は 全て ユーザ)。算術は まだ 載せない。**表現を 確定させる**
ことだけが 目的で、往復・不変条件・区間の 健全性・現行系への 射影を 検査する。

    語 = [ 溢れの 領域 h 桁 ][ 正規化された 仮数 W 桁 ][ 誤差の 領域 m 桁 ]  × 2^E

**桁は 三値 {−1, 0, +1} のみ** — 配線正規形を 変えない。フラグ用の ビットは **0 本**。

## 下端: 誤差の 領域 — **印ではなく 値**
初版は 末尾を 「連なりの パターン」= 状態の 印 として 読んでいた。**それは 情報の 捨てすぎ**
だった (2026-08-25 に 全列挙で 判明)。末尾は 半分選び定理の 通り **入れ子の 二分**で 読む:

    初期区間 = ±1 ULP。上の 桁から  −1 → 下半分 / 0 → 中央 / +1 → 上半分

  - 幅は **深さ m だけ**で 決まる (2^-m)。桁は **位置**を 決める
  - 全列挙で 「二分の 区間 ⊇ 桁の 算術値」を 確認 (深さ 5 まで 反例 0)
  - **全 3^m 個の 語が 有効**。初版で `invalid` と していた 混在パターンこそ 情報が 多い

**向きは 初版と 逆だった**: 末尾の 桁は 誤差そのものなので、全部 −1 は 誤差が 最も負
⟹ 真値 < 表示値 ⟹ **LE** (初版は これを GE と していた)。桁の 算術が 正典。

    全て −1 → 誤差区間が 全て 負 → 表示値は 上端 → LE
    全て +1 → 誤差区間が 全て 正 → 表示値は 下端 → GE
    全て  0 → 誤差区間が 0 を またぐ → **中央** (劣化ではない・幅は 他と 同一)

**`GE|LE` の 意味が 現行と 違う**: 現行は「境界なし = 最悪」。ここでは「誤差棒が 0 を
またぐ = 中央にいる」で、**幅は GE/LE の 場合と 同じ**。深さ 4 の 全 81 語で 幅は 一定
(GE 40 / LE 40 / 中央 1)。

## 上端: 溢れの 領域 — こちらは 印
正規化により 仮数の 最上位は 必ず ±1。ゆえに その 上の 桁の 0 は「溢れていない」印として
空いている。上から 最初に ±1 が 出る 位置で **何桁ぶん 溢れたか (k)**、符号で **向き**。
連なりの 端が 逆なら 向き不明。溢れの 向きは 値の 符号と 同じ (`_sat` は 符号を 保つ)
ので 別に 持たなくてよい ⟹ 3 状態。

現行の `±MAX + GE` は 溢れの **大きさを 捨てる** (真値 2·MAX も 2^50·MAX も 同じ)。
k を 残せば 区間が [MAX, ∞) から **有限幅**に なる。

**上下は 完全な 対称では ない** (初版の 主張を 訂正): 下端は 値の 精密化、上端は 印。
上端も 値 (溢れた 分の 上位桁) として 読めるかは **未検討**。

## 現行系との 関係
GE/LE/SUNK は この 語からの **射影**であって 原始的でない。SUNK も 「区間が 0 を またぐか」
から 導ける。射影は 情報を 落とすが 嘘は つかない (self_test が 検査)。
"""
import itertools

INF = float('inf')

OVF = ('none', 'over', 'unknown')          # 上端の 3 状態 (向きは 値の 符号)
F_GE, F_LE, F_SUNK = 1, 2, 4               # 現行系の ビット (cuda_total と 同じ)


# ---------------------------------------------------------------- 下端: 誤差 (値として読む)
def read_tail(tail):
    """末尾 m 桁 → 誤差の 区間 (ULP 単位)。入れ子の 二分。

    初期区間 ±1 ULP から 始め、上の 桁から 下半分/中央/上半分 と 絞る。
    戻り (lo, hi) は ULP 単位の 誤差範囲。m=0 なら (0,0) = 厳密。
    """
    m = len(tail)
    if m == 0:
        return 0.0, 0.0
    lo, hi = -1.0, 1.0
    for d in tail:
        c = (lo + hi) / 2.0
        q = (hi - lo) / 4.0
        if d == -1:
            lo, hi = lo, c
        elif d == 1:
            lo, hi = c, hi
        else:
            lo, hi = c - q, c + q
        # 上下どの枝も 幅は ちょうど 半分に なる
    return lo, hi


def tail_state(tail, sign=1):
    """誤差区間と **仮数の符号**から フラグを 導出する。

    現行の GE/LE は **絶対値の 意味** (|真| ≥ |表示| など)。ゆえに 誤差の 符号だけでは
    決まらない — 負の 値に 正の 誤差が 乗ると 大きさは **減る**。
    (2026-08-25: これを 落として self_test の 射影検査に 26 件の 嘘が 出た。)

        v > 0: 誤差 ≥ 0 → GE / 誤差 ≤ 0 → LE
        v < 0: 誤差 ≥ 0 → LE / 誤差 ≤ 0 → GE
    """
    lo, hi = read_tail(tail)
    if lo == hi == 0.0:
        return 'exact'
    if lo >= 0.0:
        return 'GE' if sign >= 0 else 'LE'
    if hi <= 0.0:
        return 'LE' if sign >= 0 else 'GE'
    return 'center'        # 0 を またぐ (劣化ではない)


# ---------------------------------------------------------------- 上端: 溢れ (印として読む)
def enc_head(state, h, k=0, s=1):
    if state == 'none':
        return [0] * h
    assert 1 <= k <= h, "溢れの 深さが 領域に 収まらない"
    if state == 'over':
        return [0] * (h - k) + [s] * k
    assert k >= 2, "向き不明は 深さ 2 以上 (印に 1 桁)"
    return [0] * (h - k) + [-s] + [s] * (k - 1)


def dec_head(head):
    h = len(head)
    if h == 0 or all(d == 0 for d in head):
        return 'none', 0, 0
    i = next(j for j, d in enumerate(head) if d != 0)
    k, s = h - i, head[i]
    body = head[i + 1:]
    if body and all(d == -s for d in body):
        return 'unknown', k, 0
    if all(d == s for d in body):
        return 'over', k, s
    return 'invalid', 0, 0


# ---------------------------------------------------------------- 語
class Word:
    """[head][mant][tail] × 2^E。mant は 正規化済み (mant[0] != 0)。"""
    __slots__ = ('head', 'mant', 'tail', 'E')

    def __init__(self, head, mant, tail, E=0):
        self.head, self.mant, self.tail, self.E = list(head), list(mant), list(tail), int(E)

    def is_ternary(self):
        return all(d in (-1, 0, 1) for d in self.head + self.mant + self.tail)

    def is_normalized(self):
        return len(self.mant) > 0 and (self.mant[0] != 0 or all(d == 0 for d in self.mant))

    def ulp(self):
        return 2.0 ** self.E

    def nominal(self):
        W = len(self.mant)
        return sum(d * (1 << (W - 1 - i)) for i, d in enumerate(self.mant)) * self.ulp()

    def ceiling(self):
        return ((1 << len(self.mant)) - 1) * self.ulp()

    def state(self):
        o, k, s = dec_head(self.head)
        v = self.nominal()
        return dict(prec=tail_state(self.tail, 1 if v >= 0 else -1), depth=len(self.tail),
                    ovf=o, over_digits=k, direction=s)

    def interval(self):
        """真値の 区間。下端が 誤差棒を、上端が 範囲を 決める。"""
        v = self.nominal()
        elo, ehi = read_tail(self.tail)
        lo, hi = v + elo * self.ulp(), v + ehi * self.ulp()
        o, k, s = dec_head(self.head)
        if o == 'over':
            C = self.ceiling()
            m_lo, m_hi = C * (2 ** (k - 1)), C * (2 ** k)
            lo, hi = (m_lo, m_hi) if s > 0 else (-m_hi, -m_lo)
        elif o == 'unknown':
            C = self.ceiling() * (2 ** k)
            lo, hi = -C, C
        return lo, hi

    def to_legacy_flag(self):
        """現行 GE/LE/SUNK への 射影。

        **末尾のパターンからでなく 区間から 導く** — 区間が 0 を またぐ 場合、大きさは
        0 まで 落ちうるので `|真| ≥ |表示|` (GE) は 成り立たない。末尾だけ 見ていると
        これを 見落とす (2026-08-25: 加算の 検証で 2000 例中 27 件の 嘘)。
        射影は 情報を 落とすが 嘘は つかない。
        """
        st = self.state()
        if st['ovf'] == 'unknown':
            return F_GE | F_LE | F_SUNK
        lo, hi = self.interval()
        v = self.nominal()
        f = 0
        straddles = (lo < 0 < hi) or (v == 0 and (lo < 0 or hi > 0))
        if straddles:
            f |= F_SUNK
        a = abs(v)
        mn = 0.0 if straddles else min(abs(lo), abs(hi))
        mx = max(abs(lo), abs(hi))
        ge = mn >= a - 1e-12          # 大きさは 常に |v| 以上
        le = mx <= a + 1e-12          # 大きさは 常に |v| 以下
        if ge and le:
            pass                      # 厳密 — 旗なし
        elif ge:
            f |= F_GE
        elif le:
            f |= F_LE
        else:
            f |= F_GE | F_LE          # 片側の 境界が 立たない
        return f

    def __repr__(self):
        st = self.state()
        return (f"Word(mant={self.mant}, tail={self.tail}, E={self.E} | "
                f"{st['prec']}/m={st['depth']} {st['ovf']}/k={st['over_digits']})")


# ---------------------------------------------------------------- 検査
def self_test(verbose=True):
    ok = True
    say = print if verbose else (lambda *a, **k: None)
    say("=" * 72)
    say("全域語 符号化 — 自己テスト (末尾は 値として 読む 版)")
    say("=" * 72)

    H = 4
    mants = [[1, 0, 1, 1], [-1, 1, 0, -1], [1, 1, 1, 1]]

    # ① 二分読みの 健全性: 区間 ⊇ 桁の 算術値 (全列挙)
    say("\n① 二分読み ⊇ 桁の 算術値 (全列挙)")
    bad = tot = 0
    for m in range(1, 6):
        for t in itertools.product((-1, 0, 1), repeat=m):
            av = sum(d * (1 << (m - 1 - i)) for i, d in enumerate(t)) / float(1 << m)
            lo, hi = read_tail(list(t))
            tot += 1
            if not (lo - 1e-9 <= av <= hi + 1e-9):
                bad += 1
    say(f"   {tot} 語中 反例 {bad}")
    ok &= (bad == 0)

    # ② 幅は 深さだけで 決まるか
    say("\n② 幅は 深さだけで 決まるか (桁は 位置だけを 変える)")
    for m in range(1, 6):
        ws = {round(hi - lo, 12) for t in itertools.product((-1, 0, 1), repeat=m)
              for lo, hi in [read_tail(list(t))]}
        say(f"   m={m}: 相異なる 幅 {len(ws)} 種 {sorted(ws)} (期待 2^-{m}·2 = {2.0/(1<<m)})")
        ok &= (len(ws) == 1)

    # ③ 導出される 状態の 内訳
    say("\n③ 導出される 状態の 内訳")
    from collections import Counter
    for m in (2, 4):
        c = Counter(tail_state(list(t), 1) for t in itertools.product((-1, 0, 1), repeat=m))
        say(f"   m={m}: {dict(c)}")
    say("   全て −1 → LE / 全て +1 → GE / 全て 0 → center  (初版とは 向きが 逆)")
    ok &= (tail_state([-1] * 4, 1) == 'LE' and tail_state([1] * 4, 1) == 'GE'
           and tail_state([0] * 4, 1) == 'center'
           and tail_state([-1] * 4, -1) == 'GE' and tail_state([1] * 4, -1) == 'LE')

    # ④ 不変条件
    say("\n④ 不変条件 (三値・正規化)")
    bad4 = 0
    for mant in mants:
        for t in itertools.product((-1, 0, 1), repeat=3):
            w = Word(enc_head('none', H), mant, list(t))
            if not (w.is_ternary() and w.is_normalized()):
                bad4 += 1
    say(f"   違反 {bad4}")
    ok &= (bad4 == 0)

    # ⑤ 区間が 表示値の 近傍を 正しく 含むか / 射影が 嘘を つかないか
    say("\n⑤ 区間の 健全性と 射影の 健全性")
    bad5 = lie = tested = 0
    for mant in mants:
        for t in itertools.product((-1, 0, 1), repeat=3):
            for o, k in (('none', 0), ('over', 2), ('unknown', 2)):
                w = Word(enc_head(o, H, k, 1), mant, list(t))
                lo, hi = w.interval()
                tested += 1
                if lo > hi:
                    bad5 += 1
                f = w.to_legacy_flag() & (F_GE | F_LE)
                v = w.nominal()
                if f == F_GE and min(abs(lo), abs(hi)) < abs(v) - 1e-9:
                    lie += 1
                if f == F_LE and max(abs(lo), abs(hi)) > abs(v) + 1e-9:
                    lie += 1
    say(f"   {tested} 通り: 区間の 破れ {bad5} / 射影の 嘘 {lie}")
    ok &= (bad5 == 0 and lie == 0)

    # ⑥ 上端の 往復
    say("\n⑥ 上端 (溢れ) の 往復")
    bad6 = n6 = 0
    for s in (1, -1):
        for k in range(1, H + 1):
            for st in ('over', 'unknown'):
                if st == 'unknown' and k < 2:
                    continue
                d = dec_head(enc_head(st, H, k, s))
                n6 += 1
                if not (d[0] == st and d[1] == k and (d[2] == s if st == 'over' else True)):
                    bad6 += 1
    say(f"   {n6} 通り中 不一致 {bad6}")
    ok &= (bad6 == 0)

    say("\n" + ("全て通過" if ok else "失敗あり"))
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if self_test() else 1)
