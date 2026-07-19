#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証
"""4値桁研究・実験1 — SD加算は ?末尾正規形を保つか + 可能性つきキャリー(利用者の案)。

結果(全列挙照合):
  ① 加算器(Parhami式・1桁先読み)の正しさ: 3^10=59,049対 全一致
  ② **?は必ず末尾にまとまる**(途中に穴 0/400) — 末尾正規形は加算で保存される
  ③ ?深さの成長: ほぼ+1(250例)・時々+2(117例・境界の先読みが?領域を覗く分)・0(16例)
  ④ 4値化の幅の損: 中央値2.03×(区間を?語で包む固有コスト≈理論下限)・最悪5×
  ⑤ **可能性つきキャリー**(各桁=集合・キャリー線も集合を運ぶ・完全に局所=ハード化可能):
     健全性 0/3830 含み損ね・真集合と 98.1% ぴったり一致・余分 0.019個/桁
  ⟹ 不確かさは 既存のキャリー機構に乗って流れる。?の特別扱いは 要らない。
"""
import itertools, random

rnd = random.Random(7)
D = (-1, 0, 1)

def value(ds):
    return sum(d << i for i, d in enumerate(ds))

def concrete_add(a, b):
    """二進SDの限定キャリー加算(1桁先読み)。桁∈{-1,0,1}保証。"""
    n = len(a); t = [a[i] + b[i] for i in range(n)]
    c = [0] * (n + 1); s = [0] * n
    for i in range(n):
        ti = t[i]; tp = t[i - 1] if i > 0 else 0
        if ti >= 2:    c[i + 1], s[i] = 1, ti - 2
        elif ti == 1:  c[i + 1], s[i] = (1, -1) if tp >= 0 else (0, 1)
        elif ti == 0:  c[i + 1], s[i] = 0, 0
        elif ti == -1: c[i + 1], s[i] = (-1, 1) if tp <= 0 else (0, -1)
        else:          c[i + 1], s[i] = -1, ti + 2
    r = [s[i] + c[i] for i in range(n)] + [c[n]]
    assert all(-1 <= x <= 1 for x in r)
    return r

def possibility_add(A, B):
    """可能性つきキャリー: 桁i は T_i, T_{i-1}, C_{i-1} しか見ない(局所)。"""
    n = len(A)
    T = [frozenset(da + db for da in A[i] for db in B[i]) for i in range(n)]
    Cout = [None] * n; S = [None] * n
    for i in range(n):
        Tp = T[i - 1] if i > 0 else frozenset([0])
        pairs = set()
        for ti in T[i]:
            if ti >= 2:    pairs.add((1, ti - 2))
            elif ti == 1:
                if any(x >= 0 for x in Tp): pairs.add((1, -1))
                if any(x < 0 for x in Tp):  pairs.add((0, 1))
            elif ti == 0:  pairs.add((0, 0))
            elif ti == -1:
                if any(x <= 0 for x in Tp): pairs.add((-1, 1))
                if any(x > 0 for x in Tp):  pairs.add((0, -1))
            else:          pairs.add((-1, ti + 2))
        Cout[i] = frozenset(c for c, s in pairs)
        S[i] = frozenset(s for c, s in pairs)
    R = []
    for i in range(n):
        Cin = Cout[i - 1] if i > 0 else frozenset([0])
        R.append(frozenset(s + c for s in S[i] for c in Cin) & frozenset(D))
    R.append(Cout[n - 1])
    return R

def main():
    for a in itertools.product(D, repeat=5):
        for b in itertools.product(D, repeat=5):
            assert value(concrete_add(list(a), list(b))) == value(a) + value(b)
    print("① 加算器: 3^10=59,049対 全一致 ✓")

    N = 9; trials = 400
    noncontig = 0; growths = {}; loss = []
    unsound = tight = total_pos = 0; extra = 0
    for _ in range(trials):
        ma, mb = rnd.randint(0, 4), rnd.randint(0, 4)
        if ma + mb > 7: continue
        a_hi = [rnd.choice(D) for _ in range(N - ma)]
        b_hi = [rnd.choice(D) for _ in range(N - mb)]
        Rex = [set() for _ in range(N + 1)]; V = set()
        for ta in itertools.product(D, repeat=ma):
            for tb in itertools.product(D, repeat=mb):
                r = concrete_add(list(ta) + a_hi, list(tb) + b_hi)
                for i, d in enumerate(r): Rex[i].add(d)
                V.add(value(r))
        unc = [i for i in range(N + 1) if len(Rex[i]) > 1]
        if unc and unc != list(range(unc[-1] + 1)): noncontig += 1
        g = (unc[-1] + 1 - max(ma, mb)) if unc else 0
        growths[g] = growths.get(g, 0) + 1
        lo4 = sum((min(Rex[i]) if len(Rex[i]) == 1 else -1) << i for i in range(N + 1))
        hi4 = sum((max(Rex[i]) if len(Rex[i]) == 1 else 1) << i for i in range(N + 1))
        Vlo, Vhi = min(V), max(V)
        assert lo4 <= Vlo and Vhi <= hi4
        assert sorted(V) == list(range(Vlo, Vhi + 1))
        if len(V) > 1: loss.append((hi4 - lo4 + 1) / (Vhi - Vlo + 1))
        A = [frozenset(D)] * ma + [frozenset([d]) for d in a_hi]
        B = [frozenset(D)] * mb + [frozenset([d]) for d in b_hi]
        Rlo = possibility_add(A, B)
        for i in range(N + 1):
            total_pos += 1
            if not Rex[i] <= Rlo[i]: unsound += 1
            if Rex[i] == set(Rlo[i]): tight += 1
            extra += len(Rlo[i]) - len(Rex[i])
    import statistics
    print(f"② ?の途中の穴: {noncontig}/{trials} ③ 深さ成長: {dict(sorted(growths.items()))}")
    print(f"④ 4値化の損: 中央値 {statistics.median(loss):.2f}× 最悪 {max(loss):.2f}×")
    print(f"⑤ 可能性つきキャリー: 含み損ね {unsound}/{total_pos}・"
          f"一致 {100*tight/total_pos:.1f}%・余分 {extra/total_pos:.3f}個/桁")
    assert noncontig == 0 and unsound == 0

if __name__ == "__main__":
    main()
