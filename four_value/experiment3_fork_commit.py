#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証
"""4値桁研究・実験3 — 利用者の2案: 符号フォーク(±2桁下げ)と上位確定の下方伝播。

A) 符号フォーク: 位kの±1不明を2分岐として運ぶ。
   実測: 分岐で変わる出力桁は**最大3桁の窓**(大半は1桁)・分岐差=常に±2·2^k
   (=±2の桁下げペイロード)・遠い不明2個は4分岐が2×2局所フォークに完全分解(1500/1500)
   ⟹ **2^nの呪いは局所性で消える**(窓が重ならない限り線形コスト)。

B) 上位確定の下方伝播(MSB先行正規化): 不確かさ区間が1つの窓(重なりあり)に
   収まる限り上から桁を確定し、再センタリングして下へ。止まったら残りは?尻尾。
   実測: 健全0違反・スラック中央1.85×最悪3.87×・**接頭辞+?尻尾の正規形が自動生成**。
   隙間つき集合{−v,+v}では停止(全0確定→幻影だらけの区間に膨張)=フォークの出番。

分業(3本柱・全て利用者発案): 区間型→B / 符号の隙間→A / 置き場→?尻尾(実験0)。
多数の重なる符号不明のみノルム多項式(実験2)へ。
"""
import random
rnd = random.Random(7)
D = (-1, 0, 1)

def value(ds): return sum(d << i for i, d in enumerate(ds))

def concrete_add(a, b):
    n = len(a); t = [a[i] + b[i] for i in range(n)]
    c = [0] * (n + 1); s = [0] * n
    for i in range(n):
        ti = t[i]; tp = t[i - 1] if i > 0 else 0
        if ti >= 2:    c[i + 1], s[i] = 1, ti - 2
        elif ti == 1:  c[i + 1], s[i] = (1, -1) if tp >= 0 else (0, 1)
        elif ti == 0:  c[i + 1], s[i] = 0, 0
        elif ti == -1: c[i + 1], s[i] = (-1, 1) if tp <= 0 else (0, -1)
        else:          c[i + 1], s[i] = -1, ti + 2
    return [s[i] + c[i] for i in range(n)] + [c[n]]

def normalize(lo, hi, k):
    digits = []; j = k
    while j >= 0:
        done = False
        for d in (0, 1, -1):
            c = d << j; r = (1 << j) - 1
            if c - r <= lo and hi <= c + r:
                digits.append(d); lo -= c; hi -= c; done = True; break
        if not done: break
        j -= 1
    return digits, j + 1, (lo, hi)

def main():
    N = 12
    # A) フォークの局所性
    maxspan = 0
    for _ in range(3000):
        k = rnd.randint(1, N - 3)
        a = [rnd.choice(D) for _ in range(N)]; b = [rnd.choice(D) for _ in range(N)]
        ap = a[:]; ap[k] = 1; am = a[:]; am[k] = -1
        rp = concrete_add(ap, b); rm = concrete_add(am, b)
        diff = [i for i in range(N + 1) if rp[i] != rm[i]]
        maxspan = max(maxspan, max(diff) - min(diff) + 1)
        assert value(rp) - value(rm) == 2 << k          # ±2 桁下げペイロード
    assert maxspan <= 3
    print(f"A① フォーク窓 ≤ {maxspan} 桁・分岐差=2·2^k 常成立 ✓")
    ok = 0
    for _ in range(1500):
        k1 = rnd.randint(1, 3); k2 = rnd.randint(k1 + 6, N - 3)
        a = [rnd.choice(D) for _ in range(N)]; b = [rnd.choice(D) for _ in range(N)]
        rs = {}
        for s1 in (1, -1):
            for s2 in (1, -1):
                aa = a[:]; aa[k1] = s1; aa[k2] = s2
                rs[(s1, s2)] = concrete_add(aa, b)
        sep = True
        for i in range(N + 1):
            dg = {ss: rs[ss][i] for ss in rs}
            dep1 = dg[(1, 1)] != dg[(-1, 1)] or dg[(1, -1)] != dg[(-1, -1)]
            dep2 = dg[(1, 1)] != dg[(1, -1)] or dg[(-1, 1)] != dg[(-1, -1)]
            if dep1 and dep2: sep = False; break
        ok += sep
    assert ok == 1500
    print(f"A② 遠い不明2個: 2×2局所フォークに分解 {ok}/1500 ✓")
    # B) 上位確定
    import statistics
    slacks = []
    for _ in range(4000):
        K = 10; w = rnd.randint(1, 200)
        lo = rnd.randint(-(1 << K) + w + 1, (1 << K) - w - 1); hi = lo + w
        digs, m, (rl, rh) = normalize(lo, hi, K)
        if m > 0:
            cap = (1 << m) - 1
            assert -cap <= rl and rh <= cap
            slacks.append((2 * cap + 1) / (w + 1))
    print(f"B① 上位確定: 健全・スラック 中央 {statistics.median(slacks):.2f}× "
          f"最悪 {max(slacks):.2f}× ✓")
    digs, m, _ = normalize(-8, 8, 10)
    assert all(d == 0 for d in digs) and m == 4
    print(f"B② 隙間つき±8では停止(?深さ{m}に膨張) = フォークの出番 ✓")

if __name__ == "__main__":
    main()
