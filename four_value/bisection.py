#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証
"""4値桁研究・実験0 — 利用者の仮説「先頭+1のとき、下の桁の−1はその桁の半分を表す」。

結果(全列挙・違反0):
  ① 先頭+1のSD語のMSB読み = 区間の半分選び(−1=下/0=中/+1=上)。冗長性=中央半分の重なり。
     特例: 真下の−1 → 値 = 先頭のちょうど半分。
  ② **?を末尾にm桁並べた語 ≡ 連続区間 [中心±(2^m−1)]**(隙間なし)。
     孤立した?は隙間を作る(非区間集合)。
  ⟹ 意味論の土台: (桁前置き, ?深さm) = 値+誤差棒 = 有効数字算術。?は末尾正規形を保つべき。
  次: Avizienis 加算は ?末尾正規形を保つか・深さは1桁しか伸びないか(実験1)。
"""
import itertools

def val(ds):
    n = len(ds)
    return sum(d * (1 << (n - 1 - i)) for i, d in enumerate(ds))

def main():
    K = 7
    viol = total = 0
    for tail in itertools.product([-1, 0, 1], repeat=K):
        ds = (1,) + tail
        v = val(ds)
        P = 0
        for i, d in enumerate(ds):
            j = K - i
            P += d * (1 << j)
            radius = (1 << j) - 1 if j > 0 else 0
            total += 1
            if not (P - radius <= v <= P + radius): viol += 1
    print(f"① 半分選び: 3^{K}={3**K:,}語×各桁 逸脱 {viol}/{total:,}")
    assert viol == 0
    for k in (3, 5, 7):
        assert (1 << k) - (1 << (k - 1)) == 1 << (k - 1)
    print("   特例: 真下の−1 = 先頭の半分 ✓")
    for prefix, m in [((1,), 3), ((1, -1), 3), ((1, 0, -1), 4), ((1, 1), 5)]:
        P = val(prefix) * (1 << m)
        vals = sorted({P + val(t) for t in itertools.product([-1, 0, 1], repeat=m)})
        assert vals == list(range(vals[0], vals[-1] + 1))
        assert vals[-1] - vals[0] + 1 == 2 * (1 << m) - 1
    print("② ?末尾 = 連続区間(幅 2·2^m−1) ✓")
    vals = sorted({val((1, d, 1)) for d in (-1, 0, 1)})
    assert vals != list(range(vals[0], vals[-1] + 1))
    print(f"   対照: 孤立? は 隙間 {vals} ✓")
    print("実験0: 仮説成立。?の意味論 = 末尾正規形で 区間/有効数字。")

if __name__ == "__main__":
    main()
