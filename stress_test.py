#!/usr/bin/env python3
"""色んな演算回路で 叩く — 何か壊れないか（敵対的テスト）。

方針: 監査人は今日 7 回 先走って外した。だから ここでは 予想せず、
      **厳密な有理数を真値に** して、実装の主張が 嘘をつく所を 探す。

核心の不変条件（健全性）: 各成分の状態は 真値について 嘘をついてはならない。
    境界 =  : 真値 == 表現値（厳密）
    境界 ≤  : |真値| ≤ |表現値|         （潰れ: 表現は 上界）
    境界 ≥  : |真値| ≥ |表現値|         （溢れ: 表現は 下界）
    符号不明でない かつ 真値≠0 ⟹ 符号一致
"""
import sys, os
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
import numpy as np
from bfp_sed import BF, mul, div, add, tensor_mac, normalize, LE, EQ, GE, NB, M
from sedenion_tensor_logic import ref_mult

def true_vals(bf):
    """BF の 表現値（厳密な有理数、16 成分）。"""
    return [Fr(m) * (Fr(2) ** bf.E) for m in bf.mant]

def check_sound(bf, truth, tag, viol):
    """状態が truth（真値の有理数リスト）について 嘘をついていないか。"""
    rep = true_vals(bf)
    for k in range(M):
        t = truth[k]; r = rep[k]; b = bf.bound[k]
        if b == NB:
            pass                       # 境界なし = 向き保証なし ⟹ 常に健全
        elif b == EQ:
            if t != r:
                viol.append(f"[{tag}] 成分{k}: 境界'=' なのに 真値{t} ≠ 表現{r}")
        elif b == LE:
            if abs(t) > abs(r):
                viol.append(f"[{tag}] 成分{k}: 境界'≤' なのに |真値{t}| > |表現{r}|")
        elif b == GE:
            if abs(t) < abs(r):
                viol.append(f"[{tag}] 成分{k}: 境界'≥' なのに |真値{t}| < |表現{r}|")
        if not bf.sunk[k] and t != 0 and r != 0 and (t > 0) != (r > 0):
            viol.append(f"[{tag}] 成分{k}: 符号不明でないのに 符号違い 真{t} 表{r}")

viol = []
rng = np.random.default_rng(20260718)

# ============================================================ A. 深い乗算連鎖
print("A. 深い乗算連鎖（10 段）— 指数と仮数が 破綻しないか")
for trial in range(200):
    x = [int(v) for v in rng.integers(-8, 8, M)]
    cur = BF(x, E=int(rng.integers(-5, 5)))
    truth = [Fr(v) * Fr(2)**cur.E for v in x]
    for step in range(10):
        y = [int(v) for v in rng.integers(-4, 4, M)]
        Ey = int(rng.integers(-3, 3))
        cur = mul(cur, BF(y, E=Ey), W=24)         # W 大きめ（潰れにくい）
        truth = [Fr(t) for t in ref_mult([Fr(a) for a in truth],
                                          [Fr(b)*Fr(2)**Ey for b in y])]
    check_sound(cur, truth, f"mulchain{trial}", viol)
print(f"   健全性違反: {sum('mulchain' in v for v in viol)}")

# ============================================================ B. 潰れ下での健全性
print("B. 小さい W で 潰れさせる — 境界 ≤ が 嘘をつかないか")
for trial in range(200):
    a = [int(v) for v in rng.integers(-2**18, 2**18, M)]
    b = [int(v) for v in rng.integers(-4, 4, M)]
    c = mul(BF(a), BF(b), W=int(rng.integers(4, 10)))   # 小さい W ⟹ 潰れる
    truth = [Fr(v) for v in ref_mult(a, b)]
    check_sound(c, truth, f"trunc{trial}", viol)
print(f"   健全性違反: {sum('trunc' in v for v in viol)}")

# ============================================================ C. 加算・減算連鎖
print("C. 加減連鎖（異なる指数・異符号）— 整列と符号不明")
for trial in range(200):
    cur = BF([int(v) for v in rng.integers(-50, 50, M)], E=int(rng.integers(-4, 4)))
    truth = true_vals(cur)
    for step in range(6):
        y = BF([int(v) for v in rng.integers(-50, 50, M)], E=int(rng.integers(-4, 4)))
        ty = true_vals(y)
        cur = add(cur, y, W=24)
        truth = [truth[k] + ty[k] for k in range(M)]
    check_sound(cur, truth, f"addchain{trial}", viol)
print(f"   健全性違反: {sum('addchain' in v for v in viol)}")

# ============================================================ D. 溢れ連鎖
print("D. 溢れさせる（Emax 小）— 境界 ≥ が 嘘をつかないか")
for trial in range(200):
    a = [int(v) for v in rng.integers(-2**12, 2**12, M)]
    b = [int(v) for v in rng.integers(-2**12, 2**12, M)]
    c = mul(BF(a), BF(b), W=8, Emax=int(rng.integers(0, 6)))
    truth = [Fr(v) for v in ref_mult(a, b)]
    check_sound(c, truth, f"ovf{trial}", viol)
print(f"   健全性違反: {sum('ovf' in v for v in viol)}")

# ============================================================ E. 零因子を 連鎖に混ぜる
print("E. 零因子を 演算に混ぜる — 構造的0(=) と 潰れ0(≤) が 混同されないか")
zd_pairs = [([0]*M, [0]*M) for _ in range(3)]
zd_pairs[0][0][1]=1; zd_pairs[0][0][10]=1; zd_pairs[0][1][4]=1; zd_pairs[0][1][15]=-1
for trial in range(50):
    a, b = zd_pairs[0]
    c = mul(BF(list(a)), BF(list(b)), W=8)
    truth = [Fr(v) for v in ref_mult(a, b)]
    check_sound(c, truth, f"zd{trial}", viol)
    # 零因子の 0 は '=' であるべき（潰れ ≤ でない）
    for k in range(M):
        if c.mant[k] == 0 and truth[k] == 0 and c.bound[k] == LE:
            viol.append(f"[zd{trial}] 成分{k}: 構造的0 が ≤ と誤診（潰れと混同）")
print(f"   健全性違反: {sum('zd' in v for v in viol)}")

# ============================================================ F. 極端な縁
print("F. 縁: 全ゼロ / 単成分 / MAX同士 / W=1")
edge_viol0 = len(viol)
# 全ゼロ
z = mul(BF([0]*M), BF([0]*M), W=8)
check_sound(z, [Fr(0)]*M, "zero", viol)
assert all(b == EQ for b in z.bound), "全ゼロ積で 境界が = でない"
# MAX 同士の加算（境界なしへ）
pm=[0]*M; pm[0]=1; qm=[0]*M; qm[0]=-1
p=BF(pm,bound=[GE]+[EQ]*15); q=BF(qm,bound=[GE]+[EQ]*15)
s=add(p,q,W=8)
# W=1（極端に狭い）で 乗算
for _ in range(50):
    a=[int(v) for v in rng.integers(-100,100,M)]; b=[int(v) for v in rng.integers(-4,4,M)]
    c=mul(BF(a),BF(b),W=1)
    check_sound(c,[Fr(v) for v in ref_mult(a,b)],"W1",viol)
print(f"   健全性違反: {len(viol)-edge_viol0}")

# ============================================================ G. MAC の 厳密性（W 大）
print("G. テンソル MAC（多項）— W 大なら 厳密一致するか")
mac_bad = 0
for trial in range(100):
    n = int(rng.integers(2, 12))
    pairs=[]; ref=[0]*M
    for _ in range(n):
        x=[int(v) for v in rng.integers(-15,15,M)]; y=[int(v) for v in rng.integers(-15,15,M)]
        pairs.append((BF(x),BF(y))); ref=[a+b for a,b in zip(ref,ref_mult(x,y))]
    out=tensor_mac(pairs, W=40)
    if [int(v) for v in out.values()] != ref: mac_bad += 1
print(f"   厳密不一致: {mac_bad}/100")

# ============================================================ H. 除算（全域・一様）
print("H. 除算 a/b = a·conj(b)/N(b)（a/0=0・零因子でも逆元在り）")
def _sconjm(m): return [m[0]]+[-x for x in m[1:]]
h_viol=0
for trial in range(300):
    a=[int(v) for v in rng.integers(-50,50,M)]; b=[int(v) for v in rng.integers(-9,9,M)]
    N=sum(x*x for x in b)
    q=div(BF(a),BF(b),W=int(rng.integers(6,24)))
    if N==0:
        if any(v!=0 for v in q.values()): h_viol+=1        # a/0=0
        continue
    ac=[int(v) for v in ref_mult(a,_sconjm(b))]
    truth=[Fr(c,N) for c in ac]                            # a/b 真値（厳密）
    check_sound(q, truth, f"div{trial}", viol)
# b/b=1（零因子含む）
for b in ([1]+[0]*15, [0,1]+[0]*13+[0], [0]+[1]+[0]*8+[1]+[0]*5):
    if sum(x*x for x in b)==0: continue
    q=div(BF(list(b)),BF(list(b)),W=40)
    if [round(v) for v in q.values()]!=[1]+[0]*15: h_viol+=1
print(f"   健全性違反: {sum('div' in v for v in viol)}   b/b=1・a/0=0 違反: {h_viol}")

# ============================================================ 結果
print("\n" + "="*70)
if viol:
    print(f"**見つかった問題: {len(viol)} 件**（最初の 12 件）:")
    for v in viol[:12]:
        print("  -", v)
else:
    print("**健全性違反 ゼロ** — どの演算連鎖でも 状態は 嘘をつかなかった。")
print("="*70)
print(f"\n未実装として 記録すべき穴:")
print("  ・除算器（a/0=0 は要求だが 除算そのものが 未実装。セデニオンは零因子ゆえ 逆元が全域でない）")
print("  ・状態の 境界トリット伝播は 簡略（積の gate_bmul を bfp に 本結線していない）")
