#!/usr/bin/env python3
# ⚠️ 生成AI使用・要検証
"""配線デザイナー — Wedderburn 分解を **逆に** 使う。

  入力: 欲しい 行列ブロック（例 [3] = 3×3 行列積、[2,2] = 独立な 2×2 を 2 本 同時）
  出力: その ブロックを 分解に 含む **最小の 群**（＝配線）と 検証。

  部品:
   1. 群カタログ: 巡回 Z/n・積 Z/a×Z/b・二面体 D_n・四元数群 Q8・二重巡回 Dic_n・交代 A4・対称 S4
      （積表は 構成的に 生成し、群公理を 検証）
   2. 既約次元の 自動検出: 正則表現の 可換子環の ランダム元を 対角化 →
      固有値の 多重度の多重集合 = {d_i が d_i 回}（複素 Wedderburn）→ 次元を 読む
   3. 発見した群の 実地検証: A4 の 3次元 標準表現で ĉ = â·b̂（フーリエで 本当に 3×3 行列積か）
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import itertools
import numpy as np

rng = np.random.default_rng(20260802)


# ---------------------------------------------------------------- 1. 群カタログ
def table_from(elems, mulg):
    idx = {g: i for i, g in enumerate(elems)}
    return [[idx[mulg(a, b)] for b in elems] for a in elems], elems

def G_cyclic(n):
    return table_from(list(range(n)), lambda a, b: (a + b) % n)

def G_prod(t1, t2):
    m1, e1 = t1; m2, e2 = t2
    elems = [(a, b) for a in range(len(e1)) for b in range(len(e2))]
    return table_from(elems, lambda x, y: (m1[x[0]][y[0]], m2[x[1]][y[1]]))

def G_dihedral(n):                                       # 位数 2n: (f,k)=s^f r^k
    elems = [(f, k) for f in (0, 1) for k in range(n)]
    def mulg(g, h):
        f1, k1 = g; f2, k2 = h
        return (f1 ^ f2, (k2 + (k1 if f2 == 0 else -k1)) % n)
    return table_from(elems, mulg)

def G_dicyclic(n):                                       # 位数 4n: a^k b^j, b²=aⁿ, bab⁻¹=a⁻¹
    elems = [(k, j) for j in (0, 1) for k in range(2 * n)]
    def mulg(g, h):
        k1, j1 = g; k2, j2 = h
        k = (k1 + (k2 if j1 == 0 else -k2)) % (2 * n)
        if j1 and j2: k = (k + n) % (2 * n)              # b² = aⁿ
        return (k, j1 ^ j2)
    return table_from(elems, mulg)

def G_perm(perms):
    return table_from(perms, lambda p, q: tuple(p[q[t]] for t in range(len(q))))

def G_A4():
    perms = [p for p in itertools.permutations(range(4))
             if int(np.sign(np.linalg.det(np.eye(4)[list(p)]))) == 1]
    return G_perm(perms)

def G_S4():
    return G_perm(list(itertools.permutations(range(4))))

def G_Q8():                                              # 既検証の 表を 再構成
    from wiring_registry import _reg_bilinear             # noqa: 依存回避のため 直接構成
    from nd_algebra import cd_omega
    OM = cd_omega(4)
    elems = [(s, e) for e in range(4) for s in (1, -1)]
    def mulg(u, v):
        (sa, ea), (sb, eb) = u, v
        return (sa * sb * int(OM[ea, eb]), ea ^ eb)
    return table_from(elems, mulg)

CATALOG = [
    ("Z/4", G_cyclic(4)), ("Z/8", G_cyclic(8)), ("Z/12", G_cyclic(12)),
    ("Z/2×Z/2", G_prod(G_cyclic(2), G_cyclic(2))),
    ("Z/2×Z/4", G_prod(G_cyclic(2), G_cyclic(4))),
    ("D3=S3", G_dihedral(3)), ("D4", G_dihedral(4)), ("D5", G_dihedral(5)), ("D6", G_dihedral(6)),
    ("Q8", G_Q8()), ("Dic3", G_dicyclic(3)),
    ("A4", G_A4()), ("S4", G_S4()),
]


def group_ok(mul):
    n = len(mul)
    e = next((i for i in range(n) if all(mul[i][j] == j and mul[j][i] == j for j in range(n))), None)
    if e is None: return False
    assoc = all(mul[mul[i][j]][k] == mul[i][mul[j][k]]
                for i in range(n) for j in range(n) for k in range(n))
    latin = all(sorted(r) == list(range(n)) for r in mul) and \
            all(sorted(mul[i][j] for i in range(n)) == list(range(n)) for j in range(n))
    return assoc and latin


# ---------------------------------------------------------------- 2. 既約次元の 検出
def irrep_dims(mul):
    """複素 Wedderburn: C[G] ≅ ⊕ M_{d_i}。可換子環の ランダム元の 固有値多重度から d_i を 読む。
       正則表現 L_g と 可換するのは 右移動 ⟹ X = Σ_g c_g R_g（ランダム係数・正規行列に 対称化）。
       各 M_{d_i} ブロックは d_i 個の 固有値を 各 d_i 重で 持つ ⟹ 多重度 m の 固有値数 = m·(dim m の 既約数)。"""
    n = len(mul)
    Rg = np.zeros((n, n, n))                             # 右移動 (R_g x)[h] = x[h·g⁻¹] ⟺ R_g[h·g, h]=1
    for g in range(n):
        for h in range(n):
            Rg[g, mul[h][g], h] = 1
    c = rng.normal(size=n) + 1j * rng.normal(size=n)
    X = sum(c[g] * Rg[g] for g in range(n))
    Xh = X + X.conj().T                                  # エルミートに（固有値 実・数値安定）
    # エルミート化で 可換子環に 留まる（R_g† = R_{g⁻¹} ⟹ 和も 可換子環の 元）
    ev = np.linalg.eigvalsh(Xh)
    tol = 1e-6 * max(1.0, np.abs(ev).max())
    mults = []
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(ev[j + 1] - ev[i]) < tol: j += 1
        mults.append(j - i + 1); i = j + 1
    from collections import Counter
    cnt = Counter(mults)                                 # 多重度 m の 固有値が いくつか
    dims = []
    for m, c_m in sorted(cnt.items()):
        assert c_m % m == 0, (m, c_m)
        dims += [m] * (c_m // m)                          # dim m の 既約が c_m/m 個
    assert sum(d * d for d in dims) == n, (dims, n)       # Σd² = |G|
    return sorted(dims)


# ---------------------------------------------------------------- 3. デザイナー
def design(target_blocks):
    """欲しい ブロック（例 [3]・[2,2]）を 含む 最小位数の 群を カタログから。"""
    hits = []
    for name, (mul, elems) in CATALOG:
        dims = irrep_dims(mul)
        pool = list(dims)
        ok = True
        for t in target_blocks:                           # 多重集合として 含むか
            if t in pool: pool.remove(t)
            else: ok = False; break
        if ok:
            hits.append((len(mul), name, dims))
    hits.sort()
    return hits


# ---------------------------------------------------------------- 検証デモ
def verify_A4_matmul3():
    """デザイナーが 3×3 に 出す A4 を 実地検証: 標準3次元表現で ĉ == â·b̂。"""
    mul, perms = G_A4()
    def std_rep(p):                                       # 置換行列を Σx=0 に 制限（整数 3×3）
        P = np.zeros((4, 4), dtype=int)
        for t in range(4): P[p[t], t] = 1
        B = np.array([[1, 0, 0], [-1, 1, 0], [0, -1, 1], [0, 0, -1]])   # e0−e1, e1−e2, e2−e3
        PB = P @ B
        # (x0,x1,x2,x3), Σ=0 の 基底係数 = (x0, x0+x1, −x3)
        return np.array([PB[0], PB[0] + PB[1], -PB[3]])
    bad = 0; trials = 1000
    for _ in range(trials):
        a = [int(v) for v in rng.integers(-9, 10, 12)]
        b = [int(v) for v in rng.integers(-9, 10, 12)]
        c = [0] * 12
        for i in range(12):
            for j in range(12):
                c[mul[i][j]] += a[i] * b[j]
        F = lambda x: sum(x[i] * std_rep(perms[i]) for i in range(12))
        if not np.array_equal(F(c), F(a) @ F(b)): bad += 1
    return bad, trials


def self_test():
    print("=" * 78)
    print("配線デザイナー — 欲しい 行列ブロック → それを 含む 最小の 群（＝配線）")
    print("=" * 78)
    print("  カタログ 群公理 検証: ", end="")
    for name, (mul, _) in CATALOG:
        assert group_ok(mul), name
    print(f"{len(CATALOG)} 群 すべて ✓（単位元・結合律・ラテン方陣）")
    print()
    print(f"  {'群':<10}{'位数':>4}   複素既約次元（自動検出）")
    for name, (mul, _) in CATALOG:
        dims = irrep_dims(mul)
        print(f"  {name:<10}{len(mul):>4}   {dims}")
    print()
    print("-" * 78)
    for target, note in [([2], "2×2 行列積 1 本"), ([3], "3×3 行列積 1 本"),
                          ([2, 2], "独立な 2×2 を 2 本 同時"), ([1, 1, 1, 1], "スカラー 4 本（可換で可）")]:
        hits = design(target)
        top = " ／ ".join(f"{nm}(位数{o}, {d})" for o, nm, d in hits[:3]) if hits else "カタログに なし"
        print(f"  欲しい {str(target):<12} {note:<22} → {top}")
    print("-" * 78)
    print()
    print("  実地検証: デザイナーが [3] に 出した **A4**（位数12・最小）で 本当に 3×3 行列積か")
    bad, trials = verify_A4_matmul3()
    print(f"    A4 標準3次元表現: ĉ == â·b̂（{trials} 試行）: 違反 **{bad}**")
    print(f"    ⟹ A4 の 群配線 1 回 = スカラー3 本 ＋ **3×3 行列積 1 本** 同時（1+1+1+9=12 ✓）")
    print()
    print("  ⟹ 設計フロー: 欲しい計算 → design() が 群を 提案 → 群表 = 配線 = レジストリに 登録可。")
    print("     注: 次元は 複素 Wedderburn。実の形は 別問題（D4 は M₂(R)・Q8 は H — 同じ 次元 [1,1,1,1,2]）。")


if __name__ == "__main__":
    self_test()
