#!/usr/bin/env python3
"""双線形ユニット: どんな双線形アルゴリズムも (U,V,W) 3 枚で 1 つの大ユニットに焼ける。

    c = W · ( (U·a) ⊙ (V·b) )        ⊙ = 要素ごと積（R 個の独立乗算 = 並列）
  前処理 U,V（線形）→ R 並列乗算セル → 後処理 W（線形）。

  ・非群の高速アルゴリズム（Strassen）も この形 ⟹ 1 ブロックに 焼ける
  ・**群代数は この特別な場合**（U,V = 単なる選択・R=M²・W = 群経路+符号）
  ・Strassen の U,V,W は 全成分 {−1,0,+1} ⟹ 前後処理は **符号つき加算網だけ（乗算器なし）**、
    真ん中の R 個の乗算だけが 本物の積 ＝ ファブリックの群ユニットに 撒ける。
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
import numpy as np
from nd_algebra import cd_omega


def bilinear_unit(U, V, W, a, b):
    """1 つの大ユニット: 前線形 → 並列乗算 → 後線形。R = U.shape[0] 個の掛け算。"""
    left  = U @ a                      # R 個の 左因子（a の線形結合）
    right = V @ b                      # R 個の 右因子
    prod  = left * right               # ⊙ R 個の 独立な積（並列）
    return W @ prod                    # 出力（積の線形結合）


# ---- Strassen 2×2（7 積・全成分 {−1,0,+1}）----
U_STR = np.array([[1,0,0,1],[0,0,1,1],[1,0,0,0],[0,0,0,1],[1,1,0,0],[-1,0,1,0],[0,1,0,-1]])
V_STR = np.array([[1,0,0,1],[1,0,0,0],[0,1,0,-1],[-1,0,1,0],[0,0,0,1],[1,1,0,0],[0,0,1,1]])
W_STR = np.array([[1,0,0,1,-1,0,1],[0,0,1,0,1,0,0],[0,1,0,1,0,0,0],[1,-1,1,0,0,1,0]])


def group_as_UVW(OM, M):
    """群代数の積を (U,V,W) 形に: R=M² 個の 積 p_{(i,j)}=a_i·b_j、W が 群経路+符号。
       ＝ U,V は 単なる選択（線形結合なし）・R 最大。**群は特別な場合**であることを示す。"""
    R = M*M
    U = np.zeros((R, M), dtype=int); V = np.zeros((R, M), dtype=int); W = np.zeros((M, R), dtype=int)
    for i in range(M):
        for j in range(M):
            r = i*M + j
            U[r, i] = 1                 # 左因子 = a_i を そのまま選ぶ
            V[r, j] = 1                 # 右因子 = b_j を そのまま選ぶ
            W[i ^ j, r] = int(OM[i, j]) # 積 a_i b_j を 経路 i⊕j へ・符号 σ(i,j)
    return U, V, W


def self_test():
    rng = np.random.default_rng(20260726)

    print("="*82)
    print("① 非群の高速アルゴリズム（Strassen 2×2）を 1 つの (U,V,W) ユニットに 焼く")
    print("="*82)
    ok = True; maxent = 0
    for _ in range(2000):
        A = rng.integers(-9, 10, (2,2)); B = rng.integers(-9, 10, (2,2))
        a = A.flatten(); b = B.flatten()
        c = bilinear_unit(U_STR, V_STR, W_STR, a, b)
        if not np.array_equal(c.reshape(2,2), A @ B): ok = False
        maxent = max(maxent, int(np.abs(np.concatenate([U_STR.ravel(),V_STR.ravel(),W_STR.ravel()])).max()))
    R = U_STR.shape[0]
    print(f"  2×2 行列積 == 焼いたユニット: {'✓' if ok else '×'}（2000 件）")
    print(f"  掛け算の数 R = **{R}**（直接は 8）  U,V,W の成分は 全て |·|≤{maxent} = {{−1,0,+1}}")
    print(f"  ⟹ 前処理 U,V・後処理 W は **符号つき加算だけ（乗算器なし）**、"
          f"真ん中の {R} 積だけが 本物の乗算 ＝ 独立 ＝ 並列にファブリックへ撒ける。")

    print()
    print("="*82)
    print("② 群代数は この形の 特別な場合（U,V=選択・R=M²・W=群経路+符号）")
    print("="*82)
    for M, name in [(2,"複素"), (4,"四元数")]:
        OM = cd_omega(M)
        U, V, W = group_as_UVW(OM, M)
        from nd_algebra import ref_mult_M
        ok = True
        for _ in range(500):
            a = rng.integers(-9, 10, M); b = rng.integers(-9, 10, M)
            c = bilinear_unit(U, V, W, a, b)
            if list(c) != ref_mult_M(list(a), list(b), OM, M): ok = False
        print(f"  {name:<6} M={M:>2}: (U,V,W) ユニット == 群積 {'✓' if ok else '×'}   "
              f"R = M² = {M*M}（U,V は 選択のみ・線形結合なし）")

    print()
    print("="*82)
    print("③ 統一像 — 「積の数 R」と「前後処理の重さ」の トレードオフ")
    print("="*82)
    print(f"  {'設計':<24}{'前処理 U,V':<18}{'積の数 R':<12}{'後処理 W':<16}{'再構成'}")
    print(f"  {'群ファブリック':<22}{'選択（自由配線）':<16}{'M²（最大）':<11}{'群経路+符号':<15}{'表1枚差替'}")
    print(f"  {'焼いた高速algo':<21}{'±1 加算網':<17}{'R<M²（少）':<12}{'±1 加算網':<16}{'固定機能'}")
    print(f"  {'一般双線形ユニット':<19}{'任意 U,V（プログラム可）':<12}{'任意 R':<13}{'任意 W':<16}{'U,V,W 差替'}")
    print("""
  ⟹ **答え: YES、非群の並列アルゴリズムも 1 つの大ユニットに 焼ける**（(U,V,W) 3 枚が その回路）。
     群ファブリックは その 特別な場合（U,V=選択・W=群経路）＝「入力を混ぜない代わりに R=M² 払う」。
     高速algo は 入力を混ぜて（U,V を密に）**R を減らす**。Strassen は U,V,W∈{−1,0,+1} ⟹
     前後は 加算網のみ・真ん中の R 積は 群ユニットに 撒ける ⟹ **無乗算器のまま 高速algo を 焼ける**。
  ⟹ トレードオフ: 焼くほど 効率（R小）だが 再構成が 重くなる（表1枚 → 加算網 → 固定機能）。
     全域算術（無NaN・±MAX/±MIN・状態フラグ）は この 3 層 全部に そのまま 乗る（直交）。""")


if __name__ == "__main__":
    self_test()
