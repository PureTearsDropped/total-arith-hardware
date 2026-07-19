#!/usr/bin/env python3
"""**数は行列かもしれない**（利用者）。ならば `x/0 = 0` は選択ではなく **Moore–Penrose の定理**である。

§7.6.7⑩ は「`x/0=0` だけが群の外にあり、導出できない。手で足すしかない」と結論した。
**利用者: 「いや、あるよ。擬似逆行列と同じだよ。数は数ではなく実は行列である可能性」**

実数 x を **1×1 行列 [x]** と読む。逆元 → **擬似逆元**。すると `x/0` は `x · 0⁺` になる。**0⁺ は何か。**
"""
import numpy as np
np.set_printoptions(precision=6, suppress=True)

print("="*92); print("① **0 の擬似逆行列**"); print("="*92)
for M in ([[0.0]], [[2.0]], [[0.0, 0.0],[0.0, 0.0]], [[1.0, 0.0],[0.0, 0.0]]):
    A = np.array(M); P = np.linalg.pinv(A)
    print(f"  pinv({A.tolist()}) = **{P.tolist()}**")
print("\n  ⟹ **`pinv(0) = 0`。実装がそう言っているだけでなく、一意に決まる。次で示す。**")

print(); print("="*92); print("② **Penrose の4条件が c = 0 を強制する**"); print("="*92)
print("""  A = [0] の擬似逆元を A⁺ = [c] と置く。Penrose の4条件:

    (1) A A⁺ A = A     →  0·c·0 = 0 = A          ✓ **どんな c でも成り立つ**
    (2) A⁺ A A⁺ = A⁺   →  c·0·c = **0**、これが A⁺ = **c** に等しくねばならない
                       →  **0 = c**   ⟹ **c = 0**   ← **ここで一意に決まる**
    (3) (A A⁺)* = A A⁺ →  0 は対称                ✓
    (4) (A⁺ A)* = A⁺ A →  0 は対称                ✓

  ⟹ **条件(2) が単独で c = 0 を強制する。**
     **`0⁺ = 0` は Moore–Penrose の定理であって、便宜でも選択でもない。**""")
# 数値で総当たり確認
def penrose_ok(a, c):
    A, P = np.array([[a]]), np.array([[c]])
    return (np.allclose(A@P@A, A) and np.allclose(P@A@P, P)
            and np.allclose((A@P).T, A@P) and np.allclose((P@A).T, P@A))
cands = [-1e9, -42.0, -1.0, -1e-9, 0.0, 1e-9, 1.0, 42.0, 1e9]
ok = [c for c in cands if penrose_ok(0.0, c)]
print(f"  数値で総当たり: A=[0] に対し4条件を全部満たす c は **{ok}** ← **0 だけ**")

print(); print("="*92); print("③ **束の四つ全部が、行列代数の定理になる**"); print("="*92)
rows = [
    ("x − x = 0",   "[x] − [x] = [0]",           "零行列",              "**定理**"),
    ("x × 0 = 0",   "[x][0] = [0]",              "行列積",              "**定理**"),
    ("**x / 0 = 0**","**[x]·[0]⁺ = [x][0] = [0]**","**Moore–Penrose**",  "**定理**（②）"),
    ("x^0 = 1",     "**A⁰ = I = [1]**",          "空の積 = 単位行列",     "**定理**"),
]
print(f"  {'関係':<14}{'行列での姿':<30}{'根拠':<22}{'種類'}")
for a,b,c,d in rows: print(f"  {a:<14}{b:<30}{c:<22}{d}")
print("""
  ⟹ **§7.6.7⑩ の「x/0 だけが自由な選択」は誤りである。**
     **行列として読めば、四つ全部が定理。自由度はゼロ。**""")

print(); print("="*92); print("④ **そして我々の「白票」は、擬似逆行列そのものだった**"); print("="*92)
y = np.array([1e12, 3.0, 5.0])
Phi_bad = np.array([[0.0],[1.0],[2.0]])          # 1行目が壊れた点（0^w = 0）
c_lstsq = np.linalg.lstsq(Phi_bad, y, rcond=None)[0]
c_pinv  = np.linalg.pinv(Phi_bad) @ y
print(f"  設計 Φ = {Phi_bad.ravel().tolist()}（1行目が壊れた点）, y = {y.tolist()}")
print(f"    lstsq が返す c  = **{c_lstsq}**")
print(f"    pinv(Φ)·y      = **{c_pinv}**")
print(f"    一致するか      = **{np.allclose(c_lstsq, c_pinv)}**")
print(f"""
  ⟹ **`np.linalg.lstsq` は擬似逆行列そのものである。**
     壊れた行 `0·c = y` は **解を持たない**。擬似逆行列は **最小ノルム最小二乗解** を返す。
     **我々が「白票」と呼んでいたものは、Moore–Penrose が「解が無いときの答え」として返すものだった。**

  ⟹ **そして §7.6.7⑤⑥⑧ の三つの発見が、一つの性質になる:**
       ⑤ **てこを持たない**（白票）      ┐
       ⑥ **尺度を持たない**             ├→ **擬似逆行列は ノルムを最小化する。**
       ⑧ **符号を持たない**             ┘   **最小ノルム ＝ 最も主張しない答え。**

     **「何も知らないとき、いちばん主張しない答えを返す」は、Moore–Penrose の定義そのものである。**""")
