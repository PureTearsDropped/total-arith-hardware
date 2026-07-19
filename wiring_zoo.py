#!/usr/bin/env python3
"""配線動物園: 「配線」= 2 枚の小表（経路 mul[i][j]=k・符号 sig[i][j]=±1）だけ。

  数学的正体 = **捻れ群代数**: 群 G（位数 M）と 符号 σ で  e_i·e_j = σ(i,j)·e_{i·j}。
  経路 = 群の積、符号 = 2-コサイクル。18 ゲート圧縮器 + accumulate 網は **固定**、
  この 2 表を 差し替えるだけで 別の計算に なる（乗算器フリー・経路も データパスも 同じ）。

  ここで やること:
   1. 色んな配線を 作る（巡回=畳み込み/DFT、ブール立方=WH、CD塔、行列、…）
   2. 各配線が **同じ 底2 乗算器フリー回路で 厳密に 動くか** を検証
   3. 各配線が **何を計算しているか**（畳み込み・WH・除算代数・行列）を 同定
   4. **繋いで大きく**: 群の直積（クロネッカー）で 2 ユニットを 1 大ブロックに
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
import numpy as np
from sedenion_tensor_logic import gate9
from sd2_core import to_sd2, from_sd2, static_tree_sd2, ripple_sd2
from nd_algebra import cd_omega
from matrix_algebra import symplectic_omega, coeff_to_mat


class Wiring:
    """配線 = (経路 mul, 符号 sig)。product = 捻れ群積。"""
    def __init__(self, M, mul, sig, name, computes):
        self.M, self.mul, self.sig, self.name, self.computes = M, mul, sig, name, computes
    def product(self, x, y):
        r = [0]*self.M
        for i in range(self.M):
            for j in range(self.M):
                r[self.mul[i][j]] += self.sig[i][j]*x[i]*y[j]
        return r
    def props(self):
        rng = np.random.default_rng(0); M = self.M
        comm = assoc = True; zdiv = False
        for _ in range(40):
            a=[int(v) for v in rng.integers(-3,3,M)]; b=[int(v) for v in rng.integers(-3,3,M)]
            c=[int(v) for v in rng.integers(-3,3,M)]
            if self.product(a,b)!=self.product(b,a): comm=False
            if self.product(self.product(a,b),c)!=self.product(a,self.product(b,c)): assoc=False
        for i in range(1,M):                             # 零因子（e_i+e_j 型を軽く探索）
            for j in range(i+1,M):
                a=[0]*M; a[i]=1; a[j]=1
                for p in range(1,M):
                    for q in range(p+1,M):
                        b=[0]*M; b[p]=1; b[q]=-1
                        if all(v==0 for v in self.product(a,b)): zdiv=True
        return comm, assoc, zdiv


# ---- 配線の 作り方 ----
def from_omega(OM, M, name, computes):
    mul = [[i ^ j for j in range(M)] for i in range(M)]          # 経路 = XOR
    sig = [[int(OM[i][j]) for j in range(M)] for i in range(M)]
    return Wiring(M, mul, sig, name, computes)

def cyclic(M):
    mul = [[(i+j) % M for j in range(M)] for i in range(M)]      # 経路 = mod 足し算
    sig = [[1]*M for _ in range(M)]
    return Wiring(M, mul, sig, f"巡回 Z/{M}", "巡回畳み込み(DFT対角化)")

def boolean_cube_trivial(n):
    M = 1 << n
    mul = [[i ^ j for j in range(M)] for i in range(M)]
    sig = [[1]*M for _ in range(M)]
    return Wiring(M, mul, sig, f"ブール立方 (Z/2)^{n}", "XOR畳み込み(WH対角化)")

def kron(W1, W2, name=None):
    """群の直積 = 2 配線を クロネッカーで 繋いで 1 大ブロック（M1·M2）。"""
    M1, M2 = W1.M, W2.M; M = M1*M2
    def idx(a, b): return a*M2 + b
    mul = [[0]*M for _ in range(M)]; sig = [[0]*M for _ in range(M)]
    for a1 in range(M1):
        for a2 in range(M2):
            for b1 in range(M1):
                for b2 in range(M2):
                    i, j = idx(a1,a2), idx(b1,b2)
                    mul[i][j] = idx(W1.mul[a1][b1], W2.mul[a2][b2])
                    sig[i][j] = W1.sig[a1][b1]*W2.sig[a2][b2]
    return Wiring(M, mul, sig, name or f"({W1.name})⊗({W2.name})",
                  f"{W1.computes} ⊗ {W2.computes}")


# ---- 同じ 乗算器フリー回路で 動くか ----
def sd2_mult_wiring(xw, yw, W, ):
    K1, K2 = len(xw[0]), len(yw[0])
    cols = [[[] for _ in range(K1+K2)] for _ in range(W.M)]
    for i in range(W.M):
        for j in range(W.M):
            k, s = W.mul[i][j], W.sig[i][j]
            for p in range(K1):
                for q in range(K2):
                    cols[k][p+q].append(s*gate9(xw[i][p], yw[j][q]))
    return [ripple_sd2(static_tree_sd2(cols[k])) for k in range(W.M)]

def runs_on_hardware(W, K=12, trials=12):
    rng = np.random.default_rng(1)
    for _ in range(trials):
        x=[int(v) for v in rng.integers(-2**(K-2), 2**(K-2), W.M)]
        y=[int(v) for v in rng.integers(-2**(K-2), 2**(K-2), W.M)]
        got=[from_sd2(w) for w in sd2_mult_wiring([to_sd2(v,K) for v in x],[to_sd2(v,K) for v in y],W)]
        if got != W.product(x, y): return False
    return True


def verify_meaning():
    """各配線が 名の通り の計算を しているか（畳み込み/WH/行列）。"""
    rng = np.random.default_rng(5); out = {}
    # 巡回 → 巡回畳み込み == numpy
    W = cyclic(8); ok = True
    for _ in range(20):
        x=[int(v) for v in rng.integers(-5,5,8)]; y=[int(v) for v in rng.integers(-5,5,8)]
        conv=[int(sum(x[i]*y[(k-i)%8] for i in range(8))) for k in range(8)]
        if W.product(x,y)!=conv: ok=False
    out['巡回=巡回畳み込み'] = ok
    # ブール立方 trivial → WH で 対角化: H(x⊛y) == (Hx)∘(Hy)
    n=3; M=1<<n; W=boolean_cube_trivial(n)
    H=np.array([[1]])
    for _ in range(n): H=np.block([[H,H],[H,-H]])
    ok=True
    for _ in range(20):
        x=np.array([int(v) for v in rng.integers(-5,5,M)]); y=np.array([int(v) for v in rng.integers(-5,5,M)])
        p=np.array(W.product(list(x),list(y)))
        if not np.array_equal(H@p, (H@x)*(H@y)): ok=False
    out['ブール立方=WH対角化'] = ok
    # symplectic → 行列積（matrix_algebra で 既検証・ここは 走るかだけ）
    OM,elems,B=symplectic_omega(2); Wm=from_omega(OM,16,"","")
    ca=[int(v) for v in rng.integers(-3,3,16)]; cb=[int(v) for v in rng.integers(-3,3,16)]
    out['symplectic=行列積'] = np.array_equal(coeff_to_mat(Wm.product(ca,cb),B),
                                              coeff_to_mat(ca,B)@coeff_to_mat(cb,B))
    return out


def self_test():
    print("="*84)
    print("配線動物園 — 「配線」= 経路表 + 符号表 の 2 枚。同じ回路で 別の計算。")
    print("="*84)
    zoo = [
        cyclic(4), cyclic(8),
        boolean_cube_trivial(2), boolean_cube_trivial(3),
        from_omega(cd_omega(2), 2, "CD 複素", "除算代数(複素)"),
        from_omega(cd_omega(4), 4, "CD 四元数", "除算代数(四元)"),
        from_omega(cd_omega(16), 16, "CD セデニオン", "除算代数(零因子)"),
        from_omega(symplectic_omega(1)[0], 4, "symplectic 2×2", "行列積 2×2"),
        from_omega(symplectic_omega(2)[0], 16, "symplectic 4×4", "行列積 4×4"),
    ]
    print(f"  {'配線':<20}{'M':>4}{'可換':>5}{'結合':>5}{'零因子':>7}{'同じ回路で厳密':>14}   計算内容")
    for W in zoo:
        comm, assoc, zdiv = W.props()
        hw = runs_on_hardware(W)
        print(f"  {W.name:<20}{W.M:>4}{('✓' if comm else '×'):>5}{('✓' if assoc else '×'):>5}"
              f"{('有' if zdiv else '無'):>7}{('✓' if hw else '×'):>14}   {W.computes}")

    print()
    print("="*84)
    print("配線の 意味（名の通り 計算しているか）")
    print("="*84)
    for k, v in verify_meaning().items():
        print(f"  {k:<28}: {'✓' if v else '×'}")

    print()
    print("="*84)
    print("繋いで大きく — 群の直積（クロネッカー）で 2 ユニットを 1 大ブロックに")
    print("="*84)
    demos = [
        kron(from_omega(symplectic_omega(1)[0],4,"2×2行列","行列2×2"),
             from_omega(cd_omega(2),2,"複素","複素"), "2×2行列 ⊗ 複素 = 複素成分の2×2行列"),
        kron(cyclic(4), from_omega(cd_omega(2),2,"複素","複素"), "巡回4 ⊗ 複素 = 複素巡回畳み込み"),
        kron(boolean_cube_trivial(2), from_omega(cd_omega(4),4,"四元","四元"), "WH4 ⊗ 四元"),
    ]
    print(f"  {'合成配線':<40}{'M':>5}{'結合':>5}{'同じ回路で厳密':>14}")
    for W in demos:
        comm, assoc, zdiv = W.props()
        print(f"  {W.name:<40}{W.M:>5}{('✓' if assoc else '×'):>5}{('✓' if runs_on_hardware(W) else '×'):>14}")
    print("""
  ⟹ 「配線」は 経路表 + 符号表 の 2 枚だけ ⟹ **差し替え = 別の計算**（回路は 固定・乗算器なし）。
     群を 選ぶ: 巡回→畳み込み(DFT) / ブール立方→WH / CD塔→除算代数 / symplectic→行列積。
     群の直積で **2 ユニットを 繋いで 1 大ブロック**（複素成分の行列・複素巡回畳み込み…）。
     全部 同じ 18 ゲート圧縮器 + XOR/群経路 + accumulate で 厳密に 動く。""")


if __name__ == "__main__":
    self_test()
