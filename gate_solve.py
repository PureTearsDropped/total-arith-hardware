#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証
"""gate_solve — solve (方程式を解く除算 L_a⁺x) を ゲートの世界へ。M/N/O・U/V/W は 保存:

    O 層: Ben-Israel 反復 X←X(2I−LX) の「テープ」(K 回 固定・判断/分岐/除算なし)
    V 層: 行列積 = bf_bilinear_unit × _matmul_uvw (matmul も ただの (U,V,W))
    U 層: L_a の構築(配線表の 重み和 = lincomb・積なし) + X₀ = Lᵀ·2^{-s}(指数シフト = ゲート0個)
    W 層: 二層検算(前向き残差→厳密解 / 正規方程式→SING / 不成立→SING|INEXACT) + フラグ
    N 層: cd_omega(M) の 差し替え(どの代数の L か)

  二部構成 (gate_series / bfp_series と 同じ 前例):
    Part A — 純ゲート EXACT モード: AND/OR/NOT/XOR まで 落とし、**厳密有理数仕様と
      Fraction 完全一致**を assert。幅の 成長則(反復ごとに ほぼ倍)も 実測 → 正規化の 必然。
    Part B — 正規化 BFP モード(実用): 共有指数 + W 桁 整数仮数・各段 正規化で コスト線形。
      最終検算は **量子化入力に対する 厳密整数演算**(誤差ゼロの 主張)で フラグを 打つ —
      cuda_fused_solve と 同じ 契約「解けたフリを しない」の ゲートモデル版。
  スカラーの 1×1 は newton_recip.py が 既に 持つ(x←x(2−ax))。本モジュールは その 行列版。
"""
from fractions import Fraction as Fr
import numpy as np

from gate_bilinear import new_counter, to_sd, neg
from gate_bfp import BF, to_bf, from_bf, bf_mul, bf_sub, bf_lincomb, bf_bilinear_unit
from wiring_registry import _matmul_uvw
from nd_algebra import cd_omega


# ================================================================ Part A: 純ゲート EXACT
def _Lmat_int(OM, a_int):
    "L[k,j] = Σ_i OM[i,j]·a_i (経路 k=i⊻j) — 整数の 重み和(配線・積なし)"
    M = len(a_int)
    L = [[0] * M for _ in range(M)]
    for i in range(M):
        for j in range(M):
            L[i ^ j][j] += int(OM[i, j]) * a_int[i]
    return L

def gate_solve_exact(OM, a_int, x_int, F, K, st):
    """ゲートで K 回の Ben-Israel。入力: 整数量子化 a·2^{-F}, x·2^{-F}。EXACT(正規化なし)。
       返り値: (y の BF リスト, 幅の推移)。"""
    M = len(a_int)
    n = M * M
    U, V, W = _matmul_uvw(M)
    L_int = _Lmat_int(OM, a_int)
    Lf = [to_bf(L_int[r][c], max(2, abs(L_int[r][c]).bit_length() + 2), -F)
          for r in range(M) for c in range(M)]
    # X₀ = Lᵀ·2^{-s}: s は ノルム上界の 2 のベキ(指数シフト = 無損失・ゲート0個)
    n1 = max(sum(abs(L_int[r][c]) for r in range(M)) for c in range(M))
    ninf = max(sum(abs(L_int[r][c]) for c in range(M)) for r in range(M))
    s = (n1 * ninf).bit_length() - 2 * F                  # 2^s ≥ ‖L‖₁‖L‖∞
    Xf = [BF(list(Lf[c * M + r].mant), Lf[c * M + r].E - s)
          for r in range(M) for c in range(M)]            # 転置 + E シフト
    widths = [max(len(b.mant) for b in Xf)]
    two = to_bf(2, 3, 0)
    for _ in range(K):
        LX = bf_bilinear_unit(U, V, W, Lf, Xf, st)
        T = []
        for e in range(n):
            if e % (M + 1) == 0:
                T.append(bf_sub(two, LX[e], st))          # 対角: 2 − (LX)_ee
            else:
                T.append(BF(neg(LX[e].mant), LX[e].E))    # 非対角: −(LX)
        Xf = bf_bilinear_unit(U, V, W, Xf, T, st)
        widths.append(max(len(b.mant) for b in Xf))
    xb = [to_bf(v, max(2, abs(v).bit_length() + 2), -F) for v in x_int]
    y = [bf_lincomb([1] * M, [bf_mul(Xf[r * M + c], xb[c], st) for c in range(M)], st)
         for r in range(M)]
    return y, widths

def spec_solve_exact(OM, a_int, x_int, F, K):
    "同じ再帰を Fraction で(①の 検証仕様)"
    M = len(a_int)
    L = [[Fr(v, 1 << F) for v in row] for row in _Lmat_int(OM, a_int)]
    Li = _Lmat_int(OM, a_int)
    n1 = max(sum(abs(Li[r][c]) for r in range(M)) for c in range(M))
    ninf = max(sum(abs(Li[r][c]) for c in range(M)) for r in range(M))
    s = (n1 * ninf).bit_length() - 2 * F                  # ゲート側と 同じ 上界規約
    X = [[L[c][r] / Fr(2) ** s for c in range(M)] for r in range(M)]
    for _ in range(K):
        LX = [[sum(L[r][k] * X[k][c] for k in range(M)) for c in range(M)] for r in range(M)]
        T = [[(2 if r == c else 0) - LX[r][c] for c in range(M)] for r in range(M)]
        X = [[sum(X[r][k] * T[k][c] for k in range(M)) for c in range(M)] for r in range(M)]
    xq = [Fr(v, 1 << F) for v in x_int]
    return [sum(X[r][c] * xq[c] for c in range(M)) for r in range(M)]


# ================================================================ Part B: 正規化 BFP(実用)
class MBF:
    "行列ブロック浮動: 整数仮数行列 + 共有指数(BFP の 哲学: ブロックに 指数 1 個)"
    __slots__ = ("m", "E")
    def __init__(self, m, E):
        # 任意精度整数(object): 幅Wの制約は _normalize が担い、中間積は厳密(int64溢れなし)
        self.m, self.E = np.asarray(m, dtype=object), int(E)

def _normalize(m, E, W):
    "最大成分が W 桁に 収まる 共有指数へ(切り捨て = 誠実さは 最終検算が 担保)"
    mx = int(np.abs(m).max())
    sh = max(0, mx.bit_length() - W)
    return (m >> sh) if sh else m, E + sh

def _matmul(A: MBF, B: MBF, W):
    return MBF(*_normalize(A.m @ B.m, A.E + B.E, W))

def bfp_solve(OM, a, x, W=24, K=25, F=None, tol=1e-3):
    """正規化 BFP solve: 量子化 → U(L構築+指数初期化) → V(K 反復・各段 W 桁) → y →
       W(量子化入力に対する 厳密整数の 二層検算 → フラグ 0/SING/SING|INEXACT)。"""
    M = len(a)
    F = F or W
    a_int = [int(round(v * (1 << F))) for v in a]
    x_int = [int(round(v * (1 << F))) for v in x]
    # --- U: L 構築(配線)+ X₀(指数シフト) ---
    L = MBF(np.array(_Lmat_int(OM, a_int)), -F)
    n1 = int(np.abs(L.m).sum(0).max()); ninf = int(np.abs(L.m).sum(1).max())
    if n1 == 0:
        return np.zeros(M), 0x07, (0.0, 0.0)              # a=0 ⟹ L⁺=0 (a/0=0 と同型)
    s = (n1 * ninf).bit_length() + 2 * L.E
    X = MBF(*_normalize(L.m.T.copy(), L.E - s, W))
    # --- V: Ben-Israel ×K (乗算のみ・各段 正規化 = コスト線形) ---
    for _ in range(K):
        LX = _matmul(L, X, W)
        assert LX.E <= 0, "正規化後の指数が正に振れた(想定外)"
        T = MBF(np.diag(np.full(M, 2 << (-LX.E))) - LX.m, LX.E)   # 2I を 同じ指数系で 厳密に
        X = _matmul(X, T, W)
    xv = MBF(np.array(x_int).reshape(M, 1), -F)
    y = _matmul(X, xv, W + 8)
    # --- W: 厳密整数の 二層検算(量子化入力に対して 誤差ゼロの 主張) ---
    Ly_m = L.m.astype(object) @ y.m.astype(object)        # 任意精度で 厳密
    Ly_E = L.E + y.E
    assert Ly_E <= -F                                     # x·2^{-F} を Ly の指数系へ 持ち上げ
    r_m = Ly_m - np.array(x_int, dtype=object).reshape(M, 1) * (1 << (-F - Ly_E))
    r1 = float(max(abs(int(v)) for v in r_m.ravel())) * 2.0 ** Ly_E
    nr_m = L.m.astype(object).T @ r_m
    r2 = float(max(abs(int(v)) for v in nr_m.ravel())) * 2.0 ** (L.E + Ly_E)
    sx = max(abs(v) for v in x) + 1e-30
    sL = float(int(np.abs(L.m).max())) * 2.0 ** L.E + 1e-30
    if r1 / sx < tol:
        flag = 0x00
    elif r2 / (sL * sx) < tol:
        flag = 0x07                                       # SING (境界なし+SUNK の 既存語彙)
    else:
        flag = 0x07 | 0x08                                # SING|INEXACT
    yf = np.array([float(int(v)) for v in y.m.ravel()]) * 2.0 ** y.E
    return yf, flag, (r1 / sx, r2 / (sL * sx))


# ================================================================ self-test
def self_test():
    rng = np.random.default_rng(0)
    print("=" * 78)
    print("Part A: 純ゲート EXACT — ゲートグラフ ≡ 厳密有理数仕様 (Fraction 一致)")
    print("=" * 78)
    M, F, K = 4, 4, 2                                     # 四元数・K=2(EXACT は 幅が 倍々)
    OM = cd_omega(M)
    a_int = [int(round(v)) for v in rng.integers(-12, 13, M)]
    x_int = [int(round(v)) for v in rng.integers(-12, 13, M)]
    st = new_counter()
    y_gate, widths = gate_solve_exact(OM, a_int, x_int, F, K, st)
    y_spec = spec_solve_exact(OM, a_int, x_int, F, K)
    exact = all(from_bf(g) == sv for g, sv in zip(y_gate, y_spec))
    print(f"  四元数 K={K}: ゲート≡仕様 {'✓厳密一致' if exact else '✗'}  "
          f"ゲート数 {sum(st.values()):,}  仮数幅の推移 {widths} ← ほぼ倍々=正規化の必然")
    assert exact

    print()
    print("=" * 78)
    print("Part B: 正規化 BFP — 実用 solve (コスト線形・厳密整数検算つき)")
    print("=" * 78)
    M16 = 16
    OM16 = cd_omega(M16)
    # ① 正則セデニオン: pinv(float64) 二証人
    a = rng.standard_normal(M16); x = rng.standard_normal(M16)
    y, flag, (q1, q2) = bfp_solve(OM16, a, x, W=24, K=25)
    Lref = np.array(_Lmat_int(OM16, [int(round(v * (1 << 24))) for v in a])) / 2.0 ** 24
    y_ref = np.linalg.pinv(Lref) @ np.array([int(round(v * (1 << 24))) for v in x]) / 2.0 ** 24
    rel = np.abs(y - y_ref).max() / np.abs(y_ref).max()
    print(f"  ① 正則(W=24): pinv 二証人 相対 {rel:.1e}・前向き残差 {q1:.1e}・flag "
          f"{'clean ✓' if flag == 0 else hex(flag)}")
    assert rel < 1e-4 and flag == 0
    # ② 零因子 e3+e10 × 一般 x: SING を 正直に・正規方程式は 成立
    z = np.zeros(M16); z[3] = 1.0; z[10] = 1.0
    y2, flag2, (p1, p2) = bfp_solve(OM16, z, x, W=24, K=30)
    print(f"  ② 零因子×一般x: 前向き {p1:.2f}(解なしを正直に)・正規方程式 {p2:.1e}・"
          f"flag {'SING ✓' if flag2 == 0x07 else hex(flag2)}")
    assert flag2 == 0x07 and p2 < 1e-3
    # ③ W のダイヤル: 幅 → 精度
    print("  ③ Wダイヤル: ", end="")
    prev = None
    for Wd in (12, 16, 24, 32):
        yw, fw, (r1w, _) = bfp_solve(OM16, a, x, W=Wd, K=25)
        relw = np.abs(yw - y_ref).max() / np.abs(y_ref).max()
        print(f"W={Wd}:{relw:.0e} ", end="")
        prev = relw
    print(" ← 幅=精度のダイヤル ✓")
    # ④ コスト会計
    Kb = 25
    mults = Kb * 2 * M16 ** 3
    print(f"  ④ コスト(セデニオン K={Kb}): 整数乗算 {mults:,} 回/solve + 検算 ~{2 * M16 ** 2}"
          f" — W桁固定でKに線形(EXACTの倍々と対照)")
    print("done — solveの意味論(二層検算・解けたフリ禁止)が ゲートモデルまで 降りた")


if __name__ == "__main__":
    self_test()
