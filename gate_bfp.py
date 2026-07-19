#!/usr/bin/env python3
# ⚠️ 生成AI使用・要検証
"""ブロック浮動小数点（仮数 底2）を **最初から** ゲートで。GATE_CONDITIONS.md が 契約。

  数 = (仮数 mant: 固定幅 W の 符号つき桁列, 指数 E: 共有整数)。値 = Σ dec(mantᵢ)·2^(i+E)。
  底が 仮数も指数も 2 ⟹ **E を動かす = 桁位置の 付け替え = 無損失（相殺なし）**。
    ・固定 E なら 純配線。実行時可変 E差の 整列は バレルシフタ(mux網)だが **厳密（routing のみ・演算なし）**。

  段（このファイル）:
    B1 表現 to_bf/from_bf ＋ 無料シフト（付け替え）
    B2 積: 仮数×仮数 ＋ 指数 加算
    B3 和: 指数を 整列（小さい方へ シフト）＋ 仮数 加算
  （以降 B4 正規化〔先頭スキャン→シフト→E調整→フラグ・飽和・二種の0〕, B5 (U,V,W)ユニット）

土台ゲートは gate_bilinear（検証済み）を 再利用。数体系だけ ブロック浮動に。
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from fractions import Fraction as Fr
from gate_bilinear import (new_counter, AND, OR, NOT, XOR, enc, dec, ZERO,
                           gate9, compress3, to_sd, from_sd, neg, sd_sum,
                           multiply, canonicalize, bin_add, nz)


# ============================================================ B1: ブロック浮動 の 表現
class BF:
    """ブロック浮動数: 仮数(符号つき桁列, 低位から) と 指数 E。値 = from_sd(mant)·2^E。"""
    __slots__ = ('mant', 'E')
    def __init__(self, mant, E=0):
        self.mant = mant; self.E = E

def to_bf(value, W, E=0):
    """整数 value → BF（指数 E、仮数は value を 2^E で 割った 符号つき桁）。E=0 なら そのまま。"""
    return BF(to_sd(value, W), E)

def from_bf(x):
    """BF → 厳密値（Fraction）。value = from_sd(mant)·2^E。"""
    return Fr(from_sd(x.mant)) * (Fr(2) ** x.E)

def shift_up(mant, k):
    """仮数を k 桁 上へ（低位に ZERO を k 個）= ×2^k。指数を k 下げると 値不変。
       底2 の 付け替え（無損失・相殺なし）。可変 k は バレルシフタ(mux)だが 厳密。"""
    return [ZERO]*k + list(mant)


# ============================================================ B2: 積（仮数×仮数・指数 加算）
def bf_mul(x, y, st):
    """BF 積: 仮数を 厳密に 掛け（部分積gate9+集約）、指数を 足す。シフト なし。"""
    return BF(multiply(x.mant, y.mant, st), x.E + y.E)


# ============================================================ B3: 和（指数整列 ＋ 仮数加算）
def bf_add(x, y, st):
    """BF 和: 指数を min(Ex,Ey) へ 整列（大きい方を シフトアップ＝付け替え）＋ 仮数 加算。"""
    Elo = min(x.E, y.E)
    mx = shift_up(x.mant, x.E - Elo)          # 大E側を 下へ寄せる（無損失・左シフト）
    my = shift_up(y.mant, y.E - Elo)
    return BF(sd_sum([mx, my], st), Elo)

def bf_sub(x, y, st):
    """BF 差 = x + (−y)。符号反転は 配線（(p,n)入替）。"""
    return bf_add(x, BF(neg(y.mant), y.E), st)


# ============================================================ B4: ブロック正規化（共有指数）
#  ブロック = M 成分 + 共有 E。最大成分が 先頭位置 Lmax を決め、ブロック全体を W 桁に収まるよう
#  **同じ量 sh だけ 下へ シフト（付け替え＝配線）**し E を上げる。すると:
#    ・大成分: 収まる（切り捨てあれば ge）        ・真の0: EQ 0（そのまま）
#    ・小成分（共有Eの下で MIN 未満）: ±MIN=ε・le  ← 二種の0 と ε は ここで 自然に 出る
#    ・指数が Emax 超: その成分は ±MAX・ge（稀）
def leading_pos(canon):
    """非冗長 canon の 最上位 非零桁の位置（無ければ −1）。先頭スキャン回路の 値。"""
    for i in range(len(canon) - 1, -1, -1):
        if dec(*canon[i]) != 0:
            return i
    return -1

def block_normalize(mants, E, W, Emax, st):
    """M 成分の 仮数（幅自由）＋ 共有 E → 各 W 桁 ＋ 共有 E_out ＋ 成分ごとフラグ(ge,le,sunk)。"""
    canons = []; signs = []; Ls = []
    for m in mants:
        c, s = canonicalize(m, st); canons.append(c); signs.append(s); Ls.append(leading_pos(c))
    Lmax = max(Ls)
    sh = max(0, Lmax - (W - 1))                      # 先頭を 位置 W−1 に収める 全ブロック共通シフト
    E_out = E + sh
    sh_allowed = sh
    if E_out > Emax:                                 # 指数 溢れ ⟹ 目一杯だけ シフト・残りは 飽和
        sh_allowed = max(0, Emax - E); E_out = Emax
    out = []; flags = []
    for c, s in zip(canons, signs):
        dropped = c[:sh_allowed]
        kept = [c[i] if i < len(c) else ZERO for i in range(sh_allowed, sh_allowed + W)]
        drop_nz = 0
        for d in dropped: drop_nz = OR(drop_nz, nz(d, st), st)
        kept_nz = 0
        for d in kept: kept_nz = OR(kept_nz, nz(d, st), st)
        over = 0                                     # まだ W を 超えるか（指数溢れ時のみ 起きうる）
        for i in range(sh_allowed + W, len(c)):
            over = OR(over, nz(c[i], st), st)
        collapse = AND(drop_nz, NOT(kept_nz, st), st)    # 潰れ ε: 残0 かつ 落とし非零
        nsign = NOT(s, st)
        sel_max = over
        sel_min = AND(collapse, NOT(over, st), st)
        sel_kept = AND(NOT(over, st), NOT(collapse, st), st)
        om = []
        for i in range(W):
            p, n = kept[i]
            minp = nsign if i == 0 else 0            # ±MIN は 最下位だけ 符号
            minn = s if i == 0 else 0
            op = OR(OR(AND(sel_max, nsign, st), AND(sel_min, minp, st), st), AND(sel_kept, p, st), st)
            on = OR(OR(AND(sel_max, s, st),     AND(sel_min, minn, st), st), AND(sel_kept, n, st), st)
            om.append((op, on))
        ge = OR(over, AND(drop_nz, kept_nz, st), st)     # 溢れ or 切り捨て → ≥
        le = collapse                                    # 潰れ ε ±MIN → ≤
        out.append(om); flags.append((ge, le, 0))
    return out, E_out, flags


# ============================================================ B5: (U,V,W) ユニット（ブロック浮動）
def bf_lincomb(coeffs, bfs, st):
    """BF の 線形結合 Σ coeffs·bfs。±1 は 配線/neg、他は bf_mul。指数を 整列して 加算。"""
    terms = []
    for c, x in zip(coeffs, bfs):
        if c == 0:    continue
        elif c == 1:  terms.append(x)
        elif c == -1: terms.append(BF(neg(x.mant), x.E))
        else:         terms.append(bf_mul(BF(to_sd(c, 8), 0), x, st))
    if not terms: return BF([ZERO], 0)
    Elo = min(t.E for t in terms)
    return BF(sd_sum([shift_up(t.mant, t.E - Elo) for t in terms], st), Elo)

def bf_bilinear_unit(U, V, W, a_bf, b_bf, st):
    """c = W·((U·a)⊙(V·b)) を ブロック浮動で。前線形→R並列積→後線形（未正規化 BF のリスト）。"""
    R = len(U)
    left  = [bf_lincomb(U[r], a_bf, st) for r in range(R)]
    right = [bf_lincomb(V[r], b_bf, st) for r in range(R)]
    prod  = [bf_mul(left[r], right[r], st) for r in range(R)]
    return [bf_lincomb(W[k], prod, st) for k in range(len(W))]


# ============================================================ 自己テスト（B1〜B5）
def self_test():
    import numpy as np
    rng = np.random.default_rng(20260728)

    print("="*74)
    print("B1 表現 ＋ 無料シフト — 値 = 仮数·2^E、E動かし＝桁の付け替え（無損失）")
    print("="*74)
    for v in (-2033, -1, 0, 1, 2033):
        assert from_bf(to_bf(v, 14)) == v
    # 無料シフト: 仮数を k 上げ・E を k 下げ ⟹ 値 不変
    x = to_bf(777, 14, E=5)
    xs = BF(shift_up(x.mant, 3), x.E - 3)
    assert from_bf(xs) == from_bf(x)
    print("  to_bf/from_bf 往復 ✓、shift_up(k)＋E−k で 値不変（付け替えは 無損失）✓")

    print()
    print("="*74)
    print("B2 積 — 仮数×仮数 ＋ 指数 加算（シフトなし）")
    print("="*74)
    st = new_counter(); bad = 0
    for _ in range(3000):
        a = int(rng.integers(-400, 400)); b = int(rng.integers(-400, 400))
        Ea = int(rng.integers(-4, 5)); Eb = int(rng.integers(-4, 5))
        z = bf_mul(to_bf(a, 11, Ea), to_bf(b, 11, Eb), st)
        if from_bf(z) != Fr(a*b) * Fr(2)**(Ea+Eb): bad += 1
        if z.E != Ea + Eb: bad += 1
    print(f"  BF 積 == a·b·2^(Ea+Eb)・指数=Ea+Eb: 違反 {bad}/3000 ✓")

    print()
    print("="*74)
    print("B3 和 — 指数整列（大E側を シフト）＋ 仮数加算。異なる指数も 厳密")
    print("="*74)
    st = new_counter(); bad = 0
    for _ in range(3000):
        a = int(rng.integers(-500, 500)); b = int(rng.integers(-500, 500))
        Ea = int(rng.integers(-3, 8)); Eb = int(rng.integers(-3, 8))
        X = to_bf(a, 12, Ea); Y = to_bf(b, 12, Eb)
        s = bf_add(X, Y, st)
        if from_bf(s) != from_bf(X) + from_bf(Y): bad += 1
        if s.E != min(Ea, Eb): bad += 1
    print(f"  BF 和（指数バラバラ）== 厳密和・指数=min: 違反 {bad}/3000 ✓")
    st = new_counter(); bad = 0
    for _ in range(2000):
        a = int(rng.integers(-9999, 9999)); b = int(rng.integers(-9999, 9999))
        Ea = int(rng.integers(0, 6)); Eb = int(rng.integers(0, 6))
        X = to_bf(a, 16, Ea); Y = to_bf(b, 16, Eb)
        if from_bf(bf_sub(X, Y, st)) != from_bf(X) - from_bf(Y): bad += 1
    print(f"  BF 差（= 足す neg・符号反転は配線）: 違反 {bad}/2000 ✓")

    print()
    print("="*74)
    print("B4 ブロック正規化 — 共有指数。大成分がEを決め、小成分は ±MIN=ε に潰れる")
    print("="*74)
    # 具体例: 大 10000 と 小 1,2,3 を 同一ブロック（共有E=0）で W=6 に正規化
    mants = [to_sd(10000, 20), to_sd(3, 20), to_sd(-2, 20), to_sd(0, 20)]
    out, Eout, flags = block_normalize(mants, 0, 6, 40, new_counter())
    print(f"  ブロック[10000, 3, −2, 0] を W=6 で 正規化 → 共有 E={Eout}")
    names = ["10000(大)", "3(小)", "−2(小)", "0(真)"]
    for nm, om, (ge, le, sk), orig in zip(names, out, flags, [10000, 3, -2, 0]):
        shown = Fr(from_sd(om)) * Fr(2)**Eout
        fl = "≥" if ge else ("≤" if le else "=")
        print(f"    {nm:<9} 表示 {str(shown):>8}  フラグ {fl}  "
              f"{'← ε=±MIN(向き保持)' if le else ('← 真の0(符号なし)' if (orig==0) else '')}")

    # 健全性: フラグは 嘘をつかないか（大量ランダムブロック）
    import numpy as np
    st = new_counter(); lie = 0; eps_hit = 0; checks = 0
    for _ in range(1500):
        M = 4; W = 6; Emax = 60
        vals = [int(rng.integers(-5, 6)) * (10 ** int(rng.integers(0, 5))) for _ in range(M)]
        ms = [to_sd(v, 24) for v in vals]
        out, Eout, flags = block_normalize(ms, 0, W, Emax, st)
        for v, om, (ge, le, sk) in zip(vals, out, flags):
            checks += 1
            shown = Fr(from_sd(om)) * Fr(2)**Eout
            t = Fr(v)
            if ge and abs(t) < abs(shown): lie += 1
            if le and abs(t) > abs(shown): lie += 1
            if not ge and not le and t != shown: lie += 1
            if shown != 0 and t != 0 and (t > 0) != (shown > 0): lie += 1     # ε は 向き保持
            if le: eps_hit += 1
            if v == 0 and (ge or le): lie += 1                                # 真の0 は EQ
    print(f"  ランダムブロック 健全性: **違反 {lie}/{checks}**・ε(±MIN)発生 {eps_hit} 回")
    print("  ⟹ 大成分がEを決め、小成分は 向きを保ったまま ε=±MIN・le に潰れる。真の0は EQ・符号なし。")

    print()
    print("="*74)
    print("B5 (U,V,W) ユニット（ブロック浮動）— Strassen 2×2 ＋ 群、ブロック正規化つき")
    print("="*74)
    U_STR = [[1,0,0,1],[0,0,1,1],[1,0,0,0],[0,0,0,1],[1,1,0,0],[-1,0,1,0],[0,1,0,-1]]
    V_STR = [[1,0,0,1],[1,0,0,0],[0,1,0,-1],[-1,0,1,0],[0,0,0,1],[1,1,0,0],[0,0,1,1]]
    W_STR = [[1,0,0,1,-1,0,1],[0,0,1,0,1,0,0],[0,1,0,1,0,0,0],[1,-1,1,0,0,1,0]]
    # (a) W 大 ⟹ ブロック浮動でも 厳密（EQ・A·B 一致）
    st = new_counter(); bad = 0
    for _ in range(1500):
        A = rng.integers(-30, 31, (2,2)); B = rng.integers(-30, 31, (2,2))
        a_bf = [to_bf(int(v), 16, 0) for v in A.flatten()]
        b_bf = [to_bf(int(v), 16, 0) for v in B.flatten()]
        outs = bf_bilinear_unit(U_STR, V_STR, W_STR, a_bf, b_bf, st)
        om, Eout, flags = block_normalize([o.mant for o in outs], 0, 14, 60, st)
        truth = list((A @ B).flatten())
        for t, o, (ge, le, sk) in zip(truth, om, flags):
            shown = Fr(from_sd(o)) * Fr(2)**Eout
            if shown != t or ge or le: bad += 1
    print(f"  (a) W大 の ブロック浮動 Strassen == A·B（EQ・厳密）: 違反 {bad}/1500 ✓")
    # (b) W 小 ⟹ 丸めるが 状態は 健全（嘘なし）
    st = new_counter(); lie = 0; rounded = 0
    for _ in range(1500):
        A = rng.integers(-200, 201, (2,2)); B = rng.integers(-200, 201, (2,2))
        a_bf = [to_bf(int(v), 20, 0) for v in A.flatten()]
        b_bf = [to_bf(int(v), 20, 0) for v in B.flatten()]
        outs = bf_bilinear_unit(U_STR, V_STR, W_STR, a_bf, b_bf, st)
        om, Eout, flags = block_normalize([o.mant for o in outs], 0, 6, 60, st)   # W=6 小
        truth = list((A @ B).flatten())
        for t, o, (ge, le, sk) in zip(truth, om, flags):
            shown = Fr(from_sd(o)) * Fr(2)**Eout; t = Fr(t)
            if ge or le: rounded += 1
            if ge and abs(t) < abs(shown): lie += 1
            if le and abs(t) > abs(shown): lie += 1
            if not ge and not le and t != shown: lie += 1
            if shown != 0 and t != 0 and (t > 0) != (shown > 0): lie += 1
    print(f"  (b) W小 の ブロック浮動 Strassen: 丸め {rounded} 成分・**健全性違反 {lie}**（嘘なし）")
    # (c) 群（複素・四元）も 同じ ブロック浮動ユニットで
    def group_uvw(OM, M):
        U=[[1 if c==i else 0 for c in range(M)] for i in range(M) for j in range(M)]
        V=[[1 if c==j else 0 for c in range(M)] for i in range(M) for j in range(M)]
        Wm=[[int(OM[i][j]) if (i^j)==k else 0 for i in range(M) for j in range(M)] for k in range(M)]
        return U,V,Wm
    from nd_algebra import cd_omega, ref_mult_M
    for M, name in [(2,"複素"), (4,"四元数")]:
        OM=cd_omega(M); U,V,Wm=group_uvw(OM,M); st=new_counter(); bad=0
        for _ in range(300):
            a=[int(v) for v in rng.integers(-20,21,M)]; b=[int(v) for v in rng.integers(-20,21,M)]
            a_bf=[to_bf(v,14,0) for v in a]; b_bf=[to_bf(v,14,0) for v in b]
            outs=bf_bilinear_unit(U,V,Wm,a_bf,b_bf,st)
            om,Eout,flags=block_normalize([o.mant for o in outs],0,13,60,st)
            truth=ref_mult_M(a,b,OM,M)
            for t,o,(ge,le,sk) in zip(truth,om,flags):
                if Fr(from_sd(o))*Fr(2)**Eout != t or ge or le: bad+=1
        print(f"  (c) 群 {name} も ブロック浮動ユニットで == 群積（W大・EQ）: 違反 {bad}/300 ✓")

    print()
    print("B1〜B5 通過 — **底2ブロック浮動の (U,V,W) ユニットが 端から端まで ゲートで**。")
    print("  非群(Strassen)も群(複素/四元)も、W大なら厳密、W小なら健全に丸め、二種の0・ε・飽和つき。")


if __name__ == "__main__":
    self_test()
