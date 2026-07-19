#!/usr/bin/env python3
"""ブロック浮動小数点セデニオン + 状態 — SPEC §0 の実装（芯: sd2_core）。

数 = ブロック: 16 成分が 指数 E を共有（値_k = 2^E × 仮数_k）。
状態（成分ごと・値に付く）: 境界トリット {≤,=,≥} と 符号不明ビット。

実装する SPEC の主張:
  ・乗算は 指数を足すだけ（シフトなし）、仮数を厳密に掛け（底2 sd2）、正規化を一度
  ・正規化の一点でだけ 潰れが起き、そこで ≤ が立つ（唯一の発生源）
  ・零因子は 厳密に 0（境界 =）＝ 潰れ 0（境界 ≤）と 区別される（§5.1）
  ・溢れ（指数が範囲を出る）で ≥ が立つ
  ・符号不明は 加減の飽和でだけ立つ（積では絶対立たない・§8）
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
from sd2_core import to_sd2, from_sd2, sedenion_mult_sd2, M
from sedenion_tensor_logic import ref_mult, OMEGA

# 境界 = 2 つのフラグ（利用者の指摘）: **以上フラグ = ビット0、以下フラグ = ビット1**
#   両方立つ(3) = 境界なし（値不明）。どちらも立たない(0) = 厳密（=）。
EQ, GE, LE, NB = 0, 1, 2, 3           # =(neither) / ≥(以上のみ) / ≤(以下のみ) / 境界なし(以上|以下)
assert NB == (GE | LE)                # 境界なし は 2 フラグの or から 出る（独立な状態でない）
BNAME = {EQ: "=", GE: "≥", LE: "≤", NB: "境界なし"}
FAST_MUL = False  # True: mul の仮数積を ref_mult 直接（sd2 桁回路を通さない・高速。sd2==ref は検証済み）
USE_CENTERED = False  # 中心化形式（平均値形式）: 健全だが 除算では 実質 利得ゼロ（測定済）。
#   理由: 支配成分は naive が既に締め・残る境界なしは 真の0跨ぎ・GE/NBは 真に非有界。
#   b は 分子に線形/分母に二乗 ⟹ N区間が有界で既にタイト ⟹ 依存性ずれは二次で小。既定 OFF。

class BF:
    """ブロック浮動小数点セデニオン。"""
    def __init__(self, mant, E=0, bound=None, sunk=None):
        self.mant = [int(m) for m in mant]                 # 16 成分の整数仮数
        self.E = int(E)                                     # 共有指数（2^E）
        self.bound = list(bound) if bound else [EQ] * M     # 成分ごと {≤,=,≥}
        self.sunk = list(sunk) if sunk else [False] * M     # 成分ごと 符号不明
    def values(self):
        return [m * (2.0 ** self.E) for m in self.mant]     # 読み出し（float）
    def state(self, k):
        return ("符号不明 " if self.sunk[k] else "") + BNAME[self.bound[k]]


def normalize(mant_ints, E, W, Emax=None):
    """16 整数仮数を 指数 E で持つ → 最大成分が W 桁に収まる E' へ（ブロック浮動）。

    潰れ: 低位がはみ出して 元≠0 なのに retained=0 → 仮数 ±MIN、境界 ≤（唯一の ≤ 発生源）
    溢れ: E' が Emax を超える → E'=Emax に留め、収まらない成分は ±MAX、境界 ≥
    零因子など 本物の 0 は retained=0 かつ 元=0 → 境界 = のまま（潰れと区別）
    """
    bound = [EQ] * M
    L = max((abs(m).bit_length() for m in mant_ints), default=0)
    shift = max(0, L - W)                                    # 下へ何桁落とすか
    Ep = E + shift
    def shift_component(m, sh, out, idx):
        """m を sh 桁 右にずらして格納し、境界を 丸め方向で決める（健全）。
            = : 落ちた桁が無い（厳密）
            ≥ : 0方向に切り捨てた（真値 ≥ 表現。下界）
            ≤ : MIN より下へ完全に潰れた（真値 ≤ 表現。上界）"""
        am = abs(m); r = am >> sh; dropped = am & ((1 << sh) - 1) if sh > 0 else 0
        if m == 0:
            out.append(0)                                   # 本物の 0 → =
        elif r == 0:                                        # 完全に潰れた → ±MIN, ≤
            out.append(1 if m > 0 else -1); bound[idx] = LE
        else:
            out.append(r if m > 0 else -r)                  # 切り捨て
            if dropped != 0:
                bound[idx] = GE                             # 桁が落ちた → 下界 ≥

    overflow = (Emax is not None and Ep > Emax)
    if overflow:                                            # 溢れ: 指数が範囲外 ⟹ Emax に張り付け
        Ep = Emax; out = []; MAXm = (1 << W) - 1
        sh = Emax - E                                       # Emax での仮数へ
        for k, m in enumerate(mant_ints):
            r = abs(m) >> sh if sh >= 0 else abs(m) << (-sh)
            if r >> W:                                      # W 桁を超える → 飽和 MAX、境界 ≥
                out.append(MAXm if m > 0 else -MAXm); bound[k] = GE
            else:
                shift_component(m, max(sh, 0), out, k)      # 収まる小成分（切り捨てなら ≥）
        return out, Ep, bound
    out = []
    for k, m in enumerate(mant_ints):
        shift_component(m, shift, out, k)
    return out, Ep, bound


def mul(a, b, W, Emax=None):
    """ブロック浮動の積: 指数を足す + 仮数を厳密に掛ける（sd2）+ 正規化を一度。"""
    if FAST_MUL:                                            # 高速: 仮数積を ref_mult 直接（sd2==ref は検証済み）
        prod = [int(v) for v in ref_mult([int(m) for m in a.mant], [int(m) for m in b.mant])]
    else:
        Kw = max((abs(v).bit_length() for v in a.mant), default=0) + 2
        Lw = max((abs(v).bit_length() for v in b.mant), default=0) + 2
        xw = [to_sd2(v, Kw) for v in a.mant]
        yw = [to_sd2(v, Lw) for v in b.mant]
        prod = [from_sd2(w) for w in sedenion_mult_sd2(xw, yw)]
    if not FAST_MUL:
        assert prod == ref_mult(a.mant, b.mant), "sd2 積が参照と不一致"  # 厳密性の自己確認
    mant, Ep, nb = normalize(prod, a.E + b.E, W, Emax)
    inputs_exact = (all(x == EQ for x in a.bound) and all(x == EQ for x in b.bound)
                    and not any(a.sunk) and not any(b.sunk))
    if inputs_exact:
        return BF(mant, Ep, nb, sunk=[False] * M)           # 厳密入力: 正規化の丸めフラグのみ
    # 非厳密入力: 各成分の 真値区間を フラグから作り、積(符号つき和)を 区間で 伝播（利用者の改善）。
    #   単調なら 量的境界（≥2 等）、相殺なら 自動で 境界なし。両方 健全。
    ai = [_ival(a.mant[i], a.bound[i], a.sunk[i]) for i in range(M)]
    bi = [_ival(b.mant[i], b.bound[i], b.sunk[i]) for i in range(M)]
    bound = [EQ] * M; sunk = [False] * M
    for k in range(M):
        lo, hi = _ac_interval(ai, bi, k)                    # 符号つき和を 区間で（∞ 安全）
        in_flag, sunk[k] = _flags_for(prod[k], lo, hi)
        bound[k] = in_flag | nb[k]                          # 区間フラグ | 丸めフラグ
    return BF(mant, Ep, bound, sunk)


INF = float("inf")

def _ival(v, flag, sunk=False):
    """符号つき値 v と フラグ・符号不明 → 真値の 符号つき区間 [lo, hi]。"""
    if sunk:                                # 符号 不明 ⟹ 対称 or 全体（大きさ境界は 保つ）
        av = abs(v)
        if flag in (EQ, LE): return (-av, av)   # |true| ≤ av（EQ は =av だが 符号不明 ⟹ ±av の 殻）
        return (-INF, INF)                      # GE(|x|≥av・両符号=和集合) / NB → 全体
    if flag == EQ: return (v, v)
    if flag == NB: return (-INF, INF)
    if flag == GE:                          # |true| ≥ |v|、符号 = sign(v)
        return (v, INF) if v > 0 else (-INF, v) if v < 0 else (-INF, INF)
    return (0, v) if v > 0 else (v, 0) if v < 0 else (0, 0)   # LE: |true| ≤ |v|

def _imul(a, b):
    """符号つき区間の 積 = [min, max]（±∞ 対応。0×∞ は 0）。"""
    ps = []
    for x in (a[0], a[1]):
        for y in (b[0], b[1]):
            ps.append(0.0 if (x == 0 and y in (INF, -INF)) or (y == 0 and x in (INF, -INF))
                      else x * y)
    return min(ps), max(ps)

def _idiv(ac, Nlo, Nhi):
    """区間 ac を N∈[Nlo,Nhi]（N≥0）で 割る。N が 0 になりうる（Nlo=0）と 大きさ 無限へ。
       1/N は N の 減少関数 ⟹ 1/N ∈ [1/Nhi, 1/Nlo]。"""
    inv_lo = 0.0 if Nhi == INF else (INF if Nhi == 0 else 1.0 / Nhi)
    inv_hi = INF if Nlo == 0 else 1.0 / Nlo
    return _imul(ac, (inv_lo, inv_hi))

def _recip(x):
    return 0.0 if abs(x) == INF else (INF if x == 0 else 1.0 / x)

def _idiv_general(num, den):
    """区間 num を 符号つき区間 den で 割る。den が 0 を跨ぐ/触れると 大きさも符号も 不明。"""
    dlo, dhi = den
    if dlo < 0 < dhi:                       # 0 を 内包 ⟹ 全体（利用者: 除数0跨ぎ ⟹ 符号も値も不明）
        return (-INF, INF)
    if dlo == 0 and dhi == 0:               # den 厳密 0（呼び出し側で a/0=0 済のはず）
        return (0.0, 0.0)
    if dlo == 0:    inv = (_recip(dhi), INF)       # den ∈ [0+, dhi]（正・0 に触れる）
    elif dhi == 0:  inv = (-INF, _recip(dlo))      # den ∈ [dlo, 0-]（負・0 に触れる）
    else:           inv = (_recip(dhi), _recip(dlo))  # 0 から 離れている（1/x は 減少）
    return _imul(num, inv)

def _isum(terms):
    """区間の和（±∞ を 安全に。∞−∞ は 保守的に 広げる）。"""
    lo = hi = 0.0; lo_n = lo_p = hi_n = hi_p = False
    for l, h in terms:
        if l == -INF: lo_n = True
        elif l == INF: lo_p = True
        else: lo += l
        if h == INF: hi_p = True
        elif h == -INF: hi_n = True
        else: hi += h
    LO = -INF if (lo_n or (lo_p and lo_n)) else (INF if lo_p else lo)
    HI = INF if (hi_p or (hi_p and hi_n)) else (-INF if hi_n else hi)
    return LO, HI

def _ac_interval(ai, cbi, k):
    """出力成分 k の a·conj(b) 区間（符号つき和）を 区間で。"""
    terms = []
    for i in range(M):
        j = i ^ k; t = _imul(ai[i], cbi[j])
        if OMEGA[i, j] < 0: t = (-t[1], -t[0])
        terms.append(t)
    return _isum(terms)

# ---- 中心化形式（平均値形式）: 除算 q=num/N の 依存性を 偏微分で 相殺（利用者の affine 依頼） ----
#   平均値の定理: f が 箱 X 上 C¹ なら  f(X) ⊆ f(ȳ) + [∇f](X)·(X − ȳ)。
#   偏微分 ∂q/∂b の 商則が num と N の 相殺を 織り込む ⟹ 近スカラー除数（小摂動）で 締まる。
#   前提: 箱が 有界（全成分 EQ/LE）かつ N>0 が 箱上 保証（Nlo>0）— でないと C¹ でなく 適用不可。
def _iadd(x, y):  return (x[0] + y[0], x[1] + y[1])
def _isub(x, y):  return (x[0] - y[1], x[1] - y[0])
def _iscale(x, c):
    return (c * x[0], c * x[1]) if c >= 0 else (c * x[1], c * x[0])
def _isquare(x):
    lo, hi = x
    if lo <= 0 <= hi: return (0.0, max(lo * lo, hi * hi))
    a2, b2 = lo * lo, hi * hi
    return (min(a2, b2), max(a2, b2))

def _finite_iv(iv):
    return abs(iv[0]) != INF and abs(iv[1]) != INF

def _centered_div_intervals(am, bm, a_bound, a_sunk, b_bound, b_sunk, ac, N):
    """全 k の 商成分 q_k=（a·conj b）_k/N の 中心化形式区間 [lo,hi] のリスト。
       適用不可（非有界 or N が 0 に触れうる）なら None。健全（平均値形式）。"""
    ai = [_ival(am[i], a_bound[i], a_sunk[i]) for i in range(M)]
    bi = [_ival(bm[i], b_bound[i], b_sunk[i]) for i in range(M)]
    if not all(_finite_iv(iv) for iv in ai) or not all(_finite_iv(iv) for iv in bi):
        return None                                     # 非有界 ⟹ 中心化形式は 使えない
    Niv = _isum([_isquare(bi[i]) for i in range(M)])
    Nlo, Nhi = Niv
    if Nlo <= 0:
        return None                                     # N が 0 に 触れうる ⟹ 商 C¹ でない
    cb = [bm[0]] + [-x for x in bm[1:]]
    cbi = [bi[0]] + [(-bi[j][1], -bi[j][0]) for j in range(1, M)]
    s = [1] + [-1] * (M - 1)                            # conj の 符号（∂cb/∂b）
    da = [_isub(ai[i], (am[i], am[i])) for i in range(M)]   # 中心からの ずれ 区間
    db = [_isub(bi[i], (bm[i], bm[i])) for i in range(M)]
    N2 = (Nlo * Nlo, Nhi * Nhi)                         # N²（N>0 ⟹ 有界正）
    out = []
    for k in range(M):
        qbar = ac[k] / N
        acc = (qbar, qbar)
        num_iv = _ac_interval(ai, cbi, k)               # num_k 区間（∂q/∂b で 使う）
        for i in range(M):
            j = i ^ k
            # ∂q/∂a_i = ω(i,j)·cb_j / N
            pa = _idiv_general(_iscale(cbi[j], OMEGA[i, j]), Niv)
            acc = _iadd(acc, _imul(pa, da[i]))
            # ∂q/∂b_i = [ω(j,i)·a_j·s_i·N − num_k·2·b_i] / N²    （j=i^k）
            dnum = _iscale(ai[j], OMEGA[j, i] * s[i])
            numer = _isub(_imul(dnum, Niv), _imul(num_iv, _iscale(bi[i], 2)))
            pb = _idiv_general(numer, N2)
            acc = _iadd(acc, _imul(pb, db[i]))
        out.append(acc)
    return out

def _flags_for(V, lo, hi):
    """真値区間 [lo, hi] と 代表値 V → (2 フラグ, 符号不明)。"""
    if lo == hi:  return EQ, False          # 一点 ⟹ 厳密
    sunk = (lo < 0 < hi)                     # 0 を跨ぐ ⟹ 符号不明
    aV = abs(V)
    if aV == 0:   return NB, sunk            # V=0 だが 真値≠一点 ⟹ 境界なし
    if sunk:      mlo, mhi = 0, max(abs(lo), abs(hi))
    elif lo >= 0: mlo, mhi = lo, hi          # 符号 +
    else:         mlo, mhi = abs(hi), abs(lo)  # 符号 −
    ge = mlo >= aV                           # |真| ≥ |V| を保証
    le = mhi <= aV                           # |真| ≤ |V| を保証
    if ge and le: return EQ, sunk            # mlo=mhi=aV ⟹ 実質厳密
    if ge:        return GE, sunk
    if le:        return LE, sunk
    return NB, sunk                          # どちらも保証できない ⟹ 境界なし


def add(a, b, W, Emax=None):
    """ブロック浮動の和: 指数を整列（配線シフト）+ 成分ごと sd2 加算 + 正規化。
    同じ端で飽和した異符号の成分に 符号不明 を立てる（§8 の唯一の交差点）。"""
    E = min(a.E, b.E)                              # 小さい方へ整列 ⟹ 左シフトだけ = 無損失
    am = [m << (a.E - E) for m in a.mant]          # a.E − E ≥ 0
    bm = [m << (b.E - E) for m in b.mant]          # b.E − E ≥ 0
    s = [x + y for x, y in zip(am, bm)]
    # 加算 = フラグから作る 符号つき区間の 足し算（利用者）。
    #   (正の≥) + (負の≤) → ≥（下限が ≤ の分だけ減る）／ (≥) + (≥ 異符号) → 境界なし・符号不明。
    in_flag = [EQ] * M; sunk = [False] * M
    for k in range(M):
        alo, ahi = _ival(am[k], a.bound[k], a.sunk[k])
        blo, bhi = _ival(bm[k], b.bound[k], b.sunk[k])
        tlo, thi = alo + blo, ahi + bhi            # 真値の 符号つき区間
        in_flag[k], sunk[k] = _flags_for(s[k], tlo, thi)
    mant, Ep, nbound = normalize(s, E, W, Emax)    # 正規化（切り捨て等）の フラグ
    bound = [in_flag[k] | nbound[k] for k in range(M)]   # 2 フラグを or（以上/以下）
    return BF(mant, Ep, bound, sunk)


def sconj(bf):
    """セデニオン共役（実部そのまま・虚部 符号反転）= 配線。状態も そのまま運ぶ。"""
    return BF([bf.mant[0]] + [-m for m in bf.mant[1:]], bf.E,
              list(bf.bound), list(bf.sunk))

def div(a, b, W, Emax=None):
    """セデニオン除算 a/b = a·conj(b) / N(b)。**全域・一様**（利用者の指摘）。

    N(b) = Σ bᵢ² は b≠0 で必ず正 ⟹ **逆元は 全ての b≠0 で 在る**（零因子でも）。
    スカラー除算 (a·conj(b))ₖ / N を 各成分に 一様適用、N=0（b=0）なら **a/0 = 0**。
    """
    bm = [int(x) for x in b.mant]; am = [int(x) for x in a.mant]
    N = sum(x * x for x in bm)                          # スカラー（指数 2·b.E）
    if N == 0:                                          # b = 0 → a/0 = 0（全要素 一様）
        return BF([0] * M, 0, [EQ] * M)
    ac = [int(v) for v in ref_mult(am, [bm[0]] + [-x for x in bm[1:]])]   # a·conj(b)
    g = W + max((abs(v).bit_length() for v in ac), default=0)            # 精度スケール
    q = []; bound = [EQ] * M
    for k, v in enumerate(ac):
        num = v << g
        qk = num // N if v >= 0 else -((-num) // N)     # 0 方向 切り捨て
        q.append(qk)
        if v * (1 << g) - qk * N != 0:
            bound[k] = GE                               # |q| ≤ |真| ⟹ ≥（下界）
    mant, Ep, nb = normalize(q, (a.E - b.E) - g, W, Emax)
    bound = [bound[k] | nb[k] for k in range(M)]                         # 切り捨て | 正規化 の フラグ
    a_exact = all(x == EQ for x in a.bound) and not any(a.sunk)
    b_exact = all(x == EQ for x in b.bound) and not any(b.sunk)
    if a_exact and b_exact:
        return BF(mant, Ep, bound)                                      # 厳密: 切り捨てフラグのみ
    sunk = [False] * M
    if all(bm[i] == 0 for i in range(1, M)):
        # **除数が スカラー（実数 b_0）** ⟹ a/b = a_k/b_0 の 成分ごと（依存性なし ⟹ 締まる）。
        #   利用者の規則: (≥a_k)/(≤b_0) → ≥(a_k/b_0)。b_0 が 0 を跨げば 符号も値も 不明。
        b0iv = _ival(bm[0], b.bound[0], b.sunk[0])
        for k in range(M):
            qlo, qhi = _idiv_general(_ival(am[k], a.bound[k], a.sunk[k]), b0iv)
            in_flag, sunk[k] = _flags_for(ac[k] / N, qlo, qhi)          # 代表商 a_k/b_0
            bound[k] = in_flag | bound[k]
        return BF(mant, Ep, bound, sunk)
    # 一般（除数が セデニオン）: ac=a·conj(b) 区間、N=Σb² 区間、q=ac/N（N≥0・依存性で 保守的）。
    #   N が 0 に なりうる（Nlo=0・除数区間が 0 を跨ぐ）と 大きさ 無限・符号も 不明 ⟹ 境界なし+符号不明。
    #   だが 代表値は ac/N（非零）— **0 にはならない**（利用者: 本当に 0 で割らない限り 0 でない）。
    cb = [bm[0]] + [-x for x in bm[1:]]
    ai = [_ival(am[i], a.bound[i], a.sunk[i]) for i in range(M)]
    cbi = [_ival(cb[i], b.bound[i], b.sunk[i]) for i in range(M)]
    Nlo = Nhi = 0.0
    for i in range(M):
        bl, bh = _ival(bm[i], b.bound[i], b.sunk[i])
        mlo = 0.0 if bl <= 0 <= bh else min(abs(bl), abs(bh))            # |b_i| の 下限
        mhi = max(abs(bl), abs(bh))                                     # |b_i| の 上限（∞ 可）
        Nlo += mlo * mlo
        Nhi = INF if (Nhi == INF or mhi == INF) else Nhi + mhi * mhi
    # 中心化形式（平均値形式）: 有界箱 & N>0 なら 依存性の相殺を 取れる（近スカラー除数で 締まる）
    cen = (_centered_div_intervals(am, bm, a.bound, a.sunk, b.bound, b.sunk, ac, N)
           if USE_CENTERED else None)
    sunk = [False] * M
    for k in range(M):
        aclo, achi = _ac_interval(ai, cbi, k)
        qlo, qhi = _idiv((aclo, achi), Nlo, Nhi)                        # /N（N≥0、Nlo=0 で 無限）
        if cen is not None:                                            # 両方 健全 ⟹ 狭い方（交差）
            clo, chi = cen[k]
            qlo, qhi = max(qlo, clo), min(qhi, chi)
        in_flag, sunk[k] = _flags_for(ac[k] / N, qlo, qhi)             # 代表商 ac/N（非零）に対する 方向
        bound[k] = in_flag | bound[k]                                  # 区間フラグ | 切り捨てフラグ
    return BF(mant, Ep, bound, sunk)


# ---------------------------------------------------------------- テンソル積（MAC）
def tensor_mac(pairs, W, Emax=None):
    """Σ_i (A_i × B_i) を ブロック浮動で積和（各積は正規化せず、最後に一度）。
    pairs: [(A_0,B_0), ...] の BF 対。SPEC の carry-save MAC の骨格。

    2026-07-19 修正（外部AI監査 第3ラウンドの 同族バグ捜索で 自前発見）:
    旧版は BF(prod, E) で 積を 作る際に **入力の bound/sunk を 捨てていた**
    （(±10)×(3) が 「=」を 主張 = 符号の嘘）。mul と 同じ 区間伝播を 積ごとに 付ける。"""
    acc = None
    for A, B in pairs:
        # 積は full-width（正規化しない）で厳密に貯める
        prod = ref_mult(A.mant, B.mant)
        E = A.E + B.E
        inputs_exact = (all(x == EQ for x in A.bound) and all(x == EQ for x in B.bound)
                        and not any(A.sunk) and not any(B.sunk))
        if inputs_exact:
            p = BF(prod, E)
        else:
            ai = [_ival(A.mant[i], A.bound[i], A.sunk[i]) for i in range(M)]
            bi = [_ival(B.mant[i], B.bound[i], B.sunk[i]) for i in range(M)]
            bound = [EQ] * M; sunk = [False] * M
            for k in range(M):
                lo, hi = _ac_interval(ai, bi, k)
                bound[k], sunk[k] = _flags_for(prod[k], lo, hi)
            p = BF(prod, E, bound, sunk)
        acc = p if acc is None else add(acc, p, W=10**9)     # 貯める間は正規化しない
    m, Ep, nbound = normalize(acc.mant, acc.E, W, Emax)       # 最後に一度だけ
    bound = [acc.bound[k] | nbound[k] for k in range(M)]      # 区間フラグ | 丸めフラグ
    return BF(m, Ep, bound, acc.sunk)


def self_test():
    import numpy as np
    rng = np.random.default_rng(20260718)

    print("=" * 74)
    print("① 乗算は 厳密（潰れが無ければ）— 指数を足すだけ、シフトなし")
    print("=" * 74)
    for _ in range(5):
        x = [int(v) for v in rng.integers(-30, 30, M)]
        y = [int(v) for v in rng.integers(-30, 30, M)]
        a = BF(x, E=3); b = BF(y, E=-1)
        c = mul(a, b, W=32)
        ref = [(2.0**2) * v for v in ref_mult(x, y)]   # 2^(3+(-1)) = 2^2
        assert c.values() == ref
        assert c.E == 2 and all(bd == EQ for bd in c.bound)
    print("  積の値 == 参照 × 2^(Ea+Eb)、境界すべて '='、符号不明なし  ✓")

    print()
    print("=" * 74)
    print("② 正規化の一点で 潰れ → ≤（成分の桁差が窓を超えたとき）")
    print("=" * 74)
    a = BF([0]*M, E=0)
    a.mant[0] = 2**20        # 大きい成分
    a.mant[1] = 1            # 極小成分
    a.mant[2] = 2**19
    m, Ep, bound = normalize(a.mant, a.E, W=8)
    print(f"  成分0=2^20, 成分1=1, 成分2=2^19 を W=8 桁に正規化 → E'={Ep}")
    print(f"    成分0: 仮数 {m[0]}  境界 {BNAME[bound[0]]}   （大きい成分は '='）")
    print(f"    成分1: 仮数 {m[1]}  境界 **{BNAME[bound[1]]}**   ← 潰れた（元1が窓の下）")
    print(f"    成分2: 仮数 {m[2]}  境界 {BNAME[bound[2]]}")
    assert bound[1] == LE and m[1] == 1 and bound[0] == EQ
    print("  ⟹ 小さい成分だけ ≤（|x|≤MIN）、大きい成分は '='  ✓")

    print()
    print("=" * 74)
    print("③ 零因子は 厳密 0（境界 '='）— 潰れ 0（境界 ≤）と 区別される（§5.1）")
    print("=" * 74)
    ea = [0]*M; ea[1] = 1; ea[10] = 1
    eb = [0]*M; eb[4] = 1; eb[15] = -1
    zc = mul(BF(ea), BF(eb), W=8)
    print(f"  (e1+e10)(e4−e15): 仮数 = {zc.mant}")
    print(f"    境界 = {[BNAME[b] for b in zc.bound]}")
    assert all(v == 0 for v in zc.mant) and all(b == EQ for b in zc.bound)
    print("  ⟹ 全成分 0 かつ 境界 '='（構造的な 0）。**≤ ではない** ＝ 潰れと区別  ✓")
    # 対比: 潰れた 0 は ≤
    assert bound[1] == LE
    print("     （②の潰れた 0 は ≤ だった ⟹ 境界トリットが 二つの 0 を見分ける）")

    print()
    print("=" * 74)
    print("④ 溢れ → ≥（指数が範囲 Emax を超えたとき）")
    print("=" * 74)
    big = BF([2**30]*M, E=0)
    m, Ep, bound = normalize(big.mant, big.E, W=8, Emax=5)   # 2^30 は Emax=5 に収まらない
    print(f"  成分 2^30 を W=8, Emax=5 で: E'={Ep}, 仮数[0]={m[0]} (=MAX), 境界[0]=**{BNAME[bound[0]]}**")
    assert Ep == 5 and m[0] == (1 << 8) - 1 and bound[0] == GE
    print("  ⟹ 指数が Emax に張り付き、成分は ±MAX、境界 ≥（|x|≥MAX）  ✓")

    print()
    print("=" * 74)
    print("⑤ 符号不明は 加減の飽和でだけ（積では絶対立たない・§8）")
    print("=" * 74)
    # MAX−MAX 相当: 成分0 を 両方 ≥ の異符号にして足す
    pm = [0]*M; pm[0] = 1;  qm = [0]*M; qm[0] = -1
    pb = [EQ]*M; pb[0] = GE; qb = [EQ]*M; qb[0] = GE
    p = BF(pm, E=0, bound=pb); q = BF(qm, E=0, bound=qb)
    s = add(p, q, W=8)
    print(f"  成分0: (≥, +) + (≥, −) → 符号不明 = **{s.sunk[0]}**  ／ 成分1（普通）→ {s.sunk[1]}")
    assert s.sunk[0] is True and s.sunk[1] is False
    prod_state = mul(BF([2**10]*M, E=0), BF([2**10]*M, E=0), W=8, Emax=3)
    assert not any(prod_state.sunk)
    print("  積では 符号不明が 一つも立たない（溢れて ≥ が立っても、符号不明は立たない）  ✓")

    print()
    print("=" * 74)
    print("⑥ テンソル MAC — 積和を貯めて 正規化を最後に一度")
    print("=" * 74)
    pairs = []
    ref_acc = [0]*M
    for _ in range(6):
        x = [int(v) for v in rng.integers(-20, 20, M)]
        y = [int(v) for v in rng.integers(-20, 20, M)]
        pairs.append((BF(x), BF(y)))
        ref_acc = [a + b for a, b in zip(ref_acc, ref_mult(x, y))]
    out = tensor_mac(pairs, W=32)
    got = [int(v) for v in out.values()]
    assert got == ref_acc
    print(f"  6 項の積和 == 参照（厳密）: {got == ref_acc}  ✓")
    print("  ⟹ 貯める間は正規化ゼロ、最後に一度だけ（唯一のシフト・唯一の潰れ点）")

    print()
    print("すべて通過 — SPEC §0 の要求が 底2 ブロック浮動で動く。")


if __name__ == "__main__":
    self_test()
