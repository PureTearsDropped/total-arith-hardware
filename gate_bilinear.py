#!/usr/bin/env python3
# ⚠️ 生成AI使用・要検証
"""(U,V,W) 双線形ユニットを **論理ゲートから 一から** 組む。GATE_CONDITIONS.md が 契約。

  段1 原始ゲート（AND/OR/NOT/XOR）＋ 符号つき桁の enc/dec        ← このファイル(ここまで)
  段2 gate9（符号つき桁 1×1 の 積）
  段3 全加算器 → 3:2 圧縮器（18 ゲート・値保存・冗長）
  （以降 段4〜9 は 続けて 足す）

方針: 小さいゲートから。全ゲートは 呼び出し回数を数える（面積の 目安）。各段に 自己テスト。
"""

# ============================================================ 段1: 原始ゲート
# 全ての回路は この 4 種（さらに NAND 1 種に 帰着可能）だけで 書く。st で 呼び出しを 数える。
def new_counter():
    return {'AND': 0, 'OR': 0, 'NOT': 0, 'XOR': 0}

def AND(a, b, st): st['AND'] += 1; return a & b
def OR (a, b, st): st['OR']  += 1; return a | b
def NOT(a,    st): st['NOT'] += 1; return a ^ 1
def XOR(a, b, st): st['XOR'] += 1; return a ^ b

# --- 符号つき桁 t ∈ {−1,0,+1} の エンコード（A1）---
#   線は ビット。1 桁 = 2 本 (p, n)。 値 = p − n。 (1,1) は 冗長ゼロ。
def enc(t):
    """t ∈ {−1,0,+1} → (p, n)。配線のみ（ゲート0）。"""
    return (1 if t == 1 else 0, 1 if t == -1 else 0)

def dec(p, n):
    """(p, n) → 値 p − n ∈ {−1,0,+1}（冗長ゼロ (1,1)→0）。読み出しのみ。"""
    return p - n


# ============================================================ 段2: gate9（桁 1×1 の 積）
def gate9(x, y, st):
    """符号つき桁の 積 t = x·y。入力・出力とも (p,n)。6 ゲート（AND4 + OR2）。
       正 = 同符号（xp∧yp ∨ xn∧yn）／ 負 = 異符号（xp∧yn ∨ xn∧yp）。"""
    xp, xn = x; yp, yn = y
    rp = OR(AND(xp, yp, st), AND(xn, yn, st), st)      # 同符号 → +
    rn = OR(AND(xp, yn, st), AND(xn, yp, st), st)      # 異符号 → −
    return (rp, rn)


# ============================================================ 段3: 全加算器 → 3:2 圧縮器
def full_adder(a, b, c, st):
    """ビットの 全加算器（a⊕b を 共有）: (sum=パリティ, carry=多数決)。5 ゲート。"""
    ab = XOR(a, b, st)
    s  = XOR(ab, c, st)                                 # sum = a⊕b⊕c
    cy = OR(AND(a, b, st), AND(c, ab, st), st)          # carry = ab ∨ c(a⊕b)
    return s, cy

def compress3(x, y, c, st):
    """符号つき桁 3 個 → low + 2·high（各 ∈{−1,0,+1}）。carry-free（Avizienis）。18 ゲート。
       正レール (xp,yp,cp) と 負レール (xn,yn,cn) に 全加算器を 1 個ずつ:
         low  = Ps0 − Ns0（和どうしの差）,  high = Pc − Nc（桁上げどうしの差）。"""
    xp, xn = x; yp, yn = y; cp, cn = c
    Ps0, Pc = full_adder(xp, yp, cp, st)               # 正レール
    Ns0, Nc = full_adder(xn, yn, cn, st)               # 負レール
    lp = AND(Ps0, NOT(Ns0, st), st); ln = AND(Ns0, NOT(Ps0, st), st)   # low = Ps0−Ns0
    hp = AND(Pc,  NOT(Nc,  st), st); hn = AND(Nc,  NOT(Pc,  st), st)   # high= Pc −Nc
    return (lp, ln), (hp, hn)


# ============================================================ 段4: 多桁の 符号つき数 と 加算
ZERO = (0, 0)

def to_sd(v, K):
    """整数 v → K 桁の 符号つき桁 [(p,n),...]（低位から）。d_i = sign(v)·bit_i(|v|)。"""
    s = 1 if v >= 0 else -1; a = abs(v)
    return [enc(s) if (a >> i) & 1 else ZERO for i in range(K)]

def from_sd(digits):
    """[(p,n),...] → 値 Σ dec(p,n)·2^i。読み出しのみ。"""
    return sum(dec(p, n) << i for i, (p, n) in enumerate(digits))

def neg(X):
    """符号反転 = 各桁の (p,n) を 入れ替える 配線のみ（ゲート0）。dec(n,p) = −dec(p,n)。"""
    return [(n, p) for (p, n) in X]

def sd_sum(nums, st):
    """符号つき数の リスト の 和 = 1 つの 符号つき数（圧縮木 + リップル）。
       U·a・W·(·)・積の桁上げ集約 は 全部 これ（係数±1は 選択/neg＝配線 ⟹ 乗算器フリー）。"""
    if not nums: return [ZERO]
    width = max(len(x) for x in nums)
    cols = [[] for _ in range(width + 1)]
    for x in nums:
        for i, d in enumerate(x):
            cols[i].append(d)
    k = 0                                              # 圧縮木: 各列 3→2、high は 次列へ
    while k < len(cols):
        while len(cols[k]) > 2:
            x = cols[k].pop(); y = cols[k].pop(); z = cols[k].pop()
            low, high = compress3(x, y, z, st)
            cols[k].append(low)
            if k + 1 >= len(cols): cols.append([])
            cols[k + 1].append(high)
        k += 1
    out = []; carry = ZERO                             # リップル: 残り≤2 + 桁上げ を 低位から
    for c in cols:
        a = c[0] if len(c) > 0 else ZERO
        b = c[1] if len(c) > 1 else ZERO
        low, carry = compress3(a, b, carry, st)
        out.append(low)
    out.append(carry)
    return out


# ============================================================ 段5: 乗算セル（多桁 × 多桁）
def multiply(X, Y, st):
    """符号つき数 X×Y = 部分積(gate9) を 位置 i+j に置いて sd_sum で 集約。
       ＝ 段2(gate9) と 段4(sd_sum) の 組み合わせ だけ。C1: x×0=0 は 自動（部分積 全0）。"""
    parts = []
    for i, xd in enumerate(X):
        for j, yd in enumerate(Y):
            pp = gate9(xd, yd, st)                     # 桁の積 ∈ {−1,0,+1}
            parts.append([ZERO] * (i + j) + [pp])      # 2^{i+j} の 重み ＝ 低位に 0 詰め
    return sd_sum(parts, st)


# ============================================================ 段6: (U,V,W) ユニット 組み立て
def lincomb(coeffs, nums, st, K=None):
    """線形結合 Σ coeffs[i]·nums[i]。±1 は 選択/neg＝配線（ゲート0）、0 は 捨てる、
       その他の 整数係数だけ multiply。U·a / V·b / W·(·) の 実体。"""
    terms = []
    for c, x in zip(coeffs, nums):
        if c == 0:      continue
        elif c == 1:    terms.append(x)                    # 配線
        elif c == -1:   terms.append(neg(x))               # 配線（(p,n)入替）
        else:           terms.append(multiply(to_sd(c, (K or 8)), x, st))
    return sd_sum(terms, st) if terms else [ZERO]

def bilinear_unit_gates(U, V, W, a_ints, b_ints, K, st):
    """c = W·((U·a)⊙(V·b)) を ゲートで。前線形→R並列積→後線形。戻りは 各出力の sd 数。"""
    A = [to_sd(v, K) for v in a_ints]
    B = [to_sd(v, K) for v in b_ints]
    R = len(U)
    left  = [lincomb(U[r], A, st, K) for r in range(R)]    # 前線形 U·a（R 個）
    right = [lincomb(V[r], B, st, K) for r in range(R)]    # 前線形 V·b（R 個）
    prod  = [multiply(left[r], right[r], st) for r in range(R)]   # ⊙ R 個の 独立積（並列）
    out   = [lincomb(W[k], prod, st, K) for k in range(len(W))]   # 後線形 W·(·)
    return out


# ============================================================ 段7a: 正規化（冗長→非冗長）
#  符号つき桁は 冗長（同じ値に 複数表現）⟹ 上位に 非零桁が あっても 値が 大きいとは 限らない
#  （相殺で 小さくなりうる）。飽和・「どちらの0か」を 健全に 見るには **まず 非冗長形に 解決**。
#  値 = P − N（P=p ビット列, N=n ビット列）を 2の補数で 引き、符号-大きさへ。＝実ハードの
#  「carry-save で 貯めて 最後に carry-propagate で 解決」そのもの。
def bin_add(A, B, cin, st):
    """2進 リップル加算（低位から・同長）: A+B+cin → (和ビット列, 桁上げ)。"""
    out = []; c = cin
    for a, b in zip(A, B):
        s, c = full_adder(a, b, c, st); out.append(s)
    return out, c

def canonicalize(digits, st):
    """冗長 sd → 非冗長 符号-大きさ sd（全 非零桁が 同符号）＋ 符号ビット。値は 不変。"""
    w = len(digits) + 1                                       # 符号ビットの 余裕
    P = [(digits[i][0] if i < len(digits) else 0) for i in range(w)]
    N = [(digits[i][1] if i < len(digits) else 0) for i in range(w)]
    notN = [NOT(b, st) for b in N]
    D, _ = bin_add(P, notN, 1, st)                           # P + ~N + 1 = P − N（2の補数）
    sign = D[-1]                                             # MSB = 符号（1=負）
    notD = [NOT(b, st) for b in D]
    magneg, _ = bin_add(notD, [0]*w, 1, st)                  # −D = ~D + 1（負のとき 大きさ）
    nsign = NOT(sign, st)
    out = []
    for i in range(w):
        mag = OR(AND(sign, magneg[i], st), AND(nsign, D[i], st), st)   # |値| のビット
        out.append((AND(mag, nsign, st), AND(mag, sign, st)))          # 正:(1,0) 負:(0,1) 0:(0,0)
    return out, sign


# ============================================================ 段7b: 飽和・二種の0・フラグ
def nz(d, st):
    """桁が 非零か（非冗長では p|n）。"""
    p, n = d; return OR(p, n, st)

def saturate(canon, sign, W, st):
    """非冗長 canon → W 桁 + (ge,le,sunk)。溢れ→±MAX/ge（C4）、収まる→そのまま、真0→EQ（D1）。
       **決して 巻き上がった 嘘値を 出さない**（±MAX に 張り付く）＝ Inf を 生まない。"""
    ov = 0                                                    # 位置≥W に 非零？（非冗長だから 健全）
    for i in range(W, len(canon)):
        ov = OR(ov, nz(canon[i], st), st)
    nsign = NOT(sign, st); nov = NOT(ov, st)
    out = []
    for i in range(W):
        p, n = canon[i] if i < len(canon) else ZERO
        op = OR(AND(ov, nsign, st), AND(nov, p, st), st)     # 溢れ: ±MAX 桁=符号 / 否: そのまま
        on = OR(AND(ov, sign,  st), AND(nov, n, st), st)
        out.append((op, on))
    return out, (ov, 0, 0)                                   # ge=ov, le=0, sunk=0

def quantize(canon, sign, sh, W, st):
    """右シフト sh（低位 sh 桁を落とす）→ W 桁 + フラグ。
       落として値が残る→切り捨て ge。**残りが全0で 落とした桁が非零→潰れ=ε=±MIN・le**（C5・D1）。"""
    drop_nz = 0
    for i in range(sh):
        drop_nz = OR(drop_nz, nz(canon[i] if i < len(canon) else ZERO, st), st)
    kept = [canon[i] if i < len(canon) else ZERO for i in range(sh, sh + W)]
    kept_nz = 0
    for d in kept: kept_nz = OR(kept_nz, nz(d, st), st)
    collapse = AND(drop_nz, NOT(kept_nz, st), st)            # ε: 残0 かつ 落とした桁 非零
    nsign = NOT(sign, st); ncol = NOT(collapse, st)
    out = []
    for i in range(W):
        p, n = kept[i]
        if i == 0:                                           # 潰れなら 最下位を ±MIN(=符号)
            out.append((OR(AND(collapse, nsign, st), AND(ncol, p, st), st),
                        OR(AND(collapse, sign,  st), AND(ncol, n, st), st)))
        else:                                                # 潰れなら 上位は 0
            out.append((AND(ncol, p, st), AND(ncol, n, st)))
    ge = AND(drop_nz, kept_nz, st)                           # 切り捨て（値残る）→ ≥
    le = collapse                                            # 潰れ ε ±MIN → ≤
    return out, (ge, le, 0)


# ============================================================ 自己テスト（段1〜7b）
def self_test():
    print("="*74)
    print("段1 原始ゲート ＋ enc/dec — 符号つき桁の 表現（A1）")
    print("="*74)
    for t in (-1, 0, 1):
        p, n = enc(t); assert dec(p, n) == t
    assert dec(1, 1) == 0                                # 冗長ゼロ
    st = new_counter()
    assert AND(1,0,st)==0 and OR(1,0,st)==1 and NOT(0,st)==1 and XOR(1,1,st)==0
    print("  enc/dec: {−1,0,+1} 往復 ✓、冗長ゼロ (1,1)→0 ✓、AND/OR/NOT/XOR ✓")

    print()
    print("="*74)
    print("段2 gate9 — 桁 1×1 の 積 t=x·y（全 9 通り・冗長ゼロ入力も）")
    print("="*74)
    st = new_counter()
    for xv in (-1, 0, 1):
        for yv in (-1, 0, 1):
            rp, rn = gate9(enc(xv), enc(yv), st)
            assert dec(rp, rn) == xv * yv, (xv, yv, rp, rn)
    # 冗長ゼロ (1,1) を 入れても 値は 0 になる（値の厳密性）
    for yv in (-1, 0, 1):
        rp, rn = gate9((1,1), enc(yv), new_counter())
        assert dec(rp, rn) == 0
    g = new_counter(); gate9(enc(1), enc(1), g)
    print(f"  9/9 通り 厳密 ✓、冗長ゼロ入力→0 ✓、1 積あたり {sum(g.values())} ゲート（AND4+OR2）")

    print()
    print("="*74)
    print("段3 3:2 圧縮器 — 27 通り すべて 値保存 low+2·high = x+y+c（C・D1 の芯）")
    print("="*74)
    st = new_counter()
    for xv in (-1, 0, 1):
        for yv in (-1, 0, 1):
            for cv in (-1, 0, 1):
                (lp, ln), (hp, hn) = compress3(enc(xv), enc(yv), enc(cv), new_counter())
                low, high = dec(lp, ln), dec(hp, hn)
                assert low + 2*high == xv + yv + cv, (xv, yv, cv, low, high)
                assert low in (-1,0,1) and high in (-1,0,1)
    gc = new_counter(); compress3(enc(1), enc(-1), enc(1), gc)
    print(f"  27/27 値保存・桁は {{−1,0,+1}} 内 ✓（冗長だが 値は 厳密）")
    print(f"  1 圧縮器あたり {sum(gc.values())} ゲート = 全加算器×2(10) + レール引き算(AND4+NOT4)")
    print(f"    内訳 AND {gc['AND']}  OR {gc['OR']}  NOT {gc['NOT']}  XOR {gc['XOR']}")

    print()
    print("="*74)
    print("段4 多桁の 符号つき数 と 加算 — U·a / W·(·) / 桁集約 の 芯（乗算器フリー）")
    print("="*74)
    import numpy as np
    rng = np.random.default_rng(20260727)
    K = 12
    for tag, (v,) in [("往復 to_sd/from_sd", (0,))]:
        for v in (-2033, -1, 0, 1, 2033):
            assert from_sd(to_sd(v, K)) == v
            assert from_sd(neg(to_sd(v, K))) == -v            # 符号反転=配線
    print("  to_sd/from_sd 往復 ✓、neg（配線のみ・ゲート0）✓")
    st = new_counter(); bad = 0
    for _ in range(3000):
        xs = [int(v) for v in rng.integers(-500, 500, int(rng.integers(2, 8)))]
        X = [to_sd(v, K) for v in xs]
        if from_sd(sd_sum(X, st)) != sum(xs): bad += 1
    print(f"  可変個の 和 sd_sum == 整数和: 違反 {bad}/3000 ✓（圧縮木＋リップル）")
    st = new_counter(); bad = 0
    for _ in range(2000):
        a = int(rng.integers(-9999, 9999)); b = int(rng.integers(-9999, 9999))
        got = from_sd(sd_sum([to_sd(a, K+2), neg(to_sd(b, K+2))], st))  # 減算 = 足す neg
        if got != a - b: bad += 1
    print(f"  減算 a−b（= a + neg b）: 違反 {bad}/2000 ✓")
    g = new_counter(); sd_sum([to_sd(123, K), to_sd(45, K), to_sd(-67, K)], g)
    print(f"  参考: 3 数×12桁 の 和 で 圧縮器 呼び {sum(g.values())} ゲート")

    print()
    print("="*74)
    print("段5 乗算セル — 部分積(gate9) + 集約(sd_sum)。(U·a)⊙(V·b) の R 個の 積")
    print("="*74)
    st = new_counter(); bad = 0
    for _ in range(3000):
        a = int(rng.integers(-400, 400)); b = int(rng.integers(-400, 400))
        if from_sd(multiply(to_sd(a, 11), to_sd(b, 11), st)) != a * b: bad += 1
    print(f"  多桁 × 多桁 == 整数積: 違反 {bad}/3000 ✓")
    for a in (-255, -1, 0, 7, 255):                       # C1: x×0=0
        assert from_sd(multiply(to_sd(a, 10), to_sd(0, 10), new_counter())) == 0
    print("  C1: x×0 = 0（部分積 全0 で 自動）✓")
    g = new_counter(); multiply(to_sd(200, 9), to_sd(200, 9), g)
    print(f"  参考: 9桁×9桁 の 1 積 で {sum(g.values())} ゲート（gate9 81個 + 集約木）")

    print()
    print("="*74)
    print("段6 (U,V,W) ユニット組み立て — Strassen 2×2 が ゲートレベルで A·B に 一致するか")
    print("="*74)
    # Strassen 2×2: R=7、U,V,W は 全成分 {−1,0,+1}（前後は 配線だけ・乗算器なし）
    U_STR = [[1,0,0,1],[0,0,1,1],[1,0,0,0],[0,0,0,1],[1,1,0,0],[-1,0,1,0],[0,1,0,-1]]
    V_STR = [[1,0,0,1],[1,0,0,0],[0,1,0,-1],[-1,0,1,0],[0,0,0,1],[1,1,0,0],[0,0,1,1]]
    W_STR = [[1,0,0,1,-1,0,1],[0,0,1,0,1,0,0],[0,1,0,1,0,0,0],[1,-1,1,0,0,1,0]]
    st = new_counter(); bad = 0
    for _ in range(2000):
        A = rng.integers(-30, 31, (2,2)); B = rng.integers(-30, 31, (2,2))
        out = bilinear_unit_gates(U_STR, V_STR, W_STR, list(A.flatten()), list(B.flatten()), 14, st)
        c = [from_sd(o) for o in out]
        if c != list((A @ B).flatten()): bad += 1
    print(f"  Strassen 2×2 ゲート版 == numpy A·B: 違反 **{bad}**/2000 ✓")
    print(f"    R=7 積・U,V,W∈{{−1,0,+1}} ⟹ 前後は 配線のみ、本物の乗算は 7 個だけ（乗算器フリー）")
    g = new_counter()
    A = rng.integers(-30,31,(2,2)); B = rng.integers(-30,31,(2,2))
    bilinear_unit_gates(U_STR, V_STR, W_STR, list(A.flatten()), list(B.flatten()), 14, g)
    print(f"    1 回の 2×2 積で {sum(g.values()):,} ゲート（AND{g['AND']} OR{g['OR']} NOT{g['NOT']} XOR{g['XOR']}）")

    # 群代数も 同じユニットで（U,V=選択・R=M²・W=経路+符号）
    print()
    def group_uvw(OM, M):
        U=[[1 if c==i else 0 for c in range(M)] for i in range(M) for j in range(M)]
        V=[[1 if c==j else 0 for c in range(M)] for i in range(M) for j in range(M)]
        W=[[int(OM[i][j]) if (i^j)==k else 0 for i in range(M) for j in range(M)] for k in range(M)]
        return U,V,W
    from nd_algebra import cd_omega, ref_mult_M
    for M,name in [(2,"複素"),(4,"四元数")]:
        OM=cd_omega(M); U,V,W=group_uvw(OM,M); st=new_counter(); bad=0
        for _ in range(400):
            a=[int(v) for v in rng.integers(-20,21,M)]; b=[int(v) for v in rng.integers(-20,21,M)]
            c=[from_sd(o) for o in bilinear_unit_gates(U,V,W,a,b,12,st)]
            if c!=ref_mult_M(a,b,OM,M): bad+=1
        print(f"  群 {name} も 同じ (U,V,W) ゲートユニットで == 群積: 違反 {bad}/400 ✓（R=M²={M*M}）")

    print()
    print("="*74)
    print("段7a 正規化（冗長→非冗長）— 飽和/0判定の 前提。値保存＋非冗長化＋冗長の罠")
    print("="*74)
    st = new_counter(); bad = 0
    for _ in range(3000):                                     # 冗長な数を sd_sum で 作って 解決
        xs = [int(v) for v in rng.integers(-800, 800, int(rng.integers(2, 6)))]
        red = sd_sum([to_sd(v, 12) for v in xs], st)          # 冗長形（相殺あり）
        canon, sign = canonicalize(red, st)
        if from_sd(canon) != sum(xs): bad += 1
        vals = [dec(p, n) for (p, n) in canon]                # 非冗長: 非零桁は 全て 同符号
        nz = [v for v in vals if v != 0]
        if nz and len(set(nz)) != 1: bad += 1
    print(f"  値保存＋非冗長（非零桁が同符号）: 違反 {bad}/3000 ✓")
    # 冗長の罠を 直接構成で: 位置0〜9 に −1、位置10 に +1 ＝ 2^10 − (2^10−1) = **値1**、
    #   でも 高位(位置10) に 非零桁が 立っている（冗長表現の 典型）。
    red = [enc(-1)]*10 + [enc(1)]                            # 値 = −1023 + 1024 = 1
    assert from_sd(red) == 1
    hi_red = any(dec(*red[i]) != 0 for i in range(10, len(red)))
    canon, _ = canonicalize(red, new_counter())
    hi_can = any(dec(*canon[i]) != 0 for i in range(10, len(canon)))
    print(f"  値1 の冗長形 [−1×10, +1]: 高位(位置10+)に 非零 {hi_red}／解決後は {hi_can}（値={from_sd(canon)}）")
    print(f"  ⟹ 冗長のまま 上位桁で overflow 判定すると **嘘（値1を 溢れと誤認）**。解決してから 判定（段7b）。")

    print()
    print("="*74)
    print("段7b 飽和・二種の0・フラグ — Inf/NaN を 生まず、状態は 嘘をつかない")
    print("="*74)
    Wc = 8; MAXv = (1 << Wc) - 1
    st = new_counter(); lie = 0; sat_hit = 0
    for _ in range(4000):                                    # 溢れ: ±MAX/ge・状態は 健全か
        v = int(rng.integers(-4*MAXv, 4*MAXv))
        canon, sign = canonicalize(to_sd(v, Wc + 4), st)
        out, (ge, le, sunk) = saturate(canon, sign, Wc, st)
        shown = from_sd(out)
        if ge and le: sat_hit += 0
        if abs(v) > MAXv: sat_hit += 1
        # 健全性: ge⟹|真|≥|表示|, le⟹|真|≤|表示|, EQ⟹一致, 符号（sunkでない）
        if ge and abs(v) < abs(shown): lie += 1
        if le and abs(v) > abs(shown): lie += 1
        if not ge and not le and v != shown: lie += 1
        if shown != 0 and v != 0 and (v > 0) != (shown > 0): lie += 1
    print(f"  溢れ 飽和: ±MAX に張り付き {sat_hit} 回・**健全性違反 {lie}/4000**（Inf を 生まない・嘘つかない）")
    # 二種の0: 真の0（EQ）と 潰れ0=ε=±MIN（le）を 区別
    z_canon, z_sign = canonicalize(to_sd(0, 8), new_counter())
    z_out, (zge, zle, _) = saturate(z_canon, z_sign, Wc, new_counter())
    eps_canon, eps_sign = canonicalize(to_sd(3, 12), new_counter())   # 値3 を sh=4 で 落とす（<2^4）
    e_out, (ege, ele, _) = quantize(eps_canon, eps_sign, 4, Wc, new_counter())
    print(f"  真の0:  表示 {from_sd(z_out)}  フラグ (ge={zge},le={zle}) ＝ EQ・符号なし ✓")
    print(f"  ε(値3を窓2^4未満へ): 表示 {from_sd(e_out)}=±MIN  フラグ (ge={ege},le={ele}) ＝ le・符号+ ✓")

    st = new_counter(); lie = 0; hit = 0
    for _ in range(2000):                                    # 潰れ ε: le・符号・健全
        v = int(rng.integers(-30, 31))                       # 小さい値
        canon, sign = canonicalize(to_sd(v, 16), st)
        out, (ge, le, sunk) = quantize(canon, sign, 6, Wc, st)   # 2^6 で 割る窓
        shown = from_sd(out)
        if 0 < abs(v) < (1 << 6): hit += 1                   # 潰れる はず
        if ge and abs(v) < abs(shown) * (1 << 6): lie += 1
        if le and abs(shown) != 1: lie += 1                  # 潰れなら ±MIN(=1)
        if le and v != 0 and (v > 0) != (shown > 0): lie += 1  # ε は 符号を 保つ
    print(f"  潰れ ε: ±MIN/le に {hit} 回・**健全性違反 {lie}/2000**（ε は 向きを 保つ）")

    print()
    print("="*74)
    print("段7b 組み込み — 飽和つき (U,V,W) ユニット: 溢れても Inf でなく ±MAX/ge（全滅しない）")
    print("="*74)
    U_STR = [[1,0,0,1],[0,0,1,1],[1,0,0,0],[0,0,0,1],[1,1,0,0],[-1,0,1,0],[0,1,0,-1]]
    V_STR = [[1,0,0,1],[1,0,0,0],[0,1,0,-1],[-1,0,1,0],[0,0,0,1],[1,1,0,0],[0,0,1,1]]
    W_STR = [[1,0,0,1,-1,0,1],[0,0,1,0,1,0,0],[0,1,0,1,0,0,0],[1,-1,1,0,0,1,0]]
    Wc2 = 10; MAX2 = (1 << Wc2) - 1
    st = new_counter(); lie = 0; ovf = 0
    for _ in range(1000):
        A = rng.integers(-60, 61, (2,2)); B = rng.integers(-60, 61, (2,2))
        outs = bilinear_unit_gates(U_STR, V_STR, W_STR, list(A.flatten()), list(B.flatten()), 16, st)
        truth = list((A @ B).flatten())
        for o, t in zip(outs, truth):
            canon, sign = canonicalize(o, st)
            sat, (ge, le, sunk) = saturate(canon, sign, Wc2, st)
            shown = from_sd(sat)
            if abs(t) > MAX2: ovf += 1
            if ge and abs(t) < abs(shown): lie += 1           # 嘘チェック
            if not ge and not le and t != shown: lie += 1
            if shown != 0 and t != 0 and (t > 0) != (shown > 0): lie += 1
    print(f"  飽和つき Strassen 2×2: 溢れ {ovf} 成分は ±MAX/ge、**健全性違反 {lie}**（Inf/NaN なし・嘘なし）")

    print()
    print("段1〜7b 通過 — 全域化の層（飽和±MAX/±MIN=ε・二種の0・フラグ）が ゲートで 動く。")
    print("残り: 段8 フラグ伝播（sunk）／段9 完全結線。核（表現→積和→飽和→状態）は 全部 ゲートに。")


if __name__ == "__main__":
    self_test()
