#!/usr/bin/env python3
# ⚠️ 生成AI使用・要検証
"""Newton 逆数 = **固定回数反復＝固定の多層配線**（算術回路）。近似商＋健全な境界。

  1/b の Newton: y ← y·(2 − b·y)、2次収束。**固定 k 段 ＝ 固定配線**（反復でなく展開）。
  固定小数（尺度 S=2^F）で: y ← ⌊ y·(2S − b·y) / S ⌋。
    ・種を 1/b の 下から取り、切り捨て（⌊⌋）も 下へ ⟹ **y は 常に S/b 以下**
      ＝ 近似は 真の 1/b の **下界** ＝ **ge（|真|≥|表示|）が 健全**。← 打ち切りが フラグに化ける
  各段の 積・差・切り捨ては gate_bfp/gate_bilinear の 検証済み ゲート演算（＝配線パターン）。
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from fractions import Fraction as Fr
from gate_bilinear import (new_counter, to_sd, from_sd, neg, sd_sum, multiply,
                           canonicalize, ZERO, OR)


def newton_recip_int(b, F, k):
    """固定小数 Newton 逆数（整数版・尺度 S=2^F・k 段）。戻り y ≈ S/b（下界）。"""
    S = 1 << F
    y = S >> b.bit_length()                       # 種 = S/2^L ≤ S/b（下界・相対誤差<2）
    if y == 0: y = 1
    for _ in range(k):                            # 固定 k 段 = 固定配線
        t = b * y                                 # 積
        u = 2 * S - t                             # 差（2 − b·y をスケール）
        y = (y * u) >> F                          # 積 ＋ 切り捨て（下へ ⟹ 下界を保つ）
    return y

def newton_recip_gates(b, F, k):
    """同じ Newton を **ゲート乗算（multiply）＋sd_sum** で。各段が 配線パターン。"""
    S = 1 << F
    def shift_right_trunc(mant, f, st):           # ⌊·/2^f⌋。**先に canonicalize**（冗長のまま
        cn, _ = canonicalize(mant, st)            #   落とすと 低位部が負のとき +1 ズレる＝段7aの罠）。
        return cn[f:] if len(cn) > f else [ZERO]   # 正の非冗長形なら 低位落とし ＝ 厳密 floor
    yv = S >> b.bit_length()
    if yv == 0: yv = 1
    y = to_sd(yv, F + 4)
    twoS = to_sd(2 * S, F + 8)
    bsd = to_sd(b, b.bit_length() + 2)
    st = new_counter()
    for _ in range(k):
        t = multiply(bsd, y, st)                          # b·y
        u = sd_sum([twoS, neg(t)], st)                    # 2S − b·y
        yu = multiply(y, u, st)                           # y·(2S−b·y)
        y = shift_right_trunc(yu, F, st)                  # ⌊·/S⌋（canonicalize→桁落とし）
    return from_sd(y), st


# ============================================================ フルゲート版（種も ゲート・監査対応）
def newton_recip_gates_full(b, F, k, Wb=22):
    """監査対応: 種生成（bit_length＋可変シフト）も ゲートに。幅 Wb は 設計定数（契約 K2）。
       種 y = 2^(F−1−L)（L=先頭位置）は **優先エンコーダの one-hot を 位置反転した 純配線**。
       **契約: 1 ≤ b < 2^F**（この域で 全中間値が 正 ⟹ trunc=floor ⟹ 整数版と 一致・下界 ge 健全。
       b ≥ 2^F の fb 経路は 配線として 純だが Newton 非収束域＝機能保証なし・再監査の指摘）。"""
    from gate_exponent import priority_encoder
    st = new_counter()
    bsd = to_sd(b, Wb)                                   # 固定幅 符号化（データ依存幅の 違反を 排除）
    _, none, onehot = priority_encoder(bsd, 6, st)       # 先頭位置の one-hot（to_sd は 非冗長）
    yp = []
    for i in range(F):                                   # y の 桁 i ← onehot[F−1−i]（配線のみ）
        pos = F - 1 - i
        yp.append(onehot[pos] if 0 <= pos < Wb else 0)
    fb = 0
    for L in range(F, Wb):                               # L ≥ F（巨大 b）→ 種 = 1
        fb = OR(fb, onehot[L], st)
    yp[0] = OR(yp[0], fb, st)
    y = [(v, 0) for v in yp] + [ZERO] * 4
    S = 1 << F
    twoS = to_sd(2 * S, F + 8)
    def srt(mant, f, st):
        cn, _ = canonicalize(mant, st)
        return cn[f:] if len(cn) > f else [ZERO]
    for _ in range(k):                                   # 固定 k 段（配線）
        t = multiply(bsd, y, st)
        u = sd_sum([twoS, neg(t)], st)
        yu = multiply(y, u, st)
        y = srt(yu, F, st)
    return from_sd(y), st


def self_test():
    import numpy as np
    rng = np.random.default_rng(20260729)
    F = 24; S = 1 << F

    print("="*74)
    print("① 固定 k 段 Newton 逆数 — 段を増やすと 誤差が 2次で 減る（固定配線）")
    print("="*74)
    for b in (3, 7, 100, 999):
        print(f"  1/{b}: ", end="")
        for k in (1, 2, 3, 5):
            y = newton_recip_int(b, F, k); err = abs(Fr(y, S) - Fr(1, b))
            print(f"k={k}:{float(err):.1e} ", end="")
        print()

    print()
    print("="*74)
    print("② 打ち切り Newton は 常に 下界 ⟹ ge（≥）が 健全（嘘をつかない）")
    print("="*74)
    bad = 0; checks = 0
    for _ in range(20000):
        b = int(rng.integers(1, 1 << 20))
        for k in (1, 2, 3, 5, 8):
            y = newton_recip_int(b, F, k); checks += 1
            if Fr(y, S) > Fr(1, b): bad += 1          # 表示 > 真 なら ge 破れ
    print(f"  y/S ≤ 1/b（表示 ≤ 真 ＝ ge健全）: 破れ **{bad}**/{checks}（0なら 打ち切りは 常に下界）")

    print()
    print("="*74)
    print("③ 完全な a/b（近似商 ＋ ge境界）— ブロック浮動の指数で レンジ縮約は タダ")
    print("="*74)
    bad = 0
    for _ in range(5000):
        a = int(rng.integers(-10**6, 10**6)); b = int(rng.integers(1, 10**6))
        y = newton_recip_int(b, F, 6)                 # 1/b の 下界
        q = Fr(a * y, S)                              # a/b ≈ a·y/S
        true = Fr(a, b)
        # a>0: 近似 ≤ 真（ge）／ a<0: 近似 ≥ 真（大きさで ge）
        if abs(q) > abs(true) + Fr(1, 100): bad += 1  # |近似| ≤ |真|（ge・小さな丸め余裕）
    print(f"  a/b 近似の 大きさ ≤ 真（ge健全・6段）: 破れ {bad}/5000 ✓")

    print()
    print("="*74)
    print("④ 同じ Newton を ゲート演算（multiply＋sd_sum）で — 整数版と 一致（配線で動く）")
    print("="*74)
    for b in (3, 7, 50, 321):
        gy, st = newton_recip_gates(b, F, 4)
        iy = newton_recip_int(b, F, 4)
        ok = (gy == iy)
        print(f"  1/{b}: ゲート版 {gy} == 整数版 {iy}  {ok}  "
              f"（4段で ゲート呼び {sum(st.values()):,}・固定配線）")

    print()
    print()
    print("=" * 74)
    print("⑤ フルゲート版（種も ゲート・監査対応）— 契約域 b < 2^F で 整数版と 一致")
    print("=" * 74)
    bad = 0; n = 0
    for b in (1, 2, 3, 7, 321, 65535, (1 << 19) + 7,
              *(int(v) for v in np.random.default_rng(9).integers(1, 1 << 20, 5))):
        for k in (1, 3):
            n += 1
            gy, _ = newton_recip_gates_full(b, F, k)
            if gy != newton_recip_int(b, F, k): bad += 1
    print(f"  gates_full == 整数版（種 = one-hot 位置反転 配線）: 違反 **{bad}**/{n}")

    print()
    print("Newton 逆数 = 固定多層配線（算術回路）が ゲートで 動く。")
    print("近似商 ＋ ge境界（下界）を 1 パスで。精度は 段数（配線の深さ）で 決まる。")


if __name__ == "__main__":
    self_test()
