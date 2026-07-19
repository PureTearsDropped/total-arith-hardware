#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証
"""bfp_sed.add の オラクル検査 — フラグ付き入力の 許容真値を 乱択し、出力ラベルの 主張と 照合。

外部AI監査(2026-07-19)が CUDA 版の 加算に 見つけた 相殺バグ（(+m,≤)+(−m,=)→(0,≤)=嘘）を、
この 区間方式の 実装が 持たないことを 確認する 常設テスト。
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bfp_sed import BF, add, EQ, GE, LE, NB, M, BNAME

W = 12
rnd = random.Random(7)
INF_CAP = 1 << 40

def sample_true(v, flag, sunk):
    """_ival と 同じ 許容集合から 整数真値を 乱択。"""
    av = abs(v)
    if sunk:
        if flag in (EQ, LE): return rnd.randint(-av, av)
        return rnd.randint(-INF_CAP, INF_CAP)
    if flag == EQ: return v
    if flag == NB: return rnd.randint(-INF_CAP, INF_CAP)
    if flag == GE:
        if v > 0: return rnd.randint(v, max(v * 8, v + 1))
        if v < 0: return rnd.randint(min(v * 8, v - 1), v)
        return rnd.randint(-INF_CAP, INF_CAP)
    if v > 0: return rnd.randint(0, v)
    if v < 0: return rnd.randint(v, 0)
    return 0

def main():
    # ① 監査の反例（CUDA 旧版は ≤ という 嘘を 出した）
    m = 64
    c = add(BF([m] + [0]*(M-1), 0, [LE] + [EQ]*(M-1)), BF([-m] + [0]*(M-1), 0), W)
    lab = ("符号不明 " if c.sunk[0] else "") + BNAME[c.bound[0]]
    print(f"① 反例 (+{m},≤)+(−{m},=): 値={c.mant[0]}·2^{c.E} ラベル=「{lab}」")
    assert c.bound[0] == NB or c.sunk[0]
    print("   ✓ 嘘なし")

    # ② オラクル
    TRIALS = 20000
    lies = checked = 0
    for _ in range(TRIALS):
        ma = [rnd.randint(-1000, 1000) for _ in range(M)]
        mb = [rnd.randint(-1000, 1000) for _ in range(M)]
        fa = [rnd.choice([EQ, EQ, GE, LE, NB]) for _ in range(M)]
        fb = [rnd.choice([EQ, EQ, GE, LE, NB]) for _ in range(M)]
        sa = [rnd.random() < 0.1 for _ in range(M)]
        sb = [rnd.random() < 0.1 for _ in range(M)]
        ta = [sample_true(ma[k], fa[k], sa[k]) for k in range(M)]
        tb = [sample_true(mb[k], fb[k], sb[k]) for k in range(M)]
        c = add(BF(ma, 0, fa, sa), BF(mb, 0, fb, sb), W)
        for k in range(M):
            ts = ta[k] + tb[k]
            V = c.mant[k] * (1 << c.E)
            fl, sk = c.bound[k], c.sunk[k]
            checked += 1
            if fl == EQ and not sk and ts != V: lies += 1
            elif fl == EQ and sk and abs(ts) != abs(V): lies += 1
            elif fl == GE and abs(ts) < abs(V): lies += 1
            elif fl == LE and abs(ts) > abs(V): lies += 1
            if not sk and V != 0 and ts != 0 and (ts > 0) != (V > 0): lies += 1
    print(f"② オラクル: {TRIALS:,}回 × {M}成分 = {checked:,} ラベル照合 → 嘘 {lies}")
    assert lies == 0
    print("結論: 区間方式の add は 相殺でも 嘘を つかない。")

if __name__ == "__main__":
    main()
