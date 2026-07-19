#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証
"""bfp_sed の オラクル検査 — add / mul / tensor_mac。

フラグ付き入力の 許容真値（_ival と 同一の 集合・表示0の 自由項 20%・SUNK 15% 込み）を
乱択し、出力ラベルの 主張と 照合する。真値側は 純Python 整数（オーバーフローなし・厳密）。

経緯: 外部AI監査(2026-07-19) が CUDA 版に 見つけた バグ族（フラグ落とし・表示0と真の0の
混同・SUNK相殺）を この 区間方式 実装でも 捜索。add/mul は 区間伝播で 元から 健全、
**tensor_mac は 積を BF(prod,E) で 作る際に 入力フラグを 捨てていた**（(±10)×(3) が
「=」を 主張）→ mul と 同じ 区間伝播を 積ごとに 付けて 修正。この テストは その 回帰。
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bfp_sed import BF, add, mul, tensor_mac, EQ, GE, LE, NB, M, BNAME
from nd_algebra import cd_omega

rnd = random.Random(7)
INF_CAP = 1 << 40
W = 26
OM = cd_omega(M)

def pyref(x, y):
    """純Python 整数の CD 積（厳密）。"""
    r = [0] * M
    for i in range(M):
        xi = int(x[i])
        if xi == 0: continue
        for j in range(M):
            yj = int(y[j])
            if yj: r[i ^ j] += int(OM[i, j]) * xi * yj
    return r

def sample_true(v, flag, sunk):
    """_ival と 同一の 許容集合から 整数真値を 乱択（v=0 の GE/NB は 全区間）。"""
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

def rand_bf():
    mant = [0 if rnd.random() < 0.2 else rnd.randint(-200, 200) for _ in range(M)]
    fl = [rnd.choice([EQ, EQ, GE, LE, NB]) for _ in range(M)]
    sk = [rnd.random() < 0.15 for _ in range(M)]
    tv = [sample_true(mant[k], fl[k], sk[k]) for k in range(M)]
    return BF(mant, 0, fl, sk), tv

def judge(c, ts, lies):
    for k in range(M):
        V = c.mant[k] * (1 << c.E)
        fl, sk, t = c.bound[k], c.sunk[k], ts[k]
        if fl == EQ and not sk and t != V: lies[0] += 1
        elif fl == EQ and sk and abs(t) != abs(V): lies[0] += 1
        elif fl == GE and abs(t) < abs(V): lies[0] += 1
        elif fl == LE and abs(t) > abs(V): lies[0] += 1
        if not sk and V != 0 and t != 0 and (t > 0) != (V > 0): lies[0] += 1
        lies[1] += 1

def main():
    total = 0
    for name, trials in [("add", 8000), ("mul", 1500), ("tensor_mac(2対)", 800)]:
        lies = [0, 0]
        for _ in range(trials):
            A, tA = rand_bf(); B, tB = rand_bf()
            if name == "add":
                c = add(A, B, W); ts = [tA[k] + tB[k] for k in range(M)]
            elif name == "mul":
                c = mul(A, B, W); ts = pyref(tA, tB)
            else:
                C, tC = rand_bf(); D, tD = rand_bf()
                c = tensor_mac([(A, B), (C, D)], W)
                p1 = pyref(tA, tB); p2 = pyref(tC, tD)
                ts = [p1[k] + p2[k] for k in range(M)]
            judge(c, ts, lies)
        print(f"  {name:<16} {trials:,}回 = {lies[1]:,} ラベル照合 → 嘘 {lies[0]}")
        assert lies[0] == 0, f"{name} で 嘘 {lies[0]}"
        total += lies[1]
    # 回帰: tensor_mac の フラグ落とし（旧版は 「=」= 符号の嘘）
    A = BF([10] + [0] * (M - 1), 0, [EQ] * M, [True] + [False] * (M - 1))
    B = BF([3] + [0] * (M - 1), 0)
    c = tensor_mac([(A, B)], W)
    lab = ("符号不明 " if c.sunk[0] else "") + BNAME[c.bound[0]]
    print(f"  回帰 (±10)×(3) MAC: ラベル=「{lab}」（旧: 「=」の嘘）")
    assert c.sunk[0]
    print(f"結論: {total:,} 照合 嘘0。add/mul=区間伝播で元から健全・tensor_mac=修正済み。")

if __name__ == "__main__":
    main()
