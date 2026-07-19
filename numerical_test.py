#!/usr/bin/env python3
"""実数値での 計算試験 — 複素・四元数回転・八元数・セデニオンを 実際の float で。

測るもの:
  ① NaN/Inf を 一度も 出さないか（全域化の 保証）
  ② 健全性: 状態(境界)が 量子化入力の 厳密値について 嘘をつかないか
  ③ 精度: float64 の 参照計算に対する 相対誤差（量子化+演算誤差）
"""
import sys, os, math, time
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sed"))
import numpy as np
import bfp_sed; bfp_sed.FAST_MUL = True
from bfp_sed import BF, mul, div, add, sconj, EQ, GE, LE, NB, M
from sedenion_tensor_logic import ref_mult

# ------------------------------------------------- 実数 ⇄ ブロック浮動
def encode(vals, W):
    """実数16成分 → BF（最大成分が W ビットに入る指数、最近接丸め）。"""
    vals = [float(v) for v in vals] + [0.0] * (M - len(vals))
    amax = max((abs(v) for v in vals), default=0.0)
    if amax == 0.0:
        return BF([0] * M, 0)
    E = math.floor(math.log2(amax)) - (W - 1)
    mant = [int(round(v / 2.0 ** E)) for v in vals]
    return BF(mant, E)

def decode(bf):
    return [m * 2.0 ** bf.E for m in bf.mant]

def has_nonfinite(bf):
    return (not all(isinstance(m, int) for m in bf.mant)) or (not isinstance(bf.E, int))

def relerr(got, ref):
    num = math.sqrt(sum((g - r) ** 2 for g, r in zip(got, ref)))
    den = math.sqrt(sum(r * r for r in ref))
    return num / den if den > 0 else num

def soundness(bf, truth_rat):
    """量子化入力の 厳密値 truth_rat（Fraction）について 状態が嘘をつかないか。"""
    rep = [Fr(m) * Fr(2) ** bf.E for m in bf.mant]
    v = 0
    for k in range(M):
        t, r, b = truth_rat[k], rep[k], bf.bound[k]
        if b == EQ and t != r: v += 1
        elif b == GE and abs(t) < abs(r): v += 1
        elif b == LE and abs(t) > abs(r): v += 1
    return v


def run():
    rng = np.random.default_rng(20260719)
    nan_total = 0; sound_total = 0; t0 = time.time()

    # =============================================== ① 複素数（M=2 に埋込）チェーン
    print("① 複素数の 乗算チェーン（単位円上の 回転） vs Python complex")
    for W in (8, 12, 16, 20):
        errs = []
        for trial in range(400):
            th = rng.uniform(0, 2 * math.pi)
            z = complex(math.cos(th), math.sin(th))                  # 単位複素数
            cur = encode([z.real, z.imag], W); zc = z
            for _ in range(15):                                      # 15 回 掛ける（回転 累積）
                th2 = rng.uniform(0, 2 * math.pi)
                w = complex(math.cos(th2), math.sin(th2))
                cur = mul(cur, encode([w.real, w.imag], W), W=W); zc = zc * w
                if has_nonfinite(cur): nan_total += 1
            d = decode(cur); errs.append(relerr([d[0], d[1]], [zc.real, zc.imag]))
        print(f"   W={W:>2}: 相対誤差 中央値 {np.median(errs):.2e}, 最大 {np.max(errs):.2e}  (2^-W={2.0**-W:.1e})")

    # =============================================== ② 四元数の 回転（M=4 に埋込）実応用
    print("\n② 四元数で 3D ベクトルを 回転 vs 厳密な 回転行列（W=20, 2000 本）")
    def qmul_bf(a, b, W): return mul(a, b, W=W)
    maxerr = []; norm_drift = []
    W = 20
    for trial in range(1000):
        # ランダム単位四元数
        u = rng.normal(size=4); u /= np.linalg.norm(u)
        q = [float(u[0]), float(u[1]), float(u[2]), float(u[3])]
        v = rng.normal(size=3); vq = [0.0, float(v[0]), float(v[1]), float(v[2])]
        # 厳密: 回転行列 R で v を回す
        w0, x, y, z = q
        R = np.array([
            [1-2*(y*y+z*z), 2*(x*y-z*w0), 2*(x*z+y*w0)],
            [2*(x*y+z*w0), 1-2*(x*x+z*z), 2*(y*z-x*w0)],
            [2*(x*z-y*w0), 2*(y*z+x*w0), 1-2*(x*x+y*y)]])
        vref = R @ v
        # bfp: v' = q v conj(q)（四元数を セデニオンの 先頭4成分に 埋込）
        Q = encode(q, W); V = encode(vq, W)
        VP = mul(mul(Q, V, W=W), sconj(Q), W=W)
        d = decode(VP)
        got = [d[1], d[2], d[3]]
        maxerr.append(relerr(got, list(vref)))
        norm_drift.append(abs(math.sqrt(sum(g*g for g in got)) - np.linalg.norm(v)) / (np.linalg.norm(v)+1e-30))
        if has_nonfinite(VP): nan_total += 1
    print(f"   回転後ベクトルの 相対誤差 中央値 {np.median(maxerr):.2e}, 最大 {np.max(maxerr):.2e}")
    print(f"   ノルム保存の ずれ 中央値 {np.median(norm_drift):.2e}（回転は 長さを 保つべき）")

    # =============================================== ③ 深い八元数チェーン + 健全性
    print("\n③ 八元数の 深い積和チェーン（実数値・W=16, 300本×40段）健全性 + NaN")
    W = 16
    for trial in range(120):
        cur = encode([float(v) for v in rng.uniform(-2, 2, 8)], W)
        ctru = [Fr(m) * Fr(2) ** cur.E for m in cur.mant]
        for step in range(25):
            y = encode([float(v) for v in rng.uniform(-2, 2, 8)], W)
            yt = [Fr(m) * Fr(2) ** y.E for m in y.mant]
            cur = mul(cur, y, W=W, Emax=600)
            ctru = ref_mult(ctru, yt)
            if has_nonfinite(cur): nan_total += 1
        sound_total += soundness(cur, ctru)

    # =============================================== ④ 極端な桁範囲（実数・実応用的）
    print("\n④ 極端な桁範囲（1e−12 と 1e+12 を 同一セデニオンに・W=24, 1000本）")
    W = 24; big_err = []
    for trial in range(250):
        vals = list(rng.uniform(-1, 1, M))
        vals[0] *= 1e12; vals[1] *= 1e-12; vals[2] *= 1e6      # 24 桁 差
        bvals = [float(v) for v in rng.uniform(-3, 3, M)]
        A = encode(vals, W); B = encode(bvals, W)
        C = mul(A, B, W=W)
        ref_f = ref_mult(vals, [float(x) for x in decode(B)])          # 精度: float64 の 積
        big_err.append(relerr(decode(C), ref_f))
        ref_e = [Fr(v) * Fr(2) ** (A.E + B.E)                          # 健全性: 量子化入力の 厳密積
                 for v in ref_mult([Fr(m) for m in A.mant], [Fr(m) for m in B.mant])]
        if has_nonfinite(C): nan_total += 1
        sound_total += soundness(C, ref_e)
    print(f"   相対誤差 中央値 {np.median(big_err):.2e}, 最大 {np.max(big_err):.2e}")

    # =============================================== ⑤ 除算の 実数精度
    print("\n⑤ 除算 a/b の 実数精度 vs 厳密（W=20, 2000本）")
    W = 20; derr = []
    for trial in range(600):
        a = [float(v) for v in rng.uniform(-5, 5, M)]
        b = [float(v) for v in rng.uniform(-5, 5, M)]
        A = encode(a, W); B = encode(b, W)
        Q = div(A, B, W=W)
        bm = [round(m) for m in B.mant]; N = sum(x*x for x in bm)
        if N == 0: continue
        ac = ref_mult([Fr(m) for m in A.mant], [bm[0]] + [-x for x in bm[1:]])
        truth = [Fr(c, N) * Fr(2) ** (A.E - B.E) for c in ac]
        derr.append(relerr(decode(Q), [float(t) for t in truth]))
        if has_nonfinite(Q): nan_total += 1
    print(f"   相対誤差 中央値 {np.median(derr):.2e}, 最大 {np.max(derr):.2e}")

    print("\n" + "=" * 60)
    print(f"**NaN/Inf 発生: {nan_total} 回**（全試験通算）")
    print(f"**健全性違反: {sound_total} 件**（状態が 嘘をついた成分）")
    print(f"経過: {time.time()-t0:.0f} 秒")
    print("=" * 60)


if __name__ == "__main__":
    run()
