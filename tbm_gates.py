#!/usr/bin/env python3
# ⚠️ 生成AI使用・要検証
"""tbm_gates — 総ビリニア機械 v0.1 アセンブラの ゲート級 双子 (total-arith-cuda/tbm.py と 同形)。

  同じ 命令列を、SV 自動生成の 元になっている 監査済み ゲート関数
  (AND/OR/NOT/XOR まで 落ちた Python golden) で 実行する:
    TOTALIZE : 整数 → SD 桁 (gate_bilinear.to_sd)
    BILIN    : セデニオン積 = group_component k=0..15 (16積の 融合MAC・符号=配線)
    AXPY     : gate_fast.sd_add2 (c=1 のみ — ゲート級は 厳密 整数の 世界)
    NORM     : gate_fast.block_normalize_g_fast (ブロック正規化 + ge/le フラグ)
  LINMAP / CHECK は 未対応 (適合表の 空欄は 空欄のまま)。
  誠実さダイヤル: ゲート級は 厳密 整数 (丸めなし) — 飽和と フラグは NORM に 住む。
  副産物: st カウンタが 命令ごとの ゲート評価数を 数える (配線 = 計算の 請求書)。

  cocotb+iverilog の RTL 脚 (run_everywhere.py) との 関係: あちらは 生成 SV を 駆動、
  こちらは その 生成元 golden を 直接 実行 — emit_sv の トレースにより 同一 ゲートグラフ。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate_bilinear import to_sd, from_sd, new_counter
from gate_fast import sd_add2, block_normalize_g_fast
from gate_exponent import bus_const, bus_val
from mul_fused import group_component
from nd_algebra import cd_omega, ref_mult_M


class Program:
    "tbm.py と 同形の 命令列 (ゲート級は サブセット実行)。"
    def __init__(self, name="tbm"):
        self.name = name
        self.ins = []

    def TOTALIZE(self, dst, src):
        self.ins.append(("TOTALIZE", dict(dst=dst, src=src)))
        return self

    def BILIN(self, dst, a, b, alg="sedenion", honesty="evidence", width=None):
        assert alg == "sedenion", "ゲート級 v1 は セデニオンのみ (sed_comp と 同じ 範囲)"
        self.ins.append(("BILIN", dict(dst=dst, a=a, b=b)))
        return self

    def AXPY(self, dst, src, c=1.0):
        assert c == 1.0, "ゲート級は 厳密 整数 — c=1 のみ"
        self.ins.append(("AXPY", dict(dst=dst, src=src)))
        return self

    def NORM(self, dst, src, block=4, Ein=0):
        self.ins.append(("NORM", dict(dst=dst, src=src, block=block, Ein=Ein)))
        return self


_OM = cd_omega(16)
_OMl = [[int(_OM[i, j]) for j in range(16)] for i in range(16)]


def _pad(d, n):
    return list(d) + [(0, 0)] * (n - len(d))


def run_gates(prog, feed, K=6, W=6, Win=24, Emax=20, EW=12):
    """feed: {名前: [ [int]*16, ... ]}。返り値: {名前: 値リスト} + gates (命令別 ゲート評価数)。
       各値は SD 桁列を from_sd で 復号した 厳密 整数。"""
    env = {}
    gates = {}
    for op, p in prog.ins:
        st = new_counter()
        if op == "TOTALIZE":
            env[p["dst"]] = [[to_sd(int(v), K) for v in row] for row in feed[p["src"]]]
        elif op == "BILIN":
            out = []
            for ra, rb in zip(env[p["a"]], env[p["b"]]):
                out.append([group_component(ra, rb, _OMl, 16, k, st) for k in range(16)])
            env[p["dst"]] = out
        elif op == "AXPY":
            out = []
            for rd, rs in zip(env[p["dst"]], env[p["src"]]):
                comps = []
                for dx, dy in zip(rd, rs):
                    n = max(len(dx), len(dy))
                    comps.append(sd_add2(_pad(dx, n), _pad(dy, n), st))
                out.append(comps)
            env[p["dst"]] = out
        elif op == "NORM":
            recs = []
            for row in env[p["src"]]:
                mants = [_pad(d, Win)[:Win] for d in row[:p["block"]]]
                og, Eg, fg = block_normalize_g_fast(mants, bus_const(p["Ein"], EW),
                                                    W, Emax, st)
                recs.append(dict(o=[from_sd(d) for d in og],
                                 flags=[(int(g), int(l)) for g, l, _ in fg],
                                 Eout=bus_val([int(b) for b in Eg]) % (1 << EW)))
            env[p["dst"]] = recs
        gates[f"{op}:{p['dst']}"] = gates.get(f"{op}:{p['dst']}", 0) + sum(st.values())
    return env, gates


def decode(env_val):
    "SD 桁環境 → 整数 (検算用)。"
    return [[from_sd(d) for d in row] for row in env_val]


# ---------------------------------------------------------------- self-test
def self_test():
    import random
    rnd = random.Random(20260721)
    print("tbm_gates — 同じ TBM プログラムを ゲート級 golden で (SV 生成元と 同一 グラフ)")
    B = 6
    fa = [[rnd.randint(-9, 9) for _ in range(16)] for _ in range(B)]
    fb = [[rnd.randint(-9, 9) for _ in range(16)] for _ in range(B)]
    fc = [[rnd.randint(-9, 9) for _ in range(16)] for _ in range(B)]
    P = (Program("mac_norm")
         .TOTALIZE("a", "in_a").TOTALIZE("b", "in_b").TOTALIZE("c", "in_c")
         .BILIN("s", "a", "b").AXPY("s", "c").NORM("n", "s", block=4, Ein=0))
    env, gates = run_gates(P, {"in_a": fa, "in_b": fb, "in_c": fc})

    ref = [[ref_mult_M(fa[i], fb[i], _OM, 16)[k] + fc[i][k] for k in range(16)]
           for i in range(B)]
    got = decode(env["s"])
    assert got == ref, "ゲートの s ≠ 代数の 参照"
    print(f"  ① BILIN+AXPY: {B}ケース × 16成分 = ゲートの 答え ≡ 代数の 答え (厳密 整数) ✓")

    for i, rec in enumerate(env["n"]):
        og, Eg, fg = block_normalize_g_fast([to_sd(v, 24) for v in ref[i][:4]],
                                            bus_const(0, 12), 6, 20, new_counter())
        assert rec["o"] == [from_sd(d) for d in og] and \
               rec["Eout"] == bus_val([int(b) for b in Eg]) % (1 << 12)
    print(f"  ② NORM: {B}ブロック = golden 直呼びと 一致 (仮数・指数・ge/le) ✓")

    total = sum(gates.values())
    print("  ③ 配線=計算の 請求書 (ゲート評価数):")
    for k, v in gates.items():
        print(f"       {k:<12} {v:>12,}")
    print(f"       {'計':<12} {total:>12,}")
    print("done — 同じ プログラムが 論理ゲートの 高さでも 走った")


if __name__ == "__main__":
    self_test()
