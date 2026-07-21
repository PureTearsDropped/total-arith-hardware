# ⚠️ 生成AI使用・要検証
"""TBM プログラムの HW 脚 — ベクタファイル駆動の cocotb ドライバ。

  total-arith-cuda/run_everywhere.py が 環境変数で ベクタ(JSON)を 渡し、
  本ドライバが 自動生成 SV (sed_comp*/sd_add2/blocknorm) を 駆動して 結果を JSON で 返す。
  検算(golden 照合)は 呼び出し側の 仕事 — ここは 運転手のみ。
    TBM_VEC : 入力 JSON  {"cases": [...]}
    TBM_OUT : 出力 JSON
"""
import os, json
import cocotb
from cocotb.triggers import Timer

TOP = os.environ.get("COCOTB_TOPLEVEL", os.environ.get("TOPLEVEL", ""))
VEC = json.load(open(os.environ["TBM_VEC"]))
OUTP = os.environ["TBM_OUT"]


def to_rails(v, width):
    s = 1 if v >= 0 else -1; a = abs(v); P = N = 0
    for i in range(width):
        if (a >> i) & 1:
            if s > 0: P |= 1 << i
            else:     N |= 1 << i
    return P, N

def rails_val(P, N, width):
    return sum((((P >> i) & 1) - ((N >> i) & 1)) * (1 << i) for i in range(width))

async def settle():
    await Timer(1, unit="ns")


if TOP.startswith("sed_comp"):
    @cocotb.test()
    async def drive(dut):
        K, M = 6, 16
        res = []
        for case in VEC["cases"]:
            for i in range(M):
                P, N = to_rails(case["a"][i], K)
                getattr(dut, f"a{i}P").value = P; getattr(dut, f"a{i}N").value = N
                P, N = to_rails(case["b"][i], K)
                getattr(dut, f"b{i}P").value = P; getattr(dut, f"b{i}N").value = N
            await settle()
            w = len(dut.zP)
            res.append(rails_val(int(dut.zP.value), int(dut.zN.value), w))
        json.dump({"z": res}, open(OUTP, "w"))

elif TOP == "sd_add2":
    @cocotb.test()
    async def drive(dut):
        N = len(dut.xP)
        res = []
        for case in VEC["cases"]:
            P, Q = to_rails(case["x"], N)
            dut.xP.value = P; dut.xN.value = Q
            P, Q = to_rails(case["y"], N)
            dut.yP.value = P; dut.yN.value = Q
            await settle()
            res.append(rails_val(int(dut.zP.value), int(dut.zN.value), N + 1))
        json.dump({"z": res}, open(OUTP, "w"))

elif TOP == "blocknorm":
    @cocotb.test()
    async def drive(dut):
        Mb, Win, W = 4, 24, 6
        res = []
        for case in VEC["cases"]:
            for i in range(Mb):
                P, N = to_rails(case["m"][i], Win)
                getattr(dut, f"m{i}P").value = P; getattr(dut, f"m{i}N").value = N
            dut.Ein.value = case["Ein"]
            await settle()
            o = [rails_val(int(getattr(dut, f"o{i}P").value),
                           int(getattr(dut, f"o{i}N").value), W) for i in range(Mb)]
            fl = [int(getattr(dut, f"flag{i}").value) for i in range(Mb)]
            res.append({"o": o, "flags": fl, "Eout": int(dut.Eout.value)})
        json.dump({"blocks": res}, open(OUTP, "w"))
