# ⚠️ 生成AI使用・要検証
"""自動生成 SV（emit_sv.py）を 監査済み Python golden と 突き合わせる cocotb テスト。"""
import os, sys, random
import cocotb
from cocotb.triggers import Timer

ARITH = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ARITH)

TOP = os.environ.get("COCOTB_TOPLEVEL", os.environ.get("TOPLEVEL", "sd_mult10"))
rnd = random.Random(20260809)


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


if TOP == "sd_mult10":
    @cocotb.test()
    async def mult_random(dut):
        K = 10; lim = (1 << (K - 1)) - 1
        for a, b in [(0, 0), (lim, lim), (-lim, lim), (1, -1)] + \
                    [(rnd.randint(-lim, lim), rnd.randint(-lim, lim)) for _ in range(300)]:
            xP, xN = to_rails(a, K); yP, yN = to_rails(b, K)
            dut.xP.value = xP; dut.xN.value = xN
            dut.yP.value = yP; dut.yN.value = yN
            await settle()
            got = rails_val(int(dut.zP.value), int(dut.zN.value), 33)
            assert got == a * b, f"mult({a},{b}) = {got} ≠ {a*b}"
        dut._log.info("sd_mult10: 304 ケース == a·b ✓")


elif TOP == "pe24":
    @cocotb.test()
    async def pe_random(dut):
        from gate_fast import priority_encoder_fast
        from gate_bilinear import to_sd, new_counter
        n, EW = 24, 6
        for v in [0, 1, -1, (1 << 22) + 5] + [rnd.randint(-(1 << 22), 1 << 22) for _ in range(300)]:
            P, N = to_rails(v, n)
            dut.dP.value = P; dut.dN.value = N
            await settle()
            Lg, ng, _ = priority_encoder_fast(to_sd(v, n), EW, new_counter())
            Lg = sum(int(b) << i for i, b in enumerate(Lg))
            assert int(dut.L.value) == Lg and int(dut.none_o.value) == int(ng), \
                f"pe24({v}): L={int(dut.L.value)}≠{Lg} none={int(dut.none_o.value)}≠{int(ng)}"
        dut._log.info("pe24: 304 ケース == golden ✓")


elif TOP == "barrel18":
    @cocotb.test()
    async def barrel_random(dut):
        n = 18
        for _ in range(400):
            v = rnd.randint(-(1 << 15), 1 << 15); s = rnd.randint(0, 20)
            P, N = to_rails(v, n)
            dut.dP.value = P; dut.dN.value = N; dut.S.value = s
            await settle()
            got = rails_val(int(dut.oP.value), int(dut.oN.value), n)
            mag = abs(v) >> s if s < 64 else 0
            want = mag if v >= 0 else -mag
            wdrop = 1 if (abs(v) & ((1 << min(s, 63)) - 1)) != 0 else 0
            assert got == want and int(dut.dropped.value) == wdrop, \
                f"barrel({v}>>{s}) = {got},{int(dut.dropped.value)} ≠ {want},{wdrop}"
        dut._log.info("barrel18: 400 ケース == 期待値 ✓")


elif TOP == "blocknorm":
    @cocotb.test()
    async def blocknorm_random(dut):
        from gate_fast import block_normalize_g_fast
        from gate_exponent import bus_const, bus_val
        from gate_bilinear import to_sd, from_sd, new_counter
        M, Win, W, Emax, EW = 4, 24, 6, 20, 12
        for _ in range(120):
            E0 = rnd.randint(0, 12)
            vals = [rnd.randint(-5, 5) * (10 ** rnd.randint(0, 5)) for _ in range(M)]
            for i, v in enumerate(vals):
                P, N = to_rails(v, Win)
                getattr(dut, f"m{i}P").value = P
                getattr(dut, f"m{i}N").value = N
            dut.Ein.value = E0
            await settle()
            og, Eg, fg = block_normalize_g_fast([to_sd(v, Win) for v in vals],
                                                bus_const(E0, EW), W, Emax, new_counter())
            assert int(dut.Eout.value) == bus_val([int(b) for b in Eg]) % (1 << EW)
            for i in range(M):
                sv = rails_val(int(getattr(dut, f"o{i}P").value),
                               int(getattr(dut, f"o{i}N").value), W)
                assert sv == from_sd(og[i]), f"成分{i}: SV={sv} py={from_sd(og[i])}"
                fl = int(getattr(dut, f"flag{i}").value)
                ge, le, _ = fg[i]
                assert fl == (int(ge) | (int(le) << 1)), f"flag{i}: {fl} ≠ ge={ge},le={le}"
        dut._log.info("blocknorm: 120 ブロック（仮数・指数・フラグ ge/le/ε）== golden ✓")


elif TOP == "sed_comp":
    @cocotb.test()
    async def sed_comp_random(dut):
        from mul_fused import group_component
        from gate_bilinear import to_sd, from_sd, new_counter
        from nd_algebra import cd_omega, ref_mult_M
        M, K, k = 16, 6, 1
        OM = cd_omega(M)
        OMl = [[int(OM[i, j]) for j in range(M)] for i in range(M)]
        for _ in range(60):
            a = [rnd.randint(-9, 9) for _ in range(M)]
            b = [rnd.randint(-9, 9) for _ in range(M)]
            for i in range(M):
                P, N = to_rails(a[i], K)
                getattr(dut, f"a{i}P").value = P; getattr(dut, f"a{i}N").value = N
                P, N = to_rails(b[i], K)
                getattr(dut, f"b{i}P").value = P; getattr(dut, f"b{i}N").value = N
            await settle()
            got = rails_val(int(dut.zP.value), int(dut.zN.value), 29)
            ref = ref_mult_M(a, b, OM, M)[k]
            assert got == ref, f"sed_comp: SV={got} ≠ 参照={ref}"
        dut._log.info("sed_comp: 60 ケース == セデニオン積成分（16積の融合MAC） ✓")
