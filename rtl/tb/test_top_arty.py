# ⚠️ 生成AI使用・要検証
"""top_arty の エンドツーエンド シミュレーション — UART 波形レベルで フレームを 流し、
実機と 同じ プロトコルで sd_add2 の 答えが 返るか（ハード到着前の 事前検証）。DIV=20 に 短縮。"""
import os, sys, random
import cocotb
from cocotb.clock import Clock
from cocotb.triggers import Timer, RisingEdge

ARITH = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ARITH)

DIV = 20
TCLK = 10  # ns（100MHz）
TBIT = DIV * TCLK

def to_rails(v, width):
    s = 1 if v >= 0 else -1; a = abs(v); P = N = 0
    for i in range(width):
        if (a >> i) & 1:
            if s > 0: P |= 1 << i
            else:     N |= 1 << i
    return P, N

def rails_val(P, N, width):
    return sum((((P >> i) & 1) - ((N >> i) & 1)) * (1 << i) for i in range(width))


async def send_byte(dut, b):
    dut.uart_txd_in.value = 0                     # start
    await Timer(TBIT, unit="ns")
    for i in range(8):
        dut.uart_txd_in.value = (b >> i) & 1
        await Timer(TBIT, unit="ns")
    dut.uart_txd_in.value = 1                     # stop
    await Timer(TBIT, unit="ns")

async def recv_byte(dut, timeout_bits=20000):
    for _ in range(timeout_bits * DIV):
        await RisingEdge(dut.CLK100MHZ)
        if int(dut.uart_rxd_out.value) == 0:
            break
    else:
        raise AssertionError("start ビットが 来ない")
    await Timer(TBIT + TBIT // 2, unit="ns")      # スタート + 半ビットで ビット0 中央
    b = 0
    for i in range(8):
        b |= int(dut.uart_rxd_out.value) << i
        await Timer(TBIT, unit="ns")
    return b


@cocotb.test()
async def hw_in_the_loop_protocol(dut):
    rnd = random.Random(20260811)
    cocotb.start_soon(Clock(dut.CLK100MHZ, TCLK, unit="ns").start())
    dut.uart_txd_in.value = 1
    await Timer(50 * TCLK, unit="ns")             # リセット解除待ち

    N = 16; lim = (1 << (N - 1)) - 1
    cases = [(0, 0), (lim, -lim), (1, -1)] + \
            [(rnd.randint(-lim, lim), rnd.randint(-lim, lim)) for _ in range(12)]
    for a, b in cases:
        xP, xN = to_rails(a, N); yP, yN = to_rails(b, N)
        payload = bytes([0xA5,
                         xP & 0xFF, xP >> 8, xN & 0xFF, xN >> 8,
                         yP & 0xFF, yP >> 8, yN & 0xFF, yN >> 8])
        for byte in payload:
            await send_byte(dut, byte)
        rx = [await recv_byte(dut) for _ in range(6)]
        zP = rx[0] | (rx[1] << 8) | (rx[2] << 16)
        zN = rx[3] | (rx[4] << 8) | (rx[5] << 16)
        got = rails_val(zP, zN, 17)
        assert got == a + b, f"UART frame ({a},{b}): got {got} ≠ {a+b}"
    dut._log.info(f"top_arty: {len(cases)} フレーム UART 往復 == a+b ✓（実機と 同一プロトコル）")
