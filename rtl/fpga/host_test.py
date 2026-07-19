#!/usr/bin/env python3
# ⚠️ 生成AI使用・要検証
"""実機テスト（hardware-in-the-loop）— Arty A7 上の sd_add2 に UART で ベクタを 流し、
Python golden と 照合。シミュレーション（test_top_arty.py）と 同一プロトコル。

  使い方: pip install pyserial して
    python3 host_test.py /dev/ttyUSB1        # Arty の UART は 2 ポート目が 多い
"""
import sys, random, time

def to_rails(v, width):
    s = 1 if v >= 0 else -1; a = abs(v); P = N = 0
    for i in range(width):
        if (a >> i) & 1:
            if s > 0: P |= 1 << i
            else:     N |= 1 << i
    return P, N

def rails_val(P, N, width):
    return sum((((P >> i) & 1) - ((N >> i) & 1)) * (1 << i) for i in range(width))

def main(port):
    import serial
    ser = serial.Serial(port, 115200, timeout=2)
    time.sleep(0.2); ser.reset_input_buffer()
    rnd = random.Random()
    Nn = 16; lim = (1 << (Nn - 1)) - 1
    ok = bad = 0
    t0 = time.time()
    trials = 1000
    for i in range(trials):
        a = rnd.randint(-lim, lim); b = rnd.randint(-lim, lim)
        xP, xN = to_rails(a, Nn); yP, yN = to_rails(b, Nn)
        ser.write(bytes([0xA5, xP & 255, xP >> 8, xN & 255, xN >> 8,
                                yP & 255, yP >> 8, yN & 255, yN >> 8]))
        rx = ser.read(6)
        if len(rx) != 6:
            print(f"[{i}] タイムアウト（受信 {len(rx)}B）— 配線/ポートを 確認"); bad += 1; break
        zP = rx[0] | (rx[1] << 8) | (rx[2] << 16)
        zN = rx[3] | (rx[4] << 8) | (rx[5] << 16)
        got = rails_val(zP, zN, 17)
        if got == a + b: ok += 1
        else:
            bad += 1
            print(f"[{i}] 不一致: {a}+{b} → FPGA {got}（期待 {a+b}）")
    dt = time.time() - t0
    print(f"\n実機 sd_add2: 一致 {ok}/{ok+bad}  ({dt:.1f}s・{(ok+bad)/dt:.0f} フレーム/s)")
    print("✓ 実シリコン上の 定数深さ SD 加算器が golden と 一致" if bad == 0 else "× 要調査")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/dev/ttyUSB1")
