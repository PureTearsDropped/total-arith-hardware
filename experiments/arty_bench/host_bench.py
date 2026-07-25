#!/usr/bin/env python3
# top_bench の ホスト側: 0xB7 → チップ内で M=1e7 ベクタ 全速実行 → 13 バイトの 集計を 検証
import sys, time, serial

M = 10_000_000
port = sys.argv[1] if len(sys.argv) > 1 else "/dev/ttyUSB1"
ser = serial.Serial(port, 115200, timeout=10)
ser.reset_input_buffer()
t0 = time.perf_counter()
ser.write(bytes([0xB7]))
resp = ser.read(13)
wall = time.perf_counter() - t0
assert len(resp) == 13 and resp[0] == 0xB8, f"応答不正: {resp.hex()}"
mism = int.from_bytes(resp[1:5], "little")
cycles = int.from_bytes(resp[5:9], "little")
lfsr_lo = int.from_bytes(resp[9:13], "little")
# LFSR の 独立複製 (ホスト側 golden)
s = 0x9E3779B97F4A7C15
MASK = (1 << 64) - 1
for _ in range(M):
    s ^= (s << 13) & MASK
    s ^= s >> 7
    s ^= (s << 17) & MASK
ok_lfsr = (s & 0xFFFFFFFF) == lfsr_lo
print(f"実機 on-chip ベンチ (M={M:,} ベクタ):")
print(f"  不一致        : {mism} / {M:,}  (チップ内 独立検算器 = 2進キャリー伝播)")
print(f"  サイクル数    : {cycles:,}  → {M/cycles:.7f} 加算/サイクル")
print(f"  実スループット: {M/(cycles*10e-9)/1e6:.1f} M 加算/s @100MHz  (壁時計 {wall*1e3:.0f}ms 送受信込み)")
print(f"  LFSR 終状態   : {'一致 ✓ (ホスト複製と bit 一致)' if ok_lfsr else f'不一致 × {lfsr_lo:08x} vs {s & 0xFFFFFFFF:08x}'}")
sys.exit(0 if (mism == 0 and ok_lfsr) else 1)
