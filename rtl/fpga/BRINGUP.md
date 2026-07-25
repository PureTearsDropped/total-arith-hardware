# Arty A7-100T bring-up — 実機火入れ手順 / board bring-up

> Status: **first light achieved 2026-07-25 on a real Arty A7-100T** — flashed this prebuilt
> `top_arty.bit` (fully open-source flow: yosys + nextpnr-xilinx + prjxray, no Vivado) via
> openFPGALoader over usbipd-win/WSL2, then `host_test.py /dev/ttyUSB1` streamed 1000 random
> vectors through the on-silicon `sd_add2`: **1000/1000 match with the Python golden, 63
> frames/s**. The verification chain is closed end-to-end: Python gates → simulated RTL →
> real silicon.
>
> **On-chip benchmark (same day):** a fabric-side harness (LFSR → `sd_add2` →
> an *independent* binary carry-propagate checker → counters; UART carries only the
> 13-byte summary) measured **100.0 M adds/s at 100 MHz (0.9999999 adds/cycle), 10,000,000
> vectors, 0 mismatches**, LFSR final state bit-matching the host replica. 823 LUT / 469 FF
> total (0.6% of the chip), Fmax 139 MHz. Throughput/latency are now *measured*, not
> tool-estimated. One contract lesson surfaced: `sd_add2` requires **canonical digits** —
> feeding redundant zeros `(1,1)` breaks it (design contract R1, confirmed on silicon).
> Harness prototype: `research-workspace/fpga_bench/`.

## 0. What the design does

`top_arty.sv` wraps `sd_add2` (the constant-depth signed-digit adder) behind a UART:
host sends `0xA5` + two 16-bit signed-digit operands (P/N rails, little-endian),
board replies 6 bytes = the 17-bit signed-digit sum's P/N rails. LEDs show liveness.
`host_test.py` streams 1000 random vectors and checks every reply against the Python golden.

## 1. One-time host setup

```bash
pip install pyserial
# a flashing tool, either:
sudo apt install openfpgaloader        # Debian/Ubuntu name: openfpgaloader
# or build from source: https://github.com/trabucayre/openFPGALoader
```

### WSL2 note (important)

WSL2 does not see USB devices by default. Two options:

- **Windows side (easiest):** run the flash + `host_test.py` from Windows Python
  (`py -m pip install pyserial`, port is `COMx` — check Device Manager; the Arty
  enumerates two ports, the UART is usually the **second** one).
- **WSL side:** install [usbipd-win](https://github.com/dorssel/usbipd-win) on Windows, then
  `usbipd bind --busid <id>` + `usbipd attach --wsl --busid <id>`; the board appears as
  `/dev/ttyUSB0/1`.

## 2. Flash

```bash
openFPGALoader -b arty_a7_100t top_arty.bit          # volatile (SRAM) load — enough for testing
```

(Re-flash after every power cycle, or later write to SPI flash with `-f`.)

## 3. Hardware-in-the-loop test

```bash
python3 host_test.py /dev/ttyUSB1     # or COMx on Windows; try the other port on timeout
```

Expected: `ok=1000 bad=0` and a vectors/sec figure. Every reply is checked against
the Python golden model — the same golden the gate layer and the RTL simulation
are checked against, closing the chain **Python gates → simulated RTL → silicon**.

## 4. Rebuilding the bitstream (optional)

Two flows, either works:

- **Open-source (what produced the committed .bit):** yosys (synth) → nextpnr-xilinx
  (place & route, `arty_nextpnr.xdc`) → prjxray fasm→frames→bit. No license required.
- **Vivado:** `vivado -mode batch -source build.tcl` (uses `arty_a7.xdc`). Requires a
  (free WebPACK) license for xc7a100t.

## 5. Troubleshooting

- **Timeout on every vector** — wrong port (Arty has two; use the other), or board not
  flashed (LEDs static), or 115200 baud mismatch.
- **Some vectors fail** — genuinely interesting; capture the failing (a, b) pair and
  compare against `rtl/tb/test_top_arty.py` in simulation, which speaks the same protocol.
- **openFPGALoader can't find the board** — on WSL2 the USB device isn't attached (see §1);
  on Linux add udev rules or run with sudo once.
