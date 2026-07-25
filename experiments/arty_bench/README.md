# arty_bench — on-chip measurement of the signed-digit adder, and the honesty-tax A/B

> ⚠️ AI-assisted; verify. / 生成AI使用・要検証。実施日 2026-07-25、実機 Arty A7-100T (xc7a100t)。

**One-line result:** on real silicon, a fully-verified totalized adder runs at **100.0 M
adds/s (0.9999999 adds/cycle @100 MHz), 10,000,000 vectors, 0 mismatches**, and the cost of
running an independent verifier alongside *every* operation is **zero cycles** — honesty is
wiring (area), not time. 日本語の詳細は下段。

## What was measured / 実験の内容

| # | experiment | result |
|---|---|---|
| 1 | **On-chip throughput** — LFSR generates a vector per cycle on the fabric; `sd_add2` (N=16, constant-depth) consumes 1/cycle; an **independent on-chip checker** (ordinary binary carry-propagate — a *different circuit family* than the DUT) verifies every result; UART carries only a 13-byte summary | **100.0 M adds/s @100 MHz**, 10 M vectors, **0 mismatches**, LFSR final state bit-matches the host replica |
| 2 | **Honesty-tax A/B** — two bitstreams differing ONLY in whether the verifier exists (`CHECK` parameter; identical pipeline control) | both finish 10 M vectors in **exactly the same cycle count** → verification costs **0 cycles**, +226 LUT (+37%) |
| 3 | **Fmax-ceiling penalty, and its removal** — same-stage verifier costs 17% of the ceiling (bare 165 MHz vs 141 MHz); *retiming the verifier one stage off the critical path* removes it (150 vs **159** MHz, parity within P&R noise). Silicon re-measured: 10,000,002 cycles / 10 M vectors, 0 mismatches | **verification is observation, not dependency** — results never wait for the check, so flag-raising logic can always be pipelined off-path; only value-changing honesty (saturate/clamp = one compare+mux) must stay in-path |
| 4 | **Unit datasheet** — every generated unit through place-and-route on the real target (auto-wrapper: shift-register inputs, XOR-folded outputs, so synthesis cannot prune) | see table below |

### Unit datasheet (xc7a100t, post-P&R, single-cycle combinational)

| unit | LUT | comb. delay | single-cycle Fmax |
|---|---|---|---|
| `pe24` priority encoder | 43 | 4.8 ns | 208 MHz |
| `barrel18` shifter | 131 | 5.9 ns | 168 MHz |
| `sd_mult10` | 936 | 12.2 ns | 82 MHz |
| `blocknorm` | 1,808 | 24.3 ns | 41 MHz |
| `sed_comp` (1 of 16 sedenion components) | 5,127 | 20.1 ns | 50 MHz |

Readings: one LUT6 absorbs ≈9.4 traced gates (sed_comp: 48,222 gates → 5,127 LUT); a full
16-component sedenion multiplier ≈ 82k LUT ≈ **65% of this chip — it fits**; the slowest
stage is **NORM, not the multiplies** (its parts pe24/barrel18 are ~5 ns each → 3-stage
pipeline reaches 100 MHz+).

### Lessons that surfaced / 表に出た教訓

1. **`sd_add2`'s contract is canonical digits** — feeding redundant zeros `(1,1)` broke
   1575/2000 vectors in simulation (design contract **R1** biting on silicon); one
   canonicalization mask (`p&~n / n&~p`) fixed it to 0.
2. **Verification is observation, not dependency** — the deepest structural fact of the
   honesty architecture: anything that only *raises flags* can be retimed off the critical
   path arbitrarily; the honesty tax on hardware is area only, and even the Fmax ceiling
   penalty is a design artifact, not a law.
3. Compare the GPU: the same per-operation verification costs **1.4–2.9× in time** there
   (measured in total-arith-cuda). Same semantics, different silicon: time-tax → wire-tax.

## Files

- `top_bench.sv` — bench core (LFSR → canonicalize → `sd_add2` → retimed independent
  checker → counters) + UART wrapper. `CHECK` parameter selects bare/verified.
- `tb_bench.sv` — iverilog testbench (M=2000; expects `mism=0`, LFSR state `7c3730d8`).
- `host_bench.py` — host side: one `0xB7` byte in, 13-byte summary out; replays the LFSR
  independently and cross-checks the final state.
- `build_ab.sh` — builds both variants with the fully open-source flow
  (yosys → nextpnr-xilinx → prjxray); tool locations overridable via env vars.
- `characterize.py` — auto-generates a characterization wrapper for ANY generated unit
  and reports LUT/FF/Fmax (`python3 characterize.py sed_comp ...`).
- `top_bench_bare.bit` / `top_bench_checked.bit` — prebuilt bitstreams (flash with
  `openFPGALoader -b arty_a7_100t <bit>`, then `python3 host_bench.py /dev/ttyUSB1`).
- `results/` — raw measured outputs (`measured.txt`) and place-and-route logs.

## これは何か（JP）

実機 Arty A7-100T 上での 4 実験の記録。①チップ内 LFSR→sd_add2→**独立検算器**（被試験回路と別方式の 2 進キャリー伝播）で 1000 万ベクタを全速検証——**100.0M 加算/s・不一致 0**。②検算器の有無だけが違う 2 枚の bitstream が**サイクル数完全一致**——検査の代金は時間ゼロ・面積 +226 LUT のみ（GPU では同じ検査に時間で 1.4〜2.9× 払った）。③クロック上限の罰金（同段 −17%）も、**検算＝観測≠依存**（結果は検算を待たない）なので 1 段の退去で消滅（150 vs 159 MHz）。経路内に残る誠実さは「値を変えるもの」（飽和＝比較+mux 1 段）だけ。④全生成ユニットの実測データシート（上表）——セデニオン積フル 16 成分はこのチップの 65% に収まり、最遅段は乗算でなく NORM。

再現: `./build_ab.sh`（要 openXC7 ツールチェーン・環境変数で位置指定）→ flash → `python3 host_bench.py /dev/ttyUSB1`。ビルド済み bitstream 同梱なのでツールチェーンなしでも実機再現可。
