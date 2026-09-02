# total-arith-hardware

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21489944.svg)](https://doi.org/10.5281/zenodo.21489944)

**Total arithmetic as a circuit — a signed-digit block-floating-point sedenion unit, built from primitive gates (AND/OR/NOT/XOR) up to synthesizable SystemVerilog, plus a "wiring = computation" fabric where swapping a wiring table turns the same unit into a different algebra.**

> ⚠️ Written with AI assistance. Independently verifiable — every claim below ships with a command that reproduces it. Verify before relying on it.

日本語の詳細は下段に。

---

## Research Use and Reproducibility / 研究利用と再現性

This software implements totalization, status propagation, singularity handling, and verification rules that differ from conventional floating-point arithmetic.

When these arithmetic semantics affect research results, it is advisable, for reproducibility, to identify this repository by URL and to record the commit ID used, together with any arithmetic settings that affect the results.

本ソフトウェアは、通常の浮動小数点演算とは異なる全域化・状態伝播・特異点処理・検証規則を実装しています。これらの算術意味論が研究結果に影響する場合は、再現性のため、本リポジトリを URL で特定し、使用したコミット ID と、結果に影響する算術設定を記録することを推奨します。

## What this is (EN)

A hardware-oriented implementation of **total arithmetic** and **"wiring = computation"**, spanning two of the four "heights" of the wider project (see *Related repositories*): the **Python gate simulation** and the **SystemVerilog (HDL) → FPGA** layers, on top of a shared integer-exact algebra core.

- **Total arithmetic** — the unit never produces `NaN` or `Inf`, and its status flags never lie. Overflow saturates to `±MAX` (flag `GE`), underflow collapses to `±MIN = ε` while preserving direction (flag `LE`), `a/0 = 0` (for a genuine zero only), and there are two kinds of zero. Adversarial sweeps produce zero `NaN`/`Inf` and zero flag lies.
- **Wiring = computation** — the multiplication is described by a *structure tensor* / wiring table. Swap the table and the same gate graph computes a different algebra: complex, quaternion, sedenion, matrix product, cyclic convolution. A registry of 19 patterns plus an inverse designer (ask for a matrix block → get the minimal group that contains it). **The wiring normal form is ternary** (coefficients ⊆ {−1, 0, +1} — TBM_SPEC §1.5 in total-arith-cuda): a permanent gatekeeper checks all 18 bilinear wirings, so the linear stages are exact adds/routes (no coefficient rounding anywhere), `R` honestly counts the true multiplies, and the gate check now covers **every** wiring (the generic `(U,V,W)` gate expansion was always universal; big wirings are spot-checked for time budget, none skipped). Custom tables with structural zeros (dual numbers, Grassmann) get **dead-product pruning** (`prune_uvw`, output-invariant — measured waste before: 25–44% of gates); merging is deliberately *not* done (order-blind merging annihilates the antisymmetric part — where non-commutativity lives — and breaks even commutative algebras; measured).
- **Signed-digit block floating point** — mantissa digits are ternary `{−1, 0, +1}` (Avizienis 1961) sharing one base-2 exponent per 16-component sedenion (a signed-ternary cousin of MXFP/microscaling). Sign symmetry means **sign flip = swapping `+`/`−` digits = pure wiring** (zero gates); the sign is carried by the mantissa, so only a "sign-unknown" bit remains as separate status.
- **From gates to HDL with no hand-transcription** — the SystemVerilog is *auto-emitted* by tracing the audited Python gate graph, so the HDL is the same gate graph by construction. Verified against the Python golden model with Icarus + cocotb.
- **Elementary functions as gates** (`remez/`, v1.2.0) — exp / expm1 / log / sqrt / rsqrt are *range reduction on the exponent bus × a Remez coefficient tape (Estrin) × a one-sided exit*; the gate graphs are bit-identical to an integer-only spec, flags included, and the spec was swept over **every finite float32 input** (5.78 × 10¹⁰ checks, 0 violations, worst one-sided error 1.0039 ulp, nearest mode ≤ 0.5001 ulp). Details in the `remez/` section below and `remez/README.md`.

### Architecture

```
Layer 1  algebra core (integer-exact, Python)     sd2_core, nd_algebra, matrix_algebra, bfp_sed
Layer 2  gate implementation (from AND/OR/NOT/XOR) gate_bilinear, gate_exponent, gate_fast, multi_add, mul_fused, remez/ ...
Layer 3  wiring (= choice of computation)          wiring_registry, wiring_designer, representation_lens, ...
Layer 4  HDL (SystemVerilog + Icarus + cocotb)     rtl/  (auto-emitted SV, testbenches, Arty A7 FPGA target)
```

### Reproduce

```bash
python3 -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt
./verify.sh                 # runs the pure-Python + gate self-tests, reports PASS/FAIL
```

- **Pure-Python / gate self-tests** — no special hardware. Reproduced here: **all 17 green, zero violations.**
- **HDL (virtual hardware / RTL simulation)** — requires `iverilog` (or `verilator`) + `cocotb`. `./verify_hdl.sh` runs every SystemVerilog module as *simulated hardware* (no FPGA) and checks it against the Python golden. Reproduced here with **iverilog 12.0 + cocotb 2.0.1: all 8 modules green** — `gate9` `compress3` `sd_add2` `sd_mult10` `pe24` `barrel18` `blocknorm` `sed_comp` (`TESTS=1 PASS=1 FAIL=0` each). Note: clear `rtl/tb/sim_build` between toplevels (`verify_hdl.sh` does this) — a stale build silently reuses the previous module.
- **FPGA — verified on real silicon (2026-07-25)**: the included prebuilt bitstream (`rtl/fpga/top_arty.bit`, Arty A7-100T, UART ⇔ `sd_add2`, fully **open-source flow**: yosys + nextpnr-xilinx + prjxray — no Vivado) was flashed to a physical board and `host_test.py` streamed 1000 random vectors through the on-silicon signed-digit adder: **1000/1000 match with the Python golden (63 frames/s)** — the chain Python gates → simulated RTL → real silicon is closed end-to-end. An on-chip benchmark then measured **100.0 M adds/s at 100 MHz (10,000,000 vectors, 0 mismatches)** and an A/B experiment showed **per-operation verification costs 0 cycles** (area only; the GPU pays 1.4–2.9× in time for the same semantics) — sources, bitstreams, raw logs and write-up in **`experiments/arty_bench/`**. A Vivado `build.tcl` is also provided for those who prefer it.

### Design contract (see `GATE_CONDITIONS.md`)

- **K1** coefficients, shift amounts and wiring tables are *synthesis-time constants* (passing a signal makes a hidden comparator).
- **K2** widths are fixed for the worst case (no data-dependent wiring width).
- **R1** canonicalize before any magnitude comparison on a redundant representation.
- **R2** truncation-as-boundary is not automatic — you buy it with monotonicity-aware directed rounding.
- **R3** re-basing base-2 (change `E` = relabel digits) is lossless.

`explorations/` holds the working notebooks behind the semantics (sign survival, `±MIN` bounds, zero conventions). They are exploratory, not part of the verified core.

---

## これは何か（JP）

**全域算術**と**「配線＝計算」**の、ハードウェア寄り実装。プロジェクト全体の「4つの高さ」（*Related repositories* 参照）のうち、**Pythonゲートシミュレーション**と **SystemVerilog(HDL)→FPGA** の層を、整数厳密の代数コアの上に載せてある。

- **全域算術** — `NaN`/`Inf` を決して作らず、状態フラグが嘘をつかない。溢れ→`±MAX`(`GE`)、潰れ→`±MIN=ε`（向き保持・`LE`）、`a/0=0`（本物の0のときだけ）、二種類の0。敵対的スイープで `NaN`/`Inf` 生成ゼロ・フラグの嘘ゼロ。
- **配線＝計算** — 乗算を*構造テンソル*（配線表）で記述。表を差し替えると同じゲートグラフが別の代数（複素・四元数・セデニオン・行列積・巡回畳み込み）になる。19パターンのレジストリ＋逆設計（欲しい行列ブロック→それを含む最小の群）。**配線正規形は三値 {−1,0,+1}**（total-arith-cuda の TBM_SPEC §1.5）——三値門番が bilinear 18/18 を恒久検査。線形段は加算と経路だけ（係数の丸めが存在しない）、`R` は真の乗算数の正直な請求書、gate 検査は全配線をカバー。表に構造的な 0 がある代数（双対数・Grassmann）は**死に積の刈り込み** `prune_uvw`（出力不変・刈る前は 25〜44% のゲートが無駄と実測）。**併合は意図的にしない**——順序無視の併合は反対称部（非可換の住処）を消し、可換代数でも壊れる（実測）。
- **符号つき3値ブロック浮動** — 仮数の桁が3値 `{−1,0,+1}`（Avizienis 1961）で、16成分セデニオンが指数を1つ共有（MXFP/microscalingの符号つき3値版）。符号対称なので **符号反転＝`+`/`−`桁の入替＝純配線**（ゲートゼロ）。符号は仮数が運び、残る状態は「符号不明」1ビットだけ。
- **ゲート→HDLを手写しゼロで** — SystemVerilogは監査済みPythonゲートグラフをトレースして*自動生成*。だからHDLは定義上同一のゲートグラフ。Icarus+cocotbでPython golden と照合。
- **初等関数をゲートで**（`remez/`、v1.2.0）— exp / expm1 / log / sqrt / rsqrt は「指数バスでの縮小 × Remez 係数テープ（Estrin）× 片側出口」。ゲート版は整数だけの spec とフラグ込みで bit 一致し、spec は **float32 の有限な全入力を総当たり**（578 億回・違反 0・片側の最悪 1.0039 ulp・最近接は 0.5001 ulp 以内）。詳細は下段の `remez/` 節と `remez/README.md`。

### 再現方法

上記 *Reproduce* のコマンド。純Python/ゲートのself-testは特別なハード不要（ここでは全緑・違反0を実測）。HDLは `iverilog`+`cocotb`。FPGA は**ビルド済み bitstream 同梱**（`rtl/fpga/top_arty.bit`、完全オープンソースフロー yosys+nextpnr-xilinx+prjxray 製・Vivado 不要、手順は `rtl/fpga/BRINGUP.md`）。

### ファイル地図

- **Layer 1 代数コア**: `sd2_core.py` `sd2_gates.py` `nd_algebra.py` `matrix_algebra.py` `bfp_sed.py`
- **Layer 2 ゲート**: `gate_bilinear.py` `gate_exponent.py` `gate_bfp.py` `gate_bfp2.py` `gate_fast.py` `multi_add.py` `mul_fused.py` `newton_recip.py`
- **Layer 3 配線**: `wiring_registry.py` `wiring_designer.py` `wiring_zoo.py` `wiring_patterns.py` `representation_lens.py` `bilinear_unit.py` `parallel_array.py`
- **Layer 4 HDL**: `rtl/`（自動生成SV・cocotb TB・Arty A7一式）
- **検証**: `interval_*.py` `numerical_test*.py` `stress_test*.py` `div_*.py`、`sed/`（セデニオン零因子ロジック）
- **初等関数**: `remez/`（Remez 係数テープ × Estrin 積和 × 片側出口 — exp/expm1/log/sqrt/rsqrt をゲートまで、sin/cos は spec まで。float32 幅は全入力総当たり済み。`remez/README.md`）
- **研究枝の記録**: `RESEARCH.md`（テーマごとの予想・測定・撤回）、`ternary_mark/`（三値語の「不明」桁と全域語）、`four_value/`（四値桁 {0,+1,−1,?} の半分選び定理と素子 2 つの差し替え）
- **文書**: `SPEC.md` `GATE_CONDITIONS.md` `PROCESS.md` `FINDINGS.md`
- **explorations/**: 意味論の作業ノート（探索的・検証済みコアではない）

---

## Related repositories

The same two ideas — *total arithmetic* and *wiring = computation* — are implemented independently at other "heights":

- **[total-arith-cuda](../../total-arith-cuda)** — the GPU (torch/CUDA) height: total arithmetic + swappable structure tensor as tensor kernels.
- **[varpro-powersum-nn](../../varpro-powersum-nn)** — the learning height: where total arithmetic pays off in training (totalized gradients keep contaminated data from killing a fit).

These three are, in effect, **three backends of one total-arithmetic contract** (CPU / GPU / hardware). A planned next step is a pluggable backend so a model can run its total arithmetic on any of them.

## 興味を持ったら / If this interests you

これは利用条件ではありません。ただの声かけです — もしこの方向性に興味を持って、議論したい・一緒に発展させたい・仕事として相談したい等があれば、この repo の Issue で気軽にどうぞ。（連絡は GitHub 経由で OK、本名は不要です。）

*Not a term of use — just an open door. If this direction interests you and you'd like to discuss it, develop it together, or talk about it as work, feel free to open an Issue. Reach me via GitHub; no real name needed.*

## License

Zero-Clause BSD (0BSD). See `LICENSE`. Do whatever you want; no attribution required.

### `gate_series.py` — transcendental functions in gates (coefficient tape × chained (U,V,W) units)

`exp`/`sin` as a FIXED WIRING — no branches, no loops, no matrix inverse: a compile-time
coefficient tape (K1: dyadic constants `round(2^P/k!)·2^-P`) drives a chain of
`bf_bilinear_unit`s; scaling `2^-s` is an exponent subtract (base-2 reattachment, **zero
gates**); swapping the tape swaps the function (exp ↔ sin) on the same skeleton. Two-level
honesty: ① the gate graph equals its exact rational (Fraction) spec **exactly** (EXACT mode,
no normalization — asserted for quaternion exp/sin and sedenion exp), ② the spec's distance
to the true function is measured (tape quantization + truncation), with a bits-in→accuracy
dial (1.9e-2 → 5.1e-3). Measured cost law of EXACT mode: 86M gates (order 4) → 263M
(order 6) → 4.3B (order 8 + 2 squarings, mantissa width 477) — squaring is the gate
monster, so practical circuits need per-stage `block_normalize` (honest rounding + ge/le
flags); connecting that to `bfp_sed`'s interval machinery is the declared next step.

### `bfp_series.py` — the normalized series unit (the declared next step, done)

Answers `gate_series.py`'s measured 4.3B-gate explosion: chain `bfp_sed`'s honest-rounding
`mul`/`add` (oracle-verified flag algebra) into the same tape-driven series — per-stage
normalization to W digits makes cost linear in order. The END-TO-END claim is then
falsified against the ideal (unrounded, exact-Fraction) computation: output flags
{=, ≥, ≤, ±} must admit the ideal value per component. Results: 0 violations across
sedenion exp (20 seeds × order 8 × 2 squarings), a W dial (W=10→28: error 6.7e-2 → 1.7e-7,
sound at every width), REAL saturation (Emax=−4 pins 15/16 components at ±MAX — still
sound), and the sin tape on the same skeleton. Honest finding: in dense long chains the
flags degrade almost entirely to "no-bound" — zero lies but thin information (the
composition-level echo of the 0.5% dense retention rule; quantitative intervals want the
4-value digit / interval representations as the next structure).

### `gate_solve.py` — the solve family in gates, M/N/O and U/V/W preserved

`nsolve` (equation-solving division `L_a⁺x`) lowered to the gate world on the same rails as
the exp units: **O** = the Ben-Israel tape `X←X(2I−LX)` (K fixed, no branches, no division —
the initial scale is an exponent shift), **V** = matrix products as `bf_bilinear_unit` ×
`_matmul_uvw` (matmul is just another (U,V,W)), **U** = `L_a` construction (pure lincomb
wiring) + `2^{-s}` init, **W** = the two-tier verification (exact-solve / SING least-squares
/ SING|INEXACT) fused into the unit — "never pretend an inconsistent system was solved" at
the gate-model level; **N** = swap `cd_omega(M)`. Two parts, same precedent as
gate_series/bfp_series: **Part A** pure-gate EXACT mode — gate graph ≡ exact rational spec
(Fraction equality asserted, quaternion K=2, 93.1M gates) with the measured width doubling
[6→34→90] that makes normalization mandatory; **Part B** normalized BFP mode — W-digit
block-float matrices, cost LINEAR in K (sedenion K=25: 204,800 integer multiplies/solve),
final verification in exact integer arithmetic (zero-error claims about the quantized
problem): regular solve matches pinv at 1.4e-7 (W=24), zero divisor honestly SING with
normal-equation residual 8e-8, and a width→accuracy dial (W=12→32: 6e-4→5e-8). The scalar
1×1 case is `newton_recip.py` — this is its matrix completion.

### `remez/` — elementary functions as gates (range reduction × Remez tape × one-sided exit)

exp / expm1 / log / sqrt / rsqrt (and sin / cos at spec level) are defined as **range reduction on
the exponent bus + a Remez polynomial evaluated by Estrin + a one-sided exit** — no Newton loop,
no truncated series. The Remez fit runs offline (mpmath); the machine carries only the coefficient
tape (K1: synthesis-time constants — a constant × signal is the wiring of the constant's set bits)
and the reduction wiring. The error bound of each quantized tape is proved with `Fraction`
arithmetic alone (`supnorm_exact`). Every function returns three modes (nearest / lower bound /
upper bound) with the total-arithmetic flags; an integer-only spec (`funcs_spec.py`) is the
reference, and the gate versions (`gate_funcs.py`, `gate_logroot.py`) are **bit-identical to it,
flags included**: exp/expm1 on 535 adversarial points (f32 and f64), log/sqrt/rsqrt on 1461 (f32)
/ 1485 (f64) points, 0 host asserts. Cost (f32 / f64 gates): exp+expm1 all six outputs 270,592 /
687,173 (depth 791 / 863); log 338,660 / 1,095,224; sqrt 387,652 / 1,366,904; rsqrt 471,319 /
1,536,011. An **exhaustive float32 sweep of the spec** (every finite float32 word: 4.28 × 10⁹ for
exp/expm1, the 2.14 × 10⁹ positive words for the other five; truth = float64 with mpmath re-checks
near the boundaries) — 1.93 × 10¹⁰ (function, input) pairs × 3 modes = 5.78 × 10¹⁰ checks — found
**0 violations** of the one-sided bounds, 0 sign errors, 0 false exactness claims and 0 missed
saturations; the worst one-sided error is 1.0005 ulp (exp, sin, cos), 1.0010 (log, expm1), 1.0009
(sqrt), 1.0039 (rsqrt) against the claimed 1.0007 / 1.0015 / 1.0015 / 1.0015 / 1.0059, and the
nearest mode stays within 0.50003–0.50010 ulp (claimed 0.5002–0.5020; 4.5 k – 65 k results per
function are not the correctly rounded float32, which is never claimed). The worst point of every
one of the seven functions is an argument whose result lands just below a power of two — the
finer binade doubles the ulp count of the same relative error — and the reduction scheme
(exponent bus, Payne–Hanek, odd/even tape) leaves no visible trace. Two disciplines came out of
it (`GATE_CONDITIONS.md`): **R4** a boundary that exists only on the gate side needs a
boundary-value unit test (two window bugs passed the 535-point differential test), and **R5**
count the sticky premise ("dropped amount < 1 unit" fails when two things drop at once). Tables,
reproduce commands and the list of withdrawn predictions: `remez/README.md`, `RESEARCH.md`
(theme 3).
