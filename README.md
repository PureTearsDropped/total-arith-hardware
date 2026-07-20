# total-arith-hardware

**Total arithmetic as a circuit — a signed-digit block-floating-point sedenion unit, built from primitive gates (AND/OR/NOT/XOR) up to synthesizable SystemVerilog, plus a "wiring = computation" fabric where swapping a wiring table turns the same unit into a different algebra.**

> ⚠️ Written with AI assistance. Independently verifiable — every claim below ships with a command that reproduces it. Verify before relying on it.

日本語の詳細は下段に。

---

## What this is (EN)

A hardware-oriented implementation of **total arithmetic** and **"wiring = computation"**, spanning two of the four "heights" of the wider project (see *Related repositories*): the **Python gate simulation** and the **SystemVerilog (HDL) → FPGA** layers, on top of a shared integer-exact algebra core.

- **Total arithmetic** — the unit never produces `NaN` or `Inf`, and its status flags never lie. Overflow saturates to `±MAX` (flag `GE`), underflow collapses to `±MIN = ε` while preserving direction (flag `LE`), `a/0 = 0` (for a genuine zero only), and there are two kinds of zero. Adversarial sweeps produce zero `NaN`/`Inf` and zero flag lies.
- **Wiring = computation** — the multiplication is described by a *structure tensor* / wiring table. Swap the table and the same gate graph computes a different algebra: complex, quaternion, sedenion, matrix product, cyclic convolution. A registry of 19 patterns plus an inverse designer (ask for a matrix block → get the minimal group that contains it).
- **Signed-digit block floating point** — mantissa digits are ternary `{−1, 0, +1}` (Avizienis 1961) sharing one base-2 exponent per 16-component sedenion (a signed-ternary cousin of MXFP/microscaling). Sign symmetry means **sign flip = swapping `+`/`−` digits = pure wiring** (zero gates); the sign is carried by the mantissa, so only a "sign-unknown" bit remains as separate status.
- **From gates to HDL with no hand-transcription** — the SystemVerilog is *auto-emitted* by tracing the audited Python gate graph, so the HDL is the same gate graph by construction. Verified against the Python golden model with Icarus + cocotb.

### Architecture

```
Layer 1  algebra core (integer-exact, Python)     sd2_core, nd_algebra, matrix_algebra, bfp_sed
Layer 2  gate implementation (from AND/OR/NOT/XOR) gate_bilinear, gate_exponent, gate_fast, multi_add, mul_fused, ...
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
- **FPGA** — requires Vivado. `cd rtl/fpga && vivado -mode batch -source build.tcl` targets an Arty A7-100T (UART ⇔ `sd_add2`). Simulated at the same protocol as the board; synthesis is left to the user.

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
- **配線＝計算** — 乗算を*構造テンソル*（配線表）で記述。表を差し替えると同じゲートグラフが別の代数（複素・四元数・セデニオン・行列積・巡回畳み込み）になる。19パターンのレジストリ＋逆設計（欲しい行列ブロック→それを含む最小の群）。
- **符号つき3値ブロック浮動** — 仮数の桁が3値 `{−1,0,+1}`（Avizienis 1961）で、16成分セデニオンが指数を1つ共有（MXFP/microscalingの符号つき3値版）。符号対称なので **符号反転＝`+`/`−`桁の入替＝純配線**（ゲートゼロ）。符号は仮数が運び、残る状態は「符号不明」1ビットだけ。
- **ゲート→HDLを手写しゼロで** — SystemVerilogは監査済みPythonゲートグラフをトレースして*自動生成*。だからHDLは定義上同一のゲートグラフ。Icarus+cocotbでPython golden と照合。

### 再現方法

上記 *Reproduce* のコマンド。純Python/ゲートのself-testは特別なハード不要（ここでは全緑・違反0を実測）。HDLは `iverilog`+`cocotb`、FPGAは Vivado。

### ファイル地図

- **Layer 1 代数コア**: `sd2_core.py` `sd2_gates.py` `nd_algebra.py` `matrix_algebra.py` `bfp_sed.py`
- **Layer 2 ゲート**: `gate_bilinear.py` `gate_exponent.py` `gate_bfp.py` `gate_bfp2.py` `gate_fast.py` `multi_add.py` `mul_fused.py` `newton_recip.py`
- **Layer 3 配線**: `wiring_registry.py` `wiring_designer.py` `wiring_zoo.py` `wiring_patterns.py` `representation_lens.py` `bilinear_unit.py` `parallel_array.py`
- **Layer 4 HDL**: `rtl/`（自動生成SV・cocotb TB・Arty A7一式）
- **検証**: `interval_*.py` `numerical_test*.py` `stress_test*.py` `div_*.py`、`sed/`（セデニオン零因子ロジック）
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
