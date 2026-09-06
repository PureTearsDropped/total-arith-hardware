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

## Basis convention and independent cross-check / 基底規約と外部検証 (2026-09-06)

**EN.** The multiplication tables here are the Cayley–Dickson algebras in the *XOR labelling*: `e_i · e_j = OMEGA[i,j] · e_{i XOR j}`,
with `OMEGA` derived from the doubling `(a,b)(c,d) = (ac − d̄b, da + bc̄)`. `sed/crosscheck_external.py` checks them against
things this repository did not write: an exact rational Cayley–Dickson product typed from the textbook definition (256/256 basis
products, 1000/1000 random sedenion products), numpy-quaternion and Quaternions.jl (quaternions, 500/500 and 300/300), and
convention-free invariants (ℍ associative; 𝕆 alternative, Moufang, norm-multiplicative, non-associative; sedenions flexible,
power-associative, norm multiplicativity broken, and exactly 336 = 84×4 zero-divisor pairs `(e_a ± e_b)(e_c ± e_d) = 0`).
**Octonions.jl labels the seven imaginary units differently** (its products are not `e_{i XOR j}`), so element-wise numbers are *not*
interchangeable with it; the two tables are the same algebra, related by 1344 signed permutations (= |Aut(Fano)| × 8), e.g.
`e_i → s_i f_{p(i)}` with `p = (0,1,2,3,4,7,6,5)`, `s = (+,−,+,−,−,−,−,−)`. Run with `JULIA=/path/to/julia` (Quaternions.jl,
Octonions.jl installed) to include the Julia references; without it the exact reference and invariants still run.

**JP.** ここの乗算表は Cayley–Dickson 代数の **XOR ラベル規約** `e_i · e_j = OMEGA[i,j] · e_{i XOR j}` で、`OMEGA` は倍加公式
`(a,b)(c,d) = (ac − d̄b, da + bc̄)` から導出している。`sed/crosscheck_external.py` は、このリポジトリが書いていないもの
（教科書の定義から別に書いた厳密有理数の CD 積、numpy-quaternion、Quaternions.jl、規約に依存しない不変量）と突き合わせる
テストで、上記の結果は全部通過。**Octonions.jl は虚数単位の番号付けが違う**（積が `e_{i XOR j}` に落ちない）ので、要素ごとの
数値はそのままでは互換でない。代数としては同一で、符号付き置換 1344 通り（Fano 平面の自己同型 168 × 符号 8）のどれかで写る。
外部ライブラリと数値を突き合わせるときはこの基底変換を挟むこと。

## Branch `binary`: the {0,1} carry-save counterpart / 2 値 carry-save 版 (2026-09-06)

**EN.** `bin2_gates.py` is the two's-complement twin of the signed-digit rails, written in the same gate style
(works on ints for the golden model, on `gate_fast.B` for depth, on `rtl/emit_sv.T` for SystemVerilog): Baugh–Wooley
partial products → Dadda reduction to **two rows** (carry-save, no carry propagation) → `add_cs` / `mac_cs` keep
accumulating in that form → one Kogge–Stone `resolve` at the boundary. Emitted and checked against the golden with
cocotb (`bin_mult11`, `bin_mult11_ks`, `bin_mac11`: `./verify_hdl.sh`). Measured on sky130 with the emitted gates
mapped 1:1 to 1× cells (cell-level STA, constants folded), same numeric range (±1024 vs the 10-digit ±1023):

| module | output | cells | area µm² | delay |
|---|---|--:|--:|--:|
| `bin_mult11` | carry-save, two rows | 563 | 3 895 | 2.6 ns |
| `bin_mult11_ks` | two's complement | 847 | 5 767 | 4.3 ns |
| `bin_mac11` (x·y + carry-save accumulator) | carry-save, two rows | 792 | 5 540 | 3.5 ns |
| `sd_mult10` (signed digits, canonical out) | canonical (p,n) digits | 4 742 | 27 659 | 7.3 ns |

Carry-free accumulation is not a property of signed digits; binary carry-save has it at one seventh of the area.
What the signed-digit rails buy is representational (negation by wiring, symmetric digits, the sign-unknown /
structural-zero flags carried in the digits). This branch is the measurement, not a decision to switch.

**JP.** `bin2_gates.py` は signed-digit レールの 2 の補数版で、同じゲート流儀（golden は int、深さは `gate_fast.B`、SV は
`rtl/emit_sv.T` で追跡）。Baugh–Wooley の部分積 → Dadda で **2 行**（carry-save、桁上げ伝播なし）→ `add_cs` / `mac_cs` は
その形のまま累積 → 境界で 1 回だけ Kogge–Stone の `resolve`。生成した SV は cocotb で golden と照合済み。上の表は同じ数の
範囲での sky130 実測（セル遅延、定数畳み込み後）。桁上げ無しの累積は signed-digit の専売ではなく、2 値の carry-save でも
面積 7 分の 1 で得られる。signed-digit が買っているのは表現（配線だけの符号反転・対称な桁・符号不明や構造的零の旗を桁で運ぶ）。
このブランチは測定であって、切り替えの決定ではない。

### Same results as the signed-digit datapath? / SD 版と同じ結果を出すか (2026-09-06)

**EN.** `bin_vs_sd_equiv.py` feeds the same integers to both goldens: scalar products agree on all 225 edge
combinations and 20 000 random pairs (signed-digit golden = signed-digit gate graph = binary rows = binary resolved
= a·b), multiply-accumulate chains of 2–16 terms agree, and the full 16-component sedenion product (OMEGA signs
folded with `neg_cs`) agrees with `sedenion_mult_sd2` and `ref_mult` on 200 random pairs.
*Flags.* The datapath flags are value-level: products never raise sign-unknown, exact inputs get flags only from
normalisation (truncation → ≥, overflow → ±MAX ≥, collapse → ±MIN ≤), flagged inputs propagate intervals — none of
that depends on the digit encoding. What does depend on it is *reading* the value: canonical signed digits show
"true zero" and the sign in the digits, carry-save rows do not (r0 + r1 ≡ 0 with both rows nonzero). The branch adds
carry-free primitives for exactly that — `is_zero_cs` (Cortadella–Llabería local test + AND tree) and `sign_cs`
(prefix carry into the MSB only), both exact (W = 4 exhaustive, W = 22 random); zero-divisor products report
structural zero on every component without resolving. Not on this branch: a binary version of the block normaliser /
BFP unit (`blocknorm`, `gate_bfp.py`) — that is the next step if the switch is made.

**JP.** `bin_vs_sd_equiv.py` で同じ整数を両方の golden に流した。スカラー積は端の全組合せ 225 と乱数 20,000 で
SD golden = SD ゲート = 2 値の 2 行 = 2 値の解決済み = a·b、2〜16 項の積和連鎖も一致、セデニオン積 16 成分（OMEGA の
符号は `neg_cs` で畳む）も `sedenion_mult_sd2` と `ref_mult` に 200 組で一致。**旗**は値の上の論理（積は符号不明を
立てない・厳密入力の旗は正規化のみ・旗付き入力は区間伝播）なので桁の符号化に依存しない。依存するのは「値を読む」
所で、正準 SD なら桁を見れば真の 0 と符号が分かるが carry-save の 2 行では分からない。そのための桁上げ無しの原始回路
`is_zero_cs`（局所条件+AND 木）と `sign_cs`（最上位への桁上げだけをプレフィックスで）を追加し厳密に検証、零因子の積は
全成分で解決せずに構造的零を報告する。未着手: 2 値版のブロック正規化器/BFP ユニット（切り替えるならそこが次）。

### Binary block normaliser / 2 値ブロック正規化器 (2026-09-06)

**EN.** `bin2_bfp.py` is the binary counterpart of `block_normalize_g_fast`: M carry-save inputs (two rows each) → sign +
W-bit magnitude per component, shared exponent, the same `ge`/`le` flags by the same rules (truncation → ≥, exponent
overflow → ±MAX ≥, collapse → ±MIN ≤, true zero unflagged). Magnitude is obtained by resolving `r0+r1` and, in parallel,
`−(r0+r1)` (rows complemented + 2), selected by the sign — the bidirectional trick of `canonicalize_fast` on rows — so the
truncation keeps the magnitude semantics of the flags (an arithmetic right shift of a negative two's-complement number
would break "dropped bits ⇒ ≥"). Golden-vs-golden: 600 random blocks × 4 components with random carry-save splits,
values, exponent and flags identical to the signed-digit normaliser; emitted `bin_blocknorm.sv` passes the cocotb test
against the signed-digit golden (120 blocks). sky130, 1× cells, constants folded: `blocknorm` (signed digits) 9,643 cells /
57,519 µm² / 21.5 ns, `bin_blocknorm` 7,587 cells / 46,168 µm² / 20.7 ns. With this the binary branch covers multiply,
multiply-accumulate, resolve, zero/sign detection and normalisation; not yet: the fused sedenion component unit
(`sed_comp`) and the exponent/ε machinery of `gate_bfp.py` on binary inputs.

**JP.** `bin2_bfp.py` は `block_normalize_g_fast` の 2 値版。carry-save の入力 M 個 → 成分ごとに符号+W ビットの大きさ、
共有指数、同じ規則の `ge`/`le` 旗。大きさは `r0+r1` と `−(r0+r1)`（行の補数+2）を並列に解いて符号で選ぶ（`canonicalize_fast`
の双方向の手を行に適用）ので、切り捨ては大きさの切り捨てになり旗の意味（落とした ⇒ ≥）が保たれる（2 の補数の算術右
シフトだと負数で破れる）。golden 同士: 乱数 600 ブロック × 4 成分（carry-save の分割も乱数）で値・指数・旗が SD 版と完全一致。
生成した `bin_blocknorm.sv` は cocotb で SD golden と一致（120 ブロック）。sky130: SD 版 9,643 セル / 57,519 µm² / 21.5 ns、
2 値版 7,587 / 46,168 / 20.7 ns。これで binary ブランチは乗算・積和・解決・零/符号判定・正規化まで。未着手は融合セデニオン
成分ユニット（`sed_comp`）と `gate_bfp.py` の指数/ε 機構の 2 値化。

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
