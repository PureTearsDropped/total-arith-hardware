#!/usr/bin/env bash
# 誠実さ税 A/B ビルド: CHECK=0(裸) / CHECK=1(全数検算並走) を 同一フローで
# ツール位置は 環境変数で 差し替え可 (既定は この計画の 開発機の 配置)
set -eu
RTL="${RTL:-../../rtl}"
YOSYS="${YOSYS:-/home/claude/yosys/yosys}"
NEXTPNR="${NEXTPNR:-/home/claude/nextpnr-xilinx/build/nextpnr-xilinx}"
CHIPDB="${CHIPDB:-/home/claude/chipdb/xc7a100t.bin}"
PRJXRAY="${PRJXRAY:-/home/claude/prjxray}"
PRJXRAY_DB="${PRJXRAY_DB:-/home/claude/prjxray-db}"
PYF2F="${PYF2F:-python3}"                 # textx + prjxray が import できる python
export PYTHONPATH="${PYTHONPATH:-}:/home/claude/openxc7/lib/python:$PRJXRAY"
for CK in 0 1; do
  $YOSYS -q -p "read_verilog -sv $RTL/sd_primitives.sv $RTL/fpga/uart.sv top_bench.sv; \
    chparam -set CHECK $CK top_bench; synth_xilinx -flatten -abc9 -arch xc7 -top top_bench; \
    write_json ab_$CK.json"
  $NEXTPNR --chipdb $CHIPDB --xdc $RTL/fpga/arty_nextpnr.xdc \
    --json ab_$CK.json --fasm ab_$CK.fasm 2> results/ab_$CK.log
  echo "== CHECK=$CK:"; grep -E "SLICE_LUTX:|Max frequency" results/ab_$CK.log | tail -2
  $PYF2F $PRJXRAY/utils/fasm2frames.py --part xc7a100tcsg324-1 \
    --db-root $PRJXRAY_DB/artix7 ab_$CK.fasm > ab_$CK.frames
  $PRJXRAY/build/tools/xc7frames2bit \
    --part_file $PRJXRAY_DB/artix7/xc7a100tcsg324-1/part.yaml \
    --part_name xc7a100tcsg324-1 --frm_file ab_$CK.frames \
    --output_file top_bench_$([ $CK = 0 ] && echo bare || echo checked).bit
done
