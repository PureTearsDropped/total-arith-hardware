#!/usr/bin/env bash
# Virtual hardware: run each SystemVerilog module in RTL simulation (iverilog)
# and check it against the Python golden via cocotb. No FPGA required.
set -u; cd "$(dirname "$0")/rtl/tb"
tmp="$(mktemp)"; pass=0; fail=0
for TL in gate9 compress3 sd_add2 sd_mult10 pe24 barrel18 blocknorm sed_comp; do
  rm -rf sim_build results.xml
  if make SIM="${SIM:-icarus}" TOPLEVEL="$TL" >"$tmp" 2>&1 && grep -q 'FAIL=0' "$tmp"; then
    echo "PASS  $TL"; pass=$((pass+1)); else echo "FAIL  $TL"; tail -4 "$tmp"; fail=$((fail+1)); fi
done
rm -f "$tmp"; echo "----"; echo "$pass passed, $fail failed"; [ "$fail" -eq 0 ]
