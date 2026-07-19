#!/usr/bin/env bash
# Gate + pure-Python self-tests (no special hardware). Reports PASS/FAIL.
set -u; PY="${PY:-python3}"; cd "$(dirname "$0")"
T="sd2_gates nd_algebra matrix_algebra bfp_sed gate_bilinear gate_bfp2 gate_fast multi_add mul_fused wiring_registry wiring_designer representation_lens bilinear_unit wiring_patterns newton_recip interval_adversarial stress_test"
pass=0; fail=0; tmp="$(mktemp)"
for t in $T; do
  if $PY -u "$t.py" >"$tmp" 2>&1; then echo "PASS  $t"; pass=$((pass+1)); else echo "FAIL  $t"; fail=$((fail+1)); tail -3 "$tmp"; fi
done
rm -f "$tmp"; echo "----"; echo "$pass passed, $fail failed"; [ "$fail" -eq 0 ]
