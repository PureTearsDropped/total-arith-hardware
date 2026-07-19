# ⚠️ 生成AI使用・要検証 — Vivado バッチ合成（プロジェクトレス）
#   vivado -mode batch -source build.tcl
#   出力: build/top_arty.bit ＋ 使用率/タイミング レポート
set part xc7a100tcsg324-1
set outdir ./build
file mkdir $outdir

read_verilog -sv ../sd_primitives.sv
read_verilog -sv uart.sv
read_verilog -sv top_arty.sv
read_xdc arty_a7.xdc

synth_design -top top_arty -part $part
report_utilization -file $outdir/util_synth.rpt

opt_design
place_design
route_design

report_utilization -file $outdir/util.rpt
report_timing_summary -file $outdir/timing.rpt

write_bitstream -force $outdir/top_arty.bit
puts "==== DONE: $outdir/top_arty.bit ===="
