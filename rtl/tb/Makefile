# ⚠️ 生成AI使用・要検証
# cocotb + Verilator で SV 実装を Python golden と 突き合わせる。
#   make TOPLEVEL=gate9      … 全16通り
#   make TOPLEVEL=compress3  … 全64通り + 値保存
#   make TOPLEVEL=sd_add2    … ランダム506 + 端（既定）
SIM ?= verilator
TOPLEVEL_LANG ?= verilog
TOPLEVEL ?= sd_add2
MODULE = test_sd_primitives

GENERATED := sd_mult10 pe24 barrel18 blocknorm sed_comp
ifneq ($(filter $(TOPLEVEL),$(GENERATED)),)
VERILOG_SOURCES = $(PWD)/../generated/$(TOPLEVEL).sv
MODULE = test_generated
else ifeq ($(TOPLEVEL),top_arty)
VERILOG_SOURCES = $(PWD)/../sd_primitives.sv $(PWD)/../fpga/uart.sv $(PWD)/../fpga/top_arty.sv
MODULE = test_top_arty
COMPILE_ARGS += -Ptop_arty.DIV=20
else
VERILOG_SOURCES = $(PWD)/../sd_primitives.sv
endif
ifeq ($(SIM),verilator)
EXTRA_ARGS += --timing -Wno-UNOPTFLAT
endif

export COCOTB_TOPLEVEL := $(TOPLEVEL)

include $(shell cocotb-config --makefiles)/Makefile.sim
