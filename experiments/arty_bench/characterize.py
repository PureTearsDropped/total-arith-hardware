#!/usr/bin/env python3
# characterize.py — 棚の 生成ユニットの 実シリコン特性表 (LUT/FF/Fmax) を 自動採取。
#   各ユニットを 汎用ラッパ (入力=シフトレジスタ・出力=XOR縮約→FF) に 包み、
#   yosys → nextpnr-xilinx (xc7a100t) を 回して 配置配線後の 数字を 拾う。
#   ラッパの 意味: 合成が 論理を 消せない (入力未知・出力観測) まま IO ピンを 3 本に 抑える。
import re, subprocess, sys, pathlib

RTL = "../../rtl"
UNITS = ["pe24", "barrel18", "sd_mult10", "blocknorm", "sed_comp"]

def ports_of(path):
    src = open(path).read()
    m = re.search(r"^module\s+(\w+)\s*\((.*?)\);", src, re.S | re.M)
    name, body = m.group(1), m.group(2)
    ins, outs = [], []
    for d, w, nm in re.findall(r"(input|output)\s+wire\s+(?:\[(\d+):0\]\s+)?(\w+)", body):
        width = int(w) + 1 if w else 1
        (ins if d == "input" else outs).append((nm, width))
    return name, ins, outs

def emit_wrapper(name, ins, outs):
    tin = sum(w for _, w in ins)
    conns, base = [], 0
    for nm, w in ins:
        conns.append(f".{nm}(sin_reg[{base + w - 1}:{base}])")
        base += w
    owires, oconns = [], []
    for nm, w in outs:
        owires.append(f"    wire [{w - 1}:0] o_{nm};")
        oconns.append(f".{nm}(o_{nm})")
    ored = ", ".join(f"o_{nm}" for nm, _ in outs)
    return f"""`default_nettype none
module chz_top (
    input  wire CLK100MHZ,
    input  wire uart_txd_in,
    output wire uart_rxd_out,
    output wire [3:0] led
);
    reg [{tin - 1}:0] sin_reg;
    always @(posedge CLK100MHZ) sin_reg <= {{sin_reg[{tin - 2}:0], uart_txd_in}};
{chr(10).join(owires)}
    {name} u ({", ".join(conns + oconns)});
    reg out_r;
    always @(posedge CLK100MHZ) out_r <= ^{{{ored}}};
    assign uart_rxd_out = out_r;
    assign led = 4'b0;
endmodule
`default_nettype wire
"""

def run(unit):
    path = f"{RTL}/generated/{unit}.sv"
    name, ins, outs = ports_of(path)
    open("chz_top.sv", "w").write(emit_wrapper(name, ins, outs))
    subprocess.run(
        f"/home/claude/yosys/yosys -q -p 'read_verilog -sv {path} chz_top.sv; "
        f"synth_xilinx -flatten -abc9 -arch xc7 -top chz_top; write_json chz.json'",
        shell=True, check=True)
    r = subprocess.run(
        f"/home/claude/nextpnr-xilinx/build/nextpnr-xilinx "
        f"--chipdb /home/claude/chipdb/xc7a100t.bin --xdc {RTL}/fpga/arty_nextpnr.xdc "
        f"--json chz.json --fasm /dev/null", shell=True, capture_output=True, text=True)
    log = r.stderr
    lut = re.findall(r"SLICE_LUTX:\s+(\d+)/", log)
    ff = re.findall(r"SLICE_FFX:\s+(\d+)/", log)
    fmax = re.findall(r"Max frequency for clock\s+'[^']+':\s+([\d.]+) MHz", log)
    nff = int(ff[-1]) if ff else 0
    tin = sum(w for _, w in ins)
    print(f"{unit:11s} LUT {int(lut[-1]) if lut else -1:>6d}  FF {nff:>5d} "
          f"(うちラッパ {tin + 1})  Fmax {float(fmax[-1]) if fmax else -1:6.1f} MHz "
          f"→ 遅延 {1000 / float(fmax[-1]) if fmax else -1:5.1f} ns", flush=True)

if __name__ == "__main__":
    print("ユニット特性表 (xc7a100t 実ターゲット・配置配線後):", flush=True)
    for u in (sys.argv[1:] or UNITS):
        try:
            run(u)
        except Exception as e:
            print(f"{u}: 失敗 {e}", flush=True)
