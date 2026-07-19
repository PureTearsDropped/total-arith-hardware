#!/usr/bin/env python3
# ⚠️ 生成AI使用・要検証
"""SV エミッタ — 監査済み Python ゲート関数を トレースして **同一ゲートグラフの SV** を 自動生成。

  仕組み: ワイヤ型 T が &,|,^ を 記録（値は 持たない）。監査済み関数を T で 走らせると
  ネットリスト（assign 列）が 溜まる → module として 書き出す。
  ⟹ 手写しゼロ・「Python と SV が 同じ回路」が 構成で 保証される。

  生成対象（1: 乗算器・優先エンコーダ・バレル・block_normalize ／ 2: 融合MAC=群積成分）:
    sd_mult10      : 10桁 × 10桁 乗算器（multiply_fast）
    pe24           : 24桁 優先エンコーダ（priority_encoder_fast）
    barrel18       : 18桁 可変右シフト（barrel_shift_right_digits）
    blocknorm      : M=4, 24桁, W=6 ブロック正規化（block_normalize_g_fast・フラグ/ε/飽和込み）
    sed_comp       : セデニオン積の 成分 k=1（group_component・16積の 融合MAC・符号=配線）
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

NET = []          # (dst, op, a, b)  op ∈ '&','|','^'
NWIRE = [0]

class T:
    """トレース ワイヤ。演算を 記録するだけ（値なし）。定数は 文字列 "1'b0"/"1'b1"。"""
    __slots__ = ('name',)
    def __init__(self, name=None):
        if name is None:
            name = f"w{NWIRE[0]}"; NWIRE[0] += 1
        self.name = name
    @staticmethod
    def _c(o):
        if isinstance(o, T): return o
        return T._const(int(o))
    @staticmethod
    def _const(v):
        t = T.__new__(T); t.name = "1'b1" if v else "1'b0"; return t
    def _op(self, o, sym):
        o = T._c(o); r = T()
        NET.append((r.name, sym, self.name, o.name))
        return r
    def __and__(self, o):  return self._op(o, '&')
    def __rand__(self, o): return T._c(o)._op(self, '&')
    def __or__(self, o):   return self._op(o, '|')
    def __ror__(self, o):  return T._c(o)._op(self, '|')
    def __xor__(self, o):  return self._op(o, '^')
    def __rxor__(self, o): return T._c(o)._op(self, '^')

def reset():
    NET.clear(); NWIRE[0] = 0

def in_bus(name, width):
    """入力バス → ビットごとの T（名前 = name[i]）。"""
    return [T(f"{name}[{i}]") for i in range(width)]

def in_digits(pname, nname, width):
    P = in_bus(pname, width); N = in_bus(nname, width)
    return [(P[i], N[i]) for i in range(width)]

def emit_module(fname, modname, inputs, outputs):
    """inputs: [(port, width)] / outputs: [(port, width, wires)]（wires は T か 定数 T）。"""
    lines = [f"// ⚠️ 生成AI使用・要検証 — emit_sv.py が 監査済み Python から 自動生成（手写しなし）",
             f"`default_nettype none", f"module {modname} ("]
    ports = [f"    input  wire [{w-1}:0] {p}" for p, w in inputs]
    ports += [f"    output wire [{w-1}:0] {p}" for p, w, _ in outputs]
    lines.append(",\n".join(ports)); lines.append(");")
    n_gate = 0
    decl = [dst for dst, _, _, _ in NET]
    for i in range(0, len(decl), 20):
        lines.append("    wire " + ", ".join(decl[i:i+20]) + ";")
    for dst, op, a, b in NET:
        lines.append(f"    assign {dst} = {a} {op} {b};"); n_gate += 1
    for p, w, wires in outputs:
        for i, t in enumerate(wires):
            nm = t.name if isinstance(t, T) else ("1'b1" if t else "1'b0")
            lines.append(f"    assign {p}[{i}] = {nm};")
    lines.append("endmodule"); lines.append("`default_nettype wire")
    with open(fname, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"  {modname:<10} → {os.path.basename(fname):<16} ゲート {n_gate:>7,} 本（assign）")

def null_st():
    from collections import defaultdict
    return defaultdict(int)

def rails(digs):
    """桁リスト → (Pワイヤ列, Nワイヤ列)。定数 ZERO は 定数線に。"""
    P = []; N = []
    for d in digs:
        p, n = d
        P.append(p if isinstance(p, T) else T._const(int(p)))
        N.append(n if isinstance(n, T) else T._const(int(n)))
    return P, N


OUT = os.path.join(os.path.dirname(__file__), "generated")
os.makedirs(OUT, exist_ok=True)

def gen_sd_mult(K=10):
    reset()
    from gate_fast import multiply_fast
    X = in_digits("xP", "xN", K); Y = in_digits("yP", "yN", K)
    Z = multiply_fast(X, Y, null_st())
    P, N = rails(Z)
    emit_module(os.path.join(OUT, "sd_mult10.sv"), "sd_mult10",
                [("xP", K), ("xN", K), ("yP", K), ("yN", K)],
                [("zP", len(Z), P), ("zN", len(Z), N)])
    return len(Z)

def gen_pe(n=24, EW=6):
    reset()
    from gate_fast import priority_encoder_fast
    D = in_digits("dP", "dN", n)
    L, none, onehot = priority_encoder_fast(D, EW, null_st())
    Lw = [x if isinstance(x, T) else T._const(int(x)) for x in L]
    emit_module(os.path.join(OUT, "pe24.sv"), "pe24",
                [("dP", n), ("dN", n)],
                [("L", EW, Lw), ("none_o", 1, [none if isinstance(none, T) else T._const(int(none))])])

def gen_barrel(n=18, SW=5):
    reset()
    from gate_exponent import barrel_shift_right_digits
    D = in_digits("dP", "dN", n)
    S = in_bus("S", SW)
    out, dropped = barrel_shift_right_digits(D, S, null_st())
    P, N = rails(out)
    emit_module(os.path.join(OUT, "barrel18.sv"), "barrel18",
                [("dP", n), ("dN", n), ("S", SW)],
                [("oP", n, P), ("oN", n, N),
                 ("dropped", 1, [dropped if isinstance(dropped, T) else T._const(int(dropped))])])

def gen_blocknorm(M=4, Win=24, W=6, Emax=20, EW=12):
    reset()
    from gate_fast import block_normalize_g_fast
    mants = [in_digits(f"m{i}P", f"m{i}N", Win) for i in range(M)]
    Ebus = in_bus("Ein", EW)
    out, E_fin, flags = block_normalize_g_fast(mants, Ebus, W, Emax, null_st())
    inputs = [(f"m{i}{r}", Win) for i in range(M) for r in ("P", "N")] + [("Ein", EW)]
    outputs = []
    for i in range(M):
        P, N = rails(out[i])
        outputs += [(f"o{i}P", W, P), (f"o{i}N", W, N)]
        ge, le, _ = flags[i]
        outputs += [(f"flag{i}", 2, [ge if isinstance(ge, T) else T._const(int(ge)),
                                     le if isinstance(le, T) else T._const(int(le))])]
    Ew = [x if isinstance(x, T) else T._const(int(x)) for x in E_fin]
    outputs.append(("Eout", EW, Ew))
    emit_module(os.path.join(OUT, "blocknorm.sv"), "blocknorm", inputs, outputs)

def gen_sed_comp(M=16, K=6, k=1):
    reset()
    from mul_fused import group_component
    from nd_algebra import cd_omega
    OM = cd_omega(M)
    OMl = [[int(OM[i, j]) for j in range(M)] for i in range(M)]
    a = [in_digits(f"a{i}P", f"a{i}N", K) for i in range(M)]
    b = [in_digits(f"b{i}P", f"b{i}N", K) for i in range(M)]
    Z = group_component(a, b, OMl, M, k, null_st())
    P, N = rails(Z)
    inputs = [(f"a{i}{r}", K) for i in range(M) for r in ("P", "N")] + \
             [(f"b{i}{r}", K) for i in range(M) for r in ("P", "N")]
    emit_module(os.path.join(OUT, "sed_comp.sv"), "sed_comp", inputs,
                [("zP", len(Z), P), ("zN", len(Z), N)])
    return len(Z)


if __name__ == "__main__":
    print("SV 自動生成（監査済み Python → 同一ゲートグラフ）:")
    wz = gen_sd_mult()
    gen_pe()
    gen_barrel()
    gen_blocknorm()
    ws = gen_sed_comp()
    print(f"  （sd_mult10 出力幅 {wz}・sed_comp 出力幅 {ws}）")
