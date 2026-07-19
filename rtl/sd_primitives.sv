// ⚠️ 生成AI使用・要検証
// 符号つき桁 (p,n) の 原始部品 — Python 実装（監査3回合格）の 構造をそのまま SV に。
// 桁 t ∈ {-1,0,+1} は 2 線 (p,n): t = p - n。(1,1) は 冗長ゼロ。
// ゲートは AND/OR/NOT/XOR 相当の 連続代入のみ（振る舞い記述の +,-,if は 使わない）。

`default_nettype none

// ---------------- gate9: 桁 1×1 の 積（6 ゲート） ----------------
module gate9 (
    input  wire xp, xn, yp, yn,
    output wire rp, rn
);
    // 同符号 → +、異符号 → −（gate_bilinear.gate9 と 同一構造）
    assign rp = (xp & yp) | (xn & yn);
    assign rn = (xp & yn) | (xn & yp);
endmodule

// ---------------- full_adder: ビット全加算器（5 ゲート） ----------------
module fa (
    input  wire a, b, c,
    output wire s, cy
);
    wire ab;
    assign ab = a ^ b;
    assign s  = ab ^ c;
    assign cy = (a & b) | (c & ab);
endmodule

// ---------------- compress3: 3:2 圧縮器（18 ゲート・値保存 low+2*high = x+y+c） ----------------
module compress3 (
    input  wire xp, xn, yp, yn, cp, cn,
    output wire lp, ln, hp, hn
);
    wire ps0, pc, ns0, nc;
    fa fap (.a(xp), .b(yp), .c(cp), .s(ps0), .cy(pc));   // 正レール
    fa fan (.a(xn), .b(yn), .c(cn), .s(ns0), .cy(nc));   // 負レール
    assign lp = ps0 & ~ns0;                              // low  = Ps0 − Ns0
    assign ln = ns0 & ~ps0;
    assign hp = pc  & ~nc;                               // high = Pc − Nc
    assign hn = nc  & ~pc;
endmodule

// ---------------- sd_add2: 定数深さ SD 加算器（Avizienis 限定キャリー） ----------------
//  入力: N 桁 × 2 数（P/N レールの パック配列）。出力: N+1 桁。
//  gate_fast.sd_add2 と 同一構造: compress3 で (l,h) → pos/neg → 転送 t・仮 w → z=w+t。
module sd_add2 #(parameter int N = 16) (
    input  wire [N-1:0] xP, xN, yP, yN,
    output wire [N:0]   zP, zN
);
    wire [N-1:0] lP, lN, hP, hN;        // s_i = l + 2h
    wire [N-1:0] pos, neg;              // s_i ≥ 1 / s_i ≤ −1
    wire [N:0]   tP, tN;                // 転送（位置 i+1 へ）
    wire [N-1:0] wP, wN;                // 仮

    genvar i;
    generate
        for (i = 0; i < N; i++) begin : g_s
            compress3 c3 (.xp(xP[i]), .xn(xN[i]), .yp(yP[i]), .yn(yN[i]),
                          .cp(1'b0), .cn(1'b0),
                          .lp(lP[i]), .ln(lN[i]), .hp(hP[i]), .hn(hN[i]));
            assign pos[i] = hP[i] | lP[i];
            assign neg[i] = hN[i] | lN[i];
        end
    endgenerate

    assign tP[0] = 1'b0;
    assign tN[0] = 1'b0;
    generate
        for (i = 0; i < N; i++) begin : g_t
            // 隣（下位）の 符号。i=0 は 定数 0 に 接地（形状分岐 = 配線）
            wire npv = (i > 0) ? neg[i-1] : 1'b0;
            wire ppv = (i > 0) ? pos[i-1] : 1'b0;
            assign tP[i+1] = hP[i] | (lP[i] & ~npv);
            assign tN[i+1] = hN[i] | (lN[i] & ~ppv);
            assign wP[i]   = (lP[i] & npv) | (lN[i] & ~ppv);
            assign wN[i]   = (lP[i] & ~npv) | (lN[i] & ppv);
        end
    endgenerate

    generate
        for (i = 0; i <= N; i++) begin : g_z
            wire wp_i = (i < N) ? wP[i] : 1'b0;
            wire wn_i = (i < N) ? wN[i] : 1'b0;
            wire zhp, zhn;    // 構成上 0（伝播 終端）— テストで 確認
            compress3 cz (.xp(wp_i), .xn(wn_i), .yp(tP[i]), .yn(tN[i]),
                          .cp(1'b0), .cn(1'b0),
                          .lp(zP[i]), .ln(zN[i]), .hp(zhp), .hn(zhn));
        end
    endgenerate
endmodule

`default_nettype wire
