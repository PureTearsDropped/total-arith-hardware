// ⚠️ 生成AI使用・要検証 — emit_sv.py が 監査済み Python から 自動生成（手写しなし）
`default_nettype none
module sd_trit (
    input  wire [0:0] tP,
    input  wire [0:0] tN,
    input  wire [7:0] xP,
    input  wire [7:0] xN,
    output wire [7:0] zP,
    output wire [7:0] zN
);
    wire w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14, w15, w16, w17, w18, w19;
    wire w20, w21, w22, w23, w24, w25, w26, w27, w28, w29, w30, w31, w32, w33, w34, w35, w36, w37, w38, w39;
    wire w40, w41, w42, w43, w44, w45, w46, w47;
    assign w0 = tP[0] & xP[0];
    assign w1 = tN[0] & xN[0];
    assign w2 = w0 | w1;
    assign w3 = tP[0] & xN[0];
    assign w4 = tN[0] & xP[0];
    assign w5 = w3 | w4;
    assign w6 = tP[0] & xP[1];
    assign w7 = tN[0] & xN[1];
    assign w8 = w6 | w7;
    assign w9 = tP[0] & xN[1];
    assign w10 = tN[0] & xP[1];
    assign w11 = w9 | w10;
    assign w12 = tP[0] & xP[2];
    assign w13 = tN[0] & xN[2];
    assign w14 = w12 | w13;
    assign w15 = tP[0] & xN[2];
    assign w16 = tN[0] & xP[2];
    assign w17 = w15 | w16;
    assign w18 = tP[0] & xP[3];
    assign w19 = tN[0] & xN[3];
    assign w20 = w18 | w19;
    assign w21 = tP[0] & xN[3];
    assign w22 = tN[0] & xP[3];
    assign w23 = w21 | w22;
    assign w24 = tP[0] & xP[4];
    assign w25 = tN[0] & xN[4];
    assign w26 = w24 | w25;
    assign w27 = tP[0] & xN[4];
    assign w28 = tN[0] & xP[4];
    assign w29 = w27 | w28;
    assign w30 = tP[0] & xP[5];
    assign w31 = tN[0] & xN[5];
    assign w32 = w30 | w31;
    assign w33 = tP[0] & xN[5];
    assign w34 = tN[0] & xP[5];
    assign w35 = w33 | w34;
    assign w36 = tP[0] & xP[6];
    assign w37 = tN[0] & xN[6];
    assign w38 = w36 | w37;
    assign w39 = tP[0] & xN[6];
    assign w40 = tN[0] & xP[6];
    assign w41 = w39 | w40;
    assign w42 = tP[0] & xP[7];
    assign w43 = tN[0] & xN[7];
    assign w44 = w42 | w43;
    assign w45 = tP[0] & xN[7];
    assign w46 = tN[0] & xP[7];
    assign w47 = w45 | w46;
    assign zP[0] = w2;
    assign zP[1] = w8;
    assign zP[2] = w14;
    assign zP[3] = w20;
    assign zP[4] = w26;
    assign zP[5] = w32;
    assign zP[6] = w38;
    assign zP[7] = w44;
    assign zN[0] = w5;
    assign zN[1] = w11;
    assign zN[2] = w17;
    assign zN[3] = w23;
    assign zN[4] = w29;
    assign zN[5] = w35;
    assign zN[6] = w41;
    assign zN[7] = w47;
endmodule
`default_nettype wire
