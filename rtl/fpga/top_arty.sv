// ⚠️ 生成AI使用・要検証
// Arty A7-100T トップ — ハードウェア版 cocotb: PC が UART で テストベクタを 流し、
// 実 FPGA 上の sd_add2(N=16) が 答えを 返す（hardware-in-the-loop）。
//
// プロトコル（115200 8N1）:
//   受信: 0xA5(同期) + xP[2B] xN[2B] yP[2B] yN[2B]（リトルエンディアン・計 9 バイト）
//   送信: zP[3B] zN[3B]（17 ビットを 3 バイトに 右詰め・計 6 バイト）
// LED: ld[0]=心拍  ld[1]=フレーム受信中  ld[2]=応答送信中  ld[3]=最後の zP の LSB
`default_nettype none

module top_arty #(parameter int DIV = 868) (   // 115200 @100MHz。シムでは 小さく 上書き
    input  wire CLK100MHZ,
    input  wire uart_txd_in,          // ホスト → FPGA
    output wire uart_rxd_out,         // FPGA → ホスト
    output wire [3:0] led
);
    localparam int N = 16;
    wire clk = CLK100MHZ;
    logic rst = 1'b1;                 // 起動後 数クロックで 解除
    logic [3:0] rstcnt = '0;
    always_ff @(posedge clk)
        if (rstcnt != 4'hF) begin rstcnt <= rstcnt + 1; rst <= 1'b1; end
        else rst <= 1'b0;

    // ---------------- UART
    logic [7:0] rx_data;  logic rx_valid;
    logic [7:0] tx_data;  logic tx_start, tx_busy;
    uart_rx #(.DIV(DIV)) u_rx (.clk(clk), .rst(rst), .rxd(uart_txd_in), .data(rx_data), .valid(rx_valid));
    uart_tx #(.DIV(DIV)) u_tx (.clk(clk), .rst(rst), .data(tx_data), .start(tx_start),
                  .txd(uart_rxd_out), .busy(tx_busy));

    // ---------------- 被試験回路: sd_add2（組合せ・監査/シム/HDL 検証済み）
    logic [N-1:0] xP, xN, yP, yN;
    wire  [N:0]   zP, zN;
    sd_add2 #(.N(N)) dut (.xP(xP), .xN(xN), .yP(yP), .yN(yN), .zP(zP), .zN(zN));

    // ---------------- フレーム FSM
    typedef enum logic [1:0] {W_SYNC, W_PAY, SEND} st_t;
    st_t st;
    logic [3:0] idx;
    logic [23:0] outP, outN;

    always_ff @(posedge clk) begin
        tx_start <= 1'b0;
        if (rst) begin st <= W_SYNC; idx <= '0; end
        else case (st)
            W_SYNC: if (rx_valid && rx_data == 8'hA5) begin st <= W_PAY; idx <= '0; end
            W_PAY: if (rx_valid) begin
                case (idx)
                    4'd0: xP[7:0]   <= rx_data;  4'd1: xP[15:8] <= rx_data;
                    4'd2: xN[7:0]   <= rx_data;  4'd3: xN[15:8] <= rx_data;
                    4'd4: yP[7:0]   <= rx_data;  4'd5: yP[15:8] <= rx_data;
                    4'd6: yN[7:0]   <= rx_data;  4'd7: yN[15:8] <= rx_data;
                    default: ;
                endcase
                if (idx == 4'd7) begin st <= SEND; idx <= '0; end
                else idx <= idx + 1;
            end
            SEND: begin
                // 1 サイクル目に 結果を ラッチ（組合せ出力は 前サイクルで 確定済み）
                if (idx == 4'd0) begin
                    outP <= {7'b0, zP}; outN <= {7'b0, zN};
                    idx <= 4'd1;
                end else if (!tx_busy && !tx_start) begin
                    case (idx)
                        4'd1: tx_data <= outP[7:0];    4'd2: tx_data <= outP[15:8];
                        4'd3: tx_data <= outP[23:16];  4'd4: tx_data <= outN[7:0];
                        4'd5: tx_data <= outN[15:8];   4'd6: tx_data <= outN[23:16];
                        default: ;
                    endcase
                    tx_start <= 1'b1;
                    if (idx == 4'd6) st <= W_SYNC; else idx <= idx + 1;
                end
            end
            default: st <= W_SYNC;
        endcase
    end

    // ---------------- LED
    logic [26:0] beat;
    always_ff @(posedge clk) beat <= beat + 1;
    assign led[0] = beat[26];                     // 心拍 ~0.75Hz
    assign led[1] = (st == W_PAY);
    assign led[2] = (st == SEND);
    assign led[3] = outP[0];
endmodule

`default_nettype wire
