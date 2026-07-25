// ⚠️ 生成AI使用・要検証
// top_bench — オンチップ・ベンチハーネス: スループット/レイテンシを「推定」から「実測」へ。
//   LFSR(xorshift64) が 毎サイクル 乱数ベクタを 生成 → sd_add2(N=16) に 1加算/サイクル →
//   独立検算器 (SD桁→2進値の キャリー伝播加算器 = 被試験回路と 別方式) が 全数照合 →
//   UART は 集計 (不一致数・サイクル数・LFSR終状態) の 13 バイトだけ 運ぶ。
// プロトコル (115200 8N1): 受信 0xB7 → M ベクタ 全速実行 → 送信 0xB8 + mism[4B] + cycles[4B] + lfsr_lo[4B] (LE)
`default_nettype none

module bench_core #(parameter int N = 16, parameter bit CHECK = 1'b1,
                    parameter logic [63:0] SEED = 64'h9E3779B97F4A7C15) (
    input  wire clk, rst, start,
    input  wire [31:0] M,
    output logic running, done,
    output logic [31:0] mism, cycles, lfsr_lo
);
    function automatic [63:0] xs64(input [63:0] v);
        reg [63:0] t;
        begin
            t = v ^ (v << 13);
            t = t ^ (t >> 7);
            t = t ^ (t << 17);
            xs64 = t;
        end
    endfunction

    logic [63:0] s;
    logic [31:0] issued, retired;
    logic        v1, v2;
    // 検算の パイプ退去: 検算は「観測」であって「依存」でない — 結果は 検算を 待たない。
    // 段1で 生レールだけ 受け、差分計算(検算器)は 段2へ → 加算器の 最長経路から 検算が 消える。
    logic [N-1:0] xP1, xN1, yP1, yN1;
    logic [N:0]   zP1, zN1;
    logic signed [17:0] sxy2, sz2;

    // 正準化: (1,1)=冗長ゼロを 排除し 桁∈{−1,0,+1} に (host_test の to_rails と 同じ 規約)
    wire [N-1:0] xP = s[15:0]  & ~s[31:16], xN = s[31:16] & ~s[15:0];
    wire [N-1:0] yP = s[47:32] & ~s[63:48], yN = s[63:48] & ~s[47:32];
    wire [N:0] zP, zN;
    sd_add2 #(.N(N)) u_add (.xP(xP), .xN(xN), .yP(yP), .yN(yN), .zP(zP), .zN(zN));

    always_ff @(posedge clk) begin
        if (rst) begin
            running <= 1'b0; done <= 1'b0; v1 <= 1'b0; v2 <= 1'b0; mism <= '0; cycles <= '0;
        end else if (start && !running) begin
            running <= 1'b1; done <= 1'b0; s <= SEED;
            issued <= '0; retired <= '0; mism <= '0; cycles <= '0; v1 <= 1'b0; v2 <= 1'b0;
        end else if (running) begin
            cycles <= cycles + 1;
            if (issued < M) begin                   // 段0→1: 発行 (生レールを 受けるだけ)
                xP1 <= xP; xN1 <= xN; yP1 <= yP; yN1 <= yN; zP1 <= zP; zN1 <= zN;
                v1 <= 1'b1; s <= xs64(s); issued <= issued + 1;
            end else v1 <= 1'b0;
            v2 <= v1;
            if (v1 && CHECK) begin                  // 段1→2: 検算器の 差分 (経路外)
                sxy2 <= ($signed({2'b00, xP1}) - $signed({2'b00, xN1}))
                      + ($signed({2'b00, yP1}) - $signed({2'b00, yN1}));
                sz2  <= $signed({1'b0, zP1}) - $signed({1'b0, zN1});
            end
            if (v2) begin                           // 段2→3: 照合と 集計
                if (CHECK && (sxy2 != sz2)) mism <= mism + 1;
                retired <= retired + 1;
                if (retired + 32'd1 == M) begin
                    running <= 1'b0; done <= 1'b1; lfsr_lo <= s[31:0];
                end
            end
        end
    end
endmodule

module top_bench #(parameter int DIV = 868, parameter bit CHECK = 1'b1,
                   parameter logic [31:0] M = 32'd10_000_000) (
    input  wire CLK100MHZ,
    input  wire uart_txd_in,
    output wire uart_rxd_out,
    output wire [3:0] led
);
    wire clk = CLK100MHZ;
    logic rst = 1'b1;
    logic [3:0] rstcnt = '0;
    always_ff @(posedge clk)
        if (rstcnt != 4'hF) begin rstcnt <= rstcnt + 1; rst <= 1'b1; end
        else rst <= 1'b0;

    logic [7:0] rx_data;  logic rx_valid;
    logic [7:0] tx_data;  logic tx_start, tx_busy;
    uart_rx #(.DIV(DIV)) u_rx (.clk(clk), .rst(rst), .rxd(uart_txd_in), .data(rx_data), .valid(rx_valid));
    uart_tx #(.DIV(DIV)) u_tx (.clk(clk), .rst(rst), .data(tx_data), .start(tx_start),
                  .txd(uart_rxd_out), .busy(tx_busy));

    logic start, running, done;
    logic [31:0] mism, cycles, lfsr_lo;
    bench_core #(.N(16), .CHECK(CHECK)) u_core (.clk(clk), .rst(rst), .start(start), .M(M),
                                 .running(running), .done(done),
                                 .mism(mism), .cycles(cycles), .lfsr_lo(lfsr_lo));

    typedef enum logic [1:0] {IDLE, RUN, SEND} st_t;
    st_t st;
    logic [3:0] idx;
    logic [103:0] resp;                            // 13 バイト (LSB から 送る)

    always_ff @(posedge clk) begin
        tx_start <= 1'b0; start <= 1'b0;
        if (rst) begin st <= IDLE; idx <= '0; end
        else case (st)
            IDLE: if (rx_valid && rx_data == 8'hB7) begin start <= 1'b1; st <= RUN; end
            RUN:  if (done) begin
                resp <= {lfsr_lo, cycles, mism, 8'hB8};
                idx <= '0; st <= SEND;
            end
            SEND: if (!tx_busy && !tx_start) begin
                tx_data <= resp[idx*8 +: 8];
                tx_start <= 1'b1;
                if (idx == 4'd12) st <= IDLE; else idx <= idx + 1;
            end
            default: st <= IDLE;
        endcase
    end

    logic [26:0] beat;
    always_ff @(posedge clk) beat <= beat + 1;
    assign led[0] = beat[26];
    assign led[1] = running;
    assign led[2] = done;
    assign led[3] = (mism != 0);
endmodule

`default_nettype wire
