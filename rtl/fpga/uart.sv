// ⚠️ 生成AI使用・要検証
// UART 8N1（既定 115200 @ 100MHz ⟹ DIV=868）。受信は ビット中央 サンプル。
`default_nettype none

module uart_rx #(parameter int DIV = 868) (
    input  wire clk, rst,
    input  wire rxd,
    output logic [7:0] data,
    output logic       valid          // 1 クロック パルス
);
    logic [1:0] sync;
    always_ff @(posedge clk) sync <= {sync[0], rxd};
    wire rx = sync[1];

    typedef enum logic [1:0] {IDLE, START, BITS, STOP} st_t;
    st_t st;
    logic [$clog2(DIV)-1:0] cnt;
    logic [2:0] bitn;
    logic [7:0] sh;

    always_ff @(posedge clk) begin
        valid <= 1'b0;
        if (rst) begin st <= IDLE; cnt <= '0; bitn <= '0; end
        else case (st)
            IDLE:  if (!rx) begin st <= START; cnt <= '0; end
            START: if (cnt == DIV/2) begin                      // スタートビット 中央
                       if (!rx) begin st <= BITS; cnt <= '0; bitn <= '0; end
                       else st <= IDLE;                          // グリッチ
                   end else cnt <= cnt + 1;
            BITS:  if (cnt == DIV-1) begin
                       cnt <= '0; sh <= {rx, sh[7:1]};           // LSB first
                       if (bitn == 3'd7) st <= STOP; else bitn <= bitn + 1;
                   end else cnt <= cnt + 1;
            STOP:  if (cnt == DIV-1) begin
                       st <= IDLE; cnt <= '0;
                       data <= sh; valid <= 1'b1;
                   end else cnt <= cnt + 1;
        endcase
    end
endmodule

module uart_tx #(parameter int DIV = 868) (
    input  wire clk, rst,
    input  wire [7:0] data,
    input  wire       start,
    output logic      txd,
    output logic      busy
);
    logic [$clog2(DIV)-1:0] cnt;
    logic [3:0] bitn;                 // start + 8 data + stop = 10
    logic [9:0] sh;

    always_ff @(posedge clk) begin
        if (rst) begin txd <= 1'b1; busy <= 1'b0; cnt <= '0; end
        else if (!busy) begin
            txd <= 1'b1;
            if (start) begin
                sh <= {1'b1, data, 1'b0};                        // stop, data(LSB first), start
                busy <= 1'b1; bitn <= 4'd0; cnt <= '0;
            end
        end else begin
            txd <= sh[0];
            if (cnt == DIV-1) begin
                cnt <= '0; sh <= {1'b1, sh[9:1]};
                if (bitn == 4'd9) busy <= 1'b0; else bitn <= bitn + 1;
            end else cnt <= cnt + 1;
        end
    end
endmodule

`default_nettype wire
