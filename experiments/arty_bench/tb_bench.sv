// bench_core の RTL シミュ (iverilog): M=2000 で mism==0 と LFSR 終状態を 検査
`timescale 1ns/1ps
module tb_bench;
    logic clk = 0, rst = 1, start = 0;
    logic running, done;
    logic [31:0] mism, cycles, lfsr_lo;
    bench_core #(.N(16)) dut (.clk(clk), .rst(rst), .start(start), .M(32'd2000),
                              .running(running), .done(done),
                              .mism(mism), .cycles(cycles), .lfsr_lo(lfsr_lo));
    always #5 clk = ~clk;
    initial begin
        repeat (4) @(posedge clk);
        rst = 0;
        @(posedge clk); start = 1; @(posedge clk); start = 0;
        wait (done);
        @(posedge clk);
        $display("SIM mism=%0d cycles=%0d lfsr_lo=%08x", mism, cycles, lfsr_lo);
        $finish;
    end
endmodule
