close all
clear
clc
%%

N01 = 0:10;
BER1 = runMain(1000, N01, 10);

figure
semilogy(10*log10(N01), BER1)
xlabel("Eb/N0 [dB]")
ylabel("BER")
