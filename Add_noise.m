function [output] = Add_noise(input, Eb_N0_ratio, M, rate)

% additive white gaussian noise
% has real and imaginairy part that are independent from each other
% N0 is noise power

fs = 2e6;
T = 1/fs;

Es = norm(input)^2;
Eb = Es/log2(M);

N0 = Eb/Eb_N0_ratio;

% n1 = wgn(length(input),1, N0);
% n2 = wgn(length(input),1, N0);


en = randn(length(input), 1)*sqrt(N0) - j*randn(length(input), 1)*sqrt(N0);
output = input + en;

figure
plot(0:T/rate:length(output)*T/rate - T/rate, abs(output))
title('noisy signal')
xlabel('Time (s)')
ylabel('Amplitude')
end