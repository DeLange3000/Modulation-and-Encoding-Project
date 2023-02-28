function [output] = Add_noise(input, Eb_N0_ratio, M)

% additive white gaussian noise
% has real and imaginairy part that are independent from each other
% N0 is noise power

Es = norm(input)^2;
Eb = Es/log2(M);

N0 = Eb_N0_ratio/Eb;

n1 = wgn(1,length(input), N0);
n2 = wgn(1,length(input), N0);


en = n1 - j*n2;
output = input; % + en;

figure
plot(output)
title('noisy signal')
xlabel('Time (s)')
ylabel('Amplitude')
end