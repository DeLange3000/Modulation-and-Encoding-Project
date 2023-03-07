function [output] = Add_noise(input, Eb_N0_ratio, M, rate)

% additive white gaussian noise
% has real and imaginairy part that are independent from each other
% N0 is noise power

fs = 2e6;
T = 1/fs;
band_width = fs*rate;

Es = norm(input)^2; %also take average over time and multilpy with bit duration
Eb =  1/(T*2)*norm(input)^2*fs/M; %Es/log2(M); %also take average over time and multilpy with bit duration and divide by 2

N0 = Eb/Eb_N0_ratio; %N0 is power spectral density
noise_power = N0*band_width % bandwidth is equal to samplingrate

% n1 = wgn(length(input),1, N0);
% n2 = wgn(length(input),1, N0);


en = randn(length(input), 1)*sqrt(noise_power) - j*randn(length(input), 1)*sqrt(noise_power);
output = input; % + en;

figure
plot(0:T/rate:length(output)*T/rate - T/rate, abs(output))
title('noisy signal')
xlabel('Time (s)')
ylabel('Amplitude')
end