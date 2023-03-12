function [output] = Add_noise(input, Eb_N0_ratio, M, rate)

% additive white gaussian noise
% has real and imaginairy part that are independent from each other
% N0 is noise power

fs = 2e6;
T = 1/fs;
band_width = fs*rate;
Tsample = 1/band_width;

%% calculate noise power

Eb = Tsample/(2*M*1e6)*norm(input)^2;

% Em = abs(input).^2*T; %also take average over time and multilpy with bit duration
% Es = 1/M*sum(Em);
% Eb = 1/log2(M)*Es;
%Eb =  1/(T*2)*norm(input)^2*fs/M; %Es/log2(M); %also take average over time and multilpy with bit duration and divide by 2

N0 = Eb/Eb_N0_ratio %N0 is power spectral density
noise_power = N0*band_width*2 % bandwidth is equal to samplingrate


%% add noise to output

en = randn(length(input), 1)*sqrt(noise_power/2) - j*randn(length(input), 1)*sqrt(noise_power/2);
output = input + en;

%% plot signal with noise in time domain

figure
hold on
plot(0:T/rate:length(output)*T/rate - T/rate, real(fftshift(output)))
%plot(0:T/rate:length(output)*T/rate - T/rate, imag(fftshift(output)))
title('noisy signal')
xlabel('Time (s)')
ylabel('Amplitude')
end