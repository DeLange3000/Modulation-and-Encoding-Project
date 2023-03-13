function [output] = Add_noise(input, Eb_N0_ratio, M, rate)

% additive white gaussian noise
% has real and imaginairy part that are independent from each other
% N0 is noise power

fs = 2e6;
T = 1/fs;
band_width = fs*rate;
Tsample = 1/band_width;

%% calculate noise power

Es = 0;
for i =1:length(input)
    Es = abs(input(i))^2;
end
Es;
Eb = Es*Tsample/(2*M);

N0 = Eb/Eb_N0_ratio; %N0 is power spectral density
noise_power = N0*band_width; % bandwidth is equal to samplingrate


%% add noise to output

en = randn(length(input), 1)*sqrt(noise_power) - j*randn(length(input), 1)*sqrt(noise_power);
output = input + en;

%% plot signal with noise in time domain
% 
% figure
% hold on
% plot(0:T/rate:length(output)*T/rate - T/rate, real(fftshift(output)))
% %plot(0:T/rate:length(output)*T/rate - T/rate, imag(fftshift(output)))
% title('noisy signal')
% xlabel('Time (s)')
% ylabel('Amplitude')
end