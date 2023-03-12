clc;
close all
clear

%% input parameters
Eb_N0_ratio_dB = 0;

Eb_N0_ratio = 10^(Eb_N0_ratio_dB/10); % noise power Eb/N0

bitstream_length = 100000*4; %length of bitstream

% modulations possible:
%   pam 1
%   pam 2
%   qam 4
%   qam 6

modulation = 'qam'; % pam or qam
number_of_bits = 4; % number of bits per symbol

upsampling_rate = 2; %rate of upsamping
Fs = 2e6; % symbol frequency rate
beta = 0.3;

%% checking compatibility
if (not(mod(bitstream_length/number_of_bits,1) == 0))
    disp('length of message is incompatible with encoding');
    return
end

%% generating bitstream

bit_stream = [];
for i = 1:bitstream_length
    bit_stream = [bit_stream; round(rand())];
end

%% encoding

encoded_signal = mapping(bit_stream, number_of_bits, modulation);

%bit_stream = [0 0 1 1 1 1]
%output = BPSK(1, bit_stream);
%output = QPSK(1, bit_stream);
%output = QAM_16(1, bit_stream);
%output = QAM_64(1, bit_stream);

if (length(encoded_signal) == 0)
    fprintf('an error has occured')
    return
end

%     figure
%     hold on
%     stem(real(output));
%     stem(imag(output));

figure
plot(real(encoded_signal), imag(encoded_signal), '*');
title('Modulated signal')
xlabel('Real axis')
ylabel('Imaginairy axis')
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')

%% creating filter

f_filter = Nyquist_filter(Fs, upsampling_rate, length(encoded_signal), beta);
f_filter = sqrt(f_filter); %make it root

%% upsample and filter signal

filtered_signal = upsampling_and_filtering(encoded_signal, upsampling_rate, f_filter);

%% add noise

noisy_signal = Add_noise(filtered_signal, Eb_N0_ratio, number_of_bits, upsampling_rate);

%% inverse nyquist filter and downsamping

filtered_signal_receiver = filtering_and_downsampling(noisy_signal, upsampling_rate, f_filter);

figure
plot(real(filtered_signal_receiver), imag(filtered_signal_receiver), '*')
title('downsampled received signal')
xlabel('Real axis')
ylabel('Imaginairy axis')
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')


%% decode

decoded = demapping(filtered_signal_receiver, number_of_bits, modulation);

%% checking BER

ERR = 0;
for i = 1:length(decoded)
    if (not(decoded(i) == bit_stream(i)))
        ERR = ERR + 1;
    end
end

BER = ERR/length(decoded)













