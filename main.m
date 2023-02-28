clc;
close all
clear

%% inputs

Eb_N0_ratio = 0; % noise power Eb/N0

bitstream_length = 16*6; %length of bitstream

% modulations possible:
%   pam 1
%   pam 2
%   qam 4
%   qam 6

modulation = 'qam'; % pam or qam
number_of_bits = 6; % number of bits per symbol

upsampling_rate = 100; %rate of upsamping
Fs = 2e6; % symbol frequency rate

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

%% upsample and filter signal

filtered_signal = Nyquist_filter(encoded_signal,upsampling_rate);

%% add noise

noisy_signal = Add_noise(filtered_signal, Eb_N0_ratio, 2^number_of_bits);

%% inverse nyquist filter and downsamping

filtered_signal_receiver = inverse_Nyquist_filter(noisy_signal, upsampling_rate);
figure
plot(real(filtered_signal_receiver), imag(filtered_signal_receiver), '*')
title('decoded received signal')
xlabel('Real axis')
ylabel('Imaginairy axis')
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')


%% decode

decoded = demapping(filtered_signal_receiver, number_of_bits, modulation);













