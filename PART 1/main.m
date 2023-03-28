clc;
close all
clear

%% input parameters


Eb_N0_ratios_dB = 0:1:30; % 0:1:0;


Eb_N0_ratio_dB = 0;

Eb_N0_ratio = 10^(Eb_N0_ratio_dB/10); % noise power Eb/N0

Eb_N0_ratios = 10.^(Eb_N0_ratios_dB/10);

bitstream_length = 1000000*6; %length of bitstream

% modulations possible:
%   pam 1
%   qam 2
%   qam 4
%   qam 6

modulation = 'qam'; % pam or qam
number_of_bits = 6; % number of bits per symbol

upsampling_rate = 4; %rate of upsamping
Fs = 2e6; % symbol frequency rate
beta = 0.3;
filter_taps = 101;

%% checking compatibility
if (not(mod(bitstream_length/number_of_bits,1) == 0))
    disp('length of message is incompatible with encoding');
    return
end

%% generating bitstream

fprintf("Generating bitstream...\n")

bit_stream = zeros(bitstream_length, 1);
for i = 1:bitstream_length
    bit_stream(i) = round(rand());
end

%% encoding

fprintf("Encoding...\n")

encoded_signal = mapping(bit_stream, number_of_bits, modulation);

%bit_stream = [0 0 1 1 1 1]
%output = BPSK(1, bit_stream);
%output = QPSK(1, bit_stream);
%output = QAM_16(1, bit_stream);
%output = QAM_64(1, bit_stream);

if (isempty(encoded_signal))
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

fprintf("Creating filter...\n")

filter = Nyquist_filter(Fs, upsampling_rate, length(encoded_signal), beta, filter_taps);

%% upsample and filter signal

fprintf("Upsampling & filtering...\n")

filtered_signal = upsampling_and_filtering(encoded_signal, upsampling_rate, filter);

%% add noise

fprintf("Adding noise...\n")

noisy_signal = Add_noise(filtered_signal, Eb_N0_ratio, number_of_bits, upsampling_rate);

noisy_signals = zeros(length(filtered_signal), length(Eb_N0_ratios));
for i = 1:length(Eb_N0_ratios)
    fprintf("("+i+")\n")
    
    noisy_signals(:, i) = Add_noise(filtered_signal, Eb_N0_ratios(i), number_of_bits, upsampling_rate);
end

%% inverse nyquist filter and downsamping

fprintf("Filtering & downsampling...\n")

filtered_signal_receiver = filtering_and_downsampling(noisy_signal, upsampling_rate, filter);

figure
plot(real(filtered_signal_receiver), imag(filtered_signal_receiver), '*')
title('downsampled received signal')
xlabel('Real axis')
ylabel('Imaginairy axis')
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')

filtered_signals_receiver = zeros(length(encoded_signal), length(Eb_N0_ratios));
for i = 1:length(Eb_N0_ratios)
    fprintf("("+i+")\n")
    filtered_signals_receiver(:, i) = filtering_and_downsampling(noisy_signals(:,i), upsampling_rate, filter);
end


%% decode

fprintf("Decoding...\n")

%decoded = demapping(filtered_signal_receiver, number_of_bits, modulation);

decodeds = zeros(length(bit_stream), length(Eb_N0_ratios));
for i = 1:length(Eb_N0_ratios)
    fprintf("("+i+")\n")
    decodeds(:, i) = demapping(filtered_signals_receiver(:,i), number_of_bits, modulation);
end

%% checking BER

fprintf("Calculating BER...\n")
    
% ERR = 0;
% for i = 1:length(decoded)
%     if (not(decoded(i) == bit_stream(i)))
%         ERR = ERR + 1;
%     end
% end
% 
% BER = ERR/length(decoded);

ERRs = zeros(length(Eb_N0_ratios), 1);
BERs = zeros(length(Eb_N0_ratios), 1);
for i = 1:length(Eb_N0_ratios)
    fprintf("("+i+")\n")
    dec = decodeds(:,i);
    for a = 1:length(dec)
        if (dec(a) ~= bit_stream(a))
            ERRs(i) = ERRs(i) + 1;
        end
    end
    BERs(i) = ERRs(i)/length(dec);
end

figure
semilogy(10*log10(Eb_N0_ratios), BERs)
xlabel("Eb/N0 [dB]")
ylabel("BER")


%% plot BER

BER_pam1 = [0.0789150000000000 0.055838000000000 0.0377490000000000 0.0226910000000000 0.0124200000000000 0.00586100000000000 0.00241900000000000 0.000763000000000000 0.000219000000000000 3.40000000000000e-05 4.00000000000000e-06];
BER_qam2 = [0.0789425000000000
0.0559850000000000
0.0377985000000000
0.0228615000000000
0.0126375000000000
0.00605500000000000
0.00233850000000000
0.000774000000000000
0.000178000000000000
3.55000000000000e-05
1.50000000000000e-06
5.00000000000000e-07];
BER_qam4 = [0.140860500000000
0.118954000000000
0.0977317500000000
0.0773932500000000
0.0585857500000000
0.0419385000000000
0.0279185000000000
0.0170302500000000
0.00930425000000000
0.00441025000000000
0.00177725000000000
0.000563000000000000
0.000143250000000000
2.45000000000000e-05
4.25000000000000e-06];
BER_qam6 = [0.199876000000000
0.177850500000000
0.157242333333333
0.137272333333333
0.118624500000000
0.100855666666667
0.0840655000000000
0.0677870000000000
0.0523825000000000
0.0386655000000000
0.0269078333333333
0.0171495000000000
0.0100190000000000
0.00514400000000000
0.00231250000000000
0.000870000000000000
0.000261333333333333
6.38333333333333e-05
9.50000000000000e-06
1.16666666666667e-06
1.66666666666667e-07];

figure
hold on
plot(0:10, BER_pam1(1:11))
plot(0:10, BER_qam2(1:11))
plot(0:10, BER_qam4(1:11))
plot(0:10, BER_qam6(1:11))
set(gca, 'YScale', 'log')
xlabel('Eb/N0')
ylabel("BER")
legend('BPSK', 'QPSK', '16QAM', '64QAM')
