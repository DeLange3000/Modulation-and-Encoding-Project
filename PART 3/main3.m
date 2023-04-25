clc;
close all
clear

%% input parameters


Eb_N0_ratios_dB = 0:1:7; % 0:1:0;


Eb_N0_ratio_dB = 0;

Eb_N0_ratio = 10^(Eb_N0_ratio_dB/10); % noise power Eb/N0

Eb_N0_ratios = 10.^(Eb_N0_ratios_dB/10);


bitstream_length = 100000; %length of bitstream

% modulations possible (for part 2):
%   pam 1

modulation = 'pam'; % pam or qam
number_of_bits = 1; % number of bits per symbol

upsampling_rate = 4; %rate of upsamping
Fs = 2e6; % symbol frequency rate
beta = 0.3;
filter_taps = 101;

% part 2
number_of_iterations_HD = 1;
number_of_iterations_SD = 21;
size_LDPC = 100;

H = generate_ldpc(size_LDPC, 2*size_LDPC, 0, 1, 2);
H = double(H); % H can sometimes be a logic operator
%% checking compatibility
if (not(mod(bitstream_length/number_of_bits,1) == 0))
    disp('length of message is incompatible with encoding');
    return
end

%% generating bitstream

fprintf("Generating bitstream...\n")

bit_stream = zeros(1, bitstream_length);
for i = 1:bitstream_length
    bit_stream(i) = round(rand());
end

%% encoding

fprintf("Encoding...\n")

[bit_stream_encoded, H] = encoding(H, bit_stream);
encoded_signal = mapping(bit_stream_encoded', number_of_bits, modulation);


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
% only add real noise!
noisy_signal = Add_noise(filtered_signal, Eb_N0_ratio, number_of_bits, upsampling_rate);

noisy_signals = zeros(length(filtered_signal), length(Eb_N0_ratios));
for i = 1:length(Eb_N0_ratios)
    fprintf("("+i+")\n")
    
    [noisy_signals(:, i), N0] = Add_noise(filtered_signal, Eb_N0_ratios(i), number_of_bits, upsampling_rate);
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
% hard decoding
decodeds_HD = zeros(length(bit_stream), length(Eb_N0_ratios));
decodeds_SD = zeros(length(bit_stream), length(Eb_N0_ratios));
for i = 1:length(Eb_N0_ratios)
    fprintf("("+i+")\n")
    decodeds = demapping(filtered_signals_receiver(:,i), number_of_bits, modulation);
    decodeds_HD(:,i) = hardDecoding(H, decodeds', number_of_iterations_HD)';
    decodeds_SD(:,i) = softDecoding(H, filtered_signals_receiver(:, i)', N0,  number_of_iterations_SD)';
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

ERRs_HD = zeros(length(Eb_N0_ratios), 1);
ERRs_SD = zeros(length(Eb_N0_ratios), 1);
BERs_HD = zeros(length(Eb_N0_ratios), 1);
BERs_SD = zeros(length(Eb_N0_ratios), 1);
for i = 1:length(Eb_N0_ratios)
    fprintf("("+i+")\n")
    dec_HD = decodeds_HD(:,i);
    dec_SD = decodeds_SD(:,i);
    for a = 1:length(dec_HD)
        if (dec_HD(a) ~= bit_stream(a))
            ERRs_HD(i) = ERRs_HD(i) + 1;
        end
        if (dec_SD(a) ~= bit_stream(a))
            ERRs_SD(i) = ERRs_SD(i) + 1;
        end
    end
    BERs_HD(i) = ERRs_HD(i)/length(dec_HD);
    BERs_SD(i) = ERRs_SD(i)/length(dec_SD);
end

figure
hold on
semilogy(10*log10(Eb_N0_ratios), BERs_HD)
semilogy(10*log10(Eb_N0_ratios), BERs_SD)
xlabel("Eb/N0 [dB]")
ylabel("BER")
legend("Hard decoding", "Soft decoding")


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



BER_HD_5_iterations = [0.0964800000000000
0.0573600000000000
0.0317000000000000
0.0133500000000000
0.00519000000000000
0.00141000000000000
0.000320000000000000
0.000130000000000000];

BER_SD_15_iterations = [0.343340000000000
0.225310000000000
0.109590000000000
0.0301700000000000
0.00838000000000000
0.00151000000000000
0.000370000000000000
0.000150000000000000];

BER_HD_1_iteration = [0.0860900000000000
0.0485300000000000
0.0260600000000000
0.0103400000000000
0.00316000000000000
0.000810000000000000
0.000120000000000000
1.00000000000000e-05];

BER_SD_1_iteration =[0.0386600000000000
0.0191700000000000
0.00759000000000000
0.00244000000000000
0.000620000000000000
9.00000000000000e-05
1.00000000000000e-05
0];

BER_HD_2_iterations = [0.0955100000000000
0.0577600000000000
0.0292700000000000
0.0125100000000000
0.00377000000000000
0.000980000000000000
0.000180000000000000
4.00000000000000e-05];

BER_SD_3_iterations = [0.0297000000000000
0.0108800000000000
0.00286000000000000
0.000370000000000000
4.00000000000000e-05
3.00000000000000e-05
0
0]

BER_HD_7_iterations = [0.0812100000000000
0.0488900000000000
0.0247300000000000
0.0103600000000000
0.00333000000000000
0.000790000000000000
0.000110000000000000
2.00000000000000e-05];
BER_SD_7_iterations = [0.129010000000000
0.0591000000000000
0.0174300000000000
0.00262000000000000
0.000680000000000000
0
0
0];

BER_HD_9_iterations = [0.0841300000000000
0.0511600000000000
0.0255000000000000
0.0112600000000000
0.00425000000000000
0.000880000000000000
0.000250000000000000
5.00000000000000e-05];

BER_SD_11_iterations = [0.270450000000000
0.163840000000000
0.0525700000000000
0.0109200000000000
0.00264000000000000
0.000210000000000000
0.000100000000000000
5.00000000000000e-05];

BER_SD_18_iterations = [0.384020000000000
0.252950000000000
0.103090000000000
0.0250600000000000
0.00862000000000000
0.00251000000000000
0.000530000000000000
0.000140000000000000];

BER_SD_21_iterations = [0.405280000000000
0.263250000000000
0.112490000000000
0.0259100000000000
0.00449000000000000
0.00238000000000000
0.000160000000000000
0];

figure
hold on
plot(0:5, BER_pam1(1:6))
% plot(0:10, BER_qam2(1:11))
% plot(0:10, BER_qam4(1:11))
% plot(0:10, BER_qam6(1:11))
%plot(0:5, BER_SD_15_iterations(1:6))
plot(0:5, BER_HD_1_iteration(1:6))
plot(0:5, BER_HD_2_iterations(1:6))
plot(0:5, BER_HD_5_iterations(1:6))
plot(0:5, BER_HD_7_iterations(1:6))
plot(0:5, BER_HD_9_iterations(1:6))
set(gca, 'YScale', 'log')
xlabel('Eb/N0')
ylabel("BER")
%legend('BPSK', 'QPSK', '16QAM', '64QAM')
legend('BPSK', '1 iteration', '2 iterations', '5 iterations', '7 iterations', '9 iterations')
title('Hard Decoding of pam')

figure
hold on
plot(0:5, BER_pam1(1:6))
% plot(0:10, BER_qam2(1:11))
% plot(0:10, BER_qam4(1:11))
% plot(0:10, BER_qam6(1:11))
%plot(0:5, BER_SD_15_iterations(1:6))
plot(0:5, BER_SD_1_iteration(1:6))
plot(0:5, BER_SD_3_iterations(1:6))
plot(0:5, BER_SD_7_iterations(1:6))
plot(0:5, BER_SD_11_iterations(1:6))
plot(0:5, BER_SD_15_iterations(1:6))
plot(0:5, BER_SD_18_iterations(1:6))
plot(0:5, BER_SD_21_iterations(1:6))
set(gca, 'YScale', 'log')
xlabel('Eb/N0')
ylabel("BER")
%legend('BPSK', 'QPSK', '16QAM', '64QAM')
legend('BPSK', '1 iteration', '3 iterations', '7 iterations', '11 iterations', '15 iterations', '18 iterations', '21 iterations')
title('Soft Decoding of pam')

