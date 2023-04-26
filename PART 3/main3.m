clc;
close all
clear

%% input parameters


Eb_N0_ratios_dB = 0:2:16; % 0:1:0;


Eb_N0_ratio_dB = 0;

Eb_N0_ratio = 10^(Eb_N0_ratio_dB/10); % noise power Eb/N0

Eb_N0_ratios = 10.^(Eb_N0_ratios_dB/10);


bitstream_length = 100000; %length of bitstream

% modulations possible (for part 2):
%   pam 1

modulation = 'qam'; % pam or qam
number_of_bits = 4; % number of bits per symbol

upsampling_rate = 50; %rate of upsamping
Fs = 2e6; % symbol frequency rate
beta = 0.3;
filter_taps = 10*upsampling_rate+1;

% part 3

ppm = 0*1e-6;
F_carrier = 2e9;
CFO = ppm*F_carrier; % carrier frequency offset
phase_offset_degrees  = 0; %[0 2 5 10 15 20 30 45 60]
phase_offset = phase_offset_degrees/360*2*pi; %between 0 and 2*pi
sample_time_shift = 0.5*1e-7; %between zero and 1/Fs = 5e-7

t = 0:1/(upsampling_rate*Fs):((bitstream_length*upsampling_rate - 1)/(upsampling_rate*Fs));

shifts = sample_time_shift/t(2);
if(mod(shifts,1) ~= 0)
    fprintf('upsample rate not compatible with sample time shift \n')
    return
end


%% checking compatibility
if (not(mod(bitstream_length/number_of_bits,1) == 0))
    disp('length of message is incompatible with encoding');
    return
end

%% generating bitstream

fprintf("Generating bitstream...\n")

bit_stream = zeros(bitstream_length,1);
for i = 1:bitstream_length
    bit_stream(i) = round(rand());
end

%% encoding

fprintf("Encoding...\n")

encoded_signal = mapping(bit_stream, number_of_bits, modulation);

if (isempty(encoded_signal))
    fprintf('an error has occured')
    return
end

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

%% adding CFO, phase offset, ...

fprintf('Adding CFO, phase offset and sampling time offset...\n')

t = (0:length(noisy_signals(:,1)) - 1)/(Fs*upsampling_rate);

for j = 1:length(Eb_N0_ratios)
    fprintf("("+j+")\n")
        noisy_signals_CFO(:, j) = noisy_signals(:,j).*exp(1i*(2*pi*CFO*t' + phase_offset));
		%noisy_signals_CFO(:, j) = noisy_signals(:,j);
%     figure
%     plot(t, noisy_signals(:,j))
end

%% inverse nyquist filter and downsamping

fprintf("Filtering & downsampling...\n")

filtered_signal_receiver = filtering_and_downsampling(noisy_signal, upsampling_rate, filter,shifts);

figure
plot(real(filtered_signal_receiver), imag(filtered_signal_receiver), '*')
title('downsampled received signal')
xlabel('Real axis')
ylabel('Imaginairy axis')
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')

filtered_signals_receiver = zeros(length(encoded_signal), length(Eb_N0_ratios));
filtered_signals_receiver_CFO = zeros(length(encoded_signal), length(Eb_N0_ratios));
for j = 1:length(Eb_N0_ratios)
    fprintf("("+j+")\n")
    filtered_signals_receiver(:, j) = filtering_and_downsampling(noisy_signals(:,j), upsampling_rate, filter, shifts);
	filtered_signals_receiver_CFO(:, j) = filtering_and_downsampling(noisy_signals_CFO(:,j), upsampling_rate, filter, shifts);
end

tnew = (0:length(filtered_signals_receiver_CFO(:, 1))-1)/Fs;

for j = 1:length(Eb_N0_ratios)
    fprintf("("+j+")\n")
        filtered_signals_receiver_CFO(:, j) = filtered_signals_receiver_CFO(:,j).*exp(-1i*(2*pi*CFO*tnew'));
%     figure
%     plot(t, noisy_signals(:,j))
end


%% decode

fprintf("Decoding...\n")
decodeds = zeros(length(bit_stream), length(Eb_N0_ratios));
decodeds_CFO = zeros(length(bit_stream), length(Eb_N0_ratios));
for i = 1:length(Eb_N0_ratios)
    fprintf("("+i+")\n")
    decodeds(:,i) = demapping(filtered_signals_receiver(:,i), number_of_bits, modulation);
	decodeds_CFO(:,i) = demapping(filtered_signals_receiver_CFO(:,i), number_of_bits, modulation);
end



%% checking BER

fprintf("Calculating BER...\n")
    
ERRs = zeros(length(Eb_N0_ratios), 1);
BERs = zeros(length(Eb_N0_ratios), 1);
ERRs_CFO = zeros(length(Eb_N0_ratios), 1);
BERs_CFO = zeros(length(Eb_N0_ratios), 1);
for i = 1:length(Eb_N0_ratios)
    fprintf("("+i+")\n")
    dec = decodeds(:,i);
	dec_CFO = decodeds_CFO(:,i);
    for a = 1:length(dec)
        if (dec(a) ~= bit_stream(a))
            ERRs(i) = ERRs(i) + 1;
		end
		if (dec_CFO(a) ~= bit_stream(a))
            ERRs_CFO(i) = ERRs_CFO(i) + 1;
        end
    end
    BERs(i) = ERRs(i)/length(dec);
	BERs_CFO(i) = ERRs_CFO(i)/length(dec);
end

figure
semilogy(10*log10(Eb_N0_ratios), BERs)
hold on
semilogy(10*log10(Eb_N0_ratios), BERs_CFO)
xlabel("Eb/N0 [dB]")
ylabel("BER")
legend("no CFO", "CFO")


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

% CFO modulation

%figure
%hold on
%plot(0:16, BER_pam1(1:11))
% plot(0:10, BER_qam2(1:11))
% plot(0:10, BER_qam4(1:11))
% plot(0:10, BER_qam6(1:11))
%set(gca, 'YScale', 'log')
%ylabel("BER")
%legend('BPSK')
% legend('BPSK', 'QPSK', '16QAM', '64QAM')

