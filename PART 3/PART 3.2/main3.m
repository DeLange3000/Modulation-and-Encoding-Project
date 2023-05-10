clc;
close all
clear

%% input parameters


Eb_N0_ratios_dB = 0:1:16; % 0:1:0;

Eb_N0_ratio_dB = 0;

Eb_N0_ratio = 10^(Eb_N0_ratio_dB/10); % noise power Eb/N0

Eb_N0_ratios = 10.^(Eb_N0_ratios_dB/10);


bitstream_length = 2000; %length of bitstream

% modulations possible (for part 2):
%   pam 1

modulation = 'qam'; % pam or qam
number_of_bits = 2; % number of bits per symbol

upsampling_rate = 50; %rate of upsamping
Fs = 2e6; % symbol frequency rate
beta = 0.3;
filter_taps = 10*upsampling_rate+1;

% part 3

ppm = 10*1e-6;
F_carrier = 2e9;
CFO = ppm*F_carrier; % carrier frequency offset
phase_offset_degrees  = 0; %[0 1 2 5 10 20 30 45 60]
phase_offset = phase_offset_degrees/360*2*pi; %between 0 and 2*pi
sample_time_shift = 0.40*1/Fs; %between zero and 1/Fs = 5e-7

t = 0:1/(upsampling_rate*Fs):((bitstream_length*upsampling_rate - 1)/(upsampling_rate*Fs));

shifts = sample_time_shift/t(2);

% garner
error_factor_K = 0.01;


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
        noisy_signals(:, j) = noisy_signals(:,j).*exp(1i*(2*pi*CFO*t' + phase_offset));
end

%% inverse nyquist filter and downsamping

fprintf("Filtering & downsampling...\n")

%filtered_signal_receiver = filtering_and_downsampling(noisy_signal, upsampling_rate, filter,shifts);

% figure
% plot(real(filtered_signal_receiver), imag(filtered_signal_receiver), '*')
% title('downsampled received signal')
% xlabel('Real axis')
% ylabel('Imaginairy axis')
% set(gca, 'XAxisLocation', 'origin')
% set(gca, 'YAxisLocation', 'origin')

filtered_signals_receiver = zeros(length(encoded_signal), length(Eb_N0_ratios));
Garner_errors = zeros(length(encoded_signal), length(Eb_N0_ratios));
for j = 1:length(Eb_N0_ratios)
    fprintf("("+j+")\n")
    [filtered_signals_receiver(:, j), Garner_errors(:,j)] = filtering_and_downsampling(noisy_signals(:,j), upsampling_rate, filter, shifts, error_factor_K);
end

tnew = (0:length(filtered_signals_receiver(:, 1))-1)/Fs;


%% decode

fprintf("Decoding...\n")
decodeds = zeros(length(bit_stream), length(Eb_N0_ratios));
for i = 1:length(Eb_N0_ratios)
    fprintf("("+i+")\n")
    decodeds(:,i) = demapping(filtered_signals_receiver(:,i), number_of_bits, modulation);
end


%% plot garner errors

for i = 1: length(Eb_N0_ratios)
    figure
    plot(sample_time_shift + Garner_errors(:,i)/Fs)
    xlabel('symbols')
    ylabel('Error of sampling time offset (expressed in symbol periods)')
    title('Eb/N0 equal to ', Eb_N0_ratios_dB(i))
end

%% checking BER

fprintf("Calculating BER...\n")
    
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
hold on
xlabel("Eb/N0 [dB]")
ylabel("BER")
legend("no CFO", "CFO")


%% plot BER


