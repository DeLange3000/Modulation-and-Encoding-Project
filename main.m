clc;
close all
clear

%% input parameters

Eb_N0_ratios_dB = -5:1:20;

Eb_N0_ratio_dB = 0;

Eb_N0_ratio = 10^(Eb_N0_ratio_dB/10); % noise power Eb/N0

Eb_N0_ratios = 10.^(Eb_N0_ratios_dB/10);

bitstream_length = 1000000*4; %length of bitstream

% modulations possible:
%   pam 1
%   pam 2
%   qam 4
%   qam 6

modulation = 'qam'; % pam or qam
number_of_bits = 16; % number of bits per symbol

upsampling_rate = 2; %rate of upsamping
Fs = 2e6; % symbol frequency rate
beta = 0.3;

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

% figure
% plot(real(encoded_signal), imag(encoded_signal), '*');
% title('Modulated signal')
% xlabel('Real axis')
% ylabel('Imaginairy axis')
% set(gca, 'XAxisLocation', 'origin')
% set(gca, 'YAxisLocation', 'origin')

%% creating filter

fprintf("Creating filter...\n")

f_filter = Nyquist_filter(Fs, upsampling_rate, length(encoded_signal), beta);
f_filter = sqrt(f_filter); %make it root

%% upsample and filter signal

fprintf("Upsampling & filtering...\n")

filtered_signal = upsampling_and_filtering(encoded_signal, upsampling_rate, f_filter);

%% add noise

fprintf("Adding noise...\n")

noisy_signal = Add_noise(filtered_signal, Eb_N0_ratio, number_of_bits, upsampling_rate);

noisy_signals = [];
for i = 1:length(Eb_N0_ratios)
    fprintf("("+i+")\n")
    noisy_signals = [noisy_signals Add_noise(filtered_signal, Eb_N0_ratios(i), number_of_bits, upsampling_rate)];
end

%% inverse nyquist filter and downsamping

fprintf("Filtering & downsampling...\n")

filtered_signal_receiver = filtering_and_downsampling(noisy_signal, upsampling_rate, f_filter);

% figure
% plot(real(filtered_signal_receiver), imag(filtered_signal_receiver), '*')
% title('downsampled received signal')
% xlabel('Real axis')
% ylabel('Imaginairy axis')
% set(gca, 'XAxisLocation', 'origin')
% set(gca, 'YAxisLocation', 'origin')

filtered_signals_receiver = [];
for i = 1:length(Eb_N0_ratios)
    fprintf("("+i+")\n")
    filtered_signals_receiver = [filtered_signals_receiver filtering_and_downsampling(noisy_signals(:,i), upsampling_rate, f_filter)];
end


%% decode

fprintf("Decoding...\n")

decoded = demapping(filtered_signal_receiver, number_of_bits, modulation);

decodeds = [];
for i = 1:length(Eb_N0_ratios)
    fprintf("("+i+")\n")
    decodeds = [decodeds demapping(filtered_signals_receiver(:,i), number_of_bits, modulation)];
end

%% checking BER

fprintf("Calculating BER...\n")
    
ERR = 0;
for i = 1:length(decoded)
    if (not(decoded(i) == bit_stream(i)))
        ERR = ERR + 1;
    end
end

BER = ERR/length(decoded);

ERRs = zeros(length(Eb_N0_ratios), 1);
BERs = zeros(length(Eb_N0_ratios), 1);
for i = 1:length(Eb_N0_ratios)
    fprintf("("+i+")\n")
    dec = decodeds(:,i);
    for j = 1:length(dec)
        if (dec(j) ~= bit_stream(j))
            ERRs(i) = ERRs(i) + 1;
        end
    end
    BERs(i) = ERRs(i)/length(dec);
end

figure
semilogy(10*log10(Eb_N0_ratios), BERs)
xlabel("Eb/N0 [dB]")
ylabel("BER")











