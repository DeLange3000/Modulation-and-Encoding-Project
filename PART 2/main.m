clc;
close all
clear

%% input parameters


Eb_N0_ratios_dB = -5:2:15; % 0:1:0;


Eb_N0_ratio_dB = 0;

Eb_N0_ratio = 10^(Eb_N0_ratio_dB/10); % noise power Eb/N0

Eb_N0_ratios = 10.^(Eb_N0_ratios_dB/10);


bitstream_length = 1280000; %length of bitstream

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
size_LDPC = 128;

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
parfor i = 1:length(Eb_N0_ratios)

    fprintf("("+i+")\n")
    filtered_signals_receiver(:, i) = filtering_and_downsampling(noisy_signals(:,i), upsampling_rate, filter);
end


%% decode

fprintf("Decoding...\n")
% hard decoding
decodeds_HD = zeros(length(bit_stream), length(Eb_N0_ratios));
decodeds_SD = zeros(length(bit_stream), length(Eb_N0_ratios));
parfor i = 1:length(Eb_N0_ratios) % parallel processing
    fprintf("("+i+")\n")
    decodeds = demapping(filtered_signals_receiver(:,i), number_of_bits, modulation);
    decodeds_HD(:,i) = hardDecoding(H, decodeds', number_of_iterations_HD)';
    decodeds_SD(:,i) = softDecoding(H, filtered_signals_receiver(:, i)', N0,  number_of_iterations_SD)';
end



%% checking BER

fprintf("Calculating BER...\n")
    

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

%%

figure
hold on
semilogy(10*log10(Eb_N0_ratios), BERs_HD)
semilogy(10*log10(Eb_N0_ratios), BERs_SD)
xlabel("Eb/N0 [dB]")
ylabel("BER")
legend("Hard decoding", "Soft decoding")




