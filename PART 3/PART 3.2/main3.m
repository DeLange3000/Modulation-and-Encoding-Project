clc;
close all
clear

%% input parameters


Eb_N0_ratios_dB = 6; %0:1:16; % 0:1:0;

Eb_N0_ratio_dB = 0;

Eb_N0_ratio = 10^(Eb_N0_ratio_dB/10); % noise power Eb/N0

Eb_N0_ratios = 10.^(Eb_N0_ratios_dB/10);


bitstream_length = 2500; %length of bitstream

% modulations possible (for part 2):
%   pam 1

modulation = 'qam'; % pam or qam
number_of_bits = 4; % number of bits per symbol

upsampling_rate = 50; %rate of upsamping
Fs = 2e6; % symbol frequency rate
beta = 0.3;
filter_taps = 10*upsampling_rate+1;

% part 3

ppms = [0 5 10 20 50 100]*1e-6;
ppm = 0*1e-6;
F_carrier = 2e9;
CFO = ppm*F_carrier; % carrier frequency offset
CFOs = ppms.*F_carrier;
phase_offset_degrees  = 0; %[0 1 2 5 10 20 30 45 60]
phase_offset = phase_offset_degrees/360*2*pi; %between 0 and 2*pi
sample_time_shift = 0.40*1/Fs; %between zero and 1/Fs = 5e-7

t = 0:1/(upsampling_rate*Fs):((bitstream_length*upsampling_rate - 1)/(upsampling_rate*Fs));

shifts = sample_time_shift/t(2);

% garner
%error_factors_K = linspace(0.001, 0.01, 6);
error_factor_K = 0.01; 
coulors = ["blue", "red", "green", "black", "cyan", "magenta"];
averaging_length = 100000;


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
figure
hold on
b = 1;
%for error_factor_K = error_factors_K
for CFO = CFOs
    Garners_mean = zeros(length(encoded_signal), length(Eb_N0_ratios), averaging_length);
    for times = 1:averaging_length

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
%     figure
%     hold on
%     plot(abs(noisy_signals(:, j)))
%     plot(angle(noisy_signals(:, j)))
        noisy_signals(:, j) = noisy_signals(:,j).*exp(1j*(2*pi*CFO*t' + phase_offset));
%     plot(abs(noisy_signals(:, j)))    
%     plot(angle(noisy_signals(:, j)))
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
        Garners_mean(:,j, times) = Garner_errors(:,j);
    end
        

    tnew = (0:length(filtered_signals_receiver(:, 1))-1)/Fs;
    
    
    %% decode
    
    fprintf("Decoding...\n")
    decodeds = zeros(length(bit_stream), length(Eb_N0_ratios));
    for i = 1:length(Eb_N0_ratios)
        fprintf("("+i+")\n")
        decodeds(:,i) = demapping(filtered_signals_receiver(:,i), number_of_bits, modulation);
    end
    
    end
    %% plot garner errors
    Garners_std = std(mean(squeeze(Garners_mean)') - (squeeze(Garners_mean)'));
    Garners_mean =mean(squeeze(Garners_mean)');
   
    for i = 1: length(Eb_N0_ratios)
        plot(sample_time_shift*Fs + Garners_mean(i,:), coulors(b))
        plot((sample_time_shift*Fs + Garners_mean(i, :)) - Garners_std(i,:), coulors(b)+'--')
        plot((sample_time_shift*Fs + Garners_mean(i, :)) + Garners_std(i,:), coulors(b)+'--')
        xlabel('symbols')
        ylabel('Time error (expressed in symbol periods) mean ± std')
        title('Eb/N0 equal to ', Eb_N0_ratios_dB(i))
    end
    b = b + 1;
end

%legend("K = "+string(error_factors_K(1)),'', '', "K = "+string(error_factors_K(2)),'', '', "K = "+string(error_factors_K(3)),'', '', "K = "+string(error_factors_K(4)),'', '', "K = "+string(error_factors_K(5)),'', '', "K = "+string(error_factors_K(6)))
legend(string(ppms(1))+" ppm",'', '', string(ppms(2))+" ppm",'', '', string(ppms(3))+" ppm",'', '', string(ppms(4))+" ppm",'', '',string(ppms(5))+" ppm",'', '', string(ppms(6))+" ppm")
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


