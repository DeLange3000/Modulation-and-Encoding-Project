clc;
close all
clear

%% input parameters


Eb_N0_ratios_dB = 0:1:20; % 0:1:0;

Eb_N0_ratio_dB = 0;

Eb_N0_ratio = 10^(Eb_N0_ratio_dB/10); % noise power Eb/N0

Eb_N0_ratios = 10.^(Eb_N0_ratios_dB/10);

len_rep = 100
bitstream_length = 10*1920; %length of bitstream

% modulations possible (for part 2):
%   pam 1

modulation = 'qam'; % pam or qam
number_of_bits = 4; % number of bits per symbol

upsampling_rate = 10; %rate of upsamping
Fs = 2e6; % symbol frequency rate
beta = 0.3;
filter_taps = 10*upsampling_rate+1;

% part 3

%ppms = [0 1 2 10 20 50]*1e-6;
ppm = 0*1e-6;
F_carrier = 2e9;
CFO = ppm*F_carrier; % carrier frequency offset
%CFOs = ppms.*F_carrier;
phase_offset_degrees  = 0; %[0 1 2 5 10 20 30 45 60]
phase_offset = phase_offset_degrees/360*2*pi; %between 0 and 2*pi
sample_time_shift = 0*1/Fs; %between zero and 1/Fs = 5e-7

t = 0:1/(upsampling_rate*Fs):((bitstream_length*upsampling_rate - 1)/(upsampling_rate*Fs));

shifts = sample_time_shift/t(2);

% garner
%error_facto	rs_K = linspace(0.001, 0.01, 6);
error_factor_K = 0.01; 
coulors = ["blue", "red", "green", "black", "cyan", "magenta"];
averaging_length = 10;

% frame sync
K = [1 8 16];
pilot_length = 20*number_of_bits;
pilot = randi(2, pilot_length,1)-1;
pilot_int = 60*number_of_bits;
prefix = round(rand(number_of_bits*round(25*rand()), 1));
toa = length(prefix)/number_of_bits;

encoded_pilot = mapping(pilot, number_of_bits, modulation);
cfo_res = zeros(len_rep, length(Eb_N0_ratios), 3);
arg = zeros(len_rep, length(Eb_N0_ratios), 3);
for kval =  1:3
for lolol = 1:len_rep
lolol
kval
%% checking compatibility
if (not(mod(bitstream_length/number_of_bits,1) == 0))
    disp('length of message is incompatible with encoding');
    return
end
if (not(mod(bitstream_length/(pilot_int-pilot_length),1) == 0))
    disp('length of message is incompatible with encoding (pilot int etc)');
    return
end
if (not(mod(pilot_int/number_of_bits,1) == 0))
    disp('pilot interval is incompatible with encoding');
    return
end
if (not(mod(pilot_length/number_of_bits,1) == 0))
    disp('pilot length is incompatible with encoding');
    return
end

%% generating bitstream

fprintf("Generating bitstream...\n")

bit_stream = zeros(bitstream_length,1);
for i = 1:bitstream_length
    bit_stream(i) = round(rand());
end

%% insert pilots

fprintf("Inserting pilot sequences...\n")

for i = 1:pilot_int:bitstream_length+pilot_length*ceil(bitstream_length/(pilot_int-pilot_length))
	begin = bit_stream(1:i-1);
	rest = bit_stream(i:end);
	bit_stream = [begin; pilot; rest];
end
bit_stream = [prefix; bit_stream];

%% encoding

fprintf("Encoding...\n")

encoded_signal = mapping(bit_stream, number_of_bits, modulation);

if (isempty(encoded_signal))
    fprintf('an error has occured')
    return
end

% figure
% plot(real(encoded_signal), imag(encoded_signal), '*');
% title('Modulated signal')
% xlabel('Real axis')
% ylabel('Imaginairy axis')
% set(gca, 'XAxisLocation', 'origin')
% set(gca, 'YAxisLocation', 'origin')

%% creating filter

fprintf("Creating filter...\n")

filter = Nyquist_filter(Fs, upsampling_rate, length(encoded_signal), beta, filter_taps);

%% upsample and filter signal

fprintf("Upsampling & filtering...\n")

filtered_signal = upsampling_and_filtering(encoded_signal, upsampling_rate, filter);

%% add noise
%figure
hold on
b = 1;
%for error_factor_K = error_factors_K
%for CFO = CFOs
    Garners_mean = zeros(length(encoded_signal), length(Eb_N0_ratios), averaging_length);
 %   for times = 1:averaging_length
times  = 1;
%fprintf("Adding noise...\n")
noisy_signal = Add_noise(filtered_signal, Eb_N0_ratio, number_of_bits, upsampling_rate);

noisy_signals = zeros(length(filtered_signal), length(Eb_N0_ratios));
for i = 1:length(Eb_N0_ratios)
 %   fprintf("("+i+")\n")
    noisy_signals(:, i) = Add_noise(filtered_signal, Eb_N0_ratios(i), number_of_bits, upsampling_rate);
end

%% adding CFO, phase offset, ...

fprintf('Adding CFO, phase offset and sampling time offset...\n')

t = (0:length(noisy_signals(:,1)) - 1)/(Fs*upsampling_rate);    

for j = 1:length(Eb_N0_ratios)
   % fprintf("("+j+")\n")
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
     %   fprintf("("+j+")\n")
        [filtered_signals_receiver(:, j), Garner_errors(:,j)] = filtering_and_downsampling(noisy_signals(:,j), upsampling_rate, filter, shifts, error_factor_K);
        Garners_mean(:,j, times) = Garner_errors(:,j);
    end
      
%	fprintf("Garnering...\n")

filtered_signals_receiver = zeros(length(encoded_signal), length(Eb_N0_ratios));
Garner_errors = zeros(length(encoded_signal), length(Eb_N0_ratios));
for j = 1:length(Eb_N0_ratios)
%    fprintf("("+j+")\n")
    [filtered_signals_receiver(:, j), Garner_errors(:,j)] = filtering_and_downsampling(noisy_signals(:,j), upsampling_rate, filter, shifts, error_factor_K);
end
%end

tnew = (0:length(filtered_signals_receiver(:, 1))-1)/Fs;

%% frame acq

%fprintf("Frame acquisition...\n")

pilot_int_enc = pilot_int/number_of_bits;
for indddex = 1:length(Eb_N0_ratios)
%fprintf("("+indddex+")\n")
[argmaxn, ~] = toa_est(filtered_signals_receiver(1:pilot_int_enc,indddex),encoded_pilot,K(kval),1/2e6);
%[arg(lolol,indddex,kval), cfo_rec(lolol, indddex, kval)] = toa_est(filtered_signals_receiver(1:pilot_int_enc,indddex),encoded_pilot,K(kval),1/2e6);
cfo_recs = [];
args_now = [];
for start = 1:pilot_int_enc:length(filtered_signals_receiver)-pilot_int_enc
	frame = filtered_signals_receiver(start+1:start+pilot_int_enc, indddex);
	[arg_now, cfo_rec_now] = toa_est(filtered_signals_receiver(start:start+pilot_int_enc, indddex), encoded_pilot, K(kval), 1/2e6);
	cfo_recs = [cfo_recs cfo_rec_now];
	args_now = [args_now mod(arg_now, pilot_int_enc)];
end
cfo_rec = mean(cfo_recs);
argmaxn = mean(args_now);
cfo_res(lolol, indddex, kval) = cfo_rec;
arg(lolol, indddex, kval) = argmaxn;

end

end
end
%% plot
avgs_cfo = zeros(length(Eb_N0_ratios), 3);
avgs_toa = zeros(length(Eb_N0_ratios), 3);

% cfo_res_noOut = rmoutliers(cfo_res);
% arg_noOut = rmoutliers(arg);

for kval = 1:3
for snr = 1:length(Eb_N0_ratios)
	cfo_res_noOut = (cfo_res(:, snr, kval));
	arg_noOut = (arg(:, snr, kval));

	avgs_cfo(snr, kval) = std(cfo_res_noOut);
	avgs_toa(snr, kval) = std(arg_noOut);
end
end

ppms = ((avgs_cfo))*1e6/F_carrier;


figure
plot(Eb_N0_ratios_dB, ppms(:,1))
hold on
plot(Eb_N0_ratios_dB, ppms(:,2))
plot(Eb_N0_ratios_dB, ppms(:,3))
title("CFO error with varying averaging window K, N = 20")
xlabel("SNR [dB]")
ylabel("CFO error [ppm]")
legend("K = 1", "K = 8", "K = 16")

% figure
% hold on
% plot(Eb_N0_ratios_dB, ppms(:,2))
% plot(Eb_N0_ratios_dB, ppms(:,3))
% title("CFO error with varying averaging window K, N = 20")
% xlabel("SNR [dB]")
% ylabel("CFO error [ppm]")
% legend("K = 8", "K = 16")

figure
plot(Eb_N0_ratios_dB, avgs_toa(:,1))
hold on
plot(Eb_N0_ratios_dB, avgs_toa(:,2))
plot(Eb_N0_ratios_dB, avgs_toa(:,3))
title("Time error with varying averaging window K, N = 20")
xlabel("SNR [dB]")
ylabel("time error [samples]")
legend("K = 1", "K = 8", "K = 16")

