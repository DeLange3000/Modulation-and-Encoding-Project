clc;
close all
clear

%% input parameters


Eb_N0_ratios_dB = 0:1:16; % 0:1:0;

Eb_N0_ratio_dB = 0;

Eb_N0_ratio = 10^(Eb_N0_ratio_dB/10); % noise power Eb/N0

Eb_N0_ratios = 10.^(Eb_N0_ratios_dB/10);


bitstream_length = 100*1920; %length of bitstream

% modulations possible (for part 2):
%   pam 1

modulation = 'qam'; % pam or qam
number_of_bits = 4; % number of bits per symbol

upsampling_rate = 50; %rate of upsamping
Fs = 2e6; % symbol frequency rate
beta = 0.3;
filter_taps = 10*upsampling_rate+1;

% part 3

%ppms = [0 1 2 10 20 50]*1e-6;
ppm = 10*1e-6;
F_carrier = 2e9;
CFO = ppm*F_carrier; % carrier frequency offset
%CFOs = ppms.*F_carrier;
phase_offset_degrees  = 0; %[0 1 2 5 10 20 30 45 60]
phase_offset = phase_offset_degrees/360*2*pi; %between 0 and 2*pi
sample_time_shift = 0*1/Fs; %between zero and 1/Fs = 5e-7

t = 0:1/(upsampling_rate*Fs):((bitstream_length*upsampling_rate - 1)/(upsampling_rate*Fs));

shifts = sample_time_shift/t(2);

% garner
%error_factors_K = linspace(0.001, 0.01, 6);
error_factor_K = 0.01; 
coulors = ["blue", "red", "green", "black", "cyan", "magenta"];
averaging_length = 10;

% frame sync
pilot_length = 40*number_of_bits;
pilot = randi(2, pilot_length,1)-1;
pilot_int = 200*number_of_bits;

encoded_pilot = mapping(pilot, number_of_bits, modulation);


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
bit_stream = [round(rand(100, 1)); bit_stream];

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
figure
hold on
b = 1;
%for error_factor_K = error_factors_K
%for CFO = CFOs
    Garners_mean = zeros(length(encoded_signal), length(Eb_N0_ratios), averaging_length);
 %   for times = 1:averaging_length
times  = 1;
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
      
	fprintf("Garnering...\n")

filtered_signals_receiver = zeros(length(encoded_signal), length(Eb_N0_ratios));
Garner_errors = zeros(length(encoded_signal), length(Eb_N0_ratios));
for j = 1:length(Eb_N0_ratios)
    fprintf("("+j+")\n")
    [filtered_signals_receiver(:, j), Garner_errors(:,j)] = filtering_and_downsampling(noisy_signals(:,j), upsampling_rate, filter, shifts, error_factor_K);
end
%end

tnew = (0:length(filtered_signals_receiver(:, 1))-1)/Fs;

%% frame acq

fprintf("Frame acquisition...\n")

pilot_int_enc = pilot_int/number_of_bits;
for indddex = 1:length(Eb_N0_ratios)
fprintf("("+indddex+")\n")
[argmaxn, ~] = toa_est(filtered_signals_receiver(1:pilot_int_enc,indddex),encoded_pilot,20,1/2e6);
cfo_recs = [];
for start = 1:pilot_int_enc:length(filtered_signals_receiver)-pilot_int_enc
	frame = filtered_signals_receiver(start+1:start+pilot_int_enc, indddex);
	cfo_recs = [cfo_recs cfo_est(filtered_signals_receiver(start:start+pilot_int_enc, indddex), encoded_pilot, argmaxn, 20, 1/2e6)];
end
cfo_rec = mean(cfo_recs);
%[~, cfo_rec] = toa_est(filtered_signals_receiver(:,indddex),encoded_pilot,20,1/2e6);
%cfo_rec = cfo_est(filtered_signals_receiver(:,indddex), encoded_pilot,0,8,1/Fs)
filtered_signals_receiver(:, indddex) = filtered_signals_receiver(:, indddex).*exp(-1j*2*pi*cfo_rec*tnew');
%toa_est(filtered_signals_receiver(:,indddex),encoded_pilot,8,1/2e6)

% interp_cfo = [];
% for i = 1:length(cfo_ests)-1
% 	range = (1:pilot_int_enc)./(pilot_int_enc);
% 	interp_range = cfo_ests(i+1)*range + cfo_ests(i)*(1.-range);
% 	interp_cfo = [interp_cfo interp_range];
% end
% lastrange = ((pilot_int_enc+1):(2*pilot_int_enc+mod(length(filtered_signals_receiver(:,indddex)), pilot_int_enc)))/pilot_int_enc;
% interp_lastrange = cfo_ests(end)*lastrange + cfo_ests(end-1)*(1.-lastrange);
% interp_cfo = [interp_cfo interp_lastrange];

% recover angles
angle_ests = [];
for pilot_start = argmaxn:pilot_int_enc:length(filtered_signals_receiver)-pilot_int_enc
	pilot_rec = filtered_signals_receiver(pilot_start+1:pilot_start+pilot_length/number_of_bits, indddex);
	angle_ests = [angle_ests median(angle(pilot_rec./encoded_pilot))];
end

figure
%plot(angle_ests, "*")
hold on
title("SNR [dB]: ", Eb_N0_ratios_dB(indddex))

% interpolate angles
interp_angles = zeros(1,argmaxn);
for i = 1:length(angle_ests)-1
%	i = 1;
%	range = (-argmaxn+1:length(filtered_signals_receiver)-argmaxn)./pilot_int_enc;
	range = (1:pilot_int_enc)./(pilot_int_enc);
	if abs(angle_ests(i+1)-angle_ests(i)) < pi/2
		interp_range = angle_ests(i+1)*range + angle_ests(i)*(1.-range);
	elseif angle_ests(i+1)-angle_ests(i) >= pi/2
		interp_range = (angle_ests(i+1)-pi)*range + angle_ests(i)*(1.-range);
	else
		interp_range = (angle_ests(i+1)+pi)*range + angle_ests(i)*(1.-range);
	end
	interp_angles = [interp_angles interp_range];
end
lastrange = ((pilot_int_enc+1):(2*pilot_int_enc+mod(length(filtered_signals_receiver(:,indddex)), pilot_int_enc)-argmaxn))/pilot_int_enc;
interp_lastrange = angle_ests(end)*lastrange + angle_ests(end-1)*(1.-lastrange);
interp_angles = [interp_angles interp_lastrange];
plot(interp_angles)

filtered_signals_receiver(:, indddex) = filtered_signals_receiver(:, indddex).*(exp(1).^(-1i*interp_angles'));


end


    %% decode
    
    fprintf("Decoding...\n")
    decodeds = zeros(length(bit_stream), length(Eb_N0_ratios));
    for i = 1:length(Eb_N0_ratios)
        fprintf("("+i+")\n")
        decodeds(:,i) = demapping(filtered_signals_receiver(:,i), number_of_bits, modulation);
    end
    
%    end
%     %% plot garner errors
%     Garners_std = std(mean(squeeze(Garners_mean)') - (squeeze(Garners_mean)'));
%     Garners_mean =mean(squeeze(Garners_mean)');
%    figure
%    hold on
%     for i = 1: length(Eb_N0_ratios)
%         plot(sample_time_shift*Fs + Garners_mean(i,:), coulors(b))
%         plot((sample_time_shift*Fs + Garners_mean(i, :)) - Garners_std(i,:), coulors(b)+'--')
%         plot((sample_time_shift*Fs + Garners_mean(i, :)) + Garners_std(i,:), coulors(b)+'--')
%         xlabel('symbols')
%         ylabel('Time error (expressed in symbol periods) mean ± std')
%         title('Eb/N0 equal to ', Eb_N0_ratios_dB(i))
%     end
%     b = b + 1;


%end


%legend("K = "+string(error_factors_K(1)),'', '', "K = "+string(error_factors_K(2)),'', '', "K = "+string(error_factors_K(3)),'', '', "K = "+string(error_factors_K(4)),'', '', "K = "+string(error_factors_K(5)),'', '', "K = "+string(error_factors_K(6)))
%legend(string(ppms(1))+" ppm",'', '', string(ppms(2))+" ppm",'', '', string(ppms(3))+" ppm",'', '', string(ppms(4))+" ppm",'', '',string(ppms(5))+" ppm",'', '', string(ppms(6))+" ppm")
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


