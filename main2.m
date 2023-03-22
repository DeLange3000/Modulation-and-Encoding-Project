clc;
close all
clear

%% input parameters


Eb_N0_ratios_dB = 0:2:18; % 0:1:0;

Eb_N0_ratio_dB = 0;

Eb_N0_ratio = 10^(Eb_N0_ratio_dB/10); % noise power Eb/N0

Eb_N0_ratios = 10.^(Eb_N0_ratios_dB/10);

bitstream_length = 1*4*5; %length of bitstream

% modulations possible:
%   pam 1
%   pam 2
%   qam 4
%   qam 6

modulation = 'qam'; % pam or qam
number_of_bits = 4; % number of bits per symbol

upsampling_rate = 4; %rate of upsamping
Fs = 2e6; % symbol frequency rate
beta = 0.3;
filter_taps = 51;

%% checking compatibility
if (not(mod(bitstream_length/number_of_bits,1) == 0))
    disp('length of message is incompatible with encoding');
    return
end

%% generating bitstream

fprintf("Generating bitstream...\n")

bitstream = zeros(1, bitstream_length);
for i = 1:bitstream_length
    bitstream(i) = round(rand());
end

%% encoding

fprintf("Encoding...\n")

H = [ 1 1 0 1 1 0 0 1 0 0;
      0 1 1 0 1 1 1 0 0 0;
      0 0 0 1 0 0 0 1 1 1;
      1 1 0 0 0 1 1 0 1 0;
      0 0 1 0 0 1 0 1 0 1 ];

encoded_bitstream = encoding(H, bitstream)

err_i = randi(length(encoded_bitstream),1,1);

for i = err_i
    encoded_bitstream(i) = mod(encoded_bitstream(i)+1, 2);
end

decoded_bitstream = hardDecoding(H, encoded_bitstream)

if decoded_bitstream == bitstream
    disp("yay!")
end



