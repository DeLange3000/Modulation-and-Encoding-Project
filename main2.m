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

bit_stream = zeros(bitstream_length, 1);
for i = 1:bitstream_length
    bit_stream(i) = round(rand());
end

%% encoding

fprintf("Encoding...\n")

H = makeLdpc(5, 10, 0, 0, 2);

encoded_bitstream = [];
for i = 1:5:bitstream_length-5
    block = bit_stream(i:i+4);
    encoded = mod(block'*H, 2);
    encoded_bitstream = [encoded_bitstream; encoded'];
end

