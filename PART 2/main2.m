clc;
close all
clear

%% input parameters

bitstream_length = 400; %length of bitstream

%% generating bitstream

fprintf("Generating bitstream...\n")

bitstream = zeros(1, bitstream_length);
for i = 1:bitstream_length
    bitstream(i) = round(rand());
end

%% encoding

fprintf("Encoding...\n")

% H = [ 1 1 0 1 1 0 0 1 0 0;
%       0 1 1 0 1 1 1 0 0 0;
%       0 0 0 1 0 0 0 1 1 1;
%       1 1 0 0 0 1 1 0 1 0;
%       0 0 1 0 0 1 0 1 0 1 ];

H = generate_ldpc(100, 200, 0, 1, 3);
H = double(H); % H can sometimes be a logic operator
[encoded_bitstream, H] = encoding(H, bitstream);


fprintf("Adding errors...\n")
errs = 1;

err_i = randi(length(encoded_bitstream),1,errs);

for i = err_i
    encoded_bitstream(i) = mod(encoded_bitstream(i)+1, 2);
end

fprintf("Decoding...\n")
bitstream;
decoded_bitstream = hardDecoding(H, encoded_bitstream);

if decoded_bitstream == bitstream
    disp("yay!")
else
	disp("Decoded ≠ original")
	disp("Errors at indexes:")
	bitstream;
	decoded_bitstream;
	disp(nonzeros((bitstream~=decoded_bitstream) .* (1:bitstream_length)))
end



