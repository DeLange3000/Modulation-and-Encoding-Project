clc;
close all
clear

%% input parameters

bitstream_length = 160; % length of bitstream
N0 = 0;
number_of_iterations_HD = 3;
number_of_iterations_SD = 50;

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

H = generate_ldpc(8, 16, 0, 1, 3);
H = double(H); % H can sometimes be a logic operator
[encoded_bitstream, H] = encoding(H, bitstream);


%% adding errors
fprintf("Adding errors...\n")
errs = 0;

err_i = randi(length(encoded_bitstream),1,errs);
err_encoded_bitstream = encoded_bitstream;

for i = err_i
    err_encoded_bitstream(i) = mod(encoded_bitstream(i)+1, 2);
end

%% hard decoding
fprintf("Hard Decoding...\n")

%----- HARD DECODING ----------------------------------
bitstream;
decoded_bitstream = hardDecoding(H, err_encoded_bitstream, number_of_iterations_HD);

if decoded_bitstream == bitstream
    disp("yay!")
else
	disp("Hard Decoded ≠ original")
	disp("Errors at indexes:")
	bitstream;
	decoded_bitstream;
	disp(nonzeros((bitstream~=decoded_bitstream) .* (1:bitstream_length)))
end

BERhard = nnz((bitstream~=decoded_bitstream) .* (1:bitstream_length))/bitstream_length

%% soft decoding	
%------ SOFT DECODING --------------------------------------

fprintf("Adding noise...\n")
err_encoded_bitstream = encoded_bitstream;

for i = 1:length(encoded_bitstream)
    if(encoded_bitstream(i) == 0)
        err_encoded_bitstream(i) = -1;
    else
        err_encoded_bitstream(i) = 1;
    end
end

err_encoded_bitstream = err_encoded_bitstream + N0/2*randn(1, length(encoded_bitstream));


fprintf("Soft Decoding...\n")
decoded_bitstream = softDecoding(H, err_encoded_bitstream, N0, number_of_iterations_SD);

% if decoded_bitstream == bitstream
%     disp("yay!")
% else
% 	disp("Soft Decoded ≠ original")
% 	disp("Errors at indexes:")
% 	bitstream;
% 	decoded_bitstream;
% 	disp(nonzeros((bitstream~=decoded_bitstream) .* (1:bitstream_length)))
% end

BERsoft = nnz((bitstream~=decoded_bitstream) .* (1:bitstream_length))/bitstream_length



