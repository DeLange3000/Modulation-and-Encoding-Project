function [BER] = runMain(bs_length, N0s, its_soft)

%% generating bitstream

fprintf("Generating bitstream...\n")

bitstream = zeros(1, bs_length);
for i = 1:bs_length
    bitstream(i) = round(rand());
end

%% encoding

fprintf("Encoding...\n")

% H = [ 1 1 0 1 1 0 0 1 0 0;
%       0 1 1 0 1 1 1 0 0 0;
%       0 0 0 1 0 0 0 1 1 1;
%       1 1 0 0 0 1 1 0 1 0;
%       0 0 1 0 0 1 0 1 0 1 ];

H = generate_ldpc(5, 10, 0, 1, 3);
H = double(H); % H can sometimes be a logic operator
[encoded_bitstream, H] = encoding(H, bitstream);

%% soft decoding	
%------ SOFT DECODING --------------------------------------

fprintf("Adding noise...\n")
err_encoded_bitstreams = zeros(length(N0s),length(encoded_bitstream));

% for i = 1:length(encoded_bitstream)
%     if(encoded_bitstream(i) == 0)
%         err_encoded_bitstreams(i) = -1;
%     else
%         err_encoded_bitstreams(i) = 1;
%     end
% end

for i = 1:length(N0s)
	err_encoded_bitstreams(i, :) = 2*(encoded_bitstream-0.5);
	err_encoded_bitstreams(i, :) = err_encoded_bitstreams(i, :) + N0s(i)/2*randn(1, length(encoded_bitstream));
end

fprintf("Soft Decoding...\n")
decoded_bitstreams = zeros(length(N0s),bs_length);
BER = zeros(1, length(N0s));
for i = 1:length(N0s)
	decoded_bitstreams(i, :) = softDecoding(H, err_encoded_bitstreams(i, :), N0s(i), its_soft);
	BER(i) = nnz((bitstream~=decoded_bitstreams(i, :)))/bs_length;
end

% if decoded_bitstream == bitstream
%     disp("yay!")
% else
% 	disp("Soft Decoded ≠ original")
% 	disp("Errors at indexes:")
% 	bitstream;
% 	decoded_bitstream;
% 	disp(nonzeros((bitstream~=decoded_bitstream) .* (1:bitstream_length)))
% end

end