function [encoded_bitstream, H] = encoding(H,bitstream)
    [M, N] = size(H);
    bitstream_length = length(bitstream);
    
    %------- for premade H -----
    mat = mod(rref(H),2);
    I = mat(:,1:M);
    P = mat(:, M+1:N)';
    G = [P I];

    
    encoded_bitstream = [];
    for i = 1:M:bitstream_length
        block = bitstream(i:i+M-1);
        [encoded, H] = encode_ldpc(block', H, 0);
        encoded_bitstream = [encoded_bitstream [encoded' block]];

        %------- for premade H ----
        %encoded = mod(block*G, 2)
        %encoded_bitstream = [encoded_bitstream encoded];
    end
end