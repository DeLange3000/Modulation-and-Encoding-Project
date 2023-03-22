function [encoded_bitstream] = hardEncoding(H,bitstream)
    [M, N] = size(H)
    bitstream_length = length(bitstream);
    
    mat = mod(rref(H),2);
    I = mat(:,1:M);
    P = mat(:, M+1:N)';
    G = [P I];
    
    encoded_bitstream = [];
    for i = 1:5:bitstream_length
        block = bitstream(i:i+M-1);
        encoded = mod(block*G, 2);
        encoded_bitstream = [encoded_bitstream encoded];
    end
end