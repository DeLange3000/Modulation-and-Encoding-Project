function [dec_bs] = hardDecoding(H, rec_bs)
    [M, N] = size(H);
    rec_bs_l = length(rec_bs);
    
    dec_bs = [];

    for block_i = 1:N:rec_bs_l
        block = rec_bs(block_i:block_i+N-1);

        v_nodes = zeros(M+1, N);
        v_nodes(1, :) = block;
        
        for i = 1:M
            % #rogier is genius
            arr = nonzeros(H(i,:) .* (1:N))';

            for j = arr
                arr_sel = arr(arr~=j);
                selected = block(arr_sel);
                if mod(sum(selected), 2) == 1
                    v_nodes(i+1, j) = 1;
                end
            end
        end
        cor_block = sum(v_nodes, 1);
        cor_block = round(cor_block/M);
        dec_block = cor_block(M+1:N);
        dec_bs = [dec_bs dec_block];
    end
    
end