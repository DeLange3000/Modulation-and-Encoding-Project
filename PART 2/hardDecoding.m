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
        %get non_zero indices of H

        for i = M+1:N
        H_nonzero = H(:,i)'.*(2:M+1); % add index due to offset in v_nodes
        H_nonzero = nonzeros(H_nonzero);
        cor_block = sum(v_nodes([1; H_nonzero],i));
        cor_block = round(cor_block/(length(H_nonzero)+1));
        dec_block = cor_block;
        dec_bs = [dec_bs dec_block];
        end
    end
    
end