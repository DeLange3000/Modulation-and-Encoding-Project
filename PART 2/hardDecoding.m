function [dec_bs] = hardDecoding(H, rec_bs, number_of_iterations)
    [M, N] = size(H);
    rec_bs_l = length(rec_bs);
    
    dec_bs = [];

    for block_i = 1:N:rec_bs_l
        for a = 1:number_of_iterations
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
    
            for i = 1:N
            H_nonzero = H(:,i)'.*(2:M+1); % add index due to offset in v_nodes
            H_nonzero = nonzeros(H_nonzero);
            cor_block = sum(v_nodes([1; H_nonzero],i));
            cor_block = round(cor_block/(length(H_nonzero)+1));
            block(i) = cor_block;
            end
            if(block*H' == 0)
                break
            end
        end
        dec_bs = [dec_bs block(M+1:N)];
    end
    
end