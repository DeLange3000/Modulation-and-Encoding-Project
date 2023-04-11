function [dec_bs] = softDecoding(H, rec_bs, sigma2_w)
    [M, N] = size(H);
    rec_bs_l = length(rec_bs);
    
    dec_bs = [];

    for block_i = 1:N:rec_bs_l
        block = rec_bs(block_i:block_i+N-1);

        v_nodes = zeros(M+1, N);
        v_nodes(1, :) = block;

		v_nodes = block;

		r = block;
		q0 = zeros(M,N);
		for i = 1:N
			q0(:,i) = 1/(1+exp(2*r(i)/sigma2_w));
		end
		q1 = 1 - q0;
		rmat0 = zeros(M, N);
		for i = 1:M
			for j = 1:N
				prod = 1;
				for p = 1:N
					if p ~= j
						prod = prod*(1-2*q1(i,p));
					end
				end
				rmat0(i,j) = 0.5*(1+prod);
			end
		end

		rmat1 = 1 - rmat0;

        for i = 1:M
            % #rogier is genius
            arr = nonzeros(H(i,:) .* (1:N))';

            for j = arr
                arr_sel = arr(arr~=j);
                selected = block(arr_sel);

				q

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