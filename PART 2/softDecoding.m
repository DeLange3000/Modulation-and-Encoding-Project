function [dec_bs] = softDecoding(H,rec_bs, N0)

[M, N] = size(H);
rec_bs_l = length(rec_bs);
var = N0/2;
dec_bs = [];

    for block_i = 1:N:rec_bs_l
        block = rec_bs(block_i:block_i+N-1);
        H_index = H.*(1:M)';
        L_q = zeros(M, N);
        L_r = zeros(M, N);
        %make L_q
        for j = 1:N
            L_q(nonzeros(H_index(:,j)),j) = -2*block(j)/var; %create L matrix
        end
        
        % loop from here 10 times

        %make L_r
        for j = 1:N
            for a = nonzeros(H_index(:,j))' %loops over non zero H indices
            y_index = 1:N;
            y_index = y_index(y_index~=j); %j itself cannot be included
                if(isempty(nonzeros(L_q(a, y_index)))) %if only zero elements are present, make L_r zero
                    L_r(a,j) = 0;
                else
                    L_r(a,j) = prod(sign(nonzeros(L_q(a, y_index))), 'all')*min(nonzeros(abs(L_q(a, y_index))), [], 'all');  %create L matrix
                end
            end
        end
        
        %add L_r to L_q
        for j = 1:N
            for a = nonzeros(H_index(:, j))'
                x_index = 1:M;
                x_index = x_index(x_index~=a); %j itself cannot be included
                L_q(a,j) =  L_q(a,j) + sum(L_r(x_index,j));
            end
        end
        
        % stop loop 10 times here

        for j = M+1:N
            L_Q = -2*block(j)/var + sum(L_r(:,j));
            if(L_Q < 0)
                dec_bs = [dec_bs 1];
            else
                dec_bs = [dec_bs 0];
            end
        end
        
       

    end


end