function [user_grid, valid_slots, integer_rate_b, integer_rate_e] = find_secure_v23(integer_rates, valid_slots, block_length)
%
t = log2(size(integer_rates,2) + 1);
N = length(valid_slots); 
N_old = block_length;
% preallocate the array that holds transmitting and fixed users
user_idx = zeros(t, length(valid_slots));
Rs = zeros(1,length(valid_slots));
sum_rate = zeros(size(Rs));
% a counter to write in user_idx
lptr = 1;
%
for i = 1:N
    % on each time slot (i.e. for each parallel subchannel)
    integer_rates(i,:);
    [A, bases_1, is_binary_1] = find_bmat([0, integer_rates(i, :)]);
    [B, bases_2, is_binary_2] = find_bmat([0, integer_rates(i + N, :)]);
    %
    if ~(is_binary_1 && is_binary_2)
        valid_slots(i,:) = [];
        integer_rates(i,:) = [];
        integer_rates(i + N,:) = [];
        fprintf('invalid matrix occurred');
        %
    else
        j = 1;
        % initialize secrecy rate with zero (minimum possible)
        rs = 0;
        secure_users = []; 
        %
        n_bases = size(bases_1,1);
        rk = size(bases_1,2);
        %
        %%%%%%%%%%%%%%%%%%%%%% without secrecy constraint, sum rate will be 
        sum_rate(i) = rk;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this is a valid system, find the secure solution
        while (rs < rk) && (j <= n_bases)
            % pick jth basis of A and construct a submatrix of B with it
            % now apply gauss jordan to a submatrix of B
            [rref,~] = b_gauss_jordan(B(:,bases_1(j,:)));
            % find rows where only a single variable is left
            ix = sum(rref, 2) == 1;
            % explicit solutions are in fact non-secure
            non_secure = find(sum(rref(ix,:),1) == 1);
            % if number of secure bits in this basis is lower than stored one
            if rs < (rk - sum(ix))
                secure_users = bases_1(j,setdiff(1:rk, non_secure));
                rs = rk - sum(ix); 
            end
            j = j + 1; 
            [valid_slots(i) rs];
        end
        % now tell us which users should transmit on this slot 
        user_idx(secure_users, lptr ) = 1;
        % sum secrecy rate on valid time slots
        Rs(lptr) = rs; 
        lptr = lptr + 1;
    end
    %
    [valid_slots(i), secure_users];
end
%
length(valid_slots);
lptr - 1;
user_idx(:,(lptr + 1):end) = [];
Rs((lptr + 1):end) = [];
%
user_grid = zeros(t,N_old);
user_grid(:, valid_slots) = user_idx; 
N = length(valid_slots); 
integer_rate_b = integer_rates(1:N,:); 
integer_rate_e = integer_rates(N+1:end,:); 
end