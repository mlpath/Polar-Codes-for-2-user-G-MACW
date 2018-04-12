%
% later this will be a function
%
% clc;
% clear all;
% close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [estimated_symbols, likelihood_vec] = sc_decoder(received_sig, frozen_ix, frozen_data, chan_type,chan_par, nlike)
%
% want to modify it so that we can tell how many of likelihoods it's going
% to compute
N = length(received_sig);%16;
%
n = log2(N) + 1;
u = zeros(1,N);
L = -1*ones(N,n);
%
count = 1;
ix_now = 1;
N_now = 1;
%
% this is the received signal samples (1 to N)
%
y = received_sig;
%
% for awgn
% y = y + sqrt(sigma2) * randn(size(y));
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is where I modify sc_decoder, instead of letting count go up to N,
% it should be working up to nlike
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% while (count <= N)
while (count <= nlike)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % it starts with computation..
    while (N_now <= N)
        % works for ix_now = 1, N_now = 1
        row = find(bit_reversed(N_now) == ix_now):N_now:N;
        %
        if (N_now == 1)
            if ~strcmp(chan_type, 'marg_gmac')
                % our old code for point-point case
                L(row,n - log2(N_now)) = lr_chan_level(chan_type,chan_par,y);% zeros(size(L(row,n - log2(N_now))));%
            else
                % when we are decoding mac signals (individually), we just
                % pass computed likelihood ratios for each user and at
                % respective receivers 
                L(row,n - log2(N_now)) = chan_par; 
            end
            ix_now = 2*ix_now - 1;
        else
            % value of count indicates next symbols being estimated..
            % so we only have them up to u(count-1)...
            %
            if mod(ix_now,2)
                % this is an upper left corner of butterfly, formula for odd indices
                L(row,n - log2(N_now)) = lr_sc_odd(L(row(1):N_now/2:N,n - log2(N_now)+1));% zeros(size(L(row,n-log2(N_now))));% 
                % stupid correction
                ix_last = ix_now;
                ix_now = 2*ix_now - 1;
%                 fprintf('\n 1 N_now = %d, ix_now = %d,count = %d \n',N_now,ix_last,count)
            else
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % figuring out the input to this function..it's last input
                % from previous stage (i.e. right-hand column)
                % should make it on the fly...
                %
                % ok, we already know which row we're on..this is the root
                tmp = bit_reversed(N);
                % for all nodes in this step..
                in_array = tmp(row) - 1;
                % 
                u_even = zeros(size(row));
                for i = 1:length(row)     
                    tmp = u(1:in_array(i)); 
                    r = row(i);
                    for j = 1:log2(N/N_now)
                        % how far away we are from the root?
                        tmp_o = tmp(1:2:end); 
                        tmp_e = tmp(2:2:end);
                        %
                        q = ((r - 1) - mod(r-1,N/(2^j)))/(N/2^j);
                        %
                        if mod(q,2)
                            tmp = tmp_e;
                        else 
                            tmp = xor(tmp_e,tmp_o);
                        end
                    end
                    % ..last element of this vector
                    u_even(i) = tmp(end);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                u_even;
                %
                % this is a lower left corner of butterfly, formula for even indices
                %
                L(row,n - log2(N_now)) = lr_sc_even(L((row(1)-N_now/2):N_now/2:N, n - log2(N_now) + 1), u_even);% zeros(size(L(row,n - log2(N_now))));%
                % stupid correction
                ix_last = ix_now;
                ix_now = 2*ix_now - 1;
%                 fprintf('\n 2 N_now = %d, ix_now = %d, count = %d \n',N_now,ix_last,count)
            end
        end
        N_now = 2*N_now;
    end

%
% because we don't have more butterfly stages in our graph!
ix_now = ix_last;
N_now = N;
% 
% updating estimated symbols
if (N_now == N)
    % store the value of estimated symbols
    ix = find(frozen_ix == count, 1);
    if ~isempty(ix)
% count
        %         find(frozen_ix == count, 1)
        % this index is frozen, discard estimated likelihood as we know it
        u(count) = frozen_data(ix);
%         fprintf('\n hi \n');
    else
        if L(row, n- log2(N_now)) > 1
            u(count) = 0;
        else
            u(count) = 1;
        end
    end
    
    % tell them we are done with updating this symbol..
    %
    count = count + 1;
end
% reset this flag
done = 1;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the second place I alter sc_decoder. We require that it only
% compute symbols and likelihoods up to nlike, whereas up to N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% while (done)  && (count <= N)
while (done) && (count <= nlike)
    % now this is node updating toward channel
    if mod(ix_now,2)
        ix_now = ix_now + 1;
        % if we are on the upper corner of a butterfly
    else
        % if we are on the lower corner of a butterfly
        ix_now = ix_now/2;
        N_now = N_now/2;
    end
    % so it seems we updated our coordinates..
    %
    N_now;
    ix_now;
    row = find(bit_reversed(N_now) == ix_now) : N_now : N;
    if L(row(1), n - log2(N_now)) == -1
        done = 0;
        % another stupid correction
        ix_last = ix_now;
    end
end
% fprintf('\n 3 N_now = %d, ix_now = %d, count = %d \n',N_now,ix_now,count);
end
estimated_symbols = u; 
% we only care about symbol level llikelihoods at the end
likelihood_vec = L(:,1);




