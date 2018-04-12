% our last idea for MAC implementation is to have four parallel tree for
% computation of each of W_i(.|00), W_i(.|01), W(.|10), W(.|11) 
% 
function [estimated_symbols, estimated_bins] = mac_decoder(received_sig, chan_id, frozen_data, frozen_map, nlike)
% here I want to write the control body of these 4 graphs
t = size(frozen_map, 1);
c1 = 0; 
c2 = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(received_sig);%16;
estimated_bins = zeros(2,N);
%
n = log2(N) + 1;
u = zeros(1,N);
L1 = -1*ones(N,n);
L2 = L1; 
L3 = L1; 
L4 = L1; 
%
count = 1;
ix_now = 1;
N_now = 1;
%
% this is the received signal samples (1 to N)
y = received_sig;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nlike is the number of symbols we aim to estimate 
while (count <= nlike)
    % it starts with computation..
    while (N_now <= N)
        % works for ix_now = 1, N_now = 1
        row = find(bit_reversed(N_now) == ix_now):N_now:N;
        %
        if (N_now == 1)
            % I don't mind writing four times, though it's inefficient and
            % not general enough to work with aribitrary number of users
            % 
            % the right thing is probably to save all values in rows and
            % then pick the appropriate one 
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % channel level probabilities for each graph at this time slot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            probs = pr_chan_level(y, chan_id, 'continuous');
            L1(row,n - log2(N_now)) = probs(1,:);% zeros(size(L(row,n - log2(N_now))));%
            L2(row,n - log2(N_now)) = probs(2,:);
            L3(row,n - log2(N_now)) = probs(3,:);
            L4(row,n - log2(N_now)) = probs(4,:);
            %
            ix_now = 2*ix_now - 1;
        else
            % value of count indicates next symbols being estimated..
            % so we only have them up to u(count-1)...
            %
            if mod(ix_now,2)
                % this is an upper left corner of butterfly, formula for odd indices
                L = [L1(row(1):N_now/2:N,n - log2(N_now)+1)';
                     L2(row(1):N_now/2:N,n - log2(N_now)+1)';
                     L3(row(1):N_now/2:N,n - log2(N_now)+1)';
                     L4(row(1):N_now/2:N,n - log2(N_now)+1)'];
                %
                probs = pr_odd_chan(L);
                if ~isempty(find(probs == 0))
                    N_now;
                end
               
                if N_now == N
                    probs;
                    c1 = c1 + 1;
                end
                % 
                L1(row,n - log2(N_now)) = probs(1,:);% zeros(size(L(row,n-log2(N_now))));% 
                L2(row,n - log2(N_now)) = probs(2,:); 
                L3(row,n - log2(N_now)) = probs(3,:); 
                L4(row,n - log2(N_now)) = probs(4,:);
                % stupid correction
                ix_last = ix_now;
                ix_now = 2*ix_now - 1;
%                 fprintf('\n 1 N_now = %d, ix_now = %d,count = %d \n',N_now,ix_last,count)
            else
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % figuring out the input to this specific likelihood 
                % function..it's last input
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
                            % the only reasonable change must be here (in
                            % using bitxor on values of u between 0 and 3)
                            tmp = bitxor(tmp_e,tmp_o);
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
                L = [L1((row(1)-N_now/2):N_now/2:N, n - log2(N_now) + 1)';
                     L2((row(1)-N_now/2):N_now/2:N, n - log2(N_now) + 1)';
                     L3((row(1)-N_now/2):N_now/2:N, n - log2(N_now) + 1)';
                     L4((row(1)-N_now/2):N_now/2:N, n - log2(N_now) + 1)'];
                %
                probs = pr_even_chan(L, u_even);
                %
                if ~isempty(find(probs == 0))
                    N_now;
                end
                if N_now == N
                    probs;
                    c2 = c2 + 1;
                end
                L1(row,n - log2(N_now)) = probs(1,:); % zeros(size(L(row,n - log2(N_now))));%
                L2(row,n - log2(N_now)) = probs(2,:);
                L3(row,n - log2(N_now)) = probs(3,:);
                L4(row,n - log2(N_now)) = probs(4,:);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                        Updating estimated symbols                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (N_now == N)
    % 
    % store all values of channel probabilities at level N in a neat vector
    isequal(L1(row,n-log2(N_now)),0);
    size(L1);
    L = [L1(row, n- log2(N_now))';
         L2(row, n- log2(N_now))';
         L3(row, n- log2(N_now))';
         L4(row, n- log2(N_now))'];
     L';
     if ~isequal(sum(L == 0),0)
         L';
         count;
     end
     if ~isequal(sum(isinf(L)),0)
         L';
         count;
     end
    % rows corresponding to frozen bits in this time slot 
    idx = find(frozen_map(:,count) == 1);
%     % places that are not frozen
%     idx_c = setdiff(1:t,idx);
%     t_prime = length(idx_c);
    % 
    % write it! no fear, no judgement
    % before I look for other problems, lets try decoding without
    % considering the value of frozen bits
    % so, we have to choose which element in L has maximum value, and
    % that's our decoded symbol. Then, what will we do for frozen bits?
    % we keep transmitting but change the frozens to their initial values.
    [~, tmp] = max(L);
%     [tmp, size(L)];
%     if tmp == 5
%         L
%     end
    % decrease row index by one to get decoded symbol- (a binary string) 
    tmp = dec2bin(tmp - 1, t); 
    % now update the value of frozen bits 
    fixed_str = num2str(frozen_data(idx,count)); 
    tmp(idx) = fixed_str; 
    u(count) = bin2dec(tmp); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% now do some preprocessing
    %
    % to find maximum likelihood message we check which of the possible
    % messages are more likely conditioned on available frozen bits
    %
    % if we have at least one non-frozen bit
%     if ~isempty(idx_c)
%         % all possible cases with length(idx_c) non-frozen bits
%         tmp = repmat(dec2bin(0,t), 2^t_prime, 1);
%         if ~isempty( idx )
%             % first convert frozen bits in this time slot to string
%             fixed_str = num2str(frozen_data(idx, count));
%             tmp(:,idx) = repmat(fixed_str',size(tmp,1),1);
%         end
%         % fill remaining places with their binary equivalents
%         tmp(:,idx_c) = dec2bin(0:2^t_prime-1);
%         % this is now the index of rows we want to compare (an integer in the
%         % set {0,1,2,...,2^t-1}
%         idx = bin2dec(tmp) + 1;
%         %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %  finally decode based on (arg)max prob on selected values of L   
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %
%         % we need the index of maximum prob element, tmp below is the index
%         [~,tmp] = max(L(idx)); 
%         % to have integers in {0,1,2,..2^t-1}
%         u(count) = idx(tmp) - 1; 
%     else
%         % if all bits are frozen, we already now the value of u(count)
%         u(count) =  sum(frozen_data(:, count) .* [2 ;1], 1); 
% %         c = c + 1; 
%     end 
    % tell them we are done with updating this symbol..
    count = count + 1;
end
% reset this flag
done = 1;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update the level and node index for the next symbol estimation phase
while (done) && (count <= nlike)
    % find a group of nodes at some level which are not yet updated
    if mod(ix_now,2)
        ix_now = ix_now + 1;
        % if we are on the upper corner of a butterfly
    else
        % if we are on the lower corner of a butterfly
        ix_now = ix_now/2;
        N_now = N_now/2;
    end
    % so it seems we updated our coordinates..
    row = find(bit_reversed(N_now) == ix_now) : N_now : N;
    % I think in all 4 tree the operations happen at the same pace, so if I
    % know that L1 is not updated, it implies the same for L2,3,4 
    %
    if L1(row(1), n - log2(N_now)) == -1
        % for debug
        L2(row(1), n - log2(N_now)) * L3(row(1), n - log2(N_now)) * L4(row(1), n - log2(N_now));
        %
        done = 0;
        % another stupid correction
        ix_last = ix_now;
    end
end
% fprintf('\n 3 N_now = %d, ix_now = %d, count = %d \n',N_now,ix_now,count);
end
% can we return the 0/1 format for this integers, too?
u_bin = dec2bin(u,t);
% because we want LSB bit to be written to first row
% good job! 
for i = 1:t
    estimated_bins(i,:) = str2num(u_bin(:,i)); 
end
%
r = sum(bsxfun(@times, estimated_bins, [2 1]'),1);
[r;u];
size(r); size(u);
isequal( r, u);
estimated_symbols = u; 
c1;
c2;


