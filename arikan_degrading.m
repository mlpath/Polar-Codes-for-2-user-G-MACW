% this is a revised version of e_arikan and e_degrading_merge
% the idea is to combine both functions in one file so that there is no
% need for saving partial parts of files and then loading them one by one
% in e_degrading. It is very likely that we can increase the speed of our
% algorithm with this modification
%
% we are also revising this function having the new control routine in
% mind (which computes all channels at a specific level of binary tree)
% 
% input to this channel are input channel id, fidelity parameter, and block
% size for memory management
% 
% output is the (degrading_merge) version of channels in next level of
% binary tree (i.e. two more .mat files are created)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function arikan_degrading(input_str, fidelity_par, max_block_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tell them where we have written the degraded channel
% input_filename = 'chan_01.mat';
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   inputs required for degrading operation               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize list of labels and list of probabilities
list_b_p = [];
list_b_m = [];
list_p_p = [];
list_p_m = [];
% number of distinct vector messages (i.e. for 2 users with binary
% alphabet we have 2^2 = 4 distinct vectors)
SZ_row = 4;
% these are used for scalarization, i.e. 
% sum{b_i * b_max ^ i}_{i = 0}^{nrow - 1} (like changing binary to decimal)
% we may reconsider this method later, but it's working for now
b_max = 2*ceil(fidelity_par /(log(2)*exp(1)));
ex_b = (0:SZ_row -1)';
% this will change after we rewrite the main control routine
% initialize the output channel where results are stored 
qt_probs = zeros(4,0);
% save output results (i.e. degraded versions of arikan plus and minus
% channels) separately
output_str_p = strcat('C:\MonaHajimomeni\mfiles\debug_mac_v2\debug_channels_14\chan_', input_str, '1.mat');
output_str_m = strcat('C:\MonaHajimomeni\mfiles\debug_mac_v2\debug_channels_14\chan_', input_str, '0.mat');
save( output_str_p, 'qt_probs', '-v7.3');
save( output_str_m, 'qt_probs', '-v7.3');
% we might not need any object as we are saving the entire file at once
% but for now I leave it as it is
out_chan_p = matfile(output_str_p, 'Writable', true);
out_chan_m = matfile(output_str_m, 'Writable', true);
% 
% define our custom function for bsxfun
bxor = @(a,b) bitxor(a,b);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       input to arikan function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty( input_str )
    input_filename = strcat('C:\MonaHajimomeni\mfiles\debug_mac_v2\debug_channels_14\chan_', input_str, '.mat');
else
    input_filename = 'raw_channel_14.mat';
end
% define an object to input file (from previous level in binary tree)
deg_chan = matfile(input_filename);
% qt_probs is the name of variable we created in degrading_merge file
[nrow, ncol] = size(deg_chan,'qt_probs');
% this buffer size must be a multiple of number of rows of probability
% table (i.e. nrow)
% if size(list_b_p,2) > 2e4
%    SZ_buff = 1e5;
% else
%    SZ_buff = 1e6;
% end
SZ_buff = max_block_size; %4e7;
% pointers
ptr2 = 0;
ntot = nrow * ncol * ncol
SZ_col = min(ntot - ptr2, SZ_buff);
% we're assuming that qt_probs has more columns that rows. i.e laying
% 4*ncol*ncol because this is our alphabet size
n_blocks = ceil(ntot/SZ_col);
%
% while ( ptr2 < ntot )
for i = 1:n_blocks
%     [i ntot/SZ_col]
    [i n_blocks]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       block-processing stuff 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each range if indices specified by ptr:ptr2 we need to obtain the
    % respective assignment in table of probability of output channel (i.e.
    % arikan plus channel, later we derive arikan minus from plus)
    %
    ptr = (i - 1) * SZ_buff + 1; 
    % 
    if SZ_col == SZ_buff
        ptr2 = i * SZ_buff;
    else
        ptr2 = ntot;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Convert Indices to Assignment (in output arikan plus channel)      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % first column assignment (i.e. variable u1)
    u1 = rem(ptr:ptr2, 4);
    % correct for zero elements
    ix = u1 == 0;
    u1(ix) = 4;
    % quotient of the division by 4
    q = floor((ptr:ptr2)/4);
    % account for those elements which give a zero remainder
    q(ix) = q(ix) - 1; 
    % assignments for second and third columns
    y1 = rem(q,ncol) + 1;
    y2 = floor(q/ncol) + 1;
    % this tells us which symbols of p(y|u) are needed to calculate this
    % block of arikan plus channel
    %
    % check which columns of the table should be read from file
    [C,~,~] = unique([unique(y1),unique(y2)]);
    % read from file, we have to keep the order for matfile to work
    ix = min(C):max(C);
    tmp = deg_chan.qt_probs(:,ix);
    % keep only those columns of table which have appeared in C. therefore,
    % this is the probability table of all unique symbols
    tmp = tmp(:,ismembc2(C,ix));
    % factor p(y2|u2) can be obtained by copying those columns of tmp where
    % y2 is a member
    B = tmp(:,ismembc2(y2,C));
    % factor p(y2|u1+u2) is obtained by copying those columns of tmp where
    % y1 is a member
    A = tmp(:,ismembc2(y1,C));
    % and then permuting the row orders to account for binary xor operation
    % u3 = u1 + u2
    u3 = bsxfun(bxor, repmat((0:(nrow - 1))',1,4),u1(1:4)-1) + 1;
    % permutation
    for j = 1:4
        % because we can't have unordered indexing with matfile objects
        % now change the order of columns- this is one of the factors
        A(:,j:4:SZ_col) = A(u3(:,j),j:4:SZ_col);
    end
    % we finally have arikan plus channel
    A = A.*B ;
    clear B;% y1 y2 u1 tmp
    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   End of calculation of arikan channels  %%%%%%%%%%%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%   Start the e_degrading_merge function %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %   
    tic
    % pointer to location of degraded tabel of plus channel   
    lp_ptr = length(list_b_p) + 1;
    % pointer to location of degraded label of minus channel
    lm_ptr = length(list_b_m) + 1; 
    % conditional probability of plus channel
    A = A / nrow; 
    % conditional probability of minus channel
    q_m = reshape(sum(A,1), 4,SZ_col / 4);
    % remove all-zero columns
    A(:, ~any(A,1)) = [];
    q_m(:, ~any(q_m,1)) = [];
    % now do the rest
    %
    % posterior of plus channel- note how we calculate marginal by sum
    post_prob_p = A./repmat(sum(A,1),SZ_row,1); 
    flag_nan_p = sum(sum(isnan(post_prob_p)));
    if flag_nan_p
        fprintf('NaN occured: %d',flag_nan_p);
    end
    % this is the output of binning, rewrite posterior variables
    post_prob_p = e_binning(post_prob_p, fidelity_par) - 1;
    % make them scalar
    post_prob_p = bsxfun(@times, post_prob_p, b_max.^ex_b);
    % update list of existing labels for plus channel
    list_b_p( lp_ptr:(size(A, 2) + lp_ptr - 1),1) = sum(post_prob_p , 1);
    clear post_prob_p; 
    % sort labels
    [list_b_p, ix] = sort(list_b_p,'ascend');
    % update list of corresponding probabilities        
    list_p_p( : , lp_ptr:(size(A, 2) + lp_ptr - 1)) = A;
    list_p_p = list_p_p(:, ix);
    % distinct labels in plus channel
    [~, ia, ic] = unique(list_b_p);
    list_b_p = list_b_p(ia,:);
    % merge (i.e. add up) the probabilities of identical labels in the list
    ix = reshape(bsxfun(@plus, repmat(ic', SZ_row, 1), (0:max(ic):(SZ_row - 1)* max(ic))')', size(list_p_p,2) * SZ_row, 1);
    list_p_p = reshape(accumarray(ix, reshape(list_p_p',size(list_p_p,2) * SZ_row, 1)), length(ia), SZ_row)';
    post_prob_m = q_m./repmat(sum(q_m,1),SZ_row,1);
    flag_nan_m = sum(sum(isnan(post_prob_m)));
    if flag_nan_m
        fprintf('NaN occured: %d',flag_nan_m);
    end   
    %
    post_prob_m = e_binning(post_prob_m, fidelity_par) - 1;
    % 
    post_prob_m = bsxfun(@times, post_prob_m, b_max.^ex_b);
    % update list of existing labels for minus channel
    list_b_m( lm_ptr:(size(q_m, 2) + lm_ptr - 1),1) = sum(post_prob_m , 1);
    clear post_prob_m;
    [list_b_m, ix] = sort(list_b_m,'ascend');
    list_p_m( : , lm_ptr:(size(q_m, 2) + lm_ptr - 1)) = q_m;
    list_p_m = list_p_m(:, ix);
    %
    clear q_m
    % repeat the same procedure for minus channel
    [~, ia, ic] = unique(list_b_m);
    list_b_m = list_b_m(ia,:);
    % now merge
    ix = reshape(bsxfun(@plus, repmat(ic', SZ_row, 1), (0:max(ic):(SZ_row - 1)* max(ic))')', size(list_p_m,2) * SZ_row, 1);
    list_p_m = reshape(accumarray(ix, reshape(list_p_m',size(list_p_m,2) * SZ_row, 1)), length(ia), SZ_row)';
    % seems like we are done
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  End of degrading_merge operation for plus and minus arikan channels
%%%
%%%%%%%%%%%%%%%%%%%%%%%  Back to arikan file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if size(list_b_p,2) > 2e4
%         SZ_buff = 1e5;
%     else
%         SZ_buff = 1e6;
%     end
        SZ_col = min(ntot - ptr2, SZ_buff);
%         i = i + 1;
    toc
end
[size(list_p_p), size(list_p_m)]
% [tot1;sum(list_p_p,2)';zeros(1,0);tot2;sum(list_p_m,2)']

% at the end we have two arrays of probability tables, list_p_p & list_p_m
out_chan_p.qt_probs(:,1:size(list_p_p,2)) = list_p_p;
out_chan_m.qt_probs(:,1:size(list_p_m,2)) = list_p_m;