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
function arikan_degrading_v2(input_str, fidelity_par, max_block_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tell them where we have written the degraded channel
% input_filename = 'chan_01.mat';
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   inputs required for degrading operation               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
qt_probs = zeros(SZ_row,0);
% save output results (i.e. degraded versions of arikan plus and minus
% channels) separately %%%%%%%%%%%%%%%%%%%%%%%%%%  Main OUTPUT Files %%%%%%
output_str_p = strcat('D:\April-queens\mfile\backup_final\debug_qpsk_4\chan_', input_str, '1.mat');
output_str_m = strcat('D:\April-queens\mfile\backup_final\debug_qpsk_4\chan_', input_str, '0.mat');
save( output_str_p, 'qt_probs', '-v7.3');
save( output_str_m, 'qt_probs', '-v7.3');
%
out_chan_p = matfile(output_str_p, 'Writable', true);
out_chan_m = matfile(output_str_m, 'Writable', true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temporary files, we create, delete, and recreate them
tmp_str_p = 'D:\April-queens\mfile\backup_final\debug_qpsk_4\temp_chan_p.mat';
tmp_str_m = 'D:\April-queens\mfile\backup_final\debug_qpsk_4\temp_chan_m.mat';
save(tmp_str_p, 'qt_probs', '-v7.3');
save(tmp_str_m, 'qt_probs', '-v7.3');
%
tmp_chan_p = matfile(tmp_str_p, 'Writable', true);
tmp_chan_m = matfile(tmp_str_m, 'Writable', true);
% initialize our big file, how much time it takes?
tic
tmp_chan_p.qt_probs(:,b_max^SZ_row) = zeros(4,1);
tmp_chan_m.qt_probs(:,b_max^SZ_row) = zeros(4,1);
toc
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define our custom function for bsxfun
bxor = @(a,b) bitxor(a,b);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       input to arikan function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty( input_str )
    input_filename = strcat('D:\April-queens\mfile\backup_final\debug_qpsk_4\chan_', input_str, '.mat');
else
    input_filename = 'raw_channel_qpsk_2_mu30.mat';
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
%
% this is a loop over blocks of newly constructed arikan channel
%
% I need a loop for eliminating zero columns from previously constructed
% degraded channel
%
% I am worried about the reading part from file, which is slow and is going
% to be repeated two times, how can I do it just once?
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
    u1 = rem(ptr:ptr2, SZ_row);
    % correct for zero elements
    ix = u1 == 0;
    u1(ix) = SZ_row;
    % quotient of the division by 4
    q = floor((ptr:ptr2)/SZ_row);
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % actually we need to find min(C)-th nonzero element to max(C)-th
    % nonzero element, but how?
    %
    % I rather not change arikan, it's working good enough
    tmp = deg_chan.qt_probs(:,ix);
    %
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
    u3 = bsxfun(bxor, repmat((0:(nrow - 1))',1,SZ_row),u1(1:SZ_row)-1) + 1;
    % permutation
    for j = 1:SZ_row
        % because we can't have unordered indexing with matfile objects
        % now change the order of columns- this is one of the factors
        A(:,j:SZ_row:SZ_col) = A(u3(:,j),j:SZ_row:SZ_col);
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
    % conditional probability of plus channel
    A = A / nrow;
    % conditional probability of minus channel
    q_m = reshape(sum(A,1), SZ_row,SZ_col / SZ_row);
    % remove all-zero columns
    A(:, ~any(A,1)) = [];
    q_m(:, ~any(q_m,1)) = [];
    %
    % now do the rest
    %
    % posterior of plus channel; note how we calculate marginal by sum %%%%
    post_prob_p = A./repmat(sum(A,1),SZ_row,1);
    flag_nan_p = sum(sum(isnan(post_prob_p)));
    if flag_nan_p
        fprintf('NaN occured: %d',flag_nan_p);
    end
    % this is the output of binning, rewrite posterior variables
    post_prob_p = e_binning(post_prob_p, fidelity_par) - 1;
    % make them scalar
    post_prob_p = bsxfun(@times, post_prob_p, b_max.^ex_b);
    % find unique labels in this block
    [u_list_p, ~, ic_p] = unique(sum(post_prob_p, 1));
    clear post_prob_p
    %
    % minus channel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    post_prob_m = q_m./repmat(sum(q_m,1),SZ_row,1);
    flag_nan_m = sum(sum(isnan(post_prob_m)));
    if flag_nan_m
        fprintf('NaN occured: %d',flag_nan_m);
    end
    % superbins (or labels)
    post_prob_m = e_binning(post_prob_m, fidelity_par) - 1;
    %
    post_prob_m = bsxfun(@times, post_prob_m, b_max.^ex_b);
    % find unique labels in this block of minus channel
    [u_list_m, ~, ic_m] = unique(sum(post_prob_m, 1));
    clear post_prob_m
    % update the temporary probability list according to unique indices
    temp_prob_p = accum_probs(ic_p, A);
    % ok, so we have a number of distinct labels
    temp_prob_m = accum_probs(ic_m, q_m);
    %
    lp_ptr = 1;
    lm_ptr = 1;
    SZ_list = 1e4;
    % until we haven't saved all labels in our list (for p and m channels)
    while ( lp_ptr <= u_list_p(end) || lm_ptr <= u_list_m(end) )
        % pointer to end of this block
        lp_ptr2 = lp_ptr + SZ_list - 1;
        lm_ptr2 = lm_ptr + SZ_list - 1;
        %
        if lp_ptr2 > u_list_p(end)
            lp_ptr2 = u_list_p(end);
        end
        %
        if lm_ptr2 > u_list_m(end)
            lm_ptr2 = u_list_m(end);
        end
        %
        ix_p = u_list_p <= lp_ptr2 & u_list_p >= lp_ptr;
        ix_m = u_list_m <= lm_ptr2 & u_list_m >= lm_ptr;
        % 
        if sum(ix_p) && (lp_ptr <= u_list_p(end))
            % read this block of file 
            tmp = tmp_chan_p.qt_probs( :, lp_ptr:lp_ptr2);
            % this is now a sum where all terms have equal size, update
            tmp(:, u_list_p(ix_p) - lp_ptr + 1) = tmp(:, u_list_p(ix_p) - lp_ptr + 1) + temp_prob_p(:, ix_p);
            % and now write it back
            tmp_chan_p.qt_probs( :, lp_ptr:lp_ptr2) = tmp; 
            % update our pointer
            lp_ptr = lp_ptr2 + 1;
        else
            % nothing to write
            lp_ptr = lp_ptr2 + 1;
        end
        %
        if sum(ix_m) && (lm_ptr <= u_list_m(end))
            % read from the file
            tmp = tmp_chan_m.qt_probs(:, lm_ptr:lm_ptr2);
            % update
            tmp(:, u_list_m(ix_m) - lm_ptr + 1) = tmp(:, u_list_m(ix_m) - lm_ptr + 1) + temp_prob_m(:, ix_m);
            % and write it back
            tmp_chan_m.qt_probs( :, lm_ptr:lm_ptr2) = tmp;
            % update our pointer
            lm_ptr = lm_ptr2 + 1;
        else
            % nothing to write
            lm_ptr = lm_ptr2 + 1;
        end
    end
    clear temp_prob_m temp_prob_p
    % don't forget to filter out possibly many zero columns from this long
    % lists which we saved!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  End of degrading_merge operation for plus and minus arikan channels
%%%
%%%%%%%%%%%%%%%%%%%%%%%  Back to arikan file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    SZ_col = min(ntot - ptr2, SZ_buff);
    %
    toc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ok, now we are done with main body of arikan and degrading_merge, we just
% need to do some polish and remove the many number of zero columns
%
%%%%%%%%%%%%%%%%% Remove the zero columns from this file %%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
tic
% reading pointers
lp_ptr = 1;
lm_ptr = 1;
% writing pointers
lP_ptr = 1;
lM_ptr = 1;
% size of current arrays
[~,SZ_p] = size(tmp_chan_p, 'qt_probs');
[~,SZ_m] = size(tmp_chan_m, 'qt_probs');
%
% until we aren't done nevigating through these files
while ( lp_ptr <= SZ_p || lm_ptr <= SZ_m )
    %
    lp_ptr2 = lp_ptr + SZ_list - 1;
    lm_ptr2 = lm_ptr + SZ_list - 1;
    if lp_ptr2 > SZ_p
        lp_ptr2 = SZ_p;
    end
    if lm_ptr2 > SZ_m
        lm_ptr2 = SZ_m;
    end
    %
    if lp_ptr <= SZ_p
        % read array from temporary file
        tmp_p = tmp_chan_p.qt_probs(:,lp_ptr:lp_ptr2);
        % remove zero columns
        tmp_p(:, ~any(tmp_p, 1)) = [];
        % write it back if it's not empty
        if ~isempty(tmp_p)
            % if there are any nonzero columns, write them to main file
            out_chan_p.qt_probs(:, lP_ptr:(lP_ptr + size(tmp_p, 2) - 1)) = tmp_p;
        end
        % update writing pointers
        lP_ptr = lP_ptr + size(tmp_p, 2);
    end
    %
    if lm_ptr <= SZ_m
        % read array from temporary file
        tmp_m = tmp_chan_m.qt_probs(:,lm_ptr:lm_ptr2);
        % remove zero columns
        tmp_m(:, ~any(tmp_m, 1)) = [];
        % write it back
        if ~isempty(tmp_m)
            % if there are any nonzero columns, write them to main
            out_chan_m.qt_probs(:, lM_ptr:(lM_ptr + size(tmp_m, 2) - 1)) = tmp_m;
        end
        % update writing pointers
        lM_ptr = lM_ptr + size(tmp_m, 2);
    end
    % update reading pointers
    lp_ptr = lp_ptr2 + 1;
    lm_ptr = lm_ptr2 + 1;
    %
end
%
% now make sure that we got different sizes at the end
[SZ_p, size(out_chan_p, 'qt_probs')]
[SZ_m, size(out_chan_m, 'qt_probs')]
% delete temporary files so that we can use them next time
delete(tmp_str_p,tmp_str_m);
% [tot1;sum(list_p_p,2)';zeros(1,0);tot2;sum(list_p_m,2)']
tic
