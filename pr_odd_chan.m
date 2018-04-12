function channel_probs = pr_odd_chan(in_seq)
%
% ok- we think that just combining the computed ratio of previous column
% together is enough for calculation of new column (i.e. without separately
% obtaining y's and u's at each node). So, I'm writing based on this idea.
% 
% in_seq argument is a matrix (rows corresponding to different hypothesis
% and columns for different branches of iterative tree)
%
% the length of in_seq
l_in = size(in_seq,2);
% our next level probabilities have the half of this length 
% this keeps the value of 4 parallel tree to compute probabilities
channel_probs = zeros(4,l_in/2); 
%
% let leave this preprocessing intact (we don't have anu divisions, so
% these are not really important, will comment them)
% ix = in_seq == 0;
% in_seq(ix) = eps;
% %
% ix = isinf(in_seq);
% in_seq(ix) = 1e10;
%
% now implement the iterative formula for individual channels 
% also remember that we are computing all hypothesis in one function
%
p = [1 2 3 4;
     2 1 4 3; 
     3 4 1 2;
     4 3 2 1];
% 
for i = 1:(l_in/2)
    % find the index of respective input sequence
    j = 2*i - 1;
    % compute the probabilities on all 4 trees
    for k = 1:4
        val_1 = in_seq(:,j+1);
        val_2 = in_seq(p(k,:),j); 
        % ith probability in kth tree
        channel_probs(k,i) =   (sum(val_1 .* val_2)/4); 
    end
    %
    % also leave this preprocessing untouched
    %
    % this won't happen too
%     if isinf(channel_probs(k,i))
%         channel_probs(k,i) = 1e10;
%     end
end

    