function channel_probs = pr_even_chan(in_seq,odd_symbol)
% ok- we think that just combining the computed ratio of previous column
% together is enough for calculation of new column (i.e. without separately
% obtaining y's and u's at each node). So, I'm writing based on this idea.
% 
% 
% because this is the even function..
%
% the length of in_seq
l_in = size(in_seq,2);
% our next level probabilities have the half of this length 
% this keeps the value of 4 parallel tree to compute probabilities
channel_probs = zeros(4,l_in/2); 
%
% some unnecessary preprocessing, as we don't have any division
% ix = in_seq == 0;
% in_seq(ix) = eps;
% %
% ix = isinf(in_seq);
% in_seq(ix) = 1e10;
%
c = 0;
size(channel_probs);
in_seq;
for i = 1:(l_in/2)
    j = 2*i - 1;
    for k = 1:4
        % odd_symbol must be between 0 and 3- this gives corresponding row
        % in in_seq matrix 
        idx = bitxor(k-1, odd_symbol(i)) + 1;
        channel_probs(k,i) = (in_seq(idx, j) * in_seq(k, j+1)/4);
    end
    % 
%     if isinf(channel_probs(k,i))
%         channel_probs(k,i) = 1e10;
%     end
end
c;


    