% want to implement a discrete channel with binary input alphabet and 
% arbitrary number of output alphabet, given the conditional probability
% table as a .mat file 
% 
% 
function output = dmc( input, cond_prob_filename )
% 
% input is a binary matrix (each realization is a 2d vector) 
% output is a vector having the channel outcome for each realization
% 
% for t = 2 users, four equiprobable case may happen:
% 00, 01, 10, 11
% 
% load channel conditional probability matrix 
load(cond_prob_filename); 
% 
% for each row of this table, we construct a partition set on the interval
% [0 1] and then drag uniformly from this interval. Depending on the
% specific partition our number fall within, we decide which output symbol
% has occured 
% 
n_in = size(qt_probs,1);
n_out = size(qt_probs,2);
partition_set = zeros(n_in, n_out + 1);
% convert 2d messages to scalars and use them to refer to rows of our table
input = sum(bsxfun(@times, input, [2 1]'), 1) + 1;  
% 
q = [zeros(n_in, 1), qt_probs];
%
for i = 1:n_in
    tmp = toeplitz(q(i,:)); 
    tmp = triu(tmp);
    % now sum on each column of this matrix
    partition_set(i,:) = sum(tmp,1);
end

% now create a vector the same length as input filled with random numbers 
% uniformly drawn from [0,1] 
tmp2 = rand(length(input),1);
%
ix = bsxfun(@ge, partition_set(input,:), tmp2);
% finally we have output alphabets of the channel
sum(ix,2);
output = size(partition_set,2) - sum(ix,2);
