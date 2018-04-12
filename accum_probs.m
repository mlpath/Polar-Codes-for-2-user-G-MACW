% this function merge the probabilities in a table based on the indices
% given by a vector 
% 
function merged_table = accum_probs(ix, prob_table)
%
SZ_r = size(prob_table,1);
M = max(ix);
ix = reshape(bsxfun(@plus, repmat(ix', SZ_r, 1), (0:M:(SZ_r - 1)* M)')', size(prob_table,2) * SZ_r, 1);
merged_table = reshape(accumarray(ix, reshape(prob_table',size(prob_table,2) * SZ_r, 1)), M, SZ_r)';
end
