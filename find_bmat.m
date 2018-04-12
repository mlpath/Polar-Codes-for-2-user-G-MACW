% this function finds a binary representation for the given rate_vector,
% and then checks if binary matrix and rate_vector are matching. If so, it
% generates a flag confirming that rate_vector is of a binary matroid,
% otherwise the flag will be 0. 
% 
% we also find a matrix including all the bases that can be derived from
% this rate_vector. We do not need to use basis enumeration algorithm
% because we already have the rank of all possible subsets in the input
% vector. Note also that all bases have equal cardinality, so we can easily
% store all of them in a matrix
%
% 
function  [B, bases, is_binary] = find_bmat( rate_vector )
% rate_vector
% find B matrix
% clc
% clear all
% % a wrong rank function
% rate_vector = [0 1 1 1 2 1 1 3];
% % a correct one
% rate_vector = [0 1 1 1 2 2 1 2];
% clc
% clear all
% close all
% 
% rate_vector = [0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3];
% rate_vector = [0, 1, 1, 1, 1, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Find a binary matrix which match with the given rank function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
is_binary = 0;
t = log2(size(rate_vector,2));
ncol = 2^t;
%
order = i_vec_order(2^t);
order_bin = fliplr(dec2bin(order, t));
%
f_rank = zeros(1,2^t);
ptr = 1;
for i = 0:t
    binom = i*ones(1,nchoosek(t,i));
    ptr2 = ptr + length(binom) - 1; 
    f_rank(ptr:ptr2) = binom; 
    ptr = ptr2 + 1;
end
% preallocation (assuming that it is full rank)
B = zeros(t, t);
% the current rate vector (i.e. only one row of v_R)
tmp =  rate_vector;
% now choose a basis of B (i.e. the index of a set with maximal rank)
J = find(order_bin(find(tmp == tmp(ncol), 1),:) == '1');
J = sort(J);
% those who aren't in our basis
Jc = setdiff(1:t, J);
% standard form: J is one basis of the matrix, so we can fill these columns
% with an identity matrix of the same size as maximum rank in rate_vector
B(1:length(J),J) = eye(length(J));
%
% now we need to fill in the other columns of the matrix
for j = 1:length(Jc)
    % find all positions that have anything except j cup J
    left_cols = setdiff(1:t, union(Jc(j),J));
    if isempty(left_cols)
        ix = 1:2^t;
    else
        ix =  order_bin(:,left_cols) == '1';
        ix = sum(ix, 2) == 0;
    end
    %
    % find the first place that we have a less than full rank position
    ix_bin = find(ix, find(tmp(ix) < f_rank(ix) , 1));
    c_set = order_bin(ix_bin(end) , :);
    % what are the elements in this set?
    c_set = find(c_set == '1');
    c_set = setdiff(c_set, Jc(j));
    B(:,Jc(j)) = mod(sum(B(:,c_set),2),2);
end
% B = B(1:tmp(ncol), :); 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               validate the result by recomputing rate_vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% in fact, this wasn't necessary if we had identified if the integer 
% rate_vector is of a binary matroid or not. Our logic is if we
% find one with the same rank, that means that the rate vector we started
% from has a binary representation and therefore is binary.
%
% find the rank function of the constructed matrix
rate_vector = zeros(1,2^t);
%
for i = 1:2^t
    J = find(order_bin(i,:) == '1');
    if ~isempty(J)
        % find the rank of this submatrix of B
        [~, r] = b_gauss_jordan(B(:,J));
        rate_vector(i) = r; 
    else
        rate_vector(i) = 0;
    end
end
% tell them is our B matches with the given (rounded) information vector
if isequal(tmp, rate_vector)
    is_binary = 1; 
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               extract all bases of the this matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% all sets with maximal rank
all_sets = tmp == tmp(ncol);
all_sets = order_bin( all_sets, : );
done = 0;
i = 1;
j = 1;
% store them in a matrix 
bases = zeros(size(all_sets, 1), tmp(ncol));
%
while (~done) && (i <= size(all_sets, 1)) 
    size(all_sets);
    ix = all_sets(i,:) == '1';
    num = sum(ix);
    if num == tmp(ncol)
        bases(j,:) = sort( find(all_sets(i,:) == '1') );
        j = j + 1;
    elseif num > tmp(ncol)
        done = 1;
    end 
    i = i + 1;
end
bases(j:end,:) = [];
% 
% 
