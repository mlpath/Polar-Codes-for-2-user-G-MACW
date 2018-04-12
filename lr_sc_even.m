function likelihood = lr_sc_even(input,odd_symbol)
% ok- we think that just combining the computed ratio of previous column
% together is enough for calculation of new column (i.e. without separately
% obtaining y's and u's at each node). So, I'm writing based on this idea.
% 
% 
% because this is the even function..
%
% the length of input
l_in = length(input);
% our next level lr's have the half of this length 
likelihood = zeros(1,l_in/2); 
%
%
ix = input == 0;
input(ix) = eps;
%
ix = isinf(input);
input(ix) = 1e10;
%
for i = 1:(l_in/2)
    j = 2*i - 1;
    likelihood(i) = input(j+1) * (input(j)^(1-2*odd_symbol(i))); 
    % 
    if isinf(likelihood(i))
        likelihood(i) = 1e10;
    end
end


    