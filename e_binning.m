function b = e_binning( x , f_par)
% the vectorized (i.e. more efficient) version of binning function as 
% defined in Equation (14) of Tal Vardy approximate MAC polar construction
% mu = 512;
alpha = 1/exp(1);
beta = 1/(log(2)*exp(1));
%
% this function accepts a vector of posterior probabilities and finds the 
% corresponding vector of super bins
% find the transformation eta
eta_x = -x.*log2(x);
eta_x(isnan(eta_x)) = 0; 
IX = find(eta_x < 0, 1);
if ~isempty( IX ) 
    disp('negative.........');
end
% find the upper rounded (i.e. quantized) levels
b = floor( 1 + f_par * eta_x ); 
% find those elements which are lower than alpha
% those which are equal to alpha,
b(x == alpha) = ceil( beta * f_par ); 
% and those which are higher than alpha
b( x > alpha ) = 2*ceil(beta * f_par) + 1 - b( x > alpha);
%
% note that the length of input and output of this function is limited by
% the size of sub-block we have in degrading_merge, so that there shouldn't
% be any need for more chunking



