function y = polynomial_fit(x, deg, coeff)
flag = 0; 
% check if x is a clumn vector
if size(x, 1) < size(x, 2)
    x = x';
    flag = 1;
end
y = bsxfun(@power, x, deg:-1:0) ;
% coeff is a row vector
y = sum(bsxfun(@times, y, coeff), 2);
% returns a column vector
if flag
    y = y';
end
end
