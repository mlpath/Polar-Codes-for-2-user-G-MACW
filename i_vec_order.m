% this function calculates a vector whose elements are in the order
% specified in Abbe paper
% 
%
function ixx = i_vec_order(in_size)
% in_size = 2^4; 
w = 2.^((log2(in_size)-1):-1:0);
%
in_size;
d_vec = bin2double(in_size);
% 
s_vec = sum(d_vec,2);
[dummy,ix] = sort(s_vec,'ascend');
%
ixx = zeros(size(ix));
for i = 1:log2(in_size)
    tmp = ix(dummy == i);
    f_vec = fliplr(d_vec(tmp,:));
    %
    f_val = f_vec * w'; 
    [d_f,ix_f] = sort(f_val,'descend');
    ixx(dummy == i) = tmp(ix_f);
end
%
ixx(2:end) = ixx(2:end) - 1; 

    
    
