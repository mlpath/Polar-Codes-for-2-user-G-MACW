% this function calculates a vector whose elements are in the order
% specified in Abbe paper
% 
function d_vec = bin2double(in_size)
%
% 
in_vec = 0:in_size - 1;
%
b_vec = dec2bin(in_vec);
d_vec = zeros(in_size,log2(in_size));
%
for j = 1:in_size
for i = 1:log2(in_size)
    if strcmp(b_vec(j,i),'1')
        d_vec(j,i) = 1; 
    else 
        d_vec(j,i) = 0; 
    end
end
end
% 
