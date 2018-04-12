% this is a simple function to calculate bit-reversed vector of input
function r_index = bit_reversed(in_size)
%
% R = zeros(in_size);
% if in_size == 1
%     R = 1;
% else
%     for i = 1:(in_size/2)
%         R(i,2*i-1) = 1;
%         R(i+(in_size/2),2*i) = 1;
%     end
% end
% % 
%
if in_size == 1
    r_index = in_size;
else
    index = 0:(in_size - 1);
    binary_val = dec2bin(index,log2(in_size));
    reversed_val = fliplr(binary_val);
    r_index = 1 + bin2dec(reversed_val);
end
