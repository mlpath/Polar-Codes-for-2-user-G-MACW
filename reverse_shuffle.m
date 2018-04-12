% this is a simple function to calculate reverse shuffle matrix
function R = reverse_shuffle(in_size)
%
R = zeros(in_size);
if in_size == 1
    R = 1;
else
    for i = 1:(in_size/2)
        R(i,2*i-1) = 1;
        R(i+(in_size/2),2*i) = 1;
    end
end
% 
R = R';
% 


