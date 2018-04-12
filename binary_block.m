% this function construct the binary form of the codewords with non-binary
% symbols. therefore, a (1 x m1) vector having symbols in GF(2^R) will become
% a (R x m1) binary block
%
function bin_codeword = binary_block(input_codeword, alph_size)
% alph_size is the size of alphabet the code symbols belong to in bits
code_width = size(input_codeword, 2); 
bin_codeword = zeros(alph_size, code_width); 
% 
tmp = dec2bin(input_codeword, alph_size)'
for j = 1:code_width
    % this is the binary equivalent of our coded symbols
    % because the size of this 'for' loop is not very large
    bin_codeword(:,j) = str2num(tmp(:,j));
end