% Ok- now we should start working around the concatenated solution, and
% first you need to review some theory about coding and of course your
% notes about the idea.
% 
% if I don't write anything right now, I wouldn't have anything at the end
% of the day
% 
% block length for polar code
N_p = 256; 
% block length for outer code 
N_o = 4; 
%
% number of users 
t = 2; 
% 
% you know the channel- right? we could use the same ones we used in our 
% previous MAC experiment, i.e. debug_channel_p10 and debug_channel_p14
% however, because the secrecy rate on that settings is not much, I rather
% run the MAC code construction for either a much better Charlie, or a
% better Alice 
%
% what is our inner polar code? it is constructed based on the
% specifications of the mentioned channels. We need to know the 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           secure structure                   
% user_grid, and respective upper/lower bound of block error rate at
% Bob/Eve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
power_id = '10'; 
filename = strcat('secure_struct_', power_id, '.mat');
load(filename); 
% user 1 index and number of secure slots
ix_1 = find(user_grid(1,:) == 1); 
k1 = length(ix_1); 
% user 2 index and number of secure slots
ix_2 = find(user_grid(2,:) == 1); 
k2 = length(ix_2); 
% 
% ok, now we should choose the parameters of outer code 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            outer code parameters (user combined/individual?)                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on the worst channel we wanna use for our secure transmission,
% choose the error correction capability of the outer code
% 
% number of errors able to correct 

% the rate of the code 

% block length (width) 
m = 4; 
% block length (height)- defined the alphabet that coded symbols are coming
% from (for instance r = 2^2 corresponds to QPSK symbols (00, 01, 10, 11)
r = 2^2; 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                raw data generation and outer block encoding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of bits in transmit sequence
nsamp1 = 1e4; 
nsamp2 = 1e5; 
s1 = mod(randi(2, 1, nsamp1), 2); 
s2 = mod(randi(2, 1, nsamp2), 2); 
% get l bits at a time and code them into an r x m block (outer = rs code)
% 
% code for user 1, s1 is the raw data
[s1, x1] = outer_coding(k1, l1, m1, r1, 'rs');
% code for user 2, s2 is the raw data  
[s2, x2] = outer_coding(k2, l2, m2, r2, 'rs');
% generate frozen bits 
F 
% we have m1 inner block for this code, and each block is t * N
X1 = zeros(t * m1, N); 
X2 =
% define n1 and n2, number of coded blocks we pich from user 1 and user 2
n1 = 1; 
n2 = 1; 
for L = 1:n_block
    % suppose we are transmitting # n_block of the final concatenated code,
    % which becomes size(x1,1)/n1 blocks for user 1 and size(x2,1)/n2 
    % blocks for user 2 
    % 
    % index of user 1 coded blocks
    i1 = 
    % index of user 2 coded blocks 
    i2 = 
% for i = 1:(size(x1,1)/n1)
    % pick n1 coded blocks of user 1, transpose to get the interleaved one
    xp1 = x1((1 + (i-1)*n1):i*n1, :)';
    % and n2 coded blocks of user 2, transpose to get the interleaved one
    xp2 = x2((1 + (i-1)*n2):i*n2, :)';
    %
    % seems like we only choose a single (i.e. common m1) such that error
    % correction capability for both users are satisfactory
    %
    for j = 1:m1
        % for each row, construct the binary block 
        xp1_bin = binary_block(xp1, log2(r1)); 
        xp2_bin = binary_block(xp2, log2(r2)); 
        % Fc_1 and Fc_2 characterize the index of extremal channels we transmit on
        X1(j, :) = reshape(xp1_bin, size(xp1_bin, 1) * size(xp2_bin, 2), 1);
        X2(j, :) = reshape(xp2_bin, size(xp2_bin, 1) * size(xp2_bin, 2), 1); 
    end 
        % 
        % because the encoding and decoding process on inner polar code is
        % transparent, we only need the error rate on these channels
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %               equivalent channel based on the error rate on
        %               extremal channels
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        e1 = bit_channel_error(Fc_1); 
        e2 = bit_channel_error(Fc_2); 
        %
        % each of the selected n1 (or n2) RS block codes for user 1 (or
        % user 2) goes through the same extremal channel (i.e. a fixed
        % index). so, each m1 block of polar code, gives the noisy version
        % of [n1 x (m1-length)] codewords. for this reason, we can assume
        % that the selected codeword blocks are transmitted through a BSC
        % with respective error probability (on extremal channels).
        % 
        % transmit on BSC channels (we have n1 blocks for user 1, and n2
        % blocks for user 2 (of course in binary form)
        %
        % X1 -> BSC_1, X2 -> BSC_2 
        % 
        X1_hat = first_level_dec(e1, X1); 
        X2_hat = first_level_dec(e2, X2); 
        % 
        % extracted (n1 x m1-length) codewords
        for i = 1:size(X1_hat, 1)/r1
            tmp(:, (1 + (i-1)*r1 : i* r1)) = X1_hat((1 + (i-1)*r1 : i*r1) , :)';
        end
        s1_hat = reshape(bin2dec(tmp'), m1, n1)';
        s2_hat = reshape(bin2dec(tmp'), m1, n1)'; 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %            RS error correction and 2nd-level decoding
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        % 
        % extract the raw data, compare and count the number of
        % (uncorrected) errors 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

    
