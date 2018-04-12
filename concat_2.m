% Ok- now we should start working around the concatenated solution, and
% first you need to review some theory about coding and of course your
% notes about the idea.
%
% if I don't write anything right now, I wouldn't have anything at the end
% of the day
%
% block length for polar code
N_p = 256;
% alphabet size (each sybol of outer code is t bits)
t = 4;
% block length for outer code- this is the maximum possible one
N_o = 2^4 - 1;
% number of users
n_users = 2;
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
% this requires us to know the error probability on polar code bit channels
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               equivalent channel based on the error rate on
%               extremal channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% error probability of user 1 on non frozen bit channels Fc_1
e1 = bit_channel_error(Fc_1);
% error probability of user 2 on non frozen bit channels Fc_2
e2 = bit_channel_error(Fc_2);
% 
% find the rate of each outer code (i.e. number of parity symbols)
% given the alphabet size t, we know which bit channels must be used in
% formula (2). we also need to pich the length of RS code: m
%
tau_1 = find_RS_rate(m, t, e1); 
tau_2 = find_RS_rate(m, t, e2); 
% 
% ok, now we know that each block of concatenated code requires:
% k1/t and k2/t outer code (therefore k1 * m and k2 * m information bits
% are transmitted) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                raw data generation and outer block encoding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for L = 1:nsamp_1
    % no that we know the structure and number of data bits in each RS
    % block, we can generate them 
    n_bits_1 = (m * k1) - 2 * t * sum(tau_1); 
    n_bits_2 = (m * k2) - 2 * t * sum(tau_2); 
    s1 = mod(randi(2,1, n_bits_1), 2); 
    s2 = mod(randi(2,1, n_bits_2), 2); 
    % now encode these into RS codewords of length m- this should be a
    % k1/t by m (resp. k2/t by m) rectangular matrix at first user (resp. 
    % at second user)
    x1 = RS_encoder(s1, t, tau_1); 
    x2 = RS_encoder(s2, t, tau_2); 
    % and next construct the interleaved codewords 
    x1_intr = x1'; 
    x2_intr = x2'; 
    % 
    xp1 = zeros(m, k1); 
    xp1_noisy = xp1; 
    x1_noisy = xp1; 
    %
    xp2 = zeros(m, k2);
    xp2_noisy = xp2; 
    x2_noisy = xp2; 
    %
    for i = 1:m
        tmp = reshape(dec2bin(x1_intr(i,:), t)', 1, m * k1 / t); 
        % ok- this is data of user 1 on i-th polar codeword 
        xp1(i,:) = str2num(tmp')'; 
        xp1_noisy(i,:) = eqv_bsc(xp1, e1);
        %
        tmp = reshape(xp1_noisy, t, k1/t); 
        tmp1 = num2str(tmp(:)); 
        x1_noisy(i,:) = bin2dec(reshape(tmp1, t, k1/t)')';    
        %
        tmp = reshape(dec2bin(x2_intr(i,:), t)', 1, m * k2 / t); 
        % and this data of user 2 on i-th polar codeword
        xp2(i,:) = str2num(tmp')';
        xp2_noisy(i,:) = eqv_bsc(xp2, e2);
        % 
        tmp = reshape(xp2_noisy, t, k2/t); 
        tmp1 = num2str(tmp(:)); 
        x2_noisy(i,:) = bin2dec(reshape(tmp1, t, k2/t)')';    
    end 
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


