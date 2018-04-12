% in this file we test the Information/Fixed bit patterns obtained by
% find_secure.m file by running actual encoding/decoding for users
%
clc
clear all
close all 
%
% this is somehow confusing to me:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Some notational conventions in interpreting order of users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in mac_awgn.mat, where we have implemented our continuous Gaussian
% probability table, we calculate values based on binary order:
% 00, 01, 10, 11
% so, rows 1 and 4 do not reveal anything about order. However, in rows 2
% and 3, we multiplay MSB by first channel gain and LSB by second channel
% gain. That is,  
%                 user 1  user 2   row #        conversion rule
%                -------------------------------------------------------
%                   0       0       1      [user 1, user 2] .* [2, 1] + 1
%                   0       1       2              
%                   1       0       3
%                   1       1       4
%
% don't know how it affects the other parts, but I satrt from pr_chan_level
%
% also, the generated data for transmission is like,
%
%               |   |   |   |   |   |   |   |   |   |   |   <- user 1
%               -----------------------------------------
%               |   |   |   |   |   |   |   |   |   |   |   <- user 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% this code loops over channels with different transmitted power 
power_id = [5 10 15 20]; %25 30 40];
% 
block_error_b = zeros(1,length(power_id));
block_error_e = block_error_b; 
%
Err_b = 0; 
Err_e = 0; 
%
for m = 1:length(power_id)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Input from find_secure_v2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
capacity_id = strcat('secrecy_graph_',num2str(power_id(m)),'.xlsx');
% by default we start from column A1 and the code length is N = 256
[rate_vectors, user_grid, ~] = find_secure_v2(capacity_id,.058, 1, 2);
% 
% for debug purpose
% user_grid = zeros(2,256);
% time_slots = [144];%64]%,75:96,98:128,250:256];
% user_grid(:,time_slots) = ones(2,length(time_slots));
%
% all secure time slots
ix_secure = any(user_grid,1); 
% sort all rate vectors
[~, all_sorted] = sort(rate_vectors, 'descend'); 
% sort ix_secure accordingly
ix_sorted = ix_secure(all_sorted); 
% time slots in increasing order of sum secrecy rate
time_slots = all_sorted(ix_sorted); 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           All other inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% number of samples generated to estimate block error probability
nsamp = 4000;
% number of users 
t = size(user_grid, 1); 
% block length
N = size(user_grid, 2);
% 
% initialize raw data matrix (three dimensional array)
U_3d = zeros(nsamp,t , N);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       MAC channel parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mac_chan_b = strcat('C:\MonaHajimomeni\mfiles\debug_mac\raw_channel_10_p',num2str(power_id(m)),'.mat');
mac_chan_e = strcat('C:\MonaHajimomeni\mfiles\debug_mac\raw_channel_4_p',num2str(power_id(m)),'.mat');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    How many time slots are we using?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% for the time being we are not using this loop here
%
% to preserve the previous setting (using all slots) start the loop from L
L = length(time_slots); 
% over all possible secure slots 
for l = L:L
    if l ~= L
        % if l has any value between 1 and L-1, make time slots after that
        % frozen, otherwise don't alter user_grid- we transmit on all slots
        new_grid = user_grid;
        new_grid(time_slots(:,(L+1):end)) = 0; 
    else 
        new_grid = user_grid;
    end
    sum(any(new_grid,1))
    % 
    % this is number of information slots given to each user
    k = sum(new_grid, 2);
    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Raw messages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% now generate information and frozen bits for users in the MAC network
for i = 1:t
    % for each user generate information bits
    U_3d(:, i, find(new_grid(i,:))) = mod(randi(2, nsamp, k(i)),2);
    % for all samples, information on these slots is fixed
    U_3d(:, i, find(~new_grid(i,:))) = repmat(mod(randi(2, 1, N - k(i)),2),nsamp,1);
end
% change the dimension and apply polar encoding on all users at once 
U_2d = reshape(U_3d, t*nsamp, N);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Polar encoding 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G_N = 1;
F = [1 0; 1 1];
for i = 1:log2(N)
    tmp = reverse_shuffle(2^i) * kron(F,G_N);
    G_N = tmp;
end
%
X_2d = mod(U_2d * G_N,2);
%
% a (t x N) matrix having only one sample of transmitted data
uc = U_2d(1:nsamp:t*nsamp,:);
% this is going to be used for extracting frozen data, together with
% the index of frozen/transmitting users (in ix_f below)
ix_f = ~new_grid;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       discrete AWGN MAC channel 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% need to repeat this experiment (coding/decoding MAC signal) nsamp times
for j = 1:nsamp
    [j,m,length(power_id)]
    % a t*N coded block X
    X = X_2d(j:nsamp:end,:);
    size(X);
    % received signals at the output of quantized (degraded) MACs
    % at bob
    y = dmc(X, mac_chan_b);
    % at eve
    z = dmc(X, mac_chan_e);
    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Joint ML successive decoding on nonfrozen bits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % decode users data at Bob
    [u_hat_b, u_bin_b] = mac_decoder_disc(y, mac_chan_b, uc, ix_f, N);
    % decode users data at Eve
    [u_hat_e, u_bin_e] = mac_decoder_disc(z,mac_chan_e, uc, ix_f, N);
    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Block error calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % from the viewpoint of vector messages MAC and nonbinary alphabets are
    % the same. Think of each vector message as a symbol between 0..3 (for
    % 2user), now a block error is defined if we do not recover the
    % transmitted symbol (i.e. sent 2 but received 1, or sent 3 but
    % received 0, etc.)
    % 
    % convert 0/1 s to integers in {0,1,..2^t-1}
    u = sum(bsxfun(@times, U_2d(j:nsamp:end,:),[2 1]'),1);
    %
    if ~isequal(u_hat_b, u)
        Err_b = Err_b + 1
    end
    %
    if ~isequal(u_hat_e, u)
        Err_e = Err_e + 1
    end
    % 
%     for i = 1:t
%     % counting block error for each user individually (test it at Eve)
% %         if ~isequal( U_2d(j + (i - 1) * nsamp,:), u_bin_e(i,:))
% %             err_e(i) = err_e(i) + 1; 
% %         end
%         %
%         f_ix = find(ix_f(i,:));
%         % finally, test if frozen indices are being decoded correctly
%         isequal( U_2d(j + (i - 1) * nsamp,f_ix), u_bin_b(i,f_ix));
%     end
%     fprintf('iteration %d, [Bob ,Eve] error %d %d, and diff %d \n',j,Err_b, Err_e, Err_e - Err_b);
end
% arrays to keep block error probability for each one more subchannel we
% use for security (to test if unreliable channels are selected or not)
%
end
block_error_b(m) = Err_b/nsamp
block_error_e(m) = Err_e/nsamp
%
Err_b = 0; 
Err_e = 0;
end
p = 10*log10(power_id/10/(.97865^2));
semilogy(p, block_error_b, '-*'); 
hold on; 
semilogy(p, block_error_e, '-r*');
