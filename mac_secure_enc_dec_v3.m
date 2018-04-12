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
%                       Input from find_secure.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by default we start from column A1 and the code length is N = 256
[~, user_grid, ~] = find_secure_v2('secrecy_test.xlsx',.05 , 1, 2);
% 
fprintf('secure bits were chosen \n');
% for debug purpose
% user_grid = zeros(2,256);
% time_slots = [144];%64]%,75:96,98:128,250:256];
% user_grid(:,time_slots) = ones(2,length(time_slots));
% before this, we transmitted over all such symbols, let's do it again
k = sum(user_grid, 2); 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       MAC channel parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
power = 2*[1 1]';
noise_var = .97865^2;
%
gains_b = [1 1.4]'; 
gains_e = [1 .4]'; 
%
mac_chan_b = struct('P', power, 'G', gains_b, 'N', noise_var);
mac_chan_e = struct('P', power, 'G', gains_e, 'N', noise_var);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           All other inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% number of samples generated to estimate block error probability
nsamp = 2000;
% number of users 
t = size(user_grid, 1); 
% block length
N = size(user_grid, 2);
%  
fprintf('Achievable sum secrecy rate is %f',k/N);
% counting block errors (i.e. a system with input symbols in {0,1, 2^t-1} 
% and thus making an error if any of users is decoded incorrectly
Err_b = 0; 
Err_e = 0; 
% counting blocks of each users individually, don't know if this is
% necessary or theoretecally valid, just want to test
err_b = zeros(t,1); 
err_e = err_b; 
% initialize raw data matrix (three dimensional array)
U_3d = zeros(nsamp,t , N);
% generate information and frozen bits for users in the MAC network
for i = 1:t
    % for each user generate information bits
    U_3d(:, i, find(user_grid(i,:))) = mod(randi(2, nsamp, k(i)),2);
    % for all samples, information on these slots is fixed
    U_3d(:, i, find(~user_grid(i,:))) = repmat(mod(randi(2, 1, N - k(i)),2),nsamp,1);
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
ix_f = ~user_grid;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       continuous AWGN MAC channel 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% need to repeat this experiment (coding/decoding MAC signal) nsamp times
for j = 1:nsamp
    j
    % a t*N coded block X
    X = X_2d(j:nsamp:end,:);
    size(X);
    % received signals at the output of quantized (degraded) MACs
    % at bob
    y = c_mac_awgn(X, mac_chan_b);
    % at eve
    z = c_mac_awgn(X, mac_chan_e);
    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Joint ML successive decoding on nonfrozen bits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % decode users data at Bob
    [u_hat_b, u_bin_b] = mac_decoder(y, mac_chan_b, uc, ix_f, N);
    % decode users data at Eve
    [u_hat_e, u_bin_e] = mac_decoder(z,mac_chan_e, uc, ix_f, N);
    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Block error calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    % I assume a block error event is occured if receiver did not recover
    % any of user's data correctly 
    
    
    % convert 0/1 s to integers in {0,1,..2^t-1}
    u = sum(bsxfun(@times, U_2d(j:nsamp:end,:),[2 1]'),1);
    %
    if ~isequal(u_hat_b, u)
        Err_b = Err_b + 1;
    end
    %
    if ~isequal(u_hat_e, u)
        Err_e = Err_e + 1;
    end
    % 
    fprintf('iteration %d Bob error %d Eve error %d', j, Err_b, Err_e);
    for i = 1:t
    % counting block error for each user individually (test it at Eve)
%         if ~isequal( U_2d(j + (i - 1) * nsamp,:), u_bin_e(i,:))
%             err_e(i) = err_e(i) + 1; 
%         end
        %
        f_ix = find(ix_f(i,:));
        % finally, test if frozen indices are being decoded correctly
        isequal( U_2d(j + (i - 1) * nsamp,f_ix), u_bin_b(i,f_ix));
    end
end
% each row is error rate of one user at respective receivers
block_error_b = Err_b/nsamp;
block_error_e = Err_e/nsamp; 
%
% block error computed separately for each user
% err_e = err_e/nsamp; 
% 
% seems like we are done with scripting
