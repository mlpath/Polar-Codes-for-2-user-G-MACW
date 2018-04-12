% in this file we test the Information/Fixed bit patterns obtained by
% find_secure.m file by running actual encoding/decoding for users
%
clc
clear all
close all 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Input from find_secure.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by default we start from column A1 and the code length is N = 256
% [~, user_grid, ~] = find_secure('secrecy_test.xlsx', 2, 1);
% 
user_grid = zeros(2,256);
time_slots = 144;%64]%,75:96,98:128,250:256];
user_grid(:,time_slots) = ones(2,length(time_slots));
k = sum(user_grid, 2); 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       MAC channel parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discrete channel filename for bob
mac_chan_b = 'raw_channel_14.mat';
mac_chan_e = 'raw_channel_10.mat';
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           All other inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% number of samples generated to estimate block error probability
nsamp = 6000;
% number of users 
t = size(user_grid, 1); 
% block length
N = size(user_grid, 2);
% counter for block error decoding 
err_b = zeros(t, 1); 
err_e = err_b; 
%
% initialize raw data matrix (three dimensional array)
U_3d = zeros(nsamp,t , N);
% 
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
for j = 1:nsamp
    j
    % a t*N coded block X
    X = X_2d(j:nsamp:end,:); 
    size(X);
    % received signals at the output of quantized (degraded) MACs
    % at bob
    y = dmc(X, mac_chan_b);
    % at eve
    z = dmc(X, mac_chan_e);
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Single user polar decoding 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % compute the likelihoods for marginal channels of each user at bob
    lr_bob = marg_mac_likelihood(mac_chan_b, y);
    % compute the likelihoods for marginal channels of each user at eve
    lr_eve = marg_mac_likelihood(mac_chan_e, z);
    %
    % this matrix together with user_grid gives the chosen value for fixed bits
    tmp = U_2d(1:nsamp:t*nsamp,:); 
    %
    % we decode signals of all users at each time sample and count errors at Bob and Eve 
    for i = 1:t
        % frozen index for user i, these are known by both receivers
        ix_f = find(~user_grid(i,:));
        % frozen bits for user i, these are announced publicly 
        uc = tmp(i, ix_f); 
        % estimate data of user i at Bob
        [u_hat_b,~] = sc_decoder(y,ix_f,uc,'marg_gmac',lr_bob(i,:),N);
        % estimate data of user i at Eve
%         [u_hat_e,~] = sc_decoder(z,ix_f,uc,'marg_gmac',lr_eve(i,:),N);
        % 
%         isequal(u_hat_e(ix_f),U_2d((i-1)*nsamp + j,ix_f));
        %
        if ~isequal(u_hat_b, U_2d((i-1) * nsamp + j,:))
            err_b(i) = err_b(i) + 1
        end
        % 
%         if ~isequal(u_hat_e, U_2d((i-1) * nsamp + j,:))
%             err_e(i) = err_e(i) + 1;
%         end
        %
    end 
end
% each row is error rate of one user at respective receivers
block_error_b = err_b/nsamp
block_error_e = err_e/nsamp; 
%
