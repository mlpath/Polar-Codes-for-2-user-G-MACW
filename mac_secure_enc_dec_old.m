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
user_grid = [zeros(2,249),ones(2,7)];
k = sum(user_grid, 2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       MAC channel parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel parameters 
sigma2 = .97865^2; 
% transmission power vector 
P = [2 2]; 
% user gains at respective receivers 
h = [1 1.4]';
g = [1 1]'; 
% save all in a structure to pass to su_decoder
chan_spec = struct('tx_power',P,'bob_gain',h,'eve_gain',g,'noise_var',sigma2); 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           All other inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% number of samples generated to estimate block error probability
nsamp = 100;
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
% just for test
X_2d_old = X_2d;
% modulate 0/1 bits to +/- sqrt(P) levels
for i = 1:t 
    tmp = sqrt(P(i)) * ones(nsamp, N); 
    lptr = (i-1)*nsamp + 1; 
    lptr2 = lptr + nsamp - 1; 
    %
    tmp( X_2d(lptr:lptr2, :) == 1 ) = -sqrt(P(i));
    X_2d(lptr:lptr2, :) = tmp;
end
%
X_3d = reshape(X_2d, nsamp, t, N);
%
for j = 1:nsamp
    j
    % a t*N coded block X
    X = X_2d(j:nsamp:end,:); 
    size(X);
    % received signals multiplied by channel gains, then sum and add noise
    % at bob
    y = sum(bsxfun(@times, X, h),1) + sqrt(sigma2) .* randn(1,N);
    size(y)
    % at eve
    z = sum(bsxfun(@times, X, g),1) + sqrt(sigma2) .* randn(1,N); 
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Single user polar decoding 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % compute the likelihoods for marginal channels of each user at bob/eve
    LK = marg_mac_likelihood(chan_spec,[y;z]);
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
        [u_hat_b,~] = sc_decoder(y,ix_f,uc,'marg_gmac',LK(i,:),N);
        % estimate data of user i at Eve
        [u_hat_e,~] = sc_decoder(z,ix_f,uc,'marg_gmac',LK(i+t,:),N);
        % 
        isequal(u_hat_e(ix_f),U_2d((i-1) * nsamp + j,ix_f))
        %
        if ~isequal(u_hat_b, U_2d((i-1) * nsamp + j,:))
            err_b(i) = err_b(i) + 1
        end
        % 
        if ~isequal(u_hat_e, U_2d((i-1) * nsamp + j,:))
            err_e(i) = err_e(i) + 1
        end
        %
    end 
end
% each row is error rate of one user at respective receivers
block_error_b = err_b/nsamp; 
block_error_e = err_e/nsamp; 
%
