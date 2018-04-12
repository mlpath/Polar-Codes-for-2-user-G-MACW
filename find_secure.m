% this code operates on rate vectors of users Bob and Eve to find a secure 
% pattern of frozen and information bits
%
% filename is the name of an excel file where we have written our results
%
% the closeness measure based on which we decide a rate vector is equal to
% the rkfunction of a binary matroid
%
% outputs of this functions are sum secrcey rate, transmitting/fixed users
% grid, and sum rate without secrecy constraint
%
clc
clear all
%
% function  [sum_secrecy_rate, user_grid, sum_rate] = find_secure( filename, bob_sheet, eve_sheet )
% 
filename = 'secrecy_test.xlsx';
bob_sheet = 2;
eve_sheet = 1;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Inputs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% length of block code
N = 256;
% number of users
t = 2;
% range of stored values in excel files. by default we start at A1, Bob is
% saved in sheet 1, Eve in sheet 2, and both have of course the same number 
% of rows equal to block length, N.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filename = 'secrecy_test.xlsx';
start_range = 'A';
end_range = char(start_range + 2^t - 1);
excel_range = strcat(start_range, '1:', end_range, num2str(N));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% this is challenging to me, how should I find closest rk function to my
% capacity vector? (except obvious exponential enumeration)
Epsilon = .4;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Stuff related to identifying binary matroid on time slots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
ncol = 2^t;
% preallocate the matrix
R = zeros(2*N, ncol);
% read from excel file (for now let's assume that there is no special
% format and we just read them independently and later mix them)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% give the [1 1.4] MAC channel to Bob, that means user 2 has better channel
% to Bob than Eve (sheet 2)
R(1:N, :) = xlsread( filename, bob_sheet, excel_range);
% give [1 1] MAC channel to Eve, user 1 have equal gain to each receiver
R((N+1):(2*N), :) = xlsread( filename, eve_sheet, excel_range);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now round it to the nearest integer values
int_R = round( R );
% ix_matroid is the index of rows that represent a rate vector which is
% close enough to one of a binary matroid
% sum of absolute error is L1 norm of error
ix_matroid = sum(abs(int_R - R), 2) <= Epsilon;
% we just care about cases where we have a binary matrix associated with
% both Eve and Bob- valid_ts is the index of time slots we are going to use
valid_ts = (find(all(reshape(ix_matroid, N, 2), 2)))';
% all rate vectors of binary matroids A and B
v_R = int_R([valid_ts, valid_ts + N],:);
% match each rate vector with a binary matroid
%
% over all valid time slot, first find these matroid and next run algorithm
%
N_old = N;
N = size(v_R, 1)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Secure design starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preallocate the array that holds transmitting and fixed users
user_idx = zeros(t, length(valid_ts));
Rs = zeros(1,length(valid_ts));
sum_rate = zeros(size(Rs));
% a counter to write in user_idx
lptr = 1;
%
for i = 1:N
    % on each time slot (i.e. for each parallel subchannel)
    [A, bases_1, is_binary_1] = find_bmat(v_R(i, :));
    [B, bases_2, is_binary_2] = find_bmat(v_R(i + N, :));
    %
    if ~(is_binary_1 && is_binary_2)
        valid_ts(i,:) = [];
        v_R(i,:) = [];
        v_R(i+N,:) = [];
        disp('invalid matrix occurred');
        %
    else
        j = 1;
        % initialize secrecy rate with zero (minimum possible)
        rs = 0;
        secure_users = []; 
        %
        n_bases = size(bases_1,1);
        rk = size(bases_1,2);
        %
        %%%%%%%%%%%%%%%%%%%%%% without secrecy constraint, sum rate will be 
        sum_rate(i) = rk;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this is a valid system, find the secure solution
        while (rs < rk) && (j <= n_bases)
            % pick jth basis of A and construct a submatrix of B with it
            % now apply gauss jordan to a submatrix of B
            [rref,~] = b_gauss_jordan(B(:,bases_1(j,:)));
            % find rows where only a single variable is left
            ix = sum(rref, 2) == 1;
            % explicit solutions are in fact non-secure
            non_secure = find(sum(rref(ix,:),1) == 1);
            % if number of secure bits in this basis is lower than stored one
            if rs < (rk - sum(ix))
                secure_users = bases_1(j,setdiff(1:rk, non_secure));
                rs = rk - sum(ix); 
            end
            j = j + 1; 
        end
        % now tell us wich users should transmit on this slot 
        user_idx(secure_users, lptr ) = 1;
        % sum secrecy rate on valid time slots
        Rs(lptr) = rs; 
        lptr = lptr + 1;
    end
    %
end
%
user_idx(:,(lptr + 1):end) = [];
Rs((lptr + 1):end) = [];
%
user_grid = zeros(t,N_old);
user_grid(:, valid_ts) = user_idx; 
sum_secrecy_rate = zeros(1,N_old); 
sum_secrecy_rate(valid_ts) = Rs; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
valid_ts;
tmp = int_R - R;
% the errors on the chosen channels
sum(abs(tmp([valid_ts, valid_ts + N_old],:)), 2);
fprintf('sum secrecy rate is %f \n',sum(sum_secrecy_rate)/256);
valid_ts((Rs>0))
[valid_ts',Rs'];
