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
% clc
% clear all
%
function  [R, user_grid, R_eve] = find_secure_v3( filename, fidelity_par, block_length,  bob_sheet, eve_sheet )
% 
% clc
% clear all
%
% filename = 'secrecy_test.xlsx';
% bob_sheet = 2;
% eve_sheet = 1;
% fidelity_par = .5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Inputs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% length of block code
N = block_length;
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
Epsilon = fidelity_par;
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
% the simplest way- we have 5 possible rank functions
rank_funs = [0 0 0;1 0 1;0 1 1;1 1 1;1 1 2];
% store l1 norm of our samples with respect to each of the these rank vecs
l1_norms = zeros(size(R,1),5);
for i = 1:5 
    % sum of differences of our samples from this candidate rank function
    l1_norms(:,i) = sum(abs(bsxfun(@minus, R(:,2:end), rank_funs(i,:))),2);
end
% one way to go is to find index of columns with error less than our
% desired value (i.e. Epsilon)
% 
% index of samples in Bob closer than epsilon to a rank function
ix1 = l1_norms(1:N,:) <= Epsilon;
sum(sum(ix1));
tmp = sum(ix1,2) > 1;
% clear such rows
ix1(tmp,:) = zeros(sum(tmp),5);
% find minimum difference element on each sample that more than one rank
% function can  be assigned
[~,ix] = min(l1_norms(tmp,:),[],2);
tmp = find(tmp);
for i = 1:length(tmp)
    ix1(tmp(i),ix(i)) = 1;
end
sum(sum(ix1));
% index of samples in Eve closer than epsilon to a rank function
ix2 = l1_norms(N+1:2*N,:) <= Epsilon; 
sum(sum(ix2));
tmp = sum(ix2,2) > 1; 
% clear such rows
ix2(tmp,:) = zeros(sum(tmp),5); 
[~,ix] = min(l1_norms(tmp,:),[],2); 
tmp = find(tmp); 
for i = 1:length(tmp)
    ix2(tmp(i),ix(i)) = 1;
end
sum(sum(ix2));
% now each sample is matched with only one matroid
% 
ix = [any(ix1,2),any(ix2,2)];
% find those rows that both Bob and Eve are close to a rank functions
valid_ts = find(all(ix,2));
length(valid_ts);
% now tell us which rank category they fall in? (make v_R)
[I1,J1] = ind2sub([length(valid_ts),5],find(ix1(valid_ts,:)));
[I2,J2] = ind2sub([length(valid_ts),5],find(ix2(valid_ts,:)));
% 
[~,u1] = sort(I1); 
[~,u2] = sort(I2); 
J1 = J1(u1);
J2 = J2(u2);
%
int_R = zeros(2*N,2^t);
int_R(valid_ts,2:end) = rank_funs(J1,:);
int_R(valid_ts + N,2:end) = rank_funs(J2,:);
v_R = int_R([valid_ts, valid_ts + N],:);
length(valid_ts);
% 
% int_R now has closes rank function to our estimated rate vector
% and valid_ts contains index of samples chosen for secure algorithm
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
        fprintf('invalid matrix occurred');
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
            [valid_ts(i) rs];
        end
        % now tell us wich users should transmit on this slot 
        user_idx(secure_users, lptr ) = 1;
        % sum secrecy rate on valid time slots
        Rs(lptr) = rs; 
        lptr = lptr + 1;
    end
    %
    [valid_ts(i), secure_users];
end
%
length(valid_ts);
lptr - 1;
user_idx(:,(lptr + 1):end) = [];
Rs((lptr + 1):end) = [];
%
user_grid = zeros(t,N_old);
user_grid(:, valid_ts) = user_idx; 
sum_secrecy_rate = zeros(1,N_old); 
sum_secrecy_rate(valid_ts) = Rs; 
%
% return the actual rate vectors of Bob samples to be used in reliability
% only return the last column (rank of this rate vector)
% 
R_eve = R(N_old + 1:end,:);
R = R(1:N_old,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_old;
fprintf('sum secrecy rate is %f \n',sum(sum(user_grid))/N_old);
user_grid;
sum(sum(user_grid))/256;
length(valid_ts);
valid_ts(Rs>0)
sum(sum_secrecy_rate)/256;
