% in this code we want to numerically demonstrate the leakage from extremal
% channels which are chosen as secure
%
% our algorithm achieves zero leakage (as well as strong secrecy) with
% infinite block length at which a11 arikan channels exactly polarize to
% one of the extremal states
%
% however, in finite block length there will be some form of leakage.
% previously, the upper bound on block error rate was used for this
% purpose. now we can use the actual mutual information for leakage.
%
%
% first identify the pattern of secure bit channels with a sufficiently
% small epsilon
%
% now we choose the channels which are close to norm-1 difference at Bob
% w.r.t accurate epsilon case, while gradually deviate from it at Eve. if
% this yields effective numerical insights it means that we have kept the
% reliability constraint and increased the leakage.
%
% we previously had observed that increasing epsilon will increase the
% secrecy rate of our algorithm at the cost of losing the accuracy. the new
% approach gives a different interpretation for this observation. in fact,
% if we stick to a fixed accuracy for reliability at Bob and increase
% epsilon only in the process of choosing the subset of secure bit chnnels
% in this set, we will see how the leakage is incresed at a specific error
% bound limit at Bob (hopefully, this is plausible).
%
%
%
% function [secrecy_rate, leakage_rate, M_leak] = leakage(Epsilon_1, Epsilon_2, eta, filename, foldername_b)
% 
% function [secrecy_rate, leakage_rate, M_leak] = leakage(Epsilon_1, Epsilon_2, eta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
%
bob_sheet = 1;
eve_sheet = 2;
t = 2;
block_length = 256;
% filename = 'secrecy_graph_4_5_15.xlsx';
filename = 'secrecy_graph_e_5.xlsx';
start_range = 'A1:';
end_range = char('A'+ 2^t - 1);
sheet_range = strcat(start_range,end_range,num2str(block_length));
%
Epsilon_1 = 0.6;
Epsilon_2 = 1;
% % %
eta = .05;
%
foldername_b = 'debug_channels_10_p11_2';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ok- let's start, first extract the integer values of Bob with epsilon 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = xlsread(filename,bob_sheet,sheet_range);
% the differences in l1-norms
l1_norms = zeros(size(R,1),5);
% the extremal MACs that Arikan channels will converge to when block_length
% goes to infinity
rank_funs = [0 0 0;1 0 1;0 1 1;1 1 1;1 1 2];
for i = 1:5
    % sum of differences of our samples from this candidate rank function
    l1_norms(:,i) = sum(abs(bsxfun(@minus, R(:,2:end), rank_funs(i,:))),2);
end
% this is the index of minimum norm-1 elements
[val_tmp,ix_tmp] = min(l1_norms,[],2);
% pick the ones which satisfy <= epsilon_1 constraint
ix_1 = find(val_tmp <= Epsilon_1);
rates_b = rank_funs(ix_tmp(ix_1),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now create the integer rates at Eve over this time indices- ok, if in
% this set of indices some bit channels are not close to integer values by
% measure Epsion_2, they will not get selected. we only work on ix_2
% indices which corresponds to the intersection of two cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_e = xlsread(filename,eve_sheet,sheet_range);
Re = tmp_e(ix_1,:);
l1_norms = zeros(size(Re,1), 5);
for i = 1:5
    % could do it inside the bob loop, anyway
    l1_norms(:,i) = sum(abs(bsxfun(@minus, Re(:,2:end), rank_funs(i,:))),2);
end
% this is the index of minimum norm-1 elements
[val_tmp,ix_tmp] = min(l1_norms,[],2);
% the intersection- pick the ones which satisfy <= epsilon_1 constraint
ix_2 = ix_1(val_tmp <= Epsilon_2);
rates_e = rank_funs(ix_tmp(val_tmp <= Epsilon_2),:);
rates_b = rates_b(val_tmp <= Epsilon_2,:);
%
fprintf('integer rate vectors with epsilon_1 and epsilon_2 extracted..\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now run the secure algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rates = [rates_b; rates_e];
[user_grid, ix_3, rates_b, rates_e] = find_secure_v23(rates, ix_2, block_length);
ix_secure = any(user_grid, 1);
secure_slots = find(ix_secure);
%
% the goal is to find rates of bob and eve at secure slots. right now we
% only have those rates on the valid slots (on which both receivers had a
% fairly close integer rate vector w.r.t Epsilon
index_4 = ismember(ix_3, secure_slots);
rates_b = rates_b(index_4,:);
rates_e = rates_e(index_4,:);
% therefore, we have the time indices and th corresponding rate vectors
% over the secure slots (with parameter Epsilon_1 at Bob and Epsilon_2 at
% Eve
fprintf('found secure time slots on extracted rate vectors at Bob & Eve..\n')
fprintf('next upper bound on block error rate on secure slots at Bob will be computed..\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we now need to compute the reliability over each of these indices and
% sort them in ascending order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
upper_b = zeros(1,length(secure_slots));
for i = 1:length(secure_slots)
    % find the binary matrix corresponding to each rate vector
    [B,~,~] = find_bmat([0, rates_b(i,:)]);
    % calculate upper bound on error probability on this channel
    B(:,~user_grid(:,secure_slots(i))) = 0;
    ix_single = sum(bsxfun(@times, B, 2.^((size(B,1)-1):-1:0)),2);
    table_id = dec2bin(secure_slots(i),log2(block_length));
    % construct the single user channel
    su_prob_v2(foldername_b, table_id,'','bob');
    % Z parameters for single user channels
    load su_pars_b
    Z1 = z_par;
    % remove zero rows
    ix_single(~ix_single) = [];
    upper_b(i) = sum(Z1(ix_single));
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finally, choose the time slots which yield a total block error rate less
% than eta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ix_cum = 0;
% sort based on the value of upper bound in Bob
[upper_bs, ix_sort] = sort(upper_b,'ascend');
ts_secure_old = secure_slots;
ts_secure = secure_slots(ix_sort);
%
for i = 1:length(ts_secure)
    tmp = sum(upper_bs(1:i));
    if tmp > eta
        break;
    end
    % the last value that agrees with eta constraint is stored in ix_cum
    ix_cum = i;
end
% ok- this is the final set of reliable and secure time slots
ts_secure(ix_cum + 1: end) = [];
% the order for mutual information vector
ix_rates = i_vec_order(2^t);
% find the equivalent decimal value
dec_id = sum(bsxfun(@times, user_grid(:,ts_secure),(2.^(0:(t-1)))'),1);
%
% ok- not very efficient maybe, but I can't think of a vectorize solution
% ofcourse, I can only find the equivalent of the unique elements in dec_id
% for example, with t = 2, the loop will only need to find 3 indices. the
% efficiency will come in large sizes. I don;t mind use the first idea here
% tic
% [dec_id_u, idec, iu] = unique(dec_id);
% ix = zeros(size(dec_id_u));
% for i = 1:length(dec_id_u)
%     ix(i) = find(ix_rates == dec_id_u(i));
% end
% ixx = ix(iu);
% toc
% %
ix = zeros(size(dec_id));
for i = 1:length(dec_id)
    ix(i) = find(ix_rates == dec_id(i));
end
% submatrix corresponding to Eve rates on secure slots
ix_leak = sub2ind([block_length, 2^t], ts_secure,ix);
leakage_bit_chans = tmp_e(ix_leak);
leakage_tot = sum(leakage_bit_chans);
% should be done now.
secrecy_rate = sum(sum(user_grid(:,ts_secure)))/block_length;
leakage_rate = leakage_tot/block_length;
%
fprintf('\n Total leakage per block is %5.6f', leakage_tot);
fprintf('\n Leakage rate is %5.6f', leakage_rate);
fprintf('\n Total number of secure bits per block is %d', sum(sum(user_grid(:,ts_secure))));
fprintf('\n Secrecy rate is %5.6f', secrecy_rate);
fprintf('\n for reference, the norm-1 parameters are %3.2f, %3.2f', Epsilon_1, Epsilon_2);
fprintf('\n Reliability parameter is %1.0e\n',eta);
fprintf('---------------------------------------------------------\n');
fprintf(' Secure and reliable time slots are \n %s \n \n', num2str(ts_secure));
%
[R(ts_secure,2:4), tmp_e(ts_secure,2:4)];
Delta = leakage_tot/sum(sum(user_grid(:,ts_secure)));
if ~isempty(leakage_bit_chans)
    M_leak = max(leakage_bit_chans);
else
    M_leak = 0; 
end
