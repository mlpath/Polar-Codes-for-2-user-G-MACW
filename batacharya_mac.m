% I'm writing this one to see if there is anyway for finding the correct
% rate of secrecy
% 
% first, we need to create the part that 
%
% find_secure_v2 finds the bitchannels that are suitable for transmission
% of secure bits given an approximation parameter epsilon
% 
% in mac_secure_enc_dec, we sort these time slots based on the value of
% joint mutual information 
% 
% we can now write a code that sorts time slots based on Batachariya
% parameters, and then:
%
%  (1) Compare this order with previous sorting based on mutual information
%  (2) See the actual value of error bounds (because when we have a
%      bitchannel with a small upper bound, we can guarantee reliable
%      channel, but when this upper bound is not small enough, we need to
%      have a lower bound)
% 
clc
clear all
close all
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Configuration parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% our value for approximation of rate vectors and matroid's rank function
Epsilon = .5;%.005;
% upper_bound limit 
u_limit = .0005;%.08;
% choose only one (pair of) channel: e.g. [1 1], [1 0.4] with p = 1
excel_id = '1024';
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   First go ahead with find_secure_v2                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% returns selected subchannels and their respective rate_vectors at Bob
filename_str = strcat('secrecy_graph_', excel_id, '.xlsx');
[rate_vectors, user_grid, out_matrix] = find_secure_v22(filename_str, Epsilon, 1024, 1, 2);
% I think we must now find the upper/lower bounds
% 
% index of all secure time slots
ix_secure = any(user_grid,1);
% time_slots of secure channels
ts_secure = find(ix_secure);
% 
if isempty(ts_secure)
    fprintf('------ secrecy outage, no time slot is available ----- \n');
    fprintf('\nhint: try changing parameters (power, epslion, etc.) \n \n');
    break;
else
% pick all rows with secure time slot (and any column)
% (rate_vector contains the information vectors of all time slots)
r_secure = rate_vectors(ts_secure,:);
% isequal(ts_secure,out_matrix(out_matrix(:,2)>0,1)');
% 
bound = zeros(length(ts_secure), 2); 
% and display the final bound as well as secure time slots and bounds for
% each of them 
fprintf('time slot\tlower bound \tupper_bound \t details \n');
fprintf('-----------------------------------------------------\n');
%
e_flag = '';
count = 0; 
% 
% over all these time slots,
for i = 1:length(ts_secure)
    % find the binary matrix - columns 3 to 6 (end) contain the rate vector
    % first column has valid time slots
    % second column accomodates achieved secrecy rate
    [B,~,~] = find_bmat(out_matrix(out_matrix(:,1) == ts_secure(i),3:end));
    % ok, now we have the binary matrix corresponding to each rate vector
    % find index of single-user channels in this bit-channel
    ix_single = sum(bsxfun(@times, B, 2.^((size(B,1)-1):-1:0)),2);
    % 
    % this is the id for initial probability table
    table_id = dec2bin(ts_secure(i),8);
    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Make sure single user channels are produced               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    test_str = strcat('su_table_', table_id, '_1.mat');
    %
%     if exist(test_str, 'file') == 0
        % means we haven't computed single-user table for this data- compute
        % this function creats all su conditional prob tables corresp. to input
        % Note that if we have log2(nrows) variables, there are nrows possible
        % subset (or rows) from which the case of all 0 values is trivial
        % because this means there's no two variables which are merging into
        % each other (i.e. it equals the current joint conditional prob table
        su_prob(table_id);
%     end 
    % this loads Z parameters
    load su_pars
    Z = z_par;
    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%       compute now upper and lower bound for this time slot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % this is rather heuristic- anyway, I am going to use sum of batacharya
    % parameters of single user channels for upperbound
    %
    upper = sum(Z(ix_single));
    %
    % and the sum of lower bounds using batacharya parameters (the largest
    % of such values among all subchannels
    %
    % argument inside summation is a vector of the same size as ix_single
    lower = sum(0.5*(1 - sqrt(1 - Z(ix_single).^2)));
    % first column is lower_bound and second is upper_bound 
    bound(i,:) = [lower, upper];
    if upper >= u_limit
        e_flag(i) = '*';
        count = count + 1; 
    else 
        e_flag(i) = ' ';
    end
    %
    fprintf('%3d \t %13.8f \t %13.8f  \t %c \n', ts_secure(i), lower, upper, e_flag(i));
end 
% 
ix_error = e_flag == '*';
% 
fprintf('\n');
fprintf('*: subchannels with upper bound errors larger than %4.3f \n\n', u_limit);  
fprintf('epsilon (approximation parameter): %4.3f \n\n', Epsilon);
fprintf('true secrecy rate is: %5.4f, initial secrecy rate was: %5.4f\n\n', (length(ts_secure) - count)/256, length(ts_secure)/256);
fprintf('final secure time slots:\n') 
disp(num2str(ts_secure(~ix_error)))
% now find error bounds for whole transmitted block
% upper bound

b_bound(2) = sum(bound(~ix_error,2));
% ix is the index of subchannels with highest error bound
[b_bound(1), ix] = max(bound(~ix_error,1));
% and print it too
fprintf('\nlower and upper bound on block error probability (adjusted): [%5.4f, %5.4f] \n',b_bound(1), b_bound(2));

% 
% what is left? 
%
% testing this code with data of secrecy_graph_10 (i.e. p = 1,[1 1],[1.4])
%
% this is what we were testing with simulations
%
% I will change epsilon, then find secure channels, and upper/lower bound
% on selected channels, then decide about whether to keep all channels or
% not. 
% 
% All of this effort is to see if we observe any meaningful result or not.
% In an ideal world, I expect that we can find some suchannels (i.e. having
% nonzero secrecy rate) with acceptable error bound for Bob
end