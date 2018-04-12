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
% clc
% clear all
% close all
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Configuration parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [secrecy_rate, ts_secure, b_bound, val_min, e_bound] = batacharya_mac_v2(epsilon, eta, excel_id)
tic
% our value for approximation of rate vectors and matroid's rank function
% epsilon = .09;%.005;
% % upper_bound limit 
% eta = .002;%.08;
% % choose only one (pair of) channel: e.g. [1 1], [1 0.4] with p = 1
% excel_id = '1024';
% returns selected subchannels and their respective rate_vectors at Bob
filename_str = strcat('secrecy_graph_', excel_id, '.xlsx');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   First go ahead with find_secure_v2                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% changed it to compute for different SNRs...

[~, user_grid, out_matrix] = find_secure_v22(filename_str, epsilon, 256, 1, 2);
% I think we must now find the upper/lower bounds
% 
% index of all secure time slots
ix_secure = any(user_grid,1);
% time_slots of secure channels
ts_secure = find(ix_secure);
bound = zeros(length(ts_secure), 2); 
L2 = zeros(length(ts_secure), 1);
% 
if isempty(ts_secure)
    fprintf('------ secrecy outage, no time slot is available ----- \n');
    fprintf('\nhint: try changing parameters (power, epslion, etc.) \n \n');
    secrecy_rate = 0; 
    b_bound = zeros(1,2);
    ts_secure = [];
    val_min = 0; 
    e_bound = 0;
else
% pick all rows with secure time slot (and any column)
% (rate_vector contains the information vectors of all time slots)
% r_secure = rate_vectors(ts_secure,:);
% isequal(ts_secure,out_matrix(out_matrix(:,2)>0,1)');
% 

% and display the final bound as well as secure time slots and bounds for
% each of them 
% fprintf('time slot\tlower bound \tupper_bound \t details \t\te_lower \t\te_upper \n');
% fprintf('-------------------------------------------------------------------------------------\n');
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
    % remove frozen indices; we are not transmitting on them
    B(:,~user_grid(:,ts_secure(i))) = 0;
    % find index of single-user channels in this bit-channel
    ix_single = sum(bsxfun(@times, B, 2.^((size(B,1)-1):-1:0)),2);
    % 
    % this is the id for initial probability table
    table_id = dec2bin(ts_secure(i),8);
    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Make sure single user channels are produced               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     test_str = strcat('su_table_', table_id, '_1.mat');
    %
%     if exist(test_str, 'file') == 0
        % means we haven't computed single-user table for this data- compute
        % this function creats all su conditional prob tables corresp. to input
        % Note that if we have log2(nrows) variables, there are nrows possible
        % subset (or rows) from which the case of all 0 values is trivial
        % because this means there's no two variables which are merging into
        % each other (i.e. it equals the current joint conditional prob table
        %
        % generate all single user channels for bob
        su_prob_v2(table_id, 'bob');
%     end 
    % this loads Z parameters
        load su_pars_b
        Z1 = z_par;
    %
    su_prob_v2(table_id, 'eve');
    load su_pars_e
    Z2 = z_par;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%       compute now upper and lower bound for this time slot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % this is rather heuristic- anyway, I am going to use sum of batacharya
    % parameters of single user channels for upperbound
    %
    % remove zero rows
    ix_single(~ix_single) = [];
    upper1 = sum(Z1(ix_single));
    upper2 = sum(Z2(ix_single));

    %
    % and the sum of lower bounds using batacharya parameters (the largest
    % of such values among all subchannels
    %
    % argument inside summation is a vector of the same size as ix_single
    lower1 = sum(0.5*(1 - sqrt(1 - Z1(ix_single).^2)));
    lower2 = sum(0.5*(1 - sqrt(1 - Z2(ix_single).^2)));
    % first column is lower_bound and second is upper_bound 
    bound(i,:) = [lower1, upper1];
    L2(i) = lower2;
    if upper1 >= eta
        e_flag(i) = '*';
        count = count + 1; 
    else 
        e_flag(i) = ' ';
    end
    %
%     fprintf('%3d \t %13.8f \t %13.8f  \t %c \t\t\t %13.8f \t %13.8f \n', ts_secure(i), lower1, upper1, e_flag(i), lower2, upper2);
end

% 
ix_error = e_flag == '*';
% 

% fprintf('\n');
% fprintf('*: subchannels with upper bound errors larger than %4.3f \n\n', eta);  
% fprintf('epsilon (approximation parameter): %4.3f \n\n', epsilon);

% now find error bounds for whole transmitted block
% upper bound

if sum(~ix_error)
%     fprintf('\nnumber of selected slots: %3d, and number of rejected ones (in percent): %4.2f \n', length(ix_error), 100*sum(ix_error)/length(ix_error));
%     fprintf('adjusted secrecy rate is: %5.4f, initial secrecy rate was: %5.4f\n\n', (length(ts_secure) - count)/str2double(excel_id), length(ts_secure)/str2double(excel_id));
%     fprintf('final secure time slots:\n')
%     secrecy_rate = (length(ts_secure) - count)/str2double(excel_id);
    secrecy_rate = (length(ts_secure) - count)/256;
    % impose both security and reliability and update ts_secure
    ts_secure = ts_secure(~ix_error);
    disp(num2str(ts_secure));
    b_bound(2) = sum(bound(~ix_error,2));
    % ix is the index of subchannels with highest error bound
    [b_bound(1), ix] = max(bound(~ix_error,1));
    [e_bound, ~] = max(L2(~ix_error));
    % and print it too
%     fprintf('\nlower and upper bound on block error probability (adjusted): [%1.2e, %1.2e] \n',b_bound(1), b_bound(2));
    %
    [val_min, ix_min] = min(L2(~ix_error));
    [val_max, ix_max] = max(L2(~ix_error));
    %
%     fprintf('\nmin and max values for lower bound of secure time slots at Eve: [%5.4f, %5.4f] time slots %3d, %3d \n', val_min, val_max, ts_secure(ix_min), ts_secure(ix_max));
else
    fprintf('no time slot with error bound less than %1.2e exists on bob channel \n', eta);
    secrecy_rate = 0; 
    b_bound = zeros(1,2);
    ts_secure = [];
    val_min = 0; 
    e_bound = 0;
end

%
% 
% All of this effort is to see if we observe any meaningful result or not.
% In an ideal world, I expect that we can find some suchannels (i.e. having
% nonzero secrecy rate) with acceptable error bound for Bob
end
toc
end