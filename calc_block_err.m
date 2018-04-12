% find block error probability for arbitrary polar channels
% function block_err = calc_block_err(ts_input, user_grid, folder_name, file_name, sheet_no)
clc
clear all
close all
%
block_length = 256; 
% varepsilon = 0.4; 
varepsilon = 0.7;
file_name = 'secrecy_graph_8_10_p20.xlsx';
sheet_no = 1;
folder_name = 'debug_channels_4_p20';
chan_id = '';
% 
R1 = xlsread(file_name, sheet_no, 'A1:D256');
R2 = xlsread(file_name, 2,'A1:D256');
% find channels in bob which are varepsilon-close to [0 1 1 1]
ts_input = find(sum(abs(bsxfun(@minus, R2, [0 1 1 1])),2) < varepsilon);
% 
% user_grid specifies which users are going to transmit on this slots
% of course, we now know that on other good channels we have to transmit,
% but to keep things consistent with other results I assume there's no err
%
% only user 2 can transmit on (111) <-> (101) to have secure communication
user_grid = zeros(2, block_length); 
% only second user transmits, and only on ts_inputs
user_grid(1,ts_input) = 1; 
% 
% initialize the output
block_err = zeros(length(ts_input), 2); 
%
% rate vectors of desired time slots 
tmp = R1(ts_input, :);
for i = 1:length(ts_input)
    % find the binary matrix - of course we suppose it is [ 0 1 0 1] ! 
    B = [1 0;0 0];
    % ok, now we have the binary matrix corresponding to each rate vector
    % remove frozen indices; we are not transmitting on them
    B(:,~user_grid(:,ts_input(i))) = 0;
    % find index of single-user channels in this bit-channel
    ix_single = sum(bsxfun(@times, B, 2.^((size(B,1)-1):-1:0)),2);
    % 
    % this is the id for initial probability table
    table_id = dec2bin(ts_input(i),log2(block_length));
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
        su_prob_v2(folder_name, table_id, chan_id, 'bob');
        %     end
        % this loads Z parameters
        load su_pars_b
        Z1 = z_par;
        %
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
    %
    % and the sum of lower bounds using batacharya parameters (the largest
    % of such values among all subchannels
    %
    % argument inside summation is a vector of the same size as ix_single
    lower1 = sum(0.5*(1 - sqrt(1 - Z1(ix_single).^2)));
    % first column is lower_bound and second is upper_bound 
    block_err(i,:) = [lower1, upper1];
end