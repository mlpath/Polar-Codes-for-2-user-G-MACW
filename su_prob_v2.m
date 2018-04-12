function su_prob_v2(data_folder, input_id, chan_id, user_id)
% this function computes 2^m - 1 possible single user conditional
% probability tables using the joint (conditional) probability table
% characterized by input_id, where m is number of variables (or inputs) of
% the initial channel
% 
if strcmp(user_id,'bob')
    % input is a probability table (i.e. one of .mat files that are generated in training) 
    input_table_str = strcat('D:\April-queens\mfile\backup_final\', data_folder, chan_id, '\chan_', input_id, '.mat');
    out_str = 'su_table_b_';
    par_str = 'su_pars_b.mat';
elseif strcmp(user_id,'eve')
    input_table_str = strcat('D:\April-queens\mfile\backup_final\', data_folder, chan_id, '\chan_', input_id, '.mat');
    out_str = 'su_table_e_';
    par_str = 'su_pars_e.mat';
end
%
% now load this table. The variable is qt_probs
load(input_table_str);
%
nrows = size(qt_probs, 1); 
% 
tmp = dec2bin(0:nrows-1);
Rows = zeros(size(tmp));
Rows(tmp == '1') = 1;
% Rows is a binary rep. of index of rows in double
% it also can be used as an id for single user channels
% remember that we generate all channels irrespective of demand for them
% 
% single user output channel 
s_probs = zeros(2, size(qt_probs,2));
z_par = zeros(1,nrows - 1);
% Rows = Rows(2:end,:);
% inital string values
sRows = tmp;
% we are using index over rows of this matrix to account for different
% single user channels, clearly all zero case is skipped (no merge)
for i = 1:(size(Rows, 1) - 1)
    % 
    su_id = find(Rows(i + 1,:) == 1);
    % if there are even number of 1s, bitxor is 0, otherwise it's 1
    ix = mod(sum(Rows(:,su_id),2), 2);
    % now ix is the index we need to add up probabilities of 0s and 1s
    s_probs(1,:) = sum(qt_probs(find(~ix),:),1)/2;
    s_probs(2,:) = sum(qt_probs(find(ix),:),1)/2; 
    output_str = strcat(out_str, num2str(bin2dec(sRows(i + 1,:))),'.mat');
    % 
    save(output_str,'s_probs','-v6');
    % sum over all outputs og the channel, return batacharya parameters
    z_par(i) = sum(sqrt(s_probs(1,:) .* s_probs(2,:)));
end
% save batacharya parameters
save(par_str, 'z_par', '-v6');
end

