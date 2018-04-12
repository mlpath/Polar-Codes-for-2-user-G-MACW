% this function computes the mutual information vector as introduced in
% Abbe and Telatar paper: m-user MAC polar codes and matroids
% this is specifically required to find the corresponding matroid to
% constructed channels in arikan procedure
% 
% for this purpose we need to compute conditional mutual information terms
% of the form: I(X[J];Y|X[Jc]) where J is a subset of users and Jc is its
% complement with respect to the set of total users E_m = {1,2,...,m}
%
clc
clear all
bin_level = 8;
%
% number of users 
t = 2; 
% alphabet size 
p = 2; 
% number of possible combinations for vector messages (of all users)
n = p^t; 
% matrix of capacities- for binary alphabet, the length of each column
% equals the number of possible subsets of a set with t elements
cap = zeros(2^bin_level, 2^t);
% the first row is always zero 
cap(:,1) = 0;    
%
for j = 0:(2^bin_level - 1)
    name_str = dec2bin(j);
    name_str = strcat('chan_', name_str, '.mat');
    load(name_str);
    % probability table
    prob_table = qt_probs;
    % remove columns with all elements less than 1e-6
    prob_table(:, all(prob_table < 1e-6,1)) = [];
    %
    % we'll need something to tell us the specific J
    % we don't start from 0 because we know that its mutual information is zero
    inf_sets = dec2bin(0:n-1);
    inf_sets = inf_sets(i_vec_order(n) + 1,:);
    %
    %
    % computation starts here..
    for i = 2:n
        i;
        % store curent set
        J = log2(n) - find(inf_sets(i,:) == '1') + 1;
        % the complement set
        Jc = log2(n) - find(inf_sets(i,:) == '0') + 1;
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % we need a function that maps each row of current conditional
        % proability table p(Y|Xj,Xjc) to the new one p(Y|Xjc) according to
        % indices that are given in Jc
        %
        % pick only those elements which belong to Jc- these are binary form of
        % row indices in new probability table
        row_ix = bin2dec(inf_sets(:,Jc)) + 1;
        %
        % now add
        if ~isempty( Jc )
            % if some user messages are given at channel output
            % first repeat it
            row_IX = row_ix(:,ones(1,size(prob_table,2)));
            row_IX = bsxfun(@plus, row_IX, (0:max(row_ix):(size(row_IX,2) - 1)*max(row_ix)));
            cond_table = reshape(accumarray(row_IX(:), prob_table(:))/(2 ^ length(J)), 2^length(Jc), size(prob_table,2));
            cap(j+1, i) = sum( sum( prob_table .* log2( prob_table ./ cond_table(row_ix,:) ) / n ) );
            
        else
            
            % if it is empty, we only have y at the output (i.e. compute p(y) )
            cond_table = sum( prob_table , 1)/ n;
            sum(cond_table,2);
            cap(j+1, i) = sum( sum( prob_table .* log2( prob_table ./ cond_table(ones(1,n),:) ) / n ) );
            
        end
        %
        % now check to which binary matroid is it closest, easy for
        % 2,3-users, but for larger t we need to first make a bank of
        % binary matroids of that size
        
    end
end
% cap(cap < 1e-8) = 0;
% xlswrite('capacity_2users.xlsx',labels,'A1');
% write in sheet 2
xlswrite('sec_qpsk_4.xlsx',cap, 1, 'A1');

    

