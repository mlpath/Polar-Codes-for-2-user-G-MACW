% this is a file to compute conditional probability density for QPSK 
% 
% function nb_cond_prob = non_binary_awgn( noise_var, npoints)
clc;
clear all; 
% Attention: You have already accounted for half energy in each dimension
% in the pdf formula. So, this parameter is N0, the total noise power 
noise_var = 2^2; 
npoints = 5e3; 
% quantization factor
fidelity_par = 30;
% file identifier 
fileid = 'qpsk_2_mu30';
% fileid = 'qpsk_047';
% output address 
out_addr = 'D:\April-queens\mfile\backup_final\mat-files\';
% noise variance 
sigma2 = noise_var; 
% number of modulation symbols determines the rows of conditional
% probability
% the most simple and inefficient way..
% m = [-1, -1; -2, 2; 1, 2; 1, -2];
% normalize it so that each symbol of the constellation has the average
% power of 1 (statistically)
m = [-1, -1; -1, 1; 1, 1; 1, -1]/sqrt(2);
% maximum value for 
m_min = min(min(m));
m_max = max(max(m)); 
% 
% this is the range of output symbols 
y_i = linspace(m_min - 5*sqrt(sigma2), m_max + 5*sqrt(sigma2), npoints);
% because of symmetry
y_q = fliplr(y_i); 
%
nb_cond_prob = zeros(length(y_i), length(y_q), 4);
for l = 1:4
for i = 1:length(y_i)
    for j = 1:length(y_q)
        tmp = exp(-((y_i(i) - m(l,1))^2 +(y_q(j) - m(l,2))^2)/sigma2)/(pi * sigma2);
        nb_cond_prob(i,j,l) = tmp;
    end
    fprintf('# %d out of %d', j+(i-1)*npoints, npoints*npoints );
    fprintf('\n');
end
end
% ok, assume we have the correct samples from conditional probability
% now compute the posterior probability
for l = 1:4
        tmp = nb_cond_prob(:,:,l) ./ sum(nb_cond_prob, 3);
        % save all the posterior probabilities in a row
        nb_post_prob(l,:) = reshape(tmp, 1, size(nb_cond_prob, 1)*size(nb_cond_prob, 2));
        nb_cond(l,:) = reshape(nb_cond_prob(:,:,l), 1, size(nb_cond_prob, 1)*size(nb_cond_prob, 2));
end
%
% 
% now what? we save this, or run e_degrading directly from here
disp('done with continuous channel');
% initialize list of labels and list of probabilities
list_b = [];
list_p = [];
% this should be  4 x something matrix 
A = nb_cond; 
% number of distinct vector messages (i.e. for 2 users with binary
% alphabet we have 2^2 = 4 distinct vectors)
SZ_row = size(A , 1);
SZ_col = size(A , 2);
l_ptr = 1;
% these are used for scalarization, i.e. 
% sum{b_i * b_max ^ i}_{i = 0}^{nrow - 1} (like changing binary to decimal)
% we may reconsider this method later, but it's working for now
b_max = 2*ceil(fidelity_par /(log(2)*exp(1)));
ex_b = (0:SZ_row -1)';
% this will change after we rewrite the main control routine
% initialize the output channel where results are stored 
qt_probs = zeros(4,0);
qt_bins = zeros(1,0);
% save output results (i.e. degraded versions of arikan plus and minus
% channels) separately
output_str = strcat(out_addr, 'raw_channel_',fileid,'.mat');
% save( output_str, 'qt_probs', 'qt_bins' , '-v7.3');
% we might not need any object as we are saving the entire file at once
% but for now I leave it as it is
% out_chan = matfile(output_str, 'Writable', true);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Start the e_degrading_merge function %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
    % conditional probability of plus channel
%     A = A / nrow;    
    % posterior of plus channel- note how we calculate marginal by sum
%     post_prob_p = A./repmat(sum(A,1),SZ_row,1);   
    post_prob_p = nb_post_prob;
    
    sum(post_prob_p,1);
    -post_prob_p.*log2(post_prob_p);
    % this is the output of binning, rewrite posterior variables
    post_prob_p = e_binning(post_prob_p, fidelity_par) - 1;
    sum(sum(isnan(post_prob_p)));
    disp('done with binning');
    % make them scalar
    post_prob_p = bsxfun(@times, post_prob_p, b_max.^ex_b);
    % update list of existing labels for plus channel
    list_b(  l_ptr:(SZ_col + l_ptr - 1), 1) = sum(post_prob_p , 1);

%     tmp_rect = reshape(list_b, npoints, npoints);
%     imshow(tmp_rect/max(max(tmp_rect)));

    % distinct labels in plus channel
    [~, ia, ic] = unique(list_b);
    % compute the accumulation of all conditional prbabilities
    Delta = y_i(end) - y_i(1);
    S = zeros(4, length(ia));
    for l = 1:4
        S(l,:) = accumarray(ic, A(l,:))*(Delta /npoints)* (Delta/npoints);
    end
    sum(S,2);
    tmp = list_b(ia);
    % 
    qt_probs = S;
    save(output_str, 'qt_probs'); 
    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for continuous channels we need to integrate over probability density
%     m = chan_spec.mean;
%     y = chan_spec.out;
%     sigma2 = chan_spec.var_n;
%     %
% this is the exact integration, using the analytical forumla. the problem
% is that I don't have any way to extract the boundaries of lables and even
% in that case I don't know how to integrate over it (quad has only
% % rectabgular upper/lower bounds)
% f = @(x, y, c1, c2)(exp(-((x - c1).^2 + (y - c2).^2)/(sigma2))/(pi*sigma2));
% for i = 1:4 
%     for j = 1:length(tmp)
%         [I(j),J(j)]
%         [y_i(1), y_q(I(j)),y_i(100), y_q(J(j))]
%         % the negative sign accounts for flipping upper and lower bounds in
%         % the integration. example: (-6,-2) to (6, -6) gives xmin = 06,
%         % ymin = -2, xmax = 6, ymax = -6, which contradicts ymax > ymin,
%         % and does the sign change happens. 
%         tmp_p(i,j) = -dblquad(@(x,y,c1,c2) f(x,y, m(i,1),m(i,2)), y_i(1), y_i(100), y_q(I(j)), y_q(J(j)));
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
