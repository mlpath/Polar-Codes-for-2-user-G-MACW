% for small files where there is no need for block processing
% wrote it specifically for discretizing the continuous input raw channel
%
clc
clear all
% function e_degrading(input_filename, fidelity_par)
fidelity_par = 20; 
input_filename = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   inputs required for degrading operation               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
if isempty(input_filename)
    Delta = 1e-4;
    sigma2 = .97865^2;
    gains = [1 1]';
    power = [2 2]';
    fileid = num2str(gains(2)*10);
    % conditional probability of 2-user mac_awgn channel
    [A,chan_spec] = mac_awgn(gains, power, sigma2, Delta);
end
disp('done with continous channel');
% initialize list of labels and list of probabilities
list_b = [];
list_p = [];
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
output_str = strcat('C:\MonaHajimomeni\mfiles\debug_mac_v2\raw_channel_',fileid,'.mat')
save( output_str, 'qt_probs', 'qt_bins' , '-v7.3');
% we might not need any object as we are saving the entire file at once
% but for now I leave it as it is
out_chan = matfile(output_str, 'Writable', true);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Start the e_degrading_merge function %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
    tic
    % conditional probability of plus channel
%     A = A / nrow;    
    % posterior of plus channel- note how we calculate marginal by sum
    post_prob_p = A./repmat(sum(A,1),SZ_row,1);   
    
    sum(post_prob_p,1)
    -post_prob_p.*log2(post_prob_p);
    % this is the output of binning, rewrite posterior variables
    post_prob_p = e_binning(post_prob_p, fidelity_par) - 1;
    sum(sum(isnan(post_prob_p)));
    disp('done with binning');
    % make them scalar
    post_prob_p = bsxfun(@times, post_prob_p, b_max.^ex_b);
    % update list of existing labels for plus channel
    list_b( l_ptr:(SZ_col + l_ptr - 1),1) = sum(post_prob_p , 1);
    list_b';
    % sort labels
    [list_b, ix] = sort(list_b,'ascend');
    % update list of corresponding probabilities        
    list_p( : , l_ptr:(SZ_col + l_ptr - 1)) = A;
    list_p = list_p(:, ix);
    % distinct labels in plus channel
    [~, ia, ic] = unique(list_b);
    list_b = list_b(ia,:);
    % for continuous channels we need to integrate over probability density
    m = chan_spec.mean;
    y = chan_spec.out;
    sigma2 = chan_spec.var_n;
    %
f = @(x, c)(exp((-(x - c).^2)/(2*sigma2))/sqrt(2*pi*sigma2));
for i = 1:SZ_row 
    for j = 1:length(ia)
        [i j length(ia)]
        % integrate from this point to next one
        if j < length(ia)
            list_p(i,j) = integral(@(x) f(x,m(i)), y(ia(j)), y(ia(j+1)));
        else
            list_p(i,j) = integral(@(x) f(x,m(i)), y(ia(j)), y(SZ_col));
        end
    end
end
list_p(:,max(size(ia))+1:end) = [];
    
%     for j = 1:length(m)
%             f = @(x)(exp((-(x - m(j)).^2)/(2*sigma2))/sqrt(2*pi*sigma2));
%             [ix_a(1),ix_a(end)];
%             if ix_a(end) < length(list_a)
%                 list_p(i,j) = quad(f,y(ix_a(1)),y(ix_a(end)+1));
%             else
%                 list_p(i,j) = quad(f,y(ix_a(1)),y(ix_a(end)));
%             end
%     end
        
    
    
    % merge (i.e. add up) the probabilities of identical labels in the list
%     ix = reshape(bsxfun(@plus, repmat(ic', SZ_row, 1), (0:max(ic):(SZ_row - 1)* max(ic))')', size(list_p,2) * SZ_row, 1);
%     list_p = reshape(accumarray(ix, reshape(list_p',size(list_p,2) * SZ_row, 1)), length(ia), SZ_row)'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  End of degrading_merge operation for small size channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
sum(list_p,2)'
% save results (including bins and probs)
out_chan.qt_probs(:,1:size(list_p,2)) = list_p;
out_chan.qt_bins(1:size(list_b,1),1) = list_b;
%%%