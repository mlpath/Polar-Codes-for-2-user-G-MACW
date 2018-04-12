% this will compute channel level likelihoods
% 
function likelihood = lr_chan_level(channel_type,channel_par_1, y)
%
% y is the vector of received signal samples
% channel_type is the specific channel type whose conditional probability
% must be computed in likelihood
% channel_par is the associated parameter of the channel
% 
% 
% length of input
l_in = length(y);
% yes! this is the channel level! 
likelihood = zeros(1,l_in);
%
channel_type = lower(channel_type);
%
if strcmp(channel_type,'erasure' )
    % erasure channel has three outputs: 0, 1, e
    e = channel_par_1;                                % erasure probability
    conditional_prob = [1 - e, 0; 0, 1 - e; e, e] ;
    % 
%     likelihood(y == 0) = conditional_prob(1,1)/conditional_prob(1,2); 
%     likelihood(y == 1) = 0; 
%     likelihood(y == -1) = conditional_prob(3,1)/conditional_prob(3,2);
    likelihood(y == 0) = 1e10;
    likelihood(y == 1) = eps; 
    likelihood(y == -1) = 1;
    %
elseif strcmp(channel_type,'bmc')
%     fprintf('hi \n');
    % binary memoryless channel has two outputs: 0 , 1
    e = channel_par_1;                                  % error probability
    conditional_prob = [1 - e, e; e, 1 - e];
    % 
    likelihood(y == 0) = conditional_prob(1,1)/conditional_prob(1,2); 
    likelihood(y == 1) = conditional_prob(2,1)/conditional_prob(2,2); 
    %
elseif strcmp(channel_type,'awgn')
    % we assume that awgn has zero mean
    sigma2 = channel_par_1;                          % variance of gaussian
    % 
    for i= 1:l_in
%         likelihood(i) = exp(-(y(i) - 0).^2/(2*sigma2))/exp(-(y(i) - 1).^2/(2*sigma2));
% we're accounting for Eb = -1 as the power of antipodal signaling
        likelihood(i) = exp(-(y(i) - 1).^2/(2*sigma2))/exp(-(y(i) + 1).^2/(2*sigma2));

    end
else
    fprintf('\n bad channel type! choose from (erasure/bmc/awgn) \n');
end


