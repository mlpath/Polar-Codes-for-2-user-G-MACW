function [mac_cond_prob,chan_spec] = mac_awgn(channel_par, P, noise_var, Delta)
%
% P is the power vector for each user transmitted signal
% P = sqrt([P1, P2, ..., Pm])
%
% channel_par is channel gains from users to the receiver
sigma2 = noise_var;
%
% number of MAC users
t = length(channel_par);
%
% the space of all possible message vectors made by t users
u_decimal = 0:(2^t - 1);
u_binary = dec2bin(u_decimal);
% mean vector for gaussian random variables
m = zeros(1,2^t);
% this accomodates conditional probability of each member of current
% output alphabet for all combinations of user messages
%
%
% binary 0 is mapped to +1 and binary 1 is mapped to -1
for i = 1:(2^t)
    % a temporary variable to keep the mapped value of message vector
    x = ones(t,1);
    % change it to -1 if it is a 1
    for j = 1:t
        % there are 2^t possible message vector
        if strcmp(u_binary(i,j),'1')
            % count the number of '1's in transmitted message vector
            x(j) = -1;
        end
    end
    % now make the mean value corresponding to this message vector
    size(P);
    size(channel_par);
    size(x);
    x;
    m(i) = (sqrt(P).*channel_par)' * x;
end
% the interval of integration is based on the most negative and the
% most positive value of mean vector
y = (min(m) - 10*sigma2) :Delta: (max(m) + 10*sigma2);
%
mac_cond_prob = zeros(2^t, length(y));
for i = 1:2^t
    mac_cond_prob(i,:) = exp(-(y - m(i)).^2/(2*sigma2))/(sqrt(2*pi*sigma2));
    size(mac_cond_prob(:,i));
end
%
chan_spec = struct('mean',m,'var_n',sigma2,'out',y);
