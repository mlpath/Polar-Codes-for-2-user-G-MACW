% this function calculates the likelihoods for marginal_gaussian_mac
% and for both users: bob and eve
%
function likelihood = marg_mac_likelihood( channel_par, rx_signals ) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Inputs to this function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% received signal at bob
y = rx_signals(1,:);
% received signal at eve 
z = rx_signals(2,:); 
% transmit power
power = channel_par.tx_power;
% number of users
n_users = length(power);
% define likelihood matrix (for t users in N time slots and both bob and eve, in vertical order)
likelihood = zeros(2* n_users, length(y));
% noise power
sigma2 = channel_par.noise_var;
% channel gains
h = channel_par.bob_gain;
g = channel_par.eve_gain;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Compute marginal conditional probabilities for users at Bob/Eve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compute mean of conditional probabilities given our MACs
tmp = dec2bin(0:2^n_users - 1);
a = ones(size(tmp));
% a is modulated message array, we need to apply power and gains, too
a(tmp == '1') = -1;
a_bob = bsxfun(@times, a, sqrt(power) .* h');
a_eve = bsxfun(@times, a, sqrt(power) .* g');
% the mean of gaussian distribution formed by MACs at Bob and Eve
m_bob = sum(a_bob, 2);
m_eve = sum(a_eve, 2);
%
% next find the marginal probabilities for all users
for i = 1:n_users
    % mariginal prob of user i conditioned on xi = 0
    ix = tmp(:,i) == '0';
    % average over all such combinations
    m1 = m_bob(ix);
    m2 = m_eve(ix);
    % cases where this user has value 1
    m3 = m_bob(~ix);
    m4 = m_eve(~ix);
    %
    s1 = 0;
    s2 = s1;
    s3 = s1;
    s4 = s1;
    for j = 1:sum(ix)
        % bob if user i send 0 
        s1 = s1 + exp(-(y - m1(j)).^2/(2*sigma2));
        % eve if user i send 0 
        s2 = s2 + exp(-(z - m2(j)).^2/(2*sigma2));
        % bob if user i send 1 
        s3 = s3 + exp(-(y - m3(j)).^2/(2*sigma2));
        % eve if user i send 1 
        s4 = s4 + exp(-(z - m4(j)).^2/(2*sigma2));
    end
    % these are marginals conditioned on 0
    s1 = s1/(2^(n_users - 1))/(sqrt(2*pi*sigma2));
    s2 = s2/(2^(n_users - 1))/(sqrt(2*pi*sigma2));
    % and conditioned on 1
    s3 = s3/(2^(n_users - 1))/(sqrt(2*pi*sigma2));
    s4 = s4/(2^(n_users - 1))/(sqrt(2*pi*sigma2));
    %
    likelihood(i,:) = s1/s3;
    likelihood(i + n_users,:) = s2/s4;
end