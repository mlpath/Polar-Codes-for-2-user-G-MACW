% this function calculates the likelihoods for marginal_gaussian_mac
% and for both users: bob and eve
%
function likelihood = marg_mac_likelihood( chan_filename, rx_signal ) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Inputs to this function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% received signal at bob

%
load(chan_filename); 
% conditional probability table for bob channel
q = qt_probs; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Compute marginal conditional probabilities for users at Bob/Eve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
tmp = dec2bin(0:size(q, 1)-1); 
n_users = log2(size(q,1));
lr_ref = zeros(n_users, size(q,2));
%
% next find the marginal probabilities for each user
for i = 1:n_users
    % mariginal prob of user i conditioned on xi = 0
    ix = tmp(:,i) == '0';
    % marginal prob of bob if 0 is transmitted for user i over marginal
    % prob if 1 is transmitted
    lr_ref(i,:) = sum(q(find(ix),:),1)./sum(q(find(~ix),:),1);
end
% define likelihood matrix (for t users in N time slots and both bob and eve, in vertical order)
likelihood = lr_ref(:,rx_signal); 
