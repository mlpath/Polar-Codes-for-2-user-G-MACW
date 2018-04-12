% in this code we try to generate some informative plots using the
% upper/lower bound obtained in batacharya_mac_v2.m
% 
% the observation we made there is that starting from very small epsilon
% values, 
% 
%       (0) effect of choosing a very large epsilon: almost all the
%       subchannels are selected. on the one hand, almost all of the
%       subchannels are errorneous, because we have rounded time to a
%       matrix that is very different from its true capability. therefore,
%       just by filtering out high upper bound errors on Bob channels, most
%       of the subchannels are rejected. we may have a low rate or zero
%       secrecy rate as a result. 
%
%       (1) number of secure time slots (i.e. secrecy rate) increases
%
%       (2) the minimum value of lower bound on secure time slots at Eve
%       decreases (secrecy challenge)
%
%       (3) the upper bound to block error probability on Bob channel
%       increases, too. (reliability challenge)
%
%       (4) what is learnt from this experiment is that we have a spectrum 
%       of options for security/reliability. as part (0) shows, it is clear
%       that we need to pick a small epsilon for our experiment to be
%       meaningful. however, it still remains to tell how small it needs to
%       be? therefore, we start from very small values, monitoring the
%       upper bound for error on Bob (reliability) and minimum value of
%       lower bound on Eve (secuirty), up to moderate values. the choice of
%       a design point depends on the level of security and reliability
%       that is desired, so, we may say that we can have .33 bps secure
%       transmission but the error rate is bounded by 0.05, and the minimum 
%       error that Eve experience is .1, whereas we may have .07 bps secure
%       rates with an error bounded by 1e-3 and minimum Eve error of .45.
%   
%       now what we need is to devise more experiment to produce figures.
%       think about how. 
% 
%       what is our output parameters? secure rate, upper bound on block
%       error probability at Bob, and minimum lower bound in Eve channels
% 
clc
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%  all input parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% block_length = {'5', '6_7', '8_5', '10', '11_2', '13_5', '15', '20', '25', '30', '40'};
%
block_length = '256';
%%%%%%%%%%%%%%%%%%%%%%% error with a fixed eta %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eta = .0005;
% epsilon = [.0001:2e-4:1e-3,3e-3:2e-3:.009,.01:.05:.9];
%
%%%%%%%%%%%%%%%%%%%%%%% error with a fixed epsilon %%%%%%%%%%%%%%%%%%%%%%%%
epsilon = 0.05;
% eta = 0.002;
% 1024
% eta = [1e-6, 2e-6, 4e-6, 8e-6, 1.6e-5, 3.2e-5, 6.4e-5, 1.28e-4, 2.56e-4, 5.12e-4, 1.024e-3, 2.048e-3, 4e-3, 8e-3, 1.6e-2, 3.2e-2, 6.4e-2, 1.28e-1];
% eta = [1e-6 2e-6 4e-6 8e-6 1e-5 2e-5 4e-5 8e-5 1e-4, 2e-4, 4e-4, 8e-4, 1.6e-3, 3.2e-3, 6.4e-3, 1.28e-2, 2.56e-2, 5.12e-2 1e-1 2e-1 4e-1 5e-1 8e-1 9e-1 1.1 1.3 2];
%
% 512
% eta = [1e-4 1.2e-4 1.7e-4 2e-4 4e-4 8e-4 1.6e-3 3.2e-3 6.4e-3 1.28e-2 2.56e-2 5.12e-2 1e-1 2e-1 4e-1];
% eta = [1.28e-4, 2.56e-4, 5.12e-4, 1.024e-3, 2.048e-3, 4e-3, 8e-3, 1.6e-2, 3.2e-2, 6.4e-2, 1.28e-1];
% eta = [1e-4, 2e-4, 4e-4, 8e-4, 1.6e-3, 3.2e-3, 6.4e-3, 1.28e-2, 2.56e-2, 5.12e-2 1e-1 2e-1 4e-1 5e-1 9e-1 1.1 1.3 2];
%
% 256
% eta = [2e-3 3e-3 1.5e-2 2e-2 4e-2 2e-1 1];
eta = [1e-5 1.28e-4, 2.56e-4, 5.12e-4, 1.024e-3, 2.048e-3, 3e-3, 5.2e-3, 5.6e-3, 7e-3, 8e-3, 1e-2, 1.6e-2, 3.2e-2, 4e-2, 5e-2, 6.4e-2, 8e-2, 1.28e-1 2e-1 4e-1 1.35 2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(epsilon, 2) == size(eta, 2)
    % this condition only happens when both eta and epsilon are scalars
    ix = 0;
else
    [~, ix] = min([size(epsilon, 2), size(eta, 2)]);
end
%
if ix == 1
    % we're now here
    epsilon_str = num2str(epsilon);
    par_id = strcat('_eps_',epsilon_str(3:end));
    % 
    M = size(eta, 2);
    file_id = strcat('btch_', block_length, par_id, '_2.mat');
elseif ix == 2
    eta_str = num2str(eta);
    par_id = strcat('_eta', eta_str(3:end));
    M = size(epsilon, 2);
    file_id = strcat('btch_', block_length, par_id, '.mat');
else
    epsilon_str = num2str(epsilon); 
    eta_str = num2str(eta);
    par_id = strcat('snr_eta_', eta_str(3:end), '_', epsilon_str(3:end));
    % 
    M = size(block_length,2);
    file_id = strcat('btch_', par_id, '.mat');
end
%
%
b_upper = zeros(1,M);
e_lower = b_upper; 
b_lower = b_upper;
e_lower_min = b_upper;
r_s = b_upper;
%
for i = 1:M
    i
%     file_id = strcat('btch_', block_length, par_id, '.mat');
    %
    if ix == 1
        [secrecy_rate, time_slots, b_bound, e_lower_min(i), e_lower(i)] = batacharya_mac_v4(epsilon, eta(i), block_length);
    elseif ix == 2
        [secrecy_rate, time_slots, b_bound, e_lower_min(i), e_lower(i)] = batacharya_mac_v2(epsilon(i), eta, block_length);
    else
        [secrecy_rate, time_slots, b_bound, e_lower_min(i), e_lower(i)] = batacharya_mac_v4(epsilon, eta, block_length{i});
    end
    % reliability parameters
    b_upper(i) = b_bound(2); 
    b_lower(i) = b_bound(1);
    % security parameters
    % e_lower_min(i);
    %
    r_s(i) = secrecy_rate;
end
r_s*256
b_upper
e_lower
e_lower_min

% P = [.5,.678, .8536, 1, 1.1253, 1.353, 1.5, 2, 2.5, 3, 4];
% sigma2 = .97865^2; 
% P = 10*log10(P/sigma2);
% plot(P, r_s, '-^b'); 
% % % 
% % reliability curve
% figure(1);
% semilogy(r_s, b_upper, '-*');
% % 
% % security curve
% figure(2);
% semilogy(r_s, e_lower_min, '-*');
% %
% save(file_id, 'r_s', 'b_upper', 'e_lower_min', 'epsilon', 'eta', '-v6'); 

%
% for i = 1:length(eta)
%     i
%     [secrecy_rate, time_slots, b_bound, e_lower_min(i), e_lower(i)] = batacharya_mac_v2(epsilon, eta(i), block_length);
%     % reliability parameters
%     b_upper(i) = b_bound(2); 
%     % security parameters
%     % e_lower_min(i);
%     %
%     r_s(i) = secrecy_rate;
% end

% ix_rs = r_s ~= 0; 
% % reliability curve
% figure(1);
% semilogy(r_s(ix_rs), b_upper(ix_rs), '-*');
% % 
% % security curve
% figure(2);
% plot(r_s(ix_rs), e_lower_min(ix_rs), '-*');
% hold on; 
% plot(r_s(ix_rs), e_lower(ix_rs),'-ok');
%
save(file_id, 'r_s', 'b_upper', 'b_lower', 'e_lower_min', 'e_lower', 'epsilon', 'eta', '-v6'); 
%