% summarize the result of our algorithm 
clc
clear all
close all
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Fig 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure; 
load btch_1024_eps_05_2
ix = b_upper > 0.5;
b_upper(ix) = [];
b_lower(ix) = [];
e_lower(ix) = [];
r_s(ix) = [];
% 
ix_rs = r_s ~= 0; 
b_upper = b_upper(ix_rs);
b_lower = b_lower(ix_rs);
e_lower = e_lower(ix_rs); 
e_lower_min = e_lower_min(ix_rs);
r_s = r_s(ix_rs);
r_s_1024 = r_s;
e_lower_1024 = e_lower;
e_lower_min_1024 = e_lower_min;
b_lower_1024 = b_lower;
b_upper_1024 = b_upper;

%
% semilogy(r_s(ix_rs), b_upper(ix_rs), '^-k'); 
% semilogy(r_s, b_upper, '-+', 'LineWidth', 2); hold on; semilogy(r_s, b_lower, '--+', r_s, e_lower, '-r+');
semilogy(r_s, b_upper, '-+k');
hold on; 
load btch_512_eps_05_2
ix = b_upper > 0.5;
b_upper(ix) = [];
b_lower(ix) = [];
e_lower(ix) = [];
r_s(ix) = [];
% 
ix_rs = r_s ~= 0; 
b_upper = b_upper(ix_rs);
b_lower = b_lower(ix_rs);
e_lower = e_lower(ix_rs); 
e_lower_min = e_lower_min(ix_rs);
r_s = r_s(ix_rs);
% semilogy(r_s(ix_rs), b_upper(ix_rs), '*-k'); 
% semilogy(r_s, b_upper, '-o', 'LineWidth', 2); hold on;semilogy(r_s, b_lower, '--o', r_s, e_lower, '-ro');
semilogy(r_s, b_upper, '-ok');%, r_s, e_lower, '-or');
r_s_512 = r_s;
e_lower_512 = e_lower;
e_lower_min_512 = e_lower_min;
b_lower_512 = b_lower;
b_upper_512 = b_upper;
hold on;
load btch_256_eps_05_2
ix = b_upper > 0.5;
b_upper(ix) = [];
b_lower(ix) = [];
e_lower(ix) = [];
r_s(ix) = [];
% 
ix_rs = r_s ~= 0; 
b_upper = b_upper(ix_rs);
b_lower = b_lower(ix_rs);
e_lower = e_lower(ix_rs);
e_lower_min = e_lower_min(ix_rs);
r_s = r_s(ix_rs);
%
r_s_256 = r_s;
e_lower_256 = e_lower;
e_lower_min_256 = e_lower_min;
b_lower_256 = b_lower;
b_upper_256 = b_upper; 

% semilogy(r_s(ix_rs), b_upper(ix_rs), '+-k'); 
% 
% semilogy(r_s, b_upper, '-*', 'LineWidth', 2); hold on; semilogy(r_s, b_lower, '--*', r_s, e_lower, '-r*');
semilogy(r_s, b_upper, '-*k');%,r_s, e_lower, '-*r');
legend('N = 1024', 'N = 512', 'N = 256');
%
hold on; 
semilogy(r_s_1024, e_lower_1024, '-+r', 'MarkerSize', 4.5);
hold on; 
semilogy(r_s_512, e_lower_512, '-or', 'MarkerSize', 4.5);
hold on; 
semilogy(r_s_256, e_lower_256, '-*r', 'MarkerSize', 4.5);
%
% hold on; 
axis([0 0.1018 10^-6 1]);
% axh = axes; 
% set(axh, 'XGrid', 'off', 'YGrid', 'on');
% 
hold on; 
semilogy(r_s_256, b_lower_256, '--*k');
hold on; 
semilogy(r_s_512, b_lower_512, '--ok');
hold on;
semilogy(r_s_1024, b_lower_1024, '--+k');
%
% hold on; 
% load btch_1024_eps_5
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), e_lower(ix_rs), '^--b', 'MarkerSize', 4.5); 
% hold on; 
% load btch_512_eps_5
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), e_lower(ix_rs), '*--b', 'MarkerSize', 4.5); 
% hold on; 
% load btch_256_eps_5
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), e_lower(ix_rs), '+--b', 'MarkerSize', 4.5); 
% % 
title('Reliability for Bob (upper bound) and Eve (lower bound) at fixed \epsilon = 0.5');
ylabel('Bounds on block error probability'); 
xlabel('Achieved secrecy rate (bps)'); 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% how can I show the performance for different SNRs?
% plot it for other SNRs..
% figure;
% load btch_6_7_eps05
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), b_upper(ix_rs), '^-k'); 
% hold on; 
% %
% load btch_8_5_eps05
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), b_upper(ix_rs), '*-k'); 
% hold on;
%
% load btch_10_eps05
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), b_upper(ix_rs), '+-k'); 
% hold on;
% %
% load btch_11_2_eps05
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), b_upper(ix_rs), '.-k'); 
% hold on;
% %
% load btch_13_5_eps05
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), b_upper(ix_rs), 's-k'); 
% hold on;
% %
% load btch_15_eps05
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), b_upper(ix_rs), 'o-k'); 
% hold on;
% %
% load btch_20_eps05
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), b_upper(ix_rs), 'd-k'); 
% hold on;
% 
% figure
% load btch_1024_eps_3_05
% ix_rs = r_s ~= 0;
% semilogy(r_s(ix_rs), b_upper(ix_rs), '.-k');
% hold on; 
% load btch_1024_eps_5
% ix_rs = r_s ~= 0;
% semilogy(r_s(ix_rs), b_upper(ix_rs), '+-b');
% hold on; 
% legend('\epsilon = 0.05','\epsilon = .5'); 
% % semilogy(r_s(ix_rs), e_lower_min(ix_rs), '.-k');
% % hold on; 
% load btch_512_eps_3_05
% ix_rs = r_s ~= 0;
% semilogy(r_s(ix_rs), b_upper(ix_rs), '.-k'); 
% hold on;
% % semilogy(r_s(ix_rs), e_lower_min(ix_rs), '.-k');
% % hold on;
% load btch_256_eps_3_05
% ix_rs = r_s ~= 0;
% semilogy(r_s(ix_rs), b_upper(ix_rs), '.-k');
% hold on; 
% % semilogy(r_s(ix_rs), e_lower_min(ix_rs), '.-k');
% % % 
% % hold on; 
% 
% % semilogy(r_s(ix_rs), e_lower_min(ix_rs), '+-b');
% % hold on; 
% load btch_512_eps_5
% ix_rs = r_s ~= 0;
% semilogy(r_s(ix_rs), b_upper(ix_rs), '+-b');
% hold on; 
% % semilogy(r_s(ix_rs), e_lower_min(ix_rs), '+-b');
% % hold on; 
% load btch_256_eps_5
% ix_rs = r_s ~= 0;
% semilogy(r_s(ix_rs), b_upper(ix_rs), '+-b');
% hold on; 
% % semilogy(r_s(ix_rs), e_lower_min(ix_rs), '+-b');
% title('Reliability vs. secure rate at \epsilon = 0.05 and \epsilon = 0.5'); 
% ylabel('Upper bound of block error probability'); 
% xlabel('Achieved secrecy rate (bps)'); 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Fig 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = [.5,.678, .8536, 1, 1.1253, 1.353, 1.5, 2, 2.5, 3, 4];
sigma2 = .97865^2; 
% P = 10*log10(P/sigma2);
P = 10*log10(P);
M = length(P);
% remove '40'
P(M) = [];
% 
load btch_snr_eta_05_1
figure;
r_s = r_s(1:M-1);
plot(P, r_s, '^k');
hold on; 
load btch_snr_eta_01_1
r_s = r_s(1:M-1);
plot(P, r_s, '*k');
hold on; 
load btch_snr_eta_005_1
r_s = r_s(1:M-1);
plot(P, r_s, 'ok');
hold on;
load btch_snr_eta_002_1
r_s = r_s(1:M-1);
plot(P, r_s, 'sk');
hold on;
legend('\eta = 0.05', '\eta = 0.01', '\eta = 0.005', '\eta = 0.002');
ylabel('Achieved secrecy rate (bps)'); 
xlabel('Tx SNR (db)');
title('Secrecy rate vs. Tx SNR as a function of reliability parameter, \eta');
load fit_256_1
x = -3.5:.1:4.7;
c = fit_01_2.coeff;
d = length(c) - 1; 
y = polynomial_fit(x, d, c);
plot(x, y, '--k', 'LineWidth', 2);
hold on;
c = fit_05_2.coeff;
d = length(c) - 1;
y = polynomial_fit(x, d, c); 
plot(x, y, '--k', 'LineWidth', 2);
hold on; 
c = fit_005_2.coeff;
d = length(c) - 1;
y = polynomial_fit(x, d, c); 
plot(x, y, '--k', 'LineWidth', 2);
hold on; 
c = fit_002_2.coeff;
d = length(c) - 1;
y = polynomial_fit(x, d, c); 
plot(x, y, '--k', 'LineWidth', 2);
hold on; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Fig 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
load btch_snr_eta_05_1
ix_rs = r_s ~= 0; 
semilogy(P(ix_rs), e_lower(ix_rs), '^k');
hold on;
load btch_snr_eta_01_1
ix_rs = r_s ~= 0; 
semilogy(P(ix_rs), e_lower(ix_rs), 'sk', 'MarkerSize', 5);
hold on; 
load btch_snr_eta_005_1
ix_rs = r_s ~= 0;
semilogy(P(ix_rs), e_lower(ix_rs), '+k', 'MarkerSize', 5);
%
legend('\eta = 0.05', '\eta = 0.01', '\eta = 0.005'); 
hold on;
%
load btch_snr_eta_05_1
ix_rs = r_s ~= 0; 
semilogy(P(ix_rs), e_lower_min(ix_rs), '.k');
hold on; 
load btch_snr_eta_01_1
ix_rs = r_s ~= 0; 
semilogy(P(ix_rs), e_lower_min(ix_rs), '*k', 'MarkerSize', 4);
hold on;
load btch_snr_eta_005_1
ix_rs = r_s ~= 0; 
semilogy(P(ix_rs), e_lower_min(ix_rs), '^k', 'MarkerSize', 4); 
% 
% fitting
load fit_256_1_eve.mat
x = -3:.1:5;
c = fit_005_eta_2.coeff;
d = length(c) - 1; 
y = polynomial_fit(x, d, c);
plot(x, y, '-k', 'LineWidth', 1);
hold on;
c = fit_01_eta_2.coeff;
d = length(c) - 1;
y = polynomial_fit(x, d, c); 
plot(x, y, '-k', 'LineWidth', 1);
hold on; 
c = fit_05_eta_2.coeff;
d = length(c) - 1;
y = polynomial_fit(x, d, c); 
plot(x, y, '-k', 'LineWidth', 1);
hold on; 
%
c = fit_min_2.coeff;
d = length(c) - 1;
y = polynomial_fit(x, d, c); 
plot(x, y, '--k', 'LineWidth', 1);
%
axis([-3 5 0.01 1]); 
ylabel('Bounds on error probability over secure time slots'); 
xlabel('Tx SNR (db)');
title('Max/min of bounds on Eve error probability on secure slots, \epsilon = 0.1'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc; close all
% figure
% % want to show the minimum value of lower bound on subchannels
% load btch_1024_eps_3_05
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), e_lower_min(ix_rs), '--*b', 'MarkerSize', 4);
% hold on; 
% load btch_1024_eps_5
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), e_lower_min(ix_rs), '--^k');
% hold on; 
% load btch_1024_eps_3_05
% semilogy(r_s(ix_rs), e_lower(ix_rs), '-*b', 'MarkerSize', 4);
% hold on; 
% load btch_1024_eps_5
% semilogy(r_s(ix_rs), e_lower(ix_rs), '-^k');
% title('Max (solid) and min (dashed) value of lower bound of error probability at secure time slots for Eve'); 
% ylabel('Bounds on block error probability');
% xlabel('Achieved secrecy rate (bps)'); 
% legend('\epsilon = 0.05', '\epsilon = 0.5');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Fig 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
figure;
load btch_10_eps_03
ix_rs = r_s ~= 0; 
semilogy(r_s(ix_rs), e_lower(ix_rs), '-*k', 'MarkerSize', 6);
hold on; 
% load btch_10_eps_05
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), e_lower(ix_rs), '->b', 'MarkerSize', 4);
% hold on; 
% load btch_10_eps_08
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), e_lower(ix_rs), '-sb', 'MarkerSize', 5);
% hold on; 
load btch_10_eps_1
ix_rs = r_s ~= 0; 
semilogy(r_s(ix_rs), e_lower(ix_rs), '-sk', 'MarkerSize', 5);
hold on; 
load btch_10_eps_3
ix_rs = r_s ~= 0; 
semilogy(r_s(ix_rs), e_lower(ix_rs), '-^k', 'MarkerSize', 5);
hold on; 
% load btch_10_eps_5
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), e_lower(ix_rs), '-ob', 'MarkerSize', 5);
% hold on; 
legend('\epsilon = 0.03', '\epsilon = 0.1', '\epsilon = 0.3');
% 
load btch_10_eps_03
ix_rs = r_s ~= 0; 
semilogy(r_s(ix_rs), e_lower_min(ix_rs), '-*b', 'MarkerSize', 6);
hold on; 
% %
% load btch_10_eps_05
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), e_lower_min(ix_rs), '-->k', 'MarkerSize', 4);
% hold on; 

% load btch_10_eps_08
% ix_rs = r_s ~= 0; 
% semilogy(r_s(ix_rs), e_lower_min(ix_rs), '--sk', 'MarkerSize', 5);
% hold on; 
%
load btch_10_eps_1
ix_rs = r_s ~= 0; 
semilogy(r_s(ix_rs), e_lower_min(ix_rs), '-sb', 'MarkerSize', 5);
hold on; 
%
load btch_10_eps_3
ix_rs = r_s ~= 0; 
hold on; 
semilogy(r_s(ix_rs), e_lower_min(ix_rs), '-^b', 'MarkerSize', 5);
% % 
% load btch_10_eps_5
% ix_rs = r_s ~= 0;
% hold on; 
% semilogy(r_s(ix_rs), e_lower_min(ix_rs), '--dk', 'MarkerSize', 5);
% %
% 
axis([0 .1018 0.01 1]); 
ylabel('Bounds on Eve error probability over secure time slots'); 
xlabel('Achieved secrecy rate (bps)');
title('Max/min of bounds on Eve error probability on secure slots, N = 256, SNR = 0.19dB'); 




