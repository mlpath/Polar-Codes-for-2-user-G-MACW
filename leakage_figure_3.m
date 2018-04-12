% in this code we need to take the samples at different Tx SNRs, set a
% limit for reliability parameter, and compute leakage instead of block
% error probability
clc
clear all 
close all
% 
e_1 = .1; 
e_2 = .1;
% eta = [1.28e-4, 2.56e-4, 5.12e-4, 7e-4, 1.024e-3, 2.048e-3, 3e-3, 5.2e-3, 7e-3, 8e-3, 1e-2, 1.6e-2, 3.2e-2, 4e-2, 5e-2, 6.4e-2, 8e-2, 1.28e-1 2e-1 4e-1 6e-1 8e-1 1];
eta = .005; 
% 
SNRs = {'5', '6_7', '8_5', '10', '11_2', '13_5', '15', '20', '25', '30', '40'};
%
r_leak = zeros(size(SNRs)); 
max_leak = r_leak;
%
for i = 1:length(SNRs)
    SNRs{i}
    filename = strcat('secrecy_graph_', SNRs{i}, '.xlsx');
    foldername = strcat('debug_channels_10_p', SNRs{i});
    % now call leakage function
    [~,r_leak(i), max_leak(i)] = leakage(0.1,0.1, eta, filename, foldername);
    r_leak(i)
end
% 
%
P = [.5,.678, .8536, 1, 1.1253, 1.353, 1.5, 2, 2.5, 3, 4];
sigma2 = .97865^2; 
P = 10*log10(P/sigma2);
%
semilogy(P, r_leak, '-k*'); 
save('leak_snr_005.mat','P','r_leak','max_leak', '-v7.3');
% hold on;
% semilogy(r_s, max_leak, '-ro'); 

