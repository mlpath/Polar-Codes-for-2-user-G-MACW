% we plot secrecy rate vs. leakage in this code
clc
clear all 
close all
% 
e_1 = .05; 
e_2 = .3;
eta = [1.28e-4, 2.56e-4, 5.12e-4, 7e-4, 1.024e-3, 2.048e-3, 3e-3, 5.2e-3, 7e-3, 8e-3, 1e-2, 1.6e-2, 3.2e-2, 4e-2, 5e-2, 6.4e-2, 8e-2, 1.28e-1 2e-1 4e-1 6e-1 8e-1 1];
%
r_s = zeros(size(e_2)); 
r_leak = r_s; 
max_leak = r_s; 
%
for i = 1:length(eta)
    i 
    [r_s(i), r_leak(i), max_leak(i)] = leakage(e_1, e_2, eta(i)); 
end
%
semilogy(r_s, r_leak, '-k*'); 
save('leak_1_3.mat','r_s','r_leak','max_leak', '-v7.3');
hold on;
semilogy(r_s, max_leak, '-ro'); 

