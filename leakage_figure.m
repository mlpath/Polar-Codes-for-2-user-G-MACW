% we plot secrecy rate vs. leakage in this code
clc
clear all 
close all
% 
e_1 = .03; 
e_2 = [e_1:.02:1];
eta = 1e-3;
%
r_s = zeros(size(e_2)); 
r_leak = r_s; 
Del = r_s; 
%
for i = 1:length(e_2)
    i 
    [r_s(i), r_leak(i), Del(i)] = leakage(e_1, e_2(i), eta); 
end
%
semilogy(r_s, r_leak, '-k*'); 
hold on;
semilogy(r_s, Del, '-ro'); 

