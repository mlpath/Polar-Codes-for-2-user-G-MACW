% this code produces a simple graph showing the number of secure time slots
% versus transmission power 
%
clc
clear all 
close all
%
id = [5, 10, 15, 20, 25, 30, 40];
P = id/10;
sum_secrecy_rate = zeros(size(id));
bob_sheet = 1; 
eve_sheet = 2; 
% 
Epsilon = .058;
% block length
N = 256; 
%
for i = 1:length(id)
    filename_str = strcat('secrecy_graph_', num2str(id(i)), '.xlsx');
    [~, user_grid, ~] = find_secure_v2(filename_str, Epsilon, N, bob_sheet, eve_sheet); 
    sum_secrecy_rate(i) = sum(sum(user_grid)) / N;
    % 
end
%
P = 10 * log10(P/(.97865^2));
plot(P, sum_secrecy_rate, '-*')
    