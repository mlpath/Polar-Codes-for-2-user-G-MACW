% load R1 and R2 rate matrices
clc
clear all
load cj_sig_code.mat
% 
M = 1; 
J = 1;
%
R1_jam = zeros(size(R1, 1),2^M - 1);
R2_jam = R1_jam;
%
for i = 1:size(R1, 1)
    R1_jam(i,:) = cj_ratevector(R1(i,:),1, 1);
    R2_jam(i,:) = cj_ratevector(R2(i,:),1, 1);
end
%
