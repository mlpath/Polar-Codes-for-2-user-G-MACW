% I need to answer the question that if a gaussian MAC with unequal channel
% gains and power levels is symmetric or not
%               /
% I(Xs;Y|Xsc) = |  sigma(
%               /
% 
function [T_1,T_2, is_t_symmetric] = T_symmetric(P1,P2,h1,h2,sigma1)
clc
% clear all
% power level
% clc; clear all; close all
% P1 = 10; 
% P2 = 10; 
% % channel coefficient
% h1 = 1;
% h2 = 1;
% % % 
% sigma1 = 0.97865;
% 
% user indicator in this problem
% S = [1, 2]; 
% L = length(S);
% a priori probability of the symbols
% p_prior = 1/2^length(S); 
% 
is_t_symmetric = 'false'; 
% conditional probability of channel output given the symbols
p_cond = @(y,x1,x2)(exp(-((y - sqrt(P1)*h1*x1 - sqrt(P2)*h2*x2).^2)./(2*(sigma1^2)))/sqrt(2*pi*(sigma1^2)));
%
p_marg12 = @(y)((p_cond(y,1,1)+p_cond(y,1,-1)+p_cond(y,-1,1)+p_cond(y,-1,-1)));
%
% 
t_1 = @(y,x1)((p_cond(y,x1,1).*log2(p_cond(y,x1,1)./p_marg12(y))) + (p_cond(y,x1,-1).*log2(p_cond(y,x1,-1)./p_marg12(y)))); 
t_2 = @(y,x2)((p_cond(y,1,x2).*log2(p_cond(y,1,x2)./p_marg12(y))) + (p_cond(y,-1,x2).*log2(p_cond(y,-1,x2)./p_marg12(y)))); 
% 
% range of integration 
y = -10*sigma1-sqrt(P1)*h1-sqrt(P2)*h2:.0005:10*sigma1+sqrt(P1)*h1+sqrt(P2)*h2;
%
x = [1 -1];
T_1 = zeros(1,2); 
T_2 = T_1; 
% 
for i = 1:2 
T_1(i) = quadgk(@(y)(t_1(y,x(i))),y(1),y(end)); 
T_2(i) = quadgk(@(y)(t_2(y,x(i))),y(1),y(end)); 
end
T_1;
T_2;
T_1 - T_2;
if sum(abs(T_1 - T_2),2) < 1e-5
    is_t_symmetric = 'true'; 
end 
% 
fprintf('t-symmetric is %s \n',is_t_symmetric);
% T_1, T_2
% is_t_symmetric
end 