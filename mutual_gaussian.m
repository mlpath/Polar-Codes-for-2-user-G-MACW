% test mutual information calculation with continuous guassian symbols 
% 
clc
clear all
%
% var_x = 1; 
% x1 = -5:.01:5;
% x2 = x1; 
% 
h1 = 1; 
h2 = 1;
% 
sigma = 1;
% y = -5*sigma + min(x1)*h1 + min(x2)*h2:.01:5*sigma + max(x1)*h1 + max(x2)*h2;
% conditional probability
p_cond = @(y,x1,x2)(exp(-((y - h1*x1 - h2*x2).^2)/(2*sigma^2))/(sqrt(2*pi*sigma^2)));
% input symbols (ind.) joint probability
p_symb = @(x1,x2)(exp(-(x1.^2 + x2.^2)/2)/(2*pi));
% 
y = -10:.01:10;
tmp = @(x1,x2)(p_cond(y,x1,x2).* p_symb(x1,x2));
itmp_12 = dblquad(tmp,-5,5,-5,5);% x1(1), x1(end), x2(1), x2(end));
p_marg12 = @(y)(itmp_12);
% 
% tmp = @(y,x1,x2)(p_cond(y,x1,x2).* log2(p_cond(y,x1,x2)./p_marg12));
% I = triplequad(tmp, x1(1),x1(end),x2(1),x2(end),y(1),y(end))
% I_theory = log2(1 + h1^2 + h2^2)/2

