% in this code I will calculate the numerical values of I(Xs;Y|Xsc) for a
% fixed power level and channel coefficient. 
%               /
% I(Xs;Y|Xsc) = |  sigma(
%               /
% 
function I = mutual_bpsk_norm(Px1, Px2, h1, h2)
% power level
% clc; clear all; close all
% P1 = 10; 
% P2 = 10; 
% % channel coefficient
% h1 = 1;
% h2 = 1;
% % 
% sigma1 = 0.97865;
% 
% user indicator in this problem
% S = [1, 2]; 
% L = length(S);
% a priori probability of the symbols
% p_prior = 1/2^length(S); 
%
x = [1 1;1 -1;-1 1;-1 -1];
x = x * diag([sqrt(Px1), sqrt(Px2)]);
% range of integration
y = -10*1+min(x(:,1))*h1+min(x(:,2))*h2:.01:max(x(:,1))*h1+max(x(:,2))*h2+10*1;
%
% conditional probability of channel output given the symbols
p_cond = @(y,x1,x2)(exp(-((y - h1*x1 - h2*x2).^2)/2)/sqrt(2*pi));
%
p_marg1 = @(y,x2)((p_cond(y,x(1,1),x2) + p_cond(y,-x(1,1),x2))/2);
p_marg2 = @(y,x1)((p_cond(y,x1,x(1,2)) + p_cond(y,x1,-x(1,2)))/2);
% 
p_marg12 = @(y)((p_cond(y,x(1,1),x(1,1))+p_cond(y,x(1,1),-x(1,1))+p_cond(y,-x(1,1),x(1,1))+p_cond(y,-x(1,1),-x(1,1)))/4);
%
% plot(y,p_cond(y,1,1),'k');
% hold on; 
% % 
% lg  = zeros(2^L - 1, L);
% str = dec2bin(1:2^L - 1, L);
% for i = 1:L
%     lg(:,i) = logical(str2num(str(:,i)));
% end
% % 
% for i = 1:2^L-1
%     % choose the appropriate term for computation
%     Lg = loical(lg(i,:));
%     % 
% end
% 
% 
itmp_1 = 0; 
itmp_2 = 0;
itmp_12 = 0; 
% 
itmp_10 = 0; 
itmp_20 = 0;
%
for i = 1:4
    x1 = x(i,1); 
    x2 = x(i,2); 
    %  
    tmp_1 = @(y)((p_cond(y,x1,x2) .* log2(p_cond(y,x1,x2)./p_marg1(y,x2)))/4);
    itmp_1 = itmp_1 + quadgk(tmp_1,y(1),y(end));
    %
    tmp_2 = @(y)((p_cond(y,x1,x2) .* log2(p_cond(y,x1,x2)./p_marg2(y,x1)))/4);
    itmp_2 = itmp_2 + quadgk(tmp_2,y(1),y(end)); 
    %
    tmp_12 = @(y)((p_cond(y,x1,x2) .* log2(p_cond(y,x1,x2)./p_marg12(y)))/4);
    itmp_12 = itmp_12 + quadgk(tmp_12,y(1),y(end)); 
    % 
end
%
for i = 2:3
    x1 = x(i,1);
    x2 = x(i,2); 
    %
    tmp_10 = @(y)(p_marg2(y,x1) .* log2(p_marg2(y,x1)./p_marg12(y))/2);
    itmp_10 = itmp_10 + quadgk(tmp_10,y(1),y(end));
    % 
    tmp_20 = @(y)(p_marg1(y,x2) .* log2(p_marg1(y,x2)./p_marg12(y))/2);
    itmp_20 = itmp_20 + quadgk(tmp_20,y(1),y(end)); 
end
I = [itmp_1,itmp_2,itmp_12,itmp_10,itmp_20];


%
% debug
%
% itmp_10 + itmp_2
% itmp_20 + itmp_1
% itmp_12
% end

    