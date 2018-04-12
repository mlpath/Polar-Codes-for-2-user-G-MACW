% in this code I will calculate the numerical values of I(Xs;Y|Xsc) for a
% fixed power level and channel coefficient. 
%               /
% I(Xs;Y|Xsc) = |  sigma(
%               /
% 
function I = mutual_bpsk(P1,P2,h1,h2,sigma1)
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
% conditional probability of channel output given the symbols
p_cond = @(y,x1,x2)(exp(-((y - sqrt(P1)*h1*x1 - sqrt(P2)*h2*x2).^2)./(2*(sigma1^2)))/sqrt(2*pi*(sigma1^2)));
%
p_marg1 = @(y,x2)((p_cond(y,1,x2) + p_cond(y,-1,x2))/2);
p_marg2 = @(y,x1)((p_cond(y,x1,1) + p_cond(y,x1,-1))/2);
% 
p_marg12 = @(y)((p_cond(y,1,1)+p_cond(y,1,-1)+p_cond(y,-1,1)+p_cond(y,-1,-1))/4);
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
x = [1 1;1 -1;-1 1;-1 -1];
% range of integration
y = -10*sigma1-sqrt(P1)*h1-sqrt(P2)*h2:.005:10*sigma1+sqrt(P1)*h1+sqrt(P2)*h2;
%
itmp_1 = 0; 
itmp_2 = 0;
itmp_12 = 0; 
% 
itmp_10 = 0; 
itmp_20 = 0;
%
pp_itmp_10 = 0; 
%
% plot(y, p_marg12(y))
s = 1e-40;
for i = 1:4
    x1 = x(i,1); 
    x2 = x(i,2); 
    %  
    tmp_1 = @(y)((p_cond(y,x1,x2) .* log2((p_cond(y,x1,x2)+s)./p_marg1(y,x2)))/4);
    itmp_1 = itmp_1 + quadgk(tmp_1,y(1),y(end));
    %
    tmp_2 = @(y)((p_cond(y,x1,x2) .* log2((p_cond(y,x1,x2)+s)./p_marg2(y,x1)))/4);
    itmp_2 = itmp_2 + quadgk(tmp_2,y(1),y(end)); 
    %
    tmp_12 = @(y)((p_cond(y,x1,x2) .* log2((p_cond(y,x1,x2)+s)./p_marg12(y)))/4);
    itmp_12 = itmp_12 + quadgk(tmp_12,y(1),y(end)); 
    % 
end
x = [-1 1];
% point to point conditional probabilities 
%
pp_cond = @(y,x1)(exp( - (y -h1*x1*sqrt(P1)).^2/(2*sigma1^2))/sqrt(2*pi*sigma1^2));
pp_marg1 = @(y)((pp_cond(y,1) + pp_cond(y,-1))/2);
% 
pp_y = -10*sigma1 - h1*sqrt(P1):.005:10*sigma1 + h1*sqrt(P1); 
for i = 1:2
    x1 = x(i);
    x2 = x1; 
    %
    tmp_10 = @(y)(p_marg2(y,x1) .* log2(p_marg2(y,x1)./p_marg12(y))/2);
    itmp_10 = itmp_10 + quadgk(tmp_10,y(1),y(end));
    % 
    tmp_20 = @(y)(p_marg1(y,x2) .* log2(p_marg1(y,x2)./p_marg12(y))/2);
    itmp_20 = itmp_20 + quadgk(tmp_20,y(1),y(end)); 
    % 
    pp_tmp_10 = @(y)(pp_cond(y,x1).*log2((pp_cond(y,x1)./pp_marg1(y)))/2);
    pp_itmp_10 = pp_itmp_10 + quadgk(pp_tmp_10,pp_y(1),pp_y(end)); 
end
I = [itmp_1,itmp_2,itmp_12,itmp_10,itmp_20];
end
% 
% pp_I1 = pp_itmp_10;

%
% debug
%
% itmp_10 + itmp_2
% itmp_20 + itmp_1
% itmp_12
% end

    