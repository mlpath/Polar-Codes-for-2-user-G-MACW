% want to calculate the tdma and superpposition region for this parameters
% in gaussian channel with gaussian symbols 
% 
clc
clear all
close all
% 
h1 = 0.6; h2 = 0.6; g1 = .5; g2 = 0.65; 
sigma12 = 0.96; sigma22 = 0.96; 
P1 = 2; 
P2 = 2; 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% 
R1_sup = max(0.5*log2(1 + (P1*h1^2)/(sigma12)) - 0.5*log2(1 + (P1*g1^2)/(sigma22+P2*g2^2)),0); 
R2_sup = max(0.5*log2(1 + (P2*h2^2)/(sigma12)) - 0.5*log2(1 + (P2*g2^2)/(sigma22+P1*g1^2)),0); 
R12_sup = max(0.5*log2(1 + (P1*h1^2 + P2*h2^2)/(sigma12)) - 0.5*log2(1 + (P1*g1^2 + P2*g2^2)/(sigma22)),0); 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for all alpha1 , alpha2 positive and ess than 1, 
a = 0.005:.01:1;
R1_tdma = zeros(size(a)); 
R2_tdma = R1_tdma; 
%
for i = 1:length(a)
    R1_tdma(i) = max(0.5* a(i) * (log2(1 + (P1*h1^2)/(sigma12*a(i))) - log2(1 + (P1*g1^2)/(sigma22 * a(i)))),0); 
    R2_tdma(i) = max(0.5* (1 - a(i)) * (log2(1 + (P2*h2^2)/(sigma12*(1 - a(i)))) - log2(1 + (P2*g2^2)/(sigma22 * (1-a(i))))),0); 
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% now how to visualize it? 
r1 = linspace(0,1,500); 
r2 = r1;
% 
n = length(r1); 
for i = 1:n
    for j = 1:n
    if (r1(i) <= R1_sup) && (r2(j) <= R2_sup) && (r1(i) + r2(j) <= R12_sup)
        plot(r1(i), r2(j), '.b')
        hold on; 
    end
    end
end
hold on; 
plot(R1_tdma, R2_tdma, '.k'); 
% 
% J = 0:0.005:10;
% p = zeros(size(J));
% for i = 1:length(J)
% p(i) = 0.5*(log2(1 + (P1*h1^2)/(sigma12 + J(i)*h2^2)) - log2(1 + (P1*g1^2)/(sigma22 + J(i)*g2^2)));
% end
% [val, I] = max(p)





