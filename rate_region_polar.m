% want to calculate the tdma and superpposition region for this parameters
% in gaussian channel with gaussian symbols 
% 
clc
clear all
close all
% 
h1 = 1.11; h2 = 1; g1 = 0.9; g2 = 0.75; 
sigma1 = 1; sigma2 = 1; 
% h1 = 1; h2 = 1; g1 = 1; g2 = 0.4; 
% sigma1 = 0.97865; sigma2 = sigma1;
P1 = 1; 
P2 = 1; 
% 
% 
Ib = mutual_bpsk(P1,P2,h1,h2,sigma1);
Ie = mutual_bpsk(P1,P2,g1,g2,sigma2); 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
R1_sup = max(Ib(1) - Ie(4),0); 
R2_sup = max(Ib(2) - Ie(5),0); 
R12_sup = max(Ib(3) - Ie(3),0); 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for polar
% 
[Ib(1),Ib(2),Ie(1),Ie(2)]
[Ib(5),Ie(1),Ie(2)]
Ie(2) - Ie(5)
if Ib(5) < Ie(2)
    fprintf('figure 1-11\n');
    R1_pol = Ib(1) - Ie(4) - (Ie(2) - Ib(5));
    R2_pol = Ib(2) - Ie(5) - (Ie(2) - Ib(5)); 
    % 
    R12_pol = Ib(3) - Ie(3) + (Ie(1) - Ie(5)) - (Ie(2) - Ib(5));
else
    fprintf('figure 1-10 or 1-9 \n');
    R1_pol = Ib(1) - Ie(4);
    R2_pol = Ib(2) - Ie(5); 
    % 
    R12_pol = Ib(3) - Ie(3) + (Ie(2) - Ie(5)) + max(Ie(1) - Ib(5), 0);
end
% 
% 
% now how to visualize it? 
r1 = linspace(0,0.9,150); 
r2 = r1;
% 
n = length(r1); 
for i = 1:n
    for j = 1:n
    if (r1(i) <= R1_pol) && (r2(j) <= R2_pol) && (r1(i) + r2(j) <= R12_pol)
        plot(r1(i), r2(j), '+k','MarkerSize', 4.5)
        hold on; 
    end
    if (r1(i) <= R1_sup) && (r2(j) <= R2_sup) && (r1(i) + r2(j) <= R12_sup)
        plot(r1(i), r2(j), 'ob','MarkerSize', 4)
        hold on; 
    end 
    end
end
figure(1)
xlabel('R_1');
ylabel('R_2');

% plot(R1_tdma, R2_tdma, '.k'); 
% 
%






