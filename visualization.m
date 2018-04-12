% visualizing the mutual information vectors 
%
clc
clear all
close all
% which element?
i = 2;
% block length 
n = 256;
% excel identifier
excel_id = '20';
% the name of excel file where mutual information vectors are stored
filename = strcat('secrecy_graph_',excel_id,'.xlsx');
img_name = strcat('visu_', excel_id,'_',num2str(i),'.png');
% image size?
sz = 3*n/8;
% range of table 
rng = strcat('A1:D', num2str(n));
% Bob
R1 = xlsread(filename,1,rng); 
% Eve
R2 = xlsread(filename,2,rng); 
% 
% figure;
% plot(R2(:,4), '.'); 
% hold on; 
% plot(R2(:,2), 'or', 'MarkerSize', 10); 
% hold on; 
% plot(R2(:,3), 'ok', 'MarkerSize', 2);
% xlabel('index of extremal channels'); 
% ylabel('mutual information (b.p.c.u)');
%
figure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r1 = repmat(R1(:,i)', sz, 1);
r2 = repmat(R2(:,i)', sz, 1); 
%
I = [r1; r2];
% zero padding for border
I = [I;zeros(1,size(I,2))];
if i == 4
    imshow(I, [0,2], 'Border', 'tight');
else
    imshow(I, 'Border', 'tight');
end
imwrite(I, img_name, 'Border', 'tight')