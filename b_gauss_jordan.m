% this function implements gauss jordan over binary matrices
% input is channel matrix in raw form
% output is rank of channel and 
function [rref, rk] = b_gauss_jordan(raw_m)
%
% clc ;
% clear all;
% examples
%
% raw_m = mod(randi(2,4,8),2);
% 
% raw_m = [ 0     0     1     0     1     1     0     0;
%           1     1     0     1     0     1     1     1;
%           1     1     1     0     1     0     1     1;
%           1     1     0     0     0     1     1     0];
% raw_m = [1 0 1 0 1 0 0;
%          0 1 1 1 0 0 1;
%          0 0 1 1 1 1 0;
%          1 1 0 0 1 0 0];
% raw_m = [0 0 0 0 1 0 1 1 ;
%          1 1 0 1 1 0 0 0 ; 
%          1 0 0 0 0 0 0 0 ;
%          1 1 0 1 1 0 0 0 ];
         

q = size(raw_m,1);
%
m = size(raw_m,2);
%
count = 0; 
%
clean = 1;
i = 1; j = 1;
while ((i <= q) && (j <= m))
    % cleaning up each column..
    if (raw_m(i,j) ~= 1)
        ix = find(raw_m(:,j) == 1);
        % this ensures that we only disturb later rows of the matrix
        tmp = ix>i;
        ix = ix(tmp);
        if isempty(ix)
            clean = 0;
            j = j+1;
        else
            tmp_m = raw_m;
            % swap the row with zero leading element with raw_m nonezero one
            tmp_m(i,:) = raw_m(ix(1),:) ;
            tmp_m(ix(1),:) = raw_m(i,:);
            raw_m = tmp_m;
        end
    end
    %
    if clean
        ix = find(raw_m(:,j) == 1);
        % this time to clean later rows
        tmp = ix~=i;
        ix = ix(tmp);
        for l=1:length(ix)
            raw_m(ix(l),:) = xor(raw_m(i,:),raw_m(ix(l),:));
        end
        %
        i = i + 1;
        j = j + 1;
    end
    clean = 1;
end
% 
% binary matrix in row-reduced echelon form
rref = raw_m;
%
for i = 1:q
    if (raw_m(i,:) == zeros(1,m));
        count = count + 1;
    end
end
rk = q - count;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
end

    