% a continuous MAC channel 
function out_seq = c_mac_awgn(in_seq, chan_spec)
% in_seq is input  0/1 symbols 
% power, gains, and noise_var are respectively transmitter power, channel
% gains, and noise variance for this receiver
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc
% clear all
% close all
% user 1 -> first row 
% user 2 -> second row 
% 
% in_seq = mod(randi(2,2,20),2);
% power = [1 1]; 
% gains = [1 1];
% noise_var = .97865^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% channel specification
powers = chan_spec.P;
gains = chan_spec.G; 
noise_var = chan_spec.N; 
%
% map zeros to 1 and ones to -1
in_seq(in_seq == 1) = -1; 
in_seq(in_seq == 0) = 1;
% number of users
t = max(size(gains));
% number of time samples 
N = size(in_seq,2);
%
if size(gains, 1) < t
    % change the dimensions
    gains = gains'; 
    powers = powers';
end
% now scale it with power and channel gains and sum it over columns to
% calculate the mean of gaussian MAC
m = sum(bsxfun(@times, in_seq, sqrt(powers).*gains),1);
% 
% now add awgn noise 
out_seq = m + sqrt( noise_var ) * randn(1,N);
%
% [in_seq; m; out_seq];
