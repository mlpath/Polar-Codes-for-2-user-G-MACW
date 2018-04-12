% this is the control structure for generating all channels up to level N
%
clc
clear all
% code block length (depth of binary tree is log2(N))
N = 256;
% fidelity parameter (number of levels in quantization/approximation)
fidelity_par = 30;
% size of block processed by arikan_degrading (affects memory)
SZ_buffer = 4e7;
% we have used a less efficient older version of degrading_merge to convert
% our continous MAC BIAWGN channel to a discrete one. these are the
% characteristics:
% spacing between samples in time = 1e-4,
% noise variance = .97865^2,
% number of users = 2,
% alphabet = binary,
%
% only run it once so that we make sure our record of raw_channel is sound
% degrading_merge(fidelity_par);
%
% initialize the binary level (starting from level 1, 2^1 nodes)
ix_old = dec2bin(0:1,1);
%
% this counter should count up to log2(N)
count = 0;
% let it continue until n = 10, then we'll check it with mutual_inf
ix_old = dec2bin(0:(2^count - 1), count );
while (count <= 7)
    for i = 1:size(ix_old,1)
    %%%%%%%%%%%%%%%%%%%%%%%%% Main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if count
            % the channel we are working on
            chan_str = ix_old(i,:);
        else
            % empty string to create chan_0 and chan_1
            chan_str = '';
        end
        % input to this function is a channel generated in previous level
        arikan_degrading_v2(chan_str, fidelity_par, SZ_buffer, t);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    count = count + 1;
    ix_old = dec2bin(0:(2^count - 1), count );
end
