% Q1 of Data Incubators challenge
%
clc
clear all
% number of input seats
N = 5e4;
n = 1:N;
% define a flag so that we know all seats are taken
occupied = 0;
% counter for random selections
l = 1;
% preallocate it
x = zeros(1,N);
% number of samples
nsamp = 100000;
% define a random variable f for fraction of occupied seats
f = zeros(1,nsamp);
% generate nsamp samples of this random variable
for i = 1:nsamp
    i
    while ~occupied
        l;
        % new travellor selects a seat
        s = randi(length(n),1);
        x(l) = n(s);
        n = n( n ~= x(l));
        % rule out seats next to s from the selection pool
        ix1 = n == (x(l)+1);
        if sum(ix1)
            n = n(~ix1);
        end
        %
        ix2 = n == (x(l)-1);
        if sum(ix2)
            n = n(~ix2);
        end
        %
        if isempty(n)
            occupied = 1;
        end
        % update the counter for next experiment
        l = l + 1;
    end
    %
    tic
    x(l :end) = 0; %[];
%     x;
%     % want to see a pattern of taken seats
%     seats = zeros(1,N);
%     seats(x) = 1;
%     % our desired random variable is fraction of occupied seats
%     f(i) = sum(seats)/N;
    f(i) = sum(x ~= 0)/N;
    % reset parameters of our random seat generator
    occupied = 0;
    l = 1;
%     x = zeros(1,N);
    n = 1:N;
end
f_5e4 = f;
save('fraction.mat','f_5e4','-v6')
% mean
m = sum(f)/nsamp;
% standard deviation
std = sqrt(sum((f - m).^2)/(nsamp-1));
%
fprintf('\n mean is %f and standard variation equal to %f \n', m, std);