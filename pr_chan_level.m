% this will compute channel level likelihoods
%
function channel_probs = pr_chan_level(in_seq, chan_filename, disc_cont)
%
chan_filename = lower( chan_filename );
if strcmpi(disc_cont, 'discrete')
    % y is the vector of received signal samples
    % chan_type is the name of conditional probability table we are going to
    % use, i.e. raw_channel_10.mat or raw_channel_14.mat
    %
    % compute a 4 x N matrix of conditional probabilities at channel level
    %
    load( chan_filename );
    % this is the conditional probabilities for 256 symbols we got at the
    % output of our underlying discrete channel
    channel_probs = 70 * qt_probs(:, in_seq);
    % channel_probs(channel_probs == 0) = 1e-300;
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % channel spec for continuous case
    powers = chan_filename.P;
    gains = chan_filename.G;
    noise_var = chan_filename.N;
    % number of users
    t = max(size(gains));
    % number of samples
    N = size(in_seq,2);
    % preallocate channel_probs for continuous case
    channel_probs = zeros(2^t,N);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         conditional probabilities for continuous awgn mac
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u = dec2bin(0:2^t - 1);
    U = ones(size(u));
    U( u == '1' ) = -1;
    % the result of this should be a 2 by 1 vector
    gains = gains .* sqrt(powers);
    % and this is giving the value of means
    size(U);
    size(gains);
    m = sum(bsxfun(@times,U',gains),1);
    %
    size(m);
    for i = 1:2^t
        % loop over possible combinations (i.e. means)
        channel_probs(i,:) = exp(-(in_seq - m(i)).^2/(2*noise_var))/(sqrt(2*pi*noise_var));
    end
    %
end
end