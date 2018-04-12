% this function finds the message set we have to check in ml decoding of
% MAC system. On each given time slot, we find the frozen data and
% calculate all possibilities of 0s and 1s given the fixed data
% 
% the output is a matrix having as many rows as the number of possible
% message sets (in string format)
% 
% if a matrix is empty, that means we don't need to find 
%
function Uvec = ml_message_set(frozen_map, frozen_data, time_slot)
% pattern of transmitters at time index count
idx = find(frozen_map(:,time_slot) == 1);
if ~isempty(idx) 
    % all possible combinations of transmitt bits at this rime slot
    tx_bits = dec2bin(0:sum(idx));
    Uvec = repmat(dec2bin(0,t), size(tx_bits,1), 1);
    Uvec(:,idx) = tx_bits; 
    idx_c = setdiff(1:t, idx); 
    Uvec(:,idx_c) = num2str(frozen_data(idx_c,time_slot));
else
    % means that all symbols are frozen 
end
end 