function rsp_onsetTimes = rspOnset_find(IBIs,rspSignal)
% This function takes interbreath intervals and a respiration signal and
% finds the respiration onset / finish point in between the two peaks that
% characterize the IBI

% ---------- Inputs -------------- %
% IBIs: Nx2 array, where the first column stores time of the later peak of
% each inter-breath interval, and the second column stores the length of
% time of the inter-breath interval
% rspSignal: Mx2 array, where the first column stores time and second
% column stores respiraiton signal values

% ---------- Outputs ------------- %
% rsp_onsetTimes: Nx1 array of onset times, one per inter-breath interval


% Initialize result
rsp_onsetTimes = zeros(size(IBIs, 1), 1);

% We need to find the onset for every IBI
for i = 1:size(IBIs, 1)
    % Store the beginning time and ending time of this IBI
    startTime = IBIs(i, 1) - IBIs(i, 2);
    endTime = IBIs(i, 1);
    
    % We isolate the signal values corresponding to this IBI
    sig_window = rspSignal(rspSignal(:, 1) >= startTime...
        & rspSignal(:, 1) <= endTime,...
        :);
    
    % We need the row index of the minimum value within this window
    [~, rowIndx] = min(sig_window(:, 2));
    
    rsp_onsetTimes(i) = sig_window(rowIndx, 1);
end



end

