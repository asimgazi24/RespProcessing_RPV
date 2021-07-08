function plausible_IBI = IBI_findplaus(rsp_peakTimes, max_IBI, min_IBI)
% This function takes as input peak times for a respiration signal and
% outputs plausible inter-breath intervals, where plausible here refers to
% physiologically plausible, i.e., if IBI > 10 s or IBI < 1.4 s, then we
% remove that inter-breath interval

% ------------ Inputs -------------- %
% rsp_peakTimes: Nx1 array of peak respiration times

% ----------- Outputs --------------- %
% plausible_IBI: ((N-1) - K)x2 array of physiologically plausible
% inter-breath intervals, where K represents the number of removed IBIs,
% the first column will store times of the second peak of each interval,
% and the second column will store the interval time

% Initialize IBI array
plausible_IBI = zeros((length(rsp_peakTimes)-1), 2);

% Store times
plausible_IBI(:, 1) = rsp_peakTimes(2:end);

% Store intervals
plausible_IBI(:, 2) = diff(rsp_peakTimes);

% Remove all implausible intervals
% max_IBI = 10;
% min_IBI = 1.4;
for i = size(plausible_IBI, 1):-1:1
    if plausible_IBI(i, 2) > max_IBI || plausible_IBI(i, 2) < min_IBI
        plausible_IBI(i, :) = [];
    end
end

end

