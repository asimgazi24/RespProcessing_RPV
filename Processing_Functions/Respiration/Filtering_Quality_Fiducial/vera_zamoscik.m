function [peakTimes, minTimeBetween] = vera_zamoscik(resp_window, Fs)
% This function applies peak and onset detection according to the methods
% detailed in V Zamoscik et al. (2018) in Psychological Medicine

% -------- Inputs --------- %
% resp_window: 1xN row vector of respiration data that has already been
% bandpass filtered
% Fs: sampling frequency of respiration data (important to convert from
% indices to times) (Hz)

% -------- Outputs ---------- %
% peakTimes: Mx1 vector of peak times detected
% minTimeBetween: scalar value storing the minimum time between parameter
% used when finding peaks


% Step 1: Normalize the window (row vector, so dimension 2) to mean 0
% stdev 1000 (typo in paper, correction based on code shared online)
sigma = 1000;
resp_window = sigma*normalize(resp_window, 2);

% Step 2: Find all peaks with a minimum time between of 1.4 seconds and a
% minimum prominence of 1/3(sigma)
minTimeBetween = 1.4;       % seconds since Fs is in Hz
%minProminence = (1/3)*sigma;
minProminence = (1/2)*sigma;        % Found that 1/2 stdev worked better
[~, peakTimes] = findpeaks(resp_window, Fs, ...
    'MinPeakDistance', minTimeBetween, 'MinPeakProminence', minProminence);


end

