function [windowed_AR_1, numWindRmvd] = ...
    AR_1Breath(data_time_array, windows, mean_IBI, scaleBool)
% This function takes as input an array that stores time and data, as well
% as another array that stores window starts and window ends. It returns a
% three column array where first two columns are window starts and ends,
% and last column contains the AR(1) computed for those windows

% ----------- Inputs --------------- %
% data_time_array: Nx2 array, where first column stores time and second
% column stores data
% windows: Mx2 array, where first column stores window starts and second
% column stores window end times
% mean_IBI: mean IBI - used to figure out how much data should be expected
% within a window
% scaleBool: normalize to make zero-lag = 1, and rest corresponding?

% ----------- Outputs ---------------- %
% windowed_AR_1: (M-K)x3 array, where first column stores window starts,
% second column stores window end times, and third column stores the
% corresponding autocorrelations at 1 breath lag.
% K here represents the windows that were removed due to a lack of data
% numWindRmvd: stores K

% Initialize result
windowed_AR_1 = [];

% Set the minimum number of datapoints required to output a value for a
% particular window
percent_rqrd = 50;          % I set 50% as my threshold
% Using the median inter-breath interval, we can figure out approx. how many
% datapoints can be expected within a window and multiply that by the
% percent_rqrd to get the minimum number of datapoints
windowTime = windows(1, 2) - windows(1, 1);     % assuming equal-sized
numData_expected = windowTime/mean_IBI;
numData_rqrd = (percent_rqrd/100)*numData_expected;

% Separate out data and time to make code below cleaner
time = data_time_array(:, 1);
data = data_time_array(:, 2);

% Counter to store number of windows with insufficient data
numWindRmvd = 0;

% Iterate through all windows
for i = 1:size(windows, 1)
    % Store this specific window's start and end times
    start_time = windows(i, 1);
    end_time = windows(i, 2);
    
   % Parse data to consider only the window specified
    data_window = data(time >= start_time & time <= end_time);
    
    % Only calculate if we have enough elements in our window
    if numel(data_window) >= numData_rqrd
        % Compute autocorrelation at one lag
        % Here we use the xcorr() function to return the autocorrelation
        % sequence and then just single out when lag = 1
        if scaleBool
            [AR_seq, lags] = xcorr(data_window, 'normalized');
        else
            [AR_seq, lags] = xcorr(data_window);
        end
        
        
        % Single out lag = 1
        AR_1 = AR_seq(lags > 0.5 & lags < 1.5);
        
        % First column will be start time, second column will be end time,
        % and third column will be CV
        windowed_AR_1 = [windowed_AR_1; start_time, end_time, AR_1];
    else
        % This window will be ignored, so increment counter
        numWindRmvd = numWindRmvd + 1;
    end
end

end

