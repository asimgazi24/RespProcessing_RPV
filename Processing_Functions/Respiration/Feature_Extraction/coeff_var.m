function [windowed_CV, numWindRmvd] = ...
    coeff_var(data_time_array, windows, mean_IBI)
% This function takes as input an array that stores time and data, as well
% as another array that stores window starts and window ends. It returns a
% three column array where first two columns are window starts and ends,
% and last column contains the coefficients of variation computed for those
% windows

% ----------- Inputs --------------- %
% data_time_array: Nx2 array, where first column stores time and second
% column stores data
% windows: Mx2 array, where first column stores window starts and second
% column stores window end times
% mean_IBI: mean IBI - used to figure out how much data should be expected
% within a window

% ----------- Outputs ---------------- %
% windowed_CV: (M-K)x3 array, where first column stores window starts,
% second column stores window end times, and third column stores the
% corresponding coefficients of variation. K here represents the windows
% that were removed due to a lack of data
% numWindRmvd: stores K

% Initialize result
windowed_CV = [];

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
        % Compute coefficient of variation (CV), which is = mean/stdev
        CV = mean(data_window)/std(data_window);
        
        % First column will be start time, second column will be end time,
        % and third column will be CV
        windowed_CV = [windowed_CV; start_time, end_time, CV];
    else
        % This window will be ignored, so increment counter
        numWindRmvd = numWindRmvd + 1;
    end
end

end

