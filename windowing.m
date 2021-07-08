function windowedData = windowing(timeWindow, windowShift, ...
    inputData)
% Function that takes as input a time window, a time shift (opposite of
% overlap) and an array of data with times and values, and returns a cell
% array of data windowed using the window length and the amount of shift 
% between each moving window, where the first cell of each cell array's row
% will be the start time of the window, the second cell in each row will be
% the end time of the window, and the third cell will contain a row vector
% of data corresponding to this window

% ----------- Inputs ----------------- %
% timeWindow: # of seconds in each moving window
% windowShift: # of seconds to shift each time before creating next window
% inputData: Nx2 array of time and data, where the first column stores time
% and the second column stores the corresponding data

% ----------- Outputs ----------------- %
% windowedData: Mx3 cell array of window starts/ends and 
% corresponding data, where first two columns store window start and end,
% respectively, and the third cell column contains row vectors (if sampling
% rate is fixed, these row vectors would be the same length; otherwise, not
% necessarily the case) of corresponding data

% Store first and last time points in the data
firstTime = inputData(1, 1);
lastTime = inputData(end, 1);

% Calculate time available
timeAvailable = lastTime - firstTime;

% Compute number of windows possible
% We take the ceiling here because we force the last window to contain the
% last datapoint, even if that disobeys the windowShift amount
numWindows = ceil((timeAvailable - timeWindow)/windowShift) + 1;

% Initialize output
windowedData = cell(numWindows, 3);

% Store the data as windows, row by row
for i = 1:(numWindows - 1)
    % Store start time and end time
    startTime = firstTime + (i-1)*windowShift;
    endTime = firstTime + (i-1)*windowShift + timeWindow;
    windowedData(i, 1:2) = {startTime, endTime};
    
    % Store corresponding data
    windowedData{i, 3} = ...
        inputData(inputData(:, 1) <= endTime ...
        & inputData(:, 1) > startTime, 2)';
end

% Store the last window
windowedData(numWindows, 1:2) = {lastTime - timeWindow, lastTime};
windowedData{numWindows, 3} = inputData(inputData(:, 1) <= lastTime & ...
    inputData(:, 1) >= lastTime - timeWindow, 2)';

end

