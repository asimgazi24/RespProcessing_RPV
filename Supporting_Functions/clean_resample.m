function [outputData, outputTime] = clean_resample(OG_data, OG_time, fs, method)

% The point of this function is to resample, but to remove the endpoint
% effects

% outputData: resultant "clean" resampled data in vector format (i.e., Mx1
% or 1xM)
% outputTime: corresponding time vector in vector format

% OG_data: original data as vector (if row, we return row; if column, we
% return column)
% OG_time: original time (nonuniform) as vector (same note as above)

% fs: desired sampling frequency

% method: desired interpolation method

data_row = false;
time_row = false;

% If in row vector form, set data_row to true and then transpose for the
% rest (we will transpose back at the end)
if size(OG_data, 2) > size(OG_data, 1)
    OG_data = OG_data';
    data_row = true;
end

% If in row vector form, transpose and make note that it was in row vector
% form
if size(OG_time, 2) > size(OG_time, 1)
    OG_time = OG_time';
    time_row = true;
end

% compute slope and offset
a(1) = (OG_data(end) - OG_data(1)) / (OG_time(end) - OG_time(1));
a(2) = OG_data(1);

% detrend the signal (OG_time may not start at zero, so that needs to be
% subtracted from the linear equation evaluation)
xdetrend = OG_data - polyval(a, OG_time - OG_time(1));

% resample using desired method at desired sampling frequency
[ydetrend, outputTime] = resample(xdetrend, OG_time, fs, method);

% add back the trend
outputData = ydetrend + polyval(a, outputTime - outputTime(1));

% If fed in as row vector, return as row vector
if data_row
    outputData = outputData';
end

% If fed in as row vector, return as row vector
if time_row
    outputTime = outputTime';
end

end
