function [respOnsetTimes, respPeakTimes] = ...
    breathDetect_RSP(resp_windowed, resp_signal, Fs, varargin)
% This function will take a respiration signal divided into windows and
% return the corresponding breath onset and peak times

% -------- Inputs --------- %
% resp_windowed: Mx3 cell array, where the first column of cells contain
% window start times, the second column of cells contain window end times,
% and the final column of cells contain row vectors of respiration data
% resp_signal: Kx2 filtered and resampled respiration signal in its entirety
% (not windowed), where first column stores time, second column stores
% signal value
% Fs: sampling frequency of respiration data (important to convert from
% indices of peaks to times)
% (Optional)
% method: 'vera_zamoscik' or 'count_adv' - which method would you like to
% use to determine peaks and onsets? (Default: 'count_adv'

% -------- Outputs ---------- %
% respOnsetTimes: Nx1 array of respiration onset times
% respPeakTimes: (N-1)x1 or Nx1 or (N+1)x1 array of respiration peak times

% Parse varargin
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'vera_zamoscik')
            vera = true;
        end
    end
end


% Set default for optional arguments
if ~exist('vera', 'var'); vera = false; end

% If we want to use Vera Zamoscik's method
if vera
    % Initialize a cell array to store peaks for each respiration window
    respPeaks_time_wndw = cell(size(resp_windowed, 1), 1);
    
    % Iterate through all respiration windows
    for i = 1:size(resp_windowed, 1)
        % Store window start for this window
        windowStart = resp_windowed{i, 1};
        
        % Store peak times
        [validPeakTimes, minTimeBetween] = ...
            vera_zamoscik(resp_windowed{i, 3}, Fs);
        
        % Shift by the window start time to revert back to Biopac time
        respPeaks_time_wndw{i} = validPeakTimes + (windowStart + (1/Fs));
    end
    
    % Now, we need to take this cell array off peaks and reduce it down to
    % a single column vector of peak times
    respPeakTimes = [];
    for i = 1:size(respPeaks_time_wndw, 1)
        for j = 1:length(respPeaks_time_wndw{i})
            respPeakTimes = [respPeakTimes; respPeaks_time_wndw{i}(j)];
        end
    end
    
    % Remove any peaks that are within a specified tolerance of each other
    % Here, we use Fs to our advantage when defining a time tolerance
    % (read uniquetol documentation to see why we divide by max when defining
    % tolerance)
    tolerance = (1/10)*(1/Fs)/max(respPeakTimes);
    respPeakTimes = uniquetol(respPeakTimes, tolerance);
    
    % Finally, before we find the corresponding onset times, we need to
    % ensure that we are still abiding by the minTimeBetween used in the
    % vera_zamoscik function - this minimum time between can get violated
    % when using overlapping windows. 
    
    % Since we likely would like to keep
    % the peak of greater value, we first need to find respPeakValues
    [~, respPeakValues] = getSampleIndices(respPeakTimes, ...
        resp_signal(:, 2), resp_signal(:, 1), Fs);
    
    for i = length(respPeakTimes):-1:2
        if (respPeakTimes(i) - respPeakTimes(i-1)) < minTimeBetween
            % If this is the case, we should only keep the time and value
            % associated with the greater peak
            if respPeakValues(i) > respPeakValues(i-1)
                respPeakTimes(i-1) = [];
                respPeakValues(i-1) = [];
            else
                respPeakTimes(i) = [];
                respPeakValues(i) = [];
            end
        end
    end
    
    % Now we can determine the breath onsets using the rest of Vera
    % Zamoscik's method, which just specifies that the minimum value in
    % between any two peaks is the respiration onset
    % Because onsets are defined in this way, we will always have one less
    % onset than the number of peaks
    respOnsetTimes = zeros(length(respPeakTimes)-1, 1);
    
    % The trick here is noticing that for the ith element of
    % respOnsetTimes, it needs to be in between the ith and (i+1)th
    % elements of respPeakTimes
    for i = 1:length(respOnsetTimes)
        % Isolate window of time and data corresponding to this onset
        onsetSearchWindow = ...
            resp_signal(resp_signal(:, 1) < respPeakTimes(i+1) & ...
            resp_signal(:, 1) > respPeakTimes(i), :);
        
        % Find the index corresponding to the minimum
        [~, indx] = min(onsetSearchWindow(:, 2));
        
        % Find the time corresponding to this index and store as our result
        respOnsetTimes(i) = onsetSearchWindow(indx, 1);
    end
    
else
    % Otherwise, use count-adv method
    
    % Our first step will be to identify all respiration onsets within each
    % window of respiration data
    % Initialize a cell array that will store vectors of varying lengths
    % containing the breath onset times (varying lengths because windows will
    % not necessarily have the same number of onsets)
    respOnsets_time_wndw = cell(size(resp_windowed, 1), 1);
    respPeaks_time_wndw = cell(size(resp_windowed, 1), 1);
    
    % For each window, we will apply the count-adv algorithm to return the set
    % of detected onsets for each window of respiration data
    for i = 1:size(resp_windowed, 1)
        % Store window start for this window
        windowStart = resp_windowed{i, 1};
        
        % Apply count-adv algorithm; we ignore peaks during segmentation
        [validPeakTimes, validOnsetTimes] = count_adv(resp_windowed{i, 3}, Fs);
        
        % Shift onset times to reference time
        respOnsets_time_wndw{i} = validOnsetTimes + (windowStart + (1/Fs));
        respPeaks_time_wndw{i} = validPeakTimes + (windowStart + (1/Fs));
        
    end
    
    % Now, we need to take these cell arrays and reduce them down to single
    % column vectors; we will then remove repeats
    respOnsetTimes = [];
    respPeakTimes = [];
    for i = 1:size(respOnsets_time_wndw, 1)
        for j = 1:length(respOnsets_time_wndw{i})
            respOnsetTimes = [respOnsetTimes; respOnsets_time_wndw{i}(j)];
        end
        for j = 1:length(respPeaks_time_wndw{i})
            respPeakTimes = [respPeakTimes; respPeaks_time_wndw{i}(j)];
        end
    end
    
    % Remove any onsets that are within a specified tolerance of each other
    % Here, we use Fs to our advantage when defining a time tolerance
    % (read uniquetol documentation to see why we divide by max when defining
    % tolerance)
    tolerance = (1/10)*(1/Fs)/max(respOnsetTimes);
    respOnsetTimes = uniquetol(respOnsetTimes, tolerance);
    respPeakTimes = uniquetol(respPeakTimes, tolerance);
    
    % For the next step, we need the corresponding respOnsetValues and
    % respPeakValues
    [~, respOnsetValues] = getSampleIndices(respOnsetTimes, ...
        resp_signal(:, 2), resp_signal(:, 1), Fs);
    [~, respPeakValues] = getSampleIndices(respPeakTimes, ...
        resp_signal(:, 2), resp_signal(:, 1), Fs);
    
    % Finally, we remove the less extreme of any extrema pairs that are of the
    % same type, i.e., if there are two mins next to each other with no max to
    % separate them, we remove the larger min (less extreme). The reason this
    % can happen even though the breath detection algorithm doesn't allow it is
    % the windowing approach. If the start and end of a window occur at just
    % the right time, sometimes a breath can get split such that two mins or
    % two maxes will end up next to each other once you put all the windowed
    % peaks and onsets together into one vector
    % First, we deal with peaks with no onset in between
    i = length(respPeakTimes);
    % I realized after implementing the while loop that a for loop
    % would've actually worked...but no point in switching over now
    while i > 1     % (No need to check first element)
        % Find the onsets between
        onsetsBetween = ...
            respOnsetTimes(respOnsetTimes < respPeakTimes(i) & ...
            respOnsetTimes > respPeakTimes(i-1));
        
        % If this array is empty, we have a problem situation
        if isempty(onsetsBetween)
            % Figure out who is larger, the ith element or the (i-1)th
            % element
            if respPeakValues(i) > respPeakValues(i-1)
                % This means that the element I am about to visit is
                % actually smaller than the one I am currently on. Note
                % that if I remove the next element, the present
                % element actually moves down one index, so I'll end up
                % in the same place at the next loop instance, which is
                % exactly what we want to make sure there's not
                % additional maxes back to back :)
                respPeakValues(i-1) = [];
                respPeakTimes(i-1) = [];
                i = i - 1;
                
            else
                % In this case, the next element we are about to visit
                % is actually greater, so we can remove the ith element
                % and just move on
                respPeakValues(i) = [];
                respPeakTimes(i) = [];
                i = i - 1;
            end
        else
            i = i - 1;       % Good for now, so move on
        end
    end
    
    % Now, we deal with onsets with no peak in between
    i = length(respOnsetTimes);
    % I realized after implementing the while loop that a for loop
    % would've actually worked...but no point in switching over now
    while i > 1     % (No need to check first element)
        % Find the peaks between
        peaksBetween = ...
            respPeakTimes(respPeakTimes < respOnsetTimes(i) & ...
            respPeakTimes > respOnsetTimes(i-1));
        
        % If this array is empty, we have a problem situation
        if isempty(peaksBetween)
            % Figure out who is smaller, the ith element or the (i-1)th
            % element
            if respOnsetValues(i) < respOnsetValues(i-1)
                % This means that the element I am about to visit is
                % actually larger than the one I am currently on. Note
                % that if I remove the next element, the present
                % element actually moves down one index, so I'll end up
                % in the same place at the next loop instance, which is
                % exactly what we want to make sure there's not
                % additional maxes back to back :)
                respOnsetValues(i-1) = [];
                respOnsetTimes(i-1) = [];
                i = i - 1;
                
            else
                % In this case, the next element we are about to visit
                % is actually smaller, so we can remove the ith element
                % and just move on
                respOnsetValues(i) = [];
                respOnsetTimes(i) = [];
                i = i - 1;
            end
        else
            i = i -1;       % Good for now, so move on
        end
    end
end

end

