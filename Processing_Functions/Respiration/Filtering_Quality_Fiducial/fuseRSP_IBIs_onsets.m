function [fused_IBIs, fused_onsetTimes] = ...
    fuseRSP_IBIs_onsets(numSignals, varargin)
% This function takes a number of respiratory signals' IBIs and onsets
% and fuses these according to the RQI values for each window of
% respiration data

% ------------ Inputs ---------------- %
% numSignals: number of signals input
% 'RQI': lets the function know that RQI arrays are coming up
% The RQI arrays are Nx3 and hold window start times in their first
% column, window end times in their second column, and RQI values across
% those time windows in the third column. You need to input a RQI for each
% of the input signals (order is important! See below)
% 'IBI': lets the function know that IBI arrays are coming up
% The IBI arrays are Mx2 and hold time in the first column and inter-breath
% intervals in the second column
% 'onset': lets the function know that onset times are coming up
% The onset time vectors are Mx1 and hold onset times corresponding to each
% of the inter-breath intervals
% 'ecg-derived': alerts the function that there are ecg-derived signals
% present (so we need to consider bad ECG times) - the bad ECG times and
% which signals are ECG-derived need to be input after this flag

% ---------- Outputs ----------------- %
% fused_IBIs: resultant array of IBIs (no way to know size a priori
% since windows can be entirely removed if all 4 signals have insufficient
% RQI within those windows)
% fused_onsetTimes: resultant vector of onset times, one per IBI

% Initalize cell array of RQI arrays
RQI_windowed = cell(1, numSignals);
IBIs = cell(1, numSignals);
onsets = cell(1, numSignals);

% Parse varargin
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'RQI')
            for i = 1:numSignals
                RQI_windowed{i} = varargin{arg + i};
            end
        end
        
        if strcmp(varargin{arg}, 'IBI')
            for i = 1:numSignals
                IBIs{i} = varargin{arg + i};
            end
        end
        
        if strcmp(varargin{arg}, 'onset')
            for i = 1:numSignals
                onsets{i} = varargin{arg + i};
            end
        end
        
        if strcmp(varargin{arg}, 'ecg-derived')
            ecg_derived = true;
            % Stores bad times for ECG
            badTimes_ecg = varargin{arg+1};
            % Stores which input arrays are ecg-derived
            indices_ecgDer = varargin{arg+2};
        end
    end
else
    disp('fuseRSP_IBIs_onsets() was used improperly')
    return
end

% Set default if user does not include ECG-derived
if ~exist('ecg_derived', 'var'); ecg_derived = false; end

% Initialize results that we will fill via appending
fused_IBIs = [];
fused_onsetTimes = [];

% Initialize RQI threshold for rejection
thresh = 0.45;      % Updated as of 02/2021
% thresh = 0.4;       % Updated as of 01/2021

% What we need to do to pull off this RQI-ranking-based fusion is to
% iterate through the windows of the RQI_windowed arrays. Since they should
% all have the same number of windows, it doesn't matter which I use. Note
% that if ECG-derived signals are included, I will need to add an
% additional step to the process, so I split the two paths here
if ecg_derived
    % ECG-derived signals are included
    for i = 1:size(RQI_windowed{1}, 1)
        % First, I will check to see if any bad ECG times are present in
        % this window. If so, I need to ignore for this iteration
        windowStart = RQI_windowed{1}(i, 1);
        windowEnd = RQI_windowed{1}(i, 2);
        
        % Check to see if any bad ECG times exist in this window
        if ~isempty(badTimes_ecg(badTimes_ecg > windowStart & ...
                badTimes_ecg < windowEnd))
            % If so, then we need to limit the current RQIs, peaks, and/or
            % onsets that we consider
            current_RQIs = zeros((numSignals - length(indices_ecgDer)), 1);
            
            % Isolate the valid signal indices to work with for this windw
            validIndices = setdiff(1:numSignals, indices_ecgDer);
            
            % If we have some valid indices left, then...
            if ~isempty(validIndices)
                for j = 1:length(validIndices)
                    % Store signal index
                    indx = validIndices(j);
                    
                    % Store RQIs for this window separately for easy access
                    current_RQIs(j) = RQI_windowed{indx}(i, 3);
                end
                
                % We need to make sure that from the remaining RQIs, at
                % least one is greater than the threshold. Otherwise, we
                % need to move onto the next window
                if max(current_RQIs) >= thresh
                    % Now that we're finally here, all that we have to do
                    % is store the fiducial points corresponding to the
                    % signal with highest RQI in this window
                    
                    % Let's first figure out which RQI is the highest in
                    % this window, i.e., which signal to draw from
                    [~, indx] = max(current_RQIs);
                    bestIndx = validIndices(indx);
                    
                    % To ensure that we are always keeping IBIs
                    % associated with onsets, we need to also use IBI row
                    % indices when getting onset times
                    
                    % Append all IBIs that are within specified
                    % window from the best signal
                    fused_IBIs = [fused_IBIs; ...
                        IBIs{bestIndx}(...
                        IBIs{bestIndx}(:, 1) >= windowStart & ...
                        IBIs{bestIndx}(:, 1) <= windowEnd, ...
                        :)];
                    
                    % Append all onset times that are within specified
                    % window from the best signal
                    fused_onsetTimes = [fused_onsetTimes; ...
                        onsets{bestIndx}(...
                        IBIs{bestIndx}(:, 1) >= windowStart & ...
                        IBIs{bestIndx}(:, 1) <= windowEnd)];
                end
            end
        else
            % First, I will make sure that at least one of the RQIs for this window
            % is greater than thresh; otherwise, we need to ignore and move onto
            % next window. To figure this out, we first find the max RQI
            current_RQIs = zeros(numSignals, 1);
            for j = 1:numSignals
                % Store RQIs for this window separately for easy access
                current_RQIs(j) = RQI_windowed{j}(i, 3);
            end
            
            % Only perform the remaining steps if a RQI exists that is greater
            % than the threshold
            if max(current_RQIs) >= thresh
                % Now that we're here, all that we have to do
                % is store the fiducial points corresponding to the
                % signal with highest RQI in this window
                
                % Let's first figure out which RQI is the highest in
                % this window, i.e., which signal to draw from
                [~, bestIndx] = max(current_RQIs);
                
                % Append all IBIs that are within specified
                % window from the best signal
                fused_IBIs = [fused_IBIs; ...
                    IBIs{bestIndx}(...
                    IBIs{bestIndx}(:, 1) >= windowStart & ...
                    IBIs{bestIndx}(:, 1) <= windowEnd, ...
                    :)];
                
                % Append all onset times that are within specified
                % window from the best signal
                fused_onsetTimes = [fused_onsetTimes; ...
                    onsets{bestIndx}(...
                    IBIs{bestIndx}(:, 1) >= windowStart & ...
                    IBIs{bestIndx}(:, 1) <= windowEnd)];
            end
        end
    end
else
    % No ECG-derived signals are present
    for i = 1:size(RQI_windowed{1}, 1)
        % store window start and end times for later
        windowStart = RQI_windowed{1}(i, 1);
        windowEnd = RQI_windowed{1}(i, 2);
        
        % First, I will make sure that at least one of the RQIs for this window
        % is greater than thresh; otherwise, we need to ignore and move onto
        % next window. To figure this out, we first find the max RQI
        current_RQIs = zeros(numSignals, 1);
        for j = 1:numSignals
            % Store RQIs for this window separately for easy access
            current_RQIs(j) = RQI_windowed{j}(i, 3);
        end
        
        % Only perform the remaining steps if a RQI exists that is greater
        % than the threshold
        if max(current_RQIs) >= thresh
            % Now that we're here, all that we have to do
            % is store the fiducial points corresponding to the
            % signal with highest RQI in this window
            
            % Let's first figure out which RQI is the highest in
            % this window, i.e., which signal to draw from
            [~, bestIndx] = max(current_RQIs);
            
            % Append all IBIs that are within specified
            % window from the best signal
            fused_IBIs = [fused_IBIs; ...
                IBIs{bestIndx}(...
                IBIs{bestIndx}(:, 1) >= windowStart & ...
                IBIs{bestIndx}(:, 1) <= windowEnd, ...
                :)];
            
            % Append all onset times that are within specified
                % window from the best signal
                fused_onsetTimes = [fused_onsetTimes; ...
                    onsets{bestIndx}(...
                    IBIs{bestIndx}(:, 1) >= windowStart & ...
                    IBIs{bestIndx}(:, 1) <= windowEnd)];
        end
    end
end

% Ensure that we did not inadvertently store repeated times
% Read uniquetol documentation to figure out why we divide by the maximum
% value when defining our tolerance of 1/10 s (since we're dealing with
% breathing here, there's no chance someone breathed twice within 100 ms,
% so anything within that tolerance is likely a repeat
tolerance = (1/10)/max(fused_onsetTimes);
[fused_onsetTimes, rowIndices , ~] = uniquetol(fused_onsetTimes, tolerance);
% Use the indices returned above to also get unique values for fused_IBIs
fused_IBIs = fused_IBIs(rowIndices, :);


end

