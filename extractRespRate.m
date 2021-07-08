function [respRate, flaggedIBIs] = extractRespRate(IBIs, varargin)
% This function takes as input inter-breath intervals and returns
% respiration rate by taking 60/IBI; option also available for flagging
% outliers

% Number of deviations used for outlier removal
% numDevs = 5;    % Updated as of 01/2021
numDevs = 4;    % Updated as of 02/2021

% ----------- Inputs ------------- %
% IBIs: Nx2 array, where first column stores time of second breath and
% second column stores interval length in time (seconds)
% (optional)
% 'flagOutlier' flag indicates that the user wants to flag outlier
% intervals

% ---------- Outputs ----------- %
% respRate: Nx2 array, where first column stores IBI time and second column
% stores 60/IBI (respiration rate in breaths per minute)
% flaggedIBIs: Row indices of IBIs flagged for having outlier respiration
% rates

% Parse varargin
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'flagOutlier')
            flagOutlier = true;
        end
    end
end


% Set default for optional arguments
if ~exist('flagOutlier', 'var'); flagOutlier = false; end

% If we do not want to flag outliers, just set flagged IBIs to an empty
% array
if ~flagOutlier; flaggedIBIs = []; end

% Initialize respiration rate
respRate = zeros(size(IBIs, 1), 2);

% Store time
respRate(:, 1) = IBIs(:, 1);

% Compute respiration rates
% To take inter-breath intervals (stored in seconds) and convert to
% respiration rate (in breaths per minute), we need to divide 60 by IBI
respRate(:, 2) = 60./IBIs(:, 2);


% If we need to flag outliers
if flagOutlier
    % Initialize flaggedIBIs array
    flaggedIBIs = [];
    
    % Initialize a t_intervals array that will help me keep track of times
    % throughout this two stage outlier flagging (essentially a dummy
    % vector)
    t_intervals = IBIs(:, 1);
    
    % Stage 1: Overall outlier flagging
    % Threshold factor determines the +- scale factor used for
    % variability (factor depends on method)
    % For example, for the median absolute deviation approach, threshold factor
    % refers to the number of MAD's used to determine outliers
    [stg1_rmvd, rmvd_bool] = rmoutliers(respRate(:, 2), ...
        'median', 'ThresholdFactor', numDevs);
    
    % Need to flag based on outliers
    % Note that t_intervals is currently the same length as respRate(:, 1)
    for i = length(t_intervals):-1:1
        if rmvd_bool(i)
            % Append this interval index to the flagged IBIs array
            flaggedIBIs = [i; flaggedIBIs];
            
            % Remove this time from the t_intervals array so that we can
            % use it in the second stage to figure out time
            t_intervals(i) = [];
        end
    end
    
    % Stage 2: Moving window outlier flagging
    movingWindow = 30;
    
    % Threshold factor determines the +- scale factor used for variability
    % (factor depends on method)
    % For example, for the median absolute deviation approach, threshold factor
    % refers to the number of MAD's used to determine outliers
    [~, rmvd_bool] = rmoutliers(stg1_rmvd, ...
        'movmedian', movingWindow, 'ThresholdFactor', numDevs);
    
    % Need to flag beats based on outliers
    % In here, we can leverage t_beats that is now the same size as
    % stg1_rmvd because we removed from it in the first stage
    for i = length(t_intervals):-1:1
        if rmvd_bool(i)
            % First, I need to find the associated interval time to flag the
            % correct interval index
            % Note that this function uses signal_values to return values; we
            % don't need any corresponding signal values so I pass time 2x
            % I also pass the value 1 as an approximate Fs since Fs is only
            % used for the tolerance (see function)
            [flagIndx, ~] = ...
                getSampleIndices(t_intervals(i), IBIs(:, 1), ...
                IBIs(:, 1), 1, 'method', 'find');
            
            % Append index of flagged interval
            flaggedIBIs = [flagIndx; flaggedIBIs];
        end
    end
end


end

