function [Rpeak_times, Rpeak_values, Rpeak_good_v_bad] = ...
    Rdetect_ECG(ecg, Fs, HRVparams, SQI_js_threshold,SQI_jw_threshold)
% Takes ECG signal as input and outputs R peak times and associated SQI
% values
% ------- Inputs -------%
% ecg: is two column array with first column as time and second column as
% signal
% Fs: sampling frequency
% HRVparams: parameters set for PhysioNet Cardiovascular Signal Toolbox
% subjectID: string variable containing subject ID (for physionet toolbox)

% --------Outputs -------%
% Rpeak_times: vector storing Biopac times associated with ECG R peaks
% Rpeak_values: vector storing ecg values at each R peak time
% Rpeak_good_v_bad: 1/0 score for good/bad R peaks

% To do this, we will first bandpass filter the ECG signal. We will then
% use the Physionet Cardiovascular Signal Toolbox for QRS detection and
% SQI-based beat scoring for signal quality purposes.


%% Band pass filter ECG
Fstop1 = 0.2;    % Left endpoint of stopband for left cutoff frequency
Fstop2 = 40.5;   % Right endpoint of stopband for right cutoff frequency
Fpass1 = 0.6;    % Beginning of passband
Fpass2 = 40;     % Ending of passband

% Default parameters used in the lab
Dstop1 = 0.001;         % First stopband attenuation
Dpass = 0.05;           % Passband ripple
Dstop2 = 0.0001;        % Second stopband attenuation (HF needs stronger attenuation)
dens = 20;              % Density factor

% Calculate the order from the parameters given for a Parks-Mclellan filter
[N, Fo, Ao, W] = firpmord(...
    [Fstop1, Fpass1, Fpass2, Fstop2]/(Fs/2), ...
    [0, 1, 0], [Dstop1, Dpass, Dstop2]);

% Calculate the coefficients using a parks-mclellan design
b = firpm(N, Fo, Ao, W, {dens});

% Create digital filter
digFilt_ecg = dfilt.dffir(b);

% Apply filter and output result
% (previous lab members were using filtfilthd, so I continue the
% tradition...though we can probably just switch to filtfilt...)
ecg_filtered = ecg;  % Initialize (and store time in one go)
ecg_filtered(:, 2) =filtfilthd(digFilt_ecg, ecg(:, 2));


%% R peak detection
% Detect R peaks and return SQI values corresponding to windows
% I don't need the rr intervals and the times associated (I can get that
% myself using diff()). I also will use SQI_js instead of SQI_jw
[~, ~, R_sampleIndx_shifted, SQI_jw, SQI_js, SQI_jw_windowStarts_shifted, SQI_js_windowStarts_shifted] = ...
    Convert_RawDataToRRIntervals(ecg_filtered(:, 2), HRVparams);

% Convert sample indices to time and then add the data start time
% We subtract one because I usually start my time vectors at t = 0
Rpeak_times = ((R_sampleIndx_shifted-1)./HRVparams.Fs) + ecg_filtered(1, 1);

% Just for debugging purposes, I will also store Rpeak_values
Rpeak_values = zeros(size(R_sampleIndx_shifted));
for i = 1:length(Rpeak_values)
    Rpeak_values(i) = ecg_filtered(R_sampleIndx_shifted(i), 2);
end

% add start time to get Biopac time
SQI_js_window_starts = SQI_js_windowStarts_shifted + ecg_filtered(1, 1);
SQI_jw_window_starts = SQI_jw_windowStarts_shifted + ecg_filtered(1, 1);


%% Assign SQI js and jw to each R peak
% js: jqrs (Pan-Tompkins) vs. sqrs algorithm
% jw: jqrs vs. wqrs algorithm

% I need to assign an SQI_js and SQI_jw score for each R peak

% SQI_js
SQI_js_Rpeaks = zeros(size(Rpeak_times));
for i = 1:length(Rpeak_times)
    % Search for the time window associated with this R peak
    windowIndx = nan;  % Initialize
    for j = 1:length(SQI_js_window_starts)
        % If we are at the end of our search, we just need to check if the
        % R peak falls after this start time
        if j == length(SQI_js_window_starts)
            % If so, then assign
            if Rpeak_times(i) >= SQI_js_window_starts(j)
                windowIndx = j;
            end
            % Otherwise, check if R peak lies in between
        else
            if Rpeak_times(i) >= SQI_js_window_starts(j) && Rpeak_times(i) < SQI_js_window_starts(j+1)
                windowIndx = j;
            end
        end
    end
    SQI_js_Rpeaks(i) = SQI_js(windowIndx);
end

% SQI_jw
SQI_jw_Rpeaks = zeros(size(Rpeak_times));
for i = 1:length(Rpeak_times)
    % Search for the time window associated with this R peak
    windowIndx = nan;  % Initialize
    for j = 1:length(SQI_jw_window_starts)
        % If we are at the end of our search, we just need to check if the
        % R peak falls after this start time
        if j == length(SQI_jw_window_starts)
            % If so, then assign
            if Rpeak_times(i) >= SQI_jw_window_starts(j)
                windowIndx = j;
            end
            % Otherwise, check if R peak lies in between
        else
            if Rpeak_times(i) >= SQI_jw_window_starts(j) && Rpeak_times(i) < SQI_jw_window_starts(j+1)
                windowIndx = j;
            end
        end
    end
    SQI_jw_Rpeaks(i) = SQI_jw(windowIndx);
end

%% Specify js and jw thresholds for later
% Thresholds found through experimentation
% SQI_js_threshold = 0.5;
% SQI_jw_threshold = 0.7;
%Edit 05-28-2022: moved these outside the function


% % To turn off ECG SQI
% SQI_js_threshold = 0;
% SQI_jw_threshold = 0;

%% Determine whether to use jw or js SQI scores
% I've found jqrs vs. sqrs to be more stable, in general
% However, sometimes it can return errant values
% This will help us use the jqrs vs. wqrs SQI labels instead when
% appropriate

% First condition for using jw SQI is if there is a time interval longer
% than 60 seconds such that js SQI returns SQI < 0.5 throughout its
% entirety. To find this, we will need to create a search mechanism

% Initialize variable to store the longest bad interval
longestBad = 0;

% Initialize variable to indicate when we're checking for interval length
checkingCurrently = false;

% Search through all R peaks for longest bad interval
for i = 1:length(Rpeak_times)
    % If we're not currently checking
    if ~checkingCurrently
        % If we fall below the threshold
        if SQI_js_Rpeaks(i) < SQI_js_threshold
            % Start keeping track of interval length
            startInterval = Rpeak_times(i);
            checkingCurrently = true;
        end
    else
        % Otherwise, we're currently checking
        % If we're no longer below threshold or we've reached the end of
        % the R peak times
        if SQI_js_Rpeaks(i) >= SQI_js_threshold || i == length(Rpeak_times)
            % Keep track of interval end
            endInterval = Rpeak_times(i);
            % Find the length of the interval
            interval_length = endInterval - startInterval;
            % If the interval length is the worst yet, store it
            if interval_length > longestBad
                longestBad = interval_length;
            end
            % Now that we're back above the threshold, we're no longer
            % checking
            checkingCurrently = false;
        end
    end
end

% With that taken care of, we now check if we should be using js or jw SQI
% Second condition for using jw instead of js SQI is if the mean < 0.85
if nanmean(SQI_js_Rpeaks) < 0.85
    useJS = false;
elseif longestBad > 60
    useJS = false;
else
    useJS = true;
end


%% Replace NaN SQI values
% With that decision made, we just need to replace NaNs and then we will
% use the thresholds specified at the beginning of this section to label R
% peaks as good/bad (1/0)

% If we're using SQI_js
if useJS
    % Create the final SQI vector we will use
    SQI_Rpeaks = zeros(size(SQI_js_Rpeaks));
    for i = 1:length(SQI_js_Rpeaks)
        % If the SQI score is nan, search for neighbors
        if isnan(SQI_js_Rpeaks(i))
            % We're going to search for nearest non-NaN value on either
            % side
            % Initialize to 1 in case forward or backward neighbors do not
            % exist
            non_NaN_forward = 1;
            non_NaN_backward = 1;
            
            % Bools to store if we found a value
            forwardFound = false;
            backwardFound = false;
            
            % Initialize searches
            forwardIndx = i;
            backwardIndx = i;
            
            % Forward search
            while ~forwardFound && forwardIndx <= length(SQI_js_Rpeaks)
                if ~isnan(SQI_js_Rpeaks(forwardIndx))
                   forwardFound = true;
                   non_NaN_forward = SQI_js_Rpeaks(forwardIndx);
                end
                forwardIndx = forwardIndx + 1;
            end
            
            % Backward search
            while ~backwardFound && backwardIndx > 0
                if ~isnan(SQI_js_Rpeaks(backwardIndx))
                   backwardFound = true;
                   non_NaN_backward = SQI_js_Rpeaks(backwardIndx);
                end
                backwardIndx = backwardIndx - 1;
            end
            
            % Now, we set the SQI value for this R peak to equal the
            % minimum of the backward and forward values
            SQI_Rpeaks(i) = min([non_NaN_backward, non_NaN_forward]);
            
        % Otherwise, just store value
        else
            SQI_Rpeaks(i) = SQI_js_Rpeaks(i);
        end
    end
    % Set threshold appropriately
    SQI_threshold = SQI_js_threshold;
    
% If we're using SQI_jw
else
    % Create the final SQI vector we will use
    SQI_Rpeaks = zeros(size(SQI_jw_Rpeaks));
    for i = 1:length(SQI_jw_Rpeaks)
        % If the SQI score is nan, search for neighbors
        if isnan(SQI_jw_Rpeaks(i))
           % We're going to search for nearest non-NaN value on either
            % side
            % Initialize to 1 in case forward or backward neighbors do not
            % exist
            non_NaN_forward = 1;
            non_NaN_backward = 1;
            
            % Bools to store if we found a value
            forwardFound = false;
            backwardFound = false;
            
            % Initialize searches
            forwardIndx = i;
            backwardIndx = i;
            
            % Forward search
            while ~forwardFound && forwardIndx <= length(SQI_jw_Rpeaks)
                if ~isnan(SQI_jw_Rpeaks(forwardIndx))
                   forwardFound = true;
                   non_NaN_forward = SQI_jw_Rpeaks(forwardIndx);
                end
                forwardIndx = forwardIndx + 1;
            end
            
            % Backward search
            while ~backwardFound && backwardIndx > 0
                if ~isnan(SQI_jw_Rpeaks(backwardIndx))
                   backwardFound = true;
                   non_NaN_backward = SQI_jw_Rpeaks(backwardIndx);
                end
                backwardIndx = backwardIndx - 1;
            end
            
            % Now, we set the SQI value for this R peak to equal the
            % minimum of the backward and forward values
            SQI_Rpeaks(i) = min([non_NaN_backward, non_NaN_forward]);
            
        % Otherwise, just store value
        else
            SQI_Rpeaks(i) = SQI_jw_Rpeaks(i);
        end
    end 
    % Set threshold appropriately
    SQI_threshold = SQI_jw_threshold;
end

%% Score the Rpeaks with binary good/bad score
% If the sqi falls below threshold, it's bad
Rpeak_good_v_bad = ones(size(Rpeak_times));
for i = 1:length(Rpeak_good_v_bad)
    if SQI_Rpeaks(i) < SQI_threshold
        Rpeak_good_v_bad(i) = 0;
    end
end

%% Label all R peaks within first and last 20 seconds as "bad" (score = 0)
% I parsed out an extra 30 seconds at the beginning and end of each dataset
% anyways, so this should still include an added 10 s. I do this because
% the SQI algorithms somewhat falter at the beginning and end of datasets
for i = 1:length(Rpeak_times)
    % If the R peak time is less than 20 s after start or greater than 20 s
    % before the end
    if Rpeak_times(i) < ecg_filtered(1, 1) + 20 || ...
            Rpeak_times(i) > ecg_filtered(end, 1) - 20
        Rpeak_good_v_bad(i) = 0;
    end
end

end

