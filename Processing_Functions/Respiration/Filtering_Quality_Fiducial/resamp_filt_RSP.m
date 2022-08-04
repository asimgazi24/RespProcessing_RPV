function rsp_result = resamp_filt_RSP(rsp_orig, newFs, varargin)
% Takes as input a respiration signal (along with corresponding time
% values) and a new sampling rate; returns a resultant respiration signal
% (along with corresponding times) sampled at the new frequency and
% filtered
% ------- Inputs ----------%
% rsp_orig: Nx2 array, where the first column contains time and the second
% column contains the associated respiration signal values
% newFs: Desired sampling frequency
% (optional)
% 'sgolay': specifies that the user would like to first apply
% Savitsky-Golay detrending before band pass filtering
% 'median': specifies that the user would like to apply a moving median to
% the data after band pass filtering

% -------- Outputs ----------%
% rsp_result: Mx2 array, where the first column contains time and the second
% column contains the associated respiration signal values, downsampled and
% filtered

% Parse varargin
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'sgolay')
            savitzkyFilt = true;
        elseif strcmp(varargin{arg}, 'median')
            medianFilt = true;
        end
    end
end

% Set default for extra filtering options
if ~exist('savitzkyFilt', 'var'); savitzkyFilt = false; end
if ~exist('medianFilt', 'var'); medianFilt = false; end

% Resample original respiration signal to new sampling frequency
[rsp_resamp(:, 2), rsp_resamp(:, 1)] = ...
    clean_resample(rsp_orig(:, 2), rsp_orig(:, 1), newFs, 'linear');


% If user wants to apply Savitzky-Golay detrending, we do so
rsp_stg1 = rsp_resamp;      % Default (and store time in one go)
if savitzkyFilt
    % In Porges-Bohrer method, they applied a 3rd order Savitzky-Golay
    % filter with window size of 10.5 s
    sgolay_ord = 3;
    sgolay_frame_len = 10.5*newFs;
    
    % Apply filter to find the polynomial best fit
    sgolay_trend = ...
        sgolayfilt(rsp_resamp(:, 2), sgolay_ord, sgolay_frame_len);
    
    % Subtract this polynomial trend from signal for detrending
    rsp_stg1(:, 2) = rsp_resamp(:, 2) - sgolay_trend;
end


% Now we bandpass filter the signal using a passband of [.1, 0.72] Hz
% 0.1 chosen because conservative passband is [0.07, 1] Hz, and I feel that
% breathing at 6 bpm (0.1 Hz) is already slow enough. 0.72 is chosen
% because in Zamoscik et al. 2018, they used minimum distance between peaks
% as 1.4 seconds, so we set our passband to end around there
Fstop1 = 0.07;    % Left endpoint of stopband for left cutoff frequency
Fstop2 = 0.75;   % Right endpoint of stopband for right cutoff frequency
Fpass1 = 0.1;    % Beginning of passband
Fpass2 = 0.72;     % Ending of passband

% Default parameters used in the lab
Dstop1 = 0.001;         % First stopband attenuation
Dpass = 0.05;           % Passband ripple
Dstop2 = 0.0001;        % Second stopband attenuation (HF needs stronger attenuation)
dens = 20;              % Density factor

% Calculate the order from the parameters given for a Parks-Mclellan filter
[N, Fo, Ao, W] = firpmord(...
    [Fstop1, Fpass1, Fpass2, Fstop2]/(newFs/2), ...
    [0, 1, 0], [Dstop1, Dpass, Dstop2]);

% Calculate the coefficients using a parks-mclellan design
b = firpm(N, Fo, Ao, W, {dens});

% Create digital filter
digFilt_rsp = dfilt.dffir(b);

% Apply filter and output result
% (previous lab members were using filtfilthd, so I continue the
% tradition...though we can probably just switch to filtfilt...)
rsp_stg2 = rsp_stg1;  % Initialize (and store time in one go)
rsp_stg2(:, 2) =filtfilthd(digFilt_rsp, rsp_stg1(:, 2));


% Finally, if the user wants to apply a median filter, we take care of that
rsp_result = rsp_stg2;      % Default (and store time in one go)
if medianFilt
    % In Zamoscik et al. 2018, they apply the median filter with a window
    % of 1 second
    rsp_result(:, 2) = movmedian(rsp_stg2(:, 2), 1*newFs);
end

end

