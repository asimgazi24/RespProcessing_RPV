function RSA_out = RSA_PBmethod(t_RR, RRintervals, windows)
% Computes RSA using the Porges-Bohrer method
% ------- Inputs ----------%
% t_RR: Mx1 time vector corresponding to RR intervals (in seconds)
% RRintervals: Mx1 vector of R peak to R peak intervals (in seconds)
% windows: Nx2 array of t_start in first column and t_end in second column
%          specifying which time windows to produce RSA values for

% --------Outputs -------%
% RSA_out: Nx1 array of RSA values corresponding to window starts and ends

% This method involves the below key steps (Lewis et al, 2012):
% 1. Convert RR intervals into milliseconds (just so that results are
% comparable)
% 2. Resample RR intervals to 2 Hz
% 3. Convolve RR interval time series with 21-point cubic polynomial filter
% 4. Subtract this convolved signal from the resampled RR interval signal
% 5. Take this residual signal and filter it using a Chebychev type 1 BPF
% with bandwidth of 0.12 - 0.4 Hz
% 6. Divide this result into 30-second epochs
% 7. Take the natural log of the variance of each of these epochs
% 8. Compute the mean of these epochs - result will be in ln(ms^2) units

%% Convert RR intervals into milliseconds
% Multiply RR intervals by 1000 since they are currently in seconds
RR_ms = 1000*RRintervals;

%% Iterate through each window provided
% To complete the steps that follow, we will first consider each window of
% RR intervals one by one in a for loop

% Initialize result
RSA_out = zeros(size(windows, 1), 1);

% Iterate
for i = 1:size(windows, 1)
    
    % Store this specific window's start and end times
    start_time = windows(i, 1);
    end_time = windows(i, 2);
    
    % Parse RR intervals to consider only the window specified
    RR_window = RR_ms(t_RR >= start_time & t_RR <= end_time);
    t_window = t_RR(t_RR >= start_time & t_RR <= end_time);
    
    %% Resample RR intervals to 2 Hz
    % MATLAB has built-in resample function, but for 0-mean time series, you
    % need to do a little extra detrending and then "retrending" (see
    % clean_resample function for more details)-->deals with endpoint effects
    
    % Since the paper does not specify the method of interpolation when
    % resampling, we will just go ahead with the MATLAB default linear method
    [RR_resampled, uniformTime] = ...
        clean_resample(RR_window, t_window, 2, 'linear');
    
    
    %% Filter result using Savitzky-Golay polynomial filter for detrending
    % Paper specifies a "21-point cubic polynomial filter"
    % This is likely a Savitzky-Golay 3rd order filter w/ 21 element window
    filt_ord = 3;
    frame_len = 21;
    sgolay_trend = sgolayfilt(RR_resampled, filt_ord, frame_len);
    
    
    %% Subtract polynomial filtered trend from resampled RR intervals
    % i.e., result = RR_resampled - sgolay_trend;
    RR_detrended = RR_resampled - sgolay_trend;
    
    
    %% Apply Chebyche type 1 BPF with pass band [0.12, 0.4] Hz
    % Since paper does not specify the filter order of their BPF, I will
    % just use the following filter design
    
    % Note that because data is sampled at 2 Hz, we do not need to multiply
    % our frequencies by anything to get it into rad/sample, i.e., usually
    % the formula is (desiredFc) / (sampling frequency) * (2) pi rad/sample
    % For us here, it's (desiredFc/2)*2 = desiredFc pi rad/sample (you ignore pi
    % when inputting into digital filter functions)
    passBand = [0.12, 0.4];   % passband of [0.12, 0.4] Hz
    stopBand = [0.11, 0.41];   % Leaving 0.01 Hz on the sides for rolloff
                               % (doubt they went stricter)
    passRipple = 1;    % 1 dB ripple (seems reasonable...)
    stopAtten = 60;    % 60 dB attenuation (.001 scale factor)
    
    % Design Chebyshev filter (for numerical reasons, we take this route of
    % returning second order sections (SOS) form, and then filtering)
    BPF_sos = cheby_BPF_sos(passBand, stopBand, passRipple, stopAtten);
    
    % Second order IIR digital filtering for use with SOS form
    RR_filtered = sosfilt(BPF_sos, RR_detrended);
    
    
    %% Divide this filtered result into 30 second epochs
    % uniformTime will inform us of each element's corresponding time
    % We will use that information to divide our result into 30 second
    % chunks
    
    % Time available is window's end time minus window's start time
    timeAvailable = uniformTime(end) - uniformTime(1);
    
    % Divide available time by 30 s to get number of epochs available
    % By using round, we ignore any leftover sections of data less than 15
    % seconds long at the end or anything like that
    % ALSO, Cheby filter introduces a delay that almost zeroes out entire
    % first 30 second epoch, so we ignore the first 30 s
    numEpochs = round(timeAvailable/30) - 1;
    
    % Initialize vector of ln(var(epoch)) values
    ln_var = zeros(numEpochs, 1);
    
    % Compute ln(variance) for each 30 s epoch, remembering to ignore first
    for j = 1:numEpochs
        epoch_start = uniformTime(1) + j*30;
        epoch_end = epoch_start + 30;
        
        % parse out RR intervals corresponding to this epoch
        RR_epoch = RR_filtered(uniformTime >= epoch_start & ...
            uniformTime <= epoch_end);
        
        % Calculate the natural log of the variance of this epoch
        ln_var(j) = log(var(RR_epoch));
    end
    
    % The RSA value for this time window is the mean of all the epochs in
    % this window
    RSA_out(i) = mean(ln_var);
    
end

end


% Function to take care of designing chebyshev filter in second order
% sections form
function sos = cheby_BPF_sos(passBand, stopBand, passRipple, stopAtten)

% Find lowest sufficient filter order for Chebychev filter
[filt_ord, ~] = cheb1ord(passBand, stopBand, passRipple, stopAtten);

% For numerical reasons, we first return zeros, poles, and gain
[zeros, poles, gain] = cheby1(filt_ord, passRipple, passBand, 'bandpass');

% Then convert to second order sections form
sos = zp2sos(zeros, poles, gain);

end

