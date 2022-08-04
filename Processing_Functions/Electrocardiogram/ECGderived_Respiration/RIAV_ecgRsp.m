function RIAV = ...
    RIAV_ecgRsp(ecg, Fs, t_Rpeaks)
% This function takes as input an ECG signal and a set of R peak times
% and returns a 2 column array containing R peak
% times in the first column and R peak amplitudes as the second column.
% This result represents respiratory-induced amplitude variation (RIAV),
% which is one of the ECG-derived respiration signals that can be extracted

% -------- Inputs ----------- %
% ecg: Nx2 array, where the first column contains time and the second
% column contains corresponding ECG datapoints
% Fs: sampling frequency of ECG signal
% t_Rpeaks: Mx1 array containing R peak times

% -------- Outputs ------------ %
% RIAV: Mx2 array, where the first column will be equivalent to t_Rpeaks,
% and the second column will contain R peak amplitudes


% We use the getSampleIndices() function to find the ECG samples associated
% with each R peak time and the R peak values
[sample_Rpeaks, max_Rpeaks] = getSampleIndices(t_Rpeaks, ...
    ecg(:, 2), ecg(:, 1), Fs, 'method', 'find');

% To find the minimum (a.k.a QRS trough) associated with each R peak, 
% I need to locate the minimum ECG value within .10 s prior to the 
% R peak
min_Rpeaks = zeros(length(sample_Rpeaks), 1);
for i = 1:length(min_Rpeaks)
    min_Rpeaks(i) = ...
        min(ecg((sample_Rpeaks(i) - .1*Fs):sample_Rpeaks(i), 2));
end


% Amplitude is given by R peak minus QRS trough
RIAV(:, 2) = max_Rpeaks - min_Rpeaks;

% First column of output will store corresponding times
RIAV(:, 1) = t_Rpeaks;


%%%% Check for inverted ECG %%%%
% In some cases, the subject was hooked up in a non-standard way, resulting
% in an inverted ECG waveform. For these ECG signals, we need to set the R
% peak as the minimum and the maximum value as the maximum ECG value within 
% .10 s prior to R peak
if nnz(RIAV(:, 2)) < 0.75*size(RIAV, 1)
    % If less than 75 % of the values are nonzero, this likely implies that
    % the ECG was inverted because the minima found were almost the same as
    % the R peak, and therefore, the amplitudes were set to 0 at least 
    % 25 % of the time
    
    % First, set the R peak values to be the minimum values instead
    min_Rpeaks = max_Rpeaks;
    
    % Now find the maximum (a.k.a inverted QRS trough) for each R peak. 
    % I need to locate the maximum ECG value within .10 s prior to the 
    % R peak
    max_Rpeaks = zeros(length(sample_Rpeaks), 1);
    for i = 1:length(max_Rpeaks)
        max_Rpeaks(i) = ...
            max(ecg((sample_Rpeaks(i) - .1*Fs):sample_Rpeaks(i), 2));
    end
    
    % Reset the RIAV values to be the new amplitudes computed
    RIAV(:, 2) = max_Rpeaks - min_Rpeaks;
    
end

end

