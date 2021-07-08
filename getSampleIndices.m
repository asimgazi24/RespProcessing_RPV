function [signalSampleIndices, corrValues] = ...
    getSampleIndices(timeSet, signal_values, signal_times, Fs, varargin)
% Given a signal vector, an associated time vector, and a separate time set
% to find indices for within the signal, this function returns signal sample
% indices corresponding to the provided time set

% ----- Inputs -----%
% timeSet: Nx1 vector of times to find corresponding sample indices for
% signal_values: Mx1 vector containing signal values
% signal_times: Mx1 vector containing signal times
% Fs: Sampling frequency of signal (can be approximate --> only used to
% determine tolerance)
% (Optional)
% method: method used to find sample indices
%   -'find': less computationally efficient, but likely more robust
%   -'initialTime': faster, but susceptible to error when datapoints are
%   missing, for example (or you're trying to find indices in a
%   nonuniformly sampled signal)


% ----- Outputs -----%
% signalSampleIndices: Nx1 vector of sample indices corresponding to
% timeSet
% corrValues: Nx1 vector of corresponding values

% Parse the various input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'method'); method = varargin{arg + 1};
        end
    end
end

% Set default method
if ~exist('method', 'var'); method = 'find'; end


% Initialize signalSampleIndices and corrValues
signalSampleIndices = zeros(size(timeSet));
corrValues = zeros(size(timeSet));

% For each time in the timeSet
for i = 1:length(timeSet)
    
    if strcmp(method, 'find')
        % Find the corresponding sample index for the given time within a
        % sample frequency-dependent tolerance
        tol = (1/Fs)/10;  % when finding floats, you need a tolerance
        signalSampleIndices(i) = find(abs(signal_times - timeSet(i)) < tol);
    else
        % Otherwise, we assume no datapoints were dropped of the signal, so
        % we can simply take each timeSet value, find its offset from the
        % initial signal time and then multiply that by Fs to get sample
        % index
        signalSampleIndices(i) = (timeSet(i) - signal_times(1))*Fs + 1;
    end
    
    % Plug in indices to signal to get corresponding values
    corrValues(i) = signal_values(signalSampleIndices(i));
    
end

end

