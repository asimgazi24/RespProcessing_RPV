function RQI_output = RQI_windowed(rsp_windowed, Fs, windowTime)
% This function takes in a windowed cell array of respiration data, the
% sampling frequency, and the time length of each window, and outputs a
% corresponding respiration quality index (RQI) for each window of
% respiration data

% ------------- Inputs -------------% 
% rsp_windowed: Nx3 cell array, where the first column of cells contain
% start time, the second column of cells contain window end time, and the
% final column of cells contains row vectors of respiration data
% Fs: sampling frequency, in Hz, of respiration data
% windowTime: number of seconds available in each window

% ------------ Outputs --------------- %
% RQI_output: Nx3 array where the first column will contain window start
% time, the second column will contain window end time, and the third
% column will contain the RQI corresponding to that window of data


% We use two RQIs here: an autocorrelation-based RQI and a FFT-based RQI
% We take the average and return that as our overall RQI

% For the sake of time, this function has only been configured to deal with
% Fs = 4 for now, but if this ever needed to be changed, this could by
% tweaking the PSD and autocorrelation parameters
if abs(Fs - 4) > .01
    disp('Note that the RQI_windowed function is configured only for Fs = 4')
end

% First, let's initialize the RQI_output array with the window times since
% that will essentially be a copy and paste from the cell array of windowed
% data
RQI_output = zeros(size(rsp_windowed, 1), 3);
for i = 1:size(RQI_output, 1)
    RQI_output(i, 1) = rsp_windowed{i, 1};
    RQI_output(i, 2) = rsp_windowed{i, 2};
end


% Now we compute the RQI based on the power spectral density of the window
% The idea behind this RQI is that if respiration is fairly clean and
% consistent over the entire window, then the maximum power in the signal
% will be highly localized to 1-2 frequency bins. If the respiration signal
% is not as clean, this max power will instead be reduced via distribution
% across other frequency bands. So, by finding the ratio of max power to
% total power, you're seeing how dominant a specific window of respiration
% frequencies is within the signal; if not too dominant, then there could
% be noise that's adding frequency content elsewhere that we don't want

% Initialize RQI_psd
RQI_psd = zeros(size(rsp_windowed, 1), 1);

% For each window of data...
for i = 1:size(rsp_windowed, 1)
    % Find the power spectral density
    psd_rsp = periodogram(rsp_windowed{i, 3}, [], windowTime*Fs);
    
    % From here on out, I will assume Fs = 4 in the sake of time (see
    % Birrenkott et al. 2018 TBME paper to figure out how to tweak it for
    % different sampling frequencies)
    
    % Compute the sums necessary for numerator
    sums = zeros(30, 1);
    for k = 3:32
        sums(k-2) = psd_rsp(k) + psd_rsp(k+1);
    end
    
    % Respiration frequencies are only present in elements 3-33
    RQI_psd(i) = max(sums)/sum(psd_rsp(3:33));
end

% We next compute the RQI based on autocorrelation of the signal at lags
% relevant to respiration frequencies. The idea here is similar in that
% periodicity is assumed to be the driving factor behind signal quality

% Initialize RQI based on autocorrelation
RQI_autocorr = zeros(size(rsp_windowed, 1), 1);

% For each window...
for i = 1:size(rsp_windowed, 1)
    % Biased keeps the scaling factor constant across all lags, as is done
    % in the paper
    [autocorrs, lags] = xcorr(rsp_windowed{i, 3}, 'biased');
    
    % Based on human respiration frequencies, we only care about lags Fs to
    % 12*Fs (1/12 is .08333, which is seen to be near the lower cutoff of
    % respiration frequencies)
    max_lag = 12*Fs;
    min_lag = Fs;
    
    RQI_autocorr(i) = ...
        max(autocorrs(lags > (min_lag - 0.1) & lags < (max_lag + 0.1)));
end

% Compute final RQI_output by taking the average of the two RQIs
RQI_output(:, 3) = (1/2)*(RQI_psd + RQI_autocorr);

end

