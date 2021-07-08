function RRintervals = RRintervals_findGood(Rpeak_times, Rpeak_good_v_bad, atrialFib, removeArrhythm,...
    movingWind,u_thresh, l_thresh,rmthreshfactor)
% Takes R peak times and associated good vs. bad binary scores and outputs
% the good RR intervals and the associated times in a single array
% ------- Inputs -------%
% Rpeak_times: vector full of R peak times
% Rpeak_good_v_bad: vector full of associated good vs. bad binary scores

% --------Outputs -------%
% RRintervals: two column array where the first column is time and second
% column contains RR intervals (note that time is associated with the
% second R peak of each RR interval)

% To do this, I will first use diff() to create RR interval array.
% From there, I will use the good vs. bad scores to filter out bad RR
% intervals.
% I will then use the physiologically sensible infimum of 30 bpm to
% remove any RR intervals that are too long.
% Finally, I will remove outliers

%% Create RR interval array using diff() and concatenate time
% I will later remove bad intervals

% Initialize to hold RR intervals
RRintervals = zeros(length(Rpeak_times) - 1, 2);

% Fill times column
RRintervals(:, 1) = Rpeak_times(2:end);

RRintervals(:, 2) = diff(Rpeak_times);

%% Use the good_v_bad scores to remove bad RR intervals
% Encoded as 1 = good; 0 = bad

% Look through all RR intervals for either RR intervals that have the
% second R peak as bad or any RR intervals where the first R peak is bad
% (we need to iterate backwards because we're removing elements)
for i = size(RRintervals, 1):-1:1
    % If second R peak is bad (and possibly first R peak too)
    if Rpeak_good_v_bad(i+1) == 0
        RRintervals(i, :) = [];
    
    % If first R peak is bad, even though second R peak is good
    elseif Rpeak_good_v_bad(i+1) - Rpeak_good_v_bad(i) == 1
        RRintervals(i, :) = [];
    end
end


%% Use 30 bpm infimum to remove any RR intervals that are too long
% HR > 30 bpm implies that RR intervals mustb be < 2

% Look through all RR intervals for intervals equal to or longer than 2 s
for i = size(RRintervals, 1):-1:1
    if RRintervals(i, 2) >= 2
        RRintervals(i, :) = [];
    end
end


%% Finally, remove outliers
% Since the HR data is generally Gaussian, I will, by default, use +- 4 standard
% deviations from the mean as my thresholds
avg_RR = mean(RRintervals(:, 2));
std_RR = std(RRintervals(:, 2));

% However, if subject has atrial fibrillations and we are instructed to 
% remove arrhythmic beats, I will use median and three scaled median
% absolute deviations (formula copied from MATLAB documentation page for
% rmoutliers function)
med_RR = median(RRintervals(:, 2));
scaled_mad_RR = (-1/(sqrt(2)*erfcinv(3/2)))*...
    median(abs(RRintervals(:, 2) - med_RR));

% If atrial fibrillation and I need to remove arrhythmic beats
if atrialFib && removeArrhythm
    % Use more stringent MAD-based outlier removal method
    upper_thresh = med_RR + u_thresh*scaled_mad_RR;
    lower_thresh = med_RR - l_thresh*scaled_mad_RR;
else
    % Otherwise, use +- 4 SD away from mean
    upper_thresh = avg_RR + u_thresh*std_RR;
    lower_thresh = avg_RR - l_thresh*std_RR;   
end

% Need to backward iterate since I will be using '= [];' to remove rows
for i = size(RRintervals, 1):-1:1
    % If above or below upper or lower thresholds, respectively
    if RRintervals(i, 2) > upper_thresh || RRintervals(i, 2) < lower_thresh
        % Remove row of data (removes time as well)
        RRintervals(i, :) = [];
    end
    
end

% I learned that for AF subjects, we need to do one more stage of outlier
% removal to truly catch most of the arrhythmic beats 
% (look at subject 126, for example). This will involve moving window 
% outlier removal
if atrialFib && removeArrhythm
%     movingWind = 30;  % number of elements in moving window
    
    % Use built-in MATLAB function
    % Threshold factor refers to the number of median absolute deviations
    % used for outlier thresholding
    [~, rmvd_bool] = rmoutliers(RRintervals(:, 2), ...
        'movmedian', movingWind, 'ThresholdFactor', rmthreshfactor);
    
    % Need to remove rows according to the boolean vector
    for i = size(RRintervals, 1):-1:1
        if rmvd_bool(i)
            % (remove time as well, so remove entire row)
            RRintervals(i, :) = [];
        end
    end
    
end


end