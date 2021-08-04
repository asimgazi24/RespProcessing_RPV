%main processing code for the GUI

function gui_process(dataFileLocation,dataFileName, removeArrhythm, isi,...
    SQI_js_threshold,SQI_jw_threshold,u_thresh, l_thresh,rmthreshfactor,...
    windowTime_rsp_peak,windowShift_rsp_peak, windowTime_rsp_RQI, windowShift_rsp_RQI,...
    RQI_thresh,windowTime_RPV, windowShift_RPV, max_IBI, min_IBI,PeakDetectMethod)

    % Create arrays to store data removal statistics
    % (ordered according to ordering in subject_IDs)
    percentRmvd_array_ecg = -1*ones(1, 1);
    percentRmvd_array_rsp = -1*ones(1, 1);

    % Create arrays to store RPV windows removed for each resp. time series
    RPVwindowsRmvd_array_RR = -1*ones(1, 1);
    RPVwindowsRmvd_array_Ti = -1*ones(1, 1);
    RPVwindowsRmvd_array_Te = -1*ones(1, 1);

    % In case we want to use a parfor loop, we will operate
    % with transparency in mind when loading the data
    data_struct = load([dataFileLocation, filesep, dataFileName]);

    ecg = data_struct.ecg;
    rsp = data_struct.rsp;
    
    %% Initialize HRVparams for PhysioNet Cardiovascular Signal Toolbox
    % Because we may use a parfor looop, we need to initialize the struct
    % within each iteration of the loop
    HRVparams = InitializeHRVparams('03-2020');
    
    % ECG SQI
    HRVparams.sqi.windowlength = 10;   % length of time considered for each SQI label
    HRVparams.sqi.increment = 10;       % number of seconds the window moves each time
    HRVparams.sqi.TimeThreshold = 0.1; % number of seconds allowed to differ between R peak labels
    HRVparams.sqi.margin = 0;          % number of seconds ignored at the beginning of window
    
    % R peak detection
    HRVparams.PeakDetect.REF_PERIOD = 0.33;  % min number of seconds between R peaks
    HRVparams.PeakDetect.THRES = 0.6;  % Energy threshold for Pan-Tompkins algorithm
    HRVparams.PeakDetect.ecgType = 'MECG';  % Default type of ECG is adult
    HRVparams.PeakDetect.windows = 10;  % window size (s) to perform QRS detection
    HRVparams.PeakDetect.SIGN_FORCE = 1; % Only detect positive peaks
    
    % Heart Rate Variability (HRV) Window Settings
    HRVparams.windowlength = 300;    % Number of seconds for each HRV window
    HRVparams.increment = 1;         % Seconds to increment the window each time
    
    % HRV normal to normal (NN) Settings
    HRVparams.preprocess.per_limit = 0.2;  % Max. percent change in neighboring NN intervals
    HRVparams.preprocess.gaplimit = 2;   % Longest RR interval allowed
    HRVparams.preprocess.lowerphysiolim = 60/180;   % minimum RR interval length
    HRVparams.preprocess.upperphysiolim = 60/30;    % maximum RR interval length
    
    % Missing data and rejection thresholds are changed below because
    % this just stops the code from being dumb since my knowledge of the
    % data informs me that there is nothing to worry about (in every 300 s
    % window, I would estimate there to be at least 200 datapoints)
    % I also remove most of the arrhythmic beats from the subjects that
    % experience atrial fibrillation (outlier removal), so I will just analyze
    % as is and then exclude those subjects later if I want
    % (if I turn AF back on, it will decide to exclude ENTIRE WINDOWS of
    % data if it classifies that window of data as arrhythmic)
    HRVparams.MissingDataThreshold = 1;  % proportion of window allowed to be missing
    HRVparams.RejectionThreshold = 1;   % proportion of window rejectable by HRV preprocessing
    HRVparams.af.on = 0;  % Turned off since I want to analyze all data, regardless of AF
    % (I have my own boolean variable for excluding most of arrhythmic data
    % anyways)
    
    % HRV time domain settings
    HRVparams.timedomain.on = 1;    % We want time domain features
    HRVparams.timedomain.alpha = 50;  % We want pNN50 (not some alpha ~= 50)
    
    % HRV frequency domain settings
    HRVparams.freq.on = 1;          % We want frequency domain features
    HRVparams.freq.method = 'lomb';  % We will use the recommended Lomb periodogram
    % Lomb periodogram allows nonuniform samp.
    HRVparams.freq.normalize_lomb = 0;  % No need to normalize since I'm only using LF/HF ratio
    
    % HRV Poincare analysis settings
    HRVparams.poincare.on = 1;
    
    % Turn everything else off to speed up the code
    % (if these are turned on, then just update HRV code as well)
    HRVparams.sd.on = 0;
    HRVparams.prsa.on = 0;
    HRVparams.MSE.on = 0;
    HRVparams.Entropy.on = 0;
    HRVparams.DFA.on = 0;
    HRVparams.HRT.on = 0;
    
    
    % Store sampling frequency as reciprocal of intersample interval
    % multiply by 1000 b/c isi_units is in milliseconds
    Fs = (1/isi)*1000;
    
    HRVparams.Fs = Fs;
    
    %% Concatenate a time vector to each biosignal for processing purposes
    % In this way, we can keep track of true times based on the sampling
    % frequency
    [ecg, rsp] = concat_time(Fs, ecg, rsp);

    
     %% Initialize time series, HRV, and respV data structures
    % Throughout the rest of the script, we will be extracting various
    % physiological features. Some will be extracted each heart beat (i.e.,
    % "beat-by-beat"), others will be extracted across windows (HRV and RPV),
    % and the respiration time series features will be extracted each breath
    physFeatures_beat = struct();
    physFeatures_HRV = struct();
    physFeatures_resp = struct();
    physFeatures_RPV = struct();
    
    %% ECG Processing
    % We will begin with ecg processing, the reason being that rsp
    % processing includes ECG-derived respiration. Note that since this
    % processing will involve SQI-based rejection for ECG, the associated
    % sample / time intervals for those thrown away beats need to be stored
    % as well so we can ignore those ECG-derived respiration windows
    
    % ------- Find Rpeaks and associated good vs. bad scores --------- %
    [Rpeak_times, Rpeak_values, Rpeak_good_v_bad] = ...
        Rdetect_ECG(ecg, Fs, HRVparams,SQI_js_threshold,SQI_jw_threshold);
    
    %Assumption: if removeArrythm, then atrialFib
    atrialFib = removeArrhythm;
    %------------ Find good RR intervals -----------------%
    % Using Rpeak times and the good vs bad scores, we can now find the good RR
    % intervals and use those not only for HR and HRV extraction
    RRintervals = RRintervals_findGood(...
        Rpeak_times, Rpeak_good_v_bad, atrialFib, removeArrhythm,...
         30,u_thresh, l_thresh,rmthreshfactor);
    
    % NOTE THAT RRintervals CONTAINS Nx2 ARRAY, WHERE FIRST COLUMN
    % IS TIME (OF 2ND PEAK) AND SECOND COLUMN CONTAINS THE ACTUAL RR INTERVALS
    
    
    %----------- Compute heart rate (HR) from RR intervals --------------%
    % HR = 60 s/min divided by the RRinterval (in s)
    % Initialize HR array within the time series struct
    physFeatures_beat.HR = zeros(size(RRintervals));
    
    % Iterate through all rows
    for j = 1:size(RRintervals, 1)
        % Copy over time column
        physFeatures_beat.HR(j, 1) = RRintervals(j, 1);
        % HR = 60/RRinterval
        physFeatures_beat.HR(j, 2) = 60/RRintervals(j, 2);
    end
    
    
    
    %------------ Extract HRV features --------------%
    % We will use PhysioNet Cardio Toolbox and then some custom code to compute
    % RSA using Porges-Bohrer method
    % First, we use the Toolbox for most of the features we want
    [HRVout, ~, HRVnames] = Main_HRV_Analysis(RRintervals(:, 2), ...
        RRintervals(:, 1), 'RRIntervals', HRVparams);
    
    % Output includes NaN for all windows that analysis could not be
    % performed (toolbox automatically starts from 1-300 s, even if
    % your data starts from later time). We need to remove all rows
    % that contain nans.
    HRVout(any(isnan(HRVout), 2), :) = [];
    
    % Next, we need to exclude the rows that correspond to time windows
    % that do not make sense, i.e., their start times come before our
    % heart rate data start or their end times fall after our heart
    % rate data end
    for j = size(HRVout, 1):-1:1
        % If the start time comes before the heart rate's start time
        % or if the end time comes before the heart rate's final time
        if HRVout(j, 1) < physFeatures_beat.HR(1, 1) ||...
                HRVout(j, 2) > physFeatures_beat.HR(end, 1)
            % Remove row
            HRVout(j, :) = [];
        end
    end
    
    % With those rows removed, we now need to extract only the columns of
    % features we want to keep. Third column of each struct field will
    % be for data; first two columns will be t_start and t_end for
    % HRV windows. (check HRVnames for column labels for HRV out)
    physFeatures_HRV.RMSSD(:, 3) = HRVout(:, 11);
    physFeatures_HRV.pNN50(:, 3) = HRVout(:, 12);
    physFeatures_HRV.lfhf(:, 3) = HRVout(:, 20);
    physFeatures_HRV.sd1sd2(:, 3) = HRVout(:, 25);
    
    % Just to initialize RSA and store the necessary window starts and
    % end times to compute RSA for, I will initialize that struct with
    % garbage values. This will allow me to do the next step in one
    % fell swoop... (see next)
    physFeatures_HRV.RSA(:, 3) = HRVout(:, 1);
    
    % Loop over all field names to store first two columns as window
    % start and end times
    labels_HRV = fieldnames(physFeatures_HRV);
    for j = 1:numel(labels_HRV)
        % t_start
        physFeatures_HRV.(labels_HRV{j})(:, 1) = HRVout(:, 1);
        % t_end
        physFeatures_HRV.(labels_HRV{j})(:, 2) = HRVout(:, 2);
    end
    
    % For the final respiratory sinus arrhythmia (RSA) feature, we use
    % custom code based on the Porges-Bohrer method, leveraging the
    % fact that we already have stored the time windows we need to
    % compute RSA for in a sliding window fashion (see above)
    physFeatures_HRV.RSA(:, 3) = ...
        RSA_PBmethod(RRintervals(:, 1), RRintervals(:, 2), ...
        physFeatures_HRV.RSA(:, 1:2));
    
    
    
    %----------Store clean R peaks ---------%
    % Since these are used later, might as well
    % store them separately for easy access
    
    % Use RRintervals to work backwards to cleanRpeak_times
    % If RRintervals is Nx2, we need a (N+1)x1 array
    % this is because RRintervals has the time of the second R peak of each
    % RR interval stored in the first column and then the actual RR
    % interval time length in the second column - this means that for 100 R
    % peaks, for example, you get 99 RR intervals
    clean_Rpeak_times = zeros(size(RRintervals,1)+1, 1);
    
    
    % For the first R peak all the way to the second to last one
    % Recall that column 1 (time) corresponds to time of second R peak
    % of each RR interval
    for j = 1:size(RRintervals, 1)
        % R_1 = t_1 - RR_1
        % R_2 = t_2 - RR_2
        % ...
        % R_N = t_N - RR_N
        clean_Rpeak_times(j) = ...
            RRintervals(j, 1) - RRintervals(j, 2);
    end
    
    % Take care of last R peak
    clean_Rpeak_times(end) = RRintervals(end, 1);
    
    
    % Use custom function to evaluate original ecg signal at clean_Rpeak_times
    % to get clean_Rpeak_values
    % Function takes as input R peak times, ecg signal, corresponding
    % sample times, and sampling frequency; outputs R peak values
    [~, clean_Rpeak_values] = getSampleIndices(clean_Rpeak_times, ...
        ecg(:, 2), ecg(:, 1), Fs);
    
    
    % Just in case I want to know how much ECG data I used/removed, I can
    % leverage the clean RR intervals by creating a new vector of clean ECG
    % times based on the RR intervals, their associated times, and the
    % sampling frequency
    cleanTimes_ecg = [];
    for j = 1:size(RRintervals, 1)
        % Concatenate as a single column vector
        % Store all times from the first R peak time to the
        % second R peak time for each RR interval
        % That's why we increment by the sample time (1/Fs)
        cleanTimes_ecg = [cleanTimes_ecg; ...
            ((RRintervals(j, 1) - RRintervals(j, 2)):...
            1/Fs:RRintervals(j, 1))'];
        
        % Eliminate repeated times
        cleanTimes_ecg = unique(cleanTimes_ecg);
    end
    
    % Compute the percent of samples removed from ECG
    percentRmvd_ecg = 100*(1 - ...
        (length(cleanTimes_ecg)/size(ecg, 1)));
    
    
    
    % ----- Store corrupt ECG times to ignore corresponding windows ---- %
    % An issue we run into later is that because we resample ECG-derived
    % respiration signals, we get datapoints for RIIV, RIFV, and RIAV
    % where we shouldn't, i.e., when the ECG data was corrupt (based on ECG bad
    % data removal scheme), we originally do not have any ecg-derived resp.
    % datapoints there; but after resampling, we do. So, by figuring out which
    % ECG times corresponded to corrupt RR intervals, we can later ignore any
    % windows of ECG-derived respiration which contain these bad data times
    
    % Find ecg{i}(:, 1)\cleanTimes_ecg{i}, i.e., find the set
    % difference to extract all elements of the original time vector
    % that are not in the clean times
    % To deal with floating point tolerance issues, we scale everything
    % by 10*Fs and then rescale back
    badTimes_ecg = setdiff(round(10*Fs*ecg(:, 1)), ...
        round(10*Fs*cleanTimes_ecg))...
        ./ (10*Fs);
    
    
    %% Respiration Processing
    % We leverage both the respiration belt and ECG-derived respiration
    % ----------- Extract ECG-Derived Respiration Signals -------------- %
    % RIFV
    % We extract RIFV first because it will be the easiest to extract and
    % should not have been altered at all by the factory setting filters
    % applied during data collection. RIFV is based on respiratory sinus
    % arrhythmia, where HR increases during inspiration and decreases during
    % expiration. By using the clean R peak times and RR intervals we extracted
    % earlier from the ECG signal, we can extract our RIFV signal by setting
    % the clean R peak times as the independent variable and the RR intervals
    % as the dependent variable, i.e., RR vs. R peaks
    % (note that this is actually what RRintervals contains)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ecgDeriveRsp_RIFV = RRintervals;
    %%%%%%%%%%%%%%%%%% We end up ignoring this one!%%%%%%%%%%%%%%%%%%
    
    % RIIV
    % Although RIIV may have been filtered out due to the factory presets on the
    % Biopac sensing system, we try to extract RIIV here by simply considering
    % the clean R peak values vs. the clean R peak times
    % (note that we already have these clean R peak values and times as their
    % own variables, so all we would need is to conctenate to form a Mx2 signal
    % array)
    % Just to make the dimensions match between RIIV and RIFV, I will
    % ignore the first R peak time and value
    ecgDeriveRsp_RIIV = ...
        [clean_Rpeak_times(2:end), clean_Rpeak_values(2:end)];
    
    % RIAV
    % RIAV will require going back to the ECG signal to find QRS troughs within
    % 0.1 s prior to R peaks (Ruangsuwana et al 2010). We can then
    % calculate R peak amplitudes and set RIAV to be these R peak
    % amplitudes vs. clean R peak times. This signal may have also been
    % attenuated by the factory preset bandpass filter
    % Just to make the dimensions match between RIIV and RIFV, I will
    % ignore the first R peak time and value
    % (note that Fs below refers to the sampling frequency of ecg)
    ecgDeriveRsp_RIAV = RIAV_ecgRsp(ecg, Fs, ...
        clean_Rpeak_times(2:end));
    
    % To keep the code clean, we will now define a cell array to store all of
    % our respiration signals, within another cell array that keeps before and
    % after lunch separate
    
    rspSignals = cell(1, 3);
    
    % Store respiration signals in their appropriate cells
    rspSignals{1} = rsp;
    rspSignals{2} = ecgDeriveRsp_RIIV;
    rspSignals{3} = ecgDeriveRsp_RIAV;
    
    
    % --------------- Resample and Filter ------------------- %
    % With ECG-Derived Respiration in hand, as well as the collected respiration
    % belt signal, we now resample the signals to 4 Hz and apply a bandpass
    % filter to attenuate frequency components outside human respiration
    % frequencies (we consider minimum time between breaths to be 1.4 s and
    % maximum time to be 10 s - see resamp_filt_RSP() function for more
    % details). If desired, the function we will use also allows for
    % Savitzky-Golay detrending and median filtering. 4 Hz is chosen to match
    % ECG-derived respiration literature from David Clifton's group (check
    % Birenkott et al. paper on respiration quality indices (RQI))
    
    % I found resampling to 4 Hz causes issues with later estimating time
    % intervals because any changes < 0.25 s are not picked up; I will instead
    % resample to 50 Hz for a separate set of signals used in peak detection
    
    % Based on RQI literature, I resample to 4 Hz
    newFs_rsp_RQI = 4;
    
    % Based on Vera Zamoscik's paper, I will resample to 50 Hz
    newFs_rsp_pk = 50;
    
    % As I will do throughout the remainder of this code, I will simply
    % loop through the respiration signals I've saved in the cell
    % array initialized above
    rspSignals_resamp_filt_RQI = cell(1, numel(rspSignals));
    rspSignals_resamp_filt_pk = cell(1, numel(rspSignals));
    
    for j = 1:numel(rspSignals)
        rspSignals_resamp_filt_RQI{j} = resamp_filt_RSP(...
            rspSignals{j}, newFs_rsp_RQI);
        rspSignals_resamp_filt_pk{j} = resamp_filt_RSP(...
            rspSignals{j}, newFs_rsp_pk);
    end
    
    
    % At this point, things get a little less linear than for previous signals,
    % where steps were taken sequentially with no branching off in the pipeline
    % or anything like that. For respiration, we will be branching off into two
    % paths from here, as will be evident by windowing the signal in two
    % separate ways. The first path involves windowing, breath
    % detection, and inter-breath interval computation. The second path involves
    % windowing and computing respiration quality indices (RQIs) to determine
    % which windows signal's IBIs to choose for each window of data
    % associated with a RQI. From that point, find the onsets associated with
    % the fused IBIs and flag outliers based on respiration rate
    % (RR), inspiration time (Ti), and expiration time (Te). By then removing
    % the associated IBIs and onsets, we arrive at a final set from which we
    % extract our final sets of features.
    
    
    % Notice that right now, the belt signal is not of the same length as the
    % ECG-derived signals, so we need to fix that first for this plan to work.
    % I will take the ECG-derived respiration signal start and end times as
    % reference since these are based on R peaks, and the first and last chunks
    % of ECG data were ignored. Hence, I will chop off the beginning and ending
    % of the respiration belt signal accordingly (to ensure the signals are
    % windowed equivalently in time)
    
    % I could use any of the 3 ECG-derived respiration signals as
    % reference since they all have the same time columns
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (may need changing if user does not use ECG-derived; check twin study
    % script to be sure)
    rsp_startTime_RQI = rspSignals_resamp_filt_RQI{2}(1, 1);
    rsp_endTime_RQI = rspSignals_resamp_filt_RQI{2}(end, 1);
    rsp_startTime_pk = rspSignals_resamp_filt_pk{2}(1, 1);
    rsp_endTime_pk = rspSignals_resamp_filt_pk{2}(end, 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Isolate the respiration belt signal that corresponds to this
    % start and end time
    rspSignals_resamp_filt_RQI{1} = ...
        rspSignals_resamp_filt_RQI{1}(...
        rspSignals_resamp_filt_RQI{1}(:, 1) >= rsp_startTime_RQI & ...
        rspSignals_resamp_filt_RQI{1}(:, 1) <= rsp_endTime_RQI, :);
    
    rspSignals_resamp_filt_pk{1} = ...
        rspSignals_resamp_filt_pk{1}(...
        rspSignals_resamp_filt_pk{1}(:, 1) >= rsp_startTime_pk & ...
        rspSignals_resamp_filt_pk{1}(:, 1) <= rsp_endTime_pk, :);
    
    
    % --------- Window the signals for peak detection ------------ %
    % For peak detection, we found that about a minute of data in each window
    % serves as a good tradeoff between making the window small enough to not
    % be affected by drastic amplitude changes, as well as large enough to not
    % factor in spikes too largely in setting thresholds. A window overlap is
    % used to ensure the boundary breaths are still detected properly.
%     windowTime_rsp_peak = 60;
%     windowShift_rsp_peak = 58;
    rspSignals_peak_wind = cell(1, numel(rspSignals));
    
    for j = 1:numel(rspSignals)
        rspSignals_peak_wind{j} = windowing(windowTime_rsp_peak, ...
            windowShift_rsp_peak, rspSignals_resamp_filt_pk{j});
    end
    
    
    % --------- Window the signals for RQI computation ------------ %
    % For RQI computation, we will first experiment with 32 s, as that is what
    % is used in literature; if that does not work, we will reduce to 16 s. No
    % overlap is necessary
%     windowTime_rsp_RQI = 16;
%     windowShift_rsp_RQI = 16;
    rspSignals_RQI_wind = cell(1, numel(rspSignals));
    
    for j = 1:numel(rspSignals)
        rspSignals_RQI_wind{j} = windowing(windowTime_rsp_RQI, ...
            windowShift_rsp_RQI, rspSignals_resamp_filt_RQI{j});
        
        % Normalize each window to zero mean, stdev = 1
        for k = 1:size(rspSignals_RQI_wind{j}, 1)
            rspSignals_RQI_wind{j}{k, 3} = ...
                normalize(rspSignals_RQI_wind{j}{k, 3}, 2);
        end
    end

    
    
    % ----------- Detect peak for each rsp signal ----------- %
    % Using the windowed cell arrays of respiration data, we will extract peaks
    % from each of the four respiration signals so we can later
    % choose which of the four sets of peaks to use in each RQI window
    % To extract the peaks, we go with the findpeaks() method
    % proposed by Vera Zamoscik et al. 2018 in their paper on respiration
    % pattern variability in depression, though the minimum prominence may
    % change depending on how I feel performance is; also note that the paper
    % had a typo --> 'variance' should be replaced by 'standard deviation.'
    % Another method that works, but just not as well with this dataset is
    % the count-adv algorithm. The option is available to use this algorithm
    % instead
    rsp_peakTimes = cell(1, numel(rspSignals));
    
    for j = 1:numel(rspSignals)
        [~, rsp_peakTimes{j}] = ...
            breathDetect_RSP(rspSignals_peak_wind{j}, ...
            rspSignals_resamp_filt_pk{j}, newFs_rsp_pk, PeakDetectMethod);
    end
    
    % Use custom function to evaluate each rsp signal at the corresponding
    % peak times, mainly for debugging purposes (e.g., by plotting)
    % Below function takes as input search times, original signal times,
    % original signal values, and sampling frequency; outputs values
    % corresponding to search times
    rsp_peakValues = cell(1, numel(rspSignals));
    
    for j = 1:numel(rspSignals)
        [~, rsp_peakValues{j}] = ...
            getSampleIndices(rsp_peakTimes{j}, ...
            rspSignals_resamp_filt_pk{j}(:, 2), ...
            rspSignals_resamp_filt_pk{j}(:, 1), newFs_rsp_pk);
    end
    
    
    % ---------- Compute inter-breath intervals ----------- %
    % Using the peaks extracted thus far, I can compute IBIs for each of the 4
    % signals, where IBI stands for inter-breath intervals, i.e., the time
    % between two breath peaks. This function finds plausible IBIs by removing
    % any IBIs > 10 s or < 1.4 s
    interBreaths = cell(1, numel(rspSignals));
    
    for j = 1:numel(rspSignals)
        interBreaths{j} = IBI_findplaus(rsp_peakTimes{j},max_IBI, min_IBI);
    end
    
    
    
    % --------- Compute associated onset times ------------ %
    % Later when we calculate inspiration time and expiration time, we will
    % need these onset times handy
    rsp_onsetTimes = cell(1, numel(rspSignals));
    
    for j = 1:numel(rspSignals)
        rsp_onsetTimes{j} = rspOnset_find(interBreaths{j}, ...
            rspSignals_resamp_filt_pk{j});
    end
    
    
    % -------- Compute RQI values for RQI windows created earlier -------- %
    % For each of the 4 signals, we will compute an RQI value for each of the
    % windows created for RQI purposes. Since this is done on a
    % window-by-window basis, the returned result will have window start,
    % window end, and the associated RQI value stored (Nx3 array, where N is
    % the number of windows)
    rsp_RQI = cell(1, numel(rspSignals));
    
    for j = 1:numel(rspSignals)
        rsp_RQI{j} = RQI_windowed(rspSignals_RQI_wind{j}, ...
            newFs_rsp_RQI, windowTime_rsp_RQI);
    end
    
    
    
    % -------- Fuse IBIs and onsets based on RQI and bad ECG times ----- %
    % The idea here is that we want to choose the inter-breath intervals (IBIs)
    % corresponding to periods where at least one of the signals was of
    % sufficient quality. Moreover, we want to choose the IBIs of the signal
    % with the highest quality, as determined by the RQI scores
    % Note that the result from this step will be a single 2 column array 
    % because we are fusing the sets of IBIs into 1 set.
    
    % Onset times will just be a vector with elements corresponding to IBIs
    % (i.e., column vector with same number of rows as fused_IBIs)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Specify the minimum RQI required to accept a window of data
    % (MAY NEED MANUAL TUNING FOR THE PARTICULAR DATASET DEPENDING ON NOISE
    % CHARACTERISTICS AND SENSITIVITY/SPECIFICITY DESIRED)
%     RQI_thresh = 0.45;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    numSignals = numel(rspSignals);
    
    % For the fusion function, we must specify which signals are
    % ecg-derived, if any
    % Specifies which signals are ECG-derived
    ecgDerived_indices = [2, 3];
    
    [fused_IBIs, fused_onsetTimes] = ...
        fuseRSP_IBIs_onsets(numSignals, RQI_thresh, 'RQI', ...
        rsp_RQI{1}, rsp_RQI{2}, rsp_RQI{3}, ...
        'IBI', interBreaths{1}, interBreaths{2}, interBreaths{3}, ...
        'onset', rsp_onsetTimes{1}, rsp_onsetTimes{2}, ...
        rsp_onsetTimes{3}, ...
        'ecg-derived', badTimes_ecg, ecgDerived_indices);
    
    
    % -------- Remove IBIs and onsets with outlying RR, Ti, or Te ------- %
    % Before we extract the final respiration rates (RR), inspiration times
    % (Ti), and expiration times (Te), we need to first remove all IBIs
    % corresponding to outlying RR, Ti, or Te
    % Flag outlier intervals and onsets based on respiration rate
    [~, flagged_respRate] = extractRespRate(fused_IBIs, ...
        'flagOutlier');
    
    % Flag outlier intervals and onsets based on inspiration time
    [~, flagged_Ti] = extractTi(fused_IBIs, ...
        fused_onsetTimes, 'flagOutlier');
    
    % Flag outlier intervals and onsets based on expiration time
    [~, flagged_Te] = extractTe(fused_IBIs, ...
        fused_onsetTimes, 'flagOutlier');
    
    % Find the union of the flagged beats for both features
    flagged_unique_rsp = ...
        unique([flagged_respRate; flagged_Ti; flagged_Te]);
    
    % Remove these intervals and onsets from fused_IBIs and rsp_onsets
    for j = size(fused_IBIs, 1):-1:1
        % If this row index is a member of the flagged ones
        if ismember(j, flagged_unique_rsp)
            % Remove IBI and onset time
            fused_IBIs(j, :) = [];
            fused_onsetTimes(j) = [];
        end
    end
    
    % -------- Extract RR, Ti, Te, and the ratio Ti/Te ----------- %
    % With the final set of IBIs and onset times determined, we can now compute
    % RR by taking 60/IBI, inspiration time (Ti) by subtracting the beginning
    % of the IBI interval from the onset time, and Te by subtracting the onset
    % time from the end of the IBI interval
    % Extract respiration rate
    [physFeatures_resp.RR, ~] = extractRespRate(fused_IBIs);
    
    % Extract inspiration time
    [physFeatures_resp.Ti, ~] = extractTi(fused_IBIs, ...
        fused_onsetTimes);
    
    % Extract expiration time
    [physFeatures_resp.Te, ~] = extractTe(fused_IBIs, ...
        fused_onsetTimes);
    
    % Store the Ti/Te ratio
    % (recall that first column is time)
    physFeatures_resp.Ti_over_Te = ...
        [physFeatures_resp.Ti(:, 1), ...
        physFeatures_resp.Ti(:, 2) ./ ...
        physFeatures_resp.Te(:, 2)];
    
    
    
    % ----------- Respiration Pattern Variability (RPV) Features --------- %
    % We will extract three variability features from each of the three
    % original time series features extracted above (RR, Ti, and Te). These
    % three variability features are coefficient of variation (CV),
    % autocorrelation at one breath lag (AR(1)), and root mean square of
    % successive differences (RMSSD). These RPV features will be
    % extracted using a rolling window, similar to the case for the HRV
    % features in the ECG section of the code.
    
%     windowTime_RPV = 300;       % 300 seconds
%     windowShift_RPV = 1;        % 299 s overlap
    
    % For the below functions, we need the mean inter-breath interval to
    % understand whether enough data exist to extract a RPV feature for that
    % window (this mean interbreath interval is used to be fair across
    % subjects - some subjects may just breathe slower, so the functions
    % below consider that when deciding whether enough data exist or not)
    mean_IBI = mean(fused_IBIs(:, 2));
    
    % Compute RPV features
    % Create window starts and ends based on fused_IBIs
    RPV_windows_cell = windowing(windowTime_RPV, ...
        windowShift_RPV, fused_IBIs);
    
    % Because a cell array is returned, let's get the first two columns
    % and store them in a separate array
    RPV_windows = zeros(size(RPV_windows_cell, 1), 2);
    for j = 1:size(RPV_windows_cell, 1)
        RPV_windows(j, :) = ...
            [RPV_windows_cell{j, 1}, RPV_windows_cell{j, 2}];
    end
    
    % Store features one by one
    % Also store number of windows removed
    % RR
    [physFeatures_RPV.RR_CV, RPV_wind_rmvd_RR] = ...
        coeff_var(physFeatures_resp.RR, RPV_windows, mean_IBI);
    
    [physFeatures_RPV.RR_AR_1, ~] = ...
        AR_1Breath(physFeatures_resp.RR, RPV_windows, mean_IBI, false);
    
    [physFeatures_RPV.RR_AR_1_norm, ~] = ...
        AR_1Breath(physFeatures_resp.RR, RPV_windows, mean_IBI, true);
    
    [physFeatures_RPV.RR_RMSSD, ~] = ...
        RMS_succ_diff(physFeatures_resp.RR, RPV_windows, mean_IBI);
    
    % Ti
    [physFeatures_RPV.Ti_CV, RPV_wind_rmvd_Ti] = ...
        coeff_var(physFeatures_resp.Ti, RPV_windows, mean_IBI);
    
    [physFeatures_RPV.Ti_AR_1, ~] = ...
        AR_1Breath(physFeatures_resp.Ti, RPV_windows, mean_IBI, false);
    
    [physFeatures_RPV.Ti_AR_1_norm, ~] = ...
        AR_1Breath(physFeatures_resp.Ti, RPV_windows, mean_IBI, true);
    
    [physFeatures_RPV.Ti_RMSSD, ~] = ...
        RMS_succ_diff(physFeatures_resp.Ti, RPV_windows, mean_IBI);
    
    % Te
    [physFeatures_RPV.Te_CV, RPV_wind_rmvd_Te] = ...
        coeff_var(physFeatures_resp.Te, RPV_windows, mean_IBI);
    
    [physFeatures_RPV.Te_AR_1, ~] = ...
        AR_1Breath(physFeatures_resp.Te, RPV_windows, mean_IBI, false);
    
    [physFeatures_RPV.Te_AR_1_norm, ~] = ...
        AR_1Breath(physFeatures_resp.Te, RPV_windows, mean_IBI, true);
    
    [physFeatures_RPV.Te_RMSSD, ~] = ...
        RMS_succ_diff(physFeatures_resp.Te, RPV_windows, mean_IBI);
    
    
    % Just in case I want to know how much respiration data I used/removed
    cleanTimes_rsp = [];
    for j = 1:size(fused_IBIs, 1)
        % Concatenate as a single column vector
        % Store all times from the first breath time to the
        % second breath time for each inter-breath interval
        % That's why we increment by the sample time
        % Sample time = (1/newFs_rsp_pk(i))
        cleanTimes_rsp = [cleanTimes_rsp; ...
            ((fused_IBIs(j, 1) - fused_IBIs(j, 2)):...
            1/newFs_rsp_pk:fused_IBIs(j, 1))'];
        
        % Eliminate repeated times
        cleanTimes_rsp = unique(cleanTimes_rsp);
    end
    
    % Compute the percent of samples removed from respiration
    % We'll just use the first respiration signal as our reference
    % (shouldn't matter since we resized them to be same length)
    percentRmvd_rsp = 100*(1 - ...
        (length(cleanTimes_rsp)/...
        size(rspSignals_resamp_filt_pk{1}, 1)));
    
    
    
    %% Save results
    % With signal processing and feature extraction complete, now it's time
    % to save everything for accessing later :)
    
    % Let's set the file name
    % Path to where processed output will be stored
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pathSaveDir = ['Processed_Features_ECG_RSP', filesep];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    saveFileName = [dataFileName(1:end-4), '_procFeat_ECG_RSP', '.mat'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Now, because we may use a for loop, we ave to implement a
    % workaround to save our workspace to the above filename
    % We will leverage my parsave function, which takes as input a filename
    % and a struct with fields populated by variables we want to save and
    % will save the file with all the variables separately labeled (i.e.,
    % won't all be stuck in a struct, so it'll be as if we were looking at
    % the workspace as we could be now)
    % To use this function, we need to save all desired variables into
    % struct fields
    save_struct = struct();
    
    % Features
    save_struct.physFeatures_beat = physFeatures_beat;
    save_struct.physFeatures_HRV = physFeatures_HRV;
    save_struct.physFeatures_resp = physFeatures_resp;
    save_struct.physFeatures_RPV = physFeatures_RPV;
%     save_struct.subjectID = subjectID;
    
    % Extra Info
    save_struct.cleanTimes_rsp = cleanTimes_rsp;
    save_struct.Fs = Fs;
    save_struct.flagged_respRate = flagged_respRate;
    save_struct.flagged_Te = flagged_Te;
    save_struct.flagged_Ti = flagged_Ti;
    save_struct.flagged_unique_rsp = flagged_unique_rsp;
    save_struct.fused_IBIs = fused_IBIs;
    save_struct.fused_onsetTimes = fused_onsetTimes;
    save_struct.HRVparams = HRVparams;
    save_struct.interBreaths = interBreaths;
    save_struct.mean_IBI = mean_IBI;
    save_struct.newFs_rsp_pk = newFs_rsp_pk;
    save_struct.newFs_rsp_RQI = newFs_rsp_RQI;
    save_struct.percentRmvd_rsp = percentRmvd_rsp;
    save_struct.RPV_wind_rmvd_RR = RPV_wind_rmvd_RR;
    save_struct.RPV_wind_rmvd_Te = RPV_wind_rmvd_Te;
    save_struct.RPV_wind_rmvd_Ti = RPV_wind_rmvd_Ti;
    save_struct.RPV_windows = RPV_windows;
    save_struct.rsp_onsetTimes = rsp_onsetTimes;
    save_struct.rsp_peakTimes = rsp_peakTimes;
    save_struct.rsp_peakValues = rsp_peakValues;
    save_struct.rsp_RQI = rsp_RQI;
    save_struct.rspSignals = rspSignals;
    save_struct.rspSignals_peak_wind = rspSignals_peak_wind;
    save_struct.rspSignals_resamp_filt_pk = rspSignals_resamp_filt_pk;
    save_struct.rspSignals_resamp_filt_RQI = rspSignals_resamp_filt_RQI;
    save_struct.rspSignals_RQI_wind = rspSignals_RQI_wind;
    save_struct.windowShift_RPV = windowShift_RPV;
    save_struct.windowShift_rsp_peak = windowShift_rsp_peak;
    save_struct.windowShift_rsp_RQI = windowShift_rsp_RQI;
    save_struct.windowTime_RPV = windowTime_RPV;
    save_struct.windowTime_rsp_peak = windowTime_rsp_peak;
    save_struct.windowTime_rsp_RQI = windowTime_rsp_RQI;
    save_struct.badTimes_ecg = badTimes_ecg;
    save_struct.clean_Rpeak_times = clean_Rpeak_times;
    save_struct.clean_Rpeak_values = clean_Rpeak_values;
    save_struct.cleanTimes_ecg = cleanTimes_ecg;
    save_struct.ecgDerived_indices = ecgDerived_indices;
    save_struct.ecgDeriveRsp_RIAV = ecgDeriveRsp_RIAV;
    save_struct.ecgDeriveRsp_RIFV = ecgDeriveRsp_RIFV;
    save_struct.ecgDeriveRsp_RIIV = ecgDeriveRsp_RIIV;
    save_struct.HRVout = HRVout;
    save_struct.percentRmvd_ecg = percentRmvd_ecg;
    save_struct.Rpeak_good_v_bad = Rpeak_good_v_bad;
    save_struct.Rpeak_times = Rpeak_times;
    save_struct.Rpeak_values = Rpeak_values;
    save_struct.RRintervals = RRintervals;
    save_struct.rsp_endTime_pk = rsp_endTime_pk;
    save_struct.rsp_endTime_RQI = rsp_endTime_RQI;
    save_struct.rsp_startTime_pk = rsp_startTime_pk;
    save_struct.rsp_startTime_RQI = rsp_startTime_RQI;
    
    % Finally, save these variables :)
  parsave([dataFileLocation, filesep, saveFileName], save_struct);
  disp(['Processing Complete. Processed data located at ' dataFileLocation, filesep, saveFileName])
    
    
    %% Store data removal statistics and such in summary arrays
    % These are updated each loop iteration so we produce a summary array
    % by the time we finish that we can store away for safe keeping so we
    % have data removal statistics handy
    
    % Percents removed
    percentRmvd_array_ecg = percentRmvd_ecg;
    
    % RSP
    percentRmvd_array_rsp = percentRmvd_rsp;
    
    % RPV windows removed
    RPVwindowsRmvd_array_RR = RPV_wind_rmvd_RR;
    RPVwindowsRmvd_array_Ti = RPV_wind_rmvd_Ti;
    RPVwindowsRmvd_array_Te = RPV_wind_rmvd_Te;
    
    



% %% Save summary info
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pathSaveDir = ['Processed_Features_ECG_RSP', filesep];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveName = [dataFileLocation,pathSaveDir, 'Summary_', date, '.mat'];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% save(saveName, 'percentRmvd_array_ecg', 'percentRmvd_array_rsp', ...
%     'RPVwindowsRmvd_array_RR', 'RPVwindowsRmvd_array_Te', ...
%     'RPVwindowsRmvd_array_Ti');

end 

% ------------------------------------------------
% Functions
% ------------------------------------------------

% Function to concatenate time to keep track for protocol info
function varargout = concat_time(Fs, varargin)

% Create Biopac time vector based on sampling frequency (Fs)
Biopac_time = [0:1/Fs:(length(varargin{1})-1)/Fs]';

% Concatenate and store
for indx = 1:length(varargin)
    varargout{indx} = [Biopac_time, varargin{indx}];
end

end


