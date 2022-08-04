# RespProcessing_RPV

UNDER CONSTRUCTION

Code associated with the following journal article on respiration signal processing and robust estimation of respiration pattern variability (RPV):

A. H. Gazi, M. T. Wittbrodt, A. B. Harrison, S. Sundararaj, N. Z. Gurel, J. A. Nye, A. J. Shah, V. Vaccarino, J. D. Bremner, and O. T. Inan, 
"Robust Estimation of Respiratory Variability Uncovers Correlates of Limbic Brain Activity and Transcutaneous Cervical Vagus Nerve Stimulation in the Context of Traumatic Stress", 
IEEE Transactions on Biomedical Engineering, in press, 2021.

If ECG processing and heart rate variability (HRV) extraction tools are used, please also cite:

A. Vest, G. Da Poian, Q. Li, C. Liu, S. Nemati, A. Shah, and G. D. Clifford, 
"An Open Source Benchmarked Toolbox for Cardiovascular Waveform and Interval Analysis", 
Physiological measurement, vol. 39, no. 10, 2018. 

Please see the original PhysioNet Cardiovascular Signal Toolbox repository for questions specifically related to ECG processing and HRV extraction: 
https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox

Loaded Variables:
(option given to allow for only one or the other, instead of both)
ecg - N x 1 array of ECG datapoints
rsp - N x 1 array of RSP datapoints
(one of the below is needed)
isi - intersample interval in milliseconds
Fs - sampling frequency in Hz
subjectID - unique identifier for this data file to keep separate from others
(below is optional)
atrialFib - true/false whether subject has known atrial fibrillation


Saved Variables:
 ----- Main Output ------
physFeatures_beat - Struct that stores beat-by-beat feature arrays in each field, where each array is 2-D with times as first column and feature values in second
physFeatures_HRV - Struct that stores windowed heart rate variability arrays in each field, where each array is 3-D with window start time as first column, window end time as second column, and feature value as third
physFeatures_resp - Struct that stores breath-by-breath feature arrays in each field, where each array is 2-D with times as first column and feature values in second
physFeatures_RPV - Struct that stores windowed respiration pattern variability arrays in each field, where each array is 3-D with window start time as first column, window end time as second column, and feature value as third
subjectID - unique identifier given to this dataset (usually corresponding to a recording for a subject)


----- Summary Arrays -----
percentRmvd_array_ecg - 1-D array storing the percent removed of ECG data from the recordings provided for each data file; ordered corresponding to ordering of data files / subject IDs
percentRmvd_array_rsp - 1-D array storing the percent removed of RSP data from the recordings provided for each data file; ordered corresponding to ordering of data files / subject IDs
RPV_wind_rmvd_RR - 1-D array storing the number of respiration rate-based RPV windows removed from the recordings provided for each data file; ordered corresponding to ordering of data files / subject IDs
RPV_wind_rmvd_Ti - 1-D array storing the number of inspiration time-based RPV windows removed from the recordings provided for each data file; ordered corresponding to ordering of data files / subject IDs
RPV_wind_rmvd_Te - 1-D array storing the number of expiration time-based RPV windows removed from the recordings provided for each data file; ordered corresponding to ordering of data files / subject IDs


------ Extra Info --------
cleanTimes_rsp - 1-D array storing approximated times corresponding to clean respiration data
Fs - sampling frequency of data in Hz
flagged_respRate - elements of fused IBIs and onsets that were flagged for outlier respiration rates
flagged_Te - elements of fused IBIs and onsets that were flagged for outlier expiration times
flagged_Ti - elements of fused IBIs and onsets that were flagged for outlier inspiration times
flagged_unique_rsp - union of flagged_respRate, flagged_Te, and flagged_Ti
fused_IBIs - 2-D array storing second breath timings and corresponding interbreath intervals after fusing
fused_onsetTimes - 1-D array storing onset timings corresponding to fused IBIs
HRVparams - PhysioNet Cardiovascular Signal Toolbox parameter struct
interBreaths - Cell array storing IBI arrays prior to fusion
mean_IBI - mean interbreath interval for this data file
newFs_rsp_pk - new sampling frequency (in Hz) used for peak detection
newFs_rsp_RQI - new sampling frequency (in Hz) used for respiration quality indexing
percentRmvd_rsp - approximate percentage of respiration data removed from the entire recording
RPV_wind_rmvd_RR - RPV windows removed for the respiration rate-based RPV features
RPV_wind_rmvd_Te - RPV windows removed for the expiration time-based RPV features
RPV_wind_rmvd_Ti - RPV windows removed for the inspiration time-based RPV features
RPV_windows - Stacked array of windows used to compute RPV values, where each window is a row
rsp_onsetTimes - Cell array of respiration onset times corresponding to each signal
rsp_peakTimes - Cell array of respiration peak times corresponding to each signal
rsp_peakValues - Cell array of respiration peak values corresponding to each signal
rsp_RQI - Cell array of respiration quality index values for each signal, for each window
rspSignals - Cell array that stores each of the respiration signals as its own 2-D array with time as first column and data as second
rspSignals_peak_wind - Cell array storing windowed respiration signals ready for peak detection (by windowing, we effectively accomplish adaptive thresholding)
rspSignals_RQI_wind - Cell array storing windowed respiration signals ready for respiration quality assessment
windowShift_RPV - Number of seconds to shift the rolling window by for computing RPV
windowShift_rsp_peak - Number of seconds to shift window by for windowing for peak detection
windowShift_rsp_RQI - Number of seconds to shift window by for windowing for RQI
windowTime_RPV - Window length, in seconds, used to compute RPV metrics
windowTime_rsp_peak - Number of seconds used for window length for peak detection windowing
windowTime_rsp_RQI - Number of seconds used for window length for RQI windowing
badTimes_ecg - 1-D array of approximate times corresponding to bad ECG data
clean_Rpeak_times - 1-D array of times corresponding to R peaks associated with "clean" RR intervals
clean_Rpeak_values - 1-D array of values corresponding to R peaks associated with "clean" RR intervals
cleanTimes_ecg - 1-D array of approximate times corresponding to clean ECG data
ecg
ecgDerived_indices - Array specifying the indices of the rspSignals cell array that correspond to ECG-derived signals; empty if no ECG-derived signals are used
ecgDeriveRsp_RIAV - 2-D array with time as first column and data as second column for ECG-derived respiration extracted using respiration-induced amplitude variations
ecgDeriveRsp_RIFV - 2-D array with time as first column and data as second column for ECG-derived respiration extracted using respiration-induced frequency variations
ecgDeriveRsp_RIIV - 2-D array with time as first column and data as second column for ECG-derived respiration extracted using respiration-induced intensity variations
HRVout - Array that stores output of PhysioNet Cardiovascular Signal Toolbox HRV calculations
percentRmvd_ecg - Approximate percent removed of ECG data from the entirety of the recording
Rpeak_good_v_bad - Array that stored 1s and 0s for good and bad R peaks, respectively, when performing signal quality assessment; corresponds to Rpeak_times (i.e., before cleaning)
Rpeak_times - Array that stores R peak times detected prior to "cleaning"
Rpeak_values - Array that stores R peak values detected prior to "cleaning"
RRintervals - 2-D array with first column storing time of second R peak of each RR interval, and second column storing the length of the clean RR intervals
rsp_endTime_pk - end time used for respiration peak detection (used for chopping and aligning all signals)
rsp_endTime_RQI - end time used for RQI (used for chopping and aligning all signals)
rsp_startTime_pk - start time used for respiration peak detection (used for chopping and aligning all signals)
rsp_startTime_RQI - start time used for RQI (used for chopping and aligning all signals)
