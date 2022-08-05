# RespProcessing_RPV

## References
Code associated with the following journal article on respiration signal processing and robust estimation of respiration pattern variability (RPV): <br/>
A. H. Gazi et al., "Robust Estimation of Respiratory Variability Uncovers Correlates of Limbic Brain Activity and Transcutaneous Cervical Vagus Nerve Stimulation in the Context of Traumatic Stress," in IEEE Transactions on Biomedical Engineering, vol. 69, no. 2, pp. 849-859, Feb. 2022, doi: 10.1109/TBME.2021.3108135.

If ECG processing and heart rate variability (HRV) extraction tools are used, please also cite: <br/>
A. Vest, et al., "An Open Source Benchmarked Toolbox for Cardiovascular Waveform and Interval Analysis", Physiological measurement, vol. 39, no. 10, 2018. <br/>
Please see the original PhysioNet Cardiovascular Signal Toolbox (PCST) repository for questions specifically related to ECG processing and HRV extraction: 
https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox


## Preparing your data

When formatting your .mat file for processing, please ensure the following variables are included: <br/>
ecg - N x 1 array of ECG datapoints<br/>
rsp - N x 1 array of RSP datapoints<br/>
subjectID - unique identifier for this data file to keep separate from others<br/>

You will also need either an intersample interval in milliseconds or sampling frequency in Hz: <br/>
isi - intersample interval in milliseconds<br/>
Fs - sampling frequency in Hz<br/>

(below is optional)<br/>
atrialFib - true/false whether subject has known atrial fibrillation<br/>

## Using the tool
Dependencies:

Signal Processing Toolbox Deep Learning Toolbox Statistics and Machine Learning Toolbox

Also make sure all folders and subfolders are added to path (i.e. Tools, Necessary Repositories, HRV_Output and folder path to files to process)

To Run: 

Double click RespGui.mlapp RespGui.mlapp calls the gui_process.m function, so if you want to understand the code better, inspect the gui_process function.

All constituent functions are organized under the folders 'Processing_Functions' or 'Supporting_Functions'. See function code for methodological details.


## Description of Saved Variables
### Main Output
physFeatures_beat - Struct that stores beat-by-beat feature arrays in each field, where each array is 2-D with times as first column and feature values in second<br/>
physFeatures_HRV - Struct that stores windowed heart rate variability arrays in each field, where each array is 3-D with window start time as first column, window end time as second column, and feature value as third<br/>
physFeatures_resp - Struct that stores breath-by-breath feature arrays in each field, where each array is 2-D with times as first column and feature values in second<br/>
physFeatures_RPV - Struct that stores windowed respiration pattern variability arrays in each field, where each array is 3-D with window start time as first column, window end time as second column, and feature value as third<br/>
subjectID - unique identifier given to this dataset (usually corresponding to a recording for a subject)<br/>


### Summary Arrays
percentRmvd_array_ecg - 1-D array storing the percent removed of ECG data from the recordings provided for each data file; ordered corresponding to ordering of data files / subject IDs<br/>
percentRmvd_array_rsp - 1-D array storing the percent removed of RSP data from the recordings provided for each data file; ordered corresponding to ordering of data files / subject IDs<br/>
RPV_wind_rmvd_RR - 1-D array storing the number of respiration rate-based RPV windows removed from the recordings provided for each data file; ordered corresponding to ordering of data files / subject IDs<br/>
RPV_wind_rmvd_Ti - 1-D array storing the number of inspiration time-based RPV windows removed from the recordings provided for each data file; ordered corresponding to ordering of data files / subject IDs<br/>
RPV_wind_rmvd_Te - 1-D array storing the number of expiration time-based RPV windows removed from the recordings provided for each data file; ordered corresponding to ordering of data files / subject IDs<br/>


### Extra Debugging Information
cleanTimes_rsp - 1-D array storing approximated times corresponding to clean respiration data<br/>
Fs - sampling frequency of data in Hz<br/>
flagged_respRate - elements of fused IBIs and onsets that were flagged for outlier respiration rates<br/>
flagged_Te - elements of fused IBIs and onsets that were flagged for outlier expiration times<br/>
flagged_Ti - elements of fused IBIs and onsets that were flagged for outlier inspiration times<br/>
flagged_unique_rsp - union of flagged_respRate, flagged_Te, and flagged_Ti<br/>
fused_IBIs - 2-D array storing second breath timings and corresponding interbreath intervals after fusing<br/>
fused_onsetTimes - 1-D array storing onset timings corresponding to fused IBIs<br/>
HRVparams - PhysioNet Cardiovascular Signal Toolbox parameter struct<br/>
interBreaths - Cell array storing IBI arrays prior to fusion<br/>
mean_IBI - mean interbreath interval for this data file<br/>
newFs_rsp_pk - new sampling frequency (in Hz) used for peak detection<br/>
newFs_rsp_RQI - new sampling frequency (in Hz) used for respiration quality indexing<br/>
percentRmvd_rsp - approximate percentage of respiration data removed from the entire recording<br/>
RPV_wind_rmvd_RR - RPV windows removed for the respiration rate-based RPV features<br/>
RPV_wind_rmvd_Te - RPV windows removed for the expiration time-based RPV features<br/>
RPV_wind_rmvd_Ti - RPV windows removed for the inspiration time-based RPV features<br/>
RPV_windows - Stacked array of windows used to compute RPV values, where each window is a row<br/>
rsp_onsetTimes - Cell array of respiration onset times corresponding to each signal<br/>
rsp_peakTimes - Cell array of respiration peak times corresponding to each signal<br/>
rsp_peakValues - Cell array of respiration peak values corresponding to each signal<br/>
rsp_RQI - Cell array of respiration quality index values for each signal, for each window<br/>
rspSignals - Cell array that stores each of the respiration signals as its own 2-D array with time as first column and data as second<br/>
rspSignals_peak_wind - Cell array storing windowed respiration signals ready for peak detection (by windowing, we effectively accomplish adaptive thresholding)<br/>
rspSignals_RQI_wind - Cell array storing windowed respiration signals ready for respiration quality assessment<br/>
windowShift_RPV - Number of seconds to shift the rolling window by for computing RPV<br/>
windowShift_rsp_peak - Number of seconds to shift window by for windowing for peak detection<br/>
windowShift_rsp_RQI - Number of seconds to shift window by for windowing for RQI<br/>
windowTime_RPV - Window length, in seconds, used to compute RPV metrics<br/>
windowTime_rsp_peak - Number of seconds used for window length for peak detection windowing<br/>
windowTime_rsp_RQI - Number of seconds used for window length for RQI windowing<br/>
badTimes_ecg - 1-D array of approximate times corresponding to bad ECG data<br/>
clean_Rpeak_times - 1-D array of times corresponding to R peaks associated with "clean" RR intervals<br/>
clean_Rpeak_values - 1-D array of values corresponding to R peaks associated with "clean" RR intervals<br/>
cleanTimes_ecg - 1-D array of approximate times corresponding to clean ECG data<br/>
ecgDerived_indices - Array specifying the indices of the rspSignals cell array that correspond to ECG-derived signals; empty if no ECG-derived signals are used<br/>
ecgDeriveRsp_RIAV - 2-D array with time as first column and data as second column for ECG-derived respiration extracted using respiration-induced amplitude variations<br/>
ecgDeriveRsp_RIFV - 2-D array with time as first column and data as second column for ECG-derived respiration extracted using respiration-induced frequency variations<br/>
ecgDeriveRsp_RIIV - 2-D array with time as first column and data as second column for ECG-derived respiration extracted using respiration-induced intensity variations<br/>
HRVout - Array that stores output of PhysioNet Cardiovascular Signal Toolbox HRV calculations<br/>
percentRmvd_ecg - Approximate percent removed of ECG data from the entirety of the recording<br/>
Rpeak_good_v_bad - Array that stored 1s and 0s for good and bad R peaks, respectively, when performing signal quality assessment; corresponds to Rpeak_times (i.e., before cleaning)<br/>
Rpeak_times - Array that stores R peak times detected prior to "cleaning"<br/>
Rpeak_values - Array that stores R peak values detected prior to "cleaning"<br/>
RRintervals - 2-D array with first column storing time of second R peak of each RR interval, and second column storing the length of the clean RR intervals<br/>
rsp_endTime_pk - end time used for respiration peak detection (used for chopping and aligning all signals)<br/>
rsp_endTime_RQI - end time used for RQI (used for chopping and aligning all signals)<br/>
rsp_startTime_pk - start time used for respiration peak detection (used for chopping and aligning all signals)<br/>
rsp_startTime_RQI - start time used for RQI (used for chopping and aligning all signals)<br/>
