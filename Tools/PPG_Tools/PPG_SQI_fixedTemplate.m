% [annot sqimatrix template valid] = PPG_SQI_buf(wave,anntime,template_ahead,windowlen,Fs)
%
% PPG_SQI.m - PPG SQI based on beat template correlation.
% (as an advice, the algorithm get 30 beats at each call and run in loop)
% by Qiao Li 30 Mar 2011
%
% input:
%     wave:       PPG data;
%     anntime:    PPG annotation time (samples), read from ple annot file,
%                 But the ann-time is the OFFSET based on wave(1)
%     template:   Last PPG beat template
%     windowlen:  length of window to calculate template(default: 30s)
%     Fs       :  sampling frequency (default Fs=125Hz)
% output:
%     annot:      ppg sqi annotation
%                     E - excellent beat;
%                     A - acceptable beat;
%                     Q - unacceptable beat
%     sqimatrix:  ppg sqi matrix
%                     [N,1]: SQI based on Direct compare
%                     [N,2]: SQI based on Linear resampling
%                     [N,3]: SQI based on Dynamic time warping
%                     [N,4]: SQI based on Clipping detection
%     template:   Current PPG beat template
%     valid:      1 or greater for valid template,
%                 0 for invalid template
%
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information
%
% 03-03-2017
% Edits by Adriana Vest
% - Changed output variable annot from numeric to cell to preserve
%   characters
% - Style changes to align loops and conditional statements
%
% 12-01-2017 Modified by Giulia Da Poian: sampling frequency as input
% parameter instead of fixed fs = 125
%
% 12-19-2017 Modified by Giulia Da Poian: replaced dp_dtw with dpfast


% sqimatrix is now a 2 column matrix; first column for direct correlation
% comparison, and second column for DTW-based correlation comparison
function sqimatrix = ...
    PPG_SQI_fixedTemplate(ppg,sampleInd_ppgBeats,template)

% Initialize SQI matrix
sqimatrix = [];

% Store template in the variable they had previously assigned to template
t=template;

% Calculate the piecewise linear approx. of template for dynamic time warping
d1=t;
% Normalize
d1=(d1-min(d1))/(max(d1)-min(d1)).*100;
[y1, pla1]=PLA(d1,1,1);

% Store template length
templatelength = length(t);

% Main Loop
for j=1:length(sampleInd_ppgBeats)
    %% SQI1: Direct compare
    % Calculate correlation coefficients based on the template
    % length
    beatbegin = sampleInd_ppgBeats(j);
    beatend = sampleInd_ppgBeats(j) + templatelength - 1;
    
    % Current beat index
    currentb = j;
    
    % Compute the correlation coefficient between ppg beat and template
    cc = corrcoef(t,ppg(beatbegin:beatend));
    c1(j) = cc(1,2);
    if (c1(j)<0)
        c1(j)=0;
    end
    
    % Store result in first column of current beat's row
    sqimatrix(currentb,1)=int8(c1(j)*100);
    
    
    %% SQI2: Dynamic Time Warping
    % Calculate correlation coefficients based on the dynamic time
    % warping
    d2=ppg(beatbegin:beatend);
    
    % if beat too long, set SQI = 0;
    if (length(d2)>length(d1)*10)
        c2(j)=0;
    else
        % Normalize
        d2=(d2-min(d2))/(max(d2)-min(d2)).*100;
        % Piecewise linear approximation
        [y2 pla2]=PLA(d2,1,1);
        
        % Check physionet cardiovascular signal toolbox functions...
        % too lazy :)
        % don't know why this is not commented better...
        [w ta tb] = simmx_dtw(y1,pla1,y2,pla2);
        try % try to use the fast version if possible
            [p,q,Dm] = dpfast(w);
        catch
            [p,q,Dm] = dp_dtw(w);
        end
        [ym1, ym2, yout1] = draw_dtw(y1,pla1,p,y2,pla2,q);
        cc=corrcoef(y1,ym2);
        c2(j)=cc(1,2);
        if (c2(j)<0)
            c2(j)=0;
        end
    end
    % Store DTW-based correlation in second column of array
    sqimatrix(currentb,2)=int8(c2(j)*100); 
end

end