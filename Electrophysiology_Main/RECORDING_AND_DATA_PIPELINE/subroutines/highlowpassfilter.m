function [data,resampt]=highlowpassfilter(data)

%% GET HIGHPASSED SIGNAL
%make filter (PASSBAND: 500 - 4000 Hz)

hpc               = 500;                          % high-pass in Hz
lpc               = 4000;
samprate          = 24414.1;            		          % in Hz
nyq               = samprate/2;			          % in Hz
hWn               = [hpc/nyq lpc/nyq];         
[hb,ha]           = butter(2,hWn); % butterworth highpass

% Hi-pass to get rid of DC offsets
hpdata           = filtfilt(hb,ha,data);          % bi-directional filter
data=hpdata;