function [ddata] = EMGpassFilt(data, samplert,resamp_rate,lowPassVal,highPassVal, n)
% LFPdownfilt applies bandpass cheby1 filter to data and then resamples to
% a lower rate.

nyq               = samplert/2;			          % in Hz
%n                 = 2;                % order of the Chebyshev Type1 filter
R                 = 0.8;                % dB of peak-to-peak ripple in the passband
Wp                = [  highPassVal/nyq  lowPassVal/nyq    ];           % define bandpass



[b,a] = cheby1(n,R,Wp);


fdata = filtfilt(b,a,data);
%resample the data


Y = fdata( 1 : round(24414.06/resamp_rate) : end );
ddata=Y;



