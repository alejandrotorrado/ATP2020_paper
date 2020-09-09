function ts = calc_tailslope(input_wf,interp_factor,sampling_rate)

% input_wf is the cell's peak-scaled waveform
[og_samples,maxwfsz] = max(size(input_wf));

if nargin == 1
    switch og_samples
        case 33
            interp_factor = 1;
        case 59
            interp_factor = 3;
        case 91
            interp_factor = 5;
        otherwise
            error('Problem with number of waveform samples!');
    end
    sampling_rate = 24414.06; % Hz
elseif nargin == 2
    sampling_rate = 24414.06; % Hz
end


% tail slope from 0.4 ms after initial minimum to end of WF
ms_past = 0.4;
samples_past = round(ms_past*interp_factor*sampling_rate*1e-3);
[~,min_peak] = min(input_wf);

sample0 =  min_peak + samples_past;
if sample0 >= (og_samples-2)  
    sample0 = round(og_samples - 10*(interp_factor/3));
end
sample1 = og_samples;
lims = [sample0 sample1];

xlims = 1:(diff(lims)+1);
[~,maxlimsz] = max(size(xlims));
if maxlimsz ~= maxwfsz
    xlims = xlims';
end
try
slope = polyfit( xlims  , input_wf( lims(1) : lims(2) ),1);
catch
    keyboard;
    disp('slope keyboard');
end
ts = slope(1);