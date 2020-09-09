function neg_pos_time = calc_negpostime(input_wf,interp_factor,sampling_rate)

% input_wf is the cell's peak-scaled waveform
og_samples = max(size(input_wf));

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

% find minimum peak of wavelength
[~,bottom] = min(input_wf);
% find maximum following that
[~,nextmax] = max(input_wf(bottom:end));
nextmax = nextmax + bottom - 1;

% find number of samples between that maximum and the bottom
np_samples = nextmax - bottom;
% number of seconds and milliseconds per sample
seconds_per_sample = 1/(interp_factor*sampling_rate);
ms_per_sample = seconds_per_sample*1e3;

% GET THE HALFWIDTH IN MILLISECONDS
neg_pos_time = np_samples * ms_per_sample;