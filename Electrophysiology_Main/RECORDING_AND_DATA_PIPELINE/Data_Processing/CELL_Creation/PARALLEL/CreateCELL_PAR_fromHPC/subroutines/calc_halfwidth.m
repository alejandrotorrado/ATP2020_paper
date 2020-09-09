function [halfwidth] = calc_halfwidth(input_wf,interp_factor,sampling_rate)

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


swf = interp(double(input_wf),100);

% find minimum of waveform
[~,bottom] = min(swf);
idx_start = 1;
while bottom == 1
    idx_start = idx_start+1;
    [~,bottom] = min(swf(idx_start:end));
end

half_max = -0.5;
% then find closest point at height -0.5
left_halfmax = find(swf(idx_start:bottom) < half_max,1,'first');
left_halfmax = left_halfmax + idx_start;
% deal with weird noisy waveforms
if left_halfmax == 1
    half_max = (-1 + swf(1))/2; % instead of -0.5 use this value  
    left_halfmax = find(swf(1:bottom) < half_max,1,'first');
end
% difference between -0.5 and value at left_halfmax
left_dist = abs(half_max - swf(left_halfmax));
% check previous point to make sure it is not closer to -0.5
try
    left_dist_b = abs(half_max - swf(left_halfmax-1));
catch
    keyboard;
end
% if it is, swap them
if left_dist_b < left_dist
    left_halfmax = left_halfmax - 1;
end

% repeat with right side of minimum
right_halfmax = find(swf(bottom:end) > half_max,1,'first');
right_halfmax = right_halfmax + bottom - 1;
% difference between -0.5 and value at left_halfmax
right_dist = abs(half_max - swf(right_halfmax));
% check previous point to make sure it is not closer to -0.5
right_dist_b = abs(half_max - swf(right_halfmax-1));
% if it is, swap them
if right_dist_b < right_dist
    right_halfmax = right_halfmax - 1;
end

% then the halfwidth is the distance between those two points
temp_HW = right_halfmax - left_halfmax;
% find how many ms correspond to each sample point
seconds_per_sample = 1/(interp_factor*sampling_rate);
ms_per_sample = seconds_per_sample*1e3;
% adjust for interpolation at start of this function
hw_samples = temp_HW/100;

% GET THE HALFWIDTH IN MILLISECONDS
halfwidth = hw_samples * ms_per_sample;