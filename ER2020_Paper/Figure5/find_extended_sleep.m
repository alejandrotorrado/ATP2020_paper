function [ext_sleep_idx] = find_extended_sleep(st,dur_thresh,short_thresh)
%
% find_extended_sleep
%
% Alejandro Torrado Pacheco - 2018
%
% This function goes through a set of statetimes and finds extended sleep
% periods based on the parameters passed in as arguments.
%
% INPUTS
%   _st: Nx2 array with state codes in col 1 and state times in col 2
%   _dur_thresh: minimum duration of extended sleep periods, in minutes
%   _short_thresh: optional; maximum interruption to consider a state as
%                  extended sleep. This is set to 60 sec by default. This
%                  means that for dur_thresh=30 for example, this code will
%                  find periods of continuous sleep without any intervening
%                  wake epochs longer than 60 sec.
%
% OUTPUTS
%   _ext_sleep_idx: Nx2 array containing indices of extended sleep periods
%                   ([start, end]). These refer to state times.

% deal with optional argument
switch nargin
    case 1
        dur_thresh = 30; % minutes
        short_thresh = 60; % seconds
    case 2
        short_thresh = 60; % seconds
end

% convert to seconds
dur_thresh_sec = dur_thresh * 60;

% get state codes and durations
codes = st(:,1);
durations = diff(st(:,2));

% convert from 4-state coding to 2-state coding (only S and W)
new_codes = zeros(size(codes));
new_codes(codes<3) = 100; % 100 = SLEEP
new_codes(codes>3) = 400; % 400 = WAKE

% make short wake periods ignorable
for cc = 2:size(new_codes,1)-1
    % if state is WAKE
    if new_codes(cc) == 400
        % and both before and after are sleep
        if new_codes(cc-1) == 100 && new_codes(cc+1) == 100
            % if the WAKE is shorter than threshold
            if durations(cc) <= short_thresh
                % do not count it as wake (give it a code of 250, in between wake and sleep)
                new_codes(cc) = 250;
            end
        end
    end
end

% find all the sleep states
sleeps = find(new_codes == 100);

% subtract state codes. This will be 0 when the state doesn't change (e.g.
% going from REM to NREM), -300 when going from W to S, and +300 when going
% from S to W
diff_new = diff(new_codes);

% based on the diff value, find when epochs of sleep begin and end.
% Interruptions shorter than the short_thresh value are ignored because
% their diff values are -150 or 150
sleep_starts = find(diff_new == -300)+1;
sleep_ends = find(diff_new == 300);

% deal with edge cases
if sleep_ends(1) < sleep_starts(1)
    sleep_starts = [1; sleep_starts];
end

if sleep_starts(end) > sleep_ends(end)
    sleep_ends = [sleep_ends; numel(new_codes)];
end

% find indices of extended sleep periods
tmp_idx = [sleep_starts, sleep_ends];

% deal with edge cases
if tmp_idx(end,2) == size(st,1)
    tmp_idx(end,2) = tmp_idx(end,2)-1;
end

if tmp_idx(end,1) == tmp_idx(end,2)
    tmp_idx(end,:) = [];
end

% now check for sleep state duration
good_idx = ones(size(tmp_idx,1),1);

% loop through indices
for xx = 1:size(tmp_idx,1)-1
    % calculate duration
    try
        duration = st(tmp_idx(xx,2)+1,2) - st(tmp_idx(xx,1),2);
    catch
        keyboard;
    end
    % if above threshold, add to array of good indices
    if duration >= dur_thresh_sec
        good_idx(xx,1) = 1;
    else
        good_idx(xx,1) = 0;
    end
end

% select the good indices
tmp_good = tmp_idx(logical(good_idx),:);

% difference between columns - this tells us about how many intermediate
% states there are within the extended sleep period
diff_idx = diff(tmp_good,[],2);

% select only extendedl sleep periods within the good ones that have at
% least 3 different states (this avoids problems with later analysis that
% would arise from ext sleep states that are just one NREM followed by one
% REM)
ext_sleep_idx = tmp_good(diff_idx>=3,:);





