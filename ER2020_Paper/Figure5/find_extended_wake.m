function [ext_wake_idx] = find_extended_wake(st,dur_thresh,short_thresh)
%
% find_extended_wake
%
% Alejandro Torrado Pacheco - 2018
%
% See also: find_extended_sleep.m
%
% This function goes through a set of statetimes and finds extended wake
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
%   _ext_wake_idx: Nx2 array containing indices of extended sleep periods
%                   ([start, end]). These refer to state times.

switch nargin
    case 1
        dur_thresh = 30; % minutes
        short_thresh = 60; % seconds
    case 2
        short_thresh = 60; % seconds
end

dur_thresh_sec = dur_thresh * 60;

codes = st(:,1);
durations = diff(st(:,2));

new_codes = zeros(size(codes));
new_codes(codes<3) = 100; % 100 = SLEEP
new_codes(codes>3) = 400; % 400 = WAKE

for cc = 2:size(new_codes,1)-1
    % if state is SLEEP
    if new_codes(cc) == 100
        % and both before and after are wake
        if new_codes(cc-1) == 400 && new_codes(cc+1) == 400
            % if the SLEEP is shorter than threshold
            if durations(cc) <= short_thresh
                % do not count it as sleep (give it a code of 250, in between wake and sleep)
                new_codes(cc) = 250;
            end
        end
    end
end

wakes = find(new_codes == 400);

diff_new = diff(new_codes);

% based on the diff value, find when epochs of sleep begin and end.
% Interruptions shorter than the short_thresh value are ignored because
% their diff values are -150 or 150
wake_starts = find(diff_new == 300)+1;
wake_ends = find(diff_new == -300);

if wake_ends(1) < wake_starts(1)
    wake_starts = [1; wake_starts];
end

if wake_starts(end) > wake_ends(end)
    wake_ends = [wake_ends; numel(new_codes)];
end

tmp_idx = [wake_starts, wake_ends];

if tmp_idx(end,2) == size(st,1)
    tmp_idx(end,2) = tmp_idx(end,2)-1;
end

if tmp_idx(end,1) == tmp_idx(end,2)
    tmp_idx(end,:) = [];
end


good_idx = ones(size(tmp_idx,1),1);

for xx = 1:size(tmp_idx,1)-1
    try
        duration = st(tmp_idx(xx,2)+1,2) - st(tmp_idx(xx,1),2);
    catch
        keyboard;
    end
    if duration >= dur_thresh_sec
        good_idx(xx,1) = 1;
    else
        good_idx(xx,1) = 0;
    end
end

tmp_good = tmp_idx(logical(good_idx),:);

% difference between columns - this tells us about how many intermediate
% states there are within the extended wake period
diff_idx = diff(tmp_good,[],2);

% select only extendedl wake periods within the good ones that have at
% least 3 different states (this avoids problems with later analysis that
% would arise from ext wake states that are just one QW followed by one
% AW for instance)
ext_wake_idx = tmp_good(diff_idx>=3,:);





