%% swdense_plot_example_states_24hr.m
%
% Alejandro Torrado Pacheco - 2019
%
% Reproduce Fig 4A - example of S-W pattern for one animal.
%
% The animal can be selected at the top of the script. For the figure in
% the paper, it was AT12.
%
% To use this script, first load CONTCELL_recov.mat

%% select animal here
anim = 'AT12';
statetimes = CONTCELL_recov.STATETIMES.(anim);

%% process statetimes

interrupt_threshold = 10; % seconds - remove states shorter than this
duration_threshold = 30; % seconds

% remove repeats
st_diff = diff(statetimes(:,1));
kill_these = find(st_diff==0);
statetimes(kill_these+1,:) = [];

% remove short states
st_timediff = diff(statetimes(:,2));
too_short = find(st_timediff <= interrupt_threshold);
statetimes(too_short,:) = [];

% remove repeats again
st_diff = diff(statetimes(:,1));
kill_these = find(st_diff==0);
statetimes(kill_these+1,:) = [];

% switch from 4-state classification (REM, NREM, QW, AW) to 2-state
% classification (Sleep, Wake).
% In statetimes, the state codes are:
% 1 = REM
% 2 = NREM
% 4 = Active Wake
% 5 = Quiet Wake
% NEW CODES:
% convert all wake to 400 and sleep to 100

% state codes are in the first column
X = statetimes(:,1);
X(X>3) = 400; % do these first otherwise it will be screwed up
X(X<3) = 100;

% state start timestamps are in column 2
T = statetimes(:,2);

% recompile new statetimes
new_ST = [X, T];

% remove repeats
st_diff = diff(new_ST(:,1));
kill_these = find(st_diff==0);
new_ST(kill_these+1,:) = [];

% separate new times and codes
new_T = new_ST(:,2);
new_X = new_ST(:,1);

%% select 24-hour period to plot

% find the first cell from this animal to get the experiment start time
all_anims = {CONTCELL_recov.MASTER.animal};
first_anim_idx = find(strcmp(all_anims,anim),1,'first');
% note: expt start time is a unix timestamp. the function unixtime converts
% it to a full date (format: [year months day hour minute second]). To get
% the experiment start time aligned to ZT0 (zeitgeber time 0, i.e. lights
% ON on the expt start day), find the unix timestamp of 7:30 am on the
% first day.
expt_start_time = unixtime(CONTCELL_recov.MASTER(first_anim_idx).EXPTSTART);
expt_start_ZT0 = [expt_start_time(1:3) 7 30 0];

% offset determines the day to plot. offset of 0 would plot the first
% experiment day (bear in mind not all states might be present on the first
% day, since experiments were often started after 7:30 am on day 1)
offset = 3;
% t0 is the beginning timestamp of the day to be plotted
t0 = expt_start_ZT0 + offset*24*3600;
% t1 is just t0 + 24 hours
t1 = t0 + 24*3600;

% find the statetimes corresponding to the selected day
ix0 = find(new_T < t0, 1, 'last');
ix1 = find(new_T > t1, 1, 'first');

% set the colors for plotting
W_col = [206 107 77]./255;
S_col = [85 119 146]./255;


%% Fig 4A - plot example sleep/wake pattern over 24 hours

% set plot limits (height and baseline) and default plotting params
h = 1;
b = 0;
setFigureDefaults;

% make figure
swd_fig = figure();
set(swd_fig,'position',[.04 .1 .9 .4]);
hold on;

% loop through states
for ee = ix0:ix1-1
    
    % for each state, find the beginning and end timestamps
    x0 = new_T(ee);
    x1 = new_T(ee+1);
    % get the rectangle width
    w = x1 - x0;
    
    % figure out if state is S or W and assign color
    s = new_X(ee);
    if s < 200
        this_col = S_col;
    elseif s > 200
        this_col = W_col;
    end
    
    % draw the rectangle
    rectangle('position',[x0, b, w, h],'facecolor',this_col);
    
end

% format axes
set(gca,'xlim',[t0 t1],'xtick',t0:3*3600:t1,'xticklabel',0:3:24,'tickdir','out');
set(gca,'ylim',[0 3]);






