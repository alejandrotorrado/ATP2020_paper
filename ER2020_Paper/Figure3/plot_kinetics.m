%% plot_kinetics
%
% Alejandro Torrado Pacheco - 2019
%
% Plot waveform kinetics for mEPSC slice recordings. This can be used to
% reproduce figure S4E, S4F, S4G (see thesis).
%
% This code compiles data for both er2 and er4 timepoints

% clera workspace
clearIDE
% clearvars -except CELLS CELL_EVENTS

% no need to change this, code processes er4 later on
timepoint = 'er2';

% load the data
if ismac
    maindir = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig3/analyzed_traces';
elseif ispc
    maindir = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig3\analyzed_traces';
end

% cell file and events file
cell_file = [maindir filesep 'CELLS_' timepoint '.mat'];
events_file = [maindir filesep 'CELL_EVENTS_TEMP_' timepoint '.mat'];


if ~exist('CELLS','var')
    fprintf('Loading cells...\n');
    tic
    loadcells = load(cell_file);
    toc
    CELLS = loadcells.CELLS;
end

if ~exist('CELL_EVENTS','var')
    fprintf('Loading events...\n');
    tic
    loadev = load(events_file);
    toc
    CELL_EVENTS = loadev.CELL_EVENTS;
end

n_cells = numel(CELLS);


%% main loop thru cells
% find control vs re-opened (dep) cells
dep = [CELLS.dep];

er_idx = find(dep == 1);
ctrl_idx = find(dep == 0);

% CONTROL
ctrl_cells = numel(ctrl_idx);
ctrl_aligned_mini = [];

% loop through control cells
for c = 1:ctrl_cells

    c_idx = ctrl_idx(c);
    
    % get aligned mini waveform
    ctrl_aligned_mini = padmat(CELL_EVENTS(c_idx).mean_aligned_mini, ctrl_aligned_mini, 2);
    % get rise and decay times
    rise_times = CELL_EVENTS(c_idx).cell_rise_times;
    decay_taus = CELL_EVENTS(c_idx).cell_tau;
    
    % get average rise and decay times for each cell
    ctrl_rise_time(c) = nanmean(rise_times);
    ctrl_tau(c) = nanmean(decay_taus);
    
end

% average across cells
mean_ctrl_rise = nanmean(ctrl_rise_time);
mean_ctrl_tau  = nanmean(ctrl_tau);
mean_ctrl_mini = nanmean(ctrl_aligned_mini,2);

% Same thing for RE-OPENED hemisphere
er_cells = numel(er_idx);
er_aligned_mini = [];

% loop through cells
for e = 1:er_cells

    e_idx = er_idx(e);
    % aligned mini WF, rise and decay times
    er_aligned_mini = padmat(CELL_EVENTS(e_idx).mean_aligned_mini, er_aligned_mini, 2);
    rise_times = CELL_EVENTS(e_idx).cell_rise_times;
    decay_taus = CELL_EVENTS(e_idx).cell_tau;
    % average rise and decay times across events for this cell
    er_rise_time(e) = nanmean(rise_times);
    er_tau(e) = nanmean(decay_taus);
    
end
% average across cells   
mean_er_rise = nanmean(er_rise_time);
mean_er_tau  = nanmean(er_tau);
mean_er_mini = nanmean(er_aligned_mini, 2);


%% AS ABOVE BUT FOR ER4 TIMEPOINT
% clear variables
clearvars -except mean_er* mean_ctrl* er_rise* er_tau* er_aligned* ctrl_aligned* ctrl_rise* ctrl_tau*

% set timepoint
timepoint = 'er4';

% load the data
if ismac
    maindir = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig3/analyzed_traces';
elseif ispc
    maindir = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig3\analyzed_traces';
end

cell_file = [maindir filesep 'CELLS_' timepoint '.mat'];
events_file = [maindir filesep 'CELL_EVENTS_' timepoint '.mat'];


if ~exist('CELLS','var')
    fprintf('Loading cells...\n');
    tic
    loadcells = load(cell_file);
    toc
    CELLS = loadcells.CELLS;
end

if ~exist('CELL_EVENTS','var')
    fprintf('Loading events...\n');
    tic
    loadev = load(events_file);
    toc
    CELL_EVENTS = loadev.CELL_EVENTS;
end

n_cells = numel(CELLS);

%% Loop through cells as above

dep = [CELLS.dep];

Ler_idx = find(dep == 1);
Lctrl_idx = find(dep == 0);

% CONTROL
Lctrl_cells = numel(Lctrl_idx);
Lctrl_aligned_mini = [];

for c = 1:Lctrl_cells

    c_idx = Lctrl_idx(c);
    
    Lctrl_aligned_mini = padmat(CELL_EVENTS(c_idx).mean_aligned_mini, Lctrl_aligned_mini, 2);
    rise_times = CELL_EVENTS(c_idx).cell_rise_times;
    decay_taus = CELL_EVENTS(c_idx).cell_tau;
    
    Lctrl_rise_time(c) = nanmean(rise_times);
    Lctrl_tau(c) = nanmean(decay_taus);
    
end

mean_Lctrl_rise = nanmean(Lctrl_rise_time);
mean_Lctrl_tau  = nanmean(Lctrl_tau);
mean_Lctrl_mini = nanmean(Lctrl_aligned_mini,2);

% RE-OPENED
Ler_cells = numel(Ler_idx);
Ler_aligned_mini = [];

for e = 1:Ler_cells

    e_idx = Ler_idx(e);
    
    Ler_aligned_mini = padmat(CELL_EVENTS(e_idx).mean_aligned_mini, Ler_aligned_mini, 2);
    rise_times = CELL_EVENTS(e_idx).cell_rise_times;
    decay_taus = CELL_EVENTS(e_idx).cell_tau;
    
    Ler_rise_time(e) = nanmean(rise_times);
    Ler_tau(e) = nanmean(decay_taus);
    
end
    
mean_Ler_rise = nanmean(Ler_rise_time);
mean_Ler_tau  = nanmean(Ler_tau);
mean_Ler_mini = nanmean(Ler_aligned_mini, 2);


%% PLOTTING for both timepoints
setFigureDefaults; % set default plotting parameters

% plotting colors
ctrl_col = [.7 .7 .7];
er_col = [112, 166, 217] ./ 255;

% figure
kinfig = figure();
set(kinfig,'position',[.04 .1 .9 .6]);

% plot average mini WFs for each condition
plot(1:150,mean_ctrl_mini(1:150).*1e12,'k','linewidth',3);
hold on;
plot(201:350,mean_er_mini(1:150).*1e12,'color',er_col,'linewidth',3);
plot(401:550,mean_Lctrl_mini(1:150).*1e12,'color','k','linewidth',3);
plot(601:750,mean_Ler_mini(1:150).*1e12,'color',er_col,'linewidth',3);
box off
line([-50 800],[min(mean_ctrl_mini)*1e12 min(mean_ctrl_mini)*1e12],...
    'linestyle','--','color','k','linewidth',1.5);
set(gca,'ylim',[-12 2],'xlim',[-50 800],'xtick',-100:50:800,'xticklabel',{''});
ylabel('Amplitude (pA)');
xlabel('Time (1 tick = 5 msec)');

%% Fig S4G - plot peak-scaled waveforms
pkfig = figure();
set(pkfig,'position',[.04 .1 .4 .6]);

pk_C2 = mean_ctrl_mini(1:150) ./ abs(min(mean_ctrl_mini(1:150)));
pk_R2 = mean_er_mini(1:150) ./ abs(min(mean_er_mini(1:150)));
pk_C4 = mean_Lctrl_mini(1:150) ./ abs(min(mean_Lctrl_mini(1:150)));
pk_R4 = mean_Ler_mini(1:150) ./ abs(min(mean_Ler_mini(1:150)));

plot(1:150,pk_C2,'k','linewidth',1.5);
hold on;
plot(1:150,pk_R2,'color',er_col,'linewidth',1.5,'linestyle','-');
plot(1:150,pk_C4,'color','k','linewidth',1.5,'linestyle','--');
plot(1:150,pk_R4,'color',er_col,'linewidth',1.5,'linestyle','--');
box off

set(gca,'ylim',[-1.2 0.2],'xlim',[-50 200],'xtick',-50:50:200,'xticklabel',{''});
ylabel('Amplitude (peak-scaled)');
xlabel('Time (1 tick = 5 msec)');

%% Fig S4E - plot rise time
ctrl_col = [.5 .5 .5];
reop_col = [7, 95, 179] ./ 255;
cset1 = [ctrl_col;ctrl_col;reop_col;reop_col];
cset2 = [ctrl_col;[1 1 1];reop_col;[1 1 1]];
labels = {'ER2 Control','ER4 Control','ER2 Re-open','ER4 Re-open'};

rise_dat = padcat(ctrl_rise_time',er_rise_time',Lctrl_rise_time',Ler_rise_time').*1e3;
timefig = plot_mini_data_4cond(rise_dat,'Rise',labels,cset1,cset2);
set(timefig,'position',[.1 .08 .3 .7]);

groupz = [repmat(1,size(ctrl_rise_time')); repmat(2,size(er_rise_time'));...
    repmat(3,size(Lctrl_rise_time')); repmat(4,size(Ler_rise_time'))];
risedat = [ctrl_rise_time'; er_rise_time'; Lctrl_rise_time'; Ler_rise_time'];

% STATS
[p_rise,~,rise_stats] = kruskalwallis(risedat,groupz,'off');
mc_rise = multcompare(rise_stats,'display','off');
disp(mc_rise);

%% Fig S4F - plot tau decay time constant
tau_dat = padcat(ctrl_tau',er_tau',Lctrl_tau',Ler_tau').*1e3;
timefig = plot_mini_data_4cond(tau_dat,'Tau',labels,cset1,cset2);
set(timefig,'position',[.1 .08 .3 .7]);


groupz = [repmat(1,size(ctrl_tau')); repmat(2,size(er_tau'));...
    repmat(3,size(Lctrl_tau')); repmat(4,size(Ler_tau'))];
kdat = [ctrl_tau'; er_tau'; Lctrl_tau'; Ler_tau'];

% STATS
[p_tau,~,tau_stats] = kruskalwallis(kdat,groupz,'off');
mc_tau = multcompare(tau_stats,'display','off');
disp(mc_tau)

