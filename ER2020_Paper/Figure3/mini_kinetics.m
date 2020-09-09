%% Mini_kinetics
%
% Alejandro Torrado Pacheco - 2019
%
% Extract kinetics of mEPSC waveforms from events

% clear workspace
clearIDE

% this code needs to be run for 2 different time points: 'er2' and 'er4'.
% In each case neurons recorded 2 or 4 days after eye re-opening will be
% analyzed. Both control and re-opened hemisphere neurons will be analyzed.
% The script goes through every event for every cell and computes rise and
% decay time. It also calculates the rise and decay time for the average WF
% for a given cell (so you can compare the two). All of the data is saved
% in a structure for later plotting.

% change this to er2 to get early ER data
timepoint = 'er4';
% saving flag. Set to 0 to NOT save data
do_save = 1;

% load the data
if ismac
    maindir = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig3/analyzed_traces';
elseif ispc
    maindir = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig3\analyzed_traces';
end

cell_file = [maindir filesep 'CELLS_' timepoint '.mat'];

if ~exist('CELLS','var')
    fprintf('Loading cells...\n');
    tic
    loadcells = load(cell_file);
    toc
    CELLS = loadcells.CELLS;
end
n_cells = numel(CELLS);

%%
% clearvars -except CELLS n_cells
% clc

% parameters
ampl_thresh = 5.0e-12; % pA - mini detection threshold
% peak_frac is a unitless variable that controls the alignment of mini
% events. For instance, if set to 0.5, minis are aligned so that the 50%
% value of the peak is aligned; if peak_frac = 1.0, minis are aligned so
% the peak is aligned.
peak_frac   = 1.0;
end_baseline_thresh = 3.0e-12;
shift_axis = 40;
sample_rate = 1e4; % 10 kHz for all my data
decay_fit_samples = 100;

test_plots = 0;
setFigureDefaults;

% loop through cells for this timepoint
for ii = 1:n_cells
    fprintf('\nCell %u of %u.\n',ii,n_cells);
    
    % get the data
    mini_data = CELLS(ii).mini_data;
    
    % find number of traces recorded for this cell
    n_traces = size(mini_data,1);
    
    % initialize variables
    cell_events = [];
    cell_starts = [];
    cell_peaks  = [];
    cell_adj_peaks = [];
    cell_adj_ampl = [];
    cell_ampl   = [];
    smooth_events_cell = [];
    adj_events_cell = [];
    shifted_events_cell = [];
    rejected_events_cell = [];
    
    % loop through traces
    for t = 1:n_traces
        
        % find events in this trace
        trace_events = mini_data{t,2}.events;
        [n_samp, n_ev] = size(trace_events);
        
        % find index of start of events, as well as peak and amplitude
        % obtained from MiniAnalysisGUI
        trace_starts = mini_data{t,2}.event_start_ind;
        trace_peaks  = mini_data{t,2}.sm_pk_ind;
        trace_ampl   = mini_data{t,2}.amp;
        
        % concatenate into master variable
        cell_events = padmat(cell_events, trace_events, 2);
        cell_starts = [cell_starts trace_starts];
        cell_peaks  = [cell_peaks trace_peaks];
        cell_ampl   = [cell_ampl trace_ampl];
        
    end
    
    % found this many events for this cell
    n_cell_events = size(cell_events,2);
    fprintf('Found %u events for this cell.\n', n_cell_events);
    
    % now loop through events
    for m = 1:n_cell_events
        % get the raw mEPSC waveform, its start and peak indices, and
        % amplitude
        mini_raw = cell_events(:,m);
        st_idx = cell_starts(m);
        pk_idx = cell_peaks(m);
        mini_amp = cell_ampl(m)*1e-12;
        
        % smooth the mini WF
        mini_smooth = smooth(mini_raw(~isnan(mini_raw)),21,'sgolay',2);
        smooth_events_cell = padmat(smooth_events_cell, mini_smooth, 2);
        
        % check if it returns to BL after peak
        end_baseline = nanmean(mini_raw(50:end));
        
        if mini_amp >= ampl_thresh && abs(end_baseline) <= end_baseline_thresh
            % find the baseline current and subtract from the mini waveform
            baseline_current = nanmedian(mini_smooth(1:st_idx));
            mini = mini_smooth - baseline_current;
            % add the smoothed mini to the adjusted events master variable
            % (adj_events_cell)
            adj_events_cell = padmat(adj_events_cell, mini, 2);
            
            % find the adjusted start and peak indices
            adj_st_amp = mini(st_idx);
            [adj_pk_amp,adj_pk_idx] = min(mini(st_idx:pk_idx+10));
            adj_pk_idx = adj_pk_idx + st_idx - 1;
            
            % align the mini according to the peak_frac parameter
            align_idx(m) = find(mini <= -(abs(adj_pk_amp) - abs(adj_st_amp))*peak_frac,1,'first');
            
            % save the adjusted peak index and amplitude
            cell_adj_peaks(m) = adj_pk_idx;
            cell_adj_ampl(m) = adj_pk_amp;
            
            % compute the rise time
            mini_rise_time(m) = (adj_pk_idx - st_idx) / sample_rate;
            
            % compute the decay tau for each mini
            % This is done by fitting an exponential function to the decay
            % phase of the mini event
            max_idx = min(adj_pk_idx+decay_fit_samples,size(mini,1));
            decay_fit_idx = adj_pk_idx:max_idx;
            decay_fit_data = mini(decay_fit_idx);
            decay_fit_data = decay_fit_data(~isnan(decay_fit_data)).*10^12;
            
            % start with an estimate of 0.005, then setup the fit
            tau_est = 0.005;
            time_pts = (0:(1/sample_rate):(length(decay_fit_data)/sample_rate)-(1/sample_rate)).';
            s = fitoptions('Method', 'NonlinearLeastSquares', ...
                'StartPoint', [mean(decay_fit_data(1:5)), tau_est],...
                'Lower', [mean(decay_fit_data(1:3))*1.1, 0.0005],...
                'Upper', [mean(decay_fit_data(1:3))*0.9, 0.010]);
            f = fittype('a*(exp(-x/b))','options',s);
            
            % do the fitting and get coefficients
            [exp_fit,gof] = fit(time_pts,decay_fit_data,f);
            cval = coeffvalues(exp_fit);
            mini_tau(m) = cval(2);
            mini_rmse(m) = gof.rmse;
            Y_mini = cval(1).*(exp(-time_pts./cval(2)));
            
            % test_plots is a flag that can be used to visualize the
            % exponential fit
            if test_plots
                figure();
                plot(mini,'-k','linewidth',2);
                hold on;
                plot(adj_pk_idx, adj_pk_amp,'rs','markerfacecolor','r','markersize',12);
                plot(st_idx, adj_st_amp,'bs','markerfacecolor','b','markersize',12);
                plot(decay_fit_idx,Y_mini.*1e-12,':m','linewidth',1.5);
                text(100,adj_pk_amp,sprintf('RMSE: %.4f',gof.rmse),...
                    'fontsize',18);
                uiwait(gcf);
            end
            
            
            % -------------------------------------------------------------
            
            % shift the minis for averaging
            shift_idx = shift_axis - align_idx(m) - 1;
            
            if shift_idx > 0
                shifted_event = [nan(shift_idx,1); mini];
            elseif shift_idx < 0
                shifted_event = mini(-shift_idx:end);
            else
                shifted_event = mini(shift_idx+1:end);
            end
            
            shifted_events_cell = padmat(shifted_events_cell, shifted_event, 2);
            
        else
            % if bad event, set everything to NaN
            adj_events_cell = padmat(adj_events_cell, nan(size(adj_events_cell,1),1), 2);
            cell_adj_peaks(m) = NaN;
            cell_adj_ampl(m) = NaN;
            mini_rise_time(m) = NaN;
            mini_tau(m) = NaN;
            mini_rmse(m) = NaN;
            
            rejected_events_cell = padmat(rejected_events_cell, mini_smooth, 1);
        end
        
    end
    
    % calculate rise time for average mini for this cell
    mean_aligned_mini = nanmean(shifted_events_cell,2);
    
    % find local max to approximate start of mini
    [~,mean_pk] = min(mean_aligned_mini);
    [~,local_max] = max(mean_aligned_mini(1:mean_pk));
    dm = diff(mean_aligned_mini(1:mean_pk));
    
    % calculate rise time
    if all(dm(local_max+1:end) < 0)
        mean_rise_time = (mean_pk - local_max) / sample_rate;
    else
        while ~all(dm(local_max+1:end) < 0)
            prev_max = local_max;
            [~,local_max] = max(mean_aligned_mini(prev_max+1:mean_pk));
            local_max = prev_max + local_max;
        end
        mean_rise_time = (mean_pk - local_max) / sample_rate;
    end
    
    
    
    % Calculate decay time constant Tau for this cell's mean aligned mini
    
    decay_fit_idx = (mean_pk):(mean_pk+decay_fit_samples);
    decay_fit_data = mean_aligned_mini(decay_fit_idx);
    decay_fit_data = decay_fit_data(~isnan(decay_fit_data)).*10^12;
    
    tau_est = 0.005;
    time_pts = (0:(1/sample_rate):(length(decay_fit_data)/sample_rate)-(1/sample_rate)).';
    s = fitoptions('Method', 'NonlinearLeastSquares', ...
        'StartPoint', [mean(decay_fit_data(1:5)), tau_est],...
        'Lower', [mean(decay_fit_data(1:3))*1.1, 0.0005],...
        'Upper', [mean(decay_fit_data(1:3))*0.9, 0.010]);
    f = fittype('a*(exp(-x/b))','options',s);
    
    [exp_fit,gof] = fit(time_pts,decay_fit_data,f);
    cval = coeffvalues(exp_fit);
    mean_tau = cval(2);
    Y = cval(1).*(exp(-time_pts./cval(2)));
    
    cell_avg_rise_time = nanmean(mini_rise_time);
    
    fprintf('\nMean WF rise: %.4f\nAvg cell rise time: %.4f\n\n',...
        mean_rise_time, cell_avg_rise_time);
    
    % save all of the data
    CELL_EVENTS(ii).raw_events = cell_events;
    CELL_EVENTS(ii).event_starts = cell_starts;
    CELL_EVENTS(ii).event_peaks = cell_peaks;
    CELL_EVENTS(ii).adj_peaks = cell_adj_peaks;
    CELL_EVENTS(ii).adj_ampl = cell_adj_ampl;
    CELL_EVENTS(ii).mini_smooth = mini_smooth;
    CELL_EVENTS(ii).mini_adj = mini;
    CELL_EVENTS(ii).align_idx = align_idx;
    CELL_EVENTS(ii).shifted_events = shifted_events_cell;
    CELL_EVENTS(ii).rejected_events = rejected_events_cell;
    CELL_EVENTS(ii).cell_rise_times = mini_rise_time;
    CELL_EVENTS(ii).cell_tau = mini_tau;
    CELL_EVENTS(ii).cell_tau_rmse = mini_rmse;
    CELL_EVENTS(ii).mean_aligned_mini = mean_aligned_mini;
    CELL_EVENTS(ii).mean_rise_time = mean_rise_time;
    CELL_EVENTS(ii).mean_tau = mean_tau;
    
    
    % more test plots
    if test_plots
        
        c1 = [.47 .67 .19];
        c2 = [0.30 .75 .93];
        c3 = [.49 .18 .56];
        
        figure(1);
        set(gcf,'position',[.06 .08 .8 .9]);
        
        subplot(1,2,1);
        plot(shifted_events_cell,'-','linewidth',0.25,'color',[0 0 0 0.1]);
        hold on;
        plot(mean_aligned_mini,'-','linewidth',3.5,'color','k');
        box off;
        
        subplot(1,2,2);
        plot(mean_aligned_mini,'-k','linewidth',2.5);
        hold on;
        plot(local_max,mean_aligned_mini(local_max),'o','color',c1,'markerfacecolor',c1,'markersize',10);
        plot(mean_pk,mean_aligned_mini(mean_pk),'o','color',c3,'markerfacecolor',c3,'markersize',10);
        plot(decay_fit_idx,Y.*1e-12,'-','color',c2,'linewidth',2.5);
        box off;
        
        uiwait(gcf);
        
    end
    
    
end


%% final save of all the data
if do_save
    save_file = [maindir filesep 'CELL_EVENTS_TEMP_' timepoint '.mat'];
    fprintf('Saving...\n')
    tic
    save(save_file,'CELL_EVENTS','-v7.3');
    toc
end


