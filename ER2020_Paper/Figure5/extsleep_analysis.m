%% extsleep_analysis.m
%
% Alejandro Torrado Pacheco - 2018
%
% Extended Sleep Analysis script
%
% This code does the analysis shown in Figure 5. It can be used to
% reproduce the results in Fig 5C, 5D, 5E, 5F.

%% SETUP
% clear workspace
clearIDE

% load contcell variable
contcell_file = 'CONTCELL_recov.mat';
if ismac
    contcell_dir = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig1';
elseif ispc
    contcell_dir = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig1';
end
load([contcell_dir filesep contcell_file]);

% separate neurons and statetimes, and get animal info
MASTER = CONTCELL_recov.MASTER;
STATES = CONTCELL_recov.STATETIMES;
DAYSTART = STATES.DAYSTART;
anim_fields = fieldnames(STATES);
anim_fields(strcmp(anim_fields,'DAYSTART')) = [];
n_anims = numel(anim_fields); % minus daystart field

cell_anims = {MASTER.animal};

% colors for plotting
c_rem   = [25 181 149]./255;
c_nrem  = [131 49 146]./255;
c_aw    = [201 28 101]./255;
c_qw    = [247 148 41]./255;

% analysis flags
save_delta_data = 0;

% period of interest (downward homeostatic recovery)
homeo_period_start = 8.0; % 8 is start of ER2
homeo_period_end = 10.0; % 10 is beginning of ER4

% cell selection parameters
 % select control (0) or re-opened (1) hemisphere
dep = 1; % change this to get Fig 5C and 5D (reopen) or 5E and 5F (control)
qthresh = 2; % only single units (quality 1 or 2)
perc_thresh = 0.75; % cells must be continuous (ON for 75% of rec time)
norm = 0; % don't normalize
bl_on = 1; % enforce cells being ON at baseline
dataset = 'recov'; % recov dataset
negpos = 0.39; % negpostime used to find RSUs (not FS cells)
tailslope = 0.005; % same as above

G_bin = 1; % seconds. Bin size for FR calculations.

% mean_t sets the normalizing mode:
%   0 for no normalization
%   1 for z-scoring to sleep epoch
%   2 for z-scoring to all recorded times for this cell
%   3 for fractional change
mean_t = 1;

% analysis parameters
dur_frac = 0.20; % only used if mean_t == 3
% extended sleep minimum duration
ext_sleep_time_thresh = 30; % minutes
% individual states duration thresholds
interrupt_threshold = 10; % seconds - get rid of states shorter than this value
duration_threshold = 30; % seconds - ignore states shorter than this (no FR calculation)

% time start and end for main FR calculation
t_start = 0;
t_stop = 12*24*3600;

% animal list selection
if dep
    anim_list = {'AT12','AT14','AT16','AT27','AT29','KH67'};
else
    anim_list = {'AT12','AT14','AT16','AT27','AT29','KH67','KH72','KH73','KH75'};
end

%% Loop through animals and find extended sleep periods

for aa = 1:n_anims
    animal = anim_fields{aa};
    
    % get statetimes for tis animal
    statetimes = STATES.(animal);
    
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
    
    % FIND EXTENDED SLEEP - see function
    ext_sleeps = find_extended_sleep(statetimes,ext_sleep_time_thresh);
    
    % save ext sleeps for this animal
    all_extsleeps{aa} = ext_sleeps;
    
    % find experiment start time for this animal
    all_anims = {CONTCELL_recov.MASTER.animal};
    first_cell = find(strcmp(animal,all_anims),1,'first');
    
    expstart_raw = unixtime(CONTCELL_recov.MASTER(first_cell).EXPTSTART);
    
    % find ZT0 on expt start day for this animal
    expstart_unix = [expstart_raw(1:3) 7 30 0];
    expstart_t = unixtime(expstart_unix);
    
    % save it
    all_expstart{aa} = expstart_t;
    
    clear  ext_sleeps expstart* statetimes
end

% initialize variables
all_NREM_zFR = [];
all_REM_zFR = [];


%% Loop through animals and do FR analysis
for aa = 1:n_anims
    animal = anim_fields{aa};
    fprintf('Animal %s.\n',animal);
    
    if any(strcmp(anim_list,animal))
        
        % define period of interest in seconds
        homeo_start = homeo_period_start*24*3600;
        homeo_end = homeo_period_end*24*3600;
        
        % find this animal's cells
        anim_rsu_idx = find(strcmp(cell_anims,animal));
        anim_cells = MASTER(anim_rsu_idx);
        
        % find the cells according to specified parameters using get_RSUs
        % function.
        rsu = get_RSUs(anim_cells,qthresh,dep,perc_thresh,norm,bl_on,dataset,...
            plot_cellsep,negpos,tailslope,0);
        n_RSU = size(rsu,2);
        
        anim_idx = find(strcmp(anim_fields,animal));
        
        % RETRIEVE STATE DATA
        
        ext_sleeps = all_extsleeps{anim_idx};
        expstart_t = all_expstart{anim_idx};
        
        statetimes = STATES.(animal);
        
        % ENSURE THIS PROCESSING IS SAME AS ABOVE - IT HAS TO MATCH BECAUSE
        % THE STATETIMES ARRAY PASSED ABOVE IS WHERE THE EXT SLEEP INDICES
        % WERE FOUND
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
        
        % get animal daystart variable
        daystart = DAYSTART.(animal);
        
        % find number of extended sleep episodes
        n_sleeps = size(ext_sleeps,1);
        
        % initialize variables
        all_sleep_dur{anim_idx} = nan(n_sleeps, n_RSU);
        deltaFR_NREM{anim_idx} = nan(n_sleeps,n_RSU);
        deltaFR_REM{anim_idx} = nan(n_sleeps,n_RSU);
        percent_NREM{anim_idx} = nan(n_sleeps,n_RSU);
        percent_REM{anim_idx} = nan(n_sleeps,n_RSU);
        num_REM_eps{anim_idx} = nan(n_sleeps, n_RSU);
        
        
        % loop through ext sleep episodes
        for tt = 1:n_sleeps
            
            fprintf('Extended sleep epoch %u of %u.\n',tt,n_sleeps);
            
            % find the states occurring during this episode
            st_sleep = statetimes(ext_sleeps(tt,1):ext_sleeps(tt,2)+1,:);
            
            % align to BL1
            st_sleep(:,2) = st_sleep(:,2) - expstart_t + daystart*24*3600;
            
            
            try
                % find duration and percent of REM and NREM states within
                % ext sleep episode
                sleep_durs = diff(st_sleep(:,2));
                total_dur = sum(sleep_durs);
                
                nrem_idx = find(st_sleep(:,1) == 2);
                if nrem_idx(end) > size(sleep_durs,1)
                    nrem_idx(end) = [];
                end
                nrem_dur = sum(sleep_durs(nrem_idx));
                nrem_percent = 100 * (nrem_dur/total_dur);
                
                rem_idx = find(st_sleep(:,1) == 1);
                if ~isempty(rem_idx)
                    if rem_idx(end) > size(sleep_durs,1)
                        rem_idx(end) = [];
                    end
                    rem_dur = sum(sleep_durs(rem_idx));
                    rem_percent = 100 * (rem_dur/total_dur);
                end     
            catch
                keyboard;
            end
            
            % get ext sleep duration
            sleep_start = st_sleep(1,2);
            sleep_end = st_sleep(end,2);
            ext_sleep_duration = sleep_end - sleep_start;
            
            % if this ext sleep episode is within the homeostatic recovery
            % window, proceed
            if sleep_start >= homeo_start && sleep_end <= homeo_end
                
                % find number and indices of NREM and REM epochs
                fprintf('good\n');
                NREM_idx = find(st_sleep(:,1) == 2);
                n_NREM = numel(NREM_idx);
                
                REM_idx = find(st_sleep(:,1) == 1);
                n_REM = numel(REM_idx);
                
                idx_NR_first = NREM_idx(1);
                idx_NR_last = NREM_idx(end);
                
                if n_REM > 0
                    idx_REM_first = REM_idx(1);
                    idx_REM_last = REM_idx(end);
                end
                
                % loop through cells
                for cc = 1:n_RSU
                    if size(rsu(cc).onTime,1) > 1
                        fprintf('Multiple on/off times.\n');
                        keyboard;
                    end
                    
                    % get spike times
                    spikes = rsu(cc).time;
                    
                    % calculate FR according to analysis parameters
                    switch mean_t
                        case 0
                            mean_long_t0 = sleep_start;
                            mean_long_t1 = sleep_end;
                            do_z = 0;
                            delta_frac = 0;
                        case 1
                            mean_long_t0 = sleep_start;
                            mean_long_t1 = sleep_end;
                            do_z = 1;
                            delta_frac = 0;
                        case 2
                            mean_long_t0 = rsu(cc).onTime(1);
                            mean_long_t1 = rsu(cc).offTime(end);
                            do_z = 1;
                            delta_frac = 0;
                        case 3
                            mean_long_t0 = sleep_start;
                            mean_long_t1 = sleep_end;
                            do_z = 0;
                            delta_frac = 1;
                    end
                    
                    % simple histogram count to get FR in bins of size G_bin
                    rate_tmp = histc(spikes,t_start:G_bin:t_stop) ./ G_bin;
                    
                    % process on/off times for this cell
                    tlist = [rsu(cc).onTime rsu(cc).offTime];
                    rate_onoff = processOnOffTimes_ATP(rate_tmp,tlist,G_bin,0,'recov',0);
                    
                    % get the FR within the normalizing episode
                    rate_long = rate_onoff(round(mean_long_t0/G_bin):round(mean_long_t1/G_bin));
                    
                    % get the mean rate and std dev within the normalizing episode
                    frmean_long = nanmean(rate_long);
                    frstd_long  = std(rate_long,0,'omitnan');
                    
                    % cycle thru NREM epochs and find mean FR in each
                    for nr = 1:n_NREM
                        
                        % find start and end time of current NREM epoch
                        NR_t0 = st_sleep(NREM_idx(nr),2);
                        NR_t1 = st_sleep(NREM_idx(nr)+1,2);
                        
                        % if cell is on, and this epoch is longer than
                        % threshold, proceed
                        if NR_t0 >=  rsu(cc).onTime(1) && NR_t1 <= rsu(cc).offTime(end) && (NR_t1-NR_t0 >= duration_threshold)
                            
                            % convert time to index according to FR bin size
                            NR_idx0 = floor( (NR_t0-mean_long_t0) / G_bin );
                            if NR_idx0 == 0, NR_idx0 = 1; end
                            NR_idx1 = ceil( (NR_t1-mean_long_t0) / G_bin );
                            
                            % find mean FR in this NREM epoch
                            NR_FR = nanmean(rate_long(NR_idx0:NR_idx1));
                            % depending on analysis parameters, z-score or not
                            % NOTE: this will either z-score to the ext
                            % sleep episode or to the whole expt time,
                            % depending on the chosen value of mean_t
                            if do_z
                                NR_zFR = (NR_FR - frmean_long) / frstd_long;
                            else
                                NR_zFR = NR_FR;
                            end
                        else
                            NR_zFR = NaN;
                        end
                        
                        % save NREM FR and time from ext sleep start
                        all_NREM_zFR{anim_idx,tt}(cc,nr) = NR_zFR;
                        all_NREM_offT{anim_idx,tt}(cc,nr) = NR_t0 - sleep_start;
                    end
                    
                    % REPEAT FOR REM
                    % cycle thru REM epochs and find mean FR in each
                    if n_REM > 0
                        for r = 1:n_REM
                            
                            R_t0 = st_sleep(REM_idx(r),2);
                            R_t1 = st_sleep(REM_idx(r)+1,2);
                            if R_t0 >=  rsu(cc).onTime(1) && R_t1 <= rsu(cc).offTime(end) && (R_t1-R_t0 >= duration_threshold)
                                R_idx0 = floor( (R_t0-mean_long_t0) / G_bin );
                                if R_idx0 == 0, R_idx0 = 1; end
                                R_idx1 = ceil( (R_t1-mean_long_t0) / G_bin );
                                
                                R_FR = nanmean(rate_long(R_idx0:R_idx1));
                                if do_z
                                    R_zFR = (R_FR - frmean_long) / frstd_long;
                                else
                                    R_zFR = R_FR;
                                end
                            else
                                R_zFR = NaN;
                            end
                            
                            all_REM_zFR{anim_idx,tt}(cc,r) = R_zFR;
                            all_REM_offT{anim_idx,tt}(cc,r) = R_t0 - sleep_start;
                            num_REM_eps{anim_idx}(tt,cc) = n_REM;
                        end
                    end
                    
                    % find difference in FR between last NREM and first NREM
                    
                    % set up here depends on chosen value of mean_t
                    if mean_t ~= 3
                        NR_first_t0 = st_sleep(idx_NR_first,2);
                        NR_first_t1 = st_sleep(idx_NR_first+1,2);
                        NR_last_t0 = st_sleep(idx_NR_last,2);
                        NR_last_t1 = st_sleep(idx_NR_last+1,2);
                        
                        NR_dur_first = NR_first_t1-NR_first_t0;
                        NR_dur_last = NR_last_t1-NR_last_t0;
                    elseif mean_t == 3
                        sleep_dur_frac = dur_frac * ext_sleep_duration;
                        
                        NR_first_t0 = st_sleep(1,2);
                        NR_first_t1 = NR_first_t0 + sleep_dur_frac;
                        NR_last_t1 = st_sleep(end,2);
                        NR_last_t0 = NR_last_t1 - sleep_dur_frac;
                    end
                    
                    % do the calculation
                    if NR_first_t0 >=  rsu(cc).onTime(1) && NR_last_t1 <= rsu(cc).offTime(end) && ...
                        NR_dur_first >= duration_threshold && NR_dur_last >= duration_threshold
                    
                        % find indices of first and last NR period
                        NR_first_idx0 = floor( (NR_first_t0-mean_long_t0) / G_bin );
                        if NR_first_idx0 == 0, NR_first_idx0 = 1; end
                        NR_first_idx1 = ceil( (NR_first_t1-mean_long_t0) / G_bin );
                        
                        NR_last_idx0 = floor( (NR_last_t0-mean_long_t0) / G_bin );
                        NR_last_idx1 = ceil( (NR_last_t1-mean_long_t0) / G_bin );
                        
                        % find men rate in each one
                        NR_first_FR = nanmean(rate_long(NR_first_idx0:NR_first_idx1));
                        NR_last_FR = nanmean(rate_long(NR_last_idx0:NR_last_idx1));
                        
                        % depending on mean_t, do z-scoring or not
                        if do_z
                            NR_first_zFR = (NR_first_FR - frmean_long) / frstd_long;
                            NR_last_zFR = (NR_last_FR - frmean_long) / frstd_long;
                        else
                            NR_first_zFR = NR_first_FR;
                            NR_last_zFR = NR_last_FR;
                        end
                        
                        % if doing the fractional method, do this
                        if delta_frac
                            NR_delta_zFR = NR_last_zFR / (NR_first_zFR);
                            if NR_delta_zFR == Inf || NR_delta_zFR == 0
                                NR_delta_zFR = NaN;
                            end
                        else
                            NR_delta_zFR = NR_last_zFR - NR_first_zFR;
                        end
                    else
                        NR_delta_zFR = NaN;
                    end
                    % compile data
                    deltaFR_NREM{anim_idx}(tt,cc) = NR_delta_zFR;
                    all_sleep_dur{anim_idx}(tt,cc) = ext_sleep_duration;
                    
                    % find difference in FR between last REM and first REM
                    if n_REM > 0 % if there are any REM periods
                        R_first_t0 = st_sleep(idx_REM_first,2);
                        R_first_t1 = st_sleep(idx_REM_first+1,2);
                        R_last_t0 = st_sleep(idx_REM_last,2);
                        R_last_t1 = st_sleep(idx_REM_last+1,2);
                        
                        R_dur_first = R_first_t1 - R_first_t0;
                        R_dur_last = R_last_t1 - R_last_t0;
                        
                        % as above...
                        if R_first_t0 >=  rsu(cc).onTime(1) && R_last_t1 <= rsu(cc).offTime(end) && ...
                           R_dur_first >= duration_threshold && R_dur_last >= duration_threshold
                    
                            R_first_idx0 = floor( (R_first_t0-mean_long_t0) / G_bin );
                            if R_first_idx0 == 0, R_first_idx0 = 1; end
                            R_first_idx1 = ceil( (R_first_t1-mean_long_t0) / G_bin );
                            
                            R_last_idx0 = floor( (R_last_t0-mean_long_t0) / G_bin );
                            R_last_idx1 = ceil( (R_last_t1-mean_long_t0) / G_bin );
                            
                            R_first_FR = nanmean(rate_long(R_first_idx0:R_first_idx1));
                            R_last_FR = nanmean(rate_long(R_last_idx0:R_last_idx1));
                            
                            if do_z
                                R_first_zFR = (R_first_FR - frmean_long) / frstd_long;
                                R_last_zFR = (R_last_FR - frmean_long) / frstd_long;
                            else
                                R_first_zFR = R_first_FR;
                                R_last_zFR = R_last_FR;
                            end
                            
                            if delta_frac
                                R_delta_zFR = R_last_zFR / R_first_zFR;
                            else
                                R_delta_zFR = R_last_zFR - R_first_zFR;
                            end
                        else
                            R_delta_zFR = NaN;
                        end
                        deltaFR_REM{anim_idx}(tt,cc) = R_delta_zFR;
                    end
                    
                    percent_NREM{anim_idx}(tt,cc) = nrem_percent;
                    if n_REM > 0
                        percent_REM{anim_idx}(tt,cc) = rem_percent;
                    end
                end
                
            else
                fprintf('bad\n');
                all_NREM_zFR{anim_idx,tt} = NaN;
                all_REM_zFR{anim_idx,tt} = NaN;
                all_NREM_offT{anim_idx,tt} = NaN;
                all_REM_offT{anim_idx,tt} = NaN;
                
            end
        end
    end
end


%% compile the data
% change in FR between first and last NREM
NREM_all_deltas = cellfun(@(x) x(:),deltaFR_NREM,'UniformOutput',false);
cat_NREM_deltas = cat(1,NREM_all_deltas{:});
cat_NREM_deltas(cat_NREM_deltas > 5) = NaN;
delta_nonans = ~isnan(cat_NREM_deltas);
NREM_delta_nonan = cat_NREM_deltas(delta_nonans);

% change in FR between first and last REM
REM_all_deltas = cellfun(@(x) x(:),deltaFR_REM,'UniformOutput',false);
cat_REM_deltas = cat(1,REM_all_deltas{:});
cat_REM_deltas(cat_REM_deltas > 5) = NaN;
rem_delta_nonans = ~isnan(cat_REM_deltas);
REM_delta_nonan = cat_REM_deltas(rem_delta_nonans);

% ext sleep durations
all_sleep_dur_cellcat = cellfun(@(x) x(:),all_sleep_dur,'UniformOutput',false);
cat_sleep_dur = cat(1,all_sleep_dur_cellcat{:});
sleepdur_nonan = cat_sleep_dur(delta_nonans);

% percent time in NREM
NREM_all_percs = cellfun(@(x) x(:),percent_NREM,'UniformOutput',false);
cat_NREM_percs = cat(1,NREM_all_percs{:});
delta_nonans = ~isnan(cat_NREM_deltas);
NREM_perc_nonan = cat_NREM_percs(delta_nonans);

% eprcent time in REM
REM_all_percs = cellfun(@(x) x(:),percent_REM,'UniformOutput',false);
cat_REM_percs = cat(1,REM_all_percs{:});
delta_nonans = ~isnan(cat_NREM_deltas);
REM_perc_nonan = cat_REM_percs(delta_nonans);

% number of REM episodes
num_REM_all = cellfun(@(x) x(:),num_REM_eps,'UniformOutput',false);
cat_REM_num = cat(1,num_REM_all{:});
delta_nonans = ~isnan(cat_NREM_deltas);
num_REM_nonan = cat_REM_num(delta_nonans);
num_REM_remnonan = cat_REM_num(rem_delta_nonans);

% sleep durations by cell - useful for later analysis
sleepdur_bycell = cellfun(@(x) nanmean(x,2),all_sleep_dur,'UniformOutput',false);
sleepdur_bycell_cat = cat(1,sleepdur_bycell{:});

% compile by cell
% NREM data averaged across epochs for a given cell
NREM_cellcat = cellfun(@(x) nanmean(x,2),deltaFR_NREM,'UniformOutput',false);
NREM_cellcat_sem = cellfun(@(x) nanstd(x,[],2) ./ sqrt(sum(~isnan(x),2)),deltaFR_NREM,'UniformOutput',false);
NREM_perccat = cellfun(@(x) nanmean(x,2),percent_NREM,'UniformOutput',false);
% collecting all epochs instead of averaging by cell
NREM_cellmean = cat(1,NREM_cellcat{:});
NREM_cellsem = cat(1,NREM_cellcat_sem{:});
NREM_cellperc = cat(1,NREM_perccat{:});
NREM_all_percent = NREM_cellperc(~isnan(NREM_cellmean));
NREM_all_deltaFR = NREM_cellmean(~isnan(NREM_cellmean));
NREM_all_deltaSEM = NREM_cellsem(~isnan(NREM_cellmean));
sleepdur_bycell_nonan = sleepdur_bycell_cat(~isnan(NREM_cellmean));
n_epochs_NREM = sum(~isnan(NREM_cellmean));


% as above but for REM
REM_cellcat = cellfun(@(x) nanmean(x,2),deltaFR_REM,'UniformOutput',false);
REM_perccat = cellfun(@(x) nanmean(x,2),percent_REM,'UniformOutput',false);
REM_cellmean = cat(1,REM_cellcat{:});
REM_cellperc = cat(1,REM_perccat{:});
REM_all_percent = REM_cellperc(~isnan(REM_cellmean));
REM_all_deltaFR = REM_cellmean(~isnan(REM_cellmean));
n_epochs_REM = sum(~isnan(REM_cellmean));

% average across cells - NREM and REM
NREM_mean_c = nanmean(NREM_cellmean);
NREM_sem_c = std(NREM_cellmean,0,'omitnan') / sqrt(n_epochs_NREM - 1);
NREM_mean = nanmean(NREM_delta_nonan);
NREM_sem = std(NREM_delta_nonan,0,'omitnan') / sqrt(sum(~isnan(NREM_delta_nonan)) - 1);

REM_mean_c = nanmean(REM_cellmean);
REM_sem_c = std(REM_cellmean,0,'omitnan') / sqrt(n_epochs_REM - 1);
REM_mean = nanmean(REM_delta_nonan);
REM_sem = std(REM_delta_nonan,0,'omitnan') / sqrt(sum(~isnan(REM_delta_nonan)) - 1);

%% z-scored FR vs time from start of ext sleep
% compile the data
cellcat_NREM_offT = cellfun(@(x) x(:),all_NREM_offT,'UniformOutput',false);
cellcat_NREM_zFR = cellfun(@(x) x(:),all_NREM_zFR,'UniformOutput',false);

cat_NREM_offT = cat(1,cellcat_NREM_offT{:});
cat_NREM_zFR = cat(1,cellcat_NREM_zFR{:});

fr_notnans = ~isnan(cat_NREM_zFR);

NREM_offT_nonan = cat_NREM_offT(fr_notnans);
NREM_zFR_nonan = cat_NREM_zFR(fr_notnans);

% this is the whole data set - lots of points!
all_data = [NREM_zFR_nonan NREM_offT_nonan];
all_data_sort = sortrows(all_data,2);

% separate into ngroups of equal size (equal # of datapoints per group)
ngroups = 10;
group_split = round(linspace(0,numel(NREM_offT_nonan),ngroups+1));
for uu = 1:ngroups
    this_group = all_data_sort(group_split(uu)+1 : group_split(uu+1), 1);
    group_means(uu) = nanmean(this_group);
    group_sem(uu) = std(this_group,0,'omitnan') / sqrt(numel(this_group)-1);
    group_means_offT(uu) = nanmean(all_data_sort(group_split(uu)+1 : group_split(uu+1), 2));
end

% do the correlation analysis
[nr_rho,nr_p_rho] = corr(NREM_zFR_nonan, NREM_offT_nonan);
nr_asterisks = get_asterisks_from_pval(nr_p_rho);
nr_lincoeff = polyfit(NREM_offT_nonan, NREM_zFR_nonan, 1);

% for mean_t different than 3 (i.e. no fractional change)
if mean_t ~= 3
    
    % DO THE PLOTS
    
    
    %% Fig 5C - NREM FR vs time correlation
    % note this because Fig 5E if dep == 0
    markersize = 4;
    mylims = 1;
    
    % this figure actually has more than is shown in the preprint's Fig 5.
    % It shows the raw data (all points) in the top subplot, and the
    % grouped data in the bottom one.
    setFigureDefaults;
    figure();
    set(gcf,'position',[.2 .1 .3 .8]);
    subplot(2,1,1);
    % all points subplot
    if dep
        scatter(NREM_offT_nonan,NREM_zFR_nonan,markersize^2,c_nrem,'filled');
    else
        scatter(NREM_offT_nonan,NREM_zFR_nonan,markersize^2,c_nrem);
    end
    hold on;
    % axes formatting
    if mylims
        xlims = [-500 6000];
        ylims = [-2 2];
        xl2 = [-250 4000];
        yl2 = [-0.2 0.4];
    else
        xlims = get(gca,'xlim');
        ylims = get(gca,'ylim');
        xl2 = get(gca,'xlim');
        yl2 = get(gca,'ylim');
    end
    X = xlims(1):10:xlims(end);
    Y = nr_lincoeff(2) + nr_lincoeff(1).*X;
    plot(X,Y,'--k','linewidth',1.5);
    set(gca,'xlim',xlims,'ylim',ylims);
    ylabel('Firing rate (z)','fontsize',21);
    
    % grouped subplot
    subplot(2,1,2);
    if dep
        mfc = c_nrem;
    else
        mfc = 'none';
    end
    errorbar(group_means_offT,group_means,group_sem,'o','capsize',12,'markersize',8,...
        'linestyle','none','markerfacecolor',mfc,'linewidth',1.5,'color',c_nrem);
    hold on;
    box off
    plot(X,Y,'--k','linewidth',1.5);
    if mylims
        set(gca,'xlim',xl2,'ylim',yl2);
        text(xl2(2)*.8,yl2(2)*.9,sprintf('r = %.4f\np = %.4f',nr_rho,nr_p_rho),'fontsize',16);
    end
    ylabel('Firing rate (z)','fontsize',21);
    xlabel('Time from start of extended sleep (s)','fontsize',20);
    
    %% Fig 5C - REM FR vs time correlation
    % everything as above, but for REM
    cellcat_REM_offT = cellfun(@(x) x(:),all_REM_offT,'UniformOutput',false);
    cellcat_REM_zFR = cellfun(@(x) x(:),all_REM_zFR,'UniformOutput',false);
    
    cat_REM_offT = cat(1,cellcat_REM_offT{:});
    cat_REM_zFR = cat(1,cellcat_REM_zFR{:});
    
    REM_fr_notnans = ~isnan(cat_REM_zFR);
    
    REM_offT_nonan = cat_REM_offT(REM_fr_notnans);
    REM_zFR_nonan = cat_REM_zFR(REM_fr_notnans);
    
    REM_all_data = [REM_zFR_nonan REM_offT_nonan];
    REM_all_data_sort = sortrows(REM_all_data,2);
    
    % ngroups is inherited from above to ensure consistency between NREM
    % and REM analyses
    REM_group_split = round(linspace(0,numel(REM_offT_nonan),ngroups+1));
    for uu = 1:ngroups
        REM_this_group = REM_all_data_sort(REM_group_split(uu)+1 : REM_group_split(uu+1), 1);
        REM_group_means(uu) = nanmean(REM_this_group);
        REM_group_sem(uu) = std(REM_this_group,0,'omitnan') / sqrt(numel(REM_this_group)-1);
        REM_group_means_offT(uu) = nanmean(REM_all_data_sort(REM_group_split(uu)+1 : REM_group_split(uu+1), 2));
    end
    
    
    % do the correlation
    [rho,p_rho] = corr(REM_offT_nonan, REM_zFR_nonan);
    asterisks = get_asterisks_from_pval(p_rho);
    lincoeff = polyfit(REM_offT_nonan, REM_zFR_nonan, 1);
    
    %% plotting
    markersize = 4;
    mylims = 1;
    
    
    setFigureDefaults;
    figure();
    set(gcf,'position',[.2 .1 .3 .8]);
    subplot(2,1,1);
    if dep
        scatter(REM_offT_nonan,REM_zFR_nonan,markersize^2,c_rem,'filled');
    else
        scatter(REM_offT_nonan,REM_zFR_nonan,markersize^2,c_rem);
    end
    hold on;
    if mylims
        xlims = [-500 6000];
        ylims = [-2 2];
        xl2 = [-250 4000];
        yl2 = [-0.2 0.4];
    else
        xlims = get(gca,'xlim');
        ylims = get(gca,'ylim');
        xl2 = get(gca,'xlim');
        yl2 = get(gca,'ylim');
    end
    X = xlims(1):10:xlims(end);
    Y = lincoeff(2) + lincoeff(1).*X;
    plot(X,Y,'--k','linewidth',1.5);
    set(gca,'xlim',xlims,'ylim',ylims);
    ylabel('Firing rate (z)','fontsize',21);
    
    if dep
        mfc = c_rem;
    else
        mfc = 'none';
    end
    subplot(2,1,2);
    errorbar(REM_group_means_offT,REM_group_means,REM_group_sem,'o','capsize',12,'markersize',8,...
        'linestyle','none','markerfacecolor',mfc,'color',c_rem,'linewidth',1.5);
    hold on;
    box off;
    plot(X,Y,'--k','linewidth',1.5);
    if mylims
        set(gca,'xlim',xl2,'ylim',yl2);
        text(xl2(2)*.8,yl2(2)*.9,sprintf('r = %.4f\np = %.4f',rho,p_rho),'fontsize',16);
    end
    ylabel('Firing rate (z)','fontsize',21);
    xlabel('Time from start of extended sleep (s)','fontsize',20);
    
    %% Fig 5D - plotting the first-last change
    % Note this becomes Fig 5F if dep == 0 (control)
    mylims = 1;
    
    % STATS - t-test versus no change (mean = 0)
    [~,p_nr] = ttest(NREM_delta_nonan,0);
    [a_nr,fsz_nr] = get_asterisks_from_pval(p_nr);
    [~,p_rem] = ttest(REM_delta_nonan);
    [a_rem,fsz_rem] = get_asterisks_from_pval(p_rem);
    
    dfig = figure();
    set(dfig,'position',[.1 .2 .25 .4]);
    box off
    hold on;
    bw = .3;
    csz = 22;
    if dep
        mfc1 = c_nrem;
        mfc2 = c_rem;
        ec1 = 'none';
        ec2 = 'none';
    else
        mfc1 = 'none';
        mfc2 = 'none';
        ec1 = c_nrem;
        ec2 = c_rem;
    end
    
    % plot the bars
    nr_bar = bar(1,NREM_mean,bw,'edgecolor',ec1,'facecolor',mfc1,'linewidth',2);
    nr_err = errorbar(1,NREM_mean,NREM_sem,'linestyle','none','capsize',csz,...
        'color',c_nrem,'linewidth',2);
    rem_bar = bar(2,REM_mean,bw,'edgecolor',ec2,'facecolor',mfc2,'linewidth',2);
    rem_err = errorbar(2,REM_mean,REM_sem,'linestyle','none','capsize',csz,...
        'color',c_rem,'linewidth',2);
    if mylims
        xl = [.5 2.5];
        yl = [-0.2 0.2];
        yt = -0.8:.1:0.8;
    else
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        yt = get(gca,'ytick');
    end
    text(.75,0.1,sprintf('p = %.4f',p_nr),'fontsize',16);
    text(1.75,0.1,sprintf('p = %.4f',p_rem),'fontsize',16);
    set(gca,'xlim',xl,'ylim',yl,'xtick',[1 2],'xticklabel',{'NREM','REM'},...
        'ytick',yt);
    ylabel('Firing rate (z)','fontsize',22);
    
    
    
else
    %% if mean_t == 3, plotting is different
    markersize = 5;
    
    all_data_d = [NREM_delta_nonan sleepdur_nonan];
    all_data_sort_d = sortrows(all_data_d,2);
    
    ngroups = 10;
    group_split_d = round(linspace(0,numel(NREM_delta_nonan),ngroups+1));
    for uu = 1:ngroups
        this_group = all_data_sort_d(group_split_d(uu)+1 : group_split_d(uu+1), 1);
        group_means_d(uu) = nanmean(this_group);
        group_sem_d(uu) = std(this_group,0,'omitnan') / sqrt(numel(this_group)-1);
        group_means_offT_d(uu) = nanmean(all_data_sort_d(group_split_d(uu)+1 : group_split_d(uu+1), 2));
        group_max_d(uu) = max(all_data_sort_d(group_split_d(uu)+1 : group_split_d(uu+1), 2));
    end

    
    [rho_d,p_rho_d] = corr(sleepdur_nonan, NREM_delta_nonan);
%     asterisks_d = get_asterisks_from_pval(p_rho_d);
    lincoeff_d = polyfit(sleepdur_nonan, NREM_delta_nonan, 1);
    
    
    mylims = 1;
    
    setFigureDefaults;
    figure();
    set(gcf,'position',[.1 .1 .6 .7]);
    subplot(2,1,1);
    if dep
        scatter(sleepdur_nonan,NREM_delta_nonan,markersize^2,'k','filled');
    else
        scatter(sleepdur_nonan,NREM_delta_nonan,markersize^2,'k');
    end
    hold on;
    if mylims
        xlims = [-500 6000];
        ylims = [0 4];
        xl2 = [-500 6000];
        yl2 = [0 2.5];
    else
        xlims = get(gca,'xlim');
        ylims = get(gca,'ylim');
        xl2 = get(gca,'xlim');
        yl2 = get(gca,'ylim');
    end
    X = xlims(1):10:xlims(end);
    Y = lincoeff_d(2) + lincoeff_d(1).*X;
    plot(X,Y,'--k','linewidth',1.5);
    text(xlims(2)*.8,ylims(2)*.9,sprintf('r = %.4f\np = %.4f',rho_d,p_rho_d),'fontsize',16);
    set(gca,'xlim',xlims,'ylim',ylims);
    ylabel('Fractional change in firing rate','fontsize',21);
    
    subplot(2,1,2);
    if dep
        mfc = 'k';
        mfc2 = [.7 .7 .7];
    else
        mfc = 'none';
        mfc2 = 'none';
    end
    
    hold on;
    errorbar(sleepdur_bycell_nonan,NREM_all_deltaFR,NREM_all_deltaSEM,'o',...
        'color',[.7 .7 .7],'markerfacecolor',mfc2,'linestyle','none',...
        'linewidth',1.5,'markersize',5);
    errorbar(group_means_offT_d,group_means_d,group_sem_d,'o','capsize',7,'markersize',7,...
        'linestyle','none','markerfacecolor',mfc,'linewidth',1.5,'color','k');
     X = xlims(1):10:xlims(end);
    Y = lincoeff_d(2) + lincoeff_d(1).*X;
    plot(X,Y,'--k','linewidth',1.5);
    plot(xlims,[1 1],'m:','linewidth',2);
    set(gca,'ylim',yl2,'xlim',xl2);
    
    
    ylabel('Fractional change in firing rate','fontsize',21);
    xlabel('Duration of extended sleep (s)','fontsize',20);
    
    
    
    if save_delta_data
        %         delta_data
        SLEEP_DELTAS.group_value                    = group_means_d;
        SLEEP_DELTAS.group_times                    = group_means_offT_d;
        SLEEP_DELTAS.params.G_bin                   = G_bin;
        SLEEP_DELTAS.params.dur_frac                = dur_frac;
        SLEEP_DELTAS.params.duration_threshold      = duration_threshold;
        SLEEP_DELTAS.params.mean_t                  = mean_t;
        SLEEP_DELTAS.params.ext_sleep_time_thresh   = ext_sleep_time_thresh;
        
        
        sf = 'SLEEP_DELTAS_v0.mat';
        sd = 'Z:\ATP_MAIN\CODE\lfp_analysis_beta\Analysis_Data';
        save([sd filesep sf],'SLEEP_DELTAS','-v7.3');
    end
end