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

%% Loop through animals and find extended wake periods

for aa = 1:n_anims
    animal = anim_fields{aa};
    
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
    
    % FIND EXTENDED WAKE
    ext_wakes = find_extended_wake(statetimes,ext_wake_time_thresh);
    
    all_extwakes{aa} = ext_wakes;
    
    expstart_raw = unixtime(statetimes(1,2));
    expstart_unix = [expstart_raw(1:3) 7 30 0];
    expstart_t = unixtime(expstart_unix);
    
    all_expstart{aa} = expstart_t;
    
    clear  ext_sleeps expstart* statetimes
end

all_AW_zFR = [];
all_QW_zFR = [];

%% Loop thru animals for FR analysis
% this is the same as in the extsleep_analysis script
for aa = 1:n_anims
    animal = anim_fields{aa};
    fprintf('Animal %s.\n',animal);
    
    if any(strcmp(anim_list,animal))
        
        anim_rsu_idx = find(strcmp(cell_anims,animal));
        
        anim_cells = MASTER(anim_rsu_idx);
        rsu = get_RSUs(anim_cells,qthresh,dep,perc_thresh,norm,bl_on,dataset,...
            plot_cellsep,negpos,tailslope,0);
        n_RSU = size(rsu,2);
        
        anim_idx = find(strcmp(anim_fields,animal));
        
        % RETRIEVE STATE DATA
        
        ext_wakes = all_extwakes{anim_idx};
        expstart_t = all_expstart{anim_idx};
        
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
        
        daystart = DAYSTART.(animal);
        
        n_wakes = size(ext_wakes,1);
        
        all_wake_dur{anim_idx} = nan(n_wakes, n_RSU);
        deltaFR_AW{anim_idx} = nan(n_wakes,n_RSU);
        deltaFR_QW{anim_idx} = nan(n_wakes,n_RSU);
        
        
        
        for tt = 1:n_wakes
            
            fprintf('Extended wake epoch %u of %u.\n',tt,n_wakes);
            
            
            st_wake = statetimes(ext_wakes(tt,1):ext_wakes(tt,2)+1,:);
            
            st_wake(:,2) = st_wake(:,2) - expstart_t + daystart*24*3600;
            
            try
                wake_durs = diff(st_wake(:,2));
                total_dur = sum(wake_durs);
                
                aw_idx = find(st_wake(:,1) == 4);
                if ~isempty(aw_idx)
                    if aw_idx(end) > size(wake_durs,1)
                        aw_idx(end) = [];
                    end
                    aw_dur = sum(wake_durs(aw_idx));
                    aw_percent = 100 * (aw_dur/total_dur);
                end
                
                qw_idx = find(st_wake(:,1) == 5);
                if ~isempty(qw_idx)
                    if qw_idx(end) > size(wake_durs,1)
                        qw_idx(end) = [];
                    end
                    qw_dur = sum(wake_durs(qw_idx));
                    qwpercent = 100 * (qw_dur/total_dur);
                end
            catch
                keyboard;
            end
            
            wake_start = st_wake(1,2);
            wake_end = st_wake(end,2);
            ext_wake_duration = wake_end - wake_start;
            
            
            
            if wake_start >= homeo_start && wake_end <= homeo_end
                
                fprintf('good\n');
                AW_idx = find(st_wake(:,1) == 4);
                n_AW = numel(AW_idx);
                QW_idx = find(st_wake(:,1) == 5);
                n_QW = numel(QW_idx);
                
                if n_AW > 0
                    idx_AW_first = AW_idx(1);
                    idx_AW_last = AW_idx(end);
                end
                
                idx_QW_first = QW_idx(1);
                idx_QW_last = QW_idx(end);
                
                for cc = 1:n_RSU
                    if size(rsu(cc).onTime,1) > 1
                        fprintf('Multiple on/off times.\n');
                        keyboard;
                    end
                    
                    spikes = rsu(cc).time;
                    switch mean_t
                        case 0
                            do_z = 0;
                            delta_frac = 0;
                        case 1
                            mean_long_t0 = wake_start;
                            mean_long_t1 = wake_end;
                            do_z = 1;
                            delta_frac = 0;
                        case 2
                            mean_long_t0 = rsu(cc).onTime(1);
                            mean_long_t1 = rsu(cc).offTime(end);
                            do_z = 1;
                            delta_frac = 0;
                        case 3
                            mean_long_t0 = wake_start;
                            mean_long_t1 = wake_end;
                            do_z = 0;
                            delta_frac = 1;
                    end
                    
                    rate_tmp = histc(spikes,t_start:G_bin:t_stop) ./ G_bin;
                    
                    tlist = [rsu(cc).onTime rsu(cc).offTime];
                    rate_onoff = processOnOffTimes_ATP(rate_tmp,tlist,G_bin,0,'recov',0);
                    
                    rate_long = rate_onoff(round(mean_long_t0/G_bin):round(mean_long_t1/G_bin));
                    
                    frmean_long = nanmean(rate_long);
                    frstd_long  = std(rate_long,0,'omitnan');
                    
                    % cycle thru NREM epochs and find mean FR in each
                    if n_AW > 0
                        for aw = 1:n_AW
                            
                            AW_t0 = st_wake(AW_idx(aw),2);
                            AW_t1 = st_wake(AW_idx(aw)+1,2);
                            if AW_t0 >=  rsu(cc).onTime(1) && AW_t1 <= rsu(cc).offTime(end) && (AW_t1-AW_t0 >= duration_threshold)
                                
                                AW_idx0 = floor( (AW_t0-mean_long_t0) / G_bin );
                                if AW_idx0 == 0, AW_idx0 = 1; end
                                AW_idx1 = ceil( (AW_t1-mean_long_t0) / G_bin );
                                
                                AW_FR = nanmean(rate_long(AW_idx0:AW_idx1));
                                if do_z
                                    AW_zFR = (AW_FR - frmean_long) / frstd_long;
                                else
                                    AW_zFR = AW_FR;
                                end
                            else
                                AW_zFR = NaN;
                            end
                            
                            all_AW_zFR{anim_idx,tt}(cc,aw) = AW_zFR;
                            all_AW_offT{anim_idx,tt}(cc,aw) = AW_t0 - wake_start;
                        end
                    end
                    % cycle thru REM epochs and find mean FR in each
                    for qw = 1:n_QW
                        
                        QW_t0 = st_wake(QW_idx(qw),2);
                        QW_t1 = st_wake(QW_idx(qw)+1,2);
                        if QW_t0 >=  rsu(cc).onTime(1) && QW_t1 <= rsu(cc).offTime(end) && (QW_t1-QW_t0 >= duration_threshold)
                            
                            QW_idx0 = floor( (QW_t0-mean_long_t0) / G_bin );
                            if QW_idx0 == 0, QW_idx0 = 1; end
                            QW_idx1 = ceil( (QW_t1-mean_long_t0) / G_bin );
                            
                            QW_FR = nanmean(rate_long(QW_idx0:QW_idx1));
                            if do_z
                                QW_zFR = (QW_FR - frmean_long) / frstd_long;
                            else
                                QW_zFR = QW_FR;
                            end
                        else
                            QW_zFR = NaN;
                        end
                        
                        all_QW_zFR{anim_idx,tt}(cc,qw) = QW_zFR;
                        all_QW_offT{anim_idx,tt}(cc,qw) = QW_t0 - wake_start;
                    end
                    
                    if n_AW > 0
                        AW_first_t0 = st_wake(idx_AW_first,2);
                        AW_first_t1 = st_wake(idx_AW_first+1,2);
                        AW_last_t0 = st_wake(idx_AW_last,2);
                        AW_last_t1 = st_wake(idx_AW_last+1,2);
                        
                        AW_dur_first = AW_first_t1 - AW_first_t0;
                        AW_dur_last = AW_last_t1 - AW_last_t0;
                        
                        
                        if AW_first_t0 >=  rsu(cc).onTime(1) && AW_last_t1 <= rsu(cc).offTime(end) && ...
                                AW_dur_first >= duration_threshold && AW_dur_last >= duration_threshold
                                                    
                            AW_first_idx0 = floor( (AW_first_t0-mean_long_t0) / G_bin );
                            if AW_first_idx0 == 0, AW_first_idx0 = 1; end
                            AW_first_idx1 = ceil( (AW_first_t1-mean_long_t0) / G_bin );
                            
                            AW_last_idx0 = floor( (AW_last_t0-mean_long_t0) / G_bin );
                            AW_last_idx1 = ceil( (AW_last_t1-mean_long_t0) / G_bin );
                            
                            AW_first_FR = nanmean(rate_long(AW_first_idx0:AW_first_idx1));
                            AW_last_FR = nanmean(rate_long(AW_last_idx0:AW_last_idx1));
                            
                            if do_z
                                AW_first_zFR = (AW_first_FR - frmean_long) / frstd_long;
                                AW_last_zFR = (AW_last_FR - frmean_long) / frstd_long;
                            else
                                AW_first_zFR = AW_first_FR;
                                AW_last_zFR = AW_last_FR;
                            end
                            
                            if delta_frac
                                AW_delta_zFR = AW_last_zFR / (AW_first_zFR);
                                if AW_delta_zFR == Inf || AW_delta_zFR == 0
                                    AW_delta_zFR = NaN;
                                end
                            else
                                AW_delta_zFR = AW_last_zFR - AW_first_zFR;
                            end
                        else
                            AW_delta_zFR = NaN;
                        end
                        deltaFR_AW{anim_idx}(tt,cc) = AW_delta_zFR;
                    end
                    all_wake_dur{anim_idx}(tt,cc) = ext_wake_duration;
                    
                    % find difference in FR between last NREM and first NREM
                    if mean_t ~= 3
                        QW_first_t0 = st_wake(idx_QW_first,2);
                        QW_first_t1 = st_wake(idx_QW_first+1,2);
                        QW_last_t0 = st_wake(idx_QW_last,2);
                        QW_last_t1 = st_wake(idx_QW_last+1,2);

                        QW_dur_first = QW_first_t1 - QW_first_t0;
                        QW_dur_last = QW_last_t1 - QW_last_t0;
                        
                    else
                        wake_dur_frac = dur_frac * ext_wake_duration;
                        QW_first_t0 = st_wake(1,2);
                        QW_first_t1 = QW_first_t0 + wake_dur_frac;
                        QW_last_t1 = st_wake(end,2);
                        QW_last_t0 = QW_last_t1 - wake_dur_frac;
                    end
                    
                    if QW_first_t0 >=  rsu(cc).onTime(1) && QW_last_t1 <= rsu(cc).offTime(end) && ...
                            QW_dur_first >= duration_threshold && QW_dur_last >= duration_threshold
                        
                        QW_first_idx0 = floor( (QW_first_t0-mean_long_t0) / G_bin );
                        if QW_first_idx0 == 0, QW_first_idx0 = 1; end
                        QW_first_idx1 = ceil( (QW_first_t1-mean_long_t0) / G_bin );
                        
                        QW_last_idx0 = floor( (QW_last_t0-mean_long_t0) / G_bin );
                        QW_last_idx1 = ceil( (QW_last_t1-mean_long_t0) / G_bin );
                        
                        QW_first_FR = nanmean(rate_long(QW_first_idx0:QW_first_idx1));
                        QW_last_FR = nanmean(rate_long(QW_last_idx0:QW_last_idx1));
                        
                        if do_z
                            QW_first_zFR = (QW_first_FR - frmean_long) / frstd_long;
                            QW_last_zFR = (QW_last_FR - frmean_long) / frstd_long;
                        else
                            QW_first_zFR = QW_first_FR;
                            QW_last_zFR = QW_last_FR;
                        end
                        
                        if delta_frac
                            QW_delta_zFR = QW_last_zFR / (QW_first_zFR);
                            if QW_delta_zFR == Inf || QW_delta_zFR == 0
                                QW_delta_zFR = NaN;
                            end
                        else
                            QW_delta_zFR = QW_last_zFR - QW_first_zFR;
                        end
                    else
                        QW_delta_zFR = NaN;
                    end
                    deltaFR_QW{anim_idx}(tt,cc) = QW_delta_zFR;
                end
                
                
            else
                fprintf('bad\n');
                all_AW_zFR{anim_idx,tt} = NaN;
                all_QW_zFR{anim_idx,tt} = NaN;
                all_AW_offT{anim_idx,tt} = NaN;
                all_QW_offT{anim_idx,tt} = NaN;
                
            end
        end
    end
end

% keyboard;
%% deltaFR compile data
QW_all_deltas = cellfun(@(x) x(:),deltaFR_QW,'UniformOutput',false);
cat_QW_deltas = cat(1,QW_all_deltas{:});
cat_QW_deltas(cat_QW_deltas > 5) = NaN;
delta_nonans = ~isnan(cat_QW_deltas);
QW_delta_nonan = cat_QW_deltas(delta_nonans);

AW_all_deltas = cellfun(@(x) x(:),deltaFR_AW,'UniformOutput',false);
cat_AW_deltas = cat(1,AW_all_deltas{:});
cat_AW_deltas(cat_AW_deltas > 5) = NaN;
delta_nonans = ~isnan(cat_AW_deltas);
AW_delta_nonan = cat_AW_deltas(delta_nonans);

all_wake_dur_cellcat = cellfun(@(x) x(:),all_wake_dur,'UniformOutput',false);
cat_wake_dur = cat(1,all_wake_dur_cellcat{:});
wakedur_nonan = cat_wake_dur(delta_nonans);


wakedur_bycell = cellfun(@(x) nanmean(x,2),all_wake_dur,'UniformOutput',false);
wakedur_bycell_cat = cat(1,wakedur_bycell{:});

AW_cellcat = cellfun(@(x) nanmean(x,2),deltaFR_AW,'UniformOutput',false);
AW_cellcat_sem = cellfun(@(x) nanstd(x,[],2) ./ sqrt(sum(~isnan(x),2)),deltaFR_AW,'UniformOutput',false);
AW_cellmean = cat(1,AW_cellcat{:});
AW_cellsem = cat(1,AW_cellcat_sem{:});
AW_all_deltaFR = AW_cellmean(~isnan(AW_cellmean));
AW_all_deltaSEM = AW_cellsem(~isnan(AW_cellmean));
n_epochs_AW = sum(~isnan(AW_cellmean));

QW_cellcat = cellfun(@(x) nanmean(x,2),deltaFR_QW,'UniformOutput',false);
QW_cellcat_sem = cellfun(@(x) nanstd(x,[],2) ./ sqrt(sum(~isnan(x),2)),deltaFR_QW,'UniformOutput',false);
QW_cellmean = cat(1,QW_cellcat{:});
QW_cellsem = cat(1,QW_cellcat_sem{:});
QW_all_deltaFR = QW_cellmean(~isnan(QW_cellmean));
QW_all_deltaSEM = QW_cellsem(~isnan(QW_cellmean));
wakedur_bycell_nonan = wakedur_bycell_cat(~isnan(QW_cellmean));
n_epochs_QW = sum(~isnan(QW_cellmean));

AW_mean_c = nanmean(AW_cellmean);
AW_sem_c = std(AW_cellmean,0,'omitnan') / sqrt(n_epochs_AW - 1);
AW_mean = nanmean(AW_delta_nonan);
AW_sem = std(AW_delta_nonan,0,'omitnan') / sqrt(sum(~isnan(AW_delta_nonan)) - 1);

QW_mean_c = nanmean(QW_cellmean);
QW_sem_c = std(QW_cellmean,0,'omitnan') / sqrt(n_epochs_QW - 1);
QW_mean = nanmean(QW_delta_nonan);
QW_sem = std(QW_delta_nonan,0,'omitnan') / sqrt(sum(~isnan(QW_delta_nonan)) - 1);

%% Fig 5G - AW FR vs time correlation
cellcat_AW_offT = cellfun(@(x) x(:),all_AW_offT,'UniformOutput',false);
cellcat_AW_zFR = cellfun(@(x) x(:),all_AW_zFR,'UniformOutput',false);

cat_AW_offT = cat(1,cellcat_AW_offT{:});
cat_AW_zFR = cat(1,cellcat_AW_zFR{:});

fr_notnans = ~isnan(cat_AW_zFR);

AW_offT_nonan = cat_AW_offT(fr_notnans);
AW_zFR_nonan = cat_AW_zFR(fr_notnans);

all_data = [AW_zFR_nonan AW_offT_nonan];
all_data_sort = sortrows(all_data,2);

ngroups = 10;
group_split = round(linspace(0,numel(AW_offT_nonan),ngroups+1));
for uu = 1:ngroups
    this_group = all_data_sort(group_split(uu)+1 : group_split(uu+1), 1);
    group_means(uu) = nanmean(this_group);
    group_sem(uu) = std(this_group,0,'omitnan') / sqrt(numel(this_group)-1);
    group_means_offT(uu) = nanmean(all_data_sort(group_split(uu)+1 : group_split(uu+1), 2));
end



[aw_rho,aw_p_rho] = corr(AW_zFR_nonan, AW_offT_nonan);
aw_asterisks = get_asterisks_from_pval(aw_p_rho);
aw_lincoeff = polyfit(AW_offT_nonan, AW_zFR_nonan, 1);


if mean_t ~= 3
    %% plotting
    markersize = 4;
    mylims = 1;
    
    
    setFigureDefaults;
    figure();
    set(gcf,'position',[.2 .1 .3 .8]);
    subplot(2,1,1);
    if dep
        scatter(AW_offT_nonan,AW_zFR_nonan,markersize^2,c_aw,'s','filled');
    else
        scatter(AW_offT_nonan,AW_zFR_nonan,markersize^2,c_aw,'s');
    end
    hold on;
    if mylims
        xlims = [-500 11000];
        ylims = [-2 2];
        xl2 = [-250 10000];
        yl2 = [-0.2 0.4];
    else
        xlims = get(gca,'xlim');
        ylims = get(gca,'ylim');
        xl2 = get(gca,'xlim');
        yl2 = get(gca,'ylim');
    end
    X = xlims(1):10:xlims(end);
    Y = aw_lincoeff(2) + aw_lincoeff(1).*X;
    plot(X,Y,'--k','linewidth',1.5);
    set(gca,'xlim',xlims,'ylim',ylims);
    ylabel('Firing rate (z)','fontsize',21);
    
    
    subplot(2,1,2);
    if dep
        mfc = c_aw;
    else
        mfc = 'none';
    end
    errorbar(group_means_offT,group_means,group_sem,'s','capsize',12,'markersize',8,...
        'linestyle','none','markerfacecolor',mfc,'linewidth',1.5,'color',c_aw);
    hold on;
    box off
    plot(X,Y,'--k','linewidth',1.5);
    if mylims
        set(gca,'xlim',xl2,'ylim',yl2);
        text(xl2(2)*.8,yl2(2)*.9,sprintf('r = %.4f\np = %.4f',aw_rho,aw_p_rho),'fontsize',16);
    end
    ylabel('Firing rate (z)','fontsize',21);
    xlabel('Time from start of extended wake (s)','fontsize',20);
    
    %% QW
    cellcat_QW_offT = cellfun(@(x) x(:),all_QW_offT,'UniformOutput',false);
    cellcat_QW_zFR = cellfun(@(x) x(:),all_QW_zFR,'UniformOutput',false);
    
    cat_QW_offT = cat(1,cellcat_QW_offT{:});
    cat_QW_zFR = cat(1,cellcat_QW_zFR{:});
    
    QW_fr_notnans = ~isnan(cat_QW_zFR);
    
    QW_offT_nonan = cat_QW_offT(QW_fr_notnans);
    QW_zFR_nonan = cat_QW_zFR(QW_fr_notnans);
    
    QW_all_data = [QW_zFR_nonan QW_offT_nonan];
    QW_all_data_sort = sortrows(QW_all_data,2);
    
    % ngroups = 5;
    QW_group_split = round(linspace(0,numel(QW_offT_nonan),ngroups+1));
    for uu = 1:ngroups
        QW_this_group = QW_all_data_sort(QW_group_split(uu)+1 : QW_group_split(uu+1), 1);
        QW_group_means(uu) = nanmean(QW_this_group);
        QW_group_sem(uu) = std(QW_this_group,0,'omitnan') / sqrt(numel(QW_this_group)-1);
        QW_group_means_offT(uu) = nanmean(QW_all_data_sort(QW_group_split(uu)+1 : QW_group_split(uu+1), 2));
    end
    
    
    
    [qw_rho,qw_p_rho] = corr(QW_offT_nonan, QW_zFR_nonan);
    qw_asterisks = get_asterisks_from_pval(qw_p_rho);
    qw_lincoeff = polyfit(QW_offT_nonan, QW_zFR_nonan, 1);
    
    %% plotting
    markersize = 4;
    mylims = 1;
    
    
    setFigureDefaults;
    figure();
    set(gcf,'position',[.2 .1 .3 .8]);
    subplot(2,1,1);
    if dep
        scatter(QW_offT_nonan,QW_zFR_nonan,markersize^2,c_qw,'s','filled');
    else
        scatter(QW_offT_nonan,QW_zFR_nonan,markersize^2,c_qw,'s');
    end
    hold on;
    if mylims
        xlims = [-500 11000];
        ylims = [-2 2];
        xl2 = [-250 10000];
        yl2 = [-0.2 0.4];
    else
        xlims = get(gca,'xlim');
        ylims = get(gca,'ylim');
        xl2 = get(gca,'xlim');
        yl2 = get(gca,'ylim');
    end
    X = xlims(1):10:xlims(end);
    Y = qw_lincoeff(2) + qw_lincoeff(1).*X;
    plot(X,Y,'--k','linewidth',1.5);
    set(gca,'xlim',xlims,'ylim',ylims);
    ylabel('Firing rate (z)','fontsize',21);
    
    if dep
        mfc = c_qw;
    else
        mfc = 'none';
    end
    subplot(2,1,2);
    errorbar(QW_group_means_offT,QW_group_means,QW_group_sem,'s','capsize',12,'markersize',8,...
        'linestyle','none','markerfacecolor',mfc,'color',c_qw,'linewidth',1.5);
    hold on;
    box off;
    plot(X,Y,'--k','linewidth',1.5);
    if mylims
        set(gca,'xlim',xl2,'ylim',yl2);
        text(xl2(2)*.8,yl2(2)*.9,sprintf('r = %.4f\np = %.4f',qw_rho,qw_p_rho),'fontsize',16);
    end
    ylabel('Firing rate (z)','fontsize',21);
    xlabel('Time from start of extended wake (s)','fontsize',20);
    
    %% Fig 5H - plotting the first-last change
    mylims = 1;
    
    [~,p_aw] = ttest(AW_delta_nonan);
    [a_aw,fsz_aw] = get_asterisks_from_pval(p_aw);
    [~,p_qw] = ttest(QW_delta_nonan);
    [a_qw,fsz_qw] = get_asterisks_from_pval(p_qw);
    
    dfig = figure();
    set(dfig,'position',[.1 .2 .25 .4]);
    box off
    hold on;
    bw = .3;
    csz = 20;
    if dep
        mfc1 = c_aw;
        mfc2 = c_qw;
        ec1 = 'none';
        ec2 = 'none';
    else
        mfc1 = 'none';
        mfc2 = 'none';
        ec1 = c_aw;
        ec2 = c_qw;
    end
    
    aw_bar = bar(1,AW_mean,bw,'edgecolor',ec1,'facecolor',mfc1,'linewidth',2);
    aw_err = errorbar(1,AW_mean,AW_sem,'linestyle','none','capsize',csz,...
        'color',c_aw,'linewidth',2);
    qw_bar = bar(2,QW_mean,bw,'edgecolor',ec2,'facecolor',mfc2,'linewidth',2);
    qw_err = errorbar(2,QW_mean,QW_sem,'linestyle','none','capsize',csz,...
        'color',c_qw,'linewidth',2);
    if mylims
        xl = [.5 2.5];
        yl = [-0.2 0.2];
        yt = -0.2:.1:0.2;
    else
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        yt = get(gca,'ytick');
    end
    text(.75,0.2,sprintf('p = %.4f',p_aw),'fontsize',16);
    text(1.75,0.2,sprintf('p = %.4f',p_qw),'fontsize',16);
    set(gca,'xlim',xl,'ylim',yl,'xtick',[1 2],'xticklabel',{'Active','Quiet'},...
        'ytick',yt);
    ylabel('Firing rate (z)','fontsize',22);
    
    
else
    %%
    markersize = 5;
    
    all_data_d = [QW_delta_nonan wakedur_nonan];
    all_data_sort_d = sortrows(all_data_d,2);
    
    ngroups = 10;
    group_split_d = round(linspace(0,numel(QW_delta_nonan),ngroups+1));
    clear group_means_d group_sem_d group_means_offT_d group_max_d
    for uu = 1:ngroups
        this_group = all_data_sort_d(group_split_d(uu)+1 : group_split_d(uu+1), 1);
        group_means_d(uu) = nanmean(this_group);
        group_sem_d(uu) = std(this_group,0,'omitnan') / sqrt(numel(this_group)-1);
        group_means_offT_d(uu) = nanmean(all_data_sort_d(group_split_d(uu)+1 : group_split_d(uu+1), 2));
        group_max_d(uu) = max(all_data_sort_d(group_split_d(uu)+1 : group_split_d(uu+1), 2));
    end

    
    [rho_d,p_rho_d] = corr(wakedur_nonan, QW_delta_nonan);
%     asterisks_d = get_asterisks_from_pval(p_rho_d);
    lincoeff_d = polyfit(wakedur_nonan, QW_delta_nonan, 1);
    
    
    mylims = 0;
    
    setFigureDefaults;
    figure();
    set(gcf,'position',[.1 .1 .6 .7]);
    subplot(2,1,1);
    if dep
        scatter(wakedur_nonan,QW_delta_nonan,markersize^2,'k','filled');
    else
        scatter(wakedur_nonan,QW_delta_nonan,markersize^2,'k');
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
    errorbar(wakedur_bycell_nonan,QW_all_deltaFR,QW_all_deltaFR,'o',...
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
    xlabel('Duration of extended wake (s)','fontsize',20);
    
    
    
    if save_delta_data
        %         delta_data
        WAKE_DELTAS.group_value                    = group_means_d;
        WAKE_DELTAS.group_times                    = group_means_offT_d;
        WAKE_DELTAS.params.G_bin                   = G_bin;
        WAKE_DELTAS.params.dur_frac                = dur_frac;
        WAKE_DELTAS.params.duration_threshold      = duration_threshold;
        WAKE_DELTAS.params.mean_t                  = mean_t;
        WAKE_DELTAS.params.ext_wake_time_thresh    = ext_wake_time_thresh;
        
        
        sf = 'WAKE_DELTAS_v0.mat';
        sd = 'Z:\ATP_MAIN\CODE\lfp_analysis_beta\Analysis_Data';
        save([sd filesep sf],'WAKE_DELTAS','-v7.3');
    end
    
    
    
end





