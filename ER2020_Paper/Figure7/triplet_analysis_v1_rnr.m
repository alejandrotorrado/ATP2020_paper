clearvars -except CONTCELL*
close all
clc

STATES = CONTCELL_recov.STATETIMES;
DAYSTART = STATES.DAYSTART;
anim_fields = fieldnames(STATES);
anim_fields(strcmp(anim_fields,'DAYSTART')) = [];
n_anims = numel(anim_fields); % minus daystart field

cell_anims = {CONTCELL_recov.MASTER.animal};


if ispc
    loadfile = 'Z:\ATP_MAIN\DATA\Eye_Reopening\recov_DATA\recov_analysis_Jul19.mat';
elseif ismac
end

recov_load = load(loadfile);
recov = recov_load.recov_analysis;

%% analysis params
save_the_data = 1;
ID_NUM = 7;

homeo_period_start = 8.0;
homeo_period_end   = 10.0;

t_start = 0;
t_stop = 12*24*3600;


dep = 1;
qthresh = 2;
perc_thresh = 0.75;
norm = 0;
bl_on = 1;
dataset = 'recov';
plot_cellsep = 0;
negpos = 0.39;
tailslope = 0.005;
big_val_thresh = 0;
small_val_thresh = 0;
toobig_ = 0;
toosmall_ = 0;

G_bin = 5; % seconds. Bin size for FR calculations.

bl_t0 = 6.5*24*3600;
bl_t1 = 7.0*24*3600;

bl_bin0 = bl_t0 / G_bin;
bl_bin1 = bl_t1 / G_bin;

use_periods = 0;
do_z = 1;

% options of mean_t:
%   0 for all time
%   1 to z-score only to triplet
%   2 to z-score to 12 hour mean
%   10 for no z-scoring (will atomatically set do_z = 0)
mean_t = 1;


duration_threshold = 10; % seconds
analysis_dur_thresh = 50; % seconds

do_spindles = 0;
do_SWA = 0;
do_power = 0;
do_thirds = 1;
do_packets = 1;
min_end_excl = 1 * 60; % seconds (minutes to exclude at state transitions)
min_fr_len = 5 * 60; % seconds (minutes to calc FR before transition)
min_end_excl_REM = 0.5 * 60; % seconds (minutes to exclude at state transitions)
min_fr_len_REM = 1.5 * 60; % seconds (minutes to calc FR before transition)
REM_dur_thresh = 3*60;


if dep
    anim_list = {'AT12','AT14','AT16','AT25','AT27','AT29','KH67'};
else
    anim_list = {'AT12','AT14','AT16','AT25','AT27','AT29','KH67','KH72','KH73','KH75'};
end

% animal order:
% AT12, KH67, AT14, AT16, AT29, KH73, KH75, AT27
rec_periods = {[8.0 9.0],[8.5 10],[8.0 9.0],[8.5 10],[8.0 9.5],...
    [8.5 10],[8.5 10],[8.5 10]};

use_idx = 1;
if dep
    RSU_idx = recov.DEPRIVED.RSU_idx;
else
    RSU_idx = recov.CONTROL.RSU_idx;
end


%% LOOP THRU AVAILABLE STATETIMES

for aa = 1:n_anims
    animal = anim_fields{aa};
    
    statetimes = STATES.(animal);
    
    % remove repeats
    st_diff = diff(statetimes(:,1));
    kill_these = find(st_diff==0);
    statetimes(kill_these+1,:) = [];
    
    % remove short states
    st_timediff = diff(statetimes(:,2));
    too_short = find(st_timediff <= duration_threshold);
    statetimes(too_short,:) = [];
    
    % remove repeats again
    st_diff = diff(statetimes(:,1));
    kill_these = find(st_diff==0);
    statetimes(kill_these+1,:) = [];
    
    % FIND NREM/REM/NREM TRIPLETS
    
    triplets = find_sleep_triplets_2(statetimes,'rnr');
    %     ext_sleeps = find_extended_sleep(statetimes);
    
    all_anim_cells = {CONTCELL_recov.MASTER.animal};
    this_anim_cells = find(strcmp(all_anim_cells,animal));
    expt_start_temp = CONTCELL_recov.MASTER(this_anim_cells(1)).EXPTSTART;
    
    expstart_raw = unixtime(expt_start_temp);
    expstart_unix = [expstart_raw(1:3) 7 30 0];
    expstart_t = unixtime(expstart_unix);
    
    all_triplets{aa} = triplets;
    %     all_extsleeps{aa} = ext_sleeps;
    all_expstart{aa} = expstart_t;
    
    
    clear triplets ext_sleeps expstart* statetimes
    
end

badcount = 0;



%% MAIN LOOP THRU CELLS

for aa = 1:n_anims
    animal = anim_fields{aa};
    fprintf('Animal %s.\n',animal);
    
    if any(strcmp(anim_list,animal))
        
        anim_rsu_idx = find(strcmp(cell_anims,animal));
        
        if use_idx
            this_rsu_idx = intersect(RSU_idx,anim_rsu_idx);
            rsu = CONTCELL_recov.MASTER(this_rsu_idx);
        else
            anim_cells = CONTCELL_recov.MASTER(anim_rsu_idx);
            rsu = get_RSUs(anim_cells,qthresh,dep,perc_thresh,norm,bl_on,dataset,...
                plot_cellsep,negpos,tailslope,0);
        end
%         keyboard;
        n_RSU = size(rsu,2);
        
        anim_idx = find(strcmp(anim_fields,animal));
        
        % RETRIEVE STATE DATA
        
        triplets = all_triplets{anim_idx};
        %         ext_sleeps = all_extsleeps{anim_idx};
        expstart_t = all_expstart{anim_idx};
        
        statetimes = STATES.(animal);
        
        % remove repeats
        st_diff = diff(statetimes(:,1));
        kill_these = find(st_diff==0);
        statetimes(kill_these+1,:) = [];
        
        % remove short states
        st_timediff = diff(statetimes(:,2));
        too_short = find(st_timediff <= duration_threshold);
        statetimes(too_short,:) = [];
        
        % remove repeats again
        st_diff = diff(statetimes(:,1));
        kill_these = find(st_diff==0);
        statetimes(kill_these+1,:) = [];
        
        daystart = DAYSTART.(animal);
        
        n_triplets = size(triplets,1);
        %         n_sleeps = size(ext_sleeps,1);
        
        if use_periods
            homeo_start = rec_periods{aa}(1)*24*3600;
            homeo_end = rec_periods{aa}(2)*24*3600;
        else
            homeo_start = homeo_period_start*24*3600;
            homeo_end = homeo_period_end*24*3600;
        end
        
        fprintf('Homeo start: %.1f | Homeo end: %.1f.\n',homeo_start/(24*3600),homeo_end/(24*3600));
        
        if n_RSU > 0 && do_power
            %% load lfp for this animal
            if ispc
                lfp_dir = 'Z:\ANIMALDATA\MLS_DATA\LFP_DATA';
            elseif ismac
                lfp_dir = '/Volumes/turrigiano-lab/ANIMALDATA/MLS_DATA/LFP_DATA';
            end
            fprintf('Loading LFP for this animal.\n');
            clear LFPinfo
            tic
            lfp_load = load([lfp_dir filesep animal filesep animal '_LFPinfo.mat']);
            LFPinfo = lfp_load.LFPinfo;
            clear lfp_load
            toc
        else
            fprintf('Skip LFP load since there are no cells for this animal.\n');
        end
        
        %% loop thru triplets
        delta_FR{anim_idx} = nan(n_triplets,n_RSU);
        allREM_dur{anim_idx} = nan(n_triplets,1);
        allNR1_dur{anim_idx} = nan(n_triplets,1);
        allNR2_dur{anim_idx} = nan(n_triplets,1);
        allREM_theta{anim_idx} = nan(n_triplets,1);
        allNR1_delta{anim_idx} = nan(n_triplets,1);
        allTD_ratio{anim_idx} = nan(n_triplets,1);
        NR1_thirds_allFR{anim_idx} = nan(n_triplets,n_RSU,3);
        NR2_thirds_allFR{anim_idx} = nan(n_triplets,n_RSU,3);
        REM_thirds_allFR{anim_idx} = nan(n_triplets,n_RSU,3);
        NR1_packets_allFR{anim_idx} = nan(n_triplets,n_RSU);
        NR2_packets_allFR{anim_idx} = nan(n_triplets,n_RSU);
        NR1_packets_allFR_z{anim_idx} = nan(n_triplets,n_RSU);
        NR2_packets_allFR_z{anim_idx} = nan(n_triplets,n_RSU);
        
        NR1_p1_allFR{anim_idx} = nan(n_triplets,n_RSU);
        NR1_p2_allFR{anim_idx} = nan(n_triplets,n_RSU);
        NR2_p1_allFR{anim_idx} = nan(n_triplets,n_RSU);
        NR2_p2_allFR{anim_idx} = nan(n_triplets,n_RSU);
        REM_p1_allFR{anim_idx} = nan(n_triplets,n_RSU);
        REM_p2_allFR{anim_idx} = nan(n_triplets,n_RSU);
        
        if do_spindles
            spindle_incidence{anim_idx} = nan(n_triplets,1);
        end
        if do_SWA
            SWA_amplitude{anim_idx} = num2cell(nan(n_triplets,1));
            if n_RSU > 0
                fprintf('Calculating mean SWA for this animal in this period (for normalization).\n\n');
                [mean_SWA{anim_idx}, ~] = get_mean_SWA(LFPinfo,homeo_start,homeo_end,expstart_t,daystart);
            else
                mean_SWA{anim_idx} = [];
            end
        end
        
        % only do all the stuff if there are cells for this animal
        if n_RSU > 0
            
            for tt = 1:n_triplets
                
                %         fprintf('Triplet %u of %u.\n',tt,n_triplets);
                
                idx_NR1_start = triplets(tt,1);
                idx_NR1_end   = triplets(tt,2);
                
                idx_REM_start = triplets(tt,2);
                idx_REM_end   = triplets(tt,3);
                
                idx_NR2_start = triplets(tt,3);
                idx_NR2_end   = triplets(tt,3) + 1;
                
                NR1_t0 = statetimes(idx_NR1_start,2);
                NR1_t1 = statetimes(idx_NR1_end,2);
                
                REM_t0 = statetimes(idx_REM_start,2);
                REM_t1 = statetimes(idx_REM_end,2);
                
                NR2_t0 = statetimes(idx_NR2_start,2);
                NR2_t1 = statetimes(idx_NR2_end,2);
                
                NR1_t0 = NR1_t0 - expstart_t + daystart*24*3600;
                NR1_t1 = NR1_t1 - expstart_t + daystart*24*3600;
                REM_t0 = REM_t0 - expstart_t + daystart*24*3600;
                REM_t1 = REM_t1 - expstart_t + daystart*24*3600;
                NR2_t0 = NR2_t0 - expstart_t + daystart*24*3600;
                NR2_t1 = NR2_t1 - expstart_t + daystart*24*3600;
                
                REM_dur = REM_t1 - REM_t0;
                NR1_dur = NR1_t1 - NR1_t0;
                NR2_dur = NR2_t1 - NR2_t0;
                
                if do_power
                    [REM_theta_power, NR1_delta_power, REM_td_ratio] = triplet_power(NR1_t0,NR2_t1,NR1_t0,NR1_t1,REM_t0,REM_t1,LFPinfo,expstart_t,daystart);
                end
                
                if do_spindles
                    [spindle_inc, n_spindles] = triplet_spindles(NR1_t0,NR2_t1,NR1_t0,NR1_t1,NR1_dur,LFPinfo,expstart_t,daystart);
                end
                
                if do_SWA
                    NR1_t = [NR1_t0, NR1_t1];
                    [swa_ampl] = triplet_SWA(NR1_t0,NR2_t1,LFPinfo,expstart_t,daystart,'first_state',NR1_t);
                end
                
                if NR1_t0 >= homeo_start && NR2_t1 <= homeo_end && all([REM_dur,NR1_dur,NR2_dur] >= analysis_dur_thresh)
                    %             fprintf('good\n');
                    
                    for cc = 1:n_RSU
                        if size(rsu(cc).onTime,1) > 1
                            fprintf('Multiple on/off times.\n');
%                             keyboard;tlist = [rsu(cc).onTime, rsu(cc).offTime];
                            this_T_idx = find(NR1_t0 >= tlist(:,1),1,'last');
                            this_onT = tlist(this_T_idx,1);
                            this_offT = tlist(this_T_idx,2);
                        else
                            this_onT = rsu(cc).onTime(1);
                            this_offT = rsu(cc).offTime(1);
                        end
                        if NR1_t0 >=  this_onT && NR2_t1 <= this_offT
                            %             fprintf('\n\tFound good on/off times!\n');
                            %             keyboard;
                            
                            % get rsu spikes
                            spikes = rsu(cc).time;
                            %                         keyboard;
                            
                            % calculate 12-hour mean
                            switch mean_t
                                case 0
                                    mean_long_t0 = rsu(cc).onTime(1);
                                    mean_long_t1 = rsu(cc).offTime(end);
                                case 1
                                    mean_long_t0 = NR1_t0;
                                    mean_long_t1 = NR2_t1;
                                case 2
                                    half_point = NR1_t0 + (NR2_t1 - NR1_t0) / 2;
                                    mean_long_t0 = half_point - 6*3600;
                                    mean_long_t1 = half_point + 6*3600;
                                case 10
                                    mean_long_t0 = rsu(cc).onTime(1);
                                    mean_long_t1 = rsu(cc).offTime(end);
                                    do_z = 0;
                            end
                            
                            rate_tmp = histc(spikes,t_start:G_bin:t_stop) ./ G_bin;
%                             keyboard;
                            
                            bl_rate = nanmean(rate_tmp(bl_bin0 : bl_bin1));
                            
                            rate_norm = rate_tmp ./ bl_rate;
                            
                            
                            if ~isnan(rate_tmp)
                                
                                tlist = [rsu(cc).onTime rsu(cc).offTime];
                                rate_onoff = processOnOffTimes_ATP(rate_tmp,tlist,G_bin,0,'recov',0);
                                
                                ix0 = floor(mean_long_t0/G_bin);
                                ix1 = ceil(mean_long_t1/G_bin);
                                %                         keyboard;
                                if ix0 <= 0, ix0 = 1; end
                                if ix1 > size(rate_onoff,1), ix1 = size(rate_onoff,1); end
                                
                                rate_long = rate_onoff(ix0:ix1);
                                
                                if sum(rate_long==0) > 0 %numel(rate_long) / 4
                                    
                                    frmean_long = nanmean(rate_long);
                                    frstd_long  = std(rate_long,0,'omitnan');
                                    
                                    NR1_idx0 = ceil( (NR1_t0-mean_long_t0) / G_bin );
                                    if NR1_idx0 == 0, NR1_idx0 = 1; end
                                    NR1_idx1 = floor( (NR1_t1-mean_long_t0) / G_bin );
                                    
                                    NR2_idx0 = ceil( (NR2_t0-mean_long_t0) / G_bin );
                                    if NR2_idx0 == 0, NR2_idx0 = 1; end
                                    NR2_idx1 = floor( (NR2_t1-mean_long_t0) / G_bin );
                                    
                                    REM_idx0 = ceil( (REM_t0-mean_long_t0) / G_bin );
                                    if REM_idx0 == 0, REM_idx0 = 1; end
                                    REM_idx1 = floor( (REM_t1-mean_long_t0) / G_bin );
                                    
                                    NR1_FR = nanmean(rate_long(NR1_idx0:NR1_idx1));
                                    NR2_FR = nanmean(rate_long(NR2_idx0:NR2_idx1));
                                    
                                    if do_thirds
                                        idx_thirds_NR1 = round(linspace(NR1_idx0,NR1_idx1,4));
                                        idx_thirds_NR2 = round(linspace(NR2_idx0,NR2_idx1,4));
                                        idx_thirds_REM = round(linspace(REM_idx0,REM_idx1,4));
                                        
                                        NR1_T1_FR = nanmean(rate_long(idx_thirds_NR1(1):idx_thirds_NR1(2)));
                                        NR1_T2_FR = nanmean(rate_long(idx_thirds_NR1(2):idx_thirds_NR1(3)));
                                        NR1_T3_FR = nanmean(rate_long(idx_thirds_NR1(3):idx_thirds_NR1(4)));
                                        
                                        NR2_T1_FR = nanmean(rate_long(idx_thirds_NR2(1):idx_thirds_NR2(2)));
                                        NR2_T2_FR = nanmean(rate_long(idx_thirds_NR2(2):idx_thirds_NR2(3)));
                                        NR2_T3_FR = nanmean(rate_long(idx_thirds_NR2(3):idx_thirds_NR2(4)));
                                        
                                        REM_T1_FR = nanmean(rate_long(idx_thirds_REM(1):idx_thirds_REM(2)));
                                        REM_T2_FR = nanmean(rate_long(idx_thirds_REM(2):idx_thirds_REM(3)));
                                        REM_T3_FR = nanmean(rate_long(idx_thirds_REM(3):idx_thirds_REM(4)));
                                    end
                                    
                                    if do_packets
                                        %                                     keyboard;
                                        NR1_packet_idx1 = ceil( (NR1_t1 - min_end_excl - mean_long_t0) / G_bin);
                                        NR1_packet_idx0 = floor( (NR1_t1 - min_end_excl - min_fr_len - mean_long_t0) / G_bin);
                                        if NR1_packet_idx0 < NR1_idx0
                                            NR1_packet_idx0 = NR1_idx0;
                                        end
                                        
                                        NR2_packet_idx0 = ceil( (NR2_t0 + min_end_excl - mean_long_t0) / G_bin);
                                        NR2_packet_idx1 = floor( (NR2_t0 + min_end_excl + min_fr_len - mean_long_t0) / G_bin);
                                        if NR2_packet_idx1 > NR2_idx1
                                            NR2_packet_idx1 = NR2_idx1;
                                        end
                                        
                                        NR1_packet_FR = nanmean(rate_long(NR1_packet_idx0:NR1_packet_idx1));
                                        NR2_packet_FR = nanmean(rate_long(NR2_packet_idx0:NR2_packet_idx1));
                                        
                                        % packets by state
                                        % REM1
                                        if NR1_dur >= REM_dur_thresh
                                            NR1_packet_1_idx0 = ceil( (NR1_t0 + min_end_excl_REM - mean_long_t0) / G_bin);
                                            NR1_packet_1_idx1 = floor( (NR1_t0 + min_end_excl_REM + min_fr_len_REM - mean_long_t0) / G_bin);
                                            NR1_packet_2_idx1 = ceil( (NR1_t1 - min_end_excl_REM - mean_long_t0) / G_bin);
                                            NR1_packet_2_idx0 = floor( (NR1_t1 - min_end_excl_REM - min_fr_len_REM - mean_long_t0) / G_bin);
                                            NR1_packet_2 = nanmean(rate_long(NR1_packet_2_idx0:NR1_packet_2_idx1));
                                            NR1_packet_1 = nanmean(rate_long(NR1_packet_1_idx0:NR1_packet_1_idx1));
                                        else
                                            NR1_packet_1 = NaN;
                                            NR1_packet_2 = NaN;                                            
                                        end
                                        % REM2 (labelled as NR2 - legacy)
                                        if NR2_dur >= REM_dur_thresh
                                            NR2_packet_1_idx0 = ceil( (NR2_t0 + min_end_excl_REM - mean_long_t0) / G_bin);
                                            NR2_packet_1_idx1 = floor( (NR2_t0 + min_end_excl_REM + min_fr_len_REM - mean_long_t0) / G_bin);
                                            NR2_packet_2_idx1 = ceil( (NR2_t1 - min_end_excl_REM - mean_long_t0) / G_bin);
                                            NR2_packet_2_idx0 = floor( (NR2_t1 - min_end_excl_REM - min_fr_len_REM - mean_long_t0) / G_bin);
                                            NR2_packet_1 = nanmean(rate_long(NR2_packet_1_idx0:NR2_packet_1_idx1));
                                            NR2_packet_2 = nanmean(rate_long(NR2_packet_2_idx0:NR2_packet_2_idx1));
                                        else
                                            NR2_packet_1 = NaN;
                                            NR2_packet_2 = NaN;
                                        end
                                        
                                        % NR (labelled as REM  -legacy)
                                        
                                        REM_packet_1_idx0 = ceil( (REM_t0 + min_end_excl - mean_long_t0) / G_bin);
                                        REM_packet_1_idx1 = floor( (REM_t0 + min_end_excl + min_fr_len - mean_long_t0) / G_bin);
                                        REM_packet_2_idx0 = floor( (REM_t1 - min_end_excl - min_fr_len - mean_long_t0) / G_bin);
                                        REM_packet_2_idx1 = ceil( (REM_t1 - min_end_excl - mean_long_t0) / G_bin);
                                        if REM_packet_1_idx1 > REM_idx1 || REM_packet_2_idx0 < REM_idx0
                                            REM_packet_1 = NaN;
                                            REM_packet_2 = NaN;
                                        else
                                            REM_packet_1 = nanmean(rate_long(REM_packet_1_idx0:REM_packet_1_idx1));
                                            REM_packet_2 = nanmean(rate_long(REM_packet_2_idx0:REM_packet_2_idx1)); 
                                        end
                                        
                                    end
                                    
                                    
                                    if do_z
                                        NR1_zFR = (NR1_FR - frmean_long) / frstd_long;
                                        NR2_zFR = (NR2_FR - frmean_long) / frstd_long;
                                        
                                        if do_thirds
                                            NR1_T1_zFR = (NR1_T1_FR - frmean_long) / frstd_long;
                                            NR1_T2_zFR = (NR1_T2_FR - frmean_long) / frstd_long;
                                            NR1_T3_zFR = (NR1_T3_FR - frmean_long) / frstd_long;
                                            
                                            NR2_T1_zFR = (NR2_T1_FR - frmean_long) / frstd_long;
                                            NR2_T2_zFR = (NR2_T2_FR - frmean_long) / frstd_long;
                                            NR2_T3_zFR = (NR2_T3_FR - frmean_long) / frstd_long;
                                            
                                            REM_T1_zFR = (REM_T1_FR - frmean_long) / frstd_long;
                                            REM_T2_zFR = (REM_T2_FR - frmean_long) / frstd_long;
                                            REM_T3_zFR = (REM_T3_FR - frmean_long) / frstd_long;
                                        end
                                        
                                        if do_packets
                                            
                                            NR1_packet_zFR = (NR1_packet_FR - frmean_long) / frstd_long;
                                            NR2_packet_zFR = (NR2_packet_FR - frmean_long) / frstd_long;
                                            
                                            NR1_packet_1_zFR = (NR1_packet_1 - frmean_long) / frstd_long;
                                            NR1_packet_2_zFR = (NR1_packet_2 - frmean_long) / frstd_long;
                                            
                                            NR2_packet_1_zFR = (NR2_packet_1 - frmean_long) / frstd_long;
                                            NR2_packet_2_zFR = (NR2_packet_2 - frmean_long) / frstd_long;
                                            
                                            REM_packet_1_zFR = (REM_packet_1 - frmean_long) / frstd_long;
                                            REM_packet_2_zFR = (REM_packet_2 - frmean_long) / frstd_long;
                                            
                                        end
                                    else
                                        NR1_zFR = NR1_FR;
                                        NR2_zFR = NR2_FR;
                                        
                                        if do_thirds
                                            NR1_T1_zFR = NR1_T1_FR;
                                            NR1_T2_zFR = NR1_T2_FR;
                                            NR1_T3_zFR = NR1_T3_FR;
                                            
                                            NR2_T1_zFR = NR2_T1_FR;
                                            NR2_T2_zFR = NR2_T2_FR;
                                            NR2_T3_zFR = NR2_T3_FR;
                                            
                                            REM_T1_zFR = REM_T1_FR;
                                            REM_T2_zFR = REM_T2_FR;
                                            REM_T3_zFR = REM_T3_FR;
                                        end
                                        
                                        if do_packets
                                            NR1_packet_zFR = NR1_packet_FR;
                                            NR2_packet_zFR = NR2_packet_FR;
                                            
                                            NR1_packet_1_zFR = NR1_packet_1;
                                            NR1_packet_2_zFR = NR1_packet_2;
                                            
                                            NR2_packet_1_zFR = NR2_packet_1;
                                            NR2_packet_2_zFR = NR2_packet_2;
                                            
                                            REM_packet_1_zFR = REM_packet_1;
                                            REM_packet_2_zFR = REM_packet_2;
                                        end
                                    end
                                    
                                    
                                    delta_FR{anim_idx}(tt,cc) = NR2_zFR - NR1_zFR;
                                    allREM_dur{anim_idx}(tt,1) = REM_dur;
                                    allNR1_dur{anim_idx}(tt,1) = NR1_dur;
                                    allNR2_dur{anim_idx}(tt,1) = NR2_dur;
                                    
                                    if do_power
                                        allREM_theta{anim_idx}(tt,1) = REM_theta_power;
                                        allNR1_delta{anim_idx}(tt,1) = NR1_delta_power;
                                        allTD_ratio{anim_idx}(tt,1) = REM_td_ratio;
                                    end
                                    
                                    if do_thirds
                                        NR1_thirds_allFR{anim_idx}(tt,cc,1) = NR1_T1_zFR;
                                        NR1_thirds_allFR{anim_idx}(tt,cc,2) = NR1_T2_zFR;
                                        NR1_thirds_allFR{anim_idx}(tt,cc,3) = NR1_T3_zFR;
                                        
                                        NR2_thirds_allFR{anim_idx}(tt,cc,1) = NR2_T1_zFR;
                                        NR2_thirds_allFR{anim_idx}(tt,cc,2) = NR2_T2_zFR;
                                        NR2_thirds_allFR{anim_idx}(tt,cc,3) = NR2_T3_zFR;
                                        
                                        REM_thirds_allFR{anim_idx}(tt,cc,1) = REM_T1_zFR;
                                        REM_thirds_allFR{anim_idx}(tt,cc,2) = REM_T2_zFR;
                                        REM_thirds_allFR{anim_idx}(tt,cc,3) = REM_T3_zFR;

                                    end
                                    
                                    if do_spindles
                                        spindle_incidence{anim_idx}(tt,1) = spindle_inc;
                                    end
                                    
                                    if do_SWA
                                        SWA_amplitude{anim_idx}{tt,1} = swa_ampl;
                                    end
                                    
                                    if do_packets
                                        NR1_packets_allFR{anim_idx}(tt,cc) = NR1_packet_zFR;
                                        NR2_packets_allFR{anim_idx}(tt,cc) = NR2_packet_zFR;
                                        
                                        NR1_p1_allFR{anim_idx}(tt,cc) = NR1_packet_1_zFR;
                                        NR1_p2_allFR{anim_idx}(tt,cc) = NR1_packet_2_zFR;
                                        NR2_p1_allFR{anim_idx}(tt,cc) = NR2_packet_1_zFR;
                                        NR2_p2_allFR{anim_idx}(tt,cc) = NR2_packet_2_zFR;
                                        REM_p1_allFR{anim_idx}(tt,cc) = REM_packet_1_zFR;
                                        REM_p2_allFR{anim_idx}(tt,cc) = REM_packet_2_zFR;
                                    end
                                else
                                    badcount = badcount+1;
                                    %                             delta_FR{anim_idx}(tt,cc) = NaN;
                                end
                                
                            else
                                %                 delta_FR(tt,cc) = NaN;
                                %             fprintf('Bad on/off times.\n');
                            end
                        end
                    end
                else
                    %             fprintf('bad\n');
                end
            end
        end
    end
end

fprintf('Total bad count: %u.\n', badcount);

%keyboard;
skip_it = 1;
if ~skip_it
    %% process data
    do_avg_by_cell = 1;
    % delta_FR_triplets = nanmean(delta_FR,2);
    
    if do_avg_by_cell
        
        delta_FR_byanim = cellfun(@(x) nanmean(x,2), delta_FR, 'UniformOutput', false);
        delta_FR_nonan_byanim_idx = cellfun(@(x) ~isnan(x), delta_FR_byanim, 'UniformOutput', false);
        for uu = 1:size(delta_FR_byanim,2)
            delta_FR_nonan_byanim{uu} = delta_FR_byanim{uu}(delta_FR_nonan_byanim_idx{uu});
        end
        delta_FR_nonan = cat(1,delta_FR_nonan_byanim{:});
        
        %     REMdur_nonan_byanim_idx = cellfun(@(x) ~isnan(x), allREM_dur, 'UniformOutput', false);
        for uu = 1:size(delta_FR_byanim,2)
            REMdur_nonan_byanim{uu} = allREM_dur{uu}(delta_FR_nonan_byanim_idx{uu});
        end
        REMdur_nonan = cat(1,REMdur_nonan_byanim{:});
        
        for uu = 1:size(delta_FR_byanim,2)
            REMtheta_nonan_byanim{uu} = allREM_theta{uu}(delta_FR_nonan_byanim_idx{uu});
        end
        REMtheta_nonan = cat(1,REMtheta_nonan_byanim{:});
        
        
        for uu = 1:size(delta_FR_byanim,2)
            NR1delta_nonan_byanim{uu} = allNR1_delta{uu}(delta_FR_nonan_byanim_idx{uu});
        end
        NR1delta_nonan = cat(1,NR1delta_nonan_byanim{:});
        
        
        for uu = 1:size(delta_FR_byanim,2)
            TDratio_nonan_byanim{uu} = allTD_ratio{uu}(delta_FR_nonan_byanim_idx{uu});
        end
        TDratio_nonan = cat(1,TDratio_nonan_byanim{:});
        
        if do_spindles
            for uu = 1:size(delta_FR_byanim,2)
                spindleinc_nonan_byanim{uu} = spindle_incidence{uu}(delta_FR_nonan_byanim_idx{uu});
            end
            spindleinc_nonan = cat(1,spindleinc_nonan_byanim{:});
        end
        
        if do_SWA
            for uu = 1:size(delta_FR_byanim,2)
                if iscell(SWA_amplitude{uu})
                    nancheck = all(cell2mat(cellfun(@(x) isnan(mean(x)),SWA_amplitude{uu},'uniformoutput',0)));
                elseif isnumeric(SWA_amplitude{uu})
                    nancheck = isnan(SWA_amplitude{uu});
                end
                if isempty(SWA_amplitude{uu}) || nancheck
                    swa_mean{uu} = [];
                else
                    swa_mean{uu} = cellfun(@nanmean, SWA_amplitude{uu});
                end
                swa_nonan_byanim{uu} = swa_mean{uu}(delta_FR_nonan_byanim_idx{uu});
                swa_nonan_byanim_norm{uu} = 100 * swa_nonan_byanim{uu} ./ mean_SWA{uu};
            end
            SWA_nonan = cat(1,swa_nonan_byanim_norm{:});
        end
        
    else
        
        delta_FR_byanim = cellfun(@(x) nanmean(x,2), delta_FR, 'UniformOutput', false);
        
    end
    
    [rho,p_rho] = corr(REMdur_nonan, delta_FR_nonan);
    asterisks = get_asterisks_from_pval(p_rho);
    lincoeff = polyfit(REMdur_nonan, delta_FR_nonan, 1);
    
    %% plotting
    markersize = 7;
    mylims = 0;
    
    
    setFigureDefaults;
    figure();
    if dep
        scatter(REMdur_nonan,delta_FR_nonan,markersize^2,'k','filled');
    else
        scatter(REMdur_nonan,delta_FR_nonan,markersize^2,'k');
    end
    hold on;
    if mylims
        xlims = [0 600];
        ylims = [-2 2];
    else
        xlims = get(gca,'xlim');
        ylims = get(gca,'ylim');
    end
    X = xlims(1):10:xlims(end);
    Y = lincoeff(2) + lincoeff(1).*X;
    plot(X,Y,'--k','linewidth',1.5);
    text(xlims(2)*.8,ylims(2)*.9,['\it r\rm' sprintf(' = %.3f%s',rho,asterisks)],'fontsize',16);
    text(xlims(2)*.8,ylims(2)*.80,['\it p\rm' sprintf(' = %.3f',p_rho)],'fontsize',16);
    set(gca,'xlim',xlims,'ylim',ylims);
    ylabel('Firing rate (z)','fontsize',21);
    xlabel('REM duration (s)','fontsize',20);
    
    %% Theta power scatter
    markersize = 7;
    mylims = 0;
    
    [rho_th,p_rho_th] = corr(REMtheta_nonan, delta_FR_nonan);
    asterisks_th = get_asterisks_from_pval(p_rho_th);
    lincoeff_th = polyfit(REMtheta_nonan, delta_FR_nonan, 1);
    
    setFigureDefaults;
    figure();
    if dep
        scatter(REMtheta_nonan,delta_FR_nonan,markersize^2,'k','filled');
    else
        scatter(REMtheta_nonan,delta_FR_nonan,markersize^2,'k');
    end
    hold on;
    if mylims
        xlims = [0 600];
        ylims = [-2 2];
    else
        xlims = get(gca,'xlim');
        ylims = get(gca,'ylim');
    end
    X = linspace(xlims(1),xlims(end),50);
    Y = lincoeff_th(2) + lincoeff_th(1).*X;
    plot(X,Y,'--k','linewidth',1.5);
    text(xlims(2)*.8,ylims(2)*.9,['\it r\rm' sprintf(' = %.3f%s',rho_th,asterisks_th)],'fontsize',16);
    text(xlims(2)*.8,ylims(2)*.80,['\it p\rm' sprintf(' = %.3f',p_rho_th)],'fontsize',16);
    set(gca,'xlim',xlims,'ylim',ylims);
    ylabel('Firing rate (z)','fontsize',21);
    xlabel('Theta power (z)','fontsize',20);
    
    
    %% Delta power scatter
    markersize = 7;
    mylims = 0;
    
    [rho_d,p_rho_d] = corr(NR1delta_nonan, delta_FR_nonan);
    asterisks_d = get_asterisks_from_pval(p_rho_d);
    lincoeff_d = polyfit(NR1delta_nonan, delta_FR_nonan, 1);
    
    setFigureDefaults;
    figure();
    if dep
        scatter(NR1delta_nonan,delta_FR_nonan,markersize^2,'k','filled');
    else
        scatter(NR1delta_nonan,delta_FR_nonan,markersize^2,'k');
    end
    hold on;
    if mylims
        xlims = [0 600];
        ylims = [-2 2];
    else
        xlims = get(gca,'xlim');
        ylims = get(gca,'ylim');
    end
    X = linspace(xlims(1),xlims(end),50);
    Y = lincoeff_d(2) + lincoeff_d(1).*X;
    plot(X,Y,'--k','linewidth',1.5);
    text(xlims(2)*.8,ylims(2)*.9,['\it r\rm' sprintf(' = %.3f%s',rho_d,asterisks_d)],'fontsize',16);
    text(xlims(2)*.8,ylims(2)*.80,['\it p\rm' sprintf(' = %.3f',p_rho_d)],'fontsize',16);
    set(gca,'xlim',xlims,'ylim',ylims);
    ylabel('Firing rate (z)','fontsize',21);
    xlabel('Delta power (z)','fontsize',20);
    
    %% Theta/delta ratio scatter
    markersize = 7;
    mylims = 0;
    
    [rho_td,p_rho_td] = corr(TDratio_nonan, delta_FR_nonan);
    asterisks_td = get_asterisks_from_pval(p_rho_td);
    lincoeff_td = polyfit(TDratio_nonan, delta_FR_nonan, 1);
    
    setFigureDefaults;
    figure();
    if dep
        scatter(TDratio_nonan,delta_FR_nonan,markersize^2,'k','filled');
    else
        scatter(TDratio_nonan,delta_FR_nonan,markersize^2,'k');
    end
    hold on;
    if mylims
        xlims = [0 600];
        ylims = [-2 2];
    else
        xlims = get(gca,'xlim');
        ylims = get(gca,'ylim');
    end
    X = linspace(xlims(1),xlims(end),50);
    Y = lincoeff_td(2) + lincoeff_td(1).*X;
    plot(X,Y,'--k','linewidth',1.5);
    text(xlims(2)*.8,ylims(2)*.9,['\it r\rm' sprintf(' = %.3f%s',rho_td,asterisks_td)],'fontsize',16);
    text(xlims(2)*.8,ylims(2)*.80,['\it p\rm' sprintf(' = %.3f',p_rho_td)],'fontsize',16);
    set(gca,'xlim',xlims,'ylim',ylims);
    ylabel('Firing rate (z)','fontsize',21);
    xlabel('Theta/delta ratio','fontsize',20);
    
    %% Spindle incidence scatter
    if do_spindles
        markersize = 7;
        mylims = 0;
        
        [rho_sp,p_rho_sp] = corr(spindleinc_nonan, delta_FR_nonan);
        asterisks_sp = get_asterisks_from_pval(p_rho_sp);
        lincoeff_sp = polyfit(spindleinc_nonan, delta_FR_nonan, 1);
        
        setFigureDefaults;
        figure();
        if dep
            scatter(spindleinc_nonan,delta_FR_nonan,markersize^2,'k','filled');
        else
            scatter(spindleinc_nonan,delta_FR_nonan,markersize^2,'k');
        end
        hold on;
        if mylims
            xlims = [0 600];
            ylims = [-2 2];
        else
            xlims = get(gca,'xlim');
            ylims = get(gca,'ylim');
        end
        X = linspace(xlims(1),xlims(end),50);
        Y = lincoeff_sp(2) + lincoeff_sp(1).*X;
        plot(X,Y,'--k','linewidth',1.5);
        text(xlims(2)*.8,ylims(2)*.9,['\it r\rm' sprintf(' = %.3f%s',rho_sp,asterisks_sp)],'fontsize',16);
        text(xlims(2)*.8,ylims(2)*.80,['\it p\rm' sprintf(' = %.3f',p_rho_sp)],'fontsize',16);
        set(gca,'xlim',xlims,'ylim',ylims);
        ylabel('Firing rate (z)','fontsize',21);
        xlabel('Spindle incidence (1/s)','fontsize',20);
    end
    
    %% SWA scatter
    if do_SWA
        markersize = 7;
        mylims = 0;
        
        [rho_swa,p_rho_swa] = corr(SWA_nonan, delta_FR_nonan);
        asterisks_swa = get_asterisks_from_pval(p_rho_swa);
        lincoeff_swa = polyfit(SWA_nonan, delta_FR_nonan, 1);
        
        setFigureDefaults;
        figure();
        if dep
            scatter(SWA_nonan,delta_FR_nonan,markersize^2,'k','filled');
        else
            scatter(SWA_nonan,delta_FR_nonan,markersize^2,'k');
        end
        hold on;
        if mylims
            xlims = [0 600];
            ylims = [-2 2];
        else
            xlims = get(gca,'xlim');
            ylims = get(gca,'ylim');
        end
        X = linspace(xlims(1),xlims(end),50);
        Y = lincoeff_swa(2) + lincoeff_swa(1).*X;
        plot(X,Y,'--k','linewidth',1.5);
        text(xlims(2)*.8,ylims(2)*.9,['\it r\rm' sprintf(' = %.3f%s',rho_swa,asterisks_swa)],'fontsize',16);
        text(xlims(2)*.8,ylims(2)*.80,['\it p\rm' sprintf(' = %.3f',p_rho_swa)],'fontsize',16);
        set(gca,'xlim',xlims,'ylim',ylims);
        ylabel('Firing rate (z)','fontsize',21);
        xlabel('SWA amplitude (% of 36-hour mean)','fontsize',20);
    end
end
%% save data

%keyboard;
if save_the_data
    fprintf('Saving data...\n');
    tic
    
    if ispc
        save_dir = 'Z:\ATP_MAIN\CODE\lfp_analysis_beta\Analysis_Data\RNR';
    elseif ismac
        save_dir = '/Volumes/turrigiano-lab/ATP_MAIN/CODE/lfp_analysis_beta/Analysis_Data/RNR';
    end
    dep_str = {'CTRL','DEP'};
    save_file = ['triplet_data_rnr_v' num2str(ID_NUM) '_'  dep_str{dep+1} '.mat'];
    
    triplet_data.delta_FR      = delta_FR;
    triplet_data.allREM_dur    = allREM_dur;
    triplet_data.allNR1_dur    = allNR1_dur;
    triplet_data.allNR2_dur    = allNR2_dur;
    
    if do_power
        triplet_data.allREM_theta  = allREM_theta;
        triplet_data.allNR1_delta  = allNR1_delta;
        triplet_data.allTD_ratio   = allTD_ratio;
    end
    if do_spindles
        triplet_data.spindle_incidence   = spindle_incidence;
    end
    if do_SWA
        triplet_data.SWA_amplitude = SWA_amplitude;
        triplet_data.mean_SWA_ampl = mean_SWA;
    end
    if do_thirds
        triplet_data.NR1_thirds_allFR = NR1_thirds_allFR;
        triplet_data.NR2_thirds_allFR = NR2_thirds_allFR;
        triplet_data.REM_thirds_allFR = REM_thirds_allFR;

    end
    if do_packets
        triplet_data.NR1_packets_allFR = NR1_packets_allFR;
        triplet_data.NR2_packets_allFR = NR2_packets_allFR;
        
        triplet_data.NR1_p1_allFR = NR1_p1_allFR;
        triplet_data.NR1_p2_allFR = NR1_p2_allFR;
        triplet_data.NR2_p1_allFR = NR2_p1_allFR;
        triplet_data.NR2_p2_allFR = NR2_p2_allFR;
        triplet_data.REM_p1_allFR = REM_p1_allFR;
        triplet_data.REM_p2_allFR = REM_p2_allFR;

    end
    
    triplet_data.params.homeo_start = homeo_period_start;
    triplet_data.params.homeo_end = homeo_period_end;
    triplet_data.params.t_start = t_start;
    triplet_data.params.t_stop = t_stop;
    triplet_data.params.dep = dep;
    triplet_data.params.qthresh = qthresh;
    triplet_data.params.perc_thresh = perc_thresh;
    triplet_data.params.norm = norm;
    triplet_data.params.bl_on = bl_on;
    triplet_data.params.dataset = dataset;
    triplet_data.params.negpos = negpos;
    triplet_data.params.tailslope = tailslope;
    triplet_data.params.G_bin = G_bin;
    triplet_data.params.use_periods = use_periods;
    triplet_data.params.mean_t = mean_t;
    triplet_data.params.duration_threshold = duration_threshold;
    triplet_data.params.anim_list = anim_list;
    triplet_data.params.rec_periods = rec_periods;
    
    save([save_dir filesep save_file],'triplet_data','-v7.3');
    toc
end







