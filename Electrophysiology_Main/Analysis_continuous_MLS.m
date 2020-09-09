%% Analysis_Continuous_MLS
% 
% Alejandro Torrado Pacheco - 2017
% 
% This is the main analysis script for electrophysiology data collected
% during eye re-opening experiments. This script can be used with all the
% different datasets (regular ER, CPP injections, Sleep Deprivation, etc.).
% 
% The script finds continuous RSUs (parameters to set threshold for this
% are set by user) and generates a plot of the average (normalized or not)
% FR for the RSU population in that dataset. It also outputs plots of the
% mean FR by animal and the mean FR by cell (organized by animal).
%
% A variant of this code can be used to generate intermediary data
% structures that can then be used to generate the figures in the paper
% (Torrado Pacheco et al., 2019; doi: doi.org/10.1101/827832 - bioRxiv).

% load CONTCELL variable manually before you start!

% clear workspace
clearvars -except CONTCELL*

%% SETUP
% select dataset to analyze. Options:
% - 'ER' : continuous eye re-opening (BL1-ER4)
% - 'recov' :  last 5 days of ER (MD4-ER4)
% - 'CPPrecov' : recov with CPP injection on ER1
% - 'CPP2recov' : recov with CPP injections on ER and ER3
% - 'SD_NEW' : sleep deprivation
dataset = 'CPP2recov';
eval(['CELL = CONTCELL_' dataset '.MASTER;']);

% define colors for plotting
cred = [0.85 0.33 0.01];
cblu = [0 0.45 0.74];
cyel = [0.93 0.69 0.13];
cpur = [0.49 0.18 0.55];
cgre = [0.47 0.67 0.19];
ccya = [0.30 0.75 0.93];
cmar = [0.64 0.08 0.18];
cblk = [0 0 0];
c_ctrl = [0 0 0];
c_dep = [112,166,217]./255;
colorpool = {cblu,cred,cyel,cpur,cgre,ccya,cmar,cblk,c_ctrl};

% parameters for analysis
% quality threshold (2 for single units)
dep = 1; % deprived (1) or control (0) hemisphere
dep_status = {'CONTROL','DEPRIVED'};

% quality threshold: should be 2
qthresh = 2;
doRSU = 1; % leave this as 1

% exclude cells based on permutation analysis of continuity (1=yes, 0=no)
permut_drop = 0;

% Analysis & plotting parameters
plotcellsep = 0; % plot RSU vs FS cell separation
normalize = 1; % normalize FRs
no_big_vals = 1; % exclude large values 
big_val_thresh = 12; % threshold for large values (e.g. 12-fold change in FR)
no_small_vals = 1; % exclude very small values
small_val_thresh = 0.01; % small value threshold (100-fold decrease)
val_thresh_perc = 0.01; % percent of value under threshold that will trigger exclusion
percThresh = 0.80; % continuity threshold - % of rec time that cells must be online to be analyzed
fprintf('Continuity threshold is: %u%%.\n\n',percThresh*100);
max_bl_change = 2; % maximum fold change during baseline for baseline to be considered stable
min_bl_change = 1/max_bl_change; % maximum downward change, same thing as above
require_BLON = 1; % require that neurons be online during baseline
blon_hour = 48; % baseline start; this is ignored for RECOV data (where blon_hour is hard-coded to 156)
if normalize && ~require_BLON
    warning('WARNING: asked to normalize but not requiring cells to be ON at baseline may result in errors.');
end


% for plotting main fig
ymax_main = 4;

% create arrays for indexing of cells
quality = [CELL.quality];
deprived = [CELL.deprived];
neg_pos_time = [CELL.neg_pos_time];
tail_slope = [CELL.tailSlope];
animalname = {CELL.animal};

% animals to examine - this depends on chosen dataset
% this section can be used to exclude animals from analysis
switch dataset
    case 'ER'
        if dep
            which_anim = {'KH67','AT14','AT16',...
                'KH73','AT29','AT22','AT25','AT29','KH75'};
        else
            which_anim = {'KH67','AT14','AT22','AT16','AT25',...
                'AT27','AT29','KH72','KH73'};
        end
    case 'recov'
        sep_dep = 1;
        if sep_dep
            if dep
                which_anim = {'KH67','AT12','AT14','AT16','AT22','AT25','AT27','AT29'};
            else
                which_anim = {'KH67','AT12','AT14','AT16','AT22','AT25','AT27','AT29','KH72','KH73','KH75'};
            end
        else
            which_anim = {'KH67','AT12','AT14','AT22','AT24','AT25','AT27','AT29','KH73','KH75'};
        end
    case 'CPPrecov'
        which_anim = {'AT19','AT20','AT23','AT32'};
    case 'CPP2recov'
        which_anim = {'AT50','AT51','AT52','AT53','AT54','AT56','AT57','AT58','AT59','AT60','AT61'};
    case 'SD_NEW'
        which_anim = {'AT63','AT64','AT65','AT66','AT67','AT68','AT69'};
end

%% Filtering cells by continuity and baseline presence
% Loop through neurons
for dd = 1:length(CELL)
    % calculate total time each cell is 'online'
    totT = [];
    for uu = 1:length(CELL(dd).onTime)
        totT(uu) = CELL(dd).offTime(uu) - CELL(dd).onTime(uu);
    end
    
    % depending on dataset, set total recording time
    switch dataset
        % 10.5 days for full continuous recordings
        case 'ER'
            if regexp(CELL(dd).animal,regexptranslate('wildcard','KH*'))
                first_number = str2double(regexp(CELL(dd).animal,'\d','match','once'));
                if first_number > 5
                    daymax = 10.5;
                else
                    daymax = 8;
                end
            elseif regexp(CELL(dd).animal,regexptranslate('wildcard','AT*'))
                daymax = 10.5;
            elseif strcmp(CELL(dd).animal,'AT09')
                daymax = 10.5;
            end
            daysRecorded = daymax - CELL(dd).dayStart;
            
        % 5 days for 'recov' data
        case {'recov','CPPrecov','CPP2recov','SD_NEW'}
            % fixed number of recorded days, from MD4 to ER4 (included)
            daysRecorded = 5;
    end
    
    % if neuron is online for long enough, set percentOn flag to 1
    if sum(totT) < (percThresh * daysRecorded*24*3600)
        percentOn(dd)        = 0;
    else
        percentOn(dd) = 1;
    end
    
    if require_BLON
        % cell must come online before start of "baseline"
        % "Baseline" is 36-60 hours for animals starting before MD1
        % "Baseline" is MD4 for AT12 (starting at MD3)
        % "Baseline" is MD4 night for all recov data
        switch dataset
            case {'ER'}
                timeon = CELL(dd).onTime(1) + CELL(dd).dayStart*24*3600;
                if timeon <= blon_hour*3600
                    onatBase(dd) = 1;
                elseif CELL(dd).onTime(1) <= 6.0*24*3600 && CELL(dd).dayStart > 2
                    onatBase(dd) = 1;
                else
                    onatBase(dd) = 0;
                end
                
            case {'recov','CPPrecov'}
                timeon = CELL(dd).onTime(1);
                if timeon <= 6.5*24*3600
                    onatBase(dd) = 1;
                else
                    onatBase(dd) = 0;
                end
            case{'SD_NEW','CPP2recov'}
                timeon = CELL(dd).onTime(1) + CELL(dd).dayStart*24*3600;
                if timeon <= 6.5*24*3600
                    onatBase(dd) = 1;
                else
                    onatBase(dd) = 0;
                end                
        end
    else
        onatBase(dd) = 1;
    end
    
end

%% INDEXING OF CELLS

% find cells corresponding to animals identified above
thisanim_idx = [];
tmp_anim = [];
for xx=1:length(which_anim)
    thisone = which_anim{xx};
    tmp_anim = [tmp_anim find(strcmp(thisone,animalname)==1)];
end
thisanim_idx = tmp_anim;
clear tmp_anim;

% initialize variables - these are all indexing flags (1 if cell is to be
% analyzed, 0 otherwise)
depidx = []; qualidx = []; percidx = []; idx0 = []; idx1 = [];
% deprived, quality, percentageOn, baselineOn - set the flags
depidx  = find(deprived == dep);
qualidx = find(quality <= qthresh);
percidx = find(percentOn == 1);
blidx = find(onatBase == 1);
% intersect the flags to find cells that satisfy all of them
idx00 = intersect(qualidx,depidx);
idx01 = intersect(blidx,percidx);
idx0 = intersect(idx00,idx01);
% final result - cross with animals to get final set of cells to analyze
idx1    = intersect(idx0,thisanim_idx);


% Cell classification into RSU/pFS cells
% Most recent negposthresh, tailslope and halfwidth should all be
% independent of interpolation factor
idx2 = []; idx3 = []; idx4 = []; idx5 = []; ns = [];

% go through cells
for cellno = idx1
    
    % set thresholds for RSU vs FS discrimination
    % these can be set separately for RSU vs FS, if desired
    negpos_thresh = 0.39;
    fs_negpos_thresh = 0.39;
    tailSlope_thresh = 0.005;
    fs_tailSlope_thresh = 0.005;
    
    % sort based on negpostime
    if CELL(cellno).neg_pos_time >= negpos_thresh
        idx2(end+1) = cellno;
    elseif CELL(cellno).neg_pos_time < fs_negpos_thresh
        idx3(end+1) = cellno;
    end
    
    % sort based on tailslope
    if CELL(cellno).tailSlope >= tailSlope_thresh
        idx4(end+1) = cellno;
    elseif CELL(cellno).tailSlope < fs_tailSlope_thresh
        idx5(end+1) = cellno;
    end
end

%% PLOT CELL SEPARATION
cellsep_count = 0;
if plotcellsep
    
    % get the negpostime, tailslope, and average peak-scaled WF
    for cx = idx1
        cellsep_count = cellsep_count + 1;
        cellsep_wf{cellsep_count} = CELL(cx).scaledWF;
        cellsep_nt(cellsep_count) = neg_pos_time(cx);
        cellsep_ts(cellsep_count) = tail_slope(cx);
        if neg_pos_time(cx) >= negpos_thresh && tail_slope(cx) >= tailSlope_thresh
            cellsep_color{cellsep_count} = cblu;
        elseif neg_pos_time(cx) < fs_negpos_thresh && tail_slope(cx) < fs_tailSlope_thresh
            cellsep_color{cellsep_count} = cyel;
        else
            cellsep_color{cellsep_count} = [.6 .6 .6];
        end
    end
    
    % make the figure
    cellsepfig = figure(222);
    set(cellsepfig,'visible','off','numbertitle','off');
    % subplot 1: plot the waveforms of each cell, color depends on RSU vs
    % FS classification
    % subplot 2: plot the negpostime vs tailslope for each WF, color by
    % cell type as above
    for cc = 1:cellsep_count
        s1 = subplot(1,2,1); hold on;
        plot(s1,cellsep_wf{cc},...
            'color',cellsep_color{cc},'linewidth',2);
        s2 = subplot(1,2,2); hold on;
        plot(s2,cellsep_nt(cc),cellsep_ts(cc),'o','color',cellsep_color{cc},...
            'MarkerFaceColor',cellsep_color{cc},'MarkerSize',10);
    end
    % threshold lines
    line([fs_negpos_thresh fs_negpos_thresh],[-0.1 fs_tailSlope_thresh],...
        'color','k','linestyle','--')
    line([0 fs_negpos_thresh],[fs_tailSlope_thresh fs_tailSlope_thresh],...
        'color','k','linestyle','--')
    line([negpos_thresh negpos_thresh],[tailSlope_thresh 0.1],...
        'color','k','linestyle','--')
    line([negpos_thresh 1],[tailSlope_thresh tailSlope_thresh],...
        'color','k','linestyle','--')
    % format axes
    set(s2,'ylim',[-0.1 0.1],'xlim',[0.2 0.60]);
    set(cellsepfig,'units','normalized','position',[0.1 0.1 0.7 0.7],'visible','on');
    set(s1,'fontsize',18,'linewidth',2);
    set(s2,'fontsize',18,'linewidth',2);
    s1.XAxis.Color = 'k';
    s2.XAxis.Color = 'k';
end


%% SETUP FOR FR ANALYSIS

% get the RSUs
RSU_idx = [];
RSU_idx = intersect(intersect(idx1,idx2),idx4);

% get the FS cells
pFS_idx = [];
pFS_idx = intersect(intersect(idx1,idx3),idx5);

% print info
fprintf('Found %u continuous RSUs that are on at BL and persist for %u%% of time.\n\n',...
    size(RSU_idx,2),percThresh*100);

fprintf('Found %u continuous pFS cells that are on at BL and persist for %u%% of time.\n\n',...
    size(pFS_idx,2),percThresh*100);

%% from permutation test analysis
% not this only works for recov data
if permut_drop && strcmp(dataset,'recov')
    if dep
        badz = [322, 384, 594, 605, 610];
    else
        badz = [157, 238, 454, 461, 485, 507];
    end
    % remove these neurons
    RSU_idx = setdiff(RSU_idx,badz);
end

% make separate cell variables for RSU & FS
% ---- NOTE: I did not analyze the FS neuron's FR activity based on MD/ER,
%            because I did not record enough of them. This code used to be
%            setup to analyze FS cells too, but I have removed that
%            section. Thus, the CELL_pFS variable is unused.
CELL_RSU = CELL(RSU_idx);
CELL_pFS = CELL(pFS_idx);

% clear variables
clear t_start t_stop max_h bin_sz bin_edges

% estimate FR using gaussian kernel. Create kernel here
G_sigma = 300; % s.dev. of Gaussian kernel (seconds) **** standard is 300! ****
G_bin = max(1,G_sigma/5); % bin size with which to estimate Gaussian kernel (seconds)
G_xvals = -5*G_sigma : G_bin : G_sigma*5; % values at which to estimate kernel
G_kernel = ( 1/(sqrt(2*pi)*G_sigma) ) *exp(-(G_xvals.^2)/(2*G_sigma^2)); % kernel creation

% total time over which to calculate FR
t_start = 0*24*3600;
t_stop  = 12*24*3600;

% normalization based on dataset
switch dataset
    case {'recov','CPPrecov','CPP2recov','SD_NEW'}
        % For normalization, get bin numbers corresponding to MD4 (night)
        % note that "recov" recordings are clustered to start on MD4 day
        bl_start = 6.5*24*3600/G_bin;
        bl_end   = 7*24*3600/G_bin;
    case {'ER'}
        % For normalization, get bin numbers corresponding to
        % BL2night/BL3day
        bl_start = 1.5*24*3600/G_bin;
        bl_end   = 2.5*24*3600/G_bin;
        
end
max_h = t_stop / 3600;


% get rid of FR artifact due to unplugging for eye re-opening procedure
switch dataset
    case {'ER','recov'}
        ER_start = (7*24*3600 + 2.0*3600)/G_bin; % 2.0 hours after start of ER1 (9:30 am)
        ER_end   = (7*24*3600 + 4.5*3600)/G_bin; % 4.5 hours after start of ER1 (12:00 pm)
    case 'CPP'
        ER_start = (5*24*3600 + 2.5*3600)/G_bin; % 2.5 hours after start of ER1 (10:00 am)
        ER_end   = (5*24*3600 + 5.0*3600)/G_bin; % 5.0 hours after start of ER1 (12:30 pm)
    case 'SD_NEW'
        ER_start = (7*24*3600 + 1.1*3600)/G_bin; % 1.1 hours after start of ER1 (~8:36 am)
        ER_end   = (7*24*3600 + 2.2*3600)/G_bin; % 2.2 hours after start of ER1 (~9:42 am)
    case {'CPPrecov','CPP2recov'}
        ER_start = (7*24*3600 + 2.0*3600)/G_bin; % 2.0 hours after start of ER1 (9:30 am)
        ER_end   = (7*24*3600 + 4.5*3600)/G_bin; % 4.5 hours after start of ER1 (12:00 pm)
end

%% Main loop for analysis of RSUs
% check for flag
if doRSU
    fprintf('\n\n      *** RSU cells. ***   \n\n');
    rsu_count = 0;
    
    % loop through cells
    for ee = 1:length(CELL_RSU)
        
        % some datasets have spiketimes aligned to their specific
        % experiment start, others already aligned to hypothetical BL1.
        % This ensures consistency: align everything to BL1.
        switch dataset
            case {'ER','CPP2recov','SD_NEW'}
                spikes  = CELL_RSU(ee).time + CELL_RSU(ee).dayStart*24*3600;
            case {'recov','CPPrecov'}
                spikes  = CELL_RSU(ee).time;
        end

        % Get spike raster in bins (G_bin sec size. This will be used to
        % estimate FR by convolving with Gaussian kernel estimated with
        % same bin size
        spike_raster = histc(spikes,t_start:G_bin:t_stop);
        % Use convolution to get mean FR estimate
        temp = conv(spike_raster,G_kernel,'same'); % make sure to have option 'same' here! 
        
        if nanmean(temp) < 50 && nanmean(temp) >= 0.001 % if rate is within "acceptable" range
            % increase counter
            
            % assign to new variable
            rate_temp = temp;
            
            % collect on/off times for cell
            ontimes = CELL_RSU(ee).onTime;
            offtimes = CELL_RSU(ee).offTime;
            timelist(:,1) = ontimes;
            timelist(:,2) = offtimes;
            
            % Process on/off times - outputs are FR arrays with 'offline'
            % values set to NaN
            try
                rate = processOnOffTimes_ATP(rate_temp,timelist,G_bin,CELL_RSU(ee).dayStart,dataset);
            catch
                keyboard;
            end
            
            % normalization
            switch dataset
                case 'ER'
                    if CELL_RSU(ee).dayStart >= 2
                        % for AT12, normalize to MD4
                        bl_start_MD = 6.5*24*3600/G_bin; % start of MD4 night
                        bl_end_MD   = 7*24*3600/G_bin; % end of MD4 night
                        bl_rate = nanmean(rate(bl_start_MD:bl_end_MD));
                        if normalize
                            rate = rate./bl_rate;
                        end
                    else
                        % otherwise normalize to pre-MD baseline
                        bl_rate = nanmean(rate(bl_start:bl_end));
                        if normalize
                            if require_BLON
                                rate = rate./bl_rate;
                            else
                                ratefrac = 0.10;
                                ratestart = floor(ratefrac * size(rate,1));
                                rateend = ceil((1-ratefrac) * size(rate,1));
                                rate = rate./nanmean(rate(ratestart:rateend));
                            end
                        end
                    end
                % for recov data, normalize to pre-ER baseline
                case {'recov','CPPrecov','SD_recov','CPP2recov','BLINJ','SD_NEW'}
                    bl_rate = nanmean(rate(bl_start:bl_end));
                    if normalize
                        rate = rate./bl_rate;
                    end
                    
            end
            
            
            if no_big_vals
                % get rid of crazy large values
                if normalize
                    % if more than 1% of values are very large, exclude
                    if sum(rate > big_val_thresh) / size(rate,1) >= val_thresh_perc
                        disp('toobig');
                        rate(:) = NaN;
                    end
                else
                    % same thing but when FRs are not normalized
                    if sum((rate./bl_rate) > big_val_thresh) / size(rate,1) >= val_thresh_perc 
                        disp('toobig');
                        rate(:) = NaN;
                    end
                end
            end
            
            if no_small_vals
                % get rid of crazy small values - similar to above
                if normalize
                    if sum(rate<small_val_thresh) / size(rate,1) >= val_thresh_perc
                        disp('teeny');
                        rate(:) = NaN;
                    end
                else
                    if sum((rate./bl_rate) < small_val_thresh) / size(rate,1)  >= val_thresh_perc
                        disp('teeny');
                        rate(:) = NaN;
                    end
                end
            end
            
            % as long as rate is not all NaNs, save it in big variable
            if ~all(isnan(rate))
                rsu_count = rsu_count + 1;
                % assign to main variable
                if exist('FRbycell_RSU','var') && size(rate,1) > size(FRbycell_RSU,2)
                    rate(size(FRbycell_RSU,2)+1:end) = [];
                end
                FRbycell_RSU(rsu_count,:) = rate;
                anims_bycell{rsu_count} = CELL_RSU(ee).animal;
                good_idxs(rsu_count) = RSU_idx(ee);
            end
            
            clear timelist
        end
        
        
    end
    
    
    %% RSU compile data
    
    % find mean and SEM
    meanFR_RSU = nanmean(FRbycell_RSU);
    semFR_RSU = nanstd(FRbycell_RSU,[],1)./sqrt(new_rsu_count-1);
    
    % get rid of ER artifact
    meanFR_RSU(ER_start:ER_end) = NaN;
    semFR_RSU(ER_start:ER_end) = NaN;

    xsize = size(meanFR_RSU,2);
     
end

%% Plotting

% flag to toggle plotting on/off
plot_main_rsu = 1;
if plot_main_rsu && doRSU
    % make figure
    main_rsu_fig = figure();
    set(main_rsu_fig,'color','white');
    hold on;
    % set x-axis limits
    switch dataset
        case {'ER'}
            xstart = 0;
            xend = 264;
        case {'recov','CPPrecov','CPP2recov','SD_NEW'}
            xstart = 156;
            xend = 264;
    end
    
    % make light/dark edges
    edgesD = (12*3600)/G_bin : (24*3600)/G_bin : xsize;
    for qq = 1:size(edgesD,2)
        q = rectangle('Position',[edgesD(qq), 0, 12*3600/G_bin, 15] );
        set(q,'facecolor',[0.9 0.9 0.9],'linestyle','none');
    end
    
    % plot shaded error bar
    switch dataset
        case {'ER','recov'}
            thiscolors = {c_ctrl,cblu};
        case {'CPPrecov'}
            thiscolors = {c_ctrl,cmar};
        case {'SD_NEW'}
            thiscolors = {c_ctrl,cpur};
        case {'CPP2recov'}
            thiscolors = {c_ctrl,cgre};
    end
    shadedErrorBar(1:xsize,meanFR_RSU,semFR_RSU,{'color',thiscolors{dep+1},'linewidth',2});
    
    % title and axes format
    title(sprintf('%s HEMISPHERE\nAverage normalized RSU firing rate (n=%u)\nContinuity threshold: %u%%',...
        dep_status{dep+1},new_rsu_count,percThresh*100),'fontsize',20);
    set(gca,'xtick',[0:(24*3600)/G_bin:xend*3600/G_bin],...
        'xticklabel',[0:24:xend],'xlim',[xstart*3600/G_bin xend*3600/G_bin]);
    set(main_rsu_fig,'units','normalized','position',[0.2 0.1 0.75 0.75]);
    
    % 'baseline' indicator line
    if normalize
        line([0 size(meanFR_RSU,2)],[1 1],'color',cyel,'linestyle','--','linewidth',4);
        set(gca,'ylim',[0 ymax_main]);
    else
        set(gca,'ylim',[0 12]);
    end
    
    set(gca,'fontsize',22,'linewidth',2,'XColor','k','YColor','k');
    xlabel('Experiment hours','fontsize',24,'color','k');
    ylabel('Mean normalized firing rate','fontsize',26,'color','k');
    
end



%% PLOT *MEAN* FR BY ANIM
plot_meanby_anim = 1;
if plot_meanby_anim
    
    meanbyanim_fig = figure(110); hold on;
    set(meanbyanim_fig,'visible','off','units','normalized',...
        'position',[0.1 0.1 0.9 0.7]);
    title(sprintf('%s HEMISPHERE. mean FR by animal.',dep_status{dep+1}),...
        'fontsize',20);
    
    anims_online = unique(anims_bycell);
    for aa = 1:size(anims_online,2)
        if ~strcmp(anims_online{aa},'NONE')
            %             anim_colors{aa} = colorpool{aa};
            animcolor_idx{aa} = find(strcmp(anims_bycell,anims_online{aa}));
        end
    end
    
    if floor(numel(anims_online)/2) == 1
        n_subrows = numel(anims_online);
        n_subcols = 1;
    else
        n_subrows = ceil(numel(anims_online)/2);
        n_subcols = ceil(numel(anims_online)/n_subrows);
        if n_subrows-n_subcols > 1
            n_subcols = n_subcols + 1;
            n_subrows = ceil(numel(anims_online)/n_subcols);
        end
    end
    
    for cc = 1:size(FRbycell_RSU,1)
        this_color_idx = find(cellfun(@(x) ismember(cc,x), animcolor_idx));
        cells_bycolor(cc,1) = this_color_idx;
        clear this_color_idx
    end
    
    
    edgesD = (12*3600)/G_bin : (24*3600)/G_bin : xsize;
    for qq = 1:size(edgesD,2)
        q = rectangle('Position',[edgesD(qq), 0, 12*3600/G_bin, 15] );
        set(q,'facecolor',[0.9 0.9 0.9],'linestyle','none');
    end
    
    
    for xx = unique(cells_bycolor)'
        thesecells = find(cells_bycolor == xx);
        meanfr_thisanim = nanmean(FRbycell_RSU(thesecells,:),1);
        subplot(n_subrows,n_subcols,xx); hold on;
        plot(meanfr_thisanim,'linewidth',1.5,'color',colorpool{xx});
        line([0 size(meanfr_thisanim,2)],[1 1],'color',c_ctrl);
        xlabel(anims_bycell{thesecells(1)});
        set(gca,'ylim',[0 max(2,1.05*max(meanfr_thisanim))],'xtick',...
            [0:24*3600/G_bin:288*3600/G_bin],'xticklabel',[0:24:288]);
        clear thesecells meanfr_thisanim
    end
    
    
    set(meanbyanim_fig,'visible','on');
end

%% PLOT EACH CELL BY ANIMAL - RSUs
plot_by_anim = 1;
if plot_by_anim
    
    figure(120); hold on;
    
    title(sprintf('%s HEMISPHERE. FR by cell and by animal.',dep_status{dep+1}),...
        'fontsize',20);
    
    anims_online = unique(anims_bycell);
    for aa = 1:size(anims_online,2)
        if ~strcmp(anims_online{aa},'NONE')
            %             anim_colors{aa} = colorpool{aa};
            animcolor_idx{aa} = find(strcmp(anims_bycell,anims_online{aa}));
        end
    end
    
    if floor(numel(anims_online)/2) == 1
        n_subrows = numel(anims_online);
        n_subcols = 1;
    else
        n_subrows = ceil(numel(anims_online)/2);
        n_subcols = ceil(numel(anims_online)/n_subrows);
        if n_subrows-n_subcols > 1
            n_subcols = n_subcols + 1;
            n_subrows = ceil(numel(anims_online)/n_subcols);
        end
    end
    
    
    edgesD = (12*3600)/G_bin : (24*3600)/G_bin : xsize;
    for qq = 1:size(edgesD,2)
        q = rectangle('Position',[edgesD(qq), 0, 12*3600/G_bin, 15] );
        set(q,'facecolor',[0.9 0.9 0.9],'linestyle','none');
    end
    
    for cc = 1:size(FRbycell_RSU,1)
        disp(cc)
        this_color_idx = find(cellfun(@(x) ismember(cc,x), animcolor_idx));
        cells_bycolor(cc,1) = this_color_idx;
        subplot(n_subrows,n_subcols,this_color_idx); hold on;
        plot(FRbycell_RSU(cc,:),'linewidth',2)
        set(gcf,'position',[.1 .1 .8 .8]);
        clear this_color_idx
    end
    
    for pp = 1:numel(anims_online)
        subplot(n_subrows,n_subcols,pp);
        title(sprintf('%s, %s HEMISPHERE',anims_online{pp},dep_status{dep+1}),'fontsize',20);
        set(gca,'xtick',[0*3600/G_bin:(12*3600)/G_bin:size(meanFR_RSU,2)],...
            'fontsize',16,'xlim',[0*3600/G_bin size(meanFR_RSU,2)],...
            'xticklabel',{'','BL1','','BL2','','BL3','','MD1','','MD2','',...
            'MD3','','MD4','','ER1','','ER2','','ER3','','ER4','','ER5'});
        %             'xticklabel',{'MD4','','ER1','','ER2','','ER3','','ER4','','ER5'});
        if normalize
            set(gca,'ylim',[0 12]);
        else
            set(gca,'ylim',[0 20]);
        end
        if pp ==2, ylabel('Normalized firing rate','fontsize',18); end
    end
    
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




