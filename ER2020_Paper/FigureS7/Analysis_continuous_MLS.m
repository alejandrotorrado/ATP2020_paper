%    for dep=1:1
% clearvars -except CONTCELL_* recov*_analysis dep
%     clearvars -except CONTCELL* ER_DATA

%% SETUP
% make variable easier to use
dataset = 'SD_NEW';
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
%     c_ctrl = [5, 83, 48]./255;
c_ctrl = [0 0 0];
c_dep = [112,166,217]./255;
colorpool = {cblu,cred,cyel,cpur,cgre,ccya,cmar,cblk,c_ctrl};

% parameters for analysis
% deprived or control hemisphere
% quality threshold (2 for single units)
dep = 1;
dep_status = {'CONTROL','DEPRIVED'};
permut_drop = 0;

qthresh = 2;
doRSU = 1;
doPFS = 0;

% Choose continuity threshold, whether to normalize to baseline, and
% parameters to define baseline
plotcellsep = 0;
normalize = 1;
no_big_vals = 1;
big_val_thresh = 12;
no_small_vals = 1;
small_val_thresh = 0.01;
val_thresh_perc = 0.01;
percThresh = 0.75;
fprintf('Continuity threshold is: %u%%.\n\n',percThresh*100);
max_bl_change = 2;
min_bl_change = 1/max_bl_change;
require_BLON = 1;
blon_hour = 48; % this is ignored for RECOV data (where blon_hour is hard-coded to 156)
if normalize && ~require_BLON
    warning('WARNING: asked to normalize but not requiring cells to be ON at baseline may result in errors.');
end
conf_filter = 0;
conf_thresh = 0.7;

% for plotting main fig
ymax_main = 4;

% create arrays for indexing of cells
quality = [CELL.quality];
deprived = [CELL.deprived];
neg_pos_time = [CELL.neg_pos_time];
tail_slope = [CELL.tailSlope];
animalname = {CELL.animal};

% animals to examine
switch dataset
    case 'ER'
        if dep
            which_anim = {'KH67','AT14','AT16',...
                'KH73','AT29','AT22','AT25','AT29','KH75'};
        else
            which_anim = {'KH67','AT14','AT22','AT16','AT25',...
                'AT27','AT29','KH72','KH73'};
        end
%                     which_anim = {'AT12'};
    case 'CPP'
        which_anim = {'AT20','AT19'};
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
        %         which_anim = {'AT25'};
    case 'SD_recov'
        which_anim = {'AT35','AT36','AT37','AT38'};
    case 'CPPrecov'
        which_anim = {'AT19','AT20','AT23','AT32'};
    case 'CPP2recov'
        which_anim = {'AT50','AT51','AT52','AT53','AT54','AT56','AT57','AT58','AT59','AT60','AT61'};
%         which_anim = {'AT47'};
    case 'BLINJ'
        which_anim = {'AT56','AT57'};
    case 'SD_NEW'
        which_anim = {'AT63','AT64','AT65'};
end

%% REVIEW THIS
for dd = 1:length(CELL)
    totT = [];
    for uu = 1:length(CELL(dd).onTime)
        totT(uu) = CELL(dd).offTime(uu) - CELL(dd).onTime(uu);
    end
    
    switch dataset
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
            
        case 'CPP'
            daymax = 11.5;
            daysRecorded = daymax - CELL(dd).dayStart;
            
        case {'recov','CPPrecov','SD_recov','CPP2recov','SD_NEW'}
            % fixed number of recorded days, from MD4 to ER4 (included)
            daysRecorded = 5;
        case 'BLINJ'
            daysRecorded = 3;
    end
    
    
    if sum(totT) < (percThresh * daysRecorded*24*3600)
        percentOn(dd)        = 0;
    else
        if CELL(dd).dayStart > 2
            if sum(totT) < (percThresh * daysRecorded*24*3600)
                percentOn(dd) = 0;
            else
                percentOn(dd) = 1;
            end
        else
            %             if strcmp(CELL(dd).animal,'AT14')
            %                 fprintf('cell %u percent on\n',dd);
            %                 ongood1 = [ongood1; dd];
            %             end
            percentOn(dd) = 1;
        end
    end
    
    % this is old - REVIEW
    if require_BLON
        % cell must come online before start of "baseline" | or half of baseline?
        % "Baseline" is 36-60 hours for animals starting before MD1
        % "Baseline" is MD4 for AT12 (starting at MD3)
        switch dataset
            case {'ER','CPP'}
                timeon = CELL(dd).onTime(1) + CELL(dd).dayStart*24*3600;
                if timeon <= blon_hour*3600
                    %             if strcmp(CELL(dd).animal,'AT14')
                    %                 fprintf('cell %u on at: %.1f\n',dd,CELL(dd).onTime(1)/3600);
                    %                 ongood2 = [ongood2; dd];
                    %             end
                    onatBase(dd) = 1;
                elseif CELL(dd).onTime(1) <= 6.0*24*3600 && CELL(dd).dayStart > 2
                    onatBase(dd) = 1;
                else
                    onatBase(dd) = 0;
                end
                
            case {'recov','CPPrecov','CPP2recov'}
                timeon = CELL(dd).onTime(1);
                if timeon <= 6.5*24*3600
                    onatBase(dd) = 1;
                else
                    onatBase(dd) = 0;
                end
            case{'SD_recov','SD_NEW'}
                timeon = CELL(dd).onTime(1) + CELL(dd).dayStart*24*3600;
                if timeon <= 6.5*24*3600
                    onatBase(dd) = 1;
                else
                    onatBase(dd) = 0;
                end
            case 'BLINJ'
                timeon = CELL(dd).onTime(1);
                if timeon <= 3*24*3600
                    onatBase(dd) = 1;
                else
                    onatBase(dd) = 0;
                end
                
        end
    else
        onatBase(dd) = 1;
    end
    
    if conf_filter
        confscore = CELL(dd).score;
        if confscore(1) == 2 && confscore(3) == 3
            if confscore(2) < conf_thresh
                quality(dd) = 2.5;
            end
        end
    end
    
end

%% INDEXING OF CELLS


thisanim_idx = [];
tmp_anim = [];
for xx=1:length(which_anim)
    thisone = which_anim{xx};
    tmp_anim = [tmp_anim find(strcmp(thisone,animalname)==1)];
end
thisanim_idx = tmp_anim;
clear tmp_anim;

depidx = []; qualidx = []; percidx = []; idx0 = []; idx1 = [];
% deprived, quality, percentageOn, baselineOn
depidx  = find(deprived == dep);
qualidx = find(quality <= qthresh);
percidx = find(percentOn == 1);
blidx = find(onatBase == 1);
% intersect them
idx00 = intersect(qualidx,depidx);
idx01 = intersect(blidx,percidx);
idx0 = intersect(idx00,idx01);
% final result - cross with animals
idx1    = intersect(idx0,thisanim_idx);


% Cell classification into RSU/pFS cells
% Most recent negposthresh, tailslope and halfwidth should all be
% independent of interpolation factor
idx2 = []; idx3 = []; idx4 = []; idx5 = []; ns = [];

for cellno = idx1
    
    negpos_thresh = 0.39;
    fs_negpos_thresh = 0.39;
    tailSlope_thresh = 0.005;
    fs_tailSlope_thresh = 0.005;
    
    if CELL(cellno).neg_pos_time >= negpos_thresh
        idx2(end+1) = cellno;
    elseif CELL(cellno).neg_pos_time < fs_negpos_thresh
        idx3(end+1) = cellno;
    end
    
    if CELL(cellno).tailSlope >= tailSlope_thresh
        idx4(end+1) = cellno;
    elseif CELL(cellno).tailSlope < fs_tailSlope_thresh
        idx5(end+1) = cellno;
    end
end

%% PLOT CELL SEPARATION
cellsep_count = 0;
if plotcellsep
    
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
    
    
    cellsepfig = figure(222);
    set(cellsepfig,'visible','off','numbertitle','off');
    for cc = 1:cellsep_count
        s1 = subplot(1,2,1); hold on;
        plot(s1,cellsep_wf{cc},...
            'color',cellsep_color{cc},'linewidth',2);
        s2 = subplot(1,2,2); hold on;
        plot(s2,cellsep_nt(cc),cellsep_ts(cc),'o','color',cellsep_color{cc},...
            'MarkerFaceColor',cellsep_color{cc},'MarkerSize',10);
    end
    line([fs_negpos_thresh fs_negpos_thresh],[-0.1 fs_tailSlope_thresh],...
        'color','k','linestyle','--')
    line([0 fs_negpos_thresh],[fs_tailSlope_thresh fs_tailSlope_thresh],...
        'color','k','linestyle','--')
    line([negpos_thresh negpos_thresh],[tailSlope_thresh 0.1],...
        'color','k','linestyle','--')
    line([negpos_thresh 1],[tailSlope_thresh tailSlope_thresh],...
        'color','k','linestyle','--')
    set(s2,'ylim',[-0.1 0.1],'xlim',[0.2 0.60]);
    set(cellsepfig,'units','normalized','position',[0.1 0.1 0.7 0.7],'visible','on');
    set(s1,'fontsize',18,'linewidth',2);
    set(s2,'fontsize',18,'linewidth',2);
    s1.XAxis.Color = 'k';
    s2.XAxis.Color = 'k';
    
    
    
end


%% SETUP FOR FR ANALYSIS

RSU_idx = [];
RSU_idx = intersect(intersect(idx1,idx2),idx4);

pFS_idx = [];
pFS_idx = intersect(intersect(idx1,idx3),idx5);

fprintf('Found %u continuous RSUs that are on at BL and persist for %u%% of time.\n\n',...
    size(RSU_idx,2),percThresh*100);

fprintf('Found %u continuous pFS cells that are on at BL and persist for %u%% of time.\n\n',...
    size(pFS_idx,2),percThresh*100);

%% from permutation test analysis
if permut_drop
    if dep
        badz = [322, 384, 594, 605, 610];
    else
        badz = [157, 238, 454, 461, 485, 507];
    end
%     keyboard;
    RSU_idx = setdiff(RSU_idx,badz);
end

CELL_RSU = CELL(RSU_idx);
CELL_pFS = CELL(pFS_idx);

clear t_start t_stop max_h bin_sz bin_edges


% estimate FR using gaussian kernel. Create kernel here
G_sigma = 300; % s.dev. of Gaussian kernel (seconds) standard is 300!
G_bin = max(1,G_sigma/5); % bin size with which to estimate Gaussian kernel (seconds)
G_xvals = -5*G_sigma : G_bin : G_sigma*5; % values at which to estimate kernel
G_kernel = ( 1/(sqrt(2*pi)*G_sigma) ) *exp(-(G_xvals.^2)/(2*G_sigma^2)); % kernel creation


t_start = 0*24*3600;
t_stop  = 12*24*3600;

switch dataset
    case {'recov','CPPrecov','SD_recov','CPP2recov','SD_NEW'}
        % For normalization, get bin numbers corresponding to MD4 (night)
        % note that "recov" recordings are clustered to start on MD4 day
        bl_start = 6.5*24*3600/G_bin;
        bl_end   = 7*24*3600/G_bin;
    case {'ER','CPP','BLINJ'}
        % For normalization, get bin numbers corresponding to
        % BL2night/BL3day
        bl_start = 1.5*24*3600/G_bin;
        bl_end   = 2.5*24*3600/G_bin;
        
end
max_h = t_stop / 3600;


% get rid of artifact due to eye re-opening
switch dataset
    case {'ER','recov','SD_recov'}
        ER_start = (7*24*3600 + 2.0*3600)/G_bin; % 2.5 hours after start of ER1 (9:00 am)
        ER_end   = (7*24*3600 + 4.5*3600)/G_bin; % 5.0 hours after start of ER1 (12:00 pm)
    case 'CPP'
        ER_start = (5*24*3600 + 2.5*3600)/G_bin; % 2.5 hours after start of ER1 (10 am)
        ER_end   = (5*24*3600 + 5.0*3600)/G_bin; % 5.0 hours after start of ER1 (12:30 pm)
    case 'SD_NEW'
        ER_start = (7*24*3600 + 1.1*3600)/G_bin; % 2.5 hours after start of ER1 (10 am)
        ER_end   = (7*24*3600 + 2.2*3600)/G_bin; % 5.0 hours after start of ER1 (12:30 pm)
    case {'CPPrecov','CPP2recov'}
        ER_start = (7*24*3600 + 2.0*3600)/G_bin; % 2.5 hours after start of ER1 (10 am)
        ER_end   = (7*24*3600 + 4.5*3600)/G_bin; % 5.0 hours after start of ER1 (12:30 pm)
end
%% RSUs

if doRSU
    fprintf('\n\n      *** RSU cells. ***   \n\n');
    rsu_count = 0;
    for ee = 1:length(CELL_RSU)
        
        switch dataset
            case {'ER','CPP','SD_recov','CPP2recov','BLINJ','SD_NEW'}
                spikes  = CELL_RSU(ee).time + CELL_RSU(ee).dayStart*24*3600;
            case {'recov','CPPrecov'}
                spikes  = CELL_RSU(ee).time;
        end
%         keyboard;
        % Get spike raster in 1 sec bins. This will be used to estimate FR by
        % convolving with Gaussian kernel estimated with same bin size
        spike_raster = histc(spikes,t_start:G_bin:t_stop);
        % Use convolution to get mean FR estimate
        temp = conv(spike_raster,G_kernel,'same');
        temp(temp==0) = NaN; % get rid of 0 values (make them NaNs)
%         figure(); plot(temp);
        
        
        
        
        if nanmean(temp) < 50 && nanmean(temp) >= 0.001 % if rate is within "acceptable" range
            % increase counter
            
            % assign to new variable
            rate_temp = temp;
            
            % collect on/off times for cell
            ontimes = CELL_RSU(ee).onTime;
            offtimes = CELL_RSU(ee).offTime;
            timelist(:,1) = ontimes;
            timelist(:,2) = offtimes;
            
            % Process on/off times - REVIEW THIS
            try
                rate = processOnOffTimes_ATP(rate_temp,timelist,G_bin,CELL_RSU(ee).dayStart,dataset);
            catch
                keyboard;
            end
            
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
                case 'CPP'
                    if require_BLON
                        bl_start = 0*24*3600/G_bin;
                        bl_end = 1*24*3600/G_bin;
                        bl_rate = nanmean(rate(bl_start:bl_end));
                        if normalize
                            rate = rate./bl_rate;
                        end
                    else
                        bl_start = 4.5*24*3600/G_bin;
                        bl_end = ceil(5.00*24*3600/G_bin);
                        bl_rate = nanmean(rate(bl_start:bl_end));
                        if normalize
                            rate = rate./bl_rate;
                        end
                    end
                case {'recov','CPPrecov','SD_recov','CPP2recov','BLINJ','SD_NEW'}
                    bl_rate = nanmean(rate(bl_start:bl_end));
                    if normalize
                        rate = rate./bl_rate;
                    end
                    
            end
            
            
            if no_big_vals
                % get rid of crazy values
                if normalize
%                     if sum(rate > big_val_thresh) >= ceil(1e2/G_bin)
                    if sum(rate > big_val_thresh) / size(rate,1) >= val_thresh_perc
                        disp('toobig');
                        rate(:) = NaN;
                    end
                else
                    if sum((rate./bl_rate) > big_val_thresh) / size(rate,1) >= val_thresh_perc % OG way: >= ceil(1e2/G_bin)
                        disp('toobig');
                        rate(:) = NaN;
                    end
                end
            end
            
            if no_small_vals
                % get rid of crazy small values
                if normalize
                    if sum(rate<small_val_thresh) / size(rate,1) >= val_thresh_perc % OG way: ceil(1e4 / G_bin)
                        disp('teeny');
                        rate(:) = NaN;
                    end
                else
                    if sum((rate./bl_rate) < small_val_thresh) / size(rate,1)  >= val_thresh_perc % OG way: ceil(1e4 / G_bin)
                        %                 keyboard;
                        disp('teeny');
                        rate(:) = NaN;
                    end
                end
            end
%             keyboard;
            
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
            %     else
            % otherwise assign dummy to anims_bycell
            %         anims_bycell{ee} = 'NONE';
        end
        
        
    end
    
    
    %% RSU compile data
    blacklist = [];
    for xc = 1:size(blacklist,2)
        badcell = blacklist(xc);
        FRbycell_RSU(badcell,:) = NaN(1,size(FRbycell_RSU,2));
    end
    [nanrow,nancol] = find(~all(isnan(FRbycell_RSU),2));
    new_rsu_count = numel(unique(nanrow));
    
    
    meanFR_RSU = nanmean(FRbycell_RSU);
    semFR_RSU = nanstd(FRbycell_RSU,[],1)./sqrt(new_rsu_count);
    % get rid of ER artifact
    if any(strcmp(dataset,{'recov','CPPrecov'}))
        MD4_start = 6*24*3600/G_bin;
        meanFR_RSU(MD4_start:MD4_start+G_sigma/10) = NaN;
        semFR_RSU(MD4_start:MD4_start+G_sigma/10) = NaN;
    end
    if ~(strcmp(dataset,'BLINJ'))
        meanFR_RSU(ER_start:ER_end) = NaN;
        semFR_RSU(ER_start:ER_end) = NaN;
    end
    % also remove first few bins (artifactual rise due to Gaussian kernel
    % smoothing) - remove G_sigma/2
%     meanFR_RSU(1:300) = NaN;
%     meanFR_RSU(end-300:end) = NaN;
%     semFR_RSU(end-300:end) = NaN;
%     semFR_RSU(1:300) = NaN;
    xsize = size(meanFR_RSU,2);
% %     
% % %     
%                 recov_analysis.(dep_status{dep+1}).RSU_FRbycell = FRbycell_RSU;
%                 recov_analysis.(dep_status{dep+1}).RSU_idx      = good_idxs;
%                 recov_analysis.(dep_status{dep+1}).RSU_count    = new_rsu_count;
%                 recov_analysis.(dep_status{dep+1}).RSU_anims    = anims_bycell;
%                 recov_analysis.(dep_status{dep+1}).percThresh   = percThresh;
%     
end

%     end
%     recov_analysis.G_bin = G_bin;
%     if ismac
%         sfile = '/Users/atorrado/Desktop/MLS_DATA/recov_analysis_SD_NEW_May2020.mat';
%     elseif ispc
%         sfile = 'Z:\ATP_MAIN\DATA\Eye_Reopening\recov_DATA\recov_analysis_SD_NEW_May2020.mat';
%     end
%     save(sfile,'recov_analysis','-v7.3');
%     keyboard;
% % %     
    
%% pFSs
if doPFS
    fprintf('\n\n      *** pFS cells. ***   \n\n');
    pfs_count = 0;
    for ee = 1:length(CELL_pFS)
        
        
        spikes  = CELL_pFS(ee).time;
        % Get spike raster in 1 sec bins. This will be used to estimate FR by
        % convolving with Gaussian kernel estimated with same bin size
        spike_raster = histc(spikes,t_start:G_bin:t_stop);
        % Use convolution to get mean FR estimate
        temp = conv(spike_raster,G_kernel,'same');
        temp(temp==0) = NaN; % get rid of 0 values (make them NaNs)
        
        if nanmean(temp) < 50 && nanmean(temp) >= 0.001 % if rate is within "acceptable" range
            % increase counter
            
            % assign to new variable
            rate_temp = temp;
            
            % collect on/off times for cell
            ontimes = CELL_pFS(ee).onTime;
            offtimes = CELL_pFS(ee).offTime;
            timelist(:,1) = ontimes;
            timelist(:,2) = offtimes;
            
            % Process on/off times - REVIEW THIS
            try
                rate = processOnOffTimes_ATP(rate_temp,timelist,G_bin,CELL_pFS(ee).dayStart,dataset);
            catch
                keyboard;
            end
            
            switch dataset
                case 'ER'
                    if CELL_pFS(ee).dayStart >= 2
                        % for AT12, normalize to MD4
                        bl_start_MD = 6.5*24*3600/G_bin; % start of MD4 night
                        bl_end_MD   = 7*24*3600/G_bin; % end of MD4 night
                        bl_rate = nanmean(rate(bl_start_MD:bl_end_MD));
                        if normalize
                            rate = rate./bl_rate;
                        end
                    else
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
                case 'CPP'
                    if require_BLON
                        bl_start = 0*24*3600/G_bin;
                        bl_end = 1*24*3600/G_bin;
                        bl_rate = nanmean(rate(bl_start:bl_end));
                        if normalize
                            rate = rate./bl_rate;
                        end
                    else
                        bl_start = 4.5*24*3600/G_bin;
                        bl_end = ceil(5.00*24*3600/G_bin);
                        bl_rate = nanmean(rate(bl_start:bl_end));
                        if normalize
                            rate = rate./bl_rate;
                        end
                    end
                case 'recov'
                    bl_rate = nanmean(rate(bl_start:bl_end));
                    if normalize
                        rate = rate./bl_rate;
                    end
                    
            end
            
            
            if no_big_vals
                % get rid of crazy values
%                 if normalize && any(rate > big_val_thresh)
                if normalize && sum(rate > big_val_thresh) >= 180
                    rate(:) = NaN;
                elseif ~normalize && any((rate./bl_rate) > big_val_thresh)
                    rate(:) = NaN;
                end
            end
            
            if no_small_vals
                % get rid of crazy small values
                if normalize
                    if sum(rate<small_val_thresh) >= 180
                        rate(:) = NaN;
                    end
                elseif ~normalize
                    if sum((rate./bl_rate) < small_val_thresh) >= 180
                        %                 keyboard;
                        disp('teeny');
                        rate(:) = NaN;
                    end
                end
            end
            
            if ~all(isnan(rate))
                pfs_count = pfs_count + 1;
                % assign to main variable
                FRbycell_pFS(pfs_count,:) = rate;
                anims_bycell_pfs{pfs_count} = CELL_pFS(ee).animal;
                good_idxs_pfs(pfs_count) = pFS_idx(ee);
                %                     keyboard;
            end
            
            clear timelist
            %     else
            % otherwise assign dummy to anims_bycell
            %         anims_bycell{ee} = 'NONE';
        end
        
        
    end
    
    
    %% pFS compile data
    blacklist = [];
    for xc = 1:size(blacklist,2)
        badcell = blacklist(xc);
        FRbycell_pFS(badcell,:) = NaN(1,size(FRbycell_pFS,2));
    end
    [nanrow,nancol] = find(~all(isnan(FRbycell_pFS),2));
    new_pfs_count = numel(unique(nanrow));
    
    
    meanFR_pFS = nanmean(FRbycell_pFS);
    semFR_pFS = nanstd(FRbycell_pFS,[],1)./sqrt(new_pfs_count);
    % get rid of ER artifact
    if strcmp(dataset,'recov')
        MD4_start = 6*24*3600/G_bin;
        meanFR_pFS(MD4_start:MD4_start+30) = NaN;
        semFR_pFS(MD4_start:MD4_start+30) = NaN;
    end
    meanFR_pFS(ER_start:ER_end) = NaN;
    semFR_pFS(ER_start:ER_end) = NaN;
    % also remove first few bins (artifactual rise due to Gaussian kernel
    % smoothing) - remove G_sigma/2
    meanFR_pFS(1:300) = NaN;
    meanFR_pFS(end-300:end) = NaN;
    semFR_pFS(end-300:end) = NaN;
    semFR_pFS(1:300) = NaN;
    xsize_fs = size(meanFR_pFS,2);
    
%             recov_analysis_FS.(dep_status{dep+1}).pFS_FRbycell = FRbycell_pFS;
%             recov_analysis_FS.(dep_status{dep+1}).pFS_idx      = good_idxs_pfs;
%             recov_analysis_FS.(dep_status{dep+1}).pFS_count    = new_pfs_count;
%             recov_analysis_FS.(dep_status{dep+1}).pFS_anims    = anims_bycell_pfs;
%             recov_analysis_FS.(dep_status{dep+1}).percThresh   = percThresh;
    
    
end
% end
% recov_analysis_FS.G_bin = G_bin;
% sfile = '/Users/atorrado/Desktop/MLS_DATA/DecData/recov_analysis_FS.mat';
% save(sfile,'recov_analysis_FS','-v7.3');






%%
%%
%%

%% MAIN FR FIGURE RSU
plot_main_rsu = 1;
if plot_main_rsu && doRSU
    main_rsu_fig = figure();
    set(main_rsu_fig,'color','white');
    hold on;
    switch dataset
        case {'ER','MD'}
            xstart = 0;
            xend = 264;
        case {'recov','CPPrecov','SD_recov','CPP2recov','SD_NEW'}
            xstart = 156;
%             xend = xstart+daysRecorded*24;
            xend = 264;
        case 'BLINJ'
            xstart = 0;
            xend = 72;
    end
    
    %     subplot(2,1,1); hold on;
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
        case {'CPP','CPPrecov'}
            thiscolors = {c_ctrl,cmar};
        case {'SD_NEW'}
            thiscolors = {c_ctrl,cpur};
        case {'CPP2recov'}
            thiscolors = {c_ctrl,cgre};
        case 'BLINJ'
            thiscolors = {c_ctrl,cpur};
    end
    shadedErrorBar(1:xsize,meanFR_RSU,semFR_RSU,{'color',thiscolors{dep+1},'linewidth',2});
    title(sprintf('%s HEMISPHERE\nAverage normalized RSU firing rate (n=%u)\nContinuity threshold: %u%%',...
        dep_status{dep+1},new_rsu_count,percThresh*100),'fontsize',20);
    set(gca,'xtick',[0:(24*3600)/G_bin:xend*3600/G_bin],...
        'xticklabel',[0:24:xend],'xlim',[xstart*3600/G_bin xend*3600/G_bin]);
    %          set(gca,'xtick',[0:(12*3600)/G_bin:size(meanFR_RSU,2)],...
    %         'xticklabel',{'','BL1','','BL2','','BL3','','MD1','','MD2','',...
    %         'MD3','','MD4','','ER1','','ER2','','ER3','','ER4','','ER5'},...
    %         'fontsize',16,'xlim',[0 size(meanFR_RSU,2)]);
    set(main_rsu_fig,'units','normalized','position',[0.2 0.1 0.75 0.75]);
    
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

%% MAIN FIGURE pFS
plot_main_pfs = 1;
if plot_main_pfs && doPFS
    main_pfs_fig = figure();
    set(main_pfs_fig,'color','white');
    hold on;
    switch dataset
        case {'ER','MD'}
            xstart = 0;
        case 'recov'
            xstart = 144;
    end
    
    %     subplot(2,1,1); hold on;
    % make light/dark edges
    edgesD = (12*3600)/G_bin : (24*3600)/G_bin : xsize_fs;
    for qq = 1:size(edgesD,2)
        q = rectangle('Position',[edgesD(qq), 0, 12*3600/G_bin, 15] );
        set(q,'facecolor',[0.9 0.9 0.9],'linestyle','none');
    end
    % plot shaded error bar
    thiscolors = {c_ctrl,cblu};
    shadedErrorBar(1:xsize_fs,meanFR_pFS,semFR_pFS,{'color',thiscolors{dep+1},'linewidth',2});
    title(sprintf('%s HEMISPHERE\nAverage normalized RSU firing rate (n=%u)\nContinuity threshold: %u%%',...
        dep_status{dep+1},new_pfs_count,percThresh*100),'fontsize',20);
    set(gca,'xtick',[0:(24*3600)/G_bin:264*3600/G_bin],...
        'xticklabel',[0:24:288],'xlim',[xstart*3600/G_bin 264*3600/G_bin]);
    %          set(gca,'xtick',[0:(12*3600)/G_bin:size(meanFR_RSU,2)],...
    %         'xticklabel',{'','BL1','','BL2','','BL3','','MD1','','MD2','',...
    %         'MD3','','MD4','','ER1','','ER2','','ER3','','ER4','','ER5'},...
    %         'fontsize',16,'xlim',[0 size(meanFR_RSU,2)]);
    set(main_pfs_fig,'units','normalized','position',[0.2 0.1 0.75 0.75]);
    
    if normalize
        line([0 size(meanFR_pFS,2)],[1 1],'color',[0.3 .3 .3],'linestyle','--','linewidth',2.5);
        set(gca,'ylim',[0 ymax_main]);
    else
        set(gca,'ylim',[0 12]);
    end
    
    set(gca,'fontsize',20);
    xlabel('Experiment hours','fontsize',22);
    ylabel('Mean normalized firing rate','fontsize',24);
    
end

keyboard


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
    
    
    %     for pp = 1:numel(anims_online)
    %         subplot(n_subrows,n_subcols,pp);
    %         title(sprintf('%s, %s HEMISPHERE',anims_online{pp},dep_status{dep+1}),'fontsize',20);
    %         set(gca,'xtick',[0:(12*3600)/G_bin:size(meanFR_RSU,2)],...
    %             'xticklabel',{'','BL1','','BL2','','BL3','','MD1','','MD2','',...
    %             'MD3','','MD4','','ER1','','ER2','','ER3','','ER4','','ER5'},...
    %             'fontsize',16,'xlim',[0 size(meanFR_RSU,2)]);
    %         if normalize
    %             set(gca,'ylim',[0 8]);
    %         else
    %             set(gca,'ylim',[0 20]);
    %         end
    %         if pp ==2, ylabel('Normalized firing rate','fontsize',18); end
    %     end
    
    
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
        plot(FRbycell_RSU(cc,:),'linewidth',2)%,'color',colorpool{this_color_idx});
        set(gcf,'position',[.1 .1 .8 .8]);
%         uiwait(gcf);
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

%
% figure(1000); hold on;
% plot(meanFR_ctrl,'color',[.3 .3 .3],'linewidth',2);
% plot(meanFR_dep,'color',cblu,'linewidth',2);
% set(gca,'ylim',[0 2.5]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if doPFS
%% PLOT *MEAN* FR BY ANIM
plot_meanby_anim_fs = 1;
if plot_meanby_anim_fs
    
    meanbyanim_fig_fs = figure(220); hold on;
    set(meanbyanim_fig_fs,'visible','off','units','normalized',...
        'position',[0.1 0.1 0.9 0.7]);
    title(sprintf('%s HEMISPHERE. mean FR (FS cells) by animal.',dep_status{dep+1}),...
        'fontsize',20);
    
    anims_online_pfs = unique(anims_bycell_pfs);
    for aa = 1:size(anims_online_pfs,2)
        if ~strcmp(anims_online_pfs{aa},'NONE')
            %             anim_colors{aa} = colorpool{aa};
            animcolor_idx_fs{aa} = find(strcmp(anims_bycell_pfs,anims_online_pfs{aa}));
        end
    end
    
    if floor(numel(anims_online_pfs)/2) == 1
        n_subrows = numel(anims_online_pfs);
        n_subcols = 1;
    else
        n_subrows = ceil(numel(anims_online_pfs)/2);
        n_subcols = ceil(numel(anims_online_pfs)/n_subrows);
        if n_subrows-n_subcols > 1
            n_subcols = n_subcols + 1;
            n_subrows = ceil(numel(anims_online_pfs)/n_subcols);
        end
    end
    
    for cc = 1:size(FRbycell_pFS,1)
        this_color_idx = find(cellfun(@(x) ismember(cc,x), animcolor_idx_fs));
        cells_bycolor_fs(cc,1) = this_color_idx;
        clear this_color_idx
    end
    
    
    edgesD = (12*3600)/G_bin : (24*3600)/G_bin : xsize;
    for qq = 1:size(edgesD,2)
        q = rectangle('Position',[edgesD(qq), 0, 12*3600/G_bin, 15] );
        set(q,'facecolor',[0.9 0.9 0.9],'linestyle','none');
    end
    
    
    for xx = unique(cells_bycolor_fs)'
        thesecells = find(cells_bycolor_fs == xx);
        meanfr_thisanim = nanmean(FRbycell_pFS(thesecells,:),1);
        subplot(n_subrows,n_subcols,xx); hold on;
        plot(meanfr_thisanim,'linewidth',1.5,'color',colorpool{xx});
        line([0 size(meanfr_thisanim,2)],[1 1],'color',c_ctrl);
        xlabel(anims_bycell_pfs{thesecells(1)});
        set(gca,'ylim',[0 max(2,1.05*max(meanfr_thisanim))],'xtick',...
            [0:24*3600/G_bin:288*3600/G_bin],'xticklabel',[0:24:288]);
        clear thesecells meanfr_thisanim
    end
    
    
    set(meanbyanim_fig_fs,'visible','on');
    
    
    %     for pp = 1:numel(anims_online)
    %         subplot(n_subrows,n_subcols,pp);
    %         title(sprintf('%s, %s HEMISPHERE',anims_online{pp},dep_status{dep+1}),'fontsize',20);
    %         set(gca,'xtick',[0:(12*3600)/G_bin:size(meanFR_RSU,2)],...
    %             'xticklabel',{'','BL1','','BL2','','BL3','','MD1','','MD2','',...
    %             'MD3','','MD4','','ER1','','ER2','','ER3','','ER4','','ER5'},...
    %             'fontsize',16,'xlim',[0 size(meanFR_RSU,2)]);
    %         if normalize
    %             set(gca,'ylim',[0 8]);
    %         else
    %             set(gca,'ylim',[0 20]);
    %         end
    %         if pp ==2, ylabel('Normalized firing rate','fontsize',18); end
    %     end
    
    
end
%% PLOT EACH CELL BY ANIMAL - RSUs
plot_by_anim_fs = 1;
if plot_by_anim_fs
    
    figure(100); hold on;
    title(sprintf('%s HEMISPHERE. FR by cell and by animal.',dep_status{dep+1}),...
        'fontsize',20);
    
    anims_online_fs = unique(anims_bycell_pfs);
    for aa = 1:size(anims_online_fs,2)
        if ~strcmp(anims_online_fs{aa},'NONE')
            %             anim_colors{aa} = colorpool{aa};
            animcolor_idx{aa} = find(strcmp(anims_bycell_pfs,anims_online_fs{aa}));
        end
    end
    
    if floor(numel(anims_online_fs)/2) == 1
        n_subrows = numel(anims_online_fs);
        n_subcols = 1;
    else
        n_subrows = ceil(numel(anims_online_fs)/2);
        n_subcols = ceil(numel(anims_online_fs)/n_subrows);
        if n_subrows-n_subcols > 1
            n_subcols = n_subcols + 1;
            n_subrows = ceil(numel(anims_online_fs)/n_subcols);
        end
    end
    
    
    edgesD = (12*3600)/G_bin : (24*3600)/G_bin : xsize;
    for qq = 1:size(edgesD,2)
        q = rectangle('Position',[edgesD(qq), 0, 12*3600/G_bin, 15] );
        set(q,'facecolor',[0.9 0.9 0.9],'linestyle','none');
    end
    
    for cc = 1:size(FRbycell_pFS,1)
        
        this_color_idx = find(cellfun(@(x) ismember(cc,x), animcolor_idx));
        cells_bycolor_fs(cc,1) = this_color_idx;
        subplot(n_subrows,n_subcols,this_color_idx); hold on;
        plot(FRbycell_pFS(cc,:),'linewidth',1,'color',colorpool{this_color_idx});
        clear this_color_idx
    end
    
    for pp = 1:numel(anims_online_fs)
        subplot(n_subrows,n_subcols,pp);
        title(sprintf('%s, %s HEMISPHERE',anims_online{pp},dep_status{dep+1}),'fontsize',20);
        set(gca,'xtick',[0:(12*3600)/G_bin:size(meanFR_RSU,2)],...
            'xticklabel',{'','BL1','','BL2','','BL3','','MD1','','MD2','',...
            'MD3','','MD4','','ER1','','ER2','','ER3','','ER4','','ER5'},...
            'fontsize',16,'xlim',[0 size(meanFR_RSU,2)]);
        if normalize
            set(gca,'ylim',[0 8]);
        else
            set(gca,'ylim',[0 20]);
        end
        if pp ==2, ylabel('Normalized firing rate','fontsize',18); end
    end
    
    
end



%
% figure(1000); hold on;
% plot(meanFR_ctrl,'color',[.3 .3 .3],'linewidth',2);
% plot(meanFR_dep,'color',cblu,'linewidth',2);
% set(gca,'ylim',[0 2.5]);
end
