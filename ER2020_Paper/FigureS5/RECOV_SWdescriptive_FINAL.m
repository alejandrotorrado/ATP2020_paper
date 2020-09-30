%% RECOV_SWdescriptive
%
% Alejandro Torrado Pacheco - 2017
%
% Plot descriptive features of sleep/wake states. This can be used to
% reproduce Figure S6 in the thesis. This is a corrected version of Fig. S5
% in the bioRxiv preprint. 


clearIDE

% load CONTCELL variable
if ismac
    contcell_file = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig1/CONTCELL_recov.mat';
elseif ispc
    contcell_file = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig1\CONTCELL_recov.mat';
end
tic;
load(contcell_file);
toc;

% comment all above and uncomment this line to debug without having to load
% CONTCELL every time
% clearvars -except CONTCELL*

% load recov data
if ismac
    loadfile = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig1/recov_analysis.mat';
%     loadfile = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Eye_Reopening/recov_DATA/recov_analysis_Jan20_nonorm_v3.mat';
elseif ispc
    loadfile = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig1\recov_analysis.mat';
end
rload = load(loadfile);
recov = rload.recov_analysis;

% define colors for plotting
cred = [0.85 0.33 0.01];
cblu = [0 0.45 0.74];
cyel = [0.93 0.69 0.13];
cpur = [0.49 0.18 0.55];
cgre = [0.47 0.67 0.19];
cg = [0, 210, 120]./255;
ccya = [0.30 0.75 0.93];
cmar = [0.64 0.08 0.18];
c_ctrl = [5, 83, 48]./255;
cgray = [.45 .45 .45];
clight = [160, 156, 0]./255;
cdark = [13, 7, 66]./255;
cblk = [0 0 0];
c_rem   = [25 181 149]./255;
c_nrem  = [131 49 146]./255;
c_aw    = [201 28 101]./255;
c_qw    = [247 148 41]./255;
s4colors = {c_nrem,c_rem,c_aw,c_qw};
s4colors_2 = {c_aw,c_qw,c_nrem,c_rem};
colorpool = {cblk,ccya,cmar,cyel,cpur,cgre,cblu,cred};
dc1 = [8, 219, 139]./255;
dc2 = [16, 164, 213]./255;
dc3 = [255, 93, 10]./255;
dc4 = [255, 217, 10]./255;
deucolorpool = {cblk,dc1,dc2,cblk,dc3,dc4};

% extract data
G_bin = recov.G_bin;

% get indices of re-opened and control cells used in main analysis
depcells = recov.DEPRIVED.RSU_idx;
ctrcells = recov.CONTROL.RSU_idx;
goodcells = [depcells';ctrcells'];
cellstat  = [ones(size(depcells'));zeros(size(ctrcells'))];
% find corresponding entries in MASTER variable
CELL_RSU = CONTCELL_recov.MASTER(goodcells);
% also get the statetimes
STATES = CONTCELL_recov.STATETIMES;

% toggles to choose which cells to examine. To reproduce Fig S6, set
% no_deps = 1 and no_ctrs = 0 (only look at control neurons)
no_deps = 1;
no_ctrs = 0;

%% ANALYSIS SETUP
% threshold for epoch rejection - only look at FR in epochs longer than
% epoch_thresh
epoch_thresh = 30; % seconds

% Days to examine - only want control/baseline days for this
% re-opened hemisphere neurons
dep_day_start = 6*24; % 6 = start of MD4
dep_start = dep_day_start.*3600;
dep_day_end = 7*24; % 7 = start of ER1
dep_end = dep_day_end.*3600;
% control hemisphere neurons
ctrl_day_start  = 6*24; % 6 = start of MD4
ctrl_start = ctrl_day_start.*3600;
ctrl_day_end    = 10*24; % 10 = start of ER4
ctrl_end = ctrl_day_end.*3600;

% settings
use_dayStart = 0; % this should be 0!
plot_states = 1; % this should be 1 for plotting

% initialize variables
REM_meanFRbyepoch = []; NREM_meanFRbyepoch = []; AW_meanFRbyepoch = []; QW_meanFRbyepoch = [];
REM_meanCVbyepoch = []; NREM_meanCVbyepoch = []; AW_meanCVbyepoch = []; QW_meanCVbyepoch = [];
REM_meanDUR_byepoch = []; NREM_meanDUR_byepoch = []; AW_meanDUR_byepoch = []; QW_meanDUR_byepoch = [];


%% MAIN ANALYSIS LOOP

% initialize counter
cc=0;
for uu = 1:length(goodcells)
    
    % loop through good cells
    fprintf('\nState analysis, cell %u of %u.\n',uu,length(CELL_RSU));
    fprintf('Animal: %s.\n',CELL_RSU(uu).animal);
    
    % get the statetimes for this animal - if animal is not in the STATES
    % structure, ignore it
    try
        st = STATES.(CELL_RSU(uu).animal);
        go_ahead = 1;
    catch
        if any(strcmp(CELL_RSU(uu).animal,{'AT22','AT25','AT24','KH72'}))
            go_ahead = 0;
            disp('NOPE!');
        else
            keyboard;
        end
    end
    
    % apply the settings choosing dep and ctrl neurons
    if no_deps
        if cellstat(uu), go_ahead = 0; end
    end
    if no_ctrs
        if ~cellstat(uu), go_ahead = 0; end
    end
    
    % if everything is OK, proceed
    if go_ahead
        
        % adjust statetimes to hypothetical BL1
        daystart = STATES.DAYSTART.(CELL_RSU(uu).animal);
        statetimes = adjust_statetimes_recov(st,CELL_RSU(uu),daystart);
        
        
        % select states during period of interest
        if cellstat(uu) == 1
            fprintf('deprived cell\n');
            for rr = 1:size(dep_start,1)
                homeo_states = statetimes( statetimes(:,2) >= dep_start(rr,1) & statetimes(:,2) < dep_end(rr,1) , : ) ;
            end
        elseif cellstat(uu) == 0
            fprintf('control cell\n');
            homeo_states = statetimes( statetimes(:,2) >= ctrl_start & statetimes(:,2) < ctrl_end , : ) ;
        end
        
        % if there are any
        if ~isempty(homeo_states)
            % updated counter
            cc = cc + 1;
            
            % use custom function to calculate mean FR, CV of ISI and epoch
            % duration, by state
            [REM_mean_cell, NREM_mean_cell, AW_mean_cell, QW_mean_cell,...
                REM_mean_epochs, NREM_mean_epochs, AW_mean_epochs, QW_mean_epochs, ...
                REM_cvcell, NREM_cvcell, AW_cvcell, QW_cvcell,...
                REM_cvepoch, NREM_cvepoch, AW_cvepoch, QW_cvepoch,...
                REM_durcell, NREM_durcell, AW_durcell, QW_durcell,...
                REM_durepoch, NREM_durepoch, AW_durepoch, QW_durepoch, zerocount] = ...
                SW_meanFR_4state( CELL_RSU(uu), homeo_states, epoch_thresh, use_dayStart);
            
            % Split into corresponding variables
            
            % FR
            REM_meanFRbycell(cc,1)   = REM_mean_cell;
            REM_meanFRbyepoch        = [REM_meanFRbyepoch; REM_mean_epochs];
            NREM_meanFRbycell(cc,1)  = NREM_mean_cell;
            NREM_meanFRbyepoch       = [NREM_meanFRbyepoch; NREM_mean_epochs];
            AW_meanFRbycell(cc,1)    = AW_mean_cell;
            AW_meanFRbyepoch         = [AW_meanFRbyepoch; AW_mean_epochs];
            QW_meanFRbycell(cc,1)    = QW_mean_cell;
            QW_meanFRbyepoch         = [QW_meanFRbyepoch; QW_mean_epochs];
            
            
            % CV
            REM_meanCVbycell(cc,1)   = REM_cvcell;
            REM_meanCVbyepoch        = [REM_meanCVbyepoch; REM_cvepoch];
            NREM_meanCVbycell(cc,1)  = NREM_cvcell;
            NREM_meanCVbyepoch       = [NREM_meanCVbyepoch; NREM_cvepoch];
            AW_meanCVbycell(cc,1)    = AW_cvcell;
            AW_meanCVbyepoch         = [AW_meanCVbyepoch; AW_cvepoch];
            QW_meanCVbycell(cc,1)    = QW_cvcell;
            QW_meanCVbyepoch         = [QW_meanCVbyepoch; QW_cvepoch];
            
            % EPOCH DURATION
            REM_meanDUR_bycell(cc,1)   = REM_durcell;
            REM_meanDUR_byepoch        = [REM_meanDUR_byepoch; REM_durepoch];
            NREM_meanDUR_bycell(cc,1)  = NREM_durcell;
            NREM_meanDUR_byepoch       = [NREM_meanDUR_byepoch; NREM_durepoch];
            AW_meanDUR_bycell(cc,1)    = AW_durcell;
            AW_meanDUR_byepoch         = [AW_meanDUR_byepoch; AW_durepoch];
            QW_meanDUR_bycell(cc,1)    = QW_durcell;
            QW_meanDUR_byepoch         = [QW_meanDUR_byepoch; QW_durepoch];
            
            
        end
        
    end
    
end

%% compile the data

% compile data and choose what to normalize to
norm_to = (AW_meanFRbycell); % usually normalize to AW
% compiled data without normalization
barFRdata_nonorm = [AW_meanFRbycell, QW_meanFRbycell, NREM_meanFRbycell, REM_meanFRbycell];
% compiled data with normalization
barFRdata = barFRdata_nonorm ./ norm_to;
% compute mean for both datasets
barFRmean = nanmean(barFRdata);
barFRmean_nonorm = nanmean(barFRdata_nonorm);
% compute sem for both datasets
barFRsem  = nanstd(barFRdata)./sqrt(sum(~isnan(barFRdata))-1);
barFRsem_nonorm  = nanstd(barFRdata_nonorm)./sqrt(sum(~isnan(barFRdata_nonorm))-1);

% compile CV data and compute mean and SEM
% not normalized
barCVdata_nonorm = [AW_meanCVbycell, QW_meanCVbycell, NREM_meanCVbycell, REM_meanCVbycell];
% normalized
barCVdata = barCVdata_nonorm ./ AW_meanCVbycell;
% mean for each
barCVmean = nanmean(barCVdata);
barCVmean_nonorm = nanmean(barCVdata_nonorm);
% SEM for each
barCVsem  = nanstd(barCVdata)./sqrt(sum(~isnan(barCVdata))-1);
barCVsem_nonorm  = nanstd(barCVdata_nonorm)./sqrt(sum(~isnan(barCVdata_nonorm))-1);

%% PLOT FR DATA

% sort the data for ladder plots
plot_dat = sortrows(barFRdata,1);
plot_dat_nonorm = sortrows(barFRdata_nonorm,1);

% STATS
% Pairwise Wilocoxon sign-rank tests with Bonferroni correction
% Using this test because data are paired, not independent. Each cell gives
% a mean FR value in each state.
% Syntax: signrank(x,y,'tail',tail)
% Compares medians of x and y distributions where samples are paired.
% Default is two-tailed (tail='both'). If set tail='right', testing that
% x>y. If tail='left', testing that x<y.

% AW v QW
p_rank(1) = signrank(plot_dat_nonorm(:,1),plot_dat_nonorm(:,2),'tail','both');
% AW v NREM
p_rank(2) = signrank(plot_dat_nonorm(:,1),plot_dat_nonorm(:,3),'tail','both');
% AW v REM
p_rank(3) = signrank(plot_dat_nonorm(:,1),plot_dat_nonorm(:,4),'tail','both');
% QW v NREM
p_rank(4) = signrank(plot_dat_nonorm(:,2),plot_dat_nonorm(:,3),'tail','both');
% QW v REM
p_rank(5) = signrank(plot_dat_nonorm(:,2),plot_dat_nonorm(:,4),'tail','both');
% NREM v REM
p_rank(6) = signrank(plot_dat_nonorm(:,3),plot_dat_nonorm(:,4),'tail','both');
 
% Bonferroni correction - multiply p-vals by # of comparisons
n_comparisons = max(size(p_rank));
p_rank = p_rank.* n_comparisons;
p_rank(p_rank>1) = 1; % p-values above 1 are just set to 1

% plotting
dotplot = figure();
set(dotplot,'color','w','units','normalized','position',[.1 .1 .9 .8]);
msz = 190; % markersize
lw = 4; % linewidth
mw = 0.2; % mean line half length
capsz = 17; % capsize for errorbars
m_alpha = 0.8; % scatter plot dot transparency
% axes 1 - no_norm
daxId(1) = axes('units','normalized','position',[.58 .1 .37 .8]);
for pp=1:4
    axes(daxId(1)); hold on;
    datacol = plot_dat_nonorm(:,pp);
    
    % choose appropriate color
    switch pp
        case 1
            this_col = c_aw;
        case 2
            this_col = c_qw;
        case 3
            this_col = c_nrem;
        case 4
            this_col = c_rem;
    end
    
    % scatter plot
    scat{pp} = scatter(repmat(pp,size(datacol)),datacol,msz,'o','MarkerFaceColor',...
        this_col,'markeredgecolor','none');
    
    % mean value line
    line([pp-mw pp+mw],[barFRmean_nonorm(pp) barFRmean_nonorm(pp)],'linewidth',lw,'color','k');
    
end
% errorbars for means
errorbar([1:size(plot_dat_nonorm,2)],barFRmean_nonorm,barFRsem_nonorm,'k',...
        'linestyle','none','linewidth',lw,'capsize',capsz);
% format axes
set(daxId(1),'xlim',[0.5 4.5],'yscale','log','fontsize',22,'xtick',[1:4],...
    'xticklabel',{'Active','Quiet','NREM','REM'},'XColor','k','YColor','k',...
    'ylim',[.01 100],'yticklabel',{'0.01','0.1','1','10','100'},'linewidth',2);
ylabel('Firing rate (Hz)','fontsize',26,'color','k');

% draw connector lines
for pt = 1:size(plot_dat_nonorm,1)
    connector{pt} = plot([1:4],plot_dat_nonorm(pt,1:4),'color',[cblk],...
        'linewidth',2);
    uistack(connector{pt},'bottom');
end


% axes 2 - normalized to AW
daxId(2) = axes('units','normalized','position',[.1 .1 .37 .8]);
for pp=1:4
    axes(daxId(2)); hold on;
    datacol = plot_dat(:,pp);
    
    switch pp
        case 1
            this_col = c_aw;
        case 2
            this_col = c_qw;
        case 3
            this_col = c_nrem;
        case 4
            this_col = c_rem;
    end
    
    scat{4+pp} = scatter(repmat(pp,size(datacol)),datacol,msz,'o','MarkerFaceColor',...
        this_col,'markeredgecolor','none');
    
    line([pp-mw pp+mw],[barFRmean(pp) barFRmean(pp)],'linewidth',lw,'color','k');

end
errorbar([1:size(plot_dat,2)],barFRmean,barFRsem,'k',...
        'linestyle','none','linewidth',lw,'capsize',capsz);
line([0 size(plot_dat,2)+1],[1 1],'linestyle','--','color',[.65 .65 .65],'linewidth',lw/2)
set(daxId(2),'xlim',[0.5 4.5],'yscale','linear','fontsize',22,'xtick',[1:4],...
    'xticklabel',{'Active','Quiet','NREM','REM'},'XColor','k','YColor','k',...
    'ylim',[0 2],'ytick',0:0.5:10,'linewidth',2);
ylabel('FR normalized to Active Wake','fontsize',26,'color','k');

% draw connector lines
for pt = 1:size(plot_dat,1)
    connector{pt} = plot([1:4],plot_dat(pt,1:4),'color',[cblk],...
        'linewidth',2);
    uistack(connector{pt},'bottom');
end

%% PLOT CV DATA

% sort the data for ladder plots
plot_cv = sortrows(barCVdata,1);
plot_cv_nonorm = sortrows(barCVdata_nonorm,1);

% STATS
% Pairwise Wilocoxon sign-rank tests with Bonferroni correction
% Using this test because data are paired, not independent. Each cell gives
% a mean FR value in each state.
% Syntax: signrank(x,y,'tail',tail)
% Compares medians of x and y distributions where samples are paired.
% Default is two-tailed (tail='both'). If set tail='right', testing that
% x>y. If tail='left', testing that x<y.

% AW v QW
p_cv_rank(1) = signrank(plot_cv_nonorm(:,1),plot_cv_nonorm(:,2),'tail','both');
% AW v NREM
p_cv_rank(2) = signrank(plot_cv_nonorm(:,1),plot_cv_nonorm(:,3),'tail','both');
% AW v REM
p_cv_rank(3) = signrank(plot_cv_nonorm(:,1),plot_cv_nonorm(:,4),'tail','both');
% QW v NREM
p_cv_rank(4) = signrank(plot_cv_nonorm(:,2),plot_cv_nonorm(:,3),'tail','both');
% QW v REM
p_cv_rank(5) = signrank(plot_cv_nonorm(:,2),plot_cv_nonorm(:,4),'tail','both');
% NREM v REM
p_cv_rank(6) = signrank(plot_cv_nonorm(:,3),plot_cv_nonorm(:,4),'tail','both');
 
% Bonferroni correction - multiply p-vals by # of comparisons
n_comparisons = max(size(p_cv_rank));
p_cv_rank = p_cv_rank.* n_comparisons;
p_cv_rank(p_cv_rank>1) = 1; % p-values above 1 are just set to 1

% plotting
dotcvplot = figure();
set(dotcvplot,'color','w','units','normalized','position',[.1 .1 .9 .8]);
msz = 190; % markersize
lw = 4; % linewidth
mw = 0.2; % mean line half length
capsz = 17; % capsize for errorbars
m_alpha = 0.8; % scatter plot dot transparency
% axes 1 - no_norm
cvaxId(1) = axes('units','normalized','position',[.58 .1 .37 .8]);
for pp=1:4
    axes(cvaxId(1)); hold on;
    datacol = plot_cv_nonorm(:,pp);
    
    % choose appropriate color
    switch pp
        case 1
            this_col = c_aw;
        case 2
            this_col = c_qw;
        case 3
            this_col = c_nrem;
        case 4
            this_col = c_rem;
    end
    
    % scatter plot
    scat{pp} = scatter(repmat(pp,size(datacol)),datacol,msz,'^','MarkerFaceColor',...
        this_col,'markeredgecolor','none');
    
    % mean value line
    line([pp-mw pp+mw],[barCVmean_nonorm(pp) barCVmean_nonorm(pp)],'linewidth',lw,'color','k');
    
end
% errorbars for means
errorbar([1:size(plot_cv_nonorm,2)],barCVmean_nonorm,barCVsem_nonorm,'k',...
        'linestyle','none','linewidth',lw,'capsize',capsz);
% format axes
set(cvaxId(1),'xlim',[0.5 4.5],'yscale','linear','fontsize',22,'xtick',[1:4],...
    'xticklabel',{'Active','Quiet','NREM','REM'},'XColor','k','YColor','k',...
    'ylim',[0.5 2],'ytick',0:.5:10,'linewidth',2);
ylabel('CV of ISI)','fontsize',26,'color','k');

% draw connector lines
for pt = 1:size(plot_cv_nonorm,1)
    connector{pt} = plot([1:4],plot_cv_nonorm(pt,1:4),'color',[cblk],...
        'linewidth',2);
    uistack(connector{pt},'bottom');
end


% axes 2 - normalized to AW
cvaxId(2) = axes('units','normalized','position',[.1 .1 .37 .8]);
for pp=1:4
    axes(cvaxId(2)); hold on;
    datacol = plot_cv(:,pp);
    
    switch pp
        case 1
            this_col = c_aw;
        case 2
            this_col = c_qw;
        case 3
            this_col = c_nrem;
        case 4
            this_col = c_rem;
    end
    
    scat{4+pp} = scatter(repmat(pp,size(datacol)),datacol,msz,'^','MarkerFaceColor',...
        this_col,'markeredgecolor','none');
    
    line([pp-mw pp+mw],[barCVmean(pp) barCVmean(pp)],'linewidth',lw,'color','k');

end
errorbar([1:size(plot_cv,2)],barCVmean,barCVsem,'k',...
        'linestyle','none','linewidth',lw,'capsize',capsz);
line([0 size(plot_cv,2)+1],[1 1],'linestyle','--','color',[.65 .65 .65],'linewidth',lw/2)
set(cvaxId(2),'xlim',[0.5 4.5],'yscale','linear','fontsize',22,'xtick',[1:4],...
    'xticklabel',{'Active','Quiet','NREM','REM'},'XColor','k','YColor','k',...
    'ylim',[0.75 1.5],'ytick',0:0.25:10,'linewidth',2);
ylabel('CV of ISI normalized to Active Wake','fontsize',26,'color','k');

% draw connector lines
for pt = 1:size(plot_cv,1)
    connector{pt} = plot([1:4],plot_cv(pt,1:4),'color',[cblk],...
        'linewidth',2);
    uistack(connector{pt},'bottom');
end



%% Bar plot for FR and CV - not used in manuscript
h1 = figure(); 
set(h1,'color','w','unit','normalized','position',[.1 .1 .6 .7]);


% axes 1 - FR
h1_a = axes(h1,'position',[.08 .1 .4 .85]); hold on;

for uu =1:4
    
    barFR{uu} = bar(uu,barFRmean(uu),.5);
    barFR{uu}.FaceColor = s4colors_2{uu};
    barFR{uu}.FaceAlpha = 1;
    barFR{uu}.EdgeColor = 'k';
    barFR{uu}.LineWidth = 2;
    barFR{uu}.LineStyle = '-';
    errFR{uu} = errorbar(uu,barFRmean(uu),barFRsem(uu),'linestyle','none',...
        'linewidth',2,'color','k','capsize',0);
end
hold off;


% Kruskal-Wallis test (not used since data violates independence)
%{
kwdata = reshape(barFRdata_nonorm,[],1);
kw_gr  = [ones(size(barFRdata,1),1);2*ones(size(barFRdata,1),1);...
    3*ones(size(barFRdata,1),1);4*ones(size(barFRdata,1),1)];
[p_kw,~,kwstats] = kruskalwallis(kwdata,kw_gr,'off');
p_mult = multcompare(kwstats,'display','off');
%}

%
h1_b = axes(h1,'position',[.58 .1 .4 .85]); hold on;
for uu =1:4
    barCV{uu} = bar(uu,barCVmean(uu),.5);
    barCV{uu}.FaceColor = s4colors_2{uu};
    barCV{uu}.FaceAlpha = 1;
    barCV{uu}.EdgeColor = 'k';
    barCV{uu}.LineWidth = 2;
    barCV{uu}.LineStyle = '-';
    errCV{uu} = errorbar(uu,barCVmean(uu),barCVsem(uu),'linestyle','none',...
        'linewidth',2,'color','k','capsize',0);
end

s_names = {'Active','Quiet','NREM','REM'};
set(h1_a,'xlim',[.5 4.5],'xtick',1:4,'xticklabel',s_names,'fontsize',14,...
    'XColor','k','YColor','k','LineWidth',1.5);
axes(h1_a); ylabel('Normalized mean firing rate','fontsize',18,'color','k');
set(h1_b,'xlim',[.5 4.5],'xtick',1:4,'xticklabel',s_names,'fontsize',14,...
    'XColor','k','YColor','k','LineWidth',1.5);
axes(h1_b); ylabel('Mean CV of ISI','fontsize',18,'color','k');


%{
kwcv = reshape(barCVdata,[],1);
kw_cvgr  = [ones(size(barCVdata,1),1);2*ones(size(barCVdata,1),1);...
    3*ones(size(barCVdata,1),1);4*ones(size(barCVdata,1),1)];
[p_kwcv,~,kwstatscv] = kruskalwallis(kwcv,kw_cvgr,'off');
p_multcv = multcompare(kwstatscv,'display','off');
%}

%% epoch duration
dur_data = {AW_meanDUR_byepoch, QW_meanDUR_byepoch, NREM_meanDUR_byepoch, REM_meanDUR_byepoch};
h2 = figure(); hold on;
set(h2,'color','w','unit','normalized','position',[.1 .1 .5 .7]);
for tt = 1:4
    semdur = nanstd(dur_data{tt})./sqrt(numel(dur_data{tt}));
    meandur = nanmean(dur_data{tt});
    durbar{tt} = bar(tt,meandur,.5);
    durbar{tt}.FaceColor = s4colors_2{tt};
    durbar{tt}.EdgeColor = 'k';
    durbar{tt}.FaceAlpha = 1;
    durbar{tt}.LineWidth = 2;
    durerr{tt} = errorbar(tt,meandur,semdur,'linestyle','none',...
        'linewidth',2,'capsize',0,'color',[0 0 0]);
end
set(gca,'linewidth',2,'fontsize',14,'XColor','k','YColor','k','xtick',1:4,...
    'xticklabel',s_names);
ylabel('Mean epoch duration (seconds)','fontsize',16);

%% duration histograms - not plotted in manuscript
gg = figure();
set(gg,'color','w','unit','normalized','position',[.1 .1 .7 .8]);


bin_edges = 0:25:2000;

nAxes = 4;
nRows = 2;
nCols = 2;

colW = .75/nCols; rowH = .7/nRows;
colX = .1 + linspace(0,0.96,nCols+1); colX = colX(1:end-1);
rowY = .1 + linspace(.96,0,nRows+1); rowY = rowY(2:end);

for dId = 1:nAxes
    rowId = ceil( dId / nCols ) ;
    colId = dId - (rowId - 1) * nCols ;
    h_ax{dId} = axes('position',[colX(colId) rowY(rowId) colW rowH]); 
    
    histdur{dId} = histogram(dur_data{dId},bin_edges);
    histdur{dId}.FaceColor = s4colors_2{dId};
    histdur{dId}.EdgeColor = 'none';
    histdur{dId}.FaceAlpha = 1;
    set(gca,'xlim',[0 2000],'fontsize',14,'XColor','k','YColor','k',...
        'linewidth',1.5,'box','off');
    ylims = get(gca,'ylim');
    set(gca,'ylim',[0 1500]);
    title(s_names{dId},'fontsize',14);
    if dId>2, xlabel('Epoch duration (seconds)','fontsize',16); end
    if mod(dId,2)==1, ylabel('Number of epochs','fontsize',16); end
    
end


%% mean percent time in each state
goodSTATES = rmfield(STATES,'DAYSTART');
% goodSTATES = STATES;
allnames = fieldnames(goodSTATES);

for ss = 1:numel(allnames)
    REMtemp = 0; NREMtemp = 0;
    AWtemp = 0; QWtemp = 0;
   thisAnimal = allnames{ss};
   fprintf('Working on animal %s (%u out of %u).\n',thisAnimal,ss,numel(allnames));
   theseStates = STATES.(thisAnimal);
   for jj = 1:size(theseStates,1)-1
       statecode = theseStates(jj,1);
       statedur = theseStates(jj+1,2) - theseStates(jj,2);
       switch statecode
           case 1
               REMtemp = REMtemp + statedur; 
           case 2
               NREMtemp = NREMtemp + statedur;
           case 4
               AWtemp = AWtemp + statedur;
           case 5
               QWtemp = QWtemp + statedur;
       end
   end
   REMtime(ss) = REMtemp;
   NREMtime(ss) = NREMtemp;
   AWtime(ss) = AWtemp;
   QWtime(ss) = QWtemp;
   totTime(ss) = theseStates(end,2) - theseStates(1,2);
   clear *temp
end

stateTime = [AWtime; QWtime; NREMtime; REMtime];
percTime = bsxfun(@rdivide,stateTime,totTime).*100;

meanPercTime = nanmean(percTime,2);
semPercTime = nanstd(percTime,[],2)./sqrt(numel(allnames));

h_time = figure(); hold on;
set(h_time,'color','w','unit','normalized','position',[.1 .1 .5 .7]);
for tt = 1:4
    tbar{tt} = bar(tt,meanPercTime(tt),.5);
    tbar{tt}.FaceColor = s4colors_2{tt};
    tbar{tt}.EdgeColor = 'k';
    tbar{tt}.FaceAlpha = 1;
    tbar{tt}.LineWidth = 2;
    terr{tt} = errorbar(tt,meanPercTime(tt),semPercTime(tt),'linestyle','none',...
        'linewidth',2,'capsize',0,'color',[0 0 0]);
end
set(gca,'linewidth',2,'fontsize',15,'XColor','k','YColor','k','xtick',1:4,...
    'xticklabel',s_names,'ytick',[0:10:100]);
ylabel('Mean percentage time in state','fontsize',17);


%% rank cells by AW 
all_FR = [AW_meanFRbycell, QW_meanFRbycell, NREM_meanFRbycell, REM_meanFRbycell];
all_FR_sorted = sortrows(all_FR,1);
all_FR_sorted = all_FR_sorted(~isnan(all_FR_sorted(:,1)),:);
n_cells = size(all_FR_sorted,1);
x_pts = 1:n_cells;

FRrank_fig = figure();
set(FRrank_fig,'color','w','unit','normalized','position',[.1 .1 .7 .8]);
nRows = 1;
nCols = 4;
colW = .7/nCols; rowH = .7/nRows;
colX = .07 + linspace(0,0.96,nCols+1); colX = colX(1:end-1);
rowY = .15 + linspace(.96,0,nRows+1); rowY = rowY(2:end);

xlabs = {'Active','Quiet','NREM','REM'};

for xx = 1:4
    % make axes
    rowId = ceil( xx / nCols ) ;
    colId = xx - (rowId - 1) * nCols ;
    rankax{xx} = axes(FRrank_fig,'position',[colX(colId) rowY(rowId) colW rowH]); hold on;
    hold on;
    set(rankax{xx},'color',[.92 .92 .92])
    
    plot(all_FR_sorted(:,xx),x_pts,'o','color',s4colors_2{xx},...
        'markerfacecolor',s4colors_2{xx},...
        'markeredgecolor','none','markersize',9);
    set(gca,'xscale','log','fontsize',14,'ytick',[],'XColor','k',...
        'YColor','k','ylim',[0 n_cells+1]);
    text(.45,25.5,xlabs{xx},'fontsize',14);
    for uu=0:n_cells
        line([0.011 100],[uu+.5 uu+.5],'color','w','linewidth',2);
    end
%     if xx>1, set(rankax{xx},'YColor','none'); end
    if xx==1
        ylabel('Cells','fontsize',16);
        text(2,105,'n = 24 cells','fontsize',14);
    elseif xx == 2
        text(50,-2.5,'Firing Rate (Hz)','fontsize',15);
    end
end
    


%% mean FR by light/dark






