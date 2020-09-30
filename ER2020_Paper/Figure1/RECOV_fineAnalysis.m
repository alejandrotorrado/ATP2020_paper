%% RECOV_fineAnalysis
%
% Alejandro Torrado Pacheco - 2017
%
% Use this script to re-create Figure 1E through 1H. 

clearIDE

% Load the data - path depends on platform
if ismac
    loadfile = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig1/recov_analysis.mat';
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
cg = [.55 .55 .55];
ccya = [0.30 0.75 0.93];
cmar = [0.64 0.08 0.18];
c_ctrl = [5, 83, 48]./255;
cgray = [.45 .45 .45];
clight = [160, 156, 0]./255;
cdark = [13, 7, 66]./255;
cblk = [0 0 0];
colorpool = {cblk,ccya,cmar,cyel,cpur,cgre,cblu,cred};
dc1 = [8, 219, 139]./255;
dc2 = [16, 164, 213]./255;
dc3 = [255, 93, 10]./255;
dc4 = [255, 217, 10]./255;
deucolorpool = {cblk,dc1,dc2,cblk,dc3,dc4};
C_DEP = [112, 166, 217] ./ 255;
C_CTR = [0, 0, 0] ./ 255;

% extract data
G_bin = recov.G_bin;
dep_status = {'CONTROL','DEPRIVED'};

% FRs of ER and CTRL RSUs
FRbycell_RSU = recov.DEPRIVED.RSU_FRbycell;
FRbycell_RSU_CTRL = recov.CONTROL.RSU_FRbycell;
% number of RSUs in each condition
DEP_RSUs = recov.DEPRIVED.RSU_count;
CTRL_RSUs = recov.CONTROL.RSU_count;


%% Setup for mean FR calculations

% define baseline and get mean FRs in that 12-hour bin
md4_0 = 6.5*24*3600/G_bin;
md4_1 = 7.0*24*3600/G_bin;
meanFR_MD4 = nanmean(FRbycell_RSU(:,floor(md4_0):ceil(md4_1)),2);
meanFR_MD4_CTRL = nanmean(FRbycell_RSU_CTRL(:,floor(md4_0):ceil(md4_1)),2);

% define Early ER (ER2)
er2_N0 = 8.0*24*3600/G_bin;
er2_N1 = 8.5*24*3600/G_bin;
meanFR_ER2_N = nanmean(FRbycell_RSU(:,floor(er2_N0):ceil(er2_N1)),2);
meanFR_ER2_N_CTRL = nanmean(FRbycell_RSU_CTRL(:,floor(er2_N0):ceil(er2_N1)),2);

% define Late ER (ER4)
er4_N0 = 10.0*24*3600/G_bin;
er4_N1 = 10.5*24*3600/G_bin;
meanFR_ER4_N = nanmean(FRbycell_RSU(:,floor(er4_N0):ceil(er4_N1)),2);
meanFR_ER4_N_CTRL = nanmean(FRbycell_RSU_CTRL(:,floor(er4_N0):ceil(er4_N1)),2);

% compile data
frdata = [meanFR_MD4 meanFR_ER2_N meanFR_ER4_N];
frdata_CTRL = [meanFR_MD4_CTRL meanFR_ER2_N_CTRL meanFR_ER4_N_CTRL];

% number of 12-hour chunks in analysis
nSeries = 3;


% eliminate NaN values
[nanrow,nancol] = find(isnan(frdata(:,1:nSeries)));
[nanrow_CTRL,nancol_CTRL] = find(isnan(frdata_CTRL(:,1:nSeries)));

frdata(nanrow,:) = [];
frdata_CTRL(nanrow_CTRL,:) = [];

% calculate change in FR in 12-hour chunks. This can be done in several
% different ways. Options (with B=FR in chunk; A=FR in baseline):
% - 'perc_real': percentage change from baseline: 100*(B-A)/A
% - 'perc': absolute percentage change from baseline: abs(100*(B-A)/A)
% - 'fold': fold change from baseline: (B/A)
% - 'pairwise': index change: (B-A)/(B+A)
% - 'raw': raw Hz change: (B-A)
%
% To re-create plots in paper, use 'perc_real'
change_mode = 'perc_real';
for ss=2:nSeries
    switch change_mode
        case 'fold'
            change_data_CTRL(:,ss) = frdata_CTRL(:,ss)./frdata_CTRL(:,1);
            change_err_CTRL(ss,1) = nanstd(change_data_CTRL(:,ss))/sqrt(numel(change_data_CTRL(:,ss))-1);
            change_data(:,ss) = frdata(:,ss)./frdata(:,1);
            change_err(ss,1) = nanstd(change_data(:,ss))/sqrt(numel(change_data(:,ss))-1);
        case 'perc'
            change_data_CTRL(:,ss) = 100*abs((frdata_CTRL(:,ss)-frdata_CTRL(:,1)))./frdata_CTRL(:,1);
            change_err_CTRL(ss,1) = nanstd(change_data_CTRL(:,ss))/sqrt(numel(change_data_CTRL(:,ss))-1);
            change_data(:,ss) = 100*abs((frdata(:,ss)-frdata(:,1)))./frdata(:,1);
            change_err(ss,1) = nanstd(change_data(:,ss))/sqrt(numel(change_data(:,ss))-1);
        case 'perc_real'
            change_data_CTRL(:,ss) = 100*((frdata_CTRL(:,ss)-frdata_CTRL(:,1)))./frdata_CTRL(:,1);
            change_err_CTRL(ss,1) = nanstd(change_data_CTRL(:,ss))/sqrt(numel(change_data_CTRL(:,ss))-1);
            change_data(:,ss) = 100*((frdata(:,ss)-frdata(:,1)))./frdata(:,1);
            change_err(ss,1) = nanstd(change_data(:,ss))/sqrt(numel(change_data(:,ss))-1);
        case 'pairwise'
            change_data_CTRL(:,ss) = (frdata_CTRL(:,ss)-frdata_CTRL(:,1))./(frdata_CTRL(:,ss)+frdata_CTRL(:,1));
            change_err_CTRL(ss,1) = nanstd(change_data_CTRL(:,ss))/sqrt(numel(change_data_CTRL(:,ss))-1);
            change_data(:,ss) = (frdata(:,ss)-frdata(:,1))./(frdata(:,ss)+frdata(:,1));
            change_err(ss,1) = nanstd(change_data(:,ss))/sqrt(numel(change_data(:,ss))-1);
        case 'raw'
            change_data_CTRL(:,ss) = frdata_CTRL(:,ss) - frdata_CTRL(:,1);
            change_err_CTRL(ss,1) = nanstd(change_data_CTRL(:,ss))/sqrt(numel(change_data_CTRL(:,ss))-1);
            change_data(:,ss) = frdata(:,ss) - frdata(:,1);
            change_err(ss,1) = nanstd(change_data(:,ss))/sqrt(numel(change_data(:,ss))-1);
    end
end

%% Fig 1E - unity plots

unityfig = figure();
set(gcf,'color','w','units','normalized','position',[.05 .1 .8 .7]);


ylabs = {'Early ER firing rate (Hz)','Late ER firing rate (Hz)'};
subtitles = {'Early eye re-opening','Late eye re-opening'};
% - Compute #rows/cols, dimensions, and positions of lower-left corners.
ax1 = axes(unityfig,'position',[.08 .11 .38 .8]);
ax2 = axes(unityfig,'position',[.59 .11 .38 .8]);
all_ax = {ax1,ax2};
% - Build subplots axes and plot data.
for dId = 1 : 2

    axes(all_ax{dId});
    hold on;
    
    % dep
    sc1(dId) = plot(frdata(:,1), frdata(:,dId+1),'o', 'markersize', 10, 'markerfacecolor', C_DEP,...
        'markeredgecolor','none');
    % ctrl
    sc2(dId) = plot(frdata_CTRL(:,1), frdata_CTRL(:,dId+1), 'o', 'markersize', 10, 'markerfacecolor', C_CTR,...
        'markeredgecolor','none');
    
    xlabel( 'Baseline firing rate (Hz)' ,'color','k','fontsize',24) ; 
    ylabel( ylabs{dId},'color','k' ,'fontsize',26) ;    
    set(gca,'xscale','log','yscale','log');
    axes_lims = [0.01 100];
    ticklabels = {'0.01','0.1','1','10','100'};
    set(gca,'xlim',axes_lims,'ylim',axes_lims,'xticklabel',ticklabels,...
        'yticklabel',ticklabels);
    
    rline{dId} = refline(1,0);
    rline{dId}.Color = [cyel 0.4];
    rline{dId}.LineWidth = 4;
    
    set(gca,'fontsize',22,'XColor','k','YColor','k','linewidth',2);
    title(subtitles{dId},'fontsize',18);
    
end




%% Fig 1G - Ladder plots

% sort the data
dotdata = sortrows(frdata,1);
dotdata_CTRL = sortrows(frdata_CTRL,1);
% exclude one cell with really low FR
[i,j] = find(dotdata_CTRL < 0.001);
dotdata_CTRL(i,:) = NaN;

% STATS
% Wilocoxon rank-sum tests
% Syntax: signrank(x,y,'tail',tail)
% Compares medians of x and y distributions where samples are paired.
% Default is two-tailed (tail='both'). If set tail='right', testing that
% x>y. If tail='left', testing that x<y.

% MD4 vs ER2
p_rank(1) = signrank(dotdata(:,1),dotdata(:,2),'tail','both');
p_rank_ctrl(1) = signrank(dotdata_CTRL(:,1),dotdata_CTRL(:,2),'tail','both');
% ER2 vs ER4
p_rank(2) = signrank(dotdata(:,2),dotdata(:,3),'tail','both');
p_rank_ctrl(2) = signrank(dotdata_CTRL(:,2),dotdata_CTRL(:,3),'tail','both');
% MD4 vs ER4
p_rank(3) = signrank(dotdata(:,1),dotdata(:,3),'tail','both');
p_rank_ctrl(3) = signrank(dotdata_CTRL(:,1),dotdata_CTRL(:,3),'tail','both');
 
% Bonferroni correction
n_comparisons = max(size(p_rank));
p_rank = p_rank.* n_comparisons;
n_comparisons_ctrl = max(size(p_rank_ctrl));
p_rank_ctrl = p_rank_ctrl.* n_comparisons_ctrl;

p_rank(p_rank>1) = 1;
p_rank_ctrl(p_rank_ctrl>1) = 1;

% plotting
clear jitter
dotplot = figure();
set(dotplot,'color','w','units','normalized','position',[.1 .1 .85 .8]);
% axes 1 - deprived
daxId(1) = axes('units','normalized','position',[.58 .1 .4 .8]);
for pp=1:3
    axes(daxId(1)); hold on;
    datacol = dotdata(:,pp);
    jitter(:,pp) = (rand(size(datacol))-0.5)/10;
    
    scat{pp} = plot(jitter(:,pp) + pp,datacol,'o','markersize',12,'MarkerFaceColor',...
        C_DEP,'markeredgecolor',C_DEP);
    
    % write p-vals on plot
    pval_x0 = pp + .1;
    pval_x1 = pp + .9;
    pval_y0 = 40;
    pval_y1 = 50;
    pval_y2 = 80;
    pval_y3 = 100;
    if pp < 3
        line([pval_x0 pval_x1],[pval_y0 pval_y0],'linewidth',2,'color',[cblk]);
        text(pval_x0+.1,pval_y1,sprintf('p = %1.4f',p_rank(pp)),'fontsize',14,...
            'color','k');
    else
        line([1 3],[pval_y2 pval_y2],'linewidth',2,'color',[cblk]);
        text(1.7,pval_y3,sprintf('p = %1.4f',p_rank(pp)),'fontsize',14,...
            'color','k');
    end
    
end
set(daxId(1),'xlim',[0.5 3.5],'yscale','log','fontsize',22,'xtick',[1:3],...
    'xticklabel',{'Baseline','Early ER','Late ER'},'XColor','k','YColor','k',...
    'ylim',[.01 100],'yticklabel',{'0.01','0.1','1','10','100'},'linewidth',2);
ylabel('Firing rate (Hz)','fontsize',26,'color','k');

% display connector lines
for pt = 1:size(dotdata,1)
    connector{pt} = plot([1:3]+jitter(pt,1:3),dotdata(pt,1:3),'color',[cblk 0.4],...
        'linewidth',2);
    uistack(connector{pt},'bottom');
end

clear jitter
% axes 2 - control
daxId(2) = axes('units','normalized','position',[.1 .1 .4 .8]);
for pp=1:3
    axes(daxId(2)); hold on;
    datacol = dotdata_CTRL(:,pp);
    jitter(:,pp) = (rand(size(datacol))-0.5)/10;
    
    scat{3+pp} = plot(jitter(:,pp) + pp,datacol,'o','markersize',12,'MarkerFaceColor',...
        C_CTR,'markeredgecolor',C_CTR);

    pval_x0 = pp + .1;
    pval_x1 = pp + .9;
    pval_y0 = 40;
    pval_y1 = 50;
    pval_y2 = 80;
    pval_y3 = 100;
    
end
set(daxId(2),'xlim',[0.5 3.5],'yscale','log','fontsize',22,'xtick',[1:3],...
    'xticklabel',{'Baseline','Early ER','Late ER'},'XColor','k','YColor','k',...
    'ylim',[.01 100],'yticklabel',{'0.01','0.1','1','10','100'},'linewidth',2);
ylabel('Firing rate (Hz)','fontsize',26,'color','k');

% draw connector lines
for pt = 1:size(dotdata_CTRL,1)
    connector{pt} = plot([1:3]+jitter(pt,1:3),dotdata_CTRL(pt,1:3),'color',[cblk 0.4],...
        'linewidth',2);
    uistack(connector{pt},'bottom');
end


%% Fig 1H - correlation to baseline FR plots
setFigureDefaults;

% choose Pearson or Spearman correlation
corr_type = 'Pearson';

% log-transform baseline data
FR_BL = log10(frdata(:,1));
ER2_change = change_data(:,2);
ER4_change = change_data(:,3);

FR_BL_C = log10(frdata_CTRL(:,1));
ER2_change_C = change_data_CTRL(:,2);
ER4_change_C = change_data_CTRL(:,3);

% axes tick formatting
logtick = log10([[0.01:0.01:0.09],[0.1:0.1:0.9],[1:9],[10:10:100]]);
logtick_label = {'10^{-2}','','','','','','','','','10^{-1}','','','','','','','','',...
    '10^{0}','','','','','','','','','10^{1}','','','','','','','','','10^{2}'};

% axes y-limits depending on choice of FR change calculation
switch change_mode
    case 'raw'
        ylims = [-5 20];
    case 'pairwise'
        ylims = [-1 1];
    case 'perc_real'
        ylims = [-200 500];
    case 'perc'
        ylims = [-100 400];
end

% make the figure
cfr_fig = figure();
set(cfr_fig,'position',[.1 .1 .6 .7]);


plot(FR_BL, ER2_change, 'o','color',C_DEP,'markerfacecolor',C_DEP,'linewidth',2,...
    'markersize',10);
hold on;

% compute correlation and plot the line
[r_dep, p_dep] = corr(FR_BL, ER2_change, 'Type', corr_type);
coeff_dep = polyfit(FR_BL,ER2_change,1);
X_d = [-1.3:.1:1.8];
Y_d = coeff_dep(1) .* X_d + coeff_dep(2);
plot(X_d, Y_d, '-', 'color', C_DEP, 'linewidth', 3);

% format axes and write r and p values
box off
set(gca,'xlim',[-2 2],'ylim',ylims,'xtick',logtick,'xticklabel',logtick_label)
xlabel('Baseline FR (Hz)');
ylabel('% change in FR at ER2');
text(0.5,ylims(2)*.8,sprintf('%s r = %.4f\nP-val = %.4f',corr_type,r_dep,p_dep),'fontsize',18);
title('Re-opened hemisphere','fontsize',20);
line([-2 2],[0 0],'linestyle','--','color',[.65 .65 .65],'linewidth',2);


%% Fig 1F - percent change swarm plot
setFigureDefaults;
% make figure
pfig = figure();
set(pfig,'position',[.1 .1 .6 .7]);

% compile data (padded with NaNs)
plot_dat = padcat(ER2_change_C, ER2_change, ER4_change_C, ER4_change);
labels = {'ER2 Control','ER2 Re-open','ER4 Control', 'ER4 Re-open'};
cset1 = [C_CTR; C_DEP; C_CTR; C_DEP];
cset2 = [C_CTR; C_DEP; [1 1 1]; [1 1 1]];
lw = 5;
mw = 0.3;
yl = ylims;
ystr = '% change in FR from MD4';

% for each column compute mean and SEM
for uu = 1:size(plot_dat,2)
    dat_mean(uu) = nanmean(plot_dat(:,uu));
    dat_sem(uu) = nanstd(plot_dat(:,uu)) ./ sqrt(sum(~isnan(plot_dat(:,uu)))-1);
end

% swarm plots using modified UnivarScatter function
UnivarScatter(plot_dat,'Label',labels,...
    'BoxType','Quart','StdColor',[1 1 1],'SEMColor',[1 1 1],'Width',1.0,...
    'Compression',25,'MeanColor',[1 1 1],...
    'MarkerEdgeColor',cset1,'MarkerFaceColor',cset2,'PointSize',80,...
    'LineWidth',1.5,'DataTransform','None');
hold on;
% show mean and SEM
for uu = 1:size(plot_dat,2)
    line([uu-mw uu+mw],[dat_mean(uu) dat_mean(uu)],'linewidth',lw,'color','k');
end
errorbar([1:size(plot_dat,2)],dat_mean,dat_sem,'k',...
        'linestyle','none','linewidth',lw,'capsize',30);
line([0 size(plot_dat,2)+1],[0 0],'linestyle','--','color',[.65 .65 .65],'linewidth',2)
    
% format axes
box off
set(gca,'xlim',[0 size(plot_dat,2)+1],'ylim',yl);
if size(plot_dat,2) >= 4
    xtickangle(45);
end
ylabel(ystr);

% STATS - Kruskal-Wallis Test with Tukey-Kramer post-hoc
[p_a, ~, s_a] = kruskalwallis(plot_dat,[],'off');
c_a = multcompare(s_a,'display','off');
disp(c_a);