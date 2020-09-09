%% RECOV_fineCPP1
%
% Alejandro Torrado Pacheco - 2018
%
% Use this script to re-create Figure 2E and 2F.

% clear workspace
clearIDE

% Load the data - path depends on platform
if ismac
    loadfile = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig2/recov_analysis_CPP1.mat';
elseif ispc
    loadfile = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig2\recov_analysis_CPP1.mat';
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
C_DEP = [0.64 0.08 0.18];
C_CTR = [0, 0, 0] ./ 255;

% extract data
G_bin = recov.G_bin;
dep_status = {'CONTROL','DEPRIVED'};
% RSU FR
FRbycell_RSU = recov.DEPRIVED.RSU_FRbycell;
% number of RSUs
DEP_RSUs = recov.DEPRIVED.RSU_count;

%% Setup for mean FR calculations

% define baseline and get mean FRs in that 12-hour bin
md4_0 = 6.5*24*3600/G_bin;
md4_1 = 7.0*24*3600/G_bin;
meanFR_MD4 = nanmean(FRbycell_RSU(:,floor(md4_0):ceil(md4_1)),2);

% define Early ER (ER2)
er2_N0 = 8.0*24*3600/G_bin;
er2_N1 = 8.5*24*3600/G_bin;
meanFR_ER2_N = nanmean(FRbycell_RSU(:,floor(er2_N0):ceil(er2_N1)),2);

% define Late ER (ER4)
er4_N0 = 10.0*24*3600/G_bin;
er4_N1 = 10.5*24*3600/G_bin;
meanFR_ER4_N = nanmean(FRbycell_RSU(:,floor(er4_N0):ceil(er4_N1)),2);

% compile data
frdata = [meanFR_MD4 meanFR_ER2_N meanFR_ER4_N];

% number of 12-hour chunks in analysis
nSeries = 3;


% eliminate NaN values
[nanrow,nancol] = find(isnan(frdata(:,1:nSeries)));
frdata(nanrow,:) = [];

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
            change_data(:,ss) = frdata(:,ss)./frdata(:,1);
            change_err(ss,1) = nanstd(change_data(:,ss))/sqrt(numel(change_data(:,ss))-1);
        case 'perc'
            change_data(:,ss) = 100*abs((frdata(:,ss)-frdata(:,1)))./frdata(:,1);
            change_err(ss,1) = nanstd(change_data(:,ss))/sqrt(numel(change_data(:,ss))-1);
        case 'perc_real'
            change_data(:,ss) = 100*((frdata(:,ss)-frdata(:,1)))./frdata(:,1);
            change_err(ss,1) = nanstd(change_data(:,ss))/sqrt(numel(change_data(:,ss))-1);
        case 'pairwise'
            change_data(:,ss) = (frdata(:,ss)-frdata(:,1))./(frdata(:,ss)+frdata(:,1));
            change_err(ss,1) = nanstd(change_data(:,ss))/sqrt(numel(change_data(:,ss))-1);
        case 'raw'
            change_data(:,ss) = frdata(:,ss) - frdata(:,1);
            change_err(ss,1) = nanstd(change_data(:,ss))/sqrt(numel(change_data(:,ss))-1);
    end
end


%% Fig 2E - Ladder plot

dotdata = sortrows(frdata,1);

% Wilocoxon rank-sum tests
% Syntax: signrank(x,y,'tail',tail)
% Compares medians of x and y distributions where samples are paired.
% Default is two-tailed (tail='both'). If set tail='right', testing that
% x>y. If tail='left', testing that x<y.

% MD4 vs ER2
p_rank(1) = signrank(dotdata(:,1),dotdata(:,2),'tail','both');

% ER2 vs ER4
p_rank(2) = signrank(dotdata(:,2),dotdata(:,3),'tail','both');

% MD4 vs ER4
p_rank(3) = signrank(dotdata(:,1),dotdata(:,3),'tail','both');

% Bonferroni correction
n_comparisons = max(size(p_rank));
p_rank = p_rank.* n_comparisons;
p_rank(p_rank>1) = 1;

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
% axes formatting
set(daxId(1),'xlim',[0.5 3.5],'yscale','log','fontsize',22,'xtick',[1:3],...
    'xticklabel',{'Baseline','Early ER','Late ER'},'XColor','k','YColor','k',...
    'ylim',[.01 100],'yticklabel',{'0.01','0.1','1','10','100'},'linewidth',2,...
    'Layer','top');
ylabel('Firing rate (Hz)','fontsize',26,'color','k');

% disply p-value lines
for pt = 1:size(dotdata,1)
    connector{pt} = plot([1:3]+jitter(pt,1:3),dotdata(pt,1:3),'color',[cblk 0.4],...
        'linewidth',2);
    uistack(connector{pt},'bottom');
end


%% Fig 2F - FR % change from MD4

% get the data
ER2_change = change_data(:,2);
ER4_change = change_data(:,3);

% start plotting
setFigureDefaults;
pfig = figure();
set(pfig,'position',[.1 .1 .4 .7]);

% compile data - pad with NaNs
plot_dat = padcat(ER2_change, ER4_change);
% set labels and colors
labels = {'CPP1-ER2','CPP1-ER4'};
cset1 = [C_DEP; C_DEP];
cset2 = [C_DEP; [1 1 1]];
% set plot parameters
lw = 5; % line width
mw = 0.3; % mean line length
ylims = [-100 400]; % y-axis limits
yl = ylims;
ystr = '% change in FR from MD4'; % y-axis label

% calculate the mean and SEM
for uu = 1:size(plot_dat,2)
    dat_mean(uu) = nanmean(plot_dat(:,uu));
    dat_sem(uu) = nanstd(plot_dat(:,uu)) ./ sqrt(sum(~isnan(plot_dat(:,uu)))-1);
end

% Use modified UnivarScatter function to make the swarm plots
UnivarScatter(plot_dat,'Label',labels,...
    'BoxType','SEM','StdColor',[1 1 1],'SEMColor',[1 1 1],'Width',1.0,...
    'Compression',25,'MeanColor',[1 1 1],...
    'MarkerEdgeColor',cset1,'MarkerFaceColor',cset2,'PointSize',80,...
    'LineWidth',1.5,'DataTransform','None');
hold on;
% plot the mean and SEM
for uu = 1:size(plot_dat,2)
    line([uu-mw uu+mw],[dat_mean(uu) dat_mean(uu)],'linewidth',lw,'color','k');
end
errorbar([1:size(plot_dat,2)],dat_mean,dat_sem,'k',...
        'linestyle','none','linewidth',lw,'capsize',30);
line([0 size(plot_dat,2)+1],[0 0],'linestyle','--','color',[.65 .65 .65],'linewidth',2)
    
% format plot
box off
set(gca,'xlim',[0.5 size(plot_dat,2)+.5],'ylim',yl);
if size(plot_dat,2) >= 4
    xtickangle(45);
end
ylabel(ystr);

% STATS
% two-sample t-test and Wilcoxon rank-sum test for comparing ER2 vs ER4
[~,p_t2] = ttest2(ER2_change, ER4_change)
[p_W] = ranksum(ER2_change, ER4_change)

% one-sample t-test to check if mean % change is different from 0
[~,t1p_er2] = ttest(ER2_change,0)
[~,t1p_er4] = ttest(ER4_change,0)
