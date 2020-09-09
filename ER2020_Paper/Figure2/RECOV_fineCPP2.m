%% RECOV_fineCPP2
%
% Alejandro Torrado Pacheco - 2019
%
% Use this script to re-create Figure 2G and 2H.

% clear workspace
% clearIDE
clearvars -except CONTCELL*

% Load the data - path depends on platform
if ismac
    loadfile = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig2/recov_analysis_CPP2.mat';
elseif ispc
    loadfile = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig2\recov_analysis_CPP2.mat';
end
rload = load(loadfile);
recov = rload.recov_analysis;

% define colors for plotting
cred = [0.85 0.33 0.01];
cblu = [0 0.45 0.74];
cyel = [0.93 0.69 0.13];
cpur = [0.49 0.18 0.55];
cgre = [0.47 0.67 0.19];
% cg = [0, 210, 120]./255;
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
C_DEP = [0.47 0.67 0.19];
C_CTR = [0, 0, 0] ./ 255;

% extract data
G_bin = recov.G_bin;
dep_status = {'CONTROL','DEPRIVED'};

% FRs of ER and CTRL RSUs
FRbycell_RSU = recov.DEPRIVED.RSU_FRbycell;
FRbycell_RSU_CTRL = recov.CONTROL.RSU_FRbycell;

if no_badz
    for uu=bad_dep
        FRbycell_RSU(uu,:) = nan(1,size(FRbycell_RSU,2));
    end
end


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
er2_N0 = 8.5*24*3600/G_bin;
er2_N1 = 9.0*24*3600/G_bin;
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



%% dot plots

dotdata = sortrows(frdata,1);
% dotstats = [frdata];
dotdata_CTRL = sortrows(frdata_CTRL,1);
% [pf,~,fstats] = friedman(dotstats);
[i,j] = find(dotdata_CTRL < 0.001)
dotdata_CTRL(i,:) = NaN;

% Wilocoxon rank-sum tests
% Syntax: signrank(x,y,'tail',tail)
% Compares medians of x and y distributions where samples are paired.
% Default is two-tailed (tail='both'). If set tail='right', testing that
% x>y. If tail='left', testing that x<y.

% BL3 vs MD2
p_rank(1) = signrank(dotdata(:,1),dotdata(:,2),'tail','both');
p_rank_ctrl(1) = signrank(dotdata_CTRL(:,1),dotdata_CTRL(:,2),'tail','both');
% MD2 vs MD4
p_rank(2) = signrank(dotdata(:,2),dotdata(:,3),'tail','both');
p_rank_ctrl(2) = signrank(dotdata_CTRL(:,2),dotdata_CTRL(:,3),'tail','both');
% BL3 vs MD4
p_rank(3) = signrank(dotdata(:,1),dotdata(:,3),'tail','both');
p_rank_ctrl(3) = signrank(dotdata_CTRL(:,1),dotdata_CTRL(:,3),'tail','both');
 
n_comparisons = max(size(p_rank));
p_rank = p_rank.* n_comparisons;
n_comparisons_ctrl = max(size(p_rank_ctrl));
p_rank_ctrl = p_rank_ctrl.* n_comparisons_ctrl;

p_rank(p_rank>1) = 1;
p_rank_ctrl(p_rank_ctrl>1) = 1;

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
%         'markerfacealpha',0.8,'MarkerEdgeColor',deucolorpool{pp},...
%         'MarkerEdgealpha',0.9);
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

% write p-vals
for pt = 1:size(dotdata,1)
    connector{pt} = plot([1:3]+jitter(pt,1:3),dotdata(pt,1:3),'color',[cblk 0.4],...
        'linewidth',2);
    uistack(connector{pt},'bottom');
end

clear jitter
% axes 1 - deprived
daxId(2) = axes('units','normalized','position',[.1 .1 .4 .8]);
for pp=1:3
    axes(daxId(2)); hold on;
    datacol = dotdata_CTRL(:,pp);
    jitter(:,pp) = (rand(size(datacol))-0.5)/10;
    
    scat{3+pp} = plot(jitter(:,pp) + pp,datacol,'o','markersize',12,'MarkerFaceColor',...
        C_CTR,'markeredgecolor',C_CTR);
%         'markerfacealpha',0.8,'MarkerEdgeColor',deucolorpool{pp},...
%         'MarkerEdgealpha',0.9);
    pval_x0 = pp + .1;
    pval_x1 = pp + .9;
    pval_y0 = 40;
    pval_y1 = 50;
    pval_y2 = 80;
    pval_y3 = 100;
%     if pp < 3
%         line([pval_x0 pval_x1],[pval_y0 pval_y0],'linewidth',2,'color',[cblk]);
%         text(pval_x0+.1,pval_y1,sprintf('p = %1.4f',p_rank_ctrl(pp)),'fontsize',14,...
%             'color','k');
%     else
%         line([1 3],[pval_y2 pval_y2],'linewidth',2,'color',[cblk]);
%         text(1.7,pval_y3,sprintf('p = %1.4f',p_rank_ctrl(pp)),'fontsize',14,...
%             'color','k');
%     end
    
end
set(daxId(2),'xlim',[0.5 3.5],'yscale','log','fontsize',22,'xtick',[1:3],...
    'xticklabel',{'Baseline','Early ER','Late ER'},'XColor','k','YColor','k',...
    'ylim',[.01 100],'yticklabel',{'0.01','0.1','1','10','100'},'linewidth',2);
ylabel('Firing rate (Hz)','fontsize',26,'color','k');
% write p-vals


for pt = 1:size(dotdata_CTRL,1)
    connector{pt} = plot([1:3]+jitter(pt,1:3),dotdata_CTRL(pt,1:3),'color',[cblk .4],...
        'linewidth',2);
    uistack(connector{pt},'bottom');
end



%% Fig 2H - % change in FR from MD4

ER2_change = change_data(:,2);
ER4_change = change_data(:,3);

setFigureDefaults;
pfig = figure();
set(pfig,'position',[.1 .1 .4 .7]);

plot_dat = padcat(ER2_change, ER4_change);
labels = {'CPP2-ER2','CPP2-ER4'};
cset1 = [C_DEP; C_DEP];
cset2 = [C_DEP; [1 1 1]];
lw = 5;
mw = 0.3;
ylims = [-100 400];
yl = ylims;
ystr = '% change in FR from MD4';

for uu = 1:size(plot_dat,2)
    dat_mean(uu) = nanmean(plot_dat(:,uu));
    dat_sem(uu) = nanstd(plot_dat(:,uu)) ./ sqrt(sum(~isnan(plot_dat(:,uu)))-1);
end


UnivarScatter(plot_dat,'Label',labels,...
    'BoxType','SEM','StdColor',[1 1 1],'SEMColor',[1 1 1],'Width',1.0,...
    'Compression',25,'MeanColor',[1 1 1],...
    'MarkerEdgeColor',cset1,'MarkerFaceColor',cset2,'PointSize',80,...
    'LineWidth',1.5,'DataTransform','None');
hold on;
for uu = 1:size(plot_dat,2)
    line([uu-mw uu+mw],[dat_mean(uu) dat_mean(uu)],'linewidth',lw,'color','k');
end
errorbar([1:size(plot_dat,2)],dat_mean,dat_sem,'k',...
        'linestyle','none','linewidth',lw,'capsize',30);
line([0 size(plot_dat,2)+1],[0 0],'linestyle','--','color',[.65 .65 .65],'linewidth',2)
    
    
box off
set(gca,'xlim',[0.5 size(plot_dat,2)+.5],'ylim',yl);
if size(plot_dat,2) >= 4
    xtickangle(45);
end
ylabel(ystr);

% stats
[~,p_t2] = ttest2(ER2_change, ER4_change)
[p_W] = ranksum(ER2_change, ER4_change)

[~,t1p_ER2] = ttest(ER2_change,0)
[~,t1p_ER4] = ttest(ER4_change,0)

