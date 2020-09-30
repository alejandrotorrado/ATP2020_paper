% further analysis, recov
clearvars -except CONTCELL*
clc

if ismac
    loadfile = '/Users/atorrado/Desktop/MLS_DATA/recov_analysis_SD_NEW_May2020.mat';
%     loadfile = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Eye_Reopening/recov_DATA/recov_analysis_SD_NEW_May2020.mat';
elseif ispc
    loadfile = 'Z:\ATP_MAIN\DATA\Eye_Reopening\recov_DATA\recov_analysis_SD_NEW_May2020.mat';
end
% loadfile = '/Users/atorrado/Desktop/MLS_DATA/recov_analysis.mat';
rload = load(loadfile);
% rload_FS = load(loadfile_FS);
recov = rload.recov_analysis;
% recov_FS = rload.recov_analysis_FS;

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


C_DEP = [0.49 0.18 0.55];
C_CTR = [0, 0, 0] ./ 255;

% extract data
G_bin = recov.G_bin;
dep_status = {'CONTROL','DEPRIVED'};




FRbycell_RSU = recov.DEPRIVED.RSU_FRbycell;


DEP_RSUs = recov.DEPRIVED.RSU_count;

%% first unity plots setup

md4_0 = 6.5*24*3600/G_bin;
md4_1 = 7.0*24*3600/G_bin;
meanFR_MD4 = nanmean(FRbycell_RSU(:,floor(md4_0):ceil(md4_1)),2);

er2_N0 = 8.0*24*3600/G_bin;
er2_N1 = 8.5*24*3600/G_bin;
meanFR_ER2_N = nanmean(FRbycell_RSU(:,floor(er2_N0):ceil(er2_N1)),2);

er4_N0 = 10.0*24*3600/G_bin;
er4_N1 = 10.5*24*3600/G_bin;
meanFR_ER4_N = nanmean(FRbycell_RSU(:,floor(er4_N0):ceil(er4_N1)),2);


frdata = [meanFR_MD4 meanFR_ER2_N meanFR_ER4_N];



nSeries = 3;
[nanrow,nancol] = find(isnan(frdata(:,1:nSeries)));


frdata(nanrow,:) = [];


change_mode = 'perc';
for ss=2:nSeries
    switch change_mode
        case 'fold'
            
            change_data(:,ss) = frdata(:,ss)./frdata(:,1);
            change_err(ss,1) = nanstd(change_data(:,ss))/sqrt(numel(change_data(:,ss))-1);
        case 'perc'
            
            change_data(:,ss) = 100*abs((frdata(:,ss)-frdata(:,1)))./frdata(:,1);
            change_err(ss,1) = nanstd(change_data(:,ss))/sqrt(numel(change_data(:,ss))-1);
            
            noabs_change_data(:,ss) = 100*((frdata(:,ss)-frdata(:,1)))./frdata(:,1);
            noabs_change_err(ss,1) = nanstd(change_data(:,ss))/sqrt(numel(change_data(:,ss))-1);
    end
end

%% unity plots
%{
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
%     rowId = ceil( dId / nCol ) ;
%     colId = dId - (rowId - 1) * nCol ;
%     axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
    %     plot( meanFR_BL3, frdata(:,dId),'o', 'markerfacecolor', cpur,...
    %         'markeredgecolor','none','markeralpha',0.7,'markersize',8) ;
%     figure(60+dId);
    axes(all_ax{dId});
    hold on;
    
    % dep
    sc1(dId) = plot(frdata(:,1), frdata(:,dId+1),'o', 'markersize', 10, 'markerfacecolor', C_DEP,...
        'markeredgecolor','none');%,'markerfacealpha',0.6);
    % ctrl
    
    
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
% lsc = legend([sc2(1) sc1(1)],{'Control','Deprived'},'fontsize',16,...
%     'location','northwest');
% lsc.Box = 'off';


%% BAR
% anova
anovadata = [change_data(:,2); change_data(:,3)];
depones = ones(size(change_data(:,2)));
anovagroups = [0.*depones; depones];
[p_anova1,~,anovastats] = anova1(anovadata,anovagroups,'off');
anovacompare = multcompare(anovastats,'ctype','tukey-kramer','display','off');



figure(916); hold on;
set(gcf,'color','w','units','normalized','position',[.15 .1 .3 .6])

barwidth = 0.3;
x_offset = .18;
xvs = [1, 2];
dep_x = xvs + x_offset;
ctrl_x = xvs - x_offset;

b0 = bar(dep_x,nanmean(change_data(:,2:end)),barwidth,'linewidth',2,...
    'facecolor',C_DEP,'edgecolor','none');
% b0.FaceAlpha = 0.5;
% b0.EdgeAlpha = 1;
e0 = errorbar(dep_x,nanmean(change_data(:,2:end)),change_err(2:end),'linestyle','none',...
    'linewidth',2,'capsize',0,'color',[0 0 0 1]);

switch change_mode
    case 'fold'
        ymax = 2.5;
        yaxlim = [0 ymax];
        yaxtick = yaxlim(1):.5:ymax;
        ystr = 'Fold change';
        p_move = .03;
        liney{1} = [2.3 2.3];
        liney{2} = [2.3 2.3];
        liney{3} = [1.7 1.7];
        textxy{1} = [.9 2.4];
        textxy{2} = [1.6 2.4];
        textxy{3} = [1.85 1.8];
        line([.5 2.5],[1 1],'linewidth',3,'color',ccya,'linestyle','--');
    case 'perc'
        ymax = 200;
        yaxlim = [0 ymax];
        yaxtick = yaxlim(1):50:ymax;
        ystr = '% change';
        p_move = .03;
        ly1 = 150;
        ly2 = 85;
        liney{1} = [ly1 ly1];
        liney{2} = [ly1 ly1];
        liney{3} = [ly2 ly2];
        textxy{1} = [.85 ly1+10];
        textxy{2} = [1.55 ly1+10];
        textxy{3} = [1.85 ly2+10];
end

%         bl_line = refline(0,0);
%         bl_line.Color = [0 0 0];

%         bl_line.LineWidth = 1.5;
%         bl_line.LineStyle = '-';

show_p = 1;
if show_p 
    
    line([ctrl_x(1)+p_move dep_x(1)-p_move],liney{1},'color',[cblk .8],'linewidth',1.5);
    text(textxy{1}(1),textxy{1}(2),sprintf('p = %.5f',anovacompare(1,6)),'fontsize',14);
%     
%     line([dep_x(1)+p_move dep_x(2)-p_move],liney{2},'color',[cblk .8],'linewidth',1.5);
%     text(textxy{2}(1),textxy{2}(2),sprintf('p = %.5f',anovacompare(1,6)),'fontsize',14);
%     
%     line([ctrl_x(2)+p_move dep_x(2)-p_move],liney{3},'color',[cblk .8],'linewidth',1.5);
%     text(textxy{3}(1),textxy{3}(2),sprintf('p = %.5f',anovacompare(5,6)),'fontsize',14);
end

set(gca,'xlim',[xvs(1) - 0.5, xvs(2) + 0.4],'ylim',yaxlim);
set(gca,'units','normalized','position',[.19 .1 .8 .8]);
xtickmarks = xvs;
xticklabs = {'Early ER','Late ER'};
set(gca,'fontsize',20,'Xcolor','k','Ycolor','k',...
    'xtick',xtickmarks,'xticklabel',xticklabs,'ytick',yaxtick,...
    'linewidth',2);
ylabel([ystr ' from baseline'],'color','k','fontsize',24);

legend([b0],{'Re-open'},'position',[.32 .8 0.1 .1],'fontsize',20);
legend boxoff
%}


%% dot plots

dotdata = sortrows(frdata,1);
% dotstats = [frdata];
% [pf,~,fstats] = friedman(dotstats);


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

 
n_comparisons = max(size(p_rank));
p_rank = p_rank.* n_comparisons;



p_rank(p_rank>1) = 1;


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


%% kernel density plots
%{
normit = 0;

kdefig = figure();
set(kdefig,'color','w','units','normalized','position',[.1 .1 .75 .8]);

myxt = [-400:200:800]';
myxl = [-100 600];
myxt_lab = [num2str(myxt) repmat('%',size(myxt))];

alpha = 0.7;
cg_alpha = cg.*alpha + (1-alpha);
cpur_alpha = cpur.*alpha + (1-alpha);

[f2,xi2] = ksdensity(change_data(:,2),'Function','pdf');
[f3,xi3] = ksdensity(change_data(:,3),'Function','pdf');
if normit
    f2=f2./max(f2);
    f3=f3./max(f3);
end
[~,pks23] = kstest2(change_data(:,2),change_data(:,3));
daxId(1) = axes('units','normalized','position',[.60 .12 .35 .75]); hold on;

a3 = area(xi3,f3);
a3.FaceColor = 'None';
a3.EdgeColor = C_DEP;
a3.LineWidth = 2;
a2 = area(xi2,f2);
a2.FaceColor = C_DEP;
a2.EdgeColor = 'None';
set(gca,'xtick',myxt','xticklabel',myxt,'xlim',myxl,...
    'fontsize',22,'XColor','k','YColor','k','linewidth',1.5,...
    'ytick',[0:0.002:0.014],'ylim',[0 0.014]);
xlabel('% change in FR from baseline','fontsize',24,'color','k');
% ylabel('Probability density','fontsize',26','color','k');
legend([a2,a3],{'Early ER','Late ER'},'fontsize',22,'box','off',...
    'location','northeast');
title('Deprived hemisphere','fontweight','normal','fontsize',24);
% 
% 
% [f4,xi4] = ksdensity(change_data_CTRL(:,2));
% [f5,xi5] = ksdensity(change_data_CTRL(:,3));
% [~,pks45] = kstest2(change_data_CTRL(:,2),change_data_CTRL(:,3));
% daxId(2) = axes('units','normalized','position',[.12 .12 .35 .75]); hold on;
% a4 = area(xi4,f4);
% a4.FaceColor = C_CTR;
% a4.EdgeColor = 'None';
% a5 = area(xi5,f5);
% a5.FaceColor = 'None';
% a5.EdgeColor = C_CTR;
% a5.LineWidth = 2;
% set(gca,'xtick',myxt','xticklabel',myxt,'xlim',myxl,...
%     'fontsize',22,'XColor','k','YColor','k','linewidth',2,...
%     'ytick',[0:0.002:0.014],'ylim',[0 0.014]);
% xlabel('% change in FR from baseline','fontsize',24,'color','k');
% ylabel('Probability density','fontsize',26','color','k');
% legend([a4,a5],{'Early ER','Late ER'},'fontsize',22,'box','off',...
%     'location','northeast');
% title('Control hemisphere','fontweight','normal','fontsize',24);
% 

%% CDF figures
lstyles = {'-','--',':'};
dep_c = {cmar,...
         [.8 .7 .4],...
         [.8 .9 .3]};
ctrl_c = {[0 0 0],...
          [0.45 0.45 0.45],...
          [0.75 0.75 0.75]};

kdx = [0.01:0.1:50];
for yy=1:3
    [dep_cdf(:,yy), dep_X(:,yy)] = ecdf(frdata(:,yy));
%     [dep_kcdf(:,yy)] = ksdensity(frdata(:,yy),kdx,'Function','cdf');
%     [ctrl_cdf(:,yy), ctrl_X(:,yy)] = ecdf(frdata_CTRL(:,yy));
end

cdf_fig = figure();
set(cdf_fig,'color','w','unit','normalized','position',[.1 .1 .7 .7]);
ax1 = axes(cdf_fig,'position',[.1 .1 .4 .8]);
ax2 = axes(cdf_fig,'position',[.58 .1 .4 .8]);



for kk = 1:3
    axes(ax1); hold on;
    depc(kk) = plot(dep_X(:,kk),dep_cdf(:,kk),'linewidth',4,'color',dep_c{kk},'linestyle',lstyles{kk});
%     depKc(kk) = plot(kdx,dep_kcdf(:,kk),'linewidth',4,'color',dep_c{kk},'linestyle',lstyles{kk});
%     axes(ax2); hold on;
%     ctrc(kk) = plot(ctrl_X(:,kk),ctrl_cdf(:,kk),'linewidth',4,'color',ctrl_c{kk},'linestyle',lstyles{kk});
end
set(ax1,'xscale','log','xlim',[0.05 40],'xtick',[0.1,1,10],...
    'fontsize',20,'linewidth',2,'ytick',[0:.5:1]);
legend([depc],{'Baseline','Early re-open','Late re-open'},'fontsize',20,...
    'box','off','location','northwest');
axes(ax1);
ylabel('Cumulative probability','fontsize',22);
title('Deprived');
set(ax2,'xscale','log','xlim',[0.05 40],'xtick',[0.1,1,10],...
    'fontsize',20,'linewidth',2,'ytick',[0:.5:1]);
% legend([ctrc],{'Baseline','Early re-open','Late re-open'},'fontsize',20,...
%     'box','off','location','northwest');
axes(ax2);
title('Control');


kw_dep_dat = [frdata(:), [1*ones(size(frdata,1),1); ...
    2*ones(size(frdata,1),1); 3*ones(size(frdata,1),1)]];

kruskalwallis(kw_dep_dat(:,1),kw_dep_dat(:,2))
%}

%% fr perc change

ER2_change = noabs_change_data(:,2);
ER4_change = noabs_change_data(:,3);

setFigureDefaults;
pfig = figure();
set(pfig,'position',[.1 .1 .4 .7]);

plot_dat = padcat(ER2_change, ER4_change);
labels = {'ER2','ER4'};
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
    'MarkerEdgeColor',cset1,'MarkerFaceColor',cset2,'PointSize',150,...
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

[~,t1p_er2] = ttest(ER2_change,0)
[~,t1p_er4] = ttest(ER4_change,0)

%% instead of multiple one-sample t-tests, use linear model with group as dummy variable

Data = [ER2_change; ER4_change];
Group1 = [ones(size(ER2_change)); zeros(size(ER4_change))];
Group2 = [zeros(size(ER2_change)); ones(size(ER4_change))];

tbl = table(Group1,Group2,Data,'VariableNames',{'Group1','Group2','Data'});

% fit a model with no intercept since mean we're comparing to is 0
model = fitlm(tbl,'Data ~ Group1 + Group2 - 1');
anova(model)








