%% This code includes figure S6 code to check parameters of S-W dense analysis

clearvars -except CONTCELL* tempSHANK CONT24_SHANK newCONTSHANK
clc
s_color = [84 119 146]./255;
w_color = [206,107,77]./255;
swcols = {w_color,s_color};
c_rem   = [25 181 149]./255;
c_nrem  = [131 49 146]./255;
c_aw    = [201 28 101]./255;
c_qw    = [247 148 41]./255;
c_fit = [255, 205, 55]./255;

% FOR FIGURE S6 USE:
% #15 for panel A
% #16 for panel B
% #17 for panel c


%% load data
ctrl_filenum    = 17;
dep_filenum     = 17;
filenum_thresh  = 20;

if ismac
%     load_dir = '/Volumes/GoogleDrive/My Drive/BrandeisNeuroscience/TurrigianoLabWork/Eye_Reopening/SW_Figures/sw_dense/Dec2017_DATA';
    load_dir = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Eye_Reopening/SWdense_data/thresh_fig_data';
elseif ispc
    load_dir = 'Z:\ATP_MAIN\DATA\Eye_Reopening\SWdense_data\thresh_fig_data';
end

load_ctrl_name = ['SWdense_CTRL_' num2str(ctrl_filenum) '.mat'];
ctrl_load = load([load_dir filesep load_ctrl_name]);
ctrl_dat = ctrl_load.swdense_data;

load_dep_name = ['SWdense_DEP_' num2str(dep_filenum) '.mat'];
dep_load = load([load_dir filesep load_dep_name]);
dep_dat = dep_load.swdense_data;

%% compile bar data
for tt = 1:size(dep_dat,2)
    
    if tt > size(ctrl_dat,2)
        ctrl_dat(tt).mean_FR_change_S = NaN;
        ctrl_dat(tt).mean_FR_change_W = NaN;
        ctrl_dat(tt).sem_FR_change_S = NaN;
        ctrl_dat(tt).sem_FR_change_W = NaN;
    end
    
bar_data(tt,:) = [ctrl_dat(tt).mean_FR_change_S, ctrl_dat(tt).mean_FR_change_W,...
    dep_dat(tt).mean_FR_change_S, dep_dat(tt).mean_FR_change_W];

sem_data(tt,:) = [ctrl_dat(tt).sem_FR_change_S, ctrl_dat(tt).sem_FR_change_W,...
    dep_dat(tt).sem_FR_change_S, dep_dat(tt).sem_FR_change_W];

end

%{
%% bar graph
barplot = figure(); hold on;
set(barplot,'color','w','unit','normalized','position',[.1 .1 .6 .7]);

bwidth = .25;
x_pos = [.8 1.2 1.8 2.2]; 

for uu=1:4
    b{uu} = bar(x_pos(uu),bar_data(1,uu), bwidth);
    b{uu}.LineWidth = 4;
    
    if uu <=2
        b{uu}.FaceColor = 'none';
        b{uu}.EdgeColor = swcols{mod(uu,2)+1};
    else
        b{uu}.EdgeColor = 'none';
        b{uu}.FaceColor = swcols{mod(uu,2)+1};
    end
    
    e{uu} = errorbar(x_pos(uu),bar_data(1,uu),sem_data(1,uu),'color','k',...
        'linestyle','none','capsize',0,'linewidth',3);
end
if dep_filenum > filenum_thresh
    bar_ylim = 1.8;
elseif dep_filenum <= filenum_thresh
    bar_ylim = 1.5;
end
no_change_line = line([.5 2.5],[1 1],'color','k','linestyle','--','linewidth',2);
set(gca,'XColor','k','YColor','k','fontsize',18,'ylim',[0 bar_ylim],...
    'xlim',[.5 2.5],'xtick',[1 2],'xticklabel',{'Control','Deprived'},...
    'linewidth',2,'ytick',0:.5:1.5);
ylabel('Fractional change in FR','fontsize',20);

%% stats
dist_data = {ctrl_dat(1).FR_change_by_cell_S, ctrl_dat(1).FR_change_by_cell_W, ...
    dep_dat(1).FR_change_by_cell_S, dep_dat(1).FR_change_by_cell_W};

anovadata = [];
anovagroups = [];
for xx=1:4
    groupdata = dist_data{xx};
    grouplabel = xx.*ones(size(groupdata));
    anovadata = [anovadata; groupdata];
    anovagroups = [anovagroups; grouplabel];
end

[p_anova,~,anovastats] = anova1(anovadata,anovagroups,'off');
% cfig = figure();
c = multcompare(anovastats,'display','off');

p_1v3 = c(2,6);
p_1v2 = c(1,6);
p_1v4 = c(3,6);
p_2v3 = c(4,6);
p_4v3 = c(6,6);
p_2v4 = c(5,6);
if dep_filenum <= filenum_thresh
    raw_pvals = [p_1v3,p_2v3,p_4v3,p_2v4];
elseif dep_filenum > filenum_thresh
    raw_pvals = [p_1v4,p_2v3,p_4v3,p_2v4];
end

for aa=1:numel(raw_pvals)
    thisval = num2str(raw_pvals(aa),'%.20f');
    zcounter = 0;
    for dd = 3:max(size(thisval))
        if str2double(thisval(dd)) == 0
            zcounter = zcounter + 1;
        else
            break;
        end
    end
%     disp(aa);
%     disp(zcounter);
    if zcounter > 4, zcounter = 4; end
    p_vals(aa) = 1*10^(-zcounter);
end
    
%% write p vals
barplot;
if dep_filenum > filenum_thresh
    p1v4 = line([x_pos(1)+.02 x_pos(4)-.02],[bar_ylim-.15 bar_ylim-.15],'color','k',...
        'linewidth',1.5);
    t1v4 = text((x_pos(1)+.02+x_pos(4)-.02)*.48,bar_ylim-.11,sprintf('p = %.4f',raw_pvals(1)),...
        'fontsize',13);
elseif dep_filenum <= filenum_thresh
    p1v3 = line([x_pos(1)+.02 x_pos(3)-.02],[bar_ylim-.15 bar_ylim-.15],'color','k',...
        'linewidth',1.5);
    t1v3 = text((x_pos(1)+.02+x_pos(3)-.02)*.48,bar_ylim-.11,sprintf('p = %.4f',raw_pvals(1)),...
        'fontsize',13);
end

p2v3 = line([x_pos(2)+.02 x_pos(3)-.02],[bar_ylim-.25 bar_ylim-.25],'color','k',...
    'linewidth',1.5);
t2v3 = text((x_pos(2)+.02+x_pos(3)-.02)*.48,bar_ylim-.21,sprintf('p = %.4f',raw_pvals(2)),...
    'fontsize',13);

p4v3 = line([x_pos(3)+.02 x_pos(4)-.02],[bar_ylim-.25 bar_ylim-.25],'color','k',...
    'linewidth',1.5);
t4v3 = text(1.9,bar_ylim-.21,sprintf('p = %.4f',raw_pvals(3)),...
    'fontsize',13);
    

p2v4 = line([x_pos(2)+.02 x_pos(4)-.02],[bar_ylim-.05 bar_ylim-.05],'color','k',...
    'linewidth',1.5);
t2v4 = text((x_pos(2)+.02+x_pos(4)-.02)*.48,bar_ylim-.01,sprintf('p = %.4f',raw_pvals(4)),...
    'fontsize',13);

%% legend
s_rec = patch([.55 .7 .7 .55],[bar_ylim-.02 bar_ylim-0.02 bar_ylim bar_ylim],'k');
s_rec.FaceColor = swcols{2};
s_rec.EdgeColor = 'none';
s_txt = text(.73,bar_ylim-.01,'Sleep-dense','fontsize',16);
w_rec = patch([.55 .7 .7 .55],[bar_ylim-.08 bar_ylim-0.08 bar_ylim-0.06 bar_ylim-.06],'k');
w_rec.FaceColor = swcols{1};
w_rec.EdgeColor = 'none';
s_txt = text(.73,bar_ylim-.07,'Wake-dense','fontsize',16);

%% CDF plot
cdfplot = figure(); hold on;
set(cdfplot,'color','w','unit','normalized','position',[.15 .15 .6 .7]);

% subplot(1,2,1); hold on;
for hh = 1:4
    [f,x] = ecdf(dist_data{hh});
    if hh<=2
        lstyle = ':';
    else
        lstyle = '-';
    end
    plot(x,f,'linewidth',5,'color',swcols{mod(hh,2)+1},'linestyle',lstyle);
    line([1 1],[0 1],'linewidth',2,'linestyle','--','color','k');
end

set(gca,'xlim',[0 2],'fontsize',18,'XColor','k','YColor','k','ytick',[0:.2:1],...
    'xtick',[0:.5:2],'linewidth',1.5);
xlabel('Fractional change in FR','fontsize',18);
ylabel('Fraction of cells','fontsize',20);

%% legend
y_s = .3;
s_c_line = line([1.88 2],[y_s y_s],'linewidth',3,'color',swcols{2},'linestyle',':');
w_c_line = line([1.88 2],[y_s-1*.05 y_s-1*.05],'linewidth',3,'color',swcols{1},'linestyle',':');
s_d_line = line([1.88 2],[y_s-2*.05 y_s-2*.05],'linewidth',3,'color',swcols{2},'linestyle','-');
w_d_line = line([1.88 2],[y_s-3*.05 y_s-3*.05],'linewidth',3,'color',swcols{1},'linestyle','-');

s_c_text = text(1.3,y_s,'Control, Sleep-dense','fontsize',14);
w_c_text = text(1.3,y_s-1*.05,'Control, Wake-dense','fontsize',14);
s_d_text = text(1.3,y_s-2*.05,'Deprived, Sleep-dense','fontsize',14);
w_d_text = text(1.3,y_s-3*.05,'Deprived, Wake-dense','fontsize',14);

%% CDF stats

[~,p_ks_ctrl] = kstest2(dist_data{1},dist_data{2});
[~,p_ks_dep] = kstest2(dist_data{3},dist_data{4});
[~,p_ks_S] = kstest2(dist_data{1},dist_data{3});
[~,p_ks_W] = kstest2(dist_data{2},dist_data{4});

all_ps = [p_ks_ctrl,p_ks_dep,p_ks_S,p_ks_W]';

[~,~,~,adj_ps] = fdr_bh(all_ps,0.05,'dep','yes');

for aa=1:numel(adj_ps)
    thisval = num2str(adj_ps(aa),'%.20f');
    zcounter = 0;
    for dd = 3:max(size(thisval))
        if str2double(thisval(dd)) == 0
            zcounter = zcounter + 1;
        else
            greater_than = str2double(['0.' repmat('0',1,dd-3) thisval(dd)]);
            less_than = str2double(['0.' repmat('0',1,dd-3) num2str(str2double(thisval(dd))+1)]);
            break;
        end
    end
%     disp(aa);
%     disp(zcounter);
    if zcounter > 4, zcounter = 4; end
    if zcounter > 2
        p_vals(aa) = 1*10^(-zcounter);
    elseif zcounter == 1
%         keyboard;
        if less_than <= .05
            p_vals(aa) = less_than;
        else
            p_vals(aa) = greater_than;
        end
    elseif zcounter < 1
        p_vals(aa) = greater_than;
    end
end


p_x = .1;
p_y = .95;

text(p_x,1,'Paired KS test BH-adjusted p-values:','fontsize',14);
text(p_x,p_y,sprintf('S-CTRL vs W-CTRL: p = %.4f',adj_ps(1)),'fontsize',14);
text(p_x,p_y-1*.05,sprintf('S-DEP vs W-DEP: p = %.4f',adj_ps(2)),'fontsize',14);
text(p_x,p_y-2*.05,sprintf('S-CTRL vs S-DEP: p = %.4f',adj_ps(3)),'fontsize',14);
text(p_x,p_y-3*.05,sprintf('W-CTRL vs W-DEP: p = %.4f',adj_ps(4)),'fontsize',14);
%}
%% checking dense threshold dependency

norm_mode = 'change'; % change, perc or norm


switch norm_mode
    case 'norm'
        yl = [0 max([2;ceil(reshape(bar_data(:,1:2)+sem_data(:,1:2),[],1))])];
        ol = [1 1];
    case 'perc'
        yl = [-100 100];
        ol = [0 0];% change these if you ever use this
    case 'change'
        yl = [-0.3 0.3];
        ol = [0 0];
end


threshfig = figure();
set(threshfig,'color','w','unit','normalized','position',[.1 .08 .30 .6],...
    'numbertitle','off');

%{
ax1 = axes(threshfig,'position',[.1 .58 .8 .35]);
hold on;
one_line1 = line([.5 3.5],ol,'color','k','linestyle',':','linewidth',2);
e_c_s = errorbar(1:3,bar_data(:,1),sem_data(:,1),'--s','color',s_color,...
    'linewidth',2,'markerfacecolor',s_color);
e_c_w = errorbar(1:3,bar_data(:,2),sem_data(:,2),'--s','color',w_color,...
    'linewidth',2,'markerfacecolor',w_color);
set(gca,'xlim',[.5 3.5],'ylim',yl,'fontsize',16,'linewidth',1.5,...
    'xtick',[1:3],'xticklabel',{'0.70','0.75','0.80'});
title(sprintf('Control hemisphere, window size = %.1f hours',dep_dat(1).t_win_hours),...
    'fontsize',20,'interpreter','none');
xlabel('Density threshold','fontsize',18);
ylabel('Fractional change in FR','fontsize',19);
%}
lw = 1.5;
% ax2 = axes(threshfig,'position',[.21 .16 .7 .45]);
hold on;
one_line2 = line([.5 3.5],ol,'color','k','linestyle',':','linewidth',lw);
e_d_s = errorbar(1:3,bar_data(:,3),sem_data(:,3),'-s','color',s_color,...
    'linewidth',lw,'markerfacecolor',s_color);
e_d_w = errorbar(1:3,bar_data(:,4),sem_data(:,4),'-s','color',w_color,...
    'linewidth',lw,'markerfacecolor',w_color);
set(gca,'xlim',[0.5 3.5],'ylim',yl,'ytick',-0.3:0.1:0.3,'fontsize',16,'linewidth',lw-0.5,...
    'xtick',[1:3],'xticklabel',{'0.70','0.75','0.80'});
title(sprintf('Re-opened hemisphere\n window size = %.1f hours',dep_dat(1).t_win_hours),...
    'fontsize',20,'interpreter','none');
xlabel('Density threshold','fontsize',18);
ylabel('Change in FR','fontsize',19);


%% Checking dependency of correlations on threshold
%% compile bar data
%{
for tt = 1:size(dep_dat,2)
    REM = dep_dat(tt).REM_cell_data;
    NREM = dep_dat(tt).NREM_cell_data;
    S = dep_dat(tt).S_cell_data;
    QW = dep_dat(tt).QW_cell_data;
    AW = dep_dat(tt).AW_cell_data;
    W = dep_dat(tt).W_cell_data;
    
    [REM_dep_corr(tt),REM_dep_pval(tt)] = corr(REM,S);
    [NREM_dep_corr(tt),NREM_dep_pval(tt)] = corr(NREM,S);
    [AW_dep_corr(tt),AW_dep_pval(tt)] = corr(AW,W);
    [QW_dep_corr(tt),QW_dep_pval(tt)] = corr(QW,W);
    
    % ctrl
    
    REM_C = ctrl_dat(tt).REM_cell_data;
    NREM_C = ctrl_dat(tt).NREM_cell_data;
    S_C = ctrl_dat(tt).S_cell_data;
    QW_C = ctrl_dat(tt).QW_cell_data;
    AW_C = ctrl_dat(tt).AW_cell_data;
    W_C = ctrl_dat(tt).W_cell_data;
    
    [REM_ctrl_corr(tt),REM_ctrl_pval(tt)] = corr(REM_C,S_C);
    [NREM_ctrl_corr(tt),NREM_ctrl_pval(tt)] = corr(NREM_C,S_C);
    [AW_ctrl_corr(tt),AW_ctrl_pval(tt)] = corr(AW_C,W_C);
    [QW_ctrl_corr(tt),QW_ctrl_pval(tt)] = corr(QW_C,W_C);

end



corrtfig = figure();
set(corrtfig,'color','w','unit','normalized','position',[.1 .08 .7 .9],...
    'numbertitle','off');


ax1 = axes(corrtfig,'position',[.1 .58 .8 .35]);
hold on;
one_line1 = line([.5 3.5],[0 0],'color','k','linestyle',':','linewidth',2);
ctrl_rem = plot(1:3,REM_ctrl_corr,'--s','color',c_rem,...
    'linewidth',3,'markerfacecolor',c_rem,'markersize',8);
ctrl_nrem = plot(1:3,NREM_ctrl_corr,'--s','color',c_nrem,...
    'linewidth',3,'markerfacecolor',c_nrem,'markersize',8);
set(gca,'xlim',[.5 3.5],'ylim',[-1 1],'fontsize',18,'linewidth',1.5,...
    'xtick',[1:3],'xticklabel',{'0.70','0.75','0.80'},'ytick',[-1:.2:1]);
title(sprintf('Control hemisphere, t_win = %.1f hours',ctrl_dat(1).t_win_hours),...
    'fontsize',22,'interpreter','none','fontweight','normal');
xlabel('Density threshold','fontsize',20);
ylabel('Pearson r','fontsize',22);
legend([ctrl_rem,ctrl_nrem],{'REM','NREM'},'fontsize',20,'box','off',...
    'Location','NorthWest');


ax2 = axes(corrtfig,'position',[.1 .09 .8 .35]);
hold on;
one_line1 = line([.5 3.5],[0 0],'color','k','linestyle',':','linewidth',2);
dep_rem = plot(1:3,REM_dep_corr,'-s','color',c_rem,...
    'linewidth',3,'markerfacecolor',c_rem,'markersize',8);
dep_nrem = plot(1:3,NREM_dep_corr,'-s','color',c_nrem,...
    'linewidth',3,'markerfacecolor',c_nrem,'markersize',8);
set(gca,'xlim',[.5 3.5],'ylim',[-1 1],'fontsize',18,'linewidth',1.5,...
    'xtick',[1:3],'xticklabel',{'0.70','0.75','0.80'},'ytick',[-1:.2:1]);
title(sprintf('Re-opened hemisphere, t_win = %.1f hours',dep_dat(1).t_win_hours),...
    'fontsize',22,'interpreter','none','fontweight','normal');
xlabel('Density threshold','fontsize',20);
ylabel('Pearson r','fontsize',22);
legend([dep_rem,dep_nrem],{'REM','NREM'},'fontsize',20,'box','off',...
    'Location','NorthWest');




%% SAME THING BUT FOR DURATION INSTEAD OF %
for tt = 1:size(dep_dat,2)
    REM = dep_dat(tt).REM_cell_dur;
    NREM = dep_dat(tt).NREM_cell_dur;
    S = dep_dat(tt).S_cell_data;
    QW = dep_dat(tt).QW_cell_dur;
    AW = dep_dat(tt).AW_cell_dur;
    W = dep_dat(tt).W_cell_data;
    
    [REM_dep_corr(tt),REM_dep_pval(tt)] = corr(REM,S);
    [NREM_dep_corr(tt),NREM_dep_pval(tt)] = corr(NREM,S);
    [AW_dep_corr(tt),AW_dep_pval(tt)] = corr(AW,W);
    [QW_dep_corr(tt),QW_dep_pval(tt)] = corr(QW,W);
    
    % ctrl
    
    REM_C = ctrl_dat(tt).REM_cell_dur;
    NREM_C = ctrl_dat(tt).NREM_cell_dur;
    S_C = ctrl_dat(tt).S_cell_data;
    QW_C = ctrl_dat(tt).QW_cell_dur;
    AW_C = ctrl_dat(tt).AW_cell_dur;
    W_C = ctrl_dat(tt).W_cell_data;
    
    [REM_ctrl_corr(tt),REM_ctrl_pval(tt)] = corr(REM_C,S_C);
    [NREM_ctrl_corr(tt),NREM_ctrl_pval(tt)] = corr(NREM_C,S_C);
    [AW_ctrl_corr(tt),AW_ctrl_pval(tt)] = corr(AW_C,W_C);
    [QW_ctrl_corr(tt),QW_ctrl_pval(tt)] = corr(QW_C,W_C);

end



corrtfig = figure();
set(corrtfig,'color','w','unit','normalized','position',[.1 .08 .7 .9],...
    'numbertitle','off');


ax1 = axes(corrtfig,'position',[.1 .58 .8 .35]);
hold on;
one_line1 = line([.5 3.5],[0 0],'color','k','linestyle',':','linewidth',2);
ctrl_rem = plot(1:3,REM_ctrl_corr,'--s','color',c_rem,...
    'linewidth',3,'markerfacecolor',c_rem,'markersize',8);
ctrl_nrem = plot(1:3,NREM_ctrl_corr,'--s','color',c_nrem,...
    'linewidth',3,'markerfacecolor',c_nrem,'markersize',8);
set(gca,'xlim',[.5 3.5],'ylim',[-1 1],'fontsize',18,'linewidth',1.5,...
    'xtick',[1:3],'xticklabel',{'0.70','0.75','0.80'},'ytick',[-1:.2:1]);
title(sprintf('Control hemisphere, t_win = %.1f hours',ctrl_dat(1).t_win_hours),...
    'fontsize',22,'interpreter','none','fontweight','normal');
xlabel('Density threshold','fontsize',20);
ylabel('Pearson r','fontsize',22);
legend([ctrl_rem,ctrl_nrem],{'REM','NREM'},'fontsize',20,'box','off',...
    'Location','NorthWest');


ax2 = axes(corrtfig,'position',[.1 .09 .8 .35]);
hold on;
one_line1 = line([.5 3.5],[0 0],'color','k','linestyle',':','linewidth',2);
dep_rem = plot(1:3,REM_dep_corr,'-s','color',c_rem,...
    'linewidth',3,'markerfacecolor',c_rem,'markersize',8);
dep_nrem = plot(1:3,NREM_dep_corr,'-s','color',c_nrem,...
    'linewidth',3,'markerfacecolor',c_nrem,'markersize',8);
set(gca,'xlim',[.5 3.5],'ylim',[-1 1],'fontsize',18,'linewidth',1.5,...
    'xtick',[1:3],'xticklabel',{'0.70','0.75','0.80'},'ytick',[-1:.2:1]);
title(sprintf('Re-opened hemisphere, t_win = %.1f hours',dep_dat(1).t_win_hours),...
    'fontsize',22,'interpreter','none','fontweight','normal');
xlabel('Density threshold','fontsize',20);
ylabel('Pearson r','fontsize',22);
legend([dep_rem,dep_nrem],{'REM','NREM'},'fontsize',20,'box','off',...
    'Location','NorthWest');
%}

