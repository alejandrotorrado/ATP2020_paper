% further analysis, recov
clearvars -except CONTCELL* REAL_CHANGE_DATA
clc

    if ismac
        loadfile = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig1/recov_analysis.mat';
    elseif ispc
        loadfile = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig1\recov_analysis.mat';
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


C_DEP = [112, 166, 217] ./ 255;
C_CTR = [0, 0, 0] ./ 255;

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

% extract data
G_bin = recov.G_bin;
dep_status = {'CONTROL','DEPRIVED'};


bad_ctrl = [4,9,16,18,21,24];
bad_dep = [54 63 71 72 73] - 34;
no_badz = 0;


FRbycell_RSU = recov.DEPRIVED.RSU_FRbycell;
FRbycell_RSU_CTRL = recov.CONTROL.RSU_FRbycell;
anim_idx = recov.DEPRIVED.RSU_anims;
anim_idx_CTRL = recov.CONTROL.RSU_anims;


if no_badz
    for uu=bad_dep
    FRbycell_RSU(uu,:) = nan(1,size(FRbycell_RSU,2));
    end
    for uu=bad_ctrl
    FRbycell_RSU(uu,:) = nan(1,size(FRbycell_RSU_CTRL,2));
    end
end

DEP_RSUs = recov.DEPRIVED.RSU_count;
CTRL_RSUs = recov.CONTROL.RSU_count;
G_bin = recov.G_bin;


% find S.D. of neurons' FR... in 12-hour bins

bins = [6.0*24: 0.5*24 : 11*24] .* (3600/G_bin);

% binned FR by cell

for b = 1:size(bins,2)-1
    b0 = bins(b);
    b1 = bins(b+1);
    
    DEP_FR(:,b) = nanmean(FRbycell_RSU(:,b0:b1),2);
    CTRL_FR(:,b) = nanmean(FRbycell_RSU_CTRL(:,b0:b1),2);
    
    DEP_SD(:,b) = nanstd(FRbycell_RSU(:,b0:b1),[],2);
    CTRL_SD(:,b) = nanstd(FRbycell_RSU_CTRL(:,b0:b1),[],2);
end

length_thresh = 60;

day0 = 8*24*3600;
day1 = 10*24*3600;


% LOOP THRU ANIMALS


for aa = 1:n_anims
    animal = anim_fields{aa};
    
    statetimes = STATES.(animal);
    
    statetimes(statetimes(:,1) > 3) = 7; % 7 for wake
    statetimes(statetimes(:,1) < 3) = 6; % 6 for sleep
    
    % remove repeats
    st_diff = diff(statetimes(:,1));
    kill_these = find(st_diff==0);
    statetimes(kill_these+1,:) = [];
    
    % remove short states
    st_timediff = diff(statetimes(:,2));
    too_short = find(st_timediff <= length_thresh);
    statetimes(too_short,:) = [];
    
    % remove repeats again
    st_diff = diff(statetimes(:,1));
    kill_these = find(st_diff==0);
    statetimes(kill_these+1,:) = [];
    
    dayst = DAYSTART.(animal);
    
    all_anims = {CONTCELL_recov.MASTER.animal};
    first_cell = find(strcmp(animal,all_anims),1,'first');
    
    expstart_raw = unixtime(CONTCELL_recov.MASTER(first_cell).EXPTSTART);
    
%     expstart_raw = unixtime(statetimes(1,2));
    expstart_unix = [expstart_raw(1:3) 7 30 0];
    expstart_t = unixtime(expstart_unix);
    
    statetimes(:,2) = statetimes(:,2) - expstart_t + dayst*24*3600;
    
    last_day0 = find(statetimes(:,2) < day0, 1, 'last');
    first_day1 = find(statetimes(:,2) > day1, 1, 'first');
    
    statetimes(1:last_day0-1, :) = NaN;
    statetimes(first_day1+1:end, :) = NaN;
    statetimes(isnan(statetimes(:,1)),:) = [];
    
    all_sleep_starts = find(statetimes(:,1) == 6);
    
    if all_sleep_starts(end) == size(statetimes,1)
        all_sleep_starts(end) = [];
    end
    
    all_sleep_durs{aa} = statetimes(all_sleep_starts+1,2) - statetimes(all_sleep_starts,2);
    
    thisanim_FR_idx = find(strcmp(anim_idx,animal));
    
    thisanim_FR{aa} = FRbycell_RSU(thisanim_FR_idx,:);
    

end

%% 

use_eq = 'NREM';

% change mode options:
% - dur: use sleep duration
% - abs: use absolute z-score change
change_mode = 'dur';


abs_change = - 0.08;

switch use_eq
    
    case 'NREM'
        
        % NREM equation
        m = -0.000036;
        b = 0;
        
    case 'REM'
        
        % REM equation
        m = -0.000043;
        b = 0;
end

x_bl0 = 6.5*24*3600/G_bin;
x_bl1 = 7.0*24*3600/G_bin;

x0 = (8*24 - 2) * (3600/G_bin);
x1 = 8*24 * (3600/G_bin);
x2 = 10*24*(3600/G_bin);

for aa = 1:n_anims
    
    sleeps = all_sleep_durs{aa};
    FRs = thisanim_FR{aa};
    
    % baseline FR
    FR_baseline{aa} = nanmean(FRs(:,x_bl0:x_bl1),2);
    
    % FR at beginning of 48-hour period
    FR_start{aa} = nanmean(FRs(:,x0:x1),2);
    
    FR_mean = nanmean(FRs(:,x1:x2),2);
    FR_std = nanstd(FRs(:,x1:x2),[],2);
    
    new_FR = FR_start{aa};
    
    switch change_mode
        
        case 'dur'
            
            for uu = 1:numel(sleeps)
                
                s_dur = sleeps(uu);
                new_FR = new_FR + ((m*s_dur + b) .* FR_std);
                
            end
    
        case 'abs'
            
            new_FR = new_FR + ( (numel(sleeps)*abs_change) .* FR_std );
            
    end
    
    new_FR(new_FR<0) = 0;
    
    final_FR{aa} = new_FR;
end


%% plot
frdata = [cat(1,FR_baseline{:}), cat(1,FR_start{:}), cat(1,final_FR{:})];

sd_thresh = 2;

switch sd_thresh
    case 0
        nochange_list = [];
    case 1
        nochange_list = [13, 16, 18, 20, 30, 36];
    case 2
        nochange_list = [3, 13, 16, 18, 19, 20, 24, 30, 32, 34, 35, 36];
end

frdata(nochange_list,3) = frdata(nochange_list,2);

dotdata = sortrows(frdata,1);

CSET = [0, 0, 0;
        0, 0, 0;
        0.9, 0.2, 0.01];

% Wilocoxon rank-sum tests
% Syntax: signrank(x,y,'tail',tail)
% Compares medians of x and y distributions where samples are paired.
% Default is two-tailed (tail='both'). If set tail='right', testing that
% x>y. If tail='left', testing that x<y.

% BL3 vs MD2
p_rank(1) = signrank(dotdata(:,1),dotdata(:,2),'tail','both');
% MD2 vs MD4
p_rank(2) = signrank(dotdata(:,2),dotdata(:,3),'tail','both');
% BL3 vs MD4
p_rank(3) = signrank(dotdata(:,1),dotdata(:,3),'tail','both');
 
n_comparisons = max(size(p_rank));
p_rank = p_rank.* n_comparisons;

p_rank(p_rank>1) = 1

clear jitter
dotplot = figure();
set(dotplot,'color','w','units','normalized','position',[.1 .1 .4 .8]);

for pp=1:3
    hold on;
    datacol = dotdata(:,pp);
    jitter(:,pp) = (rand(size(datacol))-0.5)/10;
    
    scat{pp} = plot(jitter(:,pp) + pp,datacol,'o','markersize',12,'MarkerFaceColor',...
        CSET(pp,:),'markeredgecolor','none');
%         'markerfacealpha',0.8,'MarkerEdgeColor',deucolorpool{pp},...
%         'MarkerEdgealpha',0.9);
   
    
end
set(gca,'xlim',[0.5 3.5],'yscale','log','fontsize',22,'xtick',[1:3],...
    'xticklabel',{'Baseline','Early ER','Late ER'},'XColor','k','YColor','k',...
    'ylim',[.01 100],'yticklabel',{'0.01','0.1','1','10','100'},'linewidth',2);
ylabel('Firing rate (Hz)','fontsize',26,'color','k');


for pt = 1:size(dotdata,1)
    connector{pt} = plot([1:3]+jitter(pt,1:3),dotdata(pt,1:3),'color',[.7 .7 .7],...
        'linewidth',2);
    uistack(connector{pt},'bottom');
end

%% calc change

change_method = 'perc_real';

switch change_method
    case 'index'
        change_data = (frdata(:,3) - frdata(:,1)) ./ (frdata(:,3) + frdata(:,1));
    case 'perc_real'
        change_data = 100 * (frdata(:,3) - frdata(:,1)) ./ frdata(:,1);
        change_data(change_data==-100) = NaN;
        yl = [-100 400];
end

plot_dat = padcat(REAL_CHANGE_DATA(:,3), change_data);

lw = 5;
mw = 0.25;

cset1 = [0 0 0; .9 .2 .01];

for uu = 1:size(plot_dat,2)
    dat_mean(uu) = nanmean(plot_dat(:,uu));
    dat_sem(uu) = nanstd(plot_dat(:,uu)) ./ sqrt(sum(~isnan(plot_dat(:,uu)))-1);
end

labels = {'Real data','Simulated data'};

setFigureDefaults;
f1 = figure();
set(f1,'position',[.1 .1 .4 .7]);

UnivarScatter(plot_dat,'Label',labels,...
    'BoxType','Quart','StdColor',[1 1 1],'SEMColor',[1 1 1],'Width',1.0,...
    'Compression',25,'MeanColor',[1 1 1],...
    'MarkerEdgeColor',cset1,'MarkerFaceColor',cset1,'PointSize',140,...
    'LineWidth',2,'DataTransform','None');
hold on;
for uu = 1:size(plot_dat,2)
    line([uu-mw uu+mw],[dat_mean(uu) dat_mean(uu)],'linewidth',lw,'color','k');
end
errorbar([1:size(plot_dat,2)],dat_mean,dat_sem,'k',...
        'linestyle','none','linewidth',lw,'capsize',30);
line([0 size(plot_dat,2)+1],[0 0],'linestyle','--','color',[.65 .65 .65],'linewidth',2)
    

    
box off
set(gca,'xlim',[0.5 size(plot_dat,2)+0.5],'ylim',yl);
if size(plot_dat,2) >= 4
    xtickangle(45);
end
ylabel(labels);

% stats

p_rs = ranksum(change_data,REAL_CHANGE_DATA(:,3))


