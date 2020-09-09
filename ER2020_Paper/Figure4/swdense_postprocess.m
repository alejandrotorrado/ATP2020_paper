%% swdense_postprocess
%
% Alejandro Torrado Pacheco - 2018
%
% Plot SWdense analysis data. This script can be re-used to reproduce Fig.
% 4, after running SW_Dense_MAIN.m to get the data.

% clear workspace
clearIDE

% colors for plotting
s_color = [84 119 146]./255;
w_color = [206,107,77]./255;
swcols = {w_color,s_color};
c_rem   = [25 181 149]./255;
c_nrem  = [131 49 146]./255;
c_aw    = [201 28 101]./255;
c_qw    = [247 148 41]./255;
c_fit = [255, 205, 55]./255;


if ismac
    load_dir = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig4';
elseif ispc
    load_dir = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig4';
end
% change the filenames here if needed
load_ctrl_name = 'SWdense_CTRL_Fig4.mat';
ctrl_load = load([load_dir filesep load_ctrl_name]);
ctrl_dat = ctrl_load.swdense_data;

load_dep_name = 'SWdense_DEP_Fig4.mat';
dep_load = load([load_dir filesep load_dep_name]);
dep_dat = dep_load.swdense_data;

% compile bar data - all the calculations are already done
bar_data = [ctrl_dat.mean_FR_change_S, ctrl_dat.mean_FR_change_W,...
    dep_dat.mean_FR_change_S, dep_dat.mean_FR_change_W];

sem_data = [ctrl_dat.sem_FR_change_S, ctrl_dat.sem_FR_change_W,...
    dep_dat.sem_FR_change_S, dep_dat.sem_FR_change_W];

%% Fig 4B - bar graph
barplot = figure(); hold on;
set(barplot,'color','w','unit','normalized','position',[.1 .1 .6 .7]);

% plotting params
bwidth = .25; %bar width
x_pos = [.8 1.2 1.8 2.2]; %x-axis bar position 

% make the bars and assign the colors based on condition
for uu=1:4
    b{uu} = bar(x_pos(uu),bar_data(uu), bwidth);
    b{uu}.LineWidth = 4;
    
    if uu <=2
        b{uu}.FaceColor = 'none';
        b{uu}.EdgeColor = swcols{mod(uu,2)+1};
    else
        b{uu}.EdgeColor = 'none';
        b{uu}.FaceColor = swcols{mod(uu,2)+1};
    end
    
    e{uu} = errorbar(x_pos(uu),bar_data(uu),sem_data(uu),'color','k',...
        'linestyle','none','capsize',0,'linewidth',3);
end
bar_ylim = 0.25;
no_change_line = line([.5 2.5],[0 0],'color','k','linestyle','--','linewidth',2);
set(gca,'XColor','k','YColor','k','fontsize',18,'ylim',[-bar_ylim bar_ylim],...
    'xlim',[.5 2.5],'xtick',[1 2],'xticklabel',{'Control','Deprived'},...
    'linewidth',2,'ytick',-0.2:0.1:0.2);
ylabel('Change in FR','fontsize',20);

%% STATS
% compile data
dist_data = {ctrl_dat.FR_change_by_cell_S, ctrl_dat.FR_change_by_cell_W, ...
    dep_dat.FR_change_by_cell_S, dep_dat.FR_change_by_cell_W};

% setup One-way ANOVA
anovadata = [];
anovagroups = [];
for xx=1:4
    groupdata = dist_data{xx};
    grouplabel = xx.*ones(size(groupdata));
    anovadata = [anovadata; groupdata];
    anovagroups = [anovagroups; grouplabel];
end

% do anova
[p_anova,~,anovastats] = kruskalwallis(anovadata,anovagroups,'off');
% Tukey-Kramer post-hoc
c = multcompare(anovastats,'display','off');

% get the pairwise p-values corrected for multiple comparisons
p_1v3 = c(2,6);
p_1v2 = c(1,6);
p_1v4 = c(3,6);
p_2v3 = c(4,6);
p_4v3 = c(6,6);
p_2v4 = c(5,6);

raw_pvals = [p_1v3,p_2v3,p_4v3,p_2v4];

% format p-values for writing on plot
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
    if zcounter > 4, zcounter = 4; end
    p_vals(aa) = 1*10^(-zcounter);
end
    
% write p vals on plot
barplot;


p1v3 = line([x_pos(1)+.02 x_pos(3)-.02],[bar_ylim-.10 bar_ylim-.10],'color','k',...
    'linewidth',1.5);
t1v3 = text((x_pos(1)+.02+x_pos(3)-.02)*.48,bar_ylim-.085,sprintf('p = %.4f',raw_pvals(1)),...
    'fontsize',13);

p2v3 = line([x_pos(2)+.02 x_pos(3)-.02],[bar_ylim-.15 bar_ylim-.15],'color','k',...
    'linewidth',1.5);
t2v3 = text((x_pos(2)+.02+x_pos(3)-.02)*.48,bar_ylim-.135,sprintf('p = %.4f',raw_pvals(2)),...
    'fontsize',13);

p4v3 = line([x_pos(3)+.02 x_pos(4)-.02],[bar_ylim-.15 bar_ylim-.15],'color','k',...
    'linewidth',1.5);
t4v3 = text(1.9,bar_ylim-.135,sprintf('p = %.4f',raw_pvals(3)),...
    'fontsize',13);

p2v4 = line([x_pos(2)+.02 x_pos(4)-.02],[bar_ylim-.05 bar_ylim-.05],'color','k',...
    'linewidth',1.5);
t2v4 = text((x_pos(2)+.02+x_pos(4)-.02)*.48,bar_ylim-.035,sprintf('p = %.4f',raw_pvals(4)),...
    'fontsize',13);

% make legend
s_rec = patch([.55 .7 .7 .55],[bar_ylim-.02 bar_ylim-0.02 bar_ylim bar_ylim],'k');
s_rec.FaceColor = swcols{2};
s_rec.EdgeColor = 'none';
s_txt = text(.73,bar_ylim-.01,'Sleep-dense','fontsize',16);
w_rec = patch([.55 .7 .7 .55],[bar_ylim-.08 bar_ylim-0.08 bar_ylim-0.06 bar_ylim-.06],'k');
w_rec.FaceColor = swcols{1};
w_rec.EdgeColor = 'none';
w_txt = text(.73,bar_ylim-.07,'Wake-dense','fontsize',16);

%% Fig 4C - CDF plot
cdfplot = figure(); hold on;
set(cdfplot,'color','w','unit','normalized','position',[.15 .15 .6 .7]);

% get empirical cdf from the data
for hh = 1:4
    [f,x] = ecdf(dist_data{hh});
    if hh<=2
        lstyle = ':';
    else
        lstyle = '-';
    end
    % plot the CDF
    plot(x,f,'linewidth',5,'color',swcols{mod(hh,2)+1},'linestyle',lstyle);
    line([0 0],[0 1],'linewidth',2,'linestyle','--','color','k');
end

% format axes
set(gca,'xlim',[-.6 .6],'fontsize',18,'XColor','k','YColor','k','ytick',[0:.2:1],...
    'xtick',[-.6:.2:.6],'linewidth',1.5);
xlabel('Change in FR','fontsize',18);
ylabel('Fraction of cells','fontsize',20);

% create legend
y_s = .3;
s_c_line = line([1.88 2],[y_s y_s],'linewidth',3,'color',swcols{2},'linestyle',':');
w_c_line = line([1.88 2],[y_s-1*.05 y_s-1*.05],'linewidth',3,'color',swcols{1},'linestyle',':');
s_d_line = line([1.88 2],[y_s-2*.05 y_s-2*.05],'linewidth',3,'color',swcols{2},'linestyle','-');
w_d_line = line([1.88 2],[y_s-3*.05 y_s-3*.05],'linewidth',3,'color',swcols{1},'linestyle','-');

s_c_text = text(1.3,y_s,'Control, Sleep-dense','fontsize',14);
w_c_text = text(1.3,y_s-1*.05,'Control, Wake-dense','fontsize',14);
s_d_text = text(1.3,y_s-2*.05,'Deprived, Sleep-dense','fontsize',14);
w_d_text = text(1.3,y_s-3*.05,'Deprived, Wake-dense','fontsize',14);

%% CDF STATS

% pairwise Kolmogorov-Smirnov tests
[~,p_ks_ctrl] = kstest2(dist_data{1},dist_data{2});
[~,p_ks_dep] = kstest2(dist_data{3},dist_data{4});
[~,p_ks_S] = kstest2(dist_data{1},dist_data{3});
[~,p_ks_W] = kstest2(dist_data{2},dist_data{4});

all_ps = [p_ks_ctrl,p_ks_dep,p_ks_S,p_ks_W]';

% Bonferroni correction
n_comparisons = 4;
adj_ps = all_ps .* n_comparisons;

% format p-values for writing
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

    if zcounter > 4, zcounter = 4; end
    if zcounter > 2
        p_vals(aa) = 1*10^(-zcounter);
    elseif zcounter == 1
        if less_than <= .05
            p_vals(aa) = less_than;
        else
            p_vals(aa) = greater_than;
        end
    elseif zcounter < 1
        p_vals(aa) = greater_than;
    end
end

% write p-values on plot
p_x = .1;
p_y = .95;

text(p_x,1,'Paired KS test Bonferroni-adjusted p-values:','fontsize',14);
text(p_x,p_y,sprintf('S-CTRL vs W-CTRL: p = %.4f',adj_ps(1)),'fontsize',14);
text(p_x,p_y-1*.05,sprintf('S-DEP vs W-DEP: p = %.4f',adj_ps(2)),'fontsize',14);
text(p_x,p_y-2*.05,sprintf('S-CTRL vs S-DEP: p = %.4f',adj_ps(3)),'fontsize',14);
text(p_x,p_y-3*.05,sprintf('W-CTRL vs W-DEP: p = %.4f',adj_ps(4)),'fontsize',14);
