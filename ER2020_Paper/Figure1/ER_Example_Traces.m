%% ER_Example_Traces
%
% Alejandro Torrado Pacheco - 2019
%
% Plot example FR traces for Eye Re-opening experiments (Fig 1A and 1B).


% load data - change path to file here if needed
if ismac
    load_dir = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig1';
elseif ispc
    load_dir = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig1';
end
loadfile = [load_dir filesep 'recov_analysis.mat'];
rload = load(loadfile);
recov = rload.recov_analysis;

% FRs of ER and CTRL RSUs
FRbycell_RSU = recov.DEPRIVED.RSU_FRbycell;
FRbycell_RSU_CTRL = recov.CONTROL.RSU_FRbycell;
G_bin = recov.G_bin;

% Cell indexes to reproduce paper figures
% dep: 5 and 31
% ctrl: 8 and 23

er_cell1 = 5;
er_cell2 = 31;
ctrl_cell1 = 8;
ctrl_cell2 = 23;

% ER cells FR
er_fr1 = FRbycell_RSU(er_cell1,:);
er_fr2 = FRbycell_RSU(er_cell2,:);
% CTRL cells FR
ctrl_fr1 = FRbycell_RSU_CTRL(ctrl_cell1,:);
ctrl_fr2 = FRbycell_RSU_CTRL(ctrl_cell2,:);

% set baseline period and find mean baseline FR
er_bl0 = 6.5*24*3600/G_bin;
er_bl1 = 7.0*24*3600/G_bin;
er_blfr1 = nanmean(er_fr1(1,er_bl0:er_bl1));
er_blfr2 = nanmean(er_fr2(1,er_bl0:er_bl1));

ctrl_bl0 = 6.5*24*3600/G_bin;
ctrl_bl1 = 7.0*24*3600/G_bin;
ctrl_blfr1 = nanmean(ctrl_fr1(1,ctrl_bl0:ctrl_bl1));
ctrl_blfr2 = nanmean(ctrl_fr2(1,ctrl_bl0:ctrl_bl1));

% plotting colors
er_col1 = [140,117,210]./255;
er_col2 = [98,37,105]./255;
ctrl_col1 = [6,150,104]./255;
ctrl_col2 = [2,66,74]./255;

% x-axis end
xsize = 264*3600/G_bin+1;

% plotting parameters
blw = 2; % baseline indicator line width
lw = 3; % FR line width
% x-axis limits
x0 = 6.5*24*3600 / G_bin; 
x1 = xsize;
% y-axis limits
y0 = 0.1;
y1 = 100;


% remove artifact frdue to ER unplugging
% dep: both 3.1 to 4.2
% ctrl: #11 1.5 to 3.8; #30 3.0 to 4.2
% these may need to be changed for best results

erx0_1 = (7.0*24 + 3.0)*3600 / G_bin;
erx1_1 = (7.0*24 + 4.2)*3600 / G_bin;
erx0_2 = (7.0*24 + 3.0)*3600 / G_bin;
erx1_2 = (7.0*24 + 4.2)*3600 / G_bin;

ctrx0_1 = (7.0*24 + 1.5)*3600 / G_bin;
ctrx1_1 = (7.0*24 + 3.8)*3600 / G_bin;
ctrx0_2 = (7.0*24 + 3.0)*3600 / G_bin;
ctrx1_2 = (7.0*24 + 4.2)*3600 / G_bin;

er_fr1(erx0_1:erx1_1) = NaN;
er_fr2(erx0_2:erx1_2) = NaN;
ctrl_fr1(ctrx0_1:ctrx1_1) = NaN;
ctrl_fr2(ctrx0_2:ctrx1_2) = NaN;

%% figure plotting
% ER examples
setFigureDefaults;
exfig1 = figure();
set(exfig1,'position',[.1 .1 .75 .75]);
hold on;

set(gca,'yscale','log');

edgesD = (12*3600)/G_bin : (24*3600)/G_bin : xsize;

for qq = 1:size(edgesD,2)
    q = rectangle('Position',[edgesD(qq), y0, 12*3600/G_bin, y1] );
    set(q,'facecolor',[.85 .85 .85],'linestyle','none');
end

% plot the firing rates
plot(er_fr1,'-','color',er_col1,'linewidth',lw);
plot([x0 xsize],[er_blfr1 er_blfr1],'--','color',er_col1,'linewidth',blw);
plot(er_fr2,'-','color',er_col2,'linewidth',lw);
plot([0 xsize],[er_blfr2 er_blfr2],'--','color',er_col2,'linewidth',blw);

set(gca,'xlim',[x0 xsize],'ylim',[y0 y1],'xtick',[x0 + 12*3600/G_bin : (24*3600)/G_bin : xsize],...
    'xticklabel',(x0*G_bin/3600)+12:24:264);

box off;
ylabel('Firing rate (Hz)','fontsize',28);

% CTRL examples
setFigureDefaults;
exfig2 = figure();
set(exfig2,'position',[.05 .08 .75 .75]);
hold on;

set(gca,'yscale','log');

edgesD = (12*3600)/G_bin : (24*3600)/G_bin : xsize;

for qq = 1:size(edgesD,2)
    q = rectangle('Position',[edgesD(qq), y0, 12*3600/G_bin, y1] );
    set(q,'facecolor',[.85 .85 .85],'linestyle','none');
end

% plot the firing rates
plot(ctrl_fr1,'-','color',ctrl_col1,'linewidth',lw);
plot([x0 xsize],[ctrl_blfr1 ctrl_blfr1],'--','color',ctrl_col1,'linewidth',blw);
plot(ctrl_fr2,'-','color',ctrl_col2,'linewidth',lw);
plot([0 xsize],[ctrl_blfr2 ctrl_blfr2],'--','color',ctrl_col2,'linewidth',blw);

set(gca,'xlim',[x0 xsize],'ylim',[y0 y1],'xtick',[x0 + 12*3600/G_bin : (24*3600)/G_bin : xsize],...
    'xticklabel',(x0*G_bin/3600)+12:24:264,'Layer','top');

box off;
ylabel('Firing rate (Hz)','fontsize',28);


