%% CPP2_Example_Traces
%
% Alejandro Torrado Pacheco - 2019
%
% Plot example FR traces for Eye Re-opening experiments with CPP injections
% on ER2 and ER3 (Fig 2D).

% load data - change path to file here if needed
if ismac
load_dir = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig2';
elseif ispc
    load_dir = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig2';
end
loadfile = [load_dir filesep 'recov_analysis_CPP2.mat'];
rload = load(loadfile);
recov = rload.recov_analysis;
frbycell = recov.DEPRIVED.RSU_FRbycell;
G_bin = recov.G_bin;

% example cell index
cell1 = 2;

% for c1 = 1:size(frbycell,1)

% get FR for this cell
fr1 = frbycell(cell1,:);

% baseline FR for normalization
bl0 = 6.5*24*3600/G_bin;
bl1 = 7.0*24*3600/G_bin;
blfr1 = nanmean(fr1(1,bl0:bl1));


% plotting colors
col1 = [60,179,113]./255;

% x-axis end
xsize = 264*3600/G_bin+1;

% plotting parameters
blw = 2; % baseline indicator line width
lw = 3; % FR line width
% x-axis limits
x0 = 6.5*24*3600 / G_bin; 
x1 = xsize;
% y-axis limits
y0 = 0.01;
y1 = 10;


% get rid of ER artifact
erx0_1 = (7.0*24 + 2.1)*3600 / G_bin;
erx1_1 = (7.0*24 + 3.4)*3600 / G_bin;

fr1(erx0_1:erx1_1) = NaN;
fr2(erx0_2:erx1_2) = NaN;


%% plotting
setFigureDefaults;
exfig = figure();
set(exfig,'position',[.1 .1 .7 .5]);
hold on;

% make figure and set log-scale y-axis
set(gca,'yscale','log');
xsize = 264*3600/G_bin+1;

% make L/D bars
edgesD = (12*3600)/G_bin : (24*3600)/G_bin : xsize;

for qq = 1:size(edgesD,2)
    q = rectangle('Position',[edgesD(qq), y0, 12*3600/G_bin, y1] );
    set(q,'facecolor',[.85 .85 .85],'linestyle','none');
end

% plot the FR
plot(fr1,'-','color',col1,'linewidth',lw);
plot([x0 xsize],[blfr1 blfr1],'--','color',col1,'linewidth',blw);

% axis and ticks formatting
set(gca,'xlim',[x0 xsize],'ylim',[y0 y1],'xtick',[x0+12*3600/G_bin : (24*3600)/G_bin : xsize],...
    'xticklabel',(x0*G_bin/3600)+12:24:264,'Layer','top');

box off;
ylabel('Firing rate (Hz)','fontsize',28);



