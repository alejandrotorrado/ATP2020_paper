function [fig_h] = plot_mini_data_4cond(plot_dat, parameter, labels, cset1, cset2)
% plot_mini_data_4cond
%
% This function actually makes the plots given the data. Some parameters
% are hard-coded below
%
% keyboard;
lw = 5; % hard-coded line width for means
mw = 0.3; % hard-coded length of mean lines

% Switch which measurement is being plotted
switch parameter
    case 'Ampl'
        label_str = 'mEPSC amplitude (pA)';
        yl = [6 20];
    case 'Rin'
        label_str = 'Input resistance (MOhm)';
        yl = [0 500];
    case 'Vr'
        label_str = 'Resting potential (mV)';
        yl = [-95 -55];
    case 'Cp'
        label_str = 'Capacitance (pF)';
        yl = [0 200];
    case 'Freq'
        label_str = 'mEPSC frequency (Hz)';
        yl = [0 8];
    case 'Rise'
        label_str = 'Rise time (msec)';
        yl = [1.5 2.1];
    case 'Tau'
        label_str = 'Decay time constant (msec)';
        yl = [2.0 4];
end
    

% make the figure
fig_h = figure();
set(fig_h,'position',[.06 .1 .3 .7]);

% compute mean and SEM from the data
for uu = 1:size(plot_dat,2)
    dat_mean(uu) = nanmean(plot_dat(:,uu));
    dat_sem(uu) = nanstd(plot_dat(:,uu)) ./ sqrt(sum(~isnan(plot_dat(:,uu)))-1);
end

% use a modified UnivarScatter function to make the swarm plots
UnivarScatter(plot_dat,'Label',labels,...
    'BoxType','Quart','StdColor',[1 1 1],'SEMColor',[1 1 1],'Width',1.0,...
    'Compression',25,'MeanColor',[1 1 1],...
    'MarkerEdgeColor',cset1,'MarkerFaceColor',cset2,'PointSize',100,...
    'LineWidth',1.5,'DataTransform','None');
hold on;
% plot the mean and SEM as black lines
for uu = 1:size(plot_dat,2)
    line([uu-mw uu+mw],[dat_mean(uu) dat_mean(uu)],'linewidth',lw,'color','k');
end
errorbar([1:size(plot_dat,2)],dat_mean,dat_sem,'k',...
        'linestyle','none','linewidth',lw,'capsize',30);

% format plot
box off
set(gca,'xlim',[0.5 size(plot_dat,2)+.5],'ylim',yl);
if size(plot_dat,2) >= 4
    xtickangle(45);
end
ylabel(label_str);

