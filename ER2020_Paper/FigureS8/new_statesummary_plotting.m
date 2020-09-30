
load_dir = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_FigS8';

cpp2_load = load([load_dir filesep 'state_summary_CPP2recov.mat']);
recov_load = load([load_dir filesep 'state_summary_recov.mat']);
sdnew_load = load([load_dir filesep 'state_summary_SD_NEW.mat']);

%% mean duration data

% EYE RE-OPENING

recov_NR_full = recov_load.savedata.MEAN.nrem;
recov_R_full = recov_load.savedata.MEAN.rem;

recov_NR_full(recov_NR_full==0) = NaN;
recov_R_full(recov_R_full==0) = NaN;

recov_NR_reshape = [recov_NR_full(:,13:14) recov_NR_full(:,17:18);...
    recov_NR_full(:,15:16) recov_NR_full(:,19:20)];
recov_R_reshape = [recov_R_full(:,13:14) recov_R_full(:,17:18);...
    recov_R_full(:,15:16) recov_R_full(:,19:20)];    

recov_NR = nanmean(recov_NR_reshape);
recov_NR_sem = real(nanstd(recov_NR_reshape) ./ sqrt(sum(~isnan(recov_NR_reshape))-1));
recov_R = nanmean(recov_R_reshape);
recov_R_sem = real(nanstd(recov_R_reshape) ./ sqrt(sum(~isnan(recov_R_reshape))-1));

% CPP2

cpp2_NR_full = cpp2_load.savedata.MEAN.nrem;
cpp2_R_full = cpp2_load.savedata.MEAN.rem;

cpp2_NR_full(cpp2_NR_full==0) = NaN;
cpp2_R_full(cpp2_R_full==0) = NaN;

cpp2_NR_reshape = [cpp2_NR_full(:,13:14) cpp2_NR_full(:,17:18);...
    cpp2_NR_full(:,15:16) cpp2_NR_full(:,19:20)];
cpp2_R_reshape = [cpp2_R_full(:,13:14) cpp2_R_full(:,17:18);...
    cpp2_R_full(:,15:16) cpp2_R_full(:,19:20)];    

cpp2_NR = nanmean(cpp2_NR_reshape);
cpp2_NR_sem = real(nanstd(cpp2_NR_reshape) ./ sqrt(sum(~isnan(cpp2_NR_reshape))-1));
cpp2_R = nanmean(cpp2_R_reshape);
cpp2_R_sem = real(nanstd(cpp2_R_reshape) ./ sqrt(sum(~isnan(cpp2_R_reshape))-1));

% SLEEP DEP NEW

sd_NR_full = sdnew_load.savedata.MEAN.nrem;
sd_R_full = sdnew_load.savedata.MEAN.rem;

sd_NR_reshape = [sd_NR_full(:,13:14) sd_NR_full(:,17:18);...
    sd_NR_full(:,15:16) sd_NR_full(:,19:20)];
sd_R_reshape = [sd_R_full(:,13:14) sd_R_full(:,17:18);...
    sd_R_full(:,15:16) sd_R_full(:,19:20)];    

sd_NR = nanmean(sd_NR_reshape);
sd_NR_sem = real(nanstd(sd_NR_reshape) ./ sqrt(sum(~isnan(sd_NR_reshape))-1));
sd_R = nanmean(sd_R_reshape);
sd_R_sem = real(nanstd(sd_R_reshape) ./ sqrt(sum(~isnan(sd_R_reshape))-1));


%% plotting NREM

cset = [0 0.44 0.74;
        0.91 0.765 0;
        0.83 .14 .14];
lw = 2;
csz = 0;
yls = [0 500];

nrem_dat_long = padcat(recov_NR', cpp2_NR', sd_NR');
nrem_sem_long = padcat(recov_NR_sem', cpp2_NR_sem', sd_NR_sem');

nrem_dat_L = nrem_dat_long(1:2:3,:);
nrem_sem_L = nrem_sem_long(1:2:3,:);

nrem_dat_D = nrem_dat_long(2:2:4,:);
nrem_sem_D = nrem_sem_long(2:2:4,:);


% NREM light figure

setFigureDefaults;
f1 = figure();
set(f1,'position',[.04 .06 .7 .6]);
hold on;
b = bar(nrem_dat_L',1,'grouped');
b(1).FaceColor = cset(1,:);
b(2).FaceColor = cset(2,:);
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';

ngroups = size(nrem_dat_L', 1);
nbars = size(nrem_dat_L', 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, nrem_dat_L(i,:), nrem_sem_L(i,:), 'k', 'linestyle', 'none',...
        'linewidth',lw,'capsize',csz);
end
hold off

xlabels = {'ER','CPP2','SD'};
set(gca,'xtick',1:3,'xlim',[0 4.5],'xticklabel',xlabels,'ylim',yls);
title('NREM duration - Light phase');
ylabel('Duration (seconds)');
legend({'Before','After'});
legend boxoff



% NREM dark figure

setFigureDefaults;
f1 = figure();
set(f1,'position',[.04 .06 .7 .6]);
hold on;
b = bar(nrem_dat_D',1,'grouped');
b(1).FaceColor = cset(1,:);
b(2).FaceColor = cset(2,:);
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';

ngroups = size(nrem_dat_D', 1);
nbars = size(nrem_dat_D', 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, nrem_dat_D(i,:), nrem_sem_D(i,:), 'k', 'linestyle', 'none',...
        'linewidth',lw,'capsize',csz);
end
hold off

xlabels = {'ER','CPP2','SD'};
set(gca,'xtick',1:3,'xlim',[0 4.5],'xticklabel',xlabels,'ylim',yls);
title('NREM duration - Dark phase');
ylabel('Duration (seconds)');
legend({'Before','After'});
legend boxoff

%% plotting NREM

cset = [0 0.44 0.74;
        0.91 0.765 0;
        0.83 .14 .14];
lw = 2;
csz = 0;
yls = [0 500];

rem_dat_long = padcat(recov_R', cpp2_R', sd_R');
rem_sem_long = padcat(recov_R_sem', cpp2_R_sem', sd_R_sem');

rem_dat_L = rem_dat_long(1:2:3,:);
rem_sem_L = rem_sem_long(1:2:3,:);

rem_dat_D = rem_dat_long(2:2:4,:);
rem_sem_D = rem_sem_long(2:2:4,:);


% NREM light figure

setFigureDefaults;
f1 = figure();
set(f1,'position',[.04 .06 .7 .6]);
hold on;
b = bar(rem_dat_L',1,'grouped');
b(1).FaceColor = cset(1,:);
b(2).FaceColor = cset(2,:);
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';

ngroups = size(rem_dat_L', 1);
nbars = size(rem_dat_L', 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, rem_dat_L(i,:), rem_sem_L(i,:), 'k', 'linestyle', 'none',...
        'linewidth',lw,'capsize',csz);
end
hold off

xlabels = {'ER','CPP2','SD'};
set(gca,'xtick',1:3,'xlim',[0 4.5],'xticklabel',xlabels,'ylim',yls);
title('REM duration - Light phase');
ylabel('Duration (seconds)');
legend({'Before','After'});
legend boxoff



% NREM dark figure

setFigureDefaults;
f1 = figure();
set(f1,'position',[.04 .06 .7 .6]);
hold on;
b = bar(rem_dat_D',1,'grouped');
b(1).FaceColor = cset(1,:);
b(2).FaceColor = cset(2,:);
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';

ngroups = size(rem_dat_D', 1);
nbars = size(rem_dat_D', 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, rem_dat_D(i,:), rem_sem_D(i,:), 'k', 'linestyle', 'none',...
        'linewidth',lw,'capsize',csz);
end
hold off

xlabels = {'ER','CPP2','SD'};
set(gca,'xtick',1:3,'xlim',[0 4.5],'xticklabel',xlabels,'ylim',yls);
title('REM duration - Dark phase');
ylabel('Duration (seconds)');
legend({'Before','After'});
legend boxoff

%% stats




%% average time in state data

% EYE RE-OPENING

recov_tNR_full = recov_load.savedata.DUR.nrem;
recov_tR_full = recov_load.savedata.DUR.rem;

recov_tNR_full(recov_tNR_full==0) = NaN;
recov_tR_full(recov_tR_full==0) = NaN;

recov_tNR_reshape = [recov_tNR_full(:,13:14) recov_tNR_full(:,17:18);...
    recov_tNR_full(:,15:16) recov_tNR_full(:,19:20)];
recov_tR_reshape = [recov_tR_full(:,13:14) recov_tR_full(:,17:18);...
    recov_tR_full(:,15:16) recov_tR_full(:,19:20)];    

recov_tNR = nanmean(recov_tNR_reshape);
recov_tNR_sem = real(nanstd(recov_tNR_reshape) ./ sqrt(sum(~isnan(recov_tNR_reshape))-1));
recov_tR = nanmean(recov_tR_reshape);
recov_tR_sem = real(nanstd(recov_tR_reshape) ./ sqrt(sum(~isnan(recov_tR_reshape))-1));

% CPP2

cpp2_tNR_full = cpp2_load.savedata.DUR.nrem;
cpp2_tR_full = cpp2_load.savedata.DUR.rem;

cpp2_tNR_full(cpp2_tNR_full==0) = NaN;
cpp2_tR_full(cpp2_tR_full==0) = NaN;

cpp2_tNR_reshape = [cpp2_tNR_full(:,13:14) cpp2_tNR_full(:,17:18);...
    cpp2_tNR_full(:,15:16) cpp2_tNR_full(:,19:20)];
cpp2_tR_reshape = [cpp2_tR_full(:,13:14) cpp2_tR_full(:,17:18);...
    cpp2_tR_full(:,15:16) cpp2_tR_full(:,19:20)];    

cpp2_tNR = nanmean(cpp2_tNR_reshape);
cpp2_tNR_sem = real(nanstd(cpp2_tNR_reshape) ./ sqrt(sum(~isnan(cpp2_tNR_reshape))-1));
cpp2_tR = nanmean(cpp2_tR_reshape);
cpp2_tR_sem = real(nanstd(cpp2_tR_reshape) ./ sqrt(sum(~isnan(cpp2_tR_reshape))-1));

% SLEEP DEP NEW

sd_tNR_full = sdnew_load.savedata.DUR.nrem;
sd_tR_full = sdnew_load.savedata.DUR.rem;

sd_tNR_full(sd_tNR_full==0) = NaN;
sd_tR_full(sd_tR_full==0) = NaN;

sd_tNR_reshape = [sd_tNR_full(:,13:14) sd_tNR_full(:,17:18);...
    sd_tNR_full(:,15:16) sd_tNR_full(:,19:20)];
sd_tR_reshape = [sd_tR_full(:,13:14) sd_tR_full(:,17:18);...
    sd_tR_full(:,15:16) sd_tR_full(:,19:20)];    

sd_tNR = nanmean(sd_tNR_reshape);
sd_tNR_sem = real(nanstd(sd_tNR_reshape) ./ sqrt(sum(~isnan(sd_tNR_reshape))-1));
sd_tR = nanmean(sd_tR_reshape);
sd_tR_sem = real(nanstd(sd_tR_reshape) ./ sqrt(sum(~isnan(sd_tR_reshape))-1));


%% plotting NREM

cset = [0 0.44 0.74;
        0.91 0.765 0;
        0.83 .14 .14];
lwt = 2;
cszt = 0;
ylst = [0 8];

nrem_tdat_long = padcat(recov_tNR', cpp2_tNR', sd_tNR');
nrem_tsem_long = padcat(recov_tNR_sem', cpp2_tNR_sem', sd_tNR_sem');

nrem_tdat_L = nrem_tdat_long(1:2:3,:);
nrem_tsem_L = nrem_tsem_long(1:2:3,:);

nrem_tdat_D = nrem_tdat_long(2:2:4,:);
nrem_tsem_D = nrem_tsem_long(2:2:4,:);



% NREM light figure

setFigureDefaults;
f1 = figure();
set(f1,'position',[.04 .06 .7 .6]);
hold on;
b = bar(nrem_tdat_L',1,'grouped');
b(1).FaceColor = cset(1,:);
b(2).FaceColor = cset(2,:);
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';

ngroups = size(nrem_tdat_L', 1);
nbars = size(nrem_tdat_L', 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, nrem_tdat_L(i,:), nrem_tsem_L(i,:), 'k', 'linestyle', 'none',...
        'linewidth',lw,'capsize',csz);
end
hold off

xlabels = {'ER','CPP2','SD'};
set(gca,'xtick',1:3,'xlim',[0 4.5],'xticklabel',xlabels,'ylim',ylst);
title('NREM time in state - Light phase');
ylabel('times in state(hours)');
legend({'Before','After'});
legend boxoff



% NREM dark figure

setFigureDefaults;
f1 = figure();
set(f1,'position',[.04 .06 .7 .6]);
hold on;
b = bar(nrem_tdat_D',1,'grouped');
b(1).FaceColor = cset(1,:);
b(2).FaceColor = cset(2,:);
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';

ngroups = size(nrem_tdat_D', 1);
nbars = size(nrem_tdat_D', 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, nrem_tdat_D(i,:), nrem_tsem_D(i,:), 'k', 'linestyle', 'none',...
        'linewidth',lw,'capsize',csz);
end
hold off

xlabels = {'ER','CPP2','SD'};
set(gca,'xtick',1:3,'xlim',[0 4.5],'xticklabel',xlabels,'ylim',ylst);
title('NREM time in state - Dark phase');
ylabel('times in state(hours)');
legend({'Before','After'});
legend boxoff

%% plotting NREM

cset = [0 0.44 0.74;
        0.91 0.765 0;
        0.83 .14 .14];
lwt = 2;
cszt = 0;
ylst = [0 8];

rem_tdat_long = padcat(recov_tR', cpp2_tR', sd_tR');
rem_tsem_long = padcat(recov_tR_sem', cpp2_tR_sem', sd_tR_sem');

rem_tdat_L = rem_tdat_long(1:2:3,:);
rem_tsem_L = rem_tsem_long(1:2:3,:);

rem_tdat_D = rem_tdat_long(2:2:4,:);
rem_tsem_D = rem_tsem_long(2:2:4,:);



% NREM light figure

setFigureDefaults;
f1 = figure();
set(f1,'position',[.04 .06 .7 .6]);
hold on;
b = bar(rem_tdat_L',1,'grouped');
b(1).FaceColor = cset(1,:);
b(2).FaceColor = cset(2,:);
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';

ngroups = size(rem_tdat_L', 1);
nbars = size(rem_tdat_L', 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, rem_tdat_L(i,:), rem_tsem_L(i,:), 'k', 'linestyle', 'none',...
        'linewidth',lw,'capsize',csz);
end
hold off

xlabels = {'ER','CPP2','SD'};
set(gca,'xtick',1:3,'xlim',[0 4.5],'xticklabel',xlabels,'ylim',ylst);
title('REM time in state - Light phase');
ylabel('times in state(hours)');
legend({'Before','After'});
legend boxoff



% NREM dark figure

setFigureDefaults;
f1 = figure();
set(f1,'position',[.04 .06 .7 .6]);
hold on;
b = bar(rem_tdat_D',1,'grouped');
b(1).FaceColor = cset(1,:);
b(2).FaceColor = cset(2,:);
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';

ngroups = size(rem_tdat_D', 1);
nbars = size(rem_tdat_D', 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, rem_tdat_D(i,:), rem_tsem_D(i,:), 'k', 'linestyle', 'none',...
        'linewidth',lw,'capsize',csz);
end
hold off

xlabels = {'ER','CPP2','SD'};
set(gca,'xtick',1:3,'xlim',[0 4.5],'xticklabel',xlabels,'ylim',ylst);
title('REM time in state - Dark phase');
ylabel('times in state(hours)');
legend({'Before','After'});
legend boxoff
