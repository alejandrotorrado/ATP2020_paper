%% RECOV_scalingBootstrap
%
% Alejandro Torrado Pacheco - 2019
%
% This script does a bootstrap analysis on the FR changes induced by eye
% re-opening. It can be used to reproduce Fig. S2C and S2D. 

% The principle here is simple: given that we see the % change in FR from
% baseline at ER4 has a certain distribution, can this be explained if
% instead of returning close to their baseline value neurons just assume a
% new place within the FR distribution?
%
% In other words, can we replicate the % change in FR at ER4 by NOT having
% each neuron return to its baseline FR?
%
% To do this we use a bootstrap approach, with two different methods, that
% we will call SAMPLE and SHUFFLE.
%
% ____SAMPLE METHOD
% 
% In this method, we take the distribution of FRs in ER4 and sample at
% random from this distribution to obtain an alternative distribution of
% FRs. We then compute the change in FR using the same method as we did for
% real data in Fig 1. This process is repeated a number of times (e.g.
% 10,000). We can then compute a 99% confidence interval for the average %
% change in FR at ER4 in the case of random sampling. If the real data fall
% outside this mean, it indicates that neurons do not return to randomly
% picked values within the FR distribution, but instead have an intrinsic
% set point they return to. Note that this method does NOT preserve mean
% population FR, i.e. the average pop FR in the sampled distribution is not
% constrained to be equal to the average pop FR in the experimental ER4
% data.
%
% ____SHUFFLE METHOD
%
% The shuffle method is identical to the sample one in every way, except
% that instead of randomly sampling from the ER4 distribution, the
% experimental values are shuffled to create new random distributions. This
% has the effect of preserving average population FR, but in the situation
% in which neurons do not return to their own set point. I.e. it simulates
% network-level homeostasis as opposed to individual FR set points.

% clear workspace
clearIDE

% load the data
if ismac
    loadfile = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig1/recov_analysis.mat';
elseif ispc
    loadfile = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig1\recov_analysis.mat';
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
C_DEP = [112, 166, 217] ./ 255;
C_CTR = [0, 0, 0] ./ 255;

% extract data
G_bin = recov.G_bin;
dep_status = {'CONTROL','DEPRIVED'};

% Get cell FRs
FRbycell_RSU = recov.DEPRIVED.RSU_FRbycell;
FRbycell_RSU_CTRL = recov.CONTROL.RSU_FRbycell;

DEP_RSUs = recov.DEPRIVED.RSU_count;
CTRL_RSUs = recov.CONTROL.RSU_count;



%% calculate change in FR at the different time points - baseline, early ER, late ER

md4_0 = 6.5*24*3600/G_bin;
md4_1 = 7.0*24*3600/G_bin;
meanFR_MD4 = nanmean(FRbycell_RSU(:,floor(md4_0):ceil(md4_1)),2);
meanFR_MD4_CTRL = nanmean(FRbycell_RSU_CTRL(:,floor(md4_0):ceil(md4_1)),2);

er2_N0 = 8.0*24*3600/G_bin;
er2_N1 = 8.5*24*3600/G_bin;
meanFR_ER2_N = nanmean(FRbycell_RSU(:,floor(er2_N0):ceil(er2_N1)),2);
meanFR_ER2_N_CTRL = nanmean(FRbycell_RSU_CTRL(:,floor(er2_N0):ceil(er2_N1)),2);

er4_N0 = 10.0*24*3600/G_bin;
er4_N1 = 10.5*24*3600/G_bin;
meanFR_ER4_N = nanmean(FRbycell_RSU(:,floor(er4_N0):ceil(er4_N1)),2);
meanFR_ER4_N_CTRL = nanmean(FRbycell_RSU_CTRL(:,floor(er4_N0):ceil(er4_N1)),2);

frdata = [meanFR_MD4 meanFR_ER2_N meanFR_ER4_N];
frdata_CTRL = [meanFR_MD4_CTRL meanFR_ER2_N_CTRL meanFR_ER4_N_CTRL];


nSeries = 3;
[nanrow,nancol] = find(isnan(frdata(:,1:nSeries)));
[nanrow_CTRL,nancol_CTRL] = find(isnan(frdata_CTRL(:,1:nSeries)));

frdata(nanrow,:) = [];
frdata_CTRL(nanrow_CTRL,:) = [];

change_mode = 'perc_real';
for ss=2:nSeries
    switch change_mode
        case 'raw'
            change_data_CTRL(:,ss) = frdata_CTRL(:,ss) - frdata_CTRL(:,1);
            change_err_CTRL(ss,1) = nanstd(change_data_CTRL(:,ss))/sqrt(numel(change_data_CTRL(:,ss))-1);
            change_data(:,ss) = frdata(:,ss) - frdata(:,1);
            change_err(ss,1) = nanstd(change_data(:,ss))/sqrt(numel(change_data(:,ss))-1);
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
            change_data_CTRL(:,ss) = (frdata_CTRL(:,ss)-frdata_CTRL(:,1))./(frdata_CTRL(:,ss)+frdata_CTRL(:,1))
            change_err_CTRL(ss,1) = nanstd(change_data_CTRL(:,ss))/sqrt(numel(change_data_CTRL(:,ss))-1);
            change_data(:,ss) = (frdata(:,ss)-frdata(:,1))./(frdata(:,ss)+frdata(:,1));
            change_err(ss,1) = nanstd(change_data(:,ss))/sqrt(numel(change_data(:,ss))-1);
            
    end
end


%% do the bootstrap
% CHOOSE METHOD
% 0 : SHUFFLE method
% 1 : SAMPLE method
% NOTE: Fig S2D can only be created by setting do_emp_rand = 0. However, it
% includes both SHUFFLE and SAMPLE methods.
do_emp_rand = 0;

% get the data
er4_dep = frdata(:,3);
md4_dep = frdata(:,1);
n_cells = numel(er4_dep);

% choose number of bootstrap iterations
n_bootstrap = 10000;

% do the bootstrapping
for u = 1:n_bootstrap
    if do_emp_rand
        % sample method
        % use function emprand to sample empirical distribution
        % this method is based on creating the empirical CDF, generating
        % random numbers from the uniform distribution, and then sampling
        % the CDF at those quantiles.
        er4_shuffle(:,u) = emprand(er4_dep,n_cells,1);
    else
        % shuffle method - this uses randperm
        er4_shuffle(:,u) = er4_dep(randperm(length(er4_dep)));
        % also do the sample method for Fig S2D
        er4_empshuffle(:,u) = emprand(er4_dep,n_cells,1);
    end
end

% calculate % change from MD4 for the bootstrap data
change_shuffle = 100*((er4_shuffle - md4_dep)./md4_dep);
if ~do_emp_rand
    change_emp = 100.*((er4_empshuffle - md4_dep)./md4_dep);
end

% calculate mean and std dev for bootstrap data
mean_shuffle = nanmean(change_shuffle);
std_shuffle = nanstd(change_shuffle);

%% MEAN for every bootstrap iteration
% default plotting settings
setFigureDefaults;

% make figure
bfig = figure();
set(bfig,'position',[.1 .1 .6 .7])

% plotting params
mcol = [.5 .5 .5]; % marker color
msz = 3; % marker size
lw = 2; % linewidth

% plot the mean of each bootstrap iteration as a dot
plot(mean_shuffle,'o','markerfacecolor',mcol,'color',mcol,'markersize',msz);
hold on;
% plot the overall bootstrap mean (avg across iterations) as a black line
plot(1:n_bootstrap,repmat(nanmean(mean_shuffle),1,n_bootstrap),'-k','linewidth',lw);
% plot the real mean as a red line
plot(1:n_bootstrap,repmat(nanmean(change_data(:,3)),1,n_bootstrap),'-r','linewidth',lw);
box off;
set(gca,'ylim',[-200 1200]);
xlabel('Simulation #');
ylabel('% change in FR');

%% STD DEV for every iteration
% figure
bfig1 = figure();
set(bfig1,'position',[.2 .1 .6 .7])

% params
mcol = [.5 .5 .5]; % marker color
msz = 3; % marker size
lw = 2; % linewidth

% as above but for standard deviation instead of mean
plot(std_shuffle,'o','markerfacecolor',mcol,'color',mcol,'markersize',msz);
hold on;
plot(1:n_bootstrap,repmat(nanmean(std_shuffle),1,n_bootstrap),'-k','linewidth',lw);
plot(1:n_bootstrap,repmat(nanstd(change_data(:,3)),1,n_bootstrap),'-r','linewidth',lw);
box off;
set(gca,'ylim',[-500 4000]);
xlabel('Simulation #');
ylabel('S.D. of % change in FR');


%% Fig S2C - 99% CIs for mean and STD
% compute the confidence intervals

% sort the bootstrap data
sort_mean = sort(mean_shuffle);
sort_std = sort(std_shuffle);

% choose confidence interval
alpha = 0.99;

% find the corresponding quantiles
x_d = round((n_bootstrap * (1-alpha)) / 2);
x_1 = x_d;
x_2 = n_bootstrap - x_d;

% linewidth
blw = 4;

% find the values corresponding to the CI quantiles
mean_x1 = sort_mean(x_1);
mean_x2 = sort_mean(x_2);
std_x1 = sort_std(x_1);
std_x2 = sort_std(x_2);

% make the figure
bf = figure();
set(bf,'position',[.1 .1 .3 .7]);
hold on;
% show the CI as a gray box
% mean
rectangle('position',[.8 mean_x1 .4 mean_x2-mean_x1],'facecolor',[.75 .75 .75],'edgecolor','none');
% std dev
rectangle('position',[1.8 std_x1 .4 std_x2-std_x1],'facecolor',[.75 .75 .75],'edgecolor','none')
% show the bootstrap estimate of the mean (avg across iterations)
% mean
line([.8, 1.2],[nanmean(mean_shuffle) nanmean(mean_shuffle)],'color','k','linewidth',blw);
% std dev
line([1.8, 2.2],[nanmean(std_shuffle) nanmean(std_shuffle)],'color','k','linewidth',blw);
% show the experimental data values
% mean
line([.8 1.2],[nanmean(change_data(:,3)) nanmean(change_data(:,3))],'color','r','linewidth',blw);
% std dev
line([1.8 2.2],[nanstd(change_data(:,3)) nanstd(change_data(:,3))],'color','r','linewidth',blw);
set(gca,'xlim',[.4 2.6],'ylim',[-500 3000],'xtick',[1,2],'xticklabel',{'Mean','SD'});
ylabel('% change in FR');

%% Fig S2D - compute 5 of neurons returning to within X% of their MD4 FR

% only do this if do_emp_rand == 0
if ~do_emp_rand
    
    % 10 iterations (step through from 10% to 100%)
    for u = 1:10
        
        
        % compute the % threshold
        thresh = 10*u;
        
        % find the real # of neurons coming within thresh% of their MD4 FR
        real_perc(u) = 100 * sum(abs(change_data(:,3)) <= thresh ) / numel(change_data(:,3));
        % find that same number for the shuffle method bootstrap
        shuffle_perc = 100 * sum(abs(change_shuffle) <= thresh) / size(change_shuffle,1);
        % find that same number for the sample method bootstrap
        sample_perc = 100 * sum(abs(change_emp) <= thresh) / size(change_emp,1);
        
        % average and CI for the bootstrap data
        shuffle_perc_mean(u) = nanmean(shuffle_perc);
        shuffle_perc_CI(u,:) = [prctile(shuffle_perc,0.5) prctile(shuffle_perc,99.5)];
        
        sample_perc_mean(u) = nanmean(sample_perc);
        sample_perc_CI(u,:) = [prctile(sample_perc,0.5) prctile(sample_perc,99.5)];
        
    end
    
    % now plot the figure
    bf = figure();
    set(bf,'position',[.1 .1 .8 .7]);
    hold on;
    
    % plot params
    msz = 12; % marker size
    lw = 1.5; % linewidth
    
    % extent of the CI in the plot
    e_n_sample = sample_perc_CI(:,1) - sample_perc_mean';
    e_p_sample = sample_perc_CI(:,2) - sample_perc_mean';
    
    e_n_shuffle = shuffle_perc_CI(:,1) - shuffle_perc_mean';
    e_p_shuffle = shuffle_perc_CI(:,2) - shuffle_perc_mean';
    
    % plot the data
    plot(1:10,real_perc,'or','linestyle','none','markerfacecolor','r',...
        'markersize',msz,'linewidth',lw);
    errorbar(.85:9.85,sample_perc_mean,e_n_sample,e_p_sample,...
        'ok','linestyle','none','markerfacecolor','w','markersize',msz,'linewidth',lw);
    errorbar(1.15:10.15,shuffle_perc_mean,e_n_shuffle,e_p_shuffle,...
        'ok','linestyle','none','markerfacecolor','k','markersize',msz,'linewidth',lw);
    
    set(gca,'xlim',[0.25 10.75],'ylim',[-5 100],'xtick',1:10,'xticklabel',10:10:100);
    ylabel('% neurons returning');
    xlabel('% change in FR from MD4');
    
end



