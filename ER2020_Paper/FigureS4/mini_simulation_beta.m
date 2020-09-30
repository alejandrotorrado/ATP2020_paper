%% MINI DETECTION SIMULATION
%
% Alejandro Torrado Pacheco - July 2020
%
% 
% way 1: sigma/1.8, integer prop 0.48
% way 2: sigma/2.0, integer prop 0.50
% 

% mean single channel conductance of AMPARs
single_chan_conduct = 10; % pS
% membrane potential (absolute value)
Vm = 70; 
% mean current resulting from the opening of a single AMPAR at this Vm
single_current = single_chan_conduct * Vm*1e-3; % pA
% number of AMPARs per synapse (based on data from cerebellar synapses,
% doi: 10.1523/JNEUROSCI.2861-06.2007 ), min and max
n_AMPAR_min = 2;
n_AMPAR_max = 178;
p_open = 0.5; % mean proportion of channels opened by Glu release
p_open_sigma = 1;

% make normal distribution with values above
p_dist = makedist('Normal','mu',p_open,'sigma',p_open_sigma);
Normal_P_dist = truncate(p_dist, 0, 1);
x_p = 0:0.01:1;
y_p = pdf(Normal_P_dist,x_p);
% figure();
% plot(y_p,'r','linewidth',2);

% We assume a log-normal distribution of the number of AMPARs per synapse.
% From the paper above, we find that the mean and s.d. of this distribution
% are:
mu_X = 36.4;
sigma_X = mu_X*0.92;
% The above are the expectation and s.d. of the log-normally distributed
% variable X = e^(mu+sigma*Z), where Z is a standard normal variable and mu
% and sigma are the parameters of the lognormal distribution. We then
% compute mu and sigma based on mu_X and sigma_X:
mu = log( mu_X^2 / sqrt(mu_X^2 + sigma_X^2) );
sigma = sqrt( log( sigma_X^2/mu_X^2 + 1 ) );

% create the probability distribution object from those parameters
temp_pd = makedist('Lognormal','mu',mu,'sigma',sigma/1.8);

LogNormal_PD = truncate(temp_pd, n_AMPAR_min, n_AMPAR_max);

% compute and plot the PDF in the desired range
x_vals = n_AMPAR_min : n_AMPAR_max;
y_vals = pdf(LogNormal_PD,x_vals);

figure();
plot(x_vals,y_vals,'linewidth',2);

prop_mode = 'integer';


%% Sample simulated mini events from neurons

detection_threshold = 5.0;
% number of neurons in simulated dataset
n_cells = 35;
% sample n mini events from this neuron
n_events = 100;
% generate random set of AMPAR number per synapse
n_AMPARs = random(LogNormal_PD,n_events,n_cells);
%generate random proportion of open channels
switch prop_mode
    case 'normal'
        prop_Open = random(Normal_P_dist,n_events,n_cells);
    case 'uniform'
        unifact = 1;
        prop_Open = unifact*rand(n_events,n_cells);
    case 'integer'
        pO_range = [0.40 0.60];
        prop_Open = pO_range(1) + (pO_range(2) - pO_range(1)) .* rand(1,n_cells);
end
mini_ampl = prop_Open .* n_AMPARs .* single_current; 

do_plots = 0;
if do_plots
    for uu=1:n_cells
        figure();
        histogram(mini_ampl(:,uu),0:1:120,'facecolor','none','edgecolor','k');
        text(50, 5, sprintf('Mean amplitude: %.2f', mean(mini_ampl(:,uu))),'fontsize',12);
        hold on
        newdat = mini_ampl(:,uu);
        newdat(newdat<5) = NaN;
        histogram(newdat,0:1:120,'facecolor','none','edgecolor','r');
        text(50, 7, sprintf('Mean amplitude: %.2f', nanmean(newdat)),'fontsize',12,'color','r');
        uiwait(gcf)
    end
end

mean_ampl_bycell = nanmean(mini_ampl,1);

total_mean = nanmean(mean_ampl_bycell);
total_sem = nanstd(mean_ampl_bycell) ./ sqrt(n_cells-1);

% Now sample x events by cell at random to create simulated CDF

dataload = load('/Users/atorrado/Desktop/MLS_DATA/mini_sim/C2_CDF.mat');
cdf_data = dataload.C2_CDF;
[f,x] = ecdf(cdf_data);

n_samples_range = [30 90];

sampled_events = [];
toscale_events = [];

% rng(1);

%%
pks = 1;
counter = 0;

rng('shuffle');

while pks > 0.05
    sampled_events = [];
    toscale_events = [];
    
    for cc = 1:n_cells
%         rng();
        n_samples = randi(n_samples_range,1);
        sampled_events = [sampled_events; randsample(mini_ampl(:,cc),n_samples,'false')];
    end
    
    
%     rng(randi(1));
%     n_samples = randi(n_samples_range,1);
    counter = counter + 1;
    for cc=1:n_cells
        toscale_events = [toscale_events; randsample(mini_ampl(:,cc),n_samples,'false')];
    end
    [~,pks] = kstest2(sampled_events, toscale_events);
    if mod(counter,10)==0
        disp(counter);
        disp(pks);
    end
end
disp('ended')
%%
[control_cdf,control_x] = ecdf(sampled_events);
[toscale_cdf,toscale_x] = ecdf(toscale_events);

trunc_events = sampled_events(sampled_events >= detection_threshold);
trunc_toscale = toscale_events(toscale_events >= detection_threshold);

[ctrl_trunc_cdf, ctrl_trunc_x] = ecdf(trunc_events);
[toscale_trunc_cdf, toscale_trunc_x] = ecdf(trunc_toscale);

check_control=1;
if check_control
    lw = 3;
    setFigureDefaults;
    f1=figure();
    set(f1,'position',[.05 .1 .9 .4]);
%     subplot(1,3,1);
%     plot(control_x, control_cdf, 'linewidth', 2, 'color', 'k');
%     hold on;
%     plot(x,f,'linewidth',2,'linestyle','--','color','m');
%     plot(ctrl_trunc_x,ctrl_trunc_cdf,'linewidth',2,'linestyle','--','color','b');
%     set(gca,'xlim',[0 100]);
    subplot(1,3,3);
%     plot(control_x, control_cdf, 'linewidth', 2, 'color', 'k');
    hold on; box on
    plot(x,f,'linewidth',lw,'linestyle','-','color','m');
    plot(ctrl_trunc_x,ctrl_trunc_cdf,'linewidth',lw,'linestyle','-','color','k');
    line([5 5],[0 1],'linestyle','--','color','k');
    set(gca,'xlim',[0 30]);
    subplot(1,3,1);
    plot(control_x,control_cdf, 'linewidth', lw, 'color', 'k'); hold on;
    plot(toscale_x,toscale_cdf,'linewidth',lw,'linestyle','-','color','g');
    line([5 5],[0 1],'linestyle','--','color','k');
    set(gca,'xlim',[0 30]);
    subplot(1,3,2);
    plot(ctrl_trunc_x,ctrl_trunc_cdf, 'linewidth', lw, 'color', 'k'); hold on;
    plot(toscale_trunc_x,toscale_trunc_cdf,'linewidth',lw,'linestyle','-','color','g');
    line([5 5],[0 1],'linestyle','--','color','k');
    set(gca,'xlim',[0 30]);
    
    [~,p_ks] = kstest2(trunc_events,cdf_data)
    [~,p_ku] = kuiper_2samp(trunc_events,cdf_data)
    
    [~,p_ks_x] = kstest2(sampled_events,toscale_events)
    
end

%% TTX scaling

% use the sampled_events dist to scale and simulate a TTX-treated condition
scaling_factor = 1.20;
scaled_events = toscale_events .* scaling_factor;
[f_sc,x_sc] = ecdf(scaled_events);

trunc_scaled = scaled_events(scaled_events >= detection_threshold);
[ft_sc, xt_sc] = ecdf(trunc_scaled);

f2=figure();
set(f2,'position',[.1 .1 .9 .4]);
subplot(1,3,1);
plot(control_x, control_cdf, 'linewidth', lw, 'color', 'k'); hold on;
plot(x_sc, f_sc, 'linewidth', lw, 'color', 'r');
line([detection_threshold, detection_threshold],[0 1],'linestyle','--','color','k');
set(gca,'xlim',[0 30]);
subplot(1,3,2);
plot(ctrl_trunc_x, ctrl_trunc_cdf, 'linewidth', lw, 'color', 'k'); hold on;
plot(xt_sc, ft_sc, 'linewidth', lw, 'color', 'r');
line([detection_threshold, detection_threshold],[0 1],'linestyle','--','color','k');
set(gca,'xlim',[0 30]);


%% now setup ways to recover scaling factor



interp_x = linspace(0,1,2000);
ctlr1 = interp1(control_cdf, control_x, interp_x);
ttx2 = interp1(f_sc, x_sc, interp_x);

ctrl_sorted = sort(ctlr1);
ttx_sorted = sort(ttx2);

% first is the linear fit
[coeff,~] = polyfit(ctrl_sorted,ttx_sorted,1)
x_lin = 0:0.1:100;
y_lin = x_lin.*coeff(1) + coeff(2);

interp_x = linspace(0,1,2000);
vq1 = interp1(ctrl_trunc_cdf, ctrl_trunc_x, interp_x);
vq2 = interp1(ft_sc, xt_sc, interp_x);

ctrl_trunc_sort = sort(vq1);
ttx_trunc_sort = sort(vq2);

% first is the linear fit
[coeff_trunc,~] = polyfit(ctrl_trunc_sort,ttx_trunc_sort,1)
x_lin_t = 0:0.1:100;
y_lin_t = x_lin_t.*coeff_trunc(1) + coeff_trunc(2);

ax_lim = 60;

f3=figure();
set(f3,'position',[.02 .04 .9 .4]);
% subplot(2,2,1);
% plot(ctrl_sorted,ttx_sorted,'dk','markersize',8);
% hold on;
% line([0 100],[0 100],'linestyle','--','color','k');
% plot(x_lin,y_lin,'-m','linewidth',2);
% set(gca,'xlim',[0 ax_lim],'ylim',[0 ax_lim]);
% 
% scaled_ttx = (ttx_sorted - coeff(2)) ./ coeff(1);
% 
% subplot(2,2,2);
% hold on;
% c1=cdfplot(ctrl_sorted);
% c1.LineWidth = 2;
% c1.Color = [0 0 0];
% c2 = cdfplot(ttx_sorted);
% c2.LineWidth=2;
% c2.Color = 'm';
% c3 = cdfplot(scaled_ttx);
% c3.LineWidth=2;
% c3.Color = 'm';
% c3.LineStyle = '--';
% set(gca,'xlim',[0 ax_lim]);

subplot(1,3,1);
plot(ctrl_trunc_sort,ttx_trunc_sort,'ok','markersize',8);
hold on;
line([0 100],[0 100],'linestyle','--','color','k');
plot(x_lin_t,y_lin_t,'-r','linewidth',2);
set(gca,'xlim',[0 ax_lim],'ylim',[0 ax_lim]);

scaled_ttx_trunc = (ttx_trunc_sort - coeff_trunc(2)) ./ coeff_trunc(1);
scaled_ttx_trunc = scaled_ttx_trunc(scaled_ttx_trunc >= detection_threshold);


[c1_c, c1_x] = ecdf(ctrl_trunc_sort);
[c2_c, c2_x] = ecdf(ttx_trunc_sort);
[c3_c, c3_x] = ecdf(scaled_ttx_trunc);

subplot(1,3,2);
hold on;
plot(c1_x, c1_c, 'linestyle', '-', 'color', 'k', 'linewidth', 3);
plot(c2_x, c2_c, 'linestyle', '-', 'color', 'r', 'linewidth', 3);
plot(c3_x, c3_c, 'linestyle', '--', 'color', 'r', 'linewidth', 3);
set(gca,'xlim',[0 30]);
box on

%% now do the Kim et al 2012 method
clear pval*
sc_f = 0.70 : 0.001 : 1.10;


for ii = 1:numel(sc_f)
    
    sc_fact = sc_f(ii);
    
    
    ctrl = trunc_events;
    ttx = trunc_scaled;
    ttx_sc = trunc_scaled .* sc_fact;
    ttx_sc_thresh = ttx_sc(ttx_sc >= detection_threshold);
    
    [~,pval(ii)] = kstest2(ttx_sc,ctrl);
    [~,pval_thresh(ii)] = kstest2(ttx_sc_thresh,ctrl);
    
end

[max_pval,detected_sc_idx] = max(pval);
detected_sc = 1/sc_f(detected_sc_idx);

[max_pval_th,detected_sc_idx_th] = max(pval_thresh);
detected_sc_th = 1/sc_f(detected_sc_idx_th);

f4 = figure();
set(f4,'position',[.05 .08 .5 .7]);

plot(1./sc_f,pval,'-xk','markersize',9,'linewidth',1.5)
hold on;
line([detected_sc, detected_sc],[0 max_pval],'linestyle',':','color','r','linewidth',1.5);
text(detected_sc-0.05,max_pval+0.05,...
    sprintf('Scaling factor: %.3f',detected_sc),'fontsize',12,'color','r');
set(gca,'ylim',[0 round(max_pval + 0.1,1)]);


plot(1./sc_f,pval_thresh,'-xm','markersize',9,'linewidth',1.5)
hold on;
line([detected_sc_th, detected_sc_th],[0 max_pval_th],'linestyle',':','color','r','linewidth',1.5);
text(detected_sc_th-0.05,max_pval_th+0.05,...
    sprintf('Scaling factor: %.3f',detected_sc_th),'fontsize',12,'color','r');
set(gca,'ylim',[0 round(max_pval_th + 0.1,1)]);

%% kim et al with scaling up instead of down

clear pval*
sc_f = 0.95 : 0.001 : 1.40;


for ii = 1:numel(sc_f)
    
    sc_fact = sc_f(ii);
    
    
    ctrl = trunc_events;
    ttx = trunc_scaled;
    ctrl_sc = trunc_events .* sc_fact;
    ctrl_sc_thresh = ctrl_sc(ctrl_sc >= detection_threshold);
    
    [~,pval(ii)] = kstest2(ctrl_sc,ttx);
    [~,pval_thresh(ii)] = kstest2(ctrl_sc_thresh,ttx);
    
end

[max_pval,detected_sc_idx] = max(pval);
detected_sc = sc_f(detected_sc_idx);

[max_pval_th,detected_sc_idx_th] = max(pval_thresh);
detected_sc_th = sc_f(detected_sc_idx_th);

f4 = figure();
set(f4,'position',[.05 .08 .6 .6]);
subplot(1,2,1);
plot(sc_f,pval,'-xk','markersize',8,'linewidth',2)
hold on;
line([detected_sc, detected_sc],[0 max_pval],'linestyle',':','color','r','linewidth',1.5);
text(detected_sc-0.05,max_pval+0.05,...
    sprintf('Scaling factor: %.3f',detected_sc),'fontsize',12,'color','r');
set(gca,'ylim',[0 round(max_pval + 0.1,1)]);

subplot(1,2,2);
plot(sc_f,pval_thresh,'-xk','markersize',8,'linewidth',2)
hold on;
line([detected_sc_th, detected_sc_th],[0 max_pval_th],'linestyle',':','color','r','linewidth',1.5);
text(detected_sc_th-0.05,max_pval_th+0.05,...
    sprintf('Scaling factor: %.3f',detected_sc_th),'fontsize',12,'color','r');
set(gca,'ylim',[0 round(max_pval_th + 0.1,1)]);

%% finally do the Hanes et al analysis

ratio = ttx_trunc_sort./ctrl_trunc_sort;


f5 = figure();
set(f5,'position',[.07 .12 .5 .6]);
plot(ctrl_trunc_sort,ratio,'xk','markersize',8);
set(gca,'xlim',[0 40]);
    
    

%% save data for reproducibility

minisim.DATA.control1       = sampled_events;
minisim.DATA.control2       = toscale_events;
minisim.DATA.ttx            = scaled_events;
minisim.DATA.control1_trunc = trunc_events;
minisim.DATA.control2_trunc = trunc_toscale;
minisim.DATA.ttx_trunc      = trunc_scaled;
minisim.DATA.realdata       = cdf_data;

minisim.PARAMS.single_chan_conduct  = single_chan_conduct;
minisim.PARAMS.Vm                   = Vm;
minisim.PARAMS.single_current       = single_current;
minisim.PARAMS.n_AMPAR_min          = n_AMPAR_min;
minisim.PARAMS.n_AMPAR_max          = n_AMPAR_max;
minisim.PARAMS.mu_X                 = mu_X;
minisim.PARAMS.sigma_X              = sigma_X;
minisim.PARAMS.mu                   = mu;
minisim.PARAMS.sigma                = sigma;
minisim.PARAMS.LogNormal_PD         = LogNormal_PD;
minisim.PARAMS.detection_threshold  = detection_threshold;
minisim.PARAMS.n_cells              = n_cells;
minisim.PARAMS.n_events             = n_events;
minisim.PARAMS.n_AMPARs             = n_AMPARs;
minisim.PARAMS.pO_range             = pO_range;
minisim.PARAMS.prop_Open            = prop_Open;
minisim.PARAMS.mini_ampl            = mini_ampl;
minisim.PARAMS.prop_mode            = prop_mode;
minisim.PARAMS.n_samples_range      = n_samples_range;
minisim.PARAMS.scaling_factor       = scaling_factor;

savedir = '/Volumes/turrigiano-lab/ATP_MAIN/FIGURES/Neuron_REVIEW/mini_sim/final_version';
savename = 'minisim.mat';

save([savedir filesep savename],'minisim','-v7.3');







