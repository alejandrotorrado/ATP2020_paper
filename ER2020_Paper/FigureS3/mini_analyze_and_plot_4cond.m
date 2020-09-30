clearIDE;
setFigureDefaults;

if ismac
%     local_dir = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/SliceRecordings/ANALYZED_DATA/slice_save_complete_June03_2019/June6_2019_DATA';
    local_dir = '/Users/atorrado/Desktop/MLS_DATA/slice_data';
elseif ispc
%     local_dir = 
end
ctrl2_file = [local_dir filesep 'CONTROL_DAT_er2.mat'];
ctrl4_file = [local_dir filesep 'CONTROL_DAT_er4.mat'];
reopen2_file = [local_dir filesep 'REOPEN_DAT_er2.mat'];
reopen4_file = [local_dir filesep 'REOPEN_DAT_er4.mat'];
fprintf('Loading data...\n');
CTRL2 = load(ctrl2_file);
CTRL4 = load(ctrl4_file);
REOPEN2 = load(reopen2_file);
REOPEN4 = load(reopen4_file);

%%

% colors for plotting
ctrl_col = [.5 .5 .5];
reop_col = [7, 95, 179] ./ 255;
cset1 = [ctrl_col;ctrl_col;reop_col;reop_col];
cset2 = [ctrl_col;[1 1 1];reop_col;[1 1 1]];

%% Amplitude
C2_Amp = CTRL2.CTRL_DAT.meanAmp;
C4_Amp = CTRL4.CTRL_DAT.meanAmp;
R2_Amp = REOPEN2.REOPEN_DAT.meanAmp;
R4_Amp = REOPEN4.REOPEN_DAT.meanAmp;
% keyboard;
n_C2 = sum(~isnan(C2_Amp));
n_C4 = sum(~isnan(C4_Amp));
n_R2 = sum(~isnan(R2_Amp));
n_R4 = sum(~isnan(R4_Amp));

a_groupz = [repmat(1,size(C2_Amp)); repmat(2,size(R2_Amp));...
    repmat(3,size(C4_Amp)); repmat(4,size(R4_Amp))];
a_kdat = [C2_Amp; R2_Amp; C4_Amp; R4_Amp];

[p_a,~,a_stats] = kruskalwallis(a_kdat,a_groupz,'off')
% [p_a,~,a_stats] = anova1(a_kdat,a_groupz,'off')
c_a = multcompare(a_stats,'display','off')

amp_dat = padcat(C2_Amp, C4_Amp, R2_Amp, R4_Amp);
labels = {'ER2 Control','ER4 Control','ER2 Re-open','ER4 Re-open'};

amp_fig = plot_mini_data_4cond(amp_dat,'Ampl',labels,cset1,cset2);

%% Frequency
C2_Freq = CTRL2.CTRL_DAT.meanFreq;
C4_Freq = CTRL4.CTRL_DAT.meanFreq;
R2_Freq = REOPEN2.REOPEN_DAT.meanFreq;
R4_Freq = REOPEN4.REOPEN_DAT.meanFreq;

f_groupz = [repmat(1,size(C2_Freq)); repmat(2,size(R2_Freq));...
    repmat(3,size(C4_Freq)); repmat(4,size(R4_Freq))];
f_kdat = [C2_Freq; R2_Freq; C4_Freq; R4_Freq];

[p_f,~,f_stats] = kruskalwallis(f_kdat,f_groupz,'off')
c_f = multcompare(f_stats,'display','off')

freq_dat = padcat(C2_Freq, C4_Freq, R2_Freq, R4_Freq);

freq_fig = plot_mini_data_4cond(freq_dat,'Freq',labels,cset1,cset2);

%% Input Resistance
C2_Rin = CTRL2.CTRL_DAT.Rin;
C4_Rin = CTRL4.CTRL_DAT.Rin;
R2_Rin = REOPEN2.REOPEN_DAT.Rin;
R4_Rin = REOPEN4.REOPEN_DAT.Rin;

rin_groupz = [repmat(1,size(C2_Rin)); repmat(2,size(R2_Rin));...
    repmat(3,size(C4_Rin)); repmat(4,size(R4_Rin))];
rin_kdat = [C2_Rin; R2_Rin; C4_Rin; R4_Rin];

[p_rin,~,rin_stats] = kruskalwallis(rin_kdat,rin_groupz,'off')
c_rin = multcompare(rin_stats,'display','off')


rin_dat = padcat(C2_Rin, C4_Rin, R2_Rin, R4_Rin);

rin_fig = plot_mini_data_4cond(rin_dat,'Rin',labels,cset1,cset2);

%% Capacitance
C2_Cp = CTRL2.CTRL_DAT.Cp;
C4_Cp = CTRL4.CTRL_DAT.Cp;
R2_Cp = REOPEN2.REOPEN_DAT.Cp;
R4_Cp = REOPEN4.REOPEN_DAT.Cp;

cp_groupz = [repmat(1,size(C2_Cp)); repmat(2,size(R2_Cp));...
    repmat(3,size(C4_Cp)); repmat(4,size(R4_Cp))];
cp_kdat = [C2_Cp; R2_Cp; C4_Cp; R4_Cp];

[p_cp,~,cp_stats] = kruskalwallis(cp_kdat,cp_groupz,'off')
c_rin = multcompare(cp_stats,'display','off')

cp_dat = padcat(C2_Cp, C4_Cp, R2_Cp, R4_Cp);

cp_fig = plot_mini_data_4cond(cp_dat,'Cp',labels,cset1,cset2);

%% Resting membrane potential
C2_Vr = CTRL2.CTRL_DAT.Vr;
C4_Vr = CTRL4.CTRL_DAT.Vr;
R2_Vr = REOPEN2.REOPEN_DAT.Vr;
R4_Vr = REOPEN4.REOPEN_DAT.Vr;

vr_groupz = [repmat(1,size(C2_Vr)); repmat(2,size(R2_Vr));...
    repmat(3,size(C4_Vr)); repmat(4,size(R4_Vr))];
vr_kdat = [C2_Vr; R2_Vr; C4_Vr; R4_Vr];

[p_vr,~,vr_stats] = kruskalwallis(vr_kdat,vr_groupz,'off')
c_vr = multcompare(vr_stats,'display','off')

vr_dat = padcat(C2_Vr, C4_Vr, R2_Vr, R4_Vr);

vr_fig = plot_mini_data_4cond(vr_dat,'Vr',labels,cset1,cset2);

%% CDF plots
C2_CDF = CTRL2.CTRL_DAT.CDF;
C4_CDF = CTRL4.CTRL_DAT.CDF;
R2_CDF = REOPEN2.REOPEN_DAT.CDF;
R4_CDF = REOPEN4.REOPEN_DAT.CDF;

% cdf figure
cdf_fig = figure();
set(cdf_fig,'position',[.06 .1 .5 .7]);

hold on;
% hc2 = cdfplot(C2_CDF);
% hc2.Color = 'k';
% hc2.LineWidth = 3;

hc4 = cdfplot(C4_CDF);
hc4.Color = 'k';
hc4.LineWidth = 3;
hc4.LineStyle = '--';

hr2 = cdfplot(R2_CDF);
hr2.Color = reop_col;
hr2.LineWidth = 3;

hr4 = cdfplot(R4_CDF);
hr4.Color = reop_col;
hr4.LineWidth = 3;
hr4.LineStyle = '--';

box off
grid off
legend([hc4,hr2,hr4],{'Control - ER4','Re-opened - ER2',...
    'Re-opened - ER4'},'Location','SouthEast');
legend boxoff
set(gca,'xlim',[5,30]);
xlabel('mEPSC amplitude (pA)');
ylabel('Fraction of events');
title('')

% stats for CDFs
% 1 is C2, 2 is C4, 3 is R2, 4 is R4
n_comparisons = 6;
[~, p12] = kuiper_2samp(C2_CDF, C4_CDF);
[~, p13] = kuiper_2samp(C2_CDF, R2_CDF);
[~, p14] = kuiper_2samp(C2_CDF, R4_CDF);
[~, p23] = kuiper_2samp(C4_CDF, R2_CDF);
[~, p24] = kuiper_2samp(C4_CDF, R4_CDF);
[~, p34] = kuiper_2samp(R2_CDF, R4_CDF);

% Bonferroni correction
p12 = p12 * n_comparisons;
p13 = p13 * n_comparisons;
p14 = p14 * n_comparisons;
p23 = p23 * n_comparisons;
p24 = p24 * n_comparisons;
p34 = p34 * n_comparisons;

string_table = {'C2 vs C4';num2str(p12,'%.4f');...
    'C2 vs R2';num2str(p13,'%.4f');...
    'C2 vs R4';num2str(p14,'%.4f');...
    'C4 vs R2';num2str(p23,'%.4f');...
    'C4 vs R4';num2str(p24,'%.4f');...
    'R2 vs R4';num2str(p34,'%.4f')};

fprintf('\n________________________________\n');
fprintf('| K-S test + Bonferroni p-vals |');
fprintf('\n--------------------------------\n');
fprintf('%s\t|\tp = %s\n',string_table{:});



%% reopemn cdf scaling
tt = figure();
set(tt,'position',[.06 .1 .5 .7]);

g = fittype('a+b*exp(+c*x)');
% R2_samp = randsample(R2_CDF(~isnan(R2_CDF)),numel(R4_CDF));


xq = 0:0.0005:1;

R4_nonan = R4_CDF(~isnan(R4_CDF));
R2_nonan = R2_CDF(~isnan(R2_CDF));

[cdf_R2,x_R2] = ecdf(R2_nonan);
if x_R2(1) == x_R2(2), x_R2(1) = x_R2(1) - 0.001; end
vq_R2 = interp1(cdf_R2, x_R2, xq);


[cdf_R4,x_R4] = ecdf(R4_nonan);
if x_R4(1) == x_R4(2), x_R4(1) = x_R4(1) - 0.001; end
vq_R4 = interp1(cdf_R4, x_R4, xq);


plot(sort(vq_R2),sort(vq_R4),'ko','markersize',10);
[fo,gof] = fit(sort(vq_R2)',sort(vq_R4)','poly1');
% [ef,gof_e] = fit(sort(R2_samp),sort(R4_CDF),g,'StartPoint',[0 0 0.1]);%[[ones(size(R2_samp)), -exp(-R2_samp)]\R4_CDF; 1]);
hold on;
X = 5:0.2:60;
Y = fo.p1 .* X + fo.p2;

plot(X,Y,'-k','linewidth',2);
% plot(X,ef(X),'-r','linewidth',2);
plot(1:70,1:70,'--','color',[.5 .5 .5],'linewidth',2);
text(5,45,sprintf('Linear fit: F(x) = 0.877x + 0.318\nR^2 = %.3f',gof.adjrsquare),'fontsize',20);
% text(5,55,sprintf('Exponential fit: F(x) = %.3f - %.3f*exp(-%.3fx)\nR^2 = %.3f',ef.a,ef.b,ef.c,gof_e.adjrsquare),'fontsize',20,'color','r');
title('Re-opened hemisphere');
box off
xlabel('Amplitudes on ER2 (pA)');
ylabel('Amplitudes on ER4 (pA)');
set(gca,'xlim',[0 30],'ylim',[0 30]);



%%
sc = figure();
set(sc,'position',[.06 .1 .5 .7]);



hold on;
[fr2,xr2] = ecdf(vq_R2);
hr2 = plot(xr2,fr2);
hr2.Color = [.5 .5 .5];
hr2.LineWidth = 4;

[fr4,xr4] = ecdf(vq_R4);
hr4 = plot(xr4,fr4);
hr4.Color = 'k';
hr4.LineWidth = 4;
hr4.LineStyle = '-';


R2_scaled = vq_R2 .* fo.p1 + fo.p2;

[fsr2,xsr2] = ecdf(R2_scaled);
hr4s = plot(xsr2,fsr2);
hr4s.Color = [.5 .5 .5];
hr4s.LineWidth = 4;
hr4s.LineStyle = '--';


title('');
box off
grid off
legend([hr2,hr4,hr4s],{'ER2','ER4','ER2 scaled'},'Location','SouthEast');
legend boxoff
set(gca,'xlim',[5,30]);
ylabel('Fraction of events');
xlabel('mEPSC amplitude (pA)');

[~,p01] = kstest2(vq_R2, vq_R4)
[~,p02] = kstest2(vq_R4, R2_scaled)
[~,pkui01] = kuiper_2samp(vq_R4, vq_R2)
[~,pkui02] = kuiper_2samp(R2_scaled, vq_R4)

%% CTRL vs ER2 
% n_iter = 1000;
% all_pkui2 = nan(n_iter,1);
% for k = 1:n_iter
    
%     disp(k);

xq = 0:0.0005:1;

C2_nonan = C2_CDF(~isnan(C2_CDF));
R2_nonan = R2_CDF(~isnan(R2_CDF));

[cdf_R2,x_R2] = ecdf(R2_nonan);
if x_R2(1) == x_R2(2), x_R2(1) = x_R2(1) - 0.001; end
vq_R2 = interp1(cdf_R2, x_R2, xq);


[cdf_C2,x_C2] = ecdf(C2_nonan);
if x_C2(1) == x_C2(2), x_C2(1) = x_C2(1) - 0.001; end
vq_C2 = interp1(cdf_C2, x_C2, xq);

% C2_samp = randsample(C2_nonan,numel(R2_nonan));
% [fo_c,gof_c] = fit(sort(C2_samp),sort(R2_nonan),'poly1');
% C2_scaled = C2_CDF .* fo_c.p1 + fo_c.p2;

[fo_c,gof_c] = fit(sort(vq_C2)',sort(vq_R2)','poly1');
C2_scaled = vq_C2 .* fo_c.p1 + fo_c.p2;


    %%
tt = figure();
set(tt,'position',[.06 .1 .5 .7]);
plot(sort(vq_C2),sort(vq_R2),'ko','markersize',10);

g = fittype('a-b*exp(-c*x)');
[ef,gof_e] = fit(sort(vq_C2)',sort(vq_R2)',g,'StartPoint',[0 0 1]);% [[ones(size(C2_samp)), -exp(-C2_samp)]\R2_nonan; 1]);
hold on;
X = 5:0.2:80;
Y = fo_c.p1 .* X + fo_c.p2;
plot(X,Y,'-k','linewidth',2);
plot(X,ef(X),'-r','linewidth',2);
plot(1:70,1:70,'--','color',[.5 .5 .5],'linewidth',2);
% plot(1:70,1:70,'--','color',[.4 .4 .4],'linewidth',2);
text(10,55,sprintf('Linear fit: F(x) = %.3fx + %.3f\nR^2 = %.3f',fo_c.p1,fo_c.p2,gof_c.adjrsquare),'fontsize',20);
text(10,65,sprintf('Exponential fit: F(x) = %.3f - %.3f*exp(-%.3fx)\nR^2 = %.3f',ef.a,ef.b,ef.c,gof_e.adjrsquare),'fontsize',20,'color','r');
box off
xlabel('Control - Amplitudes on ER2 (pA)');
ylabel('Re-open - Amplitudes on ER2 (pA)');
set(gca,'ylim',[0 30],'xlim',[0 30]);



%%
sc = figure();
set(sc,'position',[.06 .1 .5 .7]);



hold on;
[fc2,xc2] = ecdf(vq_C2);
hr2 = plot(xc2,fc2);
hr2.Color = 'k';
hr2.LineWidth = 4;

[fr2,xr2] = ecdf(vq_R2);
hr4 = plot(xr2,fr2);
hr4.Color = [.6 .6 .6];
hr4.LineWidth = 4;
hr4.LineStyle = '-';

[fsc2,xsc2] = ecdf(C2_scaled);
hr4s = plot(xsc2,fsc2);
hr4s.Color = 'k';
hr4s.LineWidth = 4;
hr4s.LineStyle = '--';


title('');
box off
grid off
legend([hr2,hr4,hr4s],{'Control','Re-open','Control scaled'},'Location','SouthEast');
legend boxoff
set(gca,'xlim',[5,30]);
ylabel('Fraction of events');
xlabel('mEPSC amplitude (pA)');



[~,p1] = kstest2(vq_C2, vq_R2)
[~,	p2] = kstest2(vq_R2, C2_scaled)
[~,pkui1] = kuiper_2samp(vq_C2, vq_R2)
[~,pkui2] = kuiper_2samp(C2_scaled, vq_R2)

% allpks2(k) = p2;
% all_pkui2(k) = pkui2;

% end
% 
% sorted_pkui2 = sort(all_pkui2);
% median(all_pkui2)
% 
% idx05 = 25;
% idx95 = 975;
% 
% p05 = sorted_pkui2(idx05)
% p95 = sorted_pkui2(idx95)
% 
% 
% 
% sorted_pks2 = sort(allpks2);
% median(allpks2)
% 
% pks05 = sorted_pks2(idx05)
% pks95 = sorted_pks2(idx95)
% 

%% MSE for C2 scaled vs R2

MSE_c2 = nanmean((C2_scaled - vq_R2).^2)

MSE_r4 = nanmean((R2_scaled - vq_R4).^2)

%%

R2_AMPS = REOPEN2.REOPEN_DAT.allAmps;
R4_AMPS = REOPEN4.REOPEN_DAT.allAmps;

n_quantiles = 100;
p_q = linspace(0,1,n_quantiles);

% vectorized code to create quantile-sampled mEPSC distribution

% first eliminate cells with fewer mEPSC recorded than n_quantile samples
n_minis_R2 = cellfun(@numel,R2_AMPS);
bad_idx_R2 = find(n_minis_R2 < n_quantiles);

n_minis_R4 = cellfun(@numel,R4_AMPS);
bad_idx_R4 = find(n_minis_R4 < n_quantiles);

R2_AMPS(bad_idx_R2) = [];
R4_AMPS(bad_idx_R4) = [];

n_R2 = numel(R2_AMPS);
n_R4 = numel(R4_AMPS);

% then concatenate mEPSC amplitudes for all neurons, padding with NaNs
all_R2_amps = padcat(R2_AMPS{:});
all_R4_amps = padcat(R4_AMPS{:});

% then sample using the "quantile" function (which ignores NaN). Note that
% the array of padded concatenated amplitudes needs to be transposed.
R2_samples = quantile(all_R2_amps',p_q);
R4_samples = quantile(all_R4_amps',p_q);

% estimate CDF using ecdf function
[R2_S_CDF,R2_S_X] = ecdf(R2_samples(:));
[R4_S_CDF,R4_S_X] = ecdf(R4_samples(:));

% plot CDFs
lw = 4;
scf = figure();
set(scf,'position',[0.1 0.1 0.5 0.7]);
hold on;
plot(R2_S_X,R2_S_CDF,'-','color',reop_col,'linewidth',lw);
plot(R4_S_X,R4_S_CDF,':','color','k','linewidth',lw);
set(gca,'xlim',[0 30]);

[~,p] = kstest2(R2_samples(:),R4_samples(:))
[~,pk] = kuiper_2samp(R2_samples(:), R4_samples(:))

%% Hanes et al 2020 analysis
% now use the method of Hanes et al., 2020 (https://doi.org/10.1523/JNEUROSCI.1393-19.2020)
% to check if multiplicative scaling explains the difference in
% distributions

% To do this, need to re-sample the distributions so they have the same
% number of pts. For this I will use 2xnR4 quantiles for R2 and 2xnR2
% quantiles for R4. This is because otherwise i have less than 30
% quantiles.

% get the data
R2_AMPS = REOPEN2.REOPEN_DAT.allAmps;
R4_AMPS = REOPEN4.REOPEN_DAT.allAmps;

% Eliminate neurons with less than 40 mEPSC amplitudes (set by q_threshold)
q_threshold = 100;
q_factor = 4;

n_minis_R2 = cellfun(@numel,R2_AMPS);
bad_idx_R2 = find(n_minis_R2 < q_threshold);

n_minis_R4 = cellfun(@numel,R4_AMPS);
bad_idx_R4 = find(n_minis_R4 < q_threshold);

R2_AMPS(bad_idx_R2) = [];
R4_AMPS(bad_idx_R4) = [];

% find number of remaining neurons and compute quantiles accordingly. This
% ensures that both distributions will have same number of samples. For
% example, with q_factor = 2:
% n_R2 * (2*n_R4) = 1100 for R2
% n_R4 * (2*n_R2) = 1100 for R4
n_R2 = numel(R2_AMPS);
n_R4 = numel(R4_AMPS);

q_R2 = q_factor*n_R4;
q_R4 = q_factor*n_R2;

p_R2 = linspace(0,1,q_R2);
p_R4 = linspace(0,1,q_R4);

%  concatenate mEPSC amplitudes for all neurons, padding with NaNs
all_R2_amps = padcat(R2_AMPS{:});
all_R4_amps = padcat(R4_AMPS{:});

% then sample using the "quantile" function (which ignores NaN). Note that
% the array of padded concatenated amplitudes needs to be transposed.
R2_samples = quantile(all_R2_amps',p_R2);
R4_samples = quantile(all_R4_amps',p_R4);

% rank the pooled data
R2_ranked = sort(R2_samples(:)); % this ranks from low to high
R4_ranked = sort(R4_samples(:));

% then calculate ratio
Ratio_ranked = R4_ranked./R2_ranked; % here R4 is the "treatment" distribution, putatively downscaled

% plot ratio vs R2 amplitude for each pair
ratio_fig = figure();
set(ratio_fig,'position',[0.15 0.08 0.5 0.7]);
hold on;
plot(R2_ranked,Ratio_ranked,'ok','linestyle','none');
set(gca,'xlim',[0 40]);


% also do rank-order linear fit and plot scaled CDF
[fo_c,gof_c] = fit(R2_ranked,R4_ranked,'poly1');
R2_scaled = R2_ranked .* fo_c.p1 + fo_c.p2;
R2_scaled_noint = R2_ranked .* fo_c.p1;

[R2_S_CDF,R2_S_X] = ecdf(R2_ranked(:));
[R4_S_CDF,R4_S_X] = ecdf(R4_ranked(:));
[R2_S1_CDF,R2_S1_X] = ecdf(R2_scaled(:));
[R2_S2_CDF,R2_S2_X] = ecdf(R2_scaled_noint(:));

% plot CDFs
lw = 2;
scf = figure();
set(scf,'position',[0.1 0.1 0.5 0.7]);
hold on;
plot(R2_S_X,R2_S_CDF,'-','color',reop_col,'linewidth',lw);
plot(R4_S_X,R4_S_CDF,'--','color','k','linewidth',lw);
plot(R2_S1_X,R2_S1_CDF,':','color','g','linewidth',lw);
plot(R2_S2_X,R2_S2_CDF,':','color','m','linewidth',lw);
set(gca,'xlim',[0 30]);

% p-values with k-s test
[~,p1] = kstest2(R2_ranked,R4_ranked);
[~,p2] = kstest2(R2_scaled,R4_ranked);
[~,p3] = kstest2(R2_scaled_noint,R4_ranked);


rank_fig = figure();
set(ratio_fig,'position',[0.15 0.08 0.5 0.7]);
hold on;
plot(R2_ranked,R4_ranked,'ok','linestyle','none');
line([0 40],[0 40],'linestyle','--','color','m','linewidth',2);
set(gca,'xlim',[0 40]);

%% Kim et al., 2012 (https://doi.org/10.1371/journal.pone.0037364)
% In this paper they described an iterative procedure to find the scaling
% factor. Idea is to vary the scaling factor, scale the distributions,
% compare them using a 2-sample KS test, then choose the scaling factor
% that yields the largest p-value. In addition, you correct for detection
% thresholds by removing all the data points below the threshold in the
% scaled distribution

% Choose min and max scaling factors, iteration step, and detection
% threshold. Based on plot above, "real" scaling factor is around 0.85.
scaling_min = 0.5;
scaling_max = 3.0;
scaling_step = 0.001;
detection_thresh = 5; % pA

% shorten distribution names to make the code more readable
R2 = R2_ranked;
R4 = R4_ranked;

% iterate over scaling factors
iteration = 0; % counter for iterations
for sc = scaling_min : scaling_step : scaling_max
    
    % update counter
    iteration = iteration + 1;
    
    % scaled distribution
    R2_sc = sc .* R2;
    R4_sc = R4;
    
    % remove values below detection threshold (duplicate dists to have both
    % sets of results)
    R2_det = R2_sc(R2_sc >= detection_thresh);
    R4_det = R4_sc(R4_sc >= detection_thresh);
    
    
    % do the KS test
    [~,p_temp] = kstest2(R2_sc,R4_sc);
    [~,p_det_temp] = kstest2(R2_det,R4_det);
    
    % save the p-value and the corresponding scaling factor
    p_values(iteration,1) = p_temp;
    p_values(iteration,2) = p_det_temp;
    p_values(iteration,3) = sc;
    
end

% find the scaling factor with the maximal p-value
[max_pval,max_sc_idx] = max(p_values(:,1));
[max_pval_det,max_sc_det_idx] = max(p_values(:,2));

max_sc = p_values(max_sc_idx,3);
max_sc_det = p_values(max_sc_det_idx,3);

% plot data
iterfig = figure();
set(iterfig,'position',[.05 .1 .5 .7]);
hold on
plot(p_values(:,3),p_values(:,1),'-xk','linewidth',2,'markersize',8);
plot(p_values(:,3),p_values(:,2),'-xm','linewidth',2,'markersize',8);
text(max_sc,max_pval+0.03,num2str(max_sc),'color','k','fontsize',20);
text(max_sc_det,max_pval_det+0.03,num2str(max_sc_det),'color','m','fontsize',20);
set(gca,'xlim',[0.85 1.0],'ylim',[0,0.7]);
ylabel('p-value');
xlabel('scaling factor');

% now use the two scaling factors found this way to scale the CDFs
R4_scaled_iter = R4 ./ max_sc;
R4_scaled_iter_det = R4 ./ max_sc_det;

[R2_S_CDF,R2_S_X] = ecdf(R2_ranked(:));
[R4_S_CDF,R4_S_X] = ecdf(R4_ranked(:));
[R4_S1_CDF,R4_S1_X] = ecdf(R4_scaled_iter(:));
[R4_S2_CDF,R4_S2_X] = ecdf(R4_scaled_iter_det(:));

% plot CDFs
lw = 2;
scf = figure();
set(scf,'position',[0.1 0.1 0.5 0.7]);
hold on;
plot(R2_S_X,R2_S_CDF,'-','color',reop_col,'linewidth',lw);
plot(R4_S_X,R4_S_CDF,'--','color','k','linewidth',lw);
plot(R4_S1_X,R4_S1_CDF,':','color','g','linewidth',lw);
plot(R4_S2_X,R4_S2_CDF,':','color','m','linewidth',lw);
set(gca,'xlim',[0 30]);

% p-values with k-s test
[~,p1] = kstest2(R2_ranked,R4_ranked)
[~,p2] = kstest2(R4_scaled_iter,R2_ranked)
[~,p3] = kstest2(R4_scaled_iter_det,R2_ranked)



%% Hanes et al for C2 vs R2
% As above, but for the control (ER2) vs reopen (ER2) comparison

% get the data
R2_AMPS = REOPEN2.REOPEN_DAT.allAmps;
C2_AMPS = CTRL2.CTRL_DAT.allAmps;

% Eliminate neurons with less than 40 mEPSC amplitudes (set by q_threshold)
q_threshold = 100;
q_factor = 4;

n_minis_R2 = cellfun(@numel,R2_AMPS);
bad_idx_R2 = find(n_minis_R2 < q_threshold);

n_minis_C2 = cellfun(@numel,C2_AMPS);
bad_idx_C2 = find(n_minis_C2 < q_threshold);

R2_AMPS(bad_idx_R2) = [];
C2_AMPS(bad_idx_C2) = [];

% find number of remaining neurons and compute quantiles accordingly. This
% ensures that both distributions will have same number of samples. For
% example, with q_factor = 2:
% n_R2 * (2*n_R4) = 1100 for R2
% n_R4 * (2*n_R2) = 1100 for R4
n_R2 = numel(R2_AMPS);
n_C2 = numel(C2_AMPS);

q_R2 = q_factor*n_C2;
q_C2 = q_factor*n_R2;

p_R2 = linspace(0,1,q_R2);
p_C2 = linspace(0,1,q_C2);

%  concatenate mEPSC amplitudes for all neurons, padding with NaNs
all_R2_amps = padcat(R2_AMPS{:});
all_C2_amps = padcat(C2_AMPS{:});

% then sample using the "quantile" function (which ignores NaN). Note that
% the array of padded concatenated amplitudes needs to be transposed.
R2_samples = quantile(all_R2_amps',p_R2);
C2_samples = quantile(all_C2_amps',p_C2);

% rank the pooled data
R2_ranked = sort(R2_samples(:)); % this ranks from low to high
C2_ranked = sort(C2_samples(:));

% then calculate ratio
Ratio_ranked = R2_ranked./C2_ranked; % here R2 is the "treatment" distribution, putatively potentiated

% plot ratio vs R2 amplitude for each pair
ratio_fig = figure();
set(ratio_fig,'position',[0.15 0.08 0.5 0.7]);
hold on;
plot(C2_ranked,Ratio_ranked,'ok','linestyle','none','markersize',10);
set(gca,'xlim',[0 40]);

%% Kim et al 2012 method for C2 vs R2
% form plot above, putative scaling factor is around 1.05

% Choose min and max scaling factors, iteration step, and detection
% threshold. Based on plot above, "real" scaling factor is around 0.85.
scaling_min = 0.5;
scaling_max = 3.0;
scaling_step = 0.001;
detection_thresh = 5; % pA

% shorten distribution names to make the code more readable
R2 = R2_ranked;
C2 = C2_ranked;

% iterate over scaling factors
iteration = 0; % counter for iterations
for sc = scaling_min : scaling_step : scaling_max
    
    % update counter
    iteration = iteration + 1;
    
    % scaled distribution
    R2_sc = sc .* R2;
    
    % remove values below detection threshold (duplicate dists to have both
    % sets of results)
    R2_det = R2_sc(R2_sc >= detection_thresh);
    C2_det = C2(C2 >= detection_thresh);
    
    
    % do the KS test
    try
    [~,p_temp] = kstest2(R2_sc,C2);
    [~,p_det_temp] = kstest2(R2_det,C2_det);
    catch
    keyboard
    end
    
    % save the p-value and the corresponding scaling factor
    p_values(iteration,1) = p_temp;
    p_values(iteration,2) = p_det_temp;
    p_values(iteration,3) = sc;
    
end

% find the scaling factor with the maximal p-value
[max_pval,max_sc_idx] = max(p_values(:,1));
[max_pval_det,max_sc_det_idx] = max(p_values(:,2));

max_sc = p_values(max_sc_idx,3);
max_sc_det = p_values(max_sc_det_idx,3);

% plot data
iterfig2 = figure();
set(iterfig2,'position',[.05 .14 .5 .7]);
hold on
plot(p_values(:,3),p_values(:,1),'-xk','linewidth',2,'markersize',8);
plot(p_values(:,3),p_values(:,2),'-xb','linewidth',2,'markersize',8);
text(max_sc,max_pval+0.03,num2str(max_sc),'color','k','fontsize',20);
text(max_sc_det,max_pval_det+0.03,num2str(max_sc_det),'color','m','fontsize',20);
set(gca,'xlim',[.85 1.0],'ylim',[0,1.0]);
ylabel('p-value');
xlabel('scaling factor');