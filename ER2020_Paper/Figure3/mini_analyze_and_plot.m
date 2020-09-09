%% Mini_analyze_and_plot
%
% Alejandro Torrado Pacheco - 2019
%
% This script analyzes and plots the processed mEPSC data from slice
% recordings done in early 2019. This code can be used to reproduce Fig. 3
% in the paper.
% 
% This script plots data for 4 conditions: CTRL ER2 & ER4, REOPEN ER2 &
% ER4.


% clear workspace, set default plotting params
clearIDE;
setFigureDefaults;

% load the data
if ismac
    load_dir = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig3';
elseif ispc
    local_dir = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig3';
end
ctrl2_file = [load_dir filesep 'CONTROL_DAT_er2.mat'];
ctrl4_file = [load_dir filesep 'CONTROL_DAT_er4.mat'];
reopen2_file = [load_dir filesep 'REOPEN_DAT_er2.mat'];
reopen4_file = [load_dir filesep 'REOPEN_DAT_er4.mat'];
fprintf('Loading data...\n');
CTRL2 = load(ctrl2_file);
CTRL4 = load(ctrl4_file);
REOPEN2 = load(reopen2_file);
REOPEN4 = load(reopen4_file);


% colors for plotting
ctrl_col = [.5 .5 .5];
reop_col = [7, 95, 179] ./ 255;
cset1 = [ctrl_col;ctrl_col;reop_col;reop_col];
cset2 = [ctrl_col;[1 1 1];reop_col;[1 1 1]];

%% Fig 3C - mEPSC Amplitude
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

% STATS
[p_a,~,a_stats] = kruskalwallis(a_kdat,a_groupz,'off')
% [p_a,~,a_stats] = anova1(a_kdat,a_groupz,'off')
c_a = multcompare(a_stats,'display','off')

amp_dat = padcat(C2_Amp, C4_Amp, R2_Amp, R4_Amp);
labels = {'ER2 Control','ER4 Control','ER2 Re-open','ER4 Re-open'};

amp_fig = plot_mini_data_4cond(amp_dat,'Ampl',labels,cset1,cset2);

%% Fig. S4A - Frequency
C2_Freq = CTRL2.CTRL_DAT.meanFreq;
C4_Freq = CTRL4.CTRL_DAT.meanFreq;
R2_Freq = REOPEN2.REOPEN_DAT.meanFreq;
R4_Freq = REOPEN4.REOPEN_DAT.meanFreq;

f_groupz = [repmat(1,size(C2_Freq)); repmat(2,size(R2_Freq));...
    repmat(3,size(C4_Freq)); repmat(4,size(R4_Freq))];
f_kdat = [C2_Freq; R2_Freq; C4_Freq; R4_Freq];

% STATS
[p_f,~,f_stats] = kruskalwallis(f_kdat,f_groupz,'off')
c_f = multcompare(f_stats,'display','off')

freq_dat = padcat(C2_Freq, C4_Freq, R2_Freq, R4_Freq);

freq_fig = plot_mini_data_4cond(freq_dat,'Freq',labels,cset1,cset2);

%% Fig. S4D - Input Resistance
C2_Rin = CTRL2.CTRL_DAT.Rin;
C4_Rin = CTRL4.CTRL_DAT.Rin;
R2_Rin = REOPEN2.REOPEN_DAT.Rin;
R4_Rin = REOPEN4.REOPEN_DAT.Rin;

rin_groupz = [repmat(1,size(C2_Rin)); repmat(2,size(R2_Rin));...
    repmat(3,size(C4_Rin)); repmat(4,size(R4_Rin))];
rin_kdat = [C2_Rin; R2_Rin; C4_Rin; R4_Rin];

% STATS
[p_rin,~,rin_stats] = kruskalwallis(rin_kdat,rin_groupz,'off')
c_rin = multcompare(rin_stats,'display','off')


rin_dat = padcat(C2_Rin, C4_Rin, R2_Rin, R4_Rin);

rin_fig = plot_mini_data_4cond(rin_dat,'Rin',labels,cset1,cset2);

%% Fig S4C - Capacitance
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

%% Fig S4B - Resting membrane potential
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

%% Fig 3D - CDF plots
% compile data
C2_CDF = CTRL2.CTRL_DAT.CDF;
C4_CDF = CTRL4.CTRL_DAT.CDF;
R2_CDF = REOPEN2.REOPEN_DAT.CDF;
R4_CDF = REOPEN4.REOPEN_DAT.CDF;

% cdf figure
cdf_fig = figure();
set(cdf_fig,'position',[.06 .1 .5 .7]);

hold on;

% control ER4
hc4 = cdfplot(C4_CDF);
hc4.Color = 'k';
hc4.LineWidth = 3;
hc4.LineStyle = '--';

% re-open ER2
hr2 = cdfplot(R2_CDF);
hr2.Color = reop_col;
hr2.LineWidth = 3;

% re-open ER4
hr4 = cdfplot(R4_CDF);
hr4.Color = reop_col;
hr4.LineWidth = 3;
hr4.LineStyle = '--';

% format plot
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
% The Kuiper test is more sensitive to differences away from the median
% than a Kolomgorov-Smirnov Test.
n_comparisons = 6;
[~, p12] = kuiper_2samp(C2_CDF, C4_CDF);
[~, p13] = kuiper_2samp(C2_CDF, R2_CDF);
[~, p14] = kuiper_2samp(C2_CDF, R4_CDF);
[~, p23] = kuiper_2samp(C4_CDF, R2_CDF);
[~, p24] = kuiper_2samp(C4_CDF, R4_CDF);
[~, p34] = kuiper_2samp(R2_CDF, R4_CDF);

% Bonferroni correction for multiple comparisons
p12 = p12 * n_comparisons;
p13 = p13 * n_comparisons;
p14 = p14 * n_comparisons;
p23 = p23 * n_comparisons;
p24 = p24 * n_comparisons;
p34 = p34 * n_comparisons;

% display the p-values
string_table = {'C2 vs C4';num2str(p12,'%.6f');...
    'C2 vs R2';num2str(p13,'%.6f');...
    'C2 vs R4';num2str(p14,'%.6f');...
    'C4 vs R2';num2str(p23,'%.6f');...
    'C4 vs R4';num2str(p24,'%.6f');...
    'R2 vs R4';num2str(p34,'%.6f')};

fprintf('\n________________________________\n');
fprintf('| K-S test + Bonferroni p-vals |');
fprintf('\n--------------------------------\n');
fprintf('%s\t|\tp = %s\n',string_table{:});



%% Figure S5A and S5B - Scaling test for re-opened hemisphere data
tt = figure();
set(tt,'position',[.06 .1 .5 .7]);

% x-values at which to evaluate linear fit
xq = 0:0.0005:1;

% ER2 and ER4 distributions
R4_nonan = R4_CDF(~isnan(R4_CDF));
R2_nonan = R2_CDF(~isnan(R2_CDF));

% empirical CDFs and interpolation
[cdf_R2,x_R2] = ecdf(R2_nonan);
if x_R2(1) == x_R2(2), x_R2(1) = x_R2(1) - 0.001; end
vq_R2 = interp1(cdf_R2, x_R2, xq);

[cdf_R4,x_R4] = ecdf(R4_nonan);
if x_R4(1) == x_R4(2), x_R4(1) = x_R4(1) - 0.001; end
vq_R4 = interp1(cdf_R4, x_R4, xq);

% ranked order plot
plot(sort(vq_R2),sort(vq_R4),'ko','markersize',10);
% linear fit
[fo,gof] = fit(sort(vq_R2)',sort(vq_R4)','poly1');
% plot linear fit
hold on;
X = 5:0.2:60;
Y = fo.p1 .* X + fo.p2;
plot(X,Y,'-k','linewidth',2);

% plot unity line
plot(1:70,1:70,'--','color',[.5 .5 .5],'linewidth',2);
% write linear fit equation
text(5,45,sprintf('Linear fit: F(x) = 0.877x + 0.318\nR^2 = %.3f',gof.adjrsquare),'fontsize',20);
% format plot
title('Re-opened hemisphere');
box off
xlabel('Amplitudes on ER2 (pA)');
ylabel('Amplitudes on ER4 (pA)');
set(gca,'xlim',[0 30],'ylim',[0 30]); % change the axes limits 



%% Fig 3E - Scaled Re-opened hemisphere CDFs
sc = figure();
set(sc,'position',[.06 .1 .5 .7]);

% plot the cumulative distributon for ER2
hold on;
[fr2,xr2] = ecdf(vq_R2);
hr2 = plot(xr2,fr2);
hr2.Color = [.5 .5 .5];
hr2.LineWidth = 4;

% plot CDF for ER4
[fr4,xr4] = ecdf(vq_R4);
hr4 = plot(xr4,fr4);
hr4.Color = 'k';
hr4.LineWidth = 4;
hr4.LineStyle = '-';

% scale ER2 according to fit found above 
R2_scaled = vq_R2 .* fo.p1 + fo.p2;

% plot scaled distribution
[fsr2,xsr2] = ecdf(R2_scaled);
hr4s = plot(xsr2,fsr2);
hr4s.Color = [.5 .5 .5];
hr4s.LineWidth = 4;
hr4s.LineStyle = '--';

% plot format
title('');
box off
grid off
legend([hr2,hr4,hr4s],{'ER2','ER4','ER2 scaled'},'Location','SouthEast');
legend boxoff
set(gca,'xlim',[5,30]);
ylabel('Fraction of events');
xlabel('mEPSC amplitude (pA)');

% test with KS and Kuiper test
[~,p01] = kstest2(vq_R2, vq_R4)
[~,p02] = kstest2(vq_R4, R2_scaled)
[~,pkui01] = kuiper_2samp(vq_R4, vq_R2)
[~,pkui02] = kuiper_2samp(R2_scaled, vq_R4)

%% As above but for Control ER2 vs Re-opened ER2 
xq = 0:0.0005:1;

% get CDFs
C2_nonan = C2_CDF(~isnan(C2_CDF));
R2_nonan = R2_CDF(~isnan(R2_CDF));

% interpolate
[cdf_R2,x_R2] = ecdf(R2_nonan);
if x_R2(1) == x_R2(2), x_R2(1) = x_R2(1) - 0.001; end
vq_R2 = interp1(cdf_R2, x_R2, xq);

[cdf_C2,x_C2] = ecdf(C2_nonan);
if x_C2(1) == x_C2(2), x_C2(1) = x_C2(1) - 0.001; end
vq_C2 = interp1(cdf_C2, x_C2, xq);

% linear fit and scaled control CDF
[fo_c,gof_c] = fit(sort(vq_C2)',sort(vq_R2)','poly1');
C2_scaled = vq_C2 .* fo_c.p1 + fo_c.p2;


%% Figure S5C & S5D
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
set(gca,'ylim',[0 30],'xlim',[0 30]); % change these values to get S5C or S5D



%% Fig 3F - scaled control ER2 vs re-open ER2
sc = figure();
set(sc,'position',[.06 .1 .5 .7]);


% plot distributions
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

% plot scaled control dist
[fsc2,xsc2] = ecdf(C2_scaled);
hr4s = plot(xsc2,fsc2);
hr4s.Color = 'k';
hr4s.LineWidth = 4;
hr4s.LineStyle = '--';

% axes format
title('');
box off
grid off
legend([hr2,hr4,hr4s],{'Control','Re-open','Control scaled'},'Location','SouthEast');
legend boxoff
set(gca,'xlim',[5,30]);
ylabel('Fraction of events');
xlabel('mEPSC amplitude (pA)');


% STATS - note kuiper test picks up difference in control scaled vs ER2
[~,p1] = kstest2(vq_C2, vq_R2)
[~,	p2] = kstest2(vq_R2, C2_scaled)
[~,pkui1] = kuiper_2samp(vq_C2, vq_R2)
[~,pkui2] = kuiper_2samp(C2_scaled, vq_R2)


ad_data = [C2_scaled'; vq_R2'];
ad_sample = [ones(size(C2_scaled')); 2.*ones(size(vq_R2'))];
[pad2,pad2_n] = adtest(C2_scaled, vq_R2)


%% MSE for C2 scaled vs R2 and R2_Scaled vs R4 - not plotted

MSE_c2 = nanmean((C2_scaled - vq_R2).^2);

MSE_r4 = nanmean((R2_scaled - vq_R4).^2);



