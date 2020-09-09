%% load data
clearvars -except *CONTCELL*
clc

rload = load('/Users/atorrado/Desktop/MLS_DATA/recov_bootstrap_2.mat');
recov_bstrap_small = rload.recov_bootstrap;

frload = load('/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig1/recov_analysis_ERDATA.mat');
rdat = frload.recov_analysis;
dep_rsu = rdat.DEPRIVED.RSU_idx;
ctrl_rsu = rdat.CONTROL.RSU_idx;
dep_rate = rdat.DEPRIVED.RSU_FRbycell;
ctrl_rate = rdat.CONTROL.RSU_FRbycell;
all_rsu = [ctrl_rsu'; dep_rsu'];







meanWF      = recov_bstrap_small.meanWF;
rsu_idx     = recov_bstrap_small.rsu_idx;
G_bin       = recov_bstrap_small.G_bin;
day_list    = recov_bstrap_small.day_list;


bad_cells = [4,9,16,18,21,24,54,63,71,72,73];
no_badz = 0;

if no_badz
   rsu_idx(bad_cells) = [];
   meanWF(bad_cells) = [];
end

N = numel(rsu_idx);
Nd = numel(day_list);



do_check = 0;

%% calculate data z

for j = 1:N
    
    for d = 1:(Nd-1)
        clear tmp_mse tmp_mse_p wf1* wf2*
        day1 = day_list(d) * 24*3600 / G_bin;
        day2 = day_list(d+1) * 24*3600 / G_bin;
        dw = 24*3600 / G_bin;
        
        % compare each day to the next
        wf1 = meanWF{j}(d,:);
        wf2 = meanWF{j}(d+1,:);
        wf1_p = wf1 ./ abs(min(wf1));
        wf2_p = wf2 ./ abs(min(wf2));
        
        tmp_mse = nanmean((wf1 - wf2).^2);
        tmp_mse_p = nanmean((wf1_p - wf2_p).^2);
        
%         wf_align_p = [wf1_p;wf2_p];
%         [~,minidx] = min(wf_align_p,[],2);
%         nsw = [-min(minidx)+2 : size(wf1_p,2) - min(minidx) - abs(diff(minidx)) + 1];
%         wf_aligned_p = alignrows(wf_align_p, minidx, nsw, 96);
%         wf1_a_p = wf_aligned_p(1,:);
%         wf2_a_p = wf_aligned_p(2,:);
%         
%         
%         figure(); subplot(1,2,1); hold on;
%         plot(wf1_p); plot(wf2_p)
%         subplot(1,2,2); hold on;
%         plot(wf1_a_p); plot(wf2_a_p)
%         
%         mse_a = nanmean((wf1_a_p - wf2_a_p).^2); 
%         
        mse(j,d) = tmp_mse;
        
        mse_peak(j,d) = tmp_mse_p;
        
        if do_check
            if tmp_mse_p > 0.0032
                ff1=figure();
                set(gcf,'unit','normalized','color','w','position',[.1 .1 .85 .8]);
                subplot(1,2,1); hold on;
                plot(wf1,'k','linewidth',2);
                plot(wf2,'r','linewidth',2);
                text(20,40,sprintf('Error: %.2f - cell %u, day %u.\n',tmp_mse,j,d),'fontsize',20);
                subplot(1,2,2); hold on;
                plot(wf1_p,'k','linewidth',2);
                plot(wf2_p,'r','linewidth',2);
                text(55,-0.75,sprintf('Error: %.5f - cell %u, day %u.\n',tmp_mse_p,j,d),'fontsize',20);
                
                ff2 = figure(); hold on;
                if j<=numel(ctrl_rsu)
                    title('control','fontsize',20);
                    thisrate = ctrl_rate(j,:);
                else
                    title('deprived','fontsize',20);
                    thisrate = dep_rate(j-numel(ctrl_rsu),:);
                end
                plot(thisrate);
                set(gca,'xtick',6.5*24*3600/G_bin : 24*3600/G_bin : 11.5*24*3600/G_bin,...
                    'xticklabel',{'MD4','ER1','ER2','ER3','ER4'});
                yl = get(gca,'ylim');
                ON = CONTCELL_recov.MASTER(rsu_idx(j)).onTime;
                OFF = CONTCELL_recov.MASTER(rsu_idx(j)).offTime;
                l_on = line([ON ON],yl,'color','k');
                l_off = line([ON ON],yl,'color','k');
                r1 = rectangle('position',[day1 yl(1) dw yl(2)]);
                r1.EdgeColor = 'none';
                r1.FaceColor = [0 0 0 0.5];
                r2 = rectangle('Position',[day2 yl(1) dw yl(2)]);
                r2.EdgeColor = 'none';
                r2.FaceColor = [1 0 0 0.5];
                
                uiwait(ff1);
            end
        end
        
    end
    
end


%% do the permutation test

n_this = 3;
% n_other = 3;
M = 1000;
verbose = 0;

% for each cell, pick 3 WFs and two others at random from another cell
% repeat M times per cell (see variable M above)
% note: n_this and n_other control how many of this vs other cells WFs
% are picked at random

zw = 1;
zbins = 0:zw:300;
zctr = zbins(1:end-1) + zw/2;
z = max(mse,[],2);
zc = histcounts(z,zbins);

z_peak = max(mse_peak,[],2);
zwp = 0.0001;
zbinsp = 0:zwp:0.03;
zctrp = zbinsp(1:end-1) + zwp/2;
zcp = histcounts(z_peak,zbinsp);

% loop through runs
% for every run, calculate the z and z_sh and find the 95% confidence bin
for m = 1:M
    if mod(m,20) == 0
        fprintf('Run %u of %u.\n',m,M);
    end
    % loop through cells
    for k = 1:N
%         fprintf('\nCreating shuffle dataset. Cell %u of %u.\n',k,N);
        
        clear this_ix other_ix
        this_ix = sort(randsample([1:5],n_this));
        other_ix = setdiff([1:5],this_ix);
        % set of available cells to choose from
        rand_set = setdiff([1:N],k);
        rand_cell = randsample(rand_set,1);
        
        
        meanWF_shuffle{k,m}(this_ix,:) = meanWF{k}(this_ix,:);
        meanWF_shuffle{k,m}(other_ix,:) = meanWF{rand_cell}(other_ix,:);
        
        if do_check
            figure(999);
            plot(meanWF_shuffle{k,m}(other_ix(1),:)); hold on;
            plot(meanWF_shuffle{k,m}(other_ix(2),:));
            uiwait(gcf);
        end
        
        for dd = 1:(Nd-1)
            clear tmp_mse tmp_mse_p
            % compare each day to the next
            wf1 = meanWF_shuffle{k,m}(dd,:);
            wf2 = meanWF_shuffle{k,m}(dd+1,:);
            
            wf1_p = wf1 ./ abs(min(wf1));
            wf2_p = wf2 ./ abs(min(wf2));
            
            
            
            tmp_mse = nanmean((wf1 - wf2).^2);
            tmp_mse_p = nanmean((wf1_p - wf2_p).^2);
            
            mse_shuffle{m}(k,dd) = tmp_mse;
            mse_peak_shuffle{m}(k,dd) = tmp_mse_p;
        end
    end
    
    z_sh{m} = max(mse_shuffle{m},[],2);
    zc_sh{m} = histcounts(z_sh{m},zbins);
    
    z_sh_peak{m} = max(mse_peak_shuffle{m},[],2);
    zc_shp{m} = histcounts(z_sh_peak{m},zbinsp);

    
    for bb = 1:size(zbins,2)
        conf_per_bin{m}(bb,1) = 100 * (1 - sum(z_sh{m}<zbins(bb))/sum(z<zbins(bb)));
    end

    for bb = 1:size(zbinsp,2)
        conf_per_bin_peak{m}(bb,1) = 100 * (1 - sum(z_sh_peak{m}<zbinsp(bb))/sum(z_peak<zbinsp(bb)));
    end

    % find threshold with x% confidence
    conf_thresh = 95;
    thisbin = find(conf_per_bin{m} >= conf_thresh,1,'last');
    if isempty(thisbin), thisbin = 1; end
    mse_thresh(m,1) = zbins(thisbin);
    thisbin_p = find(conf_per_bin_peak{m} >= conf_thresh,1,'last');
    if isempty(thisbin_p)
        peak_thresh(m,1) = NaN;
    else
        peak_thresh(m,1) = zbinsp(thisbin_p);
    end
    clear thisbin thisbin_p
    
    data_under = sum(z <= mse_thresh(m,1));
    shuffle_under = sum(z_sh{m} <= mse_thresh(m,1));
    conf = 100 * (1 - shuffle_under/data_under);
    
    
    datapeak_under = sum(z_peak <= peak_thresh(m,1));
    shufflepeak_under = sum(z_sh_peak{m} <= peak_thresh(m,1));
    conf_p = 100 * (1 - shufflepeak_under/datapeak_under);
    if verbose
        fprintf('Run %u. MSE <= %.1f. N(data) = %u; N(shuffle) = %u.\nConfidence: %.2f\n',...
            m, mse_thresh(m,1), data_under, shuffle_under, conf);
        fprintf('MSE_peak <= %.4f. N(data) = %u; N(shuffle) = %u.\nConfidence: %.2f\n',...
            peak_thresh(m,1), datapeak_under, shufflepeak_under, conf_p);
    end
    
end

% find bin in z_sh_peak where prob < 5% that number is from that dist
z_sh_1 = cellfun(@(x) max(x,[],2),mse_peak_shuffle,'uniformoutput',false);
z_sh_2 = cell2mat(z_sh_1);
z_sh_3 = sort(z_sh_2(:));
nsh = numel(z_sh_3);
specialbin = 0.05 * nsh;
peak_thresh_95 = z_sh_3(specialbin);
fprintf('\n\nPeak thresh from dist: %.4f\n',peak_thresh_95);



%% Fig S1D
rng(1);
cxc = copper(5);

zmean = nanmean(z_peak);
zsem  = nanstd(z_peak) / sqrt(numel(z_peak) - 1);
rand_iter = randi(1000,1,1);
zsh_iter = z_sh_peak{rand_iter};
zsh_mean = nanmean(zsh_iter);
zsh_sem = nanstd(zsh_iter) / sqrt(numel(zsh_iter) - 1);

mse_mean = mean(nanmean(mse_peak,2));
mse_sem = std(nanmean(mse_peak,2)) / sqrt(size(mse_peak,1) - 1);

mse_iter = mse_peak_shuffle{rand_iter};
mse_mean_sh = mean(nanmean(mse_iter,2));
mse_sem_sh = std(nanmean(mse_iter,2)) / sqrt(size(mse_iter,1) - 1);

[~,p_t_z,~,tstats_z] = ttest2(z_peak,zsh_iter)
[~,p_t_m,~,tstats_m] = ttest2(nanmean(mse_peak,2),nanmean(mse_iter,2))


cset = [cxc(4,:); [.6 .6 .6]];
msedat = [nanmean(mse_peak,2),nanmean(mse_iter,2)];
msegroups = [ones(size(msedat,1),1);2*ones(size(msedat,1),1)];
boxplot(msedat,msegroups); hold on;

% set(gca,'ylim',[0 .2])
set(gca,'yscale','log')

%% Fig S1 C

% change these to switch things around
my_unit = 56;
rand_unit = 49;

ms_per_sample = (1/25) / 3;

wffig = figure(); 
set(gcf,'unit','normalized','position',[.1 .1 .7 .6],'color','w');
ax1 = axes(wffig,'position',[.12 .12 .38 .8]);
ax2 = axes(wffig,'position',[.59 .12 .38 .8]);
cx = colormap(ax1,copper(5));
gx = bone(6);
colormap(ax2,gx(1:5,:));


for cc = 1:5
    axes(ax1); hold on;
    wf = meanWF{my_unit}(cc,:);
        wf = wf./abs(min(wf));
%     offset = 5*cc;
    plot(wf,'color',cx(cc,:),'linewidth',3);
    
    axes(ax2); hold on;
    wf = meanWF_shuffle{rand_unit,rand_iter}(cc,:);
    wf = wf./abs(min(wf));
%     offset = 5*cc;
    plot(wf,'color',gx(cc,:),'linewidth',3)
end


axes(ax1);
set(gca,'linewidth',2,'xtick',0:30:120,'fontsize',22,'xlim',[0 90]);
c1 = colorbar(ax1,'ticks',0.1:.2:1,'ticklabels',1:5,...
    'location','westoutside');
ylabel('Experiment day','fontsize',26);
axes(ax2);
c2 = colorbar(ax2,'ticks',0.1:.2:1,'ticklabels',[],'location','westoutside');
set(gca,'linewidth',2,'xtick',0:30:120,'fontsize',22,'xlim',[0 90]);


%% FIG S1B
my_unit = 57; % set to 57
my_unit2 = 56; % set to 56

noscalefig = figure();
set(noscalefig,'unit','normalized','position',[.1 .1 .7 .8],'color','w');
ax1 = axes(noscalefig,'position',[.05 .12 .38 .8]);
ax2 = axes(noscalefig,'position',[.59 .12 .38 .8]);


cx2 = colormap(ax1,copper(5));
gx2 = (pink(7));
colormap(ax2,gx2(1:5,:));
% colormap(ax2,gx2(end-4:end,:));

for cc = 1:5
    axes(ax1); hold on;
    wf = meanWF{my_unit}(cc,:);
%         wf = wf./abs(min(wf));
    offset = 11;
    plot(wf-cc*offset,'color',cx2(cc,:),'linewidth',3)
    
    axes(ax2); hold on;
    wf = meanWF{my_unit2}(cc,:);
%     wf = wf./abs(min(wf));
    offset = 15;
    plot(wf-cc*offset,'color',gx2(cc,:),'linewidth',3)
end

axes(ax1);
c1s = colorbar(ax1,'ticks',0.1:.2:1,'ticklabels',1:5,...
    'location','westoutside');
set(gca,'linewidth',2,'xtick',0:30:120,'fontsize',22,...
    'ylim',[-200,50],'ytick',[-150:50:50],'xlim',[0 97]);
ylabel('Experiment day','fontsize',26);
axes(ax2);
c2s = colorbar(ax2,'ticks',0.1:.2:1,'ticklabels',1:5,'location','westoutside');
set(gca,'linewidth',2,'xtick',0:30:120,'fontsize',22,...
    'ylim',[-200,50],'ytick',[-150:50:50],'xlim',[0 97]);
ylabel('Experiment day','fontsize',26);






