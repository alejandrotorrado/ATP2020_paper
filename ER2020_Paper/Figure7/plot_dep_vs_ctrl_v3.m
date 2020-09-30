close all
clearvars -except CONTCELL*
clc
% version 27 for NRN, version 7 for RNR
VERSION_NUMBER = '7';
do_rnr = 1;
%%
% plot triplet data delta FR
if ispc
    load_dir ='Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig7';
elseif ismac
    load_dir = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig7';
%     load_dir = '/Users/atorrado/Desktop/MLS_DATA/lfp_analysis_beta/Analysis_Data';
end
if do_rnr
    ctrl_data = load([load_dir filesep 'triplet_data_rnr_v' VERSION_NUMBER '_CTRL.mat']);
else
    ctrl_data = load([load_dir filesep 'triplet_data_v' VERSION_NUMBER '_CTRL.mat']);
end
triplet_CTRL = ctrl_data.triplet_data;

if isfield(triplet_CTRL,'spindle_incidence')
    ctrl_do_spindles = 1;
else
    ctrl_do_spindles = 0;
end
if isfield(triplet_CTRL,'SWA_amplitude')
    ctrl_do_swa = 1;
else
    ctrl_do_swa = 0;
end
if isfield(triplet_CTRL,'NR1_thirds_allFR')
    ctrl_do_thirds = 1;
else
    ctrl_do_thirds = 0;
end
if isfield(triplet_CTRL,'allREM_theta')
    ctrl_do_power = 1;
else
    ctrl_do_power = 0;
end
if isfield(triplet_CTRL,'NR1_packets_allFR')
    ctrl_do_packets = 1;
else
    ctrl_do_packets = 0;
end

delta_FR = triplet_CTRL.delta_FR;

allREM_dur = triplet_CTRL.allREM_dur;


if ctrl_do_power
    allREM_theta = triplet_CTRL.allREM_theta;
    allNR1_delta = triplet_CTRL.allNR1_delta;
    allTD_ratio = triplet_CTRL.allTD_ratio;
end
if ctrl_do_spindles
    spindle_incidence = triplet_CTRL.spindle_incidence;
end
if ctrl_do_swa
    SWA_amplitude = triplet_CTRL.SWA_amplitude;
    mean_SWA = triplet_CTRL.mean_SWA_ampl;
end
if ctrl_do_thirds
    NR1_third_FR = triplet_CTRL.NR1_thirds_allFR;
    NR2_third_FR = triplet_CTRL.NR2_thirds_allFR;
    REM_third_FR = triplet_CTRL.REM_thirds_allFR;
end
if ctrl_do_packets
    NR1_packet_FR = triplet_CTRL.NR1_packets_allFR;
    NR2_packet_FR = triplet_CTRL.NR2_packets_allFR;
    
    if isfield(triplet_CTRL,'NR1_p1_allFR')
        have_ps_ctrl = 1;
        NR1_p1 = triplet_CTRL.NR1_p1_allFR;
        NR1_p2 = triplet_CTRL.NR1_p2_allFR;
        NR2_p1 = triplet_CTRL.NR2_p1_allFR;
        NR2_p2 = triplet_CTRL.NR2_p2_allFR;
        REM_p1 = triplet_CTRL.REM_p1_allFR;
        REM_p2 = triplet_CTRL.REM_p2_allFR;
    else
        have_ps_ctrl = 0;
    end
end    
    %% process data
do_avg_by_cell = 1;
% delta_FR_triplets = nanmean(delta_FR,2);

if do_avg_by_cell
    
    % firing rate
    delta_FR_byanim = cellfun(@(x) nanmean(x,2), delta_FR, 'UniformOutput', false);
    delta_FR_nonan_byanim_idx = cellfun(@(x) ~isnan(x), delta_FR_byanim, 'UniformOutput', false);
    for uu = 1:size(delta_FR_byanim,2)
        delta_FR_nonan_byanim{uu} = delta_FR_byanim{uu}(delta_FR_nonan_byanim_idx{uu});
    end
    delta_FR_nonan = cat(1,delta_FR_nonan_byanim{:});
    
    % REM duration
    for uu = 1:size(delta_FR_byanim,2)
        REMdur_nonan_byanim{uu} = allREM_dur{uu}(delta_FR_nonan_byanim_idx{uu});
    end
    REMdur_nonan = cat(1,REMdur_nonan_byanim{:});
    
    if ctrl_do_power
        
        % Theta power
        for uu = 1:size(delta_FR_byanim,2)
            REMtheta_nonan_byanim{uu} = allREM_theta{uu}(delta_FR_nonan_byanim_idx{uu});
        end
        REMtheta_nonan = cat(1,REMtheta_nonan_byanim{:});
        
        % Delta power
        for uu = 1:size(delta_FR_byanim,2)
            NR1delta_nonan_byanim{uu} = allNR1_delta{uu}(delta_FR_nonan_byanim_idx{uu});
        end
        NR1delta_nonan = cat(1,NR1delta_nonan_byanim{:});
        
        % Theta-Delta ration
        for uu = 1:size(delta_FR_byanim,2)
            TDratio_nonan_byanim{uu} = allTD_ratio{uu}(delta_FR_nonan_byanim_idx{uu});
        end
        TDratio_nonan = cat(1,TDratio_nonan_byanim{:});
        
    end
    
    % Spindle incidence
    if ctrl_do_spindles
        for uu = 1:size(delta_FR_byanim,2)
            spindleinc_nonan_byanim{uu} = spindle_incidence{uu}(delta_FR_nonan_byanim_idx{uu});
        end
        spindleinc_nonan = cat(1,spindleinc_nonan_byanim{:});
    end
    
    % SWA
    if ctrl_do_swa
        for uu = 1:size(delta_FR_byanim,2)
            if iscell(SWA_amplitude{uu})
                nancheck = all(cell2mat(cellfun(@(x) isnan(mean(x)),SWA_amplitude{uu},'uniformoutput',0)));
            elseif isnumeric(SWA_amplitude{uu})
                nancheck = isnan(SWA_amplitude{uu});
            end
            if isempty(SWA_amplitude{uu}) || nancheck
                swa_mean{uu} = [];
            else
                swa_mean{uu} = cellfun(@nanmean, SWA_amplitude{uu});
            end
            swa_nonan_byanim{uu} = swa_mean{uu}(delta_FR_nonan_byanim_idx{uu});
            swa_nonan_byanim_norm{uu} = 100 * swa_nonan_byanim{uu} ./ mean_SWA{uu};
        end
        SWA_nonan = cat(1,swa_nonan_byanim_norm{:});
    end
    
    if ctrl_do_thirds
        [NR1_c] = process_sleepthird_data_new2(NR1_third_FR,'thirds');
        [NR2_c] = process_sleepthird_data_new2(NR2_third_FR,'thirds');
        [REM_c] = process_sleepthird_data_new2(REM_third_FR,'thirds');

        
%         [~,idx1_2,idx2_1] = intersect(NR1_c(:,1:3), NR2_c(:,1:3),'rows');
%         [~,idx1_r,idxr_1] = intersect(NR1_c(:,1:3), REM_c(:,1:3),'rows');
%         
%         NR1_c_2 = NR1_c(idx1_2,4:6);
%         NR2_c_1 = NR2_c(idx2_1,4:6);
%         
%         NR1_c_r = NR1_c(idx1_r,4:6);
%         REM_c_1 = REM_c(idxr_1,4:6);
    end
    
    if ctrl_do_packets
%         keyboard;
        NR1_FRpacket_c = process_sleepthird_data_new2(NR1_packet_FR,'packets');
        NR2_FRpacket_c = process_sleepthird_data_new2(NR2_packet_FR,'packets');
        
        if have_ps_ctrl
            NR1_FRp1_c = process_sleepthird_data_new2(NR1_p1,'packets');
            NR1_FRp2_c = process_sleepthird_data_new2(NR1_p2,'packets');
            NR2_FRp1_c = process_sleepthird_data_new2(NR2_p1,'packets');
            NR2_FRp2_c = process_sleepthird_data_new2(NR2_p2,'packets');
            REM_FRp1_c = process_sleepthird_data_new2(REM_p1,'packets');
            REM_FRp2_c = process_sleepthird_data_new2(REM_p2,'packets');
        end
    end
    
else
    
    delta_FR_byanim = cellfun(@(x) nanmean(x,2), delta_FR, 'UniformOutput', false);
    
end


delta_FR_CTRL = delta_FR_nonan;
REMdur_CTRL = REMdur_nonan;
if ctrl_do_power
    REMtheta_CTRL = REMtheta_nonan;
    NR1delta_CTRL = NR1delta_nonan;
    TDratio_CTRL = TDratio_nonan;
end
if ctrl_do_spindles
    spindle_CTRL = spindleinc_nonan;
end
if ctrl_do_swa
    swa_CTRL = SWA_nonan;
end

clearvars -except *_CTRL ctrl_do_* VERSION_NUMBER NR*_c* REM_c* load_dir NR*_*FRpacket* *CONTCELL* do_rnr REM_FRp* NR*_FRp* have*

%%
if do_rnr
    depr_data = load([load_dir filesep 'triplet_data_rnr_v' VERSION_NUMBER '_DEP.mat']);
else
    depr_data = load([load_dir filesep 'triplet_data_v' VERSION_NUMBER '_DEP.mat']);
end

triplet_DEP = depr_data.triplet_data;


if isfield(triplet_DEP,'spindle_incidence')
    dep_do_spindles = 1;
else
    dep_do_spindles = 0;
end
if isfield(triplet_DEP,'SWA_amplitude')
    dep_do_swa = 1;
else
    dep_do_swa = 0;
end
if isfield(triplet_DEP,'NR1_thirds_allFR')
    dep_do_thirds = 1;
else
    dep_do_thirds = 0;
end
if isfield(triplet_DEP,'allREM_theta')
    dep_do_power = 1;
else
    dep_do_power = 0;
end
if isfield(triplet_DEP,'NR1_packets_allFR')
    dep_do_packets = 1;
else
    dep_do_packets = 0;
end

delta_FR = triplet_DEP.delta_FR;
allREM_dur = triplet_DEP.allREM_dur;

if dep_do_power
    allREM_theta = triplet_DEP.allREM_theta;
    allNR1_delta = triplet_DEP.allNR1_delta;
    allTD_ratio = triplet_DEP.allTD_ratio;
end
if dep_do_spindles
    spindle_incidence = triplet_DEP.spindle_incidence;
end
if dep_do_swa
    SWA_amplitude = triplet_DEP.SWA_amplitude;
    mean_SWA = triplet_DEP.mean_SWA_ampl;
end
if dep_do_thirds
    NR1_third_FR = triplet_DEP.NR1_thirds_allFR;
    NR2_third_FR = triplet_DEP.NR2_thirds_allFR;
    REM_third_FR = triplet_DEP.REM_thirds_allFR;

    
end
if dep_do_packets
    NR1_packet_FR = triplet_DEP.NR1_packets_allFR;
    NR2_packet_FR = triplet_DEP.NR2_packets_allFR;
    
    if isfield(triplet_DEP,'NR1_p1_allFR')
        NR1_p1 = triplet_DEP.NR1_p1_allFR;
        NR1_p2 = triplet_DEP.NR1_p2_allFR;
        NR2_p1 = triplet_DEP.NR2_p1_allFR;
        NR2_p2 = triplet_DEP.NR2_p2_allFR;
        REM_p1 = triplet_DEP.REM_p1_allFR;
        REM_p2 = triplet_DEP.REM_p2_allFR;
        have_ps_dep = 1;
    else
        have_ps_dep = 0;
    end
end

%% process data
do_avg_by_cell = 1;
% delta_FR_triplets = nanmean(delta_FR,2);

if do_avg_by_cell
    
    % firing rate
    delta_FR_byanim = cellfun(@(x) nanmean(x,2), delta_FR, 'UniformOutput', false);
    delta_FR_nonan_byanim_idx = cellfun(@(x) ~isnan(x), delta_FR_byanim, 'UniformOutput', false);
    for uu = 1:size(delta_FR_byanim,2)
        delta_FR_nonan_byanim{uu} = delta_FR_byanim{uu}(delta_FR_nonan_byanim_idx{uu});
    end
    delta_FR_nonan = cat(1,delta_FR_nonan_byanim{:});
    
    % REM duration
    for uu = 1:size(delta_FR_byanim,2)
        REMdur_nonan_byanim{uu} = allREM_dur{uu}(delta_FR_nonan_byanim_idx{uu});
    end
    REMdur_nonan = cat(1,REMdur_nonan_byanim{:});
    
    if dep_do_power
        
        % Theta power
        for uu = 1:size(delta_FR_byanim,2)
            REMtheta_nonan_byanim{uu} = allREM_theta{uu}(delta_FR_nonan_byanim_idx{uu});
        end
        REMtheta_nonan = cat(1,REMtheta_nonan_byanim{:});
        
        % Delta power
        for uu = 1:size(delta_FR_byanim,2)
            NR1delta_nonan_byanim{uu} = allNR1_delta{uu}(delta_FR_nonan_byanim_idx{uu});
        end
        NR1delta_nonan = cat(1,NR1delta_nonan_byanim{:});
        
        % Theta-Delta ratio
        for uu = 1:size(delta_FR_byanim,2)
            TDratio_nonan_byanim{uu} = allTD_ratio{uu}(delta_FR_nonan_byanim_idx{uu});
        end
        TDratio_nonan = cat(1,TDratio_nonan_byanim{:});
        
    end
    
    % spindles
    if dep_do_spindles
        for uu = 1:size(delta_FR_byanim,2)
            spindleinc_nonan_byanim{uu} = spindle_incidence{uu}(delta_FR_nonan_byanim_idx{uu});
        end
        spindleinc_nonan = cat(1,spindleinc_nonan_byanim{:});
    end
    
    % SWA
    if dep_do_swa
        for uu = 1:size(delta_FR_byanim,2)
            if iscell(SWA_amplitude{uu})
                nancheck = all(cell2mat(cellfun(@(x) isnan(mean(x)),SWA_amplitude{uu},'uniformoutput',0)));
            elseif isnumeric(SWA_amplitude{uu})
                nancheck = isnan(SWA_amplitude{uu});
            end
            if isempty(SWA_amplitude{uu}) || nancheck
                swa_mean{uu} = [];
            else
                swa_mean{uu} = cellfun(@nanmean, SWA_amplitude{uu});
            end
            swa_nonan_byanim{uu} = swa_mean{uu}(delta_FR_nonan_byanim_idx{uu});
            swa_nonan_byanim_norm{uu} = 100 * swa_nonan_byanim{uu} ./ mean_SWA{uu};
        end
        SWA_nonan = cat(1,swa_nonan_byanim_norm{:});
    end
    
    if dep_do_thirds
        [NR1_d] = process_sleepthird_data_new2(NR1_third_FR,'thirds');
        [NR2_d] = process_sleepthird_data_new2(NR2_third_FR,'thirds');
        [REM_d] = process_sleepthird_data_new2(REM_third_FR,'thirds');        
    end
    
    if dep_do_packets
        NR1_FRpacket_d = process_sleepthird_data_new2(NR1_packet_FR,'packets');
        NR2_FRpacket_d = process_sleepthird_data_new2(NR2_packet_FR,'packets');
        
        if have_ps_dep
            NR1_FRp1_d = process_sleepthird_data_new2(NR1_p1,'packets');
            NR1_FRp2_d = process_sleepthird_data_new2(NR1_p2,'packets');
            NR2_FRp1_d = process_sleepthird_data_new2(NR2_p1,'packets');
            NR2_FRp2_d = process_sleepthird_data_new2(NR2_p2,'packets');
            REM_FRp1_d = process_sleepthird_data_new2(REM_p1,'packets');
            REM_FRp2_d = process_sleepthird_data_new2(REM_p2,'packets');
        end
    end

    
else
    
    delta_FR_byanim = cellfun(@(x) nanmean(x,2), delta_FR, 'UniformOutput', false);
    
end


delta_FR_DEP = delta_FR_nonan;
REMdur_DEP = REMdur_nonan;
if dep_do_power
    REMtheta_DEP = REMtheta_nonan;
    NR1delta_DEP = NR1delta_nonan;
    TDratio_DEP = TDratio_nonan;
end
if dep_do_spindles
    spindle_DEP = spindleinc_nonan;
end
if dep_do_swa
    swa_DEP = SWA_nonan;
end

clearvars -except *_CTRL *_DEP ctrl_do_* dep_do_* NR*_* REM_* *CONTCELL* do_rnr have_ps*
clc

%% General plotting params
setFigureDefaults;
ylims = [-2 2];

figpos = [.1 .1 .3 .7];

%% -------------------- Plot 1, FR vs REM duration -----------------------
xlims1 = [-50 600];

% First calculate correlation coefficients and linear fit
[rho1_c,p1_c] = corr(REMdur_CTRL, delta_FR_CTRL);
ast1_c = get_asterisks_from_pval(p1_c);
lc1_c = polyfit(REMdur_CTRL, delta_FR_CTRL, 1);
[rho1_d,p1_d] = corr(REMdur_DEP, delta_FR_DEP);
ast1_d = get_asterisks_from_pval(p1_d);
lc1_d = polyfit(REMdur_DEP, delta_FR_DEP, 1);

X1 = linspace(xlims1(1),xlims1(2),50);

f1 = figure();
set(f1,'Position',[.1 .1 .45 .8]);

% CTRL
subplot(2,1,1);
hold on;
plot(REMdur_CTRL, delta_FR_CTRL, 'ko','markerfacecolor','none',...
    'markersize',7,'linewidth',1.5);
YC1 = lc1_c(2) + lc1_c(1).*X1;
plot(X1,YC1,'--k','linewidth',1.5);
text(xlims1(2)*.7,ylims(2)*.7,sprintf('r = %.3f\np = %.3f',rho1_c,p1_c),'fontsize',20);
set(gca,'xlim',xlims1,'ylim',ylims);
ylabel('Firing rate (z)','fontsize',23);
box off

% DEP
subplot(2,1,2);
hold on;
plot(REMdur_DEP, delta_FR_DEP, 'ko','markerfacecolor','k',...
    'markersize',7,'linewidth',1.5);
YD1 = lc1_d(2) + lc1_d(1).*X1;
plot(X1,YD1,'--k','linewidth',1.5);
text(xlims1(2)*.7,ylims(2)*.7,sprintf('r = %.3f\np = %.3f',rho1_d,p1_d),'fontsize',20);
set(gca,'xlim',xlims1,'ylim',ylims);
xlabel('REM duration (s)','fontsize',22);
ylabel('Firing rate (z)','fontsize',23);
box off


%% ------------------ Plot 2, FR vs REM theta power ----------------------

if ctrl_do_power && dep_do_power
    xlims2 = [-0.7 0.7];
    
    % First calculate correlation coefficients and linear fit
    [rho2_c,p2_c] = corr(REMtheta_CTRL, delta_FR_CTRL);
    ast2_c = get_asterisks_from_pval(p2_c);
    lc2_c = polyfit(REMtheta_CTRL, delta_FR_CTRL, 1);
    [rho2_d,p2_d] = corr(REMtheta_DEP, delta_FR_DEP);
    ast2_d = get_asterisks_from_pval(p2_d);
    lc2_d = polyfit(REMtheta_DEP, delta_FR_DEP, 1);
    
    X2 = linspace(xlims2(1),xlims2(2),50);
    
    f2 = figure();
    set(f2,'Position',[.1 .1 .45 .8]);
    
    % CTRL
    subplot(2,1,1);
    hold on;
    plot(REMtheta_CTRL, delta_FR_CTRL, 'ko','markerfacecolor','none',...
        'markersize',7,'linewidth',1.5);
    Y2 = lc2_c(2) + lc2_c(1).*X2;
    plot(X2,Y2,'--k','linewidth',1.5);
    text(xlims2(2)*.7,ylims(2)*.7,sprintf('r = %.3f\np = %.3f',rho2_c,p2_c),'fontsize',20);
    set(gca,'xlim',xlims2,'ylim',ylims);
    ylabel('Firing rate (z)','fontsize',23);
    box off
    
    % DEP
    subplot(2,1,2);
    hold on;
    plot(REMtheta_DEP, delta_FR_DEP, 'ko','markerfacecolor','k',...
        'markersize',7,'linewidth',1.5);
    YD2 = lc2_d(2) + lc2_d(1).*X2;
    plot(X2,YD2,'--k','linewidth',1.5);
    text(xlims2(2)*.7,ylims(2)*.7,sprintf('r = %.3f\np = %.3f',rho2_d,p2_d),'fontsize',20);
    set(gca,'xlim',xlims2,'ylim',ylims);
    xlabel('REM theta power (z)','fontsize',22);
    ylabel('Firing rate (z)','fontsize',23);
    box off
    
    %% ------------------ Plot 3, FR vs NR1 delta power ----------------------
    xlims3 = [-1.5 1.5];
    
    % First calculate correlation coefficients and linear fit
    [rho3_c,p3_c] = corr(NR1delta_CTRL, delta_FR_CTRL);
    ast3_c = get_asterisks_from_pval(p3_c);
    lc3_c = polyfit(NR1delta_CTRL, delta_FR_CTRL, 1);
    [rho3_d,p3_d] = corr(NR1delta_DEP, delta_FR_DEP);
    ast3_d = get_asterisks_from_pval(p3_d);
    lc3_d = polyfit(NR1delta_DEP, delta_FR_DEP, 1);
    
    X3 = linspace(xlims3(1),xlims3(2),50);
    
    f3 = figure();
    set(f3,'Position',[.1 .1 .45 .8]);
    
    % CTRL
    subplot(2,1,1);
    hold on;
    plot(NR1delta_CTRL, delta_FR_CTRL, 'ko','markerfacecolor','none',...
        'markersize',7,'linewidth',1.5);
    Y3 = lc3_c(2) + lc3_c(1).*X3;
    plot(X3,Y3,'--k','linewidth',1.5);
    text(xlims3(2)*.7,ylims(2)*.7,sprintf('r = %.3f\np = %.3f',rho3_c,p3_c),'fontsize',20);
    set(gca,'xlim',xlims3,'ylim',ylims);
    ylabel('Firing rate (z)','fontsize',23);
    box off
    
    % DEP
    subplot(2,1,2);
    hold on;
    plot(NR1delta_DEP, delta_FR_DEP, 'ko','markerfacecolor','k',...
        'markersize',7,'linewidth',1.5);
    YD3 = lc3_d(2) + lc3_d(1).*X3;
    plot(X3,YD3,'--k','linewidth',1.5);
    text(xlims3(2)*.7,ylims(2)*.7,sprintf('r = %.3f\np = %.3f',rho3_d,p3_d),'fontsize',20);
    set(gca,'xlim',xlims3,'ylim',ylims);
    xlabel('NREM delta power (z)','fontsize',22);
    ylabel('Firing rate (z)','fontsize',23);
    box off
    
end

%% ------------------ Plot 4, FR vs NR1 spindle incidence ----------------------
if ctrl_do_spindles && dep_do_spindles
    
    xlims4 = [-.01 0.2];
    
    % First calculate correlation coefficients and linear fit
    [rho4_c,p4_c] = corr(spindle_CTRL, delta_FR_CTRL);
    ast4_c = get_asterisks_from_pval(p4_c);
    lc4_c = polyfit(spindle_CTRL, delta_FR_CTRL, 1);
    [rho4_d,p4_d] = corr(spindle_DEP, delta_FR_DEP);
    ast4_d = get_asterisks_from_pval(p4_d);
    lc4_d = polyfit(spindle_DEP, delta_FR_DEP, 1);
    
    X4 = linspace(xlims4(1),xlims4(2),50);
    
    f4 = figure();
    set(f4,'Position',[.1 .1 .45 .8]);
    
    % CTRL
    subplot(2,1,1);
    hold on;
    plot(spindle_CTRL, delta_FR_CTRL, 'ko','markerfacecolor','none',...
        'markersize',7,'linewidth',1.5);
    Y4 = lc4_c(2) + lc4_c(1).*X4;
    plot(X4,Y4,'--k','linewidth',1.5);
    text(xlims4(2)*.7,ylims(2)*.7,sprintf('r = %.3f\np = %.3f',rho4_c,p4_c),'fontsize',20);
    set(gca,'xlim',xlims4,'ylim',ylims);
    ylabel('Firing rate (z)','fontsize',23);
    box off
    
    % DEP
    subplot(2,1,2);
    hold on;
    plot(spindle_DEP, delta_FR_DEP, 'ko','markerfacecolor','k',...
        'markersize',7,'linewidth',1.5);
    YD4 = lc4_d(2) + lc4_d(1).*X4;
    plot(X4,YD4,'--k','linewidth',1.5);
    text(xlims4(2)*.7,ylims(2)*.7,sprintf('r = %.3f\np = %.3f',rho4_d,p4_d),'fontsize',20);
    set(gca,'xlim',xlims4,'ylim',ylims);
    xlabel('Spindle incidence (1/s)','fontsize',22);
    ylabel('Firing rate (z)','fontsize',23);
    box off
    
end


%% ----------------------- Plot 5, FR vs NR1 SWA --------------------------
if ctrl_do_swa && dep_do_swa
    
    xlims5 = [-50 350];
    
    % First calculate correlation coefficients and linear fit
    [rho5_c,p5_c] = corr(swa_CTRL, delta_FR_CTRL);
    ast5_c = get_asterisks_from_pval(p5_c);
    lc5_c = polyfit(swa_CTRL, delta_FR_CTRL, 1);
    [rho5_d,p5_d] = corr(swa_DEP, delta_FR_DEP);
    ast5_d = get_asterisks_from_pval(p5_d);
    lc5_d = polyfit(swa_DEP, delta_FR_DEP, 1);
    
    X5 = linspace(xlims5(1),xlims5(2),50);
    
    f5 = figure();
    set(f5,'Position',[.1 .1 .45 .8]);
    
    % CTRL
    subplot(2,1,1);
    hold on;
    plot(swa_CTRL, delta_FR_CTRL, 'ko','markerfacecolor','none',...
        'markersize',7,'linewidth',1.5);
    Y5 = lc5_c(2) + lc5_c(1).*X5;
    plot(X5,Y5,'--k','linewidth',1.5);
    text(xlims5(2)*.7,ylims(2)*.7,sprintf('r = %.3f\np = %.3f',rho5_c,p5_c),'fontsize',20);
    set(gca,'xlim',xlims5,'ylim',ylims);
    ylabel('Firing rate (z)','fontsize',23);
    box off
    
    % DEP
    subplot(2,1,2);
    hold on;
    plot(swa_DEP, delta_FR_DEP, 'ko','markerfacecolor','k',...
        'markersize',7,'linewidth',1.5);
    YD5 = lc5_d(2) + lc5_d(1).*X5;
    plot(X5,YD5,'--k','linewidth',1.5);
    text(xlims5(2)*.7,ylims(2)*.7,sprintf('r = %.3f\np = %.3f',rho5_d,p5_d),'fontsize',20);
    set(gca,'xlim',xlims5,'ylim',ylims);
    xlabel('SWA (% of 36-hour mean)','fontsize',22);
    ylabel('Firing rate (z)','fontsize',23);
    box off
    
end

%% thirds
if ctrl_do_thirds && dep_do_thirds
    % z-scoring parameter
    mean_t_CTRL = triplet_CTRL.params.mean_t;
    mean_t_DEP = triplet_DEP.params.mean_t;
    
    
    c_rem   = [25 181 149]./255;
    c_nrem  = [131 49 146]./255;
    if do_rnr
        c_nrem = [25 181 149]./255;
        c_rem = [131 49 146]./255;
    end
    c_aw    = [201 28 101]./255;
    c_qw    = [247 148 41]./255;
    %% data processing
       
    
    if mean_t_CTRL == 1 && mean_t_DEP == 1
        
        mean_t = 1;
        
        third_norm = 0;
        
        % compile the data
        NR1_d_plot = NR1_d;
        NR2_d_plot = NR2_d;
        REM_d_plot = REM_d;
        
        NR1_c_plot = NR1_c;
        NR2_c_plot = NR2_c;
        REM_c_plot = REM_c;
        
        % change in FR across states in triplet. For z-scoring to triplet,
        % just do the difference in z-scores
        NR1_change_d = NR1_d_plot(:,3) - NR1_d_plot(:,1);
        REM_change_d = REM_d_plot(:,3) - REM_d_plot(:,1);
        NR2_change_d = NR2_d_plot(:,3) - NR2_d_plot(:,1);
        
        NR1_change_c = NR1_c_plot(:,3) - NR1_c_plot(:,1);
        REM_change_c = REM_c_plot(:,3) - REM_c_plot(:,1);
        NR2_change_c = NR2_c_plot(:,3) - NR2_c_plot(:,1);
        
        % average FR by state in triplet
        NR1_mean_cell_d = nanmean(NR1_d_plot,2);
        NR2_mean_cell_d = nanmean(NR2_d_plot,2);
        % change between NR2 and NR1
        NR12_change_d = (NR2_mean_cell_d) - (NR1_mean_cell_d);
        % change between NR2 and NR1 - middle third only
        NR12_middle_change_d = NR2_d_plot(:,2) - NR1_d_plot(:,2);
                
        % same as above, for control hemisphere
        NR1_mean_cell_c = nanmean(NR1_c_plot,2);
        NR2_mean_cell_c = nanmean(NR2_c_plot,2);
        NR12_change_c = (NR2_mean_cell_c) - (NR1_mean_cell_c);
        NR12_middle_change_c = NR2_c_plot(:,2) - NR1_c_plot(:,2);
        
    else
        mean_t = 10;

        NR1_d_plot = NR1_d;
        NR2_d_plot = NR2_d;
        REM_d_plot = REM_d;
        
        NR1_c_plot = NR1_c;
        NR2_c_plot = NR2_c;
        REM_c_plot = REM_c;
        %     end
        
        logbins = logspace(-2,2,50);
        
        NR1_T1_d = histcounts(NR1_d_plot(:,1),logbins);
        NR1_T3_d = histcounts(NR1_d_plot(:,3),logbins);
        
        NR2_T1_d = histcounts(NR2_d_plot(:,1),logbins);
        NR2_T3_d = histcounts(NR2_d_plot(:,3),logbins);
        
        REM_T1_d = histcounts(REM_d_plot(:,1),logbins);
        REM_T3_d = histcounts(REM_d_plot(:,3),logbins);
        
        NR1_T1_c = histcounts(NR1_c_plot(:,1),logbins);
        NR1_T3_c = histcounts(NR1_c_plot(:,3),logbins);
        
        NR2_T1_c = histcounts(NR2_c_plot(:,1),logbins);
        NR2_T3_c = histcounts(NR2_c_plot(:,3),logbins);
        
        REM_T1_c = histcounts(REM_c_plot(:,1),logbins);
        REM_T3_c = histcounts(REM_c_plot(:,3),logbins);
        
        NR1_change_d = 100 .* ((NR1_d_plot(:,3) - NR1_d_plot(:,1)) ./ NR1_d_plot(:,1));
        REM_change_d = 100 .* ((REM_d_plot(:,3) - REM_d_plot(:,1)) ./ REM_d_plot(:,1));
        NR2_change_d = 100 .* ((NR2_d_plot(:,3) - NR2_d_plot(:,1)) ./ NR2_d_plot(:,1));
        
        NR1_change_c = 100 .* ((NR1_c_plot(:,3) - NR1_c_plot(:,1)) ./ NR1_c_plot(:,1));
        REM_change_c = 100 .* ((REM_c_plot(:,3) - REM_c_plot(:,1)) ./ REM_c_plot(:,1));
        NR2_change_c = 100 .* ((NR2_c_plot(:,3) - NR2_c_plot(:,1)) ./ NR2_c_plot(:,1));
        
        NR1_mean_cell_d = nanmean(NR1_d_plot,2);
        NR2_mean_cell_d = nanmean(NR2_d_plot,2);
%         NR12_change_d = 100 .* ((NR2_mean_cell_d - NR1_mean_cell_d) ./ NR1_mean_cell_d);
            NR12_change_d = (NR2_mean_cell_d) - (NR1_mean_cell_d);
        
        NR12_middle_change_d = 100 .* ((NR2_d_plot(:,2) - NR1_d_plot(:,2)) ./ NR1_d_plot(:,2));
        NR12_middle_change_c = 100 .* ((NR2_c_plot(:,2) - NR1_c_plot(:,2)) ./ NR1_c_plot(:,2));
        
        NR1_mean_cell_c = nanmean(NR1_c_plot,2);
        NR2_mean_cell_c = nanmean(NR2_c_plot,2);
%         NR12_change_c = 100 .* ((NR2_mean_cell_c - NR1_mean_cell_c) ./ NR1_mean_cell_c);
            NR12_change_c = (NR2_mean_cell_c) - (NR1_mean_cell_c);
    end
    
%     dep_NR1_means = [nanmean(NR1_T1_d), nanmean(NR1_T2_d), nanmean(NR1_T3_d)];
%     dep_REM_means = [nanmean(REM_T1_d), nanmean(REM_T2_d), nanmean(REM_T3_d)];
%     dep_NR2_means = [nanmean(NR2_T1_d), nanmean(NR2_T2_d), nanmean(NR2_T3_d)];
%     
%     dep_NR1_sem = [std(NR1_T1_d)/sqrt(numel(NR1_T1_d)-1), std(NR1_T2_d)/sqrt(numel(NR1_T2_d)-1),...
%         std(NR1_T3_d)/sqrt(numel(NR1_T3_d)-1)];
%     dep_REM_sem = [std(REM_T1_d)/sqrt(numel(REM_T1_d)-1), std(REM_T2_d)/sqrt(numel(REM_T2_d)-1),...
%         std(REM_T3_d)/sqrt(numel(REM_T3_d)-1)];
%     dep_NR2_sem = [std(NR2_T1_d)/sqrt(numel(NR2_T1_d)-1), std(NR2_T2_d)/sqrt(numel(NR2_T2_d)-1),...
%         std(NR2_T3_d)/sqrt(numel(NR2_T3_d)-1)];
%     
%     ctrl_NR1_means = [nanmean(NR1_T1_c), nanmean(NR1_T2_c), nanmean(NR1_T3_c)];
%     ctrl_REM_means = [nanmean(REM_T1_c), nanmean(REM_T2_c), nanmean(REM_T3_c)];
%     ctrl_NR2_means = [nanmean(NR2_T1_c), nanmean(NR2_T2_c), nanmean(NR2_T3_c)];
%     
%     ctrl_NR1_sem = [std(NR1_T1_c)/sqrt(numel(NR1_T1_c)-1), std(NR1_T2_c)/sqrt(numel(NR1_T2_c)-1),...
%         std(NR1_T3_c)/sqrt(numel(NR1_T3_c)-1)];
%     ctrl_REM_sem = [std(REM_T1_c)/sqrt(numel(REM_T1_c)-1), std(REM_T2_c)/sqrt(numel(REM_T2_c)-1),...
%         std(REM_T3_c)/sqrt(numel(REM_T3_c)-1)];
%     ctrl_NR2_sem = [std(NR2_T1_c)/sqrt(numel(NR2_T1_c)-1), std(NR2_T2_c)/sqrt(numel(NR2_T2_c)-1),...
%         std(NR2_T3_c)/sqrt(numel(NR2_T3_c)-1)];

    dep_NR1_means = nanmean(NR1_d_plot);
    dep_NR2_means = nanmean(NR2_d_plot);
    dep_REM_means = nanmean(REM_d_plot);
    
    dep_NR1_sem = std(NR1_d_plot,0,'omitnan') ./ sqrt(sum(~isnan(NR1_d_plot(:,1)))-1);
    dep_NR2_sem = std(NR2_d_plot,0,'omitnan') ./ sqrt(sum(~isnan(NR2_d_plot(:,1)))-1);
    dep_REM_sem = std(REM_d_plot,0,'omitnan') ./ sqrt(sum(~isnan(REM_d_plot(:,1)))-1);
    
    ctrl_NR1_means = nanmean(NR1_c_plot);
    ctrl_NR2_means = nanmean(NR2_c_plot);
    ctrl_REM_means = nanmean(REM_c_plot);
    
    ctrl_NR1_sem = std(NR1_c_plot,0,'omitnan') ./ sqrt(sum(~isnan(NR1_c_plot(:,1)))-1);
    ctrl_NR2_sem = std(NR2_c_plot,0,'omitnan') ./ sqrt(sum(~isnan(NR2_c_plot(:,1)))-1);
    ctrl_REM_sem = std(REM_c_plot,0,'omitnan') ./ sqrt(sum(~isnan(REM_c_plot(:,1)))-1);

    
    
    
    %% plotting
    setFigureDefaults;
    tfig = figure();
    set(tfig,'position',[.05 .1 .6 .8]);
    x_axis = [1,4,7,8,9,10,11,14,17];
    
    subplot(2,1,1);
    errorbar(x_axis,[ctrl_NR1_means,ctrl_REM_means,ctrl_NR2_means],...
        [ctrl_NR1_sem,ctrl_REM_sem,ctrl_NR2_sem],'k','linewidth',2);
    hold on
    errorbar(x_axis,[dep_NR1_means,dep_REM_means,dep_NR2_means],...
        [dep_NR1_sem,dep_REM_sem,dep_NR2_sem],'b','linewidth',2);
    box off
    

    
    
    %{
    lw = 2;
    hfig = figure();
    set(hfig,'position',[.05 .1 .9 .85]);
    subplot(2,3,1);
    hold on;
    plot(NR1_T1_d./sum(NR1_T1_d),'b','linewidth',lw);
    plot(NR1_T3_d./sum(NR1_T3_d),'r','linewidth',lw);
    set(gca,'xtick',[1,13,25.5,37.8,50],'xticklabel',[0.01,0.1,1.0,10,100],...
        'xlim',[1,50],'ylim',[0,0.1]);
    
    subplot(2,3,4);
    hold on;
    plot(NR1_T1_c./sum(NR1_T1_c),'k','linewidth',lw);
    plot(NR1_T3_c./sum(NR1_T3_c),'color',[.6 .6 .6],'linewidth',lw);
    set(gca,'xtick',[1,13,25.5,37.8,50],'xticklabel',[0.01,0.1,1.0,10,100],...
        'xlim',[1,50],'ylim',[0,0.1]);
    
    subplot(2,3,2);
    hold on;
    plot(NR2_T1_d./sum(NR2_T1_d),'b','linewidth',lw);
    plot(NR2_T3_d./sum(NR2_T3_d),'r','linewidth',lw);
    set(gca,'xtick',[1,13,25.5,37.8,50],'xticklabel',[0.01,0.1,1.0,10,100],...
        'xlim',[1,50],'ylim',[0,0.1]);
    subplot(2,3,5);
    hold on;
    plot(NR2_T1_c./sum(NR2_T1_c),'k','linewidth',lw);
    plot(NR2_T3_c./sum(NR2_T3_c),'color',[.6 .6 .6],'linewidth',lw);
    set(gca,'xtick',[1,13,25.5,37.8,50],'xticklabel',[0.01,0.1,1.0,10,100],...
        'xlim',[1,50],'ylim',[0,0.1]);
    
    subplot(2,3,3);
    hold on;
    plot(REM_T1_d./sum(REM_T1_d),'b','linewidth',lw);
    plot(REM_T3_d./sum(REM_T3_d),'r','linewidth',lw);
    set(gca,'xtick',[1,13,25.5,37.8,50],'xticklabel',[0.01,0.1,1.0,10,100],...
        'xlim',[1,50],'ylim',[0,0.1]);
    subplot(2,3,6);
    hold on;
    plot(REM_T1_c./sum(REM_T1_c),'k','linewidth',lw);
    plot(REM_T3_c./sum(REM_T3_c),'color',[.6 .6 .6],'linewidth',lw);
    set(gca,'xtick',[1,13,25.5,37.8,50],'xticklabel',[0.01,0.1,1.0,10,100],...
        'xlim',[1,50],'ylim',[0,0.1]);
    %}
    
    %%
    if mean_t == 1
        ylims = [-1 1];
    else
        ylims = [-40 80];
    end
    
    
    lw = 3;
    percfig = figure('position',[.1 .1 .7 .6]);
    d_data = padcat(NR1_change_d,REM_change_d,NR2_change_d);
    c_data = padcat(NR1_change_c,REM_change_c,NR2_change_c);
    cset = [c_nrem; c_rem; c_nrem];
    
    mean_data_d = nanmean(d_data);
    sem_data_d = nanstd(d_data) ./ sqrt(sum(~isnan(d_data(:,1)))-1);
    
    mean_data_c = nanmean(c_data);
    sem_data_c = nanstd(c_data) ./ sqrt(sum(~isnan(c_data(:,1)))-1);
    
    subplot(1,2,2); box off
    UnivarScatter(d_data,'Label',{'NR1','REM', 'NR2'},'BoxType','SEM',...
    'StdColor',[1 1 1],'MeanColor',[1 1 1],'SEMColor',[1 1 1],'Width',.9,'Compression',30,...
    'MarkerEdgeColor',cset,'MarkerFaceColor',cset,'PointSize',30,...
    'DataTransform','None');
    box off
    hold on
%     errorbar(1:3,mean_data_d,sem_data_d,'k','linestyle','none','capsize',15,'linewidth',lw);
    errorbar(1:3,mean_data_d,sem_data_d,'k','linestyle','none','linewidth',lw);
    for xx = 1:3
        line([xx-.25 xx+.25],[mean_data_d(xx) mean_data_d(xx)],'linewidth',lw,...
            'color','k');
    end
    line([0 4],[0 0],'linestyle','--','color',[.7 .7 .7],'linewidth',2);
    set(gca,'ylim',ylims);
    
    [p_kw_d,~,st_kw_d] = kruskalwallis(d_data,[],'off');
    c_d = multcompare(st_kw_d,'display','off');
    

    subplot(1,2,1);
    
    UnivarScatter(c_data,'Label',{'NR1','REM', 'NR2'},'BoxType','SEM',...
    'StdColor',[1 1 1],'SEMColor',[1 1 1],'MeanColor',[1 1 1],'Width',.9,'Compression',30,...
    'MarkerEdgeColor',cset,'MarkerFaceColor','none','PointSize',30,...
    'DataTransform','None');
    box off
    hold on
%     errorbar(1:3,mean_data_c,sem_data_c,'k','linestyle','none','capsize',15,'linewidth',lw);
    errorbar(1:3,mean_data_c,sem_data_c,'k','linestyle','none','linewidth',lw);
    for xx = 1:3
        line([xx-.25 xx+.25],[mean_data_c(xx) mean_data_c(xx)],'linewidth',lw,...
            'color','k');
    end
    line([0 4],[0 0],'linestyle','--','color',[.7 .7 .7],'linewidth',2);
    set(gca,'ylim',ylims);
    ylabel('% change in FR across state');
    
    [p_kw_c,~,st_kw_c] = kruskalwallis(c_data,[],'off');
    c_c = multcompare(st_kw_c,'display','off');
    
    %% paired NR thirds
    if mean_t == 1
        scale = 'linear';
    else
        scale = 'log';
    end
    
    pairfig = figure();
    set(pairfig,'position',[.04 .1 .95 .6]);
    
    subplot(1,3,1);
    pair_c = [NR1_mean_cell_c, NR2_mean_cell_c];
    plot(pair_c','k-o','markerfacecolor','k');
    box off
    set(gca,'xlim',[0.5 2.5],'yscale',scale,'xtick',[1 2],'xticklabel',{'NR1','NR2'});
    ylabel('Firing rate');
    
    subplot(1,3,2);
    pair_d = [NR1_mean_cell_d, NR2_mean_cell_d];
    plot(pair_d','b-o','markerfacecolor','b');
    box off
    set(gca,'xlim',[0.5 2.5],'yscale',scale,'xtick',[1 2],'xticklabel',{'NR1','NR2'});
    
    cset1 = [0.6 0.6 0.6; 0 0 1];
    subplot(1,3,3);
    pair_both = padcat(NR12_change_c, NR12_change_d);
    mean_pair_both = nanmean(pair_both);
    sem_pair_both = nanstd(pair_both) ./ sqrt(sum(~isnan(pair_both))-1);
    
    UnivarScatter(pair_both,'Label',{'CTRL','ER'},'BoxType','SEM',...
        'StdColor',[1 1 1],'SEMColor',[1 1 1],'MeanColor',[1 1 1],'Width',.9,'Compression',30,...
        'MarkerEdgeColor',cset1,'MarkerFaceColor',cset1,'PointSize',30,...
        'DataTransform','None');
    box off
    hold on
    errorbar(1:2,mean_pair_both,sem_pair_both,'k','linestyle','none','linewidth',lw);
%     errorbar(1:2,mean_pair_both,sem_pair_both,'k','linestyle','none','capsize',15,'linewidth',lw);
    for xx = 1:2
        line([xx-.25 xx+.25],[mean_pair_both(xx) mean_pair_both(xx)],'linewidth',lw,...
            'color','k');
    end
    line([0.5 2.5],[0 0],'linestyle','--','color',[.7 .7 .7],'linewidth',2);
    set(gca,'xlim',[.5 2.5],'xtick',[1 2],'xticklabel',{'CTRL','ER'});
    ylabel('\Delta change NR1-NR2');
    
    %% stats
    
    p_nr1c = signrank(NR1_change_c,0)
    p_remc = signrank(REM_change_c,0)
    p_nr2c = signrank(NR2_change_c,0)
    
    p_nr1d = signrank(NR1_change_d,0)
    p_remd = signrank(REM_change_d,0)
    p_nr2d = signrank(NR2_change_d,0)
    
    p_nr12c = signrank(NR1_mean_cell_c, NR2_mean_cell_c)
    p_nr12d = signrank(NR1_mean_cell_d, NR2_mean_cell_d)
    
    p_nr12c0 = signrank(NR12_change_c,0)
    p_nr12d0 = signrank(NR12_change_d,0)
    
    [~,p_nr12c0_t] = ttest(NR12_change_c,0)
    [~,p_nr12d0_t] = ttest(NR12_change_d,0)
    
    
    
    
end


%% packets
if ctrl_do_packets && dep_do_packets
    
    %% process data
    
    NR1_mean_c = nanmean(NR1_FRpacket_c);
    NR2_mean_c = nanmean(NR2_FRpacket_c);
    NR1_sem_c = std(NR1_FRpacket_c,0,'omitnan') ./ sqrt(sum(~isnan(NR1_FRpacket_c))-1);
    NR2_sem_c = std(NR2_FRpacket_c,0,'omitnan') ./ sqrt(sum(~isnan(NR2_FRpacket_c))-1);

    NR1_mean_d = nanmean(NR1_FRpacket_d);
    NR2_mean_d = nanmean(NR2_FRpacket_d);
    NR1_sem_d = std(NR1_FRpacket_d,0,'omitnan') ./ sqrt(sum(~isnan(NR1_FRpacket_d))-1);
    NR2_sem_d = std(NR2_FRpacket_d,0,'omitnan') ./ sqrt(sum(~isnan(NR2_FRpacket_d))-1);
    
    NR12_packchange_c = NR2_FRpacket_c - NR1_FRpacket_c;
    NR12_packchange_d = NR2_FRpacket_d - NR1_FRpacket_d;
    
    % z-scoring parameter
    mean_t_CTRL = triplet_CTRL.params.mean_t;
    mean_t_DEP = triplet_DEP.params.mean_t;
    
    mean_t = unique([mean_t_CTRL, mean_t_DEP]);
    if numel(mean_t) > 1
        disp('mean_t issue');
        keyboard;
    end
    
    if mean_t == 1
        scale = 'linear';
        change_mode = 'delta'; % delta or perc
        ylz = [-0.8 0.8];
    elseif mean_t== 10
        scale = 'log';
        change_mode = 'perc'; % delta or perc
        ylz = [0.01 100];
    end
    
    
    if have_ps_dep && have_ps_ctrl
        switch change_mode
            case 'delta'
                NR1_p_change_c = NR1_FRp2_c - NR1_FRp1_c;
                NR2_p_change_c = NR2_FRp2_c - NR2_FRp1_c;
                REM_p_change_c = REM_FRp2_c - REM_FRp1_c;
                
                NR1_p_change_d = NR1_FRp2_d - NR1_FRp1_d;
                NR2_p_change_d = NR2_FRp2_d - NR2_FRp1_d;
                REM_p_change_d = REM_FRp2_d - REM_FRp1_d;
                
                yls = [-1 1];
            case 'perc'
                NR1_p_change_c = 100 * (NR1_FRp2_c - NR1_FRp1_c) ./ NR1_FRp1_c;
                NR1_p_change_c = 100 * (NR2_FRp2_c - NR2_FRp1_c) ./ NR2_FRp1_c;
                NR1_p_change_c = 100 * (REM_FRp2_c - REM_FRp1_c) ./ REM_FRp1_c;
                
                NR1_p_change_d = 100 * (NR1_FRp2_d - NR1_FRp1_d) ./ NR1_FRp1_d;
                NR1_p_change_d = 100 * (NR2_FRp2_d - NR2_FRp1_d) ./ NR2_FRp1_d;
                NR1_p_change_d = 100 * (REM_FRp2_d - REM_FRp1_d) ./ REM_FRp1_d;
                
                yls = [-500 500];
        end
    end
    
    
    %% plotting
    packfig = figure();
    set(packfig,'position',[.04 .1 .95 .6]);
    
    subplot(1,3,1);
    pair_c = [NR1_FRpacket_c, NR2_FRpacket_c];
    plot(pair_c','k-o','markerfacecolor','k');
    box off
    set(gca,'ylim',ylz,'ytick',-1:0.2:1,'xlim',[0.5 2.5],'yscale',scale,'xtick',[1 2],'xticklabel',{'NR1','NR2'});
    ylabel('Firing rate');
    
    subplot(1,3,2);
    pair_d = [NR1_FRpacket_d, NR2_FRpacket_d];
    plot(pair_d','b-o','markerfacecolor','b');
    box off
    set(gca,'ylim',ylz,'ytick',-1:0.2:1,'xlim',[0.5 2.5],'yscale',scale,'xtick',[1 2],'xticklabel',{'NR1','NR2'});
    
    cset1 = [0.6 0.6 0.6; 0 0 1];
    subplot(1,3,3);
    pair_both_p = padcat(NR12_packchange_c, NR12_packchange_d);
    mean_pair_both_p = nanmean(pair_both_p);
    sem_pair_both_p = nanstd(pair_both_p) ./ sqrt(sum(~isnan(pair_both_p))-1);
    
    UnivarScatter(pair_both_p,'Label',{'CTRL','ER'},'BoxType','SEM',...
        'StdColor',[1 1 1],'SEMColor',[1 1 1],'MeanColor',[1 1 1],'Width',.9,'Compression',30,...
        'MarkerEdgeColor',cset1,'MarkerFaceColor',cset1,'PointSize',30,...
        'DataTransform','None');
    box off
    hold on
    errorbar(1:2,mean_pair_both_p,sem_pair_both_p,'k','linestyle','none','linewidth',lw);
%     errorbar(1:2,mean_pair_both_p,sem_pair_both_p,'k','linestyle','none','capsize',15,'linewidth',lw);
    for xx = 1:2
        line([xx-.25 xx+.25],[mean_pair_both_p(xx) mean_pair_both_p(xx)],'linewidth',lw,...
            'color','k');
    end
    line([0.5 2.5],[0 0],'linestyle','--','color',[.7 .7 .7],'linewidth',2);
    set(gca,'xlim',[.5 2.5],'xtick',[1 2],'xticklabel',{'CTRL','ER'},'ylim',[-0.8 0.4],'ytick',-1:.2:1);
    ylabel('\Delta change NR1-NR2');
    
    %% perc_change_plot
    
    if have_ps_dep && have_ps_ctrl
        
        
        lw = 3;
        percfig = figure('position',[.1 .1 .7 .6]);
        d_data = padcat(NR1_p_change_d,REM_p_change_d,NR2_p_change_d);
        c_data = padcat(NR1_p_change_c,REM_p_change_c,NR2_p_change_c);
        cset = [c_nrem; c_rem; c_nrem];
        
        mean_data_d = nanmean(d_data);
        sem_data_d = nanstd(d_data) ./ sqrt(sum(~isnan(d_data(:,1)))-1);
        
        mean_data_c = nanmean(c_data);
        sem_data_c = nanstd(c_data) ./ sqrt(sum(~isnan(c_data(:,1)))-1);
        
        subplot(1,2,2); box off
        UnivarScatter(d_data,'Label',{'NR1','REM', 'NR2'},'BoxType','SEM',...
            'StdColor',[1 1 1],'MeanColor',[1 1 1],'SEMColor',[1 1 1],'Width',.9,'Compression',30,...
            'MarkerEdgeColor',cset,'MarkerFaceColor',cset,'PointSize',30,...
            'DataTransform','None');
        box off
        hold on
        %     errorbar(1:3,mean_data_d,sem_data_d,'k','linestyle','none','capsize',15,'linewidth',lw);
        errorbar(1:3,mean_data_d,sem_data_d,'k','linestyle','none','linewidth',lw);
        for xx = 1:3
            line([xx-.25 xx+.25],[mean_data_d(xx) mean_data_d(xx)],'linewidth',lw,...
                'color','k');
        end
        line([0 4],[0 0],'linestyle','--','color',[.7 .7 .7],'linewidth',2);
        set(gca,'ylim',yls);
        
        [p_kw_d,~,st_kw_d] = kruskalwallis(d_data,[],'off');
        c_d = multcompare(st_kw_d,'display','off');
        
        
        subplot(1,2,1);
        
        UnivarScatter(c_data,'Label',{'NR1','REM', 'NR2'},'BoxType','SEM',...
            'StdColor',[1 1 1],'SEMColor',[1 1 1],'MeanColor',[1 1 1],'Width',.9,'Compression',30,...
            'MarkerEdgeColor',cset,'MarkerFaceColor','none','PointSize',30,...
            'DataTransform','None');
        box off
        hold on
        %     errorbar(1:3,mean_data_c,sem_data_c,'k','linestyle','none','capsize',15,'linewidth',lw);
        errorbar(1:3,mean_data_c,sem_data_c,'k','linestyle','none','linewidth',lw);
        for xx = 1:3
            line([xx-.25 xx+.25],[mean_data_c(xx) mean_data_c(xx)],'linewidth',lw,...
                'color','k');
        end
        line([0 4],[0 0],'linestyle','--','color',[.7 .7 .7],'linewidth',2);
        set(gca,'ylim',yls);
        ylabel('% change in FR across state');
        
        [p_kw_c,~,st_kw_c] = kruskalwallis(c_data,[],'off');
        c_c = multcompare(st_kw_c,'display','off');
    end
    
    %% p-values
    
    p_nr1c_pack = signrank(NR1_p_change_c,0)
    p_remc_pack = signrank(REM_p_change_c,0)
    p_nr2c_pack = signrank(NR2_p_change_c,0)
    
    p_nr1d_pack = signrank(NR1_p_change_d,0)
    p_remd_pack = signrank(REM_p_change_d,0)
    p_nr2d_pack = signrank(NR2_p_change_d,0)
    
    p_nr12c_pack = signrank(NR1_FRpacket_c, NR2_FRpacket_c)
    p_nr12d_pack = signrank(NR1_FRpacket_d, NR2_FRpacket_d)
    
    [~,p_nr12c0_t_pack] = ttest(NR12_packchange_c,0)
    [~,p_nr12d0_t_pack] = ttest(NR12_packchange_d,0)
    
end











