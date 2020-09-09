function [cell_set_RSU,rsu_idx] = get_RSUs(cell_set,qthresh,dep,...
    perc_thresh,normalize,require_BLON,dataset,plotcellsep,negposT,tailslopeT,...
    no_intermediates)
% Find RSUs in given set of cells according to parameters.
%
% INPUTS:
%  - cell_set: struct containing CELL variables from which to get RSUs
%  - qthresh: quality threshold (find all cells with qual <= qthresh)
%  - dep: 0 for control hemisphere, 1 for deprived hemisphere
%  - perc_thresh: continuity threshold (find all cells that are online for
%                 % of recording time equal to this number
%                 (perc_thresh*100)
%  - normalize: 1 or 0, normalize to baseline
%  - require_BLON: requires cells to be on during baseline period.
%                  NOTE: blon_hour hardcoded line 101.
%  - dataset: which dataset are you using. This can be 'ER', 'recov',
%             'CPP',etc. This variable is really important! It will
%             determine how long an experiment is considered and how the
%             baseline is determined.
%  - plotcellsep: 1 or 0 (default). Set to 1 to plot cell sep figure.
%  - negposT: neg-pos time threshold for RSUs (default = 0.35)
%  - tailslopeT: tail slope threshold for RSUs (default = 0)
%
% OUTPUTS:
%  - cell_set_RSU: cell structure containing all RSUs selected
%  - rsu_idx: indices to get RSUs by indexing cell_set input variable
%
%
% ATP, September 2017.

% optional arguments
if nargin < 8
    plotcellsep = 0;
end

if nargin < 9
    negposT = .35;
end

if nargin < 10
    tailslopeT = 0;
end

if nargin < 11
    no_intermediates = 0;
end

%colors
cred = [0.85 0.33 0.01];
cblu = [0 0.45 0.74];
cyel = [0.93 0.69 0.13];
cpur = [0.49 0.18 0.55];
cgre = [0.47 0.67 0.19];
ccya = [0.30 0.75 0.93];
cmar = [0.64 0.08 0.18];
cgray = [.6 .6 .6];
% analysis params

dep_status = {'CONTROL','DEPRIVED'};
fprintf('\nFinding RSUs for %s hemisphere.\n',dep_status{dep+1});
fprintf('Continuity threshold is: %u%%\n\n',perc_thresh*100);
fprintf('Negpostime threshold is: %.4f\n\n',negposT);
fprintf('Tailslope threshold is: %.4f\n\n',tailslopeT);
if normalize && ~require_BLON
    warning('WARNING: asked to normalize but not requiring cells to be ON at baseline may result in errors.');
end

% create arrays for indexing of cells
quality = [cell_set.quality];
deprived = [cell_set.deprived];
neg_pos_time = [cell_set.neg_pos_time];
tail_slope = [cell_set.tailSlope];

for dd = 1:length(cell_set)
%     fprintf('Length of cell_set: %u.\n',length(cell_set));
    totT = [];
    for uu = 1:length(cell_set(dd).onTime)
        totT(uu) = cell_set(dd).offTime(uu) - cell_set(dd).onTime(uu);
    end
    
    switch dataset
        case 'ER'
            if regexp(cell_set(dd).animal,regexptranslate('wildcard','KH*'))
                first_number = str2double(regexp(cell_set(dd).animal,'\d','match','once'));
                if first_number > 5
                    daymax = 10.5;
                else
                    daymax = 8;
                end
            elseif regexp(cell_set(dd).animal,regexptranslate('wildcard','AT*'))
                daymax = 11.5;
            elseif strcmp(cell_set(dd).animal,'AT09')
                daymax = 10.5;
            end
            daysRecorded = daymax - cell_set(dd).dayStart;
            
        case 'CPP'
            daymax = 11.5;
            daysRecorded = daymax - cell_set(dd).dayStart;
            
        case {'recov','CPPrecov','CPP2recov','SD_NEW'}
            % fixed number of recorded days, from MD4 to ER4 (included)
            daysRecorded = 5;
        case 'SHANK24'
            daysRecorded = 1;
            % 24h chunk clustering
        case 'SHANK48'
            daysRecorded = 2;
        case 'SHANK72'
            daysRecorded = 3;
        case 'SHANK_cont'
            daysRecorded = 8.5;
    end
    
    
    if sum(totT) < (perc_thresh * daysRecorded*24*3600)
        percentOn(dd)        = 0;
    else
        if cell_set(dd).dayStart > 2
            if sum(totT) < (perc_thresh * daysRecorded*24*3600)
                percentOn(dd) = 0;
            else
                percentOn(dd) = 1;
            end
        else
            %             if strcmp(cell_set(dd).animal,'AT14')
            %                 fprintf('cell %u percent on\n',dd);
            %                 ongood1 = [ongood1; dd];
            %             end
            percentOn(dd) = 1;
        end
    end
    
    % this is old - REVIEW
    if require_BLON
        % cell must come online before start of "baseline" | or half of baseline?
        % "Baseline" is 36-60 hours for animals starting before MD1
        % "Baseline" is MD4 for AT12 (starting at MD3)
        switch dataset
            case {'ER','CPP'}
                timeon = cell_set(dd).onTime(1) + cell_set(dd).dayStart*24*3600;
                blon_hour = 48;
                if timeon <= blon_hour*3600
                    %             if strcmp(cell_set(dd).animal,'AT14')
                    %                 fprintf('cell %u on at: %.1f\n',dd,cell_set(dd).onTime(1)/3600);
                    %                 ongood2 = [ongood2; dd];
                    %             end
                    onatBase(dd) = 1;
                elseif cell_set(dd).onTime(1) <= 6.0*24*3600 && cell_set(dd).dayStart > 2
                    onatBase(dd) = 1;
                else
                    onatBase(dd) = 0;
                end
                
            case {'recov','CPPrecov','CPP2recov','SD_NEW'}
                timeon = cell_set(dd).time(1);
                if timeon <= 6.5*24*3600
                    onatBase(dd) = 1;
                else
                    onatBase(dd) = 0;
                end
            case {'SHANK_cont'}
                timeon = cell_set(dd).onTime(1);
                if cell_set(dd).dayStart > 0
                    fprinf('Strange daystart value found.\n');
                    keyboard
                end
                blon_hour = 54;
                if timeon <= blon_hour*3600
                    onatBase(dd) = 1;
                else
                    onatBase(dd) = 0;
                end
                
        end
    else
        onatBase(dd) = 1;
    end
    
end

depidx = []; qualidx = []; percidx = []; idx0 = [];
% deprived, quality, percentageOn, baselineOn
depidx  = find(deprived == dep);
qualidx = find(quality <= qthresh);
percidx = find(percentOn == 1);
blidx = find(onatBase == 1);
% intersect them
idx00 = intersect(qualidx,depidx);
idx01 = intersect(blidx,percidx);
idx0 = intersect(idx00,idx01);

% Cell classification into RSU/pFS cells
% Most recent negposthresh, tailslope and halfwidth should all be
% independent of interpolation factor
idx2 = []; idx3 = []; idx4 = []; idx5 = []; ns = [];

for cellno = idx0
    
    negpos_thresh = negposT;
    fs_negpos_thresh = negposT;
    tailSlope_thresh = tailslopeT;
    fs_tailSlope_thresh = tailslopeT;
    
    if cell_set(cellno).neg_pos_time >= negpos_thresh
        idx2(end+1) = cellno;
    elseif cell_set(cellno).neg_pos_time < fs_negpos_thresh
        idx3(end+1) = cellno;
    end
    
    if cell_set(cellno).tailSlope >= tailSlope_thresh
        idx4(end+1) = cellno;
    elseif cell_set(cellno).tailSlope < fs_tailSlope_thresh
        idx5(end+1) = cellno;
    end
end

%% plot cell sep
cellsep_count = 0;
if plotcellsep
    
    for cx = idx0
        cellsep_count = cellsep_count + 1;
        cellsep_wf{cellsep_count} = cell_set(cx).scaledWF;
        cellsep_nt(cellsep_count) = neg_pos_time(cx);
        cellsep_ts(cellsep_count) = tail_slope(cx);
        if neg_pos_time(cx) >= negpos_thresh && tail_slope(cx) >= tailSlope_thresh
            cellsep_color{cellsep_count} = cblu;
        elseif neg_pos_time(cx) < fs_negpos_thresh && tail_slope(cx) < fs_tailSlope_thresh
            cellsep_color{cellsep_count} = cyel;
        else
            if no_intermediates
                cellsep_color{cellsep_count} = NaN;
            else
                cellsep_color{cellsep_count} = cgray;
            end
        end
    end
    
    
    cellsepfig = figure(222);
    set(cellsepfig,'visible','off','numbertitle','off');
    for cc = 1:cellsep_count
        s1 = subplot(1,2,1); hold on;
        try
            if ~isnan(cellsep_color{cc})
                plot(s1,cellsep_wf{cc},...
                    'color',cellsep_color{cc},'linewidth',2);
            end
        catch
            keyboard;
        end
        s2 = subplot(1,2,2); hold on;
        if ~isnan(cellsep_color{cc})
            plot(s2,cellsep_nt(cc),cellsep_ts(cc),'o','color',cellsep_color{cc},...
                'MarkerFaceColor',cellsep_color{cc},'MarkerSize',10);
        end
    end
    line([fs_negpos_thresh fs_negpos_thresh],[-0.1 fs_tailSlope_thresh],...
        'color','k','linestyle','--')
    line([0 fs_negpos_thresh],[fs_tailSlope_thresh fs_tailSlope_thresh],...
        'color','k','linestyle','--')
    line([negpos_thresh negpos_thresh],[tailSlope_thresh 0.1],...
        'color','k','linestyle','--')
    line([negpos_thresh 1],[tailSlope_thresh tailSlope_thresh],...
        'color','k','linestyle','--')
    set(s2,'ylim',[-0.1 0.1],'xlim',[0.2 0.60]);
    set(cellsepfig,'units','normalized','position',[0.1 0.1 0.8 0.9],'visible','on');
    set(s1,'fontsize',18,'linewidth',2);
    set(s2,'fontsize',18,'linewidth',2);
    s1.XAxis.Color = 'k';
    s2.XAxis.Color = 'k';
    
    
    uiwait(gcf);
end


rsu_idx = [];
rsu_idx = intersect(intersect(idx0,idx2),idx4);
fprintf('Found %u RSUs in this dataset.\n\n',numel(rsu_idx));
cell_set_RSU = cell_set(rsu_idx);
