function CH_OUT = CHUNKS_addFields(CH,animal,hem,bindir)

CH_OUT = CH; % initially make them the same

fprintf('Adding fields to CHUNKS variable for animal %s.\n',animal);

for xx = 1:size(CH,2)
    
    CELL = CH(xx).CELL;
    
    % clear current on/off times:
    for ee = 1:size(CELL,2);
        CELL(ee).onTime = [];
        CELL(ee).offTime = [];
    end
    
    
    if ~isfield(CELL,'EXPTSTART') && ~exist('tINFO','var');
        fprintf('Fetching experiment start time file for %s...\n',animal);
        load([bindir filesep animal '_expt_START.mat']);
        for cx = 1:length(CELL)
            CELL(cx).EXPTSTART = tINFO;
        end
    else
        for cx = 1:length(CELL)
            CELL(cx).EXPTSTART = tINFO;
        end
    end
    
    valchans = 1:32;
    
    tdiff   = CELL(1).EXPTSTART  - 27000;
    trem    = rem(tdiff,(24*3600));
    
    %here
    for ee = 1:length(CELL);
        CELL(ee).animal     = animal;
        % just set chunk time boundaries as on/off times for every cell
        CELL(ee).onTime     = CH(xx).START;
        CELL(ee).offTime    = CH(xx).END;
        if ~isempty(CELL(ee).time)
            if CELL(ee).time(1) > CELL(ee).EXPTSTART
                CELL(ee).time = CELL(ee).time - CELL(ee).EXPTSTART + trem;
                % removed the else statement here (used to be else: time =
                % time- trem). This is because we already align spike times to
                % 7:30 am in the binarySpikeProcessing_BETA script.
            end
        else
            fprintf('Empty times array! Cell %u, quality %u.\n',ee,CELL(ee).quality);
        end
        CELL(ee).trem       = trem;
        if ~isfield(CELL(ee),'scaledWF') || isempty(CELL(ee).scaledWF)
            CELL(ee).scaledWF = CELL(ee).meantrace / abs(min(CELL(ee).meantrace));
        end
        
        if strcmp(hem,'l') || strcmp(hem,'L') && CELL(ee).channel<17;
            CELL(ee).deprived = 1;
        elseif strcmp(hem,'l') || strcmp(hem,'L') && CELL(ee).channel>16;
            CELL(ee).deprived = 0;
        elseif strcmp(hem,'r') || strcmp(hem,'R') && CELL(ee).channel<17;
            CELL(ee).deprived = 0;
        elseif strcmp(hem,'r') || strcmp(hem,'R') && CELL(ee).channel>16;
            CELL(ee).deprived = 1;
        end
        
        % FIX THIS - - - INTERP NO LONGER EXISTS?
        % halfwidth
        CELL(ee).halfwidth = calc_halfwidth(CELL(ee).scaledWF,3);
        % neg pos time
        CELL(ee).neg_pos_time = calc_negpostime(CELL(ee).scaledWF,3);;
        
		% tail slope
        if size(CELL(ee).scaledWF,2) == 33
            lims = [20 25];
        elseif size(CELL(ee).scaledWF,2) == 161
            lims = [97 122];
        elseif size(CELL(ee).scaledWF,2) == 91
            lims = [57 70];
        elseif size(CELL(ee).scaledWF,2) == 59;
            lims = [40 50];
        end
        xlims = 1:(diff(lims)+1);
        
        slope = polyfit( xlims  , CELL(ee).scaledWF( lims(1) : lims(2) ),1);
        CELL(ee).tailSlope = slope(1);
    end
    
    % reassign
    CH_OUT(xx).CELL = [];
    CH_OUT(xx).CELL = CELL;
    
end
