%% S/W-DENSE ANALYSIS
%
% Alejandro Torrado Pacheco - 2017
%
% This code is run to obtain the SW-dense data that can then be used to
% plot Fig 4. Plot the data using: swdense_postprocess.m
%
% Note: A variant of this code can be used to produce the data for Fig S7.
%
% Break up statetimes in 4-hour blocks (30 min step?) and slide window
% looking for sleep- (S-) or wake- (W-) dense periods, defined as periods
% of time in which the animal spends at least 65%/70% of its time in S or
% W. Then find mean change in FR associated with these periods.
%
% steps:
% for each animal, figure out S-dense and W-dense states.
% Then loop through cells and get change in FR.

%% SETUP
% clear workspace
clearvars -except CONTCELL*
close all
clc

% plotting colors
s_color = [.8 .8 .8];
w_color = [89, 150, 213]./255;
c_rem   = [25 181 149]./255;
c_nrem  = [131 49 146]./255;
c_aw    = [201 28 101]./255;
c_qw    = [247 148 41]./255;

% time program
prog_start = tic;

% debugging flags
verbose = 0;
plot_verbose = 0;
rate_debug = 0;

% pick dataset - remember to load CONTCELL variable first
dataset = 'recov';
switch dataset
    case 'ER'
        CELL = CONTCELL_ER.MASTER;
        STATES = CONTCELL_ER.STATETIMES;
    case 'recov'
        CELL = CONTCELL_recov.MASTER;
        STATES = CONTCELL_recov.STATETIMES;
end
% use this to speed up the code by loading recov data
loadrecov = 1;
if loadrecov
    if ismac
        loadfile = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/recov_analysis_norm.mat';
    elseif ispc
        loadfile = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\recov_analysis_norm.mat';
    end
    rload = load(loadfile);
    recov = rload.recov_analysis;
end
G_bin = recov.G_bin;

% choose re-opened or control
dep = 0;
if dep == 1
    goodcells = recov.DEPRIVED.RSU_idx;
    FRbycell_RSU = recov.DEPRIVED.RSU_FRbycell;
    anims = recov.DEPRIVED.RSU_anims;
    n_RSU = recov.DEPRIVED.RSU_count;
elseif dep == 0
    goodcells = recov.CONTROL.RSU_idx;
    FRbycell_RSU = recov.CONTROL.RSU_FRbycell;
    anims = recov.CONTROL.RSU_anims;
    n_RSU = recov.CONTROL.RSU_count;
end

% set animals to be used in each case
if dep
    which_anim = {'AT12','KH67','AT16','AT14','AT24','AT25','AT29','AT27'};
else
    which_anim = {'AT12','KH67','AT14','AT16','AT22','AT24','AT25','AT29','AT27','KH73','KH75','KH72'};
end

all_animals = fieldnames(STATES);

% analysis parameters norm_mode specifies how the FR change is calculated.
% Call B the mean FR at the end of the window, and A the mean FR at the start
% of the window. Then:
%   _frac: FR_change = B/A;
%   _perc: FR_change = 100*(B-A)/A;
%   _change: FR_change = (B-A)/(B+A);
% To reproduce Fig. 4, select 'change'.
norm_mode = 'change'; % perc, frac or change
save_the_data = 1;

% SW dense analysis parameters
% NOTE: all of these parameters could be
% varied. The variant of this script loops through various threshold and
% window size values - see Fig S6.
dense_threshold = 70; % percentage threshold to consider epoch 'dense'
window_hours = 2.5; % size of sliding window in hours
window_step = 15; % sliding step size in minutes
% chunk of time at the start and end of each 'dense' window over which to
% calculate FR to estimate change across 'dense' window
FR_time_window = 15; % minutes


%% Loop through all animals
for aa = 1:numel(all_animals)
    animal = all_animals{aa};
    % if animal is within the list of good ones...
    if any(strcmp(animal,which_anim)) && any(strcmp(animal,anims))
        fprintf('Processing animal: %s\n\n',animal);
        
        % pull state times
        states = STATES.(animal);
        
        % clear double states
        testD   = diff(states(:,1));
        killem  = find(testD == 0);
        killem  = killem+1;
        states(killem,:) = [];
        
        % find all cells from that animal
        animal_cells = find(strcmp({CELL.animal},animal));
        
        % align to BL1 using first cell
        statetimes = adjust_statetimes_recov(states,CELL(animal_cells(1)),CONTCELL_recov.STATETIMES.DAYSTART.(animal));
        
        % Get all the cells for this animal
        these_cells = CELL(animal_cells);
        
        %% Define time parameters
        % define homeostatic recovery period
        epoch_start = 8.0; % 8 is start of ER2, 9 is start of ER3.
        epoch_end = 9.0; % 10 is start of ER4

        % edges of homeostatic period in seconds, aligned to BL1
        t0 = epoch_start*24*3600;
        t1 = epoch_end*24*3600;
        
        % time step size
        t_step_min = window_step;
        t_step = t_step_min*60;
        
        % time window size
        t_win_hours = window_hours;
        t_win = t_win_hours*3600;
        
        % threshold to consider an epoch S/W-dense
        t_dense_thresh = dense_threshold;

        % set time over which to calc FR in SW-dense window
        fr_time_win = FR_time_window*60; % seconds (minutes*60)
        
        %% Loop through states for this animal
        % initialize variables
        S_dense_counter{aa} = 0;
        W_dense_counter{aa} = 0;
        S_windows{aa} = [];
        W_windows{aa} = [];
        S_windows_slice{aa} = [];
        W_windows_slice{aa} = [];
        S_win_change_FR{aa} = [];
        W_win_change_FR{aa} = [];
        S_dense_dur{aa} = [];
        W_dense_dur{aa} = [];
        window_counter = 0;
        total_counter = 0;
        skip_to_time = 0;
        keep_going = 1;
        
        % set up the loop
        while keep_going
            % increase time window counter
            window_counter = window_counter + 1;
            total_counter = total_counter + 1;
            fprintf('Window counter: %u\n\n',window_counter);
            
            
            % define current time window
            if skip_to_time == 0
                % if not skipping ahead, just keep increasing steps
                t0_win = t0 + t_step * (window_counter - 1);
            else
                % otherwise simply skip ahead to end of last window
                t0_win = t0 + t_step * (window_counter - 1);
                skip_to_time = 0;
            end
            
            % get end of current window
            t1_win = t0_win + t_win;
            % if beyond specified time limit
            if t1_win >= t1 - t_step
                % stop the while loop
                keep_going = 0;
            end
            
            % print out
            fprintf('Processing window of time starting at t = %.4f\n',t0_win./(24*3600));
            
            
            % look at states in that window - edge states not included yet
            win_states_first = find(statetimes(:,2) > t0_win,1,'first');
            win_states_last = find(statetimes(:,2) < t1_win,1,'last');
            
            % include edge states in the slice
            slice_up_to = min(win_states_last+1,size(statetimes,1));
            state_slice = statetimes((win_states_first-1):slice_up_to,:);
            n_states = size(state_slice,1);
            
            % initialize duration counters
            sleep_dur = 0;
            wake_dur = 0;
            
            % loop through states and compute durations
            fprintf('\nCalculating state durations...\n');
            for state = 1:(n_states-1)
                if verbose, fprintf('State %u out of %u.\t-\t',state,(n_states-1)); end
                % loop through states and get duration
                % also deal with edge states
                if state > 1 && state < (n_states-1)
                    this_dur = state_slice(state+1,2) - state_slice(state,2);
                elseif state == 1
                    this_dur = state_slice(state+1,2) - t0_win;
                elseif state == (n_states-1)
                    this_dur = t1_win - state_slice(state,2);
                end
                
                % only consider states 60sec or longer
                if this_dur > 60
                    % find if i's a S or W state, add to appropriate duration counter
                    state_kind = state_slice(state,1);
                    if state_kind < 3
                        % add to sleep
                        if verbose, fprintf('Added to SLEEP.\n'); end
                        sleep_dur = sleep_dur + this_dur;
                    elseif state_kind > 3
                        % add to wake
                        if verbose, fprintf('Added to WAKE.\n'); end
                        wake_dur = wake_dur + this_dur;
                    end
                end
                clear this_dur
            end
            
            % calculate percentages of time spent in each state
            sleep_perc = sleep_dur/t_win;
            wake_perc = wake_dur/t_win;
            
            % print out
            fprintf('\nDone! Breakdown: \n');
            fprintf('SLEEP: %.1f%%\nWAKE: %.1f%%\n',sleep_perc*100,wake_perc*100);
            
            
            S_or_W = 0;
            % check if this epoch was S/W-dense
            if sleep_perc > t_dense_thresh
                fprintf('\n*** Found a SLEEP-dense epoch! ***\n\n');
                S_or_W = 1; % 1 is sleep-dense
                
                S_dense_counter{aa} = S_dense_counter{aa} + 1;
                % if S dense and last state is S state, extend window to end of
                % that state
                if state_slice(end-1,1) < 3
                    new_t1_win = state_slice(end,2);
                    if t1_win - state_slice(end-1,2) < 60
                        t1_win_diff = new_t1_win - state_slice(end-1,2);
                    else
                        t1_win_diff = new_t1_win - t1_win;
                    end
                    sleep_dur = sleep_dur + t1_win_diff;
                    t1_win = new_t1_win;
                end
                S_windows{aa} = [S_windows{aa}; [t0_win t1_win]];
                S_windows_slice{aa}{S_dense_counter{aa}} = state_slice;
                
                % skip to end of this epoch for next iteration
                if t1_win > t1-t_win
                    disp('Breaking loop!');
                    keep_going = 0;
                else
                    if verbose, fprintf('Resetting window counter and t0.\n'); end
                    skip_to_time = t1_win;
                    window_counter = 0;
                    t0 = skip_to_time;
                    S_dense_dur{aa} = [S_dense_dur{aa}; t1_win-t0_win];
                end
                
            elseif wake_perc > t_dense_thresh
                fprintf('\n*** Found a WAKE-dense epoch! ***\n\n');
                S_or_W = 2; % 2 is wake-dense
                W_dense_counter{aa} = W_dense_counter{aa} + 1;
                % if W dense and last state is W state, extend window to end of
                % that state
                if state_slice(end-1,1) > 3
                    new_t1_win = state_slice(end,2);
                    if t1_win - state_slice(end-1,2) < 60
                        t1_win_diff = new_t1_win - state_slice(end-1,2);
                    else
                        t1_win_diff = new_t1_win - t1_win;
                    end
                    wake_dur = wake_dur + t1_win_diff;
                    t1_win = new_t1_win;
                end
                W_windows{aa} = [W_windows{aa}; [t0_win t1_win]];
                W_windows_slice{aa}{W_dense_counter{aa}} = state_slice;
                
              
                
                % skip to end of this epoch for next iteration
                if t1_win > t1-t_win
                    disp('Breaking loop!');
                    keep_going = 0;
                else
                    if verbose, disp('x'), end
                    fprintf('Resetting window counter and t0.\n');
                    skip_to_time = t1_win;
                    window_counter = 0;
                    t0 = skip_to_time;
                    W_dense_dur{aa} = [W_dense_dur{aa}; t1_win-t0_win];
                end
                
            else
                fprintf('Nothing exciting in this epoch.\n\n');
                S_or_W = 0; % 0 is nothing
            end
            
            % if found a S/W-dense epoch, loop through cells and calc FR change
            if S_or_W % this will go ahead if S_or_W is > 0
                % loop through cells and calculate firing rate change
                
                %% re-opened hemisphere
                rsu_count = 0;
                for cc = 1:n_RSU
                    
                    % assign to new variable
                    rate = FRbycell_RSU(cc,:); % already normalized
                    
                    % get on/off times - these are already adjusted to be
                    % aligned to BL1 in CONTCELL_recov
                    cellidx = goodcells(cc);
                    thiscell = CONTCELL_recov.MASTER(cellidx);
                    ontimes = thiscell.onTime;
                    offtimes = thiscell.offTime;
                    timelist(:,1) = ontimes;
                    timelist(:,2) = offtimes;
                    newOnOff = timelist;
                    
                    % check for mismatches between this cell and the
                    % current animal
                    if ~strcmp(thiscell.animal,anims{cc})
                        disp('BIG PROBLEM');
                        keyboard;
                    end
                    
                    % if everything looks good
                    if ~all(isnan(rate)) && strcmp(animal,thiscell.animal)
                        rsu_count = rsu_count + 1;
                        
                        % check that cell is ON during this window
                        % find on time before this window thiscell.onTime
                        closest_on_idx = find(newOnOff(:,1) < t0_win,1,'last');
                        closest_on = newOnOff(closest_on_idx,1);
                        closest_off = newOnOff(closest_on_idx,2);
                        
                        if all(~isempty([closest_on,closest_off])) && closest_on < t0_win && closest_off > t1_win
                            % cell is on
                            if verbose, fprintf('Cell is on in the current window!\n'); end
                            % calculate FR at start and end of window
                            t00 = ceil(t0_win / G_bin);
                            t01 = floor((t0_win + fr_time_win) / G_bin);
                            t10 = ceil((t1_win - fr_time_win) / G_bin);
                            t11 = floor(t1_win / G_bin);
                            % store in main variable
                            win_start_FR{aa}(rsu_count,total_counter) = nanmean(rate(t00:t01));
                            win_end_FR{aa}(rsu_count,total_counter) = nanmean(rate(t10:t11));
                            % calculate FR change over window according to
                            % specified method
                            switch norm_mode
                                case 'frac'
                                    win_change_FR = nanmean(rate(t10:t11)) / nanmean(rate(t00:t01));
                                case 'perc'
                                    win_change_FR = (nanmean(rate(t10:t11)) - nanmean(rate(t00:t01))) / nanmean(rate(t00:t01));
                                case 'change'
                                    win_change_FR = (nanmean(rate(t10:t11)) - nanmean(rate(t00:t01))) / (nanmean(rate(t00:t01)) + nanmean(rate(t10:t11)));
                            end
                            
                            % the plot_verbose flag lets you check the
                            % identified window and the FR of this cell
                            if plot_verbose
                                figure(19218); hold on;
                                win_rate = rate(t00:t11);
                                subplot(2,1,1);
                                plot(win_rate,'linewidth',2);
                                r0 = rectangle('position',...
                                    [0 min(win_rate)-.1 t01-t00 .1+(max(win_rate)-min(win_rate))]);
                                r1 = rectangle('position',...
                                    [t10-t00 min(win_rate)-0.1 t11-t10 .1+(max(win_rate)-min(win_rate))]);
                                tfr = text(gca,40,max(win_rate),sprintf('FR change: %.2f',win_change_FR),...
                                    'visible','on');
                                title(sprintf('Animal: %s; rsu count: %u. S(1) or W(2): %u',animal,rsu_count,S_or_W));
                                subplot(2,1,2);
                                plot(rate);
                                set(gca,'xtick',[0*24*3600/G_bin:12*3600/G_bin:12*24*3600/G_bin],...
                                    'xticklabel',[0:12:12*24]./24);
                                yl = get(gca,'ylim');
                                line([t0_win t0_win]./G_bin,[0 yl(2)],'linewidth',2,'linestyle','--','color','k');
                                line([t1_win t1_win]./G_bin,[0 yl(2)],'linewidth',2,'linestyle','--','color','k');
                                uiwait(gcf);
                            end
                            if verbose, fprintf('FR change: %.3f.\n',win_change_FR); end
                            
                            % if it was a S or W dense window, save the
                            % data in a big variable
                            if S_or_W == 1
                                S_win_change_FR{aa}(rsu_count,total_counter) =...
                                    win_change_FR;
                                S_FR_trace{aa}{rsu_count,S_dense_counter{aa}} = rate(t00:t11);
                            elseif S_or_W == 2
                                W_win_change_FR{aa}(rsu_count,total_counter) =...
                                    win_change_FR;
                                W_FR_trace{aa}{rsu_count,W_dense_counter{aa}} = rate(t00:t11);
                            end
                        else
                            % cell is off
                            win_start_FR{aa}(rsu_count,total_counter) = NaN;
                            win_end_FR{aa}(rsu_count,total_counter) = NaN;
                            switch S_or_W
                                case 1
                                    S_win_change_FR{aa}(rsu_count,total_counter) =...
                                        NaN;
                                case 2
                                    W_win_change_FR{aa}(rsu_count,total_counter) =...
                                        NaN;
                            end
                            if verbose, fprintf('This cell turned off in the current window.\n'); end
                        end
                        
                    else
                        if verbose, fprintf('Bad cell.\n\n'); end
                    end
                    
                    clear timelist
                    
                    
                    % end of loop through cells for this animal and this state
                end
                
                % end of if S_or_W statement
            end
            
            % end of loop through states for this animal
        end
        
    else
        fprintf('\n\n Not doing this animal (%s)!\n\n',animal);
    end
    % end of loop through all animals
end

%% compile the data
% WAKE
FR_change_by_cell_W = [];
for aa = 1:size(W_win_change_FR,2)
    thisanim_FR = W_win_change_FR{aa};
    thisanim_FR(thisanim_FR==0) = NaN;
    ncells = size(thisanim_FR,1);
    mean_bycell = nanmean(thisanim_FR,2);
    FR_change_by_cell_W = [FR_change_by_cell_W; mean_bycell];
    clear mean_bycell thisanim_FR ncells
end

% SLEEP
FR_change_by_cell_S = [];
for aa = 1:size(S_win_change_FR,2)
    thisanim_FR = S_win_change_FR{aa};
    thisanim_FR(thisanim_FR==0) = NaN;
    ncells = size(thisanim_FR,1);
    mean_bycell = nanmean(thisanim_FR,2);
    FR_change_by_cell_S = [FR_change_by_cell_S; mean_bycell];
    clear mean_bycell thisanim_FR ncells
end

% OVERALL MEANS
mean_FR_change_S = nanmean(FR_change_by_cell_S);
sem_FR_change_S = nanstd(FR_change_by_cell_S)./sqrt(size(FR_change_by_cell_S,1)-1);

mean_FR_change_W = nanmean(FR_change_by_cell_W);
sem_FR_change_W = nanstd(FR_change_by_cell_W)./sqrt(size(FR_change_by_cell_W,1)-1);

fprintf(['\nWAKE mean: %.4f ' char(177) ' %.4f\n'],mean_FR_change_W,sem_FR_change_W);
fprintf(['\nSLEEP mean: %.4f ' char(177) ' %.4f\n'],mean_FR_change_S,sem_FR_change_S);

prog_end = toc(prog_start);

fprintf('\n\n\tThat took %.2f seconds.\n\n',prog_end);

% STATS  - one-sample t-test vs mean value (defined according to
% calculation mode)
switch norm_mode
    case {'frac'}
        mean_to_test = 1;
    case {'perc','change'}
        mean_to_test = 0;
end
% do the t-test
[h,p_val_S] = ttest(FR_change_by_cell_S,mean_to_test)
[h,p_val_W] = ttest(FR_change_by_cell_W,mean_to_test)

%% save the data

dep_status = {'CTRL','DEP'};
if save_the_data
    % this is currently an empty folder - available to check the code by
    % making a dataset and then post-processing for plotting
    if ismac
        save_dir = '/Volumes/turrigiano-lab/ATP_MAIN/DATA/Dissertation_Data/ER2020/ER_Fig4/swdense_data/';
    elseif ispc
        save_dir = 'Z:\ATP_MAIN\DATA\Dissertation_Data\ER2020\ER_Fig4\swdense_data\';
    end
    % this automatically finds how many files are in the folder and makes a
    % new file with a new integer ID (e.g. SWdense_CTRL_20.mat if there are
    % already 19 CTRL data files in the folder).
    base_save_name = ['SWdense_' dep_status{dep+1}];
    n_files = numel(dir([save_dir filesep base_save_name '*.mat']));
    save_name = [base_save_name '_' num2str(n_files+1) '.mat'];
    save_file = [save_dir filesep save_name];
    
    ss = 1;
    
    % main data outputs
    swdense_data(ss).FR_change_by_cell_S    = FR_change_by_cell_S;
    swdense_data(ss).FR_change_by_cell_W    = FR_change_by_cell_W;
    swdense_data(ss).mean_FR_change_S       = mean_FR_change_S;
    swdense_data(ss).mean_FR_change_W       = mean_FR_change_W;
    swdense_data(ss).sem_FR_change_S        = sem_FR_change_S;
    swdense_data(ss).sem_FR_change_W        = sem_FR_change_W;
    
    swdense_data(ss).frac_type              = frac_to_use;
    swdense_data(ss).S_windows              = S_windows;
    swdense_data(ss).W_windows              = W_windows;
    
    % metadata
    swdense_data(ss).S_dense_episodes       = n_S_dense;
    swdense_data(ss).W_dense_episodes       = n_W_dense;
    
    % analysis params
    swdense_data(ss).fr_time_win            = fr_time_win;
    swdense_data(ss).epoch_start            = epoch_start;
    swdense_data(ss).epoch_end              = epoch_end;
    swdense_data(ss).t_step_min             = t_step_min;
    swdense_data(ss).t_win_hours            = t_win_hours;
    swdense_data(ss).t_dense_thresh         = t_dense_thresh;
    % cell-selection params
    swdense_data(ss).dataset                = dataset;
    swdense_data(ss).G_bin                  = G_bin;
    swdense_data(ss).which_anim             = which_anim;
    
    save(save_file,'swdense_data','-v7.3');
    
end


