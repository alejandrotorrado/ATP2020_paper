%% OnOffTimes_RECOV
%
% Alejandro Torrado Pacheco - 2018
% 
% This script is used to set On/Off times for putative single units
% recorded in the last 5 days of an eye-reopening experiment (MD4 - ER4).
% The script uses an interactive window that pulls up the FR and ISI
% contamination of the neuron over the 5-day period above. The user can
% then set on/off times using one of several options:
%
%   _Press 'G': good. On/off times will just be the start and end of the
%               recording.
%   _Press 'Spacebar': custom on/off times. This will pull-up a cursor to
%                      set the ON time and then another to set the OFF time
%   _Press 'Q': change quality. Sometimes a cell's spike sorting quality
%               should be changed (e.g. high ISI contamination). Press Q to
%               do this and set the new quality manually.
%   _Press 'X': close and quit

% load the file with the spike sorting output
load_data = load('AT66_MASTER_CHUNKS.mat'); % example
% in my data, the sorting output is a structure called CHUNKS containing a
% CELL variable with all of the data for days MD4-ER4
CELL = load_data.CHUNKS.CELL; 

% clear the rest of the workspace
clearvars -except *CONTCELL* *CELL*
clc

% setup colors for plotting
blue = [0 0.44 0.75];
lightblue = [147 187 218]./255;
red = [0.85 0.33 0.10];
lightred = [236 176 146]./255;
yellow = [0.93 0.69 0.13];
lightyel = [255 232 175]./255;
purple = [0.49 0.18 0.56];

% BIN SIZE for FR calculation
binsz = 60*5; % secs

% keep CELL in case something goes wrong, so you can restart without
% reloading the whole thing
cells = CELL;
badcells = [];

% day start in this case is the expt start day relative to the 3-day
% baseline, such that BL1 = 0. For instance if you started the experiment
% on BL2, you would set daystart = 1. If you started it on BL3 you would
% set daystart = 2, etc.
% For most of my data daystart = 2
daystart = input('Day start?  ');
[cells(:).dayStart] = deal(daystart);

% This aligns everything to 7:30 am on the hypothetical BL1
mainStart = (6-daystart)*24*3600; % start of MD4

% start and end times for the periof of data we are interested in (start of
% MD4 to the end of ER4)
tstart = (6-daystart)*24*3600;
tend = (11-daystart)*24*3600;

% thresholds. These mean that if you click an ON-time that is within the
% low_thresh of the start_time, it will just set the ON-time equal to the
% start time. Same for hi_thresh/end_time
low_thresh = 3*3600/binsz;
hi_thresh = (tend-tstart-3*3600)/binsz;

% bins for isi calculation
isi_edges = 0:0.001:2;

%% Main loop that pulls up interactive window
while 1
    % loop through cells
    for cc = 1:size(cells,2)
        
        % only look at putative single and multi unit data (no noise)
        if cells(cc).quality < 4
            % if daystart is not set, get it from the cells
            if ~exist('daystart','var')
                daystart = cells(cc).dayStart;
                mainStart = (6-daystart)*24*3600; % start of MD4
                tstart = (6-daystart)*24*3600;
                tend = (11-daystart)*24*3600;
                low_thresh = 3*3600/binsz;
                hi_thresh = (tend-tstart-3*3600)/binsz;
                clear daystart
            end
            
            % bin edges for FR calculation
            binedges = tstart:binsz:tend;
            % compute cell FR
            rate = histc(cells(cc).time,binedges)./binsz;
            % because of the way the bins are set up, remove last point
            rate(end) = [];
            
            % Calculate hourly ISI contamination for each FR bin
            for xx=1:size(binedges,2)-1
                startedge = binedges(1,xx);
                endedge = binedges(1,xx+1);
                
                temptimes = cells(cc).time;
                temptimes = temptimes(temptimes >= startedge & temptimes < endedge);
                
                chunk_isi{xx} = diff(temptimes);
                
                % find isi contamination by chunk
                isi_cont(xx,1) = sum(chunk_isi{xx}<=0.002)/numel(chunk_isi{xx});
                
                clear temptimes
            end
            
            
            %% PLOTTING
            mainfig = figure();
            set(mainfig,'units','normalized','position',[0.2 0.1 0.75 0.8]);
            
            % top plot - firing rate
            s1 = subplot(2,1,1);
            % plot the light/dark rectangles
            edgesD = (12*3600)/binsz : (24*3600)/binsz : (tend-tstart)/binsz;
            for qq = 1:size(edgesD,2)
                q = rectangle('Position',[edgesD(qq), 0, 12*3600/binsz, 15] );
                set(q,'facecolor',[0.9 0.9 0.9],'linestyle','none');
            end
            hold on
            % plot the actual FR
            plot(rate,'linewidth',2);
            % plot the mean WF of this cell
            meanwf = cells(cc).scaledWF;
            WFx = linspace(ceil(size(rate,1)*.1),ceil(size(rate,1)*.5),size(meanwf,2));
            WFy = (meanwf + 1) * 0.5*max(rate) + 0.6*max(rate);
            plot(WFx,WFy,'color',[0.0 0.6 1]);
            % some text to help you know what you are looking at...
            title(sprintf('Cell %u out of %u.\n Animal: %s. Quality: %u.\n',...
                cc,size(cells,2),cells(cc).animal,cells(cc).quality));
            set(gca,'xtick',0:(12*3600)/binsz:(tend-tstart)/binsz,'xticklabel',...
                {'0','MD4','24','ER1','48','ER2','72','ER3','96','ER4','120','ER5','144'});
            set(gca,'xlim',[0 (tend-tstart)/binsz],'ylim',[0 1.4*max(rate)]);
            
            % on/off time threshold indicators
            line([low_thresh low_thresh],[0 1.4*max(rate)],...
                'color',red,'linestyle','--');
            line([hi_thresh hi_thresh],[0 1.4*max(rate)],...
                'color',red,'linestyle','--');
            
            
            % bottom plot shows the ISI contamination in each FR bin. It
            % also has text showing the average ISI contamination (mean
            % across all bins)
            s2 = subplot(2,1,2); hold on;
            plot(isi_cont,'color',yellow);
            title('ISI contamination (<= 2 ms)');
            set(gca,'ylim',[0 min(0.2,max(isi_cont))]);
            set(gca,'xtick',0:(12*3600)/binsz:(tend-tstart)/binsz,'xticklabel',0:12:24*7);
            set(gca,'xlim',[0 (tend-tstart)/binsz]);
            ylims = get(gca,'ylim'); xlims = get(gca,'xlim');
            isistr = sprintf('Mean ISI contamination (<=2ms) is %.2f%%.',100*nanmean(isi_cont));
            ampstr = sprintf('Mean WF amplitude is %.0f uV.',min(cells(cc).meantrace));
            text(s2,xlims(2)*0.4,ylims(2)*0.78,isistr,'fontsize',14);
            text(s2,xlims(2)*0.4,ylims(2)*0.7,ampstr,'fontsize',14);
            
            % Text instructions to use the interactive on/off times window
            instr_string = {'Press SPACEBAR to set On/Off times. Press G to keep current On/Off times.'; ...
                'Press Q to change cell quality. Press X to quit.'};
            annotation('textbox',[0.02 0.85 0.8 0.15],'String',instr_string,...
                'FitBoxToText','on','fontsize',12);
            
            % wait for button press. This is where the user decides what to
            % do
            goodbutton = 0;
            while ~goodbutton
                waitforbuttonpress;
                key = get(gcf,'CurrentKey');
                % spacebar - set custom on/off times
                if(strcmp (key , 'space'))
                    goodbutton = 1;
                    % NOW SET ON/OFF TIMES
                    
                    subplot(2,1,1)
                    onT = []; offT = [];
                    
                    while isempty(onT)
                        textON = text((17*3600/binsz),1.25*max(rate),'Set ON time.','visible','on','fontsize',20);
                        [onT,~] = ginput(1);
                        if onT <= low_thresh
                            onT = 0;
                        end
                        subplot(2,1,1); hold on;
                        line([onT onT],[0 1.4*max(rate)],...
                            'color',purple,'linestyle','-','linewidth',1.5);
                    end
                    
                    while isempty(offT)
                        set(textON,'visible','off');
                        textOFF = text((17*3600/binsz),1.25*max(rate),'Set OFF time.','visible','on','fontsize',20);
                        [offT,~] = ginput(1);
                        set(textOFF,'visible','off');
                        if offT >= hi_thresh
                            offT = (tend-tstart)/binsz;
                        end
                        subplot(2,1,1); hold on;
                        line([offT offT],[0 1.4*max(rate)],...
                            'color',purple,'linestyle','-','linewidth',1.5);
                    end
                    pause(0.8);
                    set(mainfig,'visible','off');
                    
                    % sometimes you want more than 1 pair of on/off times
                    notdone = 1;
                    while notdone
                        setmoretimes = input('Select more on/off times? (1=yes,0=no):  ');
                        
                        if setmoretimes
                            set(mainfig,'visible','on');
                            x_t1 = []; x_t2 = [];
                            while isempty(x_t1)
                                set(textOFF,'visible','off');
                                set(textON,'visible','on');
                                [x_t1,~] = ginput(1);
                                if x_t1 <= low_thresh
                                    x_t1 = 0;
                                end
                                subplot(2,1,1); hold on;
                                line([x_t1 x_t1],[0 1.4*max(rate)],...
                                    'color',purple,'linestyle','-','linewidth',1.5);
                            end
                            while isempty(x_t2)
                                set(textON,'visible','off');
                                set(textOFF,'visible','on');
                                [x_t2,~] = ginput(1);
                                if x_t2 >= hi_thresh
                                    x_t2 = (tend-tstart)/binsz;
                                end
                                subplot(2,1,1); hold on;
                                line([x_t2 x_t2],[0 1.4*max(rate)],...
                                    'color',purple,'linestyle','-','linewidth',1.5);
                                set(textOFF,'visible','off');
                            end
                            pause(0.8);
                            set(mainfig,'visible','off');
                            
                            onT = [onT; x_t1];
                            offT = [offT; x_t2];
                            
                            clear x_t*
                        else
                            notdone = 0;
                        end
                    end
                    
                    % save the data
                    cells(cc).onTime = [];
                    cells(cc).offTime = [];
                    cells(cc).onTime = mainStart + onT.*binsz;
                    cells(cc).offTime = mainStart + offT.*binsz;
                    
                    % when done, close figure and keep going
                    close(mainfig);
                    
                    % 'G' press : keep on/off times as they are (beginning
                    % and end of the 5-day period)
                elseif strcmp(key, 'g')
                    % if cell is good, leave on/off times as they are and keep
                    % going with the next cell
                    goodbutton = 1;
                    close(mainfig);
                    
                    % Q: change quality
                elseif strcmp(key, 'q')
                    % press q to change quality
                    goodbutton = 1;
                    updatedqual = input('New cell quality:  ');
                    cells(cc).quality = updatedqual;
                    close(mainfig);
                    
                    % X: close and quit
                elseif strcmp(key, 'x')
                    % press x to quit
                    goodbutton = 1;
                    close all;
                    break;
                else
                    disp('Try again,');
                    goodbutton = 0;
                end
            end
            
        end
        
    end
    
    % if user quits or last cell is reached, exit loop
    if strcmp(key,'x')
        disp('User quit');
        break;
    end
    
    if cc==size(cells,2)
        disp('Reached last cell! Quitting.');
        break;
    end
    
end

% The code uses a keyboard statement to stop here. Then you can manually do
% the saving part, and add the new CELL variable with ON/OFF times to a
% CONTCELL variable if you have one loaded.
keyboard;

%% saving

do_save = 1;
if do_save
    clear CELL
    CELL = cells;
    animal = CELL(1).animal;
    % saving directory and filename
    sd = ['Z:\ANIMALDATA\MLS_DATA\EyeReopen_only\' animal];
    sf = [animal '_MASTER_CHUNKS_OnOff.mat'];
    % do the saving
    tic
    save([sd filesep sf],'CELL','-v7.3');
    toc
end


%% add to CONTCELL variable
% if you have a CONTCELL variable loaded in the workspace, you can use
% these lines of code to add the new CELL variable to the main one
clc
clearvars -except CONTCELL* CELL*
qs = [CELL.quality];
badz = find(qs>2);
CELL(badz) = [];
% CELL2 = rmfield(CELL,'psort'); % sometimes you need to remove psort for consistency
CONTCELL_SD_NEW.MASTER = [CONTCELL_SD_NEW.MASTER CELL];
clearvars -except CONTCELL*
clc

