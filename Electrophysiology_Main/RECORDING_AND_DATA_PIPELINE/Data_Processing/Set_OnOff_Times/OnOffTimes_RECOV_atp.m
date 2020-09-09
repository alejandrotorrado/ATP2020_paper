% clearvars -except CONTCHUNKS_*h *CONTCELL* *CELL*

% colors
blue = [0 0.44 0.75];
lightblue = [147 187 218]./255;
red = [0.85 0.33 0.10];
lightred = [236 176 146]./255;
yellow = [0.93 0.69 0.13];
lightyel = [255 232 175]./255;
purple = [0.49 0.18 0.56];

% BIN SIZE
binsz = 60*5; % secs

% atpCONTCELL = CONTCELL_recov;
% 
% 
% cells = atpCONTCELL.MASTER;
cells = CELL;
badcells = [];
% tstart = atpCONTCELL.START; % start on MD3 (but will bump start time to MD4)
% tend = atpCONTCELL.END + 24*3600; % end on ER5

daystart = input('Day start?  ');
[cells(:).dayStart] = deal(daystart);

mainStart = (6-daystart)*24*3600; % start of MD4

tstart = (6-daystart)*24*3600;
tend = (11-daystart)*24*3600;

low_thresh = 3*3600/binsz;
hi_thresh = (tend-tstart-3*3600)/binsz;


isi_edges = 0:0.001:2;

while 1
    for cc = 1:size(cells,2)
        
        if cells(cc).quality < 3
            if ~exist('daystart','var')
                daystart = cells(cc).dayStart;
                mainStart = (6-daystart)*24*3600; % start of MD4
                tstart = (6-daystart)*24*3600;
                tend = (11-daystart)*24*3600;
                low_thresh = 3*3600/binsz;
                hi_thresh = (tend-tstart-3*3600)/binsz;
                clear daystart
            end
            
            binedges = tstart:binsz:tend;
            rate = histc(cells(cc).time,binedges)./binsz;
            rate(end) = [];
            
            
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
            
            mainfig = figure();
            set(mainfig,'units','normalized','position',[0.2 0.1 0.75 0.8]);
            
            s1 = subplot(2,1,1);
            edgesD = (12*3600)/binsz : (24*3600)/binsz : (tend-tstart)/binsz;
            for qq = 1:size(edgesD,2)
                q = rectangle('Position',[edgesD(qq), 0, 12*3600/binsz, 15] );
                set(q,'facecolor',[0.9 0.9 0.9],'linestyle','none');
            end
            hold on
            plot(rate,'linewidth',2);
            meanwf = cells(cc).scaledWF;
            WFx = linspace(ceil(size(rate,1)*.1),ceil(size(rate,1)*.5),size(meanwf,2));
            WFy = (meanwf + 1) * 0.5*max(rate) + 0.6*max(rate);
            plot(WFx,WFy,'color',[0.0 0.6 1]);
            title(sprintf('Cell %u out of %u.\n Animal: %s. Quality: %u.\n',...
                cc,size(cells,2),cells(cc).animal,cells(cc).quality));
            set(gca,'xtick',0:(12*3600)/binsz:(tend-tstart)/binsz,'xticklabel',...
                {'0','MD4','24','ER1','48','ER2','72','ER3','96','ER4','120','ER5','144'});
            set(gca,'xlim',[0 (tend-tstart)/binsz],'ylim',[0 1.4*max(rate)]);
            
            line([low_thresh low_thresh],[0 1.4*max(rate)],...
                'color',red,'linestyle','--');
            line([hi_thresh hi_thresh],[0 1.4*max(rate)],...
                'color',red,'linestyle','--');
            
            
            
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
            
            instr_string = {'Press SPACEBAR to set On/Off times. Press G to keep current On/Off times.'; ...
                'Press Q to change cell quality. Press X to quit.'};
            annotation('textbox',[0.02 0.85 0.8 0.15],'String',instr_string,...
                'FitBoxToText','on','fontsize',12);
            
            goodbutton = 0;
            while ~goodbutton
                waitforbuttonpress;
                key = get(gcf,'CurrentKey');
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
                    
                    
                    cells(cc).onTime = [];
                    cells(cc).offTime = [];
                    cells(cc).onTime = mainStart + onT.*binsz;
                    cells(cc).offTime = mainStart + offT.*binsz;
                    
                    % when done, close figure and keep going
                    close(mainfig);
                    
                elseif strcmp(key, 'g')
                    % if cell is good, leave on/off times as they are and keep
                    % going with the next cell
                    goodbutton = 1;
                    close(mainfig);
                elseif strcmp(key, 'q')
                    % press q to change quality
                    goodbutton = 1;
                    updatedqual = input('New cell quality:  ');
                    cells(cc).quality = updatedqual;
                    close(mainfig);
                    
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
    
    if strcmp(key,'x')
        disp('User quit');
        break;
    end
    
    if cc==size(cells,2)
        disp('Reached last cell! Quitting.');
        break;
    end
    
end

% for xx=1:length(cells)
%     if cells(xx).offTime < cells(xx).onTime
%         cells(xx).offTime = atpCONTCELL.START + cells(xx).offTime;
%     end
% end


keyboard;

%% saving

do_save = 1;
if do_save
    clear CELL
    CELL = cells;
    animal = CELL(1).animal;
    sd = ['Z:\ANIMALDATA\MLS_DATA\EyeReopen_only\' animal];
    sf = [animal '_MASTER_CHUNKS_OnOff.mat'];
    tic
    save([sd filesep sf],'CELL','-v7.3');
    toc
%     atpCONTCELL.CELL = [];
%     atpCONTCELL.CELL = CELL;
%     oldrecov = CONTCELL_recov;
%     clear CONTCELL_recov;
%     CONTCELL_recov = atpCONTCELL;
%     sd = '/Users/atorrado/Desktop/MLS_DATA';
%     sf = 'CONTCELL_recov_MLS_v2.mat';
%     save([sd filesep sf],'CONTCELL_recov','-v7.3');
end


%%
clc
clearvars -except CONTCELL* CELL*
qs = [CELL.quality];
badz = find(qs>2);
CELL(badz) = [];
CONTCELL_CPP2recov.MASTER = [CONTCELL_CPP2recov.MASTER CELL];
clearvars -except CONTCELL*
clc

