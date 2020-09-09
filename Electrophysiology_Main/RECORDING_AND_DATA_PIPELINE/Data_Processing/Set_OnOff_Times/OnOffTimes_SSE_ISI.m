clearvars -except *CELL*

scrsz = get(groot,'ScreenSize');

for cc=1:length(RVW1_CONTCELL)
    fprintf('Cell %u out of %u.\n',cc,length(RVW1_CONTCELL));
    ee = RVW1_CONTCELL(cc).OG_number;
    fprintf('OG number: %u.\n',ee);
    
    clear rate sse spikes mrate normrt mwfBeforeCrash mwfAfterCrash onoff
    rate    = RVW1_CONTCELL(cc).hour_rate;
    sse     = RVW1_CONTCELL(cc).sse;
    spikes  = RVW1_CONTCELL(cc).hour_spikes;
    mrate   = nanmean(rate);
    msse    = nanmean(sse);
    normrt  = rate./mrate;
    normsse = sse./msse;
    ISIcont = CELL(ee).ISIcontamination;
    
    onoff(:,1) = CELL(ee).onTime/3600;
    onoff(:,2) = CELL(ee).offTime/3600;
    
%     
%     crashStart  = find(isnan(sse),1,'first') - 1;
%     crashEnd    = find(isnan(sse),1,'last') + 1;
%     
%     mwfBeforeCrash  = mean(spikes(:,1:min(crashStart,end)),2);
%     mwfAfterCrash   = mean(spikes(:,min(crashEnd,end):end),2);
%     
%     % DOES MEAN WAVEFORM CHANGE AFTER RIG CRASH?
    figure(1), hold on
%     plot(mwfBeforeCrash,'--','linewidth',4);
%     plot(mwfAfterCrash,':','linewidth',4);
    plot(nanmean(spikes,2),'k','linewidth',2);
    
    uiwait(gcf);
    
    flag = input('Flag this cell? (y/n):  ','s');
    
    % DOES IT LOOK LIKE ON/OFF TIMES ARE WRONG?
    figure(2);
    set(gcf,'position',[scrsz(3)*0.2 scrsz(4)/2 scrsz(3)*0.8 scrsz(4)*0.8]);
    s1 = subplot(3,1,1); hold on;
    plot(rate,'color',[0 0.45 0.74],'linewidth',2)
    
    ylimit = get(gca,'ylim'); ylimit = ylimit(2);
    for uu = 1:size(onoff,1);
        l_a = line([onoff(uu,1) onoff(uu,1)],[0 ylimit]);
        set(l_a,'color','g','linewidth',2);
        l_b = line([onoff(uu,2) onoff(uu,2)],[0 ylimit]);
        set(l_b,'color','r','linewidth',2);
    end
    
    subplot(3,1,2); hold on;
    plot(sse,'color',[0.93 0.69 0.13],'linewidth',2)
    
    ylim2 = get(gca,'ylim'); ylim2=ylim2(2);
    set(gca,'ylim',[0 min(ylim2,0.5)]);
    drawnow;
    
    subplot(3,1,3); hold on;
    plot(ISIcont,'color',[0.64 0.08 0.18],'linewidth',2);
    
    pause(3);
    commandwindow;
    
    setnew = input('set new on/off times? (y/n):  ','s');
    
    if strcmp(setnew,'y')==1
    keepgoing = 0; count = 0;
        while keepgoing == 0
            count = count + 1;
            figure(2)
            set(l_a,'color',[0.25 0.25 0.25],'linewidth',2);
            set(l_b,'color',[0.75 0.75 0.75],'linewidth',2);
            
            set(gcf,'name','ON TIME','numbertitle','off');
            drawnow;
            [newontime(count),~] = ginput(1);
            set(gcf,'name','OFF TIME','numbertitle','off');
            drawnow;
            [newofftime(count),~] = ginput(1);
            
            set(gcf,'name','DONE!','numbertitle','off');
            subplot(s1);
            ylims = get(gca,'ylim');
            l_na = line([newontime(count) newontime(count)],[0 ylims(2)]);
            set(l_na,'color','g','linewidth',2);
            l_nb = line([newofftime(count) newofftime(count)],[0 ylims(2)]);
            set(l_nb,'color','r','linewidth',2);
           
            drawnow;
            pause(1);
            commandwindow;
            keepgoing = input('Set more times (0) or next cell (1) ?  ');
        end
    end
    
    close(gcf);
    
    if strcmp(flag,'y')==1
        RVW1_CONTCELL(cc).flag = 1;
    else
        RVW1_CONTCELL(cc).flag = 0;
    end
    
    if strcmp(setnew,'y')==1
        RVW1_CONTCELL(cc).newontime = newontime*3600;
        RVW1_CONTCELL(cc).newofftime = newofftime*3600;
    else
        RVW1_CONTCELL(cc).newontime = CELL(ee).onTime;
        RVW1_CONTCELL(cc).newofftime = CELL(ee).offTime;
    end
    
    clear newontime newofftime
end
