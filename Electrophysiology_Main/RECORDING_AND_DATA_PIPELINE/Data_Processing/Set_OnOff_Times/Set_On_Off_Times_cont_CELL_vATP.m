function Set_On_Off_Times_cont_CELL_vATP
animal = input('What is the animal code?','s');

[lFile, lDir] = uigetfile(cd,'Pick your CELL file.');
load([lDir lFile]);

for ee = 1:length(CELL);
    CELL(ee).animal = animal;
end

for ee = 1:length(CELL);
    
    CELL(ee).onTime     = [];
    CELL(ee).offTime    = [];
    
end

g_times = input('Are there global black out times for this animal? y or n  ','s');

if g_times == 'y';
    bot = 1; count = 0;
    while bot == 1;
        
        % ADDING different hemispheres...
        hemisphere = input('which hemisphere? L, R, or (B)oth ','s');
        
        count = count+1;
        ee = input('Which cell number?');
        
%         tdiff = CELL(ee).EXPTSTART  - (27000 + (3600*4));
        tdiff = CELL(ee).EXPTSTART;
        trem = rem(tdiff,(24*3600));
        
        A = histc(CELL(ee).time, 0 : 300 : CELL(ee).maxTime+300 );
        x = 1: 300 : CELL(ee).maxTime + 300;
        B = A/300;
        
        cellfig = figure;
        set(cellfig,'units','normalized','position',[0.01 0.25 0.98 0.5]);
        
        plot(x, B), hold on
        set(gca,'xlim',[-20000 CELL(ee).maxTime+20000]);
        
        d_n         = [];
        d_n(:,1)    = 0 : (3600*24) : (ceil(CELL(1).time(end)/(3600*24))* (3600*24) );
        d_n(:,2)    = (3600*12) : (3600*24) : (ceil(CELL(ee).time(end)/(3600*24))* (3600*24) + (3600*12));
        
        
        l_d = [];
        for rr = 1:length(d_n) - 1;
            
            l_d ( d_n(rr,1)+1 : d_n(rr,2)   ) = zeros;
            l_d ( d_n(rr,2)+1 : d_n(rr+1,1) ) = ones;
            
        end
        
        l_d(1:trem) = [];
        
        ylimit = get(gca,'ylim');
        ylimit = ylimit(2);
        l_d = l_d * ylimit;
        
        xs = 1:length(l_d);
        
        
        area( xs,l_d, 'facecolor',[0.85 0.85 0.85]);
        
        set(gca,'children',flipud(get(gca,'children')));
        
        if CELL(ee).channel<17;
            hem = 'L';
        else
            hem = 'R';
        end
        
        title([hem ' hem, deprived ' num2str(CELL(ee).deprived) ', Qual. ' num2str(CELL(ee).quality) '.']);
        
        ht = text(0, ylimit/2, 'OFF TIME','fontsize',20);
        [gx_start(count), ~]    = ginput(1);
        delete (ht);
        
        text(0, ylimit/2, 'RESUME TIME','fontsize',20);
        [gx_end(count), ~]      = ginput(1);
        
        if gx_start(count)<0;
            gx_start(count) = 0;
        end
        
        switch hemisphere;
            case ('L')
                gx_start(count,2)   = 1;
                gx_end (count,2)    = 1;
            case ('l')
                gx_start(count,2)   = 1;
                gx_end (count,2)    = 1;
            case ('R')
                gx_start(count,2)   = 2;
                gx_end (count,2)    = 2;
            case ('r')
                gx_start(count,2)   = 2;
                gx_end (count,2)    = 2;
            case('B')
                gx_start(count,2)   = 3;
                gx_end (count,2)    = 3;
            case('b')
                gx_start(count,2)   = 3;
                gx_end (count,2)    = 3;
        end
        
        close all;
        
        bot = input('Are there more blackout times? Y(1) or N(0) ');
        
    end
    
    
end



for ee = 1:length(CELL);
    
    more = 'y';
    
    count = 0;
    while more == 'y';
        
        % trem is the number of seconds between lights on and the
        % experiment start time
        count = count+1;
        
        CELL(ee) = mainBusiness_ATP(CELL(ee), count, ee);
        
        close all;
        if CELL(ee).quality<4;
            more = input('Mark more times (y) or move to the next cell (n)?','s');
        else
            more ='n';
        end
    end
    %     if ~ispc;
    %         jheapcl
    %     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%MAKE CELL TIME CHECKER  - - - SOMETIMES ON AND OFF ARE ERRORED TOO CLOSE TO ONE ANOTHER
disp('Checking for errors - will present data if it needs correction.');
keyboard
% NOTE - while you're here, FINALLY fix the code in the following loop!
for yy = 1:length(CELL);
   
    bbb = CELL(yy).offTime - CELL(yy).onTime;
    bbb = bbb/3600;
    disp(bbb);
    if min(bbb)< 5;

        CELL(yy).offTime    = [];
        CELL(yy).onTime     = [];
        
        [CELL] = mainBusiness(CELL(yy), 1, yy);
        
    end
end

if g_times == 'y';
    keyboard
    % You need to fix this. The code isn't doing what it should for putting
    % global blackout times in place. figure out what it's doing and why
    % it's failing. maybe there's just a cleaner way to do it.
    for uu = 1:length(CELL);
        for gg = 1:size(gx_end,1);
            
            if (min(CELL(uu).onTime) < gx_end(gg,1));
                    
                    if gx_end(gg,2) == 1 && CELL(uu).channel<17;
                        CELL(uu).onTime     = [CELL(uu).onTime; gx_end(gg,1)];
                        CELL(uu).offTime    = [CELL(uu).offTime; gx_start(gg,1)];
                        
                        CELL(uu).onTime     = sort(CELL(uu).onTime);
                        CELL(uu).offTime    = sort(CELL(uu).offTime);
                    end

                if gx_end(gg,2) == 2 && CELL(uu).channel>16;
                    CELL(uu).onTime     = [CELL(uu).onTime; gx_end(gg,1)];
                    CELL(uu).offTime    = [CELL(uu).offTime; gx_start(gg,1)];
                    
                    CELL(uu).onTime     = sort(CELL(uu).onTime);
                    CELL(uu).offTime    = sort(CELL(uu).offTime);
                end
                
                if gx_end(gg,2) == 3;
                    CELL(uu).onTime     = [CELL(uu).onTime; gx_end(gg,1)];
                    CELL(uu).offTime    = [CELL(uu).offTime; gx_start(gg,1)];
                    
                    CELL(uu).onTime     = sort(CELL(uu).onTime);
                    CELL(uu).offTime    = sort(CELL(uu).offTime);
                end
                
                
                
            end
            
        end
        
    end
    

    
end


for ee = 1:length(CELL);
    
    if length(CELL(ee).onTime) ~= length(CELL(ee).offTime);
        keyboard;
    end
    nspk = []; totT = [];
    
    if CELL(ee).onTime(1)/3600 < 6;
        CELL(ee).onTime(1) = 0;
    end
    
    if  (CELL(ee).time(end) - CELL(ee).offTime(end))/3600 < 8;
        CELL(ee).offTime(end) = CELL(ee).time(end);
    end
    
    for ii = 1:length(CELL(ee).onTime);
        
        nspk(ii) = sum(CELL(ee).time > CELL(ee).onTime(ii) & CELL(ee).time< CELL(ee).offTime(ii));
        totT(ii) = CELL(ee).offTime(ii) - CELL(ee).onTime(ii);
        
    end
    
    if sum(totT) < (0.7 * CELL(ee).time(end));
        CELL(ee).cont_stat = 0;
    else
        CELL(ee).cont_stat = 1;
    end
    CELL(ee).rate_CELL = sum(nspk)/sum(totT);
    
end


if isfield(CELL,'meanRate');
    CELL = rmfield(CELL,'meanRate');
end

if isfield(CELL,'ctr');
    CELL = rmfield(CELL,'ctr');
end

CELL = orderfields(CELL);
disp(['Animal is ' CELL(1).animal]);

[sFile sDir] = uiputfile(cd,'Save your updated CELL file.');

save([sDir sFile],'CELL','-v7.3');





  

end





