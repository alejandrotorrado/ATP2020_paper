% ANALYZE BEHAVIORAL STATE FOR DECLINE IN FR AFTER MD
%
% Juliet Bottorff - 2019
%
% August 2019: NEW CHANGES TO ANALYZE ONLY LIGHT PERIODS DURING EITHER MD 1 OR MD 2
% DEPENDING ON WHICH GROUP CELLS FALL INTO (CELLS THAT START DROP IN FIRST
% 24 HOURS AFTER MD AND CELLS THAT START DROP IN SECOND 24 HOURS AFTER MD) 
% ... ANALYSIS IS ONLY PERFORMED ON 12 HOURS OF LIGHT NOT ON ACTUAL
% IDENTIFIED DROP
%
% Getting rid of light/dark analysis because only analyzing light periods
%
% STEP TWO: 
% MUST RUN MD_FR_CHANGES_V4.M FIRST; INPUTS ANIMAL_STRUCTURE AND MDCELLS
% LOAD ALEJANDRO'S DATA FROM MD EXPMTS WITH STATE DATA: MDCELLS_STATES.MAT
% (CALLED CONTCELL)
%
% ALIGNS DECLINE TIME FROM EACH CELL/ANIMAL WITH ANIMALS BEHAVIORAL STATE
% AT THAT TIME
%
% OUTPUT DOWN_STATES STRUCTURE FOR USE IN MAKING GRAPHS AND PLOTS IN
% MD_STATE_TESTS.M OR MDdropFigureCode.m


clearvars -except MDCELLS animal_structure CONTCELL *CELL* *LFP*
%{
if ~exist('CONTCELL','var')
    load ('MDCELLS_STATES_JB.mat')
end
%}


% ANIMAL STRUCTURE FIELDS FROM INPUT:
%{
animal_structure.identifier
animal_structure.animal
animal_structure.ONspikes
animal_structure.normFR
animal_structure.FRslope
animal_structure.slopeThreshold
animal_structure.kde
animal_structure.decline_bins in 15 min bins from start of experiment NOT
from 7:30 am first day!!
animal_structure.MD
%}


BINTIME = animal_structure(1).bintimes(2); % Bintime carried over for conversions between FR slope sliding window array/ decline bins etc and seconds
normFRbinz = animal_structure(1).bintimes(1);


% Now analyze decline bins with statetimes
% first column in CONTCELL.VIDEO is state code, second column is time stamp
% in unix time; timestamps/states are organized by animal
AnimalNames= {'KH45', 'KH47', 'KH49', 'KH36', 'KH40', 'KH50', 'KH42', 'KH67'};

CONTCELL_NAMES = {'KH45', 'KH47', 'KH49', 'KH36', 'KH40', 'KH50', 'KH42', 'KH67'};
%State codes:
% 1 = REM
% 2 = NREM
% 3 = Gen Wake ** NOT USED **
% 4 = Active Wake
% 5 = Quiet Wake


%pull out all FR/and FR slopes in different states from each animal
all_animals = {animal_structure.animal};

random_control=0;

% initialize output structure
down_states.MD=0;

down_states.FR_slopes_sw=[];

down_states.sleepdense_sw=[];
down_states.wakedense_sw=[];
down_states.lightdense_sw=[];
down_states.darkdense_sw=[];

down_states.totalpercents=[];
down_states.total_diffFR=[];
first_states=[];
first_states_control=[];

dark_FR=[];
dark_FR_control=[];
light_FR=[];
light_FR_control=[];

AW_FR=[];
AW_FR_control=[];
QW_FR=[];
QW_FR_control=[];
REM_FR=[];
REM_FR_control=[];
NREM_FR=[];
NREM_FR_control=[];


MD_drop_light=[]; % 2 columns, first is seconds after most recent lights on that drop started, second is seconds after most recent lights on that drop ended
control_drop_light=[];

VIDEODATA = struct;

% 12/12 light dark cycle
lightson=0:(24*3600):(3600*24*8);
lightsoff=12*3600:(24*3600):(3600*24*8.5);

% clean up state time data
for names=1:length(AnimalNames)
    this_name = char(AnimalNames(names));
   
    viddata = CONTCELL.VIDEO.(this_name);
 
     
    %get rid of duplicate states
            
            testD   = diff(viddata(:,1));
            killem  = find(testD == 0);
            killem  = killem+1;
            viddata(killem,:) = [];
            if viddata(end,1) == 3
                viddata(end,:) = []; 
            end

            threes = viddata(:,1) == 3; % find the general wake code
            viddata(threes,:) = []; % just delete those entries. or see next line.
    
            % get rid of repeat button presses - second pass, clean up after converting
            % threes to subsequent value
            testD   = diff(viddata(:,1));
            killem  = find(testD == 0);
            killem  = killem+1;
            viddata(killem,:) = [];
            
            % find and eliminate transient wake episodes (e.g. the animal wakes up
            % simply to shift position)
            nsec=60; %set limit at 1 min-eliminate all states that only last 1 minute or less
            sdur = diff(viddata(:,2));
            fcks = find(viddata(1:end-1,1)>2 & sdur<=nsec ); % find instances of waking that are less than nsec seconds long

            tempf = viddata(fcks+1,1)<3; % find one fcks for which the following state is still sleeping
            fcks2 = fcks(tempf);

            viddata(fcks2,:) = [];

            % get rid of repeat button presses - third pass, clean up after deleting
            % short wake periods
            testD   = diff(viddata(:,1));
            killem  = find(testD == 0);
            killem  = killem+1;
            viddata(killem,:) = [];
            
            
      VIDEODATA.(this_name) = viddata;
            
end

cell_count=0;
% loop through each animal

for i=1:length(AnimalNames)
   
    current_animal_cells_index=ismember(all_animals, AnimalNames(i)); % index of all cells from current animal
    current_animal=find(current_animal_cells_index==1); %all the cells for the current animal being analyzed
    if isempty (current_animal)
        continue
    end
    %loop through cells for the animal
    for j=1:length(current_animal)
        cell_count=cell_count+1;
        ID = current_animal(j);
        
        % get ON spikes from  current cell
        current_spikes=animal_structure(ID).ONspikes;
        
        % get behavioral data from current animal
        videodata = VIDEODATA.(char(AnimalNames(i)));
     
            % get decline time start and end from deprived animals
            secfromstart = MDCELLS(animal_structure(ID).identifier).trem;
            
            
            if animal_structure(ID).MD == 1
                ID_drop_start=(animal_structure(ID).decline_bins(1))*BINTIME; %decline start in seconds
                ID_drop_end=(animal_structure(ID).decline_bins(2))*BINTIME; %decline end in seconds
                if ID_drop_start < (3.5*24*3600)
                    decline_time_starts = 3*24*3600;
                    decline_time_ends = 3.5*24*3600;
                else
                    decline_time_starts = 4*24*3600;
                    decline_time_ends = 4.5*24*3600;
                end
                these_off_times = MDCELLS(animal_structure(ID).identifier).offTime;
                for o = 1:length(these_off_times)
                    if these_off_times(o) > decline_time_starts && these_off_times(o) < decline_time_ends
                        disp ('OFFTIME DURING DROP ANALYSIS!')
                        disp (animal_structure(ID).identifier)
                        if decline_time_ends == 4.5*3600*24
                            disp ('second day')
                        end
                        disp (decline_time_ends-these_off_times(o))
                        decline_time_ends = these_off_times(o);
                    end

                end
            else
                these_off_times = MDCELLS(animal_structure(ID).identifier).offTime;
                decline_time_starts = [3*24*3600, 4*24*3600];
                decline_time_ends = [3.5*24*3600, 4.5*24*3600];
                
                for d = 1:2
                for o = 1:length(these_off_times)
                    if these_off_times(o) > decline_time_starts(d) && these_off_times(o) < decline_time_ends(d)
                        disp ('OFFTIME DURING DROP ANALYSIS!')
                        disp (animal_structure(ID).identifier)
                        disp ('control cell')
                        decline_time_ends(d) = these_off_times(o);
                    end

                end
                end
            end
            
            
        for dd = 1:length(decline_time_starts)
            
            decline_time_start = decline_time_starts(dd);
            decline_time_end = decline_time_ends(dd);
            
            decline_length=decline_time_end-decline_time_start;
            current_down_spikes=current_spikes(current_spikes>=decline_time_start); % pull out only current cell spikes during decline time
            current_down_spikes=current_down_spikes(current_down_spikes<=decline_time_end);
            states=videodata(:,1);
            mytime=videodata(:,2);
            
            decline_time_start_index=find(mytime>=decline_time_start,1);
            decline_time_end_index = find(mytime<=decline_time_end,1, 'last');
            decline_time=mytime(mytime>=decline_time_start); % state time data only in decline time to be analyzed
            decline_time=decline_time(decline_time<=decline_time_end);
            %find last state before the decline bin
                first_decline_time=mytime(mytime<=decline_time_start);
                first_decline_time=first_decline_time(end);
                first_state=states(find(mytime==first_decline_time));
                if animal_structure(ID).MD==1
                    first_states(end+1)=first_state;
                else 
                    first_states_control(end+1)=first_state;
                end
                
            % all (including first) state codes and time stamps to be
            % analyzed during decline
            decline_statetimes=[];
            decline_statetimes(:,2)=[decline_time_start; decline_time];
            decline_statetimes(:,1)=[first_state; states(decline_time_start_index:decline_time_end_index)];
            firstHourEnd = find(decline_statetimes(:,2)<(decline_time_start+3600), 1, 'last');
            down_states(cell_count).firstStates = decline_statetimes(1:firstHourEnd,:);
            %sliding window slope analysis: get time period and slope of FR during decline from
            %animal_structure
            sw_size = 2; % enter sliding window size to analyze FR slopes in hours
            sw_size = 3600*sw_size/BINTIME;
            x_sw=animal_structure(ID).FRslope(1,floor((decline_time_start)/BINTIME):floor((decline_time_end)/BINTIME)+sw_size);
            x_sw=x_sw.*(3600*24); %put this back in seconds for reference for statetimes
            decline_slope_sw=animal_structure(ID).FRslope(2,floor((decline_time_start)/BINTIME):floor((decline_time_end)/BINTIME));
            
            % find total change in FR over whole decline bin (first 20% to
            % last 20% in format (B-A)/(B+A)
            total_edges=linspace(decline_time_start,decline_time_end);

            [spikecounts, spikeedges]=histcounts(current_down_spikes,total_edges);
            first20percent=sum(spikecounts(1:20));
            %baseline=first20percent/(spikeedges(21)-spikeedges(1)); %FR for first 20% of segment in Hz
            
            baseline_spikes = current_spikes(current_spikes>=(decline_time_start-(10*3600)));
            baseline_spikes = baseline_spikes(baseline_spikes<decline_time_start);
            baseline = length(baseline_spikes)/(10*3600);
            
            last20percent=sum(spikecounts(79:end)); %there's only 99 bins (100 'edge' numbers)
            lastFR=last20percent/(spikeedges(end)-spikeedges(79));
         
            total_diffFR=(lastFR-baseline)/(lastFR+baseline);
            if isempty(down_states(cell_count).total_diffFR)
                down_states(cell_count).total_diffFR=total_diffFR;
            else
                down_states(cell_count).total_diffFR = mean([down_states(cell_count).total_diffFR, total_diffFR]);
            end
            down_states(cell_count).identifier = animal_structure(ID).identifier;
  
                % NOW DO SLIDING WINDOW AGAIN BUT FOR FINDING STATE-DENSE
                % CHUNKS, SO NEED TO BE ABLE TO CHANGE SIZE OF SLIDING WINDOW
                % AND WHEN STATE-DENSE SECTION IS FOUND, NEED TO SKIP TO THE
                % END OF THAT WHOLE SECTION TO KEEP SEARCHING AND NOT DOUBLE
                % COUNT

                % initialize matrices of percentages of current
                % window in each state to find state dense chunks
                percentREM_temp=NaN(1, length(decline_slope_sw));
                percentNREM_temp=NaN(1, length(decline_slope_sw));
                percentAW_temp=NaN(1, length(decline_slope_sw));
                percentQW_temp=NaN(1, length(decline_slope_sw));
             
                % initialize vectors to hold change in FR during each state
                % dense period ((B-A)/(B+A))
                REM_dense=[];
                NREM_dense=[];
                AW_dense=[];
                QW_dense=[];
                sleep_dense=[];
                wake_dense=[];
                

                w=7; %window size: number of (15 min) bins to include in each state dense window (minus 1 because for 8 bins, needs to be 8-7 to get to first bin)
                if length(decline_slope_sw)>(w+1)
                    h=((w+1)/4); % number of hours that the window covers
                    REM=0;
                    NREM=0;
                    AW=0;
                    QW=0;
                    sleep=0;
                    wake=0;
                    
               
                    % Loop through windows looking for state dense periods!
                    for m=(w+1):length(decline_slope_sw) 

                        sw_statetimes=mytime(mytime>=x_sw(m-w));
                        sw_statetimes=sw_statetimes(sw_statetimes<=x_sw(m));
                        sw_states=states(ismember(videodata(:,2),sw_statetimes));
                        REM_sw=find(sw_states==1);
                        NREM_sw=find(sw_states==2);
                        AW_sw=find(sw_states==4);
                        QW_sw=find(sw_states==5);

               
                        % get percent current window spent in each state
                        REM_time=0;
                        for R=1:length(REM_sw)
                            s=sw_statetimes(REM_sw(R));
                            if sw_statetimes(REM_sw(R))==sw_statetimes(end)
                                e=x_sw(m);
                            else
                                e=sw_statetimes(REM_sw(R)+1);
                            end
                            REM_time=REM_time+(e-s);
                        end
                        percentREM_temp(m-w)=REM_time/(x_sw(m)-x_sw(m-w));

                        NREM_time=0;
                        for R=1:length(NREM_sw)
                            s=sw_statetimes(NREM_sw(R));
                            if sw_statetimes(NREM_sw(R))==sw_statetimes(end)
                                e=x_sw(m);
                            else
                                e=sw_statetimes(NREM_sw(R)+1);
                            end
                            NREM_time=NREM_time+(e-s);
                        end
                        percentNREM_temp(m-w)=NREM_time/(x_sw(m)-x_sw(m-w));

                        AW_time=0;
                        for R=1:length(AW_sw)
                            s=sw_statetimes(AW_sw(R));
                            if sw_statetimes(AW_sw(R))==sw_statetimes(end)
                                e=x_sw(m);
                            else
                                e=sw_statetimes(AW_sw(R)+1);
                            end
                            AW_time=AW_time+(e-s);
                        end
                        percentAW_temp(m-w)=AW_time/(x_sw(m)-x_sw(m-w));

                        QW_time=0;
                        for R=1:length(QW_sw)
                            s=sw_statetimes(QW_sw(R));
                            if sw_statetimes(QW_sw(R))==sw_statetimes(end)
                                e=x_sw(m);
                            else
                                e=sw_statetimes(QW_sw(R)+1);
                            end
                            QW_time=QW_time+(e-s);
                        end
                        percentQW_temp(m-w)=QW_time/(x_sw(m)-x_sw(m-w));

                        % look for state dense periods
                        percent=.70; %percent of time to consider state dense
                        % if find a state dense period, need to wait til
                        % the END of that window to be able to find another
                        if QW>0
                            QW=QW+1;
                            if QW==(w+2)
                                QW=0;
                            end
                        end
                        
                        if AW>0
                            AW=AW+1;
                            if AW==(w+2)
                                AW=0;
                            end
                        end
                        
                        if REM>0
                            REM=REM+1;
                            if REM==(w+2)
                                REM=0;
                            end
                        end
                       
                        if NREM>0
                            NREM=NREM+1;
                            if NREM==(w+2)
                                NREM=0;
                            end
                        end
                        
                        if sleep>0
                            sleep=sleep+1;
                            if sleep==(w+2)
                                sleep=0;
                            end
                        end
                        
                        if wake>0
                            wake=wake+1;
                            if wake==(w+2)
                                wake=0;
                            end
                        end
                        

                        
                        % record change in FR during any found state dense
                        % periods
                        if percentQW_temp(m-w)>=percent
                            if QW==0
                                QW=QW+1;
                                %calculate FR change(diffFR) for this 4 hour block
                                spikes=current_spikes(current_spikes>=x_sw(m-w));
                                spikes=spikes(spikes<=x_sw(m));
                                edges=linspace(x_sw(m-w), x_sw(m));
                                [spikecounts, spikeedges]=histcounts(spikes, edges);
                                first40percent=sum(spikecounts(1:40));
                                baseline=first40percent/(spikeedges(41)-spikeedges(1)); %FR for first 20% of segment in Hz

                                last40percent=sum(spikecounts(59:end)); %there's only 99 bins (100 'edge' numbers)
                                lastFR=last40percent/(spikeedges(end)-spikeedges(59));
                                diffFR=(lastFR-baseline)/(lastFR+baseline);
                                QW_dense(end+1)=diffFR;
                            end
                        elseif percentAW_temp(m-w)>=percent
                            if AW==0
                                AW=AW+1;
                                spikes=current_spikes(current_spikes>=x_sw(m-w));
                                spikes=spikes(spikes<=x_sw(m));
                                edges=linspace(x_sw(m-w), x_sw(m));
                                [spikecounts, spikeedges]=histcounts(spikes, edges);
                                first40percent=sum(spikecounts(1:40));
                                baseline=first40percent/(spikeedges(41)-spikeedges(1)); %FR for first 20% of segment in Hz

                                last40percent=sum(spikecounts(59:end)); %there's only 99 bins (100 'edge' numbers)
                                lastFR=last40percent/(spikeedges(end)-spikeedges(59));
                                diffFR=(lastFR-baseline)/(lastFR+baseline);
                                AW_dense(end+1)=diffFR;
                            end
                        elseif percentNREM_temp(m-w)>=percent
                            if NREM==0
                                NREM=NREM+1;
                                spikes=current_spikes(current_spikes>=x_sw(m-w));
                                spikes=spikes(spikes<=x_sw(m));
                                edges=linspace(x_sw(m-w), x_sw(m));
                                [spikecounts, spikeedges]=histcounts(spikes, edges);
                                first40percent=sum(spikecounts(1:40));
                                baseline=first40percent/(spikeedges(41)-spikeedges(1)); %FR for first 20% of segment in Hz

                                last40percent=sum(spikecounts(59:end)); %there's only 99 bins (100 'edge' numbers)
                                lastFR=last40percent/(spikeedges(end)-spikeedges(59));
                                diffFR=(lastFR-baseline)/(lastFR+baseline);
                                NREM_dense(end+1)=diffFR;
                            end
                        elseif percentREM_temp(m-w)>=percent
                            if REM==0
                                REM=REM+1;
                                spikes=current_spikes(current_spikes>=x_sw(m-w));
                                spikes=spikes(spikes<=x_sw(m));
                                edges=linspace(x_sw(m-w), x_sw(m));
                                [spikecounts, spikeedges]=histcounts(spikes, edges);
                                first40percent=sum(spikecounts(1:40));
                                baseline=first40percent/(spikeedges(41)-spikeedges(1)); %FR for first 20% of segment in Hz

                                last40percent=sum(spikecounts(59:end)); %there's only 99 bins (100 'edge' numbers)
                                lastFR=last40percent/(spikeedges(end)-spikeedges(59));
                                diffFR=(lastFR-baseline)/(lastFR+baseline);
                                REM_dense(end+1)=diffFR;
                            end
                        end
                       
                        if (percentREM_temp(m-w)+percentNREM_temp(m-w))>=percent
                            if sleep==0
                                sleep=sleep+1;
                                spikes=current_spikes(current_spikes>=x_sw(m-w));
                                spikes=spikes(spikes<=x_sw(m));
                                edges=linspace(x_sw(m-w), x_sw(m));
                                [spikecounts, spikeedges]=histcounts(spikes, edges);
                                %first40percent=sum(spikecounts(1:40));
                                first20percent = sum(spikecounts(1:20));
                                baseline=first20percent/(spikeedges(21)-spikeedges(1));
                                %baseline=first40percent/(spikeedges(41)-spikeedges(1)); %FR for first 20% of segment in Hz

                                %last40percent=sum(spikecounts(59:end)); %there's only 99 bins (100 'edge' numbers)
                                %lastFR=last40percent/(spikeedges(end)-spikeedges(59));
                                last20percent = sum(spikecounts(79:end));
                                lastFR = last20percent/(spikeedges(end)-spikeedges(79));
                                diffFR=(lastFR-baseline)/(lastFR+baseline);
                                sleep_dense(end+1)=diffFR;
                            end
                        elseif (percentAW_temp(m-w)+percentQW_temp(m-w))>=percent
                            if wake==0
                                wake=wake+1;
                                spikes=current_spikes(current_spikes>=x_sw(m-w));
                                spikes=spikes(spikes<=x_sw(m));
                                edges=linspace(x_sw(m-w), x_sw(m));
                                [spikecounts, spikeedges]=histcounts(spikes, edges);
                                %first40percent=sum(spikecounts(1:40));
                                first20percent = sum(spikecounts(1:20));
                                baseline = first20percent/(spikeedges(21)-spikeedges(1));
                                %baseline=first40percent/(spikeedges(41)-spikeedges(1)); %FR for first 20% of segment in Hz

                                %last40percent=sum(spikecounts(59:end)); %there's only 99 bins (100 'edge' numbers)
                                %lastFR=last40percent/(spikeedges(end)-spikeedges(59));
                                last20percent = sum(spikecounts(79:end));
                                lastFR = last20percent/(spikeedges(end)-spikeedges(79));
                                diffFR=(lastFR-baseline)/(lastFR+baseline);
                                wake_dense(end+1)=diffFR;
                            end
                        end

                    end
                end
                % record changes in FR during state dense periods during
                % decline ((B-A)/(B+A))

                down_states(cell_count).wakedense_sw=[down_states(cell_count).wakedense_sw, wake_dense];
                down_states(cell_count).sleepdense_sw=[down_states(cell_count).sleepdense_sw, sleep_dense];


                REM_d=find(decline_statetimes(:,1)==1); %indices of REM start DURING DECLINE TIME
                NREM_d=find(decline_statetimes(:,1)==2); %indices of NREM start  DURING DECLINE TIME
                ActiveWake_d=find(decline_statetimes(:,1)==4); %indices of Active Wake start DURING DECLINE TIME
                QuietWake_d=find(decline_statetimes(:,1)==5); %indices of Quiet Wake start DURING DECLINE TIME
                
                if isempty(REM_d)
                    disp(ID)
                    disp('no decline REM')
                end
                if isempty(NREM_d)
                    disp(ID)
                    disp('no decline NREM')
                end
                if isempty(ActiveWake_d)
                    disp(ID)
                    disp('no decline Active Wake')
                end
                if isempty(QuietWake_d)
                    disp(ID)
                    disp('no decline Quiet Wake')
                end

                REM_dvec.startstop=NaN(length(REM_d), 2); %start and stop times (column 1 and 2) of each REM segment during decline
                
                NREM_dvec.startstop=NaN(length(NREM_d),2); % same but NREM
            
                ActiveWake_dvec.startstop=NaN(length(ActiveWake_d),2); %same but Active Wake
               
                QuietWake_dvec.startstop=NaN(length(QuietWake_d),2); %same but Quiet Wake
             
                % identify start stop times/spikes and average FR for each
                % epoch of each behavioral state
                for k=1:length(REM_d)
                    REM_dvec.startstop(k,1)=decline_statetimes(REM_d(k),2);

                    if decline_statetimes(REM_d(k),2)==decline_statetimes(end,2)
                        REM_dvec.startstop(k,2)=decline_time_end;
                    else
                        REM_dvec.startstop(k, 2)=decline_statetimes(REM_d(k)+1,2);
                    end
                end
                for n=1:length(NREM_d)
                    NREM_dvec.startstop(n,1)=decline_statetimes(NREM_d(n),2);
                    if decline_statetimes(NREM_d(n),2)==decline_statetimes(end,2)
                        NREM_dvec.startstop(n,2)=decline_time_end;
                    else
                        NREM_dvec.startstop(n,2)=decline_statetimes(NREM_d(n)+1,2);
                    end
            
                end
                for a=1:length(ActiveWake_d)
                    ActiveWake_dvec.startstop(a,1)=decline_statetimes(ActiveWake_d(a),2);
                    if decline_statetimes(ActiveWake_d(a),2)==decline_statetimes(end,2)
                        ActiveWake_dvec.startstop(a,2)=decline_time_end;
                    else
                        ActiveWake_dvec.startstop(a,2)=decline_statetimes(ActiveWake_d(a)+1,2);
                    end
               
                end
                for q=1:length(QuietWake_d)
                    QuietWake_dvec.startstop(q,1)=decline_statetimes(QuietWake_d(q),2);
                    if decline_statetimes(QuietWake_d(q),2)==decline_statetimes(end,2)
                        QuietWake_dvec.startstop(q,2)=decline_time_end;
                    else
                        QuietWake_dvec.startstop(q,2)=decline_statetimes(QuietWake_d(q)+1,2);
                    end
                
                end
                
                
                 % find percent of whole decline period spent in each state
                 total_REM_percent=0;
                 total_NREM_percent=0;
                 total_AW_percent=0;
                 total_QW_percent=0;
                 for r=1:length(REM_dvec.startstop(:,1))
                     t=REM_dvec.startstop(r,2)-REM_dvec.startstop(r,1);
                     total_REM_percent=total_REM_percent+t;
                 end
                 total_REM_percent=total_REM_percent/decline_length;
                 for n=1:length(NREM_dvec.startstop(:,1))
                     t=NREM_dvec.startstop(n,2)-NREM_dvec.startstop(n,1);
                     total_NREM_percent=total_NREM_percent+t;
                 end
                 total_NREM_percent=total_NREM_percent/decline_length;
                 for a=1:length(ActiveWake_dvec.startstop(:,1))
                     t=ActiveWake_dvec.startstop(a,2)-ActiveWake_dvec.startstop(a,1);
                     total_AW_percent=total_AW_percent+t;
                 end
                 total_AW_percent=total_AW_percent/decline_length;
                 for q=1:length(QuietWake_dvec.startstop(:,1))
                     t=QuietWake_dvec.startstop(q,2)-QuietWake_dvec.startstop(q,1);
                     total_QW_percent=total_QW_percent+t;
                 end
                 total_QW_percent=total_QW_percent/decline_length;
                 if isempty(down_states(cell_count).totalpercents)
                    down_states(cell_count).totalpercents=[total_REM_percent, total_NREM_percent, total_AW_percent, total_QW_percent];
                 else
                     thisperiod_percents = [total_REM_percent, total_NREM_percent, total_AW_percent, total_QW_percent];
                     down_states(cell_count).totalpercents=mean([down_states(cell_count).totalpercents;thisperiod_percents]);
                 end

            % record deprived status of each cell
            if animal_structure(ID).MD==1
                down_states(cell_count).MD=1;
            else
                down_states(cell_count).MD=0;
            end
        end
    end
end
