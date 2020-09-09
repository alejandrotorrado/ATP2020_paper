% ANALYZE BEHAVIORAL STATE FOR DECLINE IN FR AFTER MD
% 
% Juliet Bottorff - 2019
%
% STEP ONE: 
% LOAD ALEJANDRO'S DATA FROM MD EXPMTS: CALL THEM MDCELLS
%
% UPDATE TO COMBINE ORIGINAL CELLS WITH FULL BASELINE AND CELLS FROM
% SHORTENED EXPERIMENTS ; ALIGN TO START OF MD
%
% FILTERS FOR ONLY QUALITY 1,2 CELLS
% FILTERS FOR 'ON' CELLS FOR WHOLE SESSION
% FILTERS FOR RSU'S (kmeans clustering to define parameters)
%
% PLOTS EACH RSU (CONTROL AND DEPRIVED) ALONG WITH IDENTIFIED 'DECLINE BIN'
% BASED ON KERNEL DENSITY ESTIMATE AND SLOPES AND DEFINED BINS
%
% PLOTS AVERAGE OF CONTROL AND DEPRIVED WITH SHADED ERROR BARS OVER WHOLE 8
% DAY PERIOD
%
% OUTPUT ANIMAL_STRUCTURE WITH DECLINE BIN/FR INFO FOR ALL ON RSU'S USED IN
% ANALYSIS



clearvars -except MDCELLS sometjoing JMDonly *CELL* *LFP* animal_structure_ER
%{
if ~exist('MDCELLS','var')
    load ('Alejandro''s_MD_cells.mat')
end

clc
%}

%% filter to only get quality 1,2 cells
% Already done

%%

%FIRST: pick out cells that are active for >70% of total recording time to
%find continuous cells


lastcell=length([MDCELLS.quality]);
ON_cells=0; % keep track of how many continuous cells we have
ON_cells_index=NaN(1,lastcell); % index in original matrix of cells that we just found to be continuous

 
% Initialize: Rsquared values for regressions for finding decline bins
mdls_deprived=[];
mdls_control=[];
% Rsquared values for the regressions during the chosen decline bins only
Rsquared_decline_deprived=[];
Rsquared_decline_control=[];


% IF want to find decline bins using normalized FR's instead of SLOPES of
% these FR's, FRthreshold==1, ==0 otherwise
FRthreshold=0;
% CHOOSE percent of total expt time cell must be ON
ONthreshold = .7; 

myinput=input ('Are these full-length MD cells or shortened from another experiment (with shorter baseline before MD)? (Full=1, shortened=0)');


for n=1:lastcell
    this_animal=MDCELLS(n).animal;
    if myinput==0
        numDays=7-MDCELLS(n).dayStart;

    else
        numDays=8.5;
    end
    experiment_time_start=MDCELLS(n).trem; % in seconds

    experiment_time_end=experiment_time_start+(3600*24*numDays); % in seconds

    total_experiment_time=experiment_time_end-experiment_time_start;
    
        onTimes=MDCELLS(n).onTime;
        offTimes=MDCELLS(n).offTime;

        all_ontimes=offTimes(1)-onTimes(1); %length of time in seconds that cell is 'on'
        
        if length(onTimes)>1 

            for i=2:length(onTimes)

                all_ontimes=all_ontimes+(offTimes(i)-onTimes(i)); %if multiple on/off times, add all subsequent time chunks that cell is 'on'
            
            end
        end

        if all_ontimes>(ONthreshold * total_experiment_time) % analyze this cell if continuous enough
            ON_cells=ON_cells+1;
            ON_cells_index(n)=n;
        end
end

disp('ON cells found')

% now use only those ON cells for rest of the analysis
ON_cells_index=ON_cells_index(~isnan(ON_cells_index));

bintime=(12*60*60); %bin time in seconds, for baseline
bintime2=(30*60); %bin time in seconds, for calculating FR
bintime3=(15*60); %for sliding window steps

% x axis for graphs of FR over whole experiment period
full_expt_end_time=(8.5*3600*24);
firstpoint=experiment_time_start+(bintime2/2);
lastpoint=full_expt_end_time+(bintime2/2);
xaxis=firstpoint:bintime2:lastpoint;
xaxisdays=xaxis/(3600*24); %now xaxis in days not seconds

%initialize FR matrices for all cells
all_FR_normalized_RSU=NaN(length(ON_cells_index),length(xaxis));
RSU_count=0;
FR_normalized_RSU_control=NaN(length(ON_cells_index), length(xaxis));
RSU_control_count=0;

% initialize animal_structure for output info
animal_structure.identifier=[];
animal_structure.animal=[];
animal_structure.ONspikes=[];
animal_structure.normFR=[];
animal_structure.FRslope=[];
animal_structure.slopeThreshold=0;
animal_structure.posThreshold=0;
animal_structure.kde=[];
animal_structure.decline_bins=[];
animal_structure.MD=0;
animal_structure.bintimes = [bintime2, bintime3];

% loop through all continuous 'on' cells
RSU_number=0; % reference for how many RSU's are analyzed total
r=0;
remove_at_end = []; % for cells that I end up skipping because of too short baseline, need to get rid of these in animal_structure output
no_drop = [];
no_drop_ns = [];
old_animals = {'KH45', 'KH47', 'KH49', 'KH36', 'KH40', 'KH50', 'KH42'};
new_animals = {'AT14', 'AT16', 'AT29', 'KH67'};
for m=1:length(ON_cells_index)

    n=ON_cells_index(m); % index of current cell in MDCELLS structure
    
        if myinput==0
            numDays=7-MDCELLS(n).dayStart;
        else
            numDays=8.5;
        end

    
    if ismember(MDCELLS(n).animal, old_animals)
        np_thresh=10;
    else
        np_thresh=.35;
    end
    
    if MDCELLS(n).neg_pos_time>=np_thresh && MDCELLS(n).quality<3 

        RSU_number=RSU_number+1;
        animal_structure(RSU_number).identifier=n; %identifier of this cell; index in the original MDCELLS structure
        %Sort by animal too!! Necessary for state analysis!!
        animal_structure(RSU_number).animal=MDCELLS(n).animal; %identify which animal the cell is from
        animal_structure(RSU_number).channel = MDCELLS(n).channel;
        % need this in again because depends on current cell's trem
        firstpoint=0+(bintime2/2);
        full_expt_end_time= (numDays*3600*24);
        lastpoint=full_expt_end_time+(bintime2/2);
        xaxis=firstpoint:bintime2:lastpoint;
        xaxisdays=xaxis/(3600*24); %now xaxis in days not seconds
        
        sorted_time=sort(MDCELLS(n).time); % sort all spiketimes to make sure they're in order
         
        % filter out all spikes from 'off' times
        onTimes=sort(MDCELLS(n).onTime);
        offTimes=sort(MDCELLS(n).offTime);
        
        sorted_time(sorted_time<onTimes(1))=NaN;
        sorted_time(sorted_time>offTimes(end))=NaN;
        
        
        if length(onTimes)>1 
            
            for i=2:length(onTimes)
                current_working_time=sorted_time(sorted_time>offTimes(i-1));
                start_next_off=find(sorted_time==min(current_working_time));
                current_working_time=sorted_time(sorted_time<onTimes(i));
                end_next_off=find(sorted_time==max(current_working_time));
                sorted_time(start_next_off:end_next_off)=NaN;
         
            end
        end

        % spiketimes to use for analysis
        ontime_spikes=sorted_time(~isnan(sorted_time));
        animal_structure(RSU_number).ONspikes=ontime_spikes;
        

            % normalize 12 hour block #'s 3-5 (36 hours total of baseline starting after baseline day 1)
            experiment_time_start=0;
            experiment_time_end= (3600*24*numDays);
            edges=experiment_time_start:bintime:experiment_time_end;
            [bincounts, binedges]=histcounts(ontime_spikes, edges);
            
            % For MD cells taken from other experiments (some have
            % shortened baseline- total time is different: 
            % MAKE SURE basline time is long enough to trust normalization
            % general baseline 
            if myinput==0
                if MDCELLS(n).dayStart==0 % this cell has full 2.5 day baseline
                    MYbaseline_start=24*3600;
                    MYbaseline_end=3600*(24+36);
                    %baseline_FR=sum(bincounts(3:5))/(bintime*3); 
                elseif MDCELLS(n).dayStart==1 % this cell only has 1.5 day baseline max
                    MYbaseline_start=12*3600;
                    MYbaseline_end=3600*(36);
                    %baseline_FR=sum(bincounts(2:3))/(bintime*2);
                elseif MDCELLS(n).dayStart==2 % this cell has only .5 day baseline max
                    MYbaseline_start=0;
                    MYbaseline_end=12*3600;
                    %baseline_FR=bincounts(1)/bintime;
                end
                % alter 'BASELINE' if some of it is an 'offtime'- only
                % include 'ontime' for baseline measurement
                if MDCELLS(n).onTime(1)>MYbaseline_start 
                    if MDCELLS(n).offTime(1)>MYbaseline_end
                        MYbaseline_start=MDCELLS(n).onTime(1);
                    elseif MDCELLS(n).offTime(2)>MYbaseline_end
                        MYbaseline_start=MDCELLS(n).onTime(2);
                    elseif MDCELLS(n).offTime(3)>MYbaseline_end
                        MYbaseline_start=MDCELLS(n).onTime(3);
                  
                    end
                
                elseif MDCELLS(n).offTime(1)<MYbaseline_end
                    if MDCELLS(n).offTime(2)>MYbaseline_end
                        MYbaseline_start=MDCELLS(n).onTime(2);
                    elseif MDCELLS(n).offTime(3)>MYbaseline_end
                        MYbaseline_start=MDCELLS(n).onTime(3);
             
                    end
                end
                baseline_spikes=ontime_spikes(ontime_spikes>=MYbaseline_start);
                baseline_spikes=baseline_spikes(baseline_spikes<=MYbaseline_end);
                baseline_FR=length(baseline_spikes)/(MYbaseline_end-MYbaseline_start);
                if (MYbaseline_end-MYbaseline_start)<(3600*5) % make sure baseline time is long enough to mean something!
                    disp (n)
                    disp ('This cell baseline too short!')
                    remove_at_end = [remove_at_end, RSU_number];
                    
                    continue
                end
                % MAKE SURE no off time during MD period where searching for decline, because that would be misleading   
                %{
                MDperiod=[(2.5-(7-numDays))*3600*24, (5-(7-numDays))*3600*24];
                a=find(sorted_time>MDperiod(1), 1);
                b=find(sorted_time>MDperiod(2), 1);
                c=isnan(sorted_time(a:b-2));
                if length(c(c==1))>1
                    disp ('OffTime during MD time period! Skip this cell:')
                    disp (n)
                    continue
                end                    
                %}
            else % if normal MD cells with normal timeline- count full normal baseline
                baseline_FR=sum(bincounts(4:5))/(bintime*2); % 'baseline FR' is averaged FR over 36 hours right before MD
                
                if MDCELLS(n).onTime(1)>36*3600
                    baseline_sp = ontime_spikes(ontime_spikes>=MDCELLS(n).onTime(1));
                    baseline_sp = baseline_sp(baseline_sp<=60*3600);
                    if ((60*3600)-MDCELLS(n).onTime(1))<(3600*5) % make sure baseline time is long enough to mean something!
                    disp (n)
                    disp ('This cell baseline too short!')
                    remove_at_end = [remove_at_end, RSU_number];

                    continue
                    else
                        baseline_FR = length(baseline_sp)/((60*3600)-MDCELLS(n).onTime(1));
                    end
                end
                
            end
               

            % Normalize all FRs to baseline 

            edges2=experiment_time_start:bintime2:experiment_time_end; % for finding normalized FR
            edges3=experiment_time_start:bintime3:experiment_time_end; % for sliding window to find decline bins
            [bincounts2, binedges2]=histcounts(ontime_spikes, edges2);
            [bincounts3, binedges3]=histcounts(ontime_spikes, edges3);

            FR_normalized=(bincounts2/bintime2)/baseline_FR;
            animal_structure(RSU_number).baseline_FR = baseline_FR;
            for ff=1:length(MDCELLS(n).onTime)
                if MDCELLS(n).offTime(ff)<experiment_time_end
                    if length(MDCELLS(n).onTime)>ff
                        FR_normalized(round(MDCELLS(n).offTime(ff)/bintime2):round(MDCELLS(n).onTime(ff+1)/bintime2))=NaN;
                    else
                        FR_normalized(round(MDCELLS(n).offTime(ff)/bintime2):end)=NaN;
                    end
                end
            end
            FR_normalized(FR_normalized>12)=NaN; % exclude normalized FR above 12
            % normalized FR in smaller bins to find slope using regression
            % of sliding window of these FR's 
            FR_normalized_slidingwindow=(bincounts3/bintime3)/baseline_FR;
 
            FR_normalized_slidingwindow(FR_normalized_slidingwindow>12)=NaN; % filter out too high normalized FR (error)

            animal_structure(RSU_number).normFR=[xaxisdays(1:length(FR_normalized)); FR_normalized]; % now can plot this cell's normalized FR just with row 1 against row 2 of this field
            
            % find start index of MD, depending on experiment
            if myinput==0 
                if MDCELLS(n).dayStart==0
                    MD_start_bintime2=round((3600*24*2.65)/bintime2); %start analyzing MD FR's compared to baseline FR; in bin number
                    MD_start_bintime3=round((3600*24*2.65)/bintime3); 
                elseif MDCELLS(n).dayStart==1
                    MD_start_bintime2=round((3600*24*1.65)/bintime2);
                    MD_start_bintime3=round((3600*24*1.65)/bintime3); 
                elseif MDCELLS(n).dayStart==2
                    MD_start_bintime2=round((3600*24*.65)/bintime2);
                    MD_start_bintime3=round((3600*24*.65)/bintime3); 
                end
            else
                MD_start_bintime2=round((3600*24*2.65)/bintime2); 
                MD_start_bintime3=round((3600*24*2.73)/bintime3); %start analyzing MD FR's compared to baseline FR; in bin number
            end
            MD_end_bintime3=MD_start_bintime3+((3600*24*2)/bintime3); % the MD time that's actually being analyzed (with the sliding windows)
            
            slidingwindow_MD=FR_normalized_slidingwindow(MD_start_bintime3:MD_end_bintime3); %normalized FR for sliding window only during MD period
            
            %%% pick out outlier points (2 or less) that are more than 3 normalized
            %%% Hz away from all other points- just for slope finding
            %%% purposes!!   
            sorted_FRs=sort(slidingwindow_MD);
            high=find (slidingwindow_MD>=sorted_FRs(end-2)+3);
            slidingwindow_MD(slidingwindow_MD>=sorted_FRs(end-2)+3)=NaN;
            if ~isempty (high)
                high=high.*(bintime3)./(3600*24);
                disp (n)
                disp ('Has freak high FR: taken out for slope calc')
                disp (high)
            end
            
            xaxis_slope=experiment_time_start:bintime3:experiment_time_end;
            xaxis_slope_days=xaxis_slope./(3600*24);
     
            %find mean reg slope over whole experiment
            full_regression_slopes=NaN(1, length(xaxis_slope));
            
            % loop through sliding window FR's: 
            % SET SLIDING WINDOW SIZE:
            % First number in this loop;  slide in intervals of
            % bintime3
            for i=8:length(xaxis_slope)-1 
                
                x=xaxis_slope(i-7:i);
                y=FR_normalized_slidingwindow(i-7:i);
                reg_coeff=polyfit(x,y,1);
                reg_slope=reg_coeff(1,1);
                full_regression_slopes(i-7)=reg_slope;
            end
            
            animal_structure(RSU_number).FRslope=[xaxis_slope_days; full_regression_slopes];
            
            
            % x axis if plotting regression slopes
            xaxis_MD_slope=xaxis_slope(MD_start_bintime3: MD_end_bintime3);
            xaxis_MD_slope_days=xaxis_MD_slope./(3600*24);
            
            % same thing as above, but just for MD period
            % set sliding window period!
            regression_slopes=NaN(1, length(slidingwindow_MD));
            Rsquared_cell=[];
            for i=8:length(slidingwindow_MD) 
      
                x=xaxis_MD_slope(i-7:i);
                y=slidingwindow_MD(i-7:i);
                                
                [reg_coeff, S]=polyfit(x,y,1);
                mdl=fitlm(x, y);
                if MDCELLS(n).deprived==1
                    Rsquared_cell(end+1)=mdl.Rsquared.Adjusted;
                    mdls_deprived(end+1)=mdl.Rsquared.Adjusted;
                else
                    Rsquared_cell(end+1)=mdl.Rsquared.Adjusted;
                    mdls_control(end+1)=mdl.Rsquared.Adjusted;
                end
                % Only use the regression slope if Rsquared is over 0;
                % otherwise assume slope of the FRs is actually not
                % changing in any real direction
                if mdl.Rsquared.Adjusted>0
                    reg_slope=reg_coeff(1,1);
                else
                    reg_slope=0;
                end
                regression_slopes(i-7)=reg_slope;
      
            end
           
            % mean/sd for regression slopes during MD period only
            MD_slope_average=nanmean(regression_slopes);
            MD_slope_sd=nanstd(regression_slopes);
            
            % THRESHOLD of slopes to find 'decline' in FR
            threshold= MD_slope_average-MD_slope_sd;  
            pos_threshold= MD_slope_average+MD_slope_sd;
            

            % Find slopes below this threshold: cell specific
            %slopes_below_threshold=regression_slopes(regression_slopes<threshold); % for picking out especially NEGATIVE FR slopes
            slopes_above_threshold=regression_slopes(regression_slopes>pos_threshold); % for picking out especially POSITIVE FR slopes
            
            % Find slopes below low threshold but also NOT within one hour
            % of slopes above high threshold
            low_slope_indices=find (regression_slopes<threshold);
            delete= [];
            for i=1:length(low_slope_indices)
                in=low_slope_indices(i);
                if in==4
                    range=regression_slopes(in-3:in-1);
                elseif in==3
                    range=regression_slopes(in-2:in-1);
                elseif in==2
                    range=regression_slopes(in);
                elseif in==1
                    continue
                else
                    range=regression_slopes(in-4:in-1);
                end
                if sum(ismember(range, slopes_above_threshold))>=1
                    delete(end+1)=i;
         
                end
            end
            low_slope_indices(delete)=[];
            slopes_below_threshold=regression_slopes(low_slope_indices);
            
            
            animal_structure(RSU_number).slopeThreshold=threshold;
            animal_structure(RSU_number).posThreshold=pos_threshold;
 
            % KERNEL DENSITY ESTIMATE
            % finds where in the given distribution of slopes (MD period only) you are most
            % likely to find slopes of a given creteria (below set
            % threshold from above)
            xi=1:length(regression_slopes); % full range possible
            f=low_slope_indices;
            % in case no slopes below threshold
            if length(slopes_below_threshold)<1
                disp(n)
                disp('no slopes below threshold')
                disp('deprived status:')
                disp(MDCELLS(n).deprived)
                f=1;
            end
            
            bandwidth=length(xi)*.05;
            ksd=ksdensity(f, xi, 'function', 'pdf', 'bandwidth', bandwidth);  
            kde_threshold=.7; % threshold probability for first chunk of significant slopes below threshold
            %kde=ksd./max(ksd); % normalize kde distribution
            
            
            % KDE NUMBER 2: POSITIVE KERNEL DENSITY ESTIMATE
            % probability distribution of slopes ABOVE pos threshold
            xi2=1:length(regression_slopes); % full range possible
            f2=find(regression_slopes>pos_threshold); % slopes that fit criteria
            
            % in case no slopes above threshold
            if length(slopes_above_threshold)<1
                disp(n)
                disp('no slopes above threshold')
                disp('deprived status:')
                disp(MDCELLS(n).deprived)
                f=1;
            end
            
            bandwidth2=length(xi2)*.05;
            ksd2=ksdensity(f2, xi2, 'function', 'pdf', 'bandwidth', bandwidth2);  
            kde_threshold2=.5; % threshold probability for significant slopes above threshold
            %kde2=ksd2./max(ksd2); % normalize kde distribution
            
            ksd3 = ksd-ksd2;
            kde = ksd3./max(ksd3);
            
            % FIND DECLINE WINDOW! BASED ON KDE'S
      
            % using the biggest peak of kde distribution
            
            [pks, loc]=findpeaks(kde);
            %mypks=findpeaks(kde);
            %loc=mypks.loc;
            pks=kde(loc);
            max_peak_loc=loc(find(pks==max(pks)));
            decline_window_start_temp=kde(1:max_peak_loc);
            decline_window_start=find(decline_window_start_temp<=kde_threshold,1,'last');
            %decline_window_start=find(kde>=kde_threshold, 1); % first 15 min of decline just found
            %check to make sure they all have a start bin-they do!
            if decline_window_start>0
                x=1;
            else
                decline_window_start=1;
                disp(n)
                disp('No Low FR Slopes During Drop Window:  Control Cell?')
            end
   
            % find last bin of decline just found
            decline_window_end_temp=find(kde(decline_window_start:end)<kde_threshold,2);
            if length(decline_window_end_temp)<2
                decline_window_end_temp=length(kde(decline_window_start:end));
            else
                decline_window_end_temp=decline_window_end_temp(2);
            end
            if decline_window_end_temp>0 %check to see if they all have an end bin
                decline_window_end=decline_window_start+decline_window_end_temp; % add back bins up to the start of the decline, so bins are numbered from start of expt
            else
                decline_window_end=kde(end);
                disp(n)
                disp('Drop ends at the very end of Drop time period: Late or long drop?')
            end
        
            
            decline_window_start=decline_window_start+MD_start_bintime3; %add baseline time so that bins are referenced from start of the experiment instead of start of MD
            decline_window_end=decline_window_end+MD_start_bintime3;
            
            start_FR = mean(FR_normalized_slidingwindow(decline_window_start-48:decline_window_start));
            end_FR = mean(FR_normalized_slidingwindow(decline_window_end:decline_window_end+4));
            
            % check and exclude cells that don't drop significantly below
            % baseline
            if end_FR>(.8*baseline_FR) && MDCELLS(n).deprived == 1
               remove_at_end = [remove_at_end, RSU_number];
               no_drop = [no_drop, RSU_number];
               no_drop_ns = [no_drop_ns, n];

            end
            
            diff_FR = abs(end_FR-start_FR);
            new_start_FR = .2*diff_FR;
            % UPDATE DECLINE WINDOW START
            add_2_start = find(FR_normalized_slidingwindow(decline_window_start:decline_window_end)<(start_FR-new_start_FR),1);
            if ~isempty(add_2_start)
                decline_window_start = decline_window_start+add_2_start;
            end
  
            decline_bins=[decline_window_start decline_window_end]; % end decline bin index numbers; for whole experiment binned in 15 min bins
            
            min_FR_MD=find(slidingwindow_MD==min(slidingwindow_MD),1); % index of minimum FR in MD window
            animal_structure(RSU_number).kde=[xaxis_MD_slope; kde]; % To graph the kernel density estimate probability, since only calculating it over 3.5 days of MD to find decline bins
            animal_structure(RSU_number).decline_bins=decline_bins; % 15 min bin index numbers!: Found by the first chunk of negative FR slope over the .7 kde probability threshold
            
            % plot and store details of each cell in animal_structure
          
            if MDCELLS(n).deprived==1  %if cell from deprived hemisphere
                r=r+1;
                animal_structure(RSU_number).MD=1;

                % for plotting mean and sd of all cells!
                if numDays==6
                    mm=(3600*24)/(30*60);
                elseif numDays==5
                    mm=(3600*24*2)/(30*60);
                else
                    mm=1;
                end
                all_FR_normalized_RSU(m, mm:length(FR_normalized)+mm-1)= FR_normalized;
                RSU_count=RSU_count+1;
                
                    myfigname = ['Deprived' num2str(n)];
                    figure ('Name', myfigname)
                    disp(n)
                    disp('is figure number')
                    disp(r)
                    subplot (2,1,1)
                    yyaxis left
                    %line([min(xaxisdays), max(xaxisdays)], [1,1], 'Color', 'k')
                    plot(xaxisdays(1:length(FR_normalized)), FR_normalized.*baseline_FR, 'k')
                    hold on 
                    plot ([(decline_window_start*(15*60))/(3600*24), (decline_window_start*(15*60))/(3600*24)], [0,2.5], 'r', 'LineWidth', 3)
                    hold on
                    plot ([(decline_window_end*(15*60))/(3600*24), (decline_window_end*(15*60))/(3600*24)], [0,2.5], 'r', 'LineWidth', 3)
                    hold on

                    for v=1: length(MDCELLS(n).offTime)
                        plot ([MDCELLS(n).offTime(v)/(3600*24), MDCELLS(n).offTime(v)/(3600*24)], [0,3], 'b')
                    end
                    for w=1: length(MDCELLS(n).onTime)
                        plot ([MDCELLS(n).onTime(w)/(3600*24), MDCELLS(n).onTime(w)/(3600*24)], [0,3], 'y')
                    end
                        
                    xlabel('Time(days)')
                    
                    ylabel('Normalized FR')
                    
                    title(MDCELLS(n).animal)
                    set (gca, 'fontsize', 15, 'box', 'off')
                    set (gca, 'xlim', [0,8])

                    
                    subplot(2,1,2)
                    yyaxis left
                    %line([min(xaxis_slope_days), max(xaxis_slope_days)], [0,0], 'Color',  'k', 'LineWidth', .3)
                    hold on
                    plot(xaxis_slope_days, full_regression_slopes, 'b')
                    hold on
                    plot(xaxis_MD_slope_days(low_slope_indices), slopes_below_threshold, 'xr')
                    hold on
                    plot (xaxis_MD_slope_days(regression_slopes>pos_threshold), slopes_above_threshold, 'xc')
                    hold on
                    plot ([(decline_window_start*(15*60))/(3600*24), (decline_window_start*(15*60))/(3600*24)], [0,2.5], 'r', 'LineWidth', 3)
                    hold on
                    plot ([(decline_window_end*(15*60))/(3600*24), (decline_window_end*(15*60))/(3600*24)], [0,2.5], 'r', 'LineWidth', 3)
                    hold on
                    xlabel('Time(days)')
                    ylabel('Normalized FR Slope')
                    set(gca, 'ylim', [1.5*(min(full_regression_slopes)), 1.5*(max(full_regression_slopes))])
                    hold on
                    set (gca, 'fontsize', 15, 'box', 'off')
                    hold on
                    yyaxis right
                    hold on
                    line([min(xaxis_slope_days), max(xaxis_slope_days)], [.7, .7],'Color', [.1,.1,.1], 'LineStyle', '--', 'LineWidth', .3)
                    plot(xaxis_MD_slope_days, kde, 'r')
                    hold on
                    %plot (xaxis_MD_slope_days, kde2, 'c')

                    ylabel('Normalized Probability')
                    set (gca, 'ylim', [0,1], 'box', 'off')
                    hold off
                    set (gca, 'fontsize', 15)
                    set (gca, 'xlim', [0,8])
                    
            else  %cell from non-deprived hemisphere
                r=r+1;
                animal_structure(RSU_number).MD=0;
                
                if numDays==6
                    mm=(3600*24)/(30*60);
                elseif numDays==5
                    mm=(3600*24*2)/(30*60);
                else
                    mm=1;
                end
                
                FR_normalized_RSU_control(m, mm:length(FR_normalized)+mm-1)= FR_normalized;
                FR_normalized_RSU_control(m, 1:mm)=1;
                RSU_control_count=RSU_control_count+1;
               
               %{
                    figure ('Name', 'non-deprived')
                    subplot (2,1,1)
                    plot(xaxisdays(1:length(FR_normalized)), FR_normalized)
                    hold on
                    for v=1: length(MDCELLS(n).offTime)
                        plot ([MDCELLS(n).offTime(v)/(3600*24), MDCELLS(n).offTime(v)/(3600*24)], [0,3], 'b')
                    end
                    for w=1: length(MDCELLS(n).onTime)
                        plot ([MDCELLS(n).onTime(w)/(3600*24), MDCELLS(n).onTime(w)/(3600*24)], [0,3], 'y')
                    end
                    xlabel('time(days)')
                    ylabel('Normalized FR')
                    
                    subplot(2,1,2)
                    scatter(xaxis_slope, full_regression_slopes)
                    xlabel('time(days)')
                    ylabel('Normalized FR Slope')
                  %}

            end

    end
    
   
end

animal_structure_nodrop = animal_structure;
yes_drop = [1:length(animal_structure)];
yes_drop = yes_drop(~ismember(yes_drop, no_drop));
animal_structure_nodrop(yes_drop) = [];

animal_structure(remove_at_end) = [];
%set(gca,'yscale','log');
disp('number of on cells')
disp(ON_cells)

FR_means_RSU=nanmean(all_FR_normalized_RSU);
FR_means_RSU=FR_means_RSU(~isnan(FR_means_RSU));
FR_means_RSU_control=nanmean(FR_normalized_RSU_control);
FR_means_RSU_control=FR_means_RSU_control(~isnan(FR_means_RSU_control));

FR_SEMs_RSU=(nanstd(all_FR_normalized_RSU))./(sqrt(RSU_count)-1);
FR_SEMs_RSU=FR_SEMs_RSU(~isnan(FR_SEMs_RSU));
FR_SEMs_RSU_control= (nanstd(FR_normalized_RSU_control))./(sqrt(RSU_control_count));
FR_SEMs_RSU_control=FR_SEMs_RSU_control(~isnan(FR_SEMs_RSU_control));

xaxis = 0:bintime2:9.5*24*3600;
   

xaxisdays=xaxis/(3600*24); %now xaxis in days not seconds

disp('RSU_count')
disp(RSU_count)
disp('RSU_control_count')
disp(RSU_control_count)

totalcells=RSU_count+RSU_control_count;

disp ('total cells plotted:' )
disp(totalcells)


figure ()

RSU_control=plot(xaxisdays(1:length(FR_means_RSU_control)), FR_means_RSU_control, 'k', 'DisplayName', 'control');
shadedErrorBar(xaxisdays(1:length(FR_means_RSU_control)), FR_means_RSU_control, FR_SEMs_RSU_control, 'k');
hold on
plot(xaxisdays(1:length(FR_means_RSU)), FR_means_RSU, 'm', 'DisplayName', 'deprived');
shadedErrorBar(xaxisdays(1:length(FR_means_RSU)), FR_means_RSU, FR_SEMs_RSU, 'm');
%line([0,8], [1,1], 'Color', 'k')

xlabel('Time (days)')
ylabel('Normalized Firing Rate')
title('RSU''s: deprived vs. control')
set (gca, 'fontsize', 15)
set (gca, 'xlim', [0, 8], 'box', 'off')
%%

%{
% For visualizing/setting cutoff for RSU vs FS

all_neg_pos_time=[MDCELLS.neg_pos_time];
all_tailSlope=[MDCELLS.tailSlope];



%now pick out only neg_pos_time and tailSlope from active cells (determined
%above in array ON_cells_index)

%ON_cells_index=ON_cells_index(~isnan(ON_cells_index));
q3_cells_index = find([MDCELLS.quality]<3);
ON_cells_index = q3_cells_index;
ON_neg_pos_time=all_neg_pos_time(ON_cells_index);
ON_tailSlope=all_tailSlope(ON_cells_index);


X=[ON_neg_pos_time', ON_tailSlope'];
[idx, C]=kmeans(X,2);
indices_1=find(idx==1);
indices_2=find(idx==2);
indices_3=find(idx==3);

neg_pos_time_1=ON_neg_pos_time(indices_1);
tailSlope_1=ON_tailSlope(indices_1);
neg_pos_time_2=ON_neg_pos_time(indices_2);
tailSlope_2=ON_tailSlope(indices_2);
%neg_pos_time_3=ON_neg_pos_time(indices_3);
%tailSlope_3=ON_tailSlope(indices_3);



if min(neg_pos_time_1)<min(neg_pos_time_2) %then indices_1 are FS and indices_2 are RSU
    start_neg_pos_time_FS=min(neg_pos_time_1);
    end_neg_pos_time_FS=max(neg_pos_time_1);
    start_neg_pos_time_RSU=min(neg_pos_time_2);
    end_neg_pos_time_RSU=max(neg_pos_time_2);
else
    start_neg_pos_time_FS=min(neg_pos_time_2);
    end_neg_pos_time_FS=max(neg_pos_time_2);
    start_neg_pos_time_RSU=min(neg_pos_time_1);
    end_neg_pos_time_RSU=max(neg_pos_time_1);
end


figure (500)
plot (neg_pos_time_1, tailSlope_1, 'xk', 'LineWidth', 3)
hold on
%plot(neg_pos_time_3, tailSlope_3, 'xk')
plot (neg_pos_time_2, tailSlope_2, 'xm', 'LineWidth', 3)
xlabel('Neg Pos Time')
ylabel('Tail Slope')
disp(length(all_neg_pos_time))
disp(length(neg_pos_time_1)+length(neg_pos_time_2))
set (gca, 'fontsize', 20, 'box', 'off')
set(gca,'ylim',[-0.2 0.2],'xlim',[0 25]);
title ('K-Means Clustering')

%}
