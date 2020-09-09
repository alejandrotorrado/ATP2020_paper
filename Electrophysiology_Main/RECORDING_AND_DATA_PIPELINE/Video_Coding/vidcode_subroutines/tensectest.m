function [datout] = tensectest (data, tx, mvmt)
% Return a binary assessment of waking or not waking for each sample point
% of a set of polysomnographic data including the integrated delta and
% theta power of the LFP, as well as the integrated EMG signal recorded
% simultaneously. DATA should be organized as n sample points x 3 columns.
% The first column is the integral of the delta band power, the second
% column is the integral of the theta band power, and the third column is
% the EMG integral. TX is the 1 x nsample points vector containing the time
% indices of each data points in DATA. TMIN is the minimum time that an
% epoch will be considered valid. Periods of wake or not wake that have a
% duration of less than TMIN will be considered for elimination.
%
%       NOTE: Organization of DATA = [deltapower thetapower rs_emg'];
%       NOTE: REM = 1, NREM = 2, ACTIVE = 4, QUIET = 5

edges       = 0: 20 :size(data,1);
AWTHRESH    = 3.9; % SET THIS HIGH?
QWTHRESH    = 0.7;

if edges(end)~=(size(data,1)),
    edges(end + 1) = size(data,1);
end

verbose = 0;

if verbose;
    figure (3);
    plot(tx,data(:,1)),hold on % delta power
    area(tx,data(:,1)>1.5); % binary decision of NREM sleep?
    plot(mvmt(:,2),mvmt(:,1),'m')
    plot(tx,data(:,3),'--k'); % EMG trace
    plot(tx,data(:,2),'b','linewidth',2); % theta
    p = line([tx(1) tx(end)],[AWTHRESH AWTHRESH]);
    set(p,'color','r','linestyle','--');
    xlim([tx(1) tx(end)]);
    legend('Delta','Binary NREM','Movement','EMG','Theta',[num2str(AWTHRESH) ' Thresh'])
    title('All traces during search for waking');
    ylims = get(gca,'Ylim');
    ymax = ylims(2);
end

autocat = zeros(size(edges,2)-1,3);

for ee = 1:size(edges,2)-1;
    
    % data is currently 2 samples per second. break into 10 second blocks
    t0      = edges(ee)+1;
    t1      = edges(ee+1);
    temp    = data(t0:t1,:);
    delta   = mean(temp(:,1));
    theta   = mean(temp(:,2));
    emg     = mean(temp(:,3));
    emgstd  = std(temp(:,3));
    tlims   = [tx(t0) tx(t1)];
    tmpmvmt = mvmt(mvmt(:,2)>=tlims(1) & mvmt(:,2)<=tlims(2),1);
    if isempty(tmpmvmt);
        tmpmvmt = NaN;
    end
    mvmean  = mean(tmpmvmt);
    
    autocat(ee,2) = tlims(1);
    autocat(ee,3) = tlims(2);
    
    % ASSUME MOVEMENT DATA:
    if ~isnan(mvmean);
        % define quadrant of active waking:
        if mvmean > QWTHRESH;
            if mvmean >= AWTHRESH || emg > 2 || mvmean>(QWTHRESH+(AWTHRESH-QWTHRESH)/2) && emgstd > 0.5;
                autocat(ee,1) = 4;
            elseif mvmean < AWTHRESH;
                autocat(ee,1) = 5;
            elseif mvmean > QWTHRESH && emg > 1.5;
                autocat(ee,1) = 4;
            end

            % define quiet waking
        elseif delta<2.25 && theta/delta<1.5 && mvmean<AWTHRESH && (delta/theta)<2 ||...
                theta<2 && delta<2.5 && emg>0.25 && mvmean>0.05 && mvmean<AWTHRESH ||...
                theta>2 && delta>2 && delta<2.5 && mvmean>QWTHRESH/2 && mvmean<AWTHRESH || ...
                all([theta delta]<2) && delta>theta && mvmean>QWTHRESH/2 && mvmean<AWTHRESH ||...
                all([delta theta]<2) && abs(theta-delta)<0.5 && mvmean == 0 && emg > 0.75 ||...
                emg>=1.25 && delta>2 && delta<4 || mvmean<0.5 && theta <1.3 && delta < 1.0 &&...
                emg <1;
            autocat(ee,1) = 5;
            
            %define NREM
        elseif delta>2 && (theta/delta)<1.2 || delta>1.9 && theta>2 && mvmean>0.1 && mvmean<QWTHRESH ||...
                delta>2 && abs(theta/delta)<1.4 && mvmean<QWTHRESH;
            autocat(ee,1) = 2;
            
            % define REM
        elseif theta > 1.3 && (theta/delta)>1.5 || theta>2 && (theta/delta)>1.1 && mvmean<0.3 && (theta-delta)>0.2 ||...
                theta>2.5 && theta>delta && mvmean<0.15 || theta>1 && (theta/delta)>1.5 && mvmean<0.1 ... 
                || theta > 1.5 && delta < theta && mvmean < 0.15 && emg < 1.2 || delta <0.7 && (theta/delta) > 2 && mvmean < 0.05 || ...
                mvmean == 0 && (theta-delta)>0.35 && emg < 1.0;
            autocat(ee,1) = 1;
            
            % weird inbetween REM and NREM states:
        elseif delta> 2.5 && theta>2.5 && mvmean<0.1 || delta>1.9 && abs(delta-theta)<0.2 && mvmean<0.4 || ...
            all([delta,theta]<2) && delta>theta && emg < 0.5 && mvmean < 0.4 ;
            autocat(ee,1) = 2;

        
        
        elseif ~isnan(mvmean);
            % DEFAULT TO QUIET WAKING?
            %keyboard
            autocat(ee,1) = 5;
            disp('Uncategorized');
        end
    end
    
         % - - - - - - - NO MOVEMENT FROM VIDEO -  - - - - - - - - - - - - -
    if isnan(mvmean)% MVMEAN IS NAN:
    
        % quiet wake
        if (theta-delta)<1 && emg<1.05 && all([delta theta]<2.6) || ...
                emg>0.75 && emg<=1.3 && all([delta theta]<2.6) && (delta/theta)<1.9 ...
                || emg>=1.25 && delta>2 && delta<3.5 && (delta/theta)<1.5 || ...
                delta/theta < 2.5 && emg < 1.2 && all([delta theta]<2) ||...
                all([theta delta]<1.10) && emg < 1.15;
            autocat(ee,1) = 5 ;

        % REM - low EMG
        elseif emg<1.1 && (theta-delta)>0.2 && theta>1.75 || emg<1.0 && (theta-delta)>1 && theta> 1 ||...
                all([delta theta]>2.3) && (theta/delta)>1.2 || emgstd < 0.1 && all([theta delta]>2) && theta> 1 && theta/delta > 1||...
                emgstd < 0.1 && (theta/delta)>1.3 && theta > 2.5 ;
            autocat(ee,1) = 1 ;
           
            % NREM
        elseif emg<1.2 && delta>theta && delta>2 || delta>2.5 && theta<1.75 ||...
                delta>2.5 && (delta/theta)>=1.2 && emg<1.5 || delta> 3 && theta>2.5 && emg <1.0 ||...
                delta>2.5 && delta>theta && emg>0.75 || delta > 2.2 && (theta/delta) < 1.2 ||...
                delta>1.10 && theta < 1 && emg>=0.9 && emg <= 1.25;
            autocat(ee,1) = 2;
            
            % active wake
        elseif emg>=1.25 && delta<2.5 && isnan(mvmean) ||...
                emg>=1.18 && all([delta theta]<1.05);
            autocat(ee,1) = 4;

        else
            disp('Unclear')
            autocat(ee,1) = 5;
            
        end
    end

end


ticker = 1;
while ticker == 1;
    % eliminate states that only last for one 10-sec bin
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % copy the state assignment
    temp = autocat(:,1);
    %find the differential. zero when a state is constant
    temp2 = diff(temp);
    % make it binary. all changes are 1s and constants are 0s
    temp2(temp2~=0) = 1;
    % index the changes
    temp3 = find(temp2 == 1);
    % find the amount of bins between changes
    temp4 = diff(temp3);
    % look for bins changes that are only seperated by 1
    temp5 = find(temp4==1);
    % if you have staggered 1-bin states flanking one another, this will
    % need to run a couple of times to iteratively clean them out. move on
    % when finished:
    if isempty(temp5);
        ticker = 0;
    end
    % get the indices of the bins that have a 1 bin state
    kills = temp3(temp5);
    % overwrite the 1 bin state with the preceding state
    temp(kills+1) = temp(kills);
    % reassign temp into the autocat variable
    autocat(:,1) = temp;
end


if verbose;
    clr = [0.5 0.5 0.5; 0.2 0.2 0.3; NaN NaN NaN; 0.4 1.0 0.2; 1 1 0];
    disp('verbose');
    figure (3);
    for ii = 1:size(autocat,1);
        t0      = edges(ii)+1;
        t1      = edges(ii+1);
        tlims   = [tx(t0) tx(t1)];
        try
        rectangle( 'Position', [tlims(1), ymax-1, tlims(2) - tlims(1) , 1],...
            'facecolor',clr(autocat(ii),:),'linestyle','none');
        drawnow
        catch
            disp('Caught in plotting auto output');
            keyboard
        end
    end
end

if any(unique(autocat(:,1)) == 0);
    disp('Algorithm failed to ID and catch a time bin!');
    keyboard
    % find the zero entry in the first column of autocat, set ee to the row
    % number and then go through the above steps to figure out where it
    % should be ID'ing a state and why it's missing this entry. It might
    % help to turn on "verbose" and see the data to more easily understand
    % what it "should" be identifying.
end


% Set up the output variable:
datout = autocat;






