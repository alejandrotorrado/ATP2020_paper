function [datout,datraw,feat_out] = score_10sec_noMLS (data, tx, Fs)
% Return a binary assessment of waking or not waking for each sample point
% of a set of polysomnographic data including the integrated delta and
% theta power of the LFP, as well as the integrated EMG signal recorded
% simultaneously. DATA is a matrix. Each column is a variable of interest.
% Each row is a timepoint. TX is the 1 x nsample points vector containing
% the time indices of each data point in DATA.
% Fs is the sampling rate of the data.
%
%       NOTE: Organization of DATA = [deltapower, thetapower, deltatheta_difference, rs_mvt_z, emg_zscore];
%       NOTE: REM = 1, NREM = 2, ACTIVE = 4, QUIET = 5
% keyboard;


%%
nFeatures = 11; % hard-coded for now
bin_size = 10; % seconds
edges       = 0: bin_size*Fs :size(data,1);
AWTHRESH    = 3.0; % SET THIS HIGH?
QWTHRESH    = 0.7;
EMGTHRESH   = 0.9;

if edges(end)~=(size(data,1))
    edges(end + 1) = size(data,1);
end

verbose = 0;

%{
if verbose
    delta   = data(:,1);
    theta   = data(:,2);
    dtdf   = data(:,3);
    mvmt    = data(:,4);
    emg     = data(:,5);
    
    figure (3); hold on;
    yyaxis left
    area(tx,delta>.5); % binary decision of NREM sleep?
    plot(tx,delta),hold on % delta power
    plot(tx,theta,'b','linewidth',2); % theta
    
    yyaxis right
    plot(tx,mvmt,'m')
    plot(tx,emg,'--k'); % EMG trace
    
    p = line([tx(1) tx(end)],[AWTHRESH AWTHRESH]);
    set(p,'color','r','linestyle','--');
    xlim([tx(1) tx(end)]);
    legend('Delta','Binary NREM','Movement','EMG','Theta',['AW thresh: ' num2str(AWTHRESH)])
    title('All traces during search for waking');
    ylims = get(gca,'Ylim');
    ymax = ylims(2);
end
%}
%%
autocat = zeros(size(edges,2)-1,3);
features_out = zeros(size(edges,2)-1,nFeatures);


for ee = 1:size(edges,2)-1
    
    % data is currently 2 samples per second. break into 10 second blocks
    t0      = edges(ee)+1;
    t1      = edges(ee+1);
    temp    = data(t0:t1,:);
    
    % First slice input array and retrieve each column
    delta_data      = temp(:,1);
    theta_data      = temp(:,2);
    DTdiff_data     = temp(:,3);
    mvt_data        = temp(:,4);
    emg_data        = temp(:,5);
    
    % Features 1 - 5 are the average of each column of the input matrix
    delta      = nanmean(delta_data);
    theta      = nanmean(theta_data);
    dtdiff         = nanmean(DTdiff_data);
    emg           = nanmean(emg_data);
    if isempty(mvt_data)
        mvmean = NaN;
    else
        mvmean = nanmean(mvt_data);
    end
    
    % Features 6 - 11 are computed here
    emgvar          = var(emg_data);
    emg_ampl        = max(abs(diff(emg_data)));
    DT_slope        = nanmean(diff(DTdiff_data));
    delta_slope     = nanmean(diff(delta_data));
    theta_slope     = nanmean(diff(theta_data));
    mvt_slope       = nanmean(diff(mvt_data));
    
    % feature 11
    if ee>1
        prevscore = autocat(ee-1,1);
    else
        prevscore = 0;
    end
    
    % collect features in array
    features_out(ee,1)     = delta;
    features_out(ee,2)     = theta;
    features_out(ee,3)     = dtdiff;
    features_out(ee,4)     = mvmean;
    features_out(ee,5)     = emg;
    features_out(ee,6)     = emgvar;
    features_out(ee,7)     = emg_ampl;
    features_out(ee,8)     = DT_slope;
    features_out(ee,9)     = delta_slope;
    features_out(ee,10)     = theta_slope;
    features_out(ee,11)    = mvt_slope;
    features_out(ee,12)    = prevscore;
    
    % timestamp of 10-sec bin
    tlims   = [tx(t0) tx(t1)];
    autocat(ee,2) = tlims(1);
    autocat(ee,3) = tlims(2);
    
    % ASSUME MOVEMENT DATA:
    if ~isnan(mvmean)
        % define quadrant of active waking:
        if mvmean > QWTHRESH
            if mvmean >= AWTHRESH || emg > EMGTHRESH || mvmean>(QWTHRESH+(AWTHRESH-QWTHRESH)/2) && emgvar > 0.5
                autocat(ee,1) = 4;
            elseif mvmean < AWTHRESH
                autocat(ee,1) = 5;
            elseif mvmean > QWTHRESH && emg > EMGTHRESH
                autocat(ee,1) = 4;
            end

            % define quiet waking
        elseif delta<.4 && mvmean<AWTHRESH && dtdiff > .1 ||...
                theta<.2 && delta < .4 && emg>0.25 && abs(mvmean)<0.15 && mvmean<AWTHRESH ||...
                theta>.2 && delta>.3 && delta<.4 && mvmean>QWTHRESH/2 && mvmean<AWTHRESH || ...
                all([theta delta]<.3) && delta>theta && mvmean>QWTHRESH/2 && mvmean<AWTHRESH ||...
                all([delta theta]<.3) && abs(dtdiff)<0.5 && mvmean == 0 && emg > 0.75 ||...
                emg>=1.25 && delta>.3 && delta<.5 && theta < 0.3 ||...
                mvmean<0.5 && theta <.2 && delta < .2 && emg <1
            autocat(ee,1) = 5;
            
            %define NREM
        elseif delta>.5 && dtdiff > 0.3 || delta>.5 && theta>.1 && mvmean>0.05 && mvmean<QWTHRESH ||...
                delta>.5 && dtdiff > .3 && mvmean<QWTHRESH
            autocat(ee,1) = 2;
            
            % define REM
        elseif theta > .3 && dtdiff < 0 ||...
                theta>.3 && dtdiff>0 && mvmean<0 && dtdiff > .1 ||...
                theta> .3 && theta>delta && mvmean<0.15 ||...
                theta>.25 && dtdiff > .1 && mvmean<0 || ... 
                theta > .25 && delta < theta && mvmean < 0 && emg < .4 ||...
                delta < .2 && (theta/delta) > 2 && mvmean < 0 || ...
                mvmean < 0 && dtdiff > 0 && emg < 0.5
            autocat(ee,1) = 1;
            
            % weird inbetween REM and NREM states:
        elseif delta> .3 && theta>.3 && mvmean<0.05 ||...
                delta>.3 && dtdiff > 0 && dtdiff < 0.2  && mvmean<0.05 || ...
                all([delta,theta]<.3) && delta>theta && emg < 0.1 && mvmean < 0.05
            autocat(ee,1) = 2;
            
            
            
        elseif ~isnan(mvmean)
            % DEFAULT TO QUIET WAKING?
            %keyboard
            autocat(ee,1) = 5;
            disp('Uncategorized');
        end
    end
    
         % - - - - - - - NO MOVEMENT FROM VIDEO -  - - - - - - - - - - - - -
    if isnan(mvmean)% MVMEAN IS NAN:
    
        % quiet wake
        if dtdiff<.1 && emg<.4 && all([delta theta]<.3) || ...
                emg>0.3 && emg<=.7 && all([delta theta]<.35) && dtdiff < .1 || ...
                emg>=.8 && delta>.25 && delta<.4 && dtdiff< .3 || ...
                dtdiff < .3 && emg < .75 && all([delta theta]<.3) ||...
                all([theta delta]<.15) && emg < .7
            autocat(ee,1) = 5 ;

        % REM - low EMG
        elseif emg<.7 && dtdiff < 0 && theta>.25 || ...
                emg<.5 && dtdiff < -.1 && theta > .2 ||...
                all([delta theta]>.35) && dtdiff<0 || ...
                emgvar < 0.1 && all([theta delta]>.3) && theta > .15 && dtdiff<0||...
                emgvar < 0.1 && dtdiff< -0.1 && theta > .35
            autocat(ee,1) = 1 ;
           
            % NREM
        elseif emg<.4 && dtdiff > 0 && delta>.35 || ...
                delta>.5 && theta<.2 ||...
                delta>.5 && dtdiff>=0.3 && emg<.8 || ...
                delta>.5 && theta>.4 && emg <.4 ||...
                delta>.5 && dtdiff>0 && emg>0.3 ||...
                delta > .35 && dtdiff>.1 ||...
                delta>.2 && theta < .1 && emg>=0.4 && emg < .7
            autocat(ee,1) = 2;
            
            % active wake
        elseif emg>=.7 && delta<.4 && isnan(mvmean) ||...
                emg>=.65 && all([delta theta]<.12)
            autocat(ee,1) = 4;

        else
            disp('Unclear')
            autocat(ee,1) = 5;
            
        end
    end

end

autocat_out = autocat;
ticker = 1;
while ticker == 1
    % eliminate states that only last for one 10-sec bin
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % copy the state assignment
    temp = autocat_out(:,1);
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
    if isempty(temp5)
        ticker = 0;
    end
    % get the indices of the bins that have a 1 bin state
    kills = temp3(temp5);
    % overwrite the 1 bin state with the preceding state
    temp(kills+1) = temp(kills);
    % reassign temp into the autocat variable
    autocat_out(:,1) = temp;
end


if verbose
    clr = [0.5 0.5 0.5; 0.2 0.2 0.3; NaN NaN NaN; 0.4 1.0 0.2; 1 1 0];
    disp('verbose');
    figure (3);
    for ii = 1:size(autocat_out,1);
        t0      = edges(ii)+1;
        t1      = edges(ii+1);
        tlims   = [tx(t0) tx(t1)];
        try
        rectangle( 'Position', [tlims(1), ymax-1, tlims(2) - tlims(1) , 1],...
            'facecolor',clr(autocat_out(ii),:),'linestyle','none');
        drawnow
        catch
            disp('Caught in plotting auto output');
            keyboard
        end
    end
end

if any(unique(autocat_out(:,1)) == 0);
    disp('Algorithm failed to ID and catch a time bin!');
    keyboard
    % find the zero entry in the first column of autocat, set ee to the row
    % number and then go through the above steps to figure out where it
    % should be ID'ing a state and why it's missing this entry. It might
    % help to turn on "verbose" and see the data to more easily understand
    % what it "should" be identifying.
end


% Set up the output variable:
datout = autocat_out;
datraw = autocat;
feat_out = features_out;






