%% TRAIN NEW MODEL
clearvars -except Mdl CELL_BACKUP pDir animal CONTCELL* CELL*

ncol = 20;

make_new_model = 0;

if make_new_model
    
    if exist('pDir','var')
        savedpDir = pDir;
    else
        savedpDir = uigetdir(cd,'Point me to the folder with saved p variables.');
    end
    savedPs = dir([savedpDir filesep '*p.mat']);
    saved_p = [];
    for pp = 1:size(savedPs,1)
        pload = load([savedpDir filesep savedPs(pp).name]);
        saved_p = [saved_p; pload.p];
        clear pload
    end
    
    training = saved_p(:,1:ncol-1);
    scoredclass = saved_p(:,ncol);
    rng(1)
    clear Mdl;
    
    % make a logical vector that defines which columns in p are categorical
    catpre = zeros(1,ncol-1);
    catpre(14) = 1;
    catpre = logical(catpre);
    
    tic
    Mdl = TreeBagger(500,training,scoredclass,'OOBPrediction','On',...
        'Method','classification', 'CategoricalPredictors', catpre,...
        'PredictorNames',...
        {'gam1','gam2','cont','h_cont','FR','kFS','kRS','dFS','dRS','xRS',...
        'xFS','rise','decay','crv','amp','negpos','HW','tS','AC'},...
        'OOBPredictorImportance','On'); %'OOBPredictorImportance','On'
    toc
    view(Mdl.Trees{1},'Mode','graph');    
    
    oobfig = figure(); hold on;
    set(oobfig,'units','normalized','position',[.1 .1 .6 .5]);
    oobErrorBaggedEnsemble = oobError(Mdl);
    plot(oobErrorBaggedEnsemble,'color',[0 0 0 .7],'linewidth',2)
    xlabel 'Number of grown trees';
    ylabel 'Out-of-bag classification error';
    set(gca,'ylim',[0 0.4],'fontsize',12,'Xcolor','k','YColor','k');
    [minoob,minoob_idx] = min(oobErrorBaggedEnsemble);
    minline = refline(0,minoob);
    minline.Color = [0.30 .75 .93 0.7];
    minline.LineWidth = 2;
    text(40,0.35,sprintf('Best performance (OOB Error = %.2f%%) with %u trees.',100*minoob,minoob_idx),...
        'fontsize',11);
    text(40,0.31,sprintf('Final performance: OOB Error = %.2f%%',100*oobErrorBaggedEnsemble(end)),...
        'fontsize',11);
    
    
    save_model = 1;
    if save_model
        save([savedpDir filesep 'trained_Mdl.mat'],'Mdl');
    end
    
    fprintf('    \nNew model made and saved! Look in:\n%s\n\n',[savedpDir filesep 'trained_Mdl.mat']);
    fprintf('Press "continue" if you want to keep running this script.\n');
    keyboard;
else
    % load previous Mdl
    if ~exist('Mdl','var')
        % load('/Users/khengen/Google_Drive/Matlab_scripts_11_01_2012/Cell_Quality/RANDOMFOREST_TRAINED/rForest_Mdl_KH72_73_75_SHANK01_AT13_14_trained.mat');
        mdl_dir = '/Volumes/turrigiano-lab/RECORDING_AND_DATA_PIPELINE/Data_Processing/CELL_Creation/RANDOM_FOREST_TRAINED/LATEST';
%         mdl_dir = '/Users/atorrado/Google_Drive/Matlab_scripts_11_01_2012/Cell_Quality/RANDOMFOREST_TRAINED/SAVED_Ps/ATP_newCont_032817';
        mdl_name = 'trained_Mdl_KH73_75_67_45_AT12_19_20.mat';
        fprintf('Loading model: %s\n',mdl_name);
        mdl_load = load([mdl_dir filesep mdl_name]);
        Mdl = mdl_load.Mdl;
    else
        fprintf('Model already loaded.\n');
    end
end



fprintf('Ncol is set to: %u.\nModel trained on data with %u columns.\n',ncol,size(Mdl.X,2));
fprintf('\n *** Note that if number of columns is NOT ncol-1 you may be using the wrong model. ***\n\n');
pause(1.5);

% Hard-coded confidence threhsolds for quality review flagging
qconf_thresh = 0.80;
sbconf_thresh = 0.15;

interpx = 3;

if interpx == 5
    nsamp = 91;
    interp_samples = 161;
    swprng = -30:60;
    min_rng = 20:50;
elseif interpx == 3
    nsamp          = 59;
    interp_samples = 97;
    swprng = -18:40;
    min_rng = 10:30;
end

%% MAKE TEMPLATES
FSwindow =  [0.1438    0.1438    0.1440    0.1441    0.1440    0.1354    0.1125    0.0474...
    -0.0913   -0.3149   -0.6063   -0.8860   -1.0000   -0.8648   -0.5677   -0.2432...
    0.0295    0.2162    0.3171    0.3561    0.3578    0.3383    0.3067    0.2683...
    0.2267    0.1849    0.1447    0.1074    0.0736    0.0437    0.0180   -0.0028...
    -0.0181   -0.0274   -0.0308   -0.0290   -0.0252];

RSUwindow = [ 0.2884    0.2884    0.2887    0.2954    0.3040    0.3032    0.2783    0.1945...
    0.0137   -0.2757   -0.6227   -0.9030   -1.0000   -0.9124   -0.7313   -0.5368...
    -0.3646   -0.2246   -0.1126   -0.0196    0.0593    0.1253    0.1778    0.2170...
    0.2441    0.2608    0.2692    0.2709    0.2673    0.2593    0.2475    0.2320...
    0.2122    0.1873    0.1577    0.1337    0.1314];

% over sample and spline fit the windows established above (taken from the
% CONTCELL dataset at a 32 sample point per waveform rate [24414.14 Hz])
% and upsampled to match the newer approach tested and implemented summer
% 2016.
%THIS IS FOR 5x interpolation
% x0      = 1:size(FSwindow,2);
% x1      = 1:(36/160):size(FSwindow,2);
% FSwin2  = interp1(x0,FSwindow,x1,'spline');
% RSUwin2 = interp1(x0,RSUwindow,x1,'spline');


% this is for 3x interpolation
x0      = 1:size(FSwindow,2);
x1      = 1:(36/(interp_samples-1)):size(FSwindow,2);
FSwin2  = interp1(x0,FSwindow,x1,'spline');
RSUwin2 = interp1(x0,RSUwindow,x1,'spline');


% SEE NOTE BELOW to adjust the centering component for other sampling
% rates. It's commented out in the main loop, but I've hard coded it here
% for the sake of speed on this first pass.


% figure out which channels have data on them
% chan = [];
% for ee = 1:size(channel,2);
%
%     if ~isempty(channel(ee).block);
%         chan(end+1) = ee;
%     end
%
% end

% find the sample at which the windows are at minimum:
[~,rowmin0] = min([RSUwin2;FSwin2],[],2);
% check to make sure that the minimums are in the same sample:
if numel(unique(rowmin0)) == 1;
    rowmin0 = unique(rowmin0);
    %     rowmin1 = 31; % SEE NOTE ABOVE AND COMMENTED CODE BELOW IF THIS NEEDS
    %     TO BE UPDATED. - 9/19/16 updated rowmin1 calculation - it's now in
    %     the loop
else
    error('The interpolated waveform windows do no share the same minimum sample point! Stopped in beta_cluster_evaluation.m');
end

% if waveforms have been trimmed for PCA, trim down the templates to match
if nsamp ~= interp_samples
    trim_start = rowmin0 + swprng(1);
    trim_end = rowmin0 + swprng(end);
    FSwin2 = FSwin2(trim_start:trim_end);
    RSUwin2 = RSUwin2(trim_start:trim_end);
end

% then recalculate rowmin0
[~,rowmin0] = min([RSUwin2;FSwin2],[],2);
% check to make sure that the minimums are in the same sample:
if numel(unique(rowmin0)) == 1
    rowmin0 = unique(rowmin0);
    %     rowmin1 = 31; % SEE NOTE ABOVE AND COMMENTED CODE BELOW IF THIS NEEDS
    %     TO BE UPDATED. - 9/19/16 updated rowmin1 calculation - it's now in
    %     the loop
else
    error('The interpolated waveform windows do no share the same minimum sample point! Stopped in beta_cluster_evaluation.m');
end

%% GET DATA
% Determine what to use to construt p, either channel.mat files (set use_chanfiles = 1) or
% CELL_BACKUP variable (set use_chanfiles = 0)
use_chanfiles = 0;
showsavedscore = 0;
if use_chanfiles
    % if using chanfiles, point to dir and index files
    fprintf('Will use channel.mat files to construct p. Point me to the rigth directory!\n');
    pause(.4);
    sDir = uigetdir;
    chanfiles = dir([sDir filesep '*channel_*']);
    
    
    pfile = dir([sDir filesep '*sorting_p.mat']);
    if ~isempty(pfile)
        showsavedscore = 1;
        % load the saved p matrix
        load([sDir filesep pfile.name]); % this will put p_save in workspace
    end
else
    % otherwise, load CELL_BACKUP variable (if not already loaded)
    if exist('CELL_BACKUP','var') || exist('CELL','var')
        fprintf('Found CELL variable. Constructin p based on that.\n');
    else
        fprintf('CELL_BACKUP variable not found! Please point me to it.\n');
        pause(.4);
        [cF, cD] = uigetfile;
        c_loader = load([cD cF]);
        CELL_BACKUP = c_loader.CELL_BACKUP;
    end
    nclusters = 0;
    for qq = 1:length(CELL_BACKUP)
        if ~isempty(CELL_BACKUP(qq).block)
            nclusters = nclusters + size(CELL_BACKUP(qq).block.clust,2);
        end
    end
end

if ~exist('animal','var')
    animal = input('Which animal is this?  ','s');
else
    fprintf('Animal: %s\n',animal);
end


%% MAIN LOOP
p = []; score= []; scount = 0; pcount = 0; chcount = 0;
for ii =  1:31 % size(chanfiles,1)
    
    % based on use_chanfiles, figure out how to check for empty channels
    if use_chanfiles
        chanfilename = [sDir filesep 'channel_' num2str(ii) '.mat'];
        checker = exist(chanfilename,'file');
    else
        checker = ~isempty(CELL_BACKUP(ii).block);
    end
    
    if checker
        chcount = chcount +1;
        chan = ii;
        
        % based on use_chanfiles, load data
        if use_chanfiles
            % if using files, load each one
            fprintf('Loading file %s.\n',chanfilename);
            chan_loader = load(chanfilename);
            chdata = chan_loader.chdata;
        else
            % otherwise, chdata is simply the iith entry of CELL_BACKUP
            chdata = CELL_BACKUP(ii);
        end
        
        % load the block variable into chdat
        if iscell(chdata.block)
            chdat = chdata.block{1};
        else
            chdat = chdata.block(1);
        end
        
        % go through all clusters to find last spike time and infer
        % experiment end hour
        for gg = 1:size(chdat.clust,2)
            last_t(gg) = chdat.clust(gg).time(end);
        end
        max_hour = round(max(last_t./3600),0);
        
        % loop through clusters
        for ee = 1:size(chdat.clust,2)
            
            % this is data for current cluster
            this = chdat.clust(ee);
                
            
            dat = [];
            dat = this.time;
            X = diff(dat);
            h = histc(X,[0:0.001:10.0]); % previously calculated only out to 0.5 sec
            % This changed on 3/27/17 to 10 sec. This prevents over-representation
            % of noise in <3ms for slowly firing cells.
            
            if sum(h)<1000
                % if cluster has very few spikes, autmoatically assign
                % quality=4
                pcount = pcount + 1;
                p(pcount,1:(ncol-1)) = 0;
                p(pcount,ncol)       = 4;
                
                % and score reflects 100% confidence in q=4
                % move scount along to allow referencing of score to CELL
                scount = scount + 1;
                score(scount,1) = 4;
                score(scount,2) = 1.0;
                score(scount,3:4) = 0;
            else
                % otherwise go through calcuylation
                pcount = pcount + 1;

                % --------------------------------------------
                
                % check to make sure that the minimums are in the same sample:
                trace = this.meantrace';
                normtrace = trace/abs(min(trace));
                [~,rowmin1] = min(trace(min_rng,:));
                rowmin1 = rowmin1 + min_rng(1) - 1;
                if numel(unique(rowmin1))~=1
                    disp('Issue with finding minimums. Check me out.');
                    minmode = mode(rowmin1);
                    others = rowmin1(rowmin1~=minmode);
                    
                    if abs(minmode - others) == 1
                        rowmin1 = minmode;
                    else
                        % if not, one possibility is that a WF has a second minimum but
                        % that its "main" minimum is still good.
                        % check for this.
                        bad_ones_idx = find(abs(minmode-others)>1); % find "bad" WFs
                        if numel(bad_ones_idx) == 1 && bad_ones_idx == 1
                            fprintf('Bad alignment on cluster 1 - noise cluster! Ignoring.');
                            
                            rowmin1 = minmode;
                        else
                            for qq=1:length(bad_ones_idx)
                                % for each one, invert the data
                                tmp_badwf = meanwfs(:,bad_ones_idx(qq));
                                tmp_badwf_inv = 1.01*max(tmp_badwf) - tmp_badwf;
                                % then find the indices of the peaks (this will give us the
                                % indices of the minima since the data is inverted)
                                try
                                    [~,min_bad] = findpeaks(double(tmp_badwf_inv));
                                catch
                                    error('autoClusterEval failed at findpeaks.');
                                end
                                % now check if any of those peaks are in the right place
                                min_bad_diff = abs(min_bad - minmode) <= 4;
                                if any(min_bad_diff)
                                    % if there are any, substitute them in the rowmin1
                                    % array in place of the old value for that WF
                                    rowmin1(bad_ones_idx(qq)) = min_bad(min_bad_diff);
                                end
                            end
                            % go through the check again: if all WFs have restored via the
                            % loop above, the code will pass this check and rowmin1 will be
                            % ok. Otherwise, it will give a warning and hit a keyboard.
                            newminmode = mode(rowmin1);
                            newothers = rowmin1(rowmin1~=minmode);
                            if all(abs(newminmode - newothers) <= 4);
                                rowmin1 = newminmode;
                            else
                                disp('The mean normalized waveforms do no share the same minimum sample point! Warning in beta_cluster_evaluation.m');
                                rowmin1 = newminmode;
                            end
                        end
                    end
                else
                    rowmin1 = unique(rowmin1);
                end
                
                % get rangediff
                %rangediff = rowmin0-rowmin1;
                
                % - - - - - - - - - - - - - - - - - - - - -
                clear kvfs kvrsu firingrate sigdiff crv toofast isiflag autodecision
                clear rise clustvar
                
                % use this as a metric of how close the curves are to the window
                % means? is <0.2 too stringent? some combination with the
                % DISTANCE metric you're generating looks pretty good. play
                % with this and figure out a good threshold - maybe 10?
                
                [~,~,kvfs]  = kstest2(FSwin2,normtrace);
                
                [~,~,kvrsu] = kstest2(RSUwin2',normtrace);
                
                
                fstempdist  = sum( abs( FSwin2' - normtrace ));
                rsutempdist = sum( abs( RSUwin2' - normtrace  ));
                %
                %                 xcrsu       = xcorr(normtrace, RSUwin2',0,'coeff');
                %                 xcfs        = xcorr(normtrace, FSwin2',0,'coeff');
                
                chirsu      = sum(((normtrace - RSUwin2').^2)./abs(RSUwin2'));
                chifs       = sum(((normtrace - FSwin2').^2)./abs(FSwin2'));
                
                sigdiff     = min([fstempdist rsutempdist]);
                
                % fit the rising phase of the mean trace
                rising_0 = normtrace(min_rng(1):rowmin1);
                
                % now upsamples this  rising phase to 100 samples. This way,
                % for whatever size of waveform coming in, the calculated slope
                % will be comparable. This will allow us to set a single
                % threshold for distinguishing between good and bad units based
                % on the slope of the rising phase.
                X = 1:size(rising_0,1);
                try
                    Xq = linspace(X(1),X(end),100);
                catch
                    error('autoClusterEval failed at Xq calculation.');
                end
                rising = interp1(X, rising_0, Xq);
                %             D = 10:bot(uu);
                D = 1:100;
                pf = polyfit(D,rising,1);
                f = polyval(pf,D);
                rise = pf(1);
                
                % slope of decay phase
                x_decay0 = rowmin1+1;
                x_decay1 = x_decay0 + length(min_rng) - 1;
                decay_0 = normtrace(x_decay0:x_decay1);
                X2 = 1:size(decay_0,1);
                try
                    Xq2 = linspace(X2(1),X2(end),100);
                catch
                    error('autoClusterEval failed at Xq2 calculation.');
                end
                decay = interp1(X2, decay_0, Xq2);
                %             D = 10:bot(uu);
                D2 = 1:100;
                pf2 = polyfit(D2,decay,1);
                f2 = polyval(pf2,D2);
                decay_slope = pf2(1);
                
                % This calculates the changes in the second derivative of the
                % rebound of the spike (min to max). Here, we're searching for
                % (more than 1) local maxima in the rate of change, as a spike
                % from a "real" cell should only show 1 or 0 depending on the
                % curvature of the rebound.
                
                hipoint     = rowmin1 + (find( normtrace(rowmin1:end) == max(normtrace(rowmin1:end) ) )) -1;
                wftest      = normtrace(rowmin1:hipoint); % old way - looks only from min to max
                allwftest   = normtrace(rowmin1:end);
                xs          = 1:size(wftest);
                
                % Check for minima/maxima and inflection points. Assumption is that
                % there shouldn't be more than 1 local minimum/maximum oir
                % inflection points between absolute minimum of waveform and the
                %  following maximum, it it comes from a single-unit.
                
                % clear variables
                a = []; b= []; aa = [];
                a0 = []; b0 = []; aa0 = [];
                c1 = []; c2 = []; cc1 = [];
                c1_crv = []; c2_crv = []; cc1_crv = [];
                
                
                a       = diff(wftest);
                b       = diff(a);
                aa 		= diff(allwftest);
                
                % check that there are no more than one actual minima/maxima
                try
                    
                    a0(:,1) = a;
                    aa0(:,1) = aa;
                catch
                    disp('autoClusterEval failed at a0 calculation.');
                    keyboard
                end
                a0(:,2) = [a(2:end); NaN];
                c1(:,1) = sum(a0(:,1)>=0 & a0(:,2)<0);
                c1(:,2) = sum(a0(:,1)<0 & a0(:,2)>=0);
                c1_crv  = sum(c1(:)) > 1;
                
                %whole waveform
                aa0(:,2) = [aa(2:end); NaN];
                cc1(:,1) = sum(aa0(:,1)>=0 & aa0(:,2)<0);
                cc1(:,2) = sum(aa0(:,1)<0 & aa0(:,2)>=0);
                cc1_crv = sum(cc1(:)) > 1;
                
                % check that there is no more than one inflection point
                try
                    b0(:,1) = b;
                catch
                    keyboard
                end
                b0(:,2) = [b(2:end); NaN];
                c2(:,1) = sum(b0(:,1)>=0 & b0(:,2)<0);
                c2(:,2) = sum(b0(:,1)<0 & b0(:,2)>=0);
                c2_crv  = sum(c2(:)) > 1;
                
                crv = any([c1_crv c2_crv]);
                if crv == 0
                    crv = cc1_crv;
                end
                
                % GAMMA function fit to ISI dist and contamination calc
                % --------------------------------------------

                out_to = 100; % look at gamma fit for isi dist out to this many ms
                p(pcount,1:2) = gamfit(1:100,[],[],h(1:100)');
                ygam2 = gampdf(1:100,p(pcount,1),p(pcount,2));
                
                % compute contamination over whole ISI hsitogram (0:10 sec)
                contamination = 100*(sum(h(1:2))/sum(h(1:end)));
                meanfr = numel(this.time)/(this.time(end)-this.time(1));
                
                % Also use a different measure of contamination: get mean
                % noise contamination per hour in the middle 50% of
                % recording (25% to 75% of recording time) and then average
                % the hourly means together to get an overall value
                frac = 0.10; % fraction of experiment to exclude on either side
                % e.g. 0.25 excludes first and last 25%
                edges = floor(.1*max_hour) : ceil((1-frac)*max_hour); % 25% to 75% of experiment time
                for qq = 1:(size(edges,2)-1)
                    this_hour_spikes = this.time(this.time >= edges(qq)*3600 & ...
                        this.time <= edges(qq+1)*3600);
                    % only get contamination if cell is firing at more than 0.005
                    % Hz in that hour (more than 18 spikes in that hour)
                    if numel(this_hour_spikes) > 0.005*3600
                        this_hour_isi = diff(this_hour_spikes);
                        this_hour_h = histc(this_hour_isi,[0:0.001:10.0]);
                        hourly_cont(qq) = 100*sum(this_hour_h(1:2))/sum(this_hour_h(1:end));
                    else
                        hourly_cont(qq) = NaN;
                    end
                    this_hour_spikes = [];
                    this_hour_isi = [];
                    this_hour_h = [];
                end
                contamination_byhour = nanmean(hourly_cont);
             
                
                % - - - - - - - - - - - - - - - - - - - - - - -
                % - - - Check to see if there's 60 cycle noise in the spike
                % times - - - - - - - - - - - - - - - - - - - -
                Fs = 1000;           % Sampling frequency
                L = 500;             % Length of signal
                fhz = Fs*(0:(L/2))/L;
                yfft = fft(h(1:end-1));
                P2 = abs(yfft/L);
                P1 = P2(1:L/2+1);
                P1(2:end-1) = 2*P1(2:end-1);
                
                P1norm = P1/sum(P1);
                %plot(fhz,P1norm)
                peaktemplate = [0.3868 0.4064 0.4220 0.4371 0.4596 0.4859 0.5210...
                    0.5651 0.6203 0.6875 0.7944 1.0000 0.7902 0.6830 0.6138 ...
                    0.5590 0.5119 0.4745 0.4461 0.4201 0.4005]';
                
                ac_noise  = corr(peaktemplate, P1(20:40)/max(P1(20:40)) );
                
                %ac_noise = sum(P1norm(30:33));
                % - - - - - - - - - - - - - - - - - - - - - - -
                
                % get neg_pos_time as one of the parameters
                % not this assumes a sampling_rate of 24414.06 Hz
                negpostime = calc_negpostime(normtrace); % returns negpostime in milliseconds
                pfs_nptime = calc_negpostime(FSwin2);
                rsu_nptime = calc_negpostime(RSUwin2);
                
                % calc halfwidth in ms
                halfwidth = calc_halfwidth(normtrace);
                if isempty(halfwidth), halfwidth = 0; end
                pfs_hw = calc_halfwidth(FSwin2);
                rsu_hw = calc_halfwidth(RSUwin2);
                
                %tailslope
                tailslope = calc_tailslope(normtrace);
                pfs_ts = calc_tailslope(FSwin2);
                rsu_ts = calc_tailslope(RSUwin2);
                
                %% MAKE P
                p(pcount,3) = contamination;
                p(pcount,4) = contamination_byhour;
                p(pcount,5) = meanfr;
                p(pcount,6) = kvfs;
                p(pcount,7) = kvrsu;
                p(pcount,8) = fstempdist;
                p(pcount,9) = rsutempdist;
                p(pcount,10) = chirsu;
                p(pcount,11) = chifs;
                p(pcount,12) = rise;
                p(pcount,13) = decay_slope;
                p(pcount,14) = crv;
                p(pcount,15) = min(this.meantrace);
                p(pcount,16) = negpostime;
                p(pcount,17) = halfwidth;
                p(pcount,18) = tailslope;
                p(pcount,19) = ac_noise;
                
                Y = p(pcount,[1:(ncol-1)]);
                
                % Run measurements through the model and predict the
                % quality:
                
                clear label;
                scount = scount+1;
                
                model_scored = [];
                if size(Mdl.X,2) == size(Y,2)
                    model_scored = 1;
                    [label,scoretmp,~] = predict(Mdl,Y); % use model to score
                    label = str2double (label{1});
                    score(scount,1) = label;
                    score(scount,2) = scoretmp(label);
                    scoretmp2 = scoretmp;
                    scoretmp2(label) = -100;
                    [scoremax,~] = max(scoretmp2);
                    [~,scoremaxlabel] = find(scoretmp2 == scoremax);
                    if numel(scoremaxlabel)==1
                        score(scount,3) = scoremaxlabel;
                        score(scount,4) = scoremax;
                    else
                        score(scount,3) = 0;
                        score(scount,4) = 0;
                    end
                else
                    model_scored = 0;
                    fprintf('Model Mdl and scoring variable Y have different sizes. Assign quality manually.\n\n');
                    fprintf('  *** Previous quality: %u\n\n',this.quality);
                end
                
                %                 try
                %                     [label,scoretmp,~] = predict(Mdl,Y); % use model to score
                %                 catch ME
                %                     if strcmp(ME.identifier,'stats:CompactTreeBagger:predictAccum:XVarNamesSizeMismatch')
                %                         cause_msg = sprintf(['You are trying to score a variable Y that has %u columns using',...
                %                             ' a model Mdl with %u columns in its training data.\n',...
                %                             'You need to retrain this model first.'],size(Y,2),size(Mdl.X,2));
                %                         causeException = MException('MATLAB:train_new_Mdl:badModel',cause_msg);
                %                         ME = addCause(ME,causeException);
                %                     end
                %                     rethrow(ME);
                %                 end
                
                skip_4s = 0;
                if model_scored
                    assigned_qual = label;
                else
                    assigned_qual = this.quality;
                end
                
                if skip_4s && assigned_qual == 4
                    
                    fprintf('Cluster %u assigned quality is 4. Skipping.\n\n',pcount);
                    p(pcount,ncol) = assigned_qual;
                    
                    
                else
                    plot_and_score = 1;
                    
                    
                    
                    % QUALITY ASSIGNMENT:
                    if model_scored
                        p(pcount,ncol) = label;
                    elseif ~model_scored && ~plot_and_score
                        p(pcount,ncol) = this.quality;
                    end
                    
                    %% MAKE PLOTS
                    
                    
                    if plot_and_score
                        
                        g = figure(pcount);
                        
                        set(g,'units','normalized','position',[0.05 0.05 0.85 0.85]);
                        s1 = subplot(1,3,1);
                        plot(ygam2,'linewidth',2,'color',[.65 .09 .18]);
                        yvals = get(gca,'ylim');
                        if model_scored
                            txt = text(35,yvals(2)*0.75,['Predicted q: ' num2str(label) ]);
                        else
                            txt = text(35,yvals(2)*0.75,['Previous q: ' num2str(this.quality) ]);
                        end
                        if showsavedscore
                            txt2 = text(35,yvals(2)*0.65,['Saved q: ' num2str(p_save{chcount}(pcount,15))]);
                            set(txt2,'fontsize',14);
                        end
                        txt3 = text(35, yvals(2)*.9,sprintf('Cluster %u out of %u',pcount,nclusters),'fontsize',14);
                        set(txt,'fontsize',14);
                        set(gca,'fontsize',14);
                        title(['Param1: ' num2str(p(pcount,1)) '; Param2: ' num2str(p(pcount,2)) ],'fontsize',14);
                        
                        subplot(1,3,2);
                        bar(h(1:100));
                        title(sprintf('Total %% cont: %.3f\nMean hourly %% cont: %.3f',contamination,contamination_byhour),'fontsize',14);
                        yvs = get(gca,'ylim');
                        fftxt = text(40,0.5*yvs(2),['60 Hz noise: ' num2str(ac_noise) '.']);
                        set(fftxt,'fontsize',14);
                        set(gca,'xlim',[0 100],'xtick',0:20:100,'xticklabel',0:20:100,'fontsize',14);
                        xlabel('ISI (ms)','fontsize',14);
                        
                        subplot(1,3,3);
                        plot(normtrace,'k','linewidth',3), hold on
                        plot(RSUwin2,'color',[0 0.45 0.74],'linewidth',2);
                        plot(FSwin2,'color',[0.90 0.69 0.13],'linewidth',2);
                        text(61,-0.84,sprintf('pFS nptime: %.4f',pfs_nptime),'fontsize',10);
                        text(61,-0.86,sprintf('RSU nptime: %.4f',rsu_nptime),'fontsize',10);
                        text(61,-0.89,sprintf('pFS hw: %.4f',pfs_hw),'fontsize',10);
                        text(61,-0.91,sprintf('RSU hw: %.4f',rsu_hw),'fontsize',10);
                        text(61,-0.94,sprintf('pFS ts: %.4f',pfs_ts),'fontsize',10);
                        text(61,-0.96,sprintf('RSU ts: %.4f',rsu_ts),'fontsize',10);
                        text(37,-0.50,['Amp ' num2str(min(this.meantrace)) ' uV'],'fontsize',12);
                        text(37,-0.55,sprintf('Rise slope: %.4f',rise),'fontsize',12);
                        text(37,-0.60,sprintf('Decay slope: %.4f',decay_slope),'fontsize',12);
                        text(37,-0.65,sprintf('KS rsu: %.4f',kvrsu),'fontsize',12);
                        text(37,-0.70,sprintf('KS pfs: %.4f',kvfs),'fontsize',12);
                        text(37,-0.75,sprintf('dist rsu: %.4f',rsutempdist),'fontsize',12);
                        text(37,-0.80,sprintf('dist pfs: %.4f',fstempdist),'fontsize',12);
                        text(37,-0.85,sprintf('neg-pos time: %.4f ms',negpostime),'fontsize',12);
                        text(37,-0.90,sprintf('halfwidth: %.4f ms',halfwidth),'fontsize',12);
                        text(37,-0.95,sprintf('tailslope: %.4f',tailslope),'fontsize',12);
                        title(['Rate: ' num2str(meanfr) ' Hz'],'fontsize',14);
                        curr_y = get(gca,'ylim');
                        %                     set(gca,'ylim',[-1 max(curr_y(2),0.65)]);
                        set(gca,'fontsize',14);
                        
                        % The user hits [space] to accept the quality and move on,
                        % or press the corrected quality number key to update the
                        % record and then move on.
                        [x,y,button] = ginput(1);
                        
                        zztop = 0;
                        while zztop == 0
                            if strcmp(char(button),' ')
                                % go to next cell;
                                if ~model_scored
                                    p(pcount,ncol) = this.quality;
                                end
                                txt2 = text(s1,35,yvals(2)*0.65,...
                                    ['Keeping q: ' num2str(p(pcount,ncol))],...
                                    'color',[.47 .67 .20],'fontsize',14);
                                pause(0.75);
                                
                                zztop = 1;
                                close (g);
                            elseif any(strcmp(char(button),{'1', '2', '3', '4'}))
                                % Enter an edited quality and update the p entry
                                p(pcount,ncol) = str2double(char(button));
                                txt2 = text(s1,35,yvals(2)*0.65,...
                                    ['Newly assigned q: ' num2str(p(pcount,ncol))],...
                                    'color',[.47 .67 .20],'fontsize',14);
                                pause(0.75);
                                zztop = 1;
                                close (g);
                                
                            elseif ~any(strcmp(char(button),{'1', '2', '3', '4', ' '}))
                                disp('Try again ');
                                [x,y,button] = ginput(1);
                            end
                        end
                        
                    end
                end
            end
        end
    end
    
end

save_p = 1;
if save_p
    if ~exist('pDir','var')
        pDir =  uigetdir(cd,'Pick the directory where the new p variables will be saved.');
    end
    pFile = [animal '_p.mat'];
    save([pDir filesep pFile],'p');
end

