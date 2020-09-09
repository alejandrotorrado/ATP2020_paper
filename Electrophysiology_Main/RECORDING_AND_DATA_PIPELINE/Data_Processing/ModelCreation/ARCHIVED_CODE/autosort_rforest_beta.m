
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Use this is as a set up to run some of the more advanced descriptive
%   statistics on the clusters for later passing to a randomforest model
%   (treebagger) for autoclassification.
%
nsamp = 59;
interp_samples = 97;
swprng = -18:40;
min_rng = 10:30;


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
if numel(unique(rowmin0)) == 1;
    rowmin0 = unique(rowmin0);
    %     rowmin1 = 31; % SEE NOTE ABOVE AND COMMENTED CODE BELOW IF THIS NEEDS
    %     TO BE UPDATED. - 9/19/16 updated rowmin1 calculation - it's now in
    %     the loop
else
    error('The interpolated waveform windows do no share the same minimum sample point! Stopped in beta_cluster_evaluation.m');
end

% this sets up the range of points around the minimum to be compared in the
% voltage window. this assumes that rowmin1, e.g. the meantrace from each
% cluster, has fewer sample points than the RSU and FS means (as a result
% of how the clusters are generated at the moment)
% nsamp = size(block.clust(1).meantrace,2); % no need for this here -
% calculated above
% range = -(rowmin1-1):(nsamp - rowmin1); % no need for this now that we
% have updated the rowmin1 calculation
showoff = 0;
verbose = 1;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% fitting Gamma dist to ISI distribution
warning('off');
close all;
p = [];


count = 0;
for chan =  [1:7 9:15 17:23 25:31]
    
    if exist(['/Users/khengen/Desktop/PROFILESETUP/channel_' num2str(chan) '.mat'],'file') == 2
        
        load(['/Users/khengen/Desktop/PROFILESETUP/channel_' num2str(chan) '.mat'])
        
        
        for ee = 1:size(chdata.block{1}.clust,2)
            
            this = chdata.block{1}.clust(ee);
            
            dat = [];
            dat = this.time;
            X = diff(dat);
            h = histc(X,[0:0.001:0.5]);
            
            if sum(h)<1000
                % do nothing
            else
                count = count + 1;
                disp(count);
                % --------------------------------------------
                
                % check to make sure that the minimums are in the same sample:
                    trace = this.meantrace';
                    normtrace = trace/abs(min(trace));
                    [~,rowmin1] = min(trace(min_rng,:));
                    rowmin1 = rowmin1 + min_rng(1) - 1;
                    if numel(unique(rowmin1))~=1
                        disp('Issue with finding minimums. Check me out.');
                        keyboard
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
                rising_0 = normtrace(10:rowmin1);
                
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

                % --------------------------------------------
   
                
                %params(count,1:2) = gamfit(h);
                %ygam = gampdf(1:100,params(count,1),params(count,2));
                p(count,1:2) = gamfit(1:100,[],[],h(1:100)');
                ygam2 = gampdf(1:100,p(count,1),p(count,2));
                
                g = figure(count);
                set(g,'name',['Channel ' num2str(chan) ' Count ' num2str(count)]);
                set(g,'units','normalized','position',[0.1 0.1 0.9 0.9]);
                subplot(1,3,1);
                plot(ygam2);
                contamination = 100*(sum(h(1:2))/sum(h(1:500)));
                p(count,3) = contamination;
                
                title(['% Cont: ' num2str(contamination) ]);
                
                subplot(1,3,2);
                bar(h(1:100));
                title(['Param1: ' num2str(p(count,1)) '; Param2: ' num2str(p(count,2)) ]);
                subplot(1,3,3);
                plot(normtrace), hold on
                plot(RSUwin2,'g');
                plot(FSwin2,'r');
                meanfr = numel(chdata.block{1}.clust(ee).time)/(chdata.block{1}.clust(ee).time(end)-chdata.block{1}.clust(ee).time(1));
                p(count,4) = meanfr;
                
                title(['Rate: ' num2str(meanfr) ' Hz']);
   
                p(count,5) = kvfs;
                p(count,6) = kvrsu;
                p(count,7) = fstempdist;
                p(count,8) = rsutempdist;
                p(count,9) = chirsu;
                p(count,10) = chifs;
                p(count,11) = rise;
                p(count,12) = crv;
                
            
                uiwait (gcf) ;
                
                commandwindow;
                q = input('Quality? ','s');
                
                if any(strcmp(q,{'1','2','3','4'}))
                    q = str2double(q);
                end
                
                while ~any(q == [1 2 3 4])
                    q = input('try again...','s');
                    q = str2double(q);
                end
                
                p(count,13) = q;
                chdata.block{1}.clust(ee).manualQ = q;
                
            end
        end
    else
    end
    
end

[sFile sDir] = uiputfile(cd,'save p');
save([sDir sFile],'p');


%
q1 = p(:,13) == 1;
q2 = p(:,13) == 2;
q3 = p(:,13) == 3;
q4 = p(:,13) == 4;

figure(1000), hold on
scatter(p(q1,1), p(q1,2),'g')
scatter(p(q2,1), p(q2,2),'y')
scatter(p(q3,1), p(q3,2),'r')
scatter(p(q4,1), p(q4,2), [], [0.3 0.3 0.3])


% P is organized as follows:
%   Column 1 is the gamma fit value 1
%   Column 2 is the gamma fit value 2
%   Column 3 is the percent contamination under 2 msec
%   Column 4 is the mean firing rate (not taking into account on/off times)
%   Column 5 is KV FS
%   Column 6 is KV RSU
%   Column 7 is FStempdist
%   Column 8 is RSUtempdist
%   Column 9 is chi rsu
%   Column 10 is chi FS
%   Column 11 is rise
%   Column 12 is curve (inflection)
%   Column 13 is subjective quality






figure(1000), hold on
scatter3(p(q1,1), p(q1,2), p(q1,4), 'g')
scatter3(p(q2,1), p(q2,2), p(q2,4), 'y')
scatter3(p(q3,1), p(q3,2), p(q3,4),'r')
scatter3(p(q4,1), p(q4,2), p(q4,4),[], [0.3 0.3 0.3])
xlabel('Gamma fit 1')
ylabel('Gamma fit 2')
zlabel('% Contamination')