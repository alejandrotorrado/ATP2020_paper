function datOut = autoClusterEval(block,interp_samples,trim_rng,min_rng)

if nargin == 3
    min_rng = 1:size(block.clust(1).meantrace,2);
end

% add quality field:


for uu = 1:size(block.clust,2);
    
    block.clust(uu).quality = [];
    nsamples(uu) = size(block.clust(uu).meantrace,2);
    
end
nsamp = unique(nsamples);
fprintf('These waveforms have %u samples (trimmed down from %u).\n',nsamp,interp_samples);
if length(nsamp)~=1
    error('Problem with number of waveform samples in AutoClusterEval!');
end

% back up the original data for now.
OGdat = block;


% - - - - - - - - - - - - - - - - - - - - - - - -

% these are the minimum centered and normalized mean waveforms that
% correspond to FS and RSU traces, respectively. These are hard coded here
% for convienence. The code to regenerate these is found in
% 'voltage_windowstart.m', which is stored locally on KBH's laptop.
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
    trim_start = rowmin0 + trim_rng(1);
    trim_end = rowmin0 + trim_rng(end);
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
% for cc = chan;

if verbose; disp(['Starting automated cluster evaluation.']); end

for bb = 1:size(block,2);
    
    %         if verbose; disp(['Block ' num2str(bb) '.']); end
    
    
    % check to see if clusters need to be merged as a result of
    % over-clustering:
    passback = overclusteringcheck(block.clust);
    block.clust = passback;
    
    % test to see if the normalized waveform passes through the voltage
    % window:
    %
    % extract and reshape the mean traces:
    meanwfs = cell2mat({block.clust.meantrace});
    meanwfs = reshape(meanwfs,nsamp,[]);
    % normalize the mean traces:
    wfmins      = min(meanwfs,[],1);
    wfminsmtrx  = repmat(wfmins,nsamp,1);
    normwfs     = meanwfs./abs(wfminsmtrx);
    
    
    % check to make sure that the minimums are in the same sample:
    if numel(unique(rowmin0)) == 1;
        rowmin0 = unique(rowmin0);
        [~,rowmin1] = min(normwfs(min_rng,:));
        rowmin1 = rowmin1 + min_rng(1) - 1;
        if numel(unique(rowmin1))~=1;
            
            minmode = mode(rowmin1);
            others = rowmin1(rowmin1~=minmode);
            
            if abs(minmode - others) == 1;
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
        %rowmin1 = 41; % SEE NOTE ABOVE AND COMMENTED CODE BELOW IF THIS NEEDS TO BE UPDATED.
    else
        disp('The interpolated waveform windows do no share the same minimum sample point! Stopped in beta_cluster_evaluation.m');
    end
    
    % get rangediff
    rangediff = rowmin0-rowmin1;
    
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Look at the data... for debug and purposes of being critical.
    if showoff == 1;
        %             figure(1)
        %             plot(normwfs(range + rowmin1,:)); hold on
        %             plot (RSUwin2(range+rowmin0),'r--','linewidth',1.5);
        %             plot (FSwin2(range+rowmin0),'k--','linewidth',1.5);
        %             hold off
        
        figure(2), hold on
        try
            for uu = 1:size(block.clust,2);
                shadedErrorBar(1:nsamp,block.clust(uu).meantrace,...
                    block.clust(uu).tracestdev,'-k',1);
                plot(block.clust(uu).meantrace,'linewidth',2);
            end
        catch
            disp('caught')
            keyboard
        end
        hold off
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    % find indices of minimum point in spike. this should be the same
    % across all WFs as they're peak aligned, but keep it flexible for
    % other approaches
    [~,bot] = min(normwfs(min_rng,:));
    bot = bot + min_rng(1) - 1;
    
    
    % - - - - - - - - - - - - - - - - - - - - -
    clear kvfs kvrsu firingrate sigdiff crv toofast isiflag autodecision
    clear rise clustvar
    for uu = 1:size(normwfs,2);
        
        if verbose; disp(['Cluster ' num2str(uu) '.']); end
        
        % use this as a metric of how close the curves are to the window
        % means? is <0.2 too stringent? some combination with the
        % DISTANCE metric you're generating looks pretty good. play
        % with this and figure out a good threshold - maybe 10?
        try
            [~,~,kvfs(uu)]  = kstest2(FSwin2(rangediff+1:end)',normwfs(1:end-rangediff,uu));
        catch
            error('autoClusterEval failed at kvfs calculation.');
	end
        [~,~,kvrsu(uu)] = kstest2(RSUwin2(rangediff+1:end)',normwfs(1:end-rangediff,uu));
        fstempdist      = sum(abs(FSwin2(rangediff+1:end)' - normwfs(1:end-rangediff,uu)));
        rsutempdist     = sum(abs(RSUwin2(rangediff+1:end)' - normwfs(1:end-rangediff,uu)));
        
        sigdiff(uu)     = min([fstempdist rsutempdist]);
        
%         if showoff;
%             figure(100+uu);
%             plot(FSwin2(rangediff+1:end),'k--','linewidth',1.5); hold on
%             plot(RSUwin2(rangediff+1:end),'b--','linewidth',1.5);
%             plot(normwfs(1:end-rangediff,uu));
%             uiwait(gcf);
%         end
        % fit the rising phase of the mean trace
        rising_0 = normwfs(10:bot(uu),uu);
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
        %plot(D,rising,'.b')
        %             D = 10:bot(uu);
        D = 1:100;
        p = polyfit(D,rising,1);
        f = polyval(p,D);
        rise(uu) = p(1);
        
%         if showoff; plot(D,f,'--r'); uiwait(gcf); end
        
        % This calculates the changes in the second derivative of the
        % rebound of the spike (min to max). Here, we're searching for
        % (more than 1) local maxima in the rate of change, as a spike
        % from a "real" cell should only show 1 or 0 depending on the
        % curvature of the rebound.
        
        hipoint = rowmin1 + (find( normwfs(rowmin1:end,uu) == max(normwfs(rowmin1:end,uu) ) )) -1;
        wftest    = normwfs(rowmin1:hipoint,uu); % old way - looks only from min to max
        allwftest    = normwfs(rowmin1:end,uu);
        xs      = 1:size(wftest);
        
        % clear variables
        a = []; b= []; aa = [];
        a0 = []; b0 = []; aa0 = [];
        c1 = []; c2 = []; cc1 = [];
        c1_crv = []; c2_crv = []; cc1_crv = [];
        
        % Check for minima/maxima and inflection points. Assumption is that
        % there shouldn't be more than 1 local minimum/maximum oir
        % inflection points between absolute minimum of waveform and the
        %  following maximum, it it comes from a single-unit.
        a       = diff(wftest);
        b       = diff(a);
		aa 		= diff(allwftest);
        
        % check that there are no more than one actual minima/maxima 
        try
        	a0(:,1) = a;
			aa0(:,1) = aa;
        catch
	    error('autoClusterEval failed at a0 calculation.');
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

        % check that there are no more than one inflexion point 
        b0(:,1) = b;
        b0(:,2) = [b(2:end); NaN];
        c2(:,1) = sum(b0(:,1)>=0 & b0(:,2)<0);
        c2(:,2) = sum(b0(:,1)<0 & b0(:,2)>=0);
        c2_crv  = sum(c2(:)) > 1;
        
        crv(uu) = any([c1_crv c2_crv]);
		if crv(uu) == 0
			crv(uu) = cc1_crv;
		end
        % this is effectively a binary decision on whether traces can
        % qualify as single units or not. 
        % TLDR: single units are a 0 here.
        
%         if showoff;
%             %                 plot(xs+rowmin1-1, test,'r');
%             plot(xs,wftest,'r');
%             uiwait(gcf);
%             hold off
%         end
        
        % Test the ISI distribution to see if the cluster has too much
        % <3 msec contamination/error
        [~,~,toofast(uu,1)] = isidist(block.clust(uu).time,'cutoff',[]);
        
        clustvar(uu) = block.clust(uu).clustvar;
        
        if showoff;
            figure(100+uu);
            plot(FSwin2(rangediff+1:end),'b--','linewidth',1.5); hold on
            plot(RSUwin2(rangediff+1:end),'r--','linewidth',1.5);
            plot(normwfs(:,uu),'k'); hold on;
            text(32,-0.60,['FS k val is ' num2str(kvfs(uu))]);
            text(32,-0.65,['RSU k val is ' num2str(kvrsu(uu))]);
            text(32,-0.7,['dist is ' num2str(sigdiff(uu))]);
            text(32,-0.85,['Neg curvature is ' num2str(crv(uu))]);
            text(32,-0.8,['ISI error is ' num2str(toofast(uu)*100) '%']);
            text(32,-0.75,['Clustvar is ' num2str(block.clust(uu).clustvar) '.']);
            text(32,-0.5,['Spike amplitude is ' num2str(wfmins(uu))]);
            text(32,-0.55,['Rising phase best fit ' num2str(rise(uu)) '.']);
            pause(0.2);
        end
        
        if kvfs(uu) < kvrsu(uu); % dealing with an FS cell?
            
            if toofast(uu)<0.03;
                isiflag(uu) = 1;
            elseif toofast(uu)>0.03 && toofast(uu)<=0.05;
                isiflag(uu) = 2;
            elseif toofast(uu)>0.05 && toofast(uu)<=0.12;
                isiflag(uu) = 3;
			elseif toofast(uu) > 0.12
				isiflag(uu) = 4;
            end
            
        elseif kvrsu(uu) <= kvfs(uu); % dealing with RSU?
            
            if toofast(uu)<0.015;
                isiflag(uu) = 1;
            elseif toofast(uu)>0.015 && toofast(uu)<=0.03;
                isiflag(uu) = 2;
            elseif toofast(uu)>0.03 && toofast(uu) <= 0.10;
                isiflag(uu) = 3;
			elseif toofast(uu) > 0.10
				isiflag(uu) = 4;
            end
            
        end
        
        % Check the rough firing rate of the cluster:
        elpsdtime       = block.clust(uu).time(end) - block.clust(uu).time(1);
        firingrate(uu)  = numel(block.clust(uu).idx)/elpsdtime;
        
        if showoff; text(32,-0.45,['Firing rate is ' num2str(firingrate(uu))]); end
        
        try
            % Make decisions here? could do a point-total (best keep score
            % low) or do a series of ifelse statements.
            if  firingrate(uu)<0.05  && sigdiff(uu) > 20 ||  clustvar(uu)>3000 &&...
                    sigdiff(uu) > 15 || crv(uu) == 1 && sigdiff(uu)>15 ||...
                    rise(uu)>=-0.0045 || firingrate(uu)<0.05  && clustvar(uu) > 20000;
                autodecision(uu) = 4;
            else
                % - - - - - - - - - - - - - - - -
                if any([kvfs(uu) kvrsu(uu)]<=0.2) && crv(uu)==0 && sigdiff(uu)<=10 &&...
                        isiflag(uu) == 1 && clustvar(uu) < 2000;
                    
                    autodecision(uu) = 1;
                    
                elseif  any([kvfs(uu) kvrsu(uu)]<=0.2) && crv(uu)==0 && sigdiff(uu)<=15 &&...
                        isiflag(uu) <= 2 && clustvar(uu) < 10000 || isiflag(uu) == 1 &&...
                        any([kvfs(uu) kvrsu(uu)]<=0.2) && crv(uu)==0 && sigdiff(uu)<15 ||...
                        any([kvfs(uu) kvrsu(uu)]<=0.2) && isiflag(uu)<3 && rise(uu)<-0.005 && crv(uu)==0 ||...
                        kvfs(uu)<kvrsu(uu) && kvfs(uu)<0.18 && isiflag(uu)<3 && sigdiff(uu)<20 && ...
                        clustvar(uu)<10000 && crv(uu)==0 || ...
						sigdiff(uu)<5 && any([kvfs(uu) kvrsu(uu)]<=0.20) && isiflag(uu)<3 && rise(uu)<-0.005;
                    
                    autodecision(uu) = 2;
                    
                elseif any([kvfs(uu) kvrsu(uu)]<=0.2) && crv(uu)==1 && sigdiff(uu)<=10 &&...
                        isiflag(uu) <= 3 || sigdiff(uu)<20 &&  clustvar(uu)<15000 && crv(uu) == 1 ||...
                        crv(uu) == 1 && sigdiff(uu) < 15 && any([kvfs(uu) kvrsu(uu)]<0.3) && isiflag(uu)<=3 ||...
						sigdiff(uu)<20 && clustvar(uu)<20000 && any([kvfs(uu) kvrsu(uu)]<=0.3) && isiflag(uu)<=3 ||...
                        any([kvfs(uu) kvrsu(uu)]<0.3) && sigdiff(uu)<15 && crv(uu) == 0 ||...
                        any([kvfs(uu) kvrsu(uu)]<0.3) && isiflag(uu)<=3 && crv(uu) == 0 || ...
                        sigdiff(uu)<15 && rise(uu)<-0.004 && isiflag(uu)<=3 && crv(uu)==0 ||...
                        sigdiff(uu) < 20 && crv(uu) == 0 && isiflag(uu) <=3 && clustvar(uu) < 20000 ||...
						clustvar(uu)>20000 && kvfs(uu)<kvrsu(uu) && kvfs(uu)<0.20 && isiflag(uu)>=3 && sigdiff(uu)<10;
						 % FIGURE OUT THREES!!!!!
                    
                    autodecision(uu) = 3;
                    
                elseif all([kvfs(uu) kvrsu(uu)]>0.3) && crv(uu)==1 && sigdiff(uu)>20 &&...
                        isiflag(uu) >= 3 ||  all([kvfs(uu) kvrsu(uu)])>0.25 && sigdiff(uu)>20 &&...
                        isiflag(uu) >= 3 || all([kvfs(uu) kvrsu(uu)]>0.3) && ...
                        clustvar(uu) > 20000 && isiflag(uu) == 4 && rise(uu) > -0.005 || ...
						all([kvfs(uu) kvrsu(uu)]>0.3) && crv(uu) == 1 && sigdiff(uu)>12 ||...
						all([kvfs(uu) kvrsu(uu)]>0.22) && crv(uu) == 1 && isiflag(uu) == 4;
                    
                    autodecision(uu) = 4;
                    
                end
                
                % now correct for amplitude:
                if wfmins(uu) > -30 && autodecision(uu) <=2 ||...
					 all([kvfs(uu) kvrsu(uu)]<0.2) && abs(fstempdist - rsutempdist) <= 4.0 ||...
                        kvfs(uu) == kvrsu(uu);
                    autodecision(uu) = 3;
                end
                
                % - - - - - - - - - - - - - - - -
            end
            disp(['Autoquality ' num2str(autodecision(uu)) '.']);
        catch
            disp('Caught! @ 267. giving quality 5 for this cluster');
	    autodecision(uu) = 5;
            disp(['Autoquality ' num2str(autodecision(uu)) '.']);
        end
        
        try
            if showoff; text(32,-0.4,['Autosort qual is ' num2str(autodecision(uu)) '.']); end
        catch
            disp('I couldn''t make the call...');
        end
        
        
        block.clust(uu).quality = autodecision(uu);
        
        if showoff, uiwait(gcf); end
    end
    
    %         % USE THIS CODE TO ESTABLISH ROWMIN1 ALGORITHMICALLY - IT'S A
    %         WASTE OF TIME TO DO THIS IN EVERY LOOP SO IT'S HARDCODED FOR NOW.
    %         [~,rowmin1] = min(normwfs);
    %         if numel(unique(rowmin1))~=1;
    %             warning('The mean normalized waveforms do no share the same minimum sample point! Stopped in beta_cluster_evaluation.m');
    %             keyboard
    %         else
    %             rowmin1 = unique(rowmin1);
    %         end
    
end



datOut = block;

% now save the output:
% [vFile, vDir] = uiputfile('.mat','Save your data!');
% save([vDir vFile],'channel','-v7.3');




