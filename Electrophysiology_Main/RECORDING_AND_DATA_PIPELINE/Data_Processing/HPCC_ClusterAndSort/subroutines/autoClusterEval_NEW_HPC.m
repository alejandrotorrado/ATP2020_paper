function datOut = autoClusterEval_NEW_HPC(block,interp_samples,trim_rng,min_rng)

%% UPDATED 01/17/2017 by ATP
% This version now uses the up-to-date autosorting algorithm developed by
% Alejandro Torrado Pacheco, Keith B Hengen, Mara CP Rue.
% Scored against expert human cell sorters (ATP & KBH) the algorithm has a
% ~2% error rate on average. At single-/multi-unit border, the disagreement
% is ~ 15%, with equal fractions of false positives and false negatives
% (not biased towards a conservative or lax estimate of what constitutes a
% single-unit cluster).
%
% Version _HPC is for use on computing cluster.


%% SETUP


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
        fstempdist(uu)      = sum(abs(FSwin2(rangediff+1:end)' - normwfs(1:end-rangediff,uu)));
        rsutempdist(uu)     = sum(abs(RSUwin2(rangediff+1:end)' - normwfs(1:end-rangediff,uu)));
        
        xcrsu(uu) = xcorr(normwfs(1:end-rangediff,uu),RSUwin2(rangediff+1:end),0,'coeff');
        xcfs(uu) = xcorr(normwfs(1:end-rangediff,uu),FSwin2(rangediff+1:end),0,'coeff');
        
        chirsu(uu) = sum(((normwfs(1:end-rangediff,uu)-RSUwin2(rangediff+1:end)').^2)./abs(RSUwin2(rangediff+1:end)'));
        chifs(uu) = sum(((normwfs(1:end-rangediff,uu)-FSwin2(rangediff+1:end)').^2)./abs(FSwin2(rangediff+1:end)'));
        
        sigdiff(uu)     = min([fstempdist(uu) rsutempdist(uu)]);
        
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
        isi = diff(block.clust(uu).time);
        
        isi_edges = 0:0.001:2;
        [~,isiperbin,toofast3(uu,1)] = isidist(block.clust(uu).time,...
            'cutoff',[],'edges',isi_edges);
        
        toofast2(uu,1) = sum(isi<0.002)/numel(isi);
        toofast1(uu,1) = sum(isi<0.001)/numel(isi);
        
        isi1ms(uu) = sum(isi<=0.001);
        isi2ms(uu) = sum(isi<=0.002)- sum(isi<=0.001);
        isi3ms(uu) = sum(isi<=0.003)- sum(isi<=0.002);
        
        isi60Hz(uu) = max(isiperbin(15:18));%sum(isi<=0.018) - sum(isi<=0.015);
        isi30Hz(uu) = max(isiperbin(32:36));%sum(isi<=0.036) - sum(isi<=0.032);
        
        
        isi3_10(uu) = nanmean(sum_isi(isi,0.003,0.010,0.001));
        isi5_15(uu) = nanmean(sum_isi(isi,0.005,0.015,0.001));
        isi20_30(uu) = nanmean(sum_isi(isi,0.020,0.030,0.001));
        isi38_45(uu) = nanmean(sum_isi(isi,0.038,0.045,0.001));
        
        diff_isi = [0; diff(isiperbin)];
        [~,isi_peak] = max(isiperbin(4:end));
        isi_peak = isi_peak + 4 - 1;
        increases_topeak = numel(find(diff_isi(4:isi_peak) > 0));
        decreases_afterpeak = numel(find(diff_isi(isi_peak+1:isi_peak+15) < 0));
        if increases_topeak > 0.5*(isi_peak-4+1) && decreases_afterpeak > .5*15 &&...
                any(isiperbin(isi_peak) > [isi1ms(uu) isi2ms(uu) isi3ms(uu)])
            well_defined_isipeak = 1;
        else
            well_defined_isipeak = 0;
        end
        
        isiperbin_nolow = isiperbin;
        isiperbin_nolow(1:2) = [];
        
        
        
        if  all(isi1ms(uu) > 1.5.*[isi3ms(uu) isi3_10(uu) max(isiperbin_nolow)]) &&...
                abs(isi1ms(uu)-isi2ms(uu))/isi2ms(uu) < 0.10 && well_defined_isipeak == 0
            refractory(uu) = 1;
        else
            refractory(uu) = 0;
        end
        
        clustvar(uu) = block.clust(uu).clustvar;
        
        grad_pts = ceil(nsamp*.15);
        grad_15(uu) = nanmean(abs(diff(normwfs(end-grad_pts:end,uu)))); % maybe threshold at <1e-3
        
        
        
        if showoff
            figure(100+uu);
            plot(FSwin2(1:end),'b--','linewidth',1.5); hold on
            plot(RSUwin2(1:end),'r--','linewidth',1.5);
            plot(normwfs(:,uu),'k'); hold on;
            text(.55*nsamp,-0.55,['FS k val is ' num2str(kvfs(uu)) ', p-val is ' num2str(p_fs(uu))]);
            text(.55*nsamp,-0.6,['RSU k val is ' num2str(kvrsu(uu)) ', p-val is ' num2str(p_rsu(uu))]);
            text(.55*nsamp,-0.35,['RSU xcorr is ' num2str(xcrsu(uu)) ', chisq is ' num2str(chirsu(uu))]);
            text(.55*nsamp,-0.30,['FS xcorr is ' num2str(xcfs(uu)) ', chisq is ' num2str(chifs(uu))]);
            text(.55*nsamp,-0.65,['rsu dist is ' num2str(rsutempdist(uu))]);
            text(.55*nsamp,-0.7,['fs dist is ' num2str(fstempdist(uu))]);
            text(.55*nsamp,-0.90,['Neg curvature is ' num2str(crv(uu))]);
            text(.55*nsamp,-0.8,['ISI error (<3, <2): ' num2str(toofast3(uu)*100) '%; '...
                num2str(toofast2(uu)*100) '%']);
            text(.55*nsamp,-0.85,['ISI error (<1, 2/1ratio): ' num2str(toofast1(uu)*100) '%; '...
                num2str(toofast2(uu)/toofast1(uu))]);
            text(.55*nsamp,-0.75,['Clustvar is ' num2str(CONTCHUNKS_24h(ch).CELL(uu).clustvar) '.']);
            text(.55*nsamp,-0.45,['Spike amplitude is ' num2str(wfmins(uu))]);
            text(.55*nsamp,-0.5,['Rising phase best fit ' num2str(rise(uu)) '.']);
            text(.55*nsamp,-0.95,['Refractory is ' num2str(refractory(uu)) '.']);
            text(.55*nsamp,-0.40,['Mean gradient at end is ' num2str(grad_15(uu)) '.']);
            ylims = get(gca,'ylim');
            set(gca,'ylim',[-1 ylims(2)],'xlim',[1 nsamp+1]);
        end
        
        if kvfs(uu) < kvrsu(uu); % dealing with an FS cell?
            
            %set isi quality to 1-4 based on contamination
            if toofast3(uu)<0.03;
                isiflag(uu) = 1;
            elseif toofast3(uu)>0.03 && toofast3(uu)<=0.05;
                isiflag(uu) = 2;
            elseif toofast3(uu)>0.05 && toofast3(uu)<=0.12;
                isiflag(uu) = 3;
            elseif toofast3(uu)>0.12
                isiflag(uu) = 4;
            end
            
        elseif kvrsu(uu) <= kvfs(uu); % dealing with RSU?
            
            if toofast3(uu)<0.015;
                isiflag(uu) = 1;
            elseif toofast3(uu)>0.015 && toofast3(uu)<=0.03;
                isiflag(uu) = 2;
            elseif toofast3(uu)>0.03 && toofast3(uu)<=0.10;
                isiflag(uu) = 3;
            elseif toofast3(uu)>0.10
                isiflag(uu) = 4;
            end
            
        end
        
        % Check the rough firing rate of the cluster:
        elpsdtime       = block.clust(uu).time(end) - block.clust(uu).time(1);
        firingrate(uu)  = numel(block.clust(uu).idx)/elpsdtime;
        
        if showoff; text(32,-0.25,['Firing rate is ' num2str(firingrate(uu))]); end
        
        %%% this is where it changes BIGLY
        
        %---------- SET QUALITY 1-4 using nested if loops: -------------------------
        %first determine if cell is FS or RSU
        celltype{uu} = 'none';
        shape_measures = sum([kvfs(uu) < kvrsu(uu), xcfs(uu) > xcrsu(uu), chifs(uu)<chirsu(uu)]);
        if shape_measures >= 2 && fstempdist(uu) < rsutempdist(uu)
            celltype{uu} = 'fs'; % FS cell
        elseif shape_measures <= 1 && rsutempdist(uu) < fstempdist(uu)
            celltype{uu} = 'rsu';
        elseif shape_measures <= 1 && fstempdist(uu) < rsutempdist(uu) ||...
                shape_measures >= 2 && rsutempdist(uu) < fstempdist(uu)
            celltype{uu} = 'unclear';
        end
        
        
        
        % IF LOOKS LIKE FS CELL -------------------------------------------
        if strcmp(celltype{uu},'fs')
            celltype_assigned = 1;
            if verbose, disp('pFS cell'); end
            
            %now qualify 1-4 based on kvfs,sigdiff,clustvar and rising phase
            if      kvfs(uu)<=0.2 && sigdiff(uu)<=10 && clustvar(uu)<2000 && toofast3(2)<=.02 && firingrate(uu) >= 0.1
                autodecision(uu) = 1;
                
            elseif  any([kvfs(uu)<=0.20 chifs(uu) < 8]) && clustvar(uu)<10000 && sigdiff(uu)<=15 ||...
                    any([kvfs(uu)<=0.20 chifs(uu) < 8]) && sigdiff(uu)<=5 &&  rise(uu)<-0.005 ||...
                    any([kvfs(uu)<=0.25 chifs(uu) < 8]) && sigdiff(uu)<=10 && isiflag(uu)<=3 &&  rise(uu)<-0.005 ||...
                    any([kvfs(uu)<=0.20 chifs(uu) < 8]) && isiflag(uu)<=3 && rise(uu)<-0.005 ||...
                    clustvar(uu)<10000 && sigdiff(uu)<=10 || ...
                    isiflag(uu)==1 && sigdiff(uu)<=12 && clustvar(uu)<=5000 && rise(uu)<-0.006 && firingrate(uu) >= 6 || ...
                    isiflag(uu) == 1 && sigdiff(uu) <= 14 && clustvar(uu) <= 10000 && rise(uu)<-0.006 && kvfs(uu) <= 0.25
                autodecision(uu) = 2;
                
            elseif  kvfs(uu)<=0.20 && sigdiff(uu)<=10 && isiflag(uu)<=3 ||...
                    kvfs(uu)<=0.20 && sigdiff(uu)<=10 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                    kvfs(uu)<=0.21 && sigdiff(uu)<=10 && isiflag(uu)<=3 || ...
                    kvfs(uu)<=0.30 && sigdiff(uu)<=15 && isiflag(uu)<=3 ||...
                    kvfs(uu)<=0.30 && sigdiff(uu)<=20 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                    kvfs(uu)<=0.30 && sigdiff(uu)<=15 ||...
                    kvfs(uu)<=0.30 && isiflag(uu)<=3 || ...
                    sigdiff(uu)<20 && clustvar(uu)<15000||...
                    sigdiff(uu)<15 && rise(uu)<-0.004 && isiflag(uu)<=3 || ...
                    sigdiff(uu)<20 && isiflag(uu)<=3 && clustvar(uu) < 20000
                autodecision(uu) = 3;
                
            else
                autodecision(uu) = 4;
            end
            
            %neg curvature?
            if crv(uu) == 1
                if kvfs(uu)<=0.20 && sigdiff(uu)<=10 && isiflag(uu)<=3 ||...
                        kvfs(uu)<=0.20 && sigdiff(uu)<=10 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                        kvfs(uu)<=0.21 && sigdiff(uu)<=10 && isiflag(uu)<=3 || ...
                        kvfs(uu)<=0.27 && sigdiff(uu)<=15 && isiflag(uu)<=3 ||...
                        kvfs(uu)<=0.27 && sigdiff(uu)<=20 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                        kvfs(uu)<=0.27 && sigdiff(uu)<=15 ||...
                        kvfs(uu)<=0.27 && isiflag(uu)<=3 || ...
                        kvfs(uu)<=0.27 && sigdiff(uu)<20 && clustvar(uu)<15000||...
                        kvfs(uu)<=0.27 && sigdiff(uu)<15 && rise(uu)<-0.004 && isiflag(uu)<=3 || ...
                        kvfs(uu)<=0.27 && sigdiff(uu)<20 && isiflag(uu)<=3 && clustvar(uu) < 20000
                    autodecision(uu) = 3;
                else
                    autodecision(uu) = 4;
                end
            end
            
            % refractory period?
            if refractory(uu) == 1
                if kvfs(uu)<=0.20 && sigdiff(uu)<=10 && isiflag(uu)<=3 ||...
                        kvfs(uu)<=0.20 && sigdiff(uu)<=10 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                        kvfs(uu)<=0.21 && sigdiff(uu)<=10 && isiflag(uu)<=3 || ...
                        kvfs(uu)<=0.30 && sigdiff(uu)<=15 && isiflag(uu)<=3 ||...
                        kvfs(uu)<=0.30 && sigdiff(uu)<=20 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                        kvfs(uu)<=0.30 && sigdiff(uu)<=15 ||...
                        kvfs(uu)<=0.30 && isiflag(uu)<=3 || ...
                        kvfs(uu)<=0.30 && sigdiff(uu)<20 && clustvar(uu)<15000||...
                        kvfs(uu)<=0.30 && sigdiff(uu)<15 && rise(uu)<-0.004 && isiflag(uu)<=3 || ...
                        kvfs(uu)<=0.30 && sigdiff(uu)<20 && isiflag(uu)<=3 && clustvar(uu) < 20000
                    autodecision(uu) = 3;
                else
                    autodecision(uu) = 4;
                end
            end
            
            %Catch contaimination - isi is above 3% and firing rate below 4hz
            if autodecision(uu) <= 2
                contamination_flag(uu) = 0;
                if firingrate(uu)<=10
                    if toofast2(uu)>0.03 && toofast2(uu)<=0.07
                        contamination_flag(uu) = 1;
                        autodecision(uu) = 3;
                    elseif  toofast2(uu)>=0.08
                        if kvfs(uu)<0.18
                            contamination_flag(uu) = 1;
                            autodecision(uu) = 3;
                        else
                            contamination_flag(uu) = 2;
                            autodecision(uu)=4;
                        end
                    else
                        % unless very high contamination in 3ms bin,
                        % it's a 2 (and keeps cont_flag=0)
                        if toofast3(uu) >= 0.07
                            contamination_flag(uu) = 1;
                            autodecision(uu) = 3;
                        end
                        contamination_flag(uu) = 0;
                    end
                elseif firingrate(uu) > 10 && firingrate(uu) <= 20
                    if toofast2(uu)>0.04 && toofast2(uu)<=0.08
                        contamination_flag(uu) = 1;
                        autodecision(uu) = 3;
                    elseif toofast2(uu)>=0.08
                        if kvfs(uu)<0.18 || sigdiff(uu) <= 3.5
                            contamination_flag(uu) = 1;
                            autodecision(uu) = 3;
                        else
                            contamination_flag(uu) = 2;
                            autodecision(uu)=4;
                        end
                    else
                        contamination_flag(uu) = 0;
                    end
                elseif firingrate(uu) > 20
                    if toofast2(uu)>0.04 && toofast2(uu)<=0.10
                        contamination_flag(uu) = 1;
                        autodecision(uu) = 3;
                    elseif  toofast2(uu)>=0.10
                        if all([kvrsu(uu) kvfs(uu)]<0.18) || sigdiff(uu) <= 3.5
                            contamination_flag(uu) = 1;
                            autodecision(uu) = 3;
                        else
                            contamination_flag(uu) = 2;
                            autodecision(uu)=4;
                        end
                    else
                        contamination_flag(uu) = 0;
                    end
                end
            else
                contamination_flag(uu) = -1;
            end
            
            % ELSEIF LOOKS LIKE RSU ---------------------------------------
        elseif strcmp(celltype{uu},'rsu') %RSU cell
            celltype_assigned = 1;
            if verbose, disp(['RSU cell']); end
            
            %now qualify 1-4 based on kvrsu,sigdiff,clustvar and rising phase
            
                if      kvrsu(uu)<=0.12 && sigdiff(uu)<=5 && clustvar(uu)<2000 && isiflag(uu)<=2 && rise(uu) <-0.006
                    autodecision(uu) = 1;
                    
                elseif  kvrsu(uu)<=0.20 && clustvar(uu)<10000 && sigdiff(uu)<=15 && rise(uu)<-0.005||...
                        kvrsu(uu)<=0.20 && sigdiff(uu)<=5 && isiflag(uu)<=3 &&  rise(uu)<-0.005 ||...
                        kvrsu(uu)<=0.25 && sigdiff(uu)<=10 && isiflag(uu)<=3 &&  rise(uu)<-0.005 ||...
                        kvrsu(uu)<=0.20 && isiflag(uu)<=3 && rise(uu)<-0.005 ||...
                        kvrsu(uu)<=0.20 && clustvar(uu)<10000 && sigdiff(uu)<=10 && rise(uu)<-0.005 ||... % cut here
                        clustvar(uu)<10000 && sigdiff(uu)<=8 && rise(uu)<-0.005 || ...
                        isiflag(uu)<=2 && sigdiff(uu)<=10 && clustvar(uu)<=8000 && rise(uu)<-0.006 ||...
                        isiflag(uu)<=1 && sigdiff(uu)<=12 && clustvar(uu)<=6000 && rise(uu)<=-0.006 && all([toofast2(uu) toofast3(uu)]<0.01) ||...
                        isiflag(uu) <=1 && sigdiff(uu)<=8 && rise(uu)<=-0.006 && well_defined_isipeak == 1
                    
                    autodecision(uu) = 2;
                    
                elseif  kvrsu(uu)<=0.20 && sigdiff(uu)<=10 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                        kvrsu(uu)<=0.21 && sigdiff(uu)<=10 && isiflag(uu)>=3 || ...
                        kvrsu(uu)<=0.30 && sigdiff(uu)<=15 && isiflag(uu)<=3 ||...
                        kvrsu(uu)<=0.30 && sigdiff(uu)<=20 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                        kvrsu(uu)<=0.30 && sigdiff(uu)<=15 ||...
                        kvrsu(uu)<=0.30 && isiflag(uu)<=3 || ...
                        sigdiff(uu)<20 && clustvar(uu)<15000||...
                        sigdiff(uu)<15 && rise(uu)<-0.004 && isiflag(uu)<=3||...
                        sigdiff(uu)<20 && isiflag(uu)<=3 && clustvar(uu) < 20000
                    autodecision(uu) = 3;
                    
                else
                    autodecision(uu) = 4;
                end
            
            
            %neg curvature?
            if crv(uu) == 1
                if kvrsu(uu)<=0.20 && sigdiff(uu)<=10 && isiflag(uu)<=3 ||...
                        kvrsu(uu)<=0.20 && sigdiff(uu)<=10 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                        kvrsu(uu)<=0.21 && sigdiff(uu)<=10 && isiflag(uu)<=3 || ...
                        kvrsu(uu)<=0.27 && sigdiff(uu)<=15 && isiflag(uu)<=3 ||...
                        kvrsu(uu)<=0.27 && sigdiff(uu)<=20 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                        kvrsu(uu)<=0.27 && sigdiff(uu)<=15 ||...
                        kvrsu(uu)<=0.27 && isiflag(uu)<=3 || ...
                        kvrsu(uu)<=0.27 && sigdiff(uu)<20 && clustvar(uu)<15000||...
                        kvrsu(uu)<=0.27 && sigdiff(uu)<15 && rise(uu)<-0.004 && isiflag(uu)<=3 || ...
                        kvrsu(uu)<=0.27 && sigdiff(uu)<20 && isiflag(uu)<=3 && clustvar(uu) < 20000
                    autodecision(uu) = 3;
                else
                    autodecision(uu) = 4;
                end
            end
            
            % refractory period?
            if refractory(uu) == 1
                if kvrsu(uu)<=0.20 && sigdiff(uu)<=10 && isiflag(uu)<=3 ||...
                        kvrsu(uu)<=0.20 && sigdiff(uu)<=10 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                        kvrsu(uu)<=0.21 && sigdiff(uu)<=10 && isiflag(uu)<=3 || ...
                        kvrsu(uu)<=0.30 && sigdiff(uu)<=15 && isiflag(uu)<=3 ||...
                        kvrsu(uu)<=0.30 && sigdiff(uu)<=20 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                        kvrsu(uu)<=0.30 && sigdiff(uu)<=15 ||...
                        kvrsu(uu)<=0.30 && isiflag(uu)<=3 || ...
                        kvrsu(uu)<=0.30 && sigdiff(uu)<20 && clustvar(uu)<15000||...
                        kvrsu(uu)<=0.30 && sigdiff(uu)<15 && rise(uu)<-0.004 && isiflag(uu)<=3 || ...
                        kvrsu(uu)<=0.30 && sigdiff(uu)<20 && isiflag(uu)<=3 && clustvar(uu) < 20000
                    autodecision(uu) = 3;
                else
                    autodecision(uu) = 4;
                end
            end
            
            %Catch contaimination - isi is above 3% and firing rate below 4hz
            if autodecision(uu) <= 2
                contamination_flag(uu) = 0;
                if firingrate(uu)<=7
                    if toofast2(uu)>0.02 && toofast2(uu)<=0.05
                        if toofast2(uu) < 0.025 && wfmins(uu) < -40.0 && ...
                                sigdiff(uu) <= 5 && kvrsu(uu) < 0.18
                            contamination_flag(uu) = 0;
                            autodecision(uu) = 2;
                        else
                            contamination_flag(uu) = 1;
                            autodecision(uu) = 3;
                        end
                    elseif  toofast2(uu)>=0.05
                        if kvrsu(uu)<0.18 || sigdiff(uu) <= 3.5
                            contamination_flag(uu) = 1;
                            autodecision(uu) = 3;
                        else
                            contamination_flag(uu) = 2;
                            autodecision(uu)=4;
                        end
                    else
                        contamination_flag(uu) = 0;
                    end
                elseif firingrate(uu) > 7 && firingrate(uu) <= 10
                    if toofast2(uu)>0.025 && toofast2(uu)<=0.05
                        if sigdiff(uu) <= 3.5 && kvrsu(uu) <= 0.12 && ...
                                toofast2(uu)/toofast1(uu)>=7.0
                            contamination_flag(uu) = 0;
                        else
                            if toofast2(uu) < 0.028 && wfmins(uu) < -40.0 && ...
                                    sigdiff(uu) <= 5 && kvrsu(uu) < 0.18
                                contamination_flag(uu) = 0;
                                autodecision(uu) = 2;
                            else
                                contamination_flag(uu) = 1;
                                autodecision(uu) = 3;
                            end
                        end
                    elseif  toofast2(uu)>=0.05
                        if kvrsu(uu)<0.18 || sigdiff(uu) <= 3.5
                            contamination_flag(uu) = 1;
                            autodecision(uu) = 3;
                        else
                            contamination_flag(uu) = 2;
                            autodecision(uu)=4;
                        end
                    else
                        contamination_flag(uu) = 0;
                    end
                elseif firingrate(uu) > 10
                    if toofast2(uu)>0.03 && toofast2(uu)<=0.08
                        if sigdiff(uu) <= 3.5 && kvrsu(uu) <= 0.12 && ...
                                toofast2(uu)/toofast1(uu)>=7.0 && toofast2(uu) <= 0.05
                            contamination_flag(uu) = 0;
                        else
                            contamination_flag(uu) = 1;
                            autodecision(uu) = 3;
                        end
                    elseif  toofast2(uu)>=0.08
                        if kvrsu(uu)<0.18 || sigdiff(uu) <= 3.5
                            contamination_flag(uu) = 1;
                            autodecision(uu) = 3;
                        else
                            contamination_flag(uu) = 2;
                            autodecision(uu)=4;
                        end
                    else
                        if toofast2(uu)/toofast1(uu) <= 6 && sigdiff(uu) >= 4
                            if toofast2(uu) >= 0.022 && toofast3(uu) >= 0.035
                                contamination_flag(uu) = 1;
                                autodecision(uu) = 3;
                            end
                        else
                            contamination_flag(uu) = 0;
                        end
                        
                        
                    end
                end
            else
                contamination_flag(uu) = -1;
            end
            
            % ELSEIF UNCLEAR CELL TYPE --------------------------------
        elseif strcmp(celltype{uu},'unclear') %unclear cell
            celltype_assigned = 1;
            if verbose,disp(['unclassified cell type']); end
            
            
            if sigdiff(uu)<10 && any([kvrsu(uu) kvfs(uu)]<=0.20) && clustvar(uu)<17000 ||...
                    sigdiff(uu)<8 && any([kvrsu(uu) kvfs(uu)]<=0.20) || ...
                    sigdiff(uu) < 12 && any([kvrsu(uu) kvfs(uu)]<=0.15)
                autodecision(uu) = 2;
                
            elseif any([kvrsu(uu) kvfs(uu)]<=0.20) && sigdiff(uu)<=10 && isiflag(uu)<=3 ||...
                    any([kvrsu(uu) kvfs(uu)]<=0.20) && sigdiff(uu)<=10 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                    any([kvrsu(uu) kvfs(uu)]<=0.21) && sigdiff(uu)<=10 && isiflag(uu)<=3 || ...
                    any([kvrsu(uu) kvfs(uu)]<=0.30) && sigdiff(uu)<=15 && isiflag(uu)<=3 ||...
                    any([kvrsu(uu) kvfs(uu)]<=0.30) && sigdiff(uu)<=20 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                    any([kvrsu(uu) kvfs(uu)]<=0.30) && sigdiff(uu)<=15 ||...
                    any([kvrsu(uu) kvfs(uu)]<=0.30) && isiflag(uu)<=3 || ...
                    sigdiff(uu)<20 && clustvar(uu)<15000||...
                    sigdiff(uu)<15 && rise(uu)<-0.004 && isiflag(uu)<=3 || ...
                    sigdiff(uu)<20 && isiflag(uu)<=3 && clustvar(uu) < 20000
                autodecision(uu) = 3;
            else
                autodecision(uu) = 4;
            end
            
            %neg curvature?
            if crv(uu) == 1
                
                if any([kvrsu(uu) kvfs(uu)] <= 0.16) && sigdiff(uu) <= 5 && isiflag(uu) <=1 && ...
                        clustvar(uu) <= 10000 && rise(uu) <= -0.006 && autodecision(uu) == 2
                    autodecision(uu) = 2;
                elseif any([kvrsu(uu) kvfs(uu)]<=0.20) && sigdiff(uu)<=10 && isiflag(uu)<=3 ||...
                        any([kvrsu(uu) kvfs(uu)]<=0.20) && sigdiff(uu)<=10 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                        any([kvrsu(uu) kvfs(uu)]<=0.21) && sigdiff(uu)<=10 && isiflag(uu)<=3 || ...
                        any([kvrsu(uu) kvfs(uu)]<=0.27) && sigdiff(uu)<=15 && isiflag(uu)<=3 ||...
                        any([kvrsu(uu) kvfs(uu)]<=0.27) && sigdiff(uu)<=20 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                        any([kvrsu(uu) kvfs(uu)]<=0.27) && sigdiff(uu)<=15 ||...
                        any([kvrsu(uu) kvfs(uu)]<=0.27) && isiflag(uu)<=3 || ...
                        any([kvrsu(uu) kvfs(uu)]<=0.27) && sigdiff(uu)<20 && clustvar(uu)<15000||...
                        any([kvrsu(uu) kvfs(uu)]<=0.27) && sigdiff(uu)<15 && rise(uu)<-0.004 && isiflag(uu)<=3 || ...
                        any([kvrsu(uu) kvfs(uu)]<=0.27) && sigdiff(uu)<20 && isiflag(uu)<=3 && clustvar(uu) < 20000
                    autodecision(uu) = 3;
                else
                    autodecision(uu) = 4;
                end
                if firingrate(uu) < 0.05 && toofast1(uu) >= 0.012 || ...
                        firingrate(uu) < 0.05 && toofast1(uu) >= 0.12 && toofast2(uu)/toofast1(uu) <= 4 || ...
                        firingrate(uu) < 0.1 && toofast1(uu) >= 0.015 || ...
                        isiflag(uu) == 4
                    autodecision(uu) = 4;
                end
            end
            
            % refractory period?
            if refractory(uu) == 1
                if any([kvrsu(uu) kvfs(uu)]<=0.20) && sigdiff(uu)<=10 && isiflag(uu)<=3 ||...
                        any([kvrsu(uu) kvfs(uu)]<=0.20) && sigdiff(uu)<=10 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                        any([kvrsu(uu) kvfs(uu)]<=0.21) && sigdiff(uu)<=10 && isiflag(uu)<=3 || ...
                        any([kvrsu(uu) kvfs(uu)]<=0.30) && sigdiff(uu)<=15 && isiflag(uu)<=3 ||...
                        any([kvrsu(uu) kvfs(uu)]<=0.30) && sigdiff(uu)<=20 && isiflag(uu)<=3 && clustvar(uu)<20000 ||...
                        any([kvrsu(uu) kvfs(uu)]<=0.30) && sigdiff(uu)<=15 ||...
                        any([kvrsu(uu) kvfs(uu)]<=0.30) && isiflag(uu)<=3 || ...
                        any([kvrsu(uu) kvfs(uu)]<=0.30) && sigdiff(uu)<20 && clustvar(uu)<15000||...
                        any([kvrsu(uu) kvfs(uu)]<=0.30) && sigdiff(uu)<15 && rise(uu)<-0.004 && isiflag(uu)<=3 || ...
                        any([kvrsu(uu) kvfs(uu)]<=0.30) && sigdiff(uu)<20 && isiflag(uu)<=3 && clustvar(uu) < 20000
                    autodecision(uu) = 3;
                else
                    autodecision(uu) = 4;
                end
            end
            
            %Catch contaimination - isi is above 3% and firing rate below 4hz
            if autodecision(uu) == 2
                contamination_flag(uu) = 0;
                if firingrate(uu)<=10
                    if toofast2(uu)>0.03 && toofast2(uu)<=0.05
                        contamination_flag(uu) = 1;
                        autodecision(uu) = 3;
                    elseif  toofast2(uu)>=0.05
                        if any([kvrsu(uu) kvfs(uu)]<0.18)
                            contamination_flag(uu) = 1;
                            autodecision(uu) = 3;
                        else
                            contamination_flag(uu) = 2;
                            autodecision(uu)=4;
                        end
                    else
                        contamination_flag(uu) = 0;
                    end
                elseif firingrate(uu) > 10 && firingrate(uu) <= 20
                    if toofast2(uu)>0.04 && toofast2(uu)<=0.08
                        contamination_flag(uu) = 1;
                        autodecision(uu) = 3;
                    elseif  toofast2(uu)>=0.08
                        if any([kvrsu(uu) kvfs(uu)]<0.18) || sigdiff(uu) <= 3.5
                            contamination_flag(uu) = 1;
                            autodecision(uu) = 3;
                        else
                            contamination_flag(uu) = 2;
                            autodecision(uu)=4;
                        end
                    else
                        contamination_flag(uu) = 0;
                    end
                elseif firingrate(uu) > 20
                    if toofast2(uu)>0.04 && toofast2(uu)<=0.10
                        contamination_flag(uu) = 1;
                        autodecision(uu) = 3;
                    elseif  toofast2(uu)>=0.10
                        if all([kvrsu(uu) kvfs(uu)]<0.18) || sigdiff(uu) <= 3.5
                            contamination_flag(uu) = 1;
                            autodecision(uu) = 3;
                        else
                            contamination_flag(uu) = 2;
                            autodecision(uu)=4;
                        end
                    else
                        contamination_flag(uu) = 0;
                    end
                end
            else
                contamination_flag(uu) = -1;
            end
           
        elseif strcmp(celltype{uu},'none')
            disp('WARNING! CELL TYPE NOT ASSIGNED!');
            celltype_assigned = 0;
        end
        
        %%%%%%%%%%%%%%%%% ADDITIONAL QUALITY FILTERS %%%%%%%%%%%%%%%%%%%%%%
        if celltype_assigned
        %catch weird noise that should be qual 4
        if refractory(uu) == 1 && crv(uu) == 1  || ...  % both no refract and neg curb
                clustvar(uu) > 50000 && crv(uu) == 1 || ...  % crazy high clustvar
                min(kvfs(uu),kvrsu(uu)) > 0.40       ||... % k test is very very high
                toofast1(uu) > 0.1  || ... % really bad isi contamination in 1 ms bin
                toofast2(uu) > 0.12
            %abs(grad_15(uu)) < 0.001                    % too flat on end
            autodecision(uu)=4;
        end
        
        % some more corrections from 3 to 4
        if refractory(uu) == 1 || crv(uu) == 1
            if abs(grad_15(uu)) < 0.004 && autodecision(uu) == 3 || ...
                    isiflag(uu) == 4 || ...
                    isnan(toofast2(uu)/toofast1(uu)) || ...
                    isi60Hz(uu) / isi2ms(uu) > 50 && isi60Hz(uu) / isi3ms(uu) > 50
                
                autodecision(uu) = 4;
            end
            
            
        end
        
        % isi shape correction
        if autodecision(uu) == 3 && refractory(uu) == 1
            isi_lowlim = 0.005;
            isi_hilim = 0.015;
            isi_binsz = 0.001;
            toofast_later = sum_isi(isi,isi_lowlim,isi_hilim,isi_binsz);
            
            if toofast3(uu) < 0.008 || ...
                    toofast3(uu)<0.035 && isiflag(uu)<=2 && sigdiff(uu) <= 5
                % if contamination is very low or cell looks really
                % good, check if there is a peak after 3ms
                
                if nanmean(toofast_later) > 1.25*isi1ms(uu) || ...
                        nanmean(toofast_later) > 1.2*nanmean([isi1ms(uu) isi2ms(uu) isi3ms(uu)])
                    autodecision(uu) = 2;
                else
                    autodecision(uu) = 3;
                end
            end
            
            if isi1ms(uu) > 3*(isi2ms(uu)+isi3ms(uu)) && ...
                    isi1ms(uu) > 3*nanmean(toofast_later)
                autodecision(uu) = 4;
            end
        end
        
        % more isi shape correction
        % if refractory looks good, make sure it's really good
        if refractory(uu) == 0
            % this for 2's that should be 3's
            isiperbin_no2 = isiperbin;
            isiperbin_no2(2) = [];
            isiperbin_no1 = isiperbin;
            isiperbin_no1(1) = [];
            if all(isi2ms(uu) > 1.5.*[isi3ms(uu) isi1ms(uu) isi3_10(uu) max(isiperbin_no2)]) || ...
                    all(isi1ms(uu) > 1.5.*[isi3ms(uu) isi3_10(uu) max(isiperbin_nolow)])
                if autodecision(uu) == 2
                    if (toofast2(uu) > 0.01 && toofast3(uu) > 0.01 && well_defined_isipeak == 0) ||...
                            sigdiff(uu) > 5 && all([chirsu(uu) chifs(uu)]>10) && all([kvrsu(uu) kvfs(uu)]>0.10) ||...
                            sigdiff(uu) > 5 && clustvar(uu) > 20000 && firingrate(uu) < 0.5 ||...
                            wfmins(uu) > -35
                        % apply this unless contamination in 2ms and 3ms bins
                        % is very low OR shape is not very good OR amplitude
                        % is low
                        autodecision(uu) = 3;
                    end
                elseif autodecision(uu) == 3
                    if any([kvrsu(uu) kvfs(uu)]<0.18) && sigdiff(uu) <= 8
                        autodecision(uu) = 3;
                    else
                        autodecision(uu) = 4;
                    end
                end
                
            elseif all(isi60Hz(uu) > isiperbin(5:15)) &&...
                    all(isi60Hz(uu) > isiperbin(20:30)) &&...
                    all(isi30Hz(uu) > isiperbin(20:30)) &&...
                    all(isi30Hz(uu) > isiperbin(38:45))
                
                if autodecision(uu) == 2
                    if any([kvrsu(uu) kvfs(uu)] < 0.10) && sigdiff(uu) <= 6 && ...
                            (isi60Hz(uu) + isi30Hz(uu))/numel(isi) < 0.01
                        autodecision(uu) = 2;
                    else
                        autodecision(uu) = 3;
                    end
                elseif autodecision(uu) == 3
                    if any([kvrsu(uu) kvfs(uu)]<0.18) && sigdiff(uu) <= 10
                        autodecision(uu) = 3;
                    else
                        autodecision(uu) = 4;
                    end
                end
                
            else
                %this for 3's that should be 2's
                if autodecision(uu) == 3
                    % if shape looks really good and contamination is not
                    % too bad
                    if toofast2(uu) <= 0.03 && sigdiff(uu) <= 6 &&...
                            any([kvrsu(uu) kvfs(uu)]<=0.15) && toofast1(uu) <= 0.005
                        % check if there is a peak in the ISI distribution
                        % under 100 ms
                        biggerthan3 = find(isiperbin(4:100) > isiperbin(3));
                        if sum(diff(biggerthan3) == 1) > 4 % if more than 4 consectuive bins larger than 3ms bin
                            autodecision(uu) = 2;
                        end
                    end
                end
            end
        end
        
        % waveform shape correction
        if autodecision(uu) == 2
            if all([fstempdist(uu) rsutempdist(uu)] > 12) && ...
                    all([kvrsu(uu) kvfs(uu)] > 0.18) && toofast3(uu) > 0.05
                autodecision(uu) = 3;
            elseif all([fstempdist(uu) rsutempdist(uu)] > 10) && ...
                    all([chirsu(uu) chifs(uu)]>30)
                autodecision(uu) = 3;
            end
        end
        
        
        % low FR correction
        if firingrate(uu) < 0.01 && autodecision(uu) == 2
            autodecision(uu) = 3;
        elseif firingrate(uu) < 0.1 && autodecision(uu) == 2 && ...
                clustvar(uu) > 10000 && sigdiff(uu) > 8 && all([kvrsu(uu) kvfs(uu)]>0.18)
            autodecision(uu) = 3;
        end
        
        % shape corrections to convert bad 2s into 3s
        if autodecision(uu) <= 2 && abs(grad_15(uu)) <= 0.001
            if sigdiff(uu) <= 3.5 && isiflag(uu) <=2 && toofast2(uu) <= 0.02
                autodecision(uu) = 2;
            else
                autodecision(uu) = 3;
            end
        end
        
        % correct for crv=1 but beautiful otherwise
        if crv(uu) == 1 && autodecision(uu) == 3
            if sigdiff(uu) <= 3 && any([chirsu(uu) chifs(uu)]<5) && ...
                    any([xcrsu(uu) xcfs(uu)]>0.99) && ...
                    isiflag(uu) <= 3 && refractory(uu) == 0
                autodecision(uu) = 2;
            end
        end
        
        % if a few things are bad, make it a 3
        if autodecision == 2
            % voltage window for unclear cell - if it's right in between RSU
            % and FS templates becomes a 3
            % other criteria? Low FR? High clustvar?
            
            if nsamp == 59
                Vwin_3 = [38 0.18; 42 0.27];
            elseif nsamp == 91
                Vwin_3 = [62 0.18; 68 0.27];
            end
            plot_patch = 0;
            if plot_patch
                Vwin_patch = patch('Xdata',[Vwin_3(1,1) Vwin_3(2,1) Vwin_3(2,1) Vwin_3(1,1)],...
                    'Ydata',[Vwin_3(1,2) Vwin_3(1,2) Vwin_3(2,2) Vwin_3(2,2)],...
                    'facecolor','r','facealpha',0.3);
            end
            wf_section = normwfs(Vwin_3(1,1):Vwin_3(2,1),uu);
            if any(wf_section > Vwin_3(1,2) & wf_section < Vwin_3(2,2))
                window_catch = 1; % WF crosses this window. May be a multi-unit (quality 3)
            else
                window_catch = 0;
            end
            
            if ( clustvar(uu) > 20000 && isi1ms(uu) > 0.01 && isiflag(uu) > 1 &&...
                    all([kvrsu(uu) kvfs(uu)]>0.20) && all([chifs(uu) chirsu(uu)] > 10) ) ||...
                    ( all([kvrsu(uu) kvfs(uu)]>0.15) && all([chifs(uu) chirsu(uu)] > 18) && ...
                    all([xcrsu(uu) xcfs(uu)]<.97) && sigdiff(uu) > 6 ) ||...
                    window_catch == 1 && shape_measures > 0 && shape_measures < 3 && ...
                    all([chifs(uu) chirsu(uu)] > 10) && all([kvrsu(uu) kvfs(uu)]>0.10) && ...
                    sigdiff(uu) > 5 && all([xcrsu(uu) xcfs(uu)]<.97)
                
                
                
                autodecision(uu) = 3;
            end
        end
        
        
        
        
        %now correct for amplitude
        if wfmins(uu) >= -32
            if autodecision(uu) == 2 && ...
                    (all([kvrsu(uu) kvfs(uu)] > 0.2) || (toofast3(uu) > 0.05 && sigdiff(uu) >= 10) ||...
                    firingrate(uu) < 0.05)
                autodecision(uu) = 3;
            end
        end
        if wfmins(uu) >= -28 %|| wfmins < -200 && autodecision(uu) <=2
            if crv(uu) == 1 && autodecision(uu) == 3
                autodecision(uu) = 4;
            elseif autodecision(uu) <= 2
                autodecision(uu) = 3;
            end
        end
        if wfmins(uu) >= -24 && autodecision(uu) <=3
            autodecision(uu) = 4;
        end
        
	else
		autodecision(uu) = 4;
	end

        if showoff
            text(.55*nsamp,-0.2,['contamination flag is ' num2str(contamination_flag(uu)) '.']);
            text(.55*nsamp,-0.1,['Autosort qual is ' num2str(autodecision(uu)) '.']);
            uiwait(gcf);
        end
        
        if verbose, fprintf('Autoquality %u.\n\n',autodecision(uu)); end
        
        %%%%%%%% END AUTOQUALITY ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % assign Quality
        block.clust(uu).quality = autodecision(uu);
        block.clust(uu).newQuality = autodecision(uu);
        
        if showoff, uiwait(gcf); end
    end
 
end



datOut = block;

% now save the output:
% [vFile, vDir] = uiputfile('.mat','Save your data!');
% save([vDir vFile],'channel','-v7.3');




