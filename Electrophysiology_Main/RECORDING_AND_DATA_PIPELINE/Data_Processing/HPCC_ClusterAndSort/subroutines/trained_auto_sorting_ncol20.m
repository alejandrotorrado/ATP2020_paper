function [newblock,p,trimvec] = trained_auto_sorting_ncol20(Mdl,block,nsamp,interp_samples,swprng,min_rng, plotty)

% TRAINED_AUTO_SORTING is an algorithmic approach to "spike sorting", i.e.
% the process of determing whether a cluster of waveforms (produced by any
% clustering algorithm) represents a single unit, a multi unit, or a noise
% source. This script will extract 14 features (3/29/17: updated to 20)
% from the cluster and its constitutent waveforms and use these to make a
% decision. The decision making process calls upon a supervised learning
% model (machine learning) that, when built, aimed to take hundreds of
% manually coded clusters and learn the strongest set of discriminatory
% rules based on that input. Here, the same rules are applied to novel
% inputs (the clusters). This training dataset can be updated with new
% information to increase accuracy if anything is systematically
% mislabeled.
%
% KBH 2/2017

ncol = 20; %changed to 20 on 3/29/17 ATP (added hourly_cont,negpostime,halfwidth,tailslope,decay_slope)

% Hard-coded confidence threhsolds for quality review flagging
qconf_thresh = 0.80;
sbconf_thresh = 0.15;


print_check = 1;
if print_check
    fprintf('For this animal, using:\nnsamp = %u, interp_samp = %u\nswprng = -%u:%u, min_rng = %u:%u\n',...
        nsamp,interp_samples,abs(swprng(1)),swprng(end),min_rng(1),min_rng(end));
end

%nsamp = 59;
%interp_samples = 97;
%swprng = -18:40;
%min_rng = 10:30;

% FS and RSU templates

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

% this is for 3x interpolation - interpolate the template so they can be
% compare to cluster WFs
x0      = 1:size(FSwindow,2);
x1      = 1:(36/(interp_samples-1)):size(FSwindow,2);
FSwin2  = interp1(x0,FSwindow,x1,'spline');
RSUwin2 = interp1(x0,RSUwindow,x1,'spline');

% find the sample at which the templates are at minimum (neg peak):
[~,rowmin0] = min([RSUwin2;FSwin2],[],2);
% check to make sure that the minimums are in the same sample pt:
if numel(unique(rowmin0)) == 1
    rowmin0 = unique(rowmin0);
else
    error('The interpolated waveform windows do no share the same minimum sample point! ');
end

% waveforms have been trimmed for PCA, trim down the templates to match
trim_start  = rowmin0 + swprng(1);
trim_end    = rowmin0 + swprng(end);
FSwin2      = FSwin2(trim_start:trim_end);
RSUwin2     = RSUwin2(trim_start:trim_end);

% then recalculate rowmin0
[~,rowmin0] = min([RSUwin2;FSwin2],[],2);
% check to make sure that the minimums are in the same sample:
if numel(unique(rowmin0)) == 1
    rowmin0 = unique(rowmin0);
else
    error('The interpolated waveform windows do no share the same minimum sample point! Stopped in beta_cluster_evaluation.m');
end


if iscell(block)
    ogdat = block{1};
else
    ogdat = block(1);
end

trimvec = zeros(size(ogdat.clust,2),1);

% go through all cells to find last spike time and infer experiment end
% hour
for gg = 1:size(ogdat.clust,2)
    last_t(gg) = ogdat.clust(gg).time(end);
end
max_hour = round(max(last_t./3600),0);

count = 0; p = []; scount = 0; score = [];
for ee = 1:size(ogdat.clust,2)
    count = count + 1;
    
    % get clustering output
    this = ogdat.clust(ee);
    
    dat = [];
    dat = this.time;
    X = diff(dat);
    h = histc(X,[0:0.001:10.0]); % previously calculated only out to 0.5 sec
    % This changed on 3/27/17 to 10 sec. This prevents over-representation
    % of noise in <3ms for slowly firing cells.
    
    if ~isempty(h)
        if sum(h)<1000
            % if cluster contains less than 1000 spikes, it is noise
            % NOTE: this may have to be changed if chunking in small chunks
            p(count,1:(ncol-1))   = zeros;
            p(count,ncol)     = 4; % ncol is the number of entries in p, so the last entry is always quality
            ogdat.clust(ee).quality = p(count,ncol);
            
            scount = scount + 1;
            score(scount,1) = 4;
            score(scount,2) = 1.0;
            score(scount,3:4) = 0;
            ogdat.clust(ee).qflag = 0;
            ogdat.clust(ee).score = score;
        else
            
            %disp(count);
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
            
            % - - - - - - - - - - - - - - - - - - - - -
            clear kvfs kvrsu firingrate sigdiff crv toofast isiflag autodecision
            clear rise clustvar
            
            % Metrics of how close the mean WFs are to the templates
            
            % KS test metric
            [~,~,kvfs]  = kstest2(FSwin2,normtrace);
            [~,~,kvrsu] = kstest2(RSUwin2',normtrace);
            
            % Absolute distance metric
            fstempdist  = sum( abs( FSwin2' - normtrace ));
            rsutempdist = sum( abs( RSUwin2' - normtrace  ));
            
            % Chi-square gof metric
            chirsu      = sum(((normtrace - RSUwin2').^2)./abs(RSUwin2'));
            chifs       = sum(((normtrace - FSwin2').^2)./abs(FSwin2'));
            
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
                disp('Failed at b0 calculation.');
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
            p(count,1:2) = gamfit(1:out_to,[],[],h(1:out_to)'); % this adds gamma parameters to p variable
            ygam2 = gampdf(1:out_to,p(count,1),p(count,2));
            
            % ISI contamination
            contamination = 100*(sum(h(1:2))/sum(h(1:end)));
            
            % Mean FR
            meanfr = numel(this.time)/(this.time(end)-this.time(1));
            
            % Use a different measure of contamination: get mean noise
            % contamination per hour in the middle 50% of recording (25% to 75%
            % of recording time) and then average the hourly means together
            frac = 0.10; % fraction of experiment to exclude on either side
            % e.g. 0.25 excludes first and last 25%
            edges = floor(.1*max_hour) : ceil((1-frac)*max_hour);
            for ii = 1:(size(edges,2)-1)
                this_hour_spikes = this.time(this.time >= edges(ii)*3600 & ...
                    this.time <= edges(ii+1)*3600);
                % only get contamination if cell is firing at more than 0.005
                % Hz in that hour (more than 18 spikes in that hour)
                if numel(this_hour_spikes) > 0.005*3600
                    this_hour_isi = diff(this_hour_spikes);
                    this_hour_h = histc(this_hour_isi,[0:0.001:10.0]);
                    hourly_cont(ii) = 100*sum(this_hour_h(1:2))/sum(this_hour_h(1:end));
                else
                    hourly_cont(ii) = NaN;
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
            
            peaktemplate = [0.3868 0.4064 0.4220 0.4371 0.4596 0.4859 0.5210...
                0.5651 0.6203 0.6875 0.7944 1.0000 0.7902 0.6830 0.6138 ...
                0.5590 0.5119 0.4745 0.4461 0.4201 0.4005]';
            
            ac_noise  = corr(peaktemplate, P1(20:40)/max(P1(20:40)) );
            % - - - - - - - - - - - - - - - - - - - - - - -
            
            % get neg_pos_time as one of the parameters
            % not this assumes a sampling_rate of 24414.06 Hz
            negpostime = calc_negpostime(normtrace); % returns negpostime in milliseconds
            
            % calc halfwidth in ms
            halfwidth = calc_halfwidth(normtrace);
            if isempty(halfwidth), halfwidth = 0; end
            
            %tailslope
            tailslope = calc_tailslope(normtrace);
            
            %% MAKE P VARIABLE
            p(count,3) = contamination;
            p(count,4) = contamination_byhour;
            p(count,5) = meanfr;
            p(count,6) = kvfs;
            p(count,7) = kvrsu;
            p(count,8) = fstempdist;
            p(count,9) = rsutempdist;
            p(count,10) = chirsu;
            p(count,11) = chifs;
            p(count,12) = rise;
            p(count,13) = decay_slope;
            p(count,14) = crv;
            p(count,15) = min(this.meantrace);
            p(count,16) = negpostime;
            p(count,17) = halfwidth;
            p(count,18) = tailslope;
            p(count,19) = ac_noise;
            
            % Y combines all the cells
            Y = p(count,[1:(ncol-1)]);
            
            % Run measurements through the model and predict the
            % quality:
            clear label
            scount = scount+1;
            [label,scoretmp,~] = predict(Mdl,Y); % use model to score
            % create score variable. This has 4 columns, organized as:
            % assigned_quality assigned_confidence secondbest_quality secondbest_confidence
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
            
            % QUALITY ASSIGNMENT:
            p(count,ncol) = label;
            ogdat.clust(ee).quality = label;
            
            % SCORE AND QUALITY REVIEW FLAG ASSIGNMENT:
            ogdat.clust(ee).score = score(scount,:);
            if score(scount,2) <= qconf_thresh && score(scount,4) >= sbconf_thresh
                ogdat.clust(ee).qflag = 1;
            else
                ogdat.clust(ee).qflag = 0;
            end
            
            % Figure out if this should be trimmed:
            if any(p(count,ncol) == [2,3]) && any(p(count,3:4)>=3) && any(p(count,3:4)<12) &&...
                    p(count,5)>2 && p(count,15)<0.5
                trimvec(ee,1) = 1;
            end
            
            if plotty
                g = figure(count);
                set(g,'units','normalized','position',[0.1 0.1 0.9 0.9]);
                subplot(1,3,1);
                plot(ygam2);
                yvals = get(gca,'ylim');
                txt = text(35,yvals(2)*0.75,['Predicted q: ' num2str(label) ]);
                set(txt,'FontName','myriadpro','fontsize',14);
                title(['% Cont: ' num2str(contamination) ]);
                
                subplot(1,3,2);
                bar(h(1:100));
                title(['Param1: ' num2str(p(count,1)) '; Param2: ' num2str(p(count,2)) ]);
                yvs = get(gca,'ylim');
                fftxt = text(40,0.5*yvs(2),['60 Hz noise: ' num2str(ac_noise) '.']);
                set(fftxt,'FontName','myriadpro','fontsize',14);
                
                subplot(1,3,3);
                plot(normtrace,'k'), hold on
                plot(RSUwin2,'g');
                plot(FSwin2,'r');
                text(1,-0.9,['Amp ' num2str(min(this.meantrace)) ' uV']);
                title(['Rate: ' num2str(meanfr) ' Hz']);
            end
            
        end
    else
        % if empty, do nothing, assign quality 4
        p(count,1:(ncol-1))   = zeros;
        p(count,ncol)     = 4; % ncol is the number of entries in p, so the last entry is always quality
        ogdat.clust(ee).quality = p(count,ncol);
        
        scount = scount + 1;
        score(scount,1) = 4;
        score(scount,2) = 1.0;
        score(scount,3:4) = 0;
        ogdat.clust(ee).qflag = 0;
        ogdat.clust(ee).score = score(scount,:);
    end
end

newblock = ogdat;
