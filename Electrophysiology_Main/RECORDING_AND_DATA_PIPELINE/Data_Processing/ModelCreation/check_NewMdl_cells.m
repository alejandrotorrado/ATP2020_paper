 scount = 0; pcount = 0; chcount = 0;
 plist = [88,89,90,91,92,93,131,180,184,185,199];

for ii =  1:28
    
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
            pcount = pcount + 1;
            
            if any(plist == pcount);
                disp(pcount); pause(.2)
                % this is data for current cluster
                this = chdat.clust(ee);
                
                dat = [];
                dat = this.time;
                X = diff(dat);
                h = histc(X,[0:0.001:10.0]); % previously calculated only out to 0.5 sec
                % This changed on 3/27/17 to 10 sec. This prevents over-representation
                % of noise in <3ms for slowly firing cells.
                
                if sum(h)>1000
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
%                     p(pcount,1:2) = gamfit(1:100,[],[],h(1:100)');
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
                
                label = predict(Mdl,p(pcount,1:(ncol-1)));
                    
                    g = figure(pcount);
                    
                    set(g,'units','normalized','position',[0.05 0.05 0.85 0.85]);
                    s1 = subplot(1,3,1);
                    plot(ygam2,'linewidth',2,'color',[.65 .09 .18]);
                    yvals = get(gca,'ylim');
                    if model_scored
                        txt = text(35,yvals(2)*0.75,['Predicted q: ' label{1} ]);
                        txt2 = text(35,yvals(2)*0.70,['New assigned q: ' num2str(p(pcount,ncol)) ]);
                        set(txt2,'fontsize',14);
                    else
                        txt = text(35,yvals(2)*0.75,['Previous q: ' num2str(this.quality) ]);
                        txt2 = text(35,yvals(2)*0.70,['New assigned q: ' num2str(p(pcount,ncol)) ]);
                        set(txt2,'fontsize',14);
                    end
                    if showsavedscore
                        txt2 = text(35,yvals(2)*0.65,['Saved q: ' num2str(p_save{chcount}(pcount,15))]);
                        set(txt2,'fontsize',14);
                    end
                    set(txt,'fontsize',14);
                    set(gca,'fontsize',14);
                    title(sprintf('Total %% cont: %.3f\nMean hourly %% cont: %.3f',contamination,contamination_byhour),'fontsize',14);
                    
                    subplot(1,3,2);
                    bar(h(1:100));
                    title(['Param1: ' num2str(p(pcount,1)) '; Param2: ' num2str(p(pcount,2)) ],'fontsize',14);
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
                    
                    g1=figure(pcount+1);
                    set(g1,'units','normalized','position',[0.2 .1 .6 .8]);
                    subplot(2,1,1); hold on;
                    binsz = 600;
                    myrate = histc(this.time,0:binsz:264*3600)./binsz;
                    plot(myrate,'linewidth',2);
                    set(gca,'xlim',[0 264*3600/binsz]);
                    subplot(2,1,2);
                    plot(1:length(hourly_cont),hourly_cont,'color',[.93 .69 .13],...
                        'linewidth',2);
                    set(gca,'ylim',[0 8],'xlim',[0 216]);
                    uiwait(g);
                end
            end
        end
    end
end

keyboard;
modify_list = [106];
for aa = 1:length(modify_list)
    px = modify_list(aa);
    fprintf('label at pcount %u is:  %u.\n',px,p(px,ncol));
    newq = input('Enter new quality:  ');
    p(px,ncol) = newq;
    fprintf('label at pcount %u is now:  %u.\n',px,p(px,ncol));
end
