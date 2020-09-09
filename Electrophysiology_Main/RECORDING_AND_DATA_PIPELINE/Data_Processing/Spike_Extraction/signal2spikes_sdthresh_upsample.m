function [spiketimes, sweeps] = signal2spikes_sdthresh_upsample(signal, signal_starttime, signal_sampleinterval, sdthresh, numPts, upfactor)
% Extracts spike times & shapes from a signal w/ SD threshold
%
%     [SPIKETIMES, SPIKESHAPES] = SIGNAL2SPIKES(SIGNAL, SIGNAL_STARTTIME,
%          SIGNAL_INTERVAL, SDTHRESH, NUMPTS)
%
%  Examines a recording, calculates the standard deviation, and then
%  determines spike times based on negative threshold crossings of SDTHRESH
%  times the standard deviation (e.g., -4).
%
%  Inputs:
%      SIGNAL - the samples to be analyzed
%      SIGNAL_STARTTIME - the time (in seconds) of the first sample
%      SIGNAL_SAMPLEINTERVAL - The signal sample interval (in seconds)
%      SDTHRESH - Threshold will be this number multipled by standard deviation
%      NUMPTS - The number of points around each spike to extract
%
% This is the finalized version that does the upsampling (ATP, sept 2016).
% Based on signal2spikes_sdthresh_v2_KBHbeta1.m (written by KBH).

% make sure signal array is not empty
if isempty(signal)
    spiketimes = [];
    sweeps = [];
   warning('FOUND AN EMPTY CHANNEL. No signal is extracted. Processing will skip this and move on.');  
else
    
    % bandpass filter the data. The frequency bands are hard-coded in the
    % custom function below. Right now they are 500-4000 Hz.
    newsignal   = highlowpassfilter(signal);
    
    % upsample the signal
    UPsignal_sampleinterval = signal_sampleinterval/upfactor;
    upnpts      = numPts*upfactor;
    x           = 1:size(newsignal,2);
    xq          = 1:(1/upfactor):size(newsignal,2);
    newsignal   = interp1(x,newsignal,xq,'spline');
    
    
    xtrafat = round(numPts*6.25);
    
    
    thresh          = sdthresh.*std(newsignal);
    past            = newsignal<thresh;
    hits            = 0.*newsignal;
    hits(past)      = 1;
    hits            = [ 0 diff(hits) ];
    hits(hits<0)    = 0;    % finds all points where there is a crossing
    
    timeindx        = find(hits);  % the sample numbers within the chunk where there are spikes
    
    myspiketimes    = (signal_starttime + (timeindx-1) * UPsignal_sampleinterval)'; % these are the times in seconds
    
    nspikes                = length(timeindx);
    spikeshapes            = []; % not this  zeros( numPts, nspikes );
    spiketimes = [];
    
    
    %set this up to do 2/3 of num points AFTER the thresh crossing -
    % collect more of the refractory period.
    
    killtimes1          = find(timeindx<xtrafat/4 + 1);
    killtimes2          = find(timeindx>length(newsignal)-xtrafat*(3/4));
    killtimes           = unique([killtimes1 killtimes2]);
    
    
    timeindx(killtimes)     = [];
    myspiketimes(killtimes) = [];
    
    spiketimes          = myspiketimes';
    
    if ~isempty(spiketimes);
        
        tstarts(:,1)        = timeindx - xtrafat/4 - 1;
        
        var     = repmat(1:xtrafat+1,[size(myspiketimes),1]);
        var2    = repmat(tstarts(:,1),1,xtrafat+1);
        VAR     = var + var2;
        %VAR     = VAR';
        
        % this version of spikeshapes has 25% too many points - this allows for
        % centering on the minimum
        spikeshapes = newsignal(VAR)';
        
        % now find the minimum index for each spike, and define the size of the
        % desired spike around the minimum (swprng)
        [~,minidx]  = min(spikeshapes(floor(upnpts/5):floor(upnpts/2),:),[],1);
        minidx      = minidx + floor(upnpts/5)-1;
        swprng      = -(upnpts/4):(upnpts*(3/4));
        
        % align the spikeshapes on minimum and return sweeps, which is the
        % desired spike matrix with interpolation at a factor of UPFACTOR.
        sweeps = alignrows(spikeshapes', minidx, swprng, xtrafat);
        
        %  ----------------------------------------------
        % GET RID OF GIANT OUTLIERS
        % index and skip massive events
        
        [posrow,~]          = find(sweeps>3e-4 | sweeps <-7e-4); % Filter out huge traces;
        posrow              = unique(posrow);
        sweeps(posrow,:)    = [];
        spiketimes(posrow)  = [];
        %  ----------------------------------------------
        
        if size(sweeps,1) ~= length(spiketimes);
            warning('The number of spikes and times detected are not equal!!!signal2spikes_sdthresh_v2_KBHbeta1 ');
            keyboard
        end
        
    else
        sweeps      = [];
    end
end
