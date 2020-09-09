function [spiketimes_A, spikeshapes_A, info_A, spiketimes_B, spikeshapes_B, info_B] = ...
    spikesFromTDTRaw_Upsample(filepath,blocknum,animID,timeon,logfile,extract_these_chans,LGN_YES)
% SPIKESFROMTDTRAW - Read spike shapes from TDT Raw format
%
%  [SPIKETIMES, SPIKESHAPES] = SPIKESFROMTDTRAW(FILEPATH)
%
%  This function will return SPIKETIMES (returned in seconds from the
%  beginning of the recording) and SPIKESHAPES (NUMPTS around each spike time.)
%
%  SPIKETIMES and SPIKESHAPES will be cell arrays, with each cell
%  containing the data from a single channel.
% 
%  TIMEON is the unixtime time stamp of the start of the block
%
%  The data will be read in from disk in chunks of duration
%
% This version (ATP, Sept 2016) is the finalized upsampling + aligning
% version. Based on spikesFromTDTRaw_v10_kbhBeta1.m (written by KBH).

% by default, not an LGN recording
if nargin == 6
    LGN_YES = 0;
end

% this is standard number of sample points in the waveforms extracted from
% the TDT raw data at our sampling rate.
numPts = 32;
wroteToLog = 0;


% The extracted waveforms are up-sampled to improve PCA accuracy
UPSAMPLEFACTOR = 3;

% Get field names in SEV file
TTnames = SEV2mat(filepath,'JUSTNAMES',1,'VERBOSE',0);

spiketimes_A      = {};
spikeshapes_A     = {};
spiketimes_B      = {};
spikeshapes_B     = {};

% NOTE: 'RAWX' is the standard name of TDT raw data in our setup. Change
% this if needed.
% Main body of function
if strcmp(TTnames{1},'RAWX')==1
    xstr = 'RAWX';
    if animID == 1
        % doing animal 1, so output variables for animal 2 are set to 0
        disp(['Extracting channels 1:32 (animID = 1)']);
        spiketimes_B    = 0;
        spikeshapes_B   = 0;
        info_B          = 0;
        
        h0 = tic;
        % timeon is the block start time
        info_A.startTime = double(timeon);
        for ee = extract_these_chans
            % ee is the channel number in this loop
            fprintf('Channel %u\n',ee);
            TT = SEV2mat(filepath,'CHANNEL',ee,'VERBOSE',0);
            % sampling step
            if ee == 1
                dt              = 1/TT.(xstr).fs;
                info_A.OG_dt    = dt;
                info_A.dt       = dt / UPSAMPLEFACTOR;
                
            end
            % check for data in this channel and extract it
            if isfield(TT.(xstr),'data')
                dat_tmp = TT.(xstr).data;
                if ~isempty(dat_tmp)
                    signal_A{ee} = double(dat_tmp(1,:));
                else
                    signal_A{ee} = [];
                end
            else
                wroteToLog = 1;
                fprintf('%s::WARNING! No data found for channel %u, block %u.\n',filepath(5:8),ee,blocknum);
                fid = fopen(logfile,'a');
                fprintf(fid,'%s::WARNING! No data found for channel %u, block %u.\n',filepath(5:8),ee,blocknum);
                fclose(fid);
                signal_A{ee} = [];
            end
            % toggle whether it's an LGN experiment. This is the main
            % subfunction that gets the spikes
            if LGN_YES
                [spiketimes_A{ee}, spikeshapes_A{ee}] = signal2spikes_sdthresh_upsample_LGN(signal_A{ee}, 0, dt, -4, numPts, UPSAMPLEFACTOR);
            else
                [spiketimes_A{ee}, spikeshapes_A{ee}] = signal2spikes_sdthresh_upsample(signal_A{ee}, 0, dt, -4, numPts, UPSAMPLEFACTOR);
            end
            clear TT
        end
        h1 = toc(h0);
        
        disp(['I took ' num2str(h1) ' seconds to extract spikes from block ' num2str(blocknum) ' chans 1:32. (animID = 1)']);
        if wroteToLog, disp(['Something was wrong with the data (animal ' num2str(animID) '! Check the log file.']); end
    elseif animID == 2
        % doing animal 2, so output variables for animal 1 are set to 0
        disp(['Extracting channels 33:64 (animID = 2)']);
        spiketimes_A    = 0;
        spikeshapes_A   = 0;
        info_A          = 0;
        
        h2 = tic;
        % block start time
        info_B.startTime = double(timeon);
        for ee = extract_these_chans
            % ee is channel number
            fprintf('Channel %u\n',ee+32);
            TT = SEV2mat(filepath,'CHANNEL',ee+32,'VERBOSE',0);
            if ee == 1
                dt              = 1/TT.(xstr).fs;
                info_B.OG_dt    = dt;
                info_B.dt       = dt / UPSAMPLEFACTOR;
                
            end
            % find the data
            if isfield(TT.(xstr),'data')
                dat_tmp = TT.(xstr).data;
                if ~isempty(dat_tmp)
                    signal_B{ee} = double(dat_tmp(1,:));
                else
                    signal_B{ee} = [];
                end
            else
                wroteToLog = 1;
                fprintf('%s::WARNING! No data found for channel %u, block %u.\n',filepath(5:8),ee,blocknum);
                fid = fopen(logfile,'a');
                fprintf(fid,'%s::WARNING! No data found for channel %u, block %u.\n',filepath(5:8),ee,blocknum);
                fclose(fid);
                signal_B{ee} = [];
            end
            % get the spikes
            [spiketimes_B{ee}, spikeshapes_B{ee}] = signal2spikes_sdthresh_upsample(signal_B{ee}, 0, dt, -4, numPts, UPSAMPLEFACTOR);
        end
        h3 = toc(h2);
        
        disp(['I took ' num2str(h3) ' seconds to extract spikes from block ' num2str(blocknum) ' chans 33:64. (animID = 2)'])
        if wroteToLog, disp(['Something was wrong with the data (animal ' num2str(animID) '! Check the log file.']); end
    end
    
else
    % output if there is no data
    spiketimes_A    = 0;
    spikeshapes_A   = 0;
    info_A          = 0;
    spiketimes_B    = 0;
    spikeshapes_B   = 0;
    info_B          = 0;
    
end



