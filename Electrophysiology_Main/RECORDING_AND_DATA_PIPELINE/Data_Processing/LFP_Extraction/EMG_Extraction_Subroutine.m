function EMG_Extraction_Subroutine(meta)
%
% Based on ExtractEMG.m
% This version is updated to use .sev files, the new TDT data filetype. It
% is compatible with raw data obtained using the RS4 data streamer.
%
% EMG data is usually in channels 16 and 32. This script extracts filtered
% data from those channels. It also uses the .tsq file (stored in the rig
% PC) to get timing information.
%
% Modified from parent script by ATP, April 2015.
%
% Converted to a callable function by KBH, March 2016. This allows the user
% to automatically link LFP to EMG extraction without having to use
% multiple scripts and inputs.

% set up folder for writing data:
if isfield(meta,'saveDirec')
    mkdir([meta.saveDirec filesep meta.animal],['EMGDATA_' meta.animal]);
else
     mkdir(meta.bigDirec,['EMGDATA_' meta.animal]);
end

% Select the EMG channels relative to the channels recorded from the animal
if strcmp(meta.recBay,'A')
    emgch1 = 16;
    emgch2 = 32;
elseif strcmp(meta.recBay,'B')
    emgch1 = 48;
    emgch2 = 64;
end

% define the bands and rate for filtering (one band value makes it lowpass)
resamp_rate     = 500;
passband        = [10 200];

% retrieve block start times and durations
start_time  = meta.startTimes;
duration    = meta.durations;

for ii = 1:length(meta.listers)
    bn = meta.listers(ii,1);
    
    clear EMGdata
    EMGdata = struct('f_EMG_L',[],'f_EMG_R',[]);

    % directory containing data for this block:
    blockdir = dir([meta.bigDirec filesep 'Block-' num2str(bn)]);
    
    % extract data from SEV file
    %----------------------------------------------------------------------
    tic
    TT_L  = SEV2mat_update([meta.bigDirec filesep 'Block-' num2str(bn)],'CHANNEL',emgch1,'VERBOSE','false');
    f   = toc;
    fprintf('%.2f seconds to access the data on channel %u of block %u.\n',f,emgch1,bn);
        
    
    tic
    TT_R  = SEV2mat_update([meta.bigDirec filesep 'Block-' num2str(bn)],'CHANNEL',emgch2,'VERBOSE','false');
    f   = toc;
    fprintf('%.2f seconds to access the data on channel %u of block %u.\n',f,emgch2,bn);        
    
    info.duration   = duration(ii);
    info.startTime  = start_time(ii);
    info.dt         = 1/meta.samprate;
    
    try
        EMGsignal_L = double(TT_L.RAWX.data(1,:));
    catch
        disp(['unable to get signal from SEV data. Channel ' num2str(emgch1) ', block ' num2str(bn)]);
        keyboard
    end
    
    try
        EMGsignal_R = double(TT_R.RAWX.data(1,:));
    catch
        disp(['unable to get signal from SEV data. Channel ' num2str(emgch2) ', block ' num2str(bn)]);
        keyboard
    end

    %----------------------------------------------------------------------
    % From this point on, this script is identical to its original version,
    % ExtractEMG.m
    
    fEMG_L       = EMGpassFilt(EMGsignal_L, meta.samprate,resamp_rate,passband(2),passband(1),3);
    fEMG_R       = EMGpassFilt(EMGsignal_R, meta.samprate,resamp_rate,passband(2),passband(1),3);
    
    fEMG_L     = fEMG_L - (mean (fEMG_L) );
    fEMG_R     = fEMG_R - (mean (fEMG_R) );
    
    EMGdata.f_EMG_L    = abs(fEMG_L);
    EMGdata.f_EMG_R    = abs(fEMG_R);
    EMGdata.startTime  = start_time(ii);

    clear fEMG_L
    clear fEMG_R
    
    % Consolidate workspace memory.
    if mod(ii,3) == 0;
        % pack performs memory garbage collection.
        pack
    end
    
    if sum(isnan(EMGdata.f_EMG_L))>0;
        keyboard
    else
        filename =[];
        if isfield(meta,'saveDirec')
            filename = ([meta.saveDirec filesep meta.animal filesep 'EMGDATA_' meta.animal filesep [meta.animal '_EMG' num2str(bn)] ]);
        else
            filename = ([meta.bigDirec filesep 'EMGDATA_' meta.animal filesep [meta.animal '_EMG' num2str(bn)] ]);
        end
        save(filename, 'EMGdata');
    end
end


