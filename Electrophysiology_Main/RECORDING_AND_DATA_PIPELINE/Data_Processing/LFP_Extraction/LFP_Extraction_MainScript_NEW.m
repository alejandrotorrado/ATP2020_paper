% LFP_Extraction_MainScript
%
% September 2016 version (as part of data pipeline consolidation).
%
% UPDATE May 2016 (atp):
% We now have tsq timestamps for the start of each block, in the spikes
% files (under INFOA.startTime or INFOB.startTime). This script is adapted
% to extract those for each block. It also passes them to the EMG
% subroutine so they don't have to be calculated again there.
%
%
% In addition to extracting LFPs, this calls EMG_Extraction_Subfunc, which
% extracts EMG data.
%
% Extract band passed LFP (0.3 - 15 Hz) data from multichannel recordings.
% This is currently used for sleep/wake coding via the GUI
% Sleep_Wake_LFP_EMG_Vid_GUI
%
% This version of the script is adapted to be compatible with the RS4 data
% streamer, producing data files in the .SEV format. It uses the
% SEV2mat.m script (in "from./. TDT" folder) to extract data from selected
% channels.
% v10 is an update to make extraction compatible with new data chunking.
% We don't have a tsq timestamp for each block, so we won't be able to get
% duration for each block. However the folder names contain start times, so
% these are available.
%
%  Output: 'LFPinfo.mat' is a structure that is automatically saved to a
%  chosen directory. This file is read by the GUI mentioned above.
%


clearIDE;

% parameters
start_here = 1; % block on which to start LFP extraction
low_freq_band = 0.3; % low frequency bound 
high_freq_band = 15; % high frequency bound


% user inputs:
%   get EMG data too?
emggo = input('Do you want to extract EMG as well? (1=Yes, 0=No):  ');
%   get L/D transition times? NOTE: this will only work if the rig is setup
%   to record these using e.g. a light-sensitive resistor or diode
get_light_times = input('Get L/D transition times for this animal? (1=Yes, 0=No):  ');
%   save in ALLSPIKES folder or somewhere else?
save_in_other = input('Save in different directory? (1=Yes, 0=No):  ');

% pass.samprate is the sampling rate at which the signal was acquired
% resamp_rate is the goal resampling rate for the LFP (to speed things up)
pass.samprate     = 2.44140625e+04;            		          % in Hz
resamp_rate       = 200;
count=0;

% pass.bigDirec in this case should be the main folder in one of raw data
% backup drive containing all the raw SEV data
if ~exist('pass.bigDirec','var')
    pass.bigDirec   = uigetdir('C:\','Pick your backup drive containing the raw data.');
    pass.bigD            = dir(pass.bigDirec);    
end
% if saving somewhere else pick where
if save_in_other
    pass.saveDirec = uigetdir('C:\','Pick your saving directory.');
end
% if you want L/D transition times, point program to the TDT data tank
% (this data will be stored in the TEV files if you are recording it).
if get_light_times
    pass.TDT_Tank = uigetdir('E:\','Pick your TDT data tank on the Rig PC.');
end

% check if ALLSPIKES directory exists. Find it if it does.
spikeDirIdx = find(strcmp('ALLSPIKES',{pass.bigD(:).name}));
if ~isempty(spikeDirIdx)
    pass.spikeDir = [pass.bigDirec filesep pass.bigD(spikeDirIdx).name];
else
    pass.spikeDir = uigetdir('D:\','Could not find ALLSPIKES folder. Point me to it please!');
end

% Select only folders with raw data.
% This presumes the new naming convestion as of early 2016 following the
% switch to continuous datastreams that are broken into 1h blocks
% downstream of the recording. Format is folder title 'Block-N'.
pass.bigD = dir(fullfile(pass.bigDirec,'Block-*'));

% Find animal name automatically
pass.animal = pass.bigDirec(4:end);
clc;
fprintf('Current info: \n\n');
disp(pass);

% Find list of block numbers
% pass.listers is the list of block numbers, as given by folder names in the
% main data directory
count2=0; pass.listers=[];
for ii=1:length(pass.bigD);
    count2=count2+1;
    r=regexp(pass.bigD(ii).name,'-');
    pass.listers(count2,1)=str2double(pass.bigD(ii).name(r+1:end));
    pass.listers(count2,2)=ii;
end
pass.listers=sortrows(pass.listers);

% We use 3 channels of data to extract LFP. Specify which channels here.
% Confirm which bay the animal was recorded on. Request the channels to use
% for LFP extraction, then correct channels if necessary.
pass.chSwitch = input('Which animal? Animal A (ch1-32) or animal B (ch33-64):  ','s');
ch1 = input('What is the first channel you want to use for extraction? ');
ch2 = input('What is the second channel? ');
ch3 = input('What is the last channel? ');

% Make sure everything looks good (check for user input error)
inputcheck = 0;
while inputcheck == 0
    clc;
    if strcmp (pass.chSwitch,'A') || strcmp (pass.chSwitch,'a');
        pass.recBay = 'A';
        channels=[ch1  ch2  ch3];
        inputcheck = 1;
    elseif strcmp(pass.chSwitch,'B') || strcmp (pass.chSwitch,'b');
        pass.recBay = 'B';
        if max([ch1 ch2 ch3])<33;
            channels=[ch1+32  ch2+32  ch3+32];
        else
            channels=[ch1  ch2  ch3];
        end
        inputcheck = 1;
    end
    
    if inputcheck == 0;
        disp('Looks like you made a typo. Try again!');
        pass.chSwitch = input('Which animal? Animal A (ch1-32) or animal B (ch33-64):  ','s');
    end
end

%define the bands for filtering - one value makes it lowpass
LFPs    = [low_freq_band high_freq_band];
info.dt = 1/pass.samprate;
LFPinfo = struct('meanLFP',{},'spectrogram',[]);

% preallocate matrices to hold start times and durations
start_time  = zeros(size(pass.bigD,1),1);
duration    = zeros(size(pass.bigD,1),1);

% variable name for block start time in spikes file (based on animal
% selection)
infox = ['INFO' pass.recBay];

% extract LFP for each block
for tt=start_here:size(pass.listers,1) % tt is the block counter (cycles thru block numbers)
    
    bn = pass.listers(tt,1); % bn is the current block number
    
    % find directory containing data for this block
    blockdir = [pass.bigDirec filesep 'Block-' num2str(bn)];
    files = dir([blockdir filesep '*.sev']);
    filenames = char(files.name);
    dashes = regexp(filenames(1,:),'-');
    h_mark = regexp(filenames(1,:),'h.sev');
    
    % Find block start time from spikes file
    spkFile = dir([pass.spikeDir filesep pass.animal '*-' num2str(bn) '-spikes.mat']);
    infoData = load([pass.spikeDir filesep spkFile.name],infox);
    start_time(tt) = double(infoData.(infox).startTime);
    fprintf('Start time for block %u: %s\n',bn,datestr(unixtime(start_time(tt))));
    
    % Extract LFP signal from raw data
    %----------------------------------------------------------------------
    LFPsig =[];
    for ee=1:length(channels) % THIS IS ONLY 3 CHANNELS!!!!!
        
        % get raw data using TDT function SEV2mat
        % if, for whatever reason, you only want to extract a subset of
        % channels, use the function SEV2mat_update (as an example, this is
        % used in the EMG_Extraction_Subroutine)
        tic
        TT  = SEV2mat([pass.bigDirec filesep 'Block-' num2str(bn)],...
            'CHANNEL',channels(ee),'VERBOSE',0);
        f   = toc;
        fprintf('%.2f seconds to access the data on channel %u of block %u.\n',f,channels(ee),bn);
        
        % convert to double
        try
            signal = double(TT.RAWX.data(1,:));
        catch
            disp(['WARNING: Unable to get signal from SEV data. Channel ' ...
                num2str(channels(ee)) ', block ' num2str(bn)]);
            send_text_message('781-502-6988','att','LFP extraction',...
                ['Problem on ch ' num2str(channels(ee)) ', block ' num2str(bn) '.']);
            keyboard
        end
        
        % resample the LFP and filter it
        signal       = LFPdownfilt(signal,pass.samprate,resamp_rate,LFPs(2),LFPs(1),4);
        if size(signal,2) ~= size(LFPsig,2) && ee>1 && size(LFPsig,2) > 0
             %keyboard;
            if size(signal,2) > size(LFPsig,2)
                LFPsig(ee,:) = signal(1,1:size(LFPsig,2));
            elseif size(LFPsig,2) > size(signal,2)
                padded = padcat(1:size(LFPsig,2),signal);
                LFPsig(ee,:) = padded(2,:);
            end
        else
            LFPsig(ee,:) = signal;
        end
        
    end
    
    % get start_time and duration of this block
    if start_time(tt) ~= 0
        fprintf('Block %u started on: %s\n',bn,datestr(unixtime(start_time(tt))));
        duration(tt) = size(TT.RAWX.data,2)/TT.RAWX.fs;
        fprintf('Block %u duration: %4.2f sec.\n\n',bn,duration(tt));
    end

    % average LFP signal across the 3 extracted channels
    meanLFP     = mean(LFPsig);
    % calculate the actual sampling rate
    Fs          = length(meanLFP)/duration(tt);
    
    % store data in structure
    [S,F,T,P]                   = spectrogram(meanLFP,round(Fs*5),round(Fs*4.5),[LFPs(1) : 0.1 : LFPs(2)],Fs,'yaxis');
    LFPinfo(tt).meanLFP         = meanLFP;
    LFPinfo(tt).spectrogram.S   = S;
    LFPinfo(tt).spectrogram.F   = F;
    LFPinfo(tt).spectrogram.T   = T;
    LFPinfo(tt).spectrogram.P   = P;
    LFPinfo(tt).startTime       = start_time(tt);
    LFPinfo(tt).duration        = duration(tt);
    
    
    close all;
    
    
end

% legacy... i don't think we ever use this anymore? - 6/26/20
if exist('statematrix','var')
    LFPinfo.statematrix=statematrix;
end

% if saving in other directory...
if save_in_other
    mkdir([pass.saveDirec filesep pass.animal]);
    fprintf('Saving your data...\n...\n');
    ts0 = tic;
    save([pass.saveDirec filesep pass.animal filesep pass.animal '_LFPinfo.mat'],'LFPinfo','-v7.3')
    ts1 = toc(ts0);
    fprintf('That took %.2f seconds.\n\n\n',ts1);
else
    save([pass.bigDirec filesep pass.animal '_LFPinfo.mat'],'LFPinfo','-v7.3')
end


% If selected, move on to EMG extraction.
if emggo == 1
    % 'pass' contains all the arguments to pass to the EMG subroutine
    pass.startTimes = start_time;
    pass.durations = duration;
    EMG_Extraction_Subroutine(pass);
end

% if extracting LD transitions
if get_light_times
    disp('Getting L/D transition times.');
    % use the custom function to do this
    LD_data = get_LD_times(pass.TDT_Tank,pass.animal);
    
    % fix bad times - sometimes the rig would record a bunch of transition
    % times a few msec from each other. To deal with this, remove
    % transition times that are unexpectedly close to each other (<10 msec)
    on = LD_data.ON_TIMES;
    off = LD_data.OFF_TIMES;
    d_on = diff(on);
    bad_d_on = find(d_on < 0.01) + 1;
    on(bad_d_on) = [];
    d_off = diff(off);
    bad_d_off = find(d_off < 0.01) + 1;
    off(bad_d_off) = [];
    
    
    fprintf('DONE! Make sure to review LD times before saving!');
    if save_in_other
        ld_save = [pass.saveDirec filesep pass.animal filesep pass.animal '_LDdata.mat'];
    else
        ld_save = [pass.bigDirec filesep pass.animal '_LDdata.mat'];
    end
     keyboard;
    
    % replace with the list of good timestamps
    LD_data.ON_TIMES    = [];
    LD_data.ON_TIMES    = on;
    LD_data.OFF_TIMES   = [];
    LD_data.OFF_TIMES   = off;
    save(ld_save,'LD_data','-v7.3');
end
    
