function makeSpikeBinaries(blockspikes,blocktimes,blockinfo,inputdir,channels,VIPchans,animal,nsamps)
% WRITE DOC!


% make directory for saving
bindir = [inputdir filesep 'BinarySpikes'];
if ~exist(bindir,'dir');
    disp('Making BinarySpikes folder.');
    mkdir(bindir)
end
addpath(bindir);


% Create binary spikes and time files for all channels
havebinaries = dir([bindir filesep '*.bin']);
if isempty(havebinaries)
    fprintf('Creating binary spike files for the first time.\n');
    for ii=1:length(channels);
        fname = [];
        fname = [animal '_Channel_' num2str(channels(ii)) '_spikes.bin'];
        [fid, ~] = fopen( [bindir filesep fname],'w');
        fclose(fid);
        
        fnameT = [];
        fnameT = [animal '_Channel_' num2str(channels(ii)) '_times.bin'];
        [fidT, ~] = fopen( [bindir filesep fnameT],'w');
        fclose(fidT);
    end
else
    fprintf('Found binary files in BinarySpikes directory!\n');
end

% For current block, get spikes for each channel (parallelize)
parfor rr = 1 : length(VIPchans);
    
    
    fprintf('Working on channel %u.\n',VIPchans(rr));
    
    
    %% ------------------------------ SPIKES ------------------------------
    
    fname   = [];
    fname   = [animal '_Channel_' num2str(VIPchans(rr)) '_spikes.bin'];
    %         [fid, message] = fopen( [inputdir filesep 'BinarySpikes' filesep fname],'a');
    [fid, ~] = fopen( [bindir filesep fname],'a');
    
    
    spkdat = [];
    try
        % 05/04/2015: ATP added factor of 1e6 for new rig setup (data
        % was in units of V instead of uV and was messing up the
        % analysis)
        
        spkdat = 1e6.*blockspikes{VIPchans(rr)};
        
    catch
        disp('caught at 63 - makeSpikeBinaries');
        keyboard
    end
    
    % Convert to single-precision float to decrease memory space required
    % for storage
    spkdat = single(spkdat);
    
    % if spkdat is empty, just skip ahead and write an empty binary file
    if ~any(size(spkdat) == nsamps);
        warning(['Spike data matrix on channel %u does not have a size '...
            'that reflects sample size (npoints). Keyboarding.'],VIPchans(rr));
        if ~isempty(spkdat) 
            keyboard;
        end
    end
    
    if size(spkdat,2) == nsamps; % update 9/19: passed number of samples as argument. OLD NOTE: hard-coded number of samples in spikes - pass this as argument?
        spkdat = spkdat';
    end
    
    % Append to binary spike file
    if ~isempty(spkdat)
    fwrite(fid,spkdat,'single');
    fclose(fid);
    end
    
    
    %% ------------------------------ TIMES -------------------------------
    
    fnameT = [];
    fnameT = [animal '_Channel_' num2str(VIPchans(rr)) '_times.bin'];
    [fidT, ~] = fopen( [inputdir filesep 'BinarySpikes' filesep fnameT],'a');
    
    tmdat = [];
    try
        if isfield(blockinfo,'startTime');
            tmdat = blocktimes{VIPchans(rr)} + blockinfo.startTime;
        else
            warning('makeSpikeBinaries: could not find startTime of block. Keyboarding.');
            keyboard;
        end
    catch
        disp('Caught at 136')
        keyboard
    end
    
    % Append to binary times file
    fwrite (fidT,tmdat,'double');
    fclose(fidT);
    
    if size(spkdat,2)~=length(tmdat);
        warning('The number of spike and time entries are not equal.');
        keyboard
    end
    
end

