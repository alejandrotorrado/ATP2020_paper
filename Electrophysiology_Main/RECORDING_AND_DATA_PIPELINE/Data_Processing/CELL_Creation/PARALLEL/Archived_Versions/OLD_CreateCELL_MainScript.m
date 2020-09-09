% CLUSTERING AND AUTO-SORTING PIPELINE FOR BINARY SPIKE FILES
%
% Created September 2016 (ATP), as part of data pipeline consolidation
% effort.
%
% *** NOTE: to re-analyze older datasets (for which interpolation is
% required) use binarySpikeProcessing_VINTAGE.m
%
% This script contains the whole data processing pipeline from Binary spike
% files all the way to the creation of a CELL file. The main steps involved
% are clustering of extracted spikes, and automated spike-sorting.
%
% The main use of this code will be to analyze newer datasets  - spikes
% extracted at 161 sample points, with alignment and interpolation.
%
% The only other option to select is the size of the chunks of data on
% which clustering and sorting should be done. This script will allow for
% rapid re-clustering of datasets within specified time limits and in
% specific time chunks. For instance, one could recluster a dataset in
% 8-hour chunks between days 5 and 7 of the recording.

clearIDE;

% which channel to start on? (default should be 1)
start_here = 1;

% User-specified parameters
save_backup = 1;
animal      = input('Animal name?  ','s');
chanidx = input('If this animal was on ch. 1-32, type "A"; if it was on ch. 33-64, type "B":  ','s');
if chanidx == 1 || strcmp(chanidx,'a');
    chanidx = 'A';
elseif chanidx == 2 || strcmp(chanidx,'b');
    chanidx = 'B';
end
dephem = input('Which hemisphere (''L'' or ''R'') was deprived?  ','s');
binsize     = input('How many hours long should the processing chunks be (enter 0 for whole experiment)?:  ');

% Determine processing chunks
if binsize == 0
    startEdge = 0;
    houredges = [];
    nChunks = 0;
else
    startEdge   = input('On what hour do you want to start processing?:  ');
    endEdge     = input('On what hour do you want processing to end?:  ');
    houredges   = startEdge*3600 : 3600*binsize : endEdge*3600;
    nChunks = size(houredges,2) - 1;
end


% Waveform splicing parameters.
nSamples    = 97; % this is hard-coded at the moment. Could also get it from spike files?
fprintf('   *** Number of waveform samples is set to %u. If this is incorrect STOP NOW! ***\n\n',nSamples);
pause(2);
fprintf(' Waveform splicing parameters:\nmin_rng = 10:30; swpr_rng = -18:40;\n\n');
minrng = 10:30;
swprng = -18:40;


% Select the main directory containing the spike files (and binary spike
% files)
sDir = uigetdir(cd,'Pick ALLSPIKES folder');
saveDir = uigetdir(cd,'Pick your saving folder.');
bindatdir = fullfile(sDir,'BinarySpikes');

havebinaries = dir([bindatdir filesep '*_Channel_*_spikes*']);
if ~isempty(havebinaries)
    fprintf('Found %u binary spikes files! Starting processing.\n',numel(havebinaries))
else
    fprintf('     *** No binary spike files found! PROBLEM! ***\n');
    fprintf('Try using binarySpikeProcessing_VINTAGE.m for this dataset.\n');
    keyboard;
end

OGcd = cd;
cd(bindatdir);

% Find the channels for which spikes have been extracted
chantemps = dir ('*_Channel_*');
for ee = 1:size(chantemps,1);
    temp = strfind(chantemps(ee).name,'_');
    chans(ee) = str2double(chantemps(ee).name(temp(2)+1:temp(3)-1));
end
chans = unique(chans);

% Remove reference and EMG channels
skipref_Emg = 1;
if skipref_Emg == 1;
    chans(chans == 8 | chans == 16 | chans == 24 | chans == 32 |...
        chans == 40 | chans == 48 | chans == 56 | chans == 64) = [];
end

% The new spike extraction code already corrects for GMT time. This part
% gets the experiment start time and aligns all spike times to 7:30 am on
% the first day of recording.

% there should always be an experiment start file in the new datasets, but
% let's leave this here for now...
timefile = dir('*_expt_START.mat');
if ~isempty(timefile);
    load(timefile.name);
else
    cd(sDir);
    found_file = 0; file_counter = 0;
    while found_file == 0
        file_counter = file_counter + 1;
        first_spikefile = dir(['*-' num2str(file_counter) '-spikes.mat']);
        if ~isempty(first_spikefile)
            found_file = 1;
        end
    end
    infox = ['INFO' chanidx];
    first_loaded = load(first_spikefile.name,infox);
    tINFO = first_loaded.(infox).startTime;
    save([bindatdir filesep animal '_expt_START.mat'],'tINFO');
    cd(bindatdir);
end
% toffset is experiment start time in seconds from 7:30 am on day 1
toffset = mod(tINFO,24*3600) - (7.5*3600);

% make save directory
cd (saveDir);
autoclustdir = 'AutoClustData';
basename = autoclustdir;

% this little loop is to ensure that each iteration of clustering/sorting
% is saved in a separate folder without overwriting.
bcount = 0;
while exist(autoclustdir,'dir')
    bcount = bcount + 1;
    autoclustdir = [basename '_' num2str(bcount)];
end
mkdir(autoclustdir);
cd(autoclustdir);

% make data structure for saving
CELL_BACKUP = struct();

% Loop through channels
for cc = start_here:size(chans,2);
    
    clear scrubbed* spikes times sweeps minidx
    
    % select channel
    thischan = chans(cc);
    
    % don't do this channel if clu file is already there
    if exist([saveDir filesep autoclustdir filesep 'KlustTempData_' num2str(thischan) '_1.clu.1'],'file')
        fprintf('Channel %u already clustered!\n',thischan);
    else
        
        
        % load spikes
        % all of the following will make sense if you remember that
        % single-precision floats take up 4 bytes of space.
        fprintf('Loading spikes for channel %u.\n',thischan)
        spkfname = dir([bindatdir filesep '*_Channel_' num2str(thischan) '_spikes*']);
        [fid, ~] = fopen( [bindatdir filesep spkfname.name],'r');
        spkfsize = spkfname.bytes;
        if ~isequal(spkfsize/4,round(spkfsize/4))
            warning('WARNING: spike file size inconsistent with single-precision floats.');
        end
        
        % if the spikes file is very large, read in two chunks
        if spkfsize > 2.5e10
            spikes = zeros(nSamples,spkfsize/(nSamples*4));
            fprintf('file too big! Loading in 2 chunks.\n');
            half_bytes = spkfsize/2; % half number of bytes in file
            for half = 1:2
                fprintf('chunk %u\n',half);
                offset_bytes = floor((half-1)*half_bytes);
                A_start = floor((half-1)*size(spikes,2)/2 + 1);
                A_stop  = floor(half*size(spikes,2)/2);
                fseek(fid, offset_bytes, 'bof');
                tempA = fread(fid, [nSamples (A_stop-A_start+1)], 'single');
                
                spikes( : , A_start : A_stop ) = tempA;
                clear tempA
            end
            disp('Done!');
        else
            spikes = fread(fid, 'single=>single');
            spikes = reshape(spikes,nSamples,[]);
            fclose(fid);
            disp('Done!');
        end
        
        if ~isempty(spikes)
            
            % load spike timestamps
            fprintf('Loading spike timestamps for channel %u.\n',thischan)
            tmfname = dir([bindatdir filesep '*_Channel_' num2str(chans(cc)) '_times*']);
            [fid, ~] = fopen( [bindatdir filesep tmfname.name],'r');
            B = [];
            B = fread(fid,'double');
            fclose(fid);
            times = B;
            % align to experiment start time
            times = times - tINFO + toffset; % adjust all spike times to account for offset from 7:30.
            clear B
            disp('Done!');
            
            % index and skip massive events
            disp('Scrubbing the data to remove obvious outliers. ');
            sz              = size( spikes, 1 );
            origindx        = [1:length(times)]';
            [~,poscol]      = find(spikes>400 | spikes <-1200); % Filter out huge traces;
            poscol          = unique(poscol);
            scrubbedidx     = setdiff(origindx,poscol);
            scrubbedspikes  = spikes(:,scrubbedidx);
            clear spikes poscol origindx
            scrubbedtimes   = times(scrubbedidx);
            [nrow,ncol] = size(scrubbedspikes);
            if ncol>nrow;
                scrubbedspikes = scrubbedspikes';
            end
            disp('Outliers removed.');
            
            
            % Do a first alignment to minimum before interpolation.
            
            % center spikes on the minimum
            
            [~,minidx]  = min(scrubbedspikes(:,minrng),[],2);
            minidx      = minidx + (minrng(1) - 1);
            % align
            sweeps = alignrows(scrubbedspikes, minidx, swprng, nSamples - 1);
            clear scrubbedspikes
            
            % identify the spikes and indices that occur in each bin
            if ~isempty(houredges)
                subidx      = cell(nChunks,1);
                subsweeps   = cell(nChunks,1);
                for rr = 1:size(houredges,2)-1;
                    subidx{rr}      = find( scrubbedtimes>=houredges(rr) & scrubbedtimes < houredges(rr+1) );
                    subsweeps{rr}   = sweeps(subidx{rr},:); % This is WF segments from remaining spikes
                end
                clear sweeps
                
                % Parallel loop through chunks
                s0 = tic;
                for yy = 1:nChunks;
                    
                    fprintf('Block %u out of %u.\n',yy,nChunks);
                    
                    % perform clustering
                    if ~isempty(subsweeps{yy})
                    block{yy} = alignandcluster(subsweeps{yy},thischan,yy,houredges,subidx{yy},scrubbedidx,times,0,97);
                    % The last two arguments above are legacy from writing the
                    % function for data sampled at different rates. As of
                    % September 2016, the second-to-last should be a 0
                    % indicating no interpolation is required, the last one
                    % should be 97 (same as hard-coded number of samples that
                    % we get from new spike extraction).
                    
                    fprintf('Performing automated cluster evaluation.');
                    tic
                    newblock{yy} = autoClusterEval_NEW_HPC(block{yy},97,swprng,minrng);
                    % The arguments after the first above are legacy from
                    % writing the function for data sampled at different rates.
                    % As of September 2016, the second one should be 97 (same
                    % as hard-coded number of samples that we get from new
                    % spike extraction), and the third one should be the
                    % variable "swprng" (also hard-coded, see top of script).
                    % The fourth argument, which is optional, was added to
                    % avoid errors in finding the minima of traces. It is hard
                    % coded (see top of script).
                    disp('Done!');
                    toc
                    end
                    
                end
                s1 = toc(s0);
                disp(['It took ' num2str(s1) ' seconds to cluster and sort channel ' num2str(thischan) '.'])
                
                % store data
                CELL_BACKUP(thischan).block = newblock;
                CELL_BACKUP(thischan).channel = thischan;
                
                chdata.block    = newblock;
                chdata.channel  = thischan;
                
                fprintf('Saving data for channel %u...\n',thischan);
                savet0 = tic;
                save([saveDir filesep autoclustdir filesep 'channel_' num2str(thischan) '.mat'],...
                    'chdata','-v7.3');
                savet1 = toc(savet0);
                fprintf('That took %.2f seconds.\n\n',savet1);
                
            else
                % if not chunking, go through whole thing
                s0 = tic;
                subidx = [];
                
                fprintf('UsingKlustaKwik to cluster data.');
                tic
                block{1} = alignandcluster(sweeps,thischan,1,houredges,subidx,scrubbedidx,times,0,97);
                % The last two arguments above are legacy from writing the
                % function for data sampled at different rates. As of
                % September 2016, the second-to-last should be a 0
                % indicating no interpolation is required, the last one
                % should be 97 (same as hard-coded number of samples that
                % we get from new spike extraction).
                disp('Done!');
                toc
                
                
                fprintf('Performing automated cluster evaluation.');
                tic
                newblock = autoClusterEval_NEW_HPC(block{1},97,swprng,minrng);
                % The last two arguments above are legacy from writing the
                % function for data sampled at different rates. As of
                % September 2016, the second-to-last should be 97 (same as
                % hard-coded number of samples that we get from new spike
                % extraction), and the last one should be the variable
                % "swprng" (also hard-coded, see top of script).
                disp('Done!');
                toc
                
                s1 = toc(s0);
                
                disp(['It took ' num2str(s1) ' seconds to cluster and sort channel ' num2str(thischan) '.'])
                
                % store data
                CELL_BACKUP(thischan).channel = thischan;
                CELL_BACKUP(thischan).block = newblock;
                
                chdata.block    = newblock;
                chdata.channel  = thischan;
                
                fprintf('Saving data for channel %u...\n',thischan);
                savet0 = tic;
                save([saveDir filesep autoclustdir filesep 'channel_' num2str(thischan) '.mat'],...
                    'chdata','-v7.3');
                savet1 = toc(savet0);
                fprintf('That took %.2f seconds.\n\n',savet1);
            end
        else
            disp('No spikes for this channel :-( ');
        end
    end
    
end

% If clustering over whole experiment, format CELL file so that instead of
% being organized by channel it is a continuous list of all units
% extracted, excluding clusters of quality 4 (that correspond to noise). If
% clustering over whole experiment this will return a single CELL structure
% containing a list of all cells (excluding 4s) identified in the
% clustering process. If clustering in multiple chunks, it will return a
% structure called CHUNKS that will contain a list of cells (i.e. a CELL
% struct) for each chunk.
if nChunks == 0
    CELL_init = formatAutoCELL(CELL_BACKUP,nChunks);
    CELL = CELL_addFields(CELL_init,animal,dephem,tINFO);
    
    % Save
    disp('Saving your data...');
    tic;
    save([saveDir filesep autoclustdir filesep animal '_MASTER_CELL.mat'],'CELL','-v7.3');
    toc;
    if save_backup
        disp('Saving your backup data...');
        tic;
        save([saveDir filesep autoclustdir filesep animal '_backupCELL.mat'],'CELL_BACKUP','-v7.3');
        toc;
    end
elseif nChunks > 0    
    CHUNKS_init = formatAutoCELL(CELL_BACKUP,nChunks);
    
    CHUNKS = CHUNKS_addFields(CHUNKS_init,animal,dephem,tINFO);
    
    % Save
    disp('Saving your data...');
    tic;
    save([saveDir filesep autoclustdir filesep animal '_MASTER_CHUNKS.mat'],'CHUNKS','-v7.3');
    toc;
    if save_backup
        disp('Saving your backup data...');
        tic;
        save([saveDir filesep autoclustdir filesep animal '_backupCELL.mat'],'CELL_BACKUP','-v7.3');
        toc;
    end
end

% Move all temporary KlustaKwik files into an archive folder
klustarchive = [saveDir filesep autoclustdir filesep 'KlustOutput'];
mkdir(klustarchive);
movefile([saveDir filesep autoclustdir filesep 'KlustTempData*'],klustarchive);

% Copy experiment start file into the CELL output folder
copyfile([bindatdir filesep animal '_expt_START.mat'],[saveDir filesep autoclustdir]);


