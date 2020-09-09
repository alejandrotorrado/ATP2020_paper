% RE-PROCESSING PIPELINE FOR SPIKE FILES
%
%
% This script contains the whole data processing pipeline from Binary spike
% files all the way to the creation of a CELL file and automated
% spike-sorting.
%
% The main use of this code will be to re-analyze either old animals (with
% spikes extracted without alignment and interpolation) or newer datasets
% (spikes extracted at 161 sample points, with alignment and
% interpolation).
%
% There are two main options to select:
% __The first is whether the dataset is
%   old or new (see above). This information will be provided by the user by
%   specifying the number of sample points in each waveform.
%
% __The second is the size of the chunks of data on which clustering and
%   sorting should be done. This script will allow for rapid re-clustering
%   of datasets within specified time limits and in specific time chunks.
%   For instance, one could recluster a dataset in 8-hour chunks between
%   days 5 and 7 of the recording.

% clear all, close all, fclose('all'), clc;

% User-specified parameters
animal      = input('Animal name?  ','s');
chanidx = input('If this animal was on ch. 1-32, type "A"; if it was on ch. 33-64, type "B":  ','s');
if chanidx == 1 || strcmp(chanidx,'a');
    chanidx = 'A';
elseif chanidx == 2 || strcmp(chanidx,'b');
    chanidx = 'B';
end
nSamples    = input('How many sample points do your spikes have?:  ');
binsize     = input('How many hours long should the processing chunks be (enter 0 for whole experiment)?:  ');

% Determine processing chunks
if binsize == 0
    startEdge = 0;
    houredges = [];
else
    startEdge   = input('On what hour do you want to start processing?:  ');
    endEdge     = input('On what hour do you want processing to end?:  ');
    houredges   = startEdge*3600 : 3600*binsize : endEdge*3600;
end

% IMPORTANT: this parameter determines by what factor are the waveforms
% upsampled. E.g. if this is 3, the interpolated waveforms will have 3
% times as many samples as the raw data. This new quantity of sample points
% is also stored in a variable for later use.
interpFactor = 3;
interpSamples = 1 + (interpFactor * (nSamples - 1));

% Set flag to determine whether to do spline interpolation (if using
% an older dataset).
% Also set range in which to look for aligning to spike minima.
% Also set range around minimum over which to do PCA and clustering.
% *** NOTE: we may have to change how we do this for datasets that have
% a number of samples that is not 33 or 161.
if nSamples == 33
    doInterpolation = 1;
    minrng = 5:15;
    swprng = -5:15;
    if interpFactor == 3
        newminrng = 10:30;
        newswprng = -18:40;
    end
elseif nSamples == 161
    doInterpolation = 0;
    minrng = 20:50;
    swprng = -30:60;
elseif nSamples == 97
    doInterpolation = 0;
    minrng = 10:30;
    swprng = -18:40;
else
    doInterpolation = 1;
    minrng = ceil(.12*nSamples) : ceil(.31*nSamples);
    swprng = -floor(.19*nSamples) : ceil(.37*nSamples);
    % these ranges based on equivalence with nsamp = 161
end


% Select the main directory containing the spike files (and binary spike
% files)
sDir = uigetdir(cd,'Pick ALLSPIKES folder');
saveDir = uigetdir(cd,'Pick your saving folder.');
bindatdir = [sDir filesep 'BinarySpikes'];
if ~exist(bindatdir,'dir')
    mkdir(bindatdir);
end

havebinaries = dir([bindatdir filesep '*_Channel_*_spikes*']);
if ~isempty(havebinaries)
    fprintf('Found %u binary spikes files! Starting processing.\n',numel(havebinaries))
else
    allchans = 1:32;
    vipchans = setdiff(allchans,[8 16 24 32]);
    fprintf('No binary spike files found. Running script to create them. This will take a while.\n\n');
    binaryspikesbychannel_v3(sDir,animal,allchans,vipchans,chanidx,nSamples);
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

% Correct for GMT time.
% DEAL WITH THIS PER ANIMAL - GMT TIME IS A PAIN IN THE A$$.
% is there a start time file?
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
% GMTcorrect = tINFO - (3600*4);
% toffset = mod(GMTcorrect,24*3600) - (7.5*3600); % this returns the # of sec after 7:30 (lights on) that the recording started.
toffset = mod(tINFO,24*3600) - (7.5*3600);

% make save directory
cd (saveDir);
autoclustdir = 'AutoClustData';
basename = autoclustdir;
bcount = 0;
while exist(autoclustdir,'dir')
    bcount = bcount + 1;
    autoclustdir = [basename '_' num2str(bcount)];
end
mkdir(autoclustdir);
cd(autoclustdir);

% make data structure for saving
CELL = struct();

% if you're analyzing in chunks, follow this
if ~isempty(houredges)
    nChunks = size(houredges,2)-1;
    
    % Loop through channels
    for cc = 1:size(chans,2);
        
        clear scrubbed* spikes times sweeps minidx
        
        % select channel
        thischan = chans(cc);
        
        % don't do this channel if clu file is already there
        if exist([saveDir filesep autoclustdir filesep 'KlustTempData_' num2str(thischan) '_1.clu.1'],'file')
            fprintf('Channel %u already clustered! Skipping to autoeval.\n',thischan);
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
            
            subidx      = cell(size(houredges,2)-1,1);
            subsweeps   = cell(size(houredges,2)-1,1);
            for rr = 1:size(houredges,2)-1;
                subidx{rr}      = find( scrubbedtimes>=houredges(rr) & scrubbedtimes < houredges(rr+1) );
                subsweeps{rr}   = sweeps(subidx{rr},:); % This is WF segments from remaining spikes
            end
            clear sweeps
            
            % Parallel loop through chunks
            s0 = tic;
            % this below should be a PARFOR loop!
            parfor yy = 1:size(houredges,2)-1;
                
                fprintf('Block %u out of %u.\n',yy,size(houredges,2)-1);
                
                % perform clustering
                if size(subsweeps{yy},1) > 0
                    block{yy} = alignandcluster(subsweeps{yy},thischan,yy,houredges,subidx{yy},scrubbedidx,times,doInterpolation,interpSamples);
                    
                    
                    %         figure(yy), hold on
                    %         for ee = 2:size(block(yy).clust,2);
                    %             shadedErrorBar(1:size(subsweeps{yy},2),block(yy).clust(ee).meantrace,block(yy).clust(ee).tracestdev,'-k',1);
                    %             plot(block(yy).clust(ee).meantrace,'linewidth',2);
                    %         end
                    %         hold off
                    
                    
                    fprintf('Performing automated cluster evaluation.\n');
                    tic
                    newblock{yy} = autoClusterEval(block{yy},block{yy}.interpSamples,newswprng,newminrng);
                else
                    block{yy} = [];
                    newblock{yy} = [];
                end
                disp('Done!');
                toc

                
            end
            s1 = toc(s0);
            disp(['It took ' num2str(s1) ' seconds to process channel ' num2str(thischan) '.'])
            
            % store data
            CELL(thischan).block = newblock;
            CELL(thischan).channel = thischan;
            
            
            
            
            
            
        end
        
    end
else
    
    nChunks = 0;
    handles.saveDir         = saveDir;
    handles.autoclustdir    = autoclustdir;
    handles.bindatdir       = bindatdir;
    handles.nSamples        = nSamples;
    handles.toffset         = toffset;
    handles.tINFO           = tINFO;
    handles.minrng          = minrng;
    handles.swprng          = swprng;
    handles.newswprng       = newswprng;
    handles.doInterpolation = doInterpolation;
    handles.interpSamples   = interpSamples;
    % if not chunking, do the truffle shuffle - parfor across channels
    
    for cc = 1:size(chans,2)
        
        %
        %         % don't do this channel if clu file is already there
        %         if exist([saveDir filesep autoclustdir filesep 'KlustTempData_' num2str(thischan) '_1.clu.1'],'file')
        %             fprintf('Channel %u already clustered! Skipping to autoeval.\n',thischan);
        %         else
        
        
        % select channel
        thischan = chans(cc);
        newblock = notchunking (handles, thischan );
        
        channel{cc} = thischan;
        block{cc}   = newblock;
        
        
        %
        %
        %
        %             % load spikes
        %             % all of the following will make sense if you remember that
        %             % single-precision floats take up 4 bytes of space.
        %             fprintf('Loading spikes for channel %u.\n',thischan)
        %             spkfname = dir([bindatdir filesep '*_Channel_' num2str(thischan) '_spikes*']);
        %             [fid, ~] = fopen( [bindatdir filesep spkfname.name],'r');
        %             spkfsize = spkfname.bytes;
        %             if ~isequal(spkfsize/4,round(spkfsize/4))
        %                 warning('WARNING: spike file size inconsistent with single-precision floats.');
        %             end
        %             spikes = zeros(nSamples,spkfsize/(nSamples*4));
        %             % if the spikes file is very large, read in two chunks
        %             if spkfsize > 2e10
        %                 fprintf('file too big! Loading in 2 chunks.\n');
        %                 half_bytes = spkfsize/2; % half number of bytes in file
        %                 for half = 1:2
        %                     fprintf('chunk %u\n',half);
        %                     offset_bytes = floor((half-1)*half_bytes);
        %                     A_start = floor((half-1)*size(spikes,2)/2 + 1);
        %                     A_stop  = floor(half*size(spikes,2)/2);
        %                     fseek(fid, offset_bytes, 'bof');
        %                     tempA = fread(fid, [nSamples (A_stop-A_start+1)], 'single');
        %
        %                     spikes( : , A_start : A_stop ) = tempA;
        %                     clear tempA
        %                 end
        %                 disp('Done!');
        %             else
        %                 spikes = fread(fid, 'single=>single');
        %                 spikes = reshape(spikes,nSamples,[]);
        %                 fclose(fid);
        %                 disp('Done!');
        %             end
        %
        %             % load spike timestamps
        %             fprintf('Loading spike timestamps for channel %u.\n',thischan)
        %             tmfname = dir([bindatdir filesep '*_Channel_' num2str(chans(cc)) '_times*']);
        %             [fid, ~] = fopen( [bindatdir filesep tmfname.name],'r');
        %             B = [];
        %             B = fread(fid,'double');
        %             fclose(fid);
        %             times = B;
        %             % align to experiment start time
        %             times = times - tINFO + toffset; % adjust all spike times to account for offset from 7:30.
        %             clear B
        %             disp('Done!');
        %
        %             % index and skip massive events
        %             disp('Scrubbing the data to remove obvious outliers. ');
        %             sz              = size( spikes, 1 );
        %             origindx        = [1:length(times)]';
        %             [~,poscol]      = find(spikes>400 | spikes <-1200); % Filter out huge traces;
        %             poscol          = unique(poscol);
        %             scrubbedidx     = setdiff(origindx,poscol);
        %             scrubbedspikes  = spikes(:,scrubbedidx);
        %             clear spikes poscol origindx
        %             scrubbedtimes   = times(scrubbedidx);
        %             [nrow,ncol] = size(scrubbedspikes);
        %             if ncol>nrow;
        %                 scrubbedspikes = scrubbedspikes';
        %             end
        %             disp('Outliers removed.');
        %
        %
        %             % Do a first alignment to minimum before interpolation.
        %
        %             % center spikes on the minimum
        %             [~,minidx]  = min(scrubbedspikes(:,minrng),[],2);
        %             minidx      = minidx + (minrng(1) - 1);
        %             % align
        %             sweeps = alignrows(scrubbedspikes, minidx, swprng, nSamples - 1);
        %             clear scrubbedspikes
        %
        %
        %             % if not chunking, go through whole thing
        %             s0 = tic;
        %             subidx = [];
        %
        %             fprintf('UsingKlustaKwik to cluster data.');
        %             tic
        %             block{1} = alignandcluster(sweeps,thischan,1,houredges,subidx,scrubbedidx,times,doInterpolation);
        %             disp('Done!');
        %             toc
        %
        %
        %             fprintf('Performing automated cluster evaluation.');
        %             tic
        %             newblock = autoClusterEval(block{1});
        %             disp('Done!');
        %             toc
        %
        %             s1 = toc(s0);
        %
        %             disp(['It took ' num2str(s1) ' seconds to process channel ' num2str(thischan) '.'])
        %
        % store data
        
        
    end
    
    % now build your CELL structure
    for ee = 1:size(channel,2);
        CELL(ee).channel    =  channel{ee};
        CELL(ee).block      = block{ee};
    end
    
end


% Save
disp('Saving your data...');
tic;
save([saveDir filesep autoclustdir filesep animal '_autoCELL.mat'],'CELL','-v7.3');
toc;


% If clustering over whole experiment, format CELL file so that instead of
% being organized by channel it is a continuous list of all units
% extracted, excluding clusters of quality 4 (that correspond to noise). If
% clustering over whole experiment this will return a single CELL structure
% containing a list of all cells (excluding 4s) identified in the
% clustering process. If clustering in multiple chunks, it will return a
% structure called CHUNKS that will contain a list of cells (i.e. a CELL
% struct) for each chunk.
if nChunks == 0
    CELL = formatAutoCELL(CELL,nChunks);
    
    % Save
    disp('Saving your data...');
    tic;
    save([saveDir filesep autoclustdir filesep animal '_MASTER_CELL.mat'],'CELL','-v7.3');
    toc;
elseif nChunks > 0
    CHUNKS = formatAutoCELL(CELL,nChunks);
    CHUNKS = CHUNKS_addFields(CHUNKS,animal);
    % Save
    disp('Saving your data...');
    tic;
    save([saveDir filesep autoclustdir filesep animal '_MASTER_CHUNKS.mat'],'CHUNKS','-v7.3');
    toc;
end




% Move all temporary KlustaKwik files into an archive folder
klustarchive = [saveDir filesep autoclustdir filesep 'KlustOutput'];
mkdir(klustarchive);
movefile([saveDir filesep autoclustdir filesep 'KlustTempData*'],klustarchive);

% Copy experiment start file into the CELL output folder
copyfile([bindatdir filesep animal '_expt_START.mat'],[saveDir filesep autoclustdir]);



