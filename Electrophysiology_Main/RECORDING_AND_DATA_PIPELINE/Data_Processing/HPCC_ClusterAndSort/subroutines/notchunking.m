function [newblock] =  notchunking (handles, thischan)


saveDir         = handles.saveDir;
autoclustdir    = handles.autoclustdir;
bindatdir       = handles.bindatdir;
nSamples        = handles.nSamples;
toffset         = handles.toffset;
tINFO           = handles.tINFO;
minrng          = handles.minrng;
swprng          = handles.swprng;
newswprng       = handles.newswprng;
doInterpolation = handles.doInterpolation;
interpSamples   = handles.interpSamples;



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
    spikes = zeros(nSamples,spkfsize/(nSamples*4));
    % if the spikes file is very large, read in two chunks
    if spkfsize > 2e10
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
    tmfname = dir([bindatdir filesep '*_Channel_' num2str(thischan) '_times*']);
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
    
    
    % if not chunking, go through whole thing
    s0 = tic;
    subidx = [];
    
    fprintf('UsingKlustaKwik to cluster data.');
    tic
    block{1} = alignandcluster(sweeps,thischan,1,[],subidx,scrubbedidx,times,doInterpolation,interpSamples);
    disp('Done!');
    toc
    
    
    fprintf('Performing automated cluster evaluation.');
    tic
    newblock = autoClusterEval(block{1},interpSamples,newswprng);
    disp('Done!');
    toc
    
    s1 = toc(s0);
    
    disp(['It took ' num2str(s1) ' seconds to process channel ' num2str(thischan) '.'])

    
end
