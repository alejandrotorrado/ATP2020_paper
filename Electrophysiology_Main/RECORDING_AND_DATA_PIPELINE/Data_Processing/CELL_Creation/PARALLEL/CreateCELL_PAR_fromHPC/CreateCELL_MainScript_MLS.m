%% CreateCELL_MainScript_MLS
%
% MLS: Machine Learning Sorted
% This script uses a Random Forest algorithm to assign a quality score to
% clusters.
%
% This script contains the whole data processing pipeline from Binary spike
% files all the way to the creation of a CELL file and automated
% spike-sorting.
%

clearIDE;

cd('Z:\RECORDING_AND_DATA_PIPELINE\Data_Processing\CELL_Creation\PARALLEL\CreateCELL_PAR_fromHPC');
addpath(genpath(cd));



% User-specified parameters
save_backup = 1;
animal      = input('Animal name?  ','s');
chanidx = input('If this animal was on ch. 1-32, type "A"; if it was on ch. 33-64, type "B":  ','s');
if chanidx == 1 || strcmp(chanidx,'a');
    chanidx = 'A';
elseif chanidx == 2 || strcmp(chanidx,'b');
    chanidx = 'B';
end
% which channel to start on? (default should be 1)
chanstart = input('Which channel do you want to start on?  ');
dephem = input('Which hemisphere (''L'' or ''R'') was deprived?  ','s');
nSamples = 97;
fprintf('\n     Variable nSamples is set to %u.\nIf this is NOT the number of sample points in your waveforms, STOP NOW!\n\n',nSamples);
binsize     = input('How many hours long should the processing chunks be (enter 0 for whole experiment)?:  ');

% Determine processing chunks
if binsize == 0
    startEdge = 0;
    endEdge = 0;
    houredges = [];
    nChunks = 0;
else
    startEdge   = input('On what hour do you want to start processing?:  ');
    endEdge     = input('On what hour do you want processing to end?:  ');
    houredges   = startEdge*3600 : 3600*binsize : endEdge*3600;
    nChunks = size(houredges,2) - 1;
end


% IMPORTANT: interpFactor parameter determines by what factor are the waveforms
% upsampled. E.g. if this is 3, the interpolated waveforms will have 3
% times as many samples as the raw data. This new quantity of sample points
% is also stored in a variable for later use.

% Set flag to determine whether to do spline interpolation (if using
% an older dataset).
% Also set range in which to look for aligning to spike minima.
% Also set range around minimum over which to do PCA and clustering.
% *** NOTE: we may have to change how we do this for datasets that have
% a number of samples that is not 33 or 161.
interpFactor = 3;
if nSamples == 33
    doInterpolation = 1;
    minrng = 5:15;
    swprng = -5:15;
    if interpFactor == 3
        newminrng = 10:30;
        newswprng = -18:40;
    end
    interpSamples = 97;
elseif nSamples == 161
    doInterpolation = 0;
    minrng = 25:80;
    swprng = -30:60;
    newminrng = 20:50;
    newswprng = swprng;
    interpSamples = 161;
elseif nSamples == 97
    doInterpolation = 0;
    minrng = 15:50; % this is minimum range in the 1:97 WFs
    swprng = -18:40;
    newminrng = 10:30; % this is minimum range in the 1:59 WFs
    newswprng = swprng;
    interpSamples = 97;
else
    doInterpolation = 1;
    minrng = ceil(.12*nSamples) : ceil(.31*nSamples);
    swprng = -floor(.19*nSamples) : ceil(.37*nSamples);
    % these ranges based on equivalence with nsamp = 161
end

fprintf('   *** Number of waveform samples is set to %u. If this is incorrect STOP NOW! ***\n\n',nSamples);
pause(2);
fprintf(' Waveform splicing parameters:\nmin_rng = %u:%u; swprng = -%u:%u;\n\n',minrng(1),minrng(end),abs(swprng(1)),swprng(end));
fprintf(' NEW Waveform splicing parameters:\nnewmin_rng = %u:%u; newswp_rng = -%u:%u;\n\n',newminrng(1),newminrng(end),abs(newswprng(1)),newswprng(end));


% Number of clusters
nclust = 9;


% Select the main directory containing the spike files (and binary spike
% files)
sDir = uigetdir(cd,'Pick ALLSPIKES folder');
saveDir = uigetdir(cd,'Pick your saving folder.');
bindatdir = fullfile(sDir,'BinarySpikes');
if ~exist(bindatdir,'dir')
    mkdir(bindatdir);
end

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
for ee = 1:size(chantemps,1)
    temp = strfind(chantemps(ee).name,'_');
    chans(ee) = str2double(chantemps(ee).name(temp(2)+1:temp(3)-1));
end
chans = unique(chans);

% Remove reference and EMG channels
skipref_Emg = 1;
if skipref_Emg == 1
    chans(chans == 8 | chans == 16 | chans == 24 | chans == 32 |...
        chans == 40 | chans == 48 | chans == 56 | chans == 64) = [];
end

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
bcount = 0;
while exist(autoclustdir,'dir')
    bcount = bcount + 1;
    autoclustdir = [basename '_' num2str(bcount)];
end
mkdir(autoclustdir);
cd(autoclustdir);


% Load random forest model trained on previous datasets
rForest_dir = 'Z:\RECORDING_AND_DATA_PIPELINE\Data_Processing\CELL_Creation\RANDOM_FOREST_TRAINED\LATEST';
fprintf('Looking for Random Forest model in folder:\n%s\n',rForest_dir);
rForest_mdl_file = dir([rForest_dir filesep '*.mat']);
if numel(rForest_mdl_file) > 1
    % Ensure there is only one Mdl file available (latest one)
    error('PROBLEM: more than one rForest_Mdl file found!');
    
else
    Mdl_load = load([rForest_dir filesep rForest_mdl_file(1).name]);
    Mdl = Mdl_load.Mdl;
end

% make data structure for saving
CELL_BACKUP = struct();

% initialize temporary variables for parfor loop
cellarray_size = size(chans,2);
nclust_trim = cell(1,cellarray_size);
trimmers = cell(1,cellarray_size);
idx_tmp = cell(1,cellarray_size);
p_save      = cell(1,cellarray_size);

if isempty(gcp)
    poolobj = parpool('local');
end

% IF CHUNKING, PARALLELIZE AT CHUNKS
if ~isempty(houredges)
    
    % normal loop through channels
    for cc = chanstart:size(chans,2)
        
        % select channel
        thischan = chans(cc);
        
        % Loading, cleaning and aligning spikes
        [sweeps,scrubbedtimes,times,scrubbedidx,nspk_og] = ...
            loadspikes_hpc(bindatdir,thischan,nSamples,tINFO,minrng,swprng);
        
        
        fprintf('Using KlustaKwik to cluster data. CHUNKED CASE.\n');
        
        nChunks = size(houredges,2)-1;
        % identify the spikes and indices that occur in each bin
        
        subidx          = cell(nChunks,1);
        subsweeps       = cell(nChunks,1);
        p_save_chunk    = cell(nChunks,1);
        clip_combine    = cell(nChunks,1);
        
        for rr = 1:nChunks
            subidx{rr}      = find( scrubbedtimes>=houredges(rr) & scrubbedtimes < houredges(rr+1) );
            subsweeps{rr}   = sweeps(subidx{rr},:); % This is WF segments from remaining spikes
        end
        clear sweeps
        
        % Parallel loop through chunks
        s0 = tic;
        % this below should be a PARFOR loop!
        for yy = 1:nChunks
            
            fprintf('Block %u out of %u.\n',yy,size(houredges,2)-1);
            yyt0=tic;
            % perform clustering
            %             if size(subsweeps{yy},1) > 0
            if ~isempty(subsweeps{yy})
                % Use the supervised learning model for autosorting:
                
                fprintf('UsingKlustaKwik to cluster data.\n');
                tic
                
                p_out = []; trimvec = [];
                [block{yy}, p_out, trimvec] = alignandcluster_HPC(subsweeps{yy},...
                    thischan,yy,houredges,subidx{yy},scrubbedidx,times,...
                    doInterpolation,interpSamples,...
                    Mdl,newminrng,newswprng,nclust);
                % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                % Check for a trim flag and set up for further (in)action:
                trimstat = 0;
                if sum(trimvec) == 0
                    tflag = 4;
                    
                else
                    nclust_trim{yy} = nclust + 1;
                    trimmers{yy} = find(trimvec == 1);
                    tflag = 0;
                end
                
                thresh  = 0.5; % this is mostly cutting noise spikes that are _way_ out there.
                while tflag < 2 % trim twice (change this number to change number of passes)
                    % Set up for trimming the problematic clusters
                    trimstat    = trimstat + 1;
                    tflag       = tflag+1;
                    clip = []; killz = []; % Assemble all of the points to be trimmed:
                    
                    for vv = 1:sum(trimvec)
                        clust = []; dist = [];
                        clust   = block{yy}.clust(trimmers{yy}(vv));
                        dist    = clust.dists./max(clust.dists);
                        tmpkillz   = find(dist > thresh);
                        killz = clust.idx(tmpkillz);
                        clip    = [clip; killz];
                        tmpkillz = [];
                    end
                    
                    % build an idx list
                    idx_tmp_par = [];
                    for jj = 1:size(block{yy}.clust,2);
                        idx_tmp_par = [idx_tmp_par; block{yy}.clust(jj).idx];
                    end
                    idx_tmp{yy} = sortrows(idx_tmp_par);
                    
                    % build an idx list
                    idx_tmp_par = [];
                    for jj = 1:size(block{yy}.clust,2);
                        idx_tmp_par = [idx_tmp_par; block{yy}.clust(jj).idx];
                    end
                    idx_tmp{yy} = sortrows(idx_tmp_par);
                    % This needs to be built every time. The clust.idx indices
                    % refer to the binary spike file, but if some of them have been
                    % trimmed, they should be excluded from the dataset that is
                    % being trimmed again.
                    
                    clip_combine{yy} = [clip_combine{yy}; clip];% combine all WFs to be clipped
                    idx_tmp2 = idx_tmp{yy}; % make copy of array containing all indices
                    % NOTE: these correspond to indices in binary spike file!
                    % find indices of WFs to clip. Using intersect in this way:
                    % [~,IX,~] = intersect(A,B) returns IX where A(IX) = B.
                    % So I set A=idx_tmp2 (indices corresponding to spike file)
                    % B = clip_combine (indices in spike file to clip)
                    % and I get IX so that A(IX)=[] clips away the correct indices
                    % in idx_tmp2.
                    [~,clip_idx,~] = intersect(idx_tmp2,clip_combine{yy});
                    idx_tmp2(clip_idx) = []; % remove WFs to be trimmed
                    % NOTE: variable sweeps may have some spikes removed.
                    % Therefore, its indices do not cover the span of the whole
                    % binary file, i.e. instead of going 1:n_spikes they go
                    % 1:(n_spikes-x).
                    % To take this into account, we have to refer to the variable
                    % "scrubbedidx" that correctly indices sweeps. The indices
                    % contained in idx_tmp2 but corresponding to "sweeps" are
                    % contained in idx_tmp3.
                    
                    [~,idx_tmp3,~] = intersect(subidx{yy},idx_tmp2);
                    sweepsCLIPPED           = subsweeps{yy}(idx_tmp3,:);
                    subidxCLIPPED			= subidx{yy}(idx_tmp3,:);
                    
                    % Delete clustering files associated with this channel,
                    % otherwise it may kick back a message that this
                    % channel is already done.
                    delete([saveDir filesep autoclustdir filesep '*KlustTempData_' num2str(thischan) '_' num2str(yy) '.*.1']);
                    
                    fprintf('Channel %u is being trimmed!\n',thischan);
                    [block{yy}, p_out, trimvec] = alignandcluster_HPC(sweepsCLIPPED,...
                        thischan,yy,houredges,subidxCLIPPED,scrubbedidx,times,...
                        doInterpolation,interpSamples,...
                        Mdl,newminrng,newswprng,nclust_trim{yy});
                    
                    % exit the loop if nothing else needs to be trimmed
                    if sum(trimvec) == 0
                        tflag = 4;
                    else
                        trimmers{yy} = find(trimvec==1);
                    end
                    
                end
                disp('Done!');
                toc
                
                % store each chunk's p in p_save_chunk{chunk}
                p_save_chunk{yy} = p_out;
                
                
                newblock{yy} = block{yy};
                newblock{yy}.trimstat = trimstat;
                
            else
                block{yy} = [];
                newblock{yy} = [];
            end
            disp('Done!');
            yyt1 = toc(yyt0);
            fprintf('Block %u processing took %.2f seconds.\n',yy,yyt1);
            
            
        end
        s1 = toc(s0);
        disp(['It took ' num2str(s1) ' seconds to process channel ' num2str(thischan) '.'])
        
        % store this channel's p (all the chunks) in p_save{channel}
        p_save{cc} = p_save_chunk;
        
        % store data in main CELL variable
        CELL_BACKUP(cc).block = newblock;
        CELL_BACKUP(cc).channel = thischan;
        CELL_BACKUP(cc).psave = p_save{cc};
        
        % save data for each channel
        save_chdata(newblock,thischan,saveDir,autoclustdir);
        
    end
    
    % save the complete p variable
    save([saveDir filesep autoclustdir filesep animal '_sorting_p.mat'],'p_save','-v7.3');
    
    
else % IF NOT CHUNKING, PARALLELIZE ACROSS CHANNELS
    
    clip_combine = cell(1,cellarray_size);
    
    parfor cc = chanstart:size(chans,2)
        
        thischan = chans(cc);
        
        % Loading, cleaning and aligning spikes
        [sweeps,scrubbedtimes,times,scrubbedidx,nspk_og] = ...
            loadspikes_hpc(bindatdir,thischan,nSamples,tINFO,minrng,swprng);
        
        
        s0 = tic;
        subidx = [];
        
        if ~isempty(sweeps)
            fprintf('UsingKlustaKwik to cluster data. CONTINUOUS CASE.\n');
            tic
            
            p_out = []; trimvec = [];
            [block, p_out, trimvec] = alignandcluster_HPC(sweeps,thischan,1,...
                houredges,subidx,scrubbedidx,times,...
                doInterpolation,interpSamples,...
                Mdl,newminrng,newswprng,nclust);
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Check for a trim flag and set up for further (in)action:
            trimstat = 0;
            if sum(trimvec) == 0
                tflag = 4;
                
            else
                nclust_trim_tmp = nclust + 1;
                nclust_trim{cc} = nclust_trim_tmp;
                trimmers{cc} = find(trimvec == 1);
                tflag = 0;
            end
            
            thresh  = 0.5; % this is mostly cutting noise spikes that are _way_ out there.
            while tflag < 2 % trim twice (change this number to change number of passes)
                % Set up for trimming the problematic clusters
                trimstat    = trimstat + 1;
                tflag       = tflag+1;
                clip = []; killz = []; % Assemble all of the points to be trimmed:
                
                for vv = 1:sum(trimvec)
                    clust = []; dist = [];
                    clust   = block.clust(trimmers{cc}(vv));
                    dist    = clust.dists./max(clust.dists);
                    tmpkillz   = find(dist > thresh);
                    killz = clust.idx(tmpkillz);
                    clip    = [clip; killz];
                    tmpkillz = [];
                end
                
                % build an idx list
                idx_tmp_par = [];
                for jj = 1:size(block.clust,2);
                    idx_tmp_par = [idx_tmp_par; block.clust(jj).idx];
                end
                idx_tmp{cc} = sortrows(idx_tmp_par);
                % This needs to be built every time. The clust.idx indices
                % refer to the binary spike file, but if some of them have been
                % trimmed, they should be excluded from the dataset that is
                % being trimmed again.
                
                clip_combine{cc} = [clip_combine{cc}; clip]; % combine all WFs to be clipped
                idx_tmp2 = idx_tmp{cc}; % make copy of array containing all indices
                % NOTE: these correspond to indices in binary spike file!
                % find indices of WFs to clip. Using intersect in this way:
                % [~,IX,~] = intersect(A,B) returns IX where A(IX) = B.
                % So I set A=idx_tmp2 (indices corresponding to spike file)
                % B = clip_combine (indices in spike file to clip)
                % and I get IX so that A(IX)=[] clips away the correct indices
                % in idx_tmp2.
                [~,clip_idx,~] = intersect(idx_tmp2,clip_combine{cc});
                idx_tmp2(clip_idx) = []; % remove WFs to be trimmed
                % NOTE: variable sweeps may have some spikes removed.
                % Therefore, its indices do not cover the span of the whole
                % binary file, i.e. instead of going 1:n_spikes they go
                % 1:(n_spikes-x).
                % To take this into account, we have to refer to the variable
                % "scrubbedidx" that correctly indices sweeps. The indices
                % contained in idx_tmp2 but corresponding to "sweeps" are
                % contained in idx_tmp3.
                [~,idx_tmp3,~] = intersect(scrubbedidx,idx_tmp2);
                % using idx_tmp3 we can now remove the correct spikes from
                % sweeps and obtain sweepsCLIPPED (trimmed dataset).
                sweepsCLIPPED = [];
                sweepsCLIPPED = sweeps(idx_tmp3,:);
                
                %BIG NOTE:
                % Essentially the trick is:
                % If i have a set of indices x1 in array A and an array B
                % containing a subset of those indices; to find the indices x2
                % such that A(x2) = B i.e. find the positions in array A of my
                % subset of indices contained in array B, the intersect method
                % above is the best way.
                
                
                
                % Delete clustering files associated with this channel,
                % otherwise it may kick back a message that this
                % channel is already done.
                delete([saveDir filesep autoclustdir filesep '*KlustTempData_' num2str(thischan) '*']);
                
                fprintf('Channel %u is being trimmed!\n',thischan);
                [block, p_out, trimvec] = alignandcluster_HPC(sweepsCLIPPED,...
                    thischan,1,houredges,subidx,idx_tmp2,times,...
                    doInterpolation,interpSamples,...
                    Mdl,newminrng,newswprng,nclust_trim{cc});
                
                % exit the loop if nothing else needs to be trimmed
                if sum(trimvec) == 0
                    tflag = 4;
                else
                    trimmers{cc} = find(trimvec == 1);
                end
                
                %             idx_tmp2 = [];
                %             idx_tmp3 = [];
            end
            disp('Done!');
            toc
            
            
            p_save{cc} = p_out;
            
            newblock = block;
            newblock.trimstat = trimstat;
            s1 = toc(s0);
            disp(['It took ' num2str(s1) ' seconds to cluster and sort channel ' num2str(thischan) '.'])
        else
            fprintf('No spikes for this channel. Saved arrays will be empty.\n');
            block = [];
            newblock = [];
            p_save{cc} = [];
        end
        % store data
        CELL_BACKUP(cc).channel = thischan;
        CELL_BACKUP(cc).block = newblock;
        CELL_BACKUP(cc).psave = p_save{cc};
        
        save_chdata(newblock,thischan,saveDir,autoclustdir);
        
    end
    
    % save the p variable
    save([saveDir filesep autoclustdir filesep animal '_sorting_p.mat'],'p_save','-v7.3');
    
end



% Save stuff
disp('Saving ALL your data...');
saveallt0 = tic;
save([saveDir filesep autoclustdir filesep animal '_backupCELL.mat'],'CELL_BACKUP','-v7.3');
saveallt1=toc(saveallt0);
fprintf('That took %.2f seconds.\n',saveallt1);


% If clustering over whole experiment, format CELL file so that instead of
% being organized by channel it is a continuous list of all units
% extracted, excluding clusters of quality 4 (that correspond to noise). If
% clustering over whole experiment this will return a single CELL structure
% containing a list of all cells (excluding 4s) identified in the
% clustering process. If clustering in multiple chunks, it will return a
% structure called CHUNKS that will contain a list of cells (i.e. a CELL
% struct) for each chunk.
if nChunks == 0
    CELL = formatAutoCELL_HPC(CELL_BACKUP,nChunks,0);
    
    % Save
    disp('Formatting and saving your data...');
    tic;
    save([saveDir filesep autoclustdir filesep animal '_MASTER_CELL.mat'],'CELL','-v7.3');
    toc;
elseif nChunks > 0
    CHUNKS_init = formatAutoCELL_HPC(CELL_BACKUP,nChunks,0);
    
    CHUNKS = CHUNKS_addFields(CHUNKS_init,animal,dephem,bindatdir);
    % Save
    disp('Formatting and saving your data...');
    tic;
    save([saveDir filesep autoclustdir filesep animal '_MASTER_CHUNKS.mat'],'CHUNKS','-v7.3');
    toc;
end


% delete(poolobj);


% Move all temporary KlustaKwik files into an archive folder
klustarchive = [saveDir filesep autoclustdir filesep 'KlustOutput'];
mkdir(klustarchive);
movefile([saveDir filesep autoclustdir filesep 'KlustTempData*'],klustarchive);

% Copy experiment start file into the CELL output folder
copyfile([bindatdir filesep animal '_expt_START.mat'],[saveDir filesep autoclustdir]);



