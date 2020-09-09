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
chanstart = 1;

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
nSamples    = input('Number of waveform samples?  '); % this is hard-coded at the moment. Could also get it from spike files?
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
    minrng = 20:50;
    swprng = -30:60;
    newminrng = minrng;
    newswprng = swprng;
    interpSamples = 161;
elseif nSamples == 97
    doInterpolation = 0;
    minrng = 10:30;
    swprng = -18:40;
    newminrng = minrng;
    newswprng = swprng;
    interpSamples = 97;
end

fprintf('   *** Number of waveform samples is set to %u. If this is incorrect STOP NOW! ***\n\n',nSamples);
pause(2);
fprintf(' Waveform splicing parameters:\nmin_rng = %u:%u; swprng = -%u:%u;\n\n',minrng(1),minrng(end),abs(swprng(1)),swprng(end));
fprintf(' NEW Waveform splicing parameters:\nnewmin_rng = %u:%u; newswp_rng = -%u:%u;\n\n',newminrng(1),newminrng(end),abs(newswprng(1)),newswprng(end));

nclust = 9; % max number of clusters for KKwik clustering

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
if ~isempty(timefile)
    tinfo_load = load(timefile.name);
    tINFO = tinfo_load.tINFO;
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

startup;
if strcmp(myname,'marmoset')
    mdl_load = load('C:\Users\khengen.USERS\Google_Drive\Matlab_scripts_11_01_2012\Cell_Quality\RANDOMFOREST_TRAINED\rForest_Mdl_KH72_73_75_SHANK01_AT14_trained.mat');
elseif strcmp(myname, 'coypu')
    mdl_load = load('C:\Users\khengen\Google_Drive\Matlab_scripts_11_01_2012\Cell_Quality\RANDOMFOREST_TRAINED\rForest_Mdl_KH72_73_75_SHANK01_AT14_trained.mat');
end

Mdl = mdl_load.Mdl;


% initialize temporary variables for parfor loop
cellarray_size = size(chans,2) - (chanstart-1);
nclust_trim = cell(1,cellarray_size);
trimmers = cell(1,cellarray_size);
idx_tmp = cell(1,cellarray_size);
p_save      = cell(1,cellarray_size);

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
        
        subidx      = cell(nChunks,1);
        subsweeps   = cell(nChunks,1);
        clip_combine = cell(nChunks,1);
		
        for rr = 1:nChunks
            subidx{rr}      = find( scrubbedtimes>=houredges(rr) & scrubbedtimes < houredges(rr+1) );
            subsweeps{rr}   = sweeps(subidx{rr},:); % This is WF segments from remaining spikes
        end
        clear sweeps
        
        % Parallel loop through chunks
        s0 = tic;
        % this below should be a PARFOR loop!
        parfor yy = 1:nChunks
            
            fprintf('Block %u out of %u.\n',yy,size(houredges,2)-1);
            yyt0=tic;
            % perform clustering
            %             if size(subsweeps{yy},1) > 0
            if ~isempty(subsweeps{yy})
                % Use the supervised learning model for autosorting:
                
                fprintf('UsingKlustaKwik to cluster data.\n');
                tic
                
                p_out = []; trimvec = [];
                [block{yy}, p_out, trimvec] = alignandcluster_hpc(subsweeps{yy},...
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
                    
                    % build an idx list
                    idx_tmp_par = [];
                    for jj = 1:size(block{yy}.clust,2);
                        idx_tmp_par = [idx_tmp_par; block{yy}.clust(jj).idx];
                    end
                    idx_tmp{yy} = sortrows(idx_tmp_par);
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
                    
                    [~,idx_tmp2,~] = intersect(scrubbedidx,idx_tmp{yy});
                    [~,clip2,~] = intersect(idx_tmp2,clip);
                    clip_combine{yy} = [clip_combine{yy}; clip2];
                    idx_tmp2(clip_combine{yy}) = [];
                    [~,idx_tmp3,~] = intersect(subidx{yy},idx_tmp2);
                    
                    try
                        sweepsCLIPPED           = subsweeps{yy}(idx_tmp3,:);
                        subidxCLIPPED           = subidx{yy}(idx_tmp3,:);
                    catch
                        keyboard;
                    end
                    
                    
                    % this should NOT repopulate on each cycle through the while loop.
                    % This way, the indices always refer back to the
                    % ORIGINAL spike numbers.
                    
                    % Delete clustering files associated with this channel,
                    % otherwise it may kick back a message that this
                    % channel is already done.
                    delete([saveDir filesep autoclustdir filesep '*KlustTempData_' num2str(thischan) '*']);
                    
                    fprintf('Channel %u is being trimmed!\n',thischan);
                    [block{yy}, p_out, trimvec] = alignandcluster_hpc(sweepsCLIPPED,...
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
        
        
        fprintf('UsingKlustaKwik to cluster data. CONTINUOUS CASE.\n');
        tic
        
        p_out = []; trimvec = [];
        [block, p_out, trimvec] = alignandcluster_hpc(sweeps,thischan,1,...
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
            
            % build an idx list
            idx_tmp_par = [];
            for jj = 1:size(block.clust,2);
                idx_tmp_par = [idx_tmp_par; block.clust(jj).idx];
            end
            idx_tmp{cc} = sortrows(idx_tmp_par);
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
            
            
            [~,idx_tmp2,~] = intersect(scrubbedidx,idx_tmp{cc});
            [~,clip2,~] = intersect(idx_tmp2,clip);
            clip_combine{cc} = [clip_combine{cc}; clip2];
            idx_tmp2(clip_combine{cc}) = [];
			sweepsCLIPPED = [];
            sweepsCLIPPED = sweeps(idx_tmp2,:);
            
            % this should NOT repopulate on each cycle through the while loop.
            % This way, the indices always refer back to the
            % ORIGINAL spike numbers.
            
            % Delete clustering files associated with this channel,
            % otherwise it may kick back a message that this
            % channel is already done.
            delete([saveDir filesep autoclustdir filesep '*KlustTempData_' num2str(thischan) '*']);
            
            fprintf('Channel %u is being trimmed!\n',thischan);
            [block, p_out, trimvec] = alignandcluster_hpc(sweepsCLIPPED,...
                thischan,1,houredges,subidx,idx_tmp{cc},times,...
                doInterpolation,interpSamples,...
                Mdl,newminrng,newswprng,nclust_trim{cc});
            
            % exit the loop if nothing else needs to be trimmed
            if sum(trimvec) == 0
                tflag = 4;
            else
                trimmers{cc} = find(trimvec == 1);
            end
            
        end
        disp('Done!');
        toc
        
        
        p_save{cc} = p_out;
        
        newblock = block;
        newblock.trimstat = trimstat;
        s1 = toc(s0);
        disp(['It took ' num2str(s1) ' seconds to cluster and sort channel ' num2str(thischan) '.'])
        
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

% Move all temporary KlustaKwik files into an archive folder
klustarchive = [saveDir filesep autoclustdir filesep 'KlustOutput'];
mkdir(klustarchive);
movefile([saveDir filesep autoclustdir filesep 'KlustTempData*'],klustarchive);

% Copy experiment start file into the CELL output folder
copyfile([bindatdir filesep animal '_expt_START.mat'],[saveDir filesep autoclustdir]);


