%% SETUP

datadir = '/data/netapp/atorrpac/SHANK_DATA/';
datafile = [datadir filesep 'CONT24_SHANK_v7.mat'];
fprintf('Loading data from: %s\n\n',datafile);
dataload = load(datafile);

% setup other directories
bdir = '/data/netapp/atorrpac/';
addpath(bdir,'-end');

code_dir = '/home/atorrpac/Utilities/hpc_files';
addpath(genpath(code_dir),'-end');

CONT24_SHANK = dataload.CONT24_SHANK;


%clearvars -except CONT24_SHANK
verbose = 0;

nchunks = size(CONT24_SHANK,2);
fprintf('Analyzing %u chunks.\n',nchunks);

%% ANALYSIS PARAMS
q_thresh = 2;
score_thresh = 0;
dep = 1; % all DEP
use_on_off = 0;
cont_thresh = 0.33; % online for 50% of chunk
normalize = 0;
shankBL_norm = 0;
BL_ON = 0;
dataset = 'SHANK24';

% FR calculation parameters
bin_sz = 5*60; % seconds

% cell sep params
plot_cell_sep = 0;
neg_pos_threshold = 0.39;
tailslope_threshold = 0.005;

% number of cells to compare each RSU with for bootstrapping
n_compare = 20;

MU_count = 0;

exclude_anims = {'none'};

%bdir = '/Volumes/turrigiano-lab/ANIMALDATA/SHANK_DATA/MLS_SHANK/';

%% MAIN LOOP
for ch = 1:nchunks
    fprintf('\n  Chunk %u of %u.\n',ch,nchunks);
	cht0 = tic;
    
    cell        = CONT24_SHANK(ch).CELL;
    
    % get RSUs according to analysis parameters
    [shank_RSU{ch}, shank_RSU_idx{ch}] = get_RSUs_CLUSTER(cell,q_thresh,dep,...
        cont_thresh,normalize,BL_ON,dataset,plot_cell_sep,...
        neg_pos_threshold,tailslope_threshold);
    
    n_RSU = numel(shank_RSU_idx{ch});
    % find mean FR for each cell and average in a mean + sem
    kept_count = 0;
    for ii = 1:n_RSU
        fprintf('Cell %u of %u.\n',ii,n_RSU);
        tcell_0 = tic;
 
        this_cell = shank_RSU{ch}(ii);
        data_animals{ch}{ii} = this_cell.animal;
        
        if ~strcmp(this_cell.animal,exclude_anims)
            
            animdir = [bdir filesep this_cell.animal filesep 'BinarySpikes'];
            spikefile = [animdir filesep this_cell.animal '_Channel_' num2str(this_cell.channel) ...
                '_spikes.bin'];
            
            clear n_spikes n_wf first_chunk last_chunk n_toread_* offset_* bfid
            clear B0_* B1_* B0 B1 meanB0 meanB1 mse norm_mse rmse other_RSU
            n_spikes = size(this_cell.idx,1);
            n_wf = ceil(0.1*n_spikes);
            first_chunk = this_cell.idx(1:n_wf+1);
            last_chunk = this_cell.idx(end-n_wf:end);
            % only read first and last chunk to speed everything up.
            n_toread_0 = first_chunk(end) * 97;
            offset_1 = (last_chunk(1) - 1) * 97;
            n_toread_1 = last_chunk(end) * 97 - offset_1;
            % load binary waveforms
            bfid = fopen(spikefile,'rb');
            t0=tic;
            B0_raw = fread(bfid,[n_toread_0,1],'single');
            rewind = fseek(bfid,0,-1);
            seek_stat = fseek(bfid,offset_1*4,-1);
            B1_raw = fread(bfid,[n_toread_1,1],'single');
            fclose(bfid);
            t1=toc(t0);
            fprintf('Loading spike file took %.2f seconds.\n',t1);
            B0_all = reshape(B0_raw,97,[]);
            B1_all = reshape(B1_raw,97,[]);
            B0 = B0_all(:,first_chunk);
            B1_idx = last_chunk - last_chunk(1) + 1;
            B1 = B1_all(:,B1_idx);
            
            meanB0 = nanmean(B0,2);
            meanB1 = nanmean(B1,2);
            [mse,norm_mse,rmse,~] = get_WF_MSE(meanB0,meanB1,0);
            % randomly choose 10 other cells from this day
            
            other_RSU = randperm(n_RSU,n_compare);
            clear file_info new_cell fi_ix newset newidx new_animdir new_spikefile
            clear new_n_spikes new_last_chunk newfid *C1* *D1* *_mse *_norm_mse *_rmse
            
            
            %% first do this loop to have a list of cells on unique animals and channels
            for uu = 1:length(other_RSU)
                new_cell = shank_RSU{ch}(other_RSU(uu));
                file_info(uu,1) = str2double(new_cell.animal(end-1:end));
                file_info(uu,2) = new_cell.channel;
                if uu > 1
                    [~,fi_ix,~] = unique(file_info,'rows');
                    while numel(fi_ix) < size(file_info,1)
                        fprintf('Duplicate found! Getting another cell.\n');
                        newset = setdiff([1:n_RSU],other_RSU);
                        newidx = datasample(newset,1);
                        other_RSU(uu) = newidx;
                        new_cell = shank_RSU{ch}(newidx);
                        file_info(uu,1) = str2double(new_cell.animal(end-1:end));
                        file_info(uu,2) = new_cell.channel;
                        [~,fi_ix,~] = unique(file_info,'rows');
                    end
                end
            end
            
            %% then do the parallel loop
            parfor jj = 1:length(other_RSU)
                % store info about new cell
                new_cell = shank_RSU{ch}(other_RSU(jj));
                
                new_animdir = [bdir filesep new_cell.animal filesep 'BinarySpikes'];
                new_spikefile = [new_animdir filesep new_cell.animal '_Channel_'...
                    num2str(new_cell.channel) '_spikes.bin'];
                new_n_spikes = size(new_cell.idx,1);
                if n_wf >= new_n_spikes
                    new_last_chunk = new_cell.idx;
                else
                    new_last_chunk = new_cell.idx(end-n_wf:end);
                end
                % only read first and last chunk to speed everything up.
                offset_new = (new_last_chunk(1) - 1) * 97;
                n_toread_new = new_last_chunk(end) * 97 - offset_new;
                % load binary waveforms
                newfid = fopen(new_spikefile,'rb');
                t0n=tic;
                seek_stat = fseek(newfid,offset_new*4,-1);
                C1_raw = fread(newfid,[n_toread_new,1],'single');
                fclose(newfid);
                t1n=toc(t0n);
                fprintf('Loading NEW spike file took %.2f seconds.\n',t1n);
                C1_all = reshape(C1_raw,97,[]);
                C1_idx = new_last_chunk - new_last_chunk(1) + 1;
                C1 = C1_all(:,C1_idx);
                
                nc1 = floor(size(C1,2)/2);
                C1_randsamp = randperm(size(C1,2),nc1);
                B0_randsamp = setdiff(1:size(C1,2),C1_randsamp);
                D1 = [C1(:,C1_randsamp) B0(:,B0_randsamp)];
                
                meanC1 = nanmean(C1,2);
                meanD1 = nanmean(D1,2);
                
                
                [new_mse,new_norm_mse,new_rmse,~] = get_WF_MSE(meanB0,meanC1,0);
                [mix_mse,mix_norm_mse,mix_rmse,~] = get_WF_MSE(meanB0,meanD1,0);
                MSE(ii,jj,ch) = new_mse;
                normMSE(ii,jj,ch) = new_norm_mse;
                RMSE(ii,jj,ch) = new_rmse;
                IDX(ii,jj,ch) = other_RSU(jj);
                randWF(:,jj,ii,ch) = meanC1;
                
                mixMSE(ii,jj,ch) = mix_mse;
                mixnormMSE(ii,jj,ch) = mix_norm_mse;
                mixRMSE(ii,jj,ch) = mix_rmse;
                mixIDX(ii,jj,ch) = other_RSU(jj);
                mixrandWF(:,jj,ii,ch) = meanD1;
            end
            
            BOOTSTRAP.RAND.cellWF{ch}(ii,:) = meanB0;
            BOOTSTRAP.MIX.cellWF{ch}(ii,:) = meanB0;
            
            BOOTSTRAP.RAND.MSE = MSE;
            BOOTSTRAP.RAND.normMSE = normMSE;
            BOOTSTRAP.RAND.RMSE = RMSE;
            BOOTSTRAP.RAND.IDX = IDX;
            BOOTSTRAP.RAND.randWF = randWF;
            
            BOOTSTRAP.MIX.MSE = mixMSE;
            BOOTSTRAP.MIX.normMSE = mixnormMSE;
            BOOTSTRAP.MIX.RMSE = mixRMSE;
            BOOTSTRAP.MIX.IDX = mixIDX;
            BOOTSTRAP.MIX.randWF = mixrandWF;
            
            %fprintf('RAND MSE : %.2f uV.\n',nanmean(BOOTSTRAP.RAND.MSE{ch}));
            %fprintf('MIX MSE : %.2f uV.\n',nanmean(BOOTSTRAP.RAND.MSE{ch}));
            
        end
		tcell_1 = toc(tcell_0);
		fprintf('Cell %u of %u took %.2f seconds.\n',ii,n_RSU,tcell_1);
    end
	cht1 = toc(cht0);
	fprintf('\n\t *** Chunk %u took %.2f seconds. ***\n\n',ch,cht1);
end

tic
save([datadir filesep 'Bootstrap_MSE_2.mat'],'BOOTSTRAP','-v7.3');
toc

save_data = 0;
if save_data
    % data
    shankMSE.allrates   = all_rates_bychunk;
    shankMSE.keptcells  = kept_cells;
    shankMSE.startWFs   = meanB0;
    shankMSE.endWFs     = meanB1;
    shankMSE.meanWFs    = meanWF;
    shankMSE.RMSE       = RMSE;
    
    % analysis parameters
    shankMSE.params.quality_thresh          = q_thresh;
    shankMSE.params.score_thresh            = score_thresh;
    shankMSE.params.add_borderline_MUs      = add_borderline_multiunits;
    shankMSE.params.q3_thresh               = q3_thresh;
    shankMSE.params.deprived                = dep;
    shankMSE.params.use_OnOff_times         = use_on_off;
    shankMSE.params.continuity_thresh       = cont_thresh;
    shankMSE.params.normalizeFR             = normalize;
    shankMSE.params.baseline_ON             = BL_ON;
    shankMSE.params.dataset_name            = dataset;
    shankMSE.params.FR_bin_size             = bin_sz;
    shankMSE.params.neg_pos_threshold       = neg_pos_threshold;
    shankMSE.params.tailslope_threshold     = tailslope_threshold;
    shankMSE.params.exclude_animals         = exclude_anims;
    shankMSE.params.check_MSE               = check_MSE;
    
    % make file and save
    fprintf('\n\n\tSaving your data...\n');
    ts0 = tic;
    savefilename = 'SHANK24_MSE';
    files_present = dir([datadir filesep 'SHANK24_MSE*.mat']);
    nf = numel(files_present);
    new_ID = nf + 1;
    savefile = [datadir filesep savefilename '_' num2str(new_ID) '.mat'];
    save(savefile,'shankMSE','-v7.3');
    ts1 = toc(ts0);
    fprintf('That took %.2f seconds.\n',ts1);
end



