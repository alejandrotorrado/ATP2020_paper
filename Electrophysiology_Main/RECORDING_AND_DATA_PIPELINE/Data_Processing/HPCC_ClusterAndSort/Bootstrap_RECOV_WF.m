%% SETUP

datadir = '/data/netapp/atorrpac/RECOV_DATA';
datafile = [datadir filesep 'CONTCELL_recov_MLS_v13.mat'];
fprintf('Loading data from: %s\n\n',datafile);
dataload = load(datafile);
rfile = [datadir filesep 'recov_analysis_60.mat'];
rload = load(rfile);
recov = rload.recov_analysis;

% setup other directories
bdir = '/data/netapp/atorrpac/';
addpath(bdir,'-end');

code_dir = '/home/atorrpac/Utilities/hpc_files';
addpath(genpath(code_dir),'-end');

CONTCELL_recov = dataload.CONTCELL_recov;


%{
%%%%%%%%%%%%%%%%% THIS IS ONLY IF RUNNING ON ATP'S LAPTOP %%%%%%%%%%%%%%%%%
%% load data
clearvars -except *CONTCELL*
clc



if ~exist('CONTCELL_recov','var')
    cload = load('/Users/atorrado/Desktop/MLS_DATA/CONTCELL_recov_MLS_v13.mat');
    CONTCELL_recov = cload.CONTCELL_recov;
end

loadfile = '/Users/atorrado/Desktop/MLS_DATA/DecData/recov_analysis_60.mat';
% loadfile = '/Users/atorrado/Desktop/MLS_DATA/recov_analysis.mat';
rload = load(loadfile);
% rload_FS = load(loadfile_FS);
recov = rload.recov_analysis;

% local bdir
bdir = '/Volumes/turrigiano-lab/ANIMALDATA/MLS_DATA/Continuous';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%}

% extract data
G_bin = recov.G_bin;
dep_status = {'CONTROL','DEPRIVED'};

rsu_ctrl = recov.CONTROL.RSU_idx;
rsu_dep  = recov.DEPRIVED.RSU_idx;

all_RSU = sort([rsu_ctrl'; rsu_dep']);
n_RSU = numel(all_RSU);

% days to analyze
day_list = 6:10; % 6=MD4, 10=ER4 (start timepoint)

CELL = CONTCELL_recov.MASTER;

max_spike_samples = 2e5;
min_spike_samples = 5e4;


%% Loop through cells

for ii = 1:n_RSU
    this_idx = all_RSU(ii);
    this_cell = CELL(this_idx);
    this_anim = this_cell.animal;
    animdir = [bdir filesep this_cell.animal filesep 'BinarySpikes'];
    
    fprintf('Cell %u of %u (%s).\n',ii,n_RSU,this_anim);
    
    spikefile = [animdir filesep this_cell.animal '_Channel_' num2str(this_cell.channel) ...
        '_spikes.bin'];
    timefile  = [animdir filesep this_cell.animal '_Channel_' num2str(this_cell.channel) ...
        '_times.bin'];
    
    
    %% loop through days
    for dd = 1:numel(day_list)
               
        clear day_start day_end spiketimes total_spikes t0 t1 t00 t11
        clear n_spikes n_wf day_spikes wf_spikes n_toread bfid
        clear WF0* meanWF0
        
        day_start = day_list(dd) * 24 * 3600;
        day_end   = (day_list(dd) + 1) * 24 * 3600;
        
        this_anim = this_cell.animal;
        if strcmp(this_anim,'KH67')
            n_samp = 33;
        else
            n_samp = 97;
        end
        fprintf('Note: n_samp = %u!\n',n_samp);
        
        % deal with on/off times
        if numel(this_cell.onTime) == 1 && numel(this_cell.offTime == 1)
            if day_start < this_cell.onTime(1)
                day_start = this_cell.onTime(1);
            elseif day_end > this_cell.offTime(1)
                day_end = this_cell.offTime(1);
            end
        else % multiple on/off times
            %keyboard;
            if all(day_start > [this_cell.onTime(1) this_cell.offTime(1)])
                uu = 2;
            else
                uu = 1;
            end
            if day_start < this_cell.onTime(uu)
                day_start = this_cell.onTime(uu);
            elseif day_end > this_cell.offTime(uu)
                day_end = this_cell.offTime(uu);
            end
        end
        
        
        total_spikes = size(this_cell.idx,1);
        
        spiketimes = this_cell.time;
        t0 = find(spiketimes > day_start,1,'first');
        t1 = find(spiketimes < day_end,1,'last');
        
        day_spikes = this_cell.idx(t0:t1);
        n_spikes = numel(day_spikes);
        
        if n_spikes > 0
        
        
        if n_spikes < min_spike_samples
            n_wf = ceil(n_spikes/2);
        else
            n_wf = max(min_spike_samples, min(max_spike_samples, ceil(0.2*n_spikes) ) );
        end
        wf_spikes = sort(randsample(day_spikes,n_wf));
        
        % only read spikes from this day to speed things up
        offset = (min(wf_spikes) - 1) * n_samp;
        n_toread = (max(wf_spikes) - min(wf_spikes) + 1) * n_samp;
        
        % load binary waveforms
        bfid = fopen(spikefile,'rb');
        t00=tic;
        % go to section of file containing the relevant spikes
        % offset*4 because the file is written in single-precision (4-bit)
        seek_stat = fseek(bfid,offset*4,'bof'); 
        % read all spikes from wf_spikes(1) to wf_spikes(end)
        WF0_raw = fread(bfid,[n_toread,1],'single');
        % close file
        fclose(bfid);
        t11=toc(t00);
        fprintf('Loading spike file took %.2f seconds.\n',t11);
        % reshape waveforms
        WF0_all = reshape(WF0_raw,n_samp,[]);
        
        % index relevant waveforms and select them
        wf_ix = (wf_spikes - wf_spikes(1)) + 1;
        WF0 = WF0_all(:,wf_ix);
        
        if n_samp == 33
            upsamp = 3;
            upnpts = (n_samp-1)*upsamp;
            oldWF0 = WF0;
            clear WF0
            
            
            x = 1:n_samp;
            xq = 1 : 1/3 : n_samp;
            
            newWF0 = interp1(x,oldWF0,xq,'spline');
            newWF0 = newWF0';
            nx = floor(upnpts/5):floor(upnpts/2);
            nsw = -(upnpts/4): (3/4)*upnpts;
            [~,minidx2]     = min(newWF0(:,nx),[],2); % this is 20:50 if you're interpolating at 5x.
            minidx2         = minidx2 + nx(1) - 1;
            WF0        = alignrows(newWF0, minidx2, nsw, size(newWF0,2)-1);
			WF0 = WF0';
        end
        
        
        % compute mean waveform for that day
        meanWF0 = nanmean(WF0,2);
        
        else
            
            WF0 = nan(97,min_spike_samples);
            meanWF0 = nan(97,1);
        end
        
        % store in arrays
        all_waveforms{ii,dd} = WF0;
        mean_waveforms{ii}(dd,:) = meanWF0;
        
        
    end
end



save_data = 1;
if save_data
    % data
%     recov_bootstrap.allWF    = all_waveforms;
    recov_bootstrap.meanWF   = mean_waveforms;
    recov_bootstrap.rsu_idx  = all_RSU;
    recov_bootstrap.G_bin    = G_bin;
    recov_bootstrap.day_list = day_list;
    
    % make file and save
    fprintf('\n\n\tSaving your data...\n');
    ts0 = tic;
    savefilename = 'recov_bootstrap';
    files_present = dir([datadir filesep savefilename '*.mat']);
    nf = numel(files_present);
    new_ID = nf + 1;
    savefile = [datadir filesep savefilename '_' num2str(new_ID) '.mat'];
    save(savefile,'recov_bootstrap','-v7.3');
    ts1 = toc(ts0);
    fprintf('That took %.2f seconds.\n',ts1);
end



