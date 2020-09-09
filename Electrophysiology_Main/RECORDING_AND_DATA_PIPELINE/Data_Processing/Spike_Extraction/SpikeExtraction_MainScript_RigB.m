% SPIKE EXTRACTION MAIN SCRIPT - RIG B - March 2018, ATP
%
% UPDATE MARCH 2018:
% Updated script to work with new in vivo rig (SSC 1-36B).
%
% *** IMPORTANT ***
% This needs to be in the same folder as makeSpikeFiles and
% spikesFromTDTRaw_Upsample in order to work properly.
%
%
% UPDATE May 2016 (kbh + atp):
% Uses tsq file + timing channel (ch66) to get accurate block timestamps
% and save those in the spike file data (INFOX.startTime).
%
% March 2016 version.
%
%                        ****** IMPORTANT ******
%    This script should use analyzemystuff_v10 and spikesFromTDTRaw_v10
%    These are the latest versions that will do spike extraction in the
%    correct way (as of March 23, 2016).
%                        ***********************
%
% This script uses the backed up raw data from a recording and extracts:
%   a) spike shapes and times:
%         saved as [animal]-[date]-[block]-spikes.mat in ALLSPIKES folder
%   b) block numbers and start times:
%         saved as 2D array [block startTime] in [animal]_StartTimes.mat
%
%
% The script asks for user input for the following variables:
%
%  BACKUP_1: should be the 1st backup drive (e.g. F:\ ) for animal 1. Make
%            sure to select the drive and NOT the animal's folder on that
%            drive.
%            For example, BACKUP_1 = 'F:\' is OK; BACKUP_1 = 'F:\KH50\' is
%            not.
%  BACKUP_1a: as above but choose the 2nd backup drive for animal 1.
%  BACKUP_2: as above but choose the 1st backup drive for animal 2.
%  BACKUP_2a: as above but choose the 2nd backup drive for animal 2.
%  subject1: animal 1 code name (e.g. 'KH50')
%  subject2: animal 2 code name (e.g. 'KH52')
%  animalNum:   number of animals being used (1 or 2)
%
%  The variable 'expdate' needs to be hard-coded. It should be the START
%  DATE of the recording, in the format mmddyy (e.g. 062415 for June 24,
%  2015).
clearIDE
% verbose mode gate

% IMPORTANT !
% Choose starting block
start_here = 1;
last_block = 99;



%  -------------------------------------
verbose = 1;
if verbose, fprintf(['\n  *** THIS IS SPIKE EXTRACTION WITH ONLINE ALIGNMENT, '...
        'INTERPOLATION, AND CREATION OF BINARY SPIKE FILES.   *** \n\n']); end

% Basic input - animal number, names and experiment date
animalNum   = input('How many animals are you running?  ');
expdate     = input('Enter experiment start date (format: mmddyy):  ','s');


lgn_exp     = input('Is this an LGN experiment? (1=yes,0=no):  ');


% get animal names
% also check which channels should be skipped (for use with non-TDT arrays)
if animalNum == 1
    t_chan = 36;
    %t_chan = 66;
    subject1 = input('What is the animal name?','s');
    subject2 = 0;
    skipemgref(1) = input(['Skip EMG and ref for ' subject1 '? (1=yes, 0=no):  ']); % usually we skip EMG (16) and ref (8) channels
    if skipemgref(1) == 0
        % if not skipping them, choose by hemisphere
        hemtoinclude(1) = input('Choose which hemisphere you want to include channels 8 & 16 for: 0 (both), 1 (left) or 2 (right):  ');
        switch hemtoinclude(1)
            case 0  % if including both hemisphere's emg and ref, no skipchannels
                skipchans{1} = [];
            case 1  % if including only LEFT hemisphere's emg and ref, skip channels 24 and 32
                skipchans{1} = [24 32];
            case 2  % if including only RIGHT hemisphere's emg and ref, skip channels 8 and 16
                skipchans{1} = [8 16];
        end
    else % if skipping them, exclude them on both hemispheres
        skipchans{1} = [8 16 24 32];
    end
elseif animalNum == 2
    % -------- FIRST ANIMAL -------- 
    t_chan = 66;
    subject1 = input('\nAnimal 1 name? ','s');
    skipemgref(1) = input(['Skip EMG and ref for ' subject1 '? (1=yes, 0=no):  ']); % usually we skip EMG (16) and ref (8) channels
    if skipemgref(1) == 0
        % if not skipping them, choose by hemisphere
        hemtoinclude(1) = input('Choose which hemisphere you want to include channels 8 & 16 for: 0 (both), 1 (left) or 2 (right):  ');
        switch hemtoinclude(1)
            case 0  % if including both hemisphere's emg and ref, no skipchannels
                skipchans{1} = [];
            case 1  % if including only LEFT hemisphere's emg and ref, skip channels 24 and 32
                skipchans{1} = [24 32];
            case 2  % if including only RIGHT hemisphere's emg and ref, skip channels 8 and 16
                skipchans{1} = [8 16];
        end
    else % if skipping them, exclude them on both hemispheres
        skipchans{1} = [8 16 24 32];
    end
    
    % In some recordings (e.g. mice) we only implant one hemisphere, so
    % only want to extract those channels.
    skipwholehem(1) = input(['Do you want to skip a whole hemisphere for ' subject1 '? (1=yes,0=no):  ']);
    if skipwholehem(1)
        hemtoskip(1) = input('Which hemisphere do you want to skip altogether? (1=LEFT,2=RIGHT):  ');
        switch hemtoskip(1)
            case 1
                skipchans{1} = [skipchans{1} 1:16];
            case 2
                skipchans{1} = [skipchans{1} 17:32];
        end
    end
    
    % -------- SECOND ANIMAL -------- 
    subject2 = input('\nAnimal 2 name? ','s');
    skipemgref(2) = input(['Skip EMG and ref for ' subject2 '? (1=yes, 0=no):  ']); % usually we skip EMG (16) and ref (8) channels
    if skipemgref(2) == 0
        % if not skipping them, choose by hemisphere
        hemtoinclude(2) = input('Choose which hemisphere you want to include channels 8 & 16 for: 0 (both), 1 (left) or 2 (right):  ');
        switch hemtoinclude(2)
            case 0 % if including both hemisphere's emg and ref, no skipchannels
                skipchans{2} = [];
            case 1 % if including only LEFT hemisphere's emg and ref, skip channels 24 and 32
                skipchans{2} = [24 32];
            case 2 % if including only RIGHT hemisphere's emg and ref, skip channels 8 and 16
                skipchans{2} = [8 16];
        end
    else % if skipping them, exclude them on both hemispheres
        skipchans{2} = [8 16 24 32];
    end
    
    skipwholehem(2) = input(['Do you want to skip a whole hemisphere for ' subject2 '? (1=yes,0=no):  ']);
    if skipwholehem(2)
        hemtoskip(2) = input('Which hemisphere do you want to skip altogether? (1=LEFT,2=RIGHT):  ');
        switch hemtoskip(2)
            case 1
                skipchans{2} = [skipchans{2} 1:16];
            case 2
                skipchans{2} = [skipchans{2} 17:32];
        end
    end
end

tot_chans = input('How many ephys channels is the rig setup to record?  ');



% Choosing relevant directories
DATADRIVE    = uigetdir('D:\DATA\','Where do you want to copy the data?'); %Local data storage
COYPU_DRIVE  = uigetdir('V:\TDT_DATA_TANKS\','Where is the TDT data tank on Coypu that has the tsq files?'); % Rig PC tsq file storage 
BACKUP_1     = uigetdir('C:\',['backup drive 1 for ' subject1 ' (choose drive not folder)']); % backup drive 1 (animal 1)
BACKUP_1a    = uigetdir('C:\',['backup drive 2 for ' subject1 ' (choose drive not folder)']); % backup drive 2 (animal 1)
svdir_1      = [DATADRIVE filesep subject1 filesep 'ALLSPIKES'];
if ~exist(svdir_1,'dir'), mkdir(svdir_1); end

if animalNum == 2
    BACKUP_2     = uigetdir('C:\',['backup drive 1 for ' subject2 ' (choose drive not folder)']); % backup drive 1 (animal 2)
    BACKUP_2a    = uigetdir('C:\',['backup drive 2 for ' subject2 ' (choose drive not folder)']); % backup drive 2 (animal 2)
    svdir_2      = [DATADRIVE filesep subject2 filesep 'ALLSPIKES'];
    if ~exist(svdir_2,'dir'), mkdir(svdir_2); end
else
    svdir_2 = svdir_1;
end


TSQstamp = [];
last_TSQ_dir = '';
gmt_corrected = 0; % gate to tell code whether TSQ timestamp is GMT corrected
%% start here for last block
for ii = start_here:last_block % set start back to 1
    
    curr_block_1 = [BACKUP_1 filesep subject1 filesep 'Block-' num2str(ii)];
    next_block_1 = [BACKUP_1 filesep subject1 filesep 'Block-' num2str(ii+1)];
    if animalNum == 2
        curr_block_2 = [BACKUP_2 filesep subject2 filesep 'Block-' num2str(ii)];
        next_block_2 = [BACKUP_2 filesep subject2 filesep 'Block-' num2str(ii+1)];
    end
    
    %%% ----------------------
    %%% HOPEFULLY THIS WILL WORK! 
    %%% NOT SURE IF CAN READ TSQ FILES WHILE THEY'RE BEING WRITTEN...
    %
    % Check for the TSQ file - it should only exist in the first block of a
    % recording (and first after a crash). Use this to find the expt. start
    % time:
    tsqtemp = [];
    fprintf('\nLooking for TSQ timestamp...\n');
    % figure out which directory holds the correct TSQ file
    all_files = dir([curr_block_1 filesep '*.sev']);
    while isempty(all_files);
        pause(10);
        fprintf('Still waiting for files to be copied...\n');
        all_files = dir([curr_block_1 filesep '*.sev']);
    end
    fprintf('  Found files!\n\n');
    one_file = all_files(1).name;
    split_1 = regexp(one_file,'Chronic-','split');
    split_1 = split_1{end};
    split_2 = regexp(split_1,'_RAWX','split');
    split_2 = split_2{1};
    
    if isempty(regexp(last_TSQ_dir,split_2,'ONCE'))
        all_tsqdirs = dir([COYPU_DRIVE filesep 'Chronic-*']);
        cell_match = regexp({all_tsqdirs(:).name},split_2);
        newtsqdir_idx = find(~cellfun('isempty',cell_match));
        curr_TSQ_dir = [COYPU_DRIVE filesep all_tsqdirs(newtsqdir_idx).name];
        last_TSQ_dir = curr_TSQ_dir;
        tsqtemp = checkTSQ(curr_TSQ_dir);

        if isempty(TSQstamp) && not (tsqtemp) && ii > 1
            % if no previous stamp is found, no new one is found, and we are
            % past block 1, we have a problem.
            fprintf('Problem! We''re past the first block and there is no TSQ timestamp.\n');
            keyboard
        else
            fprintf('Found TSQ timestamp! This folder started on: %s.\n',datestr(unixtime(tsqtemp)));
            TSQstamp = unixtime(datevec(TimezoneConvert(datestr(unixtime(tsqtemp)),'UTC','America/New_York')));
            fprintf('GMT corrected folder start time: %s.\n\n',datestr(unixtime(TSQstamp)));
        end
        
        
    else
        % keep the current TSQstamp
        if verbose, disp('No new TSQ timestamp.'); end
        curr_TSQ_dir = last_TSQ_dir;
    end
    %%% -------------------
    
    % wait for block to finish being recorded and copied to back up
    m1 = 0;
    while m1 == 0
        if ii < last_block
            m1 = exist(next_block_1,'dir');
            fprintf('Waiting for block %u of animal 1 to become available.\n',ii);
            pause(2);
        elseif ii == last_block
            fprintf('Block %u is the last block in the recording. Assuming it is available for animal 1...\n',ii);
            m1 = 1;
        end
    end
    % wait for block to finish being recorded and copied to back up
    if animalNum == 2
        fprintf('  Animal 1 files are available!\n');
        
        m2 = 0;
        while m2 == 0
            if ii < last_block
                m2 = exist(next_block_2,'dir');
                fprintf('Waiting for block %u of animal 2 to become available.\n',ii);
                pause(2);
            elseif ii == last_block
                fprintf('Block %u is the last block in the recording. Assuming it is available for animal 2...\n',ii);
                m2 = 1;
            end
        end
        
        fprintf('  Animal 2 files are available!\nStarting spike extraction for block %u.\n\n',ii);
    else
        fprintf('  Animal 1 files are available!\nStarting spike extraction for block %u.\n\n',ii);
    end
   
    
    thisrow = [];
    % make this a parloop when done debugging
    for ee=1:animalNum
        if ee==1
            % only extract spikes if not already done
            if ~exist([DATADRIVE filesep subject1 filesep 'ALLSPIKES' filesep (subject1)...
                    '-' (expdate) '-' num2str(ii) '-' 'spikes.mat'],'file')
                disp('here1');
                % This function creates both the .mat spike files and the
                % binary spike files
                if tot_chans == 64
                    makeSpikeFiles(curr_block_1,subject1,expdate,ii,DATADRIVE,1,TSQstamp,skipchans{1},t_chan,lgn_exp);
                elseif tot_chans == 16
                    makeSpikeFiles_16ch(curr_block_1,subject1,expdate,ii,DATADRIVE,1,TSQstamp,skipchans{1},t_chan,lgn_exp);
                elseif tot_chans == 32
                    try
                        makeSpikeFiles_32ch(curr_block_1,subject1,expdate,ii,DATADRIVE,1,TSQstamp,skipchans{1},t_chan,lgn_exp);
                    catch
                        keyboard
                    end
                else
                    fprintf('You told me rig is setup to record %u channels. That''s weird.\n',tot_chans);
                    keyboard;
                end
            else
                fprintf('Spikes for block %u, animal %s already extracted! Moving on.\n',ii,subject1);
            end
        elseif ee==2
            % only extract spikes if not already done
            if ~exist([DATADRIVE filesep subject2 filesep 'ALLSPIKES' filesep (subject2)...
                    '-' (expdate) '-' num2str(ii) '-' 'spikes.mat'],'file')
                
                % This function creates both the .mat spike files and the
                % binary spike files
                if tot_chans == 64
                    makeSpikeFiles(curr_block_2,subject2,expdate,ii,DATADRIVE,2,TSQstamp,skipchans{2},t_chan,lgn_exp);
                elseif tot_chans == 16
                    makeSpikeFiles_16ch(curr_block_2,subject2,expdate,ii,DATADRIVE,2,TSQstamp,skipchans{2},t_chan,lgn_exp);
                else
                    fprintf('You told me rig is setup to record %u channels. That''s weird.\n',tot_chans);
                    keyboard;
                end
               
            else
                fprintf('Spikes for block %u, animal %s already extracted! Moving on.\n',ii,subject2);
            end
        end
    end
    
    if ii < last_block
        fprintf('Done with block %u. Moving on.\n',ii);
    elseif ii == last_block
        fprintf('Done with block %u. That was the last block in the recording!\n\n',ii);
    end
    
end

% 08/27/18
%% RUN FROM HERE TO BACKUP BINARY SPIKE FILES

% Point to backup drives
backup_spikes_1     = [BACKUP_1 filesep subject1 filesep 'ALLSPIKES'];
backup_spikes_1a    = [BACKUP_1a filesep subject1 filesep 'ALLSPIKES'];
backup_bin_1        = [BACKUP_1 filesep subject1 filesep 'ALLSPIKES' filesep 'BinarySpikes'];
backup_bin_1a       = [BACKUP_1a filesep subject1 filesep 'ALLSPIKES' filesep 'BinarySpikes'];

if animalNum == 2
    backup_spikes_2     = [BACKUP_2 filesep subject2 filesep 'ALLSPIKES'];
    backup_spikes_2a    = [BACKUP_2a filesep subject2 filesep 'ALLSPIKES'];
    backup_bin_2        = [BACKUP_2 filesep subject2 filesep 'ALLSPIKES' filesep 'BinarySpikes'];
    backup_bin_2a       = [BACKUP_2a filesep subject2 filesep 'ALLSPIKES' filesep 'BinarySpikes'];
end

% Make Experiment Start file
subjs = {subject1 subject2};
svds = {svdir_1 svdir_2};
chanID = {'A' 'B'};
for xx=1:2
    try
        animal = subjs{xx};
        svdir = svds{xx};
        chanidx = chanID{xx};
        found_file = 0; file_counter = 0;
        while found_file == 0
            file_counter = file_counter + 1;
            first_spikefile = dir([svdir filesep '*-' num2str(file_counter) '-spikes.mat']);
            if ~isempty(first_spikefile)
                found_file = 1;
            end
        end
        infox = ['INFO' chanidx];
        first_loaded = load([svdir filesep first_spikefile.name],infox);
        tINFO = first_loaded.(infox).startTime;
        save([svdir filesep 'BinarySpikes' filesep animal '_expt_START.mat'],'tINFO');
    catch
        fprintf('\n\t*** There was an error creating expt_START file for animal #%u ***\n\n',xx);
    end
end


% Backup all spike files


fprintf('\n\nBacking up the spikes file to backup drives.\nThis may take a while...\n');

try
t11 = tic;
copyfile([svdir_1 filesep '*spikes.mat'],backup_spikes_1);
copyfile([svdir_1 filesep 'BinarySpikes' filesep '*.bin'],backup_bin_1);
copyfile([svdir_1 filesep 'BinarySpikes' filesep '*.mat'],backup_bin_1);
t11e = toc(t11);
fprintf('Done backing up spike files for animal 1 (%s) on backup drive 1.\n',subject1);
fprintf('That took %.1f seconds.\n\n',t11e);

t12 = tic;
copyfile([svdir_1 filesep '*spikes.mat'],backup_spikes_1a);
copyfile([svdir_1 filesep 'BinarySpikes' filesep '*.bin'],backup_bin_1a);
copyfile([svdir_1 filesep 'BinarySpikes' filesep '*.mat'],backup_bin_1a);

t12e = toc(t12);
fprintf('Done backing up spike files for animal 1 (%s) on backup drive 2.\n',subject1);
fprintf('That took %.1f seconds.\n\n',t12e);
catch
   disp('skippppp'); 
end
t21 = tic;
copyfile([svdir_2 filesep '*spikes.mat'],backup_spikes_2);
copyfile([svdir_2 filesep 'BinarySpikes' filesep '*.bin'],backup_bin_2);
copyfile([svdir_1 filesep 'BinarySpikes' filesep '*.mat'],backup_bin_2);

t21e = toc(t21);
fprintf('Done backing up spike files for animal 2 (%s) on backup drive 1.\n',subject2);
fprintf('That took %.1f seconds.\n\n',t21e);

t22 = tic;
copyfile([svdir_2 filesep '*spikes.mat'],backup_spikes_2a);
copyfile([svdir_2 filesep 'BinarySpikes' filesep '*.bin'],backup_bin_2a);
copyfile([svdir_1 filesep 'BinarySpikes' filesep '*.mat'],backup_bin_2a);

t22e = toc(t22);
fprintf('Done backing up spike files for animal 2 (%s) on backup drive 2.\n',subject2);
fprintf('That took %.1f seconds.\n\n',t22e);

fprintf('Done backing up spike files!\n\n');

fprintf('Make sure to check block start times in the spikes files. Trust that the GMT correction is OK.\n\n');



