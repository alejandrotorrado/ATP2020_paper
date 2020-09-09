% COPY AND BACKUP MAIN SCRIPT - RIG B
%
% UPDATE MARCH 2018:
% Updated script to work with new in vivo rig (SSC 1-36B).
%
%
% ROBOCOPY syntax, for reference:
% robocopy source destination filename > NUL
% __NOTES:
%   __source and destination can be put within double quotes if they contain
%     spaces
%   __filename can contain wildcard character (*)
%   __ > NUL makes it so no output is displayed
%   __example: cmd = 'robocopy C:\ E:\Documents myfile*.txt > NUL'; % make cmd string
%              system(cmd); % matlab system call
%
%
% As of March, 2015, TDT data files are written in a .sev format. The
% original files are also written on the PC running data collection. It is
% faster and cleaner to operate with the SEV files, so this script
% coordinates the backup of the SEV files as well as spike extraction from
% each of the 1h blocks of data (as before). Downstream of this code, the
% remainining analyses should be compatible with the spikes files created
% here.
%
% NOTE: The SEV files do NOT have a timestamp for the start of the block.
% This has to be read from a .tsq file that is part of the original data
% tank output (stored on the recording PC, as mentioned above). The scripts
% coordinated here will read the TSQ file for the timestamp (in unix time)
% as well as the SEV files that are saved on the RS4.
%
%
% _________________________________________________________________________
% ************************** IMPORTANT NOTE *******************************
% **   If using RS4 workaround (i.e. saving to local) the variables      **
% **   "direc" and "tsqdir" should point to the same location - namely,  **
% **   the DataTank for the current experiment on the Rat Rig PC         **
% **   (rattisriggus.bio.brandeis.edu)                                   **
% *************************************************************************
%
%
% Based on coordinate_copy_and_extraction_v4.m (Written by KBH)
% Based on coordinate_copy_and_extraction_v8.m
%
%
%


% verbose mode gate
verbose = 1;
% backup = 0 does not copy to ext drives
% deletestuff = 0 does not delete from RS4
% both of these originally introduced for debugging purposes
backup      = 1;
deletestuff = 0;
both_hems = 0;
 both_hems = 1; % in most cases this will be =1. Only change to 0 if you want to restrict backup to channels on one hemisphere.
if ~both_hems
    one_hem_only = 1;
    fprintf(['Modify the "chanlist" variable on line 72 to choose subset of channels to backup (e.g. 1:16 or [2 4 6 8]).\n'...
        'Hit continue when you are done. Note that the channel selection will apply to both animals.\n']);
    keyboard
    chanlist = 1:16; % choose which subset of channels to backup
end



% Basic input - animal number, names and experiment date
n_animals   = input('How many animals are you running?  ');
expdate     = input('Experiment start date? (mmddyy):  ','s');
tankname    = input('Recording name (e.g. KH65_66)?  ','s');

expertype   = input('Is this experiment Chronic(1) or Crickets_BJL(0)?  '); %This line and following if statement added June 2020 to allow different rig configurations (BJL)
if expertype == 1
expname     = 'Chronic';
elseif expertype == 0
expname     = 'Crickets_BJL';
end

if n_animals == 1
%     subject1 = 'JB41';%input('What is the animal number?','s');
    subject1 = 'BL51';%input('What is the animal number?','s');
    subject2 = 0;
elseif n_animals == 2
    tanksplit = regexp(tankname,'\d','split');
    basename = tanksplit{1};
    numstartlist = regexp(tankname,'\d*');
    numstart = numstartlist(1);
    animsplt = regexp(tankname(numstart:end),'_','split');
    subject1 = [basename animsplt{1}]; % alternative: input('Animal 1 name? ','s');
    subject2 = [basename animsplt{2}]; % alternative: input('Animal 2 name? ','s');
end

fprintf('Starting copy and backup for animals %s and %s.\n',subject1,subject2);

% Choosing relevant directories
[direc]     = uigetdir('W:\','Pick your RS4 data tank OR rig pc data tank'); % RS4 data tank (OR RAT RIG DATA TANK - SEE IMPORTANT NOTE AT TOP)
[tsqdir]    = uigetdir('V:\TDT_DATA_TANKS','Select the data tank on the rig computer'); % Rig computer data tank

SDirec1     = uigetdir('D:\',['backup drive 1 for ' subject1]); % backup drive 1 (animal 1)
SDirec2     = uigetdir('D:\',['backup drive 2 for ' subject1]); % backup drive 2 (animal 1)
if n_animals == 2
    SDirec3     = uigetdir('D:\',['backup drive 1 for ' subject2]); % backup drive 3 (animal 2)
    SDirec4     = uigetdir('D:\',['backup drive 2 for ' subject2]); % backup drive 4 (animal 2)
end

% Make save directories
if exist('subject1','var');
    if verbose ==1, disp('make folders for S1'), end;
    if exist([SDirec1 filesep subject1],'dir') == 0;
        try
            mkdir( [SDirec1 filesep subject1]);
            mkdir( [SDirec2 filesep subject1]);
        catch
            mkdir( [SDirec1 filesep subject1]);
            mkdir( [SDirec2 filesep subject1]);
        end
        
    end
end

if subject2 ~= 0;
    if verbose ==1, disp('make folders for S2'), end;
    if exist([SDirec3 filesep subject2],'dir') == 0;
        try
            mkdir( [SDirec3 filesep subject2]);
            mkdir( [SDirec4 filesep subject2]);
        catch
            mkdir( [SDirec3 filesep subject2]);
            mkdir( [SDirec4 filesep subject2]);
        end
    end
end

% block folder number in which to copy TSQ files. This is an array that values will be appended to.
tsqblock = [1]; 

% Find how many data folders are in the tank name and figure out their
% names.
[dfolders,mfolders,num_folds] = get_folder_info(direc,expname);

tot_blocks = 0; % keep track of total number of blocks. Initialize as 0
new_tot_blocks = 0; % initialize as 0
curr_dir_idx = 1; % initialize index of folder to copy/extract (set to 1 initially)
lastblock = 150; % last block (don't wait for next if reach this)
% change the index above to start at a different folder than the first one


keep_going = 1; % keep going while this flag is set to "true"
while keep_going
    
    
    % get folder name
    [dfolders,mfolders,num_folds] = get_folder_info(direc,expname);
    curr_dir = mfolders(curr_dir_idx,:);
    n_files_start = size(dir([direc filesep curr_dir]),1);
    
    % check that folder exists (or raise exception if it doesn't)
    if exist([(direc) filesep curr_dir],'dir')~=0
        
        if verbose==1, disp(['Directory ' (direc) filesep curr_dir ' exists, processing...']), end;
        
        % find number of hour-long blocks in this folder
        nblocks = get_nblocks(direc,curr_dir);
        
        ii = 0; % initialize block counter (block you want to start at minus 1, e.g. set to 24 if you want to start copying from block 25)
        while ii < nblocks % should avoid copying last block (might be in acquisition)
            
            ii = ii + 1; % increase block counter
            
            %             if verbose && exist([SDirec1 filesep subject1 filesep 'Block-' num2str(tot_blocks+ii)],'dir')
            %                 fprintf('WARNING! This block directory already exists.\nCheck that you are not overwriting data.');
            %                 send_text_message('7815026988','att','COPY','Block already exists.');
            %                 keyboard;
            %             end
            
            fprintf('Directory: %s.\nNow starting block %u out of %u.\n', curr_dir, ii, nblocks);
            
            if ii == nblocks % if last available block is next
                disp('Reached last block in current directory.');
                % check if there is another folder
                [dfolders,mfolders,num_folds] = get_folder_info(direc,expname);
                if (tot_blocks + ii) >= lastblock
                    disp('This is the last block in the whole recording!');
                    keep_going = 0;
                else
                    n_files_start = size(dir([direc filesep curr_dir]),1);
                    if curr_dir_idx < num_folds
                        disp('Found more folders.');
                        curr_dir_idx = curr_dir_idx + 1;
                        new_tot_blocks = tot_blocks + ii; % update total block counter
                        fprintf('After this block is done, will move on to folder %s.\nTotal block count is now %u.\n',mfolders(curr_dir_idx,:),new_tot_blocks);
                        tsqblock = [tsqblock (new_tot_blocks + 1)];
                    else
                        % wait for next block
                        while size(dir([direc filesep curr_dir]),1) <= n_files_start
                            disp('Waiting for next block...')
                            pause(30);
                        end
%                         nblocks = nblocks + 1; % increase max number of blocks in this folder
                        n_files_start = size(dir([direc filesep curr_dir]),1);
                        fprintf('Found new block! Moving on.');
                    end
                end
            end
            
            % Copy files from RS4 and rig pc to external drives if backup = 1.
            % The code under the if statement also deletes files from the RS4
            % to make space. In general, should always have backup = 1 when
            % running experiments, and only set to 0 when debugging other parts
            % of the code (to avoid deleting files from the RS4 while we still
            % need them).
            
            
            %             Locations of backup folders on backup drives:
            location_A1 = [(SDirec1) filesep subject1 filesep 'Block-' num2str(tot_blocks + ii)]; % backup drive 1, subject 1
            location_A2 = [(SDirec2) filesep subject1 filesep 'Block-' num2str(tot_blocks + ii)]; % backup drive 2, subject 1
            if n_animals == 2
                location_B1 = [(SDirec3) filesep subject2 filesep 'Block-' num2str(tot_blocks + ii)]; % backup drive 3, subject 2
                location_B2 = [(SDirec4) filesep subject2 filesep 'Block-' num2str(tot_blocks + ii)]; % backup drive 4, subject 2
            end
            
            
            % Backup for 1 animal: copy all files in RS4 datatank to backup drives
            %
            % Backup for 2 animals: copy all files with "RAWX_ch1-32" in filename to
            %   'subject1' locations on the backup drives; copy all files with
            %   "RAWX_ch33-64" in filename to 'subject2' locations on backup drives.
            %
            % Also backup in each Block folder the corresponding '.tsq' file from
            %   the rig computer (rattisriggus). The tsq file for Block x will have
            %   timestamps valid for both animals for that particular Block.
            %
            %
            % This loop does it: change back to parfor if you mess with it
            
            
            % "block" is block hour index in TDT convention
            if ii > 1
                block = ['-' num2str(ii - 1) 'h'];
            else
                block = '';
            end
            
            if n_animals == 1
                
                if verbose == 1, disp(['Backing up data (1 animal) to backup drive 1']), end;
                % Animal 1
                %   Backing up raw data
                tic
                if both_hems
                    channels = [1:32 66];
                elseif one_hem_only
                    channels = [chanlist 66];
                end
                for cc = channels
                    % copy to backup drive 1 first (over network)
                    try
                        copy_cmd = ['robocopy ' direc filesep curr_dir ' ' location_A1 ' ' '*RAWX_ch' num2str(cc) block '.sev '...
                            ' > NUL'];
                        copystatus = system(copy_cmd);
                    catch
                        copy_cmd = ['robocopy ' direc filesep curr_dir ' ' location_A1 ' ' '*RAWX_ch' num2str(cc) block '.sev '...
                            ' > NUL'];
                        copystatus = system(copy_cmd);
                    end
                    
                    
                    % copy to drive 2 next (over USB 3.0)
                    try
                        copy_cmd = ['robocopy ' location_A1 ' ' location_A2 ' ' '*RAWX_ch' num2str(cc) block '.sev '...
                            ' > NUL'];
                        copystatus = system(copy_cmd);
                    catch
                        copy_cmd = ['robocopy ' location_A1 ' ' location_A2 ' ' '*RAWX_ch' num2str(cc) block '.sev '...
                            ' > NUL'];
                        copystatus = system(copy_cmd);
                    end
                end
                if verbose == 1, disp(['Animal 1/1: Copy to both backup drives done.']), end
                toc
                
                
            elseif n_animals == 2  % 2 animals
                
                % ___ANIMAL 1
                
                if verbose == 1, disp(['Backing up data  for animal 1.']), end;
                
                % Animal 1
                %   Backing up raw data
                % channels 1:32 are ephys data from animal 1; ch 66 is
                % experiment timer
                tic
                if both_hems
                    channels = [1:32 66];
                elseif one_hem_only
                    channels = [chanlist 66];
                end
                for cc = channels
                    % copy to backup drive 1 first (over network)
                    try
                        copy_cmd = ['robocopy ' direc filesep curr_dir ' ' location_A1 ' ' '*RAWX_ch' num2str(cc) block '.sev '...
                            ' > NUL'];
                        copystatus = system(copy_cmd);
                    catch
                        copy_cmd = ['robocopy ' direc filesep curr_dir ' ' location_A1 ' ' '*RAWX_ch' num2str(cc) block '.sev '...
                            ' > NUL'];
                        copystatus = system(copy_cmd);
                    end
                    
                    
                    % copy to drive 2 next (over USB 3.0)
                    try
                        copy_cmd = ['robocopy ' location_A1 ' ' location_A2 ' ' '*RAWX_ch' num2str(cc) block '.sev '...
                            ' > NUL'];
                        copystatus = system(copy_cmd);
                    catch
                        copy_cmd = ['robocopy ' location_A1 ' ' location_A2 ' ' '*RAWX_ch' num2str(cc) block '.sev '...
                            ' > NUL'];
                        copystatus = system(copy_cmd);
                    end
                end
                if verbose == 1, disp(['Animal 1/2: Copy to both backup drives done.']), end;
                toc
                
                
                
                
                % ___ANIMAL 2
                
                if verbose == 1, disp(['Backing up data for animal 2.']), end;
                % Animal 2
                %   Backing up raw data
                % channels 33:64 are ephys data for animal 2; ch 66 is
                % experiment timer
                tic
                if both_hems
                    channels = [33:64 66];
                elseif one_hem_only
                    channels  = [32+chanlist 66];
                end
                for cc = [33:64 66]
                    % copy to drive 1 first (over network)
                    try
                        copy_cmd = ['robocopy ' direc filesep curr_dir ' ' location_B1 ' ' '*RAWX_ch' num2str(cc) block '.sev '...
                            ' > NUL'];
                        copystatus = system(copy_cmd);
                    catch
                        copy_cmd = ['robocopy ' direc filesep curr_dir ' ' location_B1 ' ' '*RAWX_ch' num2str(cc) block '.sev '...
                            ' > NUL'];
                        copystatus = system(copy_cmd);
                    end
                    
                    % copy to drive 2 next (over USB)
                    try
                        copy_cmd = ['robocopy ' location_B1 ' ' location_B2 ' ' '*RAWX_ch' num2str(cc) block '.sev '...
                            ' > NUL'];
                        copystatus = system(copy_cmd);
                    catch
                        copy_cmd = ['robocopy ' location_B1 ' ' location_B2 ' ' '*RAWX_ch' num2str(cc) block '.sev '...
                            ' > NUL'];
                        copystatus = system(copy_cmd);
                    end
                end
                toc
                if verbose == 1, disp(['Animal 2/2: Copy to both backup drives done.']), end;
                
                
            else
                error('Do it right, son. Should be running 1 or 2 animals.');
            end
            
            
            % update number of 1h blocks in folder
            nblocks = get_nblocks(direc,curr_dir);
            
        end % end of while loop that loops through blocks in folder
        tot_blocks = new_tot_blocks;
    else
        error(['Directory ' (direc) filesep curr_dir ' does not exist.']);
    end
end

% Now backup tsq files in correct block folders
%%
%   Backing up '.tsq' file
fprintf('\n\n\n\t *** DONE BACKING UP RAW EPHYS DATA! ***\n');
fprintf('\n\t *** MOVING ON TO TSQ FILES BACKUP ***\n\n');
tsqblock = unique(tsqblock);
for cc = tsqblock
    
    % Locations of backup folders on backup drives:
    location_A1 = [(SDirec1) filesep subject1 filesep 'Block-' num2str(cc)]; % backup drive 1, subject 1
    location_A2 = [(SDirec2) filesep subject1 filesep 'Block-' num2str(cc)]; % backup drive 2, subject 1
    if n_animals == 2
        location_B1 = [(SDirec3) filesep subject2 filesep 'Block-' num2str(cc)]; % backup drive 3, subject 2
        location_B2 = [(SDirec4) filesep subject2 filesep 'Block-' num2str(cc)]; % backup drive 4, subject 2
    end
    
    all_files = dir([location_A1 filesep '*.sev']);
    one_file = all_files(1).name;
    split_1 = regexp(one_file,'_','split');
    copy_dir = fullfile(tsqdir,split_1{3});
    
    % animal 1
    if verbose
        fprintf('Copying tsq file %u out of %u to animal 1 backup drives.\n',...
            find(tsqblock==cc),numel(tsqblock));
    end
    % copy to drive 1
    tic
    try
        copy_tsq_cmd = ['robocopy ' copy_dir ' ' location_A1 ' ' '*.tsq*'...
            ' > NUL'];
        tsq_status = system(copy_tsq_cmd);
    catch
        pause(2);
        copy_tsq_cmd = ['robocopy ' copy_dir ' ' location_A1 ' ' '*.tsq*'...
            ' > NUL'];
        tsq_status = system(copy_tsq_cmd);
    end
    
    % copy to drive 2
    try
        copy_tsq_cmd = ['robocopy ' copy_dir ' ' location_A2 ' ' '*.tsq*'...
            ' > NUL'];
        tsq_status = system(copy_tsq_cmd);
    catch
        pause(2);
        copy_tsq_cmd = ['robocopy ' copy_dir ' ' location_A2 ' ' '*.tsq*'...
            ' > NUL'];
        tsq_status = system(copy_tsq_cmd);
    end
    
    if verbose, disp('Animal 1/2: copy of tsq files to both backup drives done.'); end
    toc
    
    
    if n_animals == 2
        % animal 2
        if verbose
            fprintf('Copying tsq file %u out of %u to animal 2 backup drives.\n',...
                find(tsqblock==cc),numel(tsqblock));
        end
        % copy to drive 1
        tic
        try
            copy_tsq_cmd = ['robocopy ' copy_dir ' ' location_B1 ' ' '*.tsq*'...
                ' > NUL'];
            tsq_status = system(copy_tsq_cmd);
        catch
            pause(2);
            copy_tsq_cmd = ['robocopy ' copy_dir ' ' location_B1 ' ' '*.tsq*'...
                ' > NUL'];
            tsq_status = system(copy_tsq_cmd);
        end
        
        % copy to drive 2
        try
            copy_tsq_cmd = ['robocopy ' copy_dir ' ' location_B2 ' ' '*.tsq*'...
                ' > NUL'];
            tsq_status = system(copy_tsq_cmd);
        catch
            pause(2);
            copy_tsq_cmd = ['robocopy ' copy_dir ' ' location_B2 ' ' '*.tsq*'...
                ' > NUL'];
            tsq_status = system(copy_tsq_cmd);
        end
        
        if verbose, disp('Animal 2/2: copy of tsq files to both backup drives done'); end
        toc
    end
    
end

disp(['*** DONE COPYING ***']);

% send_text_message('781-502-6988','att','MatLab Status',...
%     ['Finished processing RS4 data tank ' direc(4:end) ' ' today ]);




