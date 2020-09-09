

% createVidFile
%
% This takes the raw video files saved by the MacMini and creates a single
% avi file containing the whole video recording. This will then be used to
% perform motion tracking using blob analysis.
%
% ATP, April 2016

clearIDE;



rawVideoDir = uigetdir(cd,'Pick the folder with all the raw video files.');
saveDir = uigetdir(cd,'Pick the directory where you want to save the processed video and movement data.');
animal = input('Animal(s) name?  ','s');

clearvars -except *Dir animal vidfile framerate

framefiles = dir([rawVideoDir filesep animal '*.csv']);
N_frfiles = size(framefiles,1);
frfiledates(:,1) = datenum(cat(1,framefiles.date));
frfiledates(:,2) = 1:N_frfiles;
frfiledates = sortrows(frfiledates,1);

rawfiles = dir([rawVideoDir filesep animal '*.avi']);
N_files  = size(rawfiles,1);
filedates(:,1) = datenum(cat(1,rawfiles.date));
filedates(:,2) = 1:N_files;
filedates = sortrows(filedates,1);

if N_files ~= N_frfiles
    disp('   ******* PROBLEM: file numbers don''t match!!!!!!!!!!!!!');
    keyboard;
end

% check file quality
fprintf('\nChecking video file quality.\n\n');
bad_files = [];
for ii =1:N_files
    tic
    fprintf('Checking file %u out of %u.\n',ii,N_files);
    vidchunkidx = filedates(ii,2);
    
    vidchunkfile = [rawVideoDir filesep rawfiles(vidchunkidx).name];
    
    try
        vidchunk = VideoReader(vidchunkfile);
    catch
        fprintf('File %s (number %u out of %u) is bad.\n',rawfiles(vidchunkidx).name,ii,N_files);
        bad_files(ii,1) = vidchunkidx;
    end
    toc
end

if ~isempty(bad_files)
    for bb = 1:size(bad_files,1);
        if bad_files(bb,1) ~= 0
            vidchunkidx = bad_files(bb,1);
            vidchunkfile = [rawVideoDir filesep rawfiles(vidchunkidx).name];
            timechunkfile = [rawVideoDir filesep framefiles(vidchunkidx).name];
            movefile(vidchunkfile,[rawVideoDir filesep 'badVids' filesep]);
            movefile(timechunkfile,[rawVideoDir filesep 'badVids' filesep]);
        end
    end
    
    clear filedates frfiledates framefiles rawfiles N_frfiles N_files
    
    framefiles = dir([rawVideoDir filesep animal '*.csv']);
    N_frfiles = size(framefiles,1);
    frfiledates(:,1) = datenum(cat(1,framefiles.date));
    frfiledates(:,2) = 1:N_frfiles;
    frfiledates = sortrows(frfiledates,1);
    
    rawfiles = dir([rawVideoDir filesep animal '*.avi']);
    N_files  = size(rawfiles,1);
    filedates(:,1) = datenum(cat(1,rawfiles.date));
    filedates(:,2) = 1:N_files;
    filedates = sortrows(filedates,1);
end

vidfile = VideoWriter([saveDir filesep animal '_fullVideo.avi']);
framerate = 30;
vidfile.FrameRate = framerate;

open(vidfile);
t00 = tic;
all_frametimes = [];
bad_count = 0; bad_times = [];
allframenames = {framefiles.name};
allframe_us = regexp(allframenames,'_','split');
allframe_ending = cellfun(@(x) x{end}, allframe_us, 'UniformOutput', false);
allframe_dot = regexp(allframe_ending,'\.','split');
allframe_number = cellfun(@(x) str2double(x{1}), allframe_dot, 'UniformOutput', true)';

for ii = 1:N_files
    
    try
        t10 = tic;
        vidchunkidx = filedates(ii,2);
        vidchunkfile = [rawVideoDir filesep rawfiles(vidchunkidx).name];
        
        temp_underscore_split = regexp(vidchunkfile,'_','split');
        temp_dot_split = regexp(temp_underscore_split{end},'\.','split');
        vidchunkfile_number = str2double(temp_dot_split{1});
        
        timechunkidx = find(allframe_number == vidchunkfile_number);
        
        %         timechunkidx = frfiledates(ii,2); % this is the old way (ordering frametime files by date)
        timechunkfile = [rawVideoDir filesep framefiles(timechunkidx).name];
        
        if rawfiles(vidchunkidx).bytes > 1e6
            
            timechunk = load(timechunkfile);
            vidchunk = VideoReader(vidchunkfile);
            
            while hasFrame(vidchunk)
                frame = rgb2gray(readFrame(vidchunk));
                frame = imresize(frame,[389,519]);
                if isequal(size(frame),[389, 519])
                    writeVideo(vidfile,frame);
                end
                
            end
            
            disp('good file');
            
            all_frametimes = [all_frametimes; timechunk];
        else
            disp('bad file');
            bad_count = bad_count+1;
            bad_times(bad_count,1) = vidchunkidx;
            bad_times(bad_count,2) = timechunkidx;
        end
        t11 = toc(t10);
        fprintf('\n Took %3.2f seconds to write file %u out of %u.\n',t11,ii,N_files);
    catch
        bad_count = bad_count+1;
        bad_times(bad_count,1) = vidchunkidx;
        bad_times(bad_count,2) = timechunkidx;
        fprintf('Could not read file %u out of %u. Ignoring.\n',ii,N_files);
    end
end

close(vidfile);

% for AT12 video
need_cleanup = 1;
if need_cleanup
    % backup original frametimes
    fr_backup = all_frametimes;
    save([saveDir filesep animal '_frameTimes_BACKUP.mat'],'fr_backup');
    
    all_frametimes = sortrows(all_frametimes);
    fr_diff = diff(all_frametimes);
    fr_0 = find(fr_diff==0);
    all_frametimes(fr_0) = [];
end

putative_nframes = vidfile.FrameCount;

new_frametimes = all_frametimes;
n_removed = 0;
if ~isempty(bad_times);
    for xx = 1:size(bad_times,1);
        fprintf('bad chunk %u out of %u.\n',xx,size(bad_times,1));
        
        badchunkidx = bad_times(xx,1); % find bad video file index
        badchunkii = find(filedates(:,2)==badchunkidx); % find position of that index in filedates array
        badchunkfile = [rawVideoDir filesep rawfiles(badchunkidx).name]; % get bad file name
        temp_underscore_split = regexp(badchunkfile,'_','split');
        temp_dot_split = regexp(temp_underscore_split{end},'\.','split');
        badchunkfile_number = str2double(temp_dot_split{1}); % find number of bad video file
        badframeidx = find(allframe_number == badchunkfile_number); % match frame file number to that
        badframefile = [rawVideoDir filesep framefiles(badframeidx).name]; % get bad frametimes file
        
        lastgoodchunkii = badchunkii - 1; % get previous video file position in filedates
        if lastgoodchunkii > 0
            lastgoodidx = filedates(lastgoodchunkii,2); % retrieve previous video file index from filedates
            % make sure that index is not included in the list of bad indices
            while ismember(lastgoodidx,bad_times(:,1))
                lastgoodchunkii = lastgoodchunkii - 1; % if it is, go back one more index
                lastgoodidx = filedates(lastgoodchunkii,2);
                % repeat as necessary
            end
            
            % once we have our last good index, find the video file name it corresponds to
            lastgoodfile = [rawVideoDir filesep rawfiles(lastgoodidx).name];
            temp_underscore_split = regexp(lastgoodfile,'_','split');
            temp_dot_split = regexp(temp_underscore_split{end},'\.','split');
            lastgoodfile_number = str2double(temp_dot_split{1}); % get the video number for the last good video
            
            lastgoodframeidx = find(allframe_number == lastgoodfile_number); % match frametimes file number
            % get last good frametimes file
            lastgoodframefile = [rawVideoDir filesep framefiles(lastgoodframeidx).name];
            
            
            lastgoodframetimes = load(lastgoodframefile); % load last good frametimes
        else
            lastgoodframetimes = [];
        end
        
        badframetimes = load(badframefile);  % load bad frametimes
        remove_these = setdiff(badframetimes,lastgoodframetimes); % find frametimes to remove
        % loop through those and remove them, count removed frametimes
        for rr = 1:size(remove_these,1)
            remove_these_idx = find(new_frametimes == remove_these(rr));
            new_frametimes(remove_these_idx) = [];
            n_removed = n_removed + size(remove_these_idx,1);
            clear remove_these_idx
        end
        clear remove_these* lastgood* badchunk* badframe* temp_*
    end
end

% keyboard;
need_second_cleanup = 0;
if need_second_cleanup
    new_frametimes2 = sortrows(new_frametimes);
    fr_diff = diff(new_frametimes2);
    fr_0 = find(fr_diff==0);
    new_frametimes2(fr_0) = [];
end




t01 = toc(t00);
fprintf('\n Took %3.2f seconds to write full video (%u files, ~%u hours).\n',t01,N_files,N_files/4);
fprintf('Now saving frametimes...\n');
frametimes = new_frametimes;
save([saveDir filesep animal '_frameTimes.mat'],'frametimes');

fprintf('DONE!\n\n');

