%% MOVEMENT TRACKING ALGORITHM - VERSION BY CHUNKS
%
% INPUTS NEEDED:
%   1. AVI file path
%   2. save filename (with batch counter)
%   3. Save location
%   4. animal names
%   5. CSV file with frame timestamps
%   6. ROIs (masks) for each animal
%   7. chunk size - how many frames to read/process at once
%
% This is the only script needed for movement tracking of animals.
% The script goes through the following steps:
%   1. Look for AVI video file. If non-existent, create it from raw video
%      files acquired by the MacMini
%   2. Calculate mean intensity of each frame to determine whether each
%      frame was captured in the light or dark period
%   3. Detect moving objects in video. This relies on a few techniques: a)
%      optical flow vector calculation; b) median filtering and
%      velocity-based thresholding; c) erosion and closure of segmented
%      objects; d) blob analysis on resulting B/W image
%   4. Find the objects that were moving in the relevant ROI. Use ratio of
%      blob area to bounding box size to determine which box corresponds to
%      the rat. Then calculate movement data for that blob by comparing
%      positions across adjacent frames.
%   5. Store and save data.
%
%                      ____| IMPORTANT NOTE |___
% If you want to change this code, please use "movementTracking_play.m".
% That is a testing code that contains the main loop of this file, with
% some added testing plots that can be commented in and out. You can mess
% with that and then substitute it back in here.
%
% ATP, April 2016
%
%


% for the sake of  constructing the scripts, assume inputs are as follows:
% videofilename = path to AVI video file
% saveloc = save location
% savefile1/2 = save name for this block of video data (animal1/2)
% animal1/2 = animal1/2 name
% csvfile = frame timestamps in CSV format
% mask1 = ROI for animal 1
% mask2 = ROI for animal 2

% **** NOTE:
%    for saving, if this is going to process ~4 hrs of data, it should save
%    with a name that makes sense. The counter for the save filename should
%    be in the bash cluster script. The name should be passed as an input
%    to this script (see above)


%% LOAD FRAME TIMES

% remember these are in GMT time!
frame_times = load(csvfile);

% Store total number of frames in video
n_stamps = size(frame_times,1);

%% INITIALIZE COMPUTER VISION OBJECTS

% video file reader
vidSrc = VideoReader(videofilename);
nvid_frames = vidSrc.Duration*vidSrc.FrameRate;

% Check that number of frames is equal, as calculated from frame times and
% video file
if ~isequal(nvid_frames,n_stamps)
    fprintf('Problem: founrd %u timestamps but %u frames.\n',n_stamps,nvid_frames);
    keyboard
else
    fprintf('This video has %u frames and %u timestamps.\n\n',nvid_frames,n_stamps);
end

%% INITIALIZE MOVEMENT AND TRACK ARRAYS
% one of each for each animal

track_1 = zeros(nvid_frames,2);
track_2 = zeros(nvid_frames,2);

movement_1 = zeros(nvid_frames,2);
movement_2 = zeros(nvid_frames,2);

% initial position is middle of mask
% animal 1
trackstart_x1 = mask1(1,1) + (mask1(2,1)-mask1(1,1))/2;
trackstart_y1 = mask1(1,2) + (mask1(2,2)-mask1(1,2))/2;

% animal 2
trackstart_x2 = mask2(1,1) + (mask2(2,1)-mask2(1,1))/2;
trackstart_y2 = mask2(1,2) + (mask2(2,2)-mask2(1,2))/2;

%% SETTINGS FOR BACKGROUND SUBTRACTION
% Set "update proportion" parameter. Increasing this will make the update
% of the background more reactive, but will blend slow-moving objects into
% the background.
alpha = 0.15;

% Initialize background as first frame of video
firstframe = readFrame(vidSrc);
old_background = rgb2gray(firstframe);

%% INITIALIZE PARAMETERS FOR VIDEO ANALYSIS
% master frame counter
counter = 0;

% Make shapes for eroding/closing
SE1 = strel('line',5,45);
SE2 = strel('rectangle',[3,3]);
SE3 = strel('rectangle',[4,4]);
SE4 = strel('diamond',6);

% erode/close params
n_erode     = 2;
n_close     = 2;

% blob analysis params
connThreshold = 18;
areaThreshold = 10;

% thresholding params
thresholdFactor = 1.75; % multiply thresholdLevel by this

% chunking params
frame_chunk_sz = 2000;
nframe_chunks = ceil(nvid_frames/frame_chunk_sz);

% redefine videoreader object to avoid error
vidSrc = VideoReader(videofilename);

%% MAIN LOOP - STEP THROUGH VIDEO AND DO IMAGE PROCESSING

% loop through chunks
for f_chunk = 1:nframe_chunks
    % define chunk limits
    startpos = (f_chunk-1)*frame_chunk_sz+1;
    endpos   = min(f_chunk*frame_chunk_sz,nvid_frames);
    
    % read frames
    t10 = tic;
    framechunk = read(vidSrc,[startpos endpos]);
    t11 = toc(t10);
    fprintf('Reading chunk %u out of %u took %.1f seconds.\n',f_chunk,nframe_chunks,t11);
    
    if endpos < frame_chunk_sz
        endloop = mod(endpos,frame_chunk_sz);
    else
        endloop = frame_chunk_sz;
    end
    
    % initialize temporary movement and track chunks
    temp_track1  = zeros(endloop,2);
    temp_mvt1    = zeros(endloop,1);
    temp_track2  = zeros(endloop,2);
    temp_mvt2    = zeros(endloop,1);
    
    % start processing for current chunk
    t20 = tic;
    for frame_idx = 1:endloop
        
        % frame counter
        counter = counter + 1;
        % Load next frame
        frame = framechunk(:,:,:,frame_idx);
        frame = rgb2gray(frame);
        
        % Remove noise from image via median filtering
        %     clean_img = step(medianFilter, frame);
        clean_img = medfilt2(frame);
        
        
        % __Adaptive background subtraction
        % Calculate background
        if counter > 1
            new_background = old_background*(1-alpha) + clean_img*alpha;
        end
        
        % Remove background from image.
        if counter == 1
            foreground = abs(clean_img - old_background);
        else
            foreground = abs(clean_img - new_background);
            old_background = new_background; % update background for next frame
        end
        
        % __Thresholding image:
        % Frames should be in gray format already.
        thresholdLevel = graythresh(foreground); % find good threshold level
        
        binaryImage = foreground > thresholdFactor * thresholdLevel; % threshold image
        
        % erode/close operations on binary image
        count_erode = 0;
        count_close = 0;
        imgToProcess = binaryImage;
        
        while count_erode < n_erode
            count_erode = count_erode + 1;
            imgToProcess = imerode(imgToProcess,SE2);
        end
        while count_close < n_close
            count_close = count_close + 1;
            imgToProcess = imclose(imgToProcess,SE3);
        end

        newBinaryImage = imgToProcess;
        
        % apply both masks in parallel 
        for anim = 1:2
            eval(['binaryMask = mask' num2str(anim) ';']);
            
            % Apply binary mask
            newBinaryImage(binaryMask == 0) = 0;
            
            % FIND BLOBS USING CONNECTED COMPONENT ANALYSIS
            connComps   = bwconncomp(newBinaryImage,connThreshold);
            
            blobs       = regionprops(connComps,'Centroid','Area','BoundingBox','Extent');
            
            % store blob info in arrays
            extents         = cat(1,blobs.Extent);
            boxes           = cat(1,blobs.BoundingBox);
            centroids       = cat(1,blobs.Centroid);
            areas           = cat(1,blobs.Area);
            
            % eliminate small blobs
            largeEnough     = areas > areaThreshold;
            areas           = areas(largeEnough);
            extents         = extents(largeEnough);
            boxes           = boxes(largeEnough,:);
            centroids       = centroids(largeEnough,:);
            
            % find largest blob. this selects the rat in the vast majority of cases
            thisBlob = find(areas==max(areas));
            if numel(thisBlob) > 1
                thisBlob = thisBlob(extents(thisBlob)==max(extents(thisBlob)));
            end
            if numel(thisBlob) > 1
                thisBlob = thisBlob(1);
            end
            eval(['blobCenter' num2str(anim) ' = centroids(thisBlob,:);']);
        end
        
        % compute movement track of each animal
        if ~isempty(blobCenter)
            % animal 1
            temp_track1(frame_idx,:) = blobCenter1;
            temp_track1(frame_idx,2) = mask1(2,2) - (temp_track1(frame_idx,2)-mask1(1,2));
            
            % animal 2
            temp_track2(frame_idx,:) = blobCenter2;
            temp_track2(frame_idx,2) = mask2(2,2) - (temp_track2(frame_idx,2)-mask2(1,2));
        else
            if frame_idx > 1
                % animal 1
                temp_track1(frame_idx,:) = temp_track1(frame_idx-1,:);
                
                % animal 2
                temp_track2(frame_idx,:) = temp_track2(frame_idx-1,:);
            else
                if counter > 1
                    % animal 1
                    temp_track1(frame_idx,:) = track_1(counter-1,:);
                    
                    % animal 2
                    temp_track2(frame_idx,:) = track_2(counter-1,:);
                else
                    % animal 1
                    temp_track1(frame_idx,:) = [trackstart_x1 trackstart_y1];
                    
                    % animal 2
                    temp_track2(frame_idx,:) = [trackstart_x2 trackstart_y2];
                end
            end
        end
        
        % compute movement magnitude of animal
        if frame_idx > 1
            % animal 1
            temp_mvt1(frame_idx-1,1) = norm( temp_track1(frame_idx,:) - temp_track1(frame_idx-1,:) );
            
            % animal 2
            temp_mvt2(frame_idx-1,1) = norm( temp_track2(frame_idx,:) - temp_track2(frame_idx-1,:) );
        end
        
        if mod(counter,1000) == 0
            fprintf('Processed %u frames out of %u.\n',counter,nvid_frames);
        end
        
    end
    t21 = toc(t20);
    fprintf('Blob analysis on chunk %u out of %u took %.1f seconds.\n',f_chunk,nframe_chunks,t21);
    
    % smooth temp movement data and store temp variables in stable arrays
    t30 = tic;
    
    % animal 1
    temp_smooth1 = smooth(temp_mvt1,smooth_factor,'rloess');
    
    movement_1(start_pos:end_pos,1) = frame_times(start_pos:end_pos,1);
    movement_1(start_pos:end_pos,2) = temp_mvt1(:,1);
    
    track_1(start_pos:end_pos,:) = temp_track1(:,:);
    
    smoothed_mvt_1(start_pos:end_pos,1) = frame_times(start_pos:end_pos,1);
    smoothed_mvt_1(start_pos:end_pos,1) = temp_smooth1(:,1);
    
    % animal 2
    temp_smooth2 = smooth(temp_mvt2,smooth_factor,'rloess');
    
    movement_2(start_pos:end_pos,1) = frame_times(start_pos:end_pos,1);
    movement_2(start_pos:end_pos,2) = temp_mvt2(:,1);
    
    track_2(start_pos:end_pos,:) = temp_track2(:,:);
    
    smoothed_mvt_2(start_pos:end_pos,1) = frame_times(start_pos:end_pos,1);
    smoothed_mvt_2(start_pos:end_pos,1) = temp_smooth2(:,1);
    
    % save the data every 10 chunks (to minimize losses if it crashes)
    if mod(f_chunk,10) == 0
        % animal 1
        temp_savename1 = [animal1 '_tempsave_chunk_' num2str(f_chunk) '.mat'];
        save(temp_savename1,'movement_1','smoothed_mvt_1','track_1','-v7.3')
        
        % animal 2
        temp_savename2 = [animal2 '_tempsave_chunk_' num2str(f_chunk) '.mat'];
        save(temp_savename2,'movement_2','smoothed_mvt_2','track_2','-v7.3')
    end
    
    t31 = toc(t30);
    fprintf('Smoothing and storing data for chunk %u out of %u took %.1f seconds.\n',f_chunk,nframe_chunks,t31);
    
    % clear temp variables for next chunk
    clear temp_*
end



%% SAVE DATA AND EXIT
% make data structure for saving

% animal 1
DATA.smooth_movement        = smoothed_mvt_1;
DATA.frame_times            = frame_times;
DATA.nframes                = nvid_frames;
DATA.mask                   = mask1;
RAW.raw_movement            = movement_1;
RAW.track                   = track_1;
eval([animal1 '_out.DATA = DATA;']);
eval([animal1 '_out.RAW = RAW;']);
clear DATA RAW

% animal 2
DATA.smooth_movement        = smoothed_mvt_2;
DATA.frame_times            = frame_times;
DATA.nframes                = nvid_frames;
DATA.mask                   = mask2;
RAW.raw_movement            = movement_2;
RAW.track                   = track_2;
eval([animal2 '_out.DATA = DATA;']);
eval([animal2 '_out.RAW = RAW;']);
clear DATA RAW


% NOTE: only save params if it is the first 4hr chunk that is being
% processed.
if somethingsomething
    for aa=1:2
        params.SE1                  = SE1;
        params.SE2                  = SE2;
        params.SE3                  = SE3;
        params.SE4                  = SE4;
        params.smoothfactor         = smooth_factor;
        params.alpha                = alpha;
        params.connection_thresh    = connThreshold;
        params.area_thresh          = areaThreshold;
        params.gray_thresh_factor   = thresholdFactor;
        params.n_erode              = n_erode;
        params.n_close              = n_close;
        params.erode_first          = erode_first;
        eval(['animstring = animal' num2str(aa) ';']);
        eval([animstring '_out.params = params;']);
    end
end

% save in processed video directory
fprintf('\nSaving data...\n');
ts0 = tic;

varname1 = [animal1 '_out'];
save([saveloc filesep savefile1],varname1);

varname2 = [animal2 '_out'];
save([saveloc filesep savefile2],varname2);

ts1 = toc(ts0);
fprintf('Done! Saving both animals'' data took %.2f seconds.\n\n',ts1);


