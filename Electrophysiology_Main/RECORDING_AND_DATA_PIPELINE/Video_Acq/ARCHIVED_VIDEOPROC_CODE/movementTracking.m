%% MOVEMENT TRACKING ALGORITHM
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

%% LOOK FOR AVI VIDEO AND CREATE IF NON-EXISTENT
% select directories and animal identifier
maindir = uigetdir(cd,'Select folder with raw video data. ');
saveloc = uigetdir(cd,'Where do you want to save the processed video and movement data? ');
animal = input('Animal(s) name?  ','s');

% check for avi file
if ~exist([saveloc filesep animal '_fullVideo.avi'],'file')
    % if avi video file does not exist, create it from raw files
    fprintf('Processed video file not found.\nMaking avi video file from raw video...\nThis will run overnight.\n\n');
    framerate   = 30;
    [~, videofilename, videotimesfile, frameTimestamps] = ...
        createVidFile_fxn(maindir,saveloc,animal,framerate);
    nFrames = size(frameTimestamps,1);
    fprintf('Done!\n');
else
    % otherwise proceed
    videofilename = [saveloc filesep animal '_fullVideo.avi'];
    fprintf('Processed video file found: %s.\n',videofilename);
end

%% INITIALIZE COMPUTER VISION OBJECTS

% video file reader
vidSrc = vision.VideoFileReader(videofilename,'ImageColorSpace','Intensity');

% data converter
converter = vision.ImageDataTypeConverter('OutputDataType','uint8');

% optical flow vectors calculator
opticalFlow = vision.OpticalFlow('ReferenceFrameSource','Property',...
    'ReferenceFrameDelay',1,'Smoothness',0.8);

% shape and text inserters
lineInserter = vision.ShapeInserter('Shape','Lines','BorderColor','White');
rectInserter  = vision.ShapeInserter('Shape','Rectangles',...
    'BorderColor','White','LineWidth',2);
htextins = vision.TextInserter('Text', '%4d', 'Location',  [1 1], ...
    'Color', [1 1 1], 'FontSize', 12);

% image filter calculators
Mean1 = vision.Mean;
Mean2 = vision.Mean('RunningMean', true);
medianFilter = vision.MedianFilter;

% blob analysis object
% minimum blob area should be 250 for current recordings
% maximum blob area should be 1000 for current recordings
blobAnalysis = vision.BlobAnalysis(...
    'CentroidOutputPort', true, 'AreaOutputPort', true, ...
    'BoundingBoxOutputPort', true, 'OutputDataType', 'double', ...
    'MinimumBlobAreaSource', 'Property', 'MinimumBlobArea', 200, ...
    'MaximumBlobArea', 2000, 'MaximumCount', 10);

%% SELECT ROIs IN VIDEO

% show first frame of video
firstframe = step(vidSrc);
imshow(firstframe);
set(gcf,'name','Select corners of ROI. Close figure to continue.',...
    'numbertitle','off');
msg1 = msgbox({'Select corners of ROI, clockwise from top left corner.' ...
    'Close this figure when you are finished to continue.'},...
    'ATTENTION!','modal');

% user input selects ROI
[ROIx,ROIy] = ginput(4);
ROI_left = min(ROIx); ROI_bottom = min(ROIy);
ROI_width = max(ROIx) - ROI_left; ROI_height = max(ROIy) - ROI_bottom;
ROI_right = ROI_left + ROI_width; ROI_top = ROI_bottom + ROI_height;
rectangle('Position',[ROI_left ROI_bottom ROI_width ROI_height],'EdgeColor','r','linewidth',1);
release(vidSrc);
uiwait(gcf);
% when figure is closed, code proceeds


%% CALCULATE AND STORE MEAN INTENSITY OF EACH FRAME FOR L/D DISCRIMINATION
% get mean acquisition framerate and interval between frames
acq_rate    = nFrames/(frameTimestamps(end) - frameTimestamps(1));
frame_int   = 1/acq_rate;
jump12hour  = ceil(12*3600 / frame_int);
jump = 0; % this usually set to 1. Set to 0 if you want to look at all frames

% use the first few hours of video to find the threshold in mean intensity
% separating light and dark periods.
avgIntensity = zeros(nFrames,1);
counter0 = 0;
while counter0 <= 1000
    counter0    = counter0 + 1;
    
    % Get a sample of light and dark frames
    if mod(counter0,10) == 0 && jump
        % skip forward 12 hours
        for x = 1:jump12hour
            frame = step(vidSrc);
        end
    else
        % take frame
        frame = step(vidSrc);
    end
    
    % Get average intensity of frame
    avgIntensity(counter0,1) = mean(mean(frame));
    
    if mod(counter0,100)==0
        fprintf('Done with intensity calculation for %u frames.\n',counter0);
    end
    
end

% plot intensity of frames and use ginput to divide the two and set
% threshold manually
figure(100)
plot(avgIntensity);
[~,thresh] = ginput(1);
close(100);

fprintf('Intensity threshold set at %.3f.\n',thresh);


%% INITIALIZE VIDEO PLAYERS AND RESET VIDEO TO START
fprintf('\n\n     ____Initializing video players____     \n\n');
pos = [20 20 350 350];
videoPlayer1 = vision.VideoPlayer('Name','video1','Position',pos);
pos(1) = pos(1) + 360;
videoPlayer2 = vision.VideoPlayer('Name','video2','Position',pos);
pos(1) = pos(1) + 360;
videoPlayer3 = vision.VideoPlayer('Name','video3','Position',pos);

reset(vidSrc);


%% DISPLAY SELECTED ROI AND INITIALIZE MOVEMENT AND TRACK ARRAYS
f_roi = figure(); hold on;
imshow(firstframe);
ROIrect = rectangle('Position',[ROI_left ROI_bottom ROI_width ROI_height],...
    'EdgeColor','r','linewidth',1.5);
set(f_roi,'name','First frame with ROI');

trackstart_x = ROI_left + ROI_width/2;
trackstart_y = ROI_bottom + ROI_height/2;
track = [];
f_track = figure();
trackline = animatedline(trackstart_x,trackstart_y,'Color','g','Marker','x');

f_move = figure();
movement = [];
moveplot = animatedline('Color','b','LineWidth',2,'Marker','none');

%% MAIN LOOP - STEP THROUGH VIDEO AND DO IMAGE PROCESSING
counter = 0; checked = 0;
while ~isDone(vidSrc)
    % frame counter
    counter = counter + 1;
    
    % Load next frame
    frame = step(vidSrc);
    
    % Find optical flow vectors
    ofVectors = step(opticalFlow, frame);
    
    % Compute magnitude squared of complex OF vectors
    ofMagn = ofVectors .* conj(ofVectors);
    
    % Compute velocity threshold
    velThreshFactor = 0.3;
    velThresh = velThreshFactor * step(Mean2, step(Mean1, ofMagn));
    
    % Threshold and filter image to find objects
    segmentedObjs = step(medianFilter, ofMagn >= velThresh);
    
    % Make shapes for eroding/closing
    SE1 = strel('line',5,45);
    SE2 = strel('rectangle',[2,3]);
    SE3 = strel('rectangle',[4,4]);
    SE4 = strel('diamond',6);
    
    % apply erosion/closing to detected objects
    if mean(mean(frame)) >= thresh
        % processing for light frames
        erode1 = imerode(segmentedObjs,SE1);
        close1 = imclose(erode1,SE3);
        close2 = imclose(close1,SE4);
        
    elseif mean(mean(frame)) < thresh
        % processing for dark frames
        erode1 = imerode(segmentedObjs,SE1); 
        erode2 = imerode(erode1,SE2);
        close1 = imclose(erode1,SE3);
        close2 = imclose(close1,SE4);
        
    end
    
    % Find area, centroids and bounding boxes for detected objects
    [area, centroid, bbox] = step(blobAnalysis, close2);
   
    % --------
    % Main processing is completed at this point. The next step is to look
    % for objects that were detected within the selected ROI and find the
    % animal within those.
    
    % Select blobs with bounding boxes within ROI
    x_idx = arrayfun(@(x) iswithin(x,ROI_left,ROI_right,1),centroid(:,1));
    y_idx = arrayfun(@(x) iswithin(x,ROI_bottom,ROI_top,1),centroid(:,2));
    roiIdx = x_idx & y_idx;
    theseBlobs  = bbox(roiIdx,:);
    
    % Use ratio of areas to find rat blob
    ratio = zeros(length(roiIdx),1);
    ratio(roiIdx) = single(area(roiIdx,1))./single(bbox(roiIdx,3).*bbox(roiIdx,4));
    
    % threshold used depends on whether frame is in light or dark period
    lightThresh = 0.25;
    darkThresh  = 0.25;
    if mean(mean(frame)) > thresh
        ratioThresh = lightThresh;
    elseif mean(mean(frame)) < thresh
        ratioThresh = darkThresh;
        
    end
    ratioThIdx = ratio > ratioThresh;
    
    
    % find the right blob
    blobIdx = find(roiIdx & ratioThIdx);
    if size(blobIdx,1) > 1
        blobIdx = blobIdx(1,1);
    end
    
    % Apply object inserter for blob bounding boxes and centroids
    blobObjs = step(rectInserter,frame,bbox(blobIdx,:));
    centrObjs = insertMarker(frame,centroid(blobIdx,:));
    
    % compute movement track of animal
    if ~isempty(centroid(blobIdx,:))
        track(counter,:) = centroid(blobIdx,:);
        track(counter,2) = ROI_top - (track(counter,2)-ROI_bottom);
    else
        if counter > 1
            track(counter,:) = track(counter-1,:);
        else
            track(counter,:) = [trackstart_x trackstart_y];
        end
    end
    
    % compute movement magnitude of animal
    if counter>1
        movement(counter-1,1) = counter-1;
        movement(counter-1,2) = norm(track(counter,:)-track(counter-1,:));
        movement(counter-1,3) = frameTimestamps(counter-1,1);
    end
    
    
    % Play videos
    % only play for first few frames (set by "ntestframes")
    ntestframes = 1000;
    if counter < ntestframes && checked == 0
        out = step(rectInserter,frame,bbox);
        step(videoPlayer1, blobObjs);    % show bounding boxes
        step(videoPlayer2, close2);
        step(videoPlayer3, segmentedObjs);
    elseif checked == 0 && counter >= ntestframes
        % user can choose to abort and change parameters or proceed
        proceed = input('Continue? 1 for yes, 0 for no');
        checked = 1;
    end
    
    if counter >= ntestframes && proceed == 0
        % release objects
        disp('releasing');
        release(videoPlayer1);
        release(videoPlayer2);
        release(videoPlayer3);
        release(vidSrc);
        reset(vidSrc);
        close all
        error('user aborted operation');
    end
    
    % show Movement Figure
    if counter>1
        figure(f_move);
        addpoints(moveplot,movement(counter-1,1),movement(counter-1,2));
    end
    
    
    if mod(counter,1000) == 0
        fprintf('Processed %u frames out of %u.\n',counter,n_frames);
    end
    
    % for testing - ask to stop/proceed every 500 frames
    if mod(counter,1000) == 0
        keepgoing = input('continue? (y/n):  ','s');
        if strcmp(keepgoing,'n') == 1
            release(videoPlayer2);
            release(videoPlayer1);
            release(videoPlayer3);
            release(vidSrc);
            close all;
            break;
        end
    end
end

fprintf('Done running video.\n');

%% SMOOTH MOVEMENT DATA
% Eliminate isolated movement spikes likely coming from detection of small
% movement of metal cables while animal is essentially immobile.
% movement smoothing.
% This loop takes a movement point x and looks at movement points between
% "x-window" and "x+window". If x is non-zero and followed by a string of
% zeros of length "window" in both directions (past and future), it is
% detected as an isolated spike and set to 0.
window = 10;
for ee = 1:size(movement,1)
    if ee > window && movement(ee,2) > 0 && ee < size(movement,1)-window
        sum_previous = sum(movement(ee-window:ee-1,2));
        sum_following = sum(movement(ee+1:ee+window,2));
        if sum_previous == 0 && sum_following == 0
            movement(ee,2) = 0;
        end
    end
end

% This next loop performs a simple moving average smoothing on the movement
% data stripped of isolated spikes. The variable "span" can be increased to
% for harder smoothing.
% Both the raw and smoothed movement data will be saved.
span = 5;
smooth_mvt = smooth(movement,span);


%% SAVE DATA AND EXIT
% make data structure for saving
outdata.raw_movement        = movement;
outdata.smooth_movement     = smooth_mvt;
outdata.track               = track;
outdata.ROI                 = [ROI_left ROI_bottom ROI_width ROI_height];
outdata.firstframe          = firstframe;
outdata.meanIntensity       = avgIntensity;
outdata.intThresh           = thresh;
outdata.animal              = animal;
outdata.nframes             = nFrames;
outdata.params.SE1          = SE1;
outdata.params.SE2          = SE2;
outdata.params.SE3          = SE3;
outdata.params.SE4          = SE4;
outdata.params.vThreshFact  = velThreshFactor;
outdata.params.ratioThreshL = lightThresh;
outdata.params.ratioThreshD = darkThresh;
outdata.params.smooth_span  = span;

% save in processed video directory
fprintf('\nSaving data...\n');
tic

save([saveloc filesep animal '_movement.mat'],'outdata');

toc
fprintf('Done!\n\n');

% release objects
release(videoPlayer1);
release(videoPlayer2);
release(videoPlayer3);
release(vidSrc);
reset(vidSrc);


