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

%% SETUP FILES AND DIRECTORIES
% select avi file and save directory
[vf, vd] = uigetfile('*.avi','Pick your AVI video file.',cd);
videofilename = fullfile(vd,vf);
[ftf, ftd] = uigetfile('*.mat','Pick your frame times.',cd);
saveloc = uigetdir(cd,'Where do you want to save movement data? ');
animal = input('Which animal are you running?  ','s');

ft_data = load(fullfile(ftd,ftf));
frame_times = ft_data.frametimes;

% Store total number of frames in video
n_frames = size(frame_times,1);

%% INITIALIZE COMPUTER VISION OBJECTS

% video file reader
vidSrc = VideoReader(videofilename);
nvid_frames = vidSrc.Duration*vidSrc.FrameRate;
vidSrc.CurrentTime = 0;


if ~isequal(nvid_frames,n_frames)
    fprintf('Problem: founrd %u timestamps but %u frames.\n',n_frames,nvid_frames);
    keyboard
else
    fprintf('This video has %u frames.\n\n',nvid_frames);
end

% shape and text inserters
lineInserter = vision.ShapeInserter('Shape','Lines','BorderColor','White');
rectInserter  = vision.ShapeInserter('Shape','Rectangles',...
    'BorderColor','White','LineWidth',2);
htextins = vision.TextInserter('Text', '%4d', 'Location',  [1 1], ...
    'Color', [1 1 1], 'FontSize', 12);


%% SELECT ROIs IN VIDEO

% show first frame of video
firstframe = readFrame(vidSrc);
imshow(firstframe);
set(gcf,'name','Select corners of ROI. Close figure to continue.',...
    'numbertitle','off');
msg1 = msgbox({'Select corners of ROI (top left, then bottom right).' ...
    'Close this figure when you are finished to continue.'},...
    'ATTENTION!','modal');

for cc=1:2
    [x(cc,1),y(cc,1)]=ginput(1);
    x(x<1)=1;
    x(x>size(firstframe,2)) = size(firstframe,2);
    y(y<1)=1;
    y(y>size(firstframe,1)) = size(firstframe,1);
end
rectangle('Position',[x(1) y(1) x(2)-x(1) y(2)-y(1)],'EdgeColor','r','linewidth',1.5);
ROI = round([x y]);

% make binary image mask
binaryMask = zeros(size(firstframe,1),size(firstframe,2));
binaryMask(ROI(1,2):ROI(2,2),ROI(1,1):ROI(2,1)) = 1;

uiwait(gcf); % when figure is closed, code proceeds
vidSrc.CurrentTime = 0;

%% CALCULATE AND STORE MEAN INTENSITY OF EACH FRAME FOR L/D DISCRIMINATION

avgIntensity = [];
counter0 = 0;
% take 1 frame every "frame_interval"
frame_interval = 1500;
time_points = ceil(nvid_frames/frame_interval);
time_vector = linspace(0,vidSrc.Duration,time_points);
for xx = 1:size(time_vector,2)
    thistime = time_vector(xx);
    fprintf('Timepoint %u out of %u for intensity calculation.\n',xx,time_points);
    vidSrc.CurrentTime = thistime;
    if hasFrame(vidSrc)
        thisframe = rgb2gray(readFrame(vidSrc));
    else
        vidSrc.CurrentTime = vidSrc.CurrentTime - 1;
        thisframe = rgb2gray(readFrame(vidSrc));
    end
    % Get average intensity of frame
    avgIntensity(end+1,1) = mean(mean(thisframe));
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



%% DISPLAY SELECTED ROI AND INITIALIZE MOVEMENT AND TRACK ARRAYS
f_roi = figure(); hold on;
imshow(firstframe);
ROIrect = rectangle('Position',[ROI(1,1) ROI(1,2) ROI(2,1)-ROI(1,1) ROI(2,2)-ROI(1,2)],...
    'EdgeColor',[0.22 1.00 0.08],'linewidth',2.5);
set(f_roi,'name','First frame with ROI');

trackstart_x = ROI(1,1) + (ROI(2,1)-ROI(1,1))/2;
trackstart_y = ROI(1,2) + (ROI(2,2)-ROI(1,2))/2;
track = [];
f_track = figure();
trackline = animatedline(trackstart_x,trackstart_y,'Color',[0 .45 .74],'Marker','x','MaximumNumPoints',50);

f_move = figure();
movement = [];
moveplot = animatedline('Color',[.85 .33 .10],'LineWidth',2,'Marker','none');

%% SETTINGS FOR BACKGROUND SUBTRACTION
% Set "update proportion" parameter. Increasing this will make the update
% of the background more reactive, but will blend slow-moving objects into
% the background.
alpha = 0.15;

% Initialize background as first frame of video
old_background = rgb2gray(firstframe);

%% MAIN LOOP - STEP THROUGH VIDEO AND DO IMAGE PROCESSING
counter = 0; checked = 0;
release(videoPlayer2);
release(videoPlayer1);
release(videoPlayer3);
start_at = 0;

vidSrc = VideoReader(videofilename);

vidSrc.CurrentTime = start_at/vidSrc.FrameRate;
counter = start_at;
% show Movement Figure
plotfigs = 0;
checkmode = 1;
while hasFrame(vidSrc)
       
    % frame counter
    counter = counter + 1;
    % Load next frame
    frame = readFrame(vidSrc);
    frame = rgb2gray(frame);
    
    % Remove noise from image via median filtering
%     clean_img = step(medianFilter, frame);
    % Francesco's improved median filtering:
    clean_img = medfilt2(frame);

    
    % __Adaptive background subtraction
    % Calculate background
    if counter > start_at+1
        new_background = old_background*(1-alpha) + clean_img*alpha;
    end
    
    % Remove background from image.
    if counter == start_at+1
        foreground = abs(clean_img - old_background);
    else
        foreground = abs(clean_img - new_background);
        old_background = new_background; % update background for next frame
    end
    
    % __Thresholding image:
    % Frames should be in gray format already.
    thresholdLevel = graythresh(foreground); % find good threshold level
    thresholdFactor = 1.75; % multiply thresholdLevel by this
    binaryImage = foreground > thresholdFactor * thresholdLevel; % threshold image
    binaryImage = medfilt2(binaryImage);
    
    % Make shapes for eroding/closing
    SE1 = strel('line',5,45);
    SE2 = strel('rectangle',[3,3]);
    SE3 = strel('rectangle',[4,4]);
    SE4 = strel('diamond',6);
    
    n_erode     = 2;
    n_close     = 2;
    count_erode = 0;
    count_close = 0;
    erode_first = 1;
    imgToProcess = binaryImage;
    
    % apply erosion/closing to binary image
    if mean(mean(frame)) >= thresh
        
        % this switch determines the order of eroding/closing operations
        if erode_first
            % processing for light frames
            while count_erode < n_erode
                count_erode = count_erode + 1;
                %             disp(['light erode ' num2str(count_erode)]);
                imgToProcess = imerode(imgToProcess,SE2);
            end
            while count_close < n_close
                count_close = count_close + 1;
                %             disp(['light close ' num2str(count_close)]);
                imgToProcess = imclose(imgToProcess,SE3);
            end
        else
            while count_close < n_close
                count_close = count_close + 1;
                %             disp(['light close ' num2str(count_close)]);
                imgToProcess = imclose(imgToProcess,SE2);
            end
            while count_erode < n_erode
                count_erode = count_erode + 1;
                %             disp(['light erode ' num2str(count_erode)]);
                imgToProcess = imerode(imgToProcess,SE3);
            end
        end
        
        newBinaryImage = imgToProcess;
        
    elseif mean(mean(frame)) < thresh
        
%         imgToProcess = imcomplement(binaryImage);
        
        % processing for dark frames
        while count_erode < n_erode
            count_erode = count_erode + 1;
            %             disp(['dark erode ' num2str(count_erode)]);
            imgToProcess = imerode(imgToProcess,SE3);
        end
        while count_close < n_close
            count_close = count_close + 1;
            %             disp(['dark close ' num2str(count_close)]);
            imgToProcess = imclose(imgToProcess,SE3);
        end
        
        newBinaryImage = imgToProcess;
        
    end
    
    % Apply binary mask
    newBinaryImage(binaryMask == 0) = 0;
    
    % FIND BLOBS USING CONNECTED COMPONENT ANALYSIS
    % Steps:
    %   1. Find centroids, areas, bounding boxes and extent (filled pixel
    %   area/bounding box area)
    connThreshold = 18;
    areaThreshold = 10;
    
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
    blobCenter  = centroids(thisBlob,:);
    
    
    % Apply object inserter for blob bounding boxes and centroids
    if ~isempty(blobCenter)
%         blobObjs = step(rectInserter,frame,boxes(thisBlob,:));
        centrObjs = insertMarker(frame,blobCenter);
    else
        if counter > start_at+1
            centrObjs = insertMarker(frame,[ track(counter-start_at-1,1) ( ROI(2,2) - track(counter-start_at-1,2) + ROI(1,2) ) ] );
        else
            centrObjs = insertMarker(frame,[trackstart_x trackstart_y]);
        end
    end
    
    
    % compute movement track of animal
    if ~isempty(blobCenter)
        track(counter-start_at,:) = blobCenter;
        track(counter-start_at,2) = ROI(2,2) - (track(counter-start_at,2)-ROI(1,2));
    else
        if counter > start_at+1
            track(counter-start_at,:) = track(counter-start_at-1,:);
        else
            track(counter-start_at,:) = [trackstart_x trackstart_y];
        end
    end
    
    % compute movement magnitude of animal
    if counter>start_at+1
        movement(counter-start_at-1,1) = counter-start_at-1;
        movement(counter-start_at-1,2) = norm(track(counter-start_at,:)-track(counter-start_at-1,:));
        % PUT THIS BACK IN
        movement(counter-start_at-1,3) = frame_times(counter-start_at-1,1);
    end
    
    
    % Play videos
    % only play for first few frames (set by "ntestframes")
    ntestframes = 0;
    if ntestframes ~= 0
        if counter < ntestframes && checked == 0
            %         out = step(rectInserter,frame,bbox);
            step(videoPlayer1, centrObjs);    % show bounding boxes
            step(videoPlayer2, newBinaryImage);
            step(videoPlayer3, frame);
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
            close all
            error('user aborted operation');
        end
    end
    
    if checkmode
        step(videoPlayer1, centrObjs);    % show bounding boxes
        step(videoPlayer2, newBinaryImage);
        step(videoPlayer3, frame);
        k = waitforbuttonpress;
    end
        
    
    
    if counter>start_at+1 && plotfigs
        % plot movement
        figure(f_move);
        addpoints(moveplot,movement(counter-start_at-1,1),movement(counter-start_at-1,2));
        
        % plot tracking
        figure(f_track);
        addpoints(trackline,track(counter-start_at-1,1),track(counter-start_at-1,2));
        axis([ROI(1,1) ROI(2,1) ROI(1,2) ROI(2,2)]);
        drawnow limitrate;
        
    end
    % plot smoothed mvt every 100 frames
    plotsmooth = 0;
    if mod(counter,250) == 0 && plotsmooth
        mvt = smooth(movement(:,2),0.01,'rloess');
        figure(9999);
        plot(movement(:,2),'b')
        hold on;
        plot(mvt,'r','linewidth',2);
        pause(3.0);
        close(9999);
    end
    
    if mod(counter,100) == 0
        fprintf('Processed %u frames out of %u.\n',counter,n_frames);
    end
    
end

fprintf('Done running video.\n');

% release objects
release(videoPlayer1);
release(videoPlayer2);
release(videoPlayer3);
close all;

%% SMOOTH MOVEMENT DATA
% Smooth data using 'rloess' method. This uses a local regression using
% weighted linear least squares and a 2nd degree polynomial model. The
% robust version assigns lower weights to outliers in the regression. This
% is helpful to remove the movement spikes due to cables moving while the
% animal is still.

smooth_factor = 0.01;
smooth_chunk_size = 20000;
% smoothed_mvt = smooth(movement(:,2),smooth_factor,'rloess');
smoothed_mvt = smoothChunks(movement(:,2), smooth_factor, smooth_chunk_size);
% display smoothed movement
figure(80);
plot(movement(:,2),'b');
hold on;
plot(smoothed_mvt,'r','linewidth',2);

%% SAVE DATA AND EXIT
% make data structure for saving
outdata.raw_movement                = movement;
outdata.DATA.smooth_movement        = smoothed_mvt;
outdata.DATA.frame_times            = frame_times;
outdata.track                       = track;
outdata.ROI                         = ROI;
outdata.firstframe                  = firstframe;
outdata.intens.meanIntensity        = avgIntensity;
outdata.intens.intThresh            = thresh;
outdata.animal                      = animal;
outdata.nframes                     = n_frames;
outdata.params.SE1                  = SE1;
outdata.params.SE2                  = SE2;
outdata.params.SE3                  = SE3;
outdata.params.SE4                  = SE4;
outdata.params.mask                 = binaryMask;
outdata.params.smoothfactor         = smooth_factor;
outdata.params.alpha                = alpha;
outdata.params.connection_thresh    = connThreshold;
outdata.params.area_thresh          = areaThreshold;
outdata.params.gray_thresh_factor   = thresholdFactor;
outdata.params.n_erode              = n_erode;
outdata.params.n_close              = n_close;
outdata.params.erode_first          = erode_first;

% save in processed video directory
fprintf('\nSaving data...\n');
tic

save([saveloc filesep animal '_NEWmovement.mat'],'outdata');

toc
fprintf('Done!\n\n');


