release(videoPlayer2);
            release(videoPlayer1);
            release(videoPlayer3);
            release(vidSrc);
            reset(vidSrc);
            close all;

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
    SE1 = strel('line',3,45);
    SE2 = strel('line',5,0);
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
        erode2 = imerode(erode1,SE1);
        close1 = imclose(erode2,SE3);
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
    darkThresh  = 0.15;
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
    ntestframes = 800;
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

% 
% % movement smoothing
% window = 10;
% for ee = 1:size(movement(:,2),1)
%     if ee > window && movement(ee,2) > 0 && ee < size(movement(:,2),1)-window
%         sum_previous = sum(movement(ee-window:ee-1,2));
%         sum_following = sum(movement(ee+1:ee+window,2));
%         if sum_previous == 0 && sum_following == 0
%             movement(ee,2) = 0;
%         end
%     end
% end
            