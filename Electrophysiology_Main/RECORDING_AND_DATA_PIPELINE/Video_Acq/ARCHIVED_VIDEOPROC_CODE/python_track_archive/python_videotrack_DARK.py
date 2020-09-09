# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 15:46:28 2016

@author: atorrado
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 11:59:18 2016

@author: keithhengen
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 14:08:25 2016

@author: keithhengen
"""

import os
import numpy as np
import cv2
import math
import scipy.io as sio
import statsmodels.nonparametric.smoothers_lowess as smoo
import matplotlib.pyplot as plt


# DEFINE HELPER FUNCTIONS

def find_if_close(cnt1,cnt2,thresh):
    # find moments of two contours
    M1 = cv2.moments(cnt1)
    M2 = cv2.moments(cnt2)
    # find centers of each contour from moments
    cx1 = int(M1['m10']/M1['m00'])
    cx2 = int(M2['m10']/M2['m00'])
    cy1 = int(M1['m01']/M1['m00'])
    cy2 = int(M2['m01']/M2['m00'])
    center1 = np.array([[cx1,cy1]])
    center2 = np.array([[cx2,cy2]])
    # compute distance between contour centers
    dist = np.linalg.norm(center1-center2)
    # return true if within threshold distance
    if abs(dist) < thresh :
        return True
    else:
        return False

def selectROI(event, x, y, flags, param):
    global refPt, firstFrame
    # if the left mouse button was clicked, record the starting
    # (x, y) coordinates and indicate that cropping is being
    # performed
    if event == cv2.EVENT_LBUTTONDOWN:
        refPt = [(x, y)]
     
    # check to see if the left mouse button was released
    elif event == cv2.EVENT_LBUTTONUP:
    # record the ending (x, y) coordinates and indicate that
    # the cropping operation is finished
        refPt.append((x, y))
    
    # draw a rectangle around the region of interest
    cv2.rectangle(firstFrame, refPt[0], refPt[1], (0, 255, 0), 2)
    cv2.putText(firstFrame, "Press ESC to continue.", (int(.35*width), height - 30),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255,255,255), 2)
    cv2.imshow('firstframe',firstFrame)
    
    
    
# define path to video file name
videoDir = '/Volumes/turrigiano-lab/ATP_VideoFiles/KH63_64/testvideos/darkvid'
filename = 'KH6364test_darkVideo.avi'
animal = 'KH63'
saveloc = ''

# create video capture object pointing to video
vidSrc = cv2.VideoCapture(os.path.join(videoDir,filename))
# get frame dimensions
height, width = int(vidSrc.get(cv2.cv.CV_CAP_PROP_FRAME_HEIGHT)), int(vidSrc.get(cv2.cv.CV_CAP_PROP_FRAME_WIDTH))
half_area = 0.5*height*width
# define kernels for morphological operations
erode_kern = np.ones((2,2),np.uint8)
close_kern = np.ones((2,2),np.uint8)

# initialize parameters and set thresholds
alpha = 0.2            # blending parameter for adaptive bkg subtraction
areaThresh = 70         # minimum area for blob detection
areaFracThresh = 0.3    # minimum area fraction for blob detection (blobArea/boundingRectArea)
circThresh = 0.1       # minimum circularity for blob detection
linThresh = 2.5        # maximum linearity threshold for blob detection
binThresh = 25        # threshold for binary image thresholding
faroutThresh = 80       # maximum distance from center of largest blob
distThresh = 90        # maximum movement allowed between frames (to prevent "jumps")
alphaThresh = 30        # movement value that triggers change in alpha parameter
cableAreaMin = 50      # minimum sum of areas that turns 0 movement into some movement
minRatArea = 150
artificialMvt = 10.     # when assigning movement values arbitrarily, use this value
binZeroFactor = 2.
nErosions = 1
nDilations = 2
maxBlobsToAnalyze = 10   # maximum number of blobs to use in center of mass analysis
runningAvg = None       # averaged background

waitTime = 40           # time between frame updates
refPt = []              # container for ROI limits

# mode flags
do_distAnalysis = False # perform contour distance algorithm
debugMode = False       # debug mode flag: shows more videos
cut_after = -1
if cut_after < 0:
    cut_after = float('Inf')
start_at = 0
counter = 0           # frame counter

# initialize output arrays
raw_mvt = np.zeros((1,1))
track = np.zeros((1,2))

while True:
    # this while lo
    while counter < start_at:
        if counter % 100 == 0:
            print 'Skipped to frame {}.'.format(counter)
        vidSrc.grab()
        counter += 1
    
    
    ret, colorframe = vidSrc.read() # get frame from file
    
    if ret: # if there is a frame
        counter += 1    # increase frame counter
        firstIsGood = True
        if counter == start_at+1:
            print 'Starting at frame {}.'.format(counter)
            mcx = None
            mcy = None
            # show first frame and set the mouse callback function to selectROI
            firstFrame = colorframe.copy()
            cv2.namedWindow('firstframe')
            cv2.setMouseCallback('firstframe',selectROI)
            while True:
                cv2.putText(firstFrame, "Select ROI corners (top left to bottom right)", (int(.15*width), int(.08*height)),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255,255,255), 2)
                cv2.imshow('firstframe',firstFrame)
                key = cv2.waitKey(0)
                if key == 27:
                    break
            # if user defined an ROI, make a mask based on it
            if len(refPt) == 2:
                roi = firstFrame[refPt[0][1]:refPt[1][1], refPt[0][0]:refPt[1][0]]
                mask = np.zeros(firstFrame.shape[0:2])
                mask[refPt[0][1]:refPt[1][1], refPt[0][0]:refPt[1][0]] = 1.
                mask_half_area = mask[mask==1.].shape[0]/2

            # destroy windows    
            cv2.waitKey(250)
            cv2.destroyAllWindows()
            print 'ROI limits are {}'.format(refPt)
            
        
        # convert frame from BGR to grayscale
        grayframe = cv2.cvtColor(colorframe,cv2.COLOR_BGR2GRAY)
        # take complement image of frame
        invframe = cv2.bitwise_not(grayframe)
#        invframe = grayframe.copy()
        # smooth image to remove noise
        cleanframe = cv2.medianBlur(invframe,5)
        
        # adaptive background subtraction happens here
        if runningAvg is None:
            # if this is the firstframe, make it the background base
            runningAvg = grayframe
        else:
            # compute blended background using alpha parameter
            runningAvg = cv2.addWeighted(cleanframe,alpha,runningAvg,1-alpha,0)     
            
        # subtract background from current frame
        subframe = cv2.absdiff(cleanframe,runningAvg)
        
        # threshold subtracted frame based on light intensity levels
#        _,subframe = cv2.threshold(subframe,binThresh,255,cv2.THRESH_TOZERO)
#        _,subframe = cv2.threshold(subframe,binZeroFactor*binThresh,255,cv2.THRESH_TOZERO_INV)
#        _,subframe = cv2.threshold(subframe,binZeroFactor*binThresh,255,cv2.THRESH_TOZERO_INV)
        _,frame = cv2.threshold(subframe,binThresh,255,cv2.THRESH_BINARY)
        # apply mask to binary image
        frame[mask==0.] = 0
        # smooth binary image to remove noise
        processImg = cv2.medianBlur(frame,5)
        # perform erosion and dilation operations
        processImg = cv2.erode(processImg,erode_kern,iterations = nErosions)
        processImg = cv2.dilate(processImg,close_kern,iterations = nDilations)
        # display results of image processing
        cv2.imshow('After morphological closing',processImg)
        
        
        # At this point the image is processed and ready for blob detection.
        # here is where the dark/light switch should be.
        # if DARK -> do dark_analysis
        # elif LIGHT -> do light_analysis
        # these should be defined as functions. the only argument they need is the image
        # the parameters can be defined within each function?
        # or can pass those too
        # or can define them as globals?        
                
        
        if counter > 0:
            # find contours of white pixel blobs
            contours, hierarchy = cv2.findContours(processImg,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
            # sort contours in descending order by area. Keep only the first 10 contours
            areaM = [cv2.contourArea(i) for i in contours]
            allM = sum(areaM)
            contours = sorted(contours, key = cv2.contourArea, reverse = True)[:10]

            
            # if largest contours is greater than half the unmasked area, it is
            # most likely because of shadows or camera movement resulting in
            # substantial movement in the background. If that is the case, get 
            # rid of largest blob.
            # If not, largest blob is generally the rat.
            if len(contours) > 0:
                tmp_cnt = contours[0]
                if cv2.contourArea(tmp_cnt) > mask_half_area:
                    contours.pop(0)
                    if len(contours) > 0:
                        cnt = contours[0]
                else:
                    cnt = tmp_cnt
                
                
            else:
                cnt = np.array([[[height,width],[height,width]]])
            
            
            # if find flag, do analysis based on distance between blob centers
            if do_distAnalysis:
                LENGTH = len(contours)
                status = np.zeros((LENGTH,1))
                # loop through contours, checking distances
                for i,cnt1 in enumerate(contours):
                    x = i    
                    if i != LENGTH-1:
                        for j,cnt2 in enumerate(contours[i+1:]):
                            x = x+1
                            dist = find_if_close(cnt1,cnt2,30)
                            if dist == True:
                                val = min(status[i],status[x])
                                status[x] = status[i] = val
                            else:
                                if status[x]==status[i]:
                                    status[x] = i+1
                # if found any that are close, unify them and compute a convex hull
                if status.shape[0] > 0:                
                    unified = []
                    maximum = int(status.max())+1
                    for i in xrange(maximum):
                        pos = np.where(status==i)[0]
                        if pos.size != 0:
                            cont = np.vstack(contours[i] for i in pos)
                            hull = cv2.convexHull(cont)
                            unified.append(hull)
                
                    # find the center of the unified convex hull
                    mainhull = unified[0]
                    Mh1 = cv2.moments(mainhull)
                    cx = int(Mh1['m10']/Mh1['m00'])
                    cy = int(Mh1['m01']/Mh1['m00'])
                    
                    # draw results
                    frameclone = colorframe.copy()
                    cv2.drawContours(frameclone,unified,0,(0,0,255),-1)
                    cv2.circle(frameclone, (cx, cy), 5, (40, 120, 200), -1)
                    
                    newframeclone = cv2.addWeighted(colorframe,0.3,frameclone,0.7,0)            
                    
                    cv2.imshow('frameclone',frameclone)
                    cv2.imshow('newframeclone',newframeclone)

            # MAIN BLOB DETECTION ANALYSIS

            # initialize variables and arrays
            center = []
            areas = []
            circularity = []
            linearity = []
            faroutness = []
            w_center = []
            center.append([])
            center.append([])
            w_center.append([])
            w_center.append([])
            # look at first maxBlobsToAnalyze contours maximum (loop through them)
            irange = min(len(contours),maxBlobsToAnalyze)
            for ii in range(irange):
                thiscont = contours[ii]
                # find area, perimeter and circularity of each contour
                area = cv2.contourArea(thiscont)
                perimeter = cv2.arcLength(thiscont,True)
                circ = 4*math.pi*(area/(perimeter**2))
                # compute moments and find center of each contour
                Mii = cv2.moments(thiscont)
                x = int(Mii['m10']/Mii['m00'])
                y = int(Mii['m01']/Mii['m00'])
                center[0].append(x)
                center[1].append(y)
                # find rotated bounding rectangle for blob and get dimensions
                minrect = cv2.minAreaRect(thiscont)
                w_rect = minrect[1][0]
                h_rect = minrect[1][1]
                lin = w_rect/h_rect
                # check linearity threshold
                if lin > linThresh or lin < 1./linThresh:
                    tooLinear = True
#                    print 'too linear!'
                else:
                    tooLinear = False
#                tCorner = tuple([int(i) for i in list(minrect[0])])
#                bCorner = tuple([int(i) for i in list(tuple(np.subtract(topCorner,minrect[1])))])
                
                # find bounding rectangle to check area fraction
                _,_,w_box,h_box = cv2.boundingRect(thiscont)
                box_area = w_box*h_box
                area_fraction = area/box_area
                
                # find distance from largest blob
                # obviously, this will be 0 for ii = 0
                largest_center = np.array([center[0][0],center[1][0]])
                this_center = np.array([center[0][ii],center[1][ii]])
                dist_from_largest = np.linalg.norm(largest_center-this_center)
                
                
                faroutness.append(dist_from_largest)
                circularity.append(circ)
                linearity.append(lin)
                
                
                if circ > circThresh and area_fraction > areaFracThresh and area > areaThresh \
                    and not tooLinear:
                    # if thresholds are passed, compute a weighted coordinate
                    # for each contour. The x and y coordinates are weighted by
                    # the area of each contour
#                    print 'Blob {} passed thresholds:: Circularity was {} and Area was {}.'.format(ii,circ,area)
                    if firstIsGood:
                        if dist_from_largest < faroutThresh:
                            areas.append(area)
                            w_center[0].append(x*area)
                            w_center[1].append(y*area)
                            firstIsGood = True
                    else:
                        if dist_from_largest > faroutThresh:
                            areas.append(area)
                            w_center[0].append(x*area)
                            w_center[1].append(y*area)
                            firstIsGood = True
                else:
#                    print 'Blob {} did not pass thresholds:: Circularity was {} and Area was {}.'.format(ii,circ,area)
                    if ii == 0:
                        firstIsGood = False                    
                    continue
                
                
            # if any blobs are found, compute the coordinates of the center of 
            # mass of the top 5 largest blobs (or less if 5 are not found)
            bigM = 0.
            n_blobs = len(areas)
            if n_blobs>0:
                bigM = sum(areas)
                mcx = int(sum(w_center[0])/bigM)
                mcy = int(sum(w_center[1])/bigM)

                
            if mcx is None:
                mcx = 0
                mcy = 0
            
            if counter > start_at+1:
                print 'mcx : {}'.format(mcx)
                dist_moved = np.linalg.norm(np.array([old_mcx,old_mcy])-np.array([mcx,mcy]))
                if old_mcx != 0:
                    if (dist_moved > distThresh and bigM < minRatArea) or (dist_moved == 0 and allM > cableAreaMin) \
                        or (dist_moved > .5*distThresh and dist_moved < distThresh and allM < cableAreaMin):
                        print 'Artificial movement added. Dist: {}; bigM: {}; allM: {}.'.format(dist_moved,bigM,allM)
                        dist_moved = artificialMvt
                        if mcx != 0:    
                            if counter % 2:
                                mcx = old_mcx + 2
                                mcy = old_mcy + 2
                            else:
                                mcx = old_mcx - 2
                                mcy = old_mcy - 2
                if dist_moved > alphaThresh:
                    alpha = 0.40
                else:
                    alpha = 0.20
                    
                print 'Frame {}. Distance moved: {}.'.format(counter, dist_moved)
                raw_mvt = np.append(raw_mvt,np.array([[dist_moved]]))
                track = np.append(track,np.array([[mcx,mcy]]),axis=0)
                    

            old_mcx = mcx
            old_mcy = mcy
            
            
            
            # Draw images
            contframe = colorframe.copy()
            contframe2 = contframe.copy()
            cblu = (255,173,102)
            cv2.circle(contframe2, (mcx,mcy), 5, cblu, -1)
            cv2.putText(contframe2, "Rat", (mcx - 15, mcy - 15),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.5, cblu, 2)
                        
            # this shows the results of the center-of-mass analysis
            cv2.imshow('Rat Tracker 3: Tokyo Drift',contframe2) 


            # if in debug mode, show more intermediate steps
            if debugMode:
                cv2.imshow('Background subtracted',subframe)
                cv2.imshow('Thresholded image',frame)
                if counter == start_at+1:
                    waitTime = 600
            
            if counter == start_at+cut_after:
                break
    
            k = cv2.waitKey(waitTime)
            if k == 27:
                break


plt.plot(raw_mvt)
plt.show()
smooth_mvt = smoo.lowess(raw_mvt,range(len(raw_mvt)),it=2,frac=0.005,return_sorted = False)
plt.plot(smooth_mvt)
plt.show()
frametimes = np.zeros((counter,1))
nframes = counter
DATA = {'smooth_movement': smooth_mvt, 'frame_times': frametimes, 'nframes': nframes, 'mask': mask}
RAW = {'raw_movement': raw_mvt, 'track': track}

outdata = {'DATA': DATA, 'RAW': RAW}
outfile = animal + '_movement.mat'
sio.savemat(outfile,outdata)


cv2.destroyAllWindows()
vidSrc.release()

