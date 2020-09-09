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

# create video capture object pointing to video
vidSrc = cv2.VideoCapture(os.path.join(videoDir,filename))
# get frame dimensions
height, width = int(vidSrc.get(cv2.cv.CV_CAP_PROP_FRAME_HEIGHT)), int(vidSrc.get(cv2.cv.CV_CAP_PROP_FRAME_WIDTH))
half_area = 0.5*height*width
# define kernels for morphological operations
erode_kern = np.ones((2,2),np.uint8)
close_kern = np.ones((2,2),np.uint8)

# initialize parameters and set thresholds
alpha = 0.05            # blending parameter for adaptive bkg subtraction
areaThresh = 30         # minimum area for blob detection
circThresh = 0.10       # minimum circularity for blob detection
linThresh = 3.0         # maximum linearity threshold for blob detection
maxBlobsToAnalyze = 8   # maximum number of blobs to use in center of mass analysis
runningAvg = None       # averaged background
counter  = 0            # frame counter
refPt = []              # container for ROI limits

# mode flags
do_distAnalysis = False # perform contour distance algorithm
debugMode = False       # debug mode flag: shows more videos


while True:
    ret, colorframe = vidSrc.read() # get frame from file
    
    if ret: # if there is a frame
        counter += 1    # increase frame counter
        if counter == 1:
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
        _,frame = cv2.threshold(subframe,50,255,cv2.THRESH_BINARY)
        # apply mask to binary image
        frame[mask==0.] = 0
        # smooth binary image to remove noise
        processImg = cv2.medianBlur(frame,5)
        # perform erosion and dilation operations
        processImg = cv2.erode(processImg,erode_kern,iterations = 1)
        processImg = cv2.dilate(processImg,close_kern,iterations = 2)
        # display results of image processing
        cv2.imshow('After morphological closing',processImg)
        
        
        # At this point the image is processed and ready for blob detection.
        
        if counter > 0:
            # find contours of white pixel blobs
            contours, hierarchy = cv2.findContours(processImg,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
            # sort contours in descending order by area. Keep only the first 10 contours
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
            w_center = []
            center.append([])
            center.append([])
            w_center.append([])
            w_center.append([])
            # look at first 5 contours maximum (loop through them)
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
#                
                
                circularity.append(circ)
                linearity.append(lin)
                if circ > circThresh and area > areaThresh and not tooLinear:
                    # if thresholds are passed, compute a weighted coordinate
                    # for each contour. The x and y coordinates are weighted by
                    # the area of each contour
#                    print 'Blob {} passed thresholds:: Circularity was {} and Area was {}.'.format(ii,circ,area)
                    areas.append(area)
                    w_center[0].append(x*area)
                    w_center[1].append(y*area)
                else:
                    continue
#                    print 'Blob {} did not pass thresholds:: Circularity was {} and Area was {}.'.format(ii,circ,area)
            
            # if any blobs are found, compute the coordinates of the center of 
            # mass of the top 5 largest blobs (or less if 5 are not found)
            n_blobs = len(areas)
            if n_blobs>0:
                bigM = sum(areas)
                mcx = int(sum(w_center[0])/bigM)
                mcy = int(sum(w_center[1])/bigM)

            # ALTERNATIVE ANALYSIS
            # Find convex hull for largest blob and compute center of that. This
            # relies on the assumption that the largest blob found will be the 
            # one cooresponding to the rat. Not always a good assumption!
            M = cv2.moments(cnt)
            if M['m00'] > 0:
                hull = cv2.convexHull(cnt)
                Mh = cv2.moments(hull)
                cx = int(Mh['m10']/Mh['m00'])
                cy = int(Mh['m01']/Mh['m00'])
                 
            
            # Draw images
            contframe = colorframe.copy()
            contframe2 = contframe.copy()
            cblu = (255,173,102)
            cv2.circle(contframe2, (mcx,mcy), 5, cblu, -1)
            cv2.putText(contframe2, "Rat", (mcx - 15, mcy - 15),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.5, cblu, 2)
                        
            # this shows the results of the center-of-mass analysis
            cv2.imshow('Rat Tracker 3: Tokyo Drift',contframe2) 
            cv2.drawContours(contframe,[hull],0,(255,0,0),-1)
            cv2.drawContours(contframe,[cnt],0,(0,0,255),-1)
            cv2.circle(contframe, (cx, cy), 5, (255, 255, 255), -1)
            # this shows the results of the convex hull analysis
            # (uncomment below if needed)
#            cv2.imshow('contframe',contframe)
            
            # if in debug mode, show more intermediate steps
            if debugMode:
                cv2.imshow('Background subtracted',subframe)
                cv2.imshow('Thresholded image',frame)
             
    
            k = cv2.waitKey(20)
            if k == 27:
                break


cv2.destroyAllWindows()
vidSrc.release()

