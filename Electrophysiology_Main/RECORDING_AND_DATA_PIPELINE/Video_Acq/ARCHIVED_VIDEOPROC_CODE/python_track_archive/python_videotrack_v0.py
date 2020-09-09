# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 14:08:25 2016

@author: keithhengen
"""

import os
import numpy as np
import cv2

global refPt, firstframe

def selectROI(event, x, y, flags, params):
	# grab references to the global variables
	global refPt, firstframe
 
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
		cv2.rectangle(firstframe, refPt[0], refPt[1], (0, 255, 0), 2)
#		cv2.imshow("cropped", firstframe)

videoDir = '/Users/keithhengen/Desktop/animal_video_frames/KH67_68'
filename = 'KH67_68_200.avi'

vidSrc = cv2.VideoCapture(os.path.join(videoDir,filename))
height, width = int(vidSrc.get(cv2.cv.CV_CAP_PROP_FRAME_HEIGHT)), int(vidSrc.get(cv2.cv.CV_CAP_PROP_FRAME_WIDTH))
half_area = 0.5*height*width

refPt = []

# construct the argument parser and parse the arguments
ret,firstframe = vidSrc.read()
if ret:
    clone = firstframe.copy()
    cv2.namedWindow("firstframe")
    cv2.setMouseCallback("firstframe",selectROI)
 
    # keep looping until the 'Esc' key is pressed
    while True:
        # display the image and wait for a keypress
        cv2.imshow("firstframe", firstframe)
        key = cv2.waitKey(10) & 0xFF
         
        # if the 'r' key is pressed, reset the cropping region
        if key == ord("r"):
            firstframe = clone.copy()
         
        # if the 'c' key is pressed, break from the loop
        elif key == ord('q'):
            break
            cv2.destroyAllWindows()
 
# if there are two reference points, then crop the region of interest
# from teh image and display it
        if len(refPt) == 2:
            roi = clone[refPt[0][1]:refPt[1][1], refPt[0][0]:refPt[1][0]]
            mask = np.zeros(firstframe.shape[0:2])
            mask[refPt[0][1]:refPt[1][1], refPt[0][0]:refPt[1][0]] = 1.
            cv2.imshow("ROI", roi)

cv2.destroyAllWindows()
vidSrc.release()
vidSrc = cv2.VideoCapture(os.path.join(videoDir,filename))

maxval = 10
erode_kern = np.ones((6,4),np.uint8)
close_kern = np.ones((5,8),np.uint8)
bg2 = cv2.BackgroundSubtractorMOG2(150,maxval*maxval)
while True:
    ret, colorframe = vidSrc.read() # get frame from file
    if ret: # if there is a frame
        # convert it to grayscale
        frame = cv2.cvtColor(colorframe,cv2.COLOR_BGR2GRAY) # convert RGB to grayscale
        clean_img = cv2.medianBlur(frame,11)
        foreground = bg2.apply(clean_img)
        _,fgThresh = cv2.threshold(foreground,127,255,cv2.THRESH_BINARY)
        processImg = fgThresh
        processImg[mask==0.] = 0
        processImg = cv2.erode(processImg,erode_kern,iterations=1)
        processImg = cv2.dilate(processImg,close_kern,iterations=3)


        fgBlur = cv2.GaussianBlur(processImg,(11,11),0)
        adaptive = cv2.adaptiveThreshold(fgBlur,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,\
            cv2.THRESH_BINARY,3,1)
            
        contours,hierarchy = cv2.findContours(adaptive,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
        contours = sorted(contours, key = cv2.contourArea, reverse = True)[:10]
        tmp_cnt = contours[0]
        if cv2.contourArea(tmp_cnt) > half_area:
            contours.pop(0)
        if len(contours) > 0:
            cnt = contours[0]
        else:
            cnt = np.array([[[height,width],[height,width]]])
        reframe = cv2.cvtColor(frame,cv2.COLOR_GRAY2RGB) 
        cv2.drawContours(reframe,[cnt],0,(0,205,212),-1)
        
        cv2.imshow('grayscale',frame)
        cv2.imshow('fg',foreground)
#        cv2.imshow('thresh1',fgThresh)
#        cv2.imshow('binaryImg',binaryImg)
        cv2.imshow('process',processImg)
#        cv2.imshow('adapt',adaptive)
        cv2.imshow('contour',reframe)


        k = cv2.waitKey(20)
        if k == 27:
            break


cv2.destroyAllWindows()
vidSrc.release()