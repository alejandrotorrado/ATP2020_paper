# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 16:12:35 2016

@author: keithhengen
"""
import numpy as np
import cv2
import os
import operator
import matplotlib.pyplot as plt


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


videoDir = '/Volumes/turrigiano-lab/ATP_VideoFiles/KH63_64/test_ld'
filename = 'KH6364_LDtest_testVid.avi'
#videoDir = '/Volumes/turrigiano-lab/ATP_VideoFiles/KH63_64/testvideos/darkvid'
#filename = 'KH6364test_darkVideo.avi'
vidSrc = cv2.VideoCapture(os.path.join(videoDir,filename))
height, width = int(vidSrc.get(cv2.cv.CV_CAP_PROP_FRAME_HEIGHT)), int(vidSrc.get(cv2.cv.CV_CAP_PROP_FRAME_WIDTH))
half_area = 0.5*height*width
refPt = []
key = []
outside = []
inside = []
counter = 0
while True:

    ret, colorframe = vidSrc.read() # get frame from file
    
    if ret:
        counter+=1
        if counter==1:
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
            
            
            
            
            
    outside.append(np.mean(colorframe[mask==0.]))
    inside.append(np.mean(colorframe[mask==1.]))

    if counter % 100 == 0:
        print 'wanna close it at {}?'.format(counter)
    k = cv2.waitKey(20)
    if k == 27:
        break

cv2.destroyAllWindows()

divided = map(operator.div,outside,inside)
plt.plot(divided)
plt.show()
