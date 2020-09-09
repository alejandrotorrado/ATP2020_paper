# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 14:08:25 2016

@author: keithhengen
"""

import os
import numpy as np
import cv2


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
    cv2.imshow('firstframe',firstFrame)
    
    
    
    
videoDir = '/Volumes/turrigiano-lab/ATP_VideoFiles/KH63_64/full'
filename = 'KH63_64_fullVideo.avi'

vidSrc = cv2.VideoCapture(os.path.join(videoDir,filename))
height, width = int(vidSrc.get(cv2.cv.CV_CAP_PROP_FRAME_HEIGHT)), int(vidSrc.get(cv2.cv.CV_CAP_PROP_FRAME_WIDTH))
half_area = 0.5*height*width
maxval = 12
erode_kern = np.ones((2,2),np.uint8)
close_kern = np.ones((2,2),np.uint8)



runningAvg = None
alpha = 0.01
areaThresh = 75

counter  = 0
refPt = []

#bg2 = cv2.BackgroundSubtractorMOG2(10,maxval*maxval)
while True:
    ret, colorframe = vidSrc.read() # get frame from file
    
    if ret: # if there is a frame
        counter += 1
        if counter == 1:
            firstFrame = colorframe
            cv2.namedWindow('firstframe')
            cv2.setMouseCallback('firstframe',selectROI)
            while True:
                cv2.imshow('firstframe',firstFrame)
                key = cv2.waitKey(0)
                if key == 27:
                    break
            if len(refPt) == 2:
                roi = firstFrame[refPt[0][1]:refPt[1][1], refPt[0][0]:refPt[1][0]]
                cv2.imshow("ROI", roi)
            cv2.waitKey(0)
            cv2.destroyAllWindows()
            print 'ROI limits are {}'.format(refPt)
            
        
        # convert it to grayscale
        grayframe = cv2.cvtColor(colorframe,cv2.COLOR_BGR2GRAY) # convert RGB to grayscale
        invframe = cv2.bitwise_not(grayframe)
        cleanframe = cv2.medianBlur(invframe,5)
        
        if runningAvg is None:
            runningAvg = grayframe
        else:
            runningAvg = cv2.addWeighted(cleanframe,alpha,runningAvg,1-alpha,0)     

        subframe = cv2.absdiff(cleanframe,runningAvg)
  
        
        
#        conncomps = cv2.connectedComponents(cleanframe)

        _,fullframe = cv2.threshold(subframe,50,255,cv2.THRESH_BINARY)
#        _,thframe = cv2.threshold(frame,200,255,cv2.THRESH_BINARY)
        frame = fullframe
        frame[0:refPt[0][1],:] = 0
        frame[:,0:refPt[0][0]] = 0
        frame[refPt[1][1]:,:] = 0
        frame[:,refPt[1][0]:] = 0
        processImg = frame
        processImg = cv2.medianBlur(processImg,5)
        processImg = cv2.erode(processImg,erode_kern,iterations = 1)
        processImg = cv2.dilate(processImg,close_kern,iterations = 2)
        
       
        
        
        cv2.imshow('frame',frame)
        cv2.imshow('processImg',processImg)
        
        
        # At this point the image is processed and ready for blob detection.
        

        if counter > 0:
            base_contours, hierarchy = cv2.findContours(processImg,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
            base_contours = sorted(base_contours, key = cv2.contourArea, reverse = True)[:20]
            
            print 'Base contour array has {} elements.'.format(len(base_contours))             
            
            savethese = []
            for xx in range(len(base_contours)):
                print xx
                if cv2.contourArea(base_contours[xx]) >= areaThresh:
                    savethese.append(xx)
            
            contours = [base_contours[i] for i in savethese]
            
            print 'New contour array has {} elements.'.format(len(contours))            
            
            if len(contours) > 0:
                tmp_cnt = contours[0]
                if cv2.contourArea(tmp_cnt) > half_area:
                    contours.pop(0)
                    if len(contours) > 0:
                        cnt = contours[0]
                else:
                    cnt = tmp_cnt
                
                
            else:
                cnt = np.array([[[height,width],[height,width]]])
            
            
            center = []
            areas = []
            w_center = []
            center.append([])
            center.append([])
            w_center.append([])
            w_center.append([])
            irange = min(len(contours),3)
            for ii in range(irange):
#                print ii
                Mii = cv2.moments(contours[ii])
                x = int(Mii['m10']/Mii['m00'])
                y = int(Mii['m01']/Mii['m00'])
                area = cv2.contourArea(contours[ii])
                center[0].append(x)
                center[1].append(y)
                areas.append(area)
                w_center[0].append(x*area)
                w_center[1].append(y*area)
            
            n_blobs = len(center[0])
            if n_blobs>0:
                bigM = sum(areas)
                mcx = int(sum(w_center[0])/bigM)
                mcy = int(sum(w_center[1])/bigM)
#                mcx = sum(center[0][:])/n_blobs
#                mcy = sum(center[1][:])/n_blobs
                
            M = cv2.moments(cnt)
            if M['m00'] > 0:
                hull = cv2.convexHull(cnt)
                Mh = cv2.moments(hull)
                cx = int(Mh['m10']/Mh['m00'])
                cy = int(Mh['m01']/Mh['m00'])
                 
            
            contframe = cv2.cvtColor(grayframe,cv2.COLOR_GRAY2RGB)
            contframe2 = contframe
            cv2.circle(contframe2, (mcx,mcy), 5, (0,255,0), -1)
            cv2.putText(contframe2, "Rat", (mcx - 15, mcy - 15),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0,255,0), 2)
  
            cv2.imshow('contframe2',contframe2)
            cv2.drawContours(contframe,[hull],0,(255,0,0),-1)
            cv2.drawContours(contframe,[cnt],0,(0,0,255),-1)
            cv2.circle(contframe, (cx, cy), 5, (255, 255, 255), -1)
            
            
            
    #        cv2.imshow('thframe',thframe)
#            cv2.imshow('cleanframe',cleanframe)
            cv2.imshow('contframe',contframe)
            
#            cv2.imshow('subframe',subframe)
#            cv2.imshow('runningAvg',runningAvg)
            
    
    
            k = cv2.waitKey(20)
            if k == 27:
                break


cv2.destroyAllWindows()
vidSrc.release()

