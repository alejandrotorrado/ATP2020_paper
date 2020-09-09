# -*- coding: utf-8 -*-
"""
PYTHON VIDEO TRACKING MAIN SCRIPT



Created on Wed Jul 20 11:03:50 2016

@author: keithhengen (ATP)
"""

# import packages
import sys
sys.path.insert(0, '/Volumes/turrigiano-lab/RECORDING_AND_DATA_PIPELINE/Video_Acq/VideoProcessing/')
import os
import numpy as np
import cv2
import math
import scipy.io as sio
import Tkinter, tkFileDialog
import h5py
from PYTHON_VIDEOTRACK_FUNCTIONS import videoTrack_DARK
from PYTHON_VIDEOTRACK_FUNCTIONS import videoTrack_LIGHT
from PYTHON_VIDEOTRACK_FUNCTIONS import smoothChunks

# Define mouse callback function to select ROI in first frame
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


# File dialog to select AVI file
print('\n\n\t\t___*** Select AVI file in Finder dialog box. ***___\n\n')
root = Tkinter.Tk()
root.withdraw()
root.overrideredirect(True)
root.geometry('0x0+0+0')
root.deiconify()
root.lift()
root.focus_force()
vidFile = tkFileDialog.askopenfilename(parent=root,filetypes=[('AVI files','.avi')],initialdir="/",
            title='Choose the AVI video file.')
root.destroy()


animal = raw_input('Animal name (e.g. AT08) ?  ')
recname = raw_input('Recording name (e.g. AT07_08) ?  ')
 

# create video capture object pointing to video
vidSrc = cv2.VideoCapture(vidFile)
# get frame dimensions
height, width = int(vidSrc.get(cv2.cv.CV_CAP_PROP_FRAME_HEIGHT)), int(vidSrc.get(cv2.cv.CV_CAP_PROP_FRAME_WIDTH))
#height, width = int(vidSrc.get(cv2.CAP_PROP_FRAME_HEIGHT)), int(vidSrc.get(cv2.CAP_PROP_FRAME_WIDTH)) # if using opencv3
half_area = 0.5*height*width
nFrames = int(math.ceil(vidSrc.get(cv2.cv.CV_CAP_PROP_FRAME_COUNT)))
# nFrames = int(math.ceil(vidSrc.get(cv2.CAP_PROP_FRAME_COUNT))) # if using opencv3
# define kernels for morphological operations
erode_kern = np.ones((2,2),np.uint8)
close_kern = np.ones((2,2),np.uint8)


# SET PARAMETERS FOR LIGHT FRAMES
paramsLight = {}                        # initialize parameter dictionary
paramsLight['alpha'] = 0.05             # blending parameter for adaptive bkg subtraction
paramsLight['areaThresh'] = 30          # minimum area for blob detection
paramsLight['circThresh'] = 0.10       # minimum circularity for blob detection
paramsLight['linThresh'] = 3.0          # maximum linearity threshold for blob detection
paramsLight['binThresh'] = 50           # threshold for binary image thresholding
paramsLight['nErosions'] = 1             # number of erosion operations to perform
paramsLight['nDilations'] = 2            # number of dilation operations to perform
paramsLight['maxBlobsToAnalyze'] = 8    # maximum number of blobs to use in center of mass analysis

# SET PARAMETERS FOR DARK FRAMES
paramsDark = {}                         # initialize parameter dictionary
paramsDark['alpha'] = 0.2               # blending parameter for adaptive bkg subtraction
paramsDark['areaThresh'] = 70           # minimum area for blob detection
paramsDark['areaFracThresh'] = 0.3      # minimum area fraction for blob detection (blobArea/boundingRectArea)
paramsDark['circThresh'] = 0.1          # minimum circularity for blob detection
paramsDark['linThresh'] = 2.5           # maximum linearity threshold for blob detection
paramsDark['binThresh'] = 25            # threshold for binary image thresholding
paramsDark['faroutThresh'] = 80         # maximum distance from center of largest valid blob
paramsDark['distThresh'] = 90           # maximum movement allowed between frames (to prevent "jumps")
paramsDark['alphaThresh'] = 30          # movement value that triggers change in alpha parameter
paramsDark['cableAreaMin'] = 50         # minimum sum of areas that turns 0 movement into some movement
paramsDark['minRatArea'] = 150          # if movement but area is smaller than this, set to artificial movement
paramsDark['artificialMvt'] = 10.       # when assigning movement values arbitrarily, use this value
paramsDark['nErosions'] = 1             # number of erosion operations to perform
paramsDark['nDilations'] = 2            # number of dilation operations to perform
paramsDark['maxBlobsToAnalyze'] = 10    # maximum number of blobs to use in analysis

# General parameters
initParams = {}                 # dictionary for other general parameters, including those initialized on frame 1
runningAvg = None               # averaged background
waitTime = 2                    # minimum time between frame updates
refPt = []                      # container for ROI limits
cut_after = -1                  # stop after this amount of frames. Set to -1 to keep going until end of video.
if cut_after < 0:               
    cut_after = float('Inf')
start_at = 0                    # starting frame
counter = 0                     # frame counter

showVids = False        # display flag: show basic videos (tracking + morph)
debugMode = False       # debug mode flag: shows more videos and increases waitTime


# initialize output arrays
raw_mvt = np.zeros((1,1))
track = np.zeros((1,2))



while counter < nFrames:
    
    # this while loop skips to the chosen start frame
    while counter < start_at:
        if counter % 100 == 0:
            print('Skipped to frame {}.'.format(counter))
        vidSrc.grab()
        counter += 1
        
    if counter % 10000 == 0:
        print('Processing frame {} out of {}.'.format(counter,nFrames))
        
    ret, colorframe = vidSrc.read() # get frame from video object
    if ret: # if there is a frame
        counter += 1    # increase frame counter
        initParams['firstIsGood'] = True  # initialize first blob quality flag
        if counter == start_at+1: # if this is the first frame, select the ROI
            print('Starting at frame {}.'.format(counter))
            initParams['mcx'] = 0      # initialize center of mass coordinates as None
            initParams['mcy'] = 0
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
                initParams['mask_half_area'] = mask[mask==1.].shape[0]/2

            # destroy windows    
            cv2.waitKey(250)
            cv2.destroyAllWindows()
            print('ROI limits are {}'.format(refPt))
        
                
        intensFrac = np.mean(colorframe[mask==0.])/np.mean(colorframe[mask==1.])
        if intensFrac > 0.7:
            frameIsDark = False
        else:
            frameIsDark = True
            
        if frameIsDark:
            params = paramsDark
        else:
            params = paramsLight
        
        
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
            alpha = params['alpha']
            runningAvg = cv2.addWeighted(cleanframe,alpha,runningAvg,1-alpha,0)     
            
        # subtract background from current frame
        subframe = cv2.absdiff(cleanframe,runningAvg)
        
        # threshold subtracted frame based on light intensity levels
        binThresh = params['binThresh']
        _,frame = cv2.threshold(subframe,binThresh,255,cv2.THRESH_BINARY)
        # apply mask to binary image
        frame[mask==0.] = 0
        # smooth binary image to remove noise
        processImg = cv2.medianBlur(frame,5)
        # perform erosion and dilation operations
        nErosions = params['nErosions']
        nDilations = params['nDilations']
        processImg = cv2.erode(processImg,erode_kern,iterations = nErosions)
        processImg = cv2.dilate(processImg,close_kern,iterations = nDilations)
        
        # display results of image processing
        if showVids:
            cv2.imshow('After morphological closing',processImg)
            
        if frameIsDark:
            raw_mvt,track,new_mc,old_mc,firstIsGood,newalpha = videoTrack_DARK(processImg,params,initParams,raw_mvt,track,counter,start_at,verbose=False)
            params['alpha'] = newalpha
            initParams['firstIsGood'] = firstIsGood
        else:
            raw_mvt,track,new_mc,old_mc = videoTrack_LIGHT(processImg,params,initParams,raw_mvt,track,counter,start_at,verbose=False)
        
        # update CoM and other params
        initParams['old_mcx'] = old_mc[0]
        initParams['old_mcy'] = old_mc[1]
        initParams['mcx'] = new_mc[0]
        initParams['mcy'] = new_mc[1]
        
            
# Draw images
            # if in debug mode, show more intermediate steps
        if debugMode:
            showVids = True
            cv2.imshow('Background subtracted',subframe)
            cv2.imshow('Thresholded image',frame)
            if counter == start_at+1:
                waitTime = 600
                
        if showVids and counter > start_at+1:
            contframe = colorframe.copy()
            contframe2 = contframe.copy()
            cblu = (255,173,102)
            cv2.circle(contframe2, (initParams['mcx'],initParams['mcy']), 5, cblu, -1)
            cv2.putText(contframe2, "Rat", (initParams['mcx'] - 15, initParams['mcy'] - 15),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.5, cblu, 2)
                        
            # this shows the results of the center-of-mass analysis
            cv2.imshow('Rat Tracker 3: Tokyo Drift',contframe2) 
            
        if counter == start_at+cut_after:
            break

        k = cv2.waitKey(waitTime)
        if k == 27:
            break

# File dialog to select frametimes MAT file
videoDir = vidFile.rsplit(os.path.sep,1)[0]

ftFilename = videoDir + os.path.sep + recname + '_frameTimes.mat'
ftFile = sio.loadmat(ftFilename)
frametimes = ftFile['frametimes']

chSize = 40000
smooth_mvt = smoothChunks(raw_mvt,chSize)

print('Saving your data...')
# package data into struct and 
DATA = {'smooth_movement': smooth_mvt, 'frame_times': frametimes, 'nframes': nFrames, 'mask': mask}
RAW = {'raw_movement': raw_mvt, 'track': track}


outdata = {'DATA': DATA, 'RAW': RAW}
outfile = videoDir + os.path.sep + animal + '_pymovement.mat'
sio.savemat(outfile,outdata)
print('Done!')

cv2.destroyAllWindows()
vidSrc.release()

