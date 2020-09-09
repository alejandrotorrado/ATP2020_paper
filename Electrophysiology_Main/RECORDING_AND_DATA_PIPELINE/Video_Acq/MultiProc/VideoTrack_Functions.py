# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 12:33:16 2016

@author: keithhengen
"""

import cv2
import numpy as np
import math
import statsmodels.nonparametric.smoothers_lowess as smoo
import pdb

def videoTrack_DARK(vidframe,vidparams,genparams,raw_mvt,track,framecounter,verbose=False):
    
    
    # find contours of white pixel blobs
    (_, contours, _) = cv2.findContours(vidframe,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
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
        if cv2.contourArea(tmp_cnt) > genparams['mask_half_area']:
            contours.pop(0)
    
    
    # MAIN BLOB DETECTION ANALYSIS
    
    # initialize variables and arrays
    center = []
    areas = []
    w_center = []
    center.append([])
    center.append([])
    w_center.append([])
    w_center.append([])
    firstIsGood = genparams['firstIsGood']
    # look at first maxBlobsToAnalyze contours maximum (loop through them)
    irange = min(len(contours),vidparams['maxBlobsToAnalyze'])
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
        if lin > vidparams['linThresh'] or lin < 1./vidparams['linThresh']:
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
        
        
        if circ > vidparams['circThresh'] \
            and area_fraction > vidparams['areaFracThresh'] \
            and area > vidparams['areaThresh'] \
            and not tooLinear:
            # if thresholds are passed, compute a weighted coordinate
            # for each contour. The x and y coordinates are weighted by
            # the area of each contour
            if firstIsGood:
                if dist_from_largest < vidparams['faroutThresh']:
                    areas.append(area)
                    w_center[0].append(x*area)
                    w_center[1].append(y*area)
                    firstIsGood = True
            else:
                if dist_from_largest > vidparams['faroutThresh']:
                    areas.append(area)
                    w_center[0].append(x*area)
                    w_center[1].append(y*area)
                    firstIsGood = True
        else:
            if ii == 0:
                firstIsGood = False                    
            continue
    
    # retrieve center of mass coordinates from params
    mcx = genparams['mcx']
    mcy = genparams['mcy']
    
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
    
    if framecounter > 1:
        # retrieve previous CoM coords from params (only after 1st frame)
        old_mcx = genparams['old_mcx']
        old_mcy = genparams['old_mcy']
        
        dist_moved = np.linalg.norm(np.array([old_mcx,old_mcy])-np.array([mcx,mcy]))
        if old_mcx != 0:
            if (dist_moved > vidparams['distThresh'] and bigM < vidparams['minRatArea']) \
            or (dist_moved == 0 and allM > vidparams['cableAreaMin']) \
            or (dist_moved > .5*vidparams['distThresh'] \
            and dist_moved < vidparams['distThresh'] and allM < vidparams['cableAreaMin']):
                if verbose:
                    print('Artificial movement added. Dist: {}; bigM: {}; allM: {}.'.format(dist_moved,bigM,allM))
                dist_moved = vidparams['artificialMvt']
                if mcx != 0:    
                    if framecounter % 2:
                        mcx = old_mcx + 2
                        mcy = old_mcy + 2
                    else:
                        mcx = old_mcx - 2
                        mcy = old_mcy - 2
        if dist_moved > vidparams['alphaThresh']:
            newalpha = 0.40
        else:
            newalpha = 0.20
            
        if verbose:
            print('Frame {}. Distance moved: {}.'.format(framecounter, dist_moved))
        raw_mvt = np.append(raw_mvt,np.array([[dist_moved]]))
        track = np.append(track,np.array([[mcx,mcy]]),axis=0)
    else:
        dist_moved = np.array([[0]])
        track_mc = np.array([[0,0]])        
    
    new_mc = (mcx,mcy)
    old_mc = (mcx,mcy)
    
#    return raw_mvt,track,new_mc,old_mc,firstIsGood,newalpha
    return dist_moved,track_mc,new_mc,old_mc,firstIsGood,newalpha
    
    
    
    
    
    
def videoTrack_LIGHT(vidframe,vidparams,genparams,raw_mvt,track,framecounter,verbose=False):
    
    # find contours of white pixel blobs
    (_, contours, _) = cv2.findContours(vidframe,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
    # sort contours in descending order by area. Keep only the first 10 contours
    contours = sorted(contours, key = cv2.contourArea, reverse = True)[:10]
    
    # if largest contours is greater than half the unmasked area, it is
    # most likely because of shadows or camera movement resulting in
    # substantial movement in the background. If that is the case, get 
    # rid of largest blob.
    # If not, largest blob is generally the rat.
    if len(contours) > 0:
        tmp_cnt = contours[0]
        if cv2.contourArea(tmp_cnt) > genparams['mask_half_area']:
            contours.pop(0)
    
    
    # initialize variables and arrays
    center = []
    areas = []
    w_center = []
    center.append([])
    center.append([])
    w_center.append([])
    w_center.append([])
    # look at first 5 contours maximum (loop through them)
    irange = min(len(contours),vidparams['maxBlobsToAnalyze'])
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
        if lin > vidparams['linThresh'] or lin < 1./vidparams['linThresh']:
            tooLinear = True
        else:
            tooLinear = False
             
        
        if circ > vidparams['circThresh'] and area > vidparams['areaThresh'] and not tooLinear:
            # if thresholds are passed, compute a weighted coordinate
            # for each contour. The x and y coordinates are weighted by
            # the area of each contour
            areas.append(area)
            w_center[0].append(x*area)
            w_center[1].append(y*area)
        else:
            continue

    mcx = genparams['mcx']
    mcy = genparams['mcy']
    
    # if any blobs are found, compute the coordinates of the center of 
    # mass of the top 5 largest blobs (or less if 5 are not found)
    n_blobs = len(areas)
    if n_blobs>0:
        bigM = sum(areas)
        mcx = int(sum(w_center[0])/bigM)
        mcy = int(sum(w_center[1])/bigM)
        
#    if mcx is None:
#        mcx = 0
#        mcy = 0
    
    if framecounter > 1:
        # retrieve previous CoM coords from params (only after 1st frame)
        old_mcx = genparams['old_mcx']
        old_mcy = genparams['old_mcy']
        
        dist_moved = np.linalg.norm(np.array([old_mcx,old_mcy])-np.array([mcx,mcy]))
        if verbose:
            print('Frame {}. Distance moved: {}.'.format(framecounter, dist_moved))
#        raw_mvt = np.append(raw_mvt,np.array([[dist_moved]]))
        dist_moved = np.array([[dist_moved]])
#        track = np.append(track,np.array([[mcx,mcy]]),axis=0)
        track_mc = np.array([[mcx,mcy]])
    else:
        dist_moved = np.array([[0]])
        track_mc = np.array([[0,0]])
            
            
    new_mc = (mcx,mcy)
    old_mc = (mcx,mcy)
    
#    return raw_mvt,track,new_mc,old_mc
    return dist_moved,track_mc,new_mc,old_mc
    
    
def smoothChunks(mvt,chunkSize):
    
    smoothed = np.zeros(mvt.shape)
    
    nPts = float(mvt.shape[0])
    nChunks = int(math.ceil(nPts/chunkSize))
    
    for chunk in range(nChunks):
        print('Smoothing chunk {} of {}.'.format(chunk,nChunks))
        start_pos = (chunkSize * chunk)
        end_pos = chunkSize * (chunk + 1)
        if end_pos > nPts:
            end_pos = nPts
            
        mvtChunk = mvt[start_pos:end_pos]
        
        smoothChunk = smoo.lowess(mvtChunk,range(len(mvtChunk)),it=2,frac=0.005,return_sorted = False)
        smoothChunk[smoothChunk<0.] = 0.
        
        smoothed[start_pos:end_pos] = smoothChunk
        
    return smoothed
        