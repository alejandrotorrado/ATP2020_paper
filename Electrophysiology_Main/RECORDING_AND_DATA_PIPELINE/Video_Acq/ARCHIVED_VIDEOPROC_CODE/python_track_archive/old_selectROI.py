# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 16:32:57 2016

@author: keithhengen
"""

import cv2
 

 
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
		cv2.rectangle(image, refPt[0], refPt[1], (0, 255, 0), 2)
		cv2.imshow("cropped", firstframe)
  
  