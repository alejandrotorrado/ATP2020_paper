#! /usr/bin/env python

#import Tkinter, tkFileDialog
import cv2
import time
import numpy as np
import csv
import os
import sys
import zipfile, zlib
import datetime
import shutil
import scipy.io as sio
import multiprocessing as mpc

# _____ set up multiprocessing parameters
# first find number of physical cores available
hwlist = os.popen('sysctl hw').readlines()[1:20]
phys_str = [s for s in hwlist if "hw.physicalcpu:" in s]
if len(phys_str) == 1:
	phys_str = phys_str[0]
else:
	print('Strange... could not determine number of cores')
nCores = int(list(filter(str.isdigit,phys_str))[0])

# set up manager and queue objects
m = mpc.Manager()
q = m.JoinableQueue()
q_sig = m.JoinableQueue()

# create lock object to ensure prints are thread-safe
print_lock = mpc.Lock()

# ___ Define functions.
# I will need:
#    1) a function that does the tracking. Takes one frame (from q?) and 
# 		applies videotracking based on previous frames to it. This will
#		basically be the meat of the videotrack script.
# 	 2) a function that does the acquisition. Probably best to define the
#		capture object first and pass it as an argument. Then this function
# 		needs to get a frame from it at a given framerate and dump it in 
# 		the queue.
# 	 3) Figure out what I need for saving.

def acqFrames(queue,sigq,framerate,framecount,verbose=False):
	with print_lock:
		print('Initializing camera.')
	cap = cv2.VideoCapture(0)
	while True:
		tic1 = time.time()
		ret, rawframe = cap.read()
		if ret:
			# update framecounter
			framecount += 1
			if verbose:
				with print_lock:
					print('ACQ:: Frame {} acquired at {}.'.format(framecount,datetime.datetime.fromtimestamp(tic1).strftime('%Y-%m-%d %H:%M:%S.%f')))
			rawframe = cv2.resize(rawframe,(0,0),fx=0.25,fy=0.25)
			# put frame in queue for processing
			q.put(rawframe)
			# show frame
			cv2.imshow('acq frame',rawframe)
			# for accurate timing, measure time required to take frame
			toc1 = time.time()
			# pause for long enough to get right framerate
			paws = int( round( 1000* (1./framerate - (toc1-tic1) )))
			if verbose:
				with print_lock:
					print('ACQ:: Acquiring frame {} took {} seconds.\n'.format(framecount,toc1-tic1))
					print('paws was {}'.format(paws))
			# wait for button press
			k = cv2.waitKey(paws)
			# pressing Esc key causes infinite loop to break
			if k == 27:
				with print_lock:
					print('ACQ:: Halting acquisition on framecount {}.\n'.format(framecount))
				cv2.destroyAllWindows()
				cap.release()
				sigq.put(0)
				break
		else:
			with print_lock:
				print('ACQ:: No frame was acquired on framecount {}!\n'.format(framecount))
			continue

def vidTrack(frame):

	# tracking
	(_, contours, _) = cv2.findContours(frame,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
	contours = sorted(contours, key = cv2.contourArea, reverse = True)[:10]
	if len(contours) > 0:
		thisone = contours[0]
		area = cv2.contourArea(thisone)
	else:
		area = 0.
	return area


def runTrackProcess(queue,sig_queue,save_int,savefile):
	with print_lock:
		print('I am process ', mpc.current_process().name , ' and I just started.')
	# as long as there are frames in the queue
	t01 = time.time()
	runningAvg = None
	saved_areas = np.array([])

	while True:
		if queue.empty():
			if not sig_queue.empty():
				signal = sig_queue.get()
				if signal == 0:
					print('Found QUIT signal. Aborting tracking.')
					break
			time.sleep(0.1)
			continue

		print('Looking for frame in queue.')
		colorframe = queue.get() 		# grab a frame from queue

		# if the queue had an object but it is a "None" type, abort tracking
		if colorframe is None:
			with print_lock:
				print('Camera error. Aborting tracking...')
				queue.task_done()
				break

		print('Found frame! Processing...')
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
			alpha = 0.1
			runningAvg = cv2.addWeighted(cleanframe,alpha,runningAvg,1-alpha,0)     
		# subtract background from current frame
		subframe = cv2.absdiff(cleanframe,runningAvg)
		# threshold subtracted frame based on light intensity levels
		binThresh = 100
		_,thframe = cv2.threshold(subframe,binThresh,255,cv2.THRESH_BINARY)
			
		# run the tracking algorithm
		new_area = vidTrack(thframe)

		# store the data in an array
		saved_areas = np.append(saved_areas,np.array([new_area]))

		# if the time elapsed is longer than the set interval, save the array
		t02 = time.time()
		if t02-t01 >= save_int:
			with print_lock:
				print('TRACK:: Time elapsed: {}. Saving data.\n'.format(t02-t01))
			outdata = {'areas': saved_areas}
			sio.savemat(savefile,outdata)
			t01 = time.time()
		
		# done processing current frame
		queue.task_done()
		
	# when the while loop exits, save all
	outdata = {'areas': saved_areas}
	sio.savemat(savefile,outdata)
	sig_queue.task_done()
		
# get the relevant information about the animals and the data directory. make
# the folder if necessary.

animal = raw_input('Input animal name (e.g. KH65_66):  ')
print('You chose {} as the animal name'.format(animal))


print('\nOkay, moving on!')

datdir = '/Users/keithhengen/Desktop/animal_video_frames'
direc = os.path.join(datdir,animal)

print('The video data for {} will be saved in directory: {}'.format(animal,direc))
print('I will create that folder if it does not exist.')
try:
	os.stat(direc)
except:
	os.mkdir(direc)

# Initialize arguments for acquisition
#vidcap = cv2.VideoCapture(0) 	# video capture object
frate = 2. 						# framerate
fcount = 0 						# framecounter
verb = True 					# verbose flag (False in function by default)

# Initialize arguments for movement tracking
save_every = .5 # minutes
outfilename = animal + '_savedareas.mat'
outfile = os.path.join(direc, outfilename)

# define tracking process
proc_track = mpc.Process(target=acqFrames,args=[q,q_sig,frate,fcount,verb])

# start processes
proc_track.start()

# start acquisition
runTrackProcess(q,q_sig,save_every*60.,outfile)

# block until done
q.join()
q_sig.join()



