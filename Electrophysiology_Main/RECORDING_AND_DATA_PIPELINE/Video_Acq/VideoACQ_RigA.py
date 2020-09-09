#! /usr/bin/env python

# import packages
import cv2
import cv2.cv as cv
import time
import numpy as np
import csv
import os
import sys


# get the relevant information about the animals and the data directory. make
# the folder if necessary.
isOkAnim = 0
while isOkAnim == 0:
	animal = raw_input('Input animal name (e.g. KH65_66):  ')
	print('You chose {} as the animal name'.format(animal))
	isOkAnim = input('Is this correct? (1=yes,0=no):  ')

print('\nOkay, moving on!')

#root = Tkinter.Tk()

# Set up the directory to save the data in
datdir = '/Users/keithhengen/Desktop/animal_video_frames'
direc = os.path.join(datdir,animal)
print('The video data for {} will be saved in directory: {}'.format(animal,direc))
print('I will create that folder if it does not exist.')

try:
	os.stat(direc)
	# if is already a folder, find the highest number file and keep adding from there:
	print('I found previous video files for this animal. \nAdding to the end of the list')
	s = [ ]
	prev = os.listdir(direc)
	if not prev:
		s = 1
	else:	
		for x in prev:
			if animal in x:
				print(x)
				s.append(int(x[x.rfind('_')+1 : x.index('.')]))
		s = max(s)+1
	   
except:
	# if not, create the folder
	print('making a new directory')
	os.mkdir(direc)
	s = 1

# set up the camera and the display window to be used
print('Setting up the camera')
cv2.namedWindow('ratcam', 1)

# start the video capture object
cap = cv2.VideoCapture(0)

# make a function to update the video writer - this includes the filename
def vidsaver(direc,animal,cyclecount,dim):
	print('Inside vidsaver function')
	path = os.path.join(direc,animal+'_'+str(cyclecount)+'.avi')
	fourcc = cv2.cv.CV_FOURCC('h','2','6','3')
	video_writer = cv2.VideoWriter(path,fourcc, 30.0, dim)
	#video_writer = cv2.VideoWriter(path,0x00000021, 30.0, dim)
	print('vid saver func executed succesfully')
	return video_writer
	
# another function to save and reset the frame times array
def resetter(cyclecount,acqtime):
	print('About to write CSV file of times')
		# write the frametimes to a csv file:
	csvfile = os.path.join(direc,animal+'_'+str(cyclecount)+'.csv')
	with open(csvfile, "w") as output:
		csvwriter = csv.writer(output, lineterminator='\n')
		for val in acqtime:
			csvwriter.writerow([val])
			
	print('CSV file of times written successfully')
		# reset values, add to counter  
	frames  = 0
	acqtime = [ ]
	cyclecount +=1
	return(cyclecount,frames)

# saved video (.avi) is encoded to play back at 30 fps, but the acquisition 
# is set to collect frames at 'frate' Hz
frate = 6. # this is the number of frames collected per second

# files will turn over/save every nmin minutes
nmin        = 10 # how many minutes of video you store in each file
ncycles     = nmin*60*frate

# Deal with the Metadata
# this is a bit outdated... we used to save video metadata (start time etc.) on an
# external RaspberryPi. Here i've set it up so it's saved on a local folder
piDir = '/Users/keithhengen/Desktop/RPI'
checkerDir = os.path.join(piDir,'RECORD')
try:
	os.stat(checkerDir)
	print("Found pi folder")
except:
	#subprocess.Popen(["mount_smbfs //pi:invivo99@129.64.111.167/ExpLogs ", piDir])
	#cmd_string = "mount_smbfs //pi:invivo99@129.64.111.167/ExpLogs " + piDir
	cmd_string = "mount_smbfs //pi:invivo99@129.64.111.167/ExperimentLogs " + piDir
	os.system(cmd_string)
	print('Mounted RaspberryPi on Desktop')
	
# metadata file
metafile = os.path.join(checkerDir,animal+'_vid_metadata.csv')

# metadata to be saved
metadat = [animal,time.strftime("%d/%m/%Y"),time.strftime("%H:%M:%S"),
		   direc,str(frate),str(nmin)]

# write the metadata file
with open(metafile, "w") as output:
	csvwriter = csv.writer(output, lineterminator='\n')
	for val in metadat:
		csvwriter.writerow([val]) 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# this is a checker file. This code writes to it every so often, so that an
# external RPi can check whether the number in the file is increasing. If it 
# stops increasing, it means the code has crashed. This relies on having the
# external RPi code up and running.
macminichecker = os.path.join(checkerDir,'macministat.bin')
try:
	os.stat(macminichecker)
	os.remove(macminichecker)
	open(macminichecker, 'w+')
except:
	open(macminichecker, 'w+')

	
# initialize variables
acqtime     = [ ]
frames      = 0
cyclecount  = s
c = 0

# Video Acquisition Loop
while True:

	# get frame time and read frame
	tic1 = time.time()
	acqtime.append(tic1)
	ret, frame = cap.read()
	# ret is a Boolean variable, true if frame was acquired, false if not

	# if no frame was acquired, wait and re-try
	while not ret:
		print('ret false, re-trying in 3 seconds')
		time.sleep(1)
		cv2.destroyAllWindows()
		cap.release()
		cap = cv2.VideoCapture(0)
		ret, frame = cap.read()
	
	# convert frame to numpy array
	frame = np.asarray(frame[:,:])
	
	# set frame dimensions for re-sizing
	dim  = (704,576)
	
	# resize frame
	frame = cv2.resize(frame, dim, interpolation = cv2.INTER_AREA)
	
	
	# update the videowriter
	if frames == 0:
		vidwriter = vidsaver(direc,animal,cyclecount,dim)    
	
	# write the frame to the video file
	tic2 = time.time()
	vidwriter.write(frame)
	toc2 = time.time()
	
	# display the frame
	cv2.imshow('ratcam', frame)
	print("Frame {}".format(frames))
	
	# add to the frame counter
	frames += 1

	toc1 = time.time()
	print('That took {}'.format(toc1-tic1))
	print('{}, was writing.'.format(toc2-tic2))
	
	# if time elapsed between frames is very long, reset
	paws = int( round( 1000* (1./frate - (toc1-tic1) )))
	if paws<=0 or paws>=1000:
		print('Timing is crappy.')
		(cyclecount,frames) = resetter(cyclecount,acqtime)
	else:
		c = cv2.waitKey(paws)     # wait for keyboard input
		
	print('Exiting waitkey')
	sys.stdout.flush()

	# every 10 frames, save current time to checker file so RPi can check it   
	if frames % 10 == 0:
		timearr = np.array(time.time())
		timearr.tofile(macminichecker)
#        with open(macminichecker,'a') as f:
#            f.write(repr(time.time()) + '\n')
		print('updated PI time')
   
	# Set up so Esc key saves current video and stops acquisition   
	if c == 27: # esc key to stop video acquisition immediately
		cv2.destroyAllWindows()
		# write the frametimes to a csv file:
		csvfile = os.path.join(direc,animal+'_'+str(cyclecount)+'.csv')
		with open(csvfile, "w") as output:
			csvwriter = csv.writer(output, lineterminator='\n')
			for val in acqtime:
				csvwriter.writerow([val]) 
		try:
			os.stat(macminichecker)
			os.remove(macminichecker)
		except:
			exit(0)

		break
		
	# reset etc. after ncycles (after n minutes of video have been acquired)
	if frames == ncycles:
		(cyclecount,frames) = resetter(cyclecount,acqtime)
		acqtime = []
		
		

  
  






