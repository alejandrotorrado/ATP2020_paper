#! /usr/bin/env python

import cv2
import cv2.cv as cv
import time
import numpy as np
import csv
import os
import sys
import glob



# get the relevant information about the animals and the data directory. make
# the folder if necessary.
piDir = '/Users/keithhengen/Desktop/RPI'
checkerDir = os.path.join(piDir,'RECORD')
try:
    os.stat(checkerDir)
except:
    #subprocess.Popen(["mount_smbfs //pi:invivo99@129.64.111.167/ExpLogs ", piDir])
    cmd_string = "mount_smbfs //pi:invivo99@129.64.111.167/ExpLogs " + piDir
    os.system(cmd_string)
    
    
print('\nRetrieving metadata...\n')

metafilelist = glob.glob(os.path.join(checkerDir,'*_vid_metadata.csv'))
metafile = metafilelist[0]
metadata = []
with open(metafile, "r") as f:
    csvreader = csv.reader(f, lineterminator='\n')
    for row in csvreader:
        metadata.append(row[0])
        
animal = metadata[0]
startdate = metadata[1]
starttime = metadata[2]
direc = metadata[3]
frate = float(metadata[4])
nmin = float(metadata[5])
ncycles     = nmin*60*frate



print 'The video data for {} will be saved in directory: {}'.format(animal,direc)

try:
    os.stat(direc)
    # if is already a folder, find the highest number file:
    print 'I found previous video files for this animal. \nAdding to the end of the list'
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
    print 'making a new directory'
    os.mkdir(direc)
    s = 1

# set up the camera and the display window to be used
print 'Setting up the camera'
cv.NamedWindow('ratcam', 1)


cap = cv2.VideoCapture(0)


# make a function to update the video writer - this includes the filename
def vidsaver(direc,animal,cyclecount,dim):
    print 'Inside vidsaver function'
    
    path    = os.path.join(direc,animal+'_'+str(cyclecount)+'.avi')
    fourcc  = cv2.cv.CV_FOURCC('h', '2', '6', '3')
    
    video_writer = cv2.VideoWriter(path,fourcc, 30.0, dim)
    
    print 'vid saver func executed succesfully'                   
    return video_writer
    
    
def resetter(cyclecount,acqtime):
    print 'About to write CSV file of times'
        # write the frametimes to a csv file:
    csvfile = os.path.join(direc,animal+'_'+str(cyclecount)+'.csv')
    with open(csvfile, "w") as output:
        csvwriter = csv.writer(output, lineterminator='\n')
        for val in acqtime:
            csvwriter.writerow([val])
            
    print 'CSV file of times written successfully'
        # reset values, add to counter  
    frames  = 0
    acqtime = [ ]
    cyclecount +=1
    return(cyclecount,frames)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

macminichecker = os.path.join(checkerDir,'macministat.bin')
try:
    os.stat(macminichecker)
    os.remove(macminichecker)
    open(macminichecker, 'wb')
except:
    open(macminichecker, 'wb')

# create video restart flag
open(os.path.join(checkerDir,'vid_restarted.bin'),'wb')
timearr = np.array(time.time())
timearr.tofile(os.path.join(checkerDir,'vid_restarted.bin'))

acqtime     = [ ]
frames      = 0
cyclecount  = s
c = 0

while True:

    

    tic1 = time.time()

    acqtime.append(tic1)
       
    
    ret, frame = cap.read()
    frame = np.asarray(frame[:,:])
    r = 640.0 / frame.shape[1]    
    dim = (640, int(frame.shape[0]*r))
    frame = cv2.resize(frame, dim, interpolation = cv2.INTER_AREA)
    
    # update the videowriter
    if frames == 0:
        vidwriter = vidsaver(direc,animal,cyclecount,dim)    
    
    tic2 = time.time()
    vidwriter.write(frame)
    toc2 = time.time()
    
    cv2.imshow('ratcam', frame)
    print "Frame", frames
    
    frames += 1

    toc1 = time.time()
    print 'that took', toc1-tic1
    print toc2-tic2, ' was writing.'
    
    paws = int( round( 1000* (1./frate - (toc1-tic1) )))
    if paws<=0 or paws>=1000:
        print 'Timing is crappy.'
        (cyclecount,frames) = resetter(cyclecount,acqtime)
    else:
        c = cv.WaitKey(paws)     # wait for keyboard input
        
    print 'Exiting waitkey'
    sys.stdout.flush()

   
    if frames % 10 == 0:
        timearr = np.array(time.time())
        timearr.tofile(macminichecker)
#        with open(macminichecker,'a') as f:
#            f.write(str(time.time())+'\n')
        print 'updated PI time'
   
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
        
    # reset etc. after ncycles
    if frames == ncycles:
        (cyclecount,frames) = resetter(cyclecount,acqtime)
        acqtime = []
        
        

  
  






