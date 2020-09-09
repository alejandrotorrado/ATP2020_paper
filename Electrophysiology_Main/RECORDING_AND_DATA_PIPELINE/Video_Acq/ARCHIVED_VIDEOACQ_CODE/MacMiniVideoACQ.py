#! /usr/bin/env python

import Tkinter, tkFileDialog
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
    print 'You chose {} as the animal name'.format(animal)
    isOkAnim = input('Is this correct? (1=yes,0=no):  ')
    

print('\nOkay, moving on!')

root = Tkinter.Tk()

isOkDir = 0
while isOkDir == 0:
    datdir = tkFileDialog.askdirectory(parent=root,initialdir="/",
                title='Choose the main folder where to save your video files.')
    print 'The main video data folder will be: {}'.format(datdir)
    isOkDir = input('Is this correct? (1=yes,0=no):  ')

direc = os.path.join(datdir,animal)  
print '\nOkay. Your video data for {} will be saved in: {}'.format(animal,direc)
print 'I will create that folder if it does not exist.'
#direc = os.path.join('/Users/keithhengen/Desktop/animal_video_frames', 'KH67_68')
#animal = 'KH67_68'


try:
    os.stat(direc)
    # if is already a folder, find the highest number file:
    print 'I found previous video files for this animal. \nAdding to the end of the list'
    s = [ ]
    prev = os.listdir(direc)
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

cap     = cv.CaptureFromCAM(0)
height  = 480
width   = 640
fps     = 25

cv.SetCaptureProperty(cap, cv.CV_CAP_PROP_FRAME_HEIGHT, height)
cv.SetCaptureProperty(cap, cv.CV_CAP_PROP_FRAME_WIDTH, width)



# make a function to update the video writer - this includes the filename
def vidsaver(direc,animal,cyclecount):
    print 'Inside vidsaver function'
    path    = os.path.join(direc,animal+'_'+str(cyclecount)+'.mov')
    fourcc  = cv2.cv.CV_FOURCC('a', 'v', 'c', '1')
    writer  = cv2.VideoWriter(path, fourcc,
                         fps, (width, height),True)
    print 'vid saver func executed succesfully'                   
    return writer

# saved video (.avi) is encoded to play back at 25 fps, but the acquisition 
# is set to collect frames at 1/frate 
frate = 2. # this is the number of frames collected per second

# files will turn over/save every nmin minutes
nmin        = 15 # how many minutes you store in each file
ncycles     = nmin*60*frate

# Deal with the Metadata
piDir = '/Users/keithhengen/Desktop/RPI'
checkerDir = os.path.join(piDir,'RECORD')
try:
    os.stat(checkerDir)
except:
    #subprocess.Popen(["mount_smbfs //pi:invivo99@129.64.111.167/ExpLogs ", piDir])
    cmd_string = "mount_smbfs //pi:invivo99@129.64.111.167/ExpLogs " + piDir
    os.system(cmd_string)
    
metafile = os.path.join(checkerDir,animal+'_vid_metadata.csv')


metadat = [animal,time.strftime("%d/%m/%Y"),time.strftime("%H:%M:%S"),
           direc,str(frate),str(nmin)]

with open(metafile, "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            for val in metadat:
                writer.writerow([val]) 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

macminichecker = os.path.join(checkerDir,'macministat.bin')
try:
    os.stat(macminichecker)
    os.remove(macminichecker)
    open(macminichecker, 'w+')
except:
    open(macminichecker, 'w+')

    

acqtime     = [ ]
frames      = 0
cyclecount  = s

while True:

    # update the videowriter
    if frames == 0:
        writer = vidsaver(direc,animal,cyclecount)

    tic1 = time.time()

    acqtime.append(time.time())
       
    
    frame = cv.QueryFrame(cap)
    
    
    frame = np.asarray(frame[:,:])
    tic2 = time.time()
    writer.write(frame)
    toc2 = time.time()
    cv2.imshow('ratcam', frame)
    print "Frame", frames
    
    frames += 1

    toc1 = time.time()
    print 'that took', toc1-tic1
    print toc2-tic2, ' was writing.'
    
    paws = int( round( 1000* (1./frate - (toc1-tic1) )))
    if paws<0 or paws>1000:
        print 'Timing is crappy.'
        paws = 5

    c = cv.WaitKey(paws)     # wait for keyboard input
        
    print 'Exiting waitkey'
    sys.stdout.flush()

   
    if frames % 10 == 0:
        with open(macminichecker,'a') as f:
            f.write(repr(time.time()) + '\n')
        print 'updated PI time'
   
    if c == 27: # esc key to stop video acquisition immediately
        cv.DestroyWindow('ratcam')
        # write the frametimes to a csv file:
        csvfile = os.path.join(direc,animal+'_'+str(cyclecount)+'.csv')
        with open(csvfile, "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            for val in acqtime:
                writer.writerow([val]) 
        try:
            os.stat(macminichecker)
            os.remove(macminichecker)
        except:
            ()
        exit(0)
        
    # reset etc. after ncycles
    if frames == ncycles:
        print 'About to write CSV file of times'
        # write the frametimes to a csv file:
        csvfile = os.path.join(direc,animal+'_'+str(cyclecount)+'.csv')
        with open(csvfile, "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            for val in acqtime:
                writer.writerow([val]) 
        print 'CSV file of times written successfully'
        # reset values, add to counter  
        frames  = 0
        acqtime = [ ]
        cyclecount +=1
        

  
  






