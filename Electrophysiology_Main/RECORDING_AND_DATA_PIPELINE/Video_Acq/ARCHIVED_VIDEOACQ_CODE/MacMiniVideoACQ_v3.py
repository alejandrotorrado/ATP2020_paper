#! /usr/bin/env python

#import Tkinter, tkFileDialog
import cv2
import cv2.cv as cv
import time
import numpy as np
import csv
import os
import sys
import zipfile, zlib
import shutil



# get the relevant information about the animals and the data directory. make
# the folder if necessary.
isOkAnim = 0
while isOkAnim == 0:
    animal = raw_input('Input animal name (e.g. KH65_66):  ')
    print 'You chose {} as the animal name'.format(animal)
    isOkAnim = input('Is this correct? (1=yes,0=no):  ')
    

print('\nOkay, moving on!')

#root = Tkinter.Tk()

#isOkDir = 0
#while isOkDir == 0:
##    datdir = tkFileDialog.askdirectory(parent=root,initialdir="/",
##                title='Choose the main folder where to save your video files.')
datdir = '/Users/keithhengen/Desktop/animal_video_frames'
direc = os.path.join(datdir,animal)
archivedir = os.path.join(datdir,animal,'FRAME_ARCHIVE')

print 'The video data for {} will be saved in directory: {}'.format(animal,direc)
#isOkDir = input('Is this OK? (1=yes,0=no):  ')

#direc = os.path.join(datdir,animal)  
#print '\nOkay. Your video data for {} will be saved in: {}'.format(animal,direc)
print 'I will create that folder if it does not exist.'
#direc = os.path.join('/Users/keithhengen/Desktop/animal_video_frames', 'KH67_68')
#animal = 'KH67_68'


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
            if 'jpg' in x:
                print(x)
                s.append(int(x[x.find('_')+1 : x.rfind('_')]))
        s = max(s)+1
       
except:
    print 'Making directory {} for data storage.'.format(direc)
    os.mkdir(direc)
    os.mkdir(archivedir)
    s = 1

# set up the camera and the display window to be used
print 'Setting up the camera'
cv.NamedWindow('ratcam', 1)


cap = cv2.VideoCapture(0)


# make a function to update the video writer - this includes the filename
def saveframe(image,direc,animal,framecount,cyclecount):
    print 'Inside saveframe function'
    
    imgpath = os.path.join(direc,animal+'_'+str(cyclecount).zfill(4)+'_'+str(framecount).zfill(4)+'.jpg')
#    fourcc  = cv2.cv.CV_FOURCC('h', '2', '6', '3')
#    
#    video_writer = cv2.VideoWriter(path,fourcc, 30.0, dim)
    t1=time.time()
    cv2.imwrite(imgpath,image,(cv2.IMWRITE_JPEG_QUALITY,50))
    t2=time.time()
    print 'Frame {} saved. That took {} seconds'.format(framecount,t2-t1)                   
    
    
def resetter(direc,animal,cyclecount,acqtime,framecount):
    print 'About to write CSV file of times'
        # write the frametimes to a csv file:
    csvfile = os.path.join(direc,animal+'_'+str(cyclecount).zfill(4)+'.csv')
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in acqtime:
            writer.writerow([val])
            
    print 'CSV file of times written successfully'
    acqtime = [ ]
    framecount = 0
    return (framecount,acqtime)
    
def writezip(direc,animal,cyclecount):
    print 'Adding files to archive number {}'.format(cyclecount)
    archivepath = os.path.join(direc,animal+'_'+str(cyclecount).zfill(4)+'.zip')
    zf = zipfile.ZipFile(archivepath,mode='w')
    imgfiles = os.listdir(direc)
    tic0 = time.time()
    for f in imgfiles:
        if 'jpg' in f:
#            print(f)
            cycle = int(f[f.find('_')+1:f.rfind('_')])
        elif 'csv' in f:
#            print(f)
            cycle = int(f[f.find('_')+1:f.index('.')])
        else:
            cycle = 0
                
        if cycle==cyclecount:
            os.chdir(direc)
#            addpath = os.path.join(direc,f)
            # TODO is there a better way to compress so files take less space?
            zf.write(f,compress_type=zipfile.ZIP_DEFLATED)
            os.remove(f)
            
    zf.close()
    toc0 = time.time()
    print 'That took {} seconds.'.format(toc0-tic0)
    shutil.move(archivepath,archivedir)
    cyclecount += 1
    return cyclecount


# saved video (.avi) is encoded to play back at 25 fps, but the acquisition 
# is set to collect frames at 1/frate 
frate = 4. # this is the number of frames collected per second

# files will turn over/save every nmin minutes
nmin        = 30 # how many minutes you store in each file
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
framecount  = 0
cyclecount  = s
c = 0

while True:

    

    tic1 = time.time()

    acqtime.append(tic1)
       
    
    ret, rawframe = cap.read()
    if ret:
        grayframe = cv2.cvtColor(rawframe,cv2.COLOR_BGR2GRAY)
        frame = cv2.resize(grayframe,(0,0),fx=0.25,fy=0.25)

        # TODO while this does the saving, another thread should do the tracking
        saveframe(frame,direc,animal,framecount,cyclecount)

    
        cv2.imshow('ratcam', frame)
    
        framecount += 1

        toc1 = time.time()
        print 'Acquiring and saving frame took {} seconds.'.format(toc1-tic1)
    
        paws = int( round( 1000* (1./frate - (toc1-tic1) )))
        if paws<=0 or paws>=1000:
            print 'Timing is crappy.'
            (framecount,acqtime) = resetter(direc,animal,cyclecount,acqtime,framecount)
        else:
            c = cv.WaitKey(paws)     # wait for keyboard input
        
        print 'Exiting waitkey'
        sys.stdout.flush()
    else:
        print 'No frame was acquired on framecount {}!'.format(framecount)
        raise Exception('No frame available!')

   
    if framecount % 10 == 0:
        with open(macminichecker,'a') as f:
            f.write(repr(time.time()) + '\n')
        print 'Updated PI time'
   
    if c == 27: # esc key to stop video acquisition immediately
        cv.DestroyWindow('ratcam')
        # write the frametimes to a csv file:
        csvfile = os.path.join(direc,animal+'_'+str(cyclecount).zfill(4)+'.csv')
        with open(csvfile, "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            for val in acqtime:
                writer.writerow([val]) 
        # zip remaining frames
        _ = writezip(direc,animal,cyclecount)
        try:
            os.stat(macminichecker)
            os.remove(macminichecker)
        except:
            ()
        exit(0)
        
    # reset etc. after ncycles
    if framecount == ncycles:
        # make csv file of timestamps and update framecount
        (framecount,acqtime) = resetter(direc,animal,cyclecount,acqtime,framecount)
        # zip files and dump in long-term storage. update cyclecount
        cyclecount = writezip(direc,animal,cyclecount)
        # TODO thread call above

        


        
        

  
  






