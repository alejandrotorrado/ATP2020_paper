#! /Users/keithhengen/anaconda2/envs/py3/bin/python

#import Tkinter, tkFileDialog
#import cv2
import time
import numpy as np
import os
import sys
import cv2
import datetime
import scipy.io as sio
import multiprocessing as mpc
sys.path.insert(0, '/Users/keithhengen/Google_Drive/Matlab_scripts_11_01_2012/RECORDING_AND_DATA_PIPELINE/Video_Acq') # path where PYTHON_VIDEOTRACK_FUNCTIONS is
from VideoTrack_Functions import videoTrack_DARK
from VideoTrack_Functions import videoTrack_LIGHT
from VideoTrack_Functions import smoothChunks

print (cv2.__version__)
# _____ set up multiprocessing parameters
# first find number of physical cores available
hwlist = os.popen('sysctl hw').readlines()[1:20]
phys_str = [s for s in hwlist if "hw.physicalcpu:" in s]
if len(phys_str) == 1:
	phys_str = phys_str[0]
else:
	print('Strange... could not determine number of cores')
nCores = int(list(filter(str.isdigit,phys_str))[0])
print('Found {} cores in this machine.\n'.format(nCores))

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

def acqFrames(queue,sigq,framerate,framecount,save_interval,frame_file,savename,verbose=False):
    with print_lock:
        print('Initializing camera.')
    cap = cv2.VideoCapture(0)
    # initialize video writer file
    vidsavecount = 1
    vidname = frame_file + '_' + str(vidsavecount) + '.avi' 
    vid = cv2.VideoWriter(vidname,cv2.VideoWriter_fourcc('m','p','4','v'),30,(648,486))
    savecount = 0
    # initialize (pre-allocate) frametimes array
    frametimes = np.zeros((save_interval,2))
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
            #save frame to video file
            vid.write(rawframe)
            # show frame
            cv2.imshow('acq frame',rawframe)
            # store frame acquisition time in array
            print('framecount = {}'.format(framecount))
            frametimes[framecount - 1,0] = framecount
            frametimes[framecount - 1,1] = tic1 # this is local time in unix format!
            # if time elapsed since start of acquisition round is longer than save_interval.
            # save the frametimes array and reset it. Also reset savetime counter
            if framecount >= save_interval:
                savecount += 1
                with print_lock:
                    print('ACQ: saving frametimes (savecounter: {}).\n'.format(savecount))
                ftdata = {'frametimes': frametimes}
                sio.savemat(savename + '_' + str(savecount) + '.mat',ftdata)
                # reset frametimes array and video writer file
                frametimes = np.zeros((save_interval,2))
                framecount = 0
                vid.release()
                vidsavecount += 1
                vidname = frame_file + '_' + str(vidsavecount) + '.avi'
                vid = cv2.VideoWriter(vidname,cv2.VideoWriter_fourcc('m','p','4','v'),30,(648,486))
                
                
            # for accurate timing, measure time required to take frame
            toc1 = time.time()
            # pause for long enough to get right framerate
            paws = int( round( 1000* (1./framerate - (toc1-tic1) )))
            if verbose:
                with print_lock:
                    print('ACQ:: Acquiring frame {} took {} seconds.\n'.format(framecount,toc1-tic1))
                    print('paws was {}'.format(paws))
            # wait for button press
            if paws > 0:
                k = cv2.waitKey(paws)
            else:
                k = cv2.waitKey(1)
            # pressing Esc key causes infinite loop to break
            if k == 27:
                savecount += 1
                with print_lock:
                    print('ACQ:: Halting acquisition on framecount {}. Saving frametimes one last time (savecounter: {}).\n'.format(framecount,savecount))
                cv2.destroyAllWindows()
                cap.release()
                ftdata = {'frametimes': frametimes[0:framecount,:]}
                sio.savemat(savename + '_' + str(savecount) + '.mat',ftdata)
                vid.release()
                sigq.put(0)
                break
        else:
            with print_lock:
                print('ACQ:: No frame was acquired on framecount {}!\n'.format(framecount))
            continue


def runTrackProcess(queue,sig_queue,save_int,mask_1,mask_2,frame_file,savefile_1,savefile_2,initParams_1,initParams_2,counter):
    with print_lock:
        print('I am process ', mpc.current_process().name , ' and I just started.')
	
    # initializiation parameters
    t01 = time.time()
    savecounter = 0
    track_framecounter = 0
    runningAvg = None
    # define kernels for morphological operations
    erode_kern = np.ones((2,2),np.uint8)
    close_kern = np.ones((2,2),np.uint8)

    ###### Initialize parameters for light/dark
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

    raw_mvt_1 = np.zeros((save_int,1))
    track_1 = np.zeros((save_int,2))
    raw_mvt_2 = np.zeros((save_int,1))
    track_2 = np.zeros((save_int,2))
    

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
        # update total and internal frame counters
        counter += 1
        track_framecounter += 1
        # figure out if frame is light or dark
        intensFrac = np.mean(colorframe[mask_1==0.])/np.mean(colorframe[mask_1==1.])
        if intensFrac > 0.7:
            frameIsDark = False
            params = paramsLight
        else:
            frameIsDark = True
            params = paramsDark

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
        _,thframe = cv2.threshold(subframe,binThresh,255,cv2.THRESH_BINARY)
        
        # apply mask for animal 1
        thframe_1 = thframe.copy()
        thframe_1[mask_1==0.] = 0
        frame_1 = cv2.medianBlur(thframe_1,5)

        # apply mask for animal 2
        thframe_2 = thframe.copy()
        thframe_2[mask_2==0.] = 0
        frame_2 = cv2.medianBlur(thframe_2,5)

        # apply morphological closing
        nErosions = params['nErosions']
        nDilations = params['nDilations']
        # animal 1 frame
        frame_1 = cv2.erode(frame_1,erode_kern,iterations = nErosions)
        frame_1 = cv2.dilate(frame_1,close_kern,iterations = nDilations)
        # animal 2 frame
        frame_2 = cv2.erode(frame_2,erode_kern,iterations = nErosions)
        frame_2 = cv2.dilate(frame_2,close_kern,iterations = nDilations)

        # run the tracking algorithm
        if frameIsDark:
            # animal 1
            raw_mvt_1[track_framecounter-1,0],track_2[track_framecounter-1,:],new_mc_1,old_mc_1,firstIsGood_1,newalpha = videoTrack_DARK(frame_1,params,initParams_1,raw_mvt_1,track_1,counter,verbose=False)
            # animal 2
            raw_mvt_2[track_framecounter-1,0],track_2[track_framecounter-1,:],new_mc_2,old_mc_2,firstIsGood_2,newalpha = videoTrack_DARK(frame_2,params,initParams_2,raw_mvt_2,track_2,counter,verbose=False)
            params['alpha'] = newalpha
            initParams_1['firstIsGood'] = firstIsGood_1
            initParams_2['firstIsGood'] = firstIsGood_2
        else:
            # animal 1
            raw_mvt_1[track_framecounter-1,0],track_1[track_framecounter-1,:],new_mc_1,old_mc_1 = videoTrack_LIGHT(frame_1,params,initParams_1,raw_mvt_1,track_1,counter,verbose=False)
            # animal 2
            raw_mvt_2[track_framecounter-1,0],track_2[track_framecounter-1,:],new_mc_2,old_mc_2 = videoTrack_LIGHT(frame_2,params,initParams_2,raw_mvt_2,track_2,counter,verbose=False)
            
        # update CoM and other params
        initParams_1['old_mcx'] = old_mc_1[0]
        initParams_1['old_mcy'] = old_mc_1[1]
        initParams_1['mcx'] = new_mc_1[0]
        initParams_1['mcy'] = new_mc_1[1]
        initParams_2['old_mcx'] = old_mc_2[0]
        initParams_2['old_mcy'] = old_mc_2[1]
        initParams_2['mcx'] = new_mc_2[0]
        initParams_2['mcy'] = new_mc_2[1]
        
        # if the time elapsed is longer than the set interval, save the array
        t02 = time.time()
#        if t02-t01 >= save_int:
        if track_framecounter >= save_int:
            savecounter += 1
            with print_lock:
                print('TRACK:: Time elapsed: {}. Saving data (savecounter: {}).\n'.format(t02-t01,savecounter))
            outdata_1 = {'raw_mvt': raw_mvt_1, 'track': track_1}
            outdata_2 = {'raw_mvt': raw_mvt_2, 'track': track_2}
            sio.savemat(savefile_1 + '_' + str(savecounter) + '.mat',outdata_1)
            sio.savemat(savefile_2 + '_' + str(savecounter) + '.mat',outdata_2)
            # reset counter and arrays
            track_framecounter = 0
            
            raw_mvt_1[:,:] = 0
            track_1[:,:] = 0
            raw_mvt_2[:,:] = 0
            track_2[:,:] = 0
        
        # done processing current frame
        queue.task_done()
		
    # when the while loop exits, save all
    savecounter += 1
    with print_lock:
        print('TRACK:: quitting everything! Saving data one last time (savecounter: {}).\n'.format(savecounter))
    outdata_1 = {'raw_mvt': raw_mvt_1[0:track_framecounter,0], 'track': track_1[0:track_framecounter,:]}
    outdata_2 = {'raw_mvt': raw_mvt_2[0:track_framecounter,0], 'track': track_2[0:track_framecounter,:]}
    sio.savemat(savefile_1 + '_' + str(savecounter) + '.mat',outdata_1)
    sio.savemat(savefile_2 + '_' + str(savecounter) + '.mat',outdata_2)
    sig_queue.task_done()


# function to select ROI for tracking
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
    try:
        cv2.rectangle(firstFrame, refPt[0], refPt[1], (0, 255, 0), 2)
        cv2.putText(firstFrame, "Press ESC to continue.", (int(.35*firstFrame.shape[1]), firstFrame.shape[0] - 30),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255,255,255), 2)
        cv2.imshow('firstframe',firstFrame)
    except:
        None
            
        
            

######################################################################
############################# MAIN BODY ##############################
######################################################################

##### Initialize arguments for acquisition
# *** NOTE: this has to go here because of Tkinter bug. Call to any process.start()
# must come BEFORE import Tkinter statement
#vidcap = cv2.VideoCapture(0) 	# video capture object
frate = 2. 						# framerate
fcount = 0 						# framecounter
verb = True 					# verbose flag (False in function by default)


# get the relevant information about the animals and the data directory. make
# the folder if necessary.

animal = input('Input animal name (format: KH65_66):  ')
print('You chose {} as the animal name'.format(animal))
animal_base = animal[0:2]
animal_numbers = animal[2:]
animal_1 = animal_base + animal_numbers.split('_')[0]
animal_2 = animal_base + animal_numbers.split('_')[1]


print('\nOkay, moving on!')

datdir = '/Users/keithhengen/Desktop/animal_video_frames' 
direc = os.path.join(datdir,animal)
direc_1 = os.path.join(direc,animal_1)
direc_2 = os.path.join(direc,animal_2)

print('The video data for {} will be saved in directory: {}'.format(animal_1,direc_1))
print('The video data for {} will be saved in directory: {}'.format(animal_2,direc_2))
print('I will create those folder if they do not exist.')
# create folders if do not exist
try:
    os.stat(direc)
except:
    os.mkdir(direc)
try:
    os.stat(direc_1)
except:
    os.mkdir(direc_1)
try:
    os.stat(direc_2)
except:
    os.mkdir(direc_2)


# Initialize arguments for movement tracking
save_every = 1200 # frames - divide by frate to get minutes
framefilename = animal + '_rawframes'
framefile = os.path.join(direc, framefilename)
outfilename_1 = animal_1 + '_pymovement'
outfile_1 = os.path.join(direc_1, outfilename_1)
outfilename_2 = animal_2 + '_pymovement'
outfile_2 = os.path.join(direc_2, outfilename_2)
ft_savefile = animal + '_frametimes'
ft_save = os.path.join(direc, ft_savefile)



########### Initialize general parameters
initParams_1 = {}   # dictionary for other general parameters, including those initialized on frame 1
initParams_2 = {}   # as above, but for second animal
refPt = []          # container for ROI limits

initParams_1['mcx'] = 0
initParams_1['mcy'] = 0
initParams_1['firstIsGood'] = True
initParams_2['mcx'] = 0
initParams_2['mcy'] = 0
initParams_2['firstIsGood'] = True



########## define ROIs for each animal
vidSrc = cv2.VideoCapture(0)
ret,firstFrame = vidSrc.read()
firstFrame = cv2.resize(firstFrame,(0,0),fx=0.25,fy=0.25)
height = firstFrame.shape[0]
width = firstFrame.shape[1]
cv2.namedWindow('firstframe')
cv2.setMouseCallback('firstframe',selectROI)

# find roi for animal 1
while True:
    cv2.putText(firstFrame, "Animal 1: Select ROI corners (top left to bottom right)", (int(.05*width), int(.08*height)),
                cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255,255,255), 2)
    cv2.imshow('firstframe',firstFrame)
    key = cv2.waitKey(1500)
    if key == 27:
        cv2.destroyWindow('firstframe')
        break
# if user defined an ROI, make a mask based on it
if len(refPt) == 2:
    print('Making mask for animal 1.')
    roi_1 = firstFrame[refPt[0][1]:refPt[1][1], refPt[0][0]:refPt[1][0]]
    mask_1 = np.zeros(firstFrame.shape[0:2])
    mask_1[refPt[0][1]:refPt[1][1], refPt[0][0]:refPt[1][0]] = 1.
    initParams_1['mask_half_area'] = mask_1[mask_1==1.].shape[0]/2

# clear refPt and firstFrame and go again for animal 2
firstFrame = []
refPt = []
ret,firstFrame = vidSrc.read()
firstFrame = cv2.resize(firstFrame,(0,0),fx=0.25,fy=0.25)
cv2.namedWindow('firstframe')
cv2.setMouseCallback('firstframe',selectROI)
while True:
    cv2.putText(firstFrame, "Animal 2: Select ROI corners (top left to bottom right)", (int(.05*width), int(.08*height)),
                cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255,255,255), 2)
    cv2.imshow('firstframe',firstFrame)
    key = cv2.waitKey(0)
    if key == 27:
        cv2.destroyWindow('firstframe')
        break
# if user defined an ROI, make a mask based on it
if len(refPt) == 2:
    print('Making mask for animal 2.')
    roi_2 = firstFrame[refPt[0][1]:refPt[1][1], refPt[0][0]:refPt[1][0]]
    mask_2 = np.zeros(firstFrame.shape[0:2])
    mask_2[refPt[0][1]:refPt[1][1], refPt[0][0]:refPt[1][0]] = 1.
    initParams_2['mask_half_area'] = mask_2[mask_2==1.].shape[0]/2

vidSrc.release()

# define acquisition process
proc_track = mpc.Process(target=runTrackProcess,args=[q,q_sig,save_every,mask_1,mask_2,framefile,outfile_1,outfile_2,initParams_1,initParams_2,fcount])
#proc_acq = mpc.Process(target=acqFrames,args=[q,q_sig,frate,fcount,save_every*60.,ft_save,verb])

# start acquisition processes
proc_track.start()
#proc_acq.start()

# start tracking
#runTrackProcess(q,q_sig,save_every*60.,mask_1,mask_2,outfile_1,outfile_2,fcount)
acqFrames(q,q_sig,frate,fcount,save_every,framefile,ft_save,verb)
# block until done
q.join()
q_sig.join()



