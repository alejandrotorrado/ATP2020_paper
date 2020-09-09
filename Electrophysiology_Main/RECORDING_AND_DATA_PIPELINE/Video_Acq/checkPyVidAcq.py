#! /usr/bin/env python
"""
Created on Sat Aug  6 11:35:39 2016

@author: keithhengen
"""

import os, sys
import csv
import numpy as np
import time
from subprocess import Popen, PIPE, STDOUT, check_output
sys.path.insert(0, '/Users/keithhengen/Desktop/RPI')
import LogManager as LM
from send_text import send_all_texts
import glob

# figure out the process ID for the python process running spyder, to avoid
# terminating it when trying to restart video acquisition
pidout0 = check_output("ps aux | grep 'python' | awk 'NR==1 {print $2}'", shell=True)
pid0 = pidout0.rstrip('\n')
print('Spyder process ID: {}.'.format(pid0))

checkPi = False

def vidloopchecker(vidstat,vid,mmlast,pi,pi_base,rpitm):
    
    # check that the rpi is online by establishing that the rpi_is_up.txt file's modification time
    # is continuously increasing on each loop. this will be used to determine the potential causes
    # of a crash below - it could be that the mac mini has stopped writing data, or that the pi
    # code is crashed. this is IMPORTANT.
    
    print 'Starting buffer counter in vidloopchecker'
    sys.stdout.flush()
    vCount = 0
    vbuffwipe = 1000
    while vidstat:
        stopped = os.path.isfile(os.path.join(pi,'STOPREC.bin'))
        if stopped:
            print 'Found STOPREC flag! Recording has stopped. Quitting python.'
            sys.stdout.flush()
            # write to video log file
            vlogStr = 'Found STOPREC flag in vidloopchecker. Quitting python.'
            LM.addToLog(logs,vidlog,vlogStr,time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))
            sys.exit()
        
        if checkPi:
            rpitimetemp = os.path.getmtime(os.path.join(pi_base,'rpi_is_up.txt'));
            if rpitimetemp > rpitm:
                rpitm = rpitimetemp
                time.sleep(30)
            else:
                time.sleep(60)
            
                rpitimetemp = os.path.getmtime(os.path.join(pi_base,'rpi_is_up.txt'));
                if not rpitimetemp > rpitm:
                    tp = 1; tpcount = 0;
                    while tp:
                        tpcount = tpcount+1;
                        myMsg = 'from iDell. RPI is dead. RIP. Make a move, playa.'
                        send_all_texts(myMsg)
                        time.sleep(15)
                        if tpcount == 20:
                            tp = 0
                else:
                    rpitm = rpitimetemp

        print 'I am in vidloopchecker.'
        sys.stdout.flush()
        print 'Checking macmini matlab code...'
        sys.stdout.flush()
        
        
        mmtime_all = np.fromfile(vid, count=-1)
        while not mmtime_all:
            print('Trying to read macministat file...')
            time.sleep(2)
            mmtime_all = np.fromfile(vid, count=-1)
        try:
            mmtime_now = mmtime_all[-1]
        except:                
            time.sleep(5)
            mmtime_now = mmtime_all[-1]
        if mmtime_now > mmlast:
            print 'MacMini is OK at {}'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))

            vCount += 1
            mmlast = mmtime_now
            if vCount > vbuffwipe:
                # when there are more than vbuffwipe entries in macministat
                # wipe and leave only last one
                last_vidstamp = np.fromfile(vid)[-1]
                wiped = open(vid, 'wb')
                a = np.array(last_vidstamp,'float64')
                a.tofile(wiped)
                wiped.close()
                if len(np.fromfile(vid)) == 1:
                    print 'Successfully wiped macministat after {} writes.'.format(vbuffwipe)
                    sys.stdout.flush()
                else:
                    print 'Problems wiping macministat buffer...'
                    sys.stdout.flush()
                vCount = 0 # reset counter
	    time.sleep(30)
        else:
            # if RPi is not online, may have crashed. Halt checking if that's the case
            if not os.path.getmtime(os.path.join(pi_base,'rpi_is_up.txt')):
                print 'Looks like RPi might be dead.'
                vlogStr = 'Macmini checker found potential RPi death. Pausing vidloopchecker.'
                LM.addToLog(logs,vidlog,vlogStr,time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))
                myMsg = 'Macmini python code found potential RPi death. I will stop checking video until this is resolved.'
                send_all_texts(myMsg)
                rpidead_vcount = 0
                while not os.path.isfile(os.path.join(pi_base,'pi_exist.txt')):
                    print 'Waiting for RPi to come back online.'
                    rpidead_vcount += 1
                    myMsg = 'RPi is still down after ~2 minutes.'
                    send_all_texts(myMsg)
                    rpidead_vcount = 0
            else:
                print 'MacMini crashed. Trying to restart.'
                sys.stdout.flush()
                # create crash flag in RPi
                open(os.path.join(pi,'macminiCrash.bin'),'w')
                # write to video log file
                vlogStr = 'Found macmini crash in vidloopchecker. Create macminiCrash flag.'
                LM.addToLog(logs,vidlog,vlogStr,time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))
                # call autorestart code                
                vid_startover(vid,mmlast,pi)
         
def vid_startover(vidstat,vidtime_last,pi):
    stopped = os.path.isfile(os.path.join(pi,'STOPREC.bin'))
    if stopped:
        print 'Found STOPREC flag! Recording has stopped. Quitting python.'
        sys.stdout.flush()
        ## write to video log file
        vlogStr = 'Found STOPREC flag in vid_startover. Quitting python.'
        LM.addToLog(logs,vidlog,vlogStr,time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))
        sys.exit()
    
    # write to video log file
    vlogStr = 'Executing vid_startover. Killing matlab.'
    LM.addToLog(logs,vidlog,vlogStr,time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))
    
    #### HAVE TO CHANGE THIS PART
    print 'I\'m killing matlab on macmini at {}'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))
    try:
        # kill matlab
        pidout1 = check_output("ps aux | grep 'python' | awk 'NR==1 {print $2}'", shell=True)
        pid1 = pidout1.rstrip('\n')
        if pid1 != pid0:
            os.system('kill -9 ' + pid1)
    except:
        print 'Could not execute kill statements in vid_startover.'
    
    
    # make sure there is no vid restarted flag:
    if os.path.isfile(os.path.join(pi,'vid_restarted.bin')):        
        try:
            print 'Attempting to remove old video restart flag from RPi.'
            sys.stdout.flush()
            os.remove(os.path.join(pi,'vid_restarted.bin'))
        except:
            print 'couldn\'t delete vid restart flag.'
    
    print 'Restarting video acquisition now at {}'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))
    # write to video log file
    vlogStr = 'Executing vid_startover. Restarting matlab.'
    LM.addToLog(logs,vidlog,vlogStr,time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))
    
    output_filename = '/Users/keithhengen/Desktop/AutoVidACQ.txt'
    output_file = open(output_filename, 'w')
    Popen(['/Users/keithhengen/anaconda2/bin/python',\
    '/Users/keithhengen/Google_Drive/Matlab_scripts_11_01_2012/Recording_raw_data/newvideoacq/MacMiniVideoACQ_AutoRestart.py'],\
    shell=False, stdout=output_file, stderr=output_file)
    
    print 'Video acquisition on MacMini restarted.'
    while not os.path.isfile(os.path.join(pi,'vid_restarted.bin')):
        print 'Waiting for video restart flag...'
        sys.stdout.flush()
        time.sleep(1)
        
    time.sleep(10)
    newVt = np.fromfile(os.path.join(pi,'vid_restarted.bin'), count = -1)[-1]
    print 'Found video restart flag on RPi. Moving on.'
    sys.stdout.flush()
    # write to video log file
    vlogStr = 'Found video restart flag in vid_startover. Video acquisition restarted.'
    LM.addToLog(logs,vidlog,vlogStr,time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(newVt)))
    
    
    if newVt > vidtime_last:
        time.sleep(10)
        # write to video log file
        vlogStr = 'Resuming vidloopchecker.'
        LM.addToLog(logs,vidlog,vlogStr,time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))
        vidloopchecker(True,vidstat,vidtime_last,pi,pi_base,0)
    else:
        try:
            raise Exception('bad restart flag')
            print 'problems'
            sys.stdout.flush()
            # send text warning
        except Exception as error:
            print(repr(error))
    
    
logs = '/Users/keithhengen/Desktop/RPI/LOGS'
pi_base = '/Users/keithhengen/Desktop/RPI'
pi = '/Users/keithhengen/Desktop/RPI/RECORD'

# initialize the variable for comparing and storing the time stamp on rpi_is_up.txt, which is used
# to confirm that the code on the pi is still running. 
rpitm = 0;

# look for video metadata file
filelist=glob.glob(os.path.join(pi,'*vid_metadata.csv'))
while not os.path.isfile(filelist[0]):
    print 'Waiting for video metadata file...'
    sys.stdout.flush()
    time.sleep(1)
print 'Found video metadata file!'
sys.stdout.flush()

# retrieve video metadata
x=[]
with open(filelist[0]) as f:
    reader = csv.reader(f,lineterminator='\n')
    for row in reader:
        x.append(row)
animal = x[0][0]
print 'Starting video log file for animal: {}'.format(animal)
sys.stdout.flush()

# start video log for this animal
vidlog = 'VideoLog_' + animal + '.txt'
startVlog = 'Started ' + animal + ' video recording.'
vlogStart = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
if not os.path.isfile(os.path.join(logs,vidlog)):
    LM.createLog(logs,vidlog)
else:
    print 'Found videolog file: {}'.format(os.path.join(logs,vidlog))
    sys.stdout.flush()
LM.addToLog(logs,vidlog,startVlog,vlogStart)


vidgo = os.path.isfile(os.path.join(pi,'macministat.bin'))

# write to video log file
vlogStr = 'Starting MacMini python code to monitor video acquisition.'
LM.addToLog(logs,vidlog,vlogStr,time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))

print vidgo
while not vidgo:
    print 'No MacMini status file yet. Waiting.'
    sys.stdout.flush()
    time.sleep(2)
    vidgo = os.path.isfile(os.path.join(pi,'macministat.bin'))

print 'Found MacMini status file. Entering vidloop checker.'
sys.stdout.flush()

# write to video log file
vlogStr = 'Found macministat flag in macmini python code. Entering vid_loopchecker.'
LM.addToLog(logs,vidlog,vlogStr,time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))

vidstatfile = os.path.join(pi,'macministat.bin')
vidtime_last = 0

vidloopchecker(vidgo,vidstatfile,vidtime_last,pi,pi_base,rpitm)