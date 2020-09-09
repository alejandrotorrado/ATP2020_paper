# -*- coding: utf-8 -*-
"""
Created on Sat Aug  6 11:59:21 2016

@author: keithhengen
"""
from subprocess import Popen, PIPE, STDOUT, check_output
import time

proc = Popen(['/Users/keithhengen/anaconda2/bin/python2.7','/Users/keithhengen/Google_Drive/Matlab_scripts_11_01_2012/Recording_raw_data/newvideoacq/MacMiniVideoACQ_AutoRestart.py'], shell=False, stdout=PIPE, stderr=STDOUT)
while True:
    data = proc.stdout.readline() #block / wait
    if data:
        print data
        time.sleep(.1)