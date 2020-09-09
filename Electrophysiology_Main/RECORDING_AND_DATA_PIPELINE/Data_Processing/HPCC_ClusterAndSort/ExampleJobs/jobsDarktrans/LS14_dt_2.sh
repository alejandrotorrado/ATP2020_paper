#!/bin/bash

#SBATCH --account=neuro
#SBATCH --partition=neuro-largemem
#SBATCH --qos=medium
#SBATCH --time=48:00:00
#SBATCH --job-name=LS14_dt_2
#SBATCH --output=/home/atorrpac/job_outputs/LS14_dt2_out.txt
#SBATCH --error=/home/atorrpac/job_outputs/LS14_dt2_err.txt
#SBATCH --ntasks=1
#SBATCH -n 8
#SBATCH --mem=10000
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=atorrpac@brandeis.edu

# cd into folder containing all your matlab code
cd ~/Utilities/hpc_files

# start matlab
module load share_modules/MATLAB/R2018b

# Argument list
anim="'LS14'"	# animal name (e.g. AT30)
side="'A'"	# which side was animal recorded on (A for ch1-32, B for ch33-64)
nsamp=97	# number of samples in each waveform (default is 97)
hem="'L'"	# deprived hemisphere (L or R)
chanstart=1	# start clustering from this channel (set >1 if you want to skip channels)
subset=2000000  # size of subset of spikes to be used for clustering
chsize=4	# chunk size (in hours). 0 for whole experiment
tstart=70	# chunk start time. should be 0 if chsize=0
tstop=74	# chunk end time. should be 0 if chsize=0

# run matlab code
matlab -nodesktop -nosplash -r "AUTOSORT_HPC_func_PAR_NEW_BETA($anim,$side,$nsamp,$hem,$chanstart,$subset,$chsize,$tstart,$tstop);quit"
