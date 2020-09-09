#!/bin/bash

#SBATCH --account=neuro
#SBATCH --partition=neuro-largemem
#SBATCH --qos=medium
#SBATCH --time=72:00:00
#SBATCH --job-name=AT67r
#SBATCH --output=/home/atorrpac/job_outputs/AT67r_out.txt
#SBATCH --error=/home/atorrpac/job_outputs/AT67r_err.txt
#SBATCH --ntasks=1
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=atorrpac@brandeis.edu

# cd into folder containing all your matlab code
cd ~/Utilities/hpc_files

# start matlab
module load share_modules/MATLAB/R2018b

# Argument list
anim="'AT67'"	# animal name (e.g. AT30)
side="'B'"	# which side was animal recorded on (A for ch1-32, B for ch33-64)
nsamp=97	# number of samples in each waveform (default is 97)
hem="'L'"	# deprived hemisphere (L or R)
chanstart=22	# start clustering from this channel (set >1 if you want to skip channels)
chsize=120	# chunk size (in hours). 0 for whole experiment
tstart=96	# chunk start time. should be 0 if chsize=0
tstop=216	# chunk end time. should be 0 if chsize=0

# run matlab code
matlab -nodesktop -nosplash -r "AUTOSORT_HPC_func_PAR($anim,$side,$nsamp,$hem,$chanstart,$chsize,$tstart,$tstop);quit"
