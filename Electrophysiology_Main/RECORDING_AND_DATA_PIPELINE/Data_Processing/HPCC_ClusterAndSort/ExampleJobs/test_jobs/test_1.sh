#!/bin/bash

#SBATCH --account=neuro
#SBATCH --partition=neuro-largemem
#SBATCH --time=00:10:00
#SBATCH --job-name=test_1
#SBATCH --output=/home/atorrpac/test_1.txt
#SBATCH --error=err_1.txt
#SBATCH --ntasks=1
#SBATCH -n 7
#SBATCH --mem=12000
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=atorrpac@brandeis.edu


cd ~/Utilities/hpc_files
module load share_modules/MATLAB/R2018b

p=49

matlab -nodesktop -nosplash -r "test_func2($p);quit"
