#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=8
#SBATCH --mem=8GB
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=isescan
#SBATCH --mail-user=cmclimans@ou.edu
#SBATCH --chdir=/home/cmclimans/
#
#################################################

source activate isescan

while read line; do isescan.py --seqfile $line --output /scratch/cmclimans/isescan_out --nthread 8; done < isescan_files.txt

