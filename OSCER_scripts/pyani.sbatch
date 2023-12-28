#!/bin/bash
#
#SBATCH --partition=longjobs
#SBATCH --ntasks=12
#SBATCH --mem=24GB
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=120:00:00
#SBATCH --job-name=pyani_long
#SBATCH --mail-user=cmclimans@ou.edu
#SBATCH --chdir=/home/cmclimans/
#
#################################################




source activate pyani


average_nucleotide_identity.py -i all_genomes/ -o /scratch/cmclimans/pyani_out_long -m ANIb -g 
