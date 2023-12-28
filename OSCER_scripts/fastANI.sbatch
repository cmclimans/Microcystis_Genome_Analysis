#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=20
#SBATCH --mem=24GB
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=fastANI_all
#SBATCH --mail-user=cmclimans@ou.edu
#SBATCH --chdir=/home/cmclimans/
#
#################################################
#hostname


#source activate FastANI

module load FastANI

fastANI --ql queries.txt --rl queries.txt -o fastANI_out/fastANI_all.txt --matrix -t 20
