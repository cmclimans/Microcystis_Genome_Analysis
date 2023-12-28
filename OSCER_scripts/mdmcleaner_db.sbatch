#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=10GB
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=mdmcleaner_db
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cmclimans@ou.edu
#SBATCH --chdir=/scratch/cmclimans/pangenome/raw_fna
#
#################################################

#if anaconda 
eval "$(conda shell.bash hook)"

source activate mdmcleaner

mdmcleaner makedb -o mdmcleaner_db_node/

